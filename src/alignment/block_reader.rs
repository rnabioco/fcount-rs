//! BGZF block-level reader for parallel BAM processing.
//!
//! This module provides low-level access to decompressed BAM data chunks
//! for distribution to worker threads.

use anyhow::{Context, Result};
use noodles_bgzf as bgzf;
use noodles_sam as sam;
use std::fs::File;
use std::io::Read;
use std::num::NonZeroUsize;
use std::path::Path;

use super::minimal_parser::{get_record_size, parse_bam_record, MinimalRecord};
use crate::annotation::AnnotationIndex;

/// Size of each batch of decompressed BAM data (target ~4MB per batch)
/// Larger batches reduce channel overhead and improve cache efficiency
const BATCH_TARGET_SIZE: usize = 4 * 1024 * 1024;

/// A batch of decompressed BAM data ready for parallel parsing
#[derive(Debug)]
pub struct RecordBatch {
    /// Decompressed BAM record data (concatenated records)
    pub data: Vec<u8>,
    /// Any leftover partial record from previous batch
    pub leftover: Vec<u8>,
}

impl RecordBatch {
    pub fn new() -> Self {
        RecordBatch {
            data: Vec::with_capacity(BATCH_TARGET_SIZE + 65536),
            leftover: Vec::new(),
        }
    }

    pub fn clear(&mut self) {
        self.data.clear();
        self.leftover.clear();
    }
}

impl Default for RecordBatch {
    fn default() -> Self {
        Self::new()
    }
}

/// BAM block reader that produces batches for parallel processing
pub struct BamBlockReader {
    /// Underlying BGZF reader (multi-threaded for parallel decompression)
    reader: bgzf::MultithreadedReader<File>,
    /// BAM header
    header: sam::Header,
    /// Pre-computed ref_id -> chrom_id mapping
    ref_to_chrom: Vec<Option<u16>>,
    /// Current read buffer
    read_buf: Vec<u8>,
    /// Leftover data from previous batch (partial record)
    leftover: Vec<u8>,
    /// Whether we've reached EOF
    eof: bool,
}

impl BamBlockReader {
    /// Open a BAM file and read the header
    /// Uses multi-threaded BGZF decompression for better performance
    pub fn open(path: &Path, annotation: &AnnotationIndex) -> Result<Self> {
        Self::open_with_threads(path, annotation, 4) // Default to 4 decompression threads
    }

    /// Open a BAM file with specified number of decompression threads
    pub fn open_with_threads(
        path: &Path,
        annotation: &AnnotationIndex,
        threads: usize,
    ) -> Result<Self> {
        let file = File::open(path)
            .with_context(|| format!("Failed to open BAM file: {}", path.display()))?;

        // Use multi-threaded BGZF reader for parallel decompression
        let worker_count =
            NonZeroUsize::new(threads.max(1)).expect("thread count must be positive");
        let mut reader = bgzf::MultithreadedReader::with_worker_count(worker_count, file);

        // Read BAM magic number
        let mut magic = [0u8; 4];
        reader
            .read_exact(&mut magic)
            .context("Failed to read BAM magic")?;
        if &magic != b"BAM\x01" {
            anyhow::bail!("Invalid BAM magic number");
        }

        // Read header length and header text
        let mut header_len_buf = [0u8; 4];
        reader
            .read_exact(&mut header_len_buf)
            .context("Failed to read header length")?;
        let header_len = u32::from_le_bytes(header_len_buf) as usize;

        let mut header_text = vec![0u8; header_len];
        reader
            .read_exact(&mut header_text)
            .context("Failed to read header text")?;

        // Read reference sequences
        let mut n_ref_buf = [0u8; 4];
        reader
            .read_exact(&mut n_ref_buf)
            .context("Failed to read reference count")?;
        let n_ref = u32::from_le_bytes(n_ref_buf) as usize;

        let mut ref_names = Vec::with_capacity(n_ref);
        for _ in 0..n_ref {
            // Read name length
            let mut name_len_buf = [0u8; 4];
            reader.read_exact(&mut name_len_buf)?;
            let name_len = u32::from_le_bytes(name_len_buf) as usize;

            // Read name (including null terminator)
            let mut name = vec![0u8; name_len];
            reader.read_exact(&mut name)?;

            // Convert to string (excluding null terminator)
            let name_str = String::from_utf8_lossy(&name[..name_len.saturating_sub(1)]).to_string();
            ref_names.push(name_str);

            // Read and discard sequence length
            let mut _seq_len_buf = [0u8; 4];
            reader.read_exact(&mut _seq_len_buf)?;
        }

        // Build ref_id -> chrom_id mapping
        let ref_to_chrom: Vec<Option<u16>> = ref_names
            .iter()
            .map(|name| annotation.get_chrom_id(name))
            .collect();

        // Parse header for noodles compatibility (not strictly needed but useful)
        let header_str = String::from_utf8_lossy(&header_text);
        let header = header_str.parse::<sam::Header>().unwrap_or_default();

        Ok(BamBlockReader {
            reader,
            header,
            ref_to_chrom,
            read_buf: vec![0u8; 65536],
            leftover: Vec::new(),
            eof: false,
        })
    }

    /// Get the SAM header
    pub fn header(&self) -> &sam::Header {
        &self.header
    }

    /// Get the ref_id to chrom_id mapping
    pub fn ref_to_chrom(&self) -> &[Option<u16>] {
        &self.ref_to_chrom
    }

    /// Read the next batch of BAM records
    ///
    /// Returns None at EOF
    pub fn read_batch(&mut self) -> Result<Option<RecordBatch>> {
        if self.eof {
            return Ok(None);
        }

        let mut batch = RecordBatch::new();

        // Start with any leftover data from previous batch
        if !self.leftover.is_empty() {
            batch.data.extend_from_slice(&self.leftover);
            self.leftover.clear();
        }

        // Read until we have enough data or hit EOF
        while batch.data.len() < BATCH_TARGET_SIZE {
            match self.reader.read(&mut self.read_buf) {
                Ok(0) => {
                    self.eof = true;
                    break;
                }
                Ok(n) => {
                    batch.data.extend_from_slice(&self.read_buf[..n]);
                }
                Err(e) => return Err(e.into()),
            }
        }

        if batch.data.is_empty() {
            return Ok(None);
        }

        // Find the last complete record boundary
        // Records are: 4 bytes block_size + block_size bytes data
        let mut offset = 0;
        let mut last_complete_offset = 0;

        while offset + 4 <= batch.data.len() {
            let record_size = get_record_size(&batch.data[offset..]);
            if record_size == 0 {
                break;
            }

            let total_size = 4 + record_size;
            if offset + total_size > batch.data.len() {
                // This record is incomplete
                break;
            }

            last_complete_offset = offset + total_size;
            offset = last_complete_offset;
        }

        // Save incomplete data as leftover for next batch
        if last_complete_offset < batch.data.len() {
            self.leftover
                .extend_from_slice(&batch.data[last_complete_offset..]);
            batch.data.truncate(last_complete_offset);
        }

        if batch.data.is_empty() && !self.leftover.is_empty() {
            // We only have partial data, need to read more
            // This shouldn't happen often with 1MB batches
            return self.read_batch();
        }

        Ok(Some(batch))
    }

    /// Parse all records from a batch
    ///
    /// This is called by worker threads
    pub fn parse_batch_records<F>(
        batch: &RecordBatch,
        need_read_name: bool,
        record: &mut MinimalRecord,
        mut callback: F,
    ) where
        F: FnMut(&MinimalRecord),
    {
        let data = &batch.data;
        let mut offset = 0;

        while offset + 4 <= data.len() {
            let record_size = get_record_size(&data[offset..]);
            if record_size == 0 {
                break;
            }

            let data_start = offset + 4;
            let data_end = data_start + record_size;

            if data_end > data.len() {
                break;
            }

            // Parse this record
            if parse_bam_record(&data[data_start..data_end], record, need_read_name).is_ok() {
                callback(record);
            }

            offset = data_end;
        }
    }
}

/// Iterator over records in a batch (alternative API)
pub struct BatchRecordIter<'a> {
    data: &'a [u8],
    offset: usize,
    record: MinimalRecord,
    need_read_name: bool,
}

impl<'a> BatchRecordIter<'a> {
    pub fn new(batch: &'a RecordBatch, need_read_name: bool) -> Self {
        BatchRecordIter {
            data: &batch.data,
            offset: 0,
            record: MinimalRecord::new(),
            need_read_name,
        }
    }
}

impl<'a> Iterator for BatchRecordIter<'a> {
    type Item = MinimalRecord;

    fn next(&mut self) -> Option<Self::Item> {
        if self.offset + 4 > self.data.len() {
            return None;
        }

        let record_size = get_record_size(&self.data[self.offset..]);
        if record_size == 0 {
            return None;
        }

        let data_start = self.offset + 4;
        let data_end = data_start + record_size;

        if data_end > self.data.len() {
            return None;
        }

        // Parse this record
        if parse_bam_record(
            &self.data[data_start..data_end],
            &mut self.record,
            self.need_read_name,
        )
        .is_ok()
        {
            self.offset = data_end;
            Some(self.record.clone())
        } else {
            self.offset = data_end;
            self.next() // Skip invalid record
        }
    }
}
