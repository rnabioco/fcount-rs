use anyhow::{Context, Result};
use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_sam as sam;
use noodles_sam::alignment::record::Cigar;
use smallvec::SmallVec;
use std::fs::File;
use std::num::NonZeroUsize;
use std::path::Path;

use super::cigar::Interval;
use crate::annotation::AnnotationIndex;

/// SAM flags for read filtering
#[allow(dead_code)]
pub mod flags {
    pub const UNMAPPED: u16 = 0x4;
    pub const MATE_UNMAPPED: u16 = 0x8;
    pub const REVERSE_STRAND: u16 = 0x10;
    pub const MATE_REVERSE_STRAND: u16 = 0x20;
    pub const FIRST_IN_PAIR: u16 = 0x40;
    pub const SECOND_IN_PAIR: u16 = 0x80;
    pub const SECONDARY: u16 = 0x100;
    pub const NOT_PASSING_FILTERS: u16 = 0x200;
    pub const DUPLICATE: u16 = 0x400;
    pub const SUPPLEMENTARY: u16 = 0x800;
}

/// Processed alignment record with parsed fields
#[derive(Debug, Clone)]
pub struct AlignmentRecord {
    /// Read name (for paired-end mate matching)
    pub read_name: Vec<u8>,
    /// SAM flags
    pub flags: u16,
    /// Chromosome ID (index into annotation)
    pub chrom_id: Option<u16>,
    /// 1-based start position
    pub start: u32,
    /// Mapping quality
    pub mapq: u8,
    /// Parsed CIGAR intervals
    pub intervals: SmallVec<[Interval; 4]>,
    /// NH tag value (number of alignments, 1 if not present)
    pub nh: u8,
    /// Mate chromosome ID (for paired-end)
    pub mate_chrom_id: Option<u16>,
    /// Mate position (for paired-end)
    pub mate_start: u32,
    /// Template length (for paired-end)
    pub template_len: i32,
}

impl AlignmentRecord {
    #[inline]
    pub fn is_unmapped(&self) -> bool {
        self.flags & flags::UNMAPPED != 0
    }

    #[inline]
    pub fn is_mate_unmapped(&self) -> bool {
        self.flags & flags::MATE_UNMAPPED != 0
    }

    #[inline]
    pub fn is_reverse_strand(&self) -> bool {
        self.flags & flags::REVERSE_STRAND != 0
    }

    #[inline]
    pub fn is_mate_reverse_strand(&self) -> bool {
        self.flags & flags::MATE_REVERSE_STRAND != 0
    }

    #[inline]
    pub fn is_first_in_pair(&self) -> bool {
        self.flags & flags::FIRST_IN_PAIR != 0
    }

    #[inline]
    pub fn is_second_in_pair(&self) -> bool {
        self.flags & flags::SECOND_IN_PAIR != 0
    }

    #[inline]
    pub fn is_secondary(&self) -> bool {
        self.flags & flags::SECONDARY != 0
    }

    #[inline]
    pub fn is_duplicate(&self) -> bool {
        self.flags & flags::DUPLICATE != 0
    }

    #[inline]
    pub fn is_supplementary(&self) -> bool {
        self.flags & flags::SUPPLEMENTARY != 0
    }

    #[inline]
    pub fn is_paired(&self) -> bool {
        self.flags & 0x1 != 0
    }

    #[inline]
    pub fn is_proper_pair(&self) -> bool {
        self.flags & 0x2 != 0
    }
}

/// Alignment file reader that handles both BAM and SAM
pub struct AlignmentReader {
    inner: ReaderInner,
    header: sam::Header,
    /// Pre-computed mapping from BAM ref_id → annotation chrom_id
    /// This avoids HashMap lookups and string conversions per record
    ref_to_chrom: Vec<Option<u16>>,
}

enum ReaderInner {
    Bam(bam::io::Reader<bgzf::MultithreadedReader<File>>),
    Sam(sam::io::Reader<std::io::BufReader<File>>),
}

impl AlignmentReader {
    /// Open a BAM or SAM file with multi-threaded BGZF decompression
    /// Pre-computes ref_id → chrom_id mapping for fast lookups
    pub fn open_with_annotation(
        path: &Path,
        threads: usize,
        annotation: &AnnotationIndex,
    ) -> Result<Self> {
        // Check if it's a BAM file by extension
        let is_bam = path
            .extension()
            .map(|e| e.to_string_lossy().to_lowercase() == "bam")
            .unwrap_or(false);

        if is_bam {
            let file = File::open(path)
                .with_context(|| format!("Failed to open BAM file: {}", path.display()))?;

            // Use multi-threaded BGZF reader for parallel decompression
            let worker_count =
                NonZeroUsize::new(threads.max(1)).expect("thread count must be positive");
            let bgzf_reader = bgzf::MultithreadedReader::with_worker_count(worker_count, file);
            let mut bam_reader = bam::io::Reader::from(bgzf_reader);
            let header = bam_reader.read_header()?;

            // Pre-compute ref_id → chrom_id mapping
            let ref_to_chrom = build_ref_to_chrom_mapping(&header, annotation);

            Ok(AlignmentReader {
                inner: ReaderInner::Bam(bam_reader),
                header,
                ref_to_chrom,
            })
        } else {
            let file = File::open(path)
                .with_context(|| format!("Failed to open SAM file: {}", path.display()))?;
            let mut sam_reader = sam::io::Reader::new(std::io::BufReader::new(file));
            let header = sam_reader.read_header()?;

            // Pre-compute ref_id → chrom_id mapping
            let ref_to_chrom = build_ref_to_chrom_mapping(&header, annotation);

            Ok(AlignmentReader {
                inner: ReaderInner::Sam(sam_reader),
                header,
                ref_to_chrom,
            })
        }
    }

    /// Legacy open without annotation (for backwards compatibility)
    pub fn open(path: &Path, threads: usize) -> Result<Self> {
        // Check if it's a BAM file by extension
        let is_bam = path
            .extension()
            .map(|e| e.to_string_lossy().to_lowercase() == "bam")
            .unwrap_or(false);

        if is_bam {
            let file = File::open(path)
                .with_context(|| format!("Failed to open BAM file: {}", path.display()))?;

            let worker_count =
                NonZeroUsize::new(threads.max(1)).expect("thread count must be positive");
            let bgzf_reader = bgzf::MultithreadedReader::with_worker_count(worker_count, file);
            let mut bam_reader = bam::io::Reader::from(bgzf_reader);
            let header = bam_reader.read_header()?;

            Ok(AlignmentReader {
                inner: ReaderInner::Bam(bam_reader),
                header,
                ref_to_chrom: Vec::new(),
            })
        } else {
            let file = File::open(path)
                .with_context(|| format!("Failed to open SAM file: {}", path.display()))?;
            let mut sam_reader = sam::io::Reader::new(std::io::BufReader::new(file));
            let header = sam_reader.read_header()?;

            Ok(AlignmentReader {
                inner: ReaderInner::Sam(sam_reader),
                header,
                ref_to_chrom: Vec::new(),
            })
        }
    }

    /// Get the header
    pub fn header(&self) -> &sam::Header {
        &self.header
    }

    /// Iterate over records with fast chrom_id lookup
    pub fn records(&mut self) -> RecordIter<'_> {
        RecordIter {
            reader: self,
            record_buf: sam::alignment::RecordBuf::default(),
        }
    }
}

/// Build ref_id → chrom_id mapping from BAM header
fn build_ref_to_chrom_mapping(
    header: &sam::Header,
    annotation: &AnnotationIndex,
) -> Vec<Option<u16>> {
    header
        .reference_sequences()
        .iter()
        .map(|(name, _)| {
            std::str::from_utf8(name.as_ref())
                .ok()
                .and_then(|name_str| annotation.get_chrom_id(name_str))
        })
        .collect()
}

pub struct RecordIter<'a> {
    reader: &'a mut AlignmentReader,
    record_buf: sam::alignment::RecordBuf,
}

impl<'a> RecordIter<'a> {
    /// Read next record into the provided AlignmentRecord
    /// Returns Ok(true) if a record was read, Ok(false) at EOF
    /// Uses pre-computed ref_id → chrom_id mapping for fast lookup
    pub fn read_record(
        &mut self,
        _annotation: &AnnotationIndex, // Kept for API compatibility
        out: &mut AlignmentRecord,
    ) -> Result<bool> {
        let bytes_read = match &mut self.reader.inner {
            ReaderInner::Bam(reader) => {
                reader.read_record_buf(&self.reader.header, &mut self.record_buf)?
            }
            ReaderInner::Sam(reader) => {
                reader.read_record_buf(&self.reader.header, &mut self.record_buf)?
            }
        };

        if bytes_read == 0 {
            return Ok(false);
        }

        // Extract fields from record_buf (reuse buffer to avoid allocation)
        out.read_name.clear();
        if let Some(name) = self.record_buf.name() {
            out.read_name.extend_from_slice(name.as_ref());
        }
        out.flags = self.record_buf.flags().bits();

        // Use pre-computed ref_id → chrom_id mapping (fast array lookup instead of HashMap)
        out.chrom_id = self
            .record_buf
            .reference_sequence_id()
            .and_then(|id| self.reader.ref_to_chrom.get(id).copied().flatten());

        // Position (convert to 1-based)
        out.start = self
            .record_buf
            .alignment_start()
            .map(|p| p.get() as u32)
            .unwrap_or(0);

        // Mapping quality
        out.mapq = self
            .record_buf
            .mapping_quality()
            .map(|q| q.get())
            .unwrap_or(255);

        // Parse CIGAR
        out.intervals.clear();
        super::cigar::parse_cigar_intervals(
            self.record_buf.cigar().iter(),
            out.start,
            &mut out.intervals,
        );

        // NH tag
        out.nh = self
            .record_buf
            .data()
            .get(&sam::alignment::record::data::field::Tag::ALIGNMENT_HIT_COUNT)
            .and_then(|v| {
                use sam::alignment::record_buf::data::field::Value;
                match v {
                    Value::UInt8(n) => Some(*n),
                    Value::Int8(n) => Some(*n as u8),
                    Value::UInt16(n) => Some(*n as u8),
                    Value::Int16(n) => Some(*n as u8),
                    Value::UInt32(n) => Some(*n as u8),
                    Value::Int32(n) => Some(*n as u8),
                    _ => None,
                }
            })
            .unwrap_or(1);

        // Mate info (also uses pre-computed mapping)
        out.mate_chrom_id = self
            .record_buf
            .mate_reference_sequence_id()
            .and_then(|id| self.reader.ref_to_chrom.get(id).copied().flatten());

        out.mate_start = self
            .record_buf
            .mate_alignment_start()
            .map(|p| p.get() as u32)
            .unwrap_or(0);

        out.template_len = self.record_buf.template_length();

        Ok(true)
    }
}
