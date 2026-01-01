//! I/O utilities for annotation files.
//!
//! Shared utilities for reading annotation files with automatic gzip detection.

use anyhow::{Context, Result};
use flate2::read::GzDecoder;
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom};
use std::path::Path;

/// Buffer size for reading annotation files (1 MB)
const BUFFER_SIZE: usize = 1024 * 1024;

/// Check if a file is gzip-compressed by reading magic bytes.
///
/// Reads the first two bytes of the file and checks for the gzip magic number (0x1f 0x8b).
/// The file position is reset to the beginning after checking.
pub fn is_gzipped(file: &mut File) -> std::io::Result<bool> {
    let mut magic = [0u8; 2];
    let n = file.read(&mut magic)?;
    file.seek(SeekFrom::Start(0))?;
    Ok(n == 2 && magic == [0x1f, 0x8b])
}

/// Create a buffered reader, auto-detecting gzip compression.
///
/// Returns a boxed `BufRead` trait object that handles both plain text and gzipped files.
/// Uses a 1 MB buffer for efficient reading.
pub fn open_reader(path: &Path) -> Result<Box<dyn BufRead>> {
    let mut file =
        File::open(path).with_context(|| format!("Failed to open: {}", path.display()))?;

    if is_gzipped(&mut file)? {
        Ok(Box::new(BufReader::with_capacity(
            BUFFER_SIZE,
            GzDecoder::new(file),
        )))
    } else {
        Ok(Box::new(BufReader::with_capacity(BUFFER_SIZE, file)))
    }
}
