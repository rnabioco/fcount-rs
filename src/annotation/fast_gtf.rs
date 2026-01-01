//! Fast GTF parser optimized for feature counting.
//!
//! This parser is designed for speed over flexibility:
//! - Minimal allocations (reuses buffers)
//! - Simple field splitting (no full GFF validation)
//! - Only extracts fields we need

use anyhow::{Context, Result};
use memchr::memchr;
use rustc_hash::FxHashMap;
use std::io::BufRead;
use std::sync::Arc;

use super::feature::{Feature, Strand};
use super::index::AnnotationIndex;
use super::io::open_reader;
use crate::cli::Args;

/// Parse a GTF file quickly, extracting only what we need
pub fn load_gtf_fast(args: &Args) -> Result<AnnotationIndex> {
    let reader = open_reader(&args.annotation)?;

    let mut chrom_to_id: FxHashMap<Arc<str>, u16> = FxHashMap::default();
    let mut id_to_chrom: Vec<Arc<str>> = Vec::new();
    let mut gene_name_to_id: FxHashMap<Arc<str>, u32> = FxHashMap::default();
    let mut gene_names: Vec<Arc<str>> = Vec::new();
    let mut features: Vec<Feature> = Vec::with_capacity(500_000);

    let feature_type = args.feature_type.as_bytes();
    let gene_id_attr = "gene_id \"".to_string();
    let gene_id_attr_bytes = gene_id_attr.as_bytes();

    for line in reader.lines() {
        let line = line.context("Failed to read line")?;

        // Skip comments and empty lines
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let bytes = line.as_bytes();

        // Parse tab-separated fields: chrom, source, feature, start, end, score, strand, frame, attributes
        // We need: chrom (0), feature (2), start (3), end (4), strand (6), attributes (8)

        let mut field_starts = [0usize; 9];
        let mut field_idx = 0;
        let mut pos = 0;

        while field_idx < 8 && pos < bytes.len() {
            if let Some(tab_pos) = memchr(b'\t', &bytes[pos..]) {
                field_starts[field_idx + 1] = pos + tab_pos + 1;
                pos = pos + tab_pos + 1;
                field_idx += 1;
            } else {
                break;
            }
        }

        if field_idx < 8 {
            continue; // Malformed line
        }

        // Check feature type (field 2)
        let feat_start = field_starts[2];
        let feat_end = field_starts[3] - 1;
        if &bytes[feat_start..feat_end] != feature_type {
            continue;
        }

        // Parse chromosome (field 0)
        let chrom_end = field_starts[1] - 1;
        let chrom_str = std::str::from_utf8(&bytes[0..chrom_end]).unwrap_or("");
        let chrom_name: Arc<str> = Arc::from(chrom_str);
        let chrom_id = *chrom_to_id.entry(chrom_name.clone()).or_insert_with(|| {
            let id = id_to_chrom.len() as u16;
            id_to_chrom.push(chrom_name);
            id
        });

        // Parse start (field 3) and end (field 4)
        let start_str =
            std::str::from_utf8(&bytes[field_starts[3]..field_starts[4] - 1]).unwrap_or("0");
        let end_str =
            std::str::from_utf8(&bytes[field_starts[4]..field_starts[5] - 1]).unwrap_or("0");
        let start: u32 = start_str.parse().unwrap_or(0);
        let end: u32 = end_str.parse().unwrap_or(0);

        // Parse strand (field 6)
        let strand_byte = bytes.get(field_starts[6]).copied().unwrap_or(b'.');
        let strand = match strand_byte {
            b'+' => Strand::Forward,
            b'-' => Strand::Reverse,
            _ => Strand::Unknown,
        };

        // Parse gene_id from attributes (field 8)
        let attrs_start = field_starts[8];
        let attrs = &bytes[attrs_start..];

        // Find gene_id "value" in attributes
        let gene_id = if let Some(pos) = find_subsequence(attrs, gene_id_attr_bytes) {
            let value_start = pos + gene_id_attr_bytes.len();
            // Find closing quote
            if let Some(quote_pos) = memchr(b'"', &attrs[value_start..]) {
                std::str::from_utf8(&attrs[value_start..value_start + quote_pos]).ok()
            } else {
                None
            }
        } else {
            None
        };

        let gene_id = match gene_id {
            Some(id) => id,
            None => continue, // Skip if no gene_id
        };

        let gene_name: Arc<str> = Arc::from(gene_id);
        let gene_id = *gene_name_to_id.entry(gene_name.clone()).or_insert_with(|| {
            let id = gene_names.len() as u32;
            gene_names.push(gene_name);
            id
        });

        features.push(Feature {
            gene_id,
            chrom_id,
            start,
            end,
            strand,
        });
    }

    // Sort features by (chrom_id, start)
    features.sort_by_key(|f| (f.chrom_id, f.start));

    AnnotationIndex::new(chrom_to_id, id_to_chrom, gene_names, features)
}

/// Find a subsequence in a byte slice
#[inline]
fn find_subsequence(haystack: &[u8], needle: &[u8]) -> Option<usize> {
    if needle.is_empty() {
        return Some(0);
    }
    haystack
        .windows(needle.len())
        .position(|window| window == needle)
}
