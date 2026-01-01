//! Fast GTF parser optimized for feature counting.
//!
//! This parser is designed for speed over flexibility:
//! - Minimal allocations (reuses buffers)
//! - Simple field splitting (no full GFF validation)
//! - Only extracts fields we need
//! - Auto-detection of feature types

use anyhow::{Context, Result};
use memchr::{memchr, memmem};
use rustc_hash::{FxHashMap, FxHashSet};
use std::io::BufRead;
use std::path::Path;
use std::sync::Arc;

use super::feature::{Feature, Strand};
use super::index::AnnotationIndex;
use super::io::open_reader;
use crate::cli::Args;

/// Scan a GTF/GFF file to detect available feature types.
/// Reads up to `max_lines` non-comment lines for efficiency.
pub fn detect_feature_types(path: &Path, max_lines: usize) -> Result<Vec<String>> {
    let mut reader = open_reader(path)?;
    let mut feature_types: FxHashSet<String> = FxHashSet::default();
    let mut line_buf = String::with_capacity(1024);
    let mut lines_read = 0;

    while reader
        .read_line(&mut line_buf)
        .context("Failed to read line")?
        > 0
    {
        let line = line_buf.trim_end();

        // Skip comments and empty lines
        if line.is_empty() || line.starts_with('#') {
            line_buf.clear();
            continue;
        }

        let bytes = line.as_bytes();

        // Find field 2 (feature type) by locating tabs
        // Fields: chrom(0), source(1), feature(2), ...
        let mut tab_count = 0;
        let mut field_start = 0;

        for (i, &b) in bytes.iter().enumerate() {
            if b == b'\t' {
                tab_count += 1;
                if tab_count == 2 {
                    field_start = i + 1;
                } else if tab_count == 3 {
                    // Extract feature type
                    if let Ok(feat_type) = std::str::from_utf8(&bytes[field_start..i]) {
                        feature_types.insert(feat_type.to_string());
                    }
                    break;
                }
            }
        }

        lines_read += 1;
        if lines_read >= max_lines {
            break;
        }

        line_buf.clear();
    }

    // Return sorted for consistent output
    let mut types: Vec<String> = feature_types.into_iter().collect();
    types.sort();
    Ok(types)
}

/// Auto-select the best feature type based on detected types.
/// Priority: exon > exonic_part > CDS > gene
pub fn auto_select_feature_type(available: &[String]) -> Option<&str> {
    const PREFERRED: &[&str] = &["exon", "exonic_part", "CDS", "gene"];

    for &preferred in PREFERRED {
        if available.iter().any(|t| t == preferred) {
            return Some(preferred);
        }
    }
    available.first().map(|s| s.as_str())
}

/// Parse a GTF file quickly, extracting only what we need.
/// If `feature_type_override` is provided, it is used instead of `args.feature_type`.
pub fn load_gtf_fast(args: &Args, feature_type_override: Option<&str>) -> Result<AnnotationIndex> {
    let mut reader = open_reader(&args.annotation)?;

    let mut chrom_to_id: FxHashMap<Arc<str>, u16> = FxHashMap::default();
    let mut id_to_chrom: Vec<Arc<str>> = Vec::with_capacity(64); // ~chromosomes
    let mut gene_name_to_id: FxHashMap<Arc<str>, u32> = FxHashMap::default();
    let mut gene_names: Vec<Arc<str>> = Vec::with_capacity(30_000); // ~genes
    let mut features: Vec<Feature> = Vec::with_capacity(500_000);

    let feature_type_str = feature_type_override.unwrap_or(&args.feature_type);
    let feature_type = feature_type_str.as_bytes();
    // GTF format: gene_id "value"
    let gtf_attr = format!("{} \"", args.gene_id_attr);
    let gtf_attr_bytes = gtf_attr.as_bytes();
    // GFF3 format: gene_id=value
    let gff3_attr = format!("{}=", args.gene_id_attr);
    let gff3_attr_bytes = gff3_attr.as_bytes();

    // Reuse line buffer to avoid allocations per line
    let mut line_buf = String::with_capacity(1024);
    while reader
        .read_line(&mut line_buf)
        .context("Failed to read line")?
        > 0
    {
        // Trim trailing newline
        let line = line_buf.trim_end();

        // Skip comments and empty lines
        if line.is_empty() || line.starts_with('#') {
            line_buf.clear();
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
            line_buf.clear();
            continue; // Malformed line
        }

        // Check feature type (field 2)
        let feat_start = field_starts[2];
        let feat_end = field_starts[3] - 1;
        if &bytes[feat_start..feat_end] != feature_type {
            line_buf.clear();
            continue;
        }

        // Parse chromosome (field 0)
        let chrom_end = field_starts[1] - 1;
        let chrom_str = std::str::from_utf8(&bytes[0..chrom_end]).unwrap_or("");
        // Look up by &str first (avoids Arc allocation on hit)
        let chrom_id = if let Some(&id) = chrom_to_id.get(chrom_str) {
            id
        } else {
            let chrom_name: Arc<str> = Arc::from(chrom_str);
            let id = id_to_chrom.len() as u16;
            id_to_chrom.push(chrom_name.clone());
            chrom_to_id.insert(chrom_name, id);
            id
        };

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

        // Try GTF format first: gene_id "value"
        let gene_id = if let Some(pos) = find_subsequence(attrs, gtf_attr_bytes) {
            let value_start = pos + gtf_attr_bytes.len();
            // Find closing quote
            if let Some(quote_pos) = memchr(b'"', &attrs[value_start..]) {
                std::str::from_utf8(&attrs[value_start..value_start + quote_pos]).ok()
            } else {
                None
            }
        } else if let Some(pos) = find_subsequence(attrs, gff3_attr_bytes) {
            // Try GFF3 format: gene_id=value (terminated by ; or end of line)
            let value_start = pos + gff3_attr_bytes.len();
            let remaining = &attrs[value_start..];
            let value_end = memchr(b';', remaining).unwrap_or(remaining.len());
            std::str::from_utf8(&remaining[..value_end]).ok()
        } else {
            None
        };

        let gene_id = match gene_id {
            Some(id) => id,
            None => {
                line_buf.clear();
                continue; // Skip if no gene_id
            }
        };

        // Look up by &str first (avoids Arc allocation on hit)
        let gene_idx = if let Some(&id) = gene_name_to_id.get(gene_id) {
            id
        } else {
            let gene_name: Arc<str> = Arc::from(gene_id);
            let id = gene_names.len() as u32;
            gene_names.push(gene_name.clone());
            gene_name_to_id.insert(gene_name, id);
            id
        };

        features.push(Feature {
            start,
            end,
            gene_id: gene_idx,
            chrom_id,
            strand,
        });

        // Clear buffer for next line (keeps allocation)
        line_buf.clear();
    }

    // Sort features by (chrom_id, start) - unstable is faster and order doesn't matter for equal keys
    features.sort_unstable_by_key(|f| (f.chrom_id, f.start));

    AnnotationIndex::new(chrom_to_id, id_to_chrom, gene_names, features)
}

/// Find a subsequence in a byte slice using SIMD-accelerated search
#[inline]
fn find_subsequence(haystack: &[u8], needle: &[u8]) -> Option<usize> {
    memmem::find(haystack, needle)
}
