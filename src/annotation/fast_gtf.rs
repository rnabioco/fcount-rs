//! Fast GTF parser optimized for feature counting.
//!
//! This parser is designed for speed over flexibility:
//! - Minimal allocations (reuses buffers)
//! - Simple field splitting (no full GFF validation)
//! - Only extracts fields we need
//! - Auto-detection of feature types

use anyhow::{Context, Result};
use flate2::read::GzDecoder;
use memchr::{memchr, memmem};
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};
use std::fs::File;
use std::io::{BufRead, Read};
use std::path::Path;
use std::sync::Arc;

use super::feature::{Feature, Strand};
use super::index::AnnotationIndex;
use super::io::{is_gzipped, open_reader};
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

/// Parse an unsigned decimal integer directly from ASCII bytes.
///
/// Skips the UTF-8 validation + str::parse pipeline. Stops at the first
/// non-digit byte; returns 0 if no leading digit. Trailing whitespace and
/// other junk are ignored — adequate for GTF integer fields.
#[inline]
fn parse_u32_ascii(bytes: &[u8]) -> u32 {
    let mut n: u32 = 0;
    for &b in bytes {
        let d = b.wrapping_sub(b'0');
        if d >= 10 {
            break;
        }
        n = n.wrapping_mul(10).wrapping_add(d as u32);
    }
    n
}

/// A feature with chunk-local chrom/gene IDs, produced by `parse_chunk`.
struct LocalFeature {
    chrom_id_local: u16,
    gene_id_local: u32,
    start: u32,
    end: u32,
    strand: Strand,
}

/// Per-chunk parse output. Name slices borrow from the file buffer, so this
/// struct's lifetime is tied to it.
struct ChunkResult<'a> {
    chrom_names: Vec<&'a [u8]>,
    gene_names: Vec<&'a [u8]>,
    features: Vec<LocalFeature>,
}

/// Parse a contiguous chunk of GTF/GFF text. The chunk must start at a line
/// boundary; trailing partial lines are skipped if the chunk doesn't end at
/// `\n`. Returns features with chunk-local IDs for chrom/gene names.
fn parse_chunk<'a>(
    data: &'a [u8],
    feature_type: &[u8],
    gtf_attr_bytes: &[u8],
    gff3_attr_bytes: &[u8],
) -> ChunkResult<'a> {
    let mut chrom_to_id_local: FxHashMap<&'a [u8], u16> = FxHashMap::default();
    let mut chrom_names: Vec<&'a [u8]> = Vec::with_capacity(64);
    let mut gene_to_id_local: FxHashMap<&'a [u8], u32> = FxHashMap::default();
    let mut gene_names: Vec<&'a [u8]> = Vec::with_capacity(30_000);
    let mut features: Vec<LocalFeature> = Vec::with_capacity(data.len() / 200);

    let mut start_pos = 0usize;
    while start_pos < data.len() {
        let nl = memchr(b'\n', &data[start_pos..]).map(|n| start_pos + n);
        let line_end = nl.unwrap_or(data.len());

        // Trim trailing \r
        let line_end_trimmed = if line_end > start_pos && data[line_end - 1] == b'\r' {
            line_end - 1
        } else {
            line_end
        };

        let bytes = &data[start_pos..line_end_trimmed];

        // Advance start_pos before any continue
        start_pos = line_end + 1; // past the \n; ok if line_end == data.len()

        if bytes.is_empty() || bytes[0] == b'#' {
            continue;
        }

        // Locate tab-separated field boundaries
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
            continue;
        }

        // Cheapest filter first: feature type
        if &bytes[field_starts[2]..field_starts[3] - 1] != feature_type {
            continue;
        }

        let chrom_bytes = &bytes[0..field_starts[1] - 1];
        let chrom_id_local = if let Some(&id) = chrom_to_id_local.get(chrom_bytes) {
            id
        } else {
            let id = chrom_names.len() as u16;
            chrom_names.push(chrom_bytes);
            chrom_to_id_local.insert(chrom_bytes, id);
            id
        };

        let start = parse_u32_ascii(&bytes[field_starts[3]..field_starts[4] - 1]);
        let end = parse_u32_ascii(&bytes[field_starts[4]..field_starts[5] - 1]);

        let strand = match bytes.get(field_starts[6]).copied().unwrap_or(b'.') {
            b'+' => Strand::Forward,
            b'-' => Strand::Reverse,
            _ => Strand::Unknown,
        };

        let attrs = &bytes[field_starts[8]..];
        let gene_id_bytes: Option<&'a [u8]> =
            if let Some(p) = find_subsequence(attrs, gtf_attr_bytes) {
                let s = p + gtf_attr_bytes.len();
                memchr(b'"', &attrs[s..]).map(|q| &attrs[s..s + q])
            } else if let Some(p) = find_subsequence(attrs, gff3_attr_bytes) {
                let s = p + gff3_attr_bytes.len();
                let rem = &attrs[s..];
                let e = memchr(b';', rem).unwrap_or(rem.len());
                Some(&rem[..e])
            } else {
                None
            };

        let gene_id_bytes = match gene_id_bytes {
            Some(b) => b,
            None => continue,
        };

        let gene_id_local = if let Some(&id) = gene_to_id_local.get(gene_id_bytes) {
            id
        } else {
            let id = gene_names.len() as u32;
            gene_names.push(gene_id_bytes);
            gene_to_id_local.insert(gene_id_bytes, id);
            id
        };

        features.push(LocalFeature {
            chrom_id_local,
            gene_id_local,
            start,
            end,
            strand,
        });
    }

    ChunkResult {
        chrom_names,
        gene_names,
        features,
    }
}

/// Read the whole annotation file into a single `Vec<u8>`, transparently
/// decompressing gzip. We need a contiguous buffer so the parallel parser
/// can split it across threads.
fn read_annotation_to_vec(path: &Path) -> Result<Vec<u8>> {
    let mut file =
        File::open(path).with_context(|| format!("Failed to open: {}", path.display()))?;
    if is_gzipped(&mut file)? {
        let mut out = Vec::with_capacity(64 * 1024 * 1024);
        GzDecoder::new(file).read_to_end(&mut out)?;
        Ok(out)
    } else {
        let len = file.metadata().map(|m| m.len() as usize).unwrap_or(0);
        let mut out = Vec::with_capacity(len);
        file.read_to_end(&mut out)?;
        Ok(out)
    }
}

/// Find chunk boundaries aligned to newlines. Each returned `(start, end)`
/// range spans complete lines (start is either 0 or just after a `\n`; end is
/// either at `data.len()` or just past a `\n`).
fn line_aligned_chunks(data: &[u8], num_chunks: usize) -> Vec<(usize, usize)> {
    let n = data.len();
    if num_chunks <= 1 || n == 0 {
        return vec![(0, n)];
    }
    let target = n / num_chunks;
    let mut bounds = Vec::with_capacity(num_chunks + 1);
    bounds.push(0);
    for i in 1..num_chunks {
        let pivot = (i * target).min(n);
        // Walk forward to the next newline so each chunk starts on a line.
        let off = pivot
            + memchr(b'\n', &data[pivot..])
                .map(|o| o + 1)
                .unwrap_or(n - pivot);
        bounds.push(off.min(n));
    }
    bounds.push(n);
    bounds.windows(2).map(|w| (w[0], w[1])).collect()
}

/// Parse a GTF file quickly, extracting only what we need.
///
/// For non-trivial files this fans out across the global rayon pool: the file
/// is read into memory, split at newline boundaries, parsed in parallel, and
/// the per-chunk results are merged. For tiny files the parallel overhead
/// isn't worth it; we fall back to a single in-process call.
pub fn load_gtf_fast(args: &Args, feature_type_override: Option<&str>) -> Result<AnnotationIndex> {
    let feature_type_str = feature_type_override.unwrap_or(&args.feature_type);
    let feature_type = feature_type_str.as_bytes();
    let gtf_attr = format!("{} \"", args.gene_id_attr);
    let gff3_attr = format!("{}=", args.gene_id_attr);

    // Tiny files: just stream through the sequential reader to keep startup
    // latency low. Threshold chosen so we don't spend 1ms spinning up rayon
    // for a 100-feature test fixture.
    const PARALLEL_THRESHOLD: u64 = 1 << 20; // 1 MB
    let file_size = std::fs::metadata(&args.annotation)
        .map(|m| m.len())
        .unwrap_or(0);

    let data: Vec<u8> = if file_size < PARALLEL_THRESHOLD {
        // Small file: still read into memory so we can use the same parser.
        let mut buf = Vec::new();
        open_reader(&args.annotation)?
            .read_to_end(&mut buf)
            .context("Failed to read annotation")?;
        buf
    } else {
        read_annotation_to_vec(&args.annotation)?
    };

    // Pick chunk count from the rayon pool (which the CLI sets to args.threads
    // before counting starts). Cap each chunk at >= 256 KB so we don't shred
    // small files into useless slivers.
    let pool_threads = rayon::current_num_threads().max(1);
    let max_chunks_by_size = (data.len() / (256 * 1024)).max(1);
    let num_chunks = pool_threads.min(max_chunks_by_size).max(1);

    let bounds = line_aligned_chunks(&data, num_chunks);

    let chunk_results: Vec<ChunkResult> = bounds
        .par_iter()
        .map(|&(s, e)| {
            parse_chunk(
                &data[s..e],
                feature_type,
                gtf_attr.as_bytes(),
                gff3_attr.as_bytes(),
            )
        })
        .collect();

    // Merge chunk-local IDs into global IDs. We accumulate byte-slice keys to
    // avoid per-feature string allocation; conversion to Arc<str> happens once
    // at the end.
    let mut global_chrom_to_id: FxHashMap<&[u8], u16> = FxHashMap::default();
    let mut id_to_chrom_bytes: Vec<&[u8]> = Vec::with_capacity(64);
    let mut global_gene_to_id: FxHashMap<&[u8], u32> = FxHashMap::default();
    let mut gene_names_bytes: Vec<&[u8]> = Vec::with_capacity(30_000);

    let total_features: usize = chunk_results.iter().map(|c| c.features.len()).sum();
    let mut features: Vec<Feature> = Vec::with_capacity(total_features);

    for chunk in &chunk_results {
        // Remap chunk-local chrom IDs to global.
        let chrom_remap: Vec<u16> = chunk
            .chrom_names
            .iter()
            .map(|name| {
                if let Some(&id) = global_chrom_to_id.get(*name) {
                    id
                } else {
                    let id = id_to_chrom_bytes.len() as u16;
                    id_to_chrom_bytes.push(*name);
                    global_chrom_to_id.insert(*name, id);
                    id
                }
            })
            .collect();

        // Remap chunk-local gene IDs.
        let gene_remap: Vec<u32> = chunk
            .gene_names
            .iter()
            .map(|name| {
                if let Some(&id) = global_gene_to_id.get(*name) {
                    id
                } else {
                    let id = gene_names_bytes.len() as u32;
                    gene_names_bytes.push(*name);
                    global_gene_to_id.insert(*name, id);
                    id
                }
            })
            .collect();

        for f in &chunk.features {
            features.push(Feature {
                start: f.start,
                end: f.end,
                gene_id: gene_remap[f.gene_id_local as usize],
                chrom_id: chrom_remap[f.chrom_id_local as usize],
                strand: f.strand,
            });
        }
    }

    features.sort_unstable_by_key(|f| (f.chrom_id, f.start));

    // One-time string conversion. ~50 chrom + ~80k gene names — amortized noise.
    let id_to_chrom: Vec<Arc<str>> = id_to_chrom_bytes
        .iter()
        .map(|b| Arc::from(String::from_utf8_lossy(b).as_ref()))
        .collect();
    let chrom_to_id_str: FxHashMap<Arc<str>, u16> = global_chrom_to_id
        .into_iter()
        .map(|(k, v)| (Arc::from(String::from_utf8_lossy(k).as_ref()), v))
        .collect();
    let gene_names: Vec<Arc<str>> = gene_names_bytes
        .iter()
        .map(|b| Arc::from(String::from_utf8_lossy(b).as_ref()))
        .collect();

    AnnotationIndex::new(chrom_to_id_str, id_to_chrom, gene_names, features)
}

/// Find a subsequence in a byte slice using SIMD-accelerated search
#[inline]
fn find_subsequence(haystack: &[u8], needle: &[u8]) -> Option<usize> {
    memmem::find(haystack, needle)
}
