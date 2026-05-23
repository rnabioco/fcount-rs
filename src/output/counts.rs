use anyhow::{Context, Result};
use flate2::Compression;
use flate2::write::GzEncoder;
use rustc_hash::FxHashMap;
use std::fs::File;
use std::io::{BufWriter, Write};

use crate::annotation::AnnotationIndex;
use crate::cli::{Args, OutputFormat};
use crate::counting::CountResult;

/// Sort + merge `intervals` in place and return the collapsed length
/// (sum of merged interval lengths). Uses 1-based inclusive coordinates
/// to match featureCounts.
///
/// In-place variant lets callers reuse a single scratch `Vec` across all
/// genes instead of allocating per call.
fn collapsed_length_inplace(intervals: &mut Vec<(u32, u32)>) -> u32 {
    if intervals.is_empty() {
        return 0;
    }
    intervals.sort_unstable_by_key(|&(start, _)| start);

    let mut write_idx = 0;
    for read_idx in 1..intervals.len() {
        let (start, end) = intervals[read_idx];
        let last = intervals[write_idx];
        if start <= last.1 + 1 {
            intervals[write_idx].1 = last.1.max(end);
        } else {
            write_idx += 1;
            intervals[write_idx] = (start, end);
        }
    }
    intervals.truncate(write_idx + 1);
    intervals.iter().map(|&(s, e)| e - s + 1).sum()
}

/// Write count matrix to output file (gzipped if filename ends in .gz)
pub fn write_counts(args: &Args, annotation: &AnnotationIndex, result: &CountResult) -> Result<()> {
    let file = File::create(&args.output)
        .with_context(|| format!("Failed to create output file: {}", args.output.display()))?;

    // Check if output should be gzipped
    let use_gzip = args
        .output
        .extension()
        .map(|ext| ext == "gz")
        .unwrap_or(false);

    if use_gzip {
        let encoder = GzEncoder::new(file, Compression::default());
        let mut writer = BufWriter::new(encoder);
        write_counts_inner(args, annotation, result, &mut writer)?;
        writer.into_inner()?.finish()?;
    } else {
        let mut writer = BufWriter::new(file);
        write_counts_inner(args, annotation, result, &mut writer)?;
        writer.flush()?;
    }

    Ok(())
}

fn write_counts_inner<W: Write>(
    args: &Args,
    annotation: &AnnotationIndex,
    result: &CountResult,
    writer: &mut BufWriter<W>,
) -> Result<()> {
    match args.output_format {
        OutputFormat::Dexseq => write_dexseq_format(args, annotation, result, writer),
        OutputFormat::Featurecounts => write_featurecounts_format(args, annotation, result, writer),
    }
}

fn write_featurecounts_format<W: Write>(
    args: &Args,
    annotation: &AnnotationIndex,
    result: &CountResult,
    writer: &mut BufWriter<W>,
) -> Result<()> {
    // Write comment line with program info (like featureCounts)
    write!(
        writer,
        "# Program:fcount v{}; Command:",
        env!("CARGO_PKG_VERSION")
    )?;
    write!(writer, "\"fcount\" ")?;
    write!(writer, "\"-a\" \"{}\" ", args.annotation.display())?;
    write!(writer, "\"-o\" \"{}\" ", args.output.display())?;
    for bam in &args.bam_files {
        write!(writer, "\"{}\" ", bam.path.display())?;
    }
    writeln!(writer)?;

    // Write header (use sample names for cleaner output)
    write!(writer, "Geneid\tChr\tStart\tEnd\tStrand\tLength")?;
    for bam in &args.bam_files {
        write!(writer, "\t{}", bam.display_name())?;
    }
    writeln!(writer)?;

    if args.feature_level {
        // Feature-level output
        write_feature_level(args, annotation, result, writer)?;
    } else {
        // Gene-level output
        write_gene_level(args, annotation, result, writer)?;
    }

    Ok(())
}

/// Write DEXSeq-style output format.
///
/// Format: gene_id, exon_id, sample1, sample2, ...
/// exon_id is formatted as gene_id:E### with per-gene numbering
fn write_dexseq_format<W: Write>(
    args: &Args,
    annotation: &AnnotationIndex,
    result: &CountResult,
    writer: &mut BufWriter<W>,
) -> Result<()> {
    // Write header: gene_id, exon_id, sample1, sample2, ...
    write!(writer, "gene_id\texon_id")?;
    for bam in &args.bam_files {
        write!(writer, "\t{}", bam.display_name())?;
    }
    writeln!(writer)?;

    // Track exon number per gene for proper E001, E002, etc.
    let mut gene_exon_counts: FxHashMap<u32, u32> = FxHashMap::default();

    // For each feature, output gene_id and gene_id:E### with per-gene counter
    for (feat_idx, feature) in annotation.features.iter().enumerate() {
        // Skip features with zero counts across all samples
        let has_counts = result
            .counts_per_sample
            .iter()
            .any(|counts| counts.get(feat_idx).copied().unwrap_or(0) != 0);
        if !has_counts {
            continue;
        }

        let gene_name = annotation
            .gene_names
            .get(feature.gene_id as usize)
            .map(|s| s.as_ref())
            .unwrap_or("unknown");

        let exon_num = gene_exon_counts.entry(feature.gene_id).or_insert(0);
        *exon_num += 1;

        // Stream gene_id and exon_id (gene_name:E### with zero-padded number)
        // directly — `format!` here would allocate a String per feature.
        write!(writer, "{}\t{}:E{:03}", gene_name, gene_name, exon_num)?;

        for sample_counts in &result.counts_per_sample {
            let count = sample_counts.get(feat_idx).copied().unwrap_or(0);
            if args.fractional_counting {
                write!(writer, "\t{:.6}", count as f64 / 1_000_000.0)?;
            } else {
                write!(writer, "\t{}", count)?;
            }
        }

        writeln!(writer)?;
    }

    Ok(())
}

fn write_gene_level<W: Write>(
    args: &Args,
    annotation: &AnnotationIndex,
    result: &CountResult,
    writer: &mut BufWriter<W>,
) -> Result<()> {
    // Pre-group features by gene_id - O(n) instead of O(n*m)
    let mut gene_to_features: Vec<Vec<usize>> = vec![Vec::new(); annotation.gene_names.len()];
    for (feat_idx, feature) in annotation.features.iter().enumerate() {
        gene_to_features[feature.gene_id as usize].push(feat_idx);
    }

    // Scratch buffer reused for every gene's collapsed-length calc.
    let mut interval_scratch: Vec<(u32, u32)> = Vec::with_capacity(32);

    for (gene_idx, gene_name) in annotation.gene_names.iter().enumerate() {
        let feature_indices = &gene_to_features[gene_idx];
        if feature_indices.is_empty() {
            continue;
        }

        // Geneid
        writer.write_all(gene_name.as_bytes())?;
        writer.write_all(b"\t")?;

        // Chr — semicolon-separated, streamed directly (no Vec<&str> + join).
        for (i, &idx) in feature_indices.iter().enumerate() {
            if i > 0 {
                writer.write_all(b";")?;
            }
            let chrom_id = annotation.features[idx].chrom_id;
            if let Some(name) = annotation.id_to_chrom.get(chrom_id as usize) {
                writer.write_all(name.as_bytes())?;
            }
        }
        writer.write_all(b"\t")?;

        // Start
        for (i, &idx) in feature_indices.iter().enumerate() {
            if i > 0 {
                writer.write_all(b";")?;
            }
            write!(writer, "{}", annotation.features[idx].start)?;
        }
        writer.write_all(b"\t")?;

        // End
        for (i, &idx) in feature_indices.iter().enumerate() {
            if i > 0 {
                writer.write_all(b";")?;
            }
            write!(writer, "{}", annotation.features[idx].end)?;
        }
        writer.write_all(b"\t")?;

        // Strand
        for (i, &idx) in feature_indices.iter().enumerate() {
            if i > 0 {
                writer.write_all(b";")?;
            }
            let s: &[u8] = match annotation.features[idx].strand {
                crate::annotation::Strand::Forward => b"+",
                crate::annotation::Strand::Reverse => b"-",
                crate::annotation::Strand::Unknown => b".",
            };
            writer.write_all(s)?;
        }
        writer.write_all(b"\t")?;

        // Length (collapsed, computed via reused scratch buffer)
        interval_scratch.clear();
        for &idx in feature_indices {
            let f = &annotation.features[idx];
            interval_scratch.push((f.start, f.end));
        }
        let length = collapsed_length_inplace(&mut interval_scratch);
        write!(writer, "{}", length)?;

        // Counts — write directly, no intermediate String per cell.
        for sample_counts in &result.counts_per_sample {
            let count = sample_counts.get(gene_idx).copied().unwrap_or(0);
            if args.fractional_counting {
                write!(writer, "\t{:.6}", count as f64 / 1_000_000.0)?;
            } else {
                write!(writer, "\t{}", count)?;
            }
        }
        writer.write_all(b"\n")?;
    }

    Ok(())
}

fn write_feature_level<W: Write>(
    args: &Args,
    annotation: &AnnotationIndex,
    result: &CountResult,
    writer: &mut BufWriter<W>,
) -> Result<()> {
    for (feat_idx, feature) in annotation.features.iter().enumerate() {
        let gene_name = annotation
            .gene_names
            .get(feature.gene_id as usize)
            .map(|s| s.as_ref())
            .unwrap_or("unknown");

        let chrom = annotation
            .id_to_chrom
            .get(feature.chrom_id as usize)
            .map(|s| s.as_ref())
            .unwrap_or("unknown");

        let strand = match feature.strand {
            crate::annotation::Strand::Forward => "+",
            crate::annotation::Strand::Reverse => "-",
            crate::annotation::Strand::Unknown => ".",
        };

        // Write feature metadata
        write!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}",
            gene_name,
            chrom,
            feature.start,
            feature.end,
            strand,
            feature.len()
        )?;

        for sample_counts in &result.counts_per_sample {
            let count = sample_counts.get(feat_idx).copied().unwrap_or(0);
            if args.fractional_counting {
                write!(writer, "\t{:.6}", count as f64 / 1_000_000.0)?;
            } else {
                write!(writer, "\t{}", count)?;
            }
        }
        writeln!(writer)?;
    }

    Ok(())
}
