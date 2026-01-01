use anyhow::{Context, Result};
use flate2::write::GzEncoder;
use flate2::Compression;
use rustc_hash::FxHashMap;
use std::fs::File;
use std::io::{BufWriter, Write};

use crate::annotation::AnnotationIndex;
use crate::cli::{Args, OutputFormat};
use crate::counting::CountResult;

/// Calculate collapsed (non-overlapping) length from a set of intervals.
/// Uses 1-based inclusive coordinates to match featureCounts.
fn calculate_collapsed_length(intervals: &[(u32, u32)]) -> u32 {
    if intervals.is_empty() {
        return 0;
    }

    // Sort by start position
    let mut sorted: Vec<(u32, u32)> = intervals.to_vec();
    sorted.sort_by_key(|&(start, _)| start);

    // Merge overlapping/adjacent intervals
    let mut merged: Vec<(u32, u32)> = Vec::with_capacity(sorted.len());
    for (start, end) in sorted {
        if let Some(last) = merged.last_mut() {
            // Merge if overlapping or adjacent (start <= prev_end + 1)
            if start <= last.1 + 1 {
                last.1 = last.1.max(end);
                continue;
            }
        }
        merged.push((start, end));
    }

    // Sum lengths using 1-based inclusive: end - start + 1
    merged.iter().map(|&(s, e)| e - s + 1).sum()
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

        let exon_id = format!("{}:E{:03}", gene_name, exon_num);

        // Write gene_id and exon_id
        write!(writer, "{}\t{}", gene_name, exon_id)?;

        // Write count for each sample
        for sample_counts in &result.counts_per_sample {
            let count = sample_counts.get(feat_idx).copied().unwrap_or(0);
            let count_str = if args.fractional_counting {
                format!("{:.6}", count as f64 / 1_000_000.0)
            } else {
                count.to_string()
            };
            write!(writer, "\t{}", count_str)?;
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

    // Now iterate genes with pre-computed feature indices
    for (gene_idx, gene_name) in annotation.gene_names.iter().enumerate() {
        let feature_indices = &gene_to_features[gene_idx];
        if feature_indices.is_empty() {
            continue;
        }

        // Collect chromosome(s) - one per feature to match featureCounts format
        let chroms: Vec<_> = feature_indices
            .iter()
            .filter_map(|&idx| {
                let chrom_id = annotation.features[idx].chrom_id;
                annotation
                    .id_to_chrom
                    .get(chrom_id as usize)
                    .map(|s| s.as_ref())
            })
            .collect();
        let chrom_str = chroms.join(";");

        // Collect starts
        let starts: Vec<_> = feature_indices
            .iter()
            .map(|&idx| annotation.features[idx].start.to_string())
            .collect();
        let starts_str = starts.join(";");

        // Collect ends
        let ends: Vec<_> = feature_indices
            .iter()
            .map(|&idx| annotation.features[idx].end.to_string())
            .collect();
        let ends_str = ends.join(";");

        // Strand - one per feature to match featureCounts format
        let strands: Vec<_> = feature_indices
            .iter()
            .map(|&idx| match annotation.features[idx].strand {
                crate::annotation::Strand::Forward => "+",
                crate::annotation::Strand::Reverse => "-",
                crate::annotation::Strand::Unknown => ".",
            })
            .collect();
        let strand_str = strands.join(";");

        // Length: collapsed (non-overlapping) length to match featureCounts
        let intervals: Vec<(u32, u32)> = feature_indices
            .iter()
            .map(|&idx| {
                let f = &annotation.features[idx];
                (f.start, f.end)
            })
            .collect();
        let length = calculate_collapsed_length(&intervals);

        // Write gene metadata
        write!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}",
            gene_name, chrom_str, starts_str, ends_str, strand_str, length
        )?;

        // Write counts for each sample
        for sample_counts in &result.counts_per_sample {
            let count = sample_counts.get(gene_idx).copied().unwrap_or(0);
            let count_str = if args.fractional_counting {
                format!("{:.6}", count as f64 / 1_000_000.0)
            } else {
                count.to_string()
            };
            write!(writer, "\t{}", count_str)?;
        }
        writeln!(writer)?;
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
            gene_name, chrom, feature.start, feature.end, strand, feature.len()
        )?;

        // Write counts for each sample
        for sample_counts in &result.counts_per_sample {
            let count = sample_counts.get(feat_idx).copied().unwrap_or(0);
            let count_str = if args.fractional_counting {
                format!("{:.6}", count as f64 / 1_000_000.0)
            } else {
                count.to_string()
            };
            write!(writer, "\t{}", count_str)?;
        }
        writeln!(writer)?;
    }

    Ok(())
}
