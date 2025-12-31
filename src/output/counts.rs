use anyhow::{Context, Result};
use std::fs::File;
use std::io::{BufWriter, Write};

use crate::annotation::AnnotationIndex;
use crate::cli::Args;
use crate::counting::CountResult;

/// Write count matrix to output file
pub fn write_counts(args: &Args, annotation: &AnnotationIndex, result: &CountResult) -> Result<()> {
    let file = File::create(&args.output)
        .with_context(|| format!("Failed to create output file: {}", args.output.display()))?;
    let mut writer = BufWriter::new(file);

    // Write header
    write!(writer, "Geneid\tChr\tStart\tEnd\tStrand\tLength")?;
    for bam_path in &args.bam_files {
        let sample_name = bam_path
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("unknown");
        write!(writer, "\t{}", sample_name)?;
    }
    writeln!(writer)?;

    if args.feature_level {
        // Feature-level output
        write_feature_level(args, annotation, result, &mut writer)?;
    } else {
        // Gene-level output
        write_gene_level(args, annotation, result, &mut writer)?;
    }

    writer.flush()?;
    Ok(())
}

fn write_gene_level(
    args: &Args,
    annotation: &AnnotationIndex,
    result: &CountResult,
    writer: &mut BufWriter<File>,
) -> Result<()> {
    // For each gene, collect its features to determine coordinates
    for (gene_idx, gene_name) in annotation.gene_names.iter().enumerate() {
        // Find all features for this gene
        let gene_features: Vec<_> = annotation
            .features
            .iter()
            .filter(|f| f.gene_id as usize == gene_idx)
            .collect();

        if gene_features.is_empty() {
            continue;
        }

        // Collect chromosome(s)
        let chrom_ids: std::collections::BTreeSet<_> =
            gene_features.iter().map(|f| f.chrom_id).collect();
        let chroms: Vec<_> = chrom_ids
            .iter()
            .filter_map(|&id| annotation.id_to_chrom.get(id as usize))
            .map(|s| s.as_ref())
            .collect();
        let chrom_str = chroms.join(";");

        // Collect starts
        let starts: Vec<_> = gene_features.iter().map(|f| f.start.to_string()).collect();
        let starts_str = starts.join(";");

        // Collect ends
        let ends: Vec<_> = gene_features.iter().map(|f| f.end.to_string()).collect();
        let ends_str = ends.join(";");

        // Strand
        let strand = gene_features
            .first()
            .map(|f| match f.strand {
                crate::annotation::Strand::Forward => "+",
                crate::annotation::Strand::Reverse => "-",
                crate::annotation::Strand::Unknown => ".",
            })
            .unwrap_or(".");

        // Length (sum of all feature lengths for this gene)
        let length: u32 = gene_features.iter().map(|f| f.len()).sum();

        // Count
        let count = result.counts.get(gene_idx).copied().unwrap_or(0);
        let count_str = if args.fractional_counting {
            format!("{:.6}", count as f64 / 1_000_000.0)
        } else {
            count.to_string()
        };

        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            gene_name, chrom_str, starts_str, ends_str, strand, length, count_str
        )?;
    }

    Ok(())
}

fn write_feature_level(
    args: &Args,
    annotation: &AnnotationIndex,
    result: &CountResult,
    writer: &mut BufWriter<File>,
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

        let count = result.counts.get(feat_idx).copied().unwrap_or(0);
        let count_str = if args.fractional_counting {
            format!("{:.6}", count as f64 / 1_000_000.0)
        } else {
            count.to_string()
        };

        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            gene_name,
            chrom,
            feature.start,
            feature.end,
            strand,
            feature.len(),
            count_str
        )?;
    }

    Ok(())
}
