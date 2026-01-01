use anyhow::{Context, Result};
use bio::io::gff::{self, GffType};
use bio_types::strand::Strand as BioStrand;
use log::debug;
use rustc_hash::FxHashMap;
use std::path::Path;
use std::sync::Arc;

use super::feature::{Feature, Strand};
use super::index::AnnotationIndex;
use super::io::open_reader;
use crate::cli::Args;

/// Detect GFF type from filename (handles .gz suffix)
fn detect_gff_type(path: &Path) -> GffType {
    let filename = path
        .file_name()
        .map(|s| s.to_string_lossy().to_lowercase())
        .unwrap_or_default();

    if filename.contains(".gff3") || filename.contains(".gff") {
        GffType::GFF3
    } else {
        GffType::GTF2
    }
}

/// Load GTF/GFF file and build annotation index
pub fn load_gtf(args: &Args) -> Result<AnnotationIndex> {
    let path = &args.annotation;
    let gff_type = detect_gff_type(path);

    debug!("Parsing annotation file as {:?}", gff_type);

    let reader = open_reader(path)?;
    let mut reader = gff::Reader::new(reader, gff_type);

    // Build lookups for chromosome and gene names
    let mut chrom_to_id: FxHashMap<Arc<str>, u16> = FxHashMap::default();
    let mut id_to_chrom: Vec<Arc<str>> = Vec::new();
    let mut gene_name_to_id: FxHashMap<Arc<str>, u32> = FxHashMap::default();
    let mut gene_names: Vec<Arc<str>> = Vec::new();
    let mut features: Vec<Feature> = Vec::new();

    for result in reader.records() {
        let record = result.with_context(|| "Failed to parse GTF/GFF record")?;

        // Filter by feature type
        if record.feature_type() != args.feature_type {
            continue;
        }

        // Get chromosome ID (or create new)
        let chrom_name: Arc<str> = Arc::from(record.seqname());
        let chrom_id = *chrom_to_id.entry(chrom_name.clone()).or_insert_with(|| {
            let id = id_to_chrom.len() as u16;
            id_to_chrom.push(chrom_name);
            id
        });

        // Get gene ID from attributes
        let gene_id_str = record
            .attributes()
            .get(&args.gene_id_attr)
            .with_context(|| {
                format!(
                    "Missing '{}' attribute in GTF record at {}:{}-{}",
                    args.gene_id_attr,
                    record.seqname(),
                    record.start(),
                    record.end()
                )
            })?;

        // Get gene ID index (or create new)
        let gene_name: Arc<str> = Arc::from(gene_id_str.as_str());
        let gene_id = *gene_name_to_id.entry(gene_name.clone()).or_insert_with(|| {
            let id = gene_names.len() as u32;
            gene_names.push(gene_name);
            id
        });

        // Parse coordinates (rust-bio uses 1-based, inclusive)
        let start = *record.start() as u32;
        let end = *record.end() as u32;

        // Parse strand
        let strand = match record.strand() {
            Some(BioStrand::Forward) => Strand::Forward,
            Some(BioStrand::Reverse) => Strand::Reverse,
            _ => Strand::Unknown,
        };

        features.push(Feature {
            gene_id,
            chrom_id,
            start,
            end,
            strand,
        });
    }

    debug!(
        "Parsed {} features from {} genes on {} chromosomes",
        features.len(),
        gene_names.len(),
        id_to_chrom.len()
    );

    // Sort features by (chrom_id, start) for efficient indexing
    features.sort_by_key(|f| (f.chrom_id, f.start));

    // Build the spatial index
    AnnotationIndex::new(chrom_to_id, id_to_chrom, gene_names, features)
}
