mod fast_gtf;
mod feature;
mod gtf;
mod index;
pub mod io;

pub use feature::{Feature, Strand};
pub use index::AnnotationIndex;

use crate::cli::Args;
use anyhow::Result;
use std::path::Path;

/// Check if file path indicates a GFF3 format (vs GTF)
/// Handles both plain and gzipped files (e.g., .gff3.gz)
fn is_gff3_format(path: &Path) -> bool {
    let filename = path
        .file_name()
        .map(|s| s.to_string_lossy().to_lowercase())
        .unwrap_or_default();

    // Check for .gff3 or .gff (with optional .gz suffix)
    filename.ends_with(".gff3")
        || filename.ends_with(".gff3.gz")
        || filename.ends_with(".gff")
        || filename.ends_with(".gff.gz")
}

/// Load annotation from GTF/GFF file and build spatial index
/// Uses fast custom parser for GTF, falls back to bio for GFF3
/// Supports gzipped files (.gz extension)
pub fn load_annotation(args: &Args) -> Result<AnnotationIndex> {
    let path = &args.annotation;

    if is_gff3_format(path) {
        gtf::load_gtf(args)
    } else {
        fast_gtf::load_gtf_fast(args)
    }
}
