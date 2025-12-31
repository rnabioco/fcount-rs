mod fast_gtf;
mod feature;
mod gtf;
mod index;

pub use feature::{Feature, Strand};
pub use index::{AnnotationIndex, ChromIndex};

use crate::cli::Args;
use anyhow::Result;

/// Load annotation from GTF/GFF file and build spatial index
/// Uses fast custom parser for GTF, falls back to bio for GFF3
pub fn load_annotation(args: &Args) -> Result<AnnotationIndex> {
    let path = &args.annotation;

    // Use bio parser for GFF3, fast parser for GTF
    let is_gff3 = path
        .extension()
        .map(|e| e.to_string_lossy().to_lowercase())
        .map(|e| e == "gff3" || e == "gff")
        .unwrap_or(false);

    if is_gff3 {
        gtf::load_gtf(args)
    } else {
        fast_gtf::load_gtf_fast(args)
    }
}
