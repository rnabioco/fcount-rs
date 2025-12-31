mod feature;
mod gtf;
mod index;

pub use feature::{Feature, Strand};
pub use index::{AnnotationIndex, ChromIndex};

use crate::cli::Args;
use anyhow::Result;

/// Load annotation from GTF/GFF file and build spatial index
pub fn load_annotation(args: &Args) -> Result<AnnotationIndex> {
    gtf::load_gtf(args)
}
