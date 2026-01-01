mod fast_gtf;
mod feature;
mod index;
pub mod io;

pub use fast_gtf::{auto_select_feature_type, detect_feature_types};
pub use feature::{Feature, Strand};
pub use index::AnnotationIndex;

use crate::cli::Args;
use anyhow::{bail, Result};
use log::info;

/// Load annotation from GTF/GFF file and build spatial index.
/// Automatically detects feature type if default "exon" is not found.
/// Errors if user explicitly specifies a type that doesn't exist.
pub fn load_annotation(args: &Args) -> Result<AnnotationIndex> {
    // Scan file to detect available feature types
    let available_types = detect_feature_types(&args.annotation, 10000)?;

    if available_types.is_empty() {
        bail!(
            "No features found in annotation file: {}",
            args.annotation.display()
        );
    }

    // Determine which feature type to use
    let feature_type_override = if available_types.contains(&args.feature_type) {
        // User's choice (or default "exon") exists in file
        None
    } else if args.feature_type == "exon" {
        // User used default "exon" but it doesn't exist - auto-select
        let selected = auto_select_feature_type(&available_types)
            .ok_or_else(|| anyhow::anyhow!("No recognizable feature types in annotation file"))?;
        info!(
            "Feature type 'exon' not found, auto-detected: '{}'",
            selected
        );
        Some(selected.to_string())
    } else {
        // User explicitly specified a type that doesn't exist - error
        bail!(
            "Feature type '{}' not found in annotation file.\nAvailable types: {:?}",
            args.feature_type,
            available_types
        );
    };

    fast_gtf::load_gtf_fast(args, feature_type_override.as_deref())
}
