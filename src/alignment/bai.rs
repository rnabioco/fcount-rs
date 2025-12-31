use anyhow::{bail, Context, Result};
use noodles_bam::bai;
use noodles_sam as sam;
use std::path::Path;
use std::process::Command;

/// Represents a chromosome region in the BAM file
pub struct ChromRegion {
    /// BAM reference sequence ID
    pub ref_id: usize,
    /// Annotation chromosome ID
    pub annotation_chrom_id: u16,
    /// Chromosome name for debugging
    pub name: String,
}

/// Load or create BAI index
///
/// Looks for the BAI index file in two locations:
/// 1. `<bam_path>.bai` (e.g., `sample.bam.bai`)
/// 2. `<bam_path>` with `.bam` replaced by `.bai` (e.g., `sample.bai`)
///
/// If neither exists, auto-generates the index using `samtools index`.
pub fn load_bai_index(bam_path: &Path) -> Result<bai::Index> {
    let bai_path = bam_path.with_extension("bam.bai");
    let alt_bai_path = format!("{}.bai", bam_path.display());

    if !bai_path.exists() && !Path::new(&alt_bai_path).exists() {
        // Auto-generate with samtools
        generate_bai_index(bam_path)?;
    }

    // Load the index - try primary location first, then alternative
    let index = bai::fs::read(&bai_path)
        .or_else(|_| bai::fs::read(&alt_bai_path))
        .with_context(|| format!("Failed to load BAI index for {}", bam_path.display()))?;

    Ok(index)
}

/// Generate BAI index using samtools
fn generate_bai_index(bam_path: &Path) -> Result<()> {
    log::info!("Generating BAI index for {}...", bam_path.display());

    let output = Command::new("samtools")
        .arg("index")
        .arg(bam_path)
        .output()
        .context("Failed to run samtools. Is samtools installed and in PATH?")?;

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        bail!("samtools index failed: {}", stderr);
    }

    log::info!("BAI index generated successfully");
    Ok(())
}

/// Get chromosome regions from BAM header that have annotations
///
/// Returns a list of `ChromRegion` structs for chromosomes that exist in both
/// the BAM header and the annotation file. This allows parallel processing
/// of each chromosome region.
pub fn get_chrom_regions(
    header: &sam::Header,
    chrom_to_id: &rustc_hash::FxHashMap<std::sync::Arc<str>, u16>,
) -> Vec<ChromRegion> {
    header
        .reference_sequences()
        .iter()
        .enumerate()
        .filter_map(|(ref_id, (name, _))| {
            let name_str = std::str::from_utf8(name.as_ref()).ok()?;
            let annotation_chrom_id = chrom_to_id.get(name_str)?;
            Some(ChromRegion {
                ref_id,
                annotation_chrom_id: *annotation_chrom_id,
                name: name_str.to_string(),
            })
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use rustc_hash::FxHashMap;
    use std::sync::Arc;

    #[test]
    fn test_get_chrom_regions_empty() {
        let header = sam::Header::default();
        let chrom_to_id: FxHashMap<Arc<str>, u16> = FxHashMap::default();

        let regions = get_chrom_regions(&header, &chrom_to_id);
        assert!(regions.is_empty());
    }
}
