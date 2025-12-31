use anyhow::{Context, Result};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use crate::cli::Args;
use crate::counting::CountResult;

/// Write QC summary to file
pub fn write_summary(args: &Args, result: &CountResult) -> Result<()> {
    // Summary file is output file with .summary extension
    let summary_path = Path::new(&args.output).with_extension("summary");

    let file = File::create(&summary_path)
        .with_context(|| format!("Failed to create summary file: {}", summary_path.display()))?;
    let mut writer = BufWriter::new(file);

    // Write header
    writeln!(writer, "Status\tCount\tPercentage")?;

    let stats = &result.stats;
    let total = stats.total() as f64;

    // Helper to write a row
    let write_row = |w: &mut BufWriter<File>, name: &str, count: u64| -> Result<()> {
        let pct = if total > 0.0 {
            count as f64 / total * 100.0
        } else {
            0.0
        };
        writeln!(w, "{}\t{}\t{:.2}%", name, count, pct)?;
        Ok(())
    };

    write_row(&mut writer, "Assigned", stats.assigned)?;
    write_row(
        &mut writer,
        "Unassigned_Unmapped",
        stats.unassigned_unmapped,
    )?;
    write_row(&mut writer, "Unassigned_Read_Type", 0)?; // Not used
    write_row(
        &mut writer,
        "Unassigned_Singleton",
        stats.unassigned_singleton,
    )?;
    write_row(
        &mut writer,
        "Unassigned_MappingQuality",
        stats.unassigned_mapping_quality,
    )?;
    write_row(&mut writer, "Unassigned_Chimera", stats.unassigned_chimeric)?;
    write_row(
        &mut writer,
        "Unassigned_FragmentLength",
        stats.unassigned_fragment_length,
    )?;
    write_row(
        &mut writer,
        "Unassigned_Duplicate",
        stats.unassigned_duplicate,
    )?;
    write_row(
        &mut writer,
        "Unassigned_MultiMapping",
        stats.unassigned_multimapping,
    )?;
    write_row(
        &mut writer,
        "Unassigned_Secondary",
        stats.unassigned_secondary,
    )?;
    write_row(
        &mut writer,
        "Unassigned_NoFeatures",
        stats.unassigned_no_features,
    )?;
    write_row(
        &mut writer,
        "Unassigned_Overlapping_Length",
        stats.unassigned_overlap_length,
    )?;
    write_row(
        &mut writer,
        "Unassigned_Ambiguity",
        stats.unassigned_ambiguous,
    )?;

    writer.flush()?;

    // Also print to stderr if not quiet
    if !args.quiet {
        eprintln!();
        eprintln!("=== Summary ===");
        eprintln!("{}", stats);
        eprintln!("Assignment rate: {:.1}%", stats.assignment_rate() * 100.0);
    }

    Ok(())
}
