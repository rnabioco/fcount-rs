#![allow(dead_code)]

use anyhow::Result;
use clap::Parser;
use log::info;
use std::time::Instant;

#[cfg(not(target_env = "msvc"))]
use tikv_jemallocator::Jemalloc;

#[cfg(not(target_env = "msvc"))]
#[global_allocator]
static GLOBAL: Jemalloc = Jemalloc;

mod alignment;
mod annotation;
mod cli;
mod counting;
mod output;

use cli::Args;

/// Detect if BAM file contains paired-end reads by checking first record
fn detect_paired_end(bam_path: &std::path::Path) -> Result<bool> {
    use noodles_bam as bam;
    use noodles_sam::alignment::record::Flags;
    use std::fs::File;

    let mut reader = File::open(bam_path).map(bam::io::Reader::new)?;
    let _header = reader.read_header()?;

    // Read first record
    for result in reader.records() {
        let record: bam::Record = result?;
        let flags = record.flags();
        // Check if PAIRED flag is set
        return Ok(flags.contains(Flags::SEGMENTED));
    }

    // No records found
    Ok(false)
}

fn main() -> Result<()> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    let mut args = Args::parse();
    let total_start = Instant::now();

    // Auto-detect thread count
    let effective_threads = args.effective_threads();
    args.threads = effective_threads;

    // Auto-detect paired-end mode from first BAM file
    if !args.paired_end && !args.bam_files.is_empty() {
        match detect_paired_end(&args.bam_files[0]) {
            Ok(true) => {
                info!("Auto-detected paired-end reads, enabling paired-end mode");
                args.paired_end = true;
            }
            Ok(false) => {}
            Err(e) => {
                info!("Could not auto-detect read type: {}", e);
            }
        }
    }

    info!("fcount v{}", env!("CARGO_PKG_VERSION"));
    info!("Loading annotation from: {}", args.annotation.display());

    // Load annotation and build index
    let gtf_start = Instant::now();
    let annotation = annotation::load_annotation(&args)?;
    let gtf_elapsed = gtf_start.elapsed();
    info!(
        "Loaded {} features from {} genes across {} chromosomes in {:.2}s",
        annotation.features.len(),
        annotation.gene_names.len(),
        annotation.chrom_to_id.len(),
        gtf_elapsed.as_secs_f64()
    );

    // Process BAM files and count
    // Use parallel processing when multiple threads are available
    let bam_start = Instant::now();
    let num_files = args.bam_files.len();
    let threads_per_file = 4.min(args.threads).max(1);
    let counts = if args.threads > 1 {
        if args.paired_end {
            info!(
                "Processing {} BAM files in paired-end mode ({} threads, {} per file)",
                num_files, args.threads, threads_per_file
            );
            counting::count_reads_parallel_paired(&args, &annotation)?
        } else {
            info!(
                "Processing {} BAM files in single-end mode ({} threads, {} per file)",
                num_files, args.threads, threads_per_file
            );
            counting::count_reads_parallel(&args, &annotation)?
        }
    } else {
        counting::count_reads(&args, &annotation)?
    };
    let bam_elapsed = bam_start.elapsed();
    info!(
        "BAM processing completed in {:.2}s",
        bam_elapsed.as_secs_f64()
    );

    // Write output
    let output_start = Instant::now();
    output::write_counts(&args, &annotation, &counts)?;
    output::write_summary(&args, &counts)?;
    let output_elapsed = output_start.elapsed();

    let total_elapsed = total_start.elapsed();
    info!(
        "Done! Total: {:.2}s (GTF: {:.2}s, BAM: {:.2}s, Output: {:.2}s)",
        total_elapsed.as_secs_f64(),
        gtf_elapsed.as_secs_f64(),
        bam_elapsed.as_secs_f64(),
        output_elapsed.as_secs_f64()
    );
    Ok(())
}
