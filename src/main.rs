use anyhow::Result;
use clap::Parser;
use log::info;
use std::time::Instant;

mod alignment;
mod annotation;
mod cli;
mod counting;
mod output;

use cli::Args;

fn main() -> Result<()> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    let args = Args::parse();
    let total_start = Instant::now();

    info!("fcount-rs v{}", env!("CARGO_PKG_VERSION"));
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
    let counts = if args.threads > 1 {
        if args.paired_end {
            info!(
                "Using parallel paired-end pipeline with {} threads",
                args.threads
            );
            counting::count_reads_parallel_paired(&args, &annotation)?
        } else {
            info!(
                "Using parallel single-end pipeline with {} threads",
                args.threads
            );
            counting::count_reads_parallel(&args, &annotation)?
        }
    } else {
        counting::count_reads(&args, &annotation)?
    };
    let bam_elapsed = bam_start.elapsed();
    info!("BAM processing completed in {:.2}s", bam_elapsed.as_secs_f64());

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
