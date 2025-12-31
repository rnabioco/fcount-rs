use anyhow::Result;
use clap::Parser;
use log::info;

mod alignment;
mod annotation;
mod cli;
mod counting;
mod output;

use cli::Args;

fn main() -> Result<()> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    let args = Args::parse();

    info!("fcount-rs v{}", env!("CARGO_PKG_VERSION"));
    info!("Loading annotation from: {}", args.annotation.display());

    // Load annotation and build index
    let annotation = annotation::load_annotation(&args)?;
    info!(
        "Loaded {} features from {} genes across {} chromosomes",
        annotation.features.len(),
        annotation.gene_names.len(),
        annotation.chrom_to_id.len()
    );

    // Process BAM files and count
    // Use parallel processing for single-end mode with multiple threads
    let counts = if !args.paired_end && args.threads > 1 {
        info!(
            "Using parallel producer-consumer pipeline with {} threads",
            args.threads
        );
        counting::count_reads_parallel(&args, &annotation)?
    } else {
        counting::count_reads(&args, &annotation)?
    };

    // Write output
    output::write_counts(&args, &annotation, &counts)?;
    output::write_summary(&args, &counts)?;

    info!("Done!");
    Ok(())
}
