use anyhow::Result;
use clap::Parser;
use log::info;

mod cli;
mod annotation;
mod alignment;
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
    let counts = counting::count_reads(&args, &annotation)?;

    // Write output
    output::write_counts(&args, &annotation, &counts)?;
    output::write_summary(&args, &counts)?;

    info!("Done!");
    Ok(())
}
