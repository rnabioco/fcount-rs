use anyhow::Result;
use clap::Parser;
use log::info;
use std::time::Instant;

#[cfg(not(target_env = "msvc"))]
use tikv_jemallocator::Jemalloc;

#[cfg(not(target_env = "msvc"))]
#[global_allocator]
static GLOBAL: Jemalloc = Jemalloc;

// Consume the library crate rather than re-declaring every module, so the
// codebase compiles once (not once as lib + once as bin).
use fcount_rs::cli::{Args, OutputFormat};
use fcount_rs::{annotation, counting, output};

/// Number of leading records inspected when auto-detecting paired-end input.
///
/// Checking only the first record — as this used to — silently ran the whole
/// file single-end whenever that one record happened to be unpaired, which also
/// makes detection depend on sort order. Any segmented record in the prefix is
/// decisive, so this scans until it finds one.
const PAIRED_DETECT_RECORDS: usize = 1000;

/// Detect whether a BAM contains paired-end reads by sampling its first records.
fn detect_paired_end(bam_path: &std::path::Path) -> Result<bool> {
    use noodles_bam as bam;
    use noodles_sam::alignment::record::Flags;
    use std::fs::File;

    let mut reader = File::open(bam_path).map(bam::io::Reader::new)?;
    let _header = reader.read_header()?;

    for result in reader.records().take(PAIRED_DETECT_RECORDS) {
        let record: bam::Record = result?;
        if record.flags().contains(Flags::SEGMENTED) {
            return Ok(true);
        }
    }

    // No segmented record in the sampled prefix (or no records at all).
    Ok(false)
}

fn main() -> Result<()> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    let mut args = Args::parse();
    let total_start = Instant::now();

    // DEXSeq format requires feature-level output - auto-enable it
    if args.output_format == OutputFormat::Dexseq && !args.feature_level {
        info!("DEXSeq format requires feature-level output, enabling --feature-level");
        args.feature_level = true;
    }

    // Auto-detect thread count
    let effective_threads = args.effective_threads();
    args.threads = effective_threads;

    // Size the global rayon pool to the user's thread budget *before* anything
    // that uses rayon — including parallel GTF parsing. Otherwise rayon
    // auto-initializes to all available cores on first use and `-t` becomes a
    // suggestion rather than a budget.
    let _ = rayon::ThreadPoolBuilder::new()
        .num_threads(effective_threads)
        .build_global();

    // Auto-detect paired-end mode from first BAM file
    if !args.paired_end && !args.bam_files.is_empty() {
        match detect_paired_end(&args.bam_files[0].path) {
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

    // Process BAM files and count. A single shared pipeline handles every
    // thread count: threads_per_file_for() collapses to 1 worker at -t 1.
    let bam_start = Instant::now();
    let num_files = args.bam_files.len();
    let threads_per_file = counting::threads_per_file_for(args.threads, num_files);
    let counts = if args.paired_end {
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
