use clap::{builder::styling, Parser};
use std::path::PathBuf;

fn styles() -> styling::Styles {
    styling::Styles::styled()
        .header(styling::AnsiColor::Green.on_default().bold())
        .usage(styling::AnsiColor::Green.on_default().bold())
        .literal(styling::AnsiColor::Cyan.on_default().bold())
        .placeholder(styling::AnsiColor::Yellow.on_default())
}

#[derive(Parser, Debug, Clone)]
#[command(name = "fcount")]
#[command(author, version, about = "Ultrafast feature counting for RNA-seq data")]
#[command(styles = styles())]
#[command(
    long_about = "A high-performance Rust alternative to featureCounts for counting reads \
    that map to genomic features such as genes, exons, and transcripts."
)]
pub struct Args {
    // ============ Required Arguments ============
    /// GTF/GFF annotation file
    #[arg(short = 'a', long = "annotation", help_heading = "Input/Output")]
    pub annotation: PathBuf,

    /// Output count matrix file
    #[arg(short = 'o', long = "output", help_heading = "Input/Output")]
    pub output: PathBuf,

    /// Input BAM/SAM file(s)
    #[arg(required = true, help_heading = "Input/Output")]
    pub bam_files: Vec<PathBuf>,

    // ============ Feature Selection ============
    /// Feature type to count (e.g., exon, gene)
    #[arg(
        short = 't',
        long = "type",
        default_value = "exon",
        help_heading = "Feature Selection"
    )]
    pub feature_type: String,

    /// GTF attribute for gene ID
    #[arg(
        short = 'g',
        long = "gene-id",
        default_value = "gene_id",
        help_heading = "Feature Selection"
    )]
    pub gene_id_attr: String,

    /// Count at feature level instead of gene level
    #[arg(
        short = 'f',
        long = "feature-level",
        help_heading = "Feature Selection"
    )]
    pub feature_level: bool,

    // ============ Read Filtering ============
    /// Minimum mapping quality (0-255)
    #[arg(
        short = 'Q',
        long = "min-mapq",
        default_value = "0",
        help_heading = "Read Filtering"
    )]
    pub min_mapping_quality: u8,

    /// Count primary alignments only (skip secondary/supplementary)
    #[arg(long = "primary", help_heading = "Read Filtering")]
    pub primary_only: bool,

    /// Ignore duplicate reads (FLAG 0x400)
    #[arg(long = "ignore-dup", help_heading = "Read Filtering")]
    pub ignore_duplicates: bool,

    // ============ Paired-End Mode ============
    /// Count fragments instead of reads (paired-end mode)
    #[arg(short = 'p', long = "paired", help_heading = "Paired-End Mode")]
    pub paired_end: bool,

    /// Require both ends to be aligned
    #[arg(short = 'B', long = "both-aligned", help_heading = "Paired-End Mode")]
    pub require_both_aligned: bool,

    // ============ Strandedness ============
    /// Strand-specificity: 0=unstranded, 1=stranded, 2=reversely stranded
    #[arg(
        short = 's',
        long = "strand",
        default_value = "0",
        help_heading = "Strandedness"
    )]
    pub strand_mode: u8,

    // ============ Multi-Mapping Reads ============
    /// Count multi-mapping reads (reads with NH > 1)
    #[arg(short = 'M', long = "multi-mapping", help_heading = "Multi-Mapping")]
    pub count_multi_mapping: bool,

    /// Use fractional counting (1/NH) for multi-mappers
    #[arg(long = "fraction", help_heading = "Multi-Mapping")]
    pub fractional_counting: bool,

    // ============ Overlap Handling ============
    /// Allow reads to overlap multiple features/genes
    #[arg(short = 'O', long = "multi-overlap", help_heading = "Overlap Handling")]
    pub allow_multi_overlap: bool,

    /// Minimum number of overlapping bases required
    #[arg(
        long = "min-overlap",
        default_value = "1",
        help_heading = "Overlap Handling"
    )]
    pub min_overlap_bases: u32,

    /// Minimum fraction of read that must overlap feature (0.0-1.0)
    #[arg(
        long = "frac-overlap",
        default_value = "0.0",
        help_heading = "Overlap Handling"
    )]
    pub min_overlap_fraction: f32,

    /// Minimum fraction of feature that must be overlapped (0.0-1.0)
    #[arg(
        long = "frac-overlap-feature",
        default_value = "0.0",
        help_heading = "Overlap Handling"
    )]
    pub min_feature_overlap_fraction: f32,

    /// Assign reads to the feature with the largest overlap only
    #[arg(long = "largest-overlap", help_heading = "Overlap Handling")]
    pub largest_overlap_only: bool,

    // ============ Performance ============
    /// Number of threads for parallel processing (0 = auto-detect)
    #[arg(
        short = 'T',
        long = "threads",
        default_value = "0",
        help_heading = "Performance"
    )]
    pub threads: usize,

    // ============ Output Options ============
    /// Output detailed assignment per read (not yet implemented)
    #[arg(short = 'R', long = "details", help_heading = "Output Options")]
    pub details_file: Option<PathBuf>,

    /// Suppress progress output
    #[arg(short = 'q', long = "quiet", help_heading = "Output Options")]
    pub quiet: bool,

    // ============ Not Yet Implemented ============
    /// Do NOT count chimeric fragments (mates on different chromosomes)
    #[arg(short = 'C', long = "no-chimeric", help_heading = "Paired-End Mode")]
    pub no_chimeric: bool,
}

impl Args {
    /// Get the strand mode as an enum
    pub fn strand_mode(&self) -> StrandMode {
        match self.strand_mode {
            0 => StrandMode::Unstranded,
            1 => StrandMode::Stranded,
            2 => StrandMode::ReverselyStranded,
            _ => StrandMode::Unstranded,
        }
    }

    /// Check if we need to calculate overlap length
    pub fn need_overlap_length(&self) -> bool {
        self.min_overlap_bases > 1
            || self.min_overlap_fraction > 0.0
            || self.min_feature_overlap_fraction > 0.0
            || self.largest_overlap_only
    }

    /// Get effective thread count (auto-detect if 0)
    pub fn effective_threads(&self) -> usize {
        if self.threads == 0 {
            std::thread::available_parallelism()
                .map(|n| n.get())
                .unwrap_or(4)
        } else {
            self.threads
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StrandMode {
    Unstranded,
    Stranded,
    ReverselyStranded,
}
