use clap::Parser;
use std::path::PathBuf;

#[derive(Parser, Debug, Clone)]
#[command(name = "fcount-rs")]
#[command(author, version, about = "Ultrafast feature counting for RNA-seq data")]
#[command(
    long_about = "A high-performance Rust alternative to featureCounts for counting reads \
    that map to genomic features such as genes, exons, and transcripts."
)]
pub struct Args {
    /// GTF/GFF annotation file
    #[arg(short = 'a', long = "annotation")]
    pub annotation: PathBuf,

    /// Output count matrix file
    #[arg(short = 'o', long = "output")]
    pub output: PathBuf,

    /// Input BAM/SAM file(s)
    #[arg(required = true)]
    pub bam_files: Vec<PathBuf>,

    // === Counting options ===
    /// Feature type to count (e.g., exon, gene)
    #[arg(short = 't', long = "type", default_value = "exon")]
    pub feature_type: String,

    /// GTF attribute for gene ID
    #[arg(short = 'g', long = "gene-id", default_value = "gene_id")]
    pub gene_id_attr: String,

    /// Count at feature level instead of gene level
    #[arg(short = 'f', long = "feature-level")]
    pub feature_level: bool,

    // === Paired-end options ===
    /// Count fragments instead of reads (paired-end mode)
    #[arg(short = 'p', long = "paired")]
    pub paired_end: bool,

    /// Require both ends to be aligned
    #[arg(short = 'B', long = "both-aligned")]
    pub require_both_aligned: bool,

    /// Count chimeric fragments (mates on different chromosomes)
    #[arg(short = 'C', long = "count-chimeric")]
    pub count_chimeric: bool,

    /// Minimum fragment length for paired-end
    #[arg(long = "min-frag-len", default_value = "50")]
    pub min_fragment_length: u32,

    /// Maximum fragment length for paired-end
    #[arg(long = "max-frag-len", default_value = "600")]
    pub max_fragment_length: u32,

    // === Strand options ===
    /// Strand-specificity: 0=unstranded, 1=stranded, 2=reversely stranded
    #[arg(short = 's', long = "strand", default_value = "0")]
    pub strand_mode: u8,

    // === Multi-mapping options ===
    /// Count multi-mapping reads
    #[arg(short = 'M', long = "multi-mapping")]
    pub count_multi_mapping: bool,

    /// Use fractional counting for multi-mappers
    #[arg(long = "fraction")]
    pub fractional_counting: bool,

    /// Count primary alignments only
    #[arg(long = "primary")]
    pub primary_only: bool,

    // === Overlap options ===
    /// Allow reads to be assigned to multiple overlapping features
    #[arg(short = 'O', long = "multi-overlap")]
    pub allow_multi_overlap: bool,

    /// Minimum number of overlapping bases
    #[arg(long = "min-overlap", default_value = "1")]
    pub min_overlap_bases: u32,

    /// Minimum fraction of read that must overlap feature
    #[arg(long = "frac-overlap", default_value = "0.0")]
    pub min_overlap_fraction: f32,

    /// Minimum fraction of feature that must be overlapped
    #[arg(long = "frac-overlap-feature", default_value = "0.0")]
    pub min_feature_overlap_fraction: f32,

    /// Assign reads to the feature with the largest overlap only
    #[arg(long = "largest-overlap")]
    pub largest_overlap_only: bool,

    // === Filtering options ===
    /// Minimum mapping quality
    #[arg(short = 'Q', long = "min-mapq", default_value = "0")]
    pub min_mapping_quality: u8,

    /// Ignore duplicate reads (FLAG 0x400)
    #[arg(long = "ignore-dup")]
    pub ignore_duplicates: bool,

    // === Performance options ===
    /// Number of threads
    #[arg(short = 'T', long = "threads", default_value = "1")]
    pub threads: usize,

    // === Output options ===
    /// Output detailed assignment per read
    #[arg(short = 'R', long = "details")]
    pub details_file: Option<PathBuf>,

    /// Suppress progress output
    #[arg(short = 'q', long = "quiet")]
    pub quiet: bool,
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
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StrandMode {
    Unstranded,
    Stranded,
    ReverselyStranded,
}
