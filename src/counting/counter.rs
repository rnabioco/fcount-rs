use super::constants::FRACTION_MULTIPLIER;
use super::overlap::{Assignment, FeatureHit};
use super::stats::ReadCounters;
use crate::cli::Args;

/// Timing statistics for profiling hot paths
#[derive(Debug, Default, Clone)]
pub struct TimingStats {
    /// Time spent parsing BAM records (nanoseconds)
    pub parse_ns: u64,
    /// Time spent hashing read names (nanoseconds)
    pub hash_ns: u64,
    /// Time spent in mate tracker lookups (nanoseconds)
    pub mate_lookup_ns: u64,
    /// Time spent querying interval trees (nanoseconds)
    pub query_ns: u64,
    /// Number of records parsed
    pub records_parsed: u64,
    /// Number of interval tree queries
    pub queries_performed: u64,
}

impl TimingStats {
    pub fn merge(&mut self, other: &TimingStats) {
        self.parse_ns += other.parse_ns;
        self.hash_ns += other.hash_ns;
        self.mate_lookup_ns += other.mate_lookup_ns;
        self.query_ns += other.query_ns;
        self.records_parsed += other.records_parsed;
        self.queries_performed += other.queries_performed;
    }
}

/// Thread-local counter for accumulating counts
#[derive(Debug)]
pub struct ThreadCounter {
    /// Count table: gene_id or feature_idx -> count
    /// For fractional counting, stores fixed-point values (count * FRACTION_MULTIPLIER)
    pub counts: Vec<i64>,
    /// QC statistics
    pub stats: ReadCounters,
    /// Buffer for collecting feature hits (reused to avoid allocation)
    pub hit_buffer: Vec<FeatureHit>,
    /// Whether we're using fractional counting
    use_fractional: bool,
    /// Whether we're counting at feature level
    feature_level: bool,
    /// Timing statistics for profiling
    pub timing: TimingStats,
}

impl ThreadCounter {
    pub fn new(size: usize, args: &Args) -> Self {
        ThreadCounter {
            counts: vec![0; size],
            stats: ReadCounters::default(),
            hit_buffer: Vec::with_capacity(64),
            use_fractional: args.fractional_counting,
            feature_level: args.feature_level,
            timing: TimingStats::default(),
        }
    }

    /// Apply an assignment to the count table
    pub fn apply_assignment(&mut self, assignment: Assignment, nh: u8, args: &Args) {
        match assignment {
            Assignment::Unique(idx) => {
                let target_idx = if self.feature_level {
                    idx.feature_idx as usize
                } else {
                    idx.gene_id as usize
                };

                let count = self.calculate_count(nh, 1);
                self.counts[target_idx] += count;
                self.stats.assigned += 1;
            }

            Assignment::MultiOverlap(hits) if args.allow_multi_overlap => {
                let num_hits = hits.len();
                for hit in hits {
                    let target_idx = if self.feature_level {
                        hit.feature_idx as usize
                    } else {
                        hit.gene_id as usize
                    };

                    let count = self.calculate_count(nh, num_hits);
                    self.counts[target_idx] += count;
                }
                self.stats.assigned += 1;
            }

            Assignment::MultiOverlap(_) => {
                self.stats.unassigned_ambiguous += 1;
            }

            Assignment::NoFeature => {
                self.stats.unassigned_no_features += 1;
            }

            Assignment::Ambiguous => {
                self.stats.unassigned_ambiguous += 1;
            }
        }
    }

    /// Calculate count value considering NH tag and multi-overlap
    fn calculate_count(&self, nh: u8, num_targets: usize) -> i64 {
        if self.use_fractional {
            // Fractional count: 1 / (NH * num_targets)
            FRACTION_MULTIPLIER / (nh as i64 * num_targets as i64)
        } else {
            // Integer count: 1 for each target
            1
        }
    }

    /// Convert fixed-point counts to floating point for output
    pub fn get_float_counts(&self) -> Vec<f64> {
        if self.use_fractional {
            self.counts
                .iter()
                .map(|&c| c as f64 / FRACTION_MULTIPLIER as f64)
                .collect()
        } else {
            self.counts.iter().map(|&c| c as f64).collect()
        }
    }

    /// Merge another counter's results into this one
    pub fn merge(&mut self, other: &ThreadCounter) {
        for (i, &count) in other.counts.iter().enumerate() {
            self.counts[i] += count;
        }
        self.stats.merge(&other.stats);
    }
}
