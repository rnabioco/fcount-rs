use super::overlap::{Assignment, FeatureHit};
use super::stats::ReadCounters;
use crate::cli::Args;

/// Fixed-point multiplier for fractional counts
/// Using 1_000_000 to maintain 6 decimal places of precision
const FRACTION_MULTIPLIER: i64 = 1_000_000;

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
}

impl ThreadCounter {
    pub fn new(size: usize, args: &Args) -> Self {
        ThreadCounter {
            counts: vec![0; size],
            stats: ReadCounters::default(),
            hit_buffer: Vec::with_capacity(16),
            use_fractional: args.fractional_counting,
            feature_level: args.feature_level,
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
