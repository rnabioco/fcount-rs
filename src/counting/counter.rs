use super::constants::FRACTION_MULTIPLIER;
use super::overlap::{Assignment, FeatureHit};
use super::stats::ReadCounters;
use crate::cli::Args;

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
            hit_buffer: Vec::with_capacity(8),
            use_fractional: args.fractional_counting,
            feature_level: args.feature_level,
        }
    }

    /// Apply an assignment to the count table.
    ///
    /// `overlap_rejected` says whether any candidate feature was dropped by the
    /// `--min-overlap` / `--frac-overlap` thresholds while collecting hits. It
    /// only matters when nothing survived: such a read overlapped a feature but
    /// not by enough, which featureCounts reports as `Unassigned_Overlapping_Length`
    /// rather than `Unassigned_NoFeatures`.
    pub fn apply_assignment(
        &mut self,
        assignment: Assignment,
        nh: u16,
        args: &Args,
        overlap_rejected: bool,
    ) {
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
                if overlap_rejected {
                    self.stats.unassigned_overlap_length += 1;
                } else {
                    self.stats.unassigned_no_features += 1;
                }
            }

            Assignment::Ambiguous => {
                self.stats.unassigned_ambiguous += 1;
            }
        }
    }

    /// Calculate count value considering NH tag and multi-overlap
    fn calculate_count(&self, nh: u16, num_targets: usize) -> i64 {
        if self.use_fractional {
            // Fractional count: 1 / (NH * num_targets).
            // `parse_int_tag` already clamps NH to >= 1, but guard the divisor
            // anyway: an integer divide by zero here is a panic, not a wrong
            // number, and this runs once per assigned read on every input.
            let divisor = (nh as i64) * (num_targets as i64);
            FRACTION_MULTIPLIER / divisor.max(1)
        } else {
            // Integer count: 1 for each target
            1
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cli::OutputFormat;

    fn args_with(fractional: bool, multi_overlap: bool) -> Args {
        Args {
            annotation: std::path::PathBuf::new(),
            output: std::path::PathBuf::new(),
            bam_files: vec![],
            feature_type: "exon".to_string(),
            gene_id_attr: "gene_id".to_string(),
            feature_level: false,
            paired_end: false,
            require_both_aligned: false,
            no_chimeric: false,
            strand_mode: 0,
            count_multi_mapping: true,
            fractional_counting: fractional,
            primary_only: false,
            allow_multi_overlap: multi_overlap,
            min_overlap_bases: 1,
            min_overlap_fraction: 0.0,
            min_feature_overlap_fraction: 0.0,
            largest_overlap_only: false,
            min_mapping_quality: 0,
            ignore_duplicates: false,
            threads: 1,
            details_file: None,
            quiet: false,
            output_format: OutputFormat::Featurecounts,
        }
    }

    fn hit(gene_id: u32) -> FeatureHit {
        FeatureHit {
            feature_idx: gene_id,
            gene_id,
            overlap_len: 50,
        }
    }

    /// A read with no surviving hits is NoFeatures only when nothing was
    /// rejected by the overlap thresholds.
    #[test]
    fn test_no_feature_without_overlap_rejection() {
        let args = args_with(false, false);
        let mut c = ThreadCounter::new(4, &args);
        c.apply_assignment(Assignment::NoFeature, 1, &args, false);
        assert_eq!(c.stats.unassigned_no_features, 1);
        assert_eq!(c.stats.unassigned_overlap_length, 0);
    }

    /// Regression: a read that overlapped a feature but fell short of
    /// `--min-overlap` used to be filed under NoFeatures, leaving the
    /// Overlapping_Length row permanently at zero.
    #[test]
    fn test_overlap_rejection_reported_separately() {
        let args = args_with(false, false);
        let mut c = ThreadCounter::new(4, &args);
        c.apply_assignment(Assignment::NoFeature, 1, &args, true);
        assert_eq!(c.stats.unassigned_overlap_length, 1);
        assert_eq!(c.stats.unassigned_no_features, 0);
    }

    /// `overlap_rejected` must not disturb reads that did get assigned.
    #[test]
    fn test_overlap_rejection_ignored_when_assigned() {
        let args = args_with(false, false);
        let mut c = ThreadCounter::new(4, &args);
        c.apply_assignment(Assignment::Unique(hit(2)), 1, &args, true);
        assert_eq!(c.stats.assigned, 1);
        assert_eq!(c.counts[2], 1);
        assert_eq!(c.stats.unassigned_overlap_length, 0);
    }

    /// Fractional counting divides by NH — including for values that used to
    /// wrap in the old `u8` field.
    #[test]
    fn test_fractional_count_divides_by_nh() {
        let args = args_with(true, false);
        let mut c = ThreadCounter::new(4, &args);
        c.apply_assignment(Assignment::Unique(hit(0)), 4, &args, false);
        assert_eq!(c.counts[0], FRACTION_MULTIPLIER / 4);

        let mut c = ThreadCounter::new(4, &args);
        c.apply_assignment(Assignment::Unique(hit(0)), 300, &args, false);
        assert_eq!(c.counts[0], FRACTION_MULTIPLIER / 300);
    }

    /// NH of zero must not divide by zero, whatever produced it.
    #[test]
    fn test_fractional_count_survives_zero_nh() {
        let args = args_with(true, false);
        let mut c = ThreadCounter::new(4, &args);
        c.apply_assignment(Assignment::Unique(hit(0)), 0, &args, false);
        assert_eq!(c.counts[0], FRACTION_MULTIPLIER);
    }

    /// Multi-overlap splits across targets *and* by NH.
    #[test]
    fn test_fractional_multi_overlap_splits_both_ways() {
        let args = args_with(true, true);
        let mut c = ThreadCounter::new(4, &args);
        c.apply_assignment(
            Assignment::MultiOverlap(vec![hit(0), hit(1)]),
            2,
            &args,
            false,
        );
        assert_eq!(c.counts[0], FRACTION_MULTIPLIER / 4);
        assert_eq!(c.counts[1], FRACTION_MULTIPLIER / 4);
        assert_eq!(c.stats.assigned, 1);
    }
}
