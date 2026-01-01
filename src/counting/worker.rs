//! Worker thread for parallel BAM counting.
//!
//! Each worker has its own count table and processes batches of BAM records
//! independently. Results are merged at the end.

use std::sync::Arc;

use crate::alignment::block_reader::RecordBatch;
use crate::alignment::minimal_parser::{get_record_size, parse_bam_record, MinimalRecord};
use crate::alignment::total_overlap;
use crate::annotation::{AnnotationIndex, Strand};
use crate::cli::{Args, StrandMode};

use super::constants::FRACTION_MULTIPLIER;
use super::overlap::{Assignment, FeatureHit};
use super::stats::ReadCounters;

/// Worker context for processing BAM records
pub struct Worker {
    /// Thread-local count table
    pub counts: Vec<i64>,
    /// Thread-local statistics
    pub stats: ReadCounters,
    /// Buffer for collecting feature hits (reused per record)
    hit_buffer: Vec<FeatureHit>,
    /// Reusable record buffer
    record: MinimalRecord,
    /// ref_id -> chrom_id mapping (shared)
    ref_to_chrom: Arc<Vec<Option<u16>>>,
    /// Arguments (shared reference)
    args: Arc<Args>,
    /// Whether we're using fractional counting
    use_fractional: bool,
    /// Whether we're counting at feature level
    feature_level: bool,
}

impl Worker {
    /// Create a new worker
    pub fn new(count_size: usize, ref_to_chrom: Arc<Vec<Option<u16>>>, args: Arc<Args>) -> Self {
        let use_fractional = args.fractional_counting;
        let feature_level = args.feature_level;

        Worker {
            counts: vec![0i64; count_size],
            stats: ReadCounters::default(),
            hit_buffer: Vec::with_capacity(64),
            record: MinimalRecord::new(),
            ref_to_chrom,
            args,
            use_fractional,
            feature_level,
        }
    }

    /// Process a batch of BAM records
    pub fn process_batch(&mut self, batch: &RecordBatch, annotation: &AnnotationIndex) {
        let data = &batch.data;
        let mut offset = 0;
        let need_read_name = self.args.paired_end;

        while offset + 4 <= data.len() {
            let record_size = get_record_size(&data[offset..]);
            if record_size == 0 {
                break;
            }

            let data_start = offset + 4;
            let data_end = data_start + record_size;

            if data_end > data.len() {
                break;
            }

            // Parse this record
            // Only parse NH tag if we need it for multi-mapper filtering
            if parse_bam_record(
                &data[data_start..data_end],
                &mut self.record,
                need_read_name,
                self.args.count_multi_mapping, // Only parse NH if counting multi-mappers
            )
            .is_ok()
            {
                self.process_record(annotation);
            }

            offset = data_end;
        }
    }

    /// Process a single parsed record
    fn process_record(&mut self, annotation: &AnnotationIndex) {
        // Skip unmapped reads
        if self.record.is_unmapped() {
            self.stats.unassigned_unmapped += 1;
            return;
        }

        // Skip secondary/supplementary if primary-only mode
        if self.args.primary_only && (self.record.is_secondary() || self.record.is_supplementary())
        {
            self.stats.unassigned_secondary += 1;
            return;
        }

        // Skip low quality
        if self.record.mapq < self.args.min_mapping_quality {
            self.stats.unassigned_mapping_quality += 1;
            return;
        }

        // Skip duplicates if requested
        if self.args.ignore_duplicates && self.record.is_duplicate() {
            self.stats.unassigned_duplicate += 1;
            return;
        }

        // Skip multi-mappers if not counting them
        if !self.args.count_multi_mapping && self.record.nh > 1 {
            self.stats.unassigned_multimapping += 1;
            return;
        }

        // Map ref_id to chrom_id
        let chrom_id = match self
            .ref_to_chrom
            .get(self.record.ref_id as usize)
            .and_then(|&id| id)
        {
            Some(id) => id,
            None => {
                self.stats.unassigned_no_features += 1;
                return;
            }
        };

        // For now, process as single-end (paired-end requires mate tracking
        // which is more complex with parallel processing)
        self.process_single_record(chrom_id, annotation);
    }

    /// Process a single-end record
    fn process_single_record(&mut self, chrom_id: u16, annotation: &AnnotationIndex) {
        // Find overlapping features using callback-based query (no allocation)
        self.hit_buffer.clear();

        let is_reverse = self.record.is_reverse();
        let need_overlap = self.args.need_overlap_length();
        let min_overlap_bases = self.args.min_overlap_bases;
        let min_overlap_fraction = self.args.min_overlap_fraction;
        let min_feature_overlap_fraction = self.args.min_feature_overlap_fraction;
        let strand_mode = self.args.strand_mode();

        // Get read length once if needed
        let read_len: u32 = if min_overlap_fraction > 0.0 {
            self.record.intervals.iter().map(|i| i.len()).sum()
        } else {
            0
        };

        // Process each interval
        for i in 0..self.record.intervals.len() {
            let interval = self.record.intervals[i];

            annotation.query_overlapping(
                chrom_id,
                interval.start,
                interval.end,
                |feat_idx, feature| {
                    // Check strand if stranded mode
                    let strand_ok = match strand_mode {
                        StrandMode::Unstranded => true,
                        StrandMode::Stranded => {
                            let read_strand = if is_reverse {
                                Strand::Reverse
                            } else {
                                Strand::Forward
                            };
                            feature.strand == read_strand || feature.strand == Strand::Unknown
                        }
                        StrandMode::ReverselyStranded => {
                            let read_strand = if is_reverse {
                                Strand::Forward
                            } else {
                                Strand::Reverse
                            };
                            feature.strand == read_strand || feature.strand == Strand::Unknown
                        }
                    };
                    if !strand_ok {
                        return;
                    }

                    // Calculate overlap if needed
                    let overlap_len = if need_overlap {
                        total_overlap(&self.record.intervals, feature.start, feature.end)
                    } else {
                        1
                    };

                    // Check minimum overlap bases
                    if overlap_len < min_overlap_bases {
                        return;
                    }

                    // Check minimum overlap fraction
                    if min_overlap_fraction > 0.0 && read_len > 0 {
                        let frac = overlap_len as f32 / read_len as f32;
                        if frac < min_overlap_fraction {
                            return;
                        }
                    }

                    // Check minimum feature overlap fraction
                    if min_feature_overlap_fraction > 0.0 {
                        let feature_len = feature.len();
                        if feature_len > 0 {
                            let frac = overlap_len as f32 / feature_len as f32;
                            if frac < min_feature_overlap_fraction {
                                return;
                            }
                        }
                    }

                    self.hit_buffer.push(FeatureHit {
                        feature_idx: feat_idx,
                        gene_id: feature.gene_id,
                        overlap_len,
                    });
                },
            );
        }

        // Assign read
        let assignment = self.resolve_assignment();
        self.apply_assignment(assignment);
    }

    /// Resolve assignment from hits
    fn resolve_assignment(&self) -> Assignment {
        if self.hit_buffer.is_empty() {
            return Assignment::NoFeature;
        }

        if self.hit_buffer.len() == 1 {
            return Assignment::Unique(self.hit_buffer[0].clone());
        }

        // Group hits by gene
        let mut unique_genes = rustc_hash::FxHashSet::default();
        for hit in &self.hit_buffer {
            unique_genes.insert(hit.gene_id);
        }

        if unique_genes.len() == 1 {
            // All hits from same gene
            if self.args.largest_overlap_only {
                let best = self
                    .hit_buffer
                    .iter()
                    .max_by_key(|h| h.overlap_len)
                    .unwrap();
                return Assignment::Unique(best.clone());
            }
            return Assignment::Unique(self.hit_buffer[0].clone());
        }

        // Multiple genes
        if self.args.allow_multi_overlap {
            if self.args.largest_overlap_only {
                let max_overlap = self.hit_buffer.iter().map(|h| h.overlap_len).max().unwrap();
                let best_hits: Vec<_> = self
                    .hit_buffer
                    .iter()
                    .filter(|h| h.overlap_len == max_overlap)
                    .cloned()
                    .collect();

                if best_hits.len() == 1 {
                    return Assignment::Unique(best_hits.into_iter().next().unwrap());
                }

                // Deduplicate by gene
                let mut seen = rustc_hash::FxHashSet::default();
                let deduped: Vec<_> = best_hits
                    .into_iter()
                    .filter(|h| seen.insert(h.gene_id))
                    .collect();
                return Assignment::MultiOverlap(deduped);
            }

            let mut seen = rustc_hash::FxHashSet::default();
            let deduped: Vec<_> = self
                .hit_buffer
                .iter()
                .filter(|h| seen.insert(h.gene_id))
                .cloned()
                .collect();
            return Assignment::MultiOverlap(deduped);
        }

        Assignment::Ambiguous
    }

    /// Apply assignment to count table
    fn apply_assignment(&mut self, assignment: Assignment) {
        match assignment {
            Assignment::Unique(hit) => {
                let target_idx = if self.feature_level {
                    hit.feature_idx as usize
                } else {
                    hit.gene_id as usize
                };

                let count = self.calculate_count(1);
                self.counts[target_idx] += count;
                self.stats.assigned += 1;
            }

            Assignment::MultiOverlap(hits) if self.args.allow_multi_overlap => {
                let num_hits = hits.len();
                for hit in hits {
                    let target_idx = if self.feature_level {
                        hit.feature_idx as usize
                    } else {
                        hit.gene_id as usize
                    };

                    let count = self.calculate_count(num_hits);
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

    /// Calculate count value
    #[inline]
    fn calculate_count(&self, num_targets: usize) -> i64 {
        if self.use_fractional {
            FRACTION_MULTIPLIER / (self.record.nh as i64 * num_targets as i64)
        } else {
            1
        }
    }

    /// Get final results
    pub fn into_results(self) -> (Vec<i64>, ReadCounters) {
        (self.counts, self.stats)
    }
}
