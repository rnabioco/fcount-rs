mod counter;
mod overlap;
mod stats;

pub use counter::ThreadCounter;
pub use overlap::Assignment;
pub use stats::ReadCounters;

use anyhow::Result;
use log::{debug, info};
use rayon::prelude::*;
use std::sync::atomic::{AtomicUsize, Ordering};

use crate::alignment::{AlignmentReader, AlignmentRecord, MateTracker, PendingMate};
use crate::annotation::AnnotationIndex;
use crate::cli::Args;

/// Result of counting reads across all BAM files
pub struct CountResult {
    /// Gene-level counts (or feature-level if specified)
    pub counts: Vec<i64>,
    /// Aggregated QC statistics
    pub stats: ReadCounters,
}

/// Count reads in all BAM files
pub fn count_reads(args: &Args, annotation: &AnnotationIndex) -> Result<CountResult> {
    // Initialize thread pool
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .ok(); // Ignore if already initialized

    let num_features = annotation.features.len();
    let num_genes = annotation.gene_names.len();
    let count_size = if args.feature_level {
        num_features
    } else {
        num_genes
    };

    // Process each BAM file and collect results
    let progress = AtomicUsize::new(0);

    let results: Vec<Result<(Vec<i64>, ReadCounters)>> = args
        .bam_files
        .par_iter()
        .map(|bam_path| {
            let result = process_bam_file(bam_path, args, annotation, count_size);

            let done = progress.fetch_add(1, Ordering::Relaxed) + 1;
            info!(
                "Processed {}/{} BAM files: {}",
                done,
                args.bam_files.len(),
                bam_path.display()
            );

            result
        })
        .collect();

    // Merge results
    let mut final_counts = vec![0i64; count_size];
    let mut final_stats = ReadCounters::default();

    for result in results {
        let (counts, stats) = result?;
        for (i, count) in counts.into_iter().enumerate() {
            final_counts[i] += count;
        }
        final_stats.merge(&stats);
    }

    debug!(
        "Total: {} assigned, {} unassigned",
        final_stats.assigned,
        final_stats.total_unassigned()
    );

    Ok(CountResult {
        counts: final_counts,
        stats: final_stats,
    })
}

fn process_bam_file(
    bam_path: &std::path::Path,
    args: &Args,
    annotation: &AnnotationIndex,
    count_size: usize,
) -> Result<(Vec<i64>, ReadCounters)> {
    let mut reader = AlignmentReader::open(bam_path)?;
    let mut counter = ThreadCounter::new(count_size, args);
    let mut mate_tracker = if args.paired_end {
        Some(MateTracker::new(100_000))
    } else {
        None
    };

    let mut record = AlignmentRecord {
        read_name: Vec::new(),
        flags: 0,
        chrom_id: None,
        start: 0,
        mapq: 0,
        intervals: smallvec::SmallVec::new(),
        nh: 1,
        mate_chrom_id: None,
        mate_start: 0,
        template_len: 0,
    };

    let mut records = reader.records();

    while records.read_record(annotation, &mut record)? {
        // Skip unmapped reads
        if record.is_unmapped() {
            counter.stats.unassigned_unmapped += 1;
            continue;
        }

        // Skip secondary/supplementary if primary-only mode
        if args.primary_only && (record.is_secondary() || record.is_supplementary()) {
            counter.stats.unassigned_secondary += 1;
            continue;
        }

        // Skip low quality
        if record.mapq < args.min_mapping_quality {
            counter.stats.unassigned_mapping_quality += 1;
            continue;
        }

        // Skip duplicates if requested
        if args.ignore_duplicates && record.is_duplicate() {
            counter.stats.unassigned_duplicate += 1;
            continue;
        }

        // Skip multi-mappers if not counting them
        if !args.count_multi_mapping && record.nh > 1 {
            counter.stats.unassigned_multimapping += 1;
            continue;
        }

        // Handle paired-end vs single-end
        if args.paired_end && record.is_paired() {
            process_paired_read(&record, &mut counter, args, annotation, &mut mate_tracker);
        } else {
            process_single_read(&record, &mut counter, args, annotation);
        }
    }

    Ok((counter.counts, counter.stats))
}

fn process_single_read(
    record: &AlignmentRecord,
    counter: &mut ThreadCounter,
    args: &Args,
    annotation: &AnnotationIndex,
) {
    let chrom_id = match record.chrom_id {
        Some(id) => id,
        None => {
            counter.stats.unassigned_no_features += 1;
            return;
        }
    };

    // Find overlapping features
    counter.hit_buffer.clear();
    for interval in &record.intervals {
        for (feat_idx, feature) in annotation.find_overlapping(chrom_id, interval.start, interval.end) {
            // Check strand if stranded mode
            if !overlap::check_strand(record, feature, args) {
                continue;
            }

            // Calculate overlap if needed
            let overlap_len = if args.need_overlap_length() {
                crate::alignment::total_overlap(
                    &record.intervals,
                    feature.start,
                    feature.end,
                )
            } else {
                1
            };

            // Check minimum overlap
            if !overlap::check_overlap_thresholds(
                overlap_len,
                &record.intervals,
                feature,
                args,
            ) {
                continue;
            }

            counter.hit_buffer.push(overlap::FeatureHit {
                feature_idx: feat_idx,
                gene_id: feature.gene_id,
                overlap_len,
            });
        }
    }

    // Assign read
    let assignment = overlap::resolve_assignment(&counter.hit_buffer, args);
    counter.apply_assignment(assignment, record.nh, args);
}

fn process_paired_read(
    record: &AlignmentRecord,
    counter: &mut ThreadCounter,
    args: &Args,
    annotation: &AnnotationIndex,
    mate_tracker: &mut Option<MateTracker>,
) {
    let tracker = match mate_tracker {
        Some(t) => t,
        None => return,
    };

    // Check if mate is mapped
    if record.is_mate_unmapped() {
        if args.require_both_aligned {
            counter.stats.unassigned_singleton += 1;
            return;
        }
        // Process as single-end
        process_single_read(record, counter, args, annotation);
        return;
    }

    let chrom_id = match record.chrom_id {
        Some(id) => id,
        None => {
            counter.stats.unassigned_no_features += 1;
            return;
        }
    };

    // Check for chimeric reads
    if record.mate_chrom_id != record.chrom_id {
        if !args.count_chimeric {
            counter.stats.unassigned_chimeric += 1;
            return;
        }
    }

    // Try to find mate
    let mate = tracker.add_mate(
        &record.read_name,
        chrom_id,
        record.start,
        record.intervals.clone(),
        record.flags,
        record.mapq,
        record.nh,
    );

    if let Some(mate) = mate {
        // We have both mates - process as fragment
        process_fragment(record, &mate, counter, args, annotation);
    }
    // Otherwise, mate is now stored and we wait for the other end
}

fn process_fragment(
    record: &AlignmentRecord,
    mate: &PendingMate,
    counter: &mut ThreadCounter,
    args: &Args,
    annotation: &AnnotationIndex,
) {
    let chrom_id = match record.chrom_id {
        Some(id) => id,
        None => {
            counter.stats.unassigned_no_features += 1;
            return;
        }
    };

    // Collect hits from both mates
    counter.hit_buffer.clear();

    // Hits from current record
    for interval in &record.intervals {
        for (feat_idx, feature) in annotation.find_overlapping(chrom_id, interval.start, interval.end) {
            if !overlap::check_strand_paired(record, mate, feature, args) {
                continue;
            }

            let overlap_len = if args.need_overlap_length() {
                crate::alignment::total_overlap(
                    &record.intervals,
                    feature.start,
                    feature.end,
                ) + crate::alignment::total_overlap(
                    &mate.intervals,
                    feature.start,
                    feature.end,
                )
            } else {
                1
            };

            if !overlap::check_overlap_thresholds(
                overlap_len,
                &record.intervals,
                feature,
                args,
            ) {
                continue;
            }

            counter.hit_buffer.push(overlap::FeatureHit {
                feature_idx: feat_idx,
                gene_id: feature.gene_id,
                overlap_len,
            });
        }
    }

    // Hits from mate (if on same chromosome)
    if mate.chrom_id == chrom_id {
        for interval in &mate.intervals {
            for (feat_idx, feature) in annotation.find_overlapping(chrom_id, interval.start, interval.end) {
                // Skip if already counted
                if counter.hit_buffer.iter().any(|h| h.feature_idx == feat_idx) {
                    continue;
                }

                if !overlap::check_strand_paired(record, mate, feature, args) {
                    continue;
                }

                let overlap_len = if args.need_overlap_length() {
                    crate::alignment::total_overlap(
                        &record.intervals,
                        feature.start,
                        feature.end,
                    ) + crate::alignment::total_overlap(
                        &mate.intervals,
                        feature.start,
                        feature.end,
                    )
                } else {
                    1
                };

                if !overlap::check_overlap_thresholds(
                    overlap_len,
                    &record.intervals,
                    feature,
                    args,
                ) {
                    continue;
                }

                counter.hit_buffer.push(overlap::FeatureHit {
                    feature_idx: feat_idx,
                    gene_id: feature.gene_id,
                    overlap_len,
                });
            }
        }
    }

    // Assign fragment
    let assignment = overlap::resolve_assignment(&counter.hit_buffer, args);
    counter.apply_assignment(assignment, record.nh.max(mate.nh), args);
}
