mod constants;
mod counter;
pub mod filtering;
pub mod overlap;
mod sharded_mate_tracker;
mod stats;
mod worker;

pub use constants::{
    FRACTION_MULTIPLIER, channel_buffer_size, mate_tracker_shards, threads_per_file,
};
pub use filtering::Filterable;

pub use counter::ThreadCounter;
pub use sharded_mate_tracker::{DeferredRead, ShardedMateTracker};
pub use stats::ReadCounters;

use anyhow::Result;
use crossbeam::channel;
use indicatif::{ProgressBar, ProgressStyle};
use log::{debug, info};
use rayon::prelude::*;
use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};

use crate::alignment::block_reader::BamBlockReader;
use crate::alignment::{AlignmentReader, AlignmentRecord, MateTracker, PendingMate};
use crate::annotation::AnnotationIndex;
use crate::cli::Args;
use worker::Worker;

/// Result of counting reads across all BAM files
pub struct CountResult {
    /// Per-sample counts: counts_per_sample[sample_idx][feature_idx]
    pub counts_per_sample: Vec<Vec<i64>>,
    /// Per-sample QC statistics
    pub stats_per_sample: Vec<ReadCounters>,
}

impl CountResult {
    /// Get aggregated stats across all samples (for summary output)
    pub fn aggregated_stats(&self) -> ReadCounters {
        let mut stats = ReadCounters::default();
        for s in &self.stats_per_sample {
            stats.merge(s);
        }
        stats
    }
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
        .map(|bam_input| {
            let result = process_bam_file(&bam_input.path, args, annotation, count_size);

            let done = progress.fetch_add(1, Ordering::Relaxed) + 1;
            info!(
                "Processed {}/{} BAM files: {}",
                done,
                args.bam_files.len(),
                bam_input.display_name()
            );

            result
        })
        .collect();

    // Collect per-sample results (don't merge)
    let mut counts_per_sample = Vec::with_capacity(args.bam_files.len());
    let mut stats_per_sample = Vec::with_capacity(args.bam_files.len());

    for result in results {
        let (counts, stats) = result?;
        counts_per_sample.push(counts);
        stats_per_sample.push(stats);
    }

    let mut aggregated = ReadCounters::default();
    for s in &stats_per_sample {
        aggregated.merge(s);
    }

    debug!(
        "Total: {} assigned, {} unassigned",
        aggregated.assigned,
        aggregated.total_unassigned()
    );

    Ok(CountResult {
        counts_per_sample,
        stats_per_sample,
    })
}

fn process_bam_file(
    bam_path: &std::path::Path,
    args: &Args,
    annotation: &AnnotationIndex,
    count_size: usize,
) -> Result<(Vec<i64>, ReadCounters)> {
    // Use open_with_annotation for pre-computed ref_id → chrom_id mapping
    let mut reader = AlignmentReader::open_with_annotation(bam_path, args.threads, annotation)?;
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

    // Handle orphan mates (reads whose mate was never found)
    if let Some(mut tracker) = mate_tracker {
        let orphan_count = tracker.pending_count();
        if orphan_count > 0 {
            log::debug!("{} orphan mates remaining after processing", orphan_count);
            // Process orphans as single-end reads
            for mate in tracker.drain() {
                process_orphan_mate(&mate, &mut counter, args, annotation);
            }
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
    use crate::cli::StrandMode;

    let chrom_id = match record.chrom_id {
        Some(id) => id,
        None => {
            counter.stats.unassigned_no_features += 1;
            return;
        }
    };

    // Precompute read_len once (only needed if min_overlap_fraction > 0)
    let read_len: u32 = if args.min_overlap_fraction > 0.0 {
        record.intervals.iter().map(|i| i.len()).sum()
    } else {
        0
    };

    // Precompute expected strand once
    let expected_strand = match args.strand_mode() {
        StrandMode::Unstranded => None,
        mode => Some(overlap::apply_strand_mode(
            overlap::strand_from_reverse(record.is_reverse_strand()),
            mode,
        )),
    };

    // Find overlapping features using callback-based query (no allocation)
    counter.hit_buffer.clear();
    for interval in &record.intervals {
        annotation.query_overlapping(
            chrom_id,
            interval.start,
            interval.end,
            |feat_idx, feature| {
                // Check strand if stranded mode
                if !overlap::check_strand_fast(expected_strand, feature.strand) {
                    return;
                }

                // Calculate overlap if needed
                let overlap_len = if args.need_overlap_length() {
                    crate::alignment::total_overlap(&record.intervals, feature.start, feature.end)
                } else {
                    1
                };

                // Check minimum overlap
                if !overlap::check_overlap_thresholds(overlap_len, read_len, feature, args) {
                    return;
                }

                counter.hit_buffer.push(overlap::FeatureHit {
                    feature_idx: feat_idx,
                    gene_id: feature.gene_id,
                    overlap_len,
                });
            },
        );
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

    // Check for chimeric reads (mates on different chromosomes)
    if record.mate_chrom_id != record.chrom_id && args.no_chimeric {
        counter.stats.unassigned_chimeric += 1;
        return;
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
    use crate::cli::StrandMode;

    let chrom_id = match record.chrom_id {
        Some(id) => id,
        None => {
            counter.stats.unassigned_no_features += 1;
            return;
        }
    };

    // Precompute read_len once (only needed if min_overlap_fraction > 0)
    let read_len: u32 = if args.min_overlap_fraction > 0.0 {
        record.intervals.iter().map(|i| i.len()).sum()
    } else {
        0
    };

    // Check if unstranded mode (skip strand checks entirely)
    let is_unstranded = args.strand_mode() == StrandMode::Unstranded;

    // Collect hits from both mates using callback-based query (no allocation)
    counter.hit_buffer.clear();

    // Hits from current record
    for interval in &record.intervals {
        annotation.query_overlapping(
            chrom_id,
            interval.start,
            interval.end,
            |feat_idx, feature| {
                if !is_unstranded && !overlap::check_strand_paired(record, mate, feature, args) {
                    return;
                }

                let overlap_len = if args.need_overlap_length() {
                    crate::alignment::total_overlap(&record.intervals, feature.start, feature.end)
                        + crate::alignment::total_overlap(
                            &mate.intervals,
                            feature.start,
                            feature.end,
                        )
                } else {
                    1
                };

                if !overlap::check_overlap_thresholds(overlap_len, read_len, feature, args) {
                    return;
                }

                counter.hit_buffer.push(overlap::FeatureHit {
                    feature_idx: feat_idx,
                    gene_id: feature.gene_id,
                    overlap_len,
                });
            },
        );
    }

    // Hits from mate (if on same chromosome)
    if mate.chrom_id == chrom_id {
        // Build set of already-seen features for O(1) duplicate check (reusing allocated set)
        counter.seen_features.clear();
        counter
            .seen_features
            .extend(counter.hit_buffer.iter().map(|h| h.feature_idx));

        for interval in &mate.intervals {
            annotation.query_overlapping(
                chrom_id,
                interval.start,
                interval.end,
                |feat_idx, feature| {
                    // Skip if already counted - O(1) lookup
                    if counter.seen_features.contains(&feat_idx) {
                        return;
                    }

                    if !is_unstranded && !overlap::check_strand_paired(record, mate, feature, args)
                    {
                        return;
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

                    if !overlap::check_overlap_thresholds(overlap_len, read_len, feature, args) {
                        return;
                    }

                    counter.hit_buffer.push(overlap::FeatureHit {
                        feature_idx: feat_idx,
                        gene_id: feature.gene_id,
                        overlap_len,
                    });
                },
            );
        }
    }

    // Assign fragment
    let assignment = overlap::resolve_assignment(&counter.hit_buffer, args);
    counter.apply_assignment(assignment, record.nh.max(mate.nh), args);
}

/// Process an orphan mate (read whose mate was never found) as single-end
fn process_orphan_mate(
    mate: &PendingMate,
    counter: &mut ThreadCounter,
    args: &Args,
    annotation: &AnnotationIndex,
) {
    use crate::cli::StrandMode;

    // Precompute read_len once (only needed if min_overlap_fraction > 0)
    let read_len: u32 = if args.min_overlap_fraction > 0.0 {
        mate.intervals.iter().map(|i| i.len()).sum()
    } else {
        0
    };

    // Precompute expected strand once
    let expected_strand = match args.strand_mode() {
        StrandMode::Unstranded => None,
        mode => Some(overlap::apply_strand_mode(
            overlap::strand_from_reverse(mate.is_reverse_strand()),
            mode,
        )),
    };

    counter.hit_buffer.clear();

    for interval in &mate.intervals {
        annotation.query_overlapping(
            mate.chrom_id,
            interval.start,
            interval.end,
            |feat_idx, feature| {
                // Check strand
                if !overlap::check_strand_fast(expected_strand, feature.strand) {
                    return;
                }

                let overlap_len = if args.need_overlap_length() {
                    crate::alignment::total_overlap(&mate.intervals, feature.start, feature.end)
                } else {
                    1
                };

                if !overlap::check_overlap_thresholds(overlap_len, read_len, feature, args) {
                    return;
                }

                counter.hit_buffer.push(overlap::FeatureHit {
                    feature_idx: feat_idx,
                    gene_id: feature.gene_id,
                    overlap_len,
                });
            },
        );
    }

    let assignment = overlap::resolve_assignment(&counter.hit_buffer, args);
    counter.apply_assignment(assignment, mate.nh, args);
}

/// Process a single BAM file using parallel producer-consumer pattern
/// Uses crossbeam scoped threads for safe borrowing
fn process_bam_parallel(
    bam_path: &std::path::Path,
    args: &Args,
    annotation: &AnnotationIndex,
    count_size: usize,
    threads_per_file: usize,
) -> Result<(Vec<i64>, ReadCounters)> {
    use crate::alignment::block_reader::RecordBatch;

    // Open BAM reader with specified threads for decompression
    let mut reader = BamBlockReader::open_with_threads(bam_path, annotation, threads_per_file)?;

    // Pre-create Arc before spawning workers (avoids N clones of underlying data)
    let ref_to_chrom_arc = Arc::new(reader.ref_to_chrom().to_vec());
    let args_arc = Arc::new(args.clone());

    // Number of worker threads
    let num_workers = threads_per_file.max(1);

    // Channel for batch distribution - larger buffer to prevent worker starvation
    let (tx, rx) = channel::bounded::<RecordBatch>(channel_buffer_size(num_workers));

    // Use crossbeam scoped threads for safe borrowing
    let result = crossbeam::scope(|scope| {
        // Spawn workers - each gets its own clone of annotation for cache locality
        let worker_handles: Vec<_> = (0..num_workers)
            .map(|_| {
                let rx = rx.clone();
                let ref_to_chrom_arc = Arc::clone(&ref_to_chrom_arc);
                let args_arc = Arc::clone(&args_arc);
                let local_annotation = annotation.clone(); // Per-thread copy for cache locality

                scope.spawn(move |_| {
                    let mut worker = Worker::new(count_size, ref_to_chrom_arc, args_arc);

                    // Process batches until channel closes
                    while let Ok(batch) = rx.recv() {
                        worker.process_batch(&batch, &local_annotation);
                    }

                    worker.into_results()
                })
            })
            .collect();

        // Drop extra receiver so channel closes when producer is done
        drop(rx);

        debug!("Starting read processing: {}", bam_path.display());

        // Producer: read batches and send to workers
        while let Some(batch) = reader.read_batch().expect("Failed to read batch") {
            if tx.send(batch).is_err() {
                break; // Workers gone
            }
        }

        // Close channel to signal workers to finish
        drop(tx);

        // Collect and merge results
        let mut final_counts = vec![0i64; count_size];
        let mut final_stats = ReadCounters::default();

        for handle in worker_handles {
            let (counts, stats) = handle.join().expect("Worker thread panicked");
            for (i, &count) in counts.iter().enumerate() {
                final_counts[i] += count;
            }
            final_stats.merge(&stats);
        }

        (final_counts, final_stats)
    });

    result.map_err(|_| anyhow::anyhow!("Scoped thread panicked"))
}

/// Count reads using parallel processing (for single-end mode only for now)
pub fn count_reads_parallel(args: &Args, annotation: &AnnotationIndex) -> Result<CountResult> {
    let num_features = annotation.features.len();
    let num_genes = annotation.gene_names.len();
    let count_size = if args.feature_level {
        num_features
    } else {
        num_genes
    };

    // Process BAM files in parallel
    let num_files = args.bam_files.len();
    let tpf = threads_per_file(args.threads);

    // Create progress bar
    let pb = ProgressBar::new(num_files as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} BAM files ({eta})")
            .unwrap()
            .progress_chars("#>-"),
    );

    let results: Vec<Result<(Vec<i64>, ReadCounters)>> = args
        .bam_files
        .par_iter()
        .map(|bam_input| {
            let result = process_bam_parallel(&bam_input.path, args, annotation, count_size, tpf);
            pb.inc(1);
            pb.set_message(bam_input.display_name());
            result
        })
        .collect();

    pb.finish_with_message("Done");

    // Collect per-sample results (don't merge)
    let mut counts_per_sample = Vec::with_capacity(args.bam_files.len());
    let mut stats_per_sample = Vec::with_capacity(args.bam_files.len());

    for result in results {
        let (counts, stats) = result?;
        counts_per_sample.push(counts);
        stats_per_sample.push(stats);
    }

    let mut aggregated = ReadCounters::default();
    for s in &stats_per_sample {
        aggregated.merge(s);
    }

    debug!(
        "Total: {} assigned, {} unassigned",
        aggregated.assigned,
        aggregated.total_unassigned()
    );

    Ok(CountResult {
        counts_per_sample,
        stats_per_sample,
    })
}

/// Process a single BAM file using parallel producer-consumer pattern for paired-end reads
/// Uses a sharded mate tracker for concurrent mate matching
fn process_bam_parallel_paired(
    bam_path: &std::path::Path,
    args: &Args,
    annotation: &AnnotationIndex,
    count_size: usize,
    threads_per_file: usize,
) -> Result<(Vec<i64>, ReadCounters, counter::TimingStats)> {
    use crate::alignment::block_reader::RecordBatch;
    use crate::alignment::minimal_parser::MinimalRecord;

    // Open BAM reader with specified threads for decompression
    let mut reader = BamBlockReader::open_with_threads(bam_path, annotation, threads_per_file)?;
    let ref_to_chrom: Vec<Option<u16>> = reader.ref_to_chrom().to_vec();

    // All threads are workers
    let num_workers = threads_per_file.max(1);

    // Create sharded mate tracker with 8x shards per worker to reduce contention
    let mate_tracker = Arc::new(ShardedMateTracker::new(mate_tracker_shards(num_workers)));

    // Large channel buffer for read-ahead
    let (tx, rx) = channel::bounded::<RecordBatch>(channel_buffer_size(num_workers));

    // Use crossbeam scoped threads for safe borrowing
    let result = crossbeam::scope(|scope| {
        // Spawn workers - each gets its own clone of annotation for cache locality
        let worker_handles: Vec<_> = (0..num_workers)
            .map(|_| {
                let rx = rx.clone();
                let ref_to_chrom = ref_to_chrom.clone();
                let mate_tracker = Arc::clone(&mate_tracker);
                let local_annotation = annotation.clone(); // Per-thread copy for cache locality

                scope.spawn(move |_| {
                    let mut counter = ThreadCounter::new(count_size, args);
                    let mut record = MinimalRecord::default();

                    // Process batches until channel closes
                    while let Ok(batch) = rx.recv() {
                        process_paired_batch(
                            &batch,
                            &mut record,
                            &mut counter,
                            &ref_to_chrom,
                            &local_annotation,
                            args,
                            &mate_tracker,
                        );
                    }

                    (counter.counts, counter.stats, counter.timing)
                })
            })
            .collect();

        // Drop extra receiver so channel closes when producer is done
        drop(rx);

        debug!("Starting read processing: {}", bam_path.display());

        // Producer: read batches and send to workers
        while let Some(batch) = reader.read_batch().expect("Failed to read batch") {
            if tx.send(batch).is_err() {
                break; // Workers gone
            }
        }

        // Close channel
        drop(tx);

        // Collect and merge results
        let mut final_counts = vec![0i64; count_size];
        let mut final_stats = ReadCounters::default();
        let mut final_timing = counter::TimingStats::default();

        for handle in worker_handles {
            let (counts, stats, timing) = handle.join().expect("Worker thread panicked");
            for (i, &count) in counts.iter().enumerate() {
                final_counts[i] += count;
            }
            final_stats.merge(&stats);
            final_timing.merge(&timing);
        }

        // Handle orphan mates (reads whose mate was never found)
        let orphans = mate_tracker.drain_all();
        if !orphans.is_empty() {
            debug!("{} orphan mates remaining after processing", orphans.len());
            // These are singletons - we could count them if require_both_aligned is false
            if !args.require_both_aligned {
                // Process orphans as single-end reads (reuse hit_buffer)
                let mut hit_buffer: Vec<overlap::FeatureHit> = Vec::with_capacity(8);
                for (_hash, deferred) in orphans {
                    process_deferred_as_single(
                        &deferred,
                        &mut final_counts,
                        &mut final_stats,
                        annotation,
                        args,
                        &mut hit_buffer,
                    );
                }
            } else {
                final_stats.unassigned_singleton += orphans.len() as u64;
            }
        }

        (final_counts, final_stats, final_timing)
    });

    result.map_err(|_| anyhow::anyhow!("Scoped thread panicked"))
}

/// Process a batch of paired-end records
fn process_paired_batch(
    batch: &crate::alignment::block_reader::RecordBatch,
    record: &mut crate::alignment::minimal_parser::MinimalRecord,
    counter: &mut ThreadCounter,
    ref_to_chrom: &[Option<u16>],
    annotation: &AnnotationIndex,
    args: &Args,
    mate_tracker: &ShardedMateTracker,
) {
    use crate::alignment::minimal_parser::get_record_size;

    let mut offset = 0;
    let data = &batch.data;

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

        // Parse record - need read name for mate tracking
        // Note: parse_bam_record expects data WITHOUT the 4-byte size prefix
        if crate::alignment::minimal_parser::parse_bam_record(
            &data[data_start..data_end],
            record,
            true,                     // need_read_name for paired-end
            args.count_multi_mapping, // Only parse NH if counting multi-mappers
        )
        .is_err()
        {
            offset = data_end;
            continue;
        }

        offset = data_end;

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

        // Get chromosome ID
        let chrom_id = if record.ref_id >= 0 && (record.ref_id as usize) < ref_to_chrom.len() {
            ref_to_chrom[record.ref_id as usize]
        } else {
            None
        };

        let chrom_id = match chrom_id {
            Some(id) => id,
            None => {
                counter.stats.unassigned_no_features += 1;
                continue;
            }
        };

        // Check if mate is mapped
        if record.is_mate_unmapped() {
            if args.require_both_aligned {
                counter.stats.unassigned_singleton += 1;
                continue;
            }
            // Process as single-end
            process_minimal_single_read(record, chrom_id, counter, args, annotation);
            continue;
        }

        // Check for chimeric reads (mates on different chromosomes)
        let mate_chrom_id =
            if record.mate_ref_id >= 0 && (record.mate_ref_id as usize) < ref_to_chrom.len() {
                ref_to_chrom[record.mate_ref_id as usize]
            } else {
                None
            };

        // Only mark as chimeric if both chromosomes are in annotation AND different
        // If mate_chrom_id is None, mate is on chromosome not in annotation - not chimeric
        if let Some(mate_chrom) = mate_chrom_id {
            if mate_chrom != chrom_id && args.no_chimeric {
                counter.stats.unassigned_chimeric += 1;
                continue;
            }
        }

        // Hash the read name and check for existing mate first
        let name_hash = ShardedMateTracker::hash_name(&record.read_name);

        // Try to get existing mate (avoids cloning intervals if mate exists)
        if let Some(mate) = mate_tracker.remove_mate(name_hash) {
            // Found mate - process as fragment (no clone needed)
            process_minimal_fragment(record, chrom_id, &mate, counter, args, annotation);
        } else {
            // No mate yet - create deferred read and store it
            let deferred = DeferredRead {
                chrom_id,
                start: record.pos as u32,
                intervals: record.intervals.clone(),
                flags: record.flags,
                mapq: record.mapq,
                nh: record.nh,
            };
            mate_tracker.insert(name_hash, deferred);
        }
    }
}

/// Process a MinimalRecord as a single-end read
fn process_minimal_single_read(
    record: &crate::alignment::minimal_parser::MinimalRecord,
    chrom_id: u16,
    counter: &mut ThreadCounter,
    args: &Args,
    annotation: &AnnotationIndex,
) {
    use crate::cli::StrandMode;

    // Precompute read_len once (only needed if min_overlap_fraction > 0)
    let read_len: u32 = if args.min_overlap_fraction > 0.0 {
        record.intervals.iter().map(|i| i.len()).sum()
    } else {
        0
    };

    // Precompute expected strand once
    let expected_strand = match args.strand_mode() {
        StrandMode::Unstranded => None,
        mode => Some(overlap::apply_strand_mode(
            overlap::strand_from_reverse(record.is_reverse()),
            mode,
        )),
    };

    counter.hit_buffer.clear();

    for interval in &record.intervals {
        // Use callback-based query to avoid allocation
        annotation.query_overlapping(
            chrom_id,
            interval.start,
            interval.end,
            |feat_idx, feature| {
                if !overlap::check_strand_fast(expected_strand, feature.strand) {
                    return;
                }

                let overlap_len = if args.need_overlap_length() {
                    crate::alignment::total_overlap(&record.intervals, feature.start, feature.end)
                } else {
                    1
                };

                if !overlap::check_overlap_thresholds(overlap_len, read_len, feature, args) {
                    return;
                }

                counter.hit_buffer.push(overlap::FeatureHit {
                    feature_idx: feat_idx,
                    gene_id: feature.gene_id,
                    overlap_len,
                });
            },
        );
    }

    let assignment = overlap::resolve_assignment(&counter.hit_buffer, args);
    counter.apply_assignment(assignment, record.nh, args);
}

/// Process a fragment (both mates found) from MinimalRecord
fn process_minimal_fragment(
    record: &crate::alignment::minimal_parser::MinimalRecord,
    chrom_id: u16,
    mate: &DeferredRead,
    counter: &mut ThreadCounter,
    args: &Args,
    annotation: &AnnotationIndex,
) {
    use crate::cli::StrandMode;

    // Precompute read_len once (only needed if min_overlap_fraction > 0)
    let read_len: u32 = if args.min_overlap_fraction > 0.0 {
        record.intervals.iter().map(|i| i.len()).sum()
    } else {
        0
    };

    // Check if unstranded mode (skip strand checks entirely)
    let is_unstranded = args.strand_mode() == StrandMode::Unstranded;

    counter.hit_buffer.clear();

    let record_strand = if record.is_reverse() {
        crate::annotation::Strand::Reverse
    } else {
        crate::annotation::Strand::Forward
    };

    let mate_strand = if mate.is_reverse_strand() {
        crate::annotation::Strand::Reverse
    } else {
        crate::annotation::Strand::Forward
    };

    // Hits from current record - use callback-based query
    for interval in &record.intervals {
        annotation.query_overlapping(
            chrom_id,
            interval.start,
            interval.end,
            |feat_idx, feature| {
                if !is_unstranded
                    && !overlap::check_strand_paired_with_strands(
                        record_strand,
                        mate_strand,
                        feature,
                        args,
                    )
                {
                    return;
                }

                let overlap_len = if args.need_overlap_length() {
                    crate::alignment::total_overlap(&record.intervals, feature.start, feature.end)
                        + crate::alignment::total_overlap(
                            &mate.intervals,
                            feature.start,
                            feature.end,
                        )
                } else {
                    1
                };

                if !overlap::check_overlap_thresholds(overlap_len, read_len, feature, args) {
                    return;
                }

                counter.hit_buffer.push(overlap::FeatureHit {
                    feature_idx: feat_idx,
                    gene_id: feature.gene_id,
                    overlap_len,
                });
            },
        );
    }

    // Hits from mate (if on same chromosome)
    if mate.chrom_id == chrom_id {
        // Build set of already-seen features for O(1) duplicate check (reusing allocated set)
        counter.seen_features.clear();
        counter
            .seen_features
            .extend(counter.hit_buffer.iter().map(|h| h.feature_idx));

        for interval in &mate.intervals {
            annotation.query_overlapping(
                chrom_id,
                interval.start,
                interval.end,
                |feat_idx, feature| {
                    // Skip if already counted - O(1) lookup
                    if counter.seen_features.contains(&feat_idx) {
                        return;
                    }

                    if !is_unstranded
                        && !overlap::check_strand_paired_with_strands(
                            record_strand,
                            mate_strand,
                            feature,
                            args,
                        )
                    {
                        return;
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

                    if !overlap::check_overlap_thresholds(overlap_len, read_len, feature, args) {
                        return;
                    }

                    counter.hit_buffer.push(overlap::FeatureHit {
                        feature_idx: feat_idx,
                        gene_id: feature.gene_id,
                        overlap_len,
                    });
                },
            );
        }
    }

    let assignment = overlap::resolve_assignment(&counter.hit_buffer, args);
    counter.apply_assignment(assignment, record.nh.max(mate.nh), args);
}

/// Process a deferred read as single-end (for orphan handling)
fn process_deferred_as_single(
    deferred: &DeferredRead,
    counts: &mut [i64],
    stats: &mut ReadCounters,
    annotation: &AnnotationIndex,
    args: &Args,
    hit_buffer: &mut Vec<overlap::FeatureHit>,
) {
    use crate::cli::StrandMode;

    hit_buffer.clear();

    // Precompute read_len once (only needed if min_overlap_fraction > 0)
    let read_len: u32 = if args.min_overlap_fraction > 0.0 {
        deferred.intervals.iter().map(|i| i.len()).sum()
    } else {
        0
    };

    // Precompute expected strand once
    let expected_strand = match args.strand_mode() {
        StrandMode::Unstranded => None,
        mode => Some(overlap::apply_strand_mode(
            overlap::strand_from_reverse(deferred.is_reverse_strand()),
            mode,
        )),
    };

    for interval in &deferred.intervals {
        annotation.query_overlapping(
            deferred.chrom_id,
            interval.start,
            interval.end,
            |feat_idx, feature| {
                if !overlap::check_strand_fast(expected_strand, feature.strand) {
                    return;
                }

                let overlap_len = if args.need_overlap_length() {
                    crate::alignment::total_overlap(&deferred.intervals, feature.start, feature.end)
                } else {
                    1
                };

                if !overlap::check_overlap_thresholds(overlap_len, read_len, feature, args) {
                    return;
                }

                hit_buffer.push(overlap::FeatureHit {
                    feature_idx: feat_idx,
                    gene_id: feature.gene_id,
                    overlap_len,
                });
            },
        );
    }

    let assignment = overlap::resolve_assignment(hit_buffer, args);

    // Apply assignment directly to counts/stats
    // Use same counting as ThreadCounter (1 for non-fractional, FRACTION_MULTIPLIER for fractional)
    let count_value = if args.fractional_counting {
        FRACTION_MULTIPLIER
    } else {
        1
    };

    match assignment {
        overlap::Assignment::Unique(hit) => {
            let id = if args.feature_level {
                hit.feature_idx as usize
            } else {
                hit.gene_id as usize
            };
            counts[id] += count_value;
            stats.assigned += 1;
        }
        overlap::Assignment::Ambiguous => {
            stats.unassigned_ambiguous += 1;
        }
        overlap::Assignment::NoFeature => {
            stats.unassigned_no_features += 1;
        }
        overlap::Assignment::MultiOverlap(hits) => {
            if args.allow_multi_overlap {
                let frac = count_value / hits.len() as i64;
                for hit in hits {
                    let id = if args.feature_level {
                        hit.feature_idx as usize
                    } else {
                        hit.gene_id as usize
                    };
                    counts[id] += frac;
                }
                stats.assigned += 1;
            } else {
                stats.unassigned_ambiguous += 1;
            }
        }
    }
}

/// Count reads using parallel processing for paired-end mode
pub fn count_reads_parallel_paired(
    args: &Args,
    annotation: &AnnotationIndex,
) -> Result<CountResult> {
    let num_features = annotation.features.len();
    let num_genes = annotation.gene_names.len();
    let count_size = if args.feature_level {
        num_features
    } else {
        num_genes
    };

    // Process BAM files in parallel
    let num_files = args.bam_files.len();
    let tpf = threads_per_file(args.threads);

    // Create progress bar
    let pb = ProgressBar::new(num_files as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} BAM files ({eta})")
            .unwrap()
            .progress_chars("#>-"),
    );

    let results: Vec<Result<(Vec<i64>, ReadCounters, counter::TimingStats)>> = args
        .bam_files
        .par_iter()
        .map(|bam_input| {
            let result =
                process_bam_parallel_paired(&bam_input.path, args, annotation, count_size, tpf);
            pb.inc(1);
            pb.set_message(bam_input.display_name());
            result
        })
        .collect();

    pb.finish_with_message("Done");

    // Collect per-sample results (don't merge)
    let mut counts_per_sample = Vec::with_capacity(args.bam_files.len());
    let mut stats_per_sample = Vec::with_capacity(args.bam_files.len());

    for result in results {
        let (counts, stats, _timing) = result?;
        counts_per_sample.push(counts);
        stats_per_sample.push(stats);
    }

    let mut aggregated = ReadCounters::default();
    for s in &stats_per_sample {
        aggregated.merge(s);
    }

    debug!(
        "Parallel paired: {} assigned, {} unassigned",
        aggregated.assigned,
        aggregated.total_unassigned()
    );

    Ok(CountResult {
        counts_per_sample,
        stats_per_sample,
    })
}
