mod constants;
mod counter;
pub mod overlap;
mod sharded_mate_tracker;
mod stats;

pub use constants::{
    BGZF_THREADS_CAP, FRACTION_MULTIPLIER, channel_buffer_size, mate_tracker_shards,
    threads_per_file_for,
};

pub use counter::ThreadCounter;
pub use sharded_mate_tracker::{DeferredRead, ShardedMateTracker};
pub use stats::ReadCounters;

use anyhow::Result;
use crossbeam::channel;
use indicatif::{ProgressBar, ProgressStyle};
use log::debug;
use rayon::prelude::*;
use std::sync::Arc;

use crate::alignment::block_reader::BamBlockReader;
use crate::annotation::AnnotationIndex;
use crate::cli::Args;

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

/// Build the per-file progress bar.
///
/// Honors `--quiet`: the bar draws to stderr, so leaving it enabled made the
/// flag's promise to "suppress progress output" false. The template carries no
/// `{msg}` placeholder, so callers deliberately don't set one.
fn new_progress_bar(num_files: usize, quiet: bool) -> ProgressBar {
    if quiet {
        return ProgressBar::hidden();
    }
    let pb = ProgressBar::new(num_files as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template(
                "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} BAM files ({eta})",
            )
            .unwrap()
            .progress_chars("#>-"),
    );
    pb
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
    use crate::alignment::minimal_parser::MinimalRecord;

    // BGZF decompression saturates well below record-processing rate; cap
    // separately from record workers so -t N doesn't 2× oversubscribe the OS
    // scheduler (see BGZF_THREADS_CAP doc).
    let bgzf_threads = threads_per_file.clamp(1, BGZF_THREADS_CAP);
    let mut reader = BamBlockReader::open_with_threads(bam_path, annotation, bgzf_threads)?;

    // Share the ref_id→chrom_id table across workers (small Vec, no reason to clone).
    let ref_to_chrom: Arc<[Option<u16>]> = Arc::from(reader.ref_to_chrom().to_vec());

    // Number of worker threads
    let num_workers = threads_per_file.max(1);

    // Channel for batch distribution - larger buffer to prevent worker starvation
    let (tx, rx) = channel::bounded::<RecordBatch>(channel_buffer_size(num_workers));

    // Use crossbeam scoped threads for safe borrowing
    let result = crossbeam::scope(|scope| {
        // Spawn workers — all share the borrowed annotation index. Cloning it
        // per thread duplicates the COITree N times and pushes the working set
        // out of L3, which measured slower for single-end. Each worker
        // accumulates into its own ThreadCounter and routes every record
        // through the shared process_minimal_single_read path.
        let worker_handles: Vec<_> = (0..num_workers)
            .map(|_| {
                let rx = rx.clone();
                let ref_to_chrom = Arc::clone(&ref_to_chrom);

                scope.spawn(move |_| {
                    let mut counter = ThreadCounter::new(count_size, args);
                    let mut record = MinimalRecord::default();

                    // Process batches until channel closes
                    while let Ok(batch) = rx.recv() {
                        process_single_batch(
                            &batch,
                            &mut record,
                            &mut counter,
                            &ref_to_chrom,
                            annotation,
                            args,
                        );
                    }

                    (counter.counts, counter.stats)
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

/// Process a batch of records in single-end mode.
///
/// Mirrors the filtering preamble of `process_paired_batch`, but every
/// surviving record is counted directly as single-end — there is no mate
/// tracking. Shares `process_minimal_single_read` so the overlap/resolve/count
/// logic is identical across all code paths.
fn process_single_batch(
    batch: &crate::alignment::block_reader::RecordBatch,
    record: &mut crate::alignment::minimal_parser::MinimalRecord,
    counter: &mut ThreadCounter,
    ref_to_chrom: &[Option<u16>],
    annotation: &AnnotationIndex,
    args: &Args,
) {
    use crate::alignment::minimal_parser::get_record_size;

    // Hoist per-record args field reads out of the inner loop.
    let count_multi_mapping = args.count_multi_mapping;
    let primary_only = args.primary_only;
    let min_mapping_quality = args.min_mapping_quality;
    let ignore_duplicates = args.ignore_duplicates;
    // We need the NH tag whenever NH affects the result: to *exclude* NH>1 reads
    // when not counting multi-mappers, or to weight them when fractional. (Only
    // when counting multi-mappers non-fractionally is NH irrelevant.)
    let need_nh = !count_multi_mapping || args.fractional_counting;

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

        // Single-end never needs the read name.
        if crate::alignment::minimal_parser::parse_bam_record(
            &data[data_start..data_end],
            record,
            false,
            need_nh,
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
        if primary_only && (record.is_secondary() || record.is_supplementary()) {
            counter.stats.unassigned_secondary += 1;
            continue;
        }

        // Skip low quality
        if record.mapq < min_mapping_quality {
            counter.stats.unassigned_mapping_quality += 1;
            continue;
        }

        // Skip duplicates if requested
        if ignore_duplicates && record.is_duplicate() {
            counter.stats.unassigned_duplicate += 1;
            continue;
        }

        // Skip multi-mappers if not counting them
        if !count_multi_mapping && record.nh > 1 {
            counter.stats.unassigned_multimapping += 1;
            continue;
        }

        // Map ref_id → chrom_id
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

        process_minimal_single_read(record, chrom_id, counter, args, annotation);
    }
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
    let tpf = threads_per_file_for(args.threads, num_files);

    let pb = new_progress_bar(num_files, args.quiet);

    let results: Vec<Result<(Vec<i64>, ReadCounters)>> = args
        .bam_files
        .par_iter()
        .map(|bam_input| {
            let result = process_bam_parallel(&bam_input.path, args, annotation, count_size, tpf);
            pb.inc(1);
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
) -> Result<(Vec<i64>, ReadCounters)> {
    use crate::alignment::block_reader::RecordBatch;
    use crate::alignment::minimal_parser::MinimalRecord;

    // BGZF decompression saturates well below record-processing rate; cap
    // separately from record workers so -t N doesn't 2× oversubscribe the OS
    // scheduler (see BGZF_THREADS_CAP doc).
    let bgzf_threads = threads_per_file.clamp(1, BGZF_THREADS_CAP);
    let mut reader = BamBlockReader::open_with_threads(bam_path, annotation, bgzf_threads)?;
    // Share the ref_id→chrom_id table across workers (small Vec, no reason to clone).
    let ref_to_chrom: Arc<[Option<u16>]> = Arc::from(reader.ref_to_chrom().to_vec());

    // All threads are workers
    let num_workers = threads_per_file.max(1);

    // Create sharded mate tracker with 8x shards per worker to reduce contention
    let mate_tracker = Arc::new(ShardedMateTracker::new(mate_tracker_shards(num_workers)));

    // Large channel buffer for read-ahead
    let (tx, rx) = channel::bounded::<RecordBatch>(channel_buffer_size(num_workers));

    // Use crossbeam scoped threads for safe borrowing
    let result = crossbeam::scope(|scope| {
        // Spawn workers — each clones AnnotationIndex (COITrees + features) for
        // thread-local cache locality. Measured: sharing the index via &annotation
        // *regressed* paired-end throughput by ~5–9% at 4/8 threads on chr22,
        // because the random tree traversal across mate pairs benefits more from
        // thread-local L1/L2 residency than the clone cost hurts. (Single-end is
        // different — see process_bam_parallel for the shared-index variant.)
        let worker_handles: Vec<_> = (0..num_workers)
            .map(|_| {
                let rx = rx.clone();
                let ref_to_chrom = Arc::clone(&ref_to_chrom);
                let mate_tracker = Arc::clone(&mate_tracker);
                let local_annotation = annotation.clone();

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

                    (counter.counts, counter.stats)
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

        for handle in worker_handles {
            let (counts, stats) = handle.join().expect("Worker thread panicked");
            for (i, &count) in counts.iter().enumerate() {
                final_counts[i] += count;
            }
            final_stats.merge(&stats);
        }

        // Handle orphan mates (reads whose mate was never found)
        let orphans = mate_tracker.drain_all();
        if !orphans.is_empty() {
            debug!("{} orphan mates remaining after processing", orphans.len());
            // These are singletons - we could count them if require_both_aligned is false
            if !args.require_both_aligned {
                // Count orphans as single-end reads through the same
                // ThreadCounter the workers use, then merge like any worker's
                // result — so NH weighting and stats stay identical by
                // construction rather than by two copies agreeing.
                let mut orphan_counter = ThreadCounter::new(count_size, args);
                for (_hash, deferred) in orphans {
                    process_deferred_as_single(&deferred, &mut orphan_counter, annotation, args);
                }
                for (slot, count) in final_counts.iter_mut().zip(&orphan_counter.counts) {
                    *slot += count;
                }
                final_stats.merge(&orphan_counter.stats);
            } else {
                final_stats.unassigned_singleton += orphans.len() as u64;
            }
        }

        (final_counts, final_stats)
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

    // Hoist per-record args field reads out of the inner loop. The Args struct
    // is &-borrowed so LTO can probably do this on its own, but making the
    // hoist explicit shortens the prologue of each iteration and removes any
    // possibility the optimizer leaves a redundant load behind.
    let count_multi_mapping = args.count_multi_mapping;
    let primary_only = args.primary_only;
    let min_mapping_quality = args.min_mapping_quality;
    let ignore_duplicates = args.ignore_duplicates;
    let require_both_aligned = args.require_both_aligned;
    let no_chimeric = args.no_chimeric;
    // We need the NH tag whenever NH affects the result: to *exclude* NH>1 reads
    // when not counting multi-mappers, or to weight them when fractional. Parsing
    // it only under count_multi_mapping (the old behavior) silently left every
    // multi-mapper at nh=1, so the `nh > 1` filter never fired.
    let need_nh = !count_multi_mapping || args.fractional_counting;

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
            true, // need_read_name for paired-end
            need_nh,
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
        if primary_only && (record.is_secondary() || record.is_supplementary()) {
            counter.stats.unassigned_secondary += 1;
            continue;
        }

        // Skip low quality
        if record.mapq < min_mapping_quality {
            counter.stats.unassigned_mapping_quality += 1;
            continue;
        }

        // Skip duplicates if requested
        if ignore_duplicates && record.is_duplicate() {
            counter.stats.unassigned_duplicate += 1;
            continue;
        }

        // Skip multi-mappers if not counting them
        if !count_multi_mapping && record.nh > 1 {
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
            if require_both_aligned {
                counter.stats.unassigned_singleton += 1;
                continue;
            }
            // Process as single-end
            process_minimal_single_read(record, chrom_id, counter, args, annotation);
            continue;
        }

        // Check for chimeric reads (mates on different chromosomes).
        // Only resolve the mate chrom when the chimeric filter is actually on;
        // otherwise the lookup is pure waste.
        if no_chimeric {
            let mate_chrom_id =
                if record.mate_ref_id >= 0 && (record.mate_ref_id as usize) < ref_to_chrom.len() {
                    ref_to_chrom[record.mate_ref_id as usize]
                } else {
                    None
                };

            // Only mark as chimeric if both chromosomes are in annotation AND
            // different. If mate_chrom_id is None, mate is on a chromosome not
            // in annotation — not chimeric.
            if let Some(mate_chrom) = mate_chrom_id
                && mate_chrom != chrom_id
            {
                counter.stats.unassigned_chimeric += 1;
                continue;
            }
        }

        // Hash was precomputed inline by the parser to avoid copying the name.
        // take_or_insert_with holds the shard lock once and only constructs
        // the DeferredRead (cloning the intervals SmallVec) when this is the
        // first mate of the pair.
        let name_hash = record.read_name_hash;
        let mate = mate_tracker.take_or_insert_with(name_hash, || DeferredRead {
            chrom_id,
            intervals: record.intervals.clone(),
            flags: record.flags,
            nh: record.nh,
        });
        if let Some(mate) = mate {
            process_minimal_fragment(record, chrom_id, &mate, counter, args, annotation);
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
    count_single_end(
        &record.intervals,
        chrom_id,
        record.is_reverse(),
        record.nh,
        counter,
        args,
        annotation,
    );
}

/// Count one read's worth of intervals as a single-end read.
///
/// Shared by the single-end path, by paired records whose mate is unmapped, and
/// by orphan mates left over after a paired run — all three used to carry their
/// own copy of this loop, and the orphan copy had drifted (it ignored NH under
/// `--fraction`). Taking `intervals` + `is_reverse` + `nh` rather than a record
/// type is what lets `DeferredRead` and `MinimalRecord` share it.
fn count_single_end(
    intervals: &[crate::alignment::Interval],
    chrom_id: u16,
    is_reverse: bool,
    nh: u16,
    counter: &mut ThreadCounter,
    args: &Args,
    annotation: &AnnotationIndex,
) {
    use crate::cli::StrandMode;

    // Precompute read_len once (only needed if min_overlap_fraction > 0)
    let read_len: u32 = if args.min_overlap_fraction > 0.0 {
        intervals.iter().map(|i| i.len()).sum()
    } else {
        0
    };

    // Precompute expected strand once
    let expected_strand = match args.strand_mode() {
        StrandMode::Unstranded => None,
        mode => Some(overlap::apply_strand_mode(
            overlap::strand_from_reverse(is_reverse),
            mode,
        )),
    };

    // Hoist out of the per-feature callback — answer is constant per read.
    let needs_overlap = args.need_overlap_length();

    counter.hit_buffer.clear();
    // Set when a feature overlapped but fell short of the overlap thresholds, so
    // a read that ends up with no hits can be reported as Overlapping_Length
    // rather than NoFeatures.
    let mut overlap_rejected = false;

    for interval in intervals {
        // Use callback-based query to avoid allocation
        annotation.query_overlapping(
            chrom_id,
            interval.start,
            interval.end,
            |feat_idx, feature| {
                if !overlap::check_strand_fast(expected_strand, feature.strand) {
                    return;
                }

                // When needs_overlap is false, every overlap threshold knob is
                // at its no-op default, so the threshold check trivially passes
                // with overlap_len=1 — skip the call entirely.
                let overlap_len = if needs_overlap {
                    let ol = crate::alignment::total_overlap(intervals, feature.start, feature.end);
                    if !overlap::check_overlap_thresholds(ol, read_len, feature, args) {
                        overlap_rejected = true;
                        return;
                    }
                    ol
                } else {
                    1
                };

                counter.hit_buffer.push(overlap::FeatureHit {
                    feature_idx: feat_idx,
                    gene_id: feature.gene_id,
                    overlap_len,
                });
            },
        );
    }

    let assignment = overlap::resolve_assignment(&counter.hit_buffer, args);
    counter.apply_assignment(assignment, nh, args, overlap_rejected);
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
    let mode = args.strand_mode();
    let is_unstranded = mode == StrandMode::Unstranded;

    counter.hit_buffer.clear();

    // Precompute the expected feature strand for the record and the mate.
    // In unstranded mode these are unused, so set them to Unknown — the
    // callback gates the precomputed check on `!is_unstranded` anyway.
    let (expected_record, expected_mate) = if is_unstranded {
        (
            crate::annotation::Strand::Unknown,
            crate::annotation::Strand::Unknown,
        )
    } else {
        let record_strand = overlap::strand_from_reverse(record.is_reverse());
        let mate_strand = overlap::strand_from_reverse(mate.is_reverse_strand());
        (
            overlap::apply_strand_mode(record_strand, mode),
            overlap::apply_strand_mode(mate_strand, mode),
        )
    };

    // Hoist out of the per-feature callback — answer is constant per fragment.
    let needs_overlap = args.need_overlap_length();

    // Set when a feature overlapped but fell short of the overlap thresholds, so
    // a fragment that ends up with no hits can be reported as Overlapping_Length
    // rather than NoFeatures.
    let mut overlap_rejected = false;

    // Hits from current record - use callback-based query
    for interval in &record.intervals {
        annotation.query_overlapping(
            chrom_id,
            interval.start,
            interval.end,
            |feat_idx, feature| {
                if !is_unstranded
                    && !overlap::check_strand_paired_precomputed(
                        expected_record,
                        expected_mate,
                        feature.strand,
                    )
                {
                    return;
                }

                let overlap_len = if needs_overlap {
                    let ol = crate::alignment::total_overlap(
                        &record.intervals,
                        feature.start,
                        feature.end,
                    ) + crate::alignment::total_overlap(
                        &mate.intervals,
                        feature.start,
                        feature.end,
                    );
                    if !overlap::check_overlap_thresholds(ol, read_len, feature, args) {
                        overlap_rejected = true;
                        return;
                    }
                    ol
                } else {
                    1
                };

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
        // For dedup of features already counted from the record's intervals,
        // linear-scan the existing hit_buffer instead of (re)building an
        // FxHashSet. Typical RNA-seq fragments have 0–3 hits per mate, where
        // the scan is 0–3 cmps — cheaper than the hash + bucket lookup and
        // avoids the per-fragment clear/extend on counter.seen_features.
        // Snapshot the length BEFORE pushing mate hits so we only scan the
        // record's contributions.
        let initial_hit_count = counter.hit_buffer.len();

        for interval in &mate.intervals {
            annotation.query_overlapping(
                chrom_id,
                interval.start,
                interval.end,
                |feat_idx, feature| {
                    // Skip if already counted from the record side.
                    if counter.hit_buffer[..initial_hit_count]
                        .iter()
                        .any(|h| h.feature_idx == feat_idx)
                    {
                        return;
                    }

                    if !is_unstranded
                        && !overlap::check_strand_paired_precomputed(
                            expected_record,
                            expected_mate,
                            feature.strand,
                        )
                    {
                        return;
                    }

                    let overlap_len = if needs_overlap {
                        let ol = crate::alignment::total_overlap(
                            &record.intervals,
                            feature.start,
                            feature.end,
                        ) + crate::alignment::total_overlap(
                            &mate.intervals,
                            feature.start,
                            feature.end,
                        );
                        if !overlap::check_overlap_thresholds(ol, read_len, feature, args) {
                            overlap_rejected = true;
                            return;
                        }
                        ol
                    } else {
                        1
                    };

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
    counter.apply_assignment(assignment, record.nh.max(mate.nh), args, overlap_rejected);
}

/// Process a leftover unpaired mate as a single-end read.
///
/// Thin adapter over [`count_single_end`] — this used to be a fourth copy of the
/// hit-collection loop *plus* a hand-rolled reimplementation of
/// `ThreadCounter::apply_assignment` that omitted the NH divisor, so orphans got
/// full weight under `--fraction` while properly paired multi-mappers got 1/NH.
fn process_deferred_as_single(
    deferred: &DeferredRead,
    counter: &mut ThreadCounter,
    annotation: &AnnotationIndex,
    args: &Args,
) {
    count_single_end(
        &deferred.intervals,
        deferred.chrom_id,
        deferred.is_reverse_strand(),
        deferred.nh,
        counter,
        args,
        annotation,
    );
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
    let tpf = threads_per_file_for(args.threads, num_files);

    let pb = new_progress_bar(num_files, args.quiet);

    let results: Vec<Result<(Vec<i64>, ReadCounters)>> = args
        .bam_files
        .par_iter()
        .map(|bam_input| {
            let result =
                process_bam_parallel_paired(&bam_input.path, args, annotation, count_size, tpf);
            pb.inc(1);
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
        "Parallel paired: {} assigned, {} unassigned",
        aggregated.assigned,
        aggregated.total_unassigned()
    );

    Ok(CountResult {
        counts_per_sample,
        stats_per_sample,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `--quiet` must hide the progress bar.
    ///
    /// Checking this through redirected stderr doesn't work: indicatif already
    /// self-suppresses when stderr is not a terminal, so the bug only showed
    /// interactively — and for the same reason the non-quiet bar also reports
    /// hidden under a test harness, which is why only the quiet case is
    /// asserted here.
    #[test]
    fn test_progress_bar_hidden_when_quiet() {
        assert!(
            new_progress_bar(4, true).is_hidden(),
            "--quiet must suppress the progress bar"
        );
    }
}
