//! Shared constants for the counting module.

/// Fixed-point multiplier for fractional counts.
/// Using 1_000_000 to maintain 6 decimal places of precision.
pub const FRACTION_MULTIPLIER: i64 = 1_000_000;

// ============================================================================
// Parallel Processing Configuration
// ============================================================================

/// Soft cap on parse/count worker threads per file. We do *not* enforce this
/// when there's a single BAM (one file should be allowed to use every thread
/// the user gave us). The cap only kicks in when the user passes many BAM
/// files at once, so we don't oversubscribe the scheduler.
pub const SOFT_THREADS_PER_FILE_CAP: usize = 8;

/// Multiplier for channel buffer size relative to worker count.
/// Buffer size = num_workers * CHANNEL_BUFFER_MULTIPLIER
pub const CHANNEL_BUFFER_MULTIPLIER: usize = 4;

/// Multiplier for mate tracker shards relative to worker count.
/// Shards = num_workers * MATE_TRACKER_SHARD_MULTIPLIER
pub const MATE_TRACKER_SHARD_MULTIPLIER: usize = 8;

/// Calculate how many parse/count worker threads to spawn per BAM file.
///
/// With a single BAM, give the file the user's whole thread budget — capping
/// at 4 (the previous behavior) leaves cores idle on bigger machines.
/// With multiple BAMs, split the budget evenly (still ≥1) so we don't
/// oversubscribe.
#[inline]
pub fn threads_per_file(total_threads: usize) -> usize {
    threads_per_file_for(total_threads, 1)
}

/// Variant that knows how many BAM files we're about to process in parallel.
#[inline]
pub fn threads_per_file_for(total_threads: usize, num_files: usize) -> usize {
    let total = total_threads.max(1);
    let files = num_files.max(1);
    if files <= 1 {
        // Single file: use the whole budget, capped by the soft limit so we
        // don't oversubscribe BGZF decode threads on very-many-core boxes.
        total.min(SOFT_THREADS_PER_FILE_CAP).max(1)
    } else {
        // Many files: split the budget across them, but always at least 1.
        (total / files).clamp(1, SOFT_THREADS_PER_FILE_CAP)
    }
}

/// Calculate the channel buffer size for a given number of workers.
#[inline]
pub fn channel_buffer_size(num_workers: usize) -> usize {
    num_workers * CHANNEL_BUFFER_MULTIPLIER
}

/// Calculate the number of mate tracker shards for a given number of workers.
#[inline]
pub fn mate_tracker_shards(num_workers: usize) -> usize {
    num_workers * MATE_TRACKER_SHARD_MULTIPLIER
}
