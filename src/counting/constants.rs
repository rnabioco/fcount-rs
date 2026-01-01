//! Shared constants for the counting module.

/// Fixed-point multiplier for fractional counts.
/// Using 1_000_000 to maintain 6 decimal places of precision.
pub const FRACTION_MULTIPLIER: i64 = 1_000_000;

// ============================================================================
// Parallel Processing Configuration
// ============================================================================

/// Maximum threads to allocate per BAM file.
/// More threads show diminishing returns due to I/O saturation.
pub const MAX_THREADS_PER_FILE: usize = 4;

/// Multiplier for channel buffer size relative to worker count.
/// Buffer size = num_workers * CHANNEL_BUFFER_MULTIPLIER
pub const CHANNEL_BUFFER_MULTIPLIER: usize = 4;

/// Multiplier for mate tracker shards relative to worker count.
/// Shards = num_workers * MATE_TRACKER_SHARD_MULTIPLIER
pub const MATE_TRACKER_SHARD_MULTIPLIER: usize = 8;

/// Calculate the number of threads to use per file.
#[inline]
pub fn threads_per_file(total_threads: usize) -> usize {
    MAX_THREADS_PER_FILE.min(total_threads).max(1)
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
