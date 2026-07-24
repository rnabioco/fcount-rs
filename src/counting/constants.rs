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

/// Cap on BGZF decompression threads per BAM file. BGZF decompression with
/// libdeflate saturates well below the record-processing rate, so giving it
/// every thread the user requested just adds OS-scheduler contention without
/// improving throughput. At -t 8 with the old uncapped policy a single file
/// spawned 1 producer + 8 BGZF + 8 record workers = 17 OS threads, which is
/// >2× the user's stated budget. Capping BGZF here keeps the record-worker
/// count as the dominant thread consumer.
pub const BGZF_THREADS_CAP: usize = 4;

/// Multiplier for channel buffer size relative to worker count.
/// Buffer size = num_workers * CHANNEL_BUFFER_MULTIPLIER
pub const CHANNEL_BUFFER_MULTIPLIER: usize = 4;

/// Multiplier for mate tracker shards relative to worker count.
/// Shards = num_workers * MATE_TRACKER_SHARD_MULTIPLIER
pub const MATE_TRACKER_SHARD_MULTIPLIER: usize = 8;

/// Calculate how many parse/count worker threads to spawn per BAM file.
///
/// With a single BAM, give the file the user's *whole* thread budget. Record
/// processing — not BGZF decode — is the bottleneck (decode is separately
/// capped at [`BGZF_THREADS_CAP`]), so on a many-core box a single large BAM
/// should be allowed to scale past the soft cap. The soft cap only applies
/// when several BAMs run concurrently, to avoid oversubscribing the scheduler.
#[inline]
pub fn threads_per_file_for(total_threads: usize, num_files: usize) -> usize {
    let total = total_threads.max(1);
    let files = num_files.max(1);
    if files <= 1 {
        // Single file: use the whole budget. No soft cap — more record workers
        // directly improve throughput, and BGZF decode is capped elsewhere.
        total
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
///
/// DashMap requires the shard count to be a power of two and panics otherwise
/// (`assertion failed: shard_amount.is_power_of_two()`). `num_workers` is
/// derived from `total_threads / num_files`, so it is frequently *not* a power
/// of two (e.g. 3 workers → 24), which would make the raw `num_workers * 8`
/// product an invalid shard count. Round up to the next power of two so every
/// worker count produces a valid map.
#[inline]
pub fn mate_tracker_shards(num_workers: usize) -> usize {
    (num_workers.max(1) * MATE_TRACKER_SHARD_MULTIPLIER).next_power_of_two()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mate_tracker_shards_is_always_power_of_two() {
        // num_workers = total_threads / num_files, so every small integer is a
        // realistic input. 3 (e.g. -t 64 over 18 BAMs) is the case that used to
        // yield 24 and panic DashMap.
        for num_workers in 0..=64 {
            let shards = mate_tracker_shards(num_workers);
            assert!(
                shards.is_power_of_two(),
                "mate_tracker_shards({num_workers}) = {shards} is not a power of two",
            );
        }
    }

    #[test]
    fn threads_per_file_64_over_18_reproduces_the_panic_input() {
        // Regression guard for the exact production case: 64 threads, 18 BAMs.
        let tpf = threads_per_file_for(64, 18);
        assert_eq!(tpf, 3);
        assert!(mate_tracker_shards(tpf).is_power_of_two());
    }
}
