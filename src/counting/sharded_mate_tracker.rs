//! Concurrent mate tracker for parallel paired-end processing.
//!
//! This module provides a thread-safe mate tracker that allows multiple worker
//! threads to track and match read pairs concurrently with minimal contention.
//! Uses DashMap for lock-free concurrent access.

use dashmap::DashMap;
use dashmap::mapref::entry::Entry;
use rustc_hash::FxHasher;
use smallvec::SmallVec;
use std::hash::BuildHasherDefault;

use crate::alignment::Interval;

/// Lightweight pending read info stored during parallel processing
#[derive(Debug, Clone)]
pub struct DeferredRead {
    /// Chromosome ID
    pub chrom_id: u16,
    /// Parsed CIGAR intervals
    pub intervals: SmallVec<[Interval; 8]>,
    /// SAM flags
    pub flags: u16,
    /// NH tag value
    pub nh: u8,
}

impl DeferredRead {
    /// Check if this read is on the reverse strand
    #[inline]
    pub fn is_reverse_strand(&self) -> bool {
        self.flags & 0x10 != 0
    }
}

/// Concurrent map for mate tracking using DashMap.
///
/// DashMap provides fine-grained locking with automatic sharding,
/// significantly reducing contention compared to manual sharding.
pub struct ShardedMateTracker {
    map: DashMap<u64, DeferredRead, BuildHasherDefault<FxHasher>>,
}

impl ShardedMateTracker {
    /// Create a new mate tracker.
    ///
    /// # Arguments
    /// * `num_shards` - Number of shards for the DashMap (minimum 16)
    pub fn new(num_shards: usize) -> Self {
        ShardedMateTracker {
            map: DashMap::with_capacity_and_hasher_and_shard_amount(
                100_000,
                BuildHasherDefault::<FxHasher>::default(),
                num_shards.max(16),
            ),
        }
    }

    /// Insert a read or retrieve its mate, acquiring the shard lock once.
    ///
    /// The `make_read` closure constructs the [`DeferredRead`] lazily — it is
    /// only called when no existing mate was found, avoiding the per-record
    /// `intervals.clone()` when this is the second mate of the pair. Holding
    /// the lock across the closure is acceptable because the construction is
    /// pure data assembly (no I/O, no further locking).
    #[inline]
    pub fn take_or_insert_with<F>(
        &self,
        name_hash: u64,
        make_read: F,
    ) -> Option<DeferredRead>
    where
        F: FnOnce() -> DeferredRead,
    {
        match self.map.entry(name_hash) {
            Entry::Occupied(e) => Some(e.remove()),
            Entry::Vacant(e) => {
                e.insert(make_read());
                None
            }
        }
    }

    /// Drain all remaining pending mates.
    ///
    /// Call this after processing to handle orphan reads.
    pub fn drain_all(&self) -> Vec<(u64, DeferredRead)> {
        let mut result = Vec::with_capacity(self.map.len());
        // Use retain with false to drain all entries
        self.map.retain(|k, v| {
            result.push((*k, v.clone()));
            false
        });
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_read(chrom_id: u16, flags: u16) -> DeferredRead {
        let intervals: SmallVec<[Interval; 8]> = smallvec::smallvec![Interval {
            start: 100,
            end: 200
        }];
        DeferredRead {
            chrom_id,
            intervals,
            flags,
            nh: 1,
        }
    }

    #[test]
    fn test_mate_tracking() {
        let tracker = ShardedMateTracker::new(4);
        let hash = 0xdead_beef;

        // First call stores the read and returns None (no mate yet).
        let result = tracker.take_or_insert_with(hash, || make_read(0, 0x41));
        assert!(result.is_none());

        // Second call with the same hash returns the stored first mate.
        let mate = tracker.take_or_insert_with(hash, || make_read(0, 0x81));
        assert!(mate.is_some());
        assert_eq!(mate.unwrap().flags, 0x41);

        // Pair matched → nothing left pending.
        assert_eq!(tracker.drain_all().len(), 0);
    }

    #[test]
    fn test_different_reads() {
        let tracker = ShardedMateTracker::new(4);

        assert!(tracker.take_or_insert_with(1, || make_read(0, 0x41)).is_none());
        assert!(tracker.take_or_insert_with(2, || make_read(0, 0x41)).is_none());

        // Two distinct hashes both remain pending.
        assert_eq!(tracker.drain_all().len(), 2);
    }

    #[test]
    fn test_concurrent_access() {
        use std::sync::Arc;
        use std::thread;

        let tracker = Arc::new(ShardedMateTracker::new(8));
        let mut handles = vec![];

        for thread_id in 0..4u64 {
            let tracker = Arc::clone(&tracker);
            handles.push(thread::spawn(move || {
                for i in 0..100u64 {
                    let hash = thread_id * 1000 + i;
                    tracker.take_or_insert_with(hash, || make_read(0, 0x41));
                }
            }));
        }

        for handle in handles {
            handle.join().unwrap();
        }

        // All 400 distinct hashes should be pending (no mates matched).
        assert_eq!(tracker.drain_all().len(), 400);
    }
}
