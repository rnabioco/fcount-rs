//! Sharded concurrent mate tracker for parallel paired-end processing.
//!
//! This module provides a thread-safe mate tracker that allows multiple worker
//! threads to track and match read pairs concurrently with minimal contention.

use parking_lot::Mutex;
use rustc_hash::FxHashMap;
use smallvec::SmallVec;
use std::hash::{Hash, Hasher};

use crate::alignment::Interval;

/// Lightweight pending read info stored during parallel processing
#[derive(Debug, Clone)]
pub struct DeferredRead {
    /// Chromosome ID
    pub chrom_id: u16,
    /// 1-based start position
    pub start: u32,
    /// Parsed CIGAR intervals
    pub intervals: SmallVec<[Interval; 4]>,
    /// SAM flags
    pub flags: u16,
    /// Mapping quality
    pub mapq: u8,
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

/// Sharded concurrent map for mate tracking.
///
/// Uses multiple shards (each with its own mutex) to reduce lock contention
/// when multiple threads are inserting/retrieving mates concurrently.
pub struct ShardedMateTracker {
    shards: Vec<Mutex<FxHashMap<u64, DeferredRead>>>,
    num_shards: usize,
}

impl ShardedMateTracker {
    /// Create a new sharded mate tracker.
    ///
    /// # Arguments
    /// * `num_shards` - Number of shards (recommended: num_threads * 4)
    pub fn new(num_shards: usize) -> Self {
        let shards = (0..num_shards)
            .map(|_| Mutex::new(FxHashMap::default()))
            .collect();
        ShardedMateTracker { shards, num_shards }
    }

    /// Hash a read name for lookup
    #[inline]
    pub fn hash_name(name: &[u8]) -> u64 {
        let mut hasher = rustc_hash::FxHasher::default();
        name.hash(&mut hasher);
        hasher.finish()
    }

    /// Get the shard index for a given name hash
    #[inline]
    fn get_shard(&self, name_hash: u64) -> usize {
        (name_hash as usize) % self.num_shards
    }

    /// Insert a read or retrieve its mate.
    ///
    /// If the mate is already pending, removes and returns it.
    /// Otherwise, stores this read and returns None.
    ///
    /// # Arguments
    /// * `name_hash` - Hash of the read name (use `hash_name()`)
    /// * `read` - The deferred read information
    ///
    /// # Returns
    /// * `Some(mate)` if the mate was found and removed
    /// * `None` if this read was stored (awaiting mate)
    #[inline]
    pub fn insert_or_get(&self, name_hash: u64, read: DeferredRead) -> Option<DeferredRead> {
        let shard_idx = self.get_shard(name_hash);
        let mut shard = self.shards[shard_idx].lock();

        if let Some(mate) = shard.remove(&name_hash) {
            Some(mate)
        } else {
            shard.insert(name_hash, read);
            None
        }
    }

    /// Get the total number of pending mates across all shards.
    ///
    /// Note: This is approximate as shards can change during iteration.
    pub fn pending_count(&self) -> usize {
        self.shards.iter().map(|s| s.lock().len()).sum()
    }

    /// Drain all remaining pending mates.
    ///
    /// Call this after processing to handle orphan reads.
    pub fn drain_all(&self) -> Vec<(u64, DeferredRead)> {
        let mut result = Vec::new();
        for shard in &self.shards {
            let mut guard = shard.lock();
            result.extend(guard.drain());
        }
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sharded_mate_tracking() {
        let tracker = ShardedMateTracker::new(4);

        let intervals: SmallVec<[Interval; 4]> = smallvec::smallvec![Interval {
            start: 100,
            end: 200
        }];

        let read1 = DeferredRead {
            chrom_id: 0,
            start: 100,
            intervals: intervals.clone(),
            flags: 0x41,
            mapq: 60,
            nh: 1,
        };

        let read2 = DeferredRead {
            chrom_id: 0,
            start: 300,
            intervals: intervals.clone(),
            flags: 0x81,
            mapq: 60,
            nh: 1,
        };

        let hash = ShardedMateTracker::hash_name(b"read1");

        // First insert should return None (stored)
        let result = tracker.insert_or_get(hash, read1);
        assert!(result.is_none());
        assert_eq!(tracker.pending_count(), 1);

        // Second insert should return the first read
        let result = tracker.insert_or_get(hash, read2);
        assert!(result.is_some());
        assert_eq!(tracker.pending_count(), 0);

        let mate = result.unwrap();
        assert_eq!(mate.start, 100);
    }

    #[test]
    fn test_different_reads() {
        let tracker = ShardedMateTracker::new(4);

        let intervals: SmallVec<[Interval; 4]> = smallvec::smallvec![Interval {
            start: 100,
            end: 200
        }];

        let read = DeferredRead {
            chrom_id: 0,
            start: 100,
            intervals: intervals.clone(),
            flags: 0x41,
            mapq: 60,
            nh: 1,
        };

        let hash1 = ShardedMateTracker::hash_name(b"read1");
        let hash2 = ShardedMateTracker::hash_name(b"read2");

        tracker.insert_or_get(hash1, read.clone());
        tracker.insert_or_get(hash2, read.clone());

        assert_eq!(tracker.pending_count(), 2);
    }

    #[test]
    fn test_concurrent_access() {
        use std::sync::Arc;
        use std::thread;

        let tracker = Arc::new(ShardedMateTracker::new(8));
        let mut handles = vec![];

        // Spawn multiple threads inserting reads
        for thread_id in 0..4 {
            let tracker = Arc::clone(&tracker);
            handles.push(thread::spawn(move || {
                let intervals: SmallVec<[Interval; 4]> = smallvec::smallvec![Interval {
                    start: 100,
                    end: 200
                }];

                for i in 0..100 {
                    let read_name = format!("thread{}_read{}", thread_id, i);
                    let hash = ShardedMateTracker::hash_name(read_name.as_bytes());

                    let read = DeferredRead {
                        chrom_id: 0,
                        start: i as u32,
                        intervals: intervals.clone(),
                        flags: 0x41,
                        mapq: 60,
                        nh: 1,
                    };

                    tracker.insert_or_get(hash, read);
                }
            }));
        }

        for handle in handles {
            handle.join().unwrap();
        }

        // All reads should be pending (no mates inserted)
        assert_eq!(tracker.pending_count(), 400);
    }
}
