//! Concurrent mate tracker for parallel paired-end processing.
//!
//! This module provides a thread-safe mate tracker that allows multiple worker
//! threads to track and match read pairs concurrently with minimal contention.
//! Uses DashMap for lock-free concurrent access.

use dashmap::DashMap;
use rustc_hash::FxHasher;
use smallvec::SmallVec;
use std::hash::{BuildHasherDefault, Hash, Hasher};

use crate::alignment::Interval;

/// Lightweight pending read info stored during parallel processing
#[derive(Debug, Clone)]
pub struct DeferredRead {
    /// Chromosome ID
    pub chrom_id: u16,
    /// 1-based start position
    pub start: u32,
    /// Parsed CIGAR intervals
    pub intervals: SmallVec<[Interval; 8]>,
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

    /// Hash a read name for lookup
    #[inline]
    pub fn hash_name(name: &[u8]) -> u64 {
        let mut hasher = FxHasher::default();
        name.hash(&mut hasher);
        hasher.finish()
    }

    /// Try to remove and return an existing mate.
    ///
    /// If a mate with this hash is pending, removes and returns it.
    /// Otherwise returns None (caller should insert).
    #[inline]
    pub fn remove_mate(&self, name_hash: u64) -> Option<DeferredRead> {
        self.map.remove(&name_hash).map(|(_, v)| v)
    }

    /// Insert a read for later mate matching.
    ///
    /// Call this only after `remove_mate` returns None.
    #[inline]
    pub fn insert(&self, name_hash: u64, read: DeferredRead) {
        self.map.insert(name_hash, read);
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
        // Try to remove existing entry first (common case for second mate)
        if let Some((_, mate)) = self.map.remove(&name_hash) {
            Some(mate)
        } else {
            // Insert new entry (first mate)
            self.map.insert(name_hash, read);
            None
        }
    }

    /// Get the total number of pending mates.
    pub fn pending_count(&self) -> usize {
        self.map.len()
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

    #[test]
    fn test_mate_tracking() {
        let tracker = ShardedMateTracker::new(4);

        let intervals: SmallVec<[Interval; 8]> = smallvec::smallvec![Interval {
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

        let intervals: SmallVec<[Interval; 8]> = smallvec::smallvec![Interval {
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
                let intervals: SmallVec<[Interval; 8]> = smallvec::smallvec![Interval {
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
