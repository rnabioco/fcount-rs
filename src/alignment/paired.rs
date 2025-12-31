use rustc_hash::FxHashMap;
use smallvec::SmallVec;

use super::cigar::Interval;
use crate::annotation::Strand;

/// Pending mate information for paired-end read tracking
#[derive(Debug, Clone)]
pub struct PendingMate {
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

impl PendingMate {
    /// Check if this mate is on the reverse strand
    #[inline]
    pub fn is_reverse_strand(&self) -> bool {
        self.flags & 0x10 != 0
    }

    /// Get the strand of this read
    pub fn strand(&self) -> Strand {
        if self.is_reverse_strand() {
            Strand::Reverse
        } else {
            Strand::Forward
        }
    }
}

/// Tracks mates for paired-end read processing
///
/// When processing paired-end reads, we need to wait for both mates
/// before we can properly assign the fragment to features.
#[derive(Debug)]
pub struct MateTracker {
    /// Pending mates indexed by read name hash
    pending: FxHashMap<u64, PendingMate>,
    /// Maximum pending mates before forced eviction
    max_pending: usize,
}

impl MateTracker {
    pub fn new(max_pending: usize) -> Self {
        MateTracker {
            pending: FxHashMap::default(),
            max_pending,
        }
    }

    /// Hash a read name for lookup
    #[inline]
    fn hash_name(name: &[u8]) -> u64 {
        use std::hash::{Hash, Hasher};
        let mut hasher = rustc_hash::FxHasher::default();
        name.hash(&mut hasher);
        hasher.finish()
    }

    /// Add a mate to the tracker, returning the other mate if found
    ///
    /// If the mate's pair is already pending, returns Some(pending_mate)
    /// Otherwise, stores this mate and returns None
    pub fn add_mate(
        &mut self,
        read_name: &[u8],
        chrom_id: u16,
        start: u32,
        intervals: SmallVec<[Interval; 4]>,
        flags: u16,
        mapq: u8,
        nh: u8,
    ) -> Option<PendingMate> {
        let name_hash = Self::hash_name(read_name);

        // Check if mate is already pending
        if let Some(mate) = self.pending.remove(&name_hash) {
            return Some(mate);
        }

        // Evict old entries if too many pending
        if self.pending.len() >= self.max_pending {
            // Simple eviction: clear half the entries
            // In practice, this shouldn't happen often with sorted BAMs
            let to_remove: Vec<_> = self
                .pending
                .keys()
                .take(self.max_pending / 2)
                .copied()
                .collect();
            for key in to_remove {
                self.pending.remove(&key);
            }
        }

        // Store this mate
        self.pending.insert(
            name_hash,
            PendingMate {
                chrom_id,
                start,
                intervals,
                flags,
                mapq,
                nh,
            },
        );

        None
    }

    /// Clear all pending mates
    pub fn clear(&mut self) {
        self.pending.clear();
    }

    /// Get the number of pending mates
    pub fn pending_count(&self) -> usize {
        self.pending.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mate_tracking() {
        let mut tracker = MateTracker::new(1000);

        let intervals1: SmallVec<[Interval; 4]> = smallvec::smallvec![Interval {
            start: 100,
            end: 200
        }];
        let intervals2: SmallVec<[Interval; 4]> = smallvec::smallvec![Interval {
            start: 300,
            end: 400
        }];

        // Add first mate
        let result = tracker.add_mate(b"read1", 0, 100, intervals1.clone(), 0x41, 60, 1);
        assert!(result.is_none());
        assert_eq!(tracker.pending_count(), 1);

        // Add second mate - should return first mate
        let result = tracker.add_mate(b"read1", 0, 300, intervals2.clone(), 0x81, 60, 1);
        assert!(result.is_some());
        assert_eq!(tracker.pending_count(), 0);

        let mate = result.unwrap();
        assert_eq!(mate.start, 100);
        assert_eq!(mate.intervals.len(), 1);
    }

    #[test]
    fn test_different_reads() {
        let mut tracker = MateTracker::new(1000);

        let intervals: SmallVec<[Interval; 4]> = smallvec::smallvec![Interval {
            start: 100,
            end: 200
        }];

        // Add mates from different reads
        tracker.add_mate(b"read1", 0, 100, intervals.clone(), 0x41, 60, 1);
        tracker.add_mate(b"read2", 0, 200, intervals.clone(), 0x41, 60, 1);

        assert_eq!(tracker.pending_count(), 2);
    }
}
