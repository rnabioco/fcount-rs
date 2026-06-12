/// A genomic interval derived from CIGAR parsing
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Interval {
    /// 1-based start position (inclusive)
    pub start: u32,
    /// 1-based end position (inclusive)
    pub end: u32,
}

impl Interval {
    #[inline]
    #[allow(clippy::len_without_is_empty)] // genomic length, not a container
    pub fn len(&self) -> u32 {
        if self.end >= self.start {
            self.end - self.start + 1
        } else {
            0
        }
    }

    #[inline]
    pub fn overlaps(&self, start: u32, end: u32) -> bool {
        self.start <= end && self.end >= start
    }

    #[inline]
    pub fn overlap_len(&self, start: u32, end: u32) -> u32 {
        if !self.overlaps(start, end) {
            return 0;
        }
        let overlap_start = self.start.max(start);
        let overlap_end = self.end.min(end);
        overlap_end - overlap_start + 1
    }
}

/// Calculate total overlap between intervals and a feature
pub fn total_overlap(intervals: &[Interval], feature_start: u32, feature_end: u32) -> u32 {
    intervals
        .iter()
        .map(|i| i.overlap_len(feature_start, feature_end))
        .sum()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_interval_overlap() {
        let interval = Interval {
            start: 100,
            end: 200,
        };

        // Full overlap
        assert_eq!(interval.overlap_len(50, 250), 101);

        // Partial overlap at start
        assert_eq!(interval.overlap_len(50, 150), 51);

        // Partial overlap at end
        assert_eq!(interval.overlap_len(150, 250), 51);

        // No overlap
        assert_eq!(interval.overlap_len(300, 400), 0);

        // Contained
        assert_eq!(interval.overlap_len(120, 180), 61);
    }
}
