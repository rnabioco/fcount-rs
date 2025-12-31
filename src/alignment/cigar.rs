use noodles_sam::alignment::record::cigar::op::Kind;
use smallvec::SmallVec;

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
    pub fn len(&self) -> u32 {
        if self.end >= self.start {
            self.end - self.start + 1
        } else {
            0
        }
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.end < self.start
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

/// Parse CIGAR operations into genomic intervals
///
/// This function extracts the aligned reference intervals from a CIGAR string,
/// handling matches (M/=/X), deletions (D), and skipped regions (N).
///
/// # Arguments
/// * `cigar` - Iterator over CIGAR operations from noodles
/// * `start_pos` - 1-based start position of the alignment
/// * `out` - Output buffer for intervals (will be cleared first)
pub fn parse_cigar_intervals<'a, I>(cigar: I, start_pos: u32, out: &mut SmallVec<[Interval; 4]>)
where
    I: Iterator<Item = std::io::Result<noodles_sam::alignment::record::cigar::Op>>,
{
    out.clear();

    let mut ref_pos = start_pos;
    let mut current_interval: Option<Interval> = None;

    for op_result in cigar {
        let op = match op_result {
            Ok(op) => op,
            Err(_) => continue,
        };

        let len = op.len() as u32;

        match op.kind() {
            // Operations that consume reference and are part of alignment
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                let interval_end = ref_pos + len - 1;

                match current_interval.as_mut() {
                    Some(interval) => {
                        // Extend current interval
                        interval.end = interval_end;
                    }
                    None => {
                        // Start new interval
                        current_interval = Some(Interval {
                            start: ref_pos,
                            end: interval_end,
                        });
                    }
                }
                ref_pos += len;
            }

            // Deletion: consumes reference but extends current interval
            Kind::Deletion => {
                let interval_end = ref_pos + len - 1;
                match current_interval.as_mut() {
                    Some(interval) => {
                        interval.end = interval_end;
                    }
                    None => {
                        // Unusual: deletion without preceding match
                        // Start new interval anyway
                        current_interval = Some(Interval {
                            start: ref_pos,
                            end: interval_end,
                        });
                    }
                }
                ref_pos += len;
            }

            // Skip (intron): ends current interval, jumps reference position
            Kind::Skip => {
                // Save current interval if exists
                if let Some(interval) = current_interval.take() {
                    if !interval.is_empty() {
                        out.push(interval);
                    }
                }
                ref_pos += len;
            }

            // Operations that only consume read (soft clip, insertion)
            Kind::Insertion | Kind::SoftClip => {
                // Don't advance reference position
            }

            // Hard clip and padding don't affect coordinates
            Kind::HardClip | Kind::Pad => {}
        }
    }

    // Save final interval
    if let Some(interval) = current_interval {
        if !interval.is_empty() {
            out.push(interval);
        }
    }
}

/// Calculate total reference coverage from intervals
pub fn total_coverage(intervals: &[Interval]) -> u32 {
    intervals.iter().map(|i| i.len()).sum()
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
    use noodles_sam::alignment::record::cigar::Op;

    fn make_ops(ops: &[(Kind, usize)]) -> Vec<std::io::Result<Op>> {
        ops.iter()
            .map(|(kind, len)| Ok(Op::new(*kind, *len)))
            .collect()
    }

    #[test]
    fn test_simple_match() {
        let ops = make_ops(&[(Kind::Match, 100)]);
        let mut out = SmallVec::new();
        parse_cigar_intervals(ops.into_iter(), 1000, &mut out);

        assert_eq!(out.len(), 1);
        assert_eq!(out[0].start, 1000);
        assert_eq!(out[0].end, 1099);
    }

    #[test]
    fn test_with_intron() {
        // 50M100N50M
        let ops = make_ops(&[(Kind::Match, 50), (Kind::Skip, 100), (Kind::Match, 50)]);
        let mut out = SmallVec::new();
        parse_cigar_intervals(ops.into_iter(), 1000, &mut out);

        assert_eq!(out.len(), 2);
        assert_eq!(out[0].start, 1000);
        assert_eq!(out[0].end, 1049);
        assert_eq!(out[1].start, 1150);
        assert_eq!(out[1].end, 1199);
    }

    #[test]
    fn test_with_insertion() {
        // 50M5I50M
        let ops = make_ops(&[(Kind::Match, 50), (Kind::Insertion, 5), (Kind::Match, 50)]);
        let mut out = SmallVec::new();
        parse_cigar_intervals(ops.into_iter(), 1000, &mut out);

        assert_eq!(out.len(), 1);
        assert_eq!(out[0].start, 1000);
        assert_eq!(out[0].end, 1099); // Insertion doesn't consume reference
    }

    #[test]
    fn test_with_deletion() {
        // 50M5D50M
        let ops = make_ops(&[(Kind::Match, 50), (Kind::Deletion, 5), (Kind::Match, 50)]);
        let mut out = SmallVec::new();
        parse_cigar_intervals(ops.into_iter(), 1000, &mut out);

        assert_eq!(out.len(), 1);
        assert_eq!(out[0].start, 1000);
        assert_eq!(out[0].end, 1104); // Deletion consumes 5 extra reference bases
    }

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
