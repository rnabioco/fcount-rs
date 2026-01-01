/// Strand of a genomic feature
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Strand {
    Forward,
    Reverse,
    Unknown,
}

impl Strand {
    pub fn from_char(c: char) -> Self {
        match c {
            '+' => Strand::Forward,
            '-' => Strand::Reverse,
            _ => Strand::Unknown,
        }
    }

    pub fn from_option(opt: Option<&char>) -> Self {
        match opt {
            Some('+') => Strand::Forward,
            Some('-') => Strand::Reverse,
            _ => Strand::Unknown,
        }
    }
}

/// A genomic feature (exon, gene, etc.)
///
/// Fields are ordered for optimal memory layout (16 bytes instead of 20):
/// - 4-byte aligned fields first (start, end, gene_id)
/// - Then 2-byte field (chrom_id)
/// - Then 1-byte field (strand) with 1 byte padding
#[derive(Debug, Clone)]
pub struct Feature {
    /// 1-based start position (inclusive)
    pub start: u32,
    /// 1-based end position (inclusive)
    pub end: u32,
    /// Index into the gene_names vector
    pub gene_id: u32,
    /// Index into the chromosome name vector
    pub chrom_id: u16,
    /// Feature strand
    pub strand: Strand,
}

impl Feature {
    /// Calculate the length of this feature in bases
    #[inline]
    pub fn len(&self) -> u32 {
        self.end - self.start + 1
    }

    /// Check if the feature is empty (invalid coordinates)
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.end < self.start
    }

    /// Check if this feature overlaps with a given interval
    #[inline]
    pub fn overlaps(&self, start: u32, end: u32) -> bool {
        self.start <= end && self.end >= start
    }

    /// Calculate overlap length with a given interval
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
