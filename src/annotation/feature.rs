/// Strand of a genomic feature
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Strand {
    Forward,
    Reverse,
    Unknown,
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
    #[inline(always)]
    #[allow(clippy::len_without_is_empty)] // genomic length, not a container
    pub fn len(&self) -> u32 {
        self.end - self.start + 1
    }
}
