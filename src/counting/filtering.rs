//! Read filtering logic for feature counting.
//!
//! This module provides a unified filtering interface for all record types,
//! consolidating the filtering logic that was previously duplicated across
//! multiple processing paths.

use crate::cli::Args;

/// Result of filtering a read
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FilterResult {
    /// Read passes all filters and should be processed
    Pass,
    /// Read is unmapped
    Unmapped,
    /// Read is secondary or supplementary alignment
    Secondary,
    /// Read has low mapping quality
    MappingQuality,
    /// Read is a PCR duplicate
    Duplicate,
    /// Read is a multi-mapper (and multi-mapping not allowed)
    MultiMapping,
}

impl FilterResult {
    /// Returns true if the read passes all filters
    #[inline]
    pub fn is_pass(&self) -> bool {
        matches!(self, FilterResult::Pass)
    }
}

/// SAM flags used for filtering
pub mod flags {
    pub const UNMAPPED: u16 = 0x4;
    pub const SECONDARY: u16 = 0x100;
    pub const DUPLICATE: u16 = 0x400;
    pub const SUPPLEMENTARY: u16 = 0x800;
}

/// Trait for records that can be filtered.
///
/// This trait abstracts over the common filtering operations needed for
/// AlignmentRecord, MinimalRecord, DeferredRead, and PendingMate.
pub trait Filterable {
    /// Get the SAM flags
    fn flags(&self) -> u16;

    /// Get the mapping quality
    fn mapq(&self) -> u8;

    /// Get the NH tag value (number of alignments)
    fn nh(&self) -> u8;

    /// Check if the read is unmapped
    #[inline]
    fn is_unmapped(&self) -> bool {
        self.flags() & flags::UNMAPPED != 0
    }

    /// Check if the read is a secondary alignment
    #[inline]
    fn is_secondary(&self) -> bool {
        self.flags() & flags::SECONDARY != 0
    }

    /// Check if the read is a supplementary alignment
    #[inline]
    fn is_supplementary(&self) -> bool {
        self.flags() & flags::SUPPLEMENTARY != 0
    }

    /// Check if the read is a PCR duplicate
    #[inline]
    fn is_duplicate(&self) -> bool {
        self.flags() & flags::DUPLICATE != 0
    }

    /// Check if the read is on the reverse strand
    #[inline]
    fn is_reverse_strand(&self) -> bool {
        self.flags() & 0x10 != 0
    }
}

/// Filter a record based on the provided arguments.
///
/// This is the central filtering function that applies all read-level filters.
/// Returns a FilterResult indicating whether the read passes or why it was filtered.
#[inline]
pub fn filter_record<R: Filterable>(record: &R, args: &Args) -> FilterResult {
    // Skip unmapped reads
    if record.is_unmapped() {
        return FilterResult::Unmapped;
    }

    // Skip secondary/supplementary if primary-only mode
    if args.primary_only && (record.is_secondary() || record.is_supplementary()) {
        return FilterResult::Secondary;
    }

    // Skip low quality
    if record.mapq() < args.min_mapping_quality {
        return FilterResult::MappingQuality;
    }

    // Skip duplicates if requested
    if args.ignore_duplicates && record.is_duplicate() {
        return FilterResult::Duplicate;
    }

    // Skip multi-mappers if not counting them
    if !args.count_multi_mapping && record.nh() > 1 {
        return FilterResult::MultiMapping;
    }

    FilterResult::Pass
}

// Implement Filterable for AlignmentRecord
impl Filterable for crate::alignment::AlignmentRecord {
    #[inline]
    fn flags(&self) -> u16 {
        self.flags
    }

    #[inline]
    fn mapq(&self) -> u8 {
        self.mapq
    }

    #[inline]
    fn nh(&self) -> u8 {
        self.nh
    }
}

// Implement Filterable for MinimalRecord
impl Filterable for crate::alignment::minimal_parser::MinimalRecord {
    #[inline]
    fn flags(&self) -> u16 {
        self.flags
    }

    #[inline]
    fn mapq(&self) -> u8 {
        self.mapq
    }

    #[inline]
    fn nh(&self) -> u8 {
        self.nh
    }
}

// Implement Filterable for DeferredRead
impl Filterable for super::sharded_mate_tracker::DeferredRead {
    #[inline]
    fn flags(&self) -> u16 {
        self.flags
    }

    #[inline]
    fn mapq(&self) -> u8 {
        self.mapq
    }

    #[inline]
    fn nh(&self) -> u8 {
        self.nh
    }
}

// Implement Filterable for PendingMate
impl Filterable for crate::alignment::PendingMate {
    #[inline]
    fn flags(&self) -> u16 {
        self.flags
    }

    #[inline]
    fn mapq(&self) -> u8 {
        self.mapq
    }

    #[inline]
    fn nh(&self) -> u8 {
        self.nh
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    struct MockRecord {
        flags: u16,
        mapq: u8,
        nh: u8,
    }

    impl Filterable for MockRecord {
        fn flags(&self) -> u16 {
            self.flags
        }
        fn mapq(&self) -> u8 {
            self.mapq
        }
        fn nh(&self) -> u8 {
            self.nh
        }
    }

    fn make_args() -> Args {
        Args {
            annotation: std::path::PathBuf::new(),
            output: std::path::PathBuf::new(),
            bam_files: vec![],
            feature_type: "exon".to_string(),
            gene_id_attr: "gene_id".to_string(),
            feature_level: false,
            paired_end: false,
            require_both_aligned: false,
            no_chimeric: false,
            strand_mode: 0,
            count_multi_mapping: false,
            fractional_counting: false,
            primary_only: true,
            allow_multi_overlap: false,
            min_overlap_bases: 1,
            min_overlap_fraction: 0.0,
            min_feature_overlap_fraction: 0.0,
            largest_overlap_only: false,
            min_mapping_quality: 10,
            ignore_duplicates: true,
            threads: 1,
            details_file: None,
            quiet: false,
            output_format: crate::cli::OutputFormat::Featurecounts,
        }
    }

    #[test]
    fn test_pass_all_filters() {
        let record = MockRecord {
            flags: 0,
            mapq: 60,
            nh: 1,
        };
        let args = make_args();
        assert_eq!(filter_record(&record, &args), FilterResult::Pass);
    }

    #[test]
    fn test_filter_unmapped() {
        let record = MockRecord {
            flags: flags::UNMAPPED,
            mapq: 60,
            nh: 1,
        };
        let args = make_args();
        assert_eq!(filter_record(&record, &args), FilterResult::Unmapped);
    }

    #[test]
    fn test_filter_secondary() {
        let record = MockRecord {
            flags: flags::SECONDARY,
            mapq: 60,
            nh: 1,
        };
        let args = make_args();
        assert_eq!(filter_record(&record, &args), FilterResult::Secondary);
    }

    #[test]
    fn test_filter_supplementary() {
        let record = MockRecord {
            flags: flags::SUPPLEMENTARY,
            mapq: 60,
            nh: 1,
        };
        let args = make_args();
        assert_eq!(filter_record(&record, &args), FilterResult::Secondary);
    }

    #[test]
    fn test_filter_low_mapq() {
        let record = MockRecord {
            flags: 0,
            mapq: 5,
            nh: 1,
        };
        let args = make_args();
        assert_eq!(filter_record(&record, &args), FilterResult::MappingQuality);
    }

    #[test]
    fn test_filter_duplicate() {
        let record = MockRecord {
            flags: flags::DUPLICATE,
            mapq: 60,
            nh: 1,
        };
        let args = make_args();
        assert_eq!(filter_record(&record, &args), FilterResult::Duplicate);
    }

    #[test]
    fn test_filter_multimapper() {
        let record = MockRecord {
            flags: 0,
            mapq: 60,
            nh: 2,
        };
        let args = make_args();
        assert_eq!(filter_record(&record, &args), FilterResult::MultiMapping);
    }

    #[test]
    fn test_allow_multimapper() {
        let record = MockRecord {
            flags: 0,
            mapq: 60,
            nh: 2,
        };
        let mut args = make_args();
        args.count_multi_mapping = true;
        assert_eq!(filter_record(&record, &args), FilterResult::Pass);
    }
}
