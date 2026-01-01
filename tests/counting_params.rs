//! Integration tests for CLI parameters affecting counting behavior.
//!
//! These tests verify that CLI parameters correctly control the counting logic.

use fcount_rs::annotation::{Feature, Strand};
use fcount_rs::cli::{Args, StrandMode};
use fcount_rs::counting::overlap::{self, FeatureHit};
use std::path::PathBuf;

/// Create a minimal Args struct with defaults for testing
fn default_args() -> Args {
    Args {
        annotation: PathBuf::from("test.gtf"),
        output: PathBuf::from("test.txt"),
        bam_files: vec![PathBuf::from("test.bam")],
        feature_type: "exon".to_string(),
        gene_id_attr: "gene_id".to_string(),
        feature_level: false,
        min_mapping_quality: 0,
        primary_only: false,
        ignore_duplicates: false,
        paired_end: false,
        require_both_aligned: false,
        strand_mode: 0,
        count_multi_mapping: false,
        fractional_counting: false,
        allow_multi_overlap: false,
        min_overlap_bases: 1,
        min_overlap_fraction: 0.0,
        min_feature_overlap_fraction: 0.0,
        largest_overlap_only: false,
        threads: 1,
        details_file: None,
        quiet: false,
        no_chimeric: false,
    }
}

// =============================================================================
// Strand Mode Tests
// =============================================================================

#[test]
fn test_strand_mode_parsing() {
    let mut args = default_args();

    args.strand_mode = 0;
    assert_eq!(args.strand_mode(), StrandMode::Unstranded);

    args.strand_mode = 1;
    assert_eq!(args.strand_mode(), StrandMode::Stranded);

    args.strand_mode = 2;
    assert_eq!(args.strand_mode(), StrandMode::ReverselyStranded);

    // Invalid values should default to unstranded
    args.strand_mode = 99;
    assert_eq!(args.strand_mode(), StrandMode::Unstranded);
}

#[test]
fn test_strand_check_unstranded() {
    let args = default_args(); // strand_mode = 0

    let forward_feature = Feature {
        gene_id: 0,
        chrom_id: 0,
        start: 100,
        end: 200,
        strand: Strand::Forward,
    };

    let reverse_feature = Feature {
        gene_id: 1,
        chrom_id: 0,
        start: 100,
        end: 200,
        strand: Strand::Reverse,
    };

    // In unstranded mode, both should pass regardless of read strand
    assert!(overlap::check_strand_with_strand(
        Strand::Forward,
        &forward_feature,
        &args
    ));
    assert!(overlap::check_strand_with_strand(
        Strand::Reverse,
        &forward_feature,
        &args
    ));
    assert!(overlap::check_strand_with_strand(
        Strand::Forward,
        &reverse_feature,
        &args
    ));
    assert!(overlap::check_strand_with_strand(
        Strand::Reverse,
        &reverse_feature,
        &args
    ));
}

#[test]
fn test_strand_check_stranded() {
    let mut args = default_args();
    args.strand_mode = 1; // Stranded

    let forward_feature = Feature {
        gene_id: 0,
        chrom_id: 0,
        start: 100,
        end: 200,
        strand: Strand::Forward,
    };

    let reverse_feature = Feature {
        gene_id: 1,
        chrom_id: 0,
        start: 100,
        end: 200,
        strand: Strand::Reverse,
    };

    // In stranded mode, read strand must match feature strand
    assert!(overlap::check_strand_with_strand(
        Strand::Forward,
        &forward_feature,
        &args
    ));
    assert!(!overlap::check_strand_with_strand(
        Strand::Reverse,
        &forward_feature,
        &args
    ));
    assert!(!overlap::check_strand_with_strand(
        Strand::Forward,
        &reverse_feature,
        &args
    ));
    assert!(overlap::check_strand_with_strand(
        Strand::Reverse,
        &reverse_feature,
        &args
    ));
}

#[test]
fn test_strand_check_reversely_stranded() {
    let mut args = default_args();
    args.strand_mode = 2; // Reversely stranded

    let forward_feature = Feature {
        gene_id: 0,
        chrom_id: 0,
        start: 100,
        end: 200,
        strand: Strand::Forward,
    };

    let reverse_feature = Feature {
        gene_id: 1,
        chrom_id: 0,
        start: 100,
        end: 200,
        strand: Strand::Reverse,
    };

    // In reversely stranded mode, read strand must be opposite to feature strand
    assert!(!overlap::check_strand_with_strand(
        Strand::Forward,
        &forward_feature,
        &args
    ));
    assert!(overlap::check_strand_with_strand(
        Strand::Reverse,
        &forward_feature,
        &args
    ));
    assert!(overlap::check_strand_with_strand(
        Strand::Forward,
        &reverse_feature,
        &args
    ));
    assert!(!overlap::check_strand_with_strand(
        Strand::Reverse,
        &reverse_feature,
        &args
    ));
}

// =============================================================================
// Overlap Resolution Tests
// =============================================================================

#[test]
fn test_no_hits_returns_no_feature() {
    let args = default_args();
    let hits: Vec<FeatureHit> = vec![];

    let assignment = overlap::resolve_assignment(&hits, &args);
    assert!(matches!(assignment, overlap::Assignment::NoFeature));
}

#[test]
fn test_single_hit_returns_unique() {
    let args = default_args();
    let hits = vec![FeatureHit {
        feature_idx: 0,
        gene_id: 0,
        overlap_len: 100,
    }];

    let assignment = overlap::resolve_assignment(&hits, &args);
    assert!(matches!(assignment, overlap::Assignment::Unique(_)));
}

#[test]
fn test_multiple_hits_same_gene_returns_unique() {
    let args = default_args();
    // Two features from the same gene
    let hits = vec![
        FeatureHit {
            feature_idx: 0,
            gene_id: 0,
            overlap_len: 50,
        },
        FeatureHit {
            feature_idx: 1,
            gene_id: 0,
            overlap_len: 50,
        },
    ];

    let assignment = overlap::resolve_assignment(&hits, &args);
    assert!(matches!(assignment, overlap::Assignment::Unique(_)));
}

#[test]
fn test_multiple_genes_without_multi_overlap_returns_ambiguous() {
    let args = default_args(); // allow_multi_overlap = false
                               // Two features from different genes
    let hits = vec![
        FeatureHit {
            feature_idx: 0,
            gene_id: 0,
            overlap_len: 50,
        },
        FeatureHit {
            feature_idx: 1,
            gene_id: 1,
            overlap_len: 50,
        },
    ];

    let assignment = overlap::resolve_assignment(&hits, &args);
    assert!(matches!(assignment, overlap::Assignment::Ambiguous));
}

#[test]
fn test_multiple_genes_with_multi_overlap_returns_multi() {
    let mut args = default_args();
    args.allow_multi_overlap = true;

    let hits = vec![
        FeatureHit {
            feature_idx: 0,
            gene_id: 0,
            overlap_len: 50,
        },
        FeatureHit {
            feature_idx: 1,
            gene_id: 1,
            overlap_len: 50,
        },
    ];

    let assignment = overlap::resolve_assignment(&hits, &args);
    assert!(matches!(assignment, overlap::Assignment::MultiOverlap(_)));

    if let overlap::Assignment::MultiOverlap(multi_hits) = assignment {
        assert_eq!(multi_hits.len(), 2);
    }
}

#[test]
fn test_largest_overlap_only() {
    let mut args = default_args();
    args.allow_multi_overlap = true;
    args.largest_overlap_only = true;

    let hits = vec![
        FeatureHit {
            feature_idx: 0,
            gene_id: 0,
            overlap_len: 30,
        },
        FeatureHit {
            feature_idx: 1,
            gene_id: 1,
            overlap_len: 70,
        },
    ];

    let assignment = overlap::resolve_assignment(&hits, &args);

    // Should return only the feature with largest overlap
    if let overlap::Assignment::Unique(hit) = assignment {
        assert_eq!(hit.gene_id, 1);
        assert_eq!(hit.overlap_len, 70);
    } else {
        panic!("Expected Unique assignment with largest overlap");
    }
}

#[test]
fn test_largest_overlap_tie_returns_multi() {
    let mut args = default_args();
    args.allow_multi_overlap = true;
    args.largest_overlap_only = true;

    let hits = vec![
        FeatureHit {
            feature_idx: 0,
            gene_id: 0,
            overlap_len: 50,
        },
        FeatureHit {
            feature_idx: 1,
            gene_id: 1,
            overlap_len: 50,
        },
    ];

    let assignment = overlap::resolve_assignment(&hits, &args);

    // Ties should return MultiOverlap
    assert!(matches!(assignment, overlap::Assignment::MultiOverlap(_)));
}

// =============================================================================
// Need Overlap Length Tests
// =============================================================================

#[test]
fn test_need_overlap_length_default() {
    let args = default_args();
    // With default min_overlap_bases=1 and no fraction requirements
    assert!(!args.need_overlap_length());
}

#[test]
fn test_need_overlap_length_with_min_bases() {
    let mut args = default_args();
    args.min_overlap_bases = 10;
    assert!(args.need_overlap_length());
}

#[test]
fn test_need_overlap_length_with_frac_overlap() {
    let mut args = default_args();
    args.min_overlap_fraction = 0.5;
    assert!(args.need_overlap_length());
}

#[test]
fn test_need_overlap_length_with_frac_feature() {
    let mut args = default_args();
    args.min_feature_overlap_fraction = 0.5;
    assert!(args.need_overlap_length());
}

#[test]
fn test_need_overlap_length_with_largest_overlap() {
    let mut args = default_args();
    args.largest_overlap_only = true;
    assert!(args.need_overlap_length());
}

// =============================================================================
// Overlap Threshold Tests
// =============================================================================

#[test]
fn test_min_overlap_bases_filter() {
    use fcount_rs::alignment::Interval;
    use smallvec::SmallVec;

    let mut args = default_args();
    args.min_overlap_bases = 50;

    let feature = Feature {
        gene_id: 0,
        chrom_id: 0,
        start: 100,
        end: 200,
        strand: Strand::Forward,
    };

    let intervals: SmallVec<[Interval; 4]> = smallvec::smallvec![Interval {
        start: 100,
        end: 130
    }];

    // 30 bp overlap, less than 50 required
    assert!(!overlap::check_overlap_thresholds(
        30, &intervals, &feature, &args
    ));

    // 50 bp overlap, exactly enough
    assert!(overlap::check_overlap_thresholds(
        50, &intervals, &feature, &args
    ));

    // 100 bp overlap, more than enough
    assert!(overlap::check_overlap_thresholds(
        100, &intervals, &feature, &args
    ));
}

#[test]
fn test_min_overlap_fraction_filter() {
    use fcount_rs::alignment::Interval;
    use smallvec::SmallVec;

    let mut args = default_args();
    args.min_overlap_fraction = 0.5; // Require 50% of read to overlap

    let feature = Feature {
        gene_id: 0,
        chrom_id: 0,
        start: 100,
        end: 200,
        strand: Strand::Forward,
    };

    // 100bp read
    let intervals: SmallVec<[Interval; 4]> = smallvec::smallvec![Interval {
        start: 100,
        end: 199
    }];

    // 40bp overlap = 40% of 100bp read - should fail
    assert!(!overlap::check_overlap_thresholds(
        40, &intervals, &feature, &args
    ));

    // 50bp overlap = 50% of 100bp read - should pass
    assert!(overlap::check_overlap_thresholds(
        50, &intervals, &feature, &args
    ));

    // 80bp overlap = 80% of 100bp read - should pass
    assert!(overlap::check_overlap_thresholds(
        80, &intervals, &feature, &args
    ));
}

#[test]
fn test_min_feature_overlap_fraction_filter() {
    use fcount_rs::alignment::Interval;
    use smallvec::SmallVec;

    let mut args = default_args();
    args.min_feature_overlap_fraction = 0.5; // Require 50% of feature to be covered

    // 100bp feature
    let feature = Feature {
        gene_id: 0,
        chrom_id: 0,
        start: 100,
        end: 199,
        strand: Strand::Forward,
    };

    let intervals: SmallVec<[Interval; 4]> = smallvec::smallvec![Interval {
        start: 100,
        end: 149
    }];

    // 40bp overlap = 40% of 100bp feature - should fail
    assert!(!overlap::check_overlap_thresholds(
        40, &intervals, &feature, &args
    ));

    // 50bp overlap = 50% of 100bp feature - should pass
    assert!(overlap::check_overlap_thresholds(
        50, &intervals, &feature, &args
    ));
}
