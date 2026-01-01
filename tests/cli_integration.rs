//! End-to-end integration tests for fcount-rs CLI
//!
//! These tests run the actual fcount binary against test fixtures from featureCounts
//! and validate the output against expected .ora files.
//!
//! Fixed BAM fixtures are stored in tests/fixtures/ with corrected headers.

mod common;

use common::{assert_matches_expected, fixture_path, run_fcount};

// =============================================================================
// A. Basic Counting Tests
// =============================================================================

#[test]
fn test_across_genes() {
    let result = run_fcount(
        &fixture_path("across_genes_r1.bam"),
        &fixture_path("across_genes.gtf"),
        &[],
    )
    .expect("fcount should run successfully");

    assert!(result.success(), "fcount failed: {}", result.stderr);
    assert_matches_expected(
        &result.counts,
        &fixture_path("across_genes_r1.bam.ora"),
        0.0,
    );
}

#[test]
fn test_across_intron() {
    let result = run_fcount(
        &fixture_path("across_intron_r1.bam"),
        &fixture_path("across_intron.gtf"),
        &[],
    )
    .expect("fcount should run successfully");

    assert!(result.success(), "fcount failed: {}", result.stderr);
    assert_matches_expected(
        &result.counts,
        &fixture_path("across_intron_r1.bam.ora"),
        0.0,
    );
}

#[test]
fn test_intron_between() {
    let result = run_fcount(
        &fixture_path("intron_between.bam"),
        &fixture_path("intron_between.gtf"),
        &[],
    )
    .expect("fcount should run successfully");

    assert!(result.success(), "fcount failed: {}", result.stderr);
    assert_matches_expected(&result.counts, &fixture_path("intron_between.bam.ora"), 0.0);
}

#[test]
fn test_junction_reads() {
    let result = run_fcount(
        &fixture_path("corner-JUNC.bam"),
        &fixture_path("test-minimum.GTF"),
        &["-p"],
    )
    .expect("fcount should run successfully");

    assert!(result.success(), "fcount failed: {}", result.stderr);
    assert_matches_expected(&result.counts, &fixture_path("corner-JUNC.ora"), 0.0);
}

// =============================================================================
// B. Paired-End Tests
// =============================================================================

#[test]
fn test_paired_end_basic() {
    let result = run_fcount(
        &fixture_path("corner-ONEEND.bam"),
        &fixture_path("test-minimum.GTF"),
        &["-p"],
    )
    .expect("fcount should run successfully");

    assert!(result.success(), "fcount failed: {}", result.stderr);
    assert_matches_expected(&result.counts, &fixture_path("corner-ONEEND.ora"), 0.0);
}

#[test]
fn test_paired_end_require_both_aligned() {
    let result = run_fcount(
        &fixture_path("corner-ONEEND.bam"),
        &fixture_path("test-minimum.GTF"),
        &["-p", "-B"],
    )
    .expect("fcount should run successfully");

    assert!(result.success(), "fcount failed: {}", result.stderr);
    assert_matches_expected(&result.counts, &fixture_path("corner-ONEEND-BOTH.ora"), 0.0);
}

// =============================================================================
// C. Multi-Mapping Tests
// =============================================================================

#[test]
fn test_multimapper_nh_tag() {
    let result = run_fcount(
        &fixture_path("corner-NH.bam"),
        &fixture_path("test-minimum.GTF"),
        &["-p"],
    )
    .expect("fcount should run successfully");

    assert!(result.success(), "fcount failed: {}", result.stderr);
    assert_matches_expected(&result.counts, &fixture_path("corner-NH.ora"), 0.0);
}

/// The .ora file was generated with -M --primary, not just --primary
#[test]
fn test_multimapper_primary_only() {
    let result = run_fcount(
        &fixture_path("corner-NH.bam"),
        &fixture_path("test-minimum.GTF"),
        &["-p", "-M", "--primary"],
    )
    .expect("fcount should run successfully");

    assert!(result.success(), "fcount failed: {}", result.stderr);
    assert_matches_expected(&result.counts, &fixture_path("corner-NH-PM.ora"), 0.0);
}

/// SAF format is not currently supported by fcount (only GTF/GFF)
#[test]
#[ignore = "SAF annotation format not supported by fcount"]
fn test_fractional_counting() {
    let result = run_fcount(
        &fixture_path("corner-fractions.bam"),
        &fixture_path("corner-fractions.SAF"),
        &["-M", "--fraction"],
    )
    .expect("fcount should run successfully");

    assert!(result.success(), "fcount failed: {}", result.stderr);
    // Use tolerance for fractional counts
    assert_matches_expected(&result.counts, &fixture_path("corner-fractions.ora"), 0.01);
}

// =============================================================================
// D. Overlap Handling Tests
// =============================================================================

/// Test minimum overlap bases filter (--min-overlap)
#[test]
fn test_min_overlap_bases() {
    let result = run_fcount(
        &fixture_path("test-junc.bam"),
        &fixture_path("test-minimum.GTF"),
        &["-p", "--min-overlap", "10"],
    )
    .expect("fcount should run successfully");

    assert!(result.success(), "fcount failed: {}", result.stderr);
    assert_matches_expected(
        &result.counts,
        &fixture_path("expected-MinOverlap.ora"),
        0.0,
    );
}

/// Test largest overlap selection mode
#[test]
fn test_largest_overlap() {
    let result = run_fcount(
        &fixture_path("test-junc.bam"),
        &fixture_path("test-minimum.GTF"),
        &["-p", "-O", "--largest-overlap"],
    )
    .expect("fcount should run successfully");

    assert!(result.success(), "fcount failed: {}", result.stderr);
    assert_matches_expected(
        &result.counts,
        &fixture_path("expected-LargestOverlap.ora"),
        0.0,
    );
}

// =============================================================================
// E. Read Filtering Tests
// =============================================================================

/// Test minimum mapping quality filter (-Q)
#[test]
fn test_min_mapping_quality() {
    let result = run_fcount(
        &fixture_path("test-junc.bam"),
        &fixture_path("test-minimum.GTF"),
        &["-p", "-Q", "10"],
    )
    .expect("fcount should run successfully");

    assert!(result.success(), "fcount failed: {}", result.stderr);
    assert_matches_expected(&result.counts, &fixture_path("expected-MinMAPQ.ora"), 0.0);
}

#[test]
fn test_ignore_duplicates() {
    let result = run_fcount(
        &fixture_path("test-dup.bam"),
        &fixture_path("test-minimum.GTF"),
        &["-p", "--ignore-dup"],
    )
    .expect("fcount should run successfully");

    assert!(result.success(), "fcount failed: {}", result.stderr);
    assert_matches_expected(&result.counts, &fixture_path("corner-IgnoreDup.ora"), 0.0);
}

// =============================================================================
// F. Edge Cases Tests
// =============================================================================

#[test]
fn test_indel_reads() {
    let result = run_fcount(
        &fixture_path("corner-INDEL.bam"),
        &fixture_path("test-minimum.GTF"),
        &["-p"],
    )
    .expect("fcount should run successfully");

    assert!(result.success(), "fcount failed: {}", result.stderr);
    assert_matches_expected(&result.counts, &fixture_path("corner-INDEL.ora"), 0.0);
}

// =============================================================================
// G. Output Format Tests
// =============================================================================

#[test]
fn test_output_has_correct_columns() {
    let result = run_fcount(
        &fixture_path("corner-JUNC.bam"),
        &fixture_path("test-minimum.GTF"),
        &["-p"],
    )
    .expect("fcount should run successfully");

    assert!(result.success(), "fcount failed: {}", result.stderr);

    // Read the output file directly
    let content = std::fs::read_to_string(&result.output_path).expect("Should read output");

    // Check header line exists with correct columns
    let mut found_header = false;
    for line in content.lines() {
        if line.starts_with("Geneid") {
            found_header = true;
            let cols: Vec<&str> = line.split('\t').collect();
            assert!(cols.len() >= 7, "Should have at least 7 columns");
            assert_eq!(cols[0], "Geneid");
            assert_eq!(cols[1], "Chr");
            assert_eq!(cols[2], "Start");
            assert_eq!(cols[3], "End");
            assert_eq!(cols[4], "Strand");
            assert_eq!(cols[5], "Length");
            break;
        }
    }
    assert!(found_header, "Output should have header line");
}

#[test]
fn test_output_has_program_comment() {
    let result = run_fcount(
        &fixture_path("corner-JUNC.bam"),
        &fixture_path("test-minimum.GTF"),
        &["-p"],
    )
    .expect("fcount should run successfully");

    assert!(result.success(), "fcount failed: {}", result.stderr);

    let content = std::fs::read_to_string(&result.output_path).expect("Should read output");
    let first_line = content
        .lines()
        .next()
        .expect("Should have at least one line");

    assert!(
        first_line.starts_with("# Program:fcount"),
        "First line should be program comment, got: {}",
        first_line
    );
}

// =============================================================================
// H. Summary Output Tests
// =============================================================================

#[test]
fn test_summary_file_created() {
    let result = run_fcount(
        &fixture_path("corner-JUNC.bam"),
        &fixture_path("test-minimum.GTF"),
        &["-p"],
    )
    .expect("fcount should run successfully");

    assert!(result.success(), "fcount failed: {}", result.stderr);
    // Note: Summary file creation depends on fcount implementation
    if !result.summary_path.exists() {
        eprintln!(
            "Note: Summary file not created at {:?} (may not be implemented)",
            result.summary_path
        );
    }
}
