use rustc_hash::FxHashSet;
use smallvec::SmallVec;

use crate::annotation::{Feature, Strand};
use crate::cli::{Args, StrandMode};

/// A hit to a feature (12 bytes - trivially copyable)
#[derive(Debug, Clone, Copy)]
pub struct FeatureHit {
    /// Index into features array
    pub feature_idx: u32,
    /// Gene ID (for gene-level counting)
    pub gene_id: u32,
    /// Overlap length in bases
    pub overlap_len: u32,
}

/// Assignment result for a read/fragment
#[derive(Debug)]
pub enum Assignment {
    /// Assigned to a single feature/gene
    Unique(FeatureHit),
    /// Overlaps multiple features/genes with equal priority
    MultiOverlap(Vec<FeatureHit>),
    /// No overlapping features found
    NoFeature,
    /// Multiple features with same gene - ambiguous
    Ambiguous,
}

// ============================================================================
// Unified Strand Checking
// ============================================================================

/// Apply strand mode transformation to get the expected feature strand.
///
/// This is the core strand logic used by all strand checking functions.
/// Given a read strand and strand mode, returns the expected feature strand.
#[inline(always)]
pub fn apply_strand_mode(read_strand: Strand, mode: StrandMode) -> Strand {
    match mode {
        StrandMode::Unstranded => Strand::Unknown, // Any strand matches
        StrandMode::Stranded => read_strand,
        StrandMode::ReverselyStranded => match read_strand {
            Strand::Forward => Strand::Reverse,
            Strand::Reverse => Strand::Forward,
            Strand::Unknown => Strand::Unknown,
        },
    }
}

/// Check if a read/fragment strand is compatible with a feature strand.
///
/// This is the unified strand compatibility check used by all strand functions.
#[inline(always)]
fn is_strand_compatible(expected_strand: Strand, feature_strand: Strand) -> bool {
    // Unknown feature strand always matches
    if feature_strand == Strand::Unknown {
        return true;
    }
    // Unstranded mode (Unknown expected) always matches
    if expected_strand == Strand::Unknown {
        return true;
    }
    // Otherwise must match
    expected_strand == feature_strand
}

/// Get the strand from a boolean reverse flag.
#[inline(always)]
pub fn strand_from_reverse(is_reverse: bool) -> Strand {
    if is_reverse {
        Strand::Reverse
    } else {
        Strand::Forward
    }
}

/// Fast strand check with precomputed expected strand.
/// For unstranded mode, pass None to skip the check entirely.
#[inline(always)]
pub fn check_strand_fast(expected: Option<Strand>, feature_strand: Strand) -> bool {
    match expected {
        None => true, // Unstranded - always matches
        Some(exp) => is_strand_compatible(exp, feature_strand),
    }
}

/// Check strand compatibility using Strand enum directly (for parallel processing).
#[inline(always)]
pub fn check_strand_with_strand(read_strand: Strand, feature: &Feature, args: &Args) -> bool {
    let expected = apply_strand_mode(read_strand, args.strand_mode());
    is_strand_compatible(expected, feature.strand)
}

/// Paired-end strand check with precomputed expected strands.
///
/// Both `expected_record` and `expected_mate` come from
/// `apply_strand_mode(read_strand, mode)` evaluated once per fragment.
/// Because `strand_from_reverse` only returns Forward or Reverse, the expected
/// values are never `Strand::Unknown` in the stranded paths — so the only way
/// to accept is feature unknown, or the feature strand matches either expected.
///
/// Caller must skip the call entirely in unstranded mode (the precomputed
/// expecteds are meaningless there).
#[inline(always)]
pub fn check_strand_paired_precomputed(
    expected_record: Strand,
    expected_mate: Strand,
    feature_strand: Strand,
) -> bool {
    feature_strand == Strand::Unknown
        || feature_strand == expected_record
        || feature_strand == expected_mate
}

/// Check if overlap meets threshold requirements.
/// `read_len` should be precomputed once per read (sum of interval lengths).
/// Pass 0 if min_overlap_fraction is not used.
#[inline(always)]
pub fn check_overlap_thresholds(
    overlap_len: u32,
    read_len: u32,
    feature: &Feature,
    args: &Args,
) -> bool {
    // Check minimum absolute overlap
    if overlap_len < args.min_overlap_bases {
        return false;
    }

    // Check minimum fraction of read
    if args.min_overlap_fraction > 0.0 && read_len > 0 {
        let frac = overlap_len as f32 / read_len as f32;
        if frac < args.min_overlap_fraction {
            return false;
        }
    }

    // Check minimum fraction of feature
    if args.min_feature_overlap_fraction > 0.0 {
        let feature_len = feature.len();
        if feature_len > 0 {
            let frac = overlap_len as f32 / feature_len as f32;
            if frac < args.min_feature_overlap_fraction {
                return false;
            }
        }
    }

    true
}

/// Resolve assignment from a list of feature hits
/// Optimized to avoid allocations for common cases (1-4 hits)
pub fn resolve_assignment(hits: &[FeatureHit], args: &Args) -> Assignment {
    if hits.is_empty() {
        return Assignment::NoFeature;
    }

    if hits.len() == 1 {
        return Assignment::Unique(hits[0]); // Copy, not clone
    }

    // Quick check: are all hits from the same gene? (common case)
    let first_gene = hits[0].gene_id;
    let all_same_gene = hits.iter().skip(1).all(|h| h.gene_id == first_gene);

    if all_same_gene {
        // All hits are from the same gene - pick the one with best overlap
        if args.largest_overlap_only {
            let best = hits.iter().max_by_key(|h| h.overlap_len).unwrap();
            return Assignment::Unique(*best); // Copy, not clone
        }
        // Multiple exons of same gene - count as unique to that gene
        return Assignment::Unique(hits[0]); // Copy, not clone
    }

    // Multiple genes
    if args.allow_multi_overlap {
        if args.largest_overlap_only {
            // Find the gene(s) with largest overlap - single pass
            let max_overlap = hits.iter().map(|h| h.overlap_len).max().unwrap();

            // Use FxHashSet for O(1) gene deduplication instead of O(n) contains
            let mut best_hits: SmallVec<[FeatureHit; 8]> = SmallVec::new();
            let mut seen_genes: FxHashSet<u32> = FxHashSet::default();

            for h in hits {
                if h.overlap_len == max_overlap && seen_genes.insert(h.gene_id) {
                    best_hits.push(*h); // Copy
                }
            }

            if best_hits.len() == 1 {
                return Assignment::Unique(best_hits[0]);
            }

            return Assignment::MultiOverlap(best_hits.into_vec());
        }

        // Deduplicate by gene using FxHashSet for O(1) lookup
        let mut seen_genes: FxHashSet<u32> = FxHashSet::default();
        let mut deduped: SmallVec<[FeatureHit; 8]> = SmallVec::new();

        for h in hits {
            if seen_genes.insert(h.gene_id) {
                deduped.push(*h); // Copy
            }
        }

        return Assignment::MultiOverlap(deduped.into_vec());
    }

    // Multiple genes without multi-overlap - ambiguous
    Assignment::Ambiguous
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_args(multi_overlap: bool, largest_only: bool) -> Args {
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
            primary_only: false,
            allow_multi_overlap: multi_overlap,
            min_overlap_bases: 1,
            min_overlap_fraction: 0.0,
            min_feature_overlap_fraction: 0.0,
            largest_overlap_only: largest_only,
            min_mapping_quality: 0,
            ignore_duplicates: false,
            threads: 1,
            details_file: None,
            quiet: false,
            output_format: crate::cli::OutputFormat::Featurecounts,
        }
    }

    #[test]
    fn test_single_hit() {
        let hits = vec![FeatureHit {
            feature_idx: 0,
            gene_id: 0,
            overlap_len: 100,
        }];

        let args = make_args(false, false);
        match resolve_assignment(&hits, &args) {
            Assignment::Unique(h) => assert_eq!(h.gene_id, 0),
            _ => panic!("Expected unique assignment"),
        }
    }

    #[test]
    fn test_same_gene_multiple_exons() {
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

        let args = make_args(false, false);
        match resolve_assignment(&hits, &args) {
            Assignment::Unique(h) => assert_eq!(h.gene_id, 0),
            _ => panic!("Expected unique assignment"),
        }
    }

    #[test]
    fn test_multi_gene_no_overlap_allowed() {
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

        let args = make_args(false, false);
        match resolve_assignment(&hits, &args) {
            Assignment::Ambiguous => {}
            _ => panic!("Expected ambiguous assignment"),
        }
    }

    #[test]
    fn test_multi_gene_with_overlap_allowed() {
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

        let args = make_args(true, false);
        match resolve_assignment(&hits, &args) {
            Assignment::MultiOverlap(h) => assert_eq!(h.len(), 2),
            _ => panic!("Expected multi-overlap assignment"),
        }
    }

    #[test]
    fn test_largest_overlap_only() {
        let hits = vec![
            FeatureHit {
                feature_idx: 0,
                gene_id: 0,
                overlap_len: 100,
            },
            FeatureHit {
                feature_idx: 1,
                gene_id: 1,
                overlap_len: 50,
            },
        ];

        let args = make_args(true, true);
        match resolve_assignment(&hits, &args) {
            Assignment::Unique(h) => {
                assert_eq!(h.gene_id, 0);
                assert_eq!(h.overlap_len, 100);
            }
            _ => panic!("Expected unique assignment with largest overlap"),
        }
    }
}
