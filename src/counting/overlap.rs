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

/// The identifier a hit is counted against: the feature itself under
/// `--feature-level`, otherwise its gene.
///
/// Every grouping decision below keys on this rather than on `gene_id`. Keying
/// on the gene unconditionally — as this used to — meant that at feature level
/// two exons of one gene collapsed to a single hit, so a read spanning both
/// credited only the first even when `--multi-overlap` asked for both. That
/// lands squarely on the DEXSeq output, which forces `--feature-level` on.
#[inline(always)]
fn target_id(hit: &FeatureHit, feature_level: bool) -> u32 {
    if feature_level {
        hit.feature_idx
    } else {
        hit.gene_id
    }
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

    let feature_level = args.feature_level;

    // Quick check: do all hits resolve to the same target? (common case)
    let first_target = target_id(&hits[0], feature_level);
    let all_same_target = hits
        .iter()
        .skip(1)
        .all(|h| target_id(h, feature_level) == first_target);

    if all_same_target {
        // All hits are the same target - pick the one with best overlap
        if args.largest_overlap_only {
            let best = hits.iter().max_by_key(|h| h.overlap_len).unwrap();
            return Assignment::Unique(*best); // Copy, not clone
        }
        // e.g. multiple exons of the same gene at gene level - count once
        return Assignment::Unique(hits[0]); // Copy, not clone
    }

    // Multiple distinct targets
    if args.allow_multi_overlap {
        if args.largest_overlap_only {
            // Find the target(s) with largest overlap - single pass
            let max_overlap = hits.iter().map(|h| h.overlap_len).max().unwrap();

            // Use FxHashSet for O(1) deduplication instead of O(n) contains
            let mut best_hits: SmallVec<[FeatureHit; 8]> = SmallVec::new();
            let mut seen: FxHashSet<u32> = FxHashSet::default();

            for h in hits {
                if h.overlap_len == max_overlap && seen.insert(target_id(h, feature_level)) {
                    best_hits.push(*h); // Copy
                }
            }

            if best_hits.len() == 1 {
                return Assignment::Unique(best_hits[0]);
            }

            return Assignment::MultiOverlap(best_hits.into_vec());
        }

        // Deduplicate by target using FxHashSet for O(1) lookup
        let mut seen: FxHashSet<u32> = FxHashSet::default();
        let mut deduped: SmallVec<[FeatureHit; 8]> = SmallVec::new();

        for h in hits {
            if seen.insert(target_id(h, feature_level)) {
                deduped.push(*h); // Copy
            }
        }

        return Assignment::MultiOverlap(deduped.into_vec());
    }

    // Multiple targets without multi-overlap - ambiguous
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

    /// Two exons of the *same* gene, which is what a spliced read spanning an
    /// exon-exon junction produces.
    fn two_exons_one_gene() -> Vec<FeatureHit> {
        vec![
            FeatureHit {
                feature_idx: 0,
                gene_id: 7,
                overlap_len: 40,
            },
            FeatureHit {
                feature_idx: 1,
                gene_id: 7,
                overlap_len: 60,
            },
        ]
    }

    /// At gene level these collapse to one hit on the gene — unchanged behavior.
    #[test]
    fn test_gene_level_collapses_exons_of_same_gene() {
        let args = make_args(true, false); // -O, gene level
        match resolve_assignment(&two_exons_one_gene(), &args) {
            Assignment::Unique(h) => assert_eq!(h.gene_id, 7),
            other => panic!("expected Unique at gene level, got {:?}", other),
        }
    }

    /// Regression: at feature level with `-O`, both exons must be counted. The
    /// dedup used to key on `gene_id` regardless of level, so the second exon
    /// was dropped — silently halving per-exon counts on the DEXSeq path, which
    /// forces `--feature-level` on.
    #[test]
    fn test_feature_level_multi_overlap_keeps_both_exons() {
        let mut args = make_args(true, false); // -O
        args.feature_level = true;
        match resolve_assignment(&two_exons_one_gene(), &args) {
            Assignment::MultiOverlap(hits) => {
                assert_eq!(hits.len(), 2, "both exons should be counted");
                let mut idx: Vec<u32> = hits.iter().map(|h| h.feature_idx).collect();
                idx.sort_unstable();
                assert_eq!(idx, vec![0, 1]);
            }
            other => panic!("expected MultiOverlap at feature level, got {:?}", other),
        }
    }

    /// At feature level *without* `-O`, a read touching two features is
    /// ambiguous — matching featureCounts, which only assigns to every
    /// overlapping feature when `-O` is given.
    #[test]
    fn test_feature_level_without_multi_overlap_is_ambiguous() {
        let mut args = make_args(false, false);
        args.feature_level = true;
        match resolve_assignment(&two_exons_one_gene(), &args) {
            Assignment::Ambiguous => {}
            other => panic!("expected Ambiguous, got {:?}", other),
        }
    }

    /// `--largest-overlap` at feature level picks the single best exon rather
    /// than collapsing on gene.
    #[test]
    fn test_feature_level_largest_overlap_picks_best_exon() {
        let mut args = make_args(true, true);
        args.feature_level = true;
        match resolve_assignment(&two_exons_one_gene(), &args) {
            Assignment::Unique(h) => {
                assert_eq!(h.feature_idx, 1);
                assert_eq!(h.overlap_len, 60);
            }
            other => panic!("expected Unique, got {:?}", other),
        }
    }

    /// A single hit is still Unique at either level.
    #[test]
    fn test_feature_level_single_hit_unchanged() {
        let hits = vec![FeatureHit {
            feature_idx: 3,
            gene_id: 7,
            overlap_len: 10,
        }];
        let mut args = make_args(true, false);
        args.feature_level = true;
        match resolve_assignment(&hits, &args) {
            Assignment::Unique(h) => assert_eq!(h.feature_idx, 3),
            other => panic!("expected Unique, got {:?}", other),
        }
    }
}
