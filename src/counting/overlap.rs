use rustc_hash::FxHashSet;
use smallvec::SmallVec;

use crate::alignment::{AlignmentRecord, Interval, PendingMate};
use crate::annotation::{Feature, Strand};
use crate::cli::{Args, StrandMode};

/// A hit to a feature
#[derive(Debug, Clone)]
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

/// Check if read strand matches feature strand for single-end reads
pub fn check_strand(record: &AlignmentRecord, feature: &Feature, args: &Args) -> bool {
    match args.strand_mode() {
        StrandMode::Unstranded => true,
        StrandMode::Stranded => {
            let read_strand = if record.is_reverse_strand() {
                Strand::Reverse
            } else {
                Strand::Forward
            };
            feature.strand == read_strand || feature.strand == Strand::Unknown
        }
        StrandMode::ReverselyStranded => {
            let read_strand = if record.is_reverse_strand() {
                Strand::Forward // Reverse of the read strand
            } else {
                Strand::Reverse
            };
            feature.strand == read_strand || feature.strand == Strand::Unknown
        }
    }
}

/// Check if fragment strand matches feature strand for paired-end reads
pub fn check_strand_paired(
    record: &AlignmentRecord,
    mate: &PendingMate,
    feature: &Feature,
    args: &Args,
) -> bool {
    match args.strand_mode() {
        StrandMode::Unstranded => true,
        StrandMode::Stranded => {
            // For stranded paired-end: first read determines strand
            let first_read_reverse = if record.is_first_in_pair() {
                record.is_reverse_strand()
            } else {
                mate.is_reverse_strand()
            };

            let fragment_strand = if first_read_reverse {
                Strand::Reverse
            } else {
                Strand::Forward
            };

            feature.strand == fragment_strand || feature.strand == Strand::Unknown
        }
        StrandMode::ReverselyStranded => {
            // For reversely stranded: second read determines strand
            let second_read_reverse = if record.is_second_in_pair() {
                record.is_reverse_strand()
            } else {
                mate.is_reverse_strand()
            };

            let fragment_strand = if second_read_reverse {
                Strand::Forward
            } else {
                Strand::Reverse
            };

            feature.strand == fragment_strand || feature.strand == Strand::Unknown
        }
    }
}

/// Check strand compatibility using Strand enum directly (for parallel processing)
pub fn check_strand_with_strand(read_strand: Strand, feature: &Feature, args: &Args) -> bool {
    match args.strand_mode() {
        StrandMode::Unstranded => true,
        StrandMode::Stranded => feature.strand == read_strand || feature.strand == Strand::Unknown,
        StrandMode::ReverselyStranded => {
            let expected_strand = match read_strand {
                Strand::Forward => Strand::Reverse,
                Strand::Reverse => Strand::Forward,
                Strand::Unknown => Strand::Unknown,
            };
            feature.strand == expected_strand || feature.strand == Strand::Unknown
        }
    }
}

/// Check strand compatibility for paired-end using Strand enums directly
pub fn check_strand_paired_with_strands(
    _record_strand: Strand,
    _mate_strand: Strand,
    feature: &Feature,
    args: &Args,
) -> bool {
    // For paired-end in parallel mode, we don't have first/second in pair info
    // in a simple way, so we just check unstranded for now
    // TODO: Pass first_in_pair flag to properly handle stranded paired-end
    match args.strand_mode() {
        StrandMode::Unstranded => true,
        StrandMode::Stranded | StrandMode::ReverselyStranded => {
            // For now, accept if feature strand is unknown or if either read matches
            // This is a simplification - proper implementation needs first_in_pair info
            feature.strand == Strand::Unknown
                || feature.strand == _record_strand
                || feature.strand == _mate_strand
        }
    }
}

/// Check if overlap meets threshold requirements
pub fn check_overlap_thresholds(
    overlap_len: u32,
    intervals: &[Interval],
    feature: &Feature,
    args: &Args,
) -> bool {
    // Check minimum absolute overlap
    if overlap_len < args.min_overlap_bases {
        return false;
    }

    // Check minimum fraction of read
    if args.min_overlap_fraction > 0.0 {
        let read_len: u32 = intervals.iter().map(|i| i.len()).sum();
        if read_len > 0 {
            let frac = overlap_len as f32 / read_len as f32;
            if frac < args.min_overlap_fraction {
                return false;
            }
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
        return Assignment::Unique(hits[0].clone());
    }

    // Quick check: are all hits from the same gene? (common case)
    let first_gene = hits[0].gene_id;
    let all_same_gene = hits.iter().skip(1).all(|h| h.gene_id == first_gene);

    if all_same_gene {
        // All hits are from the same gene - pick the one with best overlap
        if args.largest_overlap_only {
            let best = hits.iter().max_by_key(|h| h.overlap_len).unwrap();
            return Assignment::Unique(best.clone());
        }
        // Multiple exons of same gene - count as unique to that gene
        return Assignment::Unique(hits[0].clone());
    }

    // Multiple genes
    if args.allow_multi_overlap {
        if args.largest_overlap_only {
            // Find the gene(s) with largest overlap - single pass
            let max_overlap = hits.iter().map(|h| h.overlap_len).max().unwrap();

            // Use SmallVec to avoid heap allocation for typical cases
            let mut best_hits: SmallVec<[FeatureHit; 8]> = SmallVec::new();
            let mut seen_genes: SmallVec<[u32; 8]> = SmallVec::new();

            for h in hits {
                if h.overlap_len == max_overlap && !seen_genes.contains(&h.gene_id) {
                    seen_genes.push(h.gene_id);
                    best_hits.push(h.clone());
                }
            }

            if best_hits.len() == 1 {
                return Assignment::Unique(best_hits.into_iter().next().unwrap());
            }

            return Assignment::MultiOverlap(best_hits.into_vec());
        }

        // Deduplicate by gene using SmallVec
        let mut seen_genes: SmallVec<[u32; 8]> = SmallVec::new();
        let mut deduped: SmallVec<[FeatureHit; 8]> = SmallVec::new();

        for h in hits {
            if !seen_genes.contains(&h.gene_id) {
                seen_genes.push(h.gene_id);
                deduped.push(h.clone());
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
            count_chimeric: false,
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
