use anyhow::Result;
use coitrees::{COITree, GenericInterval, Interval, IntervalTree};
use rustc_hash::FxHashMap;
use std::sync::Arc;

use super::feature::Feature;

/// Per-chromosome spatial index using cache-oblivious interval tree
#[derive(Clone)]
pub struct ChromIndex {
    /// COITree for fast interval queries
    /// Stores feature index as metadata
    tree: COITree<u32, u32>,
}

impl std::fmt::Debug for ChromIndex {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("ChromIndex")
            .field("num_intervals", &self.tree.len())
            .finish()
    }
}

impl ChromIndex {
    /// Create a new ChromIndex from features
    fn new(features: &[Feature], start_idx: usize) -> Self {
        let intervals: Vec<Interval<u32>> = features
            .iter()
            .enumerate()
            .map(|(i, f)| {
                // `Feature.start`/`end` come straight from GTF columns 4 and 5:
                // 1-based and end-*inclusive* (which is why `Feature::len()` is
                // `end - start + 1`). coitrees is also end-inclusive, so the
                // bounds carry over as-is. Subtracting 1 here — as this used to
                // do, under a comment claiming end-exclusive — shrank every
                // feature by a base and silently dropped reads whose only
                // overlap was the feature's last base.
                let first = f.start as i32;
                let last = f.end as i32;
                let feat_idx = (start_idx + i) as u32;
                Interval::new(first, last, feat_idx)
            })
            .collect();

        let tree: COITree<u32, u32> = COITree::new(&intervals);
        ChromIndex { tree }
    }

    /// Find features that overlap the query interval
    /// Calls the provided closure for each overlapping feature index
    #[inline]
    pub fn query<F>(&self, start: u32, end: u32, mut callback: F)
    where
        F: FnMut(u32),
    {
        // `start`/`end` are 1-based end-inclusive (CIGAR-derived read intervals),
        // matching how features were stored above and coitrees' own convention.
        let first = start as i32;
        let last = end as i32;
        self.tree
            .query(first, last, |node| callback(*node.metadata()));
    }
}

/// Main annotation index containing all features and spatial indices
#[derive(Debug, Clone)]
pub struct AnnotationIndex {
    /// Chromosome name to ID mapping
    pub chrom_to_id: FxHashMap<Arc<str>, u16>,
    /// Chromosome ID to name mapping
    pub id_to_chrom: Vec<Arc<str>>,
    /// Gene names (indexed by gene_id)
    pub gene_names: Vec<Arc<str>>,
    /// All features, sorted by (chrom_id, start)
    pub features: Vec<Feature>,
    /// Per-chromosome spatial index
    pub chrom_indices: Vec<ChromIndex>,
}

impl AnnotationIndex {
    /// Build annotation index from parsed features
    pub fn new(
        chrom_to_id: FxHashMap<Arc<str>, u16>,
        id_to_chrom: Vec<Arc<str>>,
        gene_names: Vec<Arc<str>>,
        features: Vec<Feature>,
    ) -> Result<Self> {
        let num_chroms = id_to_chrom.len();

        // Build chromosome boundaries in a single pass O(n)
        let mut chrom_boundaries: Vec<(usize, usize)> = vec![(0, 0); num_chroms];

        if !features.is_empty() {
            let mut current_chrom = features[0].chrom_id as usize;
            let mut start_idx = 0;

            for (i, f) in features.iter().enumerate() {
                let fchrom = f.chrom_id as usize;
                if fchrom != current_chrom {
                    chrom_boundaries[current_chrom] = (start_idx, i);
                    // Fill gaps for chromosomes with no features
                    chrom_boundaries[(current_chrom + 1)..fchrom].fill((i, i));
                    current_chrom = fchrom;
                    start_idx = i;
                }
            }
            chrom_boundaries[current_chrom] = (start_idx, features.len());
        }

        // Build per-chromosome COITree indices
        let chrom_indices: Vec<ChromIndex> = (0..num_chroms)
            .map(|chrom_id| {
                let (start_idx, end_idx) = chrom_boundaries[chrom_id];
                if start_idx >= end_idx {
                    ChromIndex::new(&[], 0)
                } else {
                    ChromIndex::new(&features[start_idx..end_idx], start_idx)
                }
            })
            .collect();

        Ok(AnnotationIndex {
            chrom_to_id,
            id_to_chrom,
            gene_names,
            features,
            chrom_indices,
        })
    }

    /// Get chromosome ID for a name, if it exists
    #[inline]
    pub fn get_chrom_id(&self, name: &str) -> Option<u16> {
        self.chrom_to_id.get(name).copied()
    }

    /// Find all features overlapping the given interval on a chromosome
    /// Uses callback-based API for efficiency (avoids allocations)
    #[inline]
    pub fn query_overlapping<F>(&self, chrom_id: u16, start: u32, end: u32, mut callback: F)
    where
        F: FnMut(u32, &Feature),
    {
        if let Some(chrom_index) = self.chrom_indices.get(chrom_id as usize) {
            chrom_index.query(start, end, |feat_idx| {
                if let Some(feature) = self.features.get(feat_idx as usize) {
                    callback(feat_idx, feature);
                }
            });
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::annotation::feature::Strand;

    fn make_test_index() -> AnnotationIndex {
        let mut chrom_to_id = FxHashMap::default();
        chrom_to_id.insert(Arc::from("chr1"), 0);
        chrom_to_id.insert(Arc::from("chr2"), 1);

        let id_to_chrom = vec![Arc::from("chr1"), Arc::from("chr2")];
        let gene_names = vec![Arc::from("GENE1"), Arc::from("GENE2")];

        let features = vec![
            Feature {
                gene_id: 0,
                chrom_id: 0,
                start: 100,
                end: 200,
                strand: Strand::Forward,
            },
            Feature {
                gene_id: 0,
                chrom_id: 0,
                start: 300,
                end: 400,
                strand: Strand::Forward,
            },
            Feature {
                gene_id: 1,
                chrom_id: 0,
                start: 500,
                end: 600,
                strand: Strand::Reverse,
            },
            Feature {
                gene_id: 1,
                chrom_id: 1,
                start: 1000,
                end: 2000,
                strand: Strand::Forward,
            },
        ];

        AnnotationIndex::new(chrom_to_id, id_to_chrom, gene_names, features).unwrap()
    }

    #[test]
    fn test_query_overlapping() {
        let index = make_test_index();

        let mut count = 0;
        index.query_overlapping(0, 150, 180, |_idx, _feat| {
            count += 1;
        });
        assert_eq!(count, 1);
    }

    /// Regression: features and queries are both 1-based end-inclusive, but the
    /// index used to subtract 1 from each end before handing them to coitrees.
    /// That shrank both sides by a base, so a read whose only overlap was the
    /// feature's first or last base was dropped — while featureCounts counts it
    /// (its default minimum overlap is one base).
    #[test]
    fn test_query_boundary_overlaps() {
        // A lone feature, so every assertion below is unambiguous about *which*
        // boundary it is probing.
        let mut chrom_to_id = FxHashMap::default();
        chrom_to_id.insert(Arc::from("chr1"), 0);
        let index = AnnotationIndex::new(
            chrom_to_id,
            vec![Arc::from("chr1")],
            vec![Arc::from("G1")],
            vec![Feature {
                gene_id: 0,
                chrom_id: 0,
                start: 100,
                end: 200,
                strand: Strand::Forward,
            }],
        )
        .unwrap();

        let hits = |start: u32, end: u32| {
            let mut n = 0;
            index.query_overlapping(0, start, end, |_, _| n += 1);
            n
        };

        // Single-base overlaps at each boundary must be found.
        assert_eq!(hits(50, 100), 1, "1 bp overlap at feature start");
        assert_eq!(hits(200, 300), 1, "1 bp overlap at feature end");

        // Two-base overlaps were always found; keep them pinned.
        assert_eq!(hits(50, 101), 1);
        assert_eq!(hits(199, 300), 1);

        // Fully inside, and fully spanning.
        assert_eq!(hits(150, 160), 1);
        assert_eq!(hits(1, 1000), 1);

        // Abutting reads genuinely do not overlap and must stay excluded —
        // this is the half the fix could plausibly have broken.
        assert_eq!(hits(50, 99), 0, "abuts start, no overlap");
        assert_eq!(hits(201, 300), 0, "abuts end, no overlap");
    }
}
