use anyhow::Result;
use rustc_hash::FxHashMap;
use std::sync::Arc;

use super::feature::Feature;

/// Bucket size for spatial index (128KB, matching C implementation)
const BUCKET_SIZE: u32 = 131072;

/// Per-chromosome spatial index using bucket-based lookup
#[derive(Debug)]
pub struct ChromIndex {
    /// Chromosome ID
    pub chrom_id: u16,
    /// Each bucket points to a range of feature indices: (start_idx, end_idx)
    /// Bucket i covers positions [i * BUCKET_SIZE, (i+1) * BUCKET_SIZE)
    pub buckets: Vec<(u32, u32)>,
    /// Maximum feature end position in each bucket (for early termination)
    pub bucket_max_end: Vec<u32>,
}

impl ChromIndex {
    /// Find features that potentially overlap the query interval
    /// Returns an iterator over feature indices
    pub fn find_overlapping(&self, start: u32, end: u32) -> impl Iterator<Item = u32> + '_ {
        let start_bucket = (start / BUCKET_SIZE) as usize;
        let end_bucket = (end / BUCKET_SIZE) as usize;

        let max_bucket = self.buckets.len().saturating_sub(1);

        (start_bucket..=end_bucket.min(max_bucket)).flat_map(move |bucket_idx| {
            // Early termination: if bucket's max end < query start, skip
            if self.bucket_max_end.get(bucket_idx).copied().unwrap_or(0) < start {
                return 0..0;
            }

            let (feat_start, feat_end) = self.buckets.get(bucket_idx).copied().unwrap_or((0, 0));
            feat_start..feat_end
        })
    }
}

/// Main annotation index containing all features and spatial indices
#[derive(Debug)]
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

        // Build per-chromosome indices
        let mut chrom_indices: Vec<ChromIndex> = Vec::with_capacity(num_chroms);

        for chrom_id in 0..num_chroms as u16 {
            // Find feature range for this chromosome
            let start_idx = features
                .iter()
                .position(|f| f.chrom_id == chrom_id)
                .unwrap_or(features.len());

            let end_idx = features
                .iter()
                .rposition(|f| f.chrom_id == chrom_id)
                .map(|i| i + 1)
                .unwrap_or(start_idx);

            if start_idx >= end_idx {
                // No features on this chromosome
                chrom_indices.push(ChromIndex {
                    chrom_id,
                    buckets: vec![],
                    bucket_max_end: vec![],
                });
                continue;
            }

            // Find max position to determine number of buckets
            let max_pos = features[start_idx..end_idx]
                .iter()
                .map(|f| f.end)
                .max()
                .unwrap_or(0);

            let num_buckets = (max_pos / BUCKET_SIZE + 1) as usize;

            // Build buckets
            let mut buckets: Vec<(u32, u32)> = vec![(0, 0); num_buckets];
            let mut bucket_max_end: Vec<u32> = vec![0; num_buckets];

            // Assign features to buckets based on their start position
            // Features are already sorted by start, so we can do a single pass
            let mut current_bucket = 0;
            let mut bucket_start = start_idx as u32;

            for (i, feature) in features[start_idx..end_idx].iter().enumerate() {
                let feat_idx = start_idx + i;
                let feature_bucket = (feature.start / BUCKET_SIZE) as usize;

                // Close out previous buckets
                while current_bucket < feature_bucket && current_bucket < num_buckets {
                    buckets[current_bucket] = (bucket_start, feat_idx as u32);
                    current_bucket += 1;
                    bucket_start = feat_idx as u32;
                }

                // Update max end for all buckets this feature could overlap
                // A feature starting in bucket B could overlap buckets B through B + (len/BUCKET_SIZE)
                let end_bucket = (feature.end / BUCKET_SIZE) as usize;
                for b in feature_bucket..=end_bucket.min(num_buckets - 1) {
                    bucket_max_end[b] = bucket_max_end[b].max(feature.end);
                }
            }

            // Close out remaining buckets
            while current_bucket < num_buckets {
                buckets[current_bucket] = (bucket_start, end_idx as u32);
                current_bucket += 1;
            }

            chrom_indices.push(ChromIndex {
                chrom_id,
                buckets,
                bucket_max_end,
            });
        }

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

    /// Get gene name for a gene ID
    #[inline]
    pub fn get_gene_name(&self, gene_id: u32) -> Option<&str> {
        self.gene_names.get(gene_id as usize).map(|s| s.as_ref())
    }

    /// Get feature by index
    #[inline]
    pub fn get_feature(&self, idx: u32) -> Option<&Feature> {
        self.features.get(idx as usize)
    }

    /// Find all features overlapping the given interval on a chromosome
    pub fn find_overlapping(
        &self,
        chrom_id: u16,
        start: u32,
        end: u32,
    ) -> impl Iterator<Item = (u32, &Feature)> + '_ {
        let chrom_index = self.chrom_indices.get(chrom_id as usize);

        chrom_index
            .into_iter()
            .flat_map(move |idx| idx.find_overlapping(start, end))
            .filter_map(move |feat_idx| {
                let feature = self.features.get(feat_idx as usize)?;
                if feature.overlaps(start, end) {
                    Some((feat_idx, feature))
                } else {
                    None
                }
            })
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
    fn test_find_overlapping() {
        let index = make_test_index();

        // Query that overlaps first feature
        let hits: Vec<_> = index.find_overlapping(0, 150, 180).collect();
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].1.start, 100);
        assert_eq!(hits[0].1.end, 200);

        // Query that overlaps no features
        let hits: Vec<_> = index.find_overlapping(0, 210, 290).collect();
        assert_eq!(hits.len(), 0);

        // Query that overlaps two features
        let hits: Vec<_> = index.find_overlapping(0, 180, 320).collect();
        assert_eq!(hits.len(), 2);
    }

    #[test]
    fn test_get_gene_name() {
        let index = make_test_index();
        assert_eq!(index.get_gene_name(0), Some("GENE1"));
        assert_eq!(index.get_gene_name(1), Some("GENE2"));
        assert_eq!(index.get_gene_name(99), None);
    }
}
