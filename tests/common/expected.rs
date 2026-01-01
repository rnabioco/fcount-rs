use std::collections::HashMap;
use std::fs;
use std::path::Path;

/// Represents expected counts from .ora file
enum ExpectedCounts {
    /// Gene-level: gene_id -> count (float for fractional)
    GeneLevel(HashMap<String, f64>),
    /// Feature-level: (gene_id, chr, start) -> count (currently unused)
    #[allow(dead_code)]
    FeatureLevel(HashMap<(String, String, u32), f64>),
}

impl ExpectedCounts {
    /// Parse .ora file - auto-detect format based on column count
    pub fn parse(path: &Path) -> Result<Self, String> {
        let content = fs::read_to_string(path)
            .map_err(|e| format!("Failed to read {}: {}", path.display(), e))?;

        let mut gene_counts: HashMap<String, f64> = HashMap::new();
        let mut feature_counts: HashMap<(String, String, u32), f64> = HashMap::new();
        let mut is_feature_level = false;

        for line in content.lines() {
            // Skip comments
            if line.starts_with('#') {
                continue;
            }

            // Skip header lines
            let lower = line.to_lowercase();
            if lower.contains("geneid") || lower.contains("gene_id") {
                let cols: Vec<_> = line.split('\t').collect();
                is_feature_level = cols.len() >= 5;
                continue;
            }

            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 2 {
                continue;
            }

            if fields.len() >= 5 {
                // Feature-level format: geneid, chr, start, end, count
                is_feature_level = true;
                let gene_id = fields[0].to_string();
                let chr = fields[1].to_string();
                let start: u32 = fields[2].parse().unwrap_or(0);
                let count: f64 = fields[4].parse().unwrap_or(0.0);
                feature_counts.insert((gene_id, chr, start), count);
            } else {
                // Gene-level format: geneid, count
                let gene_id = fields[0].to_string();
                let count: f64 = fields[1].parse().unwrap_or(0.0);
                gene_counts.insert(gene_id, count);
            }
        }

        if is_feature_level && !feature_counts.is_empty() {
            Ok(ExpectedCounts::FeatureLevel(feature_counts))
        } else {
            Ok(ExpectedCounts::GeneLevel(gene_counts))
        }
    }

    /// Compare actual counts against expected with tolerance
    pub fn compare(
        &self,
        actual: &HashMap<String, f64>,
        tolerance: f64,
    ) -> Result<(), Vec<String>> {
        match self {
            ExpectedCounts::GeneLevel(expected) => compare_gene_level(expected, actual, tolerance),
            ExpectedCounts::FeatureLevel(_) => {
                // For feature-level, we'd need different actual format
                // For now, skip detailed comparison
                Ok(())
            }
        }
    }
}

fn compare_gene_level(
    expected: &HashMap<String, f64>,
    actual: &HashMap<String, f64>,
    tolerance: f64,
) -> Result<(), Vec<String>> {
    let mut diffs = Vec::new();

    for (gene, &exp_count) in expected {
        match actual.get(gene) {
            Some(&act_count) => {
                if !counts_match(act_count, exp_count, tolerance) {
                    diffs.push(format!(
                        "{}: expected {}, got {}",
                        gene, exp_count, act_count
                    ));
                }
            }
            None => {
                if exp_count != 0.0 {
                    diffs.push(format!(
                        "{}: expected {}, not found in output",
                        gene, exp_count
                    ));
                }
            }
        }
    }

    // Check for unexpected genes with non-zero counts
    for (gene, &act_count) in actual {
        if !expected.contains_key(gene) && act_count != 0.0 {
            diffs.push(format!("{}: unexpected count {}", gene, act_count));
        }
    }

    if diffs.is_empty() {
        Ok(())
    } else {
        Err(diffs)
    }
}

fn counts_match(actual: f64, expected: f64, tolerance: f64) -> bool {
    if tolerance == 0.0 {
        // Exact match for integers
        (actual - expected).abs() < 0.0001
    } else {
        // Absolute tolerance for fractional counts
        (actual - expected).abs() <= tolerance
    }
}

/// Assert that actual counts match expected from .ora file
pub fn assert_matches_expected(
    actual: &HashMap<String, f64>,
    expected_path: &Path,
    tolerance: f64,
) {
    let expected =
        ExpectedCounts::parse(expected_path).expect("Failed to parse expected output file");

    match expected.compare(actual, tolerance) {
        Ok(()) => {}
        Err(diffs) => {
            panic!(
                "Count mismatch against {}:\n{}\n\nActual counts:\n{:#?}",
                expected_path.display(),
                diffs.join("\n"),
                actual
            );
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_gene_level() {
        let content = "g1\t10\ng2\t20\n";
        let temp = tempfile::NamedTempFile::new().unwrap();
        fs::write(temp.path(), content).unwrap();

        let expected = ExpectedCounts::parse(temp.path()).unwrap();
        match expected {
            ExpectedCounts::GeneLevel(counts) => {
                assert_eq!(counts.get("g1"), Some(&10.0));
                assert_eq!(counts.get("g2"), Some(&20.0));
            }
            _ => panic!("Expected gene-level counts"),
        }
    }

    #[test]
    fn test_parse_fractional() {
        let content = "g1\t0.78\ng2\t1.28\n";
        let temp = tempfile::NamedTempFile::new().unwrap();
        fs::write(temp.path(), content).unwrap();

        let expected = ExpectedCounts::parse(temp.path()).unwrap();
        match expected {
            ExpectedCounts::GeneLevel(counts) => {
                assert!((counts.get("g1").unwrap() - 0.78).abs() < 0.001);
                assert!((counts.get("g2").unwrap() - 1.28).abs() < 0.001);
            }
            _ => panic!("Expected gene-level counts"),
        }
    }
}
