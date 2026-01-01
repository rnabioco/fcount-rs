use std::path::PathBuf;

/// Local fixtures directory with corrected test files
fn local_fixtures_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures")
}

/// Original featureCounts test fixtures (symlinked from ext/)
fn ext_fixtures_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("ext/subread/test/featureCounts/data")
}

/// Get path to a specific fixture file
/// Checks local fixtures first, then falls back to ext/ directory
pub fn fixture_path(name: &str) -> PathBuf {
    let local = local_fixtures_dir().join(name);
    if local.exists() {
        local
    } else {
        ext_fixtures_dir().join(name)
    }
}
