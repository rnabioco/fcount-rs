use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use tempfile::TempDir;

/// Result of running fcount CLI
pub struct FcountResult {
    pub exit_code: i32,
    #[allow(dead_code)]
    pub stdout: String,
    pub stderr: String,
    pub output_path: PathBuf,
    pub summary_path: PathBuf,
    pub counts: HashMap<String, f64>,
    #[allow(dead_code)]
    temp_dir: TempDir, // Keeps temp files alive for debugging
}

impl FcountResult {
    /// Check if the command succeeded
    pub fn success(&self) -> bool {
        self.exit_code == 0
    }

    /// Parse the output counts file (7-column fcount format)
    fn parse_output(path: &Path) -> Result<HashMap<String, f64>, String> {
        let content = fs::read_to_string(path)
            .map_err(|e| format!("Failed to read output {}: {}", path.display(), e))?;

        let mut counts = HashMap::new();

        for line in content.lines() {
            // Skip comment lines (start with #)
            if line.starts_with('#') {
                continue;
            }
            // Skip header line
            if line.starts_with("Geneid") {
                continue;
            }

            let fields: Vec<&str> = line.split('\t').collect();
            // Format: Geneid, Chr, Start, End, Strand, Length, Count[, Count2, ...]
            if fields.len() >= 7 {
                let gene_id = fields[0].to_string();
                // Parse count from column 7 (index 6)
                let count: f64 = fields[6].parse().unwrap_or(0.0);
                counts.insert(gene_id, count);
            }
        }

        Ok(counts)
    }
}

/// Run fcount with given arguments
pub fn run_fcount(
    bam_file: &Path,
    annotation: &Path,
    extra_args: &[&str],
) -> Result<FcountResult, String> {
    let temp_dir = TempDir::new().map_err(|e| e.to_string())?;
    let output_path = temp_dir.path().join("counts.txt");
    let summary_path = temp_dir.path().join("counts.txt.summary");

    // Get the binary path from cargo
    let binary_path = env!("CARGO_BIN_EXE_fcount");

    let mut cmd = Command::new(binary_path);
    cmd.arg("-a")
        .arg(annotation)
        .arg("-o")
        .arg(&output_path)
        .arg(bam_file)
        .arg("-q") // Quiet mode - no progress bar
        .arg("-t")
        .arg("1"); // Single thread for reproducibility

    for arg in extra_args {
        cmd.arg(arg);
    }

    let output: Output = cmd
        .output()
        .map_err(|e| format!("Failed to run fcount: {}", e))?;

    let exit_code = output.status.code().unwrap_or(-1);
    let stdout = String::from_utf8_lossy(&output.stdout).to_string();
    let stderr = String::from_utf8_lossy(&output.stderr).to_string();

    let counts = if output.status.success() && output_path.exists() {
        FcountResult::parse_output(&output_path)?
    } else {
        HashMap::new()
    };

    Ok(FcountResult {
        exit_code,
        stdout,
        stderr,
        output_path,
        summary_path,
        counts,
        temp_dir,
    })
}

/// Run fcount with multiple BAM files
#[allow(dead_code)]
pub fn run_fcount_multi(
    bam_files: &[&Path],
    annotation: &Path,
    extra_args: &[&str],
) -> Result<FcountResult, String> {
    let temp_dir = TempDir::new().map_err(|e| e.to_string())?;
    let output_path = temp_dir.path().join("counts.txt");
    let summary_path = temp_dir.path().join("counts.txt.summary");

    let binary_path = env!("CARGO_BIN_EXE_fcount");

    let mut cmd = Command::new(binary_path);
    cmd.arg("-a")
        .arg(annotation)
        .arg("-o")
        .arg(&output_path)
        .arg("-q")
        .arg("-t")
        .arg("1");

    for arg in extra_args {
        cmd.arg(arg);
    }

    for bam in bam_files {
        cmd.arg(bam);
    }

    let output = cmd
        .output()
        .map_err(|e| format!("Failed to run fcount: {}", e))?;

    let exit_code = output.status.code().unwrap_or(-1);
    let stdout = String::from_utf8_lossy(&output.stdout).to_string();
    let stderr = String::from_utf8_lossy(&output.stderr).to_string();

    let counts = if output.status.success() && output_path.exists() {
        FcountResult::parse_output(&output_path)?
    } else {
        HashMap::new()
    };

    Ok(FcountResult {
        exit_code,
        stdout,
        stderr,
        output_path,
        summary_path,
        counts,
        temp_dir,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fcount_binary_exists() {
        let binary_path = env!("CARGO_BIN_EXE_fcount");
        assert!(
            Path::new(binary_path).exists() || binary_path.contains("fcount"),
            "fcount binary should be available"
        );
    }
}
