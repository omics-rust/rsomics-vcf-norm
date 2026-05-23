use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

fn ours() -> Command {
    Command::new(env!("CARGO_BIN_EXE_rsomics-vcf-norm"))
}

fn fixture() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/golden/multiallelic.vcf")
}

fn bcftools_available() -> bool {
    Command::new("bcftools")
        .arg("--version")
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .map(|s| s.success())
        .unwrap_or(false)
}

/// Data (non-header) records — headers differ (bcftools stamps its command line).
fn records(vcf: &[u8]) -> Vec<String> {
    String::from_utf8_lossy(vcf)
        .lines()
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .map(str::to_owned)
        .collect()
}

/// Compare data records from ours vs bcftools for a given VCF path.
fn split_cmp(path: &std::path::Path) {
    let ours_out = ours().arg("-m").arg(path).output().unwrap();
    assert!(
        ours_out.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&ours_out.stderr)
    );
    let theirs = Command::new("bcftools")
        .args(["norm", "-m-"])
        .arg(path)
        .output()
        .unwrap();
    assert!(theirs.status.success());
    assert_eq!(records(&ours_out.stdout), records(&theirs.stdout));
}

// Splitting multiallelics (-m) must produce the same biallelic records as
// `bcftools norm -m-` (compare data lines; headers differ between tools).
#[test]
fn split_matches_bcftools() {
    if !bcftools_available() {
        eprintln!("skipping: bcftools not found");
        return;
    }
    split_cmp(&fixture());
}

// Multiallelic split with FORMAT/sample columns: GT allele indices must be
// remapped and Number=A INFO fields split — the original fixture had no sample
// columns so that bug was hidden.
#[test]
fn split_multisample_matches_bcftools() {
    if !bcftools_available() {
        eprintln!("skipping: bcftools not found");
        return;
    }
    let multisample = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/golden/multisample.vcf");
    split_cmp(&multisample);
}
