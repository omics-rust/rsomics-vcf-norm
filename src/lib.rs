use std::io::{self, Read};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

pub struct NormStats {
    pub total: u64,
    pub split: u64,
    pub realigned: u64,
}

/// Byte index of the `n`-th tab (1-based) in `line`, or None if there are fewer.
fn nth_tab(line: &[u8], n: usize) -> Option<usize> {
    line.iter()
        .enumerate()
        .filter(|&(_, &b)| b == b'\t')
        .map(|(i, _)| i)
        .nth(n - 1)
}

/// Split multiallelic sites into biallelic records by a tab byte-scan: biallelic
/// lines pass through verbatim, multiallelic lines (comma in ALT, col 5) emit one
/// record per allele with ALT replaced and every other field kept. No full record
/// parse/re-serialization — matches `bcftools norm -m-` for the split (no -f ref,
/// so no left-align/trim).
pub fn normalize_vcf(
    input: &Path,
    output: &mut dyn io::Write,
    split_multiallelic: bool,
) -> Result<NormStats> {
    let raw = std::fs::read(input)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", input.display())))?;
    let data = if raw.starts_with(&[0x1f, 0x8b]) {
        let mut d = Vec::new();
        flate2::read::MultiGzDecoder::new(&raw[..])
            .read_to_end(&mut d)
            .map_err(RsomicsError::Io)?;
        d
    } else {
        raw
    };

    let mut stats = NormStats {
        total: 0,
        split: 0,
        realigned: 0,
    };
    for raw_line in data.split(|&b| b == b'\n') {
        let line = match raw_line.last() {
            Some(b'\r') => &raw_line[..raw_line.len() - 1],
            _ => raw_line,
        };
        if line.is_empty() {
            continue;
        }
        if line[0] == b'#' {
            output.write_all(line).map_err(RsomicsError::Io)?;
            output.write_all(b"\n").map_err(RsomicsError::Io)?;
            continue;
        }
        stats.total += 1;

        if split_multiallelic && let Some(t4) = nth_tab(line, 4) {
            let alt_start = t4 + 1;
            let alt_end = nth_tab(line, 5).unwrap_or(line.len());
            let alt = &line[alt_start..alt_end];
            if alt.contains(&b',') {
                let prefix = &line[..alt_start];
                let suffix = &line[alt_end..];
                for allele in alt.split(|&b| b == b',') {
                    output.write_all(prefix).map_err(RsomicsError::Io)?;
                    output.write_all(allele).map_err(RsomicsError::Io)?;
                    output.write_all(suffix).map_err(RsomicsError::Io)?;
                    output.write_all(b"\n").map_err(RsomicsError::Io)?;
                    stats.split += 1;
                }
                continue;
            }
        }

        output.write_all(line).map_err(RsomicsError::Io)?;
        output.write_all(b"\n").map_err(RsomicsError::Io)?;
    }

    Ok(stats)
}
