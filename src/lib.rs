use std::io::{self, Read};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

pub struct NormStats {
    pub total: u64,
    pub split: u64,
    pub realigned: u64,
}

/// Split a VCF data line for one target alt allele (1-based index into the original ALT list).
/// Rewrites the ALT field to just the target allele, remaps INFO Number=A fields by index,
/// and remaps GT allele indices so that the target allele becomes 1 and all other non-ref
/// alleles become 0 — matching `bcftools norm -m-` semantics for per-allele split.
fn split_record_for_allele(
    fields: &[&str],
    alts: &[&str],
    allele_idx: usize,
    per_allele_info_keys: &[String],
    out: &mut dyn io::Write,
) -> io::Result<()> {
    // col indices: 0=CHROM 1=POS 2=ID 3=REF 4=ALT 5=QUAL 6=FILTER 7=INFO [8=FORMAT 9..]
    let orig_alt_num = allele_idx + 1; // 1-based allele index in the original record

    for f in &fields[..4] {
        out.write_all(f.as_bytes())?;
        out.write_all(b"\t")?;
    }
    out.write_all(alts[allele_idx].as_bytes())?;
    out.write_all(b"\t")?;
    out.write_all(fields[5].as_bytes())?;
    out.write_all(b"\t")?;
    out.write_all(fields[6].as_bytes())?;

    out.write_all(b"\t")?;
    if fields[7] == "." || per_allele_info_keys.is_empty() {
        out.write_all(fields[7].as_bytes())?;
    } else {
        let info_str = fields[7];
        let mut parts: Vec<String> = Vec::new();
        for kv in info_str.split(';') {
            if let Some((key, val)) = kv.split_once('=') {
                if per_allele_info_keys.iter().any(|k| k == key) {
                    let picked = val.split(',').nth(allele_idx).unwrap_or(val);
                    parts.push(format!("{key}={picked}"));
                } else {
                    parts.push(kv.to_owned());
                }
            } else {
                parts.push(kv.to_owned());
            }
        }
        out.write_all(parts.join(";").as_bytes())?;
    }

    if fields.len() > 8 {
        let fmt = fields[8];
        out.write_all(b"\t")?;
        out.write_all(fmt.as_bytes())?;

        let gt_pos: Option<usize> = fmt
            .split(':')
            .enumerate()
            .find_map(|(i, f)| (f == "GT").then_some(i));

        for sample in &fields[9..] {
            out.write_all(b"\t")?;
            if let Some(gp) = gt_pos {
                let sub: Vec<&str> = sample.split(':').collect();
                if let Some(gt_raw) = sub.get(gp) {
                    let sep = if gt_raw.contains('|') { '|' } else { '/' };
                    let remapped: Vec<String> = gt_raw
                        .split(sep)
                        .map(|a| match a.parse::<u32>() {
                            Err(_) => a.to_owned(),
                            Ok(0) => "0".to_owned(),
                            Ok(n) if n as usize == orig_alt_num => "1".to_owned(),
                            // Other alt alleles remap to ref — bcftools norm -m- semantics.
                            Ok(_) => "0".to_owned(),
                        })
                        .collect();
                    let sep_str = if sep == '|' { "|" } else { "/" };
                    let mut full_sample = Vec::with_capacity(sub.len());
                    for pre in &sub[..gp] {
                        full_sample.push(pre.to_string());
                    }
                    full_sample.push(remapped.join(sep_str));
                    for other in &sub[gp + 1..] {
                        full_sample.push(other.to_string());
                    }
                    out.write_all(full_sample.join(":").as_bytes())?;
                } else {
                    out.write_all(sample.as_bytes())?;
                }
            } else {
                out.write_all(sample.as_bytes())?;
            }
        }
    }

    out.write_all(b"\n")
}

/// Parse `##INFO=<ID=...,Number=A,...>` lines to collect per-allele INFO keys.
fn collect_per_allele_info_keys(header: &str) -> Vec<String> {
    let mut keys = Vec::new();
    for line in header.lines() {
        if !line.starts_with("##INFO=") {
            continue;
        }
        // Quick parse: look for Number=A/R and ID=<key>
        if (line.contains("Number=A") || line.contains("Number=R"))
            && let Some(id_start) = line.find("ID=")
        {
            let rest = &line[id_start + 3..];
            let id_end = rest.find([',', '>']).unwrap_or(rest.len());
            keys.push(rest[..id_end].to_owned());
        }
    }
    keys
}

/// Split multiallelic sites into biallelic records.
/// Matches `bcftools norm -m-`: for each alt allele emit one record, reindex GT
/// allele numbers, and split Number=A INFO fields. No left-align/trim (no -f ref).
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

    let text = std::str::from_utf8(&data).map_err(|e| RsomicsError::InvalidInput(e.to_string()))?;

    // Collect header to find Number=A INFO keys before writing records
    let mut header_buf = String::new();
    let mut per_allele_keys: Vec<String> = Vec::new();
    let mut header_done = false;

    let mut stats = NormStats {
        total: 0,
        split: 0,
        realigned: 0,
    };

    for raw_line in text.split('\n') {
        let line = raw_line.trim_end_matches('\r');
        if line.is_empty() {
            continue;
        }

        if line.starts_with('#') {
            if !header_done {
                header_buf.push_str(line);
                header_buf.push('\n');
            }
            if line.starts_with("#CHROM") {
                // Finished meta-header: parse per-allele keys
                per_allele_keys = collect_per_allele_info_keys(&header_buf);
                header_done = true;
            }
            output
                .write_all(line.as_bytes())
                .map_err(RsomicsError::Io)?;
            output.write_all(b"\n").map_err(RsomicsError::Io)?;
            continue;
        }

        stats.total += 1;

        if split_multiallelic {
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() >= 5 {
                let alt_field = fields[4];
                if alt_field.contains(',') {
                    let alts: Vec<&str> = alt_field.split(',').collect();
                    let n = alts.len();
                    for (i, _) in alts.iter().enumerate() {
                        split_record_for_allele(&fields, &alts, i, &per_allele_keys, output)
                            .map_err(RsomicsError::Io)?;
                        stats.split += 1;
                    }
                    // One split call per allele counts as (n-1) extra records emitted;
                    // undo the one stat.total++ above so total stays per input record.
                    let _ = n; // n used only for clarity; stats.split += n above
                    continue;
                }
            }
        }

        output
            .write_all(line.as_bytes())
            .map_err(RsomicsError::Io)?;
        output.write_all(b"\n").map_err(RsomicsError::Io)?;
    }

    Ok(stats)
}
