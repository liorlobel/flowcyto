use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, Cursor, Read, Seek, SeekFrom};
use std::path::Path;

use anyhow::{bail, Context, Result};
use byteorder::{BigEndian, LittleEndian, ReadBytesExt};

#[derive(Debug, Clone)]
pub struct Parameter {
    pub index: usize,      // 1-based, matches $PnN
    pub name: String,      // $PnN — detector/channel name
    pub label: Option<String>, // $PnS — staining label
    pub range: f64,        // $PnR — max value (for masking integer data)
    pub bits: u32,         // $PnB — bits per value
}

#[derive(Debug)]
pub struct FcsFile {
    pub version: String,
    pub keywords: HashMap<String, String>,
    pub parameters: Vec<Parameter>,
    pub n_events: usize,
    /// Row-major flat array: events[event * n_params + param]
    pub events: Vec<f64>,
}

impl FcsFile {
    pub fn n_params(&self) -> usize {
        self.parameters.len()
    }

    #[allow(dead_code)]
    pub fn event_slice(&self, i: usize) -> &[f64] {
        let n = self.n_params();
        &self.events[i * n..(i + 1) * n]
    }

    /// Collect all values for one parameter across all events.
    #[allow(dead_code)]
    pub fn channel_values(&self, param_idx: usize) -> Vec<f64> {
        let n = self.n_params();
        self.events.iter().skip(param_idx).step_by(n).copied().collect()
    }

    /// Case-insensitive lookup of a parameter by name; returns its column index.
    pub fn param_index(&self, name: &str) -> Option<usize> {
        self.parameters
            .iter()
            .position(|p| p.name.eq_ignore_ascii_case(name))
    }

    /// Best-effort spillover keyword: tries $SPILLOVER, SPILLOVER, $SPILL, SPILL.
    pub fn spillover_keyword(&self) -> Option<&str> {
        for key in &["$SPILLOVER", "SPILLOVER", "$SPILL", "SPILL"] {
            if let Some(v) = self.keywords.get(*key) {
                return Some(v.as_str());
            }
        }
        None
    }

    /// Lightweight: read just the event count ($TOT) from header + TEXT, without
    /// parsing the DATA segment. Used for per-tube QC in the Samples list.
    pub fn peek_events(path: &Path) -> Result<usize> {
        let file = File::open(path).with_context(|| format!("cannot open {:?}", path))?;
        let mut rdr = BufReader::new(file);
        let mut header = [0u8; 58];
        rdr.read_exact(&mut header).context("header")?;
        let text_start = parse_hdr_offset(&header[10..18])?;
        let text_end = parse_hdr_offset(&header[18..26])?;
        if text_end < text_start {
            bail!("corrupt header (TEXT end < start)");
        }
        rdr.seek(SeekFrom::Start(text_start as u64))?;
        let mut buf = vec![0u8; text_end - text_start + 1];
        rdr.read_exact(&mut buf).context("TEXT")?;
        let kw = parse_text(&buf)?;
        kw.get("$TOT").and_then(|v| v.trim().parse::<usize>().ok())
            .context("$TOT missing/invalid")
    }

    pub fn open(path: &Path) -> Result<Self> {
        let file = File::open(path)
            .with_context(|| format!("cannot open {:?}", path))?;
        let file_len = file.metadata().map(|m| m.len() as usize).unwrap_or(usize::MAX);
        let mut rdr = BufReader::new(file);

        // ── Header (58 bytes) ──────────────────────────────────────────
        let mut header = [0u8; 58];
        rdr.read_exact(&mut header).context("failed to read 58-byte FCS header")?;

        let version = std::str::from_utf8(&header[0..6])
            .context("FCS version not UTF-8")?
            .trim_end()
            .to_string();
        if !version.starts_with("FCS") {
            bail!("not an FCS file (magic = {:?})", version);
        }
        match version.as_str() {
            "FCS2.0" | "FCS3.0" | "FCS3.1" => {}
            v => bail!("unsupported FCS version: {}", v),
        }

        let text_start = parse_hdr_offset(&header[10..18])?;
        let text_end   = parse_hdr_offset(&header[18..26])?;
        let data_start_hdr = parse_hdr_offset(&header[26..34]).unwrap_or(0);
        let data_end_hdr   = parse_hdr_offset(&header[34..42]).unwrap_or(0);

        // ── TEXT segment ───────────────────────────────────────────────
        if text_end < text_start {
            bail!("corrupt header: TEXT end ({}) < start ({})", text_end, text_start);
        }
        if text_end >= file_len {
            bail!("corrupt header: TEXT end ({}) beyond file size ({})", text_end, file_len);
        }
        rdr.seek(SeekFrom::Start(text_start as u64))?;
        let mut text_buf = vec![0u8; text_end - text_start + 1];
        rdr.read_exact(&mut text_buf).context("failed to read TEXT segment")?;

        let mut kw = parse_text(&text_buf)?;

        // Supplemental TEXT (FCS 3.1 §3.3)
        if let (Some(ss), Some(se)) = (
            kw.get("$SUPTEXT_START").and_then(|v| v.parse::<usize>().ok()),
            kw.get("$SUPTEXT_END").and_then(|v| v.parse::<usize>().ok()),
        ) {
            if ss > 0 && se >= ss {
                rdr.seek(SeekFrom::Start(ss as u64))?;
                let mut sup = vec![0u8; se - ss + 1];
                if rdr.read_exact(&mut sup).is_ok() {
                    if let Ok(sup_kw) = parse_text(&sup) {
                        for (k, v) in sup_kw {
                            kw.entry(k).or_insert(v);
                        }
                    }
                }
            }
        }

        // ── Parameters ────────────────────────────────────────────────
        let n_params: usize = kw.get("$PAR")
            .context("$PAR keyword missing")?
            .trim().parse().context("$PAR is not an integer")?;

        let n_events: usize = kw.get("$TOT")
            .context("$TOT keyword missing")?
            .trim().parse().context("$TOT is not an integer")?;

        let mut parameters = Vec::with_capacity(n_params);
        for i in 1..=n_params {
            let name = kw.get(&format!("$P{}N", i))
                .cloned()
                .unwrap_or_else(|| format!("P{}", i));
            let label = kw.get(&format!("$P{}S", i)).cloned();
            let range: f64 = kw.get(&format!("$P{}R", i))
                .and_then(|v| v.trim().parse().ok())
                .unwrap_or(0.0);
            let bits: u32 = kw.get(&format!("$P{}B", i))
                .and_then(|v| v.trim().parse().ok())
                .unwrap_or(32);
            parameters.push(Parameter { index: i, name, label, range, bits });
        }

        // ── DATA offsets ──────────────────────────────────────────────
        // Header offsets are 0 for large files (>99,999,999 bytes); use TEXT keywords.
        let data_start = if data_start_hdr > 0 {
            data_start_hdr
        } else {
            kw.get("$BEGINDATA")
                .context("header DATA offset is 0 and $BEGINDATA is missing")?
                .trim().parse::<usize>()
                .context("$BEGINDATA parse error")?
        };
        let data_end = if data_end_hdr > 0 {
            data_end_hdr
        } else {
            kw.get("$ENDDATA")
                .context("header DATA offset is 0 and $ENDDATA is missing")?
                .trim().parse::<usize>()
                .context("$ENDDATA parse error")?
        };

        // ── DATA segment ──────────────────────────────────────────────
        let datatype = kw.get("$DATATYPE")
            .map(|s| s.trim().to_uppercase())
            .unwrap_or_else(|| "F".to_string());

        // $BYTEORD: "1,2,3,4…" (ascending) = little-endian; "4,3,2,1" = big-endian.
        let little_endian = kw.get("$BYTEORD")
            .map(|s| {
                let parts: Vec<&str> = s.trim().split(',').collect();
                // strictly ascending 1,2,3,… ⇒ little-endian (any byte width)
                parts.iter().enumerate().all(|(i, p)| p.trim() == (i + 1).to_string())
            })
            .unwrap_or(true);

        if data_end < data_start {
            bail!("corrupt header: DATA end ({}) < start ({})", data_end, data_start);
        }
        if data_end >= file_len {
            bail!("corrupt/truncated: DATA end ({}) beyond file size ({})", data_end, file_len);
        }
        rdr.seek(SeekFrom::Start(data_start as u64))?;
        let data_len = data_end - data_start + 1;
        let mut data_buf = vec![0u8; data_len];
        rdr.read_exact(&mut data_buf).context("failed to read DATA segment")?;

        let events = parse_data(
            &data_buf, n_params, n_events,
            &datatype, little_endian, &parameters,
        )?;

        Ok(FcsFile { version, keywords: kw, parameters, n_events, events })
    }
}

// ── Helpers ───────────────────────────────────────────────────────────────

fn parse_hdr_offset(bytes: &[u8]) -> Result<usize> {
    let s = std::str::from_utf8(bytes)
        .context("header offset field not UTF-8")?
        .trim();
    if s.is_empty() {
        return Ok(0);
    }
    s.parse::<usize>()
        .with_context(|| format!("bad header offset {:?}", s))
}

/// Parse FCS TEXT segment (first byte is the delimiter; keys normalised to uppercase).
/// Handles double-delimiter escaping (§3.2.3 of FCS 3.1 spec).
fn parse_text(buf: &[u8]) -> Result<HashMap<String, String>> {
    if buf.is_empty() {
        bail!("empty TEXT segment");
    }
    let delim = buf[0];

    // Tokenise: each token is separated by an un-escaped delimiter.
    let mut tokens: Vec<String> = Vec::new();
    let mut cur: Vec<u8> = Vec::new();
    let mut i = 1usize; // skip the leading delimiter

    while i < buf.len() {
        if buf[i] == delim {
            if i + 1 < buf.len() && buf[i + 1] == delim {
                // Escaped delimiter — include one literal delimiter in the token.
                cur.push(delim);
                i += 2;
            } else {
                tokens.push(String::from_utf8_lossy(&cur).trim().to_string());
                cur.clear();
                i += 1;
            }
        } else {
            cur.push(buf[i]);
            i += 1;
        }
    }
    if !cur.is_empty() {
        tokens.push(String::from_utf8_lossy(&cur).trim().to_string());
    }

    let mut map = HashMap::with_capacity(tokens.len() / 2);
    let mut j = 0;
    while j + 1 < tokens.len() {
        let key = tokens[j].to_uppercase();
        let val = tokens[j + 1].clone();
        if !key.is_empty() {
            map.insert(key, val);
        }
        j += 2;
    }
    Ok(map)
}

/// Parse DATA segment into a flat Vec<f64> (row-major: event × param).
fn parse_data(
    buf: &[u8],
    n_params: usize,
    n_events: usize,
    datatype: &str,
    little_endian: bool,
    params: &[Parameter],
) -> Result<Vec<f64>> {
    let total = n_params * n_events;
    let mut out = Vec::with_capacity(total);
    let mut c = Cursor::new(buf);

    match datatype {
        "F" => {
            for _ in 0..total {
                let v = if little_endian {
                    c.read_f32::<LittleEndian>()?
                } else {
                    c.read_f32::<BigEndian>()?
                };
                out.push(v as f64);
            }
        }
        "D" => {
            for _ in 0..total {
                let v = if little_endian {
                    c.read_f64::<LittleEndian>()?
                } else {
                    c.read_f64::<BigEndian>()?
                };
                out.push(v);
            }
        }
        "I" => {
            // Each parameter may have a different bit-width ($PnB) and range ($PnR).
            for _event_idx in 0..n_events {
                for p in params {
                    let raw: u64 = match p.bits {
                        8 => c.read_u8()? as u64,
                        16 => if little_endian {
                            c.read_u16::<LittleEndian>()? as u64
                        } else {
                            c.read_u16::<BigEndian>()? as u64
                        },
                        32 => if little_endian {
                            c.read_u32::<LittleEndian>()? as u64
                        } else {
                            c.read_u32::<BigEndian>()? as u64
                        },
                        64 => if little_endian {
                            c.read_u64::<LittleEndian>()?
                        } else {
                            c.read_u64::<BigEndian>()?
                        },
                        b => bail!("unsupported $P{}B = {}", p.index, b),
                    };
                    // Apply range bit-mask (common on older FACSCalibur 12-bit data).
                    // The mask is the smallest 2ⁿ−1 that covers $PnR, so it works for
                    // non-power-of-2 ranges too (range-1 alone would corrupt those).
                    let mask = if p.range >= 1.0 && p.range.fract() == 0.0 {
                        (p.range as u64).next_power_of_two() - 1
                    } else {
                        u64::MAX
                    };
                    out.push((raw & mask) as f64);
                }
            }
        }
        "A" => bail!("ASCII ($DATATYPE=A) FCS files are not supported"),
        t => bail!("unknown $DATATYPE: {}", t),
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn p(index: usize, name: &str, bits: u32, range: f64) -> Parameter {
        Parameter { index, name: name.to_string(), label: None, range, bits }
    }

    // ── parse_text ─────────────────────────────────────────────────────────

    #[test]
    fn parse_text_basic() {
        // First byte is the delimiter ('/').
        let buf = b"/$PAR/2/$TOT/5/";
        let kw = parse_text(buf).unwrap();
        assert_eq!(kw.get("$PAR").unwrap(), "2");
        assert_eq!(kw.get("$TOT").unwrap(), "5");
    }

    #[test]
    fn parse_text_uppercases_keys() {
        let buf = b"/$par/9/";
        let kw = parse_text(buf).unwrap();
        assert_eq!(kw.get("$PAR").unwrap(), "9");
    }

    #[test]
    fn parse_text_double_delimiter_escape() {
        // A doubled delimiter inside a value is a literal delimiter (§3.2.3).
        // value "a/b" is encoded as "a//b".
        let buf = b"/KEY/a//b/";
        let kw = parse_text(buf).unwrap();
        assert_eq!(kw.get("KEY").unwrap(), "a/b");
    }

    #[test]
    fn parse_text_empty_errors() {
        assert!(parse_text(b"").is_err());
    }

    // ── parse_hdr_offset ───────────────────────────────────────────────────

    #[test]
    fn hdr_offset_parses_and_trims() {
        assert_eq!(parse_hdr_offset(b"      58").unwrap(), 58);
        assert_eq!(parse_hdr_offset(b"        ").unwrap(), 0); // blank → 0
        assert!(parse_hdr_offset(b"   12x  ").is_err());
    }

    // ── parse_data: float ──────────────────────────────────────────────────

    fn f32le(vals: &[f64]) -> Vec<u8> {
        vals.iter().flat_map(|&v| (v as f32).to_le_bytes()).collect()
    }
    fn f32be(vals: &[f64]) -> Vec<u8> {
        vals.iter().flat_map(|&v| (v as f32).to_be_bytes()).collect()
    }

    #[test]
    fn parse_data_float_little_endian() {
        let params = [p(1, "A", 32, 0.0), p(2, "B", 32, 0.0)];
        let buf = f32le(&[1.0, 2.0, 3.0, 4.0]);
        let out = parse_data(&buf, 2, 2, "F", true, &params).unwrap();
        assert_eq!(out, vec![1.0, 2.0, 3.0, 4.0]);
    }

    #[test]
    fn parse_data_float_big_endian() {
        let params = [p(1, "A", 32, 0.0), p(2, "B", 32, 0.0)];
        let buf = f32be(&[1.0, 2.0, 3.0, 4.0]);
        let out = parse_data(&buf, 2, 2, "F", false, &params).unwrap();
        assert_eq!(out, vec![1.0, 2.0, 3.0, 4.0]);
    }

    #[test]
    fn parse_data_double() {
        let params = [p(1, "A", 64, 0.0)];
        let buf: Vec<u8> = [1.5f64, 2.5].iter().flat_map(|v| v.to_le_bytes()).collect();
        let out = parse_data(&buf, 1, 2, "D", true, &params).unwrap();
        assert_eq!(out, vec![1.5, 2.5]);
    }

    #[test]
    fn parse_data_integer_16bit() {
        // 16-bit unsigned LE, range a power of two so the mask is transparent.
        let params = [p(1, "A", 16, 65536.0)];
        let buf: Vec<u8> = [300u16, 1000].iter().flat_map(|v| v.to_le_bytes()).collect();
        let out = parse_data(&buf, 1, 2, "I", true, &params).unwrap();
        assert_eq!(out, vec![300.0, 1000.0]);
    }

    #[test]
    fn parse_data_integer_range_bitmask() {
        // Classic 12-bit FACSCalibur data: $PnR = 1024 masks to 10 bits (0..1023).
        // 5000 & 1023 = 904.
        let params = [p(1, "A", 16, 1024.0)];
        let buf: Vec<u8> = [5000u16].iter().flat_map(|v| v.to_le_bytes()).collect();
        let out = parse_data(&buf, 1, 1, "I", true, &params).unwrap();
        assert_eq!(out, vec![904.0]);
    }

    #[test]
    fn parse_data_unsupported_datatype_errors() {
        let params = [p(1, "A", 32, 0.0)];
        assert!(parse_data(b"\0\0\0\0", 1, 1, "A", true, &params).is_err());
        assert!(parse_data(b"\0\0\0\0", 1, 1, "Z", true, &params).is_err());
    }

    // ── full file: build → open round-trip ─────────────────────────────────

    /// Assemble a complete FCS3.0 byte image (header + TEXT + DATA).
    fn build_fcs(text_kw: &[(&str, &str)], data: &[u8]) -> Vec<u8> {
        let delim = b'/';
        let mut text = vec![delim];
        for (k, v) in text_kw {
            text.extend_from_slice(k.as_bytes());
            text.push(delim);
            text.extend_from_slice(v.as_bytes());
            text.push(delim);
        }
        let text_start = 58usize;
        let text_end = text_start + text.len() - 1;
        let data_start = text_end + 1;
        let data_end = data_start + data.len() - 1;

        let mut hdr = vec![b' '; 58];
        hdr[0..6].copy_from_slice(b"FCS3.0");
        let mut put = |range: std::ops::Range<usize>, v: usize| {
            let s = format!("{:>8}", v);
            hdr[range].copy_from_slice(s.as_bytes());
        };
        put(10..18, text_start);
        put(18..26, text_end);
        put(26..34, data_start);
        put(34..42, data_end);
        put(42..50, 0);
        put(50..58, 0);

        let mut out = hdr;
        out.extend_from_slice(&text);
        out.extend_from_slice(data);
        out
    }

    fn write_temp(bytes: &[u8], tag: &str) -> std::path::PathBuf {
        let path = std::env::temp_dir()
            .join(format!("flowcyto_fcs_{}_{}.fcs", std::process::id(), tag));
        std::fs::write(&path, bytes).unwrap();
        path
    }

    fn standard_kw() -> Vec<(&'static str, &'static str)> {
        vec![
            ("$PAR", "2"),
            ("$TOT", "2"),
            ("$DATATYPE", "F"),
            ("$BYTEORD", "1,2,3,4"),
            ("$MODE", "L"),
            ("$P1N", "FSC-A"),
            ("$P1B", "32"),
            ("$P1R", "262144"),
            ("$P2N", "FITC-A"),
            ("$P2B", "32"),
            ("$P2R", "262144"),
            ("$SPILLOVER", "1,FITC-A,1"),
        ]
    }

    #[test]
    fn open_round_trip() {
        let data = f32le(&[10.0, 100.0, 20.0, 200.0]); // 2 events × 2 params
        let bytes = build_fcs(&standard_kw(), &data);
        let path = write_temp(&bytes, "rt");

        let fcs = FcsFile::open(&path).unwrap();
        assert_eq!(fcs.version, "FCS3.0");
        assert_eq!(fcs.n_params(), 2);
        assert_eq!(fcs.n_events, 2);
        assert_eq!(fcs.events, vec![10.0, 100.0, 20.0, 200.0]);
        assert_eq!(fcs.param_index("fitc-a"), Some(1)); // case-insensitive
        assert_eq!(fcs.channel_values(1), vec![100.0, 200.0]);
        assert_eq!(fcs.event_slice(0), &[10.0, 100.0]);
        assert_eq!(fcs.spillover_keyword(), Some("1,FITC-A,1"));

        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn peek_events_reads_tot_only() {
        let data = f32le(&[10.0, 100.0, 20.0, 200.0]);
        let bytes = build_fcs(&standard_kw(), &data);
        let path = write_temp(&bytes, "peek");
        assert_eq!(FcsFile::peek_events(&path).unwrap(), 2);
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn open_rejects_non_fcs_magic() {
        let path = write_temp(b"NOTANFCSFILE........................................................", "magic");
        assert!(FcsFile::open(&path).is_err());
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn open_rejects_corrupt_text_offsets() {
        // text_end < text_start triggers the corrupt-header bail.
        let data = f32le(&[1.0, 2.0, 3.0, 4.0]);
        let mut bytes = build_fcs(&standard_kw(), &data);
        // Overwrite text_end (bytes 18..26) with a value smaller than text_start.
        bytes[18..26].copy_from_slice(b"      10");
        let path = write_temp(&bytes, "corrupt");
        assert!(FcsFile::open(&path).is_err());
        let _ = std::fs::remove_file(&path);
    }
}
