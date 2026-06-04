//! Minimal FCS 3.0 writer.
//!
//! Writes a new FCS file containing the original RAW events plus a (possibly
//! overridden) `$SPILLOVER` keyword. Data is written as `$DATATYPE=F`
//! (little-endian float32), `$MODE=L`. The original file is never modified.
//!
//! Round-trip validated against flowCore (see `flowcyto rewrite-spillover`).

use std::io::Write;
use std::path::Path;

use anyhow::{bail, Context, Result};

use crate::fcs::FcsFile;

/// Fixed width for the $BEGINDATA/$ENDDATA values inside TEXT, so patching the
/// real offsets in does not change the TEXT length (avoids the offset chicken-egg).
const OFFW: usize = 12;

/// Write `orig` to `out_path` as a new FCS 3.0 file.
/// `new_spillover` is a `$SPILLOVER` *value* string; if `None`, the original's
/// embedded spillover (if any) is preserved.
pub fn write_fcs(orig: &FcsFile, new_spillover: Option<&str>, out_path: &Path) -> Result<()> {
    let n_params = orig.n_params();
    let n_events = orig.n_events;
    if n_params == 0 || n_events == 0 {
        bail!("cannot write an FCS file with 0 events or 0 parameters");
    }

    // ── Assemble keyword list ─────────────────────────────────────────
    // Controlled keys we always (re)write ourselves:
    let mut kws: Vec<(String, String)> = vec![
        ("$BEGINANALYSIS".into(), "0".into()),
        ("$ENDANALYSIS".into(),   "0".into()),
        ("$BEGINSTEXT".into(),    "0".into()),
        ("$ENDSTEXT".into(),      "0".into()),
        ("$NEXTDATA".into(),      "0".into()),
        ("$BYTEORD".into(),       "1,2,3,4".into()),
        ("$DATATYPE".into(),      "F".into()),
        ("$MODE".into(),          "L".into()),
        ("$PAR".into(),           n_params.to_string()),
        ("$TOT".into(),           n_events.to_string()),
    ];
    let begindata_idx = kws.len();
    kws.push(("$BEGINDATA".into(), format!("{:0width$}", 0, width = OFFW)));
    let enddata_idx = kws.len();
    kws.push(("$ENDDATA".into(), format!("{:0width$}", 0, width = OFFW)));

    // Per-parameter keywords for float data.
    for i in 1..=n_params {
        kws.push((format!("$P{}B", i), "32".to_string()));
        kws.push((format!("$P{}E", i), "0,0".to_string()));
    }

    // Spillover (override or preserved).
    let spill_val = match new_spillover {
        Some(s) => Some(s.to_string()),
        None => orig.spillover_keyword().map(|s| s.to_string()),
    };
    if let Some(sv) = spill_val {
        // Write under BOTH the FCS-standard key and the BD/flowCore convention key,
        // so flowCore's `keyword(ff)$SPILL`, FlowJo, and standards-compliant readers
        // all find it.
        kws.push(("$SPILLOVER".into(), sv.clone()));
        kws.push(("SPILL".into(), sv));
    }

    // Copy through all other original keywords (preserves $PnN/$PnS/$PnR/$CYT/$DATE…).
    // Skip empty values: an empty value would serialize as two adjacent delimiters,
    // which a reader (ours included) interprets as an escaped delimiter — corrupting
    // every subsequent key/value on round-trip.
    for (k, v) in &orig.keywords {
        if is_controlled(k) || k.is_empty() || v.is_empty() {
            continue;
        }
        kws.push((k.clone(), v.clone()));
    }

    // ── Choose a delimiter not present in any key or value ────────────
    let delim = pick_delimiter(&kws)
        .context("could not find a usable TEXT delimiter (all candidates occur in the data)")?;

    // ── Serialize TEXT (with placeholder offsets) to measure length ───
    let text = serialize_text(&kws, delim);
    let text_len = text.len();

    let text_start = 58usize;
    let text_end = text_start + text_len - 1;
    let data_start = text_end + 1;
    let data_len = n_events * n_params * 4;
    let data_end = data_start + data_len - 1;

    // Patch real offsets (same width → identical TEXT length).
    kws[begindata_idx].1 = format!("{:0width$}", data_start, width = OFFW);
    kws[enddata_idx].1 = format!("{:0width$}", data_end, width = OFFW);
    let text = serialize_text(&kws, delim);
    debug_assert_eq!(text.len(), text_len, "TEXT length changed after offset patch");
    if text.len() != text_len {
        bail!("internal error: TEXT length unstable");
    }

    // ── HEADER (58 bytes) ─────────────────────────────────────────────
    let mut header = [b' '; 58];
    header[0..6].copy_from_slice(b"FCS3.0");
    write_hdr_offset(&mut header[10..18], text_start);
    write_hdr_offset(&mut header[18..26], text_end);
    write_hdr_offset(&mut header[26..34], data_start);
    write_hdr_offset(&mut header[34..42], data_end);
    write_hdr_offset(&mut header[42..50], 0); // analysis start
    write_hdr_offset(&mut header[50..58], 0); // analysis end

    // ── DATA (float32 LE, list mode: event-major) ─────────────────────
    let mut data = Vec::with_capacity(data_len);
    for &v in &orig.events {
        data.extend_from_slice(&(v as f32).to_le_bytes());
    }
    if data.len() != data_len {
        bail!("internal error: data length {} != expected {}", data.len(), data_len);
    }

    // ── Write file ────────────────────────────────────────────────────
    let f = std::fs::File::create(out_path)
        .with_context(|| format!("cannot create {:?}", out_path))?;
    let mut w = std::io::BufWriter::new(f);
    w.write_all(&header)?;
    w.write_all(&text)?;
    w.write_all(&data)?;
    w.flush()?;
    Ok(())
}

// ── Helpers ─────────────────────────────────────────────────────────────────

fn is_controlled(k: &str) -> bool {
    matches!(
        k,
        "$BEGINANALYSIS" | "$ENDANALYSIS" | "$BEGINSTEXT" | "$ENDSTEXT" | "$NEXTDATA"
            | "$BYTEORD" | "$DATATYPE" | "$MODE" | "$PAR" | "$TOT"
            | "$BEGINDATA" | "$ENDDATA"
            | "$SPILLOVER" | "SPILLOVER" | "$SPILL" | "SPILL"
    ) || is_pnb_or_pne(k)
}

/// True for `$PnB` / `$PnE` keywords (we rewrite these for float data).
fn is_pnb_or_pne(k: &str) -> bool {
    if let Some(rest) = k.strip_prefix("$P") {
        if let Some(suffix) = rest.strip_suffix('B').or_else(|| rest.strip_suffix('E')) {
            return !suffix.is_empty() && suffix.bytes().all(|b| b.is_ascii_digit());
        }
    }
    false
}

fn serialize_text(kws: &[(String, String)], delim: u8) -> Vec<u8> {
    let mut t = Vec::new();
    t.push(delim);
    for (k, v) in kws {
        t.extend_from_slice(k.as_bytes());
        t.push(delim);
        t.extend_from_slice(v.as_bytes());
        t.push(delim);
    }
    t
}

/// Pick a delimiter byte that appears in no key or value (so no escaping needed).
fn pick_delimiter(kws: &[(String, String)]) -> Option<u8> {
    const CANDIDATES: [u8; 6] = [12, b'|', b'/', b'\\', 9, 30]; // FF, pipe, slash, backslash, tab, RS
    for &c in &CANDIDATES {
        let clash = kws.iter().any(|(k, v)| {
            k.as_bytes().contains(&c) || v.as_bytes().contains(&c)
        });
        if !clash {
            return Some(c);
        }
    }
    None
}

/// Write a usize into an 8-byte header field (right-justified ASCII).
/// Values too large for 8 digits are written as 0 (reader falls back to $BEGINDATA/$ENDDATA).
fn write_hdr_offset(field: &mut [u8], v: usize) {
    debug_assert_eq!(field.len(), 8);
    let s = if v <= 99_999_999 {
        format!("{:>8}", v)
    } else {
        format!("{:>8}", 0)
    };
    field.copy_from_slice(s.as_bytes());
}
