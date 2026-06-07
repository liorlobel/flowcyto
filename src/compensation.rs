use anyhow::{bail, Context, Result};
use nalgebra::DMatrix;

use crate::fcs::FcsFile;

pub struct SpilloverMatrix {
    pub channels: Vec<String>,
    /// Square matrix, row-major storage (nalgebra is column-major internally but indexed [row,col]).
    /// Convention matches flowCore: compensated = raw_row_vec × M⁻¹
    pub inv: DMatrix<f64>,
}

/// Parse the value of $SPILLOVER / SPILL into (channel names, matrix rows).
///   n,name1,...,nameN,v11,v12,...,vNN
/// where v[i][j] = fraction of fluorochrome i detected in channel j.
///
/// Does NOT invert — works even on a singular matrix, so it is safe for display.
pub fn parse_spillover(kw: &str) -> Result<(Vec<String>, Vec<Vec<f64>>)> {
    let parts: Vec<&str> = kw.split(',').collect();
    if parts.is_empty() {
        bail!("spillover keyword is empty");
    }
    let n: usize = parts[0].trim().parse()
        .context("spillover: first token (n) is not an integer")?;
    // Compute the expected token count with checked arithmetic: a crafted huge `n`
    // would otherwise overflow `1 + n + n*n` (wrapping it small in release) and then
    // slice out of bounds below.
    let expected = n.checked_mul(n)
        .and_then(|nn| nn.checked_add(n))
        .and_then(|s| s.checked_add(1))
        .context("spillover: n is implausibly large")?;
    if parts.len() < expected {
        bail!(
            "spillover has too few tokens: expected {}, got {} (n={})",
            expected, parts.len(), n
        );
    }
    let channels: Vec<String> = parts[1..=n]
        .iter().map(|s| s.trim().to_string()).collect();

    // Parse strictly: a malformed token used to become 0.0 silently (corrupting the
    // matrix), and "nan"/"inf" parse to non-finite — both poison compensation.
    let flat: Vec<f64> = parts[1 + n..1 + n + n * n]
        .iter()
        .map(|s| {
            let v: f64 = s.trim().parse()
                .with_context(|| format!("spillover: '{}' is not a number", s.trim()))?;
            if !v.is_finite() {
                bail!("spillover contains a non-finite value '{}'", s.trim());
            }
            Ok(v)
        })
        .collect::<Result<_>>()?;

    let rows: Vec<Vec<f64>> = (0..n).map(|i| flat[i * n..(i + 1) * n].to_vec()).collect();
    Ok((channels, rows))
}

/// Format (channels, rows) into a `$SPILLOVER` keyword value:
///   n,name1,...,nameN,v11,v12,...,vNN
/// Uses Rust's shortest round-trippable float formatting.
pub fn format_spillover(channels: &[String], rows: &[Vec<f64>]) -> String {
    let n = channels.len();
    let mut parts: Vec<String> = Vec::with_capacity(1 + n + n * n);
    parts.push(n.to_string());
    for c in channels {
        parts.push(c.clone());
    }
    for row in rows {
        for &v in row {
            parts.push(format!("{}", v));
        }
    }
    parts.join(",")
}

/// Load a spillover matrix from a CSV or JSON file.
///
/// CSV: header row = (corner cell), detector names…; each data row = source name, values…
/// JSON: `{"channels": ["FITC-A", …], "matrix": [[1.0, 0.21, …], …]}`
pub fn load_matrix_file(path: &std::path::Path) -> Result<(Vec<String>, Vec<Vec<f64>>)> {
    let text = std::fs::read_to_string(path)
        .with_context(|| format!("cannot read matrix file {:?}", path))?;
    let is_json = path.extension().map(|e| e.eq_ignore_ascii_case("json")).unwrap_or(false)
        || text.trim_start().starts_with('{');

    if is_json {
        #[derive(serde::Deserialize)]
        struct M { channels: Vec<String>, matrix: Vec<Vec<f64>> }
        let m: M = serde_json::from_str(&text).context("parsing matrix JSON")?;
        validate_square(&m.channels, &m.matrix)?;
        Ok((m.channels, m.matrix))
    } else {
        let mut lines = text.lines().filter(|l| !l.trim().is_empty());
        let header = lines.next().context("matrix CSV is empty")?;
        let channels: Vec<String> = header.split(',').skip(1).map(|s| s.trim().to_string()).collect();
        let mut rows: Vec<Vec<f64>> = Vec::new();
        for (li, line) in lines.enumerate() {
            let cells: Vec<&str> = line.split(',').collect();
            let vals: Vec<f64> = cells.iter().skip(1)
                .map(|s| s.trim().parse::<f64>()
                    .with_context(|| format!("matrix row {}: '{}' is not a number", li + 1, s.trim())))
                .collect::<Result<_>>()?;
            rows.push(vals);
        }
        validate_square(&channels, &rows)?;
        Ok((channels, rows))
    }
}

/// Save a spillover matrix to CSV (default) or JSON (if path ends in .json).
pub fn save_matrix_file(path: &std::path::Path, channels: &[String], rows: &[Vec<f64>]) -> Result<()> {
    // Never persist a non-finite matrix (e.g. from a degenerate control) — it would
    // poison any later `compensate`. (Defense-in-depth alongside the median fix.)
    if rows.iter().flatten().any(|v| !v.is_finite()) {
        bail!("refusing to write a spillover matrix containing non-finite (NaN/Inf) values");
    }
    let is_json = path.extension().map(|e| e.eq_ignore_ascii_case("json")).unwrap_or(false);
    if is_json {
        let obj = serde_json::json!({ "channels": channels, "matrix": rows });
        std::fs::write(path, serde_json::to_string_pretty(&obj)?)
            .with_context(|| format!("writing {:?}", path))?;
    } else {
        let mut s = String::new();
        s.push(',');
        s.push_str(&channels.join(","));
        s.push('\n');
        for (i, row) in rows.iter().enumerate() {
            s.push_str(&channels[i]);
            for &v in row {
                s.push(',');
                s.push_str(&format!("{}", v));
            }
            s.push('\n');
        }
        std::fs::write(path, s).with_context(|| format!("writing {:?}", path))?;
    }
    Ok(())
}

fn validate_square(channels: &[String], rows: &[Vec<f64>]) -> Result<()> {
    let n = channels.len();
    if n == 0 {
        bail!("matrix has no channels");
    }
    if rows.len() != n {
        bail!("matrix is not square: {} channels but {} rows", n, rows.len());
    }
    for (i, r) in rows.iter().enumerate() {
        if r.len() != n {
            bail!("row {} has {} values, expected {}", i, r.len(), n);
        }
        for (j, v) in r.iter().enumerate() {
            if !v.is_finite() {
                bail!("matrix entry [{},{}] is not finite ({})", i, j, v);
            }
        }
    }
    Ok(())
}

/// Result of computing a spillover matrix from single-stain controls.
pub struct ComputedSpillover {
    pub channels: Vec<String>,
    pub rows: Vec<Vec<f64>>,
    /// For each input control (same order as passed), the channel index it was
    /// assigned to (its primary/brightest detector).
    pub assigned: Vec<usize>,
}

/// Compute a spillover matrix from single-stain controls + an unstained control.
///
/// Algorithm (matches the standard flowCore/FACSDiva "all-events median" method):
///   1. background[j] = median of channel j over the unstained control
///   2. for each stained control: sig[j] = median_j(control) − background[j]
///   3. primary p = argmaxⱼ sig[j]  (intensity-based stain assignment)
///   4. matrix row for channel p = sig / sig[p]   (so diagonal = 1)
///
/// Computed on RAW (uncompensated) fluorescence. Scatter/Time must be excluded
/// by the caller via `fluor_channels`.
pub fn compute_spillover(
    fluor_channels: &[String],
    unstained: &FcsFile,
    controls: &[&FcsFile],
) -> Result<ComputedSpillover> {
    let n = fluor_channels.len();
    if n == 0 {
        bail!("no fluorescence channels supplied");
    }
    if controls.len() != n {
        bail!(
            "expected {} single-stain controls (one per fluorescence channel), got {}",
            n, controls.len()
        );
    }

    let background = channel_medians(unstained, fluor_channels)
        .context("computing unstained background medians")?;

    let mut rows: Vec<Option<Vec<f64>>> = vec![None; n];
    let mut assigned: Vec<usize> = vec![0; controls.len()];

    for (ci, c) in controls.iter().enumerate() {
        let med = channel_medians(c, fluor_channels)
            .with_context(|| format!("computing medians for control #{}", ci + 1))?;
        let sig: Vec<f64> = (0..n).map(|j| med[j] - background[j]).collect();

        // primary = brightest background-subtracted channel
        let p = sig
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.partial_cmp(b.1).unwrap_or(std::cmp::Ordering::Equal))
            .map(|(i, _)| i)
            .unwrap();
        let denom = sig[p];
        if denom <= 0.0 {
            bail!(
                "control #{} has non-positive signal in its brightest channel '{}' \
                 — is it really a stained control (vs unstained)?",
                ci + 1, fluor_channels[p]
            );
        }
        if rows[p].is_some() {
            bail!(
                "two controls were both assigned to channel '{}' by intensity — \
                 check for a duplicate or swapped control file",
                fluor_channels[p]
            );
        }
        rows[p] = Some(sig.iter().map(|&v| v / denom).collect());
        assigned[ci] = p;
    }

    let out_rows: Vec<Vec<f64>> = rows
        .into_iter()
        .enumerate()
        .map(|(j, r)| r.ok_or_else(|| anyhow::anyhow!(
            "no control was assigned to channel '{}' — a single stain may be missing",
            fluor_channels[j]
        )))
        .collect::<Result<_>>()?;

    Ok(ComputedSpillover { channels: fluor_channels.to_vec(), rows: out_rows, assigned })
}

fn channel_medians(fcs: &FcsFile, channels: &[String]) -> Result<Vec<f64>> {
    channels
        .iter()
        .map(|ch| {
            let idx = fcs.param_index(ch)
                .with_context(|| format!("channel '{}' not found in control file", ch))?;
            let vals = fcs.channel_values(idx);
            Ok(median(&vals))
        })
        .collect()
}

/// Longest fluorescence-channel token (channel name minus the "-A" suffix) that is a
/// substring of `filename`. Longest-match avoids the "PE" ⊂ "PE-Cy7" trap.
/// Used to cross-check intensity-based stain assignment against control filenames.
pub fn fluor_token_in_filename(channels: &[String], filename: &str) -> Option<usize> {
    let mut best: Option<(usize, usize)> = None;
    for (i, ch) in channels.iter().enumerate() {
        let token = ch.strip_suffix("-A").unwrap_or(ch);
        if filename.contains(token) {
            let len = token.len();
            if best.map(|(_, l)| len > l).unwrap_or(true) {
                best = Some((i, len));
            }
        }
    }
    best.map(|(i, _)| i)
}

/// Median (averages the two middle values for even N, matching R's `median`).
/// Non-finite values are dropped first so a NaN/Inf raw event in a control can't
/// poison the computed spillover (mirrors the popstats/stats medians).
fn median(vals: &[f64]) -> f64 {
    let mut v: Vec<f64> = vals.iter().copied().filter(|x| x.is_finite()).collect();
    let n = v.len();
    if n == 0 {
        return 0.0;
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    if n % 2 == 1 {
        v[n / 2]
    } else {
        0.5 * (v[n / 2 - 1] + v[n / 2])
    }
}

/// Largest off-diagonal spillover value (0 = identity / no compensation).
pub fn max_off_diagonal(rows: &[Vec<f64>]) -> f64 {
    let mut mx = 0.0f64;
    for (i, row) in rows.iter().enumerate() {
        for (j, &v) in row.iter().enumerate() {
            if i != j && v.abs() > mx { mx = v.abs(); }
        }
    }
    mx
}

impl SpilloverMatrix {
    pub fn from_keyword(kw: &str) -> Result<Self> {
        let (channels, rows) = parse_spillover(kw)?;
        Self::from_parts(channels, &rows)
    }

    /// Build a compensator from an explicit (channels, rows) matrix — used for
    /// user-supplied / edited override matrices.
    pub fn from_parts(channels: Vec<String>, rows: &[Vec<f64>]) -> Result<Self> {
        let n = channels.len();
        let flat: Vec<f64> = rows.iter().flatten().copied().collect();
        if flat.len() != n * n {
            bail!("spillover matrix is not {n}×{n}");
        }
        if !flat.iter().all(|v| v.is_finite()) {
            bail!("spillover matrix contains non-finite (NaN/Inf) values");
        }
        let m = DMatrix::from_row_slice(n, n, &flat);
        let inv = m.try_inverse()
            .context("spillover matrix is singular (cannot invert); check for linearly dependent channels")?;
        // Reject ill-conditioned matrices: near-linearly-dependent channels (e.g. heavy
        // tandem-dye overlap) invert to enormous entries that would silently explode
        // compensated values. nalgebra only returns None at ~exact singularity, so a
        // determinant of ~1e-13 inverts to ~1e13 with no error otherwise. A well-formed
        // compensation inverse has small entries.
        let max_inv = inv.iter().fold(0.0_f64, |mx, &x| mx.max(x.abs()));
        if !max_inv.is_finite() || max_inv > 1.0e6 {
            bail!("spillover matrix is ill-conditioned (inverse magnitude {:.3e}); \
                   check for near-linearly-dependent channels", max_inv);
        }
        Ok(SpilloverMatrix { channels, inv })
    }

    /// Apply compensation to all events in `fcs`, returning a new flat event vector.
    ///
    /// Only the fluorescence channels listed in the spillover matrix are modified;
    /// scatter (FSC, SSC) and Time channels are left untouched.
    ///
    /// Math (matching flowCore::compensate):
    ///   compensated_row = raw_row × M⁻¹
    ///   i.e.  comp[j] = Σᵢ raw[i] × inv[i,j]
    pub fn apply(&self, fcs: &FcsFile) -> Result<Vec<f64>> {
        let n_params = fcs.n_params();
        let n_spill  = self.channels.len();

        // Map each spillover channel name to its column index in the event matrix.
        let indices: Vec<usize> = self.channels
            .iter()
            .map(|name| {
                fcs.param_index(name).with_context(|| format!(
                    "spillover channel '{}' not found among FCS parameters: [{}]",
                    name,
                    fcs.parameters.iter().map(|p| p.name.as_str()).collect::<Vec<_>>().join(", ")
                ))
            })
            .collect::<Result<_>>()?;

        let mut events = fcs.events.clone();

        for ev in 0..fcs.n_events {
            let base = ev * n_params;
            // Extract raw fluorescence values for this event.
            let raw: Vec<f64> = indices.iter().map(|&i| events[base + i]).collect();
            // comp[j] = Σᵢ raw[i] * inv[i,j]
            for j in 0..n_spill {
                let mut val = 0.0f64;
                for i in 0..n_spill {
                    val += raw[i] * self.inv[(i, j)];
                }
                events[base + indices[j]] = val;
            }
        }
        Ok(events)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_util::make_fcs;

    // ── parse / format round-trip ──────────────────────────────────────────

    #[test]
    fn parse_spillover_basic() {
        let (chans, rows) = parse_spillover("2,FITC-A,PE-A,1,0.21,0.05,1").unwrap();
        assert_eq!(chans, vec!["FITC-A", "PE-A"]);
        assert_eq!(rows, vec![vec![1.0, 0.21], vec![0.05, 1.0]]);
    }

    #[test]
    fn parse_format_round_trip() {
        let chans = vec!["FITC-A".to_string(), "PE-A".to_string()];
        let rows = vec![vec![1.0, 0.21], vec![0.05, 1.0]];
        let kw = format_spillover(&chans, &rows);
        let (c2, r2) = parse_spillover(&kw).unwrap();
        assert_eq!(c2, chans);
        assert_eq!(r2, rows);
    }

    #[test]
    fn parse_spillover_too_few_tokens_errors() {
        // n=3 needs 1 + 3 + 9 = 13 tokens; supply far fewer.
        assert!(parse_spillover("3,A,B,C,1,0,0").is_err());
    }

    #[test]
    fn parse_spillover_bad_n_errors() {
        assert!(parse_spillover("notanint,A,1").is_err());
    }

    // ── compute_spillover: synthetic ground-truth recovery ─────────────────

    /// Build a single-stain control: brightest in channel `primary`, with a
    /// known spillover fraction `spill` into the *other* channel, on top of a
    /// constant background. Two events per file keeps the median exact.
    fn stain(primary: usize, signal: f64, spill: f64, bg: f64) -> FcsFile {
        let other = 1 - primary;
        let mut a = vec![bg, bg];
        let mut b = vec![bg, bg];
        a[primary] += signal;
        a[other] += signal * spill;
        b[primary] += signal;
        b[other] += signal * spill;
        make_fcs(&["FITC-A", "PE-A"], &[a, b])
    }

    #[test]
    fn compute_spillover_recovers_known_matrix() {
        let bg = 100.0;
        let unstained = make_fcs(&["FITC-A", "PE-A"], &[vec![bg, bg], vec![bg, bg]]);
        // FITC control spills 0.20 into PE; PE control spills 0.08 into FITC.
        let fitc = stain(0, 1000.0, 0.20, bg);
        let pe = stain(1, 1000.0, 0.08, bg);

        let res = compute_spillover(
            &["FITC-A".into(), "PE-A".into()],
            &unstained,
            &[&fitc, &pe],
        )
        .unwrap();

        assert_eq!(res.assigned, vec![0, 1], "intensity-based assignment");
        // Row 0 (FITC primary): [1, 0.20]; row 1 (PE primary): [0.08, 1].
        assert!((res.rows[0][0] - 1.0).abs() < 1e-9);
        assert!((res.rows[0][1] - 0.20).abs() < 1e-9);
        assert!((res.rows[1][0] - 0.08).abs() < 1e-9);
        assert!((res.rows[1][1] - 1.0).abs() < 1e-9);
    }

    #[test]
    fn compute_spillover_wrong_control_count_errors() {
        let unstained = make_fcs(&["FITC-A", "PE-A"], &[vec![0.0, 0.0]]);
        let fitc = stain(0, 1000.0, 0.2, 0.0);
        // Two channels but only one control.
        assert!(compute_spillover(
            &["FITC-A".into(), "PE-A".into()],
            &unstained,
            &[&fitc]
        )
        .is_err());
    }

    #[test]
    fn compute_spillover_two_controls_same_channel_errors() {
        let bg = 0.0;
        let unstained = make_fcs(&["FITC-A", "PE-A"], &[vec![bg, bg]]);
        let a = stain(0, 1000.0, 0.1, bg);
        let b = stain(0, 900.0, 0.1, bg); // also brightest in FITC
        assert!(compute_spillover(
            &["FITC-A".into(), "PE-A".into()],
            &unstained,
            &[&a, &b]
        )
        .is_err());
    }

    // ── SpilloverMatrix: invert + apply ────────────────────────────────────

    #[test]
    fn from_parts_singular_errors() {
        // Linearly dependent rows → singular → cannot invert.
        let chans = vec!["A".to_string(), "B".to_string()];
        let rows = vec![vec![1.0, 1.0], vec![1.0, 1.0]];
        assert!(SpilloverMatrix::from_parts(chans, &rows).is_err());
    }

    #[test]
    fn parse_spillover_rejects_non_finite_value() {
        // Audit M1b: "nan"/"inf" parse to non-finite and used to flow through silently.
        assert!(parse_spillover("2,A,B,1,nan,0,1").is_err());
        assert!(parse_spillover("2,A,B,1,inf,0,1").is_err());
        assert!(parse_spillover("2,A,B,1,xyz,0,1").is_err()); // malformed → no longer 0.0
    }

    #[test]
    fn parse_spillover_huge_n_errors_not_panics() {
        // Audit M5: a crafted huge n must error cleanly, not overflow + slice-panic.
        let kw = format!("{},A,B,1,0,0,1", usize::MAX);
        assert!(parse_spillover(&kw).is_err());
    }

    #[test]
    fn from_parts_rejects_non_finite() {
        let chans = vec!["A".to_string(), "B".to_string()];
        assert!(SpilloverMatrix::from_parts(chans, &[vec![1.0, f64::NAN], vec![0.0, 1.0]]).is_err());
    }

    #[test]
    fn from_parts_rejects_near_singular() {
        // Audit M1a: determinant ~1e-13 → inverse entries ~1e13 → must be rejected,
        // not silently used (which would explode compensated values).
        let chans = vec!["A".to_string(), "B".to_string()];
        let rows = vec![vec![1.0, 1.0], vec![1.0, 1.0 + 1e-13]];
        assert!(SpilloverMatrix::from_parts(chans, &rows).is_err());
    }

    #[test]
    fn validate_square_rejects_non_finite() {
        let chans = vec!["A".to_string(), "B".to_string()];
        let rows = vec![vec![1.0, 0.0], vec![f64::INFINITY, 1.0]];
        assert!(validate_square(&chans, &rows).is_err());
    }

    #[test]
    fn apply_identity_is_noop() {
        let fcs = make_fcs(&["FITC-A", "PE-A"], &[vec![10.0, 20.0], vec![30.0, 40.0]]);
        let m = SpilloverMatrix::from_parts(
            vec!["FITC-A".into(), "PE-A".into()],
            &[vec![1.0, 0.0], vec![0.0, 1.0]],
        )
        .unwrap();
        let out = m.apply(&fcs).unwrap();
        assert_eq!(out, fcs.events);
    }

    #[test]
    fn apply_inverts_known_spillover() {
        // Forward-mix raw signal through M, then compensate back to the truth.
        // flowCore convention: comp = raw × M⁻¹, and observed raw = true × M.
        let truth = [1000.0f64, 0.0];
        let s = [[1.0, 0.20], [0.08, 1.0]];
        // observed[j] = Σ_i true[i] * s[i][j]
        let obs0 = truth[0] * s[0][0] + truth[1] * s[1][0];
        let obs1 = truth[0] * s[0][1] + truth[1] * s[1][1];
        let fcs = make_fcs(&["FITC-A", "PE-A"], &[vec![obs0, obs1]]);

        let m = SpilloverMatrix::from_parts(
            vec!["FITC-A".into(), "PE-A".into()],
            &[s[0].to_vec(), s[1].to_vec()],
        )
        .unwrap();
        let out = m.apply(&fcs).unwrap();
        assert!((out[0] - truth[0]).abs() < 1e-6, "got {}", out[0]);
        assert!((out[1] - truth[1]).abs() < 1e-6, "got {}", out[1]);
    }

    #[test]
    fn apply_leaves_unlisted_channels_untouched() {
        // FSC-A is not in the spillover matrix; it must pass through unchanged.
        let fcs = make_fcs(
            &["FSC-A", "FITC-A", "PE-A"],
            &[vec![500.0, 10.0, 20.0]],
        );
        let m = SpilloverMatrix::from_parts(
            vec!["FITC-A".into(), "PE-A".into()],
            &[vec![1.0, 0.5], vec![0.0, 1.0]],
        )
        .unwrap();
        let out = m.apply(&fcs).unwrap();
        assert_eq!(out[0], 500.0, "FSC-A must be untouched");
    }

    #[test]
    fn apply_unknown_channel_errors() {
        let fcs = make_fcs(&["FITC-A"], &[vec![1.0]]);
        let m = SpilloverMatrix::from_parts(
            vec!["NONEXISTENT".into()],
            &[vec![1.0]],
        )
        .unwrap();
        assert!(m.apply(&fcs).is_err());
    }

    // ── matrix file IO ─────────────────────────────────────────────────────

    #[test]
    fn matrix_file_csv_round_trip() {
        let dir = std::env::temp_dir();
        let path = dir.join(format!("flowcyto_test_{}.csv", std::process::id()));
        let chans = vec!["FITC-A".to_string(), "PE-A".to_string()];
        let rows = vec![vec![1.0, 0.21], vec![0.05, 1.0]];
        save_matrix_file(&path, &chans, &rows).unwrap();
        let (c2, r2) = load_matrix_file(&path).unwrap();
        assert_eq!(c2, chans);
        assert_eq!(r2, rows);
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn matrix_file_json_round_trip() {
        let dir = std::env::temp_dir();
        let path = dir.join(format!("flowcyto_test_{}.json", std::process::id()));
        let chans = vec!["FITC-A".to_string(), "PE-A".to_string()];
        let rows = vec![vec![1.0, 0.21], vec![0.05, 1.0]];
        save_matrix_file(&path, &chans, &rows).unwrap();
        let (c2, r2) = load_matrix_file(&path).unwrap();
        assert_eq!(c2, chans);
        assert_eq!(r2, rows);
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn validate_square_rejects_non_square() {
        let chans = vec!["A".to_string(), "B".to_string()];
        let rows = vec![vec![1.0, 0.0]]; // only 1 row for 2 channels
        assert!(validate_square(&chans, &rows).is_err());
    }

    // ── small helpers ──────────────────────────────────────────────────────

    #[test]
    fn max_off_diagonal_ignores_diagonal() {
        let rows = vec![vec![1.0, 0.3], vec![-0.4, 1.0]];
        assert!((max_off_diagonal(&rows) - 0.4).abs() < 1e-12);
    }

    #[test]
    fn fluor_token_longest_match_wins() {
        let chans = vec!["PE-A".to_string(), "PE-Cy7-A".to_string()];
        // Filename contains both "PE" and "PE-Cy7"; longest must win → index 1.
        let idx = fluor_token_in_filename(&chans, "Tube_PE-Cy7_control.fcs");
        assert_eq!(idx, Some(1));
    }

    #[test]
    fn fluor_token_none_when_absent() {
        let chans = vec!["FITC-A".to_string()];
        assert_eq!(fluor_token_in_filename(&chans, "unstained.fcs"), None);
    }

    #[test]
    fn median_even_and_odd() {
        assert_eq!(median(&[3.0, 1.0, 2.0]), 2.0);
        assert_eq!(median(&[4.0, 1.0, 3.0, 2.0]), 2.5);
        assert_eq!(median(&[] as &[f64]), 0.0);
        // Non-finite values are dropped (audit N1): median over the finite {2,4,6}.
        assert_eq!(median(&[2.0, f64::NAN, 4.0, f64::INFINITY, 6.0]), 4.0);
    }
}
