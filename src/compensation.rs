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
    let expected = 1 + n + n * n;
    if parts.len() < expected {
        bail!(
            "spillover has too few tokens: expected {}, got {} (n={})",
            expected, parts.len(), n
        );
    }
    let channels: Vec<String> = parts[1..=n]
        .iter().map(|s| s.trim().to_string()).collect();

    let flat: Vec<f64> = parts[1 + n..1 + n + n * n]
        .iter()
        .map(|s| s.trim().parse::<f64>().unwrap_or(0.0))
        .collect();

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
            let mut vals = fcs.channel_values(idx);
            Ok(median(&mut vals))
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
fn median(vals: &mut [f64]) -> f64 {
    let n = vals.len();
    if n == 0 {
        return 0.0;
    }
    vals.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    if n % 2 == 1 {
        vals[n / 2]
    } else {
        0.5 * (vals[n / 2 - 1] + vals[n / 2])
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
        let m = DMatrix::from_row_slice(n, n, &flat);
        let inv = m.try_inverse()
            .context("spillover matrix is singular (cannot invert); check for linearly dependent channels")?;
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
