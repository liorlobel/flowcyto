//! `flowcyto selftest` — recompute the flowCore-validated numeric layers on a frozen
//! reference and compare to frozen flowCore golden values. It doubles as the benchmark
//! table and a CI-gated unit test, turning "validated against flowCore" into a claim
//! anyone can re-run offline.
//!
//! The reference FCS and golden values are embedded (`validation/`). Regenerate the
//! golden with `validation/gen_golden.R` (needs R + Bioconductor flowCore) whenever the
//! reference or the validated numerics change.

use crate::compensation::SpilloverMatrix;
use crate::fcs::FcsFile;
use crate::transform::AxisTransform;

const REFERENCE_FCS: &[u8] = include_bytes!("../validation/reference.fcs");
const GOLDEN_CSV: &str = include_str!("../validation/golden.csv");
const COFACTOR: f64 = 150.0;
/// Relative tolerance — the loosest the layers are validated to; actual deviations are
/// far smaller (reported in the table), bounded mostly by the golden's 10-sig-fig text.
const TOL: f64 = 1e-5;

pub struct LayerResult {
    pub layer: String,
    pub probes: usize,
    pub max_dev: f64, // max relative deviation
    pub pass: bool,
}

/// Recompute every probed quantity in flowcyto and compare to the flowCore golden.
pub fn run() -> Result<Vec<LayerResult>, String> {
    // Parse the embedded reference (via a unique temp file — FcsFile::open is path-based).
    let tmp = std::env::temp_dir().join(format!("flowcyto_selftest_{}.fcs", std::process::id()));
    std::fs::write(&tmp, REFERENCE_FCS).map_err(|e| e.to_string())?;
    let parsed = FcsFile::open(&tmp).map_err(|e| e.to_string());
    let _ = std::fs::remove_file(&tmp);
    let fcs = parsed?;
    let np = fcs.n_params();

    // flowcyto's compensation (M⁻¹ applied to the embedded $SPILLOVER).
    let kw = fcs.spillover_keyword().ok_or("reference has no $SPILLOVER")?;
    let comp = SpilloverMatrix::from_keyword(kw)
        .and_then(|m| m.apply(&fcs))
        .map_err(|e| e.to_string())?;

    let asinh = AxisTransform::Asinh { cofactor: COFACTOR }.compile();
    let logicle = AxisTransform::Logicle { t: 262144.0, w: 0.5, m: 4.5, a: 0.0 }.compile();

    // layer -> (probe count, max relative deviation)
    let mut acc: std::collections::BTreeMap<&str, (usize, f64)> = std::collections::BTreeMap::new();
    for line in GOLDEN_CSV.lines().skip(1) {
        let f: Vec<&str> = line.split(',').collect();
        if f.len() != 4 {
            continue;
        }
        let (kind, channel, key) = (f[0], f[1], f[2]);
        let golden: f64 = f[3].parse().map_err(|_| format!("bad golden {:?}", f[3]))?;
        let computed = match kind {
            "asinh" | "logicle" => {
                let v: f64 = key.parse().map_err(|_| format!("bad value {key}"))?;
                if kind == "asinh" { asinh.forward(v) } else { logicle.forward(v) }
            }
            "parse" | "comp" => {
                let ci = fcs.parameters.iter().position(|p| p.name == channel)
                    .ok_or_else(|| format!("channel {channel} not found"))?;
                let e: usize = key.parse().map_err(|_| format!("bad index {key}"))?;
                let idx = e * np + ci;
                if idx >= fcs.events.len() {
                    return Err(format!("probe index {idx} out of range"));
                }
                if kind == "parse" { fcs.events[idx] } else { comp[idx] }
            }
            _ => continue,
        };
        let rel = (computed - golden).abs() / (golden.abs() + 1.0);
        let entry = acc.entry(layer_name(kind)).or_insert((0, 0.0));
        entry.0 += 1;
        entry.1 = entry.1.max(rel);
    }
    if acc.is_empty() {
        return Err("no golden values loaded".into());
    }
    Ok(acc.into_iter()
        .map(|(layer, (probes, max_dev))| LayerResult {
            layer: layer.to_string(), probes, max_dev, pass: max_dev <= TOL,
        })
        .collect())
}

fn layer_name(kind: &str) -> &'static str {
    match kind {
        "parse" => "Parsing",
        "comp" => "Compensation",
        "asinh" => "Asinh transform",
        "logicle" => "Logicle transform",
        _ => "Other",
    }
}

/// Print the benchmark table; returns true iff every layer is within tolerance.
pub fn report(results: &[LayerResult]) -> bool {
    println!("flowcyto {} selftest — numeric layers vs frozen flowCore 2.24.0 golden values\n",
        env!("CARGO_PKG_VERSION"));
    println!("{:<20} {:>7} {:>15} {:>11} {:>8}", "Layer", "Probes", "Max rel. dev", "Tolerance", "Result");
    println!("{}", "-".repeat(64));
    let mut all = true;
    for r in results {
        all &= r.pass;
        println!("{:<20} {:>7} {:>15.2e} {:>11.0e} {:>8}",
            r.layer, r.probes, r.max_dev, TOL, if r.pass { "PASS" } else { "FAIL" });
    }
    println!();
    println!("{}", if all {
        "All numeric layers reproduce flowCore within tolerance."
    } else {
        "FAILED: a layer deviates from flowCore — see above."
    });
    all
}

#[cfg(test)]
mod tests {
    #[test]
    fn numeric_layers_match_flowcore() {
        let results = super::run().expect("selftest runs");
        assert!(!results.is_empty(), "golden values present");
        for r in &results {
            assert!(r.pass, "layer {} deviates from flowCore: max rel dev {:.2e}", r.layer, r.max_dev);
        }
    }
}
