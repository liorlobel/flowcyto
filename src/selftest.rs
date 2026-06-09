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

    // ── Gating layer: population counts + median MFI of gated populations, applied
    //    by flowcyto's own gating engine to the compensated reference and compared to
    //    flowCore's rectangleGate/polygonGate golden. Gates are in compensated/linear
    //    space (Linear transform) so both operate on identical coordinates; G2 is a
    //    child of G1 (within-parent count). Keyed (channel, key) → value.
    let gates = {
        use crate::gating::{Gate, GateShape};
        use crate::transform::AxisTransform::Linear as Lin;
        vec![
            Gate { id: 1, name: "Cells".into(), parent: None,
                x_channel: "FSC-A".into(), y_channel: "SSC-A".into(), x_transform: Lin, y_transform: Lin,
                shape: GateShape::Rect { x_min: 20000.0, x_max: 1e6, y_min: -1e6, y_max: 1e6 }, quad_group: None },
            Gate { id: 2, name: "FITCpos".into(), parent: Some(1),
                x_channel: "FITC-A".into(), y_channel: "FITC-A".into(), x_transform: Lin, y_transform: Lin,
                shape: GateShape::Range { x_min: 5000.0, x_max: 1e6 }, quad_group: None },
            Gate { id: 3, name: "PEpos".into(), parent: None,
                x_channel: "PE-A".into(), y_channel: "SSC-A".into(), x_transform: Lin, y_transform: Lin,
                shape: GateShape::Polygon { vertices: vec![[5000.0, -1e6], [50000.0, -1e6], [50000.0, 1e6], [5000.0, 1e6]] }, quad_group: None },
        ]
    };
    let own = crate::gating::compute_own_masks(&gates, &comp, &fcs.parameters, fcs.n_events);
    let by_id: std::collections::HashMap<u32, &crate::gating::Gate> = gates.iter().map(|g| (g.id, g)).collect();
    let eff = |id: u32| crate::gating::effective_mask(id, &by_id, &own, fcs.n_events);
    let count = |m: &[bool]| m.iter().filter(|&&b| b).count() as f64;
    let median_in = |m: &[bool], ch: &str| -> f64 {
        let ci = fcs.parameters.iter().position(|p| p.name == ch).unwrap();
        let mut v: Vec<f64> = (0..fcs.n_events).filter(|&e| m[e]).map(|e| comp[e * np + ci]).collect();
        v.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let n = v.len();
        if n == 0 { 0.0 } else if n % 2 == 1 { v[n / 2] } else { 0.5 * (v[n / 2 - 1] + v[n / 2]) }
    };
    let (m1, m2, m3) = (eff(1), eff(2), eff(3));
    let mut gate_results: std::collections::HashMap<(String, String), f64> = std::collections::HashMap::new();
    gate_results.insert((String::new(), "Cells_count".into()), count(&m1));
    gate_results.insert((String::new(), "FITCpos_count".into()), count(&m2));
    gate_results.insert(("FITC-A".into(), "FITCpos_median".into()), median_in(&m2, "FITC-A"));
    gate_results.insert(("FSC-A".into(), "FITCpos_median".into()), median_in(&m2, "FSC-A"));
    gate_results.insert((String::new(), "PEpos_count".into()), count(&m3));
    gate_results.insert(("PE-A".into(), "PEpos_median".into()), median_in(&m3, "PE-A"));

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
            "gate" => *gate_results
                .get(&(channel.to_string(), key.to_string()))
                .ok_or_else(|| format!("no gate result for {channel}/{key}"))?,
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
        "gate" => "Gating",
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
