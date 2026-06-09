//! Gating-ML 2.0 export — serialize flowcyto's gate tree to the ISAC Gating-ML 2.0
//! XML standard so gates interoperate with flowCore/CytoML, FlowKit, and other tools.
//!
//! Approach (the key trick): flowcyto stores gate bounds in *display* (post-transform)
//! coordinates. Each axis-aligned gate (Rect, Range) is emitted as a Gating-ML
//! `RectangleGate` with its bounds converted back to **data space** via the axis
//! transform's inverse — exact for any monotonic transform, because a rectangle in
//! display space is a rectangle in data space. So we never declare a Gating-ML
//! transform: gates are written on the (compensated) raw parameters directly. A
//! `Polygon`/`Ellipse` is emitted as a `PolygonGate` with vertices converted to data
//! space — exact on linear axes; an N-gon approximation of the display-space edges on
//! nonlinear axes (flagged as a known limitation). `Boolean` → `BooleanGate`.
//!
//! Compensation: a dimension on a channel present in the file's `$SPILLOVER` (and with
//! compensation active) carries `compensation-ref="FCS"` (FlowKit/flowCore then apply
//! the file's embedded matrix); scatter/uncompensated channels carry `"uncompensated"`.
//! A custom in-app spillover override is NOT embedded (a known limitation) — export
//! reflects the file's embedded matrix.

use crate::gating::{BoolOp, Gate, GateShape};
use crate::transform::AxisTransform;
use std::collections::{HashMap, HashSet};

const NS_G: &str = "http://www.isac-net.org/std/Gating-ML/v2.0/gating";
const NS_DT: &str = "http://www.isac-net.org/std/Gating-ML/v2.0/datatypes";
const NS_TF: &str = "http://www.isac-net.org/std/Gating-ML/v2.0/transformations";

/// Bounds beyond this (in data space) are treated as open (omitted) — the FCS data
/// lives within ±2^18, so this never clips a real boundary, only "unbounded" sentinels.
const OPEN: f64 = 1e9;

fn esc(s: &str) -> String {
    s.replace('&', "&amp;").replace('<', "&lt;").replace('>', "&gt;").replace('"', "&quot;")
}

/// A unique, XML-safe id per gate (its name, sanitized, disambiguated by numeric id on
/// collision), so `parent_id` / boolean `gateReference` always resolve.
fn build_ids(gates: &[Gate]) -> HashMap<u32, String> {
    let mut used: HashSet<String> = HashSet::new();
    let mut out: HashMap<u32, String> = HashMap::new();
    for g in gates {
        let base: String = g.name.chars()
            .map(|c| if c.is_ascii_alphanumeric() || c == '_' || c == '-' { c } else { '_' })
            .collect();
        let base = if base.is_empty() { format!("gate_{}", g.id) } else { base };
        let mut id = base.clone();
        if !used.insert(id.clone()) {
            id = format!("{}_{}", base, g.id);
            used.insert(id.clone());
        }
        out.insert(g.id, id);
    }
    out
}

fn num(v: f64) -> String {
    // Compact, locale-independent; integers stay integer-looking.
    if v == v.trunc() && v.abs() < 1e15 { format!("{}", v as i64) } else { format!("{}", v) }
}

/// One `<gating:dimension>` for an axis: min/max in DATA space (transform inverse of the
/// display bounds), omitting an open end. `comp_ref` is "FCS" or "uncompensated".
fn dimension(channel: &str, comp_ref: &str, tf: &AxisTransform, d_min: f64, d_max: f64) -> String {
    let ct = tf.compile();
    let (lo, hi) = (ct.inverse(d_min), ct.inverse(d_max));
    let mut attrs = format!("gating:compensation-ref=\"{}\"", comp_ref);
    if lo.is_finite() && lo > -OPEN { attrs.push_str(&format!(" gating:min=\"{}\"", num(lo))); }
    if hi.is_finite() && hi < OPEN { attrs.push_str(&format!(" gating:max=\"{}\"", num(hi))); }
    format!(
        "    <gating:dimension {attrs}>\n      <data-type:fcs-dimension data-type:name=\"{}\"/>\n    </gating:dimension>\n",
        esc(channel)
    )
}

/// Vertices of a shape's outline in DATA space (transform inverse per axis), for polygon
/// export. Returns `(points, exact)` — `exact` is false when a nonlinear transform makes
/// the straight-in-display edges only approximate in data space.
fn polygon_points(g: &Gate) -> (Vec<[f64; 2]>, bool) {
    let xt = g.x_transform.compile();
    let yt = g.y_transform.compile();
    let exact = matches!(g.x_transform, AxisTransform::Linear) && matches!(g.y_transform, AxisTransform::Linear);
    let pts = g.shape.outline().iter().map(|p| [xt.inverse(p[0]), yt.inverse(p[1])]).collect();
    (pts, exact)
}

/// Serialize a gate list to a Gating-ML 2.0 document. `comp_channels` are the channel
/// names present in the file's `$SPILLOVER`; `compensate` is whether compensation is
/// active. Returns `(xml, warnings)`.
pub fn to_gating_ml(gates: &[Gate], comp_channels: &[String], compensate: bool) -> (String, Vec<String>) {
    let ids = build_ids(gates);
    let mut warnings = Vec::new();
    let comp_ref = |ch: &str| -> &'static str {
        if compensate && comp_channels.iter().any(|c| c == ch) { "FCS" } else { "uncompensated" }
    };

    let mut x = String::from("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    x.push_str(&format!(
        "<gating:Gating-ML xmlns:gating=\"{NS_G}\" xmlns:data-type=\"{NS_DT}\" xmlns:transforms=\"{NS_TF}\">\n"
    ));

    for g in gates {
        let id = &ids[&g.id];
        let parent = g.parent.and_then(|p| ids.get(&p))
            .map(|pid| format!(" gating:parent_id=\"{}\"", pid)).unwrap_or_default();

        match &g.shape {
            GateShape::Rect { x_min, x_max, y_min, y_max } => {
                x.push_str(&format!("  <gating:RectangleGate gating:id=\"{id}\"{parent}>\n"));
                x.push_str(&dimension(&g.x_channel, comp_ref(&g.x_channel), &g.x_transform, *x_min, *x_max));
                x.push_str(&dimension(&g.y_channel, comp_ref(&g.y_channel), &g.y_transform, *y_min, *y_max));
                x.push_str("  </gating:RectangleGate>\n");
            }
            GateShape::Range { x_min, x_max } => {
                x.push_str(&format!("  <gating:RectangleGate gating:id=\"{id}\"{parent}>\n"));
                x.push_str(&dimension(&g.x_channel, comp_ref(&g.x_channel), &g.x_transform, *x_min, *x_max));
                x.push_str("  </gating:RectangleGate>\n");
            }
            GateShape::Polygon { .. } | GateShape::Ellipse { .. } => {
                let (pts, exact) = polygon_points(g);
                if !exact {
                    warnings.push(format!(
                        "gate '{}': {} on a nonlinear axis is exported as a straight-edge polygon approximation",
                        g.name, if matches!(g.shape, GateShape::Ellipse{..}) { "ellipse" } else { "polygon" }));
                }
                x.push_str(&format!("  <gating:PolygonGate gating:id=\"{id}\"{parent}>\n"));
                x.push_str(&format!("    <gating:dimension gating:compensation-ref=\"{}\"><data-type:fcs-dimension data-type:name=\"{}\"/></gating:dimension>\n",
                    comp_ref(&g.x_channel), esc(&g.x_channel)));
                x.push_str(&format!("    <gating:dimension gating:compensation-ref=\"{}\"><data-type:fcs-dimension data-type:name=\"{}\"/></gating:dimension>\n",
                    comp_ref(&g.y_channel), esc(&g.y_channel)));
                for p in pts {
                    x.push_str(&format!(
                        "    <gating:vertex><gating:coordinate data-type:value=\"{}\"/><gating:coordinate data-type:value=\"{}\"/></gating:vertex>\n",
                        num(p[0]), num(p[1])));
                }
                x.push_str("  </gating:PolygonGate>\n");
            }
            GateShape::Boolean { op, refs } => {
                let tag = match op { BoolOp::And => "and", BoolOp::Or => "or", BoolOp::Not => "not" };
                let valid: Vec<&String> = refs.iter().filter_map(|r| ids.get(r)).collect();
                if valid.is_empty() {
                    warnings.push(format!("boolean gate '{}' has no resolvable references — skipped", g.name));
                    continue;
                }
                x.push_str(&format!("  <gating:BooleanGate gating:id=\"{id}\"{parent}>\n    <gating:{tag}>\n"));
                for r in valid {
                    x.push_str(&format!("      <gating:gateReference gating:ref=\"{}\"/>\n", r));
                }
                x.push_str(&format!("    </gating:{tag}>\n  </gating:BooleanGate>\n"));
            }
        }
    }
    x.push_str("</gating:Gating-ML>\n");
    (x, warnings)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn lin() -> AxisTransform { AxisTransform::Linear }

    #[test]
    fn exports_rect_range_polygon_with_hierarchy_and_comp_ref() {
        let gates = vec![
            Gate { id: 1, name: "Cells".into(), parent: None, x_channel: "FSC-A".into(), y_channel: "SSC-A".into(),
                x_transform: lin(), y_transform: lin(),
                shape: GateShape::Rect { x_min: 20000.0, x_max: 1e6, y_min: -1e6, y_max: 1e6 }, quad_group: None },
            Gate { id: 2, name: "FITCpos".into(), parent: Some(1), x_channel: "FITC-A".into(), y_channel: "FITC-A".into(),
                x_transform: lin(), y_transform: lin(),
                shape: GateShape::Range { x_min: 5000.0, x_max: 1e6 }, quad_group: None },
            Gate { id: 3, name: "PEpos".into(), parent: None, x_channel: "PE-A".into(), y_channel: "SSC-A".into(),
                x_transform: lin(), y_transform: lin(),
                shape: GateShape::Polygon { vertices: vec![[5000.0, -1e6], [50000.0, -1e6], [50000.0, 1e6], [5000.0, 1e6]] }, quad_group: None },
        ];
        let (xml, warns) = to_gating_ml(&gates, &["FITC-A".into(), "PE-A".into(), "PE-Cy7-A".into()], true);
        assert!(warns.is_empty(), "linear-axis gates should export without warnings: {warns:?}");
        assert!(xml.contains("Gating-ML"));
        assert!(xml.contains("RectangleGate gating:id=\"Cells\""));
        assert!(xml.contains("gating:parent_id=\"Cells\""), "hierarchy must be emitted");
        assert!(xml.contains("PolygonGate gating:id=\"PEpos\""));
        // FSC (scatter) uncompensated; FITC (fluor in spillover) -> FCS
        assert!(xml.contains("name=\"FSC-A\"") && xml.contains("compensation-ref=\"uncompensated\""));
        assert!(xml.contains("compensation-ref=\"FCS\""), "fluor dimension should reference FCS compensation");
        assert!(xml.contains("gating:min=\"20000\""), "linear bound exported verbatim");
        assert!(xml.contains("<gating:vertex>"));
    }

    #[test]
    fn asinh_axis_aligned_gate_exports_in_data_space() {
        // Asinh-display Range -> data-space RectangleGate (exact, axis-aligned). The
        // display bound 2.0 maps to data sinh(2)*150 ≈ 543.9.
        let g = vec![Gate { id: 1, name: "pos".into(), parent: None,
            x_channel: "FITC-A".into(), y_channel: "FITC-A".into(),
            x_transform: AxisTransform::Asinh { cofactor: 150.0 }, y_transform: AxisTransform::Asinh { cofactor: 150.0 },
            shape: GateShape::Range { x_min: 2.0, x_max: 100.0 }, quad_group: None }];
        let (xml, warns) = to_gating_ml(&g, &["FITC-A".into()], true);
        assert!(warns.is_empty(), "axis-aligned asinh gate is exact (no warning)");
        let expected = (2.0_f64).sinh() * 150.0; // ≈ 543.9
        assert!(xml.contains(&format!("gating:min=\"{}\"", num(expected))), "asinh bound must be inverse-mapped to data space; xml=\n{xml}");
        assert!(!xml.contains("gating:max"), "open upper bound (asinh 100 -> huge) should be omitted");
    }
}
