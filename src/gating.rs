use std::collections::HashMap;

use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};

use crate::fcs::Parameter;
use crate::transform::AxisTransform;

/// A 2-D gate.
///
/// Geometry is stored in **display (transformed) coordinates** — the space the
/// user actually draws in. The transform context binds to the *channel*
/// (`x_channel` + `x_transform`), so a gate keeps its meaning no matter which
/// axis a channel is parked on. Membership re-applies these transforms to the
/// data before the geometric test.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Gate {
    pub id: u32,
    pub name: String,
    #[serde(default)]
    pub parent: Option<u32>,
    pub x_channel: String,
    pub y_channel: String,
    #[serde(default)]
    pub x_transform: AxisTransform,
    #[serde(default)]
    pub y_transform: AxisTransform,
    pub shape: GateShape,
    /// Quadrant linkage: the four rects of one quadrant share a group id, so their
    /// common center can be moved together. `None` for ordinary gates.
    #[serde(default)]
    pub quad_group: Option<u32>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "shape", rename_all = "snake_case")]
pub enum GateShape {
    /// Axis-aligned rectangle in display coordinates.
    Rect { x_min: f64, x_max: f64, y_min: f64, y_max: f64 },
    /// Ellipse in display coordinates (rotation in radians, 0 = axis-aligned).
    Ellipse { cx: f64, cy: f64, rx: f64, ry: f64, angle: f64 },
    /// Polygon in display coordinates (auto-closed).
    Polygon { vertices: Vec<[f64; 2]> },
    /// 1-D interval on the x channel only (drawn on a histogram). `y` is ignored;
    /// for such gates the Gate's y_channel/y_transform mirror x.
    Range { x_min: f64, x_max: f64 },
}

impl GateShape {
    /// Point-in-shape test, in display coordinates.
    pub fn contains(&self, x: f64, y: f64) -> bool {
        match self {
            GateShape::Rect { x_min, x_max, y_min, y_max } => {
                x >= *x_min && x <= *x_max && y >= *y_min && y <= *y_max
            }
            GateShape::Ellipse { cx, cy, rx, ry, angle } => {
                if *rx <= 0.0 || *ry <= 0.0 { return false; }
                let (s, c) = angle.sin_cos();
                let dx = x - cx;
                let dy = y - cy;
                // rotate point into ellipse-local frame
                let xl = dx * c + dy * s;
                let yl = -dx * s + dy * c;
                (xl / rx).powi(2) + (yl / ry).powi(2) <= 1.0
            }
            GateShape::Polygon { vertices } => point_in_polygon(x, y, vertices),
            GateShape::Range { x_min, x_max } => x >= *x_min && x <= *x_max,
        }
    }

    /// Outline as a closed list of display-space points (for rendering).
    pub fn outline(&self) -> Vec<[f64; 2]> {
        match self {
            GateShape::Rect { x_min, x_max, y_min, y_max } => vec![
                [*x_min, *y_min], [*x_max, *y_min],
                [*x_max, *y_max], [*x_min, *y_max], [*x_min, *y_min],
            ],
            GateShape::Ellipse { cx, cy, rx, ry, angle } => {
                let (s, c) = angle.sin_cos();
                let n = 64;
                (0..=n).map(|i| {
                    let th = std::f64::consts::TAU * i as f64 / n as f64;
                    let ex = rx * th.cos();
                    let ey = ry * th.sin();
                    [cx + ex * c - ey * s, cy + ex * s + ey * c]
                }).collect()
            }
            GateShape::Polygon { vertices } => {
                let mut v: Vec<[f64; 2]> = vertices.clone();
                if let Some(first) = vertices.first() {
                    v.push(*first);
                }
                v
            }
            // Range is rendered specially on the histogram; no 2-D outline.
            GateShape::Range { .. } => Vec::new(),
        }
    }

    /// A representative point for placing the gate label.
    pub fn label_anchor(&self) -> [f64; 2] {
        match self {
            GateShape::Rect { x_min, x_max, y_max, .. } => [(x_min + x_max) / 2.0, *y_max],
            GateShape::Ellipse { cx, cy, ry, .. } => [*cx, *cy + *ry],
            GateShape::Polygon { vertices } => {
                if vertices.is_empty() { return [0.0, 0.0]; }
                let n = vertices.len() as f64;
                let cx = vertices.iter().map(|p| p[0]).sum::<f64>() / n;
                let cy = vertices.iter().map(|p| p[1])
                    .max_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
                    .unwrap_or(0.0);
                [cx, cy]
            }
            GateShape::Range { x_min, x_max } => [(x_min + x_max) / 2.0, 0.0],
        }
    }
}

/// "Own" membership mask for a single gate (geometric test only, no parent).
pub fn gate_membership(
    gate: &Gate,
    events: &[f64],
    params: &[Parameter],
    n_events: usize,
    n_params: usize,
) -> Result<Vec<bool>> {
    let xi = param_idx(params, &gate.x_channel)?;
    let yi = param_idx(params, &gate.y_channel)?;
    let xt = gate.x_transform.compile();
    let yt = gate.y_transform.compile();

    let mut mask = vec![false; n_events];
    for ev in 0..n_events {
        let base = ev * n_params;
        let dx = xt.forward(events[base + xi]);
        let dy = yt.forward(events[base + yi]);
        mask[ev] = gate.shape.contains(dx, dy);
    }
    Ok(mask)
}

/// Depth-first ordering of gates as (id, depth) — parents before children.
pub fn gate_tree_order(gates: &[Gate]) -> Vec<(u32, usize)> {
    let mut out = Vec::new();
    fn visit(gates: &[Gate], parent: Option<u32>, depth: usize, out: &mut Vec<(u32, usize)>) {
        for g in gates.iter().filter(|g| g.parent == parent) {
            out.push((g.id, depth));
            visit(gates, Some(g.id), depth + 1, out);
        }
    }
    visit(gates, None, 0, &mut out);
    // Append any orphans (parent id not present) so nothing is dropped.
    for g in gates {
        if !out.iter().any(|(id, _)| *id == g.id) {
            out.push((g.id, 0));
        }
    }
    out
}

/// Effective mask of a gate = own mask AND all ancestors' masks.
/// `own` maps gate id → own membership mask.
pub fn effective_mask(
    gate_id: u32,
    gates_by_id: &HashMap<u32, &Gate>,
    own: &HashMap<u32, Vec<bool>>,
    n_events: usize,
) -> Vec<bool> {
    let mut result = vec![true; n_events];
    let mut current = Some(gate_id);
    let mut guard = 0;
    while let Some(id) = current {
        guard += 1;
        if guard > 1000 { break; } // cycle guard
        if let Some(m) = own.get(&id) {
            for i in 0..n_events {
                result[i] &= m[i];
            }
        }
        current = gates_by_id.get(&id).and_then(|g| g.parent);
    }
    result
}

#[derive(Debug)]
pub struct GateResult {
    pub name: String,
    pub n_in: usize,
    pub n_parent: usize,
    pub n_total: usize,
}

impl GateResult {
    pub fn pct_total(&self) -> f64 {
        if self.n_total == 0 { 0.0 } else { 100.0 * self.n_in as f64 / self.n_total as f64 }
    }
    pub fn pct_parent(&self) -> f64 {
        if self.n_parent == 0 { 0.0 } else { 100.0 * self.n_in as f64 / self.n_parent as f64 }
    }
}

/// Evaluate all gates (respecting hierarchy) and return counts.
pub fn apply_gates(
    gates: &[Gate],
    events: &[f64],
    params: &[Parameter],
    n_events: usize,
) -> Result<Vec<GateResult>> {
    let n_params = params.len();

    // Own masks
    let mut own: HashMap<u32, Vec<bool>> = HashMap::new();
    for g in gates {
        own.insert(g.id, gate_membership(g, events, params, n_events, n_params)?);
    }
    let by_id: HashMap<u32, &Gate> = gates.iter().map(|g| (g.id, g)).collect();

    let mut results = Vec::with_capacity(gates.len());
    for g in gates {
        let eff = effective_mask(g.id, &by_id, &own, n_events);
        let n_in = eff.iter().filter(|&&b| b).count();
        let n_parent = match g.parent {
            Some(pid) => {
                let pm = effective_mask(pid, &by_id, &own, n_events);
                pm.iter().filter(|&&b| b).count()
            }
            None => n_events,
        };
        results.push(GateResult { name: g.name.clone(), n_in, n_parent, n_total: n_events });
    }
    Ok(results)
}

// ── Helpers ───────────────────────────────────────────────────────────────

fn param_idx(params: &[Parameter], name: &str) -> Result<usize> {
    params
        .iter()
        .position(|p| p.name.eq_ignore_ascii_case(name))
        .with_context(|| {
            format!(
                "gate channel '{}' not found; available: [{}]",
                name,
                params.iter().map(|p| p.name.as_str()).collect::<Vec<_>>().join(", ")
            )
        })
}

/// Even-odd ray-casting point-in-polygon.
fn point_in_polygon(px: f64, py: f64, verts: &[[f64; 2]]) -> bool {
    let n = verts.len();
    if n < 3 { return false; }
    let mut inside = false;
    let mut j = n - 1;
    for i in 0..n {
        let (xi, yi) = (verts[i][0], verts[i][1]);
        let (xj, yj) = (verts[j][0], verts[j][1]);
        if ((yi > py) != (yj > py))
            && (px < (xj - xi) * (py - yi) / (yj - yi) + xi)
        {
            inside = !inside;
        }
        j = i;
    }
    inside
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_util::param;

    fn gate(id: u32, name: &str, parent: Option<u32>, shape: GateShape) -> Gate {
        Gate {
            id,
            name: name.to_string(),
            parent,
            x_channel: "X".to_string(),
            y_channel: "Y".to_string(),
            x_transform: AxisTransform::Linear,
            y_transform: AxisTransform::Linear,
            shape,
            quad_group: None,
        }
    }

    // ── GateShape::contains ────────────────────────────────────────────────

    #[test]
    fn rect_contains() {
        let r = GateShape::Rect { x_min: 0.0, x_max: 10.0, y_min: 0.0, y_max: 10.0 };
        assert!(r.contains(5.0, 5.0));
        assert!(r.contains(0.0, 0.0)); // boundary inclusive
        assert!(r.contains(10.0, 10.0));
        assert!(!r.contains(-0.1, 5.0));
        assert!(!r.contains(5.0, 10.1));
    }

    #[test]
    fn range_ignores_y() {
        let r = GateShape::Range { x_min: 1.0, x_max: 3.0 };
        assert!(r.contains(2.0, -9999.0));
        assert!(!r.contains(0.5, 0.0));
    }

    #[test]
    fn ellipse_axis_aligned() {
        let e = GateShape::Ellipse { cx: 0.0, cy: 0.0, rx: 2.0, ry: 1.0, angle: 0.0 };
        assert!(e.contains(0.0, 0.0));
        assert!(e.contains(2.0, 0.0)); // on boundary along major axis
        assert!(e.contains(0.0, 1.0));
        assert!(!e.contains(0.0, 1.5));
        assert!(!e.contains(2.0, 1.0));
    }

    #[test]
    fn ellipse_rotated_90_degrees() {
        // Rotating a 2×1 ellipse by 90° swaps which axis is long.
        let e = GateShape::Ellipse {
            cx: 0.0, cy: 0.0, rx: 2.0, ry: 1.0,
            angle: std::f64::consts::FRAC_PI_2,
        };
        // Now the long axis is vertical: (0, 2) inside, (2, 0) outside.
        assert!(e.contains(0.0, 1.9));
        assert!(!e.contains(1.9, 0.0));
    }

    #[test]
    fn ellipse_zero_radius_is_empty() {
        let e = GateShape::Ellipse { cx: 0.0, cy: 0.0, rx: 0.0, ry: 1.0, angle: 0.0 };
        assert!(!e.contains(0.0, 0.0));
    }

    #[test]
    fn polygon_triangle_and_concave() {
        let tri = GateShape::Polygon {
            vertices: vec![[0.0, 0.0], [4.0, 0.0], [2.0, 4.0]],
        };
        assert!(tri.contains(2.0, 1.0));
        assert!(!tri.contains(0.0, 3.0));

        // Degenerate polygon (< 3 vertices) contains nothing.
        let line = GateShape::Polygon { vertices: vec![[0.0, 0.0], [1.0, 1.0]] };
        assert!(!line.contains(0.5, 0.5));
    }

    #[test]
    fn outline_rect_is_closed() {
        let r = GateShape::Rect { x_min: 0.0, x_max: 1.0, y_min: 0.0, y_max: 1.0 };
        let o = r.outline();
        assert_eq!(o.first(), o.last());
        assert_eq!(o.len(), 5);
    }

    #[test]
    fn label_anchor_rect_top_center() {
        let r = GateShape::Rect { x_min: 0.0, x_max: 10.0, y_min: 0.0, y_max: 4.0 };
        assert_eq!(r.label_anchor(), [5.0, 4.0]);
    }

    // ── tree ordering ──────────────────────────────────────────────────────

    #[test]
    fn tree_order_parents_before_children() {
        let g = vec![
            gate(1, "root", None, GateShape::Rect { x_min: 0.0, x_max: 1.0, y_min: 0.0, y_max: 1.0 }),
            gate(2, "child", Some(1), GateShape::Rect { x_min: 0.0, x_max: 1.0, y_min: 0.0, y_max: 1.0 }),
            gate(3, "grandchild", Some(2), GateShape::Rect { x_min: 0.0, x_max: 1.0, y_min: 0.0, y_max: 1.0 }),
        ];
        let order = gate_tree_order(&g);
        assert_eq!(order, vec![(1, 0), (2, 1), (3, 2)]);
    }

    #[test]
    fn tree_order_appends_orphans() {
        // Gate 2's parent (99) is absent → it's an orphan, must still appear.
        let g = vec![
            gate(1, "root", None, GateShape::Rect { x_min: 0.0, x_max: 1.0, y_min: 0.0, y_max: 1.0 }),
            gate(2, "orphan", Some(99), GateShape::Rect { x_min: 0.0, x_max: 1.0, y_min: 0.0, y_max: 1.0 }),
        ];
        let order = gate_tree_order(&g);
        assert!(order.iter().any(|&(id, _)| id == 2), "orphan must not be dropped");
        assert_eq!(order.len(), 2);
    }

    // ── masks & counts ─────────────────────────────────────────────────────

    fn xy(params: &[Parameter], events: &[f64], n: usize) -> (Vec<Parameter>, Vec<f64>, usize) {
        (params.to_vec(), events.to_vec(), n)
    }

    #[test]
    fn effective_mask_ands_with_ancestors() {
        let params = vec![param(1, "X"), param(2, "Y")];
        // Three events at x = 1, 5, 9 (y = 0).
        let events = vec![1.0, 0.0, 5.0, 0.0, 9.0, 0.0];
        let (params, events, n) = xy(&params, &events, 3);

        let parent = gate(1, "P", None, GateShape::Rect { x_min: 0.0, x_max: 6.0, y_min: -1.0, y_max: 1.0 });
        let child = gate(2, "C", Some(1), GateShape::Rect { x_min: 4.0, x_max: 100.0, y_min: -1.0, y_max: 1.0 });

        let mut own = HashMap::new();
        own.insert(1, gate_membership(&parent, &events, &params, n, 2).unwrap());
        own.insert(2, gate_membership(&child, &events, &params, n, 2).unwrap());
        let by_id: HashMap<u32, &Gate> =
            [(1u32, &parent), (2u32, &child)].into_iter().collect();

        // Parent admits x∈[0,6] → events 0,1. Child admits x≥4 → events 1,2.
        // Effective child = AND → only event 1.
        let eff = effective_mask(2, &by_id, &own, n);
        assert_eq!(eff, vec![false, true, false]);
    }

    #[test]
    fn effective_mask_cycle_guard_terminates() {
        // a→b→a cycle; guard must stop rather than loop forever.
        let a = gate(1, "a", Some(2), GateShape::Rect { x_min: 0.0, x_max: 1.0, y_min: 0.0, y_max: 1.0 });
        let b = gate(2, "b", Some(1), GateShape::Rect { x_min: 0.0, x_max: 1.0, y_min: 0.0, y_max: 1.0 });
        let by_id: HashMap<u32, &Gate> = [(1u32, &a), (2u32, &b)].into_iter().collect();
        let mut own = HashMap::new();
        own.insert(1, vec![true]);
        own.insert(2, vec![true]);
        let _ = effective_mask(1, &by_id, &own, 1); // must return, not hang
    }

    #[test]
    fn apply_gates_counts_and_percentages() {
        let params = vec![param(1, "X"), param(2, "Y")];
        // 4 events at x = 1, 2, 3, 8.
        let events = vec![1.0, 0.0, 2.0, 0.0, 3.0, 0.0, 8.0, 0.0];

        let parent = gate(1, "P", None, GateShape::Rect { x_min: 0.0, x_max: 5.0, y_min: -1.0, y_max: 1.0 });
        let child = gate(2, "C", Some(1), GateShape::Rect { x_min: 0.0, x_max: 2.5, y_min: -1.0, y_max: 1.0 });

        let res = apply_gates(&[parent, child], &events, &params, 4).unwrap();
        // Parent: x≤5 → 3 of 4. Child: x≤2.5 ∧ parent → 2.
        assert_eq!(res[0].n_in, 3);
        assert_eq!(res[0].n_parent, 4);
        assert_eq!(res[1].n_in, 2);
        assert_eq!(res[1].n_parent, 3);
        assert!((res[1].pct_parent() - 100.0 * 2.0 / 3.0).abs() < 1e-9);
        assert!((res[1].pct_total() - 50.0).abs() < 1e-9);
    }

    #[test]
    fn gate_membership_unknown_channel_errors() {
        let params = vec![param(1, "X"), param(2, "Y")];
        let events = vec![1.0, 2.0];
        let mut g = gate(1, "g", None, GateShape::Rect { x_min: 0.0, x_max: 1.0, y_min: 0.0, y_max: 1.0 });
        g.x_channel = "NOPE".to_string();
        assert!(gate_membership(&g, &events, &params, 1, 2).is_err());
    }

    #[test]
    fn gate_result_zero_denominator_guards() {
        let r = GateResult { name: "x".into(), n_in: 0, n_parent: 0, n_total: 0 };
        assert_eq!(r.pct_parent(), 0.0);
        assert_eq!(r.pct_total(), 0.0);
    }
}
