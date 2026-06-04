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
