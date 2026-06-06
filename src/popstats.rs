//! Per-population statistics — the shared engine behind the Stats tab and the
//! multi-sample batch table.
//!
//! Pure function: takes compensated (linear) events + a gate tree, returns one
//! row per population (the ungated "All events" root plus every gate, in tree
//! order) with count, %-of-parent, %-of-total, and per-channel median (MFI),
//! mean, and CV. Statistics are computed on **compensated linear** data — the
//! space a flow cytometrist reads as MFI.

use std::collections::HashMap;

use crate::fcs::Parameter;
use crate::gating::{effective_mask, gate_membership, gate_tree_order, Gate};

#[derive(Clone)]
pub struct PopulationStat {
    pub name: String,
    pub parent_name: String,
    pub depth: usize,
    pub count: usize,
    pub pct_parent: f64,
    pub pct_total: f64,
    pub medians: Vec<f64>, // aligned to PopulationStatsTable::channels
    pub means: Vec<f64>,
    pub cvs: Vec<f64>,
}

pub struct PopulationStatsTable {
    pub channels: Vec<String>,
    pub rows: Vec<PopulationStat>,
}

/// Compute per-population statistics.
///
/// `events`        — compensated linear, row-major [n_events × n_params]
/// `stat_channels` — parameter indices to summarize (e.g. scatter + fluorescence)
pub fn population_stats(
    events: &[f64],
    params: &[Parameter],
    n_events: usize,
    gates: &[Gate],
    stat_channels: &[usize],
) -> PopulationStatsTable {
    let n_params = params.len();
    let channels: Vec<String> = stat_channels.iter().map(|&i| params[i].name.clone()).collect();

    // Own membership per gate, then effective (AND with ancestors).
    // On error (e.g. a channel this sample lacks), use an all-false mask so the
    // population reads 0 — NOT skipped, which would make `effective_mask` collapse
    // it to its parent's mask and report the parent count as a real population.
    let mut own: HashMap<u32, Vec<bool>> = HashMap::new();
    for g in gates {
        let m = gate_membership(g, events, params, n_events, n_params)
            .unwrap_or_else(|_| vec![false; n_events]);
        own.insert(g.id, m);
    }
    let by_id: HashMap<u32, &Gate> = gates.iter().map(|g| (g.id, g)).collect();

    let mut rows = Vec::with_capacity(gates.len() + 1);

    // Root: all events.
    let all_mask = vec![true; n_events];
    rows.push(stat_row("All events", "—", 0, &all_mask, n_events, n_events,
                       events, n_params, stat_channels));

    for (gid, depth) in gate_tree_order(gates) {
        let g = match by_id.get(&gid) { Some(g) => *g, None => continue };
        let eff = effective_mask(gid, &by_id, &own, n_events);
        let parent_count = match g.parent {
            Some(pid) => effective_mask(pid, &by_id, &own, n_events).iter().filter(|&&b| b).count(),
            None => n_events,
        };
        let parent_name = match g.parent {
            Some(pid) => by_id.get(&pid).map(|p| p.name.clone()).unwrap_or_else(|| "—".into()),
            None => "All events".into(),
        };
        rows.push(stat_row(&g.name, &parent_name, depth + 1, &eff, parent_count, n_events,
                           events, n_params, stat_channels));
    }

    PopulationStatsTable { channels, rows }
}

#[allow(clippy::too_many_arguments)]
fn stat_row(
    name: &str, parent_name: &str, depth: usize,
    mask: &[bool], parent_count: usize, total: usize,
    events: &[f64], n_params: usize, stat_channels: &[usize],
) -> PopulationStat {
    let count = mask.iter().filter(|&&b| b).count();
    let pct_parent = if parent_count > 0 { 100.0 * count as f64 / parent_count as f64 } else { 0.0 };
    let pct_total = if total > 0 { 100.0 * count as f64 / total as f64 } else { 0.0 };

    let mut medians = Vec::with_capacity(stat_channels.len());
    let mut means = Vec::with_capacity(stat_channels.len());
    let mut cvs = Vec::with_capacity(stat_channels.len());
    for &ci in stat_channels {
        let mut vals: Vec<f64> = Vec::with_capacity(count);
        for (ev, &inside) in mask.iter().enumerate() {
            if inside {
                vals.push(events[ev * n_params + ci]);
            }
        }
        let (med, mean, cv) = med_mean_cv(&mut vals);
        medians.push(med);
        means.push(mean);
        cvs.push(cv);
    }

    PopulationStat {
        name: name.to_string(), parent_name: parent_name.to_string(), depth,
        count, pct_parent, pct_total, medians, means, cvs,
    }
}

/// Median (R-style: average of two middle for even N), mean, CV% — single pass-ish.
fn med_mean_cv(vals: &mut [f64]) -> (f64, f64, f64) {
    let n = vals.len();
    if n == 0 {
        return (f64::NAN, f64::NAN, f64::NAN);
    }
    let mean = vals.iter().sum::<f64>() / n as f64;
    let var = if n > 1 {
        vals.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / (n - 1) as f64
    } else {
        0.0
    };
    let sd = var.sqrt();
    let cv = if mean.abs() > 1e-12 { 100.0 * sd / mean.abs() } else { f64::NAN };
    vals.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let med = if n % 2 == 1 { vals[n / 2] } else { 0.5 * (vals[n / 2 - 1] + vals[n / 2]) };
    (med, mean, cv)
}

// ── Tidy CSV (long format, R-friendly) ──────────────────────────────────────

pub const LONG_CSV_HEADER: &str =
    "sample,population,parent,depth,count,pct_parent,pct_total,channel,median,mean,cv";

/// Tidy CSV header including a leading `group` column (batch export with tags).
pub const LONG_CSV_HEADER_GROUPED: &str =
    "group,sample,population,parent,depth,count,pct_parent,pct_total,channel,median,mean,cv";

/// Like `append_long_csv`, but prepends a `group` (condition) column to every row.
pub fn append_long_csv_grouped(out: &mut String, group: &str, sample: &str, table: &PopulationStatsTable) {
    use std::fmt::Write;
    let esc = |s: &str| if s.contains(',') || s.contains('"') {
        format!("\"{}\"", s.replace('"', "\"\""))
    } else {
        s.to_string()
    };
    for r in &table.rows {
        for (k, ch) in table.channels.iter().enumerate() {
            let _ = writeln!(
                out, "{},{},{},{},{},{},{:.4},{:.4},{},{:.4},{:.4},{:.4}",
                esc(group), esc(sample), esc(&r.name), esc(&r.parent_name), r.depth,
                r.count, r.pct_parent, r.pct_total, esc(ch),
                r.medians[k], r.means[k], r.cvs[k],
            );
        }
    }
}

/// Append data rows (no header) for one sample to a long/tidy CSV buffer:
/// one row per (population × channel).
pub fn append_long_csv(out: &mut String, sample: &str, table: &PopulationStatsTable) {
    use std::fmt::Write;
    let esc = |s: &str| if s.contains(',') || s.contains('"') {
        format!("\"{}\"", s.replace('"', "\"\""))
    } else {
        s.to_string()
    };
    for r in &table.rows {
        for (k, ch) in table.channels.iter().enumerate() {
            let _ = writeln!(
                out, "{},{},{},{},{},{:.4},{:.4},{},{:.4},{:.4},{:.4}",
                esc(sample), esc(&r.name), esc(&r.parent_name), r.depth,
                r.count, r.pct_parent, r.pct_total, esc(ch),
                r.medians[k], r.means[k], r.cvs[k],
            );
        }
    }
}
