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
use crate::gating::{compute_own_masks, effective_mask, gate_tree_order, Gate};

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

    // Own membership per gate (geometric + Boolean), then effective (AND with ancestors).
    let own = compute_own_masks(gates, events, params, n_events);
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gating::GateShape;
    use crate::test_util::param;
    use crate::transform::AxisTransform;

    fn rect_gate(id: u32, parent: Option<u32>, name: &str, xmin: f64, xmax: f64) -> Gate {
        Gate {
            id,
            name: name.to_string(),
            parent,
            x_channel: "X".to_string(),
            y_channel: "Y".to_string(),
            x_transform: AxisTransform::Linear,
            y_transform: AxisTransform::Linear,
            shape: GateShape::Rect { x_min: xmin, x_max: xmax, y_min: -1.0, y_max: 1.0 },
            quad_group: None,
        }
    }

    #[test]
    fn root_row_is_all_events() {
        let params = vec![param(1, "X"), param(2, "Y")];
        let events = vec![1.0, 0.0, 2.0, 0.0, 3.0, 0.0];
        let table = population_stats(&events, &params, 3, &[], &[0]);
        assert_eq!(table.rows.len(), 1);
        assert_eq!(table.rows[0].name, "All events");
        assert_eq!(table.rows[0].count, 3);
        assert!((table.rows[0].pct_total - 100.0).abs() < 1e-9);
        // Median of X over {1,2,3} = 2.
        assert_eq!(table.rows[0].medians[0], 2.0);
    }

    #[test]
    fn counts_and_percentages_match_gate_tree() {
        let params = vec![param(1, "X"), param(2, "Y")];
        // x = 1,2,3,4,10 (y=0).
        let events = vec![1.0, 0.0, 2.0, 0.0, 3.0, 0.0, 4.0, 0.0, 10.0, 0.0];
        let parent = rect_gate(1, None, "P", 0.0, 5.0);     // 4 events
        let child = rect_gate(2, Some(1), "C", 0.0, 2.5);    // 2 events
        let table = population_stats(&events, &params, 5, &[parent, child], &[0]);

        // Rows: All events, P, C.
        assert_eq!(table.rows.len(), 3);
        let p = &table.rows[1];
        let c = &table.rows[2];
        assert_eq!(p.count, 4);
        assert!((p.pct_total - 80.0).abs() < 1e-9);
        assert!((p.pct_parent - 80.0).abs() < 1e-9); // parent is root (5 events)
        assert_eq!(c.count, 2);
        assert_eq!(c.parent_name, "P");
        assert!((c.pct_parent - 50.0).abs() < 1e-9); // 2 of P's 4
        assert!((c.pct_total - 40.0).abs() < 1e-9);
    }

    #[test]
    fn missing_channel_gate_reads_zero_not_parent() {
        // A gate on a channel this sample lacks must report 0, not collapse to parent.
        let params = vec![param(1, "X"), param(2, "Y")];
        let events = vec![1.0, 0.0, 2.0, 0.0];
        let mut g = rect_gate(1, None, "bad", 0.0, 10.0);
        g.x_channel = "ABSENT".to_string();
        let table = population_stats(&events, &params, 2, &[g], &[0]);
        assert_eq!(table.rows[1].count, 0);
    }

    #[test]
    fn med_mean_cv_known_values() {
        let (med, mean, cv) = med_mean_cv(&mut [2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0]);
        assert_eq!(med, 4.5);
        assert!((mean - 5.0).abs() < 1e-12);
        // Sample SD (n-1): sum of squared deviations = 32, /(8-1) = 32/7.
        let cv_expected = 100.0 * (32.0f64 / 7.0).sqrt() / 5.0;
        assert!((cv - cv_expected).abs() < 1e-9);
    }

    #[test]
    fn med_mean_cv_empty_is_nan() {
        let (med, mean, cv) = med_mean_cv(&mut []);
        assert!(med.is_nan() && mean.is_nan() && cv.is_nan());
    }

    #[test]
    fn long_csv_escapes_commas_in_names() {
        let table = PopulationStatsTable {
            channels: vec!["FITC-A".to_string()],
            rows: vec![PopulationStat {
                name: "CD11c+, MHCII+".to_string(), // contains a comma
                parent_name: "Live".to_string(),
                depth: 1,
                count: 10,
                pct_parent: 50.0,
                pct_total: 25.0,
                medians: vec![1234.5],
                means: vec![1200.0],
                cvs: vec![30.0],
            }],
        };
        let mut out = String::new();
        append_long_csv(&mut out, "sampleA", &table);
        assert!(out.contains("\"CD11c+, MHCII+\""), "comma'd name must be quoted: {out}");
    }

    #[test]
    fn long_csv_grouped_prepends_group() {
        let table = PopulationStatsTable {
            channels: vec!["FITC-A".to_string()],
            rows: vec![PopulationStat {
                name: "P".to_string(),
                parent_name: "All events".to_string(),
                depth: 1,
                count: 5,
                pct_parent: 100.0,
                pct_total: 100.0,
                medians: vec![1.0],
                means: vec![1.0],
                cvs: vec![0.0],
            }],
        };
        let mut out = String::new();
        append_long_csv_grouped(&mut out, "high_SAA", "tube1", &table);
        assert!(out.starts_with("high_SAA,tube1,P,"), "got: {out}");
    }

    #[test]
    fn header_constants_column_counts() {
        assert_eq!(LONG_CSV_HEADER.split(',').count(), 11);
        assert_eq!(LONG_CSV_HEADER_GROUPED.split(',').count(), 12);
    }
}
