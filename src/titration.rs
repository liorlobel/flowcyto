//! Antibody titration analysis — per-sample **stain index** across a dilution series.
//!
//! Stain index (Maecker & Trotter; the value FlowJo/FlowLogic report) =
//!     (MFI_pos − MFI_neg) / (2 × rSD_neg)
//! where `rSD` is the **robust** SD of the negative population on the antibody's channel
//! (robust to the skewed tails of compensated negatives, where parametric SD overstates
//! the spread and understates the index). This module is a pure view over the batch
//! population-stats tables — no re-streaming of events.

use crate::popstats::PopulationStatsTable;

/// Stain index `(MFI_pos − MFI_neg) / (2 × rSD_neg)`. NaN when any input is non-finite or
/// the negative spread is ≈0 (separation is then undefined).
pub fn stain_index(mfi_pos: f64, mfi_neg: f64, rsd_neg: f64) -> f64 {
    if !mfi_pos.is_finite() || !mfi_neg.is_finite() || !rsd_neg.is_finite() || rsd_neg.abs() < 1e-9 {
        return f64::NAN;
    }
    (mfi_pos - mfi_neg) / (2.0 * rsd_neg)
}

/// One titration data point: a sample, its concentration (parsed from the condition tag),
/// and the stain-index inputs/result on the chosen channel.
#[derive(Clone)]
pub struct TitrationRow {
    pub group: String,
    pub sample: String,
    pub conc: Option<f64>, // from the condition/group tag (None = non-numeric tag)
    pub pct_pos: f64,      // positive population as % of its parent
    pub mfi_pos: f64,
    pub mfi_neg: f64,
    pub rsd_neg: f64,
    pub si: f64,
}

/// Build titration rows from batch population-stats tables.
/// * `channel`  — the titrated fluorochrome's parameter name.
/// * `pos_name` — the positive population's gate name.
/// * `neg_name` — the negative population's gate name (`"All events"` = the root).
///
/// Samples whose table lacks the channel are skipped; a missing pos/neg population yields
/// NaN inputs (and therefore a NaN stain index) rather than dropping the row.
pub fn titration_rows(
    tables: &[(String, String, PopulationStatsTable)],
    channel: &str,
    pos_name: &str,
    neg_name: &str,
) -> Vec<TitrationRow> {
    let mut out = Vec::with_capacity(tables.len());
    for (group, sample, table) in tables {
        let ci = match table.channels.iter().position(|c| c.eq_ignore_ascii_case(channel)) {
            Some(i) => i,
            None => continue,
        };
        let pos = table.rows.iter().find(|r| r.name == pos_name);
        let neg = table.rows.iter().find(|r| r.name == neg_name);
        let (mfi_pos, pct_pos) = pos.map_or((f64::NAN, f64::NAN), |p| (p.medians[ci], p.pct_parent));
        let (mfi_neg, rsd_neg) = neg.map_or((f64::NAN, f64::NAN), |n| (n.medians[ci], n.rsds[ci]));
        out.push(TitrationRow {
            group: group.clone(),
            sample: sample.clone(),
            conc: group.trim().parse::<f64>().ok().filter(|c| c.is_finite()),
            pct_pos,
            mfi_pos,
            mfi_neg,
            rsd_neg,
            si: stain_index(mfi_pos, mfi_neg, rsd_neg),
        });
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::popstats::{PopulationStat, PopulationStatsTable};

    fn stat(name: &str, pct_parent: f64, median: f64, rsd: f64) -> PopulationStat {
        PopulationStat {
            name: name.into(),
            parent_name: "All events".into(),
            depth: 1,
            count: 100,
            pct_parent,
            pct_total: pct_parent,
            medians: vec![median],
            means: vec![median],
            cvs: vec![0.0],
            rsds: vec![rsd],
        }
    }

    #[test]
    fn stain_index_known() {
        // (1000 − 100) / (2 × 50) = 9.0
        assert_eq!(stain_index(1000.0, 100.0, 50.0), 9.0);
    }

    #[test]
    fn stain_index_degenerate_is_nan() {
        assert!(stain_index(1000.0, 100.0, 0.0).is_nan()); // zero negative spread
        assert!(stain_index(f64::NAN, 100.0, 50.0).is_nan()); // missing positive
        assert!(stain_index(1000.0, 100.0, f64::NAN).is_nan()); // missing rSD
    }

    #[test]
    fn titration_rows_compute_si_and_parse_conc() {
        let table = PopulationStatsTable {
            channels: vec!["FITC-A".into()],
            rows: vec![
                stat("All events", 100.0, 100.0, 50.0), // negative proxy
                stat("CD3+", 80.0, 1000.0, 30.0),        // positive
            ],
        };
        let tables = vec![("10".to_string(), "tube1".to_string(), table)];
        let rows = titration_rows(&tables, "FITC-A", "CD3+", "All events");
        assert_eq!(rows.len(), 1);
        let r = &rows[0];
        assert_eq!(r.conc, Some(10.0));
        assert!((r.mfi_pos - 1000.0).abs() < 1e-9);
        assert!((r.mfi_neg - 100.0).abs() < 1e-9);
        assert!((r.si - 9.0).abs() < 1e-9); // (1000−100)/(2×50)
        assert!((r.pct_pos - 80.0).abs() < 1e-9);
    }

    #[test]
    fn titration_rows_nonnumeric_tag_is_none() {
        let table = PopulationStatsTable {
            channels: vec!["FITC-A".into()],
            rows: vec![stat("All events", 100.0, 100.0, 50.0), stat("P", 50.0, 500.0, 25.0)],
        };
        let tables = vec![("unstained".to_string(), "t".to_string(), table)];
        let rows = titration_rows(&tables, "FITC-A", "P", "All events");
        assert_eq!(rows[0].conc, None);
    }
}
