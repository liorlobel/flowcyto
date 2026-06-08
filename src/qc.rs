//! Acquisition-quality QC metrics: flow-rate (clog/bubble) stability and
//! margin/saturation events. Pure functions — the per-population %viable is computed
//! separately from the (flowCore-validated) gating engine.
//!
//! Flow-rate stability is the real clog detector: a clog or bubble interrupts the
//! stream, so *events per unit time* drops in that interval while the cells that do get
//! through still look normal. Binning a channel's median over time would miss that — so
//! we bin **event counts over the Time channel**.

use crate::fcs::Parameter;

#[derive(Clone, Debug)]
pub struct FlowRate {
    /// Events per equal-width Time interval (for the sparkline).
    pub bins: Vec<usize>,
    /// Worst bin's absolute deviation from the median bin count, as a percentage.
    pub max_dev_pct: f64,
    pub flagged: bool,
}

/// Flow-rate stability over the acquisition Time channel. Splits the run into equal-width
/// **time** intervals and counts events per interval; a clog shows up as an interval far
/// below the median rate (a bubble/resume as a spike above it).
///
/// Returns `None` when there's no usable Time signal (channel absent, constant, or all
/// non-finite) — the caller should then report flow-rate as N/A rather than fake it from
/// event order (equal-size order chunks carry no rate information).
///
/// `n_bins` is clamped so each bin holds ~100 events on average (avoids flagging Poisson
/// noise on sparse bins). `flag_dev_pct` is the deviation that trips the flag.
pub fn flow_rate(time: &[f64], n_bins: usize, flag_dev_pct: f64) -> Option<FlowRate> {
    let n = time.len();
    if n == 0 || n_bins == 0 {
        return None;
    }
    let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
    for &t in time {
        if t.is_finite() {
            lo = lo.min(t);
            hi = hi.max(t);
        }
    }
    if !lo.is_finite() || !hi.is_finite() || hi <= lo {
        return None; // constant / absent / non-finite Time → no flow-rate signal
    }
    let nb = n_bins.min((n / 100).max(1)); // keep ~>=100 events/bin
    let span = hi - lo;
    let mut bins = vec![0usize; nb];
    for &t in time {
        if !t.is_finite() {
            continue;
        }
        let mut b = ((t - lo) / span * nb as f64) as usize;
        if b >= nb {
            b = nb - 1; // the max value lands in the last bin
        }
        bins[b] += 1;
    }
    let mut sorted = bins.clone();
    sorted.sort_unstable();
    let med = sorted[nb / 2].max(1) as f64; // robust baseline rate
    let max_dev = bins
        .iter()
        .map(|&c| (c as f64 - med).abs() / med * 100.0)
        .fold(0.0_f64, f64::max);
    Some(FlowRate {
        bins,
        max_dev_pct: max_dev,
        flagged: max_dev > flag_dev_pct,
    })
}

/// Per-channel percentage of events pegged at the top of the channel's dynamic range
/// (`$PnR`) — off-scale / saturated acquisition. Computed on **raw** data (the
/// acquisition range, before compensation). Channels without a meaningful range are
/// skipped. Returns `(channel_name, percent)` for each kept channel.
pub fn margin_events(raw: &[f64], params: &[Parameter], n_events: usize) -> Vec<(String, f64)> {
    let n_params = params.len();
    if n_params == 0 || n_events == 0 || raw.len() < n_events * n_params {
        return Vec::new();
    }
    params
        .iter()
        .enumerate()
        .filter(|(_, p)| p.range >= 2.0) // need a real ADC ceiling
        .map(|(ci, p)| {
            let ceil = p.range * 0.999;
            let mut hits = 0usize;
            for ev in 0..n_events {
                if raw[ev * n_params + ci] >= ceil {
                    hits += 1;
                }
            }
            (p.name.clone(), 100.0 * hits as f64 / n_events as f64)
        })
        .collect()
}

/// Worst (largest) margin percentage across channels, for a single at-a-glance flag.
pub fn worst_margin(margins: &[(String, f64)]) -> Option<(String, f64)> {
    margins
        .iter()
        .cloned()
        .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_util::param;

    #[test]
    fn flow_rate_uniform_is_not_flagged() {
        // Even acquisition: time increases linearly with event index → flat rate.
        let time: Vec<f64> = (0..1000).map(|i| i as f64).collect();
        let fr = flow_rate(&time, 10, 40.0).unwrap();
        assert!(!fr.flagged, "uniform flow should not flag (max_dev {:.1}%)", fr.max_dev_pct);
        assert!(fr.max_dev_pct < 5.0);
    }

    #[test]
    fn flow_rate_clog_gap_is_flagged() {
        // A clog: no events acquired in the time window [400,600) → that bin is empty.
        let mut time: Vec<f64> = Vec::new();
        for i in 0..1000 {
            let t = i as f64;
            if (400.0..600.0).contains(&t) { continue; } // the clog gap
            time.push(t);
        }
        let fr = flow_rate(&time, 10, 40.0).unwrap();
        assert!(fr.flagged, "a clog gap must flag (max_dev {:.1}%)", fr.max_dev_pct);
        // The clog interval is a bin well below the busy bins.
        let mn = *fr.bins.iter().min().unwrap();
        let mx = *fr.bins.iter().max().unwrap();
        assert!((mn as f64) < 0.5 * mx as f64, "clog bin should be sparse: min={mn} max={mx}");
    }

    #[test]
    fn flow_rate_constant_or_absent_time_is_none() {
        assert!(flow_rate(&[5.0; 500], 10, 40.0).is_none(), "constant Time → no signal");
        assert!(flow_rate(&[], 10, 40.0).is_none(), "no events → None");
        assert!(flow_rate(&[f64::NAN; 100], 10, 40.0).is_none(), "all-NaN Time → None");
    }

    #[test]
    fn margin_events_counts_saturated() {
        // 2 channels, range 1000; channel 0 has 1 of 4 events at the ceiling.
        let params = vec![param(1, "FSC-A"), param(2, "PE-A")];
        // events row-major [ev*2 + ch]
        let raw = vec![
            10.0, 20.0,
            999.9, 30.0,   // FSC-A saturated (>= 999.0)
            40.0, 50.0,
            60.0, 70.0,
        ];
        // give both params range 1000
        let params: Vec<_> = params.into_iter().map(|mut p| { p.range = 1000.0; p }).collect();
        let m = margin_events(&raw, &params, 4);
        let fsc = m.iter().find(|(n, _)| n == "FSC-A").unwrap().1;
        assert!((fsc - 25.0).abs() < 1e-9, "1 of 4 saturated = 25% (got {fsc})");
        let pe = m.iter().find(|(n, _)| n == "PE-A").unwrap().1;
        assert_eq!(pe, 0.0);
    }
}
