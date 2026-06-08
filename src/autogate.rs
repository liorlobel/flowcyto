//! Auto-gating *suggestions* (the user reviews and adjusts them):
//!  - `valley_threshold`: a 1-D split at the density valley between two peaks
//!    (openCyto `mindensity` style — the field-standard for flow, and unlike Otsu it
//!    handles rare-but-separated populations: a small positive peak still has a deep
//!    valley below it). The valley *depth* is an honest bimodality measure.
//!  - `singlet_polygon`: a diagonal band on FSC area×height, centered on the density
//!    **mode** of the area/height ratio (robust to the doublets it should exclude,
//!    which sit above the diagonal and would bias a mean/median upward).
//!
//! Pure + unit-tested on synthetic known-answer data (flowCore has no equivalent to
//! validate against — these are heuristics, deliberately labeled as suggestions).

/// Histogram `vals` into `n_bins` over [min,max]; returns (bins, lo, hi) or None.
fn histogram(vals: &[f64], n_bins: usize) -> Option<(Vec<f64>, f64, f64)> {
    let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
    for &v in vals {
        if v.is_finite() {
            lo = lo.min(v);
            hi = hi.max(v);
        }
    }
    if !lo.is_finite() || hi <= lo || n_bins < 2 {
        return None;
    }
    let span = hi - lo;
    let mut h = vec![0.0f64; n_bins];
    for &v in vals {
        if v.is_finite() {
            let mut b = ((v - lo) / span * n_bins as f64) as usize;
            if b >= n_bins {
                b = n_bins - 1;
            }
            h[b] += 1.0;
        }
    }
    Some((h, lo, hi))
}

/// Moving-average smoothing, `passes` times with half-window `w`.
fn smooth(h: &[f64], w: usize, passes: usize) -> Vec<f64> {
    let mut cur = h.to_vec();
    for _ in 0..passes {
        let prev = cur.clone();
        for i in 0..prev.len() {
            let lo = i.saturating_sub(w);
            let hi = (i + w + 1).min(prev.len());
            cur[i] = prev[lo..hi].iter().sum::<f64>() / (hi - lo) as f64;
        }
    }
    cur
}

/// Local maxima at least `min_h` tall (with a strict rise on one side, so noise
/// plateaus don't all register).
fn find_peaks(s: &[f64], min_h: f64) -> Vec<usize> {
    let mut peaks = Vec::new();
    for i in 1..s.len().saturating_sub(1) {
        if s[i] >= s[i - 1] && s[i] >= s[i + 1] && s[i] >= min_h && (s[i] > s[i - 1] || s[i] > s[i + 1]) {
            peaks.push(i);
        }
    }
    peaks
}

/// 1-D split at the density valley between the two tallest peaks. Returns
/// `(threshold, depth)` where depth ∈ [0,1] is how deep the valley is relative to the
/// shorter peak (≈1 = cleanly bimodal, ≈0 = shallow). `None` when the channel isn't
/// bimodal (fewer than two peaks) — caller should then ask the user to gate manually
/// rather than place a meaningless cut.
pub fn valley_threshold(values: &[f64], n_bins: usize) -> Option<(f64, f64)> {
    let finite = values.iter().filter(|v| v.is_finite()).count();
    if finite < 50 {
        return None;
    }
    let (hist, lo, hi) = histogram(values, n_bins)?;
    let w = (n_bins / 48).max(1);
    let sm = smooth(&hist, w, 3);
    let max_h = sm.iter().cloned().fold(0.0_f64, f64::max);
    if max_h <= 0.0 {
        return None;
    }
    // A peak must be ≥2% of the tallest — low enough to catch a rare positive
    // population, high enough to ignore noise.
    let peaks = find_peaks(&sm, max_h * 0.02);
    if peaks.len() < 2 {
        return None;
    }
    let mut pk = peaks;
    pk.sort_by(|&a, &b| sm[b].partial_cmp(&sm[a]).unwrap_or(std::cmp::Ordering::Equal));
    let (mut p0, mut p1) = (pk[0], pk[1]);
    if p0 > p1 {
        std::mem::swap(&mut p0, &mut p1);
    }
    let (mut vbin, mut vval) = (p0, sm[p0]);
    for b in p0..=p1 {
        if sm[b] < vval {
            vval = sm[b];
            vbin = b;
        }
    }
    let shorter = sm[p0].min(sm[p1]);
    let depth = if shorter > 0.0 { ((shorter - vval) / shorter).clamp(0.0, 1.0) } else { 0.0 };
    let span = hi - lo;
    let threshold = lo + (vbin as f64 + 0.5) / n_bins as f64 * span;
    Some((threshold, depth))
}

/// Densest value (mode) of `vals` via a histogram — robust to a minority tail.
fn mode(vals: &[f64], n_bins: usize) -> Option<f64> {
    let (hist, lo, hi) = histogram(vals, n_bins)?;
    let best = hist.iter().enumerate()
        .max_by(|a, b| a.1.partial_cmp(b.1).unwrap_or(std::cmp::Ordering::Equal))
        .map(|(i, _)| i)?;
    Some(lo + (best as f64 + 0.5) / n_bins as f64 * (hi - lo))
}

/// Suggested singlet gate on (area, height): a diagonal band `area ≈ slope·height`,
/// where `slope` is the **mode** of the area/height ratio over a robust height core
/// (so the doublets — which raise the ratio — neither bias the center nor leak in).
/// Returns parallelogram vertices in `[area, height]` coords (X=area, Y=height), or
/// `None` if degenerate. `tol` is the relative half-width of the band (e.g. 0.15).
pub fn singlet_polygon(area: &[f64], height: &[f64], tol: f64) -> Option<Vec<[f64; 2]>> {
    let n = area.len().min(height.len());
    let mut hs: Vec<f64> = (0..n).map(|i| height[i]).filter(|h| h.is_finite() && *h > 0.0).collect();
    if hs.len() < 50 {
        return None;
    }
    hs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let pct = |p: f64| hs[((p * (hs.len() - 1) as f64) as usize).min(hs.len() - 1)];
    let (h_lo, h_hi) = (pct(0.02), pct(0.98)); // trim debris / saturating tails
    if h_hi <= h_lo {
        return None;
    }
    let ratios: Vec<f64> = (0..n)
        .filter_map(|i| {
            let (a, h) = (area[i], height[i]);
            if a.is_finite() && h >= h_lo && h <= h_hi {
                Some(a / h)
            } else {
                None
            }
        })
        .collect();
    if ratios.len() < 50 {
        return None;
    }
    let slope = mode(&ratios, 64)?;
    if !(slope.is_finite() && slope > 0.0) {
        return None;
    }
    let edge = |h: f64, sign: f64| [slope * h * (1.0 + sign * tol), h];
    Some(vec![
        edge(h_lo, -1.0),
        edge(h_hi, -1.0),
        edge(h_hi, 1.0),
        edge(h_lo, 1.0),
    ])
}

#[cfg(test)]
mod tests {
    use super::*;

    // Deterministic pseudo-random in [0,1) (no rng dependency in tests).
    fn lcg(seed: &mut u64) -> f64 {
        *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        (*seed >> 33) as f64 / (1u64 << 31) as f64
    }
    fn gauss(seed: &mut u64, mu: f64, sd: f64) -> f64 {
        // Box-Muller
        let (u1, u2) = (lcg(seed).max(1e-12), lcg(seed));
        mu + sd * (-2.0 * u1.ln()).sqrt() * (std::f64::consts::TAU * u2).cos()
    }

    #[test]
    fn valley_splits_balanced_bimodal() {
        let mut s = 1;
        let mut v = Vec::new();
        for _ in 0..3000 { v.push(gauss(&mut s, 1.0, 0.4)); }
        for _ in 0..3000 { v.push(gauss(&mut s, 5.0, 0.4)); }
        let (t, depth) = valley_threshold(&v, 128).expect("bimodal → Some");
        assert!((2.0..4.0).contains(&t), "threshold in the valley ~3, got {t}");
        assert!(depth > 0.7, "clean split → deep valley, got {depth}");
    }

    #[test]
    fn valley_finds_rare_positive_population() {
        // The case Otsu/η fails: 95% negative, 5% positive but cleanly separated.
        let mut s = 42;
        let mut v = Vec::new();
        for _ in 0..9500 { v.push(gauss(&mut s, 1.0, 0.4)); }
        for _ in 0..500  { v.push(gauss(&mut s, 5.0, 0.4)); }
        let (t, _) = valley_threshold(&v, 128).expect("rare-but-separated → Some");
        assert!((2.0..4.5).contains(&t), "valley between the modes, got {t}");
    }

    #[test]
    fn valley_none_on_unimodal() {
        let mut s = 7;
        let v: Vec<f64> = (0..4000).map(|_| gauss(&mut s, 3.0, 0.8)).collect();
        assert!(valley_threshold(&v, 128).is_none(), "single population → no split");
    }

    #[test]
    fn singlet_band_excludes_doublets() {
        // Singlets on the diagonal (ratio≈1) + a 15% doublet cluster (ratio≈2, above).
        let mut s = 3;
        let mut area = Vec::new();
        let mut height = Vec::new();
        let mut is_doublet = Vec::new();
        for _ in 0..8500 {
            let h = 5000.0 + lcg(&mut s) * 55000.0;
            area.push(h + gauss(&mut s, 0.0, 1500.0)); // ratio ≈ 1
            height.push(h);
            is_doublet.push(false);
        }
        for _ in 0..1500 {
            let h = 5000.0 + lcg(&mut s) * 55000.0;
            area.push(2.0 * h + gauss(&mut s, 0.0, 2000.0)); // ratio ≈ 2
            height.push(h);
            is_doublet.push(true);
        }
        let poly = singlet_polygon(&area, &height, 0.15).expect("Some");
        let inside = |x: f64, y: f64| point_in_poly(x, y, &poly);
        let mut s_in = 0usize;
        let mut d_in = 0usize;
        let (mut s_tot, mut d_tot) = (0usize, 0usize);
        for i in 0..area.len() {
            let hit = inside(area[i], height[i]);
            if is_doublet[i] { d_tot += 1; if hit { d_in += 1; } }
            else { s_tot += 1; if hit { s_in += 1; } }
        }
        let s_frac = s_in as f64 / s_tot as f64;
        let d_frac = d_in as f64 / d_tot as f64;
        assert!(s_frac > 0.7, "keep most singlets, got {:.2}", s_frac);
        assert!(d_frac < 0.1, "exclude doublets, got {:.2}", d_frac);
    }

    // Local even-odd point-in-polygon (mirrors gating::point_in_polygon for the test).
    fn point_in_poly(px: f64, py: f64, v: &[[f64; 2]]) -> bool {
        let n = v.len();
        let mut inside = false;
        let mut j = n - 1;
        for i in 0..n {
            let (xi, yi) = (v[i][0], v[i][1]);
            let (xj, yj) = (v[j][0], v[j][1]);
            if ((yi > py) != (yj > py)) && (px < (xj - xi) * (py - yi) / (yj - yi) + xi) {
                inside = !inside;
            }
            j = i;
        }
        inside
    }
}
