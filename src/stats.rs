/// Per-channel summary statistics.
#[derive(Debug)]
#[allow(dead_code)]
pub struct Stats {
    pub name: String,
    pub n: usize,
    pub mean: f64,
    pub geo_mean: f64,  // geometric mean (positive values only)
    pub median: f64,
    pub std: f64,
    pub cv: f64,        // coefficient of variation (%)
    pub min: f64,
    pub max: f64,
    pub p05: f64,
    pub p25: f64,
    pub p75: f64,
    pub p95: f64,
}

impl Stats {
    pub fn compute(name: &str, values: &[f64]) -> Self {
        let n = values.len();
        if n == 0 {
            return Self::empty(name);
        }

        let mean = values.iter().sum::<f64>() / n as f64;

        // Geometric mean: exp(mean of ln(x)) for x > 0
        let geo_mean = {
            let pos: Vec<f64> = values.iter().filter(|&&v| v > 0.0).copied().collect();
            if pos.is_empty() {
                f64::NAN
            } else {
                let ln_mean = pos.iter().map(|v| v.ln()).sum::<f64>() / pos.len() as f64;
                ln_mean.exp()
            }
        };

        let variance = if n > 1 {
            values.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / (n - 1) as f64
        } else {
            0.0
        };
        let std = variance.sqrt();
        let cv = if mean.abs() > 1e-12 { 100.0 * std / mean.abs() } else { f64::NAN };

        let mut sorted = values.to_vec();
        sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

        let min = sorted[0];
        let max = sorted[n - 1];
        let median = percentile_sorted(&sorted, 50.0);
        let p05   = percentile_sorted(&sorted, 5.0);
        let p25   = percentile_sorted(&sorted, 25.0);
        let p75   = percentile_sorted(&sorted, 75.0);
        let p95   = percentile_sorted(&sorted, 95.0);

        Stats { name: name.to_string(), n, mean, geo_mean, median, std, cv, min, max, p05, p25, p75, p95 }
    }

    fn empty(name: &str) -> Self {
        Stats {
            name: name.to_string(), n: 0,
            mean: f64::NAN, geo_mean: f64::NAN, median: f64::NAN,
            std: f64::NAN, cv: f64::NAN, min: f64::NAN, max: f64::NAN,
            p05: f64::NAN, p25: f64::NAN, p75: f64::NAN, p95: f64::NAN,
        }
    }
}

/// Linear interpolation percentile on a pre-sorted slice.
fn percentile_sorted(sorted: &[f64], p: f64) -> f64 {
    let n = sorted.len();
    if n == 1 {
        return sorted[0];
    }
    let rank = p / 100.0 * (n - 1) as f64;
    let lo = rank.floor() as usize;
    let hi = (lo + 1).min(n - 1);
    let frac = rank - lo as f64;
    sorted[lo] + frac * (sorted[hi] - sorted[lo])
}

/// Print a stats table to stdout.
pub fn print_stats_table(stats: &[Stats]) {
    println!(
        "{:<20} {:>8} {:>12} {:>12} {:>12} {:>8} {:>8} {:>8} {:>8} {:>8}",
        "Channel", "N", "Mean", "GeoMean", "Median", "SD", "CV%", "p05", "p75", "p95"
    );
    println!("{}", "─".repeat(112));
    for s in stats {
        println!(
            "{:<20} {:>8} {:>12.1} {:>12.1} {:>12.1} {:>8.1} {:>8.1} {:>8.1} {:>8.1} {:>8.1}",
            s.name, s.n, s.mean, s.geo_mean, s.median, s.std, s.cv, s.p05, s.p75, s.p95
        );
    }
}
