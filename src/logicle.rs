//! Logicle (biexponential) data scale.
//!
//! Direct port of the Moore & Parks reference implementation
//! ("Update for the logicle data scale including operational code
//! implementations", Cytometry A, 2012) — the same C++ that flowCore wraps.
//! `R_zeroin` (Brent) is replaced by bisection: the root function is
//! monotonic on (0, b], so bisection is robust and exact to machine epsilon.
//!
//! `scale(value)`   : data  → display coordinate (≈ [0,1] for data in [0,T])
//! `inverse(scale)` : display coordinate → data
//!
//! `scale` is the expensive (root-found) direction; `inverse` is closed-form.

const LN_10: f64 = std::f64::consts::LN_10; // ln(10)
const TAYLOR_LENGTH: usize = 16;

#[derive(Clone, Debug)]
pub struct Logicle {
    // input parameters (kept for reference / serialization upstream)
    #[allow(dead_code)]
    pub t: f64,
    #[allow(dead_code)]
    pub w: f64,
    #[allow(dead_code)]
    pub m: f64,
    #[allow(dead_code)]
    pub a: f64,

    // biexponential coefficients
    a_c: f64,
    b_c: f64,
    c_c: f64,
    d_c: f64,
    f_c: f64,

    // breakpoints (normalized display coords)
    x1: f64,
    x_taylor: f64,

    taylor: [f64; TAYLOR_LENGTH],
}

impl Logicle {
    /// Build a logicle scale. Typical BD parameters: T=262144, W=0.5, M=4.5, A=0.
    pub fn new(t: f64, w: f64, m: f64, a: f64) -> Result<Logicle, String> {
        if t <= 0.0 { return Err("logicle: T must be positive".into()); }
        if w < 0.0 { return Err("logicle: W must be non-negative".into()); }
        if m <= 0.0 { return Err("logicle: M must be positive".into()); }
        if 2.0 * w > m { return Err("logicle: W is too large (need 2W ≤ M)".into()); }
        if -a > w || a + w > m - w {
            return Err("logicle: A is out of range".into());
        }

        let w_n = w / (m + a);
        let x2 = a / (m + a);
        let x1 = x2 + w_n;
        let x0 = x2 + 2.0 * w_n;
        let b = (m + a) * LN_10;
        let d = solve(b, w_n);

        let c_a = (x0 * (b + d)).exp();
        let mf_a = (b * x1).exp() - c_a / (d * x1).exp();
        let a_c = t / ((b.exp() - mf_a) - c_a / d.exp());
        let c_c = c_a * a_c;
        let f_c = -mf_a * a_c;

        // Taylor series about x1 (used near the linearization point for accuracy)
        let x_taylor = x1 + w_n / 4.0;
        let mut pos_coef = a_c * (b * x1).exp();
        let mut neg_coef = -c_c / (d * x1).exp();
        let mut taylor = [0.0f64; TAYLOR_LENGTH];
        for i in 0..TAYLOR_LENGTH {
            pos_coef *= b / (i as f64 + 1.0);
            neg_coef *= -d / (i as f64 + 1.0);
            taylor[i] = pos_coef + neg_coef;
        }
        taylor[1] = 0.0; // exact: first-order term vanishes at x1 by construction

        Ok(Logicle {
            t, w, m, a,
            a_c, b_c: b, c_c, d_c: d, f_c,
            x1, x_taylor, taylor,
        })
    }

    fn series_biexponential(&self, scale: f64) -> f64 {
        // Horner evaluation of the Taylor series (skips taylor[1] == 0).
        let x = scale - self.x1;
        let mut sum = self.taylor[TAYLOR_LENGTH - 1] * x;
        let mut i = TAYLOR_LENGTH as isize - 2;
        while i >= 2 {
            sum = (sum + self.taylor[i as usize]) * x;
            i -= 1;
        }
        (sum * x + self.taylor[0]) * x
    }

    /// data value → display coordinate.
    pub fn scale(&self, value: f64) -> f64 {
        if value == 0.0 {
            return self.x1;
        }
        let negative = value < 0.0;
        let value = if negative { -value } else { value };

        // initial guess
        let mut x = if value < self.f_c {
            self.x1 + value / self.taylor[0]
        } else {
            (value / self.a_c).ln() / self.b_c
        };

        let mut tolerance = 3.0 * f64::EPSILON;
        if x > 1.0 {
            tolerance = 3.0 * x * f64::EPSILON;
        }

        // Halley's method (cubic convergence)
        for _ in 0..20 {
            let ae2bx = self.a_c * (self.b_c * x).exp();
            let ce2mdx = self.c_c / (self.d_c * x).exp();
            let y = if x < self.x_taylor {
                self.series_biexponential(x) - value
            } else {
                (ae2bx + self.f_c) - (ce2mdx + value)
            };
            let abe2bx = self.b_c * ae2bx;
            let cde2mdx = self.d_c * ce2mdx;
            let dy = abe2bx + cde2mdx;
            let ddy = self.b_c * abe2bx - self.d_c * cde2mdx;
            let delta = y / (dy * (1.0 - y * ddy / (2.0 * dy * dy)));
            x -= delta;
            if delta.abs() < tolerance {
                return if negative { 2.0 * self.x1 - x } else { x };
            }
        }
        // Did not fully converge — return best estimate rather than panic.
        if negative { 2.0 * self.x1 - x } else { x }
    }

    /// display coordinate → data value (closed form).
    pub fn inverse(&self, scale: f64) -> f64 {
        let negative = scale < self.x1;
        let scale = if negative { 2.0 * self.x1 - scale } else { scale };
        let inverse = if scale < self.x_taylor {
            self.series_biexponential(scale)
        } else {
            (self.a_c * (self.b_c * scale).exp() + self.f_c)
                - self.c_c / (self.d_c * scale).exp()
        };
        if negative { -inverse } else { inverse }
    }
}

/// Solve  2·(ln(x) − ln(b)) + w·(b + x) = 0  for x in (0, b].
/// The function is strictly increasing in x, so bisection converges reliably.
fn solve(b: f64, w: f64) -> f64 {
    if w == 0.0 {
        return b; // logicle degenerates to asinh
    }
    let f = |x: f64| 2.0 * (x.ln() - b.ln()) + w * (b + x);
    let mut lo = f64::MIN_POSITIVE;
    let mut hi = b;
    let tol = 2.0 * b * f64::EPSILON;
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if f(mid) > 0.0 { hi = mid; } else { lo = mid; }
        if hi - lo <= tol { break; }
    }
    0.5 * (lo + hi)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_trip() {
        let lg = Logicle::new(262144.0, 0.5, 4.5, 0.0).unwrap();
        for &v in &[-1000.0, -10.0, 0.0, 10.0, 1000.0, 100000.0, 262144.0] {
            let s = lg.scale(v);
            let back = lg.inverse(s);
            assert!((back - v).abs() < 1e-3 * (v.abs() + 1.0),
                "round-trip failed for {}: got {}", v, back);
        }
    }

    #[test]
    fn monotonic() {
        let lg = Logicle::new(262144.0, 0.5, 4.5, 0.0).unwrap();
        let mut last = f64::NEG_INFINITY;
        let mut v = -10000.0;
        while v <= 262144.0 {
            let s = lg.scale(v);
            assert!(s > last, "not monotonic at {}", v);
            last = s;
            v += 137.0;
        }
    }
}
