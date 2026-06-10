use serde::{Deserialize, Serialize};

use crate::logicle::Logicle;

/// Default logicle parameters for a BD instrument (24-bit linear range).
pub const DEFAULT_LOGICLE_T: f64 = 262144.0;
pub const DEFAULT_LOGICLE_W: f64 = 0.5;
pub const DEFAULT_LOGICLE_M: f64 = 4.5;
pub const DEFAULT_LOGICLE_A: f64 = 0.0;

/// A per-axis display transform. Serializable so it can be stored with gates.
///
/// `forward`  : data value → display coordinate (what the plot shows)
/// `inverse`  : display coordinate → data value (for gate clicks, tick labels)
#[derive(Clone, Debug, Default, Serialize, Deserialize, PartialEq)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum AxisTransform {
    /// Identity — for FSC/SSC/Time.
    #[default]
    Linear,
    /// log10, with values ≤ `floor` clamped (compensated data has negatives).
    Log { floor: f64 },
    /// asinh(x / cofactor).
    Asinh { cofactor: f64 },
    /// Logicle / biexponential (Moore & Parks).
    Logicle { t: f64, w: f64, m: f64, a: f64 },
}

impl AxisTransform {
    pub fn default_logicle() -> Self {
        AxisTransform::Logicle {
            t: DEFAULT_LOGICLE_T,
            w: DEFAULT_LOGICLE_W,
            m: DEFAULT_LOGICLE_M,
            a: DEFAULT_LOGICLE_A,
        }
    }

    pub fn default_log() -> Self {
        AxisTransform::Log { floor: 1.0 }
    }

    pub fn short_label(&self) -> &'static str {
        match self {
            AxisTransform::Linear => "Linear",
            AxisTransform::Log { .. } => "Log",
            AxisTransform::Asinh { .. } => "Asinh",
            AxisTransform::Logicle { .. } => "Logicle",
        }
    }

    /// Cache key that DISTINGUISHES transform parameters (unlike `short_label`, which is
    /// param-free). Plot/stat caches key on this so editing a cofactor / logicle param /
    /// log floor invalidates them — `short_label` alone would leave them silently stale.
    pub fn key(&self) -> String {
        match self {
            AxisTransform::Linear => "Linear".into(),
            AxisTransform::Log { floor } => format!("Log:{floor}"),
            AxisTransform::Asinh { cofactor } => format!("Asinh:{cofactor}"),
            AxisTransform::Logicle { t, w, m, a } => format!("Logicle:{t},{w},{m},{a}"),
        }
    }

    /// Compile into a runtime engine (builds the Logicle scale once).
    /// Falls back to Linear if logicle parameters are invalid.
    pub fn compile(&self) -> CompiledTransform {
        match self {
            AxisTransform::Linear => CompiledTransform::Linear,
            // Clamp floor positive so log10 stays finite even if a deserialized session
            // carries a non-positive floor (the UI only ever sets floor = 1.0).
            AxisTransform::Log { floor } => CompiledTransform::Log { floor: floor.max(f64::MIN_POSITIVE) },
            AxisTransform::Asinh { cofactor } => {
                CompiledTransform::Asinh { cofactor: *cofactor }
            }
            AxisTransform::Logicle { t, w, m, a } => match Logicle::new(*t, *w, *m, *a) {
                Ok(lg) => CompiledTransform::Logicle(Box::new(lg)),
                Err(_) => CompiledTransform::Linear,
            },
        }
    }
}

/// Runtime form of a transform — holds the expensive-to-build Logicle engine.
#[derive(Clone)]
pub enum CompiledTransform {
    Linear,
    Log { floor: f64 },
    Asinh { cofactor: f64 },
    Logicle(Box<Logicle>),
}

impl CompiledTransform {
    /// data value → display coordinate.
    ///
    /// Logicle is scaled by M so output is in [0, M] decades, matching
    /// flowCore::logicleTransform / FlowJo (the reference `Logicle::scale`
    /// itself returns the normalized [0,1] coordinate).
    pub fn forward(&self, x: f64) -> f64 {
        match self {
            CompiledTransform::Linear => x,
            CompiledTransform::Log { floor } => x.max(*floor).log10(),
            CompiledTransform::Asinh { cofactor } => (x / cofactor).asinh(),
            CompiledTransform::Logicle(lg) => lg.scale(x) * lg.m,
        }
    }

    /// display coordinate → data value.
    pub fn inverse(&self, y: f64) -> f64 {
        match self {
            CompiledTransform::Linear => y,
            CompiledTransform::Log { .. } => 10f64.powf(y),
            CompiledTransform::Asinh { cofactor } => y.sinh() * cofactor,
            CompiledTransform::Logicle(lg) => lg.inverse(y / lg.m),
        }
    }
}

/// Asinh transform: the de-facto standard for flow cytometry visualisation
/// when a full logicle implementation isn't required.
///
///   f(x) = asinh(x / cofactor)
///
/// Cofactor guidance (choose based on instrument noise floor):
///   150   — common for older cytometers or pre-gated data
///   5      — for mass cytometry (CyTOF)
///
/// This matches the `asinh` transform in flowCore / Bioconductor.
pub fn asinh(x: f64, cofactor: f64) -> f64 {
    (x / cofactor).asinh()
}

/// Apply asinh in-place to selected column indices across all events.
///
/// `events` is a flat row-major array [n_events × n_params].
/// Only the channels at `indices` are transformed; scatter channels are left as-is.
pub fn apply_asinh(events: &mut [f64], n_params: usize, indices: &[usize], cofactor: f64) {
    for ev_base in (0..events.len()).step_by(n_params) {
        for &idx in indices {
            events[ev_base + idx] = asinh(events[ev_base + idx], cofactor);
        }
    }
}

/// Determine which parameter indices should be transformed (i.e., fluorescence channels).
///
/// Heuristic: exclude parameters whose names start with FSC, SSC, or Time/time,
/// since those are scatter / time channels that should stay in linear scale.
pub fn fluorescence_indices(params: &[crate::fcs::Parameter]) -> Vec<usize> {
    params
        .iter()
        .enumerate()
        .filter(|(_, p)| {
            let n = p.name.to_uppercase();
            !n.starts_with("FSC")
                && !n.starts_with("SSC")
                && !n.eq_ignore_ascii_case("TIME")
        })
        .map(|(i, _)| i)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_util::param;

    fn round_trip(t: &AxisTransform, x: f64) {
        let c = t.compile();
        let back = c.inverse(c.forward(x));
        assert!(
            (back - x).abs() < 1e-6 * (x.abs() + 1.0),
            "{:?} round-trip failed: {} → {}",
            t, x, back
        );
    }

    #[test]
    fn linear_round_trip() {
        for &x in &[-100.0, 0.0, 42.0, 1e5] {
            round_trip(&AxisTransform::Linear, x);
        }
    }

    #[test]
    fn log_round_trip_and_floor() {
        let t = AxisTransform::Log { floor: 1.0 };
        for &x in &[1.0, 10.0, 1000.0, 262144.0] {
            round_trip(&t, x);
        }
        // Values below the floor clamp to log10(floor) = 0.
        assert_eq!(t.compile().forward(0.5), 0.0);
        assert_eq!(t.compile().forward(-50.0), 0.0);
    }

    #[test]
    fn log_floor_non_positive_is_clamped_finite() {
        // Audit L5: a deserialized Log with floor <= 0 must not emit -inf/NaN.
        let c = AxisTransform::Log { floor: 0.0 }.compile();
        assert!(c.forward(-100.0).is_finite());
        assert!(c.forward(0.0).is_finite());
        assert!(AxisTransform::Log { floor: -5.0 }.compile().forward(1.0).is_finite());
    }

    #[test]
    fn asinh_round_trip() {
        let t = AxisTransform::Asinh { cofactor: 150.0 };
        for &x in &[-1000.0, -10.0, 0.0, 10.0, 1000.0, 262144.0] {
            round_trip(&t, x);
        }
    }

    #[test]
    fn asinh_matches_free_function() {
        let c = CompiledTransform::Asinh { cofactor: 5.0 };
        assert!((c.forward(37.0) - asinh(37.0, 5.0)).abs() < 1e-12);
    }

    #[test]
    fn logicle_round_trip_and_m_scaling() {
        let t = AxisTransform::default_logicle();
        for &x in &[-1000.0, 0.0, 100.0, 10000.0, 262144.0] {
            round_trip(&t, x);
        }
        // forward output should land in the [0, M] decade range.
        let c = t.compile();
        let top = c.forward(262144.0);
        assert!((top - DEFAULT_LOGICLE_M).abs() < 1e-3, "top decade ≈ M, got {}", top);
    }

    #[test]
    fn logicle_invalid_params_fall_back_to_linear() {
        // W too large for M makes Logicle::new fail → compile() falls back to Linear.
        let bad = AxisTransform::Logicle { t: 262144.0, w: 100.0, m: 4.5, a: 0.0 };
        let c = bad.compile();
        assert!(matches!(c, CompiledTransform::Linear));
        // And forward is then identity (no panic).
        assert_eq!(c.forward(123.0), 123.0);
    }

    #[test]
    fn apply_asinh_only_touches_selected_indices() {
        // 2 events × 3 params; transform only column 1.
        let mut events = vec![10.0, 100.0, 1000.0, 20.0, 200.0, 2000.0];
        apply_asinh(&mut events, 3, &[1], 150.0);
        assert_eq!(events[0], 10.0, "col 0 untouched");
        assert_eq!(events[2], 1000.0, "col 2 untouched");
        assert!((events[1] - asinh(100.0, 150.0)).abs() < 1e-12);
        assert!((events[4] - asinh(200.0, 150.0)).abs() < 1e-12);
    }

    #[test]
    fn fluorescence_indices_excludes_scatter_and_time() {
        let params = vec![
            param(1, "FSC-A"),
            param(2, "SSC-H"),
            param(3, "FITC-A"),
            param(4, "PE-A"),
            param(5, "Time"),
        ];
        assert_eq!(fluorescence_indices(&params), vec![2, 3]);
    }

    #[test]
    fn axis_transform_serde_round_trip() {
        for t in [
            AxisTransform::Linear,
            AxisTransform::default_log(),
            AxisTransform::Asinh { cofactor: 150.0 },
            AxisTransform::default_logicle(),
        ] {
            let json = serde_json::to_string(&t).unwrap();
            let back: AxisTransform = serde_json::from_str(&json).unwrap();
            assert_eq!(back, t);
        }
    }

    #[test]
    fn short_label_is_stable() {
        assert_eq!(AxisTransform::Linear.short_label(), "Linear");
        assert_eq!(AxisTransform::default_log().short_label(), "Log");
        assert_eq!(AxisTransform::Asinh { cofactor: 1.0 }.short_label(), "Asinh");
        assert_eq!(AxisTransform::default_logicle().short_label(), "Logicle");
    }

    #[test]
    fn key_distinguishes_parameters() {
        // The cache key MUST change with parameters (unlike short_label) — else plot/stat
        // caches go stale on a cofactor/logicle edit (audit M1).
        assert_ne!(
            AxisTransform::Asinh { cofactor: 150.0 }.key(),
            AxisTransform::Asinh { cofactor: 200.0 }.key()
        );
        assert_eq!(
            AxisTransform::Asinh { cofactor: 150.0 }.key(),
            AxisTransform::Asinh { cofactor: 150.0 }.key()
        );
        let lg = |w: f64| AxisTransform::Logicle { t: 262144.0, w, m: 4.5, a: 0.0 }.key();
        assert_ne!(lg(0.5), lg(1.0));
        // short_label stays param-free (those two would collide under it).
        assert_eq!(
            AxisTransform::Asinh { cofactor: 150.0 }.short_label(),
            AxisTransform::Asinh { cofactor: 200.0 }.short_label()
        );
    }
}
