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

    /// Compile into a runtime engine (builds the Logicle scale once).
    /// Falls back to Linear if logicle parameters are invalid.
    pub fn compile(&self) -> CompiledTransform {
        match self {
            AxisTransform::Linear => CompiledTransform::Linear,
            AxisTransform::Log { floor } => CompiledTransform::Log { floor: *floor },
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
