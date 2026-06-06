//! Shared constructors for unit tests across modules.
//!
//! Compiled only under `cfg(test)`. Lets every module build an in-memory
//! `FcsFile` from a row-major event matrix without touching the disk.

use std::collections::HashMap;

use crate::fcs::{FcsFile, Parameter};

/// Build a `Parameter` with sensible defaults (32-bit float, 262144 range).
pub fn param(index: usize, name: &str) -> Parameter {
    Parameter {
        index,
        name: name.to_string(),
        label: None,
        range: 262144.0,
        bits: 32,
    }
}

/// Build an in-memory `FcsFile` from channel names and a row-major event matrix
/// (`rows[event][param]`). Panics if the rows are ragged.
pub fn make_fcs(names: &[&str], rows: &[Vec<f64>]) -> FcsFile {
    let n_params = names.len();
    let n_events = rows.len();
    let parameters: Vec<Parameter> =
        names.iter().enumerate().map(|(i, n)| param(i + 1, n)).collect();

    let mut events = Vec::with_capacity(n_params * n_events);
    for r in rows {
        assert_eq!(r.len(), n_params, "ragged event row");
        events.extend_from_slice(r);
    }

    // Populate the keywords a real parsed file would carry, so code that reads
    // them back (e.g. the FCS writer copying through $PnN) behaves realistically.
    let mut keywords = HashMap::new();
    keywords.insert("$PAR".to_string(), n_params.to_string());
    keywords.insert("$TOT".to_string(), n_events.to_string());
    keywords.insert("$DATATYPE".to_string(), "F".to_string());
    keywords.insert("$BYTEORD".to_string(), "1,2,3,4".to_string());
    keywords.insert("$MODE".to_string(), "L".to_string());
    for (i, n) in names.iter().enumerate() {
        keywords.insert(format!("$P{}N", i + 1), n.to_string());
        keywords.insert(format!("$P{}B", i + 1), "32".to_string());
        keywords.insert(format!("$P{}R", i + 1), "262144".to_string());
    }

    FcsFile {
        version: "FCS3.0".to_string(),
        keywords,
        parameters,
        n_events,
        events,
    }
}
