use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::mpsc::{Receiver, TryRecvError};
use std::sync::Arc;

use egui::{Color32, RichText, Stroke};
use egui_extras::{Column, TableBuilder};
use egui_phosphor::regular as icon;
use egui_plot::{
    GridMark, Legend, Line, Plot, PlotBounds, PlotImage, PlotPoint, PlotPoints, Points,
    Polygon as PlotPolygon, Text as PlotText,
};

use crate::compensation::{
    compute_spillover, fluor_token_in_filename, format_spillover, load_matrix_file,
    max_off_diagonal, parse_spillover, save_matrix_file, SpilloverMatrix,
};
use crate::fcs_write;
use crate::fcs::FcsFile;
use crate::qc;
use crate::gating::{compute_own_masks, effective_mask, BoolOp, Gate, GateShape};
use crate::popstats::{
    append_long_csv, append_long_csv_grouped, population_stats, PopulationStatsTable,
    LONG_CSV_HEADER, LONG_CSV_HEADER_GROUPED,
};
use crate::transform::{AxisTransform, CompiledTransform};

const MAX_SCATTER: usize = 20_000;
const DENSITY_BINS: usize = 160;
const N_BUCKETS: usize = 14;
/// Squared normalized distance (fraction of plot size) within which a click on the
/// first polygon vertex closes the polygon — generous so it's easy to hit.
const POLY_CLOSE_SQ: f64 = 0.0036; // ~6% of the plot's width/height
/// Rect/ellipse drags smaller than this (squared, normalized) are ignored, so a
/// stray click-with-jitter doesn't create a tiny accidental gate.
const MIN_DRAG_SQ: f64 = 0.00035; // ~1.9% of the plot's width/height
const QC_MIN_EVENTS: usize = 5_000;   // tubes below this are flagged in the Samples list
const QC_MIN_VIABLE: f64 = 50.0;      // %viable below this is flagged
const QC_FLOW_DEV: f64 = 40.0;        // flow-rate bin deviation % above this is flagged
const QC_MARGIN_PCT: f64 = 1.0;       // saturated/margin % above this is flagged
const QC_FLOW_BINS: usize = 20;       // time bins for the flow-rate check
const REF_OVERLAY_MAX: usize = 8_000; // points drawn for a reference-overlay sample

/// Compact event count: 141 → "141", 13897 → "13.9k", 354295 → "354.3k".
fn fmt_count(n: usize) -> String {
    if n >= 1_000_000 { format!("{:.1}M", n as f64 / 1e6) }
    else if n >= 1000 { format!("{:.1}k", n as f64 / 1e3) }
    else { n.to_string() }
}

/// Tiny inline bar sparkline of per-time-bin event counts — a clog shows as a dip.
fn flow_sparkline(ui: &mut egui::Ui, bins: &[usize], flagged: bool) -> egui::Response {
    let (w, h) = (66.0_f32, 16.0_f32);
    let (rect, resp) = ui.allocate_exact_size(egui::vec2(w, h), egui::Sense::hover());
    if !bins.is_empty() {
        let mx = bins.iter().copied().max().unwrap_or(1).max(1) as f32;
        let bw = w / bins.len() as f32;
        let col = if flagged { Color32::from_rgb(220, 150, 50) } else { Color32::from_rgb(120, 165, 120) };
        let p = ui.painter();
        for (i, &c) in bins.iter().enumerate() {
            let bh = (c as f32 / mx) * h;
            let x = rect.left() + i as f32 * bw;
            p.rect_filled(
                egui::Rect::from_min_size(egui::pos2(x, rect.bottom() - bh), egui::vec2((bw - 1.0).max(1.0), bh)),
                0.0, col,
            );
        }
    }
    resp
}

// ── Colors ────────────────────────────────────────────────────────────────

/// Density-scatter colormap. Viridis/Magma/Cividis/Turbo are perceptually-uniform,
/// colorblind-safe maps (Viridis is the publication default); Jet is the legacy
/// rainbow, kept for familiarity but perceptually misleading.
#[derive(Clone, Copy, PartialEq)]
enum ColorMap { Viridis, Magma, Turbo, Cividis, Jet }

impl ColorMap {
    const ALL: [ColorMap; 5] =
        [ColorMap::Viridis, ColorMap::Magma, ColorMap::Turbo, ColorMap::Cividis, ColorMap::Jet];
    fn label(self) -> &'static str {
        match self {
            ColorMap::Viridis => "Viridis", ColorMap::Magma => "Magma", ColorMap::Turbo => "Turbo",
            ColorMap::Cividis => "Cividis", ColorMap::Jet => "Jet",
        }
    }
    fn from_name(s: &str) -> Option<Self> {
        ColorMap::ALL.into_iter().find(|c| c.label() == s)
    }
}

fn density_color(bucket: usize, n: usize, dark: bool, cmap: ColorMap) -> Color32 {
    let t = bucket as f32 / (n.saturating_sub(1).max(1)) as f32;
    density_color_t(t, dark, cmap)
}

/// Colormap lookup for a continuous density value t ∈ [0,1] (used by the filled heatmap).
fn density_color_t(t: f32, dark: bool, cmap: ColorMap) -> Color32 {
    let (r, g, b) = match cmap {
        ColorMap::Jet => return jet_color(t, dark),
        ColorMap::Viridis => lerp_anchors(t, &VIRIDIS),
        ColorMap::Magma => lerp_anchors(t, &MAGMA),
        ColorMap::Cividis => lerp_anchors(t, &CIVIDIS),
        ColorMap::Turbo => turbo(t),
    };
    Color32::from_rgba_unmultiplied(r, g, b, 220)
}

/// Piecewise-linear interpolation across an anchor table.
fn lerp_anchors(t: f32, a: &[(u8, u8, u8)]) -> (u8, u8, u8) {
    let n = a.len();
    let x = t.clamp(0.0, 1.0) * (n - 1) as f32;
    let i = (x.floor() as usize).min(n - 2);
    let f = x - i as f32;
    let m = |p: u8, q: u8| (p as f32 + (q as f32 - p as f32) * f) as u8;
    (m(a[i].0, a[i + 1].0), m(a[i].1, a[i + 1].1), m(a[i].2, a[i + 1].2))
}

// 9-anchor samples of the matplotlib maps, sampled at t = k/8 from the reference
// 256-entry tables (verified against matplotlib's _cm_listed.py: max Δ ≤ 1/255).
const VIRIDIS: [(u8, u8, u8); 9] = [
    (68, 1, 84), (71, 45, 123), (59, 82, 139), (44, 114, 142), (33, 145, 140),
    (39, 173, 129), (92, 200, 99), (170, 220, 50), (253, 231, 37),
];
const MAGMA: [(u8, u8, u8); 9] = [
    (0, 0, 4), (29, 17, 71), (81, 18, 124), (131, 38, 129), (183, 55, 121),
    (229, 80, 100), (251, 135, 97), (254, 194, 135), (252, 253, 191),
];
const CIVIDIS: [(u8, u8, u8); 9] = [
    (0, 34, 78), (26, 56, 111), (67, 78, 108), (97, 101, 111), (125, 124, 120),
    (154, 147, 118), (187, 173, 109), (221, 200, 88), (254, 232, 56),
];

/// Turbo via the standard degree-5 polynomial approximation (Mikhailov/Google).
/// Evaluated in f64 so the published coefficients keep full precision.
fn turbo(t: f32) -> (u8, u8, u8) {
    let t = t.clamp(0.0, 1.0) as f64;
    let (t2, t3, t4, t5) = (t * t, t * t * t, t * t * t * t, t * t * t * t * t);
    let r = 0.13572138 + 4.61539260 * t - 42.66032258 * t2 + 132.13108234 * t3 - 152.94239396 * t4 + 59.28637943 * t5;
    let g = 0.09140261 + 2.19418839 * t + 4.84296658 * t2 - 14.18503333 * t3 + 4.27729857 * t4 + 2.82956604 * t5;
    let b = 0.10667330 + 12.64194608 * t - 60.58204836 * t2 + 110.36276771 * t3 - 89.90310912 * t4 + 27.34824973 * t5;
    let u = |x: f64| (x.clamp(0.0, 1.0) * 255.0) as u8;
    (u(r), u(g), u(b))
}

fn jet_color(t: f32, dark: bool) -> Color32 {
    let (r, g, b) = if t < 0.25 {
        let s = t / 0.25; (0.0, s, 1.0)
    } else if t < 0.5 {
        let s = (t - 0.25) / 0.25; (0.0, 1.0, 1.0 - s)
    } else if t < 0.75 {
        let s = (t - 0.5) / 0.25; (s, 1.0, 0.0)
    } else {
        let s = (t - 0.75) / 0.25; (1.0, 1.0 - s * 0.8, 0.0)
    };
    // On a light background, lift the floor so the lowest-density blue isn't invisible.
    let floor = if dark { 0.0 } else { 0.15 };
    let scale = 1.0 - floor;
    Color32::from_rgba_unmultiplied(
        ((floor + r * scale) * 255.0) as u8,
        ((floor + g * scale) * 255.0) as u8,
        ((floor + b * scale) * 255.0) as u8,
        220,
    )
}

/// Heat tint for a spillover-matrix cell: neutral gray on the diagonal,
/// orange→red proportional to spillover magnitude off-diagonal (saturating at ~25%).
fn spill_cell_color(v: f64, is_diag: bool, dark: bool) -> Color32 {
    if is_diag {
        return if dark { Color32::from_rgb(55, 62, 72) } else { Color32::from_rgb(205, 212, 222) };
    }
    let a = (v.abs() / 0.25).clamp(0.0, 1.0);
    if a < 0.01 { return Color32::TRANSPARENT; }
    Color32::from_rgba_unmultiplied(230, 90, 40, (a * 210.0) as u8)
}

/// Stable color for a gate, keyed by its (immutable) `id` so a population keeps
/// its color across deletes/reorders — not by list position, which would
/// reshuffle every other gate's color whenever one is removed.
fn gate_color(id: u32) -> (Color32, Color32) {
    const BASES: [(u8, u8, u8); 8] = [
        (220, 40, 40), (30, 160, 30), (40, 90, 230),
        (210, 130, 0), (170, 0, 170), (0, 150, 160),
        (230, 90, 0), (120, 0, 220),
    ];
    let (r, g, b) = BASES[(id as usize) % 8];
    (Color32::from_rgb(r, g, b), Color32::from_rgba_unmultiplied(r, g, b, 26))
}

// ── State ─────────────────────────────────────────────────────────────────

#[derive(PartialEq, Clone, Copy)]
enum DrawMode { Navigate, Rect, Ellipse, Polygon, Quadrant, Edit }

#[derive(PartialEq, Clone, Copy)]
enum ActiveTab { Plot, Histogram, Stats, Batch, Spillover, QC }

#[derive(PartialEq, Clone, Copy)]
enum BatchMetric { PctParent, PctTotal, Mfi }

impl BatchMetric {
    fn label(self) -> &'static str {
        match self { BatchMetric::PctParent => "% parent", BatchMetric::PctTotal => "% total", BatchMetric::Mfi => "MFI" }
    }
}

/// Cached density buckets for one cell of the 2×2 grid view.
struct GridCell {
    xi: usize,
    yi: usize,
    x_label: String,
    y_label: String,
    pop: Option<u32>,               // active "gate from here" population this was built for
    gen: u64,                       // data generation this was built for
    buckets: Vec<Vec<[f64; 2]>>,
}

struct ScatterCache {
    buckets: Vec<Vec<[f64; 2]>>,
    x_ch: usize,
    y_ch: usize,
    x_label: String,  // transform short label, to detect changes
    y_label: String,
    pop: Option<u32>, // active population the plot is restricted to (gate-from-here)
    back_pts: Vec<[f64; 2]>, // parent population events (greyed) when backgating
    contours: Vec<[[f64; 2]; 2]>, // iso-density contour line segments (display coords)
    heat: Option<HeatGrid>, // density grid for the filled-heatmap render mode
}

/// Per-bin normalized density (sqrt-scaled, 0 = empty) for the filled heatmap,
/// with the display-coord extent it spans.
struct HeatGrid {
    t: Vec<f32>,      // len n*n, index = by*n + bx
    n: usize,
    extent: [f64; 4], // [xmin, xmax, ymin, ymax] in display coords
}

/// A user-supplied / edited spillover matrix that overrides the embedded one.
#[derive(Clone, serde::Serialize, serde::Deserialize)]
struct SpillOverride {
    channels: Vec<String>,
    rows: Vec<Vec<f64>>,
}

/// A saved analysis session — enough to reopen the workspace and resume.
#[derive(serde::Serialize, serde::Deserialize)]
struct Session {
    sample_paths: Vec<PathBuf>,
    #[serde(default)]
    groups: Vec<String>,
    #[serde(default)]
    active_sample: usize,
    do_compensate: bool,
    #[serde(default)]
    dark_mode: bool,
    #[serde(default)]
    viridis: bool, // legacy (pre-0.1.9); superseded by `colormap`
    #[serde(default)]
    colormap: Option<String>,
    channel_tf: Vec<AxisTransform>,
    x_ch: usize,
    y_ch: usize,
    #[serde(default)]
    hist_ch: usize,
    gates: Vec<Gate>,
    #[serde(default)]
    spill_override: Option<SpillOverride>,
}

#[derive(Clone, Copy, PartialEq)]
enum HistNorm { Modal, Count }

#[derive(Clone, Copy, PartialEq)]
enum HistMode { Populations, Samples }

struct HistSeries {
    name: String,
    color: Color32,
    values: Vec<f64>, // per-bin, normalized per `HistNorm`
}

struct HistogramData {
    x_ch: usize,
    x_label: String, // transform short-label, for staleness detection
    norm: HistNorm,
    mode: HistMode,
    sample_pop: Option<u32>,
    centers: Vec<f64>, // bin centers in display coords
    series: Vec<HistSeries>,
}

/// A sample in the workspace. Only the *active* sample's events are held in
/// memory (`self.fcs`/`self.compensated`); the rest are just path + name.
struct SampleRef {
    path: PathBuf,
    name: String,
    n_events: Option<usize>, // from a lightweight $TOT read, for QC display
    group: String,           // user condition/group tag, carried into the batch CSV
}

/// Cached points of a reference sample overlaid behind the active scatter.
struct RefScatter {
    ref_idx: usize,
    x_ch: usize,
    y_ch: usize,
    x_label: String,
    y_label: String,
    points: Vec<[f64; 2]>,
}

/// Result of a streamed batch run — per-sample population stat tables + skips.
struct BatchResult {
    tables: Vec<(String, String, PopulationStatsTable)>, // (group, sample, table)
    skipped: Vec<(String, String)>, // (sample, reason)
}

/// Messages streamed from the background batch worker to the UI thread.
enum BatchMsg {
    Progress { done: usize, total: usize, name: String },
    Table(String, String, PopulationStatsTable), // (group, sample, table)
    Skip(String, String),                        // (sample, reason)
    Done,
}

/// One row of the QC scan — per-tube acquisition-quality summary.
#[derive(Clone)]
struct QcRow {
    group: String,
    name: String,
    n_events: usize,
    viable: Option<f64>,              // designated live gate's % of parent (None = not set)
    flow_dev: Option<f64>,            // flow-rate max deviation %, None = no usable Time
    flow_bins: Vec<usize>,            // per-bin event counts (sparkline)
    worst_margin: Option<(String, f64)>, // (channel, % saturated)
    flags: Vec<String>,              // human-readable reasons, empty = pass
}

enum QcMsg {
    Progress { done: usize, total: usize, name: String },
    Row(Box<QcRow>),
    Skip(String, String),
    Done,
}

pub struct FlowCytoApp {
    fcs: Option<FcsFile>,            // the ACTIVE sample's parsed file
    file_path: Option<PathBuf>,
    compensated: Vec<f64>,           // active sample: raw → optional compensate (DATA space)

    recent_files: Vec<PathBuf>,      // recently opened FCS paths (most-recent first)

    // Workspace of samples (active one is loaded into `fcs`/`compensated`).
    samples: Vec<SampleRef>,
    active_sample: usize,
    batch: Option<BatchResult>,
    batch_channel: usize,            // channel index whose MFI shows in the batch table
    batch_plot_pop: Option<u32>,     // population charted across samples in the Batch tab (None = none)
    batch_plot_metric: BatchMetric,  // which stat the batch chart shows
    batch_rx: Option<Receiver<BatchMsg>>,        // streamed results from the worker thread
    batch_cancel: Option<Arc<AtomicBool>>,       // set to stop the worker early
    batch_progress: Option<(usize, usize)>,      // (done, total) while a batch runs
    // QC suite — per-tube acquisition quality, streamed like the batch.
    qc_live_gate: Option<u32>,                   // gate designated as the viability population
    qc_rows: Vec<QcRow>,
    qc_rx: Option<Receiver<QcMsg>>,
    qc_cancel: Option<Arc<AtomicBool>>,
    qc_progress: Option<(usize, usize)>,
    // Manual update check (the app's only network path; fires only on explicit click).
    update_rx: Option<Receiver<Result<crate::update::UpdateInfo, String>>>,
    update_msg: Option<String>,                  // last update-check result message
    update_found: Option<crate::update::UpdateInfo>, // Some when a newer release exists
    ref_sample: Option<usize>,       // reference sample overlaid behind the active scatter
    ref_scatter: Option<RefScatter>,

    do_compensate: bool,
    dark_mode: bool,
    colormap: ColorMap,
    cursor_label: Option<String>, // live data-coords readout under the plot cursor
    last_plot_rect: Option<egui::Rect>, // screen rect of the most recent plot (for PNG crop)
    screenshot_pending: bool,           // a "Save plot PNG" request is awaiting the captured frame
    screenshot_sent: bool,              // the viewport screenshot command has been dispatched
    #[cfg(target_os = "macos")]
    dock_icon_ticks: u8,                // frames over which to (re)assert the Dock icon

    /// When Some, this matrix replaces the embedded $SPILLOVER everywhere.
    spill_override: Option<SpillOverride>,

    channel_tf: Vec<AxisTransform>,  // per-channel display transform
    x_ch: usize,
    y_ch: usize,

    // manual axis limits (DATA units)
    x_manual: bool,
    x_lo: f64,
    x_hi: f64,
    y_manual: bool,
    y_lo: f64,
    y_hi: f64,

    scatter: Option<ScatterCache>,
    grid_mode: bool,                  // Plot tab: show a grid of plots
    grid_cols: usize,                 // grid columns (1..=6)
    grid_rows: usize,                 // grid rows (1..=6)
    grid_channels: Vec<(usize, usize)>, // per-cell (x, y) channel indices
    grid_cache: Vec<Option<GridCell>>,  // per-cell density-bucket cache
    active_grid_cell: Option<usize>,    // cell index owning the in-progress draw gesture
    data_gen: u64,                    // bumped whenever compensated data changes

    gates: Vec<Gate>,
    undo_stack: Vec<Vec<Gate>>, // gate-tree snapshots for undo (pre-mutation states)
    redo_stack: Vec<Vec<Gate>>,
    next_gate_id: u32,
    gate_counts: HashMap<u32, (usize, usize)>, // id → (n_in_effective, n_parent)
    new_gate_parent: Option<u32>,
    bool_op: BoolOp,                          // op for the boolean-gate builder
    bool_refs: std::collections::HashSet<u32>, // gates selected in the boolean-gate builder

    draw_mode: DrawMode,
    drag_start: Option<[f64; 2]>,
    drag_current: Option<[f64; 2]>,
    poly_vertices: Vec<[f64; 2]>,
    selected_gate: Option<u32>,   // gate selected in the tree → shown in the numeric inspector
    active_pop: Option<u32>,      // "gate from here": plot restricted to this population's events
    backgate: bool,               // show the active population's PARENT events greyed behind it
    show_contours: bool,          // overlay iso-density contour lines (Plot tab)
    density_fill: bool,           // render a filled density heatmap instead of dots
    lock_view: bool,              // freeze plot pan/zoom so the view holds still while gating
    fit_view: bool,              // one-shot: reframe the next plot to the data extent
    scatter_hidden: std::collections::HashSet<u32>, // gate ids hidden on the scatter overlay
    grab_handle: Option<usize>,   // which handle of the selected gate is being dragged (Edit mode)
    gate_move_last: Option<[f64; 2]>, // last cursor pos (gate-display coords) while dragging a gate body

    pop_stats: Option<PopulationStatsTable>, // cached per-population table (Stats tab)

    // Histogram tab state
    hist_ch: usize,                    // channel histogrammed (independent of the Plot X axis)
    hist_norm: HistNorm,
    hist_mode: HistMode,               // overlay populations (1 sample) or samples
    hist_sample_pop: Option<u32>,      // in Samples mode: which population to histogram (None = all)
    hist_hidden: std::collections::HashSet<u32>, // gate ids hidden in overlay
    hist_all_hidden: bool,
    hist_draw_interval: bool,
    hist_cache: Option<HistogramData>,

    active_tab: ActiveTab,
    status: String,

    needs_reprocess: bool,
    needs_rescatter: bool,
    needs_regate: bool,

    #[cfg(target_os = "macos")]
    menu: Option<MenuState>, // native macOS menu bar (built once the NSApp exists)
}

impl Default for FlowCytoApp {
    fn default() -> Self {
        FlowCytoApp {
            fcs: None, file_path: None, compensated: Vec::new(),
            recent_files: Vec::new(),
            samples: Vec::new(), active_sample: 0, batch: None, batch_channel: 0,
            batch_plot_pop: None, batch_plot_metric: BatchMetric::PctParent,
            batch_rx: None, batch_cancel: None, batch_progress: None,
            qc_live_gate: None, qc_rows: Vec::new(), qc_rx: None, qc_cancel: None, qc_progress: None,
            update_rx: None, update_msg: None, update_found: None,
            ref_sample: None, ref_scatter: None,
            do_compensate: false, dark_mode: true,
            colormap: ColorMap::Viridis,
            cursor_label: None,
            last_plot_rect: None,
            screenshot_pending: false,
            screenshot_sent: false,
            #[cfg(target_os = "macos")]
            dock_icon_ticks: 0,
            spill_override: None,
            channel_tf: Vec::new(), x_ch: 0, y_ch: 1,
            x_manual: false, x_lo: 0.0, x_hi: 262144.0,
            y_manual: false, y_lo: 0.0, y_hi: 262144.0,
            scatter: None,
            grid_mode: false,
            grid_cols: 2,
            grid_rows: 2,
            grid_channels: vec![(0, 1), (0, 2), (0, 3), (1, 2)],
            grid_cache: vec![None, None, None, None],
            active_grid_cell: None,
            data_gen: 0,
            gates: Vec::new(), undo_stack: Vec::new(), redo_stack: Vec::new(), next_gate_id: 1,
            gate_counts: HashMap::new(), new_gate_parent: None,
            bool_op: BoolOp::And, bool_refs: std::collections::HashSet::new(),
            draw_mode: DrawMode::Navigate,
            drag_start: None, drag_current: None, poly_vertices: Vec::new(),
            selected_gate: None,
            active_pop: None,
            backgate: false,
            show_contours: false,
            density_fill: false,
            lock_view: true,
            fit_view: false,
            scatter_hidden: std::collections::HashSet::new(),
            grab_handle: None,
            gate_move_last: None,
            pop_stats: None,
            hist_ch: 0,
            hist_norm: HistNorm::Modal,
            hist_mode: HistMode::Populations,
            hist_sample_pop: None,
            hist_hidden: std::collections::HashSet::new(),
            hist_all_hidden: false,
            hist_draw_interval: false,
            hist_cache: None,
            active_tab: ActiveTab::Plot,
            status: "Open an FCS file to begin.".into(),
            needs_reprocess: false, needs_rescatter: false, needs_regate: false,
            #[cfg(target_os = "macos")]
            menu: None,
        }
    }
}

// ── Data loading & pipeline ───────────────────────────────────────────────

impl FlowCytoApp {
    /// Convenience for launching with one file (CLI initial file).
    pub fn load_file(&mut self, path: &Path) {
        self.add_files(vec![path.to_path_buf()]);
    }

    /// Add files to the workspace. If the workspace was empty, the first becomes
    /// active with a fresh panel setup; otherwise they are appended (analysis kept).
    pub fn add_files(&mut self, paths: Vec<PathBuf>) {
        self.push_recent(&paths);
        let was_empty = self.samples.is_empty();
        for p in paths {
            let name = p.file_stem().map(|s| s.to_string_lossy().to_string())
                .unwrap_or_else(|| "sample".into());
            let n_events = FcsFile::peek_events(&p).ok();
            self.samples.push(SampleRef { path: p, name, n_events, group: String::new() });
        }
        self.batch = None;
        if was_empty && !self.samples.is_empty() {
            self.activate_sample(0, true);
        } else {
            self.status = format!("{} samples in workspace", self.samples.len());
        }
    }

    /// Record opened files at the front of the recent list (deduped, capped) and
    /// persist them.
    fn push_recent(&mut self, paths: &[PathBuf]) {
        for p in paths {
            self.recent_files.retain(|x| x != p);
            self.recent_files.insert(0, p.clone());
        }
        self.recent_files.truncate(12);
        self.save_recent();
    }

    fn save_recent(&self) {
        if let Some(store) = recent_store_path() {
            if let Some(dir) = store.parent() { let _ = std::fs::create_dir_all(dir); }
            let list: Vec<String> = self.recent_files.iter().map(|p| p.to_string_lossy().to_string()).collect();
            if let Ok(j) = serde_json::to_string_pretty(&list) { let _ = std::fs::write(&store, j); }
        }
    }

    /// Load a sample's events into the active slot.
    /// `fresh` = brand-new workspace (reset panel + gates); otherwise keep the
    /// gating tree/compensation and re-key transforms by channel name.
    fn activate_sample(&mut self, i: usize, fresh: bool) {
        if i >= self.samples.len() { return; }
        let path = self.samples[i].path.clone();
        let fcs = match FcsFile::open(&path) {
            Ok(f) => f,
            Err(e) => { self.status = format!("Error: {}", e); return; }
        };

        if fresh {
            self.do_compensate = fcs.spillover_keyword().is_some();
            self.channel_tf = default_transforms(&fcs);
            self.x_ch = 0;
            self.y_ch = 3.min(fcs.n_params().saturating_sub(1));
            self.gates.clear();
            self.gate_counts.clear();
            self.next_gate_id = 1;
            self.new_gate_parent = None;
            self.active_pop = None;
            self.selected_gate = None;
            self.spill_override = None;
        } else {
            // Carry transforms across by channel NAME (panels may differ in order).
            self.channel_tf = rekey_transforms(self.fcs.as_ref(), &self.channel_tf, &fcs);
            self.x_ch = self.x_ch.min(fcs.n_params().saturating_sub(1));
            self.y_ch = self.y_ch.min(fcs.n_params().saturating_sub(1));
            // Warn if gates reference channels this sample lacks.
            let missing = missing_gate_channels(&self.gates, &fcs);
            if !missing.is_empty() {
                self.status = format!("{} {} missing channel(s): {} — gates may not apply",
                    icon::WARNING, self.samples[i].name, missing.join(", "));
            } else {
                self.status = format!("Active: {} ({} events)", self.samples[i].name, fcs.n_events);
            }
        }

        self.x_manual = false;
        self.y_manual = false;
        self.scatter = None;
        self.file_path = Some(path);
        self.fcs = Some(fcs);
        self.active_sample = i;
        self.reprocess();
        self.needs_rescatter = true;
        self.needs_regate = true;
    }

    /// Apply the current compensation settings (override > embedded) to a file's
    /// raw events, returning a compensated-linear event buffer. Used for the active
    /// sample (reprocess) and for each streamed sample in a batch run.
    /// Returns Ok(raw) when no matrix is present (not an error), Err only when a
    /// matrix exists but fails to build/invert/apply — so callers never silently
    /// pass off raw data as "compensated".
    fn compensate_events(&self, fcs: &FcsFile) -> Result<Vec<f64>, String> {
        compensate_for(fcs, self.do_compensate, self.spill_override.as_ref())
    }

    fn reprocess(&mut self) {
        let fcs = match &self.fcs { Some(f) => f, None => return };
        match self.compensate_events(fcs) {
            Ok(ev) => self.compensated = ev,
            Err(e) => {
                self.status = format!("{} Compensation failed — showing RAW data: {}", icon::WARNING, e);
                self.compensated = fcs.events.clone();
            }
        }
        self.pop_stats = None; // data changed → population stats stale
        self.hist_cache = None;
        self.ref_scatter = None; // compensation changed → reference overlay stale
        self.data_gen = self.data_gen.wrapping_add(1); // invalidate grid caches
    }

    fn cur_tf(&self, ch: usize) -> AxisTransform {
        self.channel_tf.get(ch).cloned().unwrap_or(AxisTransform::Linear)
    }

    fn rebuild_scatter(&mut self) {
        let (n_events, n_params) = match &self.fcs {
            Some(f) => (f.n_events, f.n_params()), None => return,
        };
        if n_events == 0 || n_params == 0 { return; }
        // Guard: data not yet (re)processed for this file — avoid indexing panic.
        if self.compensated.len() < n_events * n_params { return; }

        let xi = self.x_ch.min(n_params - 1);
        let yi = self.y_ch.min(n_params - 1);
        let xt = self.cur_tf(xi).compile();
        let yt = self.cur_tf(yi).compile();

        // "Gate from here": restrict to the active population's events.
        let kept: Vec<usize> = match self.active_pop.map(|p| self.pop_mask(p)) {
            Some(m) => (0..n_events).filter(|&e| m.get(e).copied().unwrap_or(false)).collect(),
            None => (0..n_events).collect(),
        };
        let nk = kept.len();

        // Display-space coords of the kept events.
        let dx: Vec<f64> = kept.iter().map(|&e| xt.forward(self.compensated[e * n_params + xi])).collect();
        let dy: Vec<f64> = kept.iter().map(|&e| yt.forward(self.compensated[e * n_params + yi])).collect();

        let mut buckets: Vec<Vec<[f64; 2]>> = vec![Vec::new(); N_BUCKETS];
        let mut contours: Vec<[[f64; 2]; 2]> = Vec::new();
        let mut heat: Option<HeatGrid> = None;
        if nk > 0 {
            let (xmin, xmax) = data_range(&dx);
            let (ymin, ymax) = data_range(&dy);
            let hist = density_hist(&dx, &dy, DENSITY_BINS, xmin, xmax, ymin, ymax);
            let max_d = hist.iter().flat_map(|r| r.iter()).copied().max().unwrap_or(1).max(1);
            let n_sample = MAX_SCATTER.min(nk);
            let step = (nk / n_sample).max(1);
            for k in (0..nk).step_by(step) {
                let (x, y) = (dx[k], dy[k]);
                let bx = bin_of(x, xmin, xmax, DENSITY_BINS);
                let by = bin_of(y, ymin, ymax, DENSITY_BINS);
                let t = (hist[bx][by] as f64 / max_d as f64).sqrt();
                let b = ((t * (N_BUCKETS - 1) as f64) as usize).min(N_BUCKETS - 1);
                buckets[b].push([x, y]);
            }
            if self.show_contours {
                let levels: Vec<f64> = [0.04, 0.10, 0.20, 0.35, 0.55, 0.80]
                    .iter().map(|f| f * max_d as f64).collect();
                contours = contour_segments(&hist, DENSITY_BINS, xmin, xmax, ymin, ymax, &levels);
            }
            // Normalized density grid for the filled-heatmap mode.
            let mut t = vec![0.0f32; DENSITY_BINS * DENSITY_BINS];
            for bx in 0..DENSITY_BINS {
                for by in 0..DENSITY_BINS {
                    let c = hist[bx][by];
                    if c > 0 { t[by * DENSITY_BINS + bx] = ((c as f64 / max_d as f64).sqrt()) as f32; }
                }
            }
            heat = Some(HeatGrid { t, n: DENSITY_BINS, extent: [xmin, xmax, ymin, ymax] });
        }

        // Backgating: parent population's events (greyed) for context behind the child.
        let back_pts: Vec<[f64; 2]> = if self.backgate {
            if let Some(ap) = self.active_pop {
                let parent = self.gates.iter().find(|g| g.id == ap).and_then(|g| g.parent);
                let pmask = match parent { Some(pid) => self.pop_mask(pid), None => vec![true; n_events] };
                let idxs: Vec<usize> = (0..n_events).filter(|&e| pmask.get(e).copied().unwrap_or(false)).collect();
                let step = (idxs.len() / REF_OVERLAY_MAX.min(idxs.len()).max(1)).max(1);
                idxs.iter().step_by(step)
                    .map(|&e| [xt.forward(self.compensated[e * n_params + xi]), yt.forward(self.compensated[e * n_params + yi])])
                    .collect()
            } else { Vec::new() }
        } else { Vec::new() };

        self.scatter = Some(ScatterCache {
            buckets, x_ch: xi, y_ch: yi,
            x_label: self.cur_tf(xi).short_label().to_string(),
            y_label: self.cur_tf(yi).short_label().to_string(),
            pop: self.active_pop,
            back_pts,
            contours,
            heat,
        });
    }

    /// Density buckets for an arbitrary channel pair (grid cells). Respects the
    /// active "gate from here" population, same as the main scatter.
    /// Event indices shown in the grid: the active "gate from here" population, or
    /// all events. Computed once per frame and shared across cells (it's the only
    /// expensive shared step, so we avoid recomputing it per cell).
    fn grid_kept(&self) -> Vec<usize> {
        let n_events = self.fcs.as_ref().map(|f| f.n_events).unwrap_or(0);
        match self.active_pop.map(|p| self.pop_mask(p)) {
            Some(m) => (0..n_events).filter(|&e| m.get(e).copied().unwrap_or(false)).collect(),
            None => (0..n_events).collect(),
        }
    }

    fn compute_cell_buckets(&self, xi: usize, yi: usize, kept: &[usize], cap: usize) -> Vec<Vec<[f64; 2]>> {
        let (n_events, n_params) = match &self.fcs {
            Some(f) => (f.n_events, f.n_params()), None => return Vec::new(),
        };
        if n_events == 0 || n_params == 0 || self.compensated.len() < n_events * n_params {
            return Vec::new();
        }
        let xi = xi.min(n_params - 1);
        let yi = yi.min(n_params - 1);
        let xt = self.cur_tf(xi).compile();
        let yt = self.cur_tf(yi).compile();
        let dx: Vec<f64> = kept.iter().map(|&e| xt.forward(self.compensated[e * n_params + xi])).collect();
        let dy: Vec<f64> = kept.iter().map(|&e| yt.forward(self.compensated[e * n_params + yi])).collect();
        bucketize(&dx, &dy, cap)
    }

    /// The gate drawn on channels (x_base, y_base) whose shape contains the plot
    /// point `p` (in the plot's display coords). Prefers the deepest (most specific)
    /// matching gate, so double-clicking nested gates drills to the innermost.
    fn gate_at_point(&self, x_base: &str, y_base: &str, xt: &CompiledTransform, yt: &CompiledTransform, p: [f64; 2]) -> Option<u32> {
        let depths: HashMap<u32, usize> = crate::gating::gate_tree_order(&self.gates).into_iter().collect();
        let mut best: Option<(u32, usize)> = None;
        for g in &self.gates {
            if g.x_channel.eq_ignore_ascii_case(x_base) && g.y_channel.eq_ignore_ascii_case(y_base) {
                let gxt = g.x_transform.compile();
                let gyt = g.y_transform.compile();
                let gx = gxt.forward(xt.inverse(p[0]));
                let gy = gyt.forward(yt.inverse(p[1]));
                if g.shape.contains(gx, gy) {
                    let d = depths.get(&g.id).copied().unwrap_or(0);
                    if best.map(|(_, bd)| d > bd).unwrap_or(true) { best = Some((g.id, d)); }
                }
            }
        }
        best.map(|(id, _)| id)
    }

    /// Write the events in a population to a new FCS file (raw events + the sample's
    /// spillover, so it can be re-compensated downstream).
    fn export_gated_fcs(&mut self, gid: u32) {
        let mask = self.pop_mask(gid);
        let (subset, gname) = {
            let fcs = match &self.fcs { Some(f) => f, None => return };
            let n = fcs.n_params();
            let mut sub: Vec<f64> = Vec::new();
            let mut count = 0usize;
            for e in 0..fcs.n_events {
                if mask.get(e).copied().unwrap_or(false) {
                    sub.extend_from_slice(&fcs.events[e * n..(e + 1) * n]);
                    count += 1;
                }
            }
            if count == 0 { self.status = "That population has 0 events — nothing to export.".into(); return; }
            let subset = FcsFile {
                version: fcs.version.clone(),
                keywords: fcs.keywords.clone(),
                parameters: fcs.parameters.clone(),
                n_events: count,
                events: sub,
            };
            let gname = self.gates.iter().find(|g| g.id == gid).map(|g| g.name.clone()).unwrap_or_else(|| "gated".into());
            (subset, gname)
        };
        let default = self.file_path.as_ref().and_then(|p| p.file_stem())
            .map(|s| format!("{}_{}.fcs", s.to_string_lossy(), sanitize(&gname)))
            .unwrap_or_else(|| "gated.fcs".into());
        if let Some(path) = rfd::FileDialog::new().add_filter("FCS", &["fcs"]).set_file_name(default).save_file() {
            self.status = match fcs_write::write_fcs(&subset, None, &path) {
                Ok(_) => format!("Exported {} events ({}) → {}", subset.n_events, gname, path.display()),
                Err(e) => format!("Export error: {}", e),
            };
        }
    }

    /// Frame the plot on a gate: switch X/Y to the gate's channels and set manual
    /// axis limits to its extent (+margin), in data units.
    fn zoom_to_gate(&mut self, gid: u32) {
        let g = match self.gates.iter().find(|g| g.id == gid) { Some(g) => g.clone(), None => return };
        let bbox = match shape_display_bbox(&g.shape) {
            Some(b) => b,
            None => { self.status = "Interval gates have no 2-D extent to zoom to.".into(); return; }
        };
        let (xi, yi) = match &self.fcs {
            Some(f) => match (f.param_index(&g.x_channel), f.param_index(&g.y_channel)) {
                (Some(a), Some(b)) => (a, b),
                _ => { self.status = "Gate channels not in this sample.".into(); return; }
            },
            None => return,
        };
        let (dxmin, dxmax, dymin, dymax) = bbox;
        let xt = g.x_transform.compile();
        let yt = g.y_transform.compile();
        let mx = (dxmax - dxmin).abs().max(1e-9) * 0.15;
        let my = (dymax - dymin).abs().max(1e-9) * 0.15;
        let (x0, x1) = (xt.inverse(dxmin - mx), xt.inverse(dxmax + mx));
        let (y0, y1) = (yt.inverse(dymin - my), yt.inverse(dymax + my));
        self.x_ch = xi; self.y_ch = yi;
        self.x_lo = x0.min(x1); self.x_hi = x0.max(x1); self.x_manual = true;
        self.y_lo = y0.min(y1); self.y_hi = y0.max(y1); self.y_manual = true;
        self.scatter = None;
        self.needs_rescatter = true;
        self.status = format!("Zoomed to “{}”", g.name);
    }

    /// Drill into the gate under a double-clicked point ("gate from here"). Returns
    /// true if a gate was hit.
    fn drill_at(&mut self, x_name: &str, y_name: &str, xt: &CompiledTransform, yt: &CompiledTransform, p: [f64; 2]) -> bool {
        match self.gate_at_point(&x_name_base(x_name), &x_name_base(y_name), xt, yt, p) {
            Some(id) => {
                self.active_pop = Some(id);
                self.new_gate_parent = Some(id);
                self.selected_gate = Some(id);
                self.scatter = None;
                if let Some(g) = self.gates.iter().find(|g| g.id == id) {
                    self.status = format!("Gating from “{}” (double-click a child to go deeper)", g.name);
                }
                true
            }
            None => false,
        }
    }

    /// Ancestor path (root → `pop`) as (id, name), for the gate-from-here breadcrumb.
    fn population_path(&self, pop: Option<u32>) -> Vec<(u32, String)> {
        let mut chain = Vec::new();
        let mut cur = pop;
        let mut guard = 0;
        while let Some(id) = cur {
            guard += 1; if guard > 1000 { break; }
            match self.gates.iter().find(|g| g.id == id) {
                Some(g) => { chain.push((id, g.name.clone())); cur = g.parent; }
                None => break,
            }
        }
        chain.reverse();
        chain
    }

    /// Effective membership mask of a population (AND of the gate with its ancestors).
    fn pop_mask(&self, pop: u32) -> Vec<bool> {
        let fcs = match &self.fcs { Some(f) => f, None => return Vec::new() };
        // Match the buffer guard every other engine consumer uses (audit N4).
        if self.compensated.len() < fcs.n_events * fcs.n_params() { return Vec::new(); }
        let own = compute_own_masks(&self.gates, &self.compensated, &fcs.parameters, fcs.n_events);
        let by_id: HashMap<u32, &Gate> = self.gates.iter().map(|g| (g.id, g)).collect();
        effective_mask(pop, &by_id, &own, fcs.n_events)
    }

    /// Build the reference-overlay point cloud (a faded background sample), resolving
    /// channels by NAME in the reference and using the active sample's transforms.
    fn rebuild_ref_scatter(&mut self) {
        let ref_idx = match self.ref_sample {
            Some(i) if i < self.samples.len() && Some(i) != Some(self.active_sample) => i,
            _ => { self.ref_scatter = None; return; }
        };
        let n_params = match &self.fcs { Some(f) => f.n_params(), None => { self.ref_scatter = None; return; } };
        if n_params == 0 { self.ref_scatter = None; return; }
        let xi = self.x_ch.min(n_params - 1);
        let yi = self.y_ch.min(n_params - 1);
        let (x_name, y_name) = {
            let f = self.fcs.as_ref().unwrap();
            (f.parameters[xi].name.clone(), f.parameters[yi].name.clone())
        };
        let xt = self.cur_tf(xi).compile();
        let yt = self.cur_tf(yi).compile();

        let rf = match FcsFile::open(&self.samples[ref_idx].path) {
            Ok(f) => f, Err(_) => { self.ref_scatter = None; return; }
        };
        let (rxi, ryi) = match (rf.param_index(&x_name), rf.param_index(&y_name)) {
            (Some(a), Some(b)) => (a, b),
            _ => { self.ref_scatter = None; return; } // panel mismatch
        };
        let ev = match self.compensate_events(&rf) { Ok(e) => e, Err(_) => rf.events.clone() };
        let (ne, np) = (rf.n_events, rf.n_params());
        if ne == 0 || ev.len() < ne * np { self.ref_scatter = None; return; }
        let step = (ne / REF_OVERLAY_MAX.min(ne).max(1)).max(1);
        let points: Vec<[f64; 2]> = (0..ne).step_by(step)
            .map(|e| [xt.forward(ev[e * np + rxi]), yt.forward(ev[e * np + ryi])])
            .collect();
        self.ref_scatter = Some(RefScatter {
            ref_idx, x_ch: xi, y_ch: yi,
            x_label: self.cur_tf(xi).short_label().to_string(),
            y_label: self.cur_tf(yi).short_label().to_string(),
            points,
        });
    }

    fn regate(&mut self) {
        let fcs = match &self.fcs { Some(f) => f, None => return };
        let n = fcs.n_params();
        if self.compensated.len() < fcs.n_events * n { return; } // not yet processed

        let own = compute_own_masks(&self.gates, &self.compensated, &fcs.parameters, fcs.n_events);
        let by_id: HashMap<u32, &Gate> = self.gates.iter().map(|g| (g.id, g)).collect();

        let mut counts = HashMap::new();
        for g in &self.gates {
            let eff = effective_mask(g.id, &by_id, &own, fcs.n_events);
            let n_in = eff.iter().filter(|&&b| b).count();
            let n_parent = match g.parent {
                Some(pid) => effective_mask(pid, &by_id, &own, fcs.n_events).iter().filter(|&&b| b).count(),
                None => fcs.n_events,
            };
            counts.insert(g.id, (n_in, n_parent));
        }
        self.gate_counts = counts;
        self.pop_stats = None; // gates changed → population stats stale
        self.hist_cache = None;
        if self.active_pop.is_some() {
            self.scatter = None; // filtered view depends on masks
            self.data_gen = self.data_gen.wrapping_add(1); // grid cells depend on masks too
        }
    }

    // ── Undo / redo (gate-tree snapshots) ─────────────────────────────

    /// Snapshot the current gate tree before a mutation. Clears the redo stack.
    fn push_undo(&mut self) {
        let snap = self.gates.clone();
        self.push_undo_state(snap);
    }

    /// Push an explicit pre-mutation snapshot (used when the snapshot must be
    /// taken before an in-place widget edit). Clears the redo stack.
    fn push_undo_state(&mut self, snap: Vec<Gate>) {
        const MAX: usize = 100;
        self.undo_stack.push(snap);
        if self.undo_stack.len() > MAX { self.undo_stack.remove(0); }
        self.redo_stack.clear();
    }

    fn undo(&mut self) {
        if let Some(prev) = self.undo_stack.pop() {
            self.redo_stack.push(std::mem::replace(&mut self.gates, prev));
            self.after_gate_restore("Undo");
        } else {
            self.status = "Nothing to undo.".into();
        }
    }

    fn redo(&mut self) {
        if let Some(next) = self.redo_stack.pop() {
            self.undo_stack.push(std::mem::replace(&mut self.gates, next));
            self.after_gate_restore("Redo");
        } else {
            self.status = "Nothing to redo.".into();
        }
    }

    /// Re-derive ids/selections and mark caches stale after restoring a snapshot.
    fn after_gate_restore(&mut self, what: &str) {
        self.next_gate_id = self.gates.iter().map(|g| g.id).max().unwrap_or(0) + 1;
        if let Some(s) = self.selected_gate { if !self.gates.iter().any(|g| g.id == s) { self.selected_gate = None; } }
        if let Some(p) = self.new_gate_parent { if !self.gates.iter().any(|g| g.id == p) { self.new_gate_parent = None; } }
        if let Some(a) = self.active_pop { if !self.gates.iter().any(|g| g.id == a) { self.active_pop = None; } }
        if let Some(h) = self.hist_sample_pop { if !self.gates.iter().any(|g| g.id == h) { self.hist_sample_pop = None; } }
        if let Some(q) = self.qc_live_gate { if !self.gates.iter().any(|g| g.id == q) { self.qc_live_gate = None; } }
        self.hist_hidden.retain(|id| self.gates.iter().any(|g| g.id == *id));
        self.scatter = None;
        self.hist_cache = None;
        self.needs_regate = true;
        self.status = format!("{} ({} undo / {} redo)", what, self.undo_stack.len(), self.redo_stack.len());
    }

    fn commit_gate(&mut self, shape: GateShape) {
        self.commit_gate_on(self.x_ch, self.y_ch, shape);
    }

    /// Commit a 2-D gate on an explicit channel pair (used by grid cells, which
    /// each have their own X/Y independent of the single-plot axes).
    fn commit_gate_on(&mut self, xi: usize, yi: usize, shape: GateShape) {
        if self.fcs.is_none() { return; }
        self.push_undo();
        let fcs = match &self.fcs { Some(f) => f, None => return };
        let n = fcs.n_params();
        if n == 0 { return; }
        let xi = xi.min(n - 1);
        let yi = yi.min(n - 1);
        let id = self.next_gate_id;
        self.next_gate_id += 1;

        let gate = Gate {
            id,
            name: format!("Gate {}", id),
            parent: self.new_gate_parent,
            x_channel: fcs.parameters[xi].name.clone(),
            y_channel: fcs.parameters[yi].name.clone(),
            x_transform: self.cur_tf(xi),
            y_transform: self.cur_tf(yi),
            shape,
            quad_group: None,
        };
        self.gates.push(gate);
        self.needs_regate = true;
        self.status = format!("Added Gate {}", id);
    }

    /// Commit a quadrant gate: 4 rectangle populations split at (cx, cy) in display
    /// coords on the current X/Y channels — the natural tool for e.g. CD103×CD11b.
    fn commit_quadrant(&mut self, cx: f64, cy: f64) {
        self.commit_quadrant_on(self.x_ch, self.y_ch, cx, cy);
    }

    fn commit_quadrant_on(&mut self, xi: usize, yi: usize, cx: f64, cy: f64) {
        if self.fcs.is_none() { return; }
        self.push_undo();
        let fcs = match &self.fcs { Some(f) => f, None => return };
        let n = fcs.n_params();
        if n == 0 { return; }
        let xi = xi.min(n - 1);
        let yi = yi.min(n - 1);
        let xn = fcs.parameters[xi].name.clone();
        let yn = fcs.parameters[yi].name.clone();
        let (xt, yt) = (self.cur_tf(xi), self.cur_tf(yi));
        let xs = fcs.parameters[xi].label.clone().filter(|l| !l.is_empty()).unwrap_or_else(|| short_chan(&xn));
        let ys = fcs.parameters[yi].label.clone().filter(|l| !l.is_empty()).unwrap_or_else(|| short_chan(&yn));
        const BIG: f64 = 1.0e12;
        // Half-open split: the '-' side's upper edge is one ulp below the center so a
        // point exactly on a center line belongs to the '+' side only (no double-count).
        let cxm = cx.next_down();
        let cym = cy.next_down();
        // (name, x_range, y_range): UR, UL, LL, LR
        let quads = [
            (format!("{}+ {}+", xs, ys), (cx, BIG), (cy, BIG)),
            (format!("{}- {}+", xs, ys), (-BIG, cxm), (cy, BIG)),
            (format!("{}- {}-", xs, ys), (-BIG, cxm), (-BIG, cym)),
            (format!("{}+ {}-", xs, ys), (cx, BIG), (-BIG, cym)),
        ];
        let parent = self.new_gate_parent;
        let group = self.next_gate_id; // the first member's id doubles as the group id
        for (name, (x0, x1), (y0, y1)) in quads {
            let id = self.next_gate_id;
            self.next_gate_id += 1;
            self.gates.push(Gate {
                id, name, parent,
                x_channel: xn.clone(), y_channel: yn.clone(),
                x_transform: xt.clone(), y_transform: yt.clone(),
                shape: GateShape::Rect { x_min: x0, x_max: x1, y_min: y0, y_max: y1 },
                quad_group: Some(group),
            });
        }
        self.needs_regate = true;
        self.status = format!("Added linked quadrant on {}×{}", xs, ys);
    }

    /// Compensated (linear) values of channel `ci` over the active "gate-from-here"
    /// population (or all events). Used by the auto-gate suggestions.
    fn channel_for_pop(&self, ci: usize) -> Vec<f64> {
        let fcs = match &self.fcs { Some(f) => f, None => return Vec::new() };
        let n = fcs.n_params();
        if ci >= n || self.compensated.len() < fcs.n_events * n { return Vec::new(); }
        let mask = self.active_pop.map(|p| self.pop_mask(p));
        (0..fcs.n_events)
            .filter(|&e| mask.as_ref().map(|m| m.get(e).copied().unwrap_or(false)).unwrap_or(true))
            .map(|e| self.compensated[e * n + ci])
            .collect()
    }

    /// Auto-threshold the current X channel at the density valley between its two
    /// peaks (a *suggestion* — the user adjusts it). Computed in the channel's current
    /// display space, so it works on logicle/log fluorescence.
    fn auto_threshold(&mut self) {
        let (xi, name, np) = match &self.fcs {
            Some(f) => {
                let xi = self.x_ch.min(f.n_params().saturating_sub(1));
                (xi, f.parameters.get(xi).map(|p| p.name.clone()).unwrap_or_default(), f.n_params())
            }
            None => return,
        };
        if np == 0 { return; }
        let xt = self.cur_tf(xi);
        let ct = xt.compile();
        let vals: Vec<f64> = self.channel_for_pop(xi).iter().map(|&v| ct.forward(v)).collect();
        match crate::autogate::valley_threshold(&vals, 128) {
            Some((thr, depth)) => {
                let dmax = vals.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
                let xmax = if dmax.is_finite() && dmax > thr { dmax } else { thr };
                self.push_undo();
                let id = self.next_gate_id; self.next_gate_id += 1;
                self.gates.push(Gate {
                    id, name: format!("{}+ (auto)", short_chan(&name)),
                    parent: self.new_gate_parent,
                    x_channel: name.clone(), y_channel: name,
                    x_transform: xt.clone(), y_transform: xt,
                    shape: GateShape::Range { x_min: thr, x_max: xmax },
                    quad_group: None,
                });
                self.needs_regate = true;
                self.status = format!("Auto-threshold (valley depth {:.0}%) — review & adjust", depth * 100.0);
            }
            None => self.status = "No clear two populations on this channel — draw a gate manually.".into(),
        }
    }

    /// Suggest a singlet gate (FSC-A × FSC-H) as a diagonal band — a starting point.
    fn suggest_singlets(&mut self) {
        let found = self.fcs.as_ref().and_then(|fcs| {
            let pick = |f: &dyn Fn(&str) -> bool| fcs.parameters.iter().position(|p| f(&p.name));
            let ai = pick(&|n| n.eq_ignore_ascii_case("FSC-A"))
                .or_else(|| pick(&|n| { let u = n.to_uppercase(); u.starts_with("FSC") && u.ends_with("-A") }));
            let hi = pick(&|n| n.eq_ignore_ascii_case("FSC-H"))
                .or_else(|| pick(&|n| { let u = n.to_uppercase(); u.starts_with("FSC") && u.ends_with("-H") }));
            match (ai, hi) {
                (Some(a), Some(h)) => Some((a, h, fcs.parameters[a].name.clone(), fcs.parameters[h].name.clone())),
                _ => None,
            }
        });
        let (ai, hi, an, hn) = match found {
            Some(x) => x,
            None => { self.status = "Need FSC-A and FSC-H channels for a singlet gate.".into(); return; }
        };
        let area = self.channel_for_pop(ai);
        let height = self.channel_for_pop(hi);
        match crate::autogate::singlet_polygon(&area, &height, 0.15) {
            Some(vertices) => {
                self.push_undo();
                let id = self.next_gate_id; self.next_gate_id += 1;
                let lin = AxisTransform::Linear; // scatter: display coords == data coords
                self.gates.push(Gate {
                    id, name: "Singlets (auto)".into(),
                    parent: self.new_gate_parent,
                    x_channel: an, y_channel: hn,
                    x_transform: lin.clone(), y_transform: lin,
                    shape: GateShape::Polygon { vertices },
                    quad_group: None,
                });
                self.x_ch = ai; self.y_ch = hi; // show FSC-A × FSC-H so the gate is visible
                self.needs_rescatter = true;
                self.needs_regate = true;
                self.status = "Suggested singlet gate on FSC-A × FSC-H — review & adjust.".into();
            }
            None => self.status = "Couldn't fit a singlet gate from the scatter data.".into(),
        }
    }

    /// Commit a 1-D interval gate on the current X channel (drawn on a histogram).
    fn commit_range_gate(&mut self, x_min: f64, x_max: f64) {
        if self.fcs.is_none() { return; }
        self.push_undo();
        let fcs = match &self.fcs { Some(f) => f, None => return };
        let n = fcs.n_params();
        if n == 0 { return; }
        let xi = self.hist_ch.min(n - 1);
        let id = self.next_gate_id;
        self.next_gate_id += 1;
        let ch = fcs.parameters[xi].name.clone();
        let tf = self.cur_tf(xi);
        let gate = Gate {
            id,
            name: format!("Interval {}", id),
            parent: self.new_gate_parent,
            x_channel: ch.clone(),
            y_channel: ch,
            x_transform: tf.clone(),
            y_transform: tf,
            shape: GateShape::Range { x_min: x_min.min(x_max), x_max: x_min.max(x_max) },
            quad_group: None,
        };
        self.gates.push(gate);
        self.needs_regate = true;
        self.hist_cache = None;
        self.status = format!("Added interval gate {}", id);
    }

    /// Create a Boolean (AND/OR/NOT) gate combining the selected populations.
    fn commit_boolean(&mut self, op: BoolOp, refs: Vec<u32>) {
        if refs.is_empty() { self.status = "Select at least one gate for the boolean.".into(); return; }
        if op == BoolOp::Not && refs.len() != 1 { self.status = "NOT combines exactly one gate.".into(); return; }
        self.push_undo();
        let ref_names: Vec<String> = refs.iter()
            .filter_map(|r| self.gates.iter().find(|g| g.id == *r).map(|g| g.name.clone())).collect();
        let name = match op {
            BoolOp::And => ref_names.join(" & "),
            BoolOp::Or => ref_names.join(" | "),
            BoolOp::Not => format!("NOT {}", ref_names.first().cloned().unwrap_or_default()),
        };
        let id = self.next_gate_id;
        self.next_gate_id += 1;
        self.gates.push(Gate {
            id, name, parent: self.new_gate_parent,
            x_channel: String::new(), y_channel: String::new(),
            x_transform: AxisTransform::Linear, y_transform: AxisTransform::Linear,
            shape: GateShape::Boolean { op, refs },
            quad_group: None,
        });
        self.needs_regate = true;
        self.status = "Added boolean gate".into();
    }

    /// Finish the in-progress polygon, committing it to the right channel pair —
    /// the active grid cell's channels in grid mode, else the single-plot axes.
    fn finish_polygon(&mut self) {
        let verts = std::mem::take(&mut self.poly_vertices);
        if verts.len() >= 3 {
            if self.grid_mode {
                if let Some(c) = self.active_grid_cell {
                    let (xi, yi) = self.grid_channels.get(c).copied().unwrap_or((self.x_ch, self.y_ch));
                    self.commit_gate_on(xi, yi, GateShape::Polygon { vertices: verts });
                }
            } else {
                self.commit_gate(GateShape::Polygon { vertices: verts });
            }
        }
        self.draw_mode = DrawMode::Navigate;
        self.active_grid_cell = None;
    }

    /// Build per-population binned histogram of the current X channel (display space).
    fn rebuild_histogram(&mut self) {
        let (n_events, n_params) = match &self.fcs {
            Some(f) => (f.n_events, f.n_params()), None => return,
        };
        if n_events == 0 || n_params == 0 || self.compensated.len() < n_events * n_params {
            return;
        }
        let xi = self.hist_ch.min(n_params - 1);
        let xt = self.cur_tf(xi).compile();
        let norm = self.hist_norm;
        let x_name = self.fcs.as_ref().map(|f| f.parameters[xi].name.clone()).unwrap_or_default();
        const B: usize = 256;

        let (centers, series) = if self.hist_mode == HistMode::Samples && self.samples.len() > 1 {
            // ── Overlay SAMPLES: stream each file once, collecting the chosen
            //    population's X display-values, then bin into a COMMON range so the
            //    curves are aligned (binning each by its own range would misalign them).
            let pop = self.hist_sample_pop;
            let mut per_sample: Vec<(String, Color32, Vec<f64>)> = Vec::new();
            let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
            for (si, s) in self.samples.iter().enumerate() {
                let fcs = match FcsFile::open(&s.path) { Ok(f) => f, Err(_) => continue };
                let cxi = match fcs.param_index(&x_name) { Some(i) => i, None => continue };
                let (ne, np) = (fcs.n_events, fcs.n_params());
                let ev = match self.compensate_events(&fcs) { Ok(e) => e, Err(_) => continue };
                if ev.len() < ne * np { continue; }

                let mask: Vec<bool> = if let Some(gid) = pop {
                    if !missing_gate_channels(&self.gates, &fcs).is_empty() { continue; }
                    let own = compute_own_masks(&self.gates, &ev, &fcs.parameters, ne);
                    let by_id: HashMap<u32, &Gate> = self.gates.iter().map(|g| (g.id, g)).collect();
                    effective_mask(gid, &by_id, &own, ne)
                } else {
                    vec![true; ne]
                };

                let vals: Vec<f64> = (0..ne).filter(|&e| mask[e])
                    .map(|e| xt.forward(ev[e * np + cxi])).collect();
                for &v in &vals {
                    if v.is_finite() { lo = lo.min(v); hi = hi.max(v); }
                }
                per_sample.push((s.name.clone(), sample_color(si), vals));
                // fcs + ev drop here → only one channel's values retained per sample
            }
            if !lo.is_finite() || !hi.is_finite() || lo >= hi { lo = 0.0; hi = 1.0; }
            let span = (hi - lo).max(1e-9);
            let centers: Vec<f64> = (0..B).map(|b| lo + (b as f64 + 0.5) * span / B as f64).collect();
            let binof = |x: f64| (((x - lo) / span * B as f64) as isize).clamp(0, B as isize - 1) as usize;
            let series = per_sample.into_iter().map(|(name, color, vals)| {
                let mut c = vec![0.0; B];
                for v in vals { c[binof(v)] += 1.0; }
                HistSeries { name, color, values: normalize_hist(c, norm) }
            }).collect();
            (centers, series)
        } else {
            // ── Overlay POPULATIONS within the active sample ──
            let dx: Vec<f64> = (0..n_events)
                .map(|e| xt.forward(self.compensated[e * n_params + xi])).collect();
            let (lo, hi) = data_range(&dx);
            let span = (hi - lo).max(1e-9);
            let centers: Vec<f64> = (0..B).map(|b| lo + (b as f64 + 0.5) * span / B as f64).collect();
            let binof = |x: f64| (((x - lo) / span * B as f64) as isize).clamp(0, B as isize - 1) as usize;

            let own = match &self.fcs {
                Some(fcs) => compute_own_masks(&self.gates, &self.compensated, &fcs.parameters, n_events),
                None => HashMap::new(),
            };
            let by_id: HashMap<u32, &Gate> = self.gates.iter().map(|g| (g.id, g)).collect();
            let mut series = Vec::new();
            if !self.hist_all_hidden {
                let mut c = vec![0.0; B];
                for &x in &dx { c[binof(x)] += 1.0; }
                series.push(HistSeries { name: "All events".into(), color: Color32::GRAY, values: normalize_hist(c, norm) });
            }
            for g in self.gates.iter() {
                if self.hist_hidden.contains(&g.id) { continue; }
                let eff = effective_mask(g.id, &by_id, &own, n_events);
                let mut c = vec![0.0; B];
                for e in 0..n_events {
                    if eff[e] { c[binof(dx[e])] += 1.0; }
                }
                series.push(HistSeries { name: g.name.clone(), color: gate_color(g.id).0, values: normalize_hist(c, norm) });
            }
            (centers, series)
        };

        let x_label = self.cur_tf(xi).short_label().to_string();
        self.hist_cache = Some(HistogramData {
            x_ch: xi, x_label, norm, mode: self.hist_mode, sample_pop: self.hist_sample_pop, centers, series,
        });
    }
}

// ── eframe::App ───────────────────────────────────────────────────────────

impl eframe::App for FlowCytoApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        ctx.set_visuals(themed_visuals(self.dark_mode));

        // macOS Dock icon: assert it over the first few frames. Setting it in the
        // run_native closure does NOT stick — eframe/winit finishes launching after
        // the closure and macOS installs the default icon during that, clobbering an
        // early set. Re-applying for the first frames (after launch completes) wins.
        #[cfg(target_os = "macos")]
        if self.dock_icon_ticks < 8 {
            set_dock_icon();
            self.dock_icon_ticks += 1;
        }

        if self.needs_reprocess {
            self.reprocess();
            self.needs_reprocess = false;
            self.needs_rescatter = true;
            self.needs_regate = true;
        }
        if self.needs_regate {
            self.regate();
            self.needs_regate = false;
        }
        if self.needs_rescatter {
            self.rebuild_scatter();
            self.needs_rescatter = false;
        }

        if self.draw_mode != DrawMode::Navigate {
            ctx.set_cursor_icon(egui::CursorIcon::Crosshair);
        }

        #[cfg(target_os = "macos")]
        self.handle_menu_events(ctx);
        self.handle_dropped_files(ctx);
        self.poll_batch(ctx);
        self.poll_qc(ctx);
        self.poll_updates(ctx);
        self.handle_keys(ctx);

        self.panel_top(ctx);
        self.panel_left(ctx);
        self.panel_status(ctx);
        self.panel_central(ctx);

        self.poll_screenshot(ctx);
    }
}

// ── UI panels ─────────────────────────────────────────────────────────────

impl FlowCytoApp {
    /// Open .fcs files dropped onto the window; show a hint while hovering.
    fn handle_dropped_files(&mut self, ctx: &egui::Context) {
        let hovering = ctx.input(|i| !i.raw.hovered_files.is_empty());
        if hovering {
            egui::Area::new(egui::Id::new("drop_hint")).fixed_pos(egui::pos2(20.0, 40.0)).show(ctx, |ui| {
                egui::Frame::popup(ui.style()).show(ui, |ui| {
                    ui.label(RichText::new(format!("{} Drop .fcs files to open", icon::ARROW_LINE_DOWN)).size(16.0));
                });
            });
        }
        let dropped: Vec<PathBuf> = ctx.input(|i| {
            i.raw.dropped_files.iter()
                .filter_map(|f| f.path.clone())
                .filter(|p| p.extension().map(|e| e.eq_ignore_ascii_case("fcs")).unwrap_or(false))
                .collect()
        });
        if !dropped.is_empty() { self.add_files(dropped); }
    }

    /// Keyboard shortcuts. Skipped while a text field has focus so typing names
    /// or numbers is never hijacked.
    fn handle_keys(&mut self, ctx: &egui::Context) {
        if ctx.wants_keyboard_input() { return; }
        use egui::Key;
        struct Keys {
            undo: bool, redo: bool, save: bool, esc: bool,
            rect: bool, ellipse: bool, poly: bool, quad: bool, edit: bool, nav: bool,
            tabs: [bool; 6],
        }
        let k = ctx.input(|i| {
            let cmd = i.modifiers.command;
            let plain = !cmd && !i.modifiers.shift && !i.modifiers.alt;
            Keys {
                undo: cmd && !i.modifiers.shift && i.key_pressed(Key::Z),
                redo: cmd && i.modifiers.shift && i.key_pressed(Key::Z),
                save: cmd && i.key_pressed(Key::S),
                esc: i.key_pressed(Key::Escape),
                rect: plain && i.key_pressed(Key::R),
                ellipse: plain && i.key_pressed(Key::E),
                poly: plain && i.key_pressed(Key::P),
                quad: plain && i.key_pressed(Key::Q),
                edit: plain && i.key_pressed(Key::G),
                nav: plain && i.key_pressed(Key::V),
                tabs: [
                    plain && i.key_pressed(Key::Num1), plain && i.key_pressed(Key::Num2),
                    plain && i.key_pressed(Key::Num3), plain && i.key_pressed(Key::Num4),
                    plain && i.key_pressed(Key::Num5), plain && i.key_pressed(Key::Num6),
                ],
            }
        });

        if k.save { self.save_gates(); }
        if k.undo { self.undo(); }
        if k.redo { self.redo(); }

        let set_mode = |app: &mut Self, m: DrawMode| {
            app.draw_mode = if app.draw_mode == m { DrawMode::Navigate } else { m };
            app.drag_start = None; app.drag_current = None; app.poly_vertices.clear();
        };
        if k.rect { set_mode(self, DrawMode::Rect); }
        if k.ellipse { set_mode(self, DrawMode::Ellipse); }
        if k.poly { set_mode(self, DrawMode::Polygon); }
        if k.quad { set_mode(self, DrawMode::Quadrant); }
        if k.edit { set_mode(self, DrawMode::Edit); }
        if k.nav || k.esc {
            self.draw_mode = DrawMode::Navigate;
            self.drag_start = None; self.drag_current = None; self.poly_vertices.clear();
        }

        const TAB_ORDER: [ActiveTab; 6] = [
            ActiveTab::Plot, ActiveTab::Histogram, ActiveTab::Stats, ActiveTab::Batch, ActiveTab::Spillover, ActiveTab::QC,
        ];
        for (t, &on) in TAB_ORDER.iter().zip(k.tabs.iter()) {
            if on { self.active_tab = *t; }
        }
    }

    /// Dispatch native-menu clicks. Actions mirror the toolbar/keyboard exactly,
    /// so the menu never diverges from the in-app controls.
    #[cfg(target_os = "macos")]
    fn handle_menu_events(&mut self, ctx: &egui::Context) {
        enum A { Open, SaveGates, SaveSession, LoadSession, SavePlot, CheckUpdates, Undo, Redo, Theme, Tab(usize), ZoomIn, ZoomOut, ZoomReset }
        let mut acts: Vec<A> = Vec::new();
        if let Some(st) = &self.menu {
            while let Ok(ev) = muda::MenuEvent::receiver().try_recv() {
                let id = ev.id;
                if id == st.open { acts.push(A::Open); }
                else if id == st.save_gates { acts.push(A::SaveGates); }
                else if id == st.save_session { acts.push(A::SaveSession); }
                else if id == st.load_session { acts.push(A::LoadSession); }
                else if id == st.save_plot { acts.push(A::SavePlot); }
                else if id == st.check_updates { acts.push(A::CheckUpdates); }
                else if id == st.undo { acts.push(A::Undo); }
                else if id == st.redo { acts.push(A::Redo); }
                else if id == st.theme { acts.push(A::Theme); }
                else if id == st.zoom_in { acts.push(A::ZoomIn); }
                else if id == st.zoom_out { acts.push(A::ZoomOut); }
                else if id == st.zoom_reset { acts.push(A::ZoomReset); }
                else if let Some(i) = st.tabs.iter().position(|t| *t == id) { acts.push(A::Tab(i)); }
            }
        }
        const TABS: [ActiveTab; 6] = [
            ActiveTab::Plot, ActiveTab::Histogram, ActiveTab::Stats, ActiveTab::Batch, ActiveTab::Spillover, ActiveTab::QC,
        ];
        for a in acts {
            match a {
                A::Open => if let Some(paths) = rfd::FileDialog::new()
                    .add_filter("FCS files", &["fcs", "FCS"]).pick_files() { self.add_files(paths); },
                A::SaveGates => self.save_gates(),
                A::SaveSession => self.save_session(),
                A::LoadSession => self.load_session(),
                A::SavePlot => self.request_plot_png(),
                A::CheckUpdates => self.start_update_check(),
                A::Undo => self.undo(),
                A::Redo => self.redo(),
                A::Theme => { self.dark_mode = !self.dark_mode; self.needs_rescatter = true; }
                A::ZoomIn => ctx.set_zoom_factor((ctx.zoom_factor() * 1.1).min(3.0)),
                A::ZoomOut => ctx.set_zoom_factor((ctx.zoom_factor() / 1.1).max(0.5)),
                A::ZoomReset => ctx.set_zoom_factor(1.0),
                A::Tab(i) => self.active_tab = TABS[i],
            }
        }
    }

    fn panel_top(&mut self, ctx: &egui::Context) {
        egui::TopBottomPanel::top("toolbar").show(ctx, |ui| {
            ui.horizontal(|ui| {
                if ui.button(format!("{}  Open FCS", icon::FILE_PLUS)).clicked() {
                    if let Some(paths) = rfd::FileDialog::new()
                        .add_filter("FCS files", &["fcs", "FCS"]).pick_files()
                    {
                        self.add_files(paths);
                    }
                }
                // Recent files
                {
                    let recent = self.recent_files.clone();
                    let mut pick: Option<PathBuf> = None;
                    let mut clear = false;
                    ui.menu_button(format!("{} Recent", icon::CLOCK_COUNTER_CLOCKWISE), |ui| {
                        if recent.is_empty() { ui.label(RichText::new("(none)").weak()); }
                        for p in &recent {
                            let name = p.file_name().map(|s| s.to_string_lossy().to_string()).unwrap_or_default();
                            if ui.button(name).on_hover_text(p.to_string_lossy()).clicked() {
                                pick = Some(p.clone()); ui.close_menu();
                            }
                        }
                        if !recent.is_empty() {
                            ui.separator();
                            if ui.button("Clear recent").clicked() { clear = true; ui.close_menu(); }
                        }
                    });
                    if clear { self.recent_files.clear(); self.save_recent(); }
                    if let Some(p) = pick { self.add_files(vec![p]); }
                }
                if ui.button(format!("{} Save session", icon::FLOPPY_DISK)).on_hover_text("Save workspace + gates + transforms + compensation").clicked() {
                    self.save_session();
                }
                if ui.button(format!("{} Load session", icon::FOLDER_OPEN)).clicked() { self.load_session(); }
                if ui.button(format!("{} Report", icon::FILE_TEXT)).on_hover_text("Export a self-contained HTML reproducibility report (provenance + gates + stats)").clicked() { self.export_report(); }
                ui.separator();
                if ui.checkbox(&mut self.do_compensate, "Compensate").changed() {
                    self.needs_reprocess = true;
                }
                ui.separator();
                // Theme toggle
                let theme_label = if self.dark_mode { format!("{} Light", icon::SUN) } else { format!("{} Dark", icon::MOON) };
                if ui.button(theme_label).clicked() {
                    self.dark_mode = !self.dark_mode;
                    self.needs_rescatter = true; // recolor density floor
                }
                ui.separator();
                ui.selectable_value(&mut self.active_tab, ActiveTab::Plot, format!("{} Plot", icon::CHART_SCATTER));
                ui.selectable_value(&mut self.active_tab, ActiveTab::Histogram, format!("{} Histogram", icon::CHART_BAR));
                ui.selectable_value(&mut self.active_tab, ActiveTab::Stats, format!("{} Stats", icon::TABLE));
                ui.selectable_value(&mut self.active_tab, ActiveTab::Batch, format!("{} Batch", icon::STACK));
                ui.selectable_value(&mut self.active_tab, ActiveTab::Spillover, format!("{} Spillover", icon::GRID_NINE));
                ui.selectable_value(&mut self.active_tab, ActiveTab::QC, format!("{} QC", icon::HEARTBEAT));

                // Right-aligned: manual update check (the app's only network call).
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    self.ui_update_check(ui);
                });
            });
        });
    }

    /// Toolbar "Check for Updates" control. Rendered inside a right-to-left layout, so
    /// the first widget added is the right-most. Network happens only on click.
    fn ui_update_check(&mut self, ui: &mut egui::Ui) {
        if self.update_rx.is_some() {
            ui.add(egui::Spinner::new());
            ui.label(RichText::new("checking…").weak());
            return;
        }
        if ui.button(format!("{} Updates", icon::ARROWS_CLOCKWISE))
            .on_hover_text("Check GitHub for a newer version — the app's only network call, made only when you click this")
            .clicked()
        {
            self.start_update_check();
        }
        if let Some(info) = self.update_found.clone() {
            if ui.button(RichText::new(format!("{} Download v{}", icon::DOWNLOAD_SIMPLE, info.latest))
                    .strong().color(Color32::from_rgb(90, 200, 90)))
                .on_hover_text(format!("Open {}", info.url))
                .clicked()
            {
                let _ = open::that(&info.url);
            }
        } else if let Some(msg) = &self.update_msg {
            ui.label(RichText::new(msg).weak());
        }
    }

    /// Spawn the (blocking) GitHub Releases check on a background thread.
    fn start_update_check(&mut self) {
        let (tx, rx) = std::sync::mpsc::channel();
        std::thread::spawn(move || { let _ = tx.send(crate::update::check_latest()); });
        self.update_rx = Some(rx);
        self.update_msg = Some("Checking for updates…".into());
        self.update_found = None;
    }

    /// Drain the update-check result (called each frame while a check is in flight).
    fn poll_updates(&mut self, ctx: &egui::Context) {
        let done = match &self.update_rx {
            Some(rx) => match rx.try_recv() {
                Ok(result) => Some(result),
                Err(TryRecvError::Empty) => { ctx.request_repaint(); None }
                Err(TryRecvError::Disconnected) => Some(Err("update check stopped unexpectedly".to_string())),
            },
            None => return,
        };
        if let Some(result) = done {
            self.update_rx = None;
            match result {
                Ok(info) => {
                    self.update_msg = Some(if info.newer {
                        format!("v{} available (you have {})", info.latest, info.current)
                    } else {
                        format!("Up to date (v{})", info.current)
                    });
                    self.update_found = if info.newer { Some(info) } else { None };
                }
                Err(e) => {
                    self.update_msg = Some(format!("Update check failed: {e}"));
                    self.update_found = None;
                }
            }
        }
    }

    fn panel_left(&mut self, ctx: &egui::Context) {
        egui::SidePanel::left("left").resizable(true).default_width(260.0).show(ctx, |ui| {
            egui::ScrollArea::vertical().show(ui, |ui| {
                self.ui_samples(ui);
                self.ui_channels(ui);
                ui.separator();
                self.ui_axis_limits(ui);
                ui.separator();
                self.ui_gates(ui);
            });
        });
    }

    fn ui_channels(&mut self, ui: &mut egui::Ui) {
        ui.heading("Channels");
        // Extract owned values so no immutable borrow of `self` is held across
        // the `&mut self` transform-combo calls below.
        let (names, n_events, has_spill) = match &self.fcs {
            Some(f) => (
                f.parameters.iter().map(|p| match &p.label {
                    Some(l) if !l.is_empty() => format!("{} ({})", p.name, l),
                    _ => p.name.clone(),
                }).collect::<Vec<String>>(),
                f.n_events,
                f.spillover_keyword().is_some(),
            ),
            None => { ui.label(RichText::new("No file loaded.").color(Color32::GRAY)); return; }
        };
        let n = names.len();

        let mut ch_changed = false;
        let mut tf_changed = false;

        // X axis
        ui.horizontal(|ui| {
            ui.label("X:");
            egui::ComboBox::from_id_salt("xch")
                .selected_text(names.get(self.x_ch).cloned().unwrap_or_default())
                .width(150.0)
                .show_ui(ui, |ui| {
                    for (i, nm) in names.iter().enumerate() {
                        ch_changed |= ui.selectable_value(&mut self.x_ch, i, nm).changed();
                    }
                });
        });
        tf_changed |= self.ui_transform_combo(ui, "xtf", self.x_ch);

        // Y axis
        ui.horizontal(|ui| {
            ui.label("Y:");
            egui::ComboBox::from_id_salt("ych")
                .selected_text(names.get(self.y_ch).cloned().unwrap_or_default())
                .width(150.0)
                .show_ui(ui, |ui| {
                    for (i, nm) in names.iter().enumerate() {
                        ch_changed |= ui.selectable_value(&mut self.y_ch, i, nm).changed();
                    }
                });
        });
        tf_changed |= self.ui_transform_combo(ui, "ytf", self.y_ch);

        // Apply the X channel's transform to every fluorescence channel at once —
        // saves setting Logicle/asinh one-by-one across a big panel.
        if ui.add(egui::Button::new(RichText::new(format!("{} X scale → all fluorescence", icon::ARROW_LINE_DOWN)).small()))
            .on_hover_text("Set every fluorescence channel's scale to match the current X channel")
            .clicked()
        {
            let tf = self.cur_tf(self.x_ch);
            if let Some(fcs) = &self.fcs {
                for i in crate::transform::fluorescence_indices(&fcs.parameters) {
                    if i < self.channel_tf.len() { self.channel_tf[i] = tf.clone(); }
                }
            }
            tf_changed = true;
        }

        ui.add_space(4.0);
        if let Some(p) = &self.file_path {
            ui.label(RichText::new(p.file_name().unwrap_or_default().to_string_lossy()).small());
        }
        ui.label(format!("{} events · {} ch", n_events, n));
        if has_spill {
            ui.label(RichText::new(format!("$SPILLOVER {}", icon::CHECK)).small().color(Color32::from_rgb(80, 180, 80)));
        }

        if ch_changed || tf_changed { self.needs_rescatter = true; }
    }

    /// Transform picker for one channel. Returns true if changed.
    fn ui_transform_combo(&mut self, ui: &mut egui::Ui, id: &str, ch: usize) -> bool {
        let mut changed = false;
        let cur = self.cur_tf(ch);
        ui.horizontal(|ui| {
            ui.add_space(18.0);
            ui.label(RichText::new("scale:").small());
            egui::ComboBox::from_id_salt(id)
                .selected_text(cur.short_label())
                .width(110.0)
                .show_ui(ui, |ui| {
                    let opts = [
                        ("Linear", AxisTransform::Linear),
                        ("Log", AxisTransform::default_log()),
                        ("Asinh", AxisTransform::Asinh { cofactor: 150.0 }),
                        ("Logicle", AxisTransform::default_logicle()),
                    ];
                    for (lbl, tf) in opts {
                        let selected = cur.short_label() == lbl;
                        if ui.selectable_label(selected, lbl).clicked() && !selected {
                            if ch < self.channel_tf.len() { self.channel_tf[ch] = tf; }
                            changed = true;
                        }
                    }
                });
            // Cofactor editor for asinh
            if let AxisTransform::Asinh { cofactor } = self.cur_tf(ch) {
                let mut cf = cofactor;
                if ui.add(egui::DragValue::new(&mut cf).speed(5.0).range(1.0..=100000.0)).changed() {
                    if ch < self.channel_tf.len() { self.channel_tf[ch] = AxisTransform::Asinh { cofactor: cf }; }
                    changed = true;
                }
            }
        });
        changed
    }

    /// Compact per-axis scale picker for a grid cell. Edits the channel's shared
    /// transform (same one the left panel uses), so it updates everywhere.
    fn grid_axis_scale(&mut self, ui: &mut egui::Ui, id: &str, ch: usize) {
        let cur = self.cur_tf(ch);
        egui::ComboBox::from_id_salt(id).selected_text(cur.short_label()).width(74.0).show_ui(ui, |ui| {
            let opts = [
                ("Linear", AxisTransform::Linear),
                ("Log", AxisTransform::default_log()),
                ("Asinh", AxisTransform::Asinh { cofactor: 150.0 }),
                ("Logicle", AxisTransform::default_logicle()),
            ];
            for (lbl, tf) in opts {
                let selected = cur.short_label() == lbl;
                if ui.selectable_label(selected, lbl).clicked() && !selected && ch < self.channel_tf.len() {
                    self.channel_tf[ch] = tf;
                }
            }
        });
    }

    fn ui_axis_limits(&mut self, ui: &mut egui::Ui) {
        ui.heading("Axis limits");
        ui.label(RichText::new("(data units; off = auto-fit)").small().color(Color32::GRAY));
        let mut changed = false;
        ui.horizontal(|ui| {
            changed |= ui.checkbox(&mut self.x_manual, "X").changed();
            ui.add_enabled(self.x_manual, egui::DragValue::new(&mut self.x_lo).prefix("min ").speed(100.0));
            ui.add_enabled(self.x_manual, egui::DragValue::new(&mut self.x_hi).prefix("max ").speed(100.0));
        });
        ui.horizontal(|ui| {
            changed |= ui.checkbox(&mut self.y_manual, "Y").changed();
            ui.add_enabled(self.y_manual, egui::DragValue::new(&mut self.y_lo).prefix("min ").speed(100.0));
            ui.add_enabled(self.y_manual, egui::DragValue::new(&mut self.y_hi).prefix("max ").speed(100.0));
        });
        let _ = changed; // bounds applied each frame in scatter_plot
    }

    fn ui_gates(&mut self, ui: &mut egui::Ui) {
        ui.heading("Gates");
        // "Gate from here" breadcrumb: which population the plot is restricted to.
        {
            let path = self.population_path(self.active_pop);
            ui.horizontal_wrapped(|ui| {
                ui.label(RichText::new("Viewing:").small());
                if ui.small_button("All events").clicked() {
                    self.active_pop = None; self.scatter = None;
                }
                for (id, nm) in &path {
                    ui.label(RichText::new("›").small().color(Color32::GRAY));
                    if ui.small_button(nm).clicked() {
                        self.active_pop = Some(*id); self.new_gate_parent = Some(*id); self.scatter = None;
                    }
                }
                if self.active_pop.is_some() {
                    ui.separator();
                    if ui.checkbox(&mut self.backgate, "backgate").changed() { self.scatter = None; }
                }
            });
        }
        if self.fcs.is_none() {
            ui.label(RichText::new("Load a file to gate.").color(Color32::GRAY));
            return;
        }

        // Draw-mode buttons
        ui.horizontal(|ui| {
            self.draw_btn(ui, DrawMode::Rect, &format!("{} Rect", icon::RECTANGLE));
            self.draw_btn(ui, DrawMode::Ellipse, &format!("{} Ellipse", icon::CIRCLE));
            self.draw_btn(ui, DrawMode::Polygon, &format!("{} Polygon", icon::POLYGON));
            self.draw_btn(ui, DrawMode::Quadrant, &format!("{} Quad", icon::CROSSHAIR_SIMPLE));
            self.draw_btn(ui, DrawMode::Edit, &format!("{} Edit", icon::PENCIL_SIMPLE));
        });
        ui.horizontal(|ui| {
            if ui.button(format!("{} Suggest singlets", icon::MAGIC_WAND))
                .on_hover_text("Fit a diagonal FSC-A × FSC-H singlet band (a starting suggestion to review)")
                .clicked() { self.suggest_singlets(); }
            if ui.button(format!("{} Auto-threshold X", icon::MAGIC_WAND))
                .on_hover_text("Split the current X channel at its density valley (a suggestion to review)")
                .clicked() { self.auto_threshold(); }
        });
        ui.label(RichText::new("keys: R/E/P/Q draw · G edit · V/Esc nav · ⌘Z undo")
            .small().color(Color32::GRAY))
            .on_hover_text("1–5 switch tabs · ⌘S save gates · ⌘⇧Z redo");
        if self.draw_mode != DrawMode::Navigate {
            ui.horizontal(|ui| {
                let hint = match self.draw_mode {
                    DrawMode::Polygon => "click to add points · close: click the green dot, double-click, or Finish",
                    DrawMode::Quadrant => "click to set the quadrant center",
                    DrawMode::Edit => "select a gate (■): drag handles to resize, body to move, outer handle to rotate",
                    _ => "drag on plot",
                };
                ui.label(RichText::new(format!("{} {}", icon::PENCIL_SIMPLE, hint)).color(Color32::from_rgb(220, 170, 0)).small());
                if self.draw_mode == DrawMode::Polygon && self.poly_vertices.len() >= 3
                    && ui.small_button("Finish").clicked()
                {
                    self.finish_polygon();
                }
                if ui.small_button("Cancel").clicked() {
                    self.draw_mode = DrawMode::Navigate;
                    self.drag_start = None; self.drag_current = None; self.poly_vertices.clear();
                }
            });
        }

        // Parent selector for next gate
        let gate_names: Vec<(u32, String)> = self.gates.iter().map(|g| (g.id, g.name.clone())).collect();
        ui.horizontal(|ui| {
            ui.label(RichText::new("New gate parent:").small());
            let sel = self.new_gate_parent
                .and_then(|pid| gate_names.iter().find(|(id, _)| *id == pid).map(|(_, n)| n.clone()))
                .unwrap_or_else(|| "All events".into());
            egui::ComboBox::from_id_salt("parentsel").selected_text(sel).width(140.0).show_ui(ui, |ui| {
                ui.selectable_value(&mut self.new_gate_parent, None, "All events");
                for (id, nm) in &gate_names {
                    ui.selectable_value(&mut self.new_gate_parent, Some(*id), nm);
                }
            });
        });

        // Undo / redo
        ui.horizontal(|ui| {
            if ui.add_enabled(!self.undo_stack.is_empty(), egui::Button::new(format!("{} Undo", icon::ARROW_ARC_LEFT)))
                .on_hover_text("Undo gate change (Ctrl/Cmd+Z)").clicked() { self.undo(); }
            if ui.add_enabled(!self.redo_stack.is_empty(), egui::Button::new(format!("{} Redo", icon::ARROW_ARC_RIGHT)))
                .on_hover_text("Redo (Ctrl/Cmd+Shift+Z)").clicked() { self.redo(); }
        });

        // Save / load
        ui.horizontal(|ui| {
            if ui.button(format!("{} Save", icon::FLOPPY_DISK)).clicked() { self.save_gates(); }
            if ui.button(format!("{} Load", icon::FOLDER_OPEN)).clicked() { self.load_gates(); }
        });

        // Boolean (AND/OR/NOT) gate builder.
        if !self.gates.is_empty() {
            let gate_list: Vec<(u32, String)> = self.gates.iter().map(|g| (g.id, g.name.clone())).collect();
            let mut create_bool = false;
            egui::CollapsingHeader::new(format!("{} Boolean gate", icon::PLUS)).id_salt("bool_builder").show(ui, |ui| {
                ui.horizontal(|ui| {
                    ui.label("Op:");
                    ui.selectable_value(&mut self.bool_op, BoolOp::And, "AND");
                    ui.selectable_value(&mut self.bool_op, BoolOp::Or, "OR");
                    ui.selectable_value(&mut self.bool_op, BoolOp::Not, "NOT");
                });
                for (id, name) in &gate_list {
                    let mut on = self.bool_refs.contains(id);
                    if ui.checkbox(&mut on, name).changed() {
                        if on { self.bool_refs.insert(*id); } else { self.bool_refs.remove(id); }
                    }
                }
                let hint = if self.bool_op == BoolOp::Not { "(pick exactly 1)" } else { "(pick 2+)" };
                ui.horizontal(|ui| {
                    if ui.button("Create").clicked() { create_bool = true; }
                    ui.label(RichText::new(hint).small().color(Color32::GRAY));
                });
            });
            if create_bool {
                let refs: Vec<u32> = gate_list.iter().map(|(id, _)| *id).filter(|id| self.bool_refs.contains(id)).collect();
                self.commit_boolean(self.bool_op, refs);
                self.bool_refs.clear();
            }
        }

        ui.separator();

        // Hierarchical gate list (depth-first).
        let order = crate::gating::gate_tree_order(&self.gates);
        let mut to_delete: Option<u32> = None;
        let mut reparent: Option<(u32, Option<u32>)> = None;
        let mut zoom_to: Option<u32> = None;
        let mut export_gid: Option<u32> = None;
        for (gid, depth) in order {
            let idx = match self.gates.iter().position(|g| g.id == gid) { Some(i) => i, None => continue };
            let (color, _) = gate_color(gid);
            let (n_in, n_parent) = self.gate_counts.get(&gid).copied().unwrap_or((0, 0));
            let pct_par = if n_parent > 0 { 100.0 * n_in as f64 / n_parent as f64 } else { 0.0 };
            let (name, shape_lbl, xch, ych, is_bool) = {
                let g = &self.gates[idx];
                (g.name.clone(), shape_label(&g.shape), g.x_channel.clone(), g.y_channel.clone(),
                 matches!(g.shape, GateShape::Boolean { .. }))
            };

            ui.horizontal(|ui| {
                ui.add_space(depth as f32 * 14.0);
                let sel = self.selected_gate == Some(gid);
                if ui.add(egui::SelectableLabel::new(sel, RichText::new("■").color(color))).clicked() {
                    self.selected_gate = if sel { None } else { Some(gid) };
                }
                let mut nm = name.clone();
                if ui.add(egui::TextEdit::singleline(&mut nm).desired_width(96.0)).changed() {
                    self.gates[idx].name = nm;
                }
                if ui.small_button("▶").on_hover_text("gate from here (plot only this population)").clicked() {
                    self.active_pop = Some(gid); self.new_gate_parent = Some(gid); self.scatter = None;
                }
                let hidden = self.scatter_hidden.contains(&gid);
                if ui.small_button(if hidden { icon::EYE_SLASH } else { icon::EYE })
                    .on_hover_text("show / hide this gate's outline on the plot").clicked()
                {
                    if hidden { self.scatter_hidden.remove(&gid); } else { self.scatter_hidden.insert(gid); }
                }
                if ui.small_button(icon::MAGNIFYING_GLASS_PLUS).on_hover_text("zoom the plot to this gate").clicked() { zoom_to = Some(gid); }
                if ui.small_button(icon::TRASH).clicked() { to_delete = Some(gid); }
            });
            ui.horizontal(|ui| {
                ui.add_space(depth as f32 * 14.0 + 16.0);
                let detail = if is_bool {
                    format!("{} · {}/{} ({:.1}%)", shape_lbl, n_in, n_parent, pct_par)
                } else {
                    format!("{} · {}×{} · {}/{} ({:.1}%)", shape_lbl, xch, ych, n_in, n_parent, pct_par)
                };
                ui.label(RichText::new(detail).small().color(Color32::GRAY));
            });
            // Reparent combo
            ui.horizontal(|ui| {
                ui.add_space(depth as f32 * 14.0 + 16.0);
                ui.label(RichText::new("parent:").small());
                let cur_parent_name = self.gates[idx].parent
                    .and_then(|pid| self.gates.iter().find(|g| g.id == pid).map(|g| g.name.clone()))
                    .unwrap_or_else(|| "All".into());
                egui::ComboBox::from_id_salt(format!("rep{}", gid))
                    .selected_text(cur_parent_name).width(120.0).show_ui(ui, |ui| {
                        if ui.selectable_label(false, "All events").clicked() {
                            reparent = Some((gid, None));
                        }
                        for g in &self.gates {
                            if g.id == gid { continue; }
                            if ui.selectable_label(false, &g.name).clicked() {
                                reparent = Some((gid, Some(g.id)));
                            }
                        }
                    });
            });

            // Numeric inspector for the selected gate (bounds edited in DATA units).
            if self.selected_gate == Some(gid) {
                // Snapshot the pre-edit tree so a drag on any field is undoable.
                let before = self.gates.clone();
                let xt = self.gates[idx].x_transform.compile();
                let yt = self.gates[idx].y_transform.compile();
                let ind = depth as f32 * 14.0 + 16.0;
                let quad_group = self.gates[idx].quad_group;
                if let Some(group) = quad_group {
                    // Linked quadrant: one center editor drives all four rectangles.
                    if let Some((cx, cy)) = self.quadrant_center(group) {
                        let mut dx = xt.inverse(cx);
                        let mut dy = yt.inverse(cy);
                        let (mut q_changed, mut q_started) = (false, false);
                        ui.horizontal(|ui| { ui.add_space(ind); ui.label(RichText::new("quad center x").small());
                            let r = ui.add(egui::DragValue::new(&mut dx).speed(10.0));
                            if r.drag_started() { q_started = true; }
                            if r.changed() { q_changed = true; } });
                        ui.horizontal(|ui| { ui.add_space(ind); ui.label(RichText::new("quad center y").small());
                            let r = ui.add(egui::DragValue::new(&mut dy).speed(10.0));
                            if r.drag_started() { q_started = true; }
                            if r.changed() { q_changed = true; } });
                        ui.horizontal(|ui| { ui.add_space(ind);
                            ui.label(RichText::new("linked quadrant — moves all 4").small().color(Color32::GRAY)); });
                        if q_started { self.push_undo_state(before); }
                        if q_changed {
                            self.set_quadrant_center(group, xt.forward(dx), yt.forward(dy));
                            self.needs_regate = true;
                            if self.active_pop.is_some() { self.scatter = None; }
                        }
                    }
                } else {
                    let mut changed = false;
                    let mut started = false;
                    // (changed, drag_started) for a DragValue editing a data-unit field.
                    let row = |ui: &mut egui::Ui, lab: &str, disp: &mut f64, t: &CompiledTransform| -> (bool, bool) {
                        let mut v = t.inverse(*disp);
                        let (mut ch, mut st) = (false, false);
                        ui.horizontal(|ui| {
                            ui.add_space(ind);
                            ui.label(RichText::new(lab).small());
                            let resp = ui.add(egui::DragValue::new(&mut v).speed(10.0));
                            if resp.drag_started() { st = true; }
                            if resp.changed() { *disp = t.forward(v); ch = true; }
                        });
                        (ch, st)
                    };
                    macro_rules! apply { ($e:expr) => {{ let (c, s) = $e; changed |= c; started |= s; }} }
                    match &mut self.gates[idx].shape {
                        GateShape::Rect { x_min, x_max, y_min, y_max } => {
                            apply!(row(ui, "x min", x_min, &xt)); apply!(row(ui, "x max", x_max, &xt));
                            apply!(row(ui, "y min", y_min, &yt)); apply!(row(ui, "y max", y_max, &yt));
                        }
                        GateShape::Range { x_min, x_max } => {
                            apply!(row(ui, "min", x_min, &xt)); apply!(row(ui, "max", x_max, &xt));
                        }
                        GateShape::Ellipse { cx, cy, rx, ry, angle } => {
                            apply!(row(ui, "center x", cx, &xt)); apply!(row(ui, "center y", cy, &yt));
                            // radii are in display units; edit directly
                            ui.horizontal(|ui| { ui.add_space(ind); ui.label(RichText::new("radius x (disp)").small());
                                let r = ui.add(egui::DragValue::new(rx).speed(0.01));
                                if r.drag_started() { started = true; }
                                if r.changed() { changed = true; } });
                            ui.horizontal(|ui| { ui.add_space(ind); ui.label(RichText::new("radius y (disp)").small());
                                let r = ui.add(egui::DragValue::new(ry).speed(0.01));
                                if r.drag_started() { started = true; }
                                if r.changed() { changed = true; } });
                            // rotation in degrees (stored as radians)
                            let mut deg = angle.to_degrees();
                            ui.horizontal(|ui| { ui.add_space(ind); ui.label(RichText::new("rotation (°)").small());
                                let r = ui.add(egui::DragValue::new(&mut deg).speed(1.0).range(-180.0..=180.0));
                                if r.drag_started() { started = true; }
                                if r.changed() { *angle = deg.to_radians(); changed = true; } });
                        }
                        GateShape::Polygon { .. } => {
                            ui.horizontal(|ui| { ui.add_space(ind);
                                ui.label(RichText::new("polygon — redraw to change").small().color(Color32::GRAY)); });
                        }
                        GateShape::Boolean { .. } => {
                            ui.horizontal(|ui| { ui.add_space(ind);
                                ui.label(RichText::new("boolean gate — delete & recreate to change").small().color(Color32::GRAY)); });
                        }
                    }
                    if started { self.push_undo_state(before); }
                    if changed { self.needs_regate = true; }
                }
                ui.horizontal(|ui| {
                    ui.add_space(ind);
                    if ui.button(format!("{} Export events → .fcs", icon::EXPORT))
                        .on_hover_text("Write this population's events to a new FCS file").clicked()
                    { export_gid = Some(gid); }
                });
            }
            ui.separator();
        }

        if let Some(gid) = zoom_to { self.zoom_to_gate(gid); }
        if let Some(gid) = export_gid { self.export_gated_fcs(gid); }
        if let Some(gid) = to_delete {
            self.push_undo();
            // Re-parent children of the deleted gate to its parent (avoid orphans).
            let parent_of = self.gates.iter().find(|g| g.id == gid).and_then(|g| g.parent);
            for g in &mut self.gates {
                if g.parent == Some(gid) { g.parent = parent_of; }
                // Drop the deleted gate from any boolean combination that referenced it.
                if let GateShape::Boolean { refs, .. } = &mut g.shape { refs.retain(|r| *r != gid); }
            }
            self.gates.retain(|g| g.id != gid);
            self.bool_refs.remove(&gid);
            if self.new_gate_parent == Some(gid) { self.new_gate_parent = parent_of; }
            if self.hist_sample_pop == Some(gid) { self.hist_sample_pop = None; self.hist_cache = None; }
            if self.qc_live_gate == Some(gid) { self.qc_live_gate = None; }
            if self.selected_gate == Some(gid) { self.selected_gate = None; }
            if self.active_pop == Some(gid) { self.active_pop = parent_of; self.scatter = None; }
            self.hist_hidden.remove(&gid);
            self.needs_regate = true;
        }
        if let Some((gid, new_parent)) = reparent {
            if !self.would_cycle(gid, new_parent) {
                self.push_undo();
                if let Some(g) = self.gates.iter_mut().find(|g| g.id == gid) { g.parent = new_parent; }
                self.needs_regate = true;
            } else {
                self.status = "Cannot set parent: would create a cycle".into();
            }
        }
    }

    fn draw_btn(&mut self, ui: &mut egui::Ui, mode: DrawMode, label: &str) {
        let active = self.draw_mode == mode;
        let txt = if active { RichText::new(label).color(Color32::from_rgb(220, 170, 0)) } else { RichText::new(label) };
        if ui.selectable_label(active, txt).clicked() {
            self.draw_mode = if active { DrawMode::Navigate } else { mode };
            self.drag_start = None; self.drag_current = None; self.poly_vertices.clear();
        }
    }

    /// Shared (display-coord) center of a quadrant group, read from any member.
    fn quadrant_center(&self, group: u32) -> Option<(f64, f64)> {
        const HALF: f64 = 5.0e11;
        for g in &self.gates {
            if g.quad_group == Some(group) {
                if let GateShape::Rect { x_min, x_max, y_min, y_max } = &g.shape {
                    let cx = if *x_min > -HALF { *x_min } else { *x_max };
                    let cy = if *y_min > -HALF { *y_min } else { *y_max };
                    return Some((cx, cy));
                }
            }
        }
        None
    }

    /// Move a whole quadrant group's shared center (display coords) — all four
    /// rectangles update together.
    fn set_quadrant_center(&mut self, group: u32, cx: f64, cy: f64) {
        const HALF: f64 = 5.0e11;
        for g in &mut self.gates {
            if g.quad_group == Some(group) {
                if let GateShape::Rect { x_min, x_max, y_min, y_max } = &mut g.shape {
                    if *x_min > -HALF { *x_min = cx; } else if *x_max < HALF { *x_max = cx; }
                    if *y_min > -HALF { *y_min = cy; } else if *y_max < HALF { *y_max = cy; }
                }
            }
        }
    }

    fn would_cycle(&self, gid: u32, new_parent: Option<u32>) -> bool {
        let mut cur = new_parent;
        let mut guard = 0;
        while let Some(id) = cur {
            if id == gid { return true; }
            guard += 1; if guard > 1000 { return true; }
            cur = self.gates.iter().find(|g| g.id == id).and_then(|g| g.parent);
        }
        false
    }

    fn save_gates(&mut self) {
        if self.gates.is_empty() { self.status = "No gates to save.".into(); return; }
        if let Some(path) = rfd::FileDialog::new()
            .add_filter("JSON", &["json"]).set_file_name("gates.json").save_file()
        {
            match serde_json::to_string_pretty(&self.gates).map_err(|e| e.to_string())
                .and_then(|s| std::fs::write(&path, s).map_err(|e| e.to_string()))
            {
                Ok(_) => self.status = format!("Saved {} gates → {}", self.gates.len(), path.display()),
                Err(e) => self.status = format!("Save error: {}", e),
            }
        }
    }

    fn load_gates(&mut self) {
        if let Some(path) = rfd::FileDialog::new().add_filter("JSON", &["json"]).pick_file() {
            match std::fs::read_to_string(&path).map_err(|e| e.to_string())
                .and_then(|s| serde_json::from_str::<Vec<Gate>>(&s).map_err(|e| e.to_string()))
            {
                Ok(gates) => {
                    self.push_undo();
                    self.next_gate_id = gates.iter().map(|g| g.id).max().unwrap_or(0) + 1;
                    self.gates = gates;
                    self.new_gate_parent = None;
                    self.needs_regate = true;
                    self.status = format!("Loaded {} gates from {}", self.gates.len(), path.display());
                }
                Err(e) => self.status = format!("Load error: {}", e),
            }
        }
    }

    // ── Session save / load ───────────────────────────────────────────

    fn save_session(&mut self) {
        if self.samples.is_empty() { self.status = "No workspace to save.".into(); return; }
        let session = Session {
            sample_paths: self.samples.iter().map(|s| s.path.clone()).collect(),
            groups: self.samples.iter().map(|s| s.group.clone()).collect(),
            active_sample: self.active_sample,
            do_compensate: self.do_compensate,
            dark_mode: self.dark_mode,
            viridis: self.colormap == ColorMap::Viridis,
            colormap: Some(self.colormap.label().to_string()),
            channel_tf: self.channel_tf.clone(),
            x_ch: self.x_ch,
            y_ch: self.y_ch,
            hist_ch: self.hist_ch,
            gates: self.gates.clone(),
            spill_override: self.spill_override.clone(),
        };
        if let Some(path) = rfd::FileDialog::new()
            .add_filter("JSON", &["json"]).set_file_name("session.json").save_file()
        {
            match serde_json::to_string_pretty(&session).map_err(|e| e.to_string())
                .and_then(|s| std::fs::write(&path, s).map_err(|e| e.to_string()))
            {
                Ok(_) => self.status = format!("Saved session ({} samples) → {}", self.samples.len(), path.display()),
                Err(e) => self.status = format!("Session save error: {}", e),
            }
        }
    }

    fn load_session(&mut self) {
        let path = match rfd::FileDialog::new().add_filter("JSON", &["json"]).pick_file() {
            Some(p) => p, None => return,
        };
        let session: Session = match std::fs::read_to_string(&path).map_err(|e| e.to_string())
            .and_then(|s| serde_json::from_str(&s).map_err(|e| e.to_string()))
        {
            Ok(s) => s,
            Err(e) => { self.status = format!("Session load error: {}", e); return; }
        };

        // Reset the workspace, then reopen the saved files (sample 0 fresh).
        self.samples.clear();
        self.fcs = None; self.compensated.clear();
        self.ref_sample = None; self.ref_scatter = None; self.batch = None;
        self.undo_stack.clear(); self.redo_stack.clear();
        // Keep (path, group) paired while dropping missing files, so tags stay
        // aligned, and track how the active index shifts past dropped predecessors.
        let mut surviving: Vec<(PathBuf, String)> = Vec::new();
        let mut missing = 0usize;
        let mut active_shift = 0usize; // survivors before the originally-active sample
        for (i, p) in session.sample_paths.iter().enumerate() {
            if p.exists() {
                if i < session.active_sample { active_shift += 1; }
                let g = session.groups.get(i).cloned().unwrap_or_default();
                surviving.push((p.clone(), g));
            } else {
                missing += 1;
            }
        }
        if surviving.is_empty() { self.status = "Session files not found on disk.".into(); return; }
        let paths: Vec<PathBuf> = surviving.iter().map(|(p, _)| p.clone()).collect();
        self.add_files(paths);

        // Restore tags + analysis state.
        for (s, (_, g)) in self.samples.iter_mut().zip(surviving.iter()) { s.group = g.clone(); }
        self.do_compensate = session.do_compensate;
        self.dark_mode = session.dark_mode;
        self.colormap = session.colormap.as_deref().and_then(ColorMap::from_name)
            .unwrap_or(if session.viridis { ColorMap::Viridis } else { ColorMap::Jet });
        self.spill_override = session.spill_override;
        // Validate a session-supplied override is square before any panel indexes it
        // (every other override entry point validates; a hand-edited session might not).
        if let Some(ov) = &self.spill_override {
            let n = ov.channels.len();
            if n == 0 || ov.rows.len() != n || ov.rows.iter().any(|r| r.len() != n) {
                self.spill_override = None;
            }
        }
        self.gates = session.gates;
        self.next_gate_id = self.gates.iter().map(|g| g.id).max().unwrap_or(0) + 1;
        // Reconcile gate-id-keyed UI state against the loaded gates so a stale id can't
        // dangle (effective_mask of an unknown id returns all-true → wrong histogram).
        if let Some(h) = self.hist_sample_pop { if !self.gates.iter().any(|g| g.id == h) { self.hist_sample_pop = None; } }
        if let Some(q) = self.qc_live_gate { if !self.gates.iter().any(|g| g.id == q) { self.qc_live_gate = None; } }
        self.hist_hidden.retain(|id| self.gates.iter().any(|g| g.id == *id));
        self.scatter_hidden.retain(|id| self.gates.iter().any(|g| g.id == *id));
        self.bool_refs.retain(|id| self.gates.iter().any(|g| g.id == *id));

        // Activate the saved sample (loads its events, keeps the gate tree), then
        // override the transforms/axes with the saved ones (same files → same panel).
        let active = active_shift.min(self.samples.len().saturating_sub(1));
        self.activate_sample(active, false);
        if session.channel_tf.len() == self.channel_tf.len() {
            self.channel_tf = session.channel_tf;
        }
        let np = self.fcs.as_ref().map(|f| f.n_params()).unwrap_or(1).max(1);
        self.x_ch = session.x_ch.min(np - 1);
        self.y_ch = session.y_ch.min(np - 1);
        self.hist_ch = session.hist_ch.min(np - 1);
        self.needs_reprocess = true; self.needs_rescatter = true; self.needs_regate = true;
        self.scatter = None; self.hist_cache = None;
        self.status = if missing == 0 {
            format!("Loaded session ({} samples) from {}", self.samples.len(), path.display())
        } else {
            format!("Loaded session — {} {} file(s) missing and skipped", icon::WARNING, missing)
        };
    }

    fn panel_status(&mut self, ctx: &egui::Context) {
        egui::TopBottomPanel::bottom("status").show(ctx, |ui| {
            ui.horizontal(|ui| {
                ui.label(RichText::new(&self.status).small());
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    if let Some(sc) = &self.scatter {
                        let shown: usize = sc.buckets.iter().map(|b| b.len()).sum();
                        let total = self.fcs.as_ref().map(|f| f.n_events).unwrap_or(0);
                        ui.label(RichText::new(format!("showing {}/{} events", shown, total))
                            .small().color(Color32::GRAY));
                    }
                    if let Some(c) = &self.cursor_label {
                        ui.separator();
                        ui.label(RichText::new(format!("⌖ {}", c)).small().monospace()
                            .color(Color32::from_rgb(120, 160, 210)));
                    }
                });
            });
        });
    }

    fn panel_central(&mut self, ctx: &egui::Context) {
        egui::CentralPanel::default().show(ctx, |ui| {
            match self.active_tab {
                ActiveTab::Plot => self.scatter_plot(ui),
                ActiveTab::Histogram => self.histogram_view(ui),
                ActiveTab::Stats => self.stats_table(ui),
                ActiveTab::Batch => self.batch_view(ui),
                ActiveTab::Spillover => self.spillover_view(ui),
                ActiveTab::QC => self.qc_view(ui),
            }
        });
    }

    // ── Adjustable grid of plots (up to 6×6) ──────────────────────────

    fn grid_view(&mut self, ui: &mut egui::Ui) {
        // Clear any stuck gesture owner when not actively drawing (Esc/Cancel/tool-
        // switch don't go through a normal drag-stop, which would otherwise leave
        // `active_grid_cell` set and lock out new gestures in every cell).
        if self.draw_mode == DrawMode::Navigate { self.active_grid_cell = None; }
        ui.label(RichText::new("Each cell has its own X/Y. Pick a draw tool (left), then draw/edit gates in any cell; “gate from here” drills all cells.")
            .small().color(Color32::GRAY));

        let cols = self.grid_cols.clamp(1, 6);
        let rows = self.grid_rows.clamp(1, 6);
        let n = cols * rows;
        let n_params = self.fcs.as_ref().map(|f| f.n_params()).unwrap_or(1).max(1);
        // Grow per-cell state to fit the grid; new cells get varied default pairs.
        while self.grid_channels.len() < n {
            let k = self.grid_channels.len();
            self.grid_channels.push((0, (k + 1).min(n_params - 1)));
        }
        while self.grid_cache.len() < n { self.grid_cache.push(None); }

        // Point budget per cell: keep the whole grid's drawn points bounded so
        // large grids stay smooth (cells are small, so fewer points still read well).
        let cap = (30_000 / n).clamp(800, MAX_SCATTER);

        let avail = ui.available_size();
        let gap = 8.0;
        let cw = (avail.x / cols as f32 - gap).max(110.0);
        // Reserve room for the two header rows (channels + scale pickers) per cell.
        const HEADER: f32 = 56.0;
        let ch = (avail.y / rows as f32 - HEADER).max(70.0);

        // Compute the shared kept-event list at most once per frame, only if a cell
        // actually rebuilds (avoids per-frame mask recompute when everything's cached).
        let mut kept: Option<Vec<usize>> = None;
        for row in 0..rows {
            ui.horizontal(|ui| {
                for col in 0..cols {
                    let idx = row * cols + col;
                    ui.allocate_ui(egui::vec2(cw, ch + HEADER), |ui| {
                        ui.vertical(|ui| { self.grid_cell(ui, idx, cw, ch, cap, &mut kept); });
                    });
                }
            });
        }
    }

    fn grid_cell(&mut self, ui: &mut egui::Ui, idx: usize, cw: f32, ch: f32, cap: usize, kept: &mut Option<Vec<usize>>) {
        let n_params = self.fcs.as_ref().map(|f| f.n_params()).unwrap_or(1).max(1);
        if idx >= self.grid_channels.len() { return; }
        let (mut xi, mut yi) = self.grid_channels[idx];
        xi = xi.min(n_params - 1);
        yi = yi.min(n_params - 1);

        // channel names for the cell pickers
        let names: Vec<String> = self.fcs.as_ref().map(|f|
            (0..f.n_params()).map(|i| param_full_label(f, i)).collect()).unwrap_or_default();

        ui.horizontal(|ui| {
            egui::ComboBox::from_id_salt(format!("gx{}", idx))
                .selected_text(names.get(xi).cloned().unwrap_or_default()).width(cw * 0.42)
                .show_ui(ui, |ui| { for (i, nm) in names.iter().enumerate() {
                    ui.selectable_value(&mut xi, i, nm); } });
            ui.label("×");
            egui::ComboBox::from_id_salt(format!("gy{}", idx))
                .selected_text(names.get(yi).cloned().unwrap_or_default()).width(cw * 0.42)
                .show_ui(ui, |ui| { for (i, nm) in names.iter().enumerate() {
                    ui.selectable_value(&mut yi, i, nm); } });
        });
        self.grid_channels[idx] = (xi, yi);
        // per-axis scale pickers (Linear/Log/Asinh/Logicle) right in the cell
        ui.horizontal(|ui| {
            ui.label(RichText::new("scale").small().color(Color32::GRAY));
            self.grid_axis_scale(ui, &format!("gxt{}", idx), xi);
            self.grid_axis_scale(ui, &format!("gyt{}", idx), yi);
        });

        // (re)build this cell's buckets if channels/transforms/data changed
        let cur_xt = self.cur_tf(xi);
        let cur_yt = self.cur_tf(yi);
        let stale = self.grid_cache[idx].as_ref().map(|c| {
            c.xi != xi || c.yi != yi || c.gen != self.data_gen || c.pop != self.active_pop
                || c.x_label != cur_xt.short_label() || c.y_label != cur_yt.short_label()
        }).unwrap_or(true);
        if stale {
            let k = kept.get_or_insert_with(|| self.grid_kept());
            let buckets = self.compute_cell_buckets(xi, yi, k, cap);
            self.grid_cache[idx] = Some(GridCell {
                xi, yi, x_label: cur_xt.short_label().to_string(),
                y_label: cur_yt.short_label().to_string(),
                pop: self.active_pop, gen: self.data_gen, buckets,
            });
        }

        let dark = self.dark_mode;
        let cmap = self.colormap;
        let xt = cur_xt.compile();
        let yt = cur_yt.compile();
        let buckets = self.grid_cache[idx].as_ref().map(|c| c.buckets.clone()).unwrap_or_default();
        let (x_name, y_name) = (names.get(xi).cloned().unwrap_or_default(), names.get(yi).cloned().unwrap_or_default());

        // gates on this channel pair (outline + label), remapped to current transforms
        let gate_draws: Vec<(u32, String, GateShape, CompiledTransform, CompiledTransform)> =
            self.gates.iter().filter_map(|g| {
                if !self.scatter_hidden.contains(&g.id)
                    && g.x_channel.eq_ignore_ascii_case(&x_name_base(&x_name))
                    && g.y_channel.eq_ignore_ascii_case(&x_name_base(&y_name))
                {
                    Some((g.id, g.name.clone(), g.shape.clone(), g.x_transform.compile(), g.y_transform.compile()))
                } else { None }
            }).collect();
        let gate_clamp = {
            let (cxr, cyr) = scatter_display_range(&buckets);
            let mx = (cxr.1 - cxr.0).abs().max(1e-9) * 0.05;
            let my = (cyr.1 - cyr.0).abs().max(1e-9) * 0.05;
            (cxr.0 - mx, cxr.1 + mx, cyr.0 - my, cyr.1 + my)
        };

        let xt_fmt = xt.clone(); let yt_fmt = yt.clone();
        let xt_grid = xt.clone(); let yt_grid = yt.clone();
        let lin_x = matches!(cur_xt, AxisTransform::Linear);
        let lin_y = matches!(cur_yt, AxisTransform::Linear);

        // ── interaction setup (mirrors the single plot, scoped to this cell) ──
        let drawing = self.draw_mode != DrawMode::Navigate;
        let mode = self.draw_mode;
        let is_active = self.active_grid_cell == Some(idx);
        let can_start = self.active_grid_cell.is_none();
        let edit_gate: Option<(usize, GateShape, CompiledTransform, CompiledTransform)> =
            if mode == DrawMode::Edit {
                self.selected_gate.and_then(|sid| self.gates.iter().position(|g| g.id == sid))
                    .and_then(|gi| {
                        let g = &self.gates[gi];
                        if g.x_channel.eq_ignore_ascii_case(&x_name_base(&x_name))
                            && g.y_channel.eq_ignore_ascii_case(&x_name_base(&y_name)) {
                            Some((gi, g.shape.clone(), g.x_transform.compile(), g.y_transform.compile()))
                        } else { None }
                    })
            } else { None };

        let cur_ds = self.drag_start;
        let cur_dc = self.drag_current;
        let cur_poly = self.poly_vertices.clone();
        let cur_grab = self.grab_handle;
        let cur_move_last = self.gate_move_last;
        let mut next_ds = cur_ds;
        let mut next_dc = cur_dc;
        let mut next_grab = cur_grab;
        let mut next_move_last = cur_move_last;
        let mut next_active = self.active_grid_cell;
        let mut new_shape: Option<GateShape> = None;
        let mut new_quadrant: Option<[f64; 2]> = None;
        let mut poly_add: Option<[f64; 2]> = None;
        let mut poly_finish = false;
        let mut exit_draw = false;
        let mut begin_edit = false;
        let mut handle_update: Option<(usize, usize, f64, f64)> = None;
        let mut gate_translate: Option<(usize, f64, f64)> = None;
        let mut hover_disp: Option<[f64; 2]> = None;
        let mut dbl_drill: Option<[f64; 2]> = None;
        // Locked view: freeze bounds to the data extent so the cell can't auto-fit/drift.
        let lock_view = self.lock_view;
        let view_bounds = if lock_view {
            let (dxr, dyr) = scatter_display_range(&buckets);
            let mx = (dxr.1 - dxr.0).abs().max(1e-9) * 0.05;
            let my = (dyr.1 - dyr.0).abs().max(1e-9) * 0.05;
            Some(PlotBounds::from_min_max([dxr.0 - mx, dyr.0 - my], [dxr.1 + mx, dyr.1 + my]))
        } else { None };

        let resp = Plot::new(format!("grid_{}_{}_{}_{}_{}_{}", idx, xi, yi, cur_xt.short_label(), cur_yt.short_label(), lock_view))
            .width(cw).height(ch).allow_scroll(false)
            .allow_drag(!drawing && !lock_view).allow_zoom(!drawing && !lock_view)
            .allow_double_click_reset(false) // double-click = drill into a gate
            .x_axis_label(&x_name).y_axis_label(&y_name)
            .x_axis_formatter(move |gm: GridMark, _r| fmt_data_tick(xt_fmt.inverse(gm.value)))
            .y_axis_formatter(move |gm: GridMark, _r| fmt_data_tick(yt_fmt.inverse(gm.value)))
            .x_grid_spacer(move |inp| nonlinear_grid(&xt_grid, lin_x, inp))
            .y_grid_spacer(move |inp| nonlinear_grid(&yt_grid, lin_y, inp))
            .show(ui, |pu| {
                if let Some(b) = view_bounds {
                    pu.set_plot_bounds(b);
                    pu.set_auto_bounds(egui::Vec2b::new(false, false));
                }
                for (k, pts) in buckets.iter().enumerate() {
                    if pts.is_empty() { continue; }
                    pu.points(Points::new(PlotPoints::new(pts.clone()))
                        .radius(1.2).color(density_color(k, N_BUCKETS, dark, cmap)));
                }
                for (gid, name, shape, gxt, gyt) in &gate_draws {
                    let (outline, fill) = gate_color(*gid);
                    let pts: Vec<[f64; 2]> = shape.outline().iter().map(|p| [
                        xt.forward(gxt.inverse(p[0])).clamp(gate_clamp.0, gate_clamp.1),
                        yt.forward(gyt.inverse(p[1])).clamp(gate_clamp.2, gate_clamp.3),
                    ]).collect();
                    if pts.len() >= 2 {
                        pu.polygon(PlotPolygon::new(PlotPoints::new(pts))
                            .stroke(Stroke::new(1.4, outline)).fill_color(fill));
                        let a = shape.label_anchor();
                        let ax = xt.forward(gxt.inverse(a[0])).clamp(gate_clamp.0, gate_clamp.1);
                        let ay = yt.forward(gyt.inverse(a[1])).clamp(gate_clamp.2, gate_clamp.3);
                        pu.text(PlotText::new(PlotPoint::new(ax, ay), RichText::new(name).color(outline).size(11.0)));
                    }
                }

                hover_disp = pu.pointer_coordinate().map(|p| [p.x, p.y]);
                if !drawing && pu.response().double_clicked() {
                    dbl_drill = pu.pointer_coordinate().map(|p| [p.x, p.y]);
                }
                if drawing {
                    let ptr = pu.pointer_coordinate();
                    let bounds = pu.plot_bounds();
                    let (r_started, r_dragged, r_stopped, r_clicked, r_dbl) = {
                        let r = pu.response();
                        (r.drag_started(), r.dragged(), r.drag_stopped(), r.clicked(), r.double_clicked())
                    };
                    match mode {
                        DrawMode::Rect | DrawMode::Ellipse => {
                            if is_active {
                                if let (Some(s), Some(c)) = (cur_ds, cur_dc) {
                                    pu.line(Line::new(PlotPoints::new(rubber_band(mode, s, c)))
                                        .color(Color32::from_rgb(240, 190, 40)).width(1.5));
                                }
                            }
                            if r_started && can_start { next_active = Some(idx); next_ds = ptr.map(|p| [p.x, p.y]); next_dc = next_ds; }
                            if r_dragged && is_active { if let Some(p) = ptr { next_dc = Some([p.x, p.y]); } }
                            if r_stopped && is_active {
                                let end = ptr.map(|p| [p.x, p.y]).or(next_dc);
                                if let (Some(s), Some(c)) = (next_ds, end) {
                                    let (w, h) = (bounds.width().max(1e-9), bounds.height().max(1e-9));
                                    if ((c[0]-s[0])/w).powi(2) + ((c[1]-s[1])/h).powi(2) > MIN_DRAG_SQ {
                                        new_shape = Some(shape_from_drag(mode, s, c)); exit_draw = true;
                                    }
                                }
                                next_ds = None; next_dc = None; next_active = None;
                            }
                        }
                        DrawMode::Polygon => {
                            let (w, h) = (bounds.width().max(1e-9), bounds.height().max(1e-9));
                            let near_first = |p: PlotPoint| cur_poly.first()
                                .map(|f| ((p.x - f[0])/w).powi(2) + ((p.y - f[1])/h).powi(2) < POLY_CLOSE_SQ)
                                .unwrap_or(false);
                            if is_active && !cur_poly.is_empty() {
                                let mut line = cur_poly.clone();
                                if let Some(p) = ptr { line.push([p.x, p.y]); }
                                pu.line(Line::new(PlotPoints::new(line)).color(Color32::from_rgb(240, 190, 40)).width(1.5));
                                pu.points(Points::new(PlotPoints::new(cur_poly.clone())).radius(3.0).color(Color32::from_rgb(240, 190, 40)));
                                // Highlight the first vertex as a close target (green when the cursor is in range).
                                if cur_poly.len() >= 3 {
                                    let hot = ptr.map(near_first).unwrap_or(false);
                                    let col = if hot { Color32::from_rgb(80, 220, 120) } else { Color32::from_rgb(240, 190, 40) };
                                    pu.points(Points::new(PlotPoints::new(vec![cur_poly[0]])).radius(if hot {7.0} else {5.0}).color(col));
                                }
                            }
                            // Place a vertex on click OR release (forgiving if the mouse moves while clicking).
                            if (r_clicked || r_stopped) && (is_active || (can_start && cur_poly.is_empty())) {
                                if let Some(p) = ptr {
                                    if cur_poly.len() >= 3 && near_first(p) { poly_finish = true; }
                                    else { if cur_poly.is_empty() { next_active = Some(idx); } poly_add = Some([p.x, p.y]); }
                                }
                            }
                            if r_dbl && is_active && cur_poly.len() >= 3 { poly_finish = true; }
                        }
                        DrawMode::Quadrant => {
                            if let Some(p) = ptr {
                                let (x0, x1, y0, y1) = (bounds.min()[0], bounds.max()[0], bounds.min()[1], bounds.max()[1]);
                                pu.line(Line::new(PlotPoints::new(vec![[p.x, y0], [p.x, y1]])).color(Color32::from_rgb(240, 190, 40)).width(1.2));
                                pu.line(Line::new(PlotPoints::new(vec![[x0, p.y], [x1, p.y]])).color(Color32::from_rgb(240, 190, 40)).width(1.2));
                            }
                            if r_clicked { if let Some(p) = ptr { new_quadrant = Some([p.x, p.y]); } }
                        }
                        DrawMode::Edit => {
                            if let Some((gi, shape, gxt, gyt)) = &edit_gate {
                                let handles: Vec<[f64; 2]> = gate_handles(shape).iter()
                                    .map(|p| [xt.forward(gxt.inverse(p[0])).clamp(gate_clamp.0, gate_clamp.1),
                                              yt.forward(gyt.inverse(p[1])).clamp(gate_clamp.2, gate_clamp.3)]).collect();
                                pu.points(Points::new(PlotPoints::new(handles.clone())).radius(5.0).color(Color32::from_rgb(240, 190, 40)));
                                let w = bounds.width().max(1e-9); let h = bounds.height().max(1e-9);
                                let to_gate = |p: PlotPoint| [gxt.forward(xt.inverse(p.x)), gyt.forward(yt.inverse(p.y))];
                                if r_started && can_start {
                                    if let Some(p) = ptr {
                                        let mut best: Option<(usize, f64)> = None;
                                        for (i, hp) in handles.iter().enumerate() {
                                            let d = ((p.x - hp[0]) / w).powi(2) + ((p.y - hp[1]) / h).powi(2);
                                            if d < 0.0025 && best.map(|(_, bd)| d < bd).unwrap_or(true) { best = Some((i, d)); }
                                        }
                                        next_grab = best.map(|(i, _)| i);
                                        if next_grab.is_none() {
                                            let g = to_gate(p);
                                            next_move_last = if shape.contains(g[0], g[1]) { Some(g) } else { None };
                                        } else { next_move_last = None; }
                                        if next_grab.is_some() || next_move_last.is_some() { begin_edit = true; next_active = Some(idx); }
                                    }
                                }
                                if r_dragged && is_active {
                                    if let (Some(p), Some(hg)) = (ptr, next_grab) {
                                        let g = to_gate(p); handle_update = Some((*gi, hg, g[0], g[1]));
                                    } else if let (Some(p), Some(last)) = (ptr, next_move_last) {
                                        let g = to_gate(p); gate_translate = Some((*gi, g[0] - last[0], g[1] - last[1])); next_move_last = Some(g);
                                    }
                                }
                                if r_stopped && is_active { next_grab = None; next_move_last = None; next_active = None; }
                            }
                        }
                        DrawMode::Navigate => {}
                    }
                }
            });

        // ── apply interaction (only the active cell yields commits) ──
        self.drag_start = next_ds;
        self.drag_current = next_dc;
        self.active_grid_cell = next_active;
        if begin_edit { self.push_undo(); }
        self.grab_handle = next_grab;
        self.gate_move_last = next_move_last;
        if let Some(v) = poly_add { self.poly_vertices.push(v); }
        if poly_finish { self.finish_polygon(); }
        if let Some(shape) = new_shape { self.commit_gate_on(xi, yi, shape); self.draw_mode = DrawMode::Navigate; }
        if let Some(c) = new_quadrant { self.commit_quadrant_on(xi, yi, c[0], c[1]); self.draw_mode = DrawMode::Navigate; }
        if exit_draw { self.draw_mode = DrawMode::Navigate; }
        if let Some((gi, hh, gx, gy)) = handle_update {
            if gi < self.gates.len() { apply_gate_handle(&mut self.gates[gi].shape, hh, gx, gy); self.needs_regate = true; }
        }
        if let Some((gi, dx, dy)) = gate_translate {
            if gi < self.gates.len() { translate_shape(&mut self.gates[gi].shape, dx, dy); self.needs_regate = true; }
        }
        if let Some(p) = dbl_drill { self.drill_at(&x_name, &y_name, &xt, &yt, p); }
        if resp.response.hovered() {
            if let Some(d) = hover_disp {
                self.cursor_label = Some(format!("{} {} · {} {}",
                    short_chan(&x_name_base(&x_name)), fmt_data_tick(xt.inverse(d[0])),
                    short_chan(&x_name_base(&y_name)), fmt_data_tick(yt.inverse(d[1]))));
            }
        }
    }

    /// Inline spillover adjuster for the current X↔Y pair, with live plot update.
    /// Targets the classic "this channel under-corrects into that one" problem
    /// (e.g. MHCII → CD11b) without leaving the Plot tab.
    fn ui_comp_preview(&mut self, ui: &mut egui::Ui, x_name: &str, y_name: &str) {
        let xb = x_name_base(x_name);
        let yb = x_name_base(y_name);
        egui::CollapsingHeader::new(format!("{} Compensation (current X↔Y)", icon::SCALES)).id_salt("comp_preview").show(ui, |ui| {
            let mat = self.active_matrix();
            let channels = match &mat { Some((c, _)) => c.clone(), None => {
                ui.label(RichText::new("No spillover matrix for this file — create one on the Spillover tab.")
                    .small().color(Color32::GRAY));
                return;
            }};
            if !self.do_compensate {
                ui.label(RichText::new("Tip: turn on “Compensate” in the toolbar to see the effect live.")
                    .small().color(Color32::from_rgb(220, 170, 60)));
            }
            let xi = channels.iter().position(|c| c.eq_ignore_ascii_case(&xb));
            let yi = channels.iter().position(|c| c.eq_ignore_ascii_case(&yb));
            let (xi, yi) = match (xi, yi) {
                (Some(a), Some(b)) if a != b => (a, b),
                _ => {
                    ui.label(RichText::new("Put two different fluorescence channels on X and Y to adjust their spillover.")
                        .small().color(Color32::GRAY));
                    return;
                }
            };
            let rows = mat.as_ref().map(|(_, r)| r.clone()).unwrap_or_default();
            let (mut sxy, mut syx) = (rows[xi][yi], rows[yi][xi]);
            let mut changed = false;
            ui.horizontal(|ui| {
                ui.label(RichText::new(format!("{} → {}", short_chan(&xb), short_chan(&yb))).small());
                if ui.add(egui::DragValue::new(&mut sxy).speed(0.001).fixed_decimals(4).range(-2.0..=2.0)).changed() { changed = true; }
                ui.label(RichText::new(format!("   {} → {}", short_chan(&yb), short_chan(&xb))).small());
                if ui.add(egui::DragValue::new(&mut syx).speed(0.001).fixed_decimals(4).range(-2.0..=2.0)).changed() { changed = true; }
            });
            if changed {
                if self.spill_override.is_none() { self.start_override(); }
                if let Some(ov) = &mut self.spill_override {
                    let a = ov.channels.iter().position(|c| c.eq_ignore_ascii_case(&xb));
                    let b = ov.channels.iter().position(|c| c.eq_ignore_ascii_case(&yb));
                    if let (Some(a), Some(b)) = (a, b) { ov.rows[a][b] = sxy; ov.rows[b][a] = syx; }
                }
                if !self.do_compensate { self.do_compensate = true; }
                self.needs_reprocess = true;
            }
            if self.spill_override.is_some() {
                ui.horizontal(|ui| {
                    ui.label(RichText::new("● override active").small().color(Color32::from_rgb(230, 140, 40)));
                    if ui.small_button(format!("{} reset to embedded", icon::ARROW_COUNTER_CLOCKWISE)).clicked() {
                        self.spill_override = None; self.needs_reprocess = true;
                    }
                });
            }
        });
    }

    // ── Scatter plot ──────────────────────────────────────────────────

    /// Centered empty-state shown in the plot area before any file is loaded.
    fn empty_state(&mut self, ui: &mut egui::Ui) {
        ui.vertical_centered(|ui| {
            ui.add_space((ui.available_height() * 0.30).max(20.0));
            ui.label(RichText::new(icon::CHART_SCATTER).size(58.0).color(Color32::from_rgb(38, 162, 156)));
            ui.add_space(10.0);
            ui.label(RichText::new("No data loaded").heading());
            ui.add_space(2.0);
            ui.label(RichText::new("Open an FCS file to plot and gate your samples.").color(Color32::GRAY));
            ui.add_space(16.0);
            if ui.button(RichText::new(format!("{}  Open FCS…", icon::FILE_PLUS)).size(15.0)).clicked() {
                if let Some(paths) = rfd::FileDialog::new()
                    .add_filter("FCS files", &["fcs", "FCS"]).pick_files()
                {
                    self.add_files(paths);
                }
            }
            ui.add_space(6.0);
            ui.label(RichText::new("…or drag a .fcs file onto the window").small().weak());
        });
    }

    fn scatter_plot(&mut self, ui: &mut egui::Ui) {
        if self.fcs.is_none() {
            self.empty_state(ui);
            return;
        }

        // plot controls
        ui.horizontal(|ui| {
            ui.label("Layout:");
            let was_grid = self.grid_mode;
            ui.selectable_value(&mut self.grid_mode, false, "Single");
            ui.selectable_value(&mut self.grid_mode, true, "Grid");
            if self.grid_mode {
                let (c0, r0) = (self.grid_cols, self.grid_rows);
                ui.add(egui::DragValue::new(&mut self.grid_cols).range(1..=6).speed(0.05).prefix("cols "));
                ui.label("×");
                ui.add(egui::DragValue::new(&mut self.grid_rows).range(1..=6).speed(0.05).prefix("rows "));
                // Cell point-budgets change with grid size → rebuild all cells.
                if was_grid != self.grid_mode || self.grid_cols != c0 || self.grid_rows != r0 {
                    self.grid_cache.iter_mut().for_each(|c| *c = None);
                }
            }
            ui.separator();
            if !self.grid_mode {
                if ui.checkbox(&mut self.show_contours, "Contours").changed() { self.scatter = None; }
                ui.checkbox(&mut self.density_fill, "Fill")
                    .on_hover_text("Filled density heatmap instead of dots");
                ui.separator();
            }
            ui.checkbox(&mut self.lock_view, format!("{} Lock view", icon::LOCK_SIMPLE))
                .on_hover_text("Freeze pan/zoom so the plot holds still while gating (uncheck to drag-pan / pinch-zoom)");
            if !self.grid_mode && !self.lock_view
                && ui.button(format!("{} Fit", icon::ARROWS_OUT)).on_hover_text("Reset zoom to fit all data").clicked()
            {
                self.x_manual = false; self.y_manual = false; self.fit_view = true; self.scatter = None;
            }
            ui.separator();
            ui.label("Colormap:");
            egui::ComboBox::from_id_salt("cmap").selected_text(self.colormap.label()).show_ui(ui, |ui| {
                for cm in ColorMap::ALL {
                    ui.selectable_value(&mut self.colormap, cm, cm.label());
                }
            });
            // PNG export captures one plot; not meaningful for the grid overview.
            if !self.grid_mode {
                ui.separator();
                if ui.button(format!("{} Save plot…", icon::CAMERA)).clicked() { self.request_plot_png(); }
            }
        });

        if self.grid_mode { self.grid_view(ui); return; }

        let n_params = self.fcs.as_ref().map(|f| f.n_params()).unwrap_or(1).max(1);
        let xi = self.x_ch.min(n_params - 1);
        let yi = self.y_ch.min(n_params - 1);
        let cur_xt = self.cur_tf(xi);
        let cur_yt = self.cur_tf(yi);
        let xt = cur_xt.compile();
        let yt = cur_yt.compile();
        let dark = self.dark_mode;
        let cmap = self.colormap;

        let (x_name, y_name) = {
            let f = self.fcs.as_ref().unwrap();
            (param_full_label(f, xi), param_full_label(f, yi))
        };

        self.ui_comp_preview(ui, &x_name, &y_name);

        // Stale-scatter guard: rebuild if channel/transform changed without a flag.
        let stale = self.scatter.as_ref().map(|s| {
            s.x_ch != xi || s.y_ch != yi
                || s.x_label != cur_xt.short_label()
                || s.y_label != cur_yt.short_label()
                || s.pop != self.active_pop
        }).unwrap_or(true);
        if stale { self.rebuild_scatter(); }

        // Reference-overlay: rebuild if the chosen sample / channels / transforms changed.
        if self.ref_sample.is_some() {
            let ref_stale = self.ref_scatter.as_ref().map(|r| {
                Some(r.ref_idx) != self.ref_sample || r.x_ch != xi || r.y_ch != yi
                    || r.x_label != cur_xt.short_label() || r.y_label != cur_yt.short_label()
            }).unwrap_or(true);
            if ref_stale { self.rebuild_ref_scatter(); }
        } else if self.ref_scatter.is_some() {
            self.ref_scatter = None;
        }
        let ref_pts: Vec<[f64; 2]> = self.ref_scatter.as_ref().map(|r| r.points.clone()).unwrap_or_default();
        let ref_name: String = self.ref_sample.and_then(|i| self.samples.get(i)).map(|s| s.name.clone()).unwrap_or_default();

        let buckets: Vec<Vec<[f64; 2]>> = self.scatter.as_ref().map(|s| s.buckets.clone()).unwrap_or_default();
        let back_pts: Vec<[f64; 2]> = self.scatter.as_ref().map(|s| s.back_pts.clone()).unwrap_or_default();
        let contours: Vec<[[f64; 2]; 2]> = self.scatter.as_ref().map(|s| s.contours.clone()).unwrap_or_default();

        // Filled-density heatmap: build a texture from the cached density grid.
        let heat_image: Option<(egui::TextureHandle, PlotPoint, egui::Vec2)> = if self.density_fill {
            self.scatter.as_ref().and_then(|s| s.heat.as_ref()).map(|h| {
                let n = h.n;
                let mut px = vec![Color32::TRANSPARENT; n * n];
                for r in 0..n {            // image row 0 = top = highest-y bin
                    let by = n - 1 - r;
                    for c in 0..n {
                        let t = h.t[by * n + c];
                        if t > 0.0 { px[r * n + c] = density_color_t(t, dark, cmap); }
                    }
                }
                let img = egui::ColorImage { size: [n, n], pixels: px };
                let tex = ui.ctx().load_texture("density_heat", img, egui::TextureOptions::LINEAR);
                let [x0, x1, y0, y1] = h.extent;
                let center = PlotPoint::new((x0 + x1) / 2.0, (y0 + y1) / 2.0);
                let size = egui::vec2((x1 - x0) as f32, (y1 - y0) as f32);
                (tex, center, size)
            })
        } else { None };

        // Gate render data for the CURRENT channel pair (any transform — remap below).
        let gate_draws: Vec<(u32, String, GateShape, CompiledTransform, CompiledTransform)> =
            self.gates.iter().filter_map(|g| {
                if !self.scatter_hidden.contains(&g.id)
                    && g.x_channel.eq_ignore_ascii_case(&x_name_base(&x_name))
                    && g.y_channel.eq_ignore_ascii_case(&x_name_base(&y_name))
                {
                    // label = name + "%parent (count)" from the gate-count cache
                    let label = match self.gate_counts.get(&g.id) {
                        Some((n_in, n_par)) => {
                            let pct = if *n_par > 0 { 100.0 * *n_in as f64 / *n_par as f64 } else { 0.0 };
                            format!("{}\n{:.1}% ({})", g.name, pct, n_in)
                        }
                        None => g.name.clone(),
                    };
                    Some((g.id, label, g.shape.clone(),
                          g.x_transform.compile(), g.y_transform.compile()))
                } else { None }
            }).collect();

        // Manual bounds in display coords.
        let manual_bounds: Option<PlotBounds> = {
            if self.x_manual || self.y_manual {
                // need both axes; fall back to auto on the unset axis via current scatter range
                let (dxr, dyr) = scatter_display_range(&buckets);
                let xlo = if self.x_manual { xt.forward(self.x_lo) } else { dxr.0 };
                let xhi = if self.x_manual { xt.forward(self.x_hi) } else { dxr.1 };
                let ylo = if self.y_manual { yt.forward(self.y_lo) } else { dyr.0 };
                let yhi = if self.y_manual { yt.forward(self.y_hi) } else { dyr.1 };
                Some(PlotBounds::from_min_max([xlo, ylo], [xhi, yhi]))
            } else { None }
        };
        // "Fit" reframes once to the data extent (one-shot); locked view freezes
        // bounds to the data extent (+margin) so the plot can't auto-fit/drift as the
        // rubber-band, gates, or cursor preview change.
        let data_extent = || {
            let (dxr, dyr) = scatter_display_range(&buckets);
            let mx = (dxr.1 - dxr.0).abs().max(1e-9) * 0.05;
            let my = (dyr.1 - dyr.0).abs().max(1e-9) * 0.05;
            PlotBounds::from_min_max([dxr.0 - mx, dyr.0 - my], [dxr.1 + mx, dyr.1 + my])
        };
        let view_bounds = if self.fit_view {
            self.fit_view = false; // consume the one-shot
            Some(data_extent())
        } else {
            manual_bounds.or_else(|| if self.lock_view { Some(data_extent()) } else { None })
        };

        // Clamp box for gate outlines so open/±∞ bounds (e.g. quadrant gates) don't
        // blow up the plot's auto-fit. Gate *membership* still uses the real bounds.
        let gate_clamp = {
            let (cxr, cyr) = scatter_display_range(&buckets);
            let mx = (cxr.1 - cxr.0).abs().max(1e-9) * 0.05;
            let my = (cyr.1 - cyr.0).abs().max(1e-9) * 0.05;
            (cxr.0 - mx, cxr.1 + mx, cyr.0 - my, cyr.1 + my)
        };

        let drawing = self.draw_mode != DrawMode::Navigate;
        let lock_view = self.lock_view;
        let mode = self.draw_mode;

        // Edit mode: the selected gate (if it's on the current axes), for handle dragging.
        let edit_gate: Option<(usize, GateShape, CompiledTransform, CompiledTransform)> =
            if mode == DrawMode::Edit {
                self.selected_gate.and_then(|sid| self.gates.iter().position(|g| g.id == sid))
                    .and_then(|idx| {
                        let g = &self.gates[idx];
                        if g.x_channel.eq_ignore_ascii_case(&x_name_base(&x_name))
                            && g.y_channel.eq_ignore_ascii_case(&x_name_base(&y_name)) {
                            Some((idx, g.shape.clone(), g.x_transform.compile(), g.y_transform.compile()))
                        } else { None }
                    })
            } else { None };

        // Interaction scratch
        let cur_ds = self.drag_start;
        let cur_dc = self.drag_current;
        let cur_poly = self.poly_vertices.clone();
        let mut next_ds = cur_ds;
        let mut next_dc = cur_dc;
        let mut new_shape: Option<GateShape> = None;
        let mut poly_add: Option<[f64; 2]> = None;
        let mut poly_finish = false;
        let mut exit_draw = false;
        let mut new_quadrant: Option<[f64; 2]> = None;
        let cur_grab = self.grab_handle;
        let mut next_grab = cur_grab;
        let mut handle_update: Option<(usize, usize, f64, f64)> = None; // (gate_idx, handle, gx, gy)
        let mut hover_disp: Option<[f64; 2]> = None; // cursor position in display coords
        let cur_move_last = self.gate_move_last;
        let mut next_move_last = cur_move_last;
        let mut gate_translate: Option<(usize, f64, f64)> = None; // (gate_idx, dx, dy) in gate-display coords
        let mut begin_edit = false; // an Edit gesture (resize/move/rotate) started this frame → snapshot undo

        // Axis formatters (data units) — clone compiled transforms into closures.
        let xt_fmt = xt.clone();
        let yt_fmt = yt.clone();
        let xt_grid = xt.clone();
        let yt_grid = yt.clone();
        let lin_x = matches!(cur_xt, AxisTransform::Linear);
        let lin_y = matches!(cur_yt, AxisTransform::Linear);

        // Manual-limit flags are part of the plot id: toggling them OFF spawns a
        // fresh plot whose auto_bounds defaults back to true (re-fits to data),
        // instead of staying frozen at the bounds we locked with set_auto_bounds(false).
        let plot = Plot::new(format!(
            "scatter_{}_{}_{}_{}_{}_{}_{}",
            xi, yi, cur_xt.short_label(), cur_yt.short_label(), self.x_manual, self.y_manual, lock_view
        ))
            .legend(Legend::default())
            .x_axis_label(&x_name)
            .y_axis_label(&y_name)
            .allow_drag(!drawing && !lock_view)
            .allow_zoom(!drawing && !lock_view)
            .allow_scroll(false)
            .allow_double_click_reset(false) // double-click = drill into a gate
            .x_axis_formatter(move |gm: GridMark, _r| fmt_data_tick(xt_fmt.inverse(gm.value)))
            .y_axis_formatter(move |gm: GridMark, _r| fmt_data_tick(yt_fmt.inverse(gm.value)))
            .x_grid_spacer(move |inp| nonlinear_grid(&xt_grid, lin_x, inp))
            .y_grid_spacer(move |inp| nonlinear_grid(&yt_grid, lin_y, inp));

        let mut dbl_drill: Option<[f64; 2]> = None;
        let plot_response = plot.show(ui, |pu| {
            if let Some(b) = view_bounds {
                pu.set_plot_bounds(b);
                pu.set_auto_bounds(egui::Vec2b::new(false, false));
            }
            hover_disp = pu.pointer_coordinate().map(|p| [p.x, p.y]);
            if !drawing && pu.response().double_clicked() {
                dbl_drill = pu.pointer_coordinate().map(|p| [p.x, p.y]);
            }

            // backgate: parent population events (faint blue-grey, for context)
            if !back_pts.is_empty() {
                let c = if dark { Color32::from_rgba_unmultiplied(120, 140, 175, 55) }
                        else { Color32::from_rgba_unmultiplied(90, 110, 150, 55) };
                pu.points(Points::new(PlotPoints::new(back_pts.clone())).radius(1.1).color(c)
                    .name("parent (backgate)"));
            }

            // reference overlay (faded grey, behind the active sample)
            if !ref_pts.is_empty() {
                let grey = if dark { Color32::from_rgba_unmultiplied(170, 170, 170, 70) }
                           else { Color32::from_rgba_unmultiplied(110, 110, 110, 70) };
                pu.points(Points::new(PlotPoints::new(ref_pts.clone())).radius(1.2).color(grey)
                    .name(format!("ref: {}", ref_name)));
            }

            // density: filled heatmap, or sampled dots
            if let Some((tex, center, size)) = &heat_image {
                pu.image(PlotImage::new(tex.id(), *center, *size));
            } else {
                for (k, pts) in buckets.iter().enumerate() {
                    if pts.is_empty() { continue; }
                    pu.points(Points::new(PlotPoints::new(pts.clone()))
                        .radius(1.4).color(density_color(k, N_BUCKETS, dark, cmap)));
                }
            }

            // iso-density contour lines (on top of the dots)
            if !contours.is_empty() {
                let cc = if dark { Color32::from_rgba_unmultiplied(235, 235, 235, 130) }
                         else { Color32::from_rgba_unmultiplied(40, 40, 40, 130) };
                for seg in &contours {
                    pu.line(Line::new(PlotPoints::new(vec![seg[0], seg[1]])).color(cc).width(0.8));
                }
            }

            // gate overlays (remap stored→data→current display)
            for (gid, name, shape, gxt, gyt) in &gate_draws {
                let (outline, fill) = gate_color(*gid);
                let pts: Vec<[f64; 2]> = shape.outline().iter()
                    .map(|p| [
                        xt.forward(gxt.inverse(p[0])).clamp(gate_clamp.0, gate_clamp.1),
                        yt.forward(gyt.inverse(p[1])).clamp(gate_clamp.2, gate_clamp.3),
                    ])
                    .collect();
                if pts.len() >= 2 {
                    pu.polygon(PlotPolygon::new(PlotPoints::new(pts.clone()))
                        .stroke(Stroke::new(1.6, outline)).fill_color(fill));
                    let anchor = shape.label_anchor();
                    let ax = xt.forward(gxt.inverse(anchor[0])).clamp(gate_clamp.0, gate_clamp.1);
                    let ay = yt.forward(gyt.inverse(anchor[1])).clamp(gate_clamp.2, gate_clamp.3);
                    pu.text(PlotText::new(PlotPoint::new(ax, ay), RichText::new(name).color(outline).size(12.0)));
                }
            }

            // in-progress drawing — capture response state as owned values up
            // front so no borrow of `pu` is held across the draw (pu.line/points) calls.
            if drawing {
                let ptr = pu.pointer_coordinate();
                let bounds = pu.plot_bounds();
                let (r_started, r_dragged, r_stopped, r_clicked, r_dbl) = {
                    let resp = pu.response();
                    (resp.drag_started(), resp.dragged(), resp.drag_stopped(),
                     resp.clicked(), resp.double_clicked())
                };

                match mode {
                    DrawMode::Rect | DrawMode::Ellipse => {
                        if let (Some(s), Some(c)) = (cur_ds, cur_dc) {
                            let band = rubber_band(mode, s, c);
                            pu.line(Line::new(PlotPoints::new(band))
                                .color(Color32::from_rgb(240, 190, 40)).width(1.5));
                        }
                        if r_started { next_ds = ptr.map(|p| [p.x, p.y]); next_dc = next_ds; }
                        if r_dragged { if let Some(p) = ptr { next_dc = Some([p.x, p.y]); } }
                        if r_stopped {
                            let end = ptr.map(|p| [p.x, p.y]).or(next_dc);
                            if let (Some(s), Some(c)) = (next_ds, end) {
                                let (w, h) = (bounds.width().max(1e-9), bounds.height().max(1e-9));
                                if ((c[0]-s[0])/w).powi(2) + ((c[1]-s[1])/h).powi(2) > MIN_DRAG_SQ {
                                    new_shape = Some(shape_from_drag(mode, s, c));
                                    exit_draw = true;
                                }
                            }
                            next_ds = None; next_dc = None;
                        }
                    }
                    DrawMode::Polygon => {
                        let (w, h) = (bounds.width().max(1e-9), bounds.height().max(1e-9));
                        let near_first = |p: PlotPoint| cur_poly.first()
                            .map(|f| ((p.x - f[0])/w).powi(2) + ((p.y - f[1])/h).powi(2) < POLY_CLOSE_SQ)
                            .unwrap_or(false);
                        if !cur_poly.is_empty() {
                            let mut line: Vec<[f64; 2]> = cur_poly.clone();
                            if let Some(p) = ptr { line.push([p.x, p.y]); }
                            pu.line(Line::new(PlotPoints::new(line))
                                .color(Color32::from_rgb(240, 190, 40)).width(1.5));
                            pu.points(Points::new(PlotPoints::new(cur_poly.clone()))
                                .radius(3.0).color(Color32::from_rgb(240, 190, 40)));
                            // Highlight the first vertex as a close target (green when the cursor is in range).
                            if cur_poly.len() >= 3 {
                                let hot = ptr.map(near_first).unwrap_or(false);
                                let col = if hot { Color32::from_rgb(80, 220, 120) } else { Color32::from_rgb(240, 190, 40) };
                                pu.points(Points::new(PlotPoints::new(vec![cur_poly[0]]))
                                    .radius(if hot { 7.0 } else { 5.0 }).color(col));
                            }
                        }
                        // Place a vertex on click OR release (forgiving if the mouse moves while clicking).
                        if r_clicked || r_stopped {
                            if let Some(p) = ptr {
                                if cur_poly.len() >= 3 && near_first(p) { poly_finish = true; }
                                else { poly_add = Some([p.x, p.y]); }
                            }
                        }
                        if r_dbl && cur_poly.len() >= 3 { poly_finish = true; }
                    }
                    DrawMode::Edit => {
                        if let Some((gidx, shape, gxt, gyt)) = &edit_gate {
                            let handles: Vec<[f64; 2]> = gate_handles(shape).iter()
                                .map(|p| [xt.forward(gxt.inverse(p[0])).clamp(gate_clamp.0, gate_clamp.1),
                                          yt.forward(gyt.inverse(p[1])).clamp(gate_clamp.2, gate_clamp.3)])
                                .collect();
                            pu.points(Points::new(PlotPoints::new(handles.clone()))
                                .radius(5.0).color(Color32::from_rgb(240, 190, 40)));
                            let w = bounds.width().max(1e-9);
                            let h = bounds.height().max(1e-9);
                            // gate-display coords of the cursor (this gate's own transforms)
                            let to_gate = |p: PlotPoint| [gxt.forward(xt.inverse(p.x)), gyt.forward(yt.inverse(p.y))];
                            if r_started {
                                if let Some(p) = ptr {
                                    let mut best: Option<(usize, f64)> = None;
                                    for (i, hp) in handles.iter().enumerate() {
                                        let d = ((p.x - hp[0]) / w).powi(2) + ((p.y - hp[1]) / h).powi(2);
                                        if d < 0.0025 && best.map(|(_, bd)| d < bd).unwrap_or(true) { best = Some((i, d)); }
                                    }
                                    next_grab = best.map(|(i, _)| i);
                                    // No handle grabbed but cursor is inside the gate → drag the whole body.
                                    if next_grab.is_none() {
                                        let g = to_gate(p);
                                        next_move_last = if shape.contains(g[0], g[1]) { Some(g) } else { None };
                                    } else {
                                        next_move_last = None;
                                    }
                                    if next_grab.is_some() || next_move_last.is_some() { begin_edit = true; }
                                }
                            }
                            if r_dragged {
                                if let (Some(p), Some(hg)) = (ptr, next_grab) {
                                    let g = to_gate(p);
                                    handle_update = Some((*gidx, hg, g[0], g[1]));
                                } else if let (Some(p), Some(last)) = (ptr, next_move_last) {
                                    let g = to_gate(p);
                                    gate_translate = Some((*gidx, g[0] - last[0], g[1] - last[1]));
                                    next_move_last = Some(g);
                                }
                            }
                            if r_stopped { next_grab = None; next_move_last = None; }
                        }
                    }
                    DrawMode::Quadrant => {
                        // crosshair preview at cursor; click sets the split point
                        if let Some(p) = ptr {
                            let (x0, x1, y0, y1) = (bounds.min()[0], bounds.max()[0], bounds.min()[1], bounds.max()[1]);
                            pu.line(Line::new(PlotPoints::new(vec![[p.x, y0], [p.x, y1]]))
                                .color(Color32::from_rgb(240, 190, 40)).width(1.2));
                            pu.line(Line::new(PlotPoints::new(vec![[x0, p.y], [x1, p.y]]))
                                .color(Color32::from_rgb(240, 190, 40)).width(1.2));
                        }
                        if r_clicked {
                            if let Some(p) = ptr { new_quadrant = Some([p.x, p.y]); }
                        }
                    }
                    DrawMode::Navigate => {}
                }
            }
        });

        // commit interaction
        self.drag_start = next_ds;
        self.drag_current = next_dc;
        if let Some(v) = poly_add { self.poly_vertices.push(v); }
        if poly_finish { self.finish_polygon(); }
        if let Some(shape) = new_shape { self.commit_gate(shape); }
        if let Some(c) = new_quadrant { self.commit_quadrant(c[0], c[1]); self.draw_mode = DrawMode::Navigate; }
        if exit_draw { self.draw_mode = DrawMode::Navigate; }
        if begin_edit { self.push_undo(); }
        self.grab_handle = next_grab;
        self.gate_move_last = next_move_last;
        if let Some((gidx, h, gx, gy)) = handle_update {
            if gidx < self.gates.len() {
                apply_gate_handle(&mut self.gates[gidx].shape, h, gx, gy);
                self.needs_regate = true;
                if self.active_pop.is_some() { self.scatter = None; }
            }
        }
        if let Some((gidx, dx, dy)) = gate_translate {
            if gidx < self.gates.len() {
                translate_shape(&mut self.gates[gidx].shape, dx, dy);
                self.needs_regate = true;
                if self.active_pop.is_some() { self.scatter = None; }
            }
        }
        if let Some(p) = dbl_drill { self.drill_at(&x_name, &y_name, &xt, &yt, p); }

        // Live cursor readout (display → data units), shown in the status bar.
        self.cursor_label = hover_disp.map(|d| {
            let (xb, yb) = (x_name_base(&x_name), x_name_base(&y_name));
            format!("{} {} · {} {}",
                short_chan(&xb), fmt_data_tick(xt.inverse(d[0])),
                short_chan(&yb), fmt_data_tick(yt.inverse(d[1])))
        });

        // Remember the plot's screen rect so a screenshot can be cropped to it.
        self.last_plot_rect = Some(plot_response.response.rect);
    }

    // ── Histogram ─────────────────────────────────────────────────────

    fn histogram_view(&mut self, ui: &mut egui::Ui) {
        if self.fcs.is_none() {
            ui.centered_and_justified(|ui| {
                ui.label(RichText::new("Load a file to display histograms.").color(Color32::GRAY));
            });
            return;
        }
        let n_params = self.fcs.as_ref().map(|f| f.n_params()).unwrap_or(1).max(1);
        if self.hist_ch >= n_params { self.hist_ch = self.x_ch.min(n_params - 1); }
        let xi = self.hist_ch.min(n_params - 1);
        let cur_xt = self.cur_tf(xi);
        let x_name = self.fcs.as_ref().map(|f| param_full_label(f, xi)).unwrap_or_default();

        let stale = self.hist_cache.as_ref().map(|h| {
            h.x_ch != xi || h.x_label != cur_xt.short_label() || h.norm != self.hist_norm
                || h.mode != self.hist_mode || h.sample_pop != self.hist_sample_pop
        }).unwrap_or(true);
        if stale { self.rebuild_histogram(); }

        // channel names for the picker
        let ch_names: Vec<String> = self.fcs.as_ref().map(|f|
            (0..f.n_params()).map(|i| param_full_label(f, i)).collect()).unwrap_or_default();

        // controls
        ui.horizontal(|ui| {
            ui.label("Channel:");
            egui::ComboBox::from_id_salt("histch").selected_text(&x_name).width(180.0).show_ui(ui, |ui| {
                for (i, nm) in ch_names.iter().enumerate() {
                    if ui.selectable_label(self.hist_ch == i, nm).clicked() {
                        self.hist_ch = i; self.hist_cache = None;
                    }
                }
            });
            ui.separator();
            ui.label("Y:");
            if ui.selectable_value(&mut self.hist_norm, HistNorm::Modal, "Modal %").clicked() { self.hist_cache = None; }
            if ui.selectable_value(&mut self.hist_norm, HistNorm::Count, "Count").clicked() { self.hist_cache = None; }
            ui.separator();
            let drawing = self.hist_draw_interval;
            let lbl = if drawing {
                RichText::new(format!("{} Interval — drag on plot", icon::PENCIL_SIMPLE)).color(Color32::from_rgb(220, 170, 0))
            } else {
                RichText::new("+ Interval gate")
            };
            if ui.selectable_label(drawing, lbl).clicked() {
                self.hist_draw_interval = !drawing;
                self.drag_start = None; self.drag_current = None;
            }
            ui.separator();
            if ui.button(format!("{} Save plot…", icon::CAMERA)).clicked() { self.request_plot_png(); }
        });

        // overlay mode: populations (1 sample) vs samples (workspace)
        ui.horizontal(|ui| {
            ui.label("Overlay:");
            if ui.selectable_value(&mut self.hist_mode, HistMode::Populations, "Populations").clicked() { self.hist_cache = None; }
            let samples_enabled = self.samples.len() > 1;
            if ui.add_enabled(samples_enabled,
                egui::SelectableLabel::new(self.hist_mode == HistMode::Samples, "Samples")).clicked()
            {
                self.hist_mode = HistMode::Samples; self.hist_cache = None;
            }
            if !samples_enabled {
                ui.label(RichText::new("(open >1 sample for sample overlay)").small().color(Color32::GRAY));
            }
        });

        if self.hist_mode == HistMode::Samples && self.samples.len() > 1 {
            // pick which population to histogram across samples
            let pop_name = match self.hist_sample_pop {
                None => "All events".to_string(),
                Some(gid) => self.gates.iter().find(|g| g.id == gid).map(|g| g.name.clone()).unwrap_or_else(|| "All events".into()),
            };
            let gate_opts: Vec<(u32, String)> = self.gates.iter().map(|g| (g.id, g.name.clone())).collect();
            ui.horizontal(|ui| {
                ui.label("Population:");
                egui::ComboBox::from_id_salt("histpop").selected_text(pop_name).show_ui(ui, |ui| {
                    if ui.selectable_label(self.hist_sample_pop.is_none(), "All events").clicked() {
                        self.hist_sample_pop = None; self.hist_cache = None;
                    }
                    for (id, name) in &gate_opts {
                        if ui.selectable_label(self.hist_sample_pop == Some(*id), name).clicked() {
                            self.hist_sample_pop = Some(*id); self.hist_cache = None;
                        }
                    }
                });
            });
        } else {
            // population visibility toggles (populations mode)
            let gate_list: Vec<(u32, String, Color32)> = self.gates.iter()
                .map(|g| (g.id, g.name.clone(), gate_color(g.id).0)).collect();
            ui.horizontal_wrapped(|ui| {
                ui.label("Show:");
                let mut all_vis = !self.hist_all_hidden;
                if ui.checkbox(&mut all_vis, "All events").changed() {
                    self.hist_all_hidden = !all_vis; self.hist_cache = None;
                }
                for (id, name, color) in &gate_list {
                    let mut vis = !self.hist_hidden.contains(id);
                    if ui.checkbox(&mut vis, RichText::new(name).color(*color)).changed() {
                        if vis { self.hist_hidden.remove(id); } else { self.hist_hidden.insert(*id); }
                        self.hist_cache = None;
                    }
                }
            });
        }
        ui.separator();

        if self.hist_cache.is_none() { self.rebuild_histogram(); }
        let data = match &self.hist_cache { Some(d) => d, None => return };

        // clone render data out before the plot closure
        let series: Vec<(String, Color32, Vec<[f64; 2]>)> = data.series.iter().map(|s| {
            let pts: Vec<[f64; 2]> = data.centers.iter().zip(s.values.iter())
                .map(|(&x, &y)| [x, y]).collect();
            (s.name.clone(), s.color, pts)
        }).collect();

        let xt = cur_xt.compile();
        let x_base = x_name_base(&x_name);
        let range_gates: Vec<(String, Color32, f64, f64)> = self.gates.iter()
            .filter_map(|g| {
                if let GateShape::Range { x_min, x_max } = &g.shape {
                    if g.x_channel.eq_ignore_ascii_case(&x_base) {
                        let gxt = g.x_transform.compile();
                        return Some((g.name.clone(), gate_color(g.id).0,
                                     xt.forward(gxt.inverse(*x_min)), xt.forward(gxt.inverse(*x_max))));
                    }
                }
                None
            }).collect();

        let drawing = self.hist_draw_interval;
        let lock_view = self.lock_view;
        let view_bounds: Option<PlotBounds> = if lock_view {
            let x0 = data.centers.first().copied().unwrap_or(0.0);
            let x1 = data.centers.last().copied().unwrap_or(1.0);
            let mx = (x1 - x0).abs().max(1e-9) * 0.02;
            let ymax_d = match data.norm {
                HistNorm::Modal => 1.0,
                HistNorm::Count => data.series.iter().flat_map(|s| s.values.iter())
                    .cloned().fold(0.0_f64, f64::max).max(1.0),
            };
            Some(PlotBounds::from_min_max([x0 - mx, -0.02], [x1 + mx, ymax_d * 1.08]))
        } else { None };
        let cur_ds = self.drag_start;
        let cur_dc = self.drag_current;
        let mut next_ds = cur_ds;
        let mut next_dc = cur_dc;
        let mut new_interval: Option<(f64, f64)> = None;

        let xt_fmt = xt.clone();
        let xt_grid = xt.clone();
        let lin_x = matches!(cur_xt, AxisTransform::Linear);
        let norm_tag = match self.hist_norm { HistNorm::Modal => "m", HistNorm::Count => "c" };
        let y_label = match self.hist_norm { HistNorm::Modal => "Modal (peak = 1)", HistNorm::Count => "Count" };

        let plot = Plot::new(format!("hist_{}_{}_{}", xi, cur_xt.short_label(), norm_tag))
            .legend(Legend::default())
            .x_axis_label(&x_name)
            .y_axis_label(y_label)
            .allow_drag(!drawing && !lock_view).allow_zoom(!drawing && !lock_view).allow_scroll(false)
            .x_axis_formatter(move |gm: GridMark, _r| fmt_data_tick(xt_fmt.inverse(gm.value)))
            .x_grid_spacer(move |inp| nonlinear_grid(&xt_grid, lin_x, inp));

        let hist_response = plot.show(ui, |pu| {
            if let Some(b) = view_bounds {
                pu.set_plot_bounds(b);
                pu.set_auto_bounds(egui::Vec2b::new(false, false));
            }
            let bounds = pu.plot_bounds();
            let (ymin, ymax) = (bounds.min()[1], bounds.max()[1]);

            for (name, color, pts) in &series {
                pu.line(Line::new(PlotPoints::new(pts.clone())).color(*color).width(1.6).name(name.as_str()));
            }
            for (name, color, a, b) in &range_gates {
                pu.line(Line::new(PlotPoints::new(vec![[*a, ymin], [*a, ymax]])).color(*color).width(1.2));
                pu.line(Line::new(PlotPoints::new(vec![[*b, ymin], [*b, ymax]])).color(*color).width(1.2));
                pu.text(PlotText::new(PlotPoint::new((a + b) / 2.0, ymax), RichText::new(name).color(*color).size(11.0)));
            }

            if drawing {
                let ptr = pu.pointer_coordinate();
                let (started, dragged, stopped) = {
                    let r = pu.response();
                    (r.drag_started(), r.dragged(), r.drag_stopped())
                };
                if let (Some(s), Some(c)) = (cur_ds, cur_dc) {
                    let (lo, hi) = (s[0].min(c[0]), s[0].max(c[0]));
                    pu.polygon(PlotPolygon::new(PlotPoints::new(vec![[lo, ymin], [hi, ymin], [hi, ymax], [lo, ymax]]))
                        .stroke(Stroke::new(1.2, Color32::from_rgb(240, 190, 40)))
                        .fill_color(Color32::from_rgba_unmultiplied(240, 190, 40, 30)));
                }
                if started { next_ds = ptr.map(|p| [p.x, 0.0]); next_dc = next_ds; }
                if dragged { if let Some(p) = ptr { next_dc = Some([p.x, 0.0]); } }
                if stopped {
                    let end = ptr.map(|p| [p.x, 0.0]).or(next_dc);
                    if let (Some(s), Some(e)) = (next_ds, end) {
                        if (e[0] - s[0]).abs() > 1e-9 { new_interval = Some((s[0], e[0])); }
                    }
                    next_ds = None; next_dc = None;
                }
            }
        });

        self.drag_start = next_ds;
        self.drag_current = next_dc;
        self.last_plot_rect = Some(hist_response.response.rect);
        if let Some((a, b)) = new_interval {
            self.commit_range_gate(a, b);
            self.hist_draw_interval = false;
        }
    }

    // ── Stats table ───────────────────────────────────────────────────

    fn stats_table(&mut self, ui: &mut egui::Ui) {
        if self.fcs.is_none() {
            ui.centered_and_justified(|ui| {
                ui.label(RichText::new("Load a file to see statistics.").color(Color32::GRAY));
            });
            return;
        }

        // Stat channels = all except Time. Compute + cache lazily.
        let stat_channels: Vec<usize> = match &self.fcs {
            Some(fcs) => fcs.parameters.iter().enumerate()
                .filter(|(_, p)| !p.name.eq_ignore_ascii_case("Time"))
                .map(|(i, _)| i).collect(),
            None => return,
        };
        if self.pop_stats.is_none() {
            if let Some(fcs) = &self.fcs {
                // Match the buffer-length guard every other rebuild path uses, so a
                // short `compensated` (e.g. mid-load) can't index out of bounds.
                if self.compensated.len() >= fcs.n_events * fcs.n_params() {
                    self.pop_stats = Some(population_stats(
                        &self.compensated, &fcs.parameters, fcs.n_events, &self.gates, &stat_channels,
                    ));
                }
            }
        }
        let table = match &self.pop_stats { Some(t) => t, None => return };

        let stage = if self.do_compensate { "compensated" } else { "raw" };
        ui.horizontal(|ui| {
            ui.heading("Population statistics");
            ui.label(RichText::new(format!("· MFI = median on {} (linear) data", stage))
                .small().color(Color32::GRAY));
        });
        let mut do_export = false;
        let mut do_copy = false;
        ui.horizontal(|ui| {
            if ui.button(format!("{} Export CSV (tidy)", icon::EXPORT)).clicked() { do_export = true; }
            if ui.button(format!("{} Copy (TSV)", icon::COPY)).on_hover_text("Copy the table as tab-separated text (paste into R/Excel)").clicked() { do_copy = true; }
            ui.label(RichText::new(format!("{} populations × {} channels", table.rows.len(), table.channels.len()))
                .small().color(Color32::GRAY));
        });
        if do_copy {
            let mut s = String::from("Population\tCount\t%Parent\t%Total");
            for ch in &table.channels { s.push_str(&format!("\tMFI {}", ch)); }
            s.push('\n');
            for r in &table.rows {
                s.push_str(&format!("{}\t{}\t{:.2}\t{:.2}", r.name, r.count, r.pct_parent, r.pct_total));
                for m in &r.medians { s.push_str(&format!("\t{}", fmt(*m))); }
                s.push('\n');
            }
            ui.output_mut(|o| o.copied_text = s);
            self.status = "Copied population stats to clipboard (TSV).".into();
        }
        ui.separator();

        egui::ScrollArea::both().show(ui, |ui| {
            egui::Grid::new("popstats_grid").striped(true).spacing([14.0, 4.0]).show(ui, |ui| {
                ui.label(RichText::new("Population").strong());
                ui.label(RichText::new("Count").strong());
                ui.label(RichText::new("%Parent").strong());
                ui.label(RichText::new("%Total").strong());
                for ch in &table.channels {
                    ui.label(RichText::new(format!("MFI {}", ch)).strong());
                }
                ui.end_row();

                for r in &table.rows {
                    let indent = "    ".repeat(r.depth);
                    ui.label(RichText::new(format!("{}{}", indent, r.name)).monospace());
                    ui.label(r.count.to_string());
                    ui.label(format!("{:.2}%", r.pct_parent));
                    ui.label(format!("{:.2}%", r.pct_total));
                    for &m in &r.medians {
                        ui.label(fmt(m));
                    }
                    ui.end_row();
                }
            });
        });

        if do_export {
            self.export_popstats_csv();
        }
    }

    fn export_popstats_csv(&mut self) {
        let table = match &self.pop_stats { Some(t) => t, None => return };
        let default = self.file_path.as_ref()
            .and_then(|p| p.file_stem())
            .map(|s| format!("{}_popstats.csv", s.to_string_lossy()))
            .unwrap_or_else(|| "popstats.csv".into());
        let sample = self.file_path.as_ref()
            .and_then(|p| p.file_stem()).map(|s| s.to_string_lossy().to_string())
            .unwrap_or_else(|| "sample".into());
        if let Some(path) = rfd::FileDialog::new()
            .add_filter("CSV", &["csv"]).set_file_name(default).save_file()
        {
            let mut s = String::new();
            s.push_str(LONG_CSV_HEADER);
            s.push('\n');
            append_long_csv(&mut s, &sample, table);
            self.status = match std::fs::write(&path, s) {
                Ok(_) => format!("Exported population stats → {}", path.display()),
                Err(e) => format!("Export error: {}", e),
            };
        }
    }

    // ── Plot PNG export ───────────────────────────────────────────────

    /// Request a full-viewport screenshot; the captured frame is cropped to the
    /// last plot's rect and saved as PNG once it arrives (see `poll_screenshot`).
    fn request_plot_png(&mut self) {
        if self.last_plot_rect.is_none() {
            self.status = "Open the Plot or Histogram tab first.".into();
            return;
        }
        self.screenshot_pending = true;
        self.screenshot_sent = false;
    }

    /// Drive the async screenshot: dispatch the command once, then catch the
    /// delivered image on a later frame and save it.
    fn poll_screenshot(&mut self, ctx: &egui::Context) {
        if !self.screenshot_pending { return; }
        let shot = ctx.input(|i| i.events.iter().find_map(|e| {
            if let egui::Event::Screenshot { image, .. } = e { Some(image.clone()) } else { None }
        }));
        if let Some(img) = shot {
            let ppp = ctx.pixels_per_point();
            self.screenshot_pending = false;
            self.screenshot_sent = false;
            self.save_screenshot_png(&img, ppp);
        } else if !self.screenshot_sent {
            ctx.send_viewport_cmd(egui::ViewportCommand::Screenshot);
            self.screenshot_sent = true;
            ctx.request_repaint();
        }
    }

    /// Crop a captured viewport image to the plot rect and write it as a PNG.
    fn save_screenshot_png(&mut self, image: &egui::ColorImage, ppp: f32) {
        let rect = match self.last_plot_rect { Some(r) => r, None => return };
        let [iw, ih] = image.size;
        // Plot rect is in points; the captured image is in physical pixels.
        let x0 = ((rect.min.x * ppp).floor() as usize).min(iw);
        let y0 = ((rect.min.y * ppp).floor() as usize).min(ih);
        let x1 = ((rect.max.x * ppp).ceil() as usize).min(iw);
        let y1 = ((rect.max.y * ppp).ceil() as usize).min(ih);
        let (cw, ch) = (x1.saturating_sub(x0), y1.saturating_sub(y0));
        if cw == 0 || ch == 0 { self.status = "Plot rect empty — nothing to save.".into(); return; }

        let mut buf: Vec<u8> = Vec::with_capacity(cw * ch * 4);
        for y in y0..y1 {
            for x in x0..x1 {
                let p = image.pixels[y * iw + x];
                buf.extend_from_slice(&[p.r(), p.g(), p.b(), p.a()]);
            }
        }
        let default = self.file_path.as_ref().and_then(|p| p.file_stem())
            .map(|s| format!("{}_plot.png", s.to_string_lossy()))
            .unwrap_or_else(|| "plot.png".into());
        if let Some(path) = rfd::FileDialog::new()
            .add_filter("PNG", &["png"]).set_file_name(default).save_file()
        {
            let res = image::RgbaImage::from_raw(cw as u32, ch as u32, buf)
                .ok_or_else(|| "buffer size mismatch".to_string())
                .and_then(|img| img.save(&path).map_err(|e| e.to_string()));
            self.status = match res {
                Ok(_) => format!("Saved plot → {}", path.display()),
                Err(e) => format!("PNG save error: {}", e),
            };
        }
    }

    // ── QC suite ──────────────────────────────────────────────────────

    /// Stream an acquisition-QC pass over every workspace sample (open → metrics →
    /// drop, on a worker thread, same pattern as the batch). Per tube: event count,
    /// %viable (the designated live gate's % of parent, via the validated gating
    /// engine), flow-rate stability (clog detection on the Time channel), and
    /// margin/saturation events.
    fn run_qc(&mut self) {
        if self.samples.is_empty() {
            self.status = "No samples to QC.".into();
            return;
        }
        self.cancel_qc();
        let samples: Vec<(PathBuf, String, String)> = self.samples.iter()
            .map(|s| (s.path.clone(), s.name.clone(), s.group.clone())).collect();
        let total = samples.len();
        let gates = self.gates.clone();
        let do_comp = self.do_compensate;
        let ov = self.spill_override.clone();
        let live = self.qc_live_gate;
        let cancel = Arc::new(AtomicBool::new(false));
        let cw = cancel.clone();
        let (tx, rx) = std::sync::mpsc::channel::<QcMsg>();
        std::thread::spawn(move || {
            for (i, (path, name, group)) in samples.into_iter().enumerate() {
                if cw.load(Ordering::Relaxed) { break; }
                let _ = tx.send(QcMsg::Progress { done: i, total, name: name.clone() });
                let fcs = match FcsFile::open(&path) {
                    Ok(f) => f,
                    Err(e) => { let _ = tx.send(QcMsg::Skip(name, e.to_string())); continue; }
                };
                // Margin/saturation on RAW data (the acquisition range).
                let worst_margin = qc::worst_margin(&qc::margin_events(&fcs.events, &fcs.parameters, fcs.n_events));
                // Flow-rate on the Time channel (None if absent/constant).
                let flow = fcs.param_index("Time")
                    .map(|ti| fcs.channel_values(ti))
                    .and_then(|t| qc::flow_rate(&t, QC_FLOW_BINS, QC_FLOW_DEV));
                // %viable: the designated gate's effective count / its parent's.
                let viable = live.and_then(|gid| {
                    let events = compensate_for(&fcs, do_comp, ov.as_ref()).ok()?;
                    let own = compute_own_masks(&gates, &events, &fcs.parameters, fcs.n_events);
                    let by_id: HashMap<u32, &Gate> = gates.iter().map(|g| (g.id, g)).collect();
                    let g = by_id.get(&gid)?;
                    let count = |id: u32| effective_mask(id, &by_id, &own, fcs.n_events).iter().filter(|&&b| b).count();
                    let cnt = count(gid);
                    let par = match g.parent { Some(p) => count(p), None => fcs.n_events };
                    Some(if par > 0 { 100.0 * cnt as f64 / par as f64 } else { 0.0 })
                });
                let mut flags = Vec::new();
                if fcs.n_events < QC_MIN_EVENTS { flags.push(format!("low events ({})", fcs.n_events)); }
                if let Some(v) = viable { if v < QC_MIN_VIABLE { flags.push(format!("low viability ({:.0}%)", v)); } }
                if let Some(fr) = &flow { if fr.flagged { flags.push(format!("flow-rate Δ{:.0}%", fr.max_dev_pct)); } }
                if let Some((ch, m)) = &worst_margin { if *m > QC_MARGIN_PCT { flags.push(format!("{:.1}% off-scale on {}", m, ch)); } }
                let row = QcRow {
                    group, name, n_events: fcs.n_events, viable,
                    flow_dev: flow.as_ref().map(|f| f.max_dev_pct),
                    flow_bins: flow.map(|f| f.bins).unwrap_or_default(),
                    worst_margin, flags,
                };
                let _ = tx.send(QcMsg::Row(Box::new(row)));
            }
            let _ = tx.send(QcMsg::Done);
        });
        self.qc_rows.clear();
        self.qc_rx = Some(rx);
        self.qc_cancel = Some(cancel);
        self.qc_progress = Some((0, total));
        self.status = format!("QC scan started: 0/{}", total);
    }

    fn cancel_qc(&mut self) {
        if let Some(c) = &self.qc_cancel { c.store(true, Ordering::Relaxed); }
        self.qc_rx = None; self.qc_cancel = None; self.qc_progress = None;
    }

    /// Drain QC results each frame while a scan runs (mirrors poll_batch, incl. the
    /// disconnect-without-Done = worker-panic detection).
    fn poll_qc(&mut self, ctx: &egui::Context) {
        if self.qc_rx.is_none() { return; }
        let mut msgs: Vec<QcMsg> = Vec::new();
        let mut disconnected = false;
        if let Some(rx) = &self.qc_rx {
            loop {
                match rx.try_recv() {
                    Ok(m) => msgs.push(m),
                    Err(TryRecvError::Empty) => break,
                    Err(TryRecvError::Disconnected) => { disconnected = true; break; }
                }
            }
        }
        let mut got_done = false;
        for m in msgs {
            match m {
                QcMsg::Progress { done, total, name } => {
                    self.qc_progress = Some((done, total));
                    self.status = format!("QC {}/{}: {}", done, total, name);
                }
                QcMsg::Row(r) => self.qc_rows.push(*r),
                QcMsg::Skip(s, reason) => self.qc_rows.push(QcRow {
                    group: String::new(), name: s, n_events: 0, viable: None, flow_dev: None,
                    flow_bins: Vec::new(), worst_margin: None, flags: vec![format!("skipped: {}", reason)],
                }),
                QcMsg::Done => got_done = true,
            }
        }
        if got_done || disconnected {
            let n = self.qc_rows.len();
            let bad = self.qc_rows.iter().filter(|r| !r.flags.is_empty()).count();
            self.qc_rx = None; self.qc_cancel = None; self.qc_progress = None;
            self.status = if got_done {
                format!("QC done: {} tube(s), {} flagged", n, bad)
            } else {
                format!("{} QC ended early — a sample may have failed", icon::WARNING)
            };
        } else {
            ctx.request_repaint();
        }
    }

    fn qc_view(&mut self, ui: &mut egui::Ui) {
        ui.heading(format!("{} Acquisition QC", icon::HEARTBEAT));
        if self.samples.is_empty() {
            ui.label(RichText::new("Open one or more FCS files, then run a QC scan.").color(Color32::GRAY));
            return;
        }
        let running = self.qc_rx.is_some();
        ui.horizontal(|ui| {
            ui.label("Viability gate:");
            let cur = self.qc_live_gate
                .and_then(|id| self.gates.iter().find(|g| g.id == id))
                .map(|g| g.name.clone()).unwrap_or_else(|| "(none)".into());
            egui::ComboBox::from_id_salt("qc_live").selected_text(cur).show_ui(ui, |ui| {
                ui.selectable_value(&mut self.qc_live_gate, None, "(none)");
                let gates: Vec<(u32, String)> = self.gates.iter().map(|g| (g.id, g.name.clone())).collect();
                for (id, name) in gates { ui.selectable_value(&mut self.qc_live_gate, Some(id), name); }
            });
            ui.separator();
            if running {
                ui.add(egui::Spinner::new());
                if let Some((d, t)) = self.qc_progress { ui.label(format!("{}/{}", d, t)); }
                if ui.button(format!("{} Cancel", icon::X)).clicked() { self.cancel_qc(); }
            } else {
                if ui.button(format!("{} Run QC scan", icon::HEARTBEAT)).clicked() { self.run_qc(); }
                if !self.qc_rows.is_empty() && ui.button(format!("{} Export CSV", icon::EXPORT)).clicked() {
                    self.export_qc_csv();
                }
            }
        });
        ui.label(RichText::new(format!(
            "Flags: events <{} · viability <{:.0}% · flow-rate Δ>{:.0}% (clog) · off-scale >{:.0}%. \
             Pick a live/viable gate above for %viable.",
            QC_MIN_EVENTS, QC_MIN_VIABLE, QC_FLOW_DEV, QC_MARGIN_PCT)).small().weak());
        ui.separator();
        if self.qc_rows.is_empty() {
            ui.label(RichText::new("Run a scan to see per-tube quality.").color(Color32::GRAY));
            return;
        }
        let amber = Color32::from_rgb(220, 150, 50);
        let green = Color32::from_rgb(90, 180, 90);
        TableBuilder::new(ui).striped(true)
            .column(Column::auto().at_least(150.0))
            .column(Column::auto().at_least(64.0))
            .column(Column::auto().at_least(64.0))
            .column(Column::auto().at_least(74.0))
            .column(Column::auto().at_least(80.0))
            .column(Column::remainder().at_least(160.0))
            .header(20.0, |mut h| {
                for t in ["Sample", "Events", "%Viable", "Flow", "Off-scale", "Status"] {
                    h.col(|ui| { ui.strong(t); });
                }
            })
            .body(|mut body| {
                for r in &self.qc_rows {
                    body.row(20.0, |mut row| {
                        let bad = !r.flags.is_empty();
                        row.col(|ui| { ui.label(&r.name); });
                        row.col(|ui| {
                            let t = RichText::new(fmt_count(r.n_events));
                            ui.label(if r.n_events < QC_MIN_EVENTS { t.color(amber) } else { t });
                        });
                        row.col(|ui| {
                            ui.label(r.viable.map(|v| format!("{:.1}%", v)).unwrap_or_else(|| "—".into()));
                        });
                        row.col(|ui| {
                            if r.flow_bins.is_empty() {
                                ui.label(RichText::new("n/a").weak()).on_hover_text("no usable Time channel");
                            } else {
                                let flagged = r.flow_dev.map(|d| d > QC_FLOW_DEV).unwrap_or(false);
                                flow_sparkline(ui, &r.flow_bins, flagged).on_hover_text(format!(
                                    "flow-rate Δ {:.0}% over {} time bins (a dip = clog)",
                                    r.flow_dev.unwrap_or(0.0), r.flow_bins.len()));
                            }
                        });
                        row.col(|ui| {
                            ui.label(r.worst_margin.as_ref().map(|(_, m)| format!("{:.1}%", m)).unwrap_or_else(|| "—".into()));
                        });
                        row.col(|ui| {
                            if bad {
                                ui.label(RichText::new(format!("{} {}", icon::WARNING, r.flags.join("; "))).color(amber));
                            } else {
                                ui.label(RichText::new(format!("{} pass", icon::CHECK)).color(green));
                            }
                        });
                    });
                }
            });
    }

    fn export_qc_csv(&mut self) {
        let Some(path) = rfd::FileDialog::new().add_filter("CSV", &["csv"])
            .set_file_name("qc.csv").save_file() else { return; };
        let esc = |x: &str| if x.contains(',') || x.contains('"') {
            format!("\"{}\"", x.replace('"', "\"\""))
        } else { x.to_string() };
        let mut s = String::from("group,sample,events,pct_viable,flow_dev_pct,offscale_pct,offscale_channel,flags\n");
        for r in &self.qc_rows {
            let v = r.viable.map(|v| format!("{:.2}", v)).unwrap_or_default();
            let fd = r.flow_dev.map(|d| format!("{:.2}", d)).unwrap_or_default();
            let (mc, mp) = r.worst_margin.clone()
                .map(|(c, m)| (c, format!("{:.2}", m))).unwrap_or_default();
            s.push_str(&format!("{},{},{},{},{},{},{},{}\n",
                esc(&r.group), esc(&r.name), r.n_events, v, fd, mp, esc(&mc), esc(&r.flags.join("; "))));
        }
        match std::fs::write(&path, s) {
            Ok(_) => self.status = format!("Wrote QC table → {}", path.display()),
            Err(e) => self.status = format!("QC export error: {}", e),
        }
    }

    // ── Batch ─────────────────────────────────────────────────────────

    /// Start a background batch over every workspace sample. Each file is loaded,
    /// summarized, and dropped on a worker thread (memory stays flat), with results
    /// streamed back so the UI never freezes. Re-runs cancel any in-flight batch.
    fn run_batch(&mut self) {
        if self.samples.is_empty() {
            self.status = "No samples in the workspace.".into();
            return;
        }
        self.cancel_batch();

        let samples: Vec<(PathBuf, String, String)> = self.samples.iter()
            .map(|s| (s.path.clone(), s.name.clone(), s.group.clone())).collect();
        let total = samples.len();
        let gates = self.gates.clone();
        let do_comp = self.do_compensate;
        let ov = self.spill_override.clone();

        let cancel = Arc::new(AtomicBool::new(false));
        let cancel_worker = cancel.clone();
        let (tx, rx) = std::sync::mpsc::channel::<BatchMsg>();

        std::thread::spawn(move || {
            for (i, (path, name, group)) in samples.into_iter().enumerate() {
                if cancel_worker.load(Ordering::Relaxed) { break; }
                let _ = tx.send(BatchMsg::Progress { done: i, total, name: name.clone() });
                let fcs = match FcsFile::open(&path) {
                    Ok(f) => f,
                    Err(e) => { let _ = tx.send(BatchMsg::Skip(name, e.to_string())); continue; }
                };
                let missing = missing_gate_channels(&gates, &fcs);
                if !missing.is_empty() {
                    let _ = tx.send(BatchMsg::Skip(name, format!("missing channels: {}", missing.join(", "))));
                    continue;
                }
                let events = match compensate_for(&fcs, do_comp, ov.as_ref()) {
                    Ok(e) => e,
                    Err(e) => { let _ = tx.send(BatchMsg::Skip(name, format!("compensation failed: {}", e))); continue; }
                };
                let stat_channels: Vec<usize> = fcs.parameters.iter().enumerate()
                    .filter(|(_, p)| !p.name.eq_ignore_ascii_case("Time"))
                    .map(|(i, _)| i).collect();
                let table = population_stats(&events, &fcs.parameters, fcs.n_events, &gates, &stat_channels);
                let _ = tx.send(BatchMsg::Table(group, name, table));
                // `fcs` and `events` drop here → flat memory
            }
            let _ = tx.send(BatchMsg::Done);
        });

        self.batch = Some(BatchResult { tables: Vec::new(), skipped: Vec::new() });
        self.batch_rx = Some(rx);
        self.batch_cancel = Some(cancel);
        self.batch_progress = Some((0, total));
        self.status = format!("Batch started: 0/{}", total);
    }

    /// Signal the worker to stop and detach from it.
    fn cancel_batch(&mut self) {
        if let Some(c) = &self.batch_cancel { c.store(true, Ordering::Relaxed); }
        self.batch_rx = None;
        self.batch_cancel = None;
        self.batch_progress = None;
    }

    /// Drain streamed batch results into `self.batch`; called every frame while running.
    fn poll_batch(&mut self, ctx: &egui::Context) {
        if self.batch_rx.is_none() { return; }
        let mut msgs: Vec<BatchMsg> = Vec::new();
        let mut disconnected = false;
        if let Some(rx) = &self.batch_rx {
            loop {
                match rx.try_recv() {
                    Ok(m) => msgs.push(m),
                    Err(TryRecvError::Empty) => break,
                    Err(TryRecvError::Disconnected) => { disconnected = true; break; }
                }
            }
        }
        let mut got_done = false;
        for m in msgs {
            match m {
                BatchMsg::Progress { done, total, name } => {
                    self.batch_progress = Some((done, total));
                    self.status = format!("Batch {}/{}: {}", done, total, name);
                }
                BatchMsg::Table(g, s, t) => { if let Some(b) = &mut self.batch { b.tables.push((g, s, t)); } }
                BatchMsg::Skip(s, r) => { if let Some(b) = &mut self.batch { b.skipped.push((s, r)); } }
                BatchMsg::Done => { got_done = true; }
            }
        }
        // The worker sends Done on the normal path; a disconnect *without* a prior Done
        // means it panicked mid-batch — surface that instead of reporting clean success.
        if got_done || disconnected {
            let n = self.batch.as_ref().map(|b| b.tables.len()).unwrap_or(0);
            let ns = self.batch.as_ref().map(|b| b.skipped.len()).unwrap_or(0);
            self.batch_rx = None;
            self.batch_cancel = None;
            self.batch_progress = None;
            self.status = if got_done {
                format!("Batch done: {} processed, {} skipped", n, ns)
            } else {
                format!("{} Batch ended early ({} processed, {} skipped) — a sample may have failed", icon::WARNING, n, ns)
            };
        } else {
            ctx.request_repaint(); // keep ticking while the worker runs
        }
    }

    fn export_batch_csv(&mut self) {
        let batch = match &self.batch { Some(b) => b, None => return };
        if let Some(path) = rfd::FileDialog::new()
            .add_filter("CSV", &["csv"]).set_file_name("batch_popstats.csv").save_file()
        {
            let mut s = String::new();
            s.push_str(LONG_CSV_HEADER_GROUPED);
            s.push('\n');
            for (group, sample, table) in &batch.tables {
                append_long_csv_grouped(&mut s, group, sample, table);
            }
            self.status = match std::fs::write(&path, s) {
                Ok(_) => format!("Exported batch ({} samples) → {}", batch.tables.len(), path.display()),
                Err(e) => format!("Export error: {}", e),
            };
        }
    }

    /// Export the tidy popstats CSV **plus** a companion starter `.R` script wired to
    /// the group/condition tags — so the analysis bridges straight into R.
    fn export_r_bundle(&mut self) {
        let (csv, groups, pops) = match &self.batch {
            Some(b) if !b.tables.is_empty() => {
                let mut csv = String::new();
                csv.push_str(LONG_CSV_HEADER_GROUPED);
                csv.push('\n');
                let mut groups: Vec<String> = Vec::new();
                for (g, s, t) in &b.tables {
                    append_long_csv_grouped(&mut csv, g, s, t);
                    if !g.is_empty() && !groups.contains(g) { groups.push(g.clone()); }
                }
                let pops: Vec<String> = b.tables.first()
                    .map(|(_, _, t)| t.rows.iter().map(|r| r.name.clone()).collect())
                    .unwrap_or_default();
                (csv, groups, pops)
            }
            _ => { self.status = "Run a batch first, then export the R bundle.".into(); return; }
        };
        let Some(path) = rfd::FileDialog::new().add_filter("CSV", &["csv"])
            .set_file_name("flowcyto_popstats.csv").save_file() else { return; };
        let csv_name = path.file_name().map(|s| s.to_string_lossy().to_string())
            .unwrap_or_else(|| "popstats.csv".into());
        let r_path = path.with_extension("R");
        let r = crate::r_bridge::r_script(&csv_name, &groups, &pops, env!("CARGO_PKG_VERSION"));
        self.status = match std::fs::write(&path, csv).and_then(|_| std::fs::write(&r_path, r)) {
            Ok(_) => format!("Wrote {} + {}", csv_name,
                r_path.file_name().map(|s| s.to_string_lossy().to_string()).unwrap_or_default()),
            Err(e) => format!("R bundle export error: {}", e),
        };
    }

    /// Gather the current analysis into a provenance bundle for the HTML report.
    fn build_report_data(&self) -> Option<crate::report::ReportData> {
        use crate::report::{ChannelRow, GateRow, ReportData};
        let fcs = self.fcs.as_ref()?;
        let np = fcs.n_params();
        let kwv = |k: &str| fcs.keywords.get(k).filter(|v| !v.is_empty()).cloned();
        let file_name = self.samples.get(self.active_sample)
            .map(|s| s.path.file_name().map(|x| x.to_string_lossy().to_string())
                .unwrap_or_else(|| format!("{}.fcs", s.name)))
            .unwrap_or_else(|| "sample.fcs".into());
        let cytometer = kwv("$CYT").or_else(|| kwv("$CYTOMETER")).unwrap_or_else(|| "—".into());
        let acq_date = kwv("$DATE").unwrap_or_else(|| "—".into());

        let channels: Vec<ChannelRow> = fcs.parameters.iter().enumerate().map(|(i, p)| ChannelRow {
            name: p.name.clone(),
            stain: p.label.clone().filter(|l| !l.is_empty()).unwrap_or_else(|| "—".into()),
            transform: describe_tf(&self.cur_tf(i)),
        }).collect();

        let (mut spill_channels, mut spill_rows, mut spill_source) = (Vec::new(), Vec::new(), String::new());
        if self.do_compensate {
            if let Some(ov) = &self.spill_override {
                spill_channels = ov.channels.clone();
                spill_rows = ov.rows.clone();
                spill_source = "manual override".into();
            } else if let Some(kw) = fcs.spillover_keyword() {
                if let Ok((ch, rows)) = crate::compensation::parse_spillover(kw) {
                    spill_channels = ch;
                    spill_rows = rows;
                    let creator = ["CREATOR", "$CREATOR", "APPLICATION"].iter().find_map(|&k| kwv(k));
                    let mut prov: Vec<String> = creator.into_iter().collect();
                    if let Some(c) = kwv("$CYT").or_else(|| kwv("$CYTOMETER")) { prov.push(c); }
                    if let Some(d) = kwv("$DATE") { prov.push(d); }
                    spill_source = if prov.is_empty() { "embedded $SPILLOVER".into() }
                        else { format!("embedded $SPILLOVER ({})", prov.join(" · ")) };
                }
            }
        }

        let mut gates = Vec::new();
        for (gid, depth) in crate::gating::gate_tree_order(&self.gates) {
            if let Some(g) = self.gates.iter().find(|x| x.id == gid) {
                let parent = g.parent.and_then(|p| self.gates.iter().find(|x| x.id == p))
                    .map(|x| x.name.clone()).unwrap_or_else(|| "All events".into());
                let channels = match &g.shape {
                    GateShape::Boolean { .. } => "—".into(),
                    GateShape::Range { .. } => format!("{} ({})", g.x_channel, g.x_transform.short_label()),
                    _ => format!("{} ({}) × {} ({})", g.x_channel, g.x_transform.short_label(),
                                 g.y_channel, g.y_transform.short_label()),
                };
                gates.push(GateRow { depth, name: g.name.clone(), parent, channels, shape: shape_desc(&g.shape) });
            }
        }

        let stats = if self.compensated.len() >= fcs.n_events * np {
            let stat_channels: Vec<usize> = fcs.parameters.iter().enumerate()
                .filter(|(_, p)| !p.name.eq_ignore_ascii_case("Time")).map(|(i, _)| i).collect();
            Some(population_stats(&self.compensated, &fcs.parameters, fcs.n_events, &self.gates, &stat_channels))
        } else { None };

        Some(ReportData {
            version: env!("CARGO_PKG_VERSION").into(),
            generated: crate::report::now_utc(),
            file_name, n_events: fcs.n_events, cytometer, acq_date,
            channels, compensated: self.do_compensate, spill_channels, spill_rows, spill_source,
            gates, stats,
            gates_json: serde_json::to_string_pretty(&self.gates).unwrap_or_default(),
        })
    }

    /// Write a self-contained HTML reproducibility report and open it in the browser.
    fn export_report(&mut self) {
        let data = match self.build_report_data() {
            Some(d) => d,
            None => { self.status = "Open a file first.".into(); return; }
        };
        let html = crate::report::html(&data);
        let Some(path) = rfd::FileDialog::new().add_filter("HTML", &["html"])
            .set_file_name("flowcyto_report.html").save_file() else { return; };
        self.status = match std::fs::write(&path, html) {
            Ok(_) => { let _ = open::that(&path); format!("Wrote report → {} (opened in browser)", path.display()) }
            Err(e) => format!("Report export error: {}", e),
        };
    }

    fn batch_view(&mut self, ui: &mut egui::Ui) {
        let mut do_run = false;
        let mut do_export = false;
        let mut do_copy = false;
        let mut do_cancel = false;
        let mut do_r = false;
        let running = self.batch_rx.is_some();
        ui.horizontal(|ui| {
            ui.heading("Batch");
            if ui.add_enabled(!running, egui::Button::new("▶ Run over all samples")).clicked() { do_run = true; }
            if running && ui.button(format!("{} Cancel", icon::X)).clicked() { do_cancel = true; }
            if self.batch.is_some() && !running && ui.button(format!("{} Export combined CSV", icon::EXPORT)).clicked() { do_export = true; }
            if self.batch.is_some() && !running && ui.button(format!("{} Copy (TSV)", icon::COPY)).on_hover_text("Copy the visible table as tab-separated text").clicked() { do_copy = true; }
            if self.batch.is_some() && !running && ui.button(format!("{} R bundle", icon::EXPORT)).on_hover_text("Export the tidy CSV + a starter .R analysis script wired to your group tags").clicked() { do_r = true; }
        });
        ui.label(RichText::new(format!(
            "{} samples · {} gates · streamed on a worker thread (UI stays responsive)", self.samples.len(), self.gates.len()
        )).small().color(Color32::GRAY));
        if let Some((done, total)) = self.batch_progress {
            let frac = if total > 0 { done as f32 / total as f32 } else { 0.0 };
            ui.add(egui::ProgressBar::new(frac).show_percentage()
                .text(format!("{}/{}", done, total)));
        }
        ui.separator();

        if do_run { self.run_batch(); }
        if do_cancel { self.cancel_batch(); self.status = "Batch cancelled.".into(); }
        if do_export { self.export_batch_csv(); }
        if do_r { self.export_r_bundle(); }

        // Clone display data out of the batch borrow. Each row carries its OWN
        // table's channel list so the MFI is looked up by NAME (panels may differ
        // in channel order across tubes), not by a shared positional index.
        #[allow(clippy::type_complexity)]
        let (channels, rows, skipped): (Vec<String>, Vec<(String, String, String, usize, usize, f64, f64, Vec<f64>, Vec<String>)>, Vec<(String, String)>) =
            match &self.batch {
                None => {
                    ui.label(RichText::new("Click “Run over all samples” to compute population stats across the workspace.")
                        .color(Color32::GRAY));
                    return;
                }
                Some(b) => {
                    // Union of channel names across all tables (first table's order, then extras).
                    let mut channels: Vec<String> = Vec::new();
                    for (_, _, t) in &b.tables {
                        for c in &t.channels {
                            if !channels.iter().any(|x| x.eq_ignore_ascii_case(c)) { channels.push(c.clone()); }
                        }
                    }
                    let mut rows = Vec::new();
                    for (group, sample, table) in &b.tables {
                        for r in &table.rows {
                            rows.push((group.clone(), sample.clone(), r.name.clone(), r.depth, r.count, r.pct_parent, r.pct_total, r.medians.clone(), table.channels.clone()));
                        }
                    }
                    (channels, rows, b.skipped.clone())
                }
            };

        if channels.is_empty() {
            ui.label(RichText::new("No samples produced stats (all skipped?).").color(Color32::from_rgb(220, 150, 50)));
        } else {
            let ci = self.batch_channel.min(channels.len() - 1);
            let sel_name = channels[ci].clone();
            if do_copy {
                let mut s = format!("Group\tSample\tPopulation\tCount\t%Parent\t%Total\tMFI {}\n", sel_name);
                for (group, sample, name, depth, count, pp, pt, medians, ch_names) in &rows {
                    let mfi = ch_names.iter().position(|c| c.eq_ignore_ascii_case(&sel_name))
                        .and_then(|k| medians.get(k).copied()).unwrap_or(f64::NAN);
                    s.push_str(&format!("{}\t{}\t{}{}\t{}\t{:.2}\t{:.2}\t{}\n",
                        group, sample, "  ".repeat(*depth), name, count, pp, pt, fmt(mfi)));
                }
                ui.output_mut(|o| o.copied_text = s);
                self.status = "Copied batch table to clipboard (TSV).".into();
            }
            ui.horizontal(|ui| {
                ui.label("MFI channel:");
                egui::ComboBox::from_id_salt("batchch").selected_text(&channels[ci]).show_ui(ui, |ui| {
                    for (i, c) in channels.iter().enumerate() {
                        ui.selectable_value(&mut self.batch_channel, i, c);
                    }
                });
            });
            let num = |ui: &mut egui::Ui, s: String| {
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| ui.monospace(s));
            };
            TableBuilder::new(ui)
                .striped(true)
                .resizable(true)
                .cell_layout(egui::Layout::left_to_right(egui::Align::Center))
                .column(Column::auto().at_least(80.0))                 // Group
                .column(Column::auto().at_least(120.0).clip(true))     // Sample
                .column(Column::auto().at_least(140.0).clip(true))     // Population
                .column(Column::auto().at_least(60.0))                 // Count
                .column(Column::auto().at_least(64.0))                 // %Parent
                .column(Column::auto().at_least(64.0))                 // %Total
                .column(Column::remainder().at_least(80.0))            // MFI
                .header(22.0, |mut h| {
                    for t in ["Group", "Sample", "Population", "Count", "%Parent", "%Total"] {
                        h.col(|ui| { ui.strong(t); });
                    }
                    h.col(|ui| { ui.strong(format!("MFI {}", sel_name)); });
                })
                .body(|body| {
                    body.rows(20.0, rows.len(), |mut row| {
                        let (group, sample, name, depth, count, pp, pt, medians, ch_names) = &rows[row.index()];
                        row.col(|ui| { ui.label(RichText::new(group).color(Color32::GRAY)); });
                        row.col(|ui| { ui.label(sample); });
                        row.col(|ui| { ui.label(format!("{}{}", "  ".repeat(*depth), name)); });
                        row.col(|ui| { num(ui, count.to_string()); });
                        row.col(|ui| { num(ui, format!("{:.2}%", pp)); });
                        row.col(|ui| { num(ui, format!("{:.2}%", pt)); });
                        row.col(|ui| {
                            let mfi = ch_names.iter().position(|c| c.eq_ignore_ascii_case(&sel_name))
                                .and_then(|k| medians.get(k).copied()).unwrap_or(f64::NAN);
                            num(ui, fmt(mfi));
                        });
                    });
                });
        }

        // ── Quick chart: one population's metric across samples, grouped by tag ──
        if !rows.is_empty() {
            ui.separator();
            egui::CollapsingHeader::new(format!("{} Chart across samples", icon::CHART_BAR)).id_salt("batch_chart_hdr").show(ui, |ui| {
                let cur_name = match self.batch_plot_pop {
                    None => "All events".to_string(),
                    Some(id) => self.gates.iter().find(|g| g.id == id).map(|g| g.name.clone()).unwrap_or_else(|| "All events".into()),
                };
                let gate_opts: Vec<(Option<u32>, String)> = std::iter::once((None, "All events".to_string()))
                    .chain(self.gates.iter().map(|g| (Some(g.id), g.name.clone()))).collect();
                let mfi_name = channels.get(self.batch_channel.min(channels.len().saturating_sub(1))).cloned().unwrap_or_default();
                ui.horizontal(|ui| {
                    ui.label("Population:");
                    egui::ComboBox::from_id_salt("batchplotpop").selected_text(&cur_name).show_ui(ui, |ui| {
                        for (id, name) in &gate_opts { ui.selectable_value(&mut self.batch_plot_pop, *id, name); }
                    });
                    ui.separator();
                    ui.label("Y:");
                    ui.selectable_value(&mut self.batch_plot_metric, BatchMetric::PctParent, "% parent");
                    ui.selectable_value(&mut self.batch_plot_metric, BatchMetric::PctTotal, "% total");
                    ui.selectable_value(&mut self.batch_plot_metric, BatchMetric::Mfi, "MFI");
                    if self.batch_plot_metric == BatchMetric::Mfi {
                        ui.label(RichText::new(format!("({})", mfi_name)).small().color(Color32::GRAY));
                    }
                });

                let metric = self.batch_plot_metric;
                let sel = cur_name.clone();
                let mut groups: Vec<String> = Vec::new();
                for (g, _s, name, _d, _c, _pp, _pt, _m, _cn) in &rows {
                    if name == &sel && !groups.contains(g) { groups.push(g.clone()); }
                }
                let mut series: Vec<Vec<[f64; 2]>> = vec![Vec::new(); groups.len()];
                for (g, _s, name, _d, _c, pp, pt, medians, ch_names) in &rows {
                    if name != &sel { continue; }
                    let val = match metric {
                        BatchMetric::PctParent => *pp,
                        BatchMetric::PctTotal => *pt,
                        BatchMetric::Mfi => ch_names.iter().position(|c| c.eq_ignore_ascii_case(&mfi_name))
                            .and_then(|k| medians.get(k).copied()).unwrap_or(f64::NAN),
                    };
                    if !val.is_finite() { continue; }
                    if let Some(gi) = groups.iter().position(|x| x == g) {
                        let k = series[gi].len();
                        let jitter = (k as f64 % 5.0) * 0.06 - 0.12; // spread points within a group
                        series[gi].push([gi as f64 + jitter, val]);
                    }
                }
                let ylab = metric.label();
                let gnames = groups.clone();
                if groups.is_empty() {
                    ui.label(RichText::new("Tag samples with a group (Samples panel) and pick a population.").small().color(Color32::GRAY));
                } else {
                    Plot::new("batch_chart").height(260.0).allow_scroll(false)
                        .legend(Legend::default())
                        .y_axis_label(ylab)
                        .x_axis_formatter(move |gm: GridMark, _r| {
                            let i = gm.value.round() as isize;
                            if i >= 0 && (i as usize) < gnames.len() && (gm.value - i as f64).abs() < 0.01 {
                                gnames[i as usize].clone()
                            } else { String::new() }
                        })
                        .show(ui, |pu| {
                            for (gi, pts) in series.iter().enumerate() {
                                pu.points(Points::new(PlotPoints::new(pts.clone()))
                                    .radius(4.0).color(sample_color(gi)).name(&groups[gi]));
                            }
                        });
                }
            });
        }

        if !skipped.is_empty() {
            ui.separator();
            ui.label(RichText::new("Skipped:").strong().color(Color32::from_rgb(220, 150, 50)));
            for (s, reason) in &skipped {
                ui.label(RichText::new(format!("  {} — {}", s, reason)).small());
            }
        }
    }

    fn ui_samples(&mut self, ui: &mut egui::Ui) {
        // Show the panel for any non-empty workspace so the group-tag editor and Clear
        // button are reachable even with a single sample (the switch/overlay controls
        // simply have nothing to act on yet).
        if self.samples.is_empty() {
            return;
        }
        let n_low = self.samples.iter().filter(|s| s.n_events.map(|n| n < QC_MIN_EVENTS).unwrap_or(false)).count();
        ui.heading(format!("Samples ({})", self.samples.len()));
        if n_low > 0 {
            ui.label(RichText::new(format!("{} {} low-event tube(s)", icon::WARNING, n_low))
                .small().color(Color32::from_rgb(220, 150, 50)));
        }
        let active = self.active_sample;
        let mut switch_to: Option<usize> = None;
        let mut clear = false;
        let mut overlay_pick: Option<Option<usize>> = None;
        egui::ScrollArea::vertical().max_height(220.0).id_salt("samples_scroll").show(ui, |ui| {
            for (i, s) in self.samples.iter_mut().enumerate() {
                let low = s.n_events.map(|n| n < QC_MIN_EVENTS).unwrap_or(false);
                let ev = s.n_events.map(fmt_count).unwrap_or_else(|| "?".into());
                let lbl = format!("{}{}  · {} ev", if low { format!("{} ", icon::WARNING) } else { String::new() }, s.name, ev);
                let txt = if low { RichText::new(lbl).color(Color32::from_rgb(220, 150, 50)) } else { RichText::new(lbl) };
                ui.horizontal(|ui| {
                    if ui.selectable_label(i == active, txt).clicked() { switch_to = Some(i); }
                    // reference-overlay toggle (👁) — overlay this sample behind the active one
                    let is_ref = self.ref_sample == Some(i);
                    if i != active && ui.selectable_label(is_ref, icon::EYE).on_hover_text("overlay behind active sample").clicked() {
                        overlay_pick = Some(if is_ref { None } else { Some(i) });
                    }
                });
                ui.horizontal(|ui| {
                    ui.add_space(14.0);
                    ui.label(RichText::new("group:").small().color(Color32::GRAY));
                    ui.add(egui::TextEdit::singleline(&mut s.group)
                        .desired_width(130.0).hint_text("(condition)"));
                });
            }
        });
        ui.horizontal(|ui| {
            if ui.small_button(format!("{} Clear", icon::X)).clicked() { clear = true; }
            if self.ref_sample.is_some() && ui.small_button(format!("{} overlay", icon::EYE_SLASH)).clicked() { overlay_pick = Some(None); }
        });
        ui.separator();

        if let Some(r) = overlay_pick { self.ref_sample = r; self.ref_scatter = None; }

        if let Some(i) = switch_to {
            if i != active {
                if self.ref_sample == Some(i) { self.ref_sample = None; self.ref_scatter = None; }
                self.activate_sample(i, false);
            }
        }
        if clear {
            self.samples.clear();
            self.fcs = None;
            self.compensated.clear();
            self.batch = None;
            self.active_sample = 0;
            self.ref_sample = None;
            self.ref_scatter = None;
            self.gates.clear();
            self.active_pop = None;
            self.selected_gate = None;
            self.status = "Workspace cleared.".into();
        }
    }

    // ── Spillover matrix viewer ───────────────────────────────────────

    fn spillover_view(&mut self, ui: &mut egui::Ui) {
        if self.fcs.is_none() {
            ui.centered_and_justified(|ui| {
                ui.label(RichText::new("Load a file to view its spillover matrix.").color(Color32::GRAY));
            });
            return;
        }
        let dark = self.dark_mode;

        // ── Owned extraction (drop the fcs borrow before mutating self) ──
        let (stain_map, prov_line, has_embedded) = {
            let fcs = self.fcs.as_ref().unwrap();
            let stain_map: std::collections::HashMap<String, String> = fcs.parameters.iter()
                .filter_map(|p| p.label.as_ref().filter(|l| !l.is_empty())
                    .map(|l| (p.name.to_uppercase(), format!(" ({})", l))))
                .collect();
            let get = |k: &str| fcs.keywords.get(k).map(String::as_str);
            let prov: Vec<String> = [
                get("CREATOR").or(get("$CREATOR")).or(get("APPLICATION")),
                get("$CYT").or(get("$CYTOMETER")),
                get("$DATE"),
            ].into_iter().flatten().map(String::from).collect();
            let prov_line = if prov.is_empty() { None } else { Some(prov.join(" · ")) };
            (stain_map, prov_line, fcs.spillover_keyword().is_some())
        };
        let stain = |name: &str| stain_map.get(&name.to_uppercase()).cloned().unwrap_or_default();

        let override_active = self.spill_override.is_some();
        let active = self.active_matrix(); // Option<(channels, rows)>

        ui.heading("Spillover (compensation) matrix");

        // Status line
        if override_active {
            ui.label(RichText::new("● OVERRIDE ACTIVE — using your custom matrix (not the embedded one)")
                .color(Color32::from_rgb(230, 140, 40)).strong());
            if !self.do_compensate {
                ui.label(RichText::new("Enable “Compensate” in the toolbar to apply it to the plot/stats.")
                    .small().color(Color32::from_rgb(220, 170, 60)));
            }
        } else if has_embedded {
            if let Some(p) = &prov_line {
                ui.label(RichText::new(format!("Embedded matrix · created at acquisition by: {}", p))
                    .small().color(Color32::GRAY));
            } else {
                ui.label(RichText::new("Embedded matrix from the file.").small().color(Color32::GRAY));
            }
        } else {
            ui.label(RichText::new("No embedded matrix in this file (uncompensated).")
                .color(Color32::from_rgb(210, 150, 60)));
        }

        // ── Toolbar ──────────────────────────────────────────────────
        let (mut act_edit, mut act_reset, mut act_load, mut act_save, mut act_write, mut act_compute) =
            (false, false, false, false, false, false);
        ui.horizontal(|ui| {
            if override_active {
                if ui.button(format!("{} Reset to embedded", icon::ARROW_COUNTER_CLOCKWISE)).clicked() { act_reset = true; }
            } else {
                let lbl = if has_embedded { format!("{} Edit / Override", icon::PENCIL_SIMPLE) } else { format!("{} Create override (identity)", icon::PENCIL_SIMPLE) };
                if ui.button(lbl).clicked() { act_edit = true; }
            }
            if ui.button(format!("{} Compute from controls…", icon::FLASK)).clicked() { act_compute = true; }
            if ui.button(format!("{} Load matrix…", icon::FOLDER_OPEN)).clicked() { act_load = true; }
            if active.is_some() && ui.button(format!("{} Save matrix…", icon::FLOPPY_DISK)).clicked() { act_save = true; }
            if active.is_some() && ui.button(format!("{} Write new .fcs…", icon::FILE_PLUS)).clicked() { act_write = true; }
        });
        if override_active {
            ui.label(RichText::new("Drag any cell to edit. Diagonal is normally 1.0.").small().color(Color32::GRAY));
        }
        ui.separator();

        // ── Matrix display / edit ────────────────────────────────────
        let mut edited_rows: Option<Vec<Vec<f64>>> = None;
        if let Some((channels, rows)) = &active {
            let n = channels.len();
            let mx = max_off_diagonal(rows);
            if mx < 1e-9 {
                ui.label(RichText::new(format!("{} Identity matrix — NO real compensation encoded.", icon::WARNING))
                    .color(Color32::from_rgb(220, 150, 50)).strong());
            } else {
                ui.label(RichText::new(format!(
                    "Max off-diagonal spillover: {:.4} ({:.1}%)", mx, mx * 100.0
                )).color(Color32::from_rgb(80, 180, 80)));
            }
            ui.add_space(4.0);

            let mut local = rows.clone();
            let mut changed = false;
            egui::ScrollArea::both().show(ui, |ui| {
                egui::Grid::new("spill_grid").striped(false).spacing([4.0, 3.0]).show(ui, |ui| {
                    ui.label("");
                    for c in channels {
                        ui.label(RichText::new(c).strong().monospace());
                    }
                    ui.end_row();
                    for i in 0..n {
                        ui.label(RichText::new(format!("{}{}", channels[i], stain(&channels[i])))
                            .strong().monospace());
                        for j in 0..n {
                            if override_active {
                                let mut v = local[i][j];
                                if ui.add(egui::DragValue::new(&mut v).speed(0.001)
                                    .range(-5.0..=5.0).fixed_decimals(4)).changed()
                                {
                                    local[i][j] = v;
                                    changed = true;
                                }
                            } else {
                                let v = local[i][j];
                                let is_diag = i == j;
                                let txt = if is_diag { format!("{:>6}", "1.000") } else { format!("{:>6.3}", v) };
                                ui.label(RichText::new(txt).monospace()
                                    .background_color(spill_cell_color(v, is_diag, dark)));
                            }
                        }
                        ui.end_row();
                    }
                });
            });
            if changed { edited_rows = Some(local); }
        }

        ui.add_space(10.0);
        ui.label(RichText::new(
            "flowcyto applies the inverse of the active matrix to the fluorescence channels \
             (compensated = raw × M⁻¹). “Write new .fcs” saves a fresh file with the original raw \
             events + this matrix — your original file is never modified."
        ).small().color(Color32::GRAY));

        // ── Apply actions (self mutations after rendering) ───────────
        if let Some(er) = edited_rows {
            if let Some(ov) = &mut self.spill_override { ov.rows = er; }
            self.needs_reprocess = true;
        }
        if act_edit { self.start_override(); self.needs_reprocess = true; }
        if act_reset { self.spill_override = None; self.needs_reprocess = true; }
        if act_load { self.load_override(); }
        if act_compute { self.compute_from_controls(); }
        if act_save || act_write {
            if let Some((ch, rows)) = self.active_matrix() {
                if act_save { self.save_active_matrix(&ch, &rows); }
                if act_write { self.write_fcs_override(&ch, &rows); }
            }
        }
    }

    /// The matrix currently in effect: the override if set, else the embedded one.
    fn active_matrix(&self) -> Option<(Vec<String>, Vec<Vec<f64>>)> {
        if let Some(ov) = &self.spill_override {
            return Some((ov.channels.clone(), ov.rows.clone()));
        }
        let fcs = self.fcs.as_ref()?;
        let kw = fcs.spillover_keyword()?;
        parse_spillover(kw).ok()
    }

    /// Begin an override: copy the embedded matrix if present, else an identity
    /// matrix over the fluorescence channels.
    fn start_override(&mut self) {
        if let Some((ch, rows)) = self.active_matrix() {
            self.spill_override = Some(SpillOverride { channels: ch, rows });
            return;
        }
        if let Some(fcs) = &self.fcs {
            let idx = crate::transform::fluorescence_indices(&fcs.parameters);
            let ch: Vec<String> = idx.iter().map(|&i| fcs.parameters[i].name.clone()).collect();
            let n = ch.len();
            let rows: Vec<Vec<f64>> = (0..n)
                .map(|i| (0..n).map(|j| if i == j { 1.0 } else { 0.0 }).collect())
                .collect();
            self.spill_override = Some(SpillOverride { channels: ch, rows });
        }
    }

    fn load_override(&mut self) {
        if let Some(path) = rfd::FileDialog::new()
            .add_filter("matrix", &["csv", "json", "CSV", "JSON"]).pick_file()
        {
            match load_matrix_file(&path) {
                Ok((ch, rows)) => {
                    let missing: Vec<String> = match &self.fcs {
                        Some(fcs) => ch.iter().filter(|c| fcs.param_index(c).is_none()).cloned().collect(),
                        None => vec![],
                    };
                    if !missing.is_empty() {
                        self.status = format!("Matrix channel(s) not in this file: {}", missing.join(", "));
                        return;
                    }
                    self.spill_override = Some(SpillOverride { channels: ch, rows });
                    self.needs_reprocess = true;
                    self.status = format!("Loaded override matrix from {}", path.display());
                }
                Err(e) => self.status = format!("Load error: {}", e),
            }
        }
    }

    /// Compute a fresh spillover matrix from single-stain controls (+ unstained)
    /// and install it as the override.
    fn compute_from_controls(&mut self) {
        // Fluorescence channels from the currently loaded data file.
        let fluor: Vec<String> = match &self.fcs {
            Some(fcs) => crate::transform::fluorescence_indices(&fcs.parameters)
                .iter().map(|&i| fcs.parameters[i].name.clone()).collect(),
            None => return,
        };
        if fluor.is_empty() {
            self.status = "No fluorescence channels detected in the loaded file.".into();
            return;
        }

        let unstained = match rfd::FileDialog::new()
            .set_title("Select the UNSTAINED control")
            .add_filter("FCS", &["fcs", "FCS"]).pick_file()
        { Some(p) => p, None => return };

        let ctrls = match rfd::FileDialog::new()
            .set_title(format!("Select the {} single-stain controls (one per fluorophore)", fluor.len()))
            .add_filter("FCS", &["fcs", "FCS"]).pick_files()
        { Some(v) => v, None => return };

        if ctrls.len() != fluor.len() {
            self.status = format!(
                "Need {} single-stain controls (one per fluorescence channel: {}), got {}",
                fluor.len(), fluor.join(", "), ctrls.len()
            );
            return;
        }

        let unst = match FcsFile::open(&unstained) {
            Ok(f) => f,
            Err(e) => { self.status = format!("Unstained control: {}", e); return; }
        };
        let mut ctrl_files: Vec<(PathBuf, FcsFile)> = Vec::new();
        for p in &ctrls {
            match FcsFile::open(p) {
                Ok(f) => ctrl_files.push((p.clone(), f)),
                Err(e) => { self.status = format!("{}: {}", p.display(), e); return; }
            }
        }
        let refs: Vec<&FcsFile> = ctrl_files.iter().map(|(_, f)| f).collect();

        match compute_spillover(&fluor, &unst, &refs) {
            Ok(comp) => {
                // Filename correctness guard.
                let mut mism: Vec<String> = Vec::new();
                for (ci, &p) in comp.assigned.iter().enumerate() {
                    let fname = ctrl_files[ci].0.file_name()
                        .map(|s| s.to_string_lossy().to_string()).unwrap_or_default();
                    if fluor_token_in_filename(&fluor, &fname) != Some(p) {
                        mism.push(format!("{} → {}", fname, fluor[p]));
                    }
                }
                self.spill_override = Some(SpillOverride { channels: comp.channels, rows: comp.rows });
                self.needs_reprocess = true;
                if !self.do_compensate { self.do_compensate = true; }
                self.status = if mism.is_empty() {
                    format!("Computed spillover from controls {} (all stains matched filenames). Override active.", icon::CHECK)
                } else {
                    format!("Computed spillover — {} {} stain(s) disagree with filename: {}", icon::WARNING,
                        mism.len(), mism.join("; "))
                };
            }
            Err(e) => self.status = format!("Compute failed: {}", e),
        }
    }

    fn save_active_matrix(&mut self, channels: &[String], rows: &[Vec<f64>]) {
        if let Some(path) = rfd::FileDialog::new()
            .add_filter("CSV/JSON", &["csv", "json"]).set_file_name("spillover.csv").save_file()
        {
            match save_matrix_file(&path, channels, rows) {
                Ok(_) => self.status = format!("Saved matrix → {}", path.display()),
                Err(e) => self.status = format!("Save error: {}", e),
            }
        }
    }

    fn write_fcs_override(&mut self, channels: &[String], rows: &[Vec<f64>]) {
        if let Err(e) = SpilloverMatrix::from_parts(channels.to_vec(), rows) {
            self.status = format!("Matrix not invertible — fix before writing: {}", e);
            return;
        }
        let spill = format_spillover(channels, rows);
        let default_name = self.file_path.as_ref()
            .and_then(|p| p.file_stem())
            .map(|s| format!("{}_recomp.fcs", s.to_string_lossy()))
            .unwrap_or_else(|| "out.fcs".into());
        if let Some(path) = rfd::FileDialog::new()
            .add_filter("FCS", &["fcs"]).set_file_name(default_name).save_file()
        {
            let result = match &self.fcs {
                Some(fcs) => fcs_write::write_fcs(fcs, Some(&spill), &path),
                None => return,
            };
            self.status = match result {
                Ok(_) => format!("Wrote new FCS with this matrix → {}", path.display()),
                Err(e) => format!("Write error: {}", e),
            };
        }
    }
}

// ── Free helpers ────────────────────────────────────────────────────────────

/// Compensation that doesn't borrow the app — usable from a background thread.
/// Returns Ok(raw) when no matrix applies; Err only when a matrix exists but fails.
fn compensate_for(fcs: &FcsFile, do_compensate: bool, ov: Option<&SpillOverride>) -> Result<Vec<f64>, String> {
    if !do_compensate {
        return Ok(fcs.events.clone());
    }
    if let Some(ov) = ov {
        return SpilloverMatrix::from_parts(ov.channels.clone(), &ov.rows)
            .and_then(|s| s.apply(fcs)).map_err(|e| e.to_string());
    }
    if let Some(kw) = fcs.spillover_keyword() {
        return SpilloverMatrix::from_keyword(kw)
            .and_then(|s| s.apply(fcs)).map_err(|e| e.to_string());
    }
    Ok(fcs.events.clone()) // no matrix present → uncompensated (legitimate)
}

/// Smooth + normalize raw bin counts per the chosen mode.
fn normalize_hist(counts: Vec<f64>, norm: HistNorm) -> Vec<f64> {
    let counts = smooth_hist(&counts, 2);
    match norm {
        HistNorm::Modal => {
            let mx = counts.iter().cloned().fold(0.0_f64, f64::max).max(1e-9);
            counts.iter().map(|c| c / mx).collect()
        }
        HistNorm::Count => counts,
    }
}

/// Distinct color per overlaid sample.
fn sample_color(i: usize) -> Color32 {
    const P: [(u8, u8, u8); 8] = [
        (60, 120, 220), (220, 80, 60), (40, 170, 90), (210, 140, 0),
        (150, 60, 200), (0, 160, 170), (230, 100, 160), (130, 130, 130),
    ];
    let (r, g, b) = P[i % 8];
    Color32::from_rgb(r, g, b)
}

/// Moving-average histogram smoothing (two passes of a (2·radius+1)-bin window),
/// so sparse populations read as smooth curves rather than spikes.
fn smooth_hist(counts: &[f64], radius: usize) -> Vec<f64> {
    let pass = |c: &[f64]| -> Vec<f64> {
        let n = c.len();
        (0..n).map(|i| {
            let lo = i.saturating_sub(radius);
            let hi = (i + radius + 1).min(n);
            c[lo..hi].iter().sum::<f64>() / (hi - lo) as f64
        }).collect()
    };
    pass(&pass(counts))
}

/// Make a string safe to embed in a filename.
fn sanitize(s: &str) -> String {
    s.chars().map(|c| if c.is_alphanumeric() || c == '-' || c == '_' { c } else { '_' }).collect()
}

/// Per-user path where the recent-files list is stored.
fn recent_store_path() -> Option<PathBuf> {
    let dir = if cfg!(target_os = "macos") {
        std::env::var_os("HOME").map(|h| PathBuf::from(h).join("Library/Application Support/flowcyto"))
    } else if cfg!(windows) {
        std::env::var_os("APPDATA").map(|a| PathBuf::from(a).join("flowcyto"))
    } else {
        std::env::var_os("XDG_CONFIG_HOME").map(PathBuf::from)
            .or_else(|| std::env::var_os("HOME").map(|h| PathBuf::from(h).join(".config")))
            .map(|c| c.join("flowcyto"))
    };
    dir.map(|d| d.join("recent.json"))
}

/// Load the recent-files list, dropping entries that no longer exist on disk.
fn load_recent() -> Vec<PathBuf> {
    let p = match recent_store_path() { Some(p) => p, None => return Vec::new() };
    let txt = match std::fs::read_to_string(&p) { Ok(t) => t, Err(_) => return Vec::new() };
    serde_json::from_str::<Vec<String>>(&txt).unwrap_or_default()
        .into_iter().map(PathBuf::from).filter(|p| p.exists()).collect()
}

/// Axis-aligned bounding box of a 2-D gate in display coords (xmin,xmax,ymin,ymax).
/// `None` for 1-D interval gates.
fn shape_display_bbox(s: &GateShape) -> Option<(f64, f64, f64, f64)> {
    match s {
        GateShape::Rect { x_min, x_max, y_min, y_max } => Some((*x_min, *x_max, *y_min, *y_max)),
        GateShape::Ellipse { cx, cy, rx, ry, angle } => {
            // Half-extents of a rotated ellipse (axis-aligned box of the rotated shape).
            let (s, c) = angle.sin_cos();
            let hx = (rx * c).hypot(ry * s);
            let hy = (rx * s).hypot(ry * c);
            Some((cx - hx, cx + hx, cy - hy, cy + hy))
        }
        GateShape::Polygon { vertices } => {
            if vertices.is_empty() { return None; }
            let (mut xmn, mut xmx, mut ymn, mut ymx) = (f64::MAX, f64::MIN, f64::MAX, f64::MIN);
            for v in vertices {
                xmn = xmn.min(v[0]); xmx = xmx.max(v[0]);
                ymn = ymn.min(v[1]); ymx = ymx.max(v[1]);
            }
            Some((xmn, xmx, ymn, ymx))
        }
        GateShape::Range { .. } => None,
        GateShape::Boolean { .. } => None,
    }
}

fn shape_label(s: &GateShape) -> &'static str {
    match s {
        GateShape::Rect { .. } => "rect",
        GateShape::Ellipse { .. } => "ellipse",
        GateShape::Polygon { .. } => "polygon",
        GateShape::Range { .. } => "interval",
        GateShape::Boolean { .. } => "boolean",
    }
}

/// Default per-channel transforms: scatter/Time linear, fluorescence logicle.
fn default_transforms(fcs: &FcsFile) -> Vec<AxisTransform> {
    fcs.parameters.iter().map(|p| {
        let n = p.name.to_uppercase();
        if n.starts_with("FSC") || n.starts_with("SSC") || n.eq_ignore_ascii_case("TIME") {
            AxisTransform::Linear
        } else {
            AxisTransform::default_logicle()
        }
    }).collect()
}

/// Carry transforms to a new sample by channel NAME (panels may differ in order),
/// falling back to defaults for channels the previous sample didn't have.
fn rekey_transforms(
    old_fcs: Option<&FcsFile>, old_tf: &[AxisTransform], new_fcs: &FcsFile,
) -> Vec<AxisTransform> {
    let mut by_name: HashMap<String, AxisTransform> = HashMap::new();
    if let Some(of) = old_fcs {
        for (i, p) in of.parameters.iter().enumerate() {
            if let Some(t) = old_tf.get(i) {
                by_name.insert(p.name.to_uppercase(), t.clone());
            }
        }
    }
    let defaults = default_transforms(new_fcs);
    new_fcs.parameters.iter().enumerate().map(|(i, p)| {
        by_name.get(&p.name.to_uppercase()).cloned().unwrap_or_else(|| defaults[i].clone())
    }).collect()
}

/// Channels referenced by any gate that are absent from `fcs` (panel mismatch).
fn missing_gate_channels(gates: &[crate::gating::Gate], fcs: &FcsFile) -> Vec<String> {
    let mut miss: Vec<String> = Vec::new();
    for g in gates {
        for ch in [&g.x_channel, &g.y_channel] {
            if fcs.param_index(ch).is_none() && !miss.iter().any(|m| m.eq_ignore_ascii_case(ch)) {
                miss.push(ch.clone());
            }
        }
    }
    miss
}

fn param_full_label(f: &FcsFile, i: usize) -> String {
    let p = &f.parameters[i];
    match &p.label { Some(l) if !l.is_empty() => format!("{} ({})", p.name, l), _ => p.name.clone() }
}

/// Strip the " (stain)" suffix to recover the bare $PnN channel name.
fn x_name_base(full: &str) -> String {
    match full.find(" (") { Some(i) => full[..i].to_string(), None => full.to_string() }
}

/// Compact channel name (drop the "-A" area suffix) for gate labels.
fn short_chan(name: &str) -> String {
    name.strip_suffix("-A").unwrap_or(name).to_string()
}

/// Human description of a display transform (for the reproducibility report).
fn describe_tf(t: &AxisTransform) -> String {
    match t {
        AxisTransform::Linear => "Linear".into(),
        AxisTransform::Log { .. } => "Log".into(),
        AxisTransform::Asinh { cofactor } => format!("Asinh (cofactor {})", cofactor),
        AxisTransform::Logicle { t, w, m, a } => format!("Logicle (t={}, w={}, m={}, a={})", t, w, m, a),
    }
}

/// Format a gate bound, collapsing the ±1e12 quadrant sentinels to ±∞.
fn rb(v: f64) -> String {
    if v >= 1e11 { "+∞".into() } else if v <= -1e11 { "−∞".into() } else { format!("{:.1}", v) }
}

/// One-line gate shape description in display coordinates (for the report).
fn shape_desc(s: &GateShape) -> String {
    match s {
        GateShape::Rect { x_min, x_max, y_min, y_max } =>
            format!("rect: X∈[{}, {}], Y∈[{}, {}]", rb(*x_min), rb(*x_max), rb(*y_min), rb(*y_max)),
        GateShape::Ellipse { cx, cy, rx, ry, angle } =>
            format!("ellipse: c=({:.1}, {:.1}), r=({:.1}, {:.1}), θ={:.0}°", cx, cy, rx, ry, angle.to_degrees()),
        GateShape::Polygon { vertices } => format!("polygon: {} vertices", vertices.len()),
        GateShape::Range { x_min, x_max } => format!("interval: [{}, {}]", rb(*x_min), rb(*x_max)),
        GateShape::Boolean { op, refs } => format!("{:?} of {} gate(s)", op, refs.len()),
    }
}

/// Iso-density contour line segments via marching squares over the density grid.
/// `hist[i][j]` is the count in bin (i,j); returns line segments in display coords.
fn contour_segments(hist: &[Vec<u32>], n: usize, xmin: f64, xmax: f64, ymin: f64, ymax: f64, levels: &[f64]) -> Vec<[[f64; 2]; 2]> {
    let mut segs = Vec::new();
    if n < 2 || hist.len() < n { return segs; }
    let cx = |i: usize| xmin + (i as f64 + 0.5) * (xmax - xmin) / n as f64;
    let cy = |j: usize| ymin + (j as f64 + 0.5) * (ymax - ymin) / n as f64;
    let lerp = |t0: f64, t1: f64, p0: [f64; 2], p1: [f64; 2], lev: f64| -> [f64; 2] {
        let d = t1 - t0;
        let f = if d.abs() < 1e-12 { 0.5 } else { ((lev - t0) / d).clamp(0.0, 1.0) };
        [p0[0] + f * (p1[0] - p0[0]), p0[1] + f * (p1[1] - p0[1])]
    };
    for &lev in levels {
        for i in 0..n - 1 {
            for j in 0..n - 1 {
                let (va, vb, vc, vd) = (hist[i][j] as f64, hist[i + 1][j] as f64,
                                        hist[i + 1][j + 1] as f64, hist[i][j + 1] as f64);
                let mut case = 0u8;
                if va > lev { case |= 1; }
                if vb > lev { case |= 2; }
                if vc > lev { case |= 4; }
                if vd > lev { case |= 8; }
                if case == 0 || case == 15 { continue; }
                let (pa, pb, pc, pd) = ([cx(i), cy(j)], [cx(i + 1), cy(j)],
                                        [cx(i + 1), cy(j + 1)], [cx(i), cy(j + 1)]);
                let eb = lerp(va, vb, pa, pb, lev); // bottom edge
                let er = lerp(vb, vc, pb, pc, lev); // right
                let et = lerp(vc, vd, pc, pd, lev); // top
                let el = lerp(vd, va, pd, pa, lev); // left
                match case {
                    1 | 14 => segs.push([el, eb]),
                    2 | 13 => segs.push([eb, er]),
                    3 | 12 => segs.push([el, er]),
                    4 | 11 => segs.push([er, et]),
                    6 | 9  => segs.push([eb, et]),
                    7 | 8  => segs.push([el, et]),
                    5      => { segs.push([el, eb]); segs.push([er, et]); }
                    10     => { segs.push([eb, er]); segs.push([et, el]); }
                    _ => {}
                }
            }
        }
    }
    segs
}

/// Draggable handle points of a gate, in the gate's display coords.
fn gate_handles(shape: &GateShape) -> Vec<[f64; 2]> {
    match shape {
        GateShape::Rect { x_min, x_max, y_min, y_max } =>
            vec![[*x_min, *y_min], [*x_max, *y_min], [*x_max, *y_max], [*x_min, *y_max]],
        GateShape::Ellipse { cx, cy, rx, ry, angle } => {
            // Handles live on the ellipse's own (rotated) axes:
            //   0 = +major, 1 = +minor, 2 = −major, 3 = −minor, 4 = rotation.
            let (s, c) = angle.sin_cos();
            let major = |k: f64| [cx + k * rx * c, cy + k * rx * s];
            let minor = |k: f64| [cx - k * ry * s, cy + k * ry * c];
            vec![major(1.0), minor(1.0), major(-1.0), minor(-1.0), major(1.35)]
        }
        GateShape::Polygon { vertices } => vertices.clone(),
        GateShape::Range { x_min, x_max } => vec![[*x_min, 0.0], [*x_max, 0.0]],
        GateShape::Boolean { .. } => Vec::new(),
    }
}

/// Move handle `h` of `shape` to gate-display coords (gx, gy).
fn apply_gate_handle(shape: &mut GateShape, h: usize, gx: f64, gy: f64) {
    match shape {
        GateShape::Rect { x_min, x_max, y_min, y_max } => {
            match h { 0 => { *x_min = gx; *y_min = gy; } 1 => { *x_max = gx; *y_min = gy; }
                      2 => { *x_max = gx; *y_max = gy; } 3 => { *x_min = gx; *y_max = gy; } _ => {} }
            if x_min > x_max { std::mem::swap(x_min, x_max); }
            if y_min > y_max { std::mem::swap(y_min, y_max); }
        }
        GateShape::Ellipse { cx, cy, rx, ry, angle } => {
            // Center stays fixed; project the dragged point onto the ellipse's axes.
            let (dx, dy) = (gx - *cx, gy - *cy);
            let (s, c) = angle.sin_cos();
            match h {
                0 | 2 => { *rx = (dx * c + dy * s).abs().max(1e-9); }       // major axis
                1 | 3 => { *ry = (-dx * s + dy * c).abs().max(1e-9); }      // minor axis
                4 => { *angle = dy.atan2(dx); }                            // rotation handle
                _ => {}
            }
        }
        GateShape::Polygon { vertices } => { if h < vertices.len() { vertices[h] = [gx, gy]; } }
        GateShape::Range { x_min, x_max } => { if h == 0 { *x_min = gx; } else { *x_max = gx; } }
        GateShape::Boolean { .. } => {}
    }
}

/// Translate an entire gate by (dx, dy) in its display coordinates.
fn translate_shape(shape: &mut GateShape, dx: f64, dy: f64) {
    match shape {
        GateShape::Rect { x_min, x_max, y_min, y_max } => {
            *x_min += dx; *x_max += dx; *y_min += dy; *y_max += dy;
        }
        GateShape::Ellipse { cx, cy, .. } => { *cx += dx; *cy += dy; }
        GateShape::Polygon { vertices } => { for v in vertices { v[0] += dx; v[1] += dy; } }
        GateShape::Range { x_min, x_max } => { *x_min += dx; *x_max += dx; }
        GateShape::Boolean { .. } => {}
    }
}

fn rubber_band(mode: DrawMode, s: [f64; 2], c: [f64; 2]) -> Vec<[f64; 2]> {
    match mode {
        DrawMode::Ellipse => {
            let cx = (s[0] + c[0]) / 2.0; let cy = (s[1] + c[1]) / 2.0;
            let rx = (c[0] - s[0]).abs() / 2.0; let ry = (c[1] - s[1]).abs() / 2.0;
            (0..=48).map(|i| {
                let th = std::f64::consts::TAU * i as f64 / 48.0;
                [cx + rx * th.cos(), cy + ry * th.sin()]
            }).collect()
        }
        _ => vec![[s[0], s[1]], [c[0], s[1]], [c[0], c[1]], [s[0], c[1]], [s[0], s[1]]],
    }
}

fn shape_from_drag(mode: DrawMode, s: [f64; 2], c: [f64; 2]) -> GateShape {
    match mode {
        DrawMode::Ellipse => GateShape::Ellipse {
            cx: (s[0] + c[0]) / 2.0, cy: (s[1] + c[1]) / 2.0,
            rx: (c[0] - s[0]).abs() / 2.0, ry: (c[1] - s[1]).abs() / 2.0, angle: 0.0,
        },
        _ => GateShape::Rect {
            x_min: s[0].min(c[0]), x_max: s[0].max(c[0]),
            y_min: s[1].min(c[1]), y_max: s[1].max(c[1]),
        },
    }
}

fn fmt(v: f64) -> String {
    if !v.is_finite() { return "—".into(); }
    if v.abs() >= 1000.0 { format!("{:.0}", v) } else { format!("{:.2}", v) }
}

/// Tick label in data units with K/M suffixes (0, 30K, 262K, -10K, 0.5).
fn fmt_data_tick(v: f64) -> String {
    let a = v.abs();
    if a < 0.5 { return "0".into(); }
    let sign = if v < 0.0 { "-" } else { "" };
    if a >= 1e6 {
        format!("{}{:.1}M", sign, a / 1e6)
    } else if a >= 1e3 {
        format!("{}{:.0}K", sign, a / 1e3)
    } else if a >= 1.0 {
        format!("{}{:.0}", sign, a)
    } else {
        format!("{:.1}", v)
    }
}

fn data_range(v: &[f64]) -> (f64, f64) {
    let mut lo = f64::INFINITY; let mut hi = f64::NEG_INFINITY;
    for &x in v {
        if x.is_finite() { lo = lo.min(x); hi = hi.max(x); }
    }
    if !lo.is_finite() || !hi.is_finite() { return (0.0, 1.0); }
    if lo >= hi { hi = lo + 1.0; }
    (lo, hi)
}

fn scatter_display_range(buckets: &[Vec<[f64; 2]>]) -> ((f64, f64), (f64, f64)) {
    let mut xlo = f64::INFINITY; let mut xhi = f64::NEG_INFINITY;
    let mut ylo = f64::INFINITY; let mut yhi = f64::NEG_INFINITY;
    for b in buckets { for p in b {
        xlo = xlo.min(p[0]); xhi = xhi.max(p[0]);
        ylo = ylo.min(p[1]); yhi = yhi.max(p[1]);
    }}
    if !xlo.is_finite() { xlo = 0.0; xhi = 1.0; }
    if !ylo.is_finite() { ylo = 0.0; yhi = 1.0; }
    ((xlo, xhi), (ylo, yhi))
}

fn bin_of(x: f64, lo: f64, hi: f64, n: usize) -> usize {
    if hi <= lo { return 0; }
    ((x - lo) / (hi - lo) * n as f64).clamp(0.0, (n - 1) as f64) as usize
}

/// Density-sample display-space points into N_BUCKETS color buckets (shared by
/// the single scatter and the grid cells). `cap` bounds the points drawn so dense
/// grids (up to 36 cells) stay responsive.
fn bucketize(dx: &[f64], dy: &[f64], cap: usize) -> Vec<Vec<[f64; 2]>> {
    let mut buckets: Vec<Vec<[f64; 2]>> = vec![Vec::new(); N_BUCKETS];
    let nk = dx.len();
    if nk == 0 { return buckets; }
    let (xmin, xmax) = data_range(dx);
    let (ymin, ymax) = data_range(dy);
    let hist = density_hist(dx, dy, DENSITY_BINS, xmin, xmax, ymin, ymax);
    let max_d = hist.iter().flat_map(|r| r.iter()).copied().max().unwrap_or(1).max(1);
    let n_sample = cap.max(1).min(nk);
    let step = (nk / n_sample).max(1);
    for k in (0..nk).step_by(step) {
        let (x, y) = (dx[k], dy[k]);
        let bx = bin_of(x, xmin, xmax, DENSITY_BINS);
        let by = bin_of(y, ymin, ymax, DENSITY_BINS);
        let t = (hist[bx][by] as f64 / max_d as f64).sqrt();
        let b = ((t * (N_BUCKETS - 1) as f64) as usize).min(N_BUCKETS - 1);
        buckets[b].push([x, y]);
    }
    buckets
}

fn density_hist(xs: &[f64], ys: &[f64], n: usize, xlo: f64, xhi: f64, ylo: f64, yhi: f64) -> Vec<Vec<u32>> {
    let mut h = vec![vec![0u32; n]; n];
    for (&x, &y) in xs.iter().zip(ys.iter()) {
        let bx = bin_of(x, xlo, xhi, n);
        let by = bin_of(y, ylo, yhi, n);
        h[bx][by] = h[bx][by].saturating_add(1);
    }
    h
}

/// Grid marks: decades for nonlinear axes (data 0, ±10^k mapped to display);
/// for linear axes, fall back to a simple uniform grid.
fn nonlinear_grid(ct: &CompiledTransform, linear: bool, inp: egui_plot::GridInput) -> Vec<GridMark> {
    let (lo, hi) = inp.bounds;
    if linear {
        // default-ish: ~10 ticks
        let span = (hi - lo).abs().max(1e-9);
        let raw_step = span / 8.0;
        let mag = 10f64.powf(raw_step.log10().floor());
        let step = (raw_step / mag).round().max(1.0) * mag;
        let mut marks = Vec::new();
        let start = (lo / step).floor() as i64;
        let end = (hi / step).ceil() as i64;
        for k in start..=end {
            marks.push(GridMark { value: k as f64 * step, step_size: step });
        }
        return marks;
    }

    // nonlinear: place at 0 and ±10^k within data range
    let d_lo = ct.inverse(lo);
    let d_hi = ct.inverse(hi);
    let mut marks = Vec::new();
    // zero
    let z = ct.forward(0.0);
    if z >= lo && z <= hi { marks.push(GridMark { value: z, step_size: (hi - lo) / 5.0 }); }
    // Decade marks only (clean 10²/10³/10⁴ labels, the flow-cytometry convention).
    for sign in [1.0f64, -1.0] {
        for k in 0..7 {
            let data = sign * 10f64.powi(k);
            if data >= d_lo.min(d_hi) && data <= d_lo.max(d_hi) {
                let disp = ct.forward(data);
                if disp >= lo && disp <= hi {
                    marks.push(GridMark { value: disp, step_size: (hi - lo) / 5.0 });
                }
            }
        }
    }
    marks
}

// ── Entry point ────────────────────────────────────────────────────────────

// ── Native macOS menu bar (muda) ─────────────────────────────────────────────

/// Menu item ids we dispatch on. Holds the live `Menu` so the native menu isn't
/// torn down. Actions mirror the in-app toolbar/keyboard so nothing diverges.
#[cfg(target_os = "macos")]
struct MenuState {
    _menu: muda::Menu,
    open: muda::MenuId,
    save_gates: muda::MenuId,
    save_session: muda::MenuId,
    load_session: muda::MenuId,
    save_plot: muda::MenuId,
    check_updates: muda::MenuId,
    undo: muda::MenuId,
    redo: muda::MenuId,
    theme: muda::MenuId,
    zoom_in: muda::MenuId,
    zoom_out: muda::MenuId,
    zoom_reset: muda::MenuId,
    tabs: [muda::MenuId; 6], // Plot, Histogram, Stats, Batch, Spillover, QC
}

#[cfg(target_os = "macos")]
fn build_menu() -> MenuState {
    use muda::{accelerator::Accelerator, AboutMetadata, IsMenuItem, Menu, MenuItem,
               PredefinedMenuItem, Submenu};
    let acc = |s: &str| s.parse::<Accelerator>().ok();
    let menu = Menu::new();

    // Application menu (About / Hide / Quit).
    let app_m = Submenu::new("flowcyto", true);
    let about = PredefinedMenuItem::about(Some("About flowcyto"), Some(AboutMetadata {
        name: Some("flowcyto".into()),
        version: Some(env!("CARGO_PKG_VERSION").to_string()),
        ..Default::default()
    }));
    let check_updates = MenuItem::new("Check for Updates…", true, None);
    let _ = app_m.append_items(&[
        &about, &check_updates, &PredefinedMenuItem::separator(),
        &PredefinedMenuItem::hide(None), &PredefinedMenuItem::quit(None),
    ]);

    // File.
    let file_m = Submenu::new("File", true);
    let open = MenuItem::new("Open FCS…", true, acc("CmdOrCtrl+O"));
    let save_gates = MenuItem::new("Save Gates…", true, acc("CmdOrCtrl+S"));
    let save_session = MenuItem::new("Save Session…", true, acc("CmdOrCtrl+Shift+S"));
    let load_session = MenuItem::new("Load Session…", true, None);
    let save_plot = MenuItem::new("Save Plot as PNG…", true, None);
    let _ = file_m.append_items(&[
        &open, &PredefinedMenuItem::separator(),
        &save_gates, &save_session, &load_session,
        &PredefinedMenuItem::separator(), &save_plot,
    ]);

    // Edit.
    let edit_m = Submenu::new("Edit", true);
    let undo = MenuItem::new("Undo", true, acc("CmdOrCtrl+Z"));
    let redo = MenuItem::new("Redo", true, acc("CmdOrCtrl+Shift+Z"));
    let _ = edit_m.append_items(&[&undo, &redo]);

    // View (tab switches + theme).
    let view_m = Submenu::new("View", true);
    let names = ["Plot", "Histogram", "Stats", "Batch", "Spillover", "QC"];
    let tab_items: Vec<MenuItem> = names.iter().enumerate()
        .map(|(i, n)| MenuItem::new(*n, true, acc(&format!("CmdOrCtrl+{}", i + 1)))).collect();
    let theme = MenuItem::new("Toggle Light/Dark", true, None);
    let zoom_in = MenuItem::new("Zoom In", true, acc("CmdOrCtrl+Plus"));
    let zoom_out = MenuItem::new("Zoom Out", true, acc("CmdOrCtrl+-"));
    let zoom_reset = MenuItem::new("Actual Size", true, acc("CmdOrCtrl+0"));
    let tab_refs: Vec<&dyn IsMenuItem> = tab_items.iter().map(|m| m as &dyn IsMenuItem).collect();
    let _ = view_m.append_items(&tab_refs);
    let _ = view_m.append_items(&[
        &PredefinedMenuItem::separator() as &dyn IsMenuItem, &theme,
        &PredefinedMenuItem::separator(), &zoom_in, &zoom_out, &zoom_reset,
    ]);

    let _ = menu.append_items(&[&app_m, &file_m, &edit_m, &view_m]);
    menu.init_for_nsapp();

    MenuState {
        open: open.id().clone(),
        save_gates: save_gates.id().clone(),
        save_session: save_session.id().clone(),
        load_session: load_session.id().clone(),
        save_plot: save_plot.id().clone(),
        check_updates: check_updates.id().clone(),
        undo: undo.id().clone(),
        redo: redo.id().clone(),
        theme: theme.id().clone(),
        zoom_in: zoom_in.id().clone(),
        zoom_out: zoom_out.id().clone(),
        zoom_reset: zoom_reset.id().clone(),
        tabs: std::array::from_fn(|i| tab_items[i].id().clone()),
        _menu: menu,
    }
}

/// Bundle the Inter UI font (Regular + SemiBold for headings) and the Phosphor icon
/// font, so the whole app renders with crisp, consistent typography and line-icons
/// instead of the built-in font + platform emoji.
fn setup_fonts(ctx: &egui::Context) {
    use egui::{FontData, FontDefinitions, FontFamily};
    let mut fonts = FontDefinitions::default();
    fonts.font_data.insert(
        "Inter".to_owned(),
        FontData::from_static(include_bytes!("../assets/fonts/Inter-Regular.ttf")),
    );
    fonts.font_data.insert(
        "Inter-SemiBold".to_owned(),
        FontData::from_static(include_bytes!("../assets/fonts/Inter-SemiBold.ttf")),
    );
    // Inter as the primary proportional face; keep egui's bundled font as a fallback
    // for any glyph Inter lacks.
    fonts.families.entry(FontFamily::Proportional).or_default().insert(0, "Inter".to_owned());
    // A dedicated SemiBold family for headings.
    fonts.families.insert(
        FontFamily::Name("heading".into()),
        vec!["Inter-SemiBold".to_owned(), "Inter".to_owned()],
    );
    // Phosphor icon glyphs. `add_to_fonts` registers the "phosphor" font and inserts
    // it at index 1 of the Proportional family (right after Inter, ahead of egui's
    // default fallback fonts) so PUA icon codepoints resolve to Phosphor.
    //
    // NOTE: the bundled Inter TTFs are subset to DROP their Private-Use-Area glyphs.
    // Upstream Inter maps ~745 PUA codepoints; since Inter sits at index 0 it would
    // otherwise shadow the Phosphor icons (which live in the same PUA range) with
    // blank/placeholder glyphs — the cause of the "weird icon" rendering. With those
    // stripped, the index-1 Phosphor face is reached for every icon. Regenerate the
    // subset with `assets/fonts/strip-pua.py` if the Inter TTFs are ever updated.
    egui_phosphor::add_to_fonts(&mut fonts, egui_phosphor::Variant::Regular);
    // `add_to_fonts` only touches Proportional; add Phosphor to the heading and
    // monospace families too, since headings (e.g. the QC tab) and some labels embed
    // inline icons.
    fonts.families.entry(FontFamily::Name("heading".into())).or_default().insert(1, "phosphor".to_owned());
    fonts.families.entry(FontFamily::Monospace).or_default().insert(1, "phosphor".to_owned());
    ctx.set_fonts(fonts);
}

/// One-time style setup: spacing, padding, rounded widgets, readable text sizes.
/// (Spacing/text persist across the per-frame `set_visuals` calls.)
fn setup_style(ctx: &egui::Context) {
    use egui::{FontFamily::{self, Monospace, Proportional}, FontId, TextStyle};
    let mut style = (*ctx.style()).clone();
    style.spacing.item_spacing = egui::vec2(8.0, 6.0);
    style.spacing.button_padding = egui::vec2(8.0, 4.0);
    style.spacing.interact_size.y = 24.0;
    style.spacing.menu_margin = egui::Margin::same(6.0);
    style.text_styles = [
        (TextStyle::Heading, FontId::new(16.5, FontFamily::Name("heading".into()))),
        (TextStyle::Body, FontId::new(14.0, Proportional)),
        (TextStyle::Button, FontId::new(14.0, Proportional)),
        (TextStyle::Monospace, FontId::new(13.0, Monospace)),
        (TextStyle::Small, FontId::new(11.0, Proportional)),
    ].into();
    ctx.set_style(style);
}

/// Dark/light visuals with a teal accent (matching the app icon), cohesive panel
/// tones, soft shadows, and accent-tinted widget states. Rebuilt each frame because
/// `set_visuals` resets visuals.
fn themed_visuals(dark: bool) -> egui::Visuals {
    let mut v = if dark { egui::Visuals::dark() } else { egui::Visuals::light() };
    let accent = Color32::from_rgb(38, 162, 156);
    v.selection.bg_fill = accent.linear_multiply(if dark { 0.55 } else { 0.35 });
    v.selection.stroke = Stroke::new(1.0, accent);
    v.hyperlink_color = accent;

    if dark {
        let panel = Color32::from_rgb(30, 34, 41);
        let widget = Color32::from_rgb(42, 47, 56);
        let widget_hi = Color32::from_rgb(52, 58, 69);
        v.window_fill = Color32::from_rgb(24, 27, 33);
        v.panel_fill = panel;
        v.faint_bg_color = Color32::from_rgb(37, 42, 50);     // table stripe / faint fills
        v.extreme_bg_color = Color32::from_rgb(18, 20, 25);   // text edits, plot background
        v.widgets.noninteractive.bg_fill = panel;
        v.widgets.noninteractive.weak_bg_fill = panel;
        v.widgets.noninteractive.bg_stroke = Stroke::new(1.0, Color32::from_rgb(46, 51, 60));
        v.widgets.noninteractive.fg_stroke = Stroke::new(1.0, Color32::from_rgb(192, 198, 207));
        v.widgets.inactive.bg_fill = widget;
        v.widgets.inactive.weak_bg_fill = widget;
        v.widgets.inactive.bg_stroke = Stroke::new(1.0, Color32::from_rgb(58, 64, 75));
        v.widgets.inactive.fg_stroke = Stroke::new(1.0, Color32::from_rgb(208, 213, 221));
        v.widgets.hovered.bg_fill = widget_hi;
        v.widgets.hovered.weak_bg_fill = widget_hi;
        v.widgets.hovered.bg_stroke = Stroke::new(1.0, accent.linear_multiply(0.75));
        v.widgets.active.bg_fill = accent.linear_multiply(0.9);
        v.widgets.active.weak_bg_fill = accent.linear_multiply(0.9);
        v.widgets.active.bg_stroke = Stroke::new(1.0, accent);
    } else {
        let panel = Color32::from_rgb(252, 252, 253);
        let widget = Color32::from_rgb(255, 255, 255);
        let widget_hi = Color32::from_rgb(236, 245, 244);
        v.window_fill = Color32::from_rgb(246, 247, 250);
        v.panel_fill = panel;
        v.faint_bg_color = Color32::from_rgb(241, 243, 246);
        v.extreme_bg_color = Color32::from_rgb(255, 255, 255);
        v.widgets.noninteractive.bg_fill = panel;
        v.widgets.noninteractive.weak_bg_fill = panel;
        v.widgets.noninteractive.bg_stroke = Stroke::new(1.0, Color32::from_rgb(223, 227, 233));
        v.widgets.noninteractive.fg_stroke = Stroke::new(1.0, Color32::from_rgb(58, 64, 72));
        v.widgets.inactive.bg_fill = widget;
        v.widgets.inactive.weak_bg_fill = Color32::from_rgb(244, 246, 248);
        v.widgets.inactive.bg_stroke = Stroke::new(1.0, Color32::from_rgb(214, 219, 226));
        v.widgets.inactive.fg_stroke = Stroke::new(1.0, Color32::from_rgb(44, 49, 57));
        v.widgets.hovered.bg_fill = widget_hi;
        v.widgets.hovered.weak_bg_fill = widget_hi;
        v.widgets.hovered.bg_stroke = Stroke::new(1.0, accent.linear_multiply(0.8));
        v.widgets.active.bg_fill = accent.linear_multiply(0.22);
        v.widgets.active.weak_bg_fill = accent.linear_multiply(0.22);
        v.widgets.active.bg_stroke = Stroke::new(1.0, accent);
    }

    let r = egui::Rounding::same(6.0);
    for w in [
        &mut v.widgets.noninteractive, &mut v.widgets.inactive,
        &mut v.widgets.hovered, &mut v.widgets.active, &mut v.widgets.open,
    ] {
        w.rounding = r;
    }
    v.window_rounding = egui::Rounding::same(10.0);
    v.menu_rounding = egui::Rounding::same(8.0);
    v.window_shadow = egui::epaint::Shadow {
        offset: egui::vec2(0.0, 4.0),
        blur: 18.0,
        spread: 0.0,
        color: Color32::from_black_alpha(if dark { 96 } else { 38 }),
    };
    v.popup_shadow = egui::epaint::Shadow {
        offset: egui::vec2(0.0, 3.0),
        blur: 12.0,
        spread: 0.0,
        color: Color32::from_black_alpha(if dark { 90 } else { 30 }),
    };
    v
}

/// The application icon, embedded once and shared by every platform's icon path
/// (the eframe window icon below + the macOS Dock setter). Same source PNG the
/// `.app` bundle's `.icns` and the Windows `.exe` resource are generated from.
const APP_ICON_PNG: &[u8] = include_bytes!("../packaging/icon.png");

/// Decode the embedded icon into eframe's `IconData` for the window/taskbar icon.
/// Used by winit's `with_icon` (effective on **Windows** and Linux/X11; a no-op on
/// macOS, where `set_dock_icon` handles it). Returns `None` if decoding fails so a
/// bad icon never blocks the GUI.
fn app_icon() -> Option<egui::IconData> {
    let img = image::load_from_memory(APP_ICON_PNG).ok()?.to_rgba8();
    let (width, height) = img.dimensions();
    Some(egui::IconData { rgba: img.into_raw(), width, height })
}

/// Set the macOS Dock / application icon at runtime from the embedded PNG.
///
/// winit does not support setting the dock icon on macOS (`with_icon` is a no-op
/// there), and the bundle's `.icns` only covers the installed `.app` — so the bare
/// binary (e.g. `cargo run` / `flowcyto gui file.fcs` from a terminal) otherwise
/// shows a generic icon. `NSApplication setApplicationIconImage:` fixes both: it
/// applies to the bare binary and harmlessly re-asserts the bundle's icon.
///
/// Called from `update()` (not the `run_native` closure — a closure-time set is
/// clobbered when eframe finishes launching and macOS installs the default icon).
/// Every step is guarded and bails silently on failure — a missing icon must never
/// panic the GUI.
#[cfg(target_os = "macos")]
fn set_dock_icon() {
    use objc2::ClassType; // brings `alloc` into scope
    use objc2_app_kit::{NSApplication, NSImage};
    use objc2_foundation::{MainThreadMarker, NSData};

    let Some(mtm) = MainThreadMarker::new() else { return };
    let data = NSData::with_bytes(APP_ICON_PNG);
    let Some(image) = NSImage::initWithData(NSImage::alloc(), &data) else { return };
    let app = NSApplication::sharedApplication(mtm);
    // SAFETY: called on the main thread with a valid NSImage we just constructed;
    // setting the application icon image has no other invariants.
    unsafe { app.setApplicationIconImage(Some(&image)); }
}

pub fn run_gui(initial_file: Option<&Path>) {
    // Window/taskbar icon: effective on Windows + Linux/X11. On Windows this sets it
    // via WM_SETICON, independent of the .exe's embedded resource icon (which winit
    // does not install on its window class — it registers with hIcon = 0). No-op on
    // macOS; `set_dock_icon` covers that from `update()`.
    let mut viewport = egui::ViewportBuilder::default()
        .with_inner_size([1500.0_f32, 950.0_f32]).with_title("flowcyto");
    if let Some(icon) = app_icon() {
        viewport = viewport.with_icon(icon);
    }
    let options = eframe::NativeOptions { viewport, ..Default::default() };
    eframe::run_native("flowcyto", options, Box::new(|cc| {
        setup_fonts(&cc.egui_ctx);
        setup_style(&cc.egui_ctx);
        let mut app = FlowCytoApp { recent_files: load_recent(), ..Default::default() };
        if let Some(p) = initial_file { app.load_file(p); }
        #[cfg(target_os = "macos")]
        { app.menu = Some(build_menu()); }
        Ok(Box::new(app))
    })).expect("failed to launch GUI");
}
