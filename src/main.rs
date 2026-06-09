// Numeric/matrix code indexes intentionally (matrix multiply, event offsets,
// 2-D grids); a few CLI handlers and display tuples are inherently wide.
#![allow(clippy::needless_range_loop, clippy::too_many_arguments, clippy::type_complexity)]
// On Windows release builds, use the GUI subsystem so no console window flashes
// when the app is launched from a shortcut. (CLI use re-attaches the parent
// console at startup — see `maybe_attach_console`.)
#![cfg_attr(all(target_os = "windows", not(debug_assertions)), windows_subsystem = "windows")]

mod autogate;
mod compensation;
mod fcs;
mod fcs_write;
mod gating;
mod gatingml;
mod gui;
mod logicle;
mod popstats;
mod qc;
mod r_bridge;
mod report;
mod selftest;
mod stats;
#[cfg(test)]
mod test_util;
mod transform;
mod update;

use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use clap::{Parser, Subcommand, ValueEnum};

use compensation::{
    compute_spillover, fluor_token_in_filename, format_spillover, load_matrix_file,
    max_off_diagonal, parse_spillover, save_matrix_file, SpilloverMatrix,
};
use fcs::FcsFile;
use gating::{apply_gates, Gate};
use stats::{print_stats_table, Stats};
use transform::{apply_asinh, fluorescence_indices};

// ── CLI definition ────────────────────────────────────────────────────────

#[derive(Parser)]
#[command(
    name = "flowcyto",
    about = "Flow cytometry FCS analyzer — compensation · asinh transform · gating · stats · GUI",
    version
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Launch the graphical interface (default if no subcommand given).
    Gui {
        /// Optional FCS file to open immediately.
        file: Option<PathBuf>,
    },

    /// Print FCS metadata (version, keywords, channel list).
    Info {
        file: PathBuf,
    },

    /// Check GitHub for a newer release (the only command that touches the network).
    Update,

    /// Verify the numeric layers against frozen flowCore golden values (offline).
    Selftest,

    /// Print per-channel summary statistics.
    Stats {
        file: PathBuf,
        #[arg(long)]
        compensate: bool,
        #[arg(long, value_enum, default_value = "none")]
        transform: TransformArg,
        #[arg(long, default_value_t = 150.0, value_parser = parse_positive)]
        cofactor: f64,
    },

    /// Export event data to CSV.
    Export {
        file: PathBuf,
        #[arg(short, long, help = "Output CSV path (default: <file>.csv)")]
        output: Option<PathBuf>,
        #[arg(long)]
        compensate: bool,
        #[arg(long, value_enum, default_value = "none")]
        transform: TransformArg,
        #[arg(long, default_value_t = 150.0, value_parser = parse_positive)]
        cofactor: f64,
    },

    /// Apply gate definitions from a JSON file and report population counts.
    Gate {
        file: PathBuf,
        #[arg(short, long)]
        gates: PathBuf,
        #[arg(long)]
        compensate: bool,
        #[arg(long, value_enum, default_value = "none")]
        transform: TransformArg,
        #[arg(long, default_value_t = 150.0, value_parser = parse_positive)]
        cofactor: f64,
    },

    /// Per-population statistics (count, %parent, %total, per-channel MFI) → tidy CSV.
    Popstats {
        file: PathBuf,
        #[arg(short, long)]
        gates: PathBuf,
        #[arg(long)]
        compensate: bool,
        #[arg(short, long)]
        output: Option<PathBuf>,
    },

    /// Show the $SPILLOVER (compensation) matrix embedded in the FCS file.
    Spillover {
        file: PathBuf,
    },

    /// Export gate definitions (JSON) to Gating-ML 2.0 XML for interoperability with
    /// flowCore/CytoML, FlowKit, etc. The sample FCS supplies the $SPILLOVER channel
    /// names so compensated fluorescence dimensions reference the file's matrix.
    GatingMl {
        file: PathBuf,
        #[arg(short, long)]
        gates: PathBuf,
        #[arg(long)]
        compensate: bool,
        #[arg(short, long, help = "Output .xml path")]
        output: PathBuf,
    },

    /// Compute a spillover matrix from single-stain controls + an unstained control.
    ComputeSpillover {
        /// Unstained control FCS file.
        #[arg(long)]
        unstained: PathBuf,
        /// Single-stain control FCS files (repeat --control once per fluorophore).
        #[arg(long = "control", required = true)]
        controls: Vec<PathBuf>,
        /// Write the computed matrix to this CSV/JSON file.
        #[arg(short, long)]
        output: Option<PathBuf>,
    },

    /// Write a NEW .fcs file with an overridden $SPILLOVER matrix (original untouched).
    RewriteSpillover {
        file: PathBuf,
        /// CSV/JSON matrix file to write in. If omitted, the embedded matrix is preserved.
        #[arg(long)]
        matrix: Option<PathBuf>,
        /// Output .fcs path.
        #[arg(short, long)]
        output: PathBuf,
    },

    /// Dump raw + transformed values for one channel (for validation vs flowCore).
    TransformDump {
        file: PathBuf,
        #[arg(long)]
        channel: String,
        #[arg(long, value_enum, default_value = "logicle")]
        method: MethodArg,
        #[arg(long)]
        compensate: bool,
        #[arg(long, default_value_t = 150.0, value_parser = parse_positive)]
        cofactor: f64,
        // logicle params
        #[arg(long, default_value_t = 262144.0)]
        t: f64,
        #[arg(long, default_value_t = 0.5)]
        w: f64,
        #[arg(long, default_value_t = 4.5)]
        m: f64,
        #[arg(long, default_value_t = 0.0)]
        a: f64,
        #[arg(short, long)]
        output: Option<PathBuf>,
    },
}

#[derive(Clone, ValueEnum)]
enum MethodArg {
    Linear,
    Log,
    Asinh,
    Logicle,
}

#[derive(Clone, ValueEnum)]
enum TransformArg {
    None,
    Asinh,
}

/// clap value parser: the asinh cofactor divides the data, so it must be positive and
/// finite (0 → inf/NaN, negative → silently sign-flipped).
fn parse_positive(s: &str) -> Result<f64, String> {
    let v: f64 = s.parse().map_err(|_| format!("`{s}` is not a number"))?;
    if v.is_finite() && v > 0.0 { Ok(v) } else { Err("must be a positive number".into()) }
}

// ── Entry point ───────────────────────────────────────────────────────────

/// On Windows the binary is built for the GUI subsystem (no console), so when it
/// is run from a terminal with CLI arguments we reattach to the parent console
/// so stdout/stderr are visible. No-op elsewhere / for the GUI.
#[cfg(windows)]
fn maybe_attach_console() {
    if std::env::args().len() > 1 {
        use windows_sys::Win32::System::Console::{AttachConsole, ATTACH_PARENT_PROCESS};
        let _ = unsafe { AttachConsole(ATTACH_PARENT_PROCESS) };
    }
}

fn main() -> Result<()> {
    #[cfg(windows)]
    maybe_attach_console();

    // If no subcommand was given, open the GUI directly.
    let args: Vec<String> = std::env::args().collect();
    if args.len() <= 1 {
        gui::run_gui(None);
        return Ok(());
    }

    let cli = Cli::parse();
    match cli.command {
        Command::Gui { file } => {
            gui::run_gui(file.as_deref());
            Ok(())
        }
        Command::Info { file } => cmd_info(&file),
        Command::Update => cmd_update(),
        Command::Selftest => cmd_selftest(),
        Command::Stats { file, compensate, transform, cofactor } =>
            cmd_stats(&file, compensate, &transform, cofactor),
        Command::Export { file, output, compensate, transform, cofactor } =>
            cmd_export(&file, output.as_deref(), compensate, &transform, cofactor),
        Command::Gate { file, gates, compensate, transform, cofactor } =>
            cmd_gate(&file, &gates, compensate, &transform, cofactor),
        Command::Popstats { file, gates, compensate, output } =>
            cmd_popstats(&file, &gates, compensate, output.as_deref()),
        Command::Spillover { file } => cmd_spillover(&file),
        Command::GatingMl { file, gates, compensate, output } =>
            cmd_gatingml(&file, &gates, compensate, &output),
        Command::ComputeSpillover { unstained, controls, output } =>
            cmd_compute_spillover(&unstained, &controls, output.as_deref()),
        Command::RewriteSpillover { file, matrix, output } =>
            cmd_rewrite_spillover(&file, matrix.as_deref(), &output),
        Command::TransformDump { file, channel, method, compensate, cofactor, t, w, m, a, output } =>
            cmd_transform_dump(&file, &channel, &method, compensate, cofactor, t, w, m, a, output.as_deref()),
    }
}

fn cmd_transform_dump(
    path: &Path,
    channel: &str,
    method: &MethodArg,
    do_comp: bool,
    cofactor: f64,
    t: f64, w: f64, m: f64, a: f64,
    output: Option<&std::path::Path>,
) -> Result<()> {
    use transform::AxisTransform;

    let fcs = FcsFile::open(path)?;
    let (events, _label) = prepare_events(&fcs, do_comp, &TransformArg::None, cofactor)?;

    let ci = fcs.param_index(channel)
        .with_context(|| format!("channel '{}' not found", channel))?;

    let axis = match method {
        MethodArg::Linear  => AxisTransform::Linear,
        MethodArg::Log     => AxisTransform::default_log(),
        MethodArg::Asinh   => AxisTransform::Asinh { cofactor },
        MethodArg::Logicle => AxisTransform::Logicle { t, w, m, a },
    };
    let ct = axis.compile();

    let n = fcs.n_params();
    let raw: Vec<f64> = events.iter().skip(ci).step_by(n).copied().collect();

    let mut out: Box<dyn std::io::Write> = match output {
        Some(p) => Box::new(std::io::BufWriter::new(std::fs::File::create(p)?)),
        None => Box::new(std::io::BufWriter::new(std::io::stdout())),
    };
    use std::io::Write;
    // With --compensate the first column holds post-compensation values, not raw —
    // label it honestly so a flowCore comparison isn't fooled.
    let in_label = if do_comp { "compensated" } else { "raw" };
    writeln!(out, "{in_label},transformed")?;
    for v in raw {
        writeln!(out, "{:.6},{:.10}", v, ct.forward(v))?;
    }
    Ok(())
}

// ── Commands ──────────────────────────────────────────────────────────────

fn cmd_selftest() -> Result<()> {
    let results = selftest::run().map_err(|e| anyhow::anyhow!(e))?;
    if selftest::report(&results) {
        Ok(())
    } else {
        anyhow::bail!("selftest failed — a numeric layer deviates from flowCore");
    }
}

fn cmd_update() -> Result<()> {
    match update::check_latest() {
        Ok(info) if info.newer => {
            println!("A newer version is available: v{} (you have v{})", info.latest, info.current);
            println!("Download: {}", info.url);
        }
        Ok(info) => println!("You're on the latest version (v{}).", info.current),
        Err(e) => anyhow::bail!("update check failed: {}", e),
    }
    Ok(())
}

fn cmd_info(path: &Path) -> Result<()> {
    let fcs = FcsFile::open(path)?;
    println!("File    : {}", path.display());
    println!("Version : {}", fcs.version);
    println!("Events  : {}", fcs.n_events);
    println!("Params  : {}", fcs.n_params());

    let instrument = fcs.keywords.get("$CYT").or(fcs.keywords.get("$CYTOMETER"))
        .map(String::as_str).unwrap_or("—");
    println!("Cytometer : {}", instrument);

    let date = fcs.keywords.get("$DATE").map(String::as_str).unwrap_or("—");
    println!("Date    : {}", date);

    println!("\nChannels:");
    for p in &fcs.parameters {
        let label = p.label.as_deref().unwrap_or("");
        println!(
            "  {:2}  {:<20}  stain={:<20}  bits={:2}  range={}",
            p.index, p.name, label, p.bits, p.range
        );
    }

    if let Some(sp) = fcs.spillover_keyword() {
        let preview: String = sp.chars().take(80).collect();
        println!("\n$SPILLOVER (truncated): {}", preview);
        println!("(run `flowcyto spillover {}` for the full matrix)", path.display());
    } else {
        println!("\n$SPILLOVER: not present");
    }
    Ok(())
}

fn cmd_gatingml(path: &Path, gates_path: &Path, compensate: bool, output: &Path) -> Result<()> {
    let fcs = FcsFile::open(path)?;
    // Fluorescence channels (from the embedded $SPILLOVER) get compensation-ref="FCS".
    let comp_channels: Vec<String> = match fcs.spillover_keyword() {
        Some(kw) => parse_spillover(kw).map(|(ch, _)| ch).unwrap_or_default(),
        None => Vec::new(),
    };
    let gates: Vec<Gate> = serde_json::from_str(
        &std::fs::read_to_string(gates_path)
            .with_context(|| format!("reading gates {:?}", gates_path))?,
    )
    .context("parsing gates JSON")?;

    let (xml, warnings) = gatingml::to_gating_ml(&gates, &comp_channels, compensate);
    std::fs::write(output, &xml).with_context(|| format!("writing {:?}", output))?;
    println!("Wrote {} gate(s) → {} (Gating-ML 2.0)", gates.len(), output.display());
    if compensate && comp_channels.is_empty() {
        println!("  note: --compensate set but the file has no $SPILLOVER; dimensions are uncompensated");
    }
    for w in &warnings {
        println!("  ⚠ {}", w);
    }
    Ok(())
}

fn cmd_spillover(path: &Path) -> Result<()> {
    let fcs = FcsFile::open(path)?;
    println!("File : {}", path.display());

    let kw = match fcs.spillover_keyword() {
        Some(k) => k,
        None => {
            println!("\nNo $SPILLOVER / $SPILL keyword present.");
            println!("→ This file is uncompensated; the Compensate option would do nothing.");
            return Ok(());
        }
    };

    let (channels, rows) = parse_spillover(kw)?;
    let n = channels.len();

    // stain label per channel (from $PnS)
    let stain = |name: &str| -> String {
        fcs.parameters.iter()
            .find(|p| p.name.eq_ignore_ascii_case(name))
            .and_then(|p| p.label.clone())
            .filter(|l| !l.is_empty())
            .map(|l| format!(" ({})", l))
            .unwrap_or_default()
    };

    println!("\n$SPILLOVER matrix — {} fluorescence channels", n);
    println!("Row = source dye, column = detector. Diagonal = 1.0 by definition.\n");

    println!("Detector legend:");
    for (i, c) in channels.iter().enumerate() {
        println!("  [{}] {}{}", i + 1, c, stain(c));
    }

    // header
    print!("\n{:<22}", "Source \\ Detector");
    for i in 0..n { print!("{:>9}", format!("[{}]", i + 1)); }
    println!();

    for (i, row) in rows.iter().enumerate() {
        let label = format!("[{}] {}", i + 1, channels[i]);
        print!("{:<22}", truncate(&label, 22));
        for (j, &v) in row.iter().enumerate() {
            if i == j {
                print!("{:>9}", "1");
            } else {
                print!("{:>9.4}", v);
            }
        }
        println!();
    }

    let mx = max_off_diagonal(&rows);
    println!();
    if mx < 1e-9 {
        println!("Max off-diagonal spillover: 0.0000 → identity matrix (NO real compensation).");
    } else {
        println!("Max off-diagonal spillover: {:.4} ({:.1}%) → real compensation matrix.", mx, mx * 100.0);
    }

    // provenance
    let get = |k: &str| fcs.keywords.get(k).map(String::as_str);
    let creator = get("CREATOR").or(get("$CREATOR")).or(get("APPLICATION"));
    let cyt = get("$CYT").or(get("$CYTOMETER"));
    let date = get("$DATE");
    let mut prov: Vec<String> = Vec::new();
    if let Some(c) = creator { prov.push(c.to_string()); }
    if let Some(c) = cyt { prov.push(c.to_string()); }
    if let Some(d) = date { prov.push(d.to_string()); }
    if !prov.is_empty() {
        println!("Created by: {}", prov.join(" · "));
    }
    Ok(())
}

fn truncate(s: &str, n: usize) -> String {
    if s.chars().count() <= n { s.to_string() }
    else { s.chars().take(n.saturating_sub(1)).collect::<String>() + "…" }
}

fn cmd_popstats(
    path: &Path,
    gates_path: &Path,
    do_comp: bool,
    output: Option<&std::path::Path>,
) -> Result<()> {
    let fcs = FcsFile::open(path)?;
    // Compensated linear data — the MFI space (no display transform).
    let (events, label) = prepare_events(&fcs, do_comp, &TransformArg::None, 0.0)?;

    let gates: Vec<gating::Gate> = serde_json::from_str(
        &std::fs::read_to_string(gates_path)
            .with_context(|| format!("reading gates {:?}", gates_path))?,
    ).context("parsing gates JSON")?;

    // Summarize all channels except Time.
    let stat_channels: Vec<usize> = fcs.parameters.iter().enumerate()
        .filter(|(_, p)| !p.name.eq_ignore_ascii_case("Time"))
        .map(|(i, _)| i)
        .collect();

    let table = popstats::population_stats(&events, &fcs.parameters, fcs.n_events, &gates, &stat_channels);

    println!("── {} ─── {}", path.display(), label);
    println!("{:<28} {:>9} {:>9} {:>9}", "Population", "Count", "%Parent", "%Total");
    println!("{}", "─".repeat(58));
    for r in &table.rows {
        let indent = "  ".repeat(r.depth);
        let name = format!("{}{}", indent, r.name);
        println!("{:<28} {:>9} {:>8.2}% {:>8.2}%", truncate(&name, 28), r.count, r.pct_parent, r.pct_total);
    }

    if let Some(o) = output {
        let mut s = String::new();
        s.push_str(popstats::LONG_CSV_HEADER);
        s.push('\n');
        let sample = path.file_stem().map(|x| x.to_string_lossy().to_string())
            .unwrap_or_else(|| "sample".into());
        popstats::append_long_csv(&mut s, &sample, &table);
        std::fs::write(o, s).with_context(|| format!("writing {:?}", o))?;
        println!("\nWrote tidy per-population stats → {}", o.display());
    }
    Ok(())
}

fn cmd_compute_spillover(
    unstained: &std::path::Path,
    controls: &[PathBuf],
    output: Option<&std::path::Path>,
) -> Result<()> {
    let unst = FcsFile::open(unstained)?;
    let ctrl_files: Vec<FcsFile> = controls.iter()
        .map(|p| FcsFile::open(p))
        .collect::<Result<_>>()?;

    // Fluorescence channels (exclude FSC/SSC/Time), from the unstained control.
    let fluor_idx = transform::fluorescence_indices(&unst.parameters);
    let fluor: Vec<String> = fluor_idx.iter().map(|&i| unst.parameters[i].name.clone()).collect();

    let refs: Vec<&FcsFile> = ctrl_files.iter().collect();
    let comp = compute_spillover(&fluor, &unst, &refs)?;

    println!("Computed spillover from {} single-stain controls + unstained", controls.len());
    println!("(positive-gated median, background-subtracted, row-normalized by primary)\n");

    // matrix
    print!("{:<22}", "Source \\ Detector");
    for i in 0..fluor.len() { print!("{:>9}", format!("[{}]", i + 1)); }
    println!();
    for (i, row) in comp.rows.iter().enumerate() {
        print!("{:<22}", truncate(&format!("[{}] {}", i + 1, fluor[i]), 22));
        for (j, &v) in row.iter().enumerate() {
            if i == j { print!("{:>9}", "1"); } else { print!("{:>9.4}", v); }
        }
        println!();
    }

    // Stain assignment + filename correctness guard
    println!("\nStain assignment (intensity-based) vs filename:");
    let mut all_ok = true;
    for (ci, &p) in comp.assigned.iter().enumerate() {
        let fname = controls[ci].file_name().map(|s| s.to_string_lossy().to_string()).unwrap_or_default();
        let by_name = fluor_token_in_filename(&fluor, &fname);
        let ok = by_name == Some(p);
        if !ok { all_ok = false; }
        println!("  {:<55} → {:<16} {}",
            truncate(&fname, 55), fluor[p],
            if ok { "✓" } else { "⚠ MISMATCH vs filename" });
    }
    if all_ok {
        println!("✓ All controls' brightest channel matches their filename fluorophore.");
    } else {
        println!("⚠ At least one control's brightest channel disagrees with its filename — \
                  check for swapped/mislabeled controls before trusting this matrix.");
    }

    let mx = max_off_diagonal(&comp.rows);
    println!("\nMax off-diagonal spillover: {:.4} ({:.1}%)", mx, mx * 100.0);

    if let Some(o) = output {
        save_matrix_file(o, &comp.channels, &comp.rows)?;
        println!("\nWrote matrix → {}", o.display());
        println!("Apply it with:  flowcyto rewrite-spillover <data.fcs> --matrix {} -o <out.fcs>", o.display());
    }
    Ok(())
}

fn cmd_rewrite_spillover(
    path: &Path,
    matrix: Option<&std::path::Path>,
    output: &std::path::Path,
) -> Result<()> {
    let fcs = FcsFile::open(path)?;

    let new_spill: Option<String> = match matrix {
        Some(mp) => {
            let (channels, rows) = load_matrix_file(mp)?;
            // Validate every matrix channel exists in the file.
            for c in &channels {
                if fcs.param_index(c).is_none() {
                    anyhow::bail!(
                        "matrix channel '{}' not found in FCS parameters: [{}]",
                        c,
                        fcs.parameters.iter().map(|p| p.name.as_str()).collect::<Vec<_>>().join(", ")
                    );
                }
            }
            // Sanity: ensure it's invertible (else compensation downstream fails).
            SpilloverMatrix::from_parts(channels.clone(), &rows)
                .context("supplied matrix cannot be inverted")?;
            println!("Override matrix: {} channels from {}", channels.len(), mp.display());
            Some(format_spillover(&channels, &rows))
        }
        None => {
            println!("No matrix supplied — preserving the embedded $SPILLOVER.");
            None
        }
    };

    fcs_write::write_fcs(&fcs, new_spill.as_deref(), output)?;
    println!("Wrote {} events × {} channels → {}", fcs.n_events, fcs.n_params(), output.display());
    println!("(original file unchanged)");
    Ok(())
}

fn cmd_stats(
    path: &Path,
    do_comp: bool,
    transform: &TransformArg,
    cofactor: f64,
) -> Result<()> {
    let fcs = FcsFile::open(path)?;
    let (events, label) = prepare_events(&fcs, do_comp, transform, cofactor)?;

    println!("── {} ─── {}", path.display(), label);
    let stats: Vec<Stats> = fcs
        .parameters
        .iter()
        .enumerate()
        .map(|(i, p)| {
            let vals: Vec<f64> = events
                .iter()
                .skip(i)
                .step_by(fcs.n_params())
                .copied()
                .collect();
            Stats::compute(&p.name, &vals)
        })
        .collect();

    print_stats_table(&stats);
    Ok(())
}

fn cmd_export(
    path: &Path,
    output: Option<&std::path::Path>,
    do_comp: bool,
    transform: &TransformArg,
    cofactor: f64,
) -> Result<()> {
    let fcs = FcsFile::open(path)?;
    let (events, label) = prepare_events(&fcs, do_comp, transform, cofactor)?;

    let out_path = output.map(|p| p.to_path_buf()).unwrap_or_else(|| {
        let mut p = path.to_path_buf();
        p.set_extension("csv");
        p
    });

    // Never write over the input — e.g. an FCS file named `foo.csv`, where
    // set_extension("csv") is a no-op and the default output would clobber the source.
    let same_as_input = out_path.as_path() == path
        || matches!((out_path.canonicalize(), path.canonicalize()), (Ok(a), Ok(b)) if a == b);
    if same_as_input {
        anyhow::bail!(
            "refusing to overwrite the input file {:?} — pass -o to choose a different output path",
            path
        );
    }

    let mut wtr = csv::Writer::from_path(&out_path)
        .with_context(|| format!("cannot create {:?}", out_path))?;

    let headers: Vec<&str> = fcs.parameters.iter().map(|p| p.name.as_str()).collect();
    wtr.write_record(&headers)?;

    let n = fcs.n_params();
    for ev in 0..fcs.n_events {
        let base = ev * n;
        let row: Vec<String> = events[base..base + n]
            .iter()
            .map(|v| if v.is_finite() { format!("{:.4}", v) } else { "NA".to_string() })
            .collect();
        wtr.write_record(&row)?;
    }
    wtr.flush()?;

    println!("Exported {} events × {} channels → {:?}  ({})",
        fcs.n_events, fcs.n_params(), out_path, label);
    Ok(())
}

fn cmd_gate(
    path: &Path,
    gates_path: &Path,
    do_comp: bool,
    transform: &TransformArg,
    cofactor: f64,
) -> Result<()> {
    let fcs = FcsFile::open(path)?;
    let (events, label) = prepare_events(&fcs, do_comp, transform, cofactor)?;

    let gates_json = std::fs::read_to_string(gates_path)
        .with_context(|| format!("cannot read gates file {:?}", gates_path))?;
    let gates: Vec<Gate> = serde_json::from_str(&gates_json)
        .context("failed to parse gates JSON")?;

    let results = apply_gates(&gates, &events, &fcs.parameters, fcs.n_events)?;

    println!("── {} ─── {} ─── gates: {}", path.display(), label, gates_path.display());
    println!("{:<24} {:>10} {:>10} {:>10} {:>9}", "Gate", "In", "Parent", "%Parent", "%Total");
    println!("{}", "─".repeat(67));
    for r in &results {
        println!("{:<24} {:>10} {:>10} {:>9.2}% {:>8.2}%",
            r.name, r.n_in, r.n_parent, r.pct_parent(), r.pct_total());
    }
    Ok(())
}

// ── Shared pipeline ───────────────────────────────────────────────────────

fn prepare_events(
    fcs: &FcsFile,
    do_comp: bool,
    transform: &TransformArg,
    cofactor: f64,
) -> Result<(Vec<f64>, String)> {
    let mut events = fcs.events.clone();
    let mut stages: Vec<&str> = vec!["raw"];

    if do_comp {
        let kw = fcs.spillover_keyword()
            .context("--compensate requested but no $SPILLOVER/$SPILL keyword in FCS file")?;
        let spill = SpilloverMatrix::from_keyword(kw)?;
        events = spill.apply(fcs)?;
        stages.push("compensated");
    }

    if matches!(transform, TransformArg::Asinh) {
        let fluor_idx = fluorescence_indices(&fcs.parameters);
        apply_asinh(&mut events, fcs.n_params(), &fluor_idx, cofactor);
        stages.push("asinh");
    }

    Ok((events, stages.join(" → ")))
}
