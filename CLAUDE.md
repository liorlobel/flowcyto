# flowcyto — developer guide & session handoff

A Rust **CLI + GUI** for analyzing BD flow-cytometry **.fcs** files (FCS 2.0/3.0/3.1).
Built incrementally and **validated against R/flowCore at every numeric layer**. ~5k LOC.

## Build & run  ⚠️ cargo is NOT on the normal PATH
The `~/.cargo/bin/cargo` symlink is broken. Use one of:
```bash
export PATH="/Users/liorlobel/.rustup/toolchains/stable-aarch64-apple-darwin/bin:$PATH"
cd /Users/liorlobel/flowcyto && /opt/homebrew/bin/rustup run stable cargo build --release
```
- Binary: `target/release/flowcyto`
- GUI: `flowcyto gui <file.fcs>`  (or `flowcyto` with no args → GUI)
- **macOS installer:** `./packaging/make-macos-app.sh` → `dist/flowcyto.app` + `dist/flowcyto-<version>.dmg` (drag-to-Applications). Builds host-arch (Apple Silicon), generates the `.icns` from `packaging/icon.png`, writes Info.plist with the Cargo version, ad-hoc code-signs. Not notarized → first launch is right-click → Open. `dist/` is git-ignored.
- Always finish with: `cargo build --release` clean, `cargo clippy --release --all-targets` = **0 warnings**, `cargo test --release` = **88 tests pass** (84 unit + 4 CLI integration). Unit tests live inline (`#[cfg(test)] mod tests`) in each module; `src/test_util.rs` is a `cfg(test)`-only in-memory `FcsFile` builder; `tests/cli.rs` drives the real binary against `tests/fixtures/tiny.fcs`. Add a regression test alongside any numeric change.

## Architecture (src/)
| file | role |
|---|---|
| `fcs.rs` | FCS parser (DATATYPE F/D/I, BYTEORD, offsets, `peek_events` for QC) |
| `fcs_write.rs` | FCS 3.0 writer (raw events as F-LE; writes `$SPILLOVER`+`SPILL`) |
| `compensation.rs` | spillover: parse / apply (M⁻¹) / `compute_spillover` from single stains / matrix CSV·JSON IO |
| `transform.rs` | `AxisTransform` (Linear/Log/Asinh/Logicle) + `CompiledTransform` (forward/inverse) |
| `logicle.rs` | Moore & Parks logicle (faithful port; `scale`/`inverse`) |
| `gating.rs` | `Gate` + `GateShape` (Rect/Ellipse/Polygon/Range), hierarchical `effective_mask`, `gate_tree_order` |
| `popstats.rs` | **pure** per-population stats engine (count/%parent/%total/median-MFI/mean/CV) — also the batch engine |
| `stats.rs` | per-channel whole-file stats (CLI `stats`) |
| `gui.rs` | egui GUI (2.6k LOC) — tabs Plot/Histogram/Stats/Batch/Spillover |
| `main.rs` | clap CLI |

**CLI:** `info stats export gate popstats spillover compute-spillover rewrite-spillover transform-dump gui`

**GUI:** left panel = Samples (QC event counts, 👁 overlay) · Channels (X/Y + per-axis transform) · Axis limits · Gates (draw ▭⬭⬠ ✛Quad ✎Edit, tree, ▶ gate-from-here, numeric inspector, save/load). Tabs: Plot (density + contours + gates + control overlay + backgate), Histogram (overlays: populations or samples), Stats (per-population table + CSV), Batch (streamed multi-sample → CSV), Spillover (view/edit/compute/write matrix).

## Validation discipline (THE most important habit — keep it)
R + **flowCore 2.24.0** are installed; flowCore is the oracle. Validate every numeric change before building on it:
- counts must match **EXACTLY**; medians/transforms to ~1e-5..1e-11 (float rounding).
- Harnesses: `flowcyto transform-dump`, `flowcyto popstats -o csv`, `flowcyto spillover`, `flowcyto compute-spillover`.
- Proven equal to flowCore: parsing, compensation, asinh, logicle, per-population counts+MFI, spillover-from-controls (also synthetic ground-truth recovery), FCS writer round-trip (0e+00), quadrant partition.

## GUI dev workflow (how features were validated)
- Screenshot: launch GUI in background → `osascript -e 'tell application "System Events" to set frontmost of (first process whose name contains "flowcyto") to true'` → `screencapture -x /tmp/x.png` → Read the PNG.
- **Temp-screenshot pattern:** temporarily edit the `// TEMP` line in `run_gui` (bottom of gui.rs) to preset state (channels, gates, toggles), screenshot, then **revert** — always `grep -n TEMP src/gui.rs` before declaring done.
- Borrow checker: clone render data OUT before `plot.show(...)` closures; capture egui `Response` booleans as owned values before any `pu.line/points/polygon` call.
- Caches (`scatter`/`pop_stats`/`hist_cache`/`ref_scatter`/`gate_counts`) invalidate via `None` + `needs_reprocess/regate/rescatter`. There was a frame-ordering class of bug — keep `compensated` consistent with `fcs` before any panel renders.

## Status — feature roadmap
FlowJo-parity: ✅ per-population stats · 1D histograms+overlays · multi-sample batch. UI bundles: ✅ A (quadrant + numeric + drag-resize gates), ✅ B (%/count labels, gate-from-here, backgating), ✅ C (control overlay + per-tube QC), ◑ D (✅ contours+legend; **TODO: multi-plot layout grid**, zebra plots).
Other possible: boolean gates, subset-FCS export, per-tube %viable QC scan (needs a Live gate), tSNE/UMAP/FlowSOM.

## The real analysis done with it (cDC in cAPC/SAA-diet experiments)
Data: `…/cAPC_SAA_Diets/*_cAPC_mice_myeloid_stain/` (4 usable experiments; 02_12_19 excluded — no controls). Panel: FITC=CD11c, PE=CD103, PerCP-Cy5-5=CD11b, PE-Cy7=MHCII, PacBlue=CD45, AmCyan=Live/Dead. cDC = CD11c⁺MHCII⁺; cDC1 = CD103⁺CD11b⁻; cDC2 = CD103⁻CD11b⁺.
**Gotchas that bit us:** (1) compensation control tubes are stored **uncompensated (identity matrix)** — compensate samples with the **sample's own** embedded `$SPILLOVER`, not the unstained's. (2) cDC1/cDC2 is **compensation-sensitive** (embedded under-corrects MHCII→CD11b: 0.145 vs 0.325 from single stains) — only 06_20_18 has single stains. (3) **mice are the replicates** (pool them, experiment as a fixed block); a binomial GLMM on cell counts pseudoreplicates → spurious p<0.0001 (use sample-level / OLRE). (4) flag bad tubes (one had 141 CD45 events, 4.5% viable).
**Findings:** tumor cDC1-depletion robust (4/4 experiments); high-SAA diet ↑ colonic cDC (p≈0.007, mice pooled, blocked); no MLN diet effect. Results CSVs saved in the experiment folders (`cDC_*_results.csv`, `cDC_QC_table.csv`).

> More detail in the user's auto-memory `project_flowcyto.md`. NOTE: the multi-experiment cDC gating + all statistics were run in **R/flowCore** (flowcyto spot-checked == flowCore); flowcyto did the compensation-from-single-stains + validated single-experiment gating.
