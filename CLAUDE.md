# flowcyto — developer guide & session handoff

A Rust **CLI + GUI** for analyzing BD flow-cytometry **.fcs** files (FCS 2.0/3.0/3.1).
Built incrementally and **validated against R/flowCore at every numeric layer**. ~6.5k LOC.

## Build & run  ⚠️ cargo is NOT on the normal PATH
The `~/.cargo/bin/cargo` symlink is broken. Use one of:
```bash
export PATH="/Users/liorlobel/.rustup/toolchains/stable-aarch64-apple-darwin/bin:$PATH"
cd /Users/liorlobel/flowcyto && /opt/homebrew/bin/rustup run stable cargo build --release
```
- Binary: `target/release/flowcyto`
- GUI: `flowcyto gui <file.fcs>`  (or `flowcyto` with no args → GUI)
- **macOS installer:** `./packaging/make-macos-app.sh` → `dist/flowcyto.app` + `dist/flowcyto-<version>.dmg` (drag-to-Applications). Builds host-arch (Apple Silicon), generates the `.icns` from `packaging/icon.png` (itself regenerable via `python3 packaging/make-icon.py` — a procedural viridis density-plot + gate-ring mark), writes Info.plist with the Cargo version, ad-hoc code-signs, and bundles an "Open Me First.txt" Gatekeeper guide into the DMG. Not notarized (needs paid Apple Developer Program) → recipients clear Gatekeeper once via `xattr -dr com.apple.quarantine /Applications/flowcyto.app` or "Open Anyway"; see `INSTALL.md`. `dist/` is git-ignored.
- **Releasing:** bump `version` in `Cargo.toml`, commit, then push a `vX.Y.Z` tag. Two GitHub Actions workflows (`.github/workflows/macos-installer.yml`, `windows-installer.yml`) build on `macos-latest` (clippy+test gate → `.dmg`) and `windows-latest` (test → Inno Setup `.exe`) and attach both installers to the release (auto-created with `--generate-notes`; edit notes after if desired). Both also run on manual dispatch (artifact only). No need to build installers locally anymore.
- Always finish with: `cargo build --release` clean, `cargo clippy --release --all-targets` = **0 warnings**, `cargo test --release` = **109 tests pass** (105 unit + 4 CLI integration). Unit tests live inline (`#[cfg(test)] mod tests`) in each module; `src/test_util.rs` is a `cfg(test)`-only in-memory `FcsFile` builder; `tests/cli.rs` drives the real binary against `tests/fixtures/tiny.fcs`. Add a regression test alongside any numeric change.
- Current released version: **0.1.11** (latest GitHub release; macOS `.dmg` + Windows `.exe`, unsigned).

## Architecture (src/)
| file | role |
|---|---|
| `fcs.rs` | FCS parser (DATATYPE F/D/I, BYTEORD, offsets, `peek_events` for QC) |
| `fcs_write.rs` | FCS 3.0 writer (raw events as F-LE; writes `$SPILLOVER`+`SPILL`) |
| `compensation.rs` | spillover: parse / apply (M⁻¹) / `compute_spillover` from single stains / matrix CSV·JSON IO |
| `transform.rs` | `AxisTransform` (Linear/Log/Asinh/Logicle) + `CompiledTransform` (forward/inverse) |
| `logicle.rs` | Moore & Parks logicle (faithful port; `scale`/`inverse`) |
| `gating.rs` | `Gate` + `GateShape` (Rect/Ellipse/Polygon/Range/**Boolean**), hierarchical `effective_mask`, `gate_tree_order`, **`compute_own_masks`** (the one place that builds every gate's own mask — geometric gates then Boolean AND/OR/NOT combos in dependency order; all views + CLI go through it) |
| `popstats.rs` | **pure** per-population stats engine (count/%parent/%total/median-MFI/mean/CV) — also the batch engine |
| `stats.rs` | per-channel whole-file stats (CLI `stats`) |
| `gui.rs` | egui GUI (~4.5k LOC) — tabs Plot/Histogram/Stats/Batch/Spillover; native macOS menu bar via `muda` (cfg-gated) |
| `update.rs` | manual "check for updates" vs the GitHub Releases API (version compare; the app's **only** network call, fired only on explicit user action) |
| `qc.rs` | **pure** acquisition-QC metrics — flow-rate (clog) stability over the Time channel + margin/saturation events (unit-tested on synthetic known-answer data); %viable comes from the gating engine |
| `main.rs` | clap CLI |

**CLI:** `info update stats export gate popstats spillover compute-spillover rewrite-spillover transform-dump gui`

**GUI:** left panel = Samples (QC counts, 👁 overlay, **group/condition tags**) · Channels (X/Y + per-axis scale, "apply X scale to all fluorescence") · Axis limits · Gates. Toolbar: Open/**Recent**/**Save+Load session**/Compensate/theme/tabs. Gating: draw ▭⬭⬠ ✛Quad ✎Edit (drag body to move, rotate ellipse), **double-click a gate to drill in**, per-gate **👁 hide** + **⊕ zoom-to-gate**, **➕ Boolean (AND/OR/NOT)** builder, **undo/redo**, numeric inspector, save/load JSON, **export a population → .fcs**. Tabs: Plot (density dots or **filled heatmap "Fill"**, contours, gates, control overlay, backgate, **🔒 Lock view** = frozen pan/zoom, **adjustable Single / cols×rows grid up to 6×6**, Viridis/Jet colormap, **📷 Save plot PNG**, inline ⚖ compensation preview), Histogram (overlays + interval gates), Stats (table + CSV + **📋 Copy TSV**), Batch (threaded multi-sample → CSV + **📋 Copy** + **📊 chart across samples**), Spillover (view/edit/compute/write matrix). **Drag-and-drop .fcs** to open; keyboard shortcuts (R/E/P/Q/G/V/Esc, ⌘Z/⌘S/⌘1–5) + ⌘+/− UI zoom.

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

## QC suite + %viable (2026-06-08, v0.1.11) — roadmap item #1–2
New **QC tab** (6th tab; `ActiveTab::QC`, ⌘6). Streams every workspace sample on a worker thread (same pattern as `run_batch` — `run_qc`/`poll_qc`/`QcMsg`) and reports per tube: event count; **%viable** = a user-designated **live gate** (`qc_live_gate: Option<u32>`, reconciled in `after_gate_restore`/`load_session`/on-delete like `hist_sample_pop`) computed as its effective count ÷ parent count via the validated gating engine; **flow-rate stability** (`qc::flow_rate` — bins events over the **Time** channel; a clog/bubble = a dip below the median bin rate; rendered as an inline `flow_sparkline`; **None** when Time is absent/constant, shown as "n/a" — never faked from event order); **margin/off-scale** events (`qc::margin_events` — % at `$PnR` ceiling, on RAW data). Flagged tubes highlighted amber; CSV export. `qc.rs` is pure + unit-tested (clog/uniform/constant-Time/margin known-answer). Thresholds: `QC_MIN_VIABLE`/`QC_FLOW_DEV`/`QC_MARGIN_PCT`. NOTE (design, per review): channel-median drift would NOT catch clogs (the cells that pass look normal) — flow-rate (events/time) is the real clog detector. `gen_fcs.py` (screenshot helper, in /tmp) now emits a Time channel with an injected clog gap.

## Re-audit fixes (2026-06-07, v0.1.10) — see `AUDIT.md` "v0.1.9 re-audit"
Fresh full audit of v0.1.9 → all 7 Lows fixed + 2 follow-ups. **N1:** `compensation::median` (the `compute_spillover` path) now drops non-finite — the one median site the v0.1.7 M1d fix missed; `save_matrix_file` also refuses non-finite. **N2:** `update::check_latest` no longer trusts the response `html_url` — always opens the constant `RELEASES_PAGE` (don't hand a network-controlled string to `open::that`). **N3:** Windows CI passes the tag via `env:`/`$env:TAG` (was `${{ github.ref_name }}` in PowerShell source — injection sink). **N4:** `pop_mask` got the `compensated.len()` buffer guard its siblings have. **N5:** integer bitmask uses `checked_next_power_of_two` ($PnR>2^63). **N6:** finished the Phosphor migration (~13 leftover emoji → icons; the "every emoji" claim is true again). **N7:** CI actions SHA-pinned + Windows now runs the clippy `-D warnings` gate. **Follow-ups:** added `assets/fonts/Phosphor-MIT.txt`; replaced the hand-typed Viridis/Cividis LUTs (off by up to 61/255) with reference matplotlib values (Δ≤1, verified vs `_cm_listed.py`). All prior fixes confirmed intact; `cargo audit` = 0 vulns.

## UI redesign (2026-06-07, v0.1.9)
Visual pass over the egui GUI (no numeric changes). **Typography:** bundles **Inter** (Regular + SemiBold for headings) via `setup_fonts` + `assets/fonts/` (`include_bytes!`; OFL license alongside). **Icons:** `egui-phosphor` (0.7.3, the egui-0.29 line) replaces every UI emoji — aliased `use egui_phosphor::regular as icon;`, used as `format!("{} Label", icon::NAME)`. **Theme:** `themed_visuals` deepened — cohesive panel/window fills, soft `window_shadow`/`popup_shadow`, accent-tinted hover/active widget states, teal active tab (both dark + light). **Colormaps:** `ColorMap` now Viridis/Magma/Turbo/Cividis/Jet (`lerp_anchors` over 9-anchor LUTs + `turbo()` polynomial in f64 to avoid the excessive-precision lint); persisted by name in `Session.colormap` (back-compat with the old `viridis` bool). **Empty-state:** `empty_state()` (centered icon + heading + Open button). New deps: `egui-phosphor`, plus the two Inter TTFs (~650 KB) in `assets/`. Kept the neutral ■/▶/› glyphs. Gate-label "chips" intentionally skipped (egui_plot text has no background). Verified: build/clippy(-D)/109 tests green + screenshots on a synthetic 5-population FCS.

## Audit & hardening pass (2026-06-07) — see `AUDIT.md` for the full record
Full-codebase audit (6 parallel review passes, every finding re-verified at its `file:line`), then **all High/Medium + 10/12 Low findings fixed** with 16 new regression tests. Numerics re-validated unchanged (the exact-count + median tests still pass). Highlights:
- **(HIGH) Boolean `NOT`/`OR` of a missing-channel gate** reported the *parent's* count (~100%) / silently dropped the ref — reintroducing the geometric "reports parent count" bug via the Boolean path. Fixed in `gating::compute_own_masks` by tracking **unevaluable** gate ids (missing channel) distinctly from *empty*: a Boolean referencing an unevaluable population is itself unevaluable (all-false). NOT of a *legitimately empty* gate still correctly returns all-true.
- **(MED) Silent non-finite / ill-conditioned compensation** → poisoned MFI. `parse_spillover` now errors on malformed/`nan`/`inf` tokens (was `unwrap_or(0.0)`); `from_parts`/`validate_square` reject non-finite; `from_parts` rejects ill-conditioned inverses (max|inv|>1e6, catches near-singular that nalgebra inverts to ~1e13); `med_mean_cv`/`Stats::compute` drop non-finite before reducing.
- **(MED) FCS parser DoS** — `$PAR×$TOT` is now `checked_mul` + bounded by the DATA segment size (was an unbounded `Vec::with_capacity` → OOM/overflow); `$PAR` capped; `parse_spillover` `n` bounded (was overflow→slice-panic).
- **(MED) `export` data loss** — refuses to overwrite the input (e.g. an FCS named `foo.csv`).
- **Lows:** CLI `--cofactor` must be >0; `transform-dump` labels compensated input honestly; CSV non-finite → `NA`; `Log` floor clamped >0; supplemental TEXT reads the standard `$BEGINSTEXT`/`$ENDSTEXT` (+ bound); exotic `$BYTEORD` permutations rejected; `peek_events` bounded; Stats-tab buffer guard; `load_session` validates the override matrix + reconciles stale gate-ids; batch worker panic no longer reported as success; quadrant gates half-open (no center double-count); rotated-ellipse zoom bbox; Samples panel shown for a single sample.
- **Intentionally left:** L4 (compensate CLI-error vs GUI-no-op divergence — each reasonable) and L6 (logicle→Linear fallback is correct/tested).

## Status — feature roadmap
FlowJo-parity ✅: per-population stats · 1D histograms+overlays · multi-sample batch · quadrant/numeric/drag-resize/rotate gates · %/count labels · gate-from-here + double-click drill · backgating · control overlay + per-tube QC · contours+legend · **multi-plot grid (up to 6×6)** · **boolean (AND/OR/NOT) gates** · **subset-FCS export** · filled-density heatmap · undo/redo · sessions · clipboard/recent/drag-drop · native menu bar · cross-platform CI installers (macOS + Windows).
Still open / ideas: per-tube %viable QC scan (needs a Live gate), zebra plots, code-signing+notarization (needs paid Apple + Windows certs), universal/Intel macOS build, tSNE/UMAP/FlowSOM.

## The real analysis done with it (cDC in cAPC/SAA-diet experiments)
Data: `…/cAPC_SAA_Diets/*_cAPC_mice_myeloid_stain/` (4 usable experiments; 02_12_19 excluded — no controls). Panel: FITC=CD11c, PE=CD103, PerCP-Cy5-5=CD11b, PE-Cy7=MHCII, PacBlue=CD45, AmCyan=Live/Dead. cDC = CD11c⁺MHCII⁺; cDC1 = CD103⁺CD11b⁻; cDC2 = CD103⁻CD11b⁺.
**Gotchas that bit us:** (1) compensation control tubes are stored **uncompensated (identity matrix)** — compensate samples with the **sample's own** embedded `$SPILLOVER`, not the unstained's. (2) cDC1/cDC2 is **compensation-sensitive** (embedded under-corrects MHCII→CD11b: 0.145 vs 0.325 from single stains) — only 06_20_18 has single stains. (3) **mice are the replicates** (pool them, experiment as a fixed block); a binomial GLMM on cell counts pseudoreplicates → spurious p<0.0001 (use sample-level / OLRE). (4) flag bad tubes (one had 141 CD45 events, 4.5% viable).
**Findings:** tumor cDC1-depletion robust (4/4 experiments); high-SAA diet ↑ colonic cDC (p≈0.007, mice pooled, blocked); no MLN diet effect. Results CSVs saved in the experiment folders (`cDC_*_results.csv`, `cDC_QC_table.csv`).

> More detail in the user's auto-memory `project_flowcyto.md`. NOTE: the multi-experiment cDC gating + all statistics were run in **R/flowCore** (flowcyto spot-checked == flowCore); flowcyto did the compensation-from-single-stains + validated single-experiment gating.
