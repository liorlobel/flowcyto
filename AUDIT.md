# flowcyto — code audit (v0.1.6)

> **Status (2026-06-07): all High/Medium findings and ten of twelve Lows are FIXED**
> in the working tree (16 new regression tests; 102 unit + 4 CLI tests pass, clippy
> clean). The two exceptions are deliberate, documented non-changes: **L4** (CLI errors
> vs GUI no-ops on `--compensate` with no matrix — each behavior is individually
> reasonable; the GUI already comments its no-op) and **L6** (the logicle→Linear
> fallback is correct and tested; surfacing it cleanly needs a UI change not worth the
> risk for a Low). This file is retained as the audit record.

---

# v0.1.9 re-audit (2026-06-07) — addendum

> **Status: all of N1–N7 FIXED in v0.1.10**, plus the two follow-ups (Phosphor MIT
> license added next to the Inter OFL; the hand-typed Viridis/Cividis colormap LUTs —
> which were off by up to 61/255 — replaced with reference matplotlib values, Δ≤1).
> 109 tests pass, clippy `-D warnings` clean. N6 also makes the "replaces every UI
> emoji" claim true again.

Fresh full audit of the current tree (**v0.1.9, commit `3a848e2`, ~8.9k LOC**) after
three releases of change since the v0.1.6 audit above: the hardening fixes (v0.1.7),
the network-touching in-app updater (`update.rs`, v0.1.8), and a heavy `gui.rs`
redesign (+480 lines: fonts/icons/theme/colormaps, v0.1.9). Five parallel review
passes — update/network · numeric-core edges · gui state/regression · gui
rendering · CLI/deps/CI — **every finding re-verified at its `file:line`.** Static
analysis (no exploit run). The agents cross-checked against the fixed list above so
this is *new since v0.1.6 + regressions*, not a rehash.

**Headline: the codebase is in good shape.** All prior fixes survived the churn —
**5/5 GUI fixes confirmed intact**, **5/5 v0.1.7 numeric fixes verified sound** on the
untested edges, and **`cargo audit` reports 0 vulnerabilities** (the new pure-Rust
network stack is clean). No High/Medium. Seven Lows (one Low–Medium), summarized:

| # | Severity | Finding | Location |
|---|---|---|---|
| N1 | Low–Med | `compute_spillover`'s `median()` doesn't drop non-finite — the one median site the v0.1.7 M1d fix **missed**; CLI `compute-spillover -o` can silently write a NaN/wrong matrix to disk | `compensation.rs:271-282` |
| N2 | Low | `open::that(&info.url)` opens a **response-controlled** URL (`html_url` from the API) — the app's only attacker-influenceable sink | `update.rs:53-57` → `gui.rs:1443` |
| N3 | Low | **Windows** CI interpolates `${{ github.ref_name }}` into PowerShell **source** (script-injection sink with a `contents:write` token); macOS uses the safe env-var form | `windows-installer.yml:58` |
| N4 | Low (latent) | `pop_mask` / `grid_kept` omit the `compensated.len()` guard their siblings have (same class as the fixed L8) — latent OOB panic, not reachable on current frame ordering | `gui.rs:824`, `:677` |
| N5 | Low | integer range bitmask `next_power_of_two() - 1` overflows for `$PnR > 2^63` (debug panic; release wraps to a benign no-op mask) | `fcs.rs:364` |
| N6 | Low (cosmetic) | Phosphor-icon migration **incomplete** — ~13 buttons still use raw emoji/symbols (`💾`×5, `📁`,`📝`,`✖`,`↶`,`↷`,`⇊`,`⤢`,`⤓`); gate Save/Load (`💾`/`📁`) sit right beside session Save/Load that already use Phosphor | `gui.rs:1734-1735` + 11 more |
| N7 | Low (hardening) | CI actions pinned to mutable tags (`@v6`/`@v7`) not SHAs; Windows job skips the clippy `-D warnings` gate the macOS job runs | both workflows |

### Detail on the two worth acting on first

**N1 — the missed non-finite median (the most actionable).** `compensation.rs:271-282`:
the standalone `median()` (used by `channel_medians` → `compute_spillover`) sorts with
`partial_cmp().unwrap_or(Equal)` and returns `vals[n/2]` **without** the
`retain(is_finite)` that v0.1.7 added to the popstats and stats medians. A NaN raw
event in a control file survives the sort, becomes the channel median, bypasses the
`denom <= 0.0` guard (false for NaN), and yields a NaN matrix row — which the CLI
`compute-spillover -o` then writes to disk via `save_matrix_file` (no finiteness check
on **save**; `validate_square` only runs on load). The GUI path is safe (it rebuilds
via `from_parts`, which rejects non-finite). **Fix:** filter non-finite at the top of
`median()` (mirror `med_mean_cv`), and/or add a finiteness check to the matrix writer.

**N3 — Windows CI script injection.** `windows-installer.yml:58` is `$tag = "${{ github.ref_name }}"`
— GitHub substitutes the context into the PowerShell source *before* execution, and
PowerShell evaluates `$(...)` inside double-quoted strings, while `git check-ref-format`
permits `$()"` in tag names. The macOS workflow (`:50`, `tag="${GITHUB_REF_NAME}"`) is
the correct pattern. **Fix:** pass via `env:` and read `$env:TAG`. (Low — only someone
who can push tags can trigger it, but trivial to close.)

### Confirmed intact / clean (so they aren't re-litigated)
- **5/5 prior GUI fixes intact**, code re-read: `load_session` square-validates the
  override + reconciles `hist_sample_pop`/`hist_hidden`/`scatter_hidden`/`bool_refs`;
  `poll_batch` distinguishes `Done` vs disconnect (`got_done` checked first); Stats-tab
  `compensated.len()` guard present; quadrant gates half-open via `next_down()`;
  `shape_display_bbox` accounts for ellipse rotation.
- **5/5 v0.1.7 numeric fixes sound** on untested edges: gating unevaluable-tracking
  (transitive via the `unevaluable` set; cycle/self-ref terminate), compensation
  ill-conditioning gate (`max|inv|>1e6` rejects no plausible separable matrix — keep
  it), `parse_spillover` checked arithmetic + finiteness, fcs `$PAR×$TOT` bound (rejects
  no valid file), popstats/stats non-finite filtering (true no-op on finite data).
- **New v0.1.9 logic clean**: colormap session round-trip + legacy `viridis`-bool
  back-compat correct; `lerp_anchors` index math in-bounds (NaN saturates to 0, no
  panic); `turbo` clamped; font/style ordering correct; `themed_visuals` correct for
  dark **and** light (`linear_multiply` scales opacity, verified).
- **Dependencies clean**: `cargo audit` = **0 vulnerabilities** (10 informational
  unmaintained/unsound, all pre-existing — the Linux-only GTK cluster via `muda`, and
  `paste` via `nalgebra`; none in the new network surface). Pure-Rust `rustls` stack,
  **no `openssl`/`native-tls`**, no async runtime, no process-spawn beyond `open`'s
  platform opener (URL passed as a process arg, not via a shell). All deps from
  crates.io; no git/path/patched deps.

### Benign notes (no fix needed)
- The native "Check for Updates…" **menu** path lacks the double-click guard the
  toolbar button has — rapid menu clicks orphan harmless 15 s-bounded threads. One-line
  guard if desired.
- Still no `[profile.release]` → overflow-checks off in release. The one overflow that
  mattered (v0.1.6 M2) is now guarded by explicit `checked_mul`, so this is no longer a
  live crash vector; `overflow-checks = true` would be belt-and-suspenders.

### Suggested fix order
**N1** (real silent-bad-output through the CLI) → **N3** (CI injection, one-line) →
**N2** (open the constant `RELEASES_PAGE`, drop `html_url`) → **N4/N5** (cheap guards) →
**N6** (cosmetic emoji sweep) → **N7** (CI hardening). All Low; none blocks a release.

---

Full-codebase audit of the ~8.2k-LOC Rust tree at commit `5c8ad67`. Six parallel
review passes (FCS binary I/O · compensation+transforms · gating+stats · CLI ·
GUI state/cache · GUI gating geometry); **every finding below was independently
re-verified by reading the cited `file:line`.** Numeric parity vs R/flowCore
(parsing, compensation, asinh, logicle, per-population counts/MFI) was taken as
given and not re-derived — the focus is the edges validation can't reach:
non-finite propagation, malformed-input robustness, silent mis-counts, and GUI
state/cache correctness.

**Headline:** the engine is well-built and defensively coded; most historically
risky areas (channel-index clamping, cache invalidation, logicle convergence,
the integer bitmask) are correctly handled (see *Verified correct* below). The
findings cluster into three real themes: **(1)** one genuine correctness
regression — Boolean `NOT` of a missing-channel gate reports ~100% of its parent;
**(2)** a family of *silent non-finite / ill-conditioned* values that can poison
reported MFI with no error; **(3)** robustness gaps on malformed FCS/matrix input.

| Severity | Count |
|---|---|
| High | 1 |
| Medium | 5 |
| Low | 12 |
| Verified-correct (non-bugs) | 9 |

---

## HIGH

### H1 — Boolean `NOT` of a missing-channel gate reports the *parent's* count (silent ~100%)
- **Location:** `src/gating.rs:248-251` (the `Not` arm), root cause at `src/gating.rs:189-190`
- **Confidence:** High · reproducible by trace, no build needed
- **Mechanism:** A geometric gate on a channel the sample lacks correctly collapses
  to an all-false *own* mask (`gate_membership` errs → line 190
  `unwrap_or_else(|_| vec![false; n_events])`). But "channel absent" and
  "legitimately empty" are now indistinguishable. `boolean_mask`'s `Not`
  (line 249) inverts the referenced gate's **effective** mask:
  missing-channel gate `G` → `effective_mask(G)` all-false → `NOT` → **all-true**
  → `effective_mask(NOT-gate)` ANDs with its parent chain → count = **parent count**.
- **Why it matters:** This re-introduces the *exact* "reports the parent's count as
  a population" bug the project already fixed for geometric gates — now via the
  Boolean path. A natural `Live = NOT Dead` gate, evaluated on a tube missing the
  viability channel, silently reports ~100% live, and the error propagates to
  every child population and CSV row. Silent wrong counts in a scientific tool.
- **Reachability (verified):** Live via the **CLI** (`flowcyto popstats`/`gate`
  call `population_stats`/`apply_gates` directly — `main.rs:402`, `main.rs:608` —
  with no missing-channel guard) and via the **interactive Stats tab**
  (`gui.rs:3143`). **The Batch tab is protected**: `missing_gate_channels`
  (`gui.rs:4132`) skips any sample lacking a gate's channel before
  `population_stats` runs (`gui.rs:3329-3333`). `GateShape::Boolean` is
  serde-derived, so a hand-written/session gate triggers it headlessly.
- **Sibling issue → M4** below: `OR` has the same root cause, opposite direction.
- **Fix:** Distinguish *unevaluable* from *empty*. Track the gate ids whose channel
  was missing (the `Err` arm at 189-190) and force any Boolean that transitively
  references an unevaluable gate to all-false (or surface NaN/"n/a"). Cleaner:
  make own masks tri-state (`Option<Vec<bool>>`) so `NOT`/`OR` of an unevaluable
  ref stays unevaluable rather than inverting an all-false placeholder to all-true.
  **Key on *unevaluable*, not *empty*:** `NOT` of a legitimately empty gate (present
  channel, zero events) must still return all-true — the complement of an empty
  population really is everything — so a naive "NOT-of-empty → empty" fix is wrong.
- **Tested?** No. `missing_channel_gate_reads_zero_not_parent` (popstats.rs:234)
  covers a *geometric* missing-channel gate; `boolean_gates_and_or_not`
  (gating.rs:535) covers `NOT` on a *present* channel — the cross of the two
  (the bug) is untested.

---

## MEDIUM

### M1 — Silent non-finite / ill-conditioned compensation → poisoned MFI (a chained family)
Four related gaps let `NaN`/`Inf`/exploded values enter the event matrix and
silently become the reported median/mean/CV. No single one is high-likelihood on
a pristine FACSDiva file, but they share one cheap fix (finiteness/conditioning
guards) and the failure mode — *silent wrong numbers* — is the worst class here.

- **M1a — near-singular spillover inverts silently to huge values.**
  `src/compensation.rs:290-291`. `m.try_inverse()` is trusted to reject
  non-invertible matrices, but nalgebra returns `None` only at ~exact singularity;
  a determinant ~1e-13 (two nearly-collinear tandem dyes) inverts to entries
  ~1e13 with no condition-number/magnitude check, so `apply` (`:329`) yields
  compensated values ~1e16 or ±Inf. **Fix:** after inversion, reject if
  `!inv.iter().all(|x| x.is_finite())` or if `max|inv|` / reconstruction residual
  exceeds a sane bound. *(Confidence High that no guard exists; the 1e13 magnitude
  is a reviewer repro.)*
- **M1b — `parse_spillover` swallows bad tokens as 0.0 and admits `nan`/`inf`.**
  `src/compensation.rs:37` — `.unwrap_or(0.0)`. This is the **default** path
  (the sample's own embedded `$SPILLOVER`). A mangled token silently becomes 0.0
  (corrupting the matrix, pushing toward singular); `"nan"`/`"inf"` parse to
  `Ok(non-finite)` and pass straight through to the inverse. **Fix:** parse with
  context-bearing error + reject non-finite (mirror the CSV path).
- **M1c — CSV matrix import is only half-hardened.** `src/compensation.rs:85-91`
  + `validate_square` (`:121-135`). The CSV value parse correctly errors on
  `"abc"`, but `validate_square` checks **only dimensions**, so `nan`/`inf` cells
  import silently. (JSON path is safe — serde rejects `NaN`/`Infinity`.) **Fix:**
  add `if !rows.iter().flatten().all(|v| v.is_finite()) { bail! }` to
  `validate_square` (closes CSV *and* M5's session path).
- **M1d — no finiteness guard before median.** `src/popstats.rs:122-123`
  (also `compensation.rs:256`, `stats.rs:49`). `sort_by(partial_cmp().unwrap_or(Equal))`
  leaves a `NaN` in place (it compares `Equal` to everything), so `vals[n/2]` can
  *return* NaN, and with non-finite values present the sort isn't a valid total
  order (median becomes arbitrary). This is the amplifier: any upstream non-finite
  (M1a-c, or H?/L cofactor/log) silently becomes the reported MFI; counts/% are
  unaffected so it's easy to miss. **Fix:** `vals.retain(|v| v.is_finite())` before
  reducing (and/or guard at the compensation boundary).
- **Tested?** No test exercises a non-finite value *inside* a non-empty vector or a
  near-singular matrix; only exact-singular (`from_parts_singular_errors`) and
  empty-input NaN are covered.

### M2 — FCS parser: unbounded/overflowing allocation from `$PAR`×`$TOT` (crash on a malformed file)
- **Location:** `src/fcs.rs:286-287` (and the sibling `Vec::with_capacity(n_params)` at `:150`)
- **Confidence:** High
- **Mechanism:** `let total = n_params * n_events;` then `Vec::with_capacity(total)`,
  with `$PAR`/`$TOT` parsed from untrusted TEXT (`:142-148`) and **no bound against
  the (correctly bounded) DATA segment length**. `$PAR=2,$TOT=10^18` → `with_capacity`
  requests ~16 EB → process abort. Release builds ship with overflow-checks off
  (no `[profile.release]`), so an overflowing product instead wraps to a small
  capacity, the short read succeeds, and `FcsFile.n_events` retains the huge value
  → later `events[ev*n + ci]` indexing panics (see **L8**, which is the same wrap
  feeding the one unguarded Stats path).
- **Why it matters:** A single malformed/hostile tube crashes the app — and the GUI
  does directory-scale batch loading. No OOB write (memory-safe); this is DoS, not RCE.
- **Fix:** `n_params.checked_mul(n_events)` and reject if `total * bytes_per_value`
  can't fit `data_len` (already known) before allocating; same guard before `:150`.
- **Tested?** No.

### M3 — `export` silently overwrites — and destroys the input if it ends in `.csv`
- **Location:** `src/main.rs:565-571`
- **Confidence:** High
- **Mechanism:** With no `-o`, the output path is the input with extension set to
  `csv`; `csv::Writer::from_path` truncates it with no existence check. `set_extension("csv")`
  on a `.csv` input is a **no-op → output == input**, so running `export` on a valid
  FCS file named `data.csv` overwrites it in place (irreversible). Independently, any
  pre-existing `tube1.csv` beside `tube1.fcs` is clobbered without warning.
- **Fix:** `bail!` if `out_path == path`; optionally refuse an existing target
  unless `--force`.
- **Tested?** No (`export` has no CLI test).

### M4 — Boolean `OR` silently drops a missing-channel reference
- **Location:** `src/gating.rs:243-247`
- **Confidence:** High
- **Mechanism:** Same root cause as **H1**: a ref whose channel is absent contributes
  an all-false effective mask, so `A OR B` with `B` missing silently degrades to `A`
  — no error, no indication the OR ran over fewer inputs. Cannot inflate to the
  parent (so < H1), but under-counts silently. Same reachability (CLI + interactive;
  Batch skips). **Fix:** the unevaluable-tracking from H1.
- **Tested?** No.

### M5 — `parse_spillover` panics on a crafted/huge `n` (integer overflow → OOB slice)
- **Location:** `src/compensation.rs:25-26`, `:32`
- **Confidence:** High
- **Mechanism:** `let expected = 1 + n + n*n;` with `n` from the first `$SPILLOVER`
  token. For `n ≈ usize::MAX`, `n*n` wraps so `expected` collapses to ~1, the
  `parts.len() < expected` guard passes, and `parts[1..=n]` (`:32`) builds an
  out-of-bounds slice → panic. Reachable from a real file's `$SPILLOVER` via
  `from_keyword`, and from the GUI Spillover tab / CLI `spillover` (both call
  `parse_spillover` to display).
- **Fix:** checked arithmetic + bound `n` (e.g. `> parts.len()` or a channel cap) before slicing.
- **Tested?** No (overflow/large-`n` case uncovered).

---

## LOW

- **L1 — `transform-dump` mislabels compensated data as `raw`.** `src/main.rs:238,259-262`.
  With `--compensate`, `prepare_events` returns compensated values that are then
  written under the `raw` column header. In a harness whose stated purpose is
  flowCore validation, the label is misleading (values are correct *untransformed*,
  but post-compensation). Fix: keep the true raw value in the `raw` column, or rename
  it `input`.
- **L2 — CLI `--cofactor 0` (and `=-150`) → silent `inf`/`NaN`/sign-flip.**
  `src/main.rs:68-69,81-82,94-95,147-148` → `src/transform.rs:88,115`. The GUI clamps
  cofactor to `1.0..=100000.0` (`gui.rs:1493`); the CLI has no `value_parser`, so a
  typo yields non-finite output (feeds M1d) with exit 0. Fix: enforce `cofactor > 0`.
- **L3 — CSV/TSV float cells render `NaN`/`inf` literally.** `src/main.rs:582`,
  `src/popstats.rs:147-151,168-172`. An empty population (normal in deep trees)
  writes literal `NaN` for median/mean/cv; R `read.csv` then coerces the whole
  column to character. Fix: map non-finite → empty/`NA` via one helper.
- **L4 — `--compensate` with no matrix: CLI errors, GUI silently no-ops.**
  `src/main.rs:631-633` (bails) vs `src/gui.rs:4004` (`Ok(events.clone())`). Behavioral
  divergence; pick one policy or document it.
- **L5 — `Log { floor }` with `floor ≤ 0` → `-inf`/`NaN`.** `src/transform.rs:87`.
  `floor` is hardcoded `1.0` with no UI control, but the field is `Deserialize`, so a
  hand-edited session/gate JSON with `floor ≤ 0` produces non-finite display coords.
  Latent. Fix: clamp `floor` positive on `compile()`.
- **L6 — invalid logicle params silently fall back to Linear.** `src/transform.rs:61-64`.
  By design (prevents a panic) and tested, but a fluorescence axis silently switching
  biexponential→linear with no banner can be mistaken for real structure. Fix: surface
  the fallback in the UI.
- **L7 — supplemental TEXT uses a non-standard keyword (dead for real files) + unbounded alloc.**
  `src/fcs.rs:124-130`. The parser reads `$SUPTEXT_START`/`$SUPTEXT_END`, but the FCS
  3.x standard (and this project's own writer, `fcs_write.rs:35-36`) use
  `$BEGINSTEXT`/`$ENDSTEXT` — so any real file's supplemental TEXT is silently
  ignored, and the `vec![0u8; se-ss+1]` (`:130`) lacks the `file_len` bound the main
  TEXT path has (`:114`). Fix: read the standard keywords; bound `se` by `file_len`.
- **L8 — Stats tab omits the buffer-length guard its siblings have (latent panic).**
  `src/gui.rs:3141-3146` → `popstats.rs:93`. Every other consumer guards
  `compensated.len() >= n_events*n_params` (e.g. `gui.rs:562`); the Stats tab doesn't.
  Not reachable today (reprocess runs synchronously before render) but it's the exact
  buffer a future async load — or the M2 wrap — would underflow. Fix: add the one-line guard.
- **L9 — `load_session` installs `spill_override` from JSON with no shape check.**
  `src/gui.rs:2005` → indexed unchecked at `gui.rs:2426` (`rows[xi][yi]`). Every other
  override entry validates squareness (`load_matrix_file`→`validate_square`); a
  hand-edited/corrupt session with a ragged override panics when the Spillover/preview
  panel renders. Fix: validate (reuse `from_parts`) and drop with a status message.
- **L10 — `load_session` doesn't reconcile gate-id-keyed UI state.** `src/gui.rs:2006`.
  The reset block (1978-1981) and `self.gates = session.gates` leave `hist_sample_pop`
  (and `scatter_hidden`/`hist_hidden`/`bool_refs`) dangling against the new gate set.
  A stale `hist_sample_pop` makes `effective_mask(unknown_id)` return all-true
  (`gating.rs:269/274`), so the Samples-overlay histogram shows *all events* labeled as
  the missing population. `after_gate_restore` (`gui.rs:905`) already does the right
  reconciliation for undo/redo — reuse it here.
- **L11 — batch worker panic is reported as success.** `src/gui.rs:3341` (unwrapped
  `population_stats`) + `poll_batch:3373,3394`. A panic drops `tx`; `poll_batch` reads
  `Disconnected` identically to the normal `Done` and prints "Batch done: N processed".
  Memory-safe (no Mutex to poison) but a crashed sample silently vanishes from the
  exported CSV. Fix: only treat an explicit `BatchMsg::Done` as success; flag
  disconnect-without-Done.
- **L12 — minor geometry/robustness nits.**
  (a) **Quadrant gates double-count exact center lines** — `gui.rs:967-972` build four
  rects sharing `cx`/`cy`, and `Rect::contains` (`gating.rs:60-61`) is inclusive on
  both ends, so a point exactly on a center line lands in two quadrants. Measure-zero
  for continuous float data (why flowCore validation passed) but real for integer/
  floor-clamped ties; FlowJo assigns each boundary to one quadrant. (b) **`shape_display_bbox`
  ignores ellipse rotation** — `gui.rs:4075` drops `angle`, so "⊕ zoom to gate" crops a
  rotated ellipse's tips (membership/rendering unaffected). (c) **Exotic `$BYTEORD`
  permutations silently decode as big-endian** — `fcs.rs:190-196` treats anything but
  strictly-ascending as BE; `4,3,2,1` is correct, but a true permutation like `2,1,4,3`
  mis-decodes (rare instruments; no memory issue). (d) **`peek_events` lacks the
  `file_len` alloc bound** `open` has — `fcs.rs:76` vs `:114` (worst case a ~100 MB
  transient alloc; header offsets cap at 8 digits). (e) **`ui_samples` early-returns with
  one sample** — `gui.rs:3617` hides the Clear button and per-sample group-tag editor
  until a 2nd file is added (UX, not correctness).

---

## Verified correct (non-bugs — checked so they aren't re-investigated)

- **Integer range bitmask** `next_power_of_two()-1` (`fcs.rs:337-342`) — correct for
  non-power-of-2 ranges; the prior `range-1` corruption is genuinely gone. (Tested.)
- **"Fewer-channel sample → index panic" hypothesis — DISPROVEN.** `x_ch`/`y_ch`/`hist_ch`/
  grid/batch channel indices are clamped with `.min(n_params-1)` at activation *and* every
  use site (`gui.rs:504-505,564-565,997,2113-2114,3477,…`). Residual risk is the documented
  *wrong-channel* (silent, no crash), not a panic.
- **Logicle convergence** (`logicle.rs:99-171`) — Halley capped at 20 iters / bisection at
  200, returns best estimate rather than panicking; `scale(0)` and `w==0` short-circuit. No
  infinite loop, no `ln(0)`.
- **Boolean dependency ordering / cycles** (`gating.rs:198-219`) — terminates on cycle/
  self-ref/missing-ref via the `progressed` flag; unresolvable → all-false. (Tested.)
- **point_in_polygon** (`gating.rs:339-355`), **ellipse rotation + zero-radius guard**
  (`:63-72`), **inverted/empty rect/range** — all correct.
- **Percentage div-by-zero** is guarded everywhere (returns 0.0, not NaN) —
  `popstats.rs:83-84`, `gating.rs:288-293`.
- **GUI coordinate-transform backbone** — every remap is the correct
  `current.forward(gate.inverse(·))` / `gate.forward(current.inverse(·))` direction; gates
  store their own per-channel transform; `gate_clamp` is render-only (membership uses real
  bounds). Nonlinear (Log/Logicle) cases checked for space-mixing — none found.
- **Cache invalidation / `compensated`↔`fcs` frame-ordering** — `reprocess()` is the only
  event-data writer of `compensated` and runs synchronously before any render; caches carry
  per-frame staleness keys. The historically-buggy invariant holds. (L8/L10 are the two
  small gaps.)
- **FCS writer** (`fcs_write.rs`) — offset back-patch at fixed width with a length assert,
  delimiter chosen absent from all data, empty keys/values skipped, `$SPILLOVER`+`SPILL`
  both written. Round-trip clean.
- **CLI error handling** — no `unwrap`/`expect`/`panic`/raw-index in `main.rs`; every
  subcommand threads `anyhow` with context; `rewrite-spillover` validates channels +
  invertibility before writing.

---

## Test-coverage gaps (the seams the bugs live in)
Engine unit tests are strong, but none cover: Boolean × missing-channel (H1/M4);
non-finite values inside compensation/stats or near-singular inversion (M1);
oversized/overflowing `$PAR`/`$TOT` or `$SPILLOVER` `n` (M2/M5); and **every CLI
subcommand that writes files or parses gates/matrix JSON** (`export`, `gate`,
`popstats`, `compute-spillover`, `rewrite-spillover`, `transform-dump`) — i.e.
M3, L1, L2, L3 are all unguarded by tests.

## Suggested fix order
1. **H1 + M4** — one change (track unevaluable gate ids) fixes the only true
   correctness regression; add the Boolean×missing-channel regression test.
2. **M1** — one shared finiteness guard (`validate_square` + `parse_spillover` +
   a `retain(is_finite)` before median) closes M1b/c/d; add a conditioning check for M1a.
3. **M2** (+ L8) — `checked_mul` + `data_len` bound at the parser; add the Stats-tab guard.
4. **M3** — the `out_path == path` guard (cheap, prevents data loss).
