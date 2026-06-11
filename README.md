<p align="center">
  <img src="packaging/icon.png" alt="flowcyto logo" width="160">
</p>

<h1 align="center">flowcyto</h1>

A fast, native **macOS, Windows &amp; Linux app (and CLI)** for analyzing BD
flow-cytometry `.fcs` files — compensation, Linear/Log/Asinh/Logicle transforms,
hierarchical gating, per-population statistics, antibody-titration stain
indices, and multi-sample batch export (CSV/XLSX).
Every numeric layer is cross-validated against
R/[flowCore](https://bioconductor.org/packages/flowCore/).

## Install

Grab the latest build from the
[Releases page](https://github.com/liorlobel/flowcyto/releases).

**macOS (Apple Silicon)** — download **`flowcyto-<version>.dmg`**, drag **flowcyto**
into **Applications**, then clear Gatekeeper once on first launch (the app is
un-notarized):

```bash
xattr -dr com.apple.quarantine /Applications/flowcyto.app
```

(or right-click → Open on macOS ≤14, or "Open Anyway" in System Settings →
Privacy & Security on macOS 15+). Full details: [INSTALL.md](INSTALL.md).

**Windows (x64)** — download **`flowcyto-<version>-setup.exe`** and run it
(Start-Menu shortcut + uninstaller). It's unsigned, so SmartScreen may warn:
click **More info → Run anyway**.

**Linux (x64)** — two options:
- **AppImage** (**`flowcyto-<version>-x86_64.AppImage`**) — `chmod +x` it and run; no
  install needed. It expects FUSE; if that's unavailable, run it with
  `./flowcyto-<version>-x86_64.AppImage --appimage-extract-and-run`.
- **Debian/Ubuntu `.deb`** (**`flowcyto_<version>-1_amd64.deb`**) —
  `sudo apt install ./flowcyto_<version>-1_amd64.deb`, which pulls in the few
  GL/xkb libraries it needs and adds a desktop entry. Built on Ubuntu 22.04 for
  broad glibc compatibility.

## Tutorial: gate a sample in 5 steps

1. **Open data.** Launch flowcyto and click **📂 Open FCS** (select one or more
   files). The first sample loads automatically. If the file carries an embedded
   `$SPILLOVER` matrix, **Compensate** turns on by default.

2. **Choose axes & scale.** In the left panel under **Channels**, set **X** and
   **Y**, then pick each axis's **scale** — keep `Linear` for FSC/SSC, use
   `Logicle` (or `Asinh`) for fluorescence. The **⇊ X scale → all fluorescence**
   button applies the X scale to every fluorescence channel at once.

3. **Draw a gate.** Under **Gates**, pick a tool — **▭ Rect**, **⬭ Ellipse**,
   **⬠ Polygon**, or **✛ Quad** — then drag (or click vertices for a polygon) on
   the plot. The gate appears in the tree with its **% of parent** and count. Use
   **✎ Edit** to drag handles to resize, drag the body to move, or rotate an
   ellipse. `Ctrl/Cmd+Z` undoes.

4. **Drill down.** **Double-click inside a gate** to "gate from here" — the plot
   restricts to that population so you can draw child gates on it (e.g. singlets →
   live → CD45⁺ → lineage). The **Viewing:** breadcrumb shows your path; click any
   level to jump back. Switch the left X/Y to gate on different channels, or use
   the **Grid** layout (top of the Plot tab) — an adjustable **cols × rows** grid
   from 1×1 up to **6×6** where each cell has its own X/Y and scale, so you can
   see (and gate on) the whole gating sequence at once.

5. **Read the numbers.** The **Stats** tab shows a per-population table
   (count, %parent, %total, median MFI per channel) — export it as a tidy **CSV**
   or a formatted **XLSX**. With multiple files open, the **Batch** tab runs your
   whole gate tree over every sample and exports one combined table (CSV or XLSX) —
   tag each sample with a **group** (condition) in the Samples list and it becomes a
   column in the output.

**Save your work:** **💾 Save** (in Gates) writes the gate tree to JSON; the
toolbar's **🖫 Save session** stores everything — samples, gates, transforms, and
compensation — to reopen later.

### Handy extras
- **Reproducibility report** (toolbar → **Report**): one self-contained `.html`
  capturing the full provenance of an analysis — file + instrument, the
  compensation matrix and its source, display transforms, the gating strategy,
  per-population statistics, and the gate JSON for re-import. Drop it into a
  paper's supplement so a result is reproducible.
- **R bridge** (Batch tab → **R bundle**): exports the tidy per-population CSV
  **plus** a ready-to-run starter `.R` (base R, no packages) wired to your
  group/condition tags — reads the CSV, reshapes to population frequencies,
  summarizes and plots by group, and scaffolds a comparison with the right
  flow-cytometry caveats (the tube is the unit of replication, block by
  experiment, mind rare populations and multiple comparisons).
- **Auto-gate suggestions** (in the Gates panel — starting points you review):
  **Suggest singlets** fits a diagonal FSC-A × FSC-H band, and **Auto-threshold X**
  splits the current channel at its density valley (handles a small/rare positive
  population, where a simple Otsu cut would not).
- **QC** tab: a per-tube acquisition-quality scan across the whole workspace —
  event count, **%viable** (pick any gate as the live/viable population; its % of
  parent comes from the same validated gating engine), **flow-rate stability**
  (a clog/bubble shows as a dip in the per-tube event-rate sparkline), and
  **off-scale/saturation** events. Flagged tubes are highlighted; export the table
  as CSV.
- **Histogram** tab: 1-D overlays of populations or samples; drag to add an
  interval gate.
- **Titration** tab: for an antibody dilution series, the per-sample **stain index**
  — `(MFI⁺ − MFI⁻) / (2·rSD⁻)`, the robust-SD definition FlowJo uses — on the chosen
  channel, with the optimal (highest-index) concentration highlighted and a
  stain-index-vs-dose curve. Tag each sample with its concentration, pick the
  positive and negative populations, and export CSV/XLSX. (Runs over the Batch
  results, so it's a view across your whole dilution series.)
- **Spillover** tab: view/edit the compensation matrix, or compute one from
  single-stain controls; the inline **⚖ Compensation** panel on the Plot tab
  adjusts the current X↔Y spillover with a live preview; and the **N×N cross-check**
  lays out every fluorochrome pair on compensated data with a residual-spillover
  flag, so over-/under-compensation is visible at a glance (the flag is a
  spillover-coefficient estimate — quantitative on single-stain controls, a visual
  check on mixed samples).
- **Menu bar (macOS):** native **File / Edit / View** menus mirror the in-app
  controls — Open FCS (⌘O), Save Gates (⌘S) / Session (⇧⌘S), Undo/Redo (⌘Z / ⇧⌘Z),
  switch tabs (⌘1–⌘6), toggle light/dark.
- **Keyboard:** `R`/`E`/`P`/`Q` draw tools, `G` edit, `V`/`Esc` navigate,
  `Ctrl/Cmd+Z` undo, `1`–`6` switch tabs.
- **Appearance:** light/dark themes and five density colormaps — **Viridis**
  (default), **Magma**, **Turbo**, **Cividis** (all perceptually-uniform and
  colorblind-safe), and **Jet** (legacy). The chosen colormap is saved with the
  session.
- **Save plot** exports the current plot as a PNG.
- **Check for updates:** the toolbar's **⟳ Updates** button (or **flowcyto →
  Check for Updates…** on macOS) compares your version against the latest GitHub
  release and opens the download page if there's a newer one. It's the app's only
  network access and runs *only* when you click it — never on launch.

## Command line

The same engine is scriptable:

```bash
flowcyto info    sample.fcs                         # metadata + channel list
flowcyto stats   sample.fcs --compensate            # per-channel summary
flowcyto popstats sample.fcs -g gates.json -o csv   # per-population stats → tidy CSV
flowcyto gate    sample.fcs -g gates.json           # gate counts
flowcyto gui     sample.fcs                          # open in the GUI
flowcyto update                                      # check GitHub for a newer release
flowcyto selftest                                    # verify the numerics vs flowCore
```

### Validated against flowCore — and you can check it yourself

`flowcyto selftest` recomputes parsing, compensation, asinh, logicle, and gating
(population counts + median MFI) on a bundled
reference and compares to **frozen R/[flowCore](https://bioconductor.org/packages/flowCore/)
golden values** — offline, no R needed. It prints a benchmark table and exits non-zero
on any deviation (it's also part of the test suite, so CI guards it on every change):

```
Layer                 Probes    Max rel. dev   Tolerance   Result
Asinh transform           15        4.73e-11        1e-5     PASS
Compensation             120        4.47e-10        1e-5     PASS
Gating                     6        1.90e-10        1e-5     PASS
Logicle transform         15        1.07e-10        1e-5     PASS
Parsing                  120        3.82e-10        1e-5     PASS
```

(Regenerate the golden with `validation/gen_golden.R` if you have R + flowCore.)

## Building from source

Requires the Rust toolchain.

```bash
cargo build --release          # binary at target/release/flowcyto
cargo test --release           # 134 tests
./packaging/make-macos-app.sh  # macOS: build the .app + .dmg
./packaging/make-appimage.sh   # Linux: build the AppImage (.deb via `cargo deb`)
```

The installers are also built by CI: pushing a `vX.Y.Z` tag runs the
[macOS](.github/workflows/macos-installer.yml),
[Windows](.github/workflows/windows-installer.yml), and
[Linux](.github/workflows/linux-installer.yml) workflows, which build the
`.dmg` (on a Mac runner), the Inno Setup `.exe` (on a Windows runner), and the
`.deb` + AppImage (on Ubuntu) and attach them all to the GitHub release.
