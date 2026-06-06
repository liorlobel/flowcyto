# flowcyto

A fast, native **macOS app (and CLI)** for analyzing BD flow-cytometry `.fcs`
files — compensation, Linear/Log/Asinh/Logicle transforms, hierarchical gating,
per-population statistics, and multi-sample batch export. Every numeric layer is
cross-validated against R/[flowCore](https://bioconductor.org/packages/flowCore/).

## Install

Download **`flowcyto-<version>.dmg`** from the
[Releases page](https://github.com/liorlobel/flowcyto/releases), drag **flowcyto**
into **Applications**, then clear Gatekeeper once on first launch (the app is
un-notarized):

```bash
xattr -dr com.apple.quarantine /Applications/flowcyto.app
```

(or right-click → Open on macOS ≤14, or "Open Anyway" in System Settings →
Privacy & Security on macOS 15+). Full details: [INSTALL.md](INSTALL.md).

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
   the **2×2 grid** layout (top of the Plot tab) to see several gates at once.

5. **Read the numbers.** The **Stats** tab shows a per-population table
   (count, %parent, %total, median MFI per channel) with **💾 Export CSV (tidy)**.
   With multiple files open, the **Batch** tab runs your whole gate tree over every
   sample and exports one combined tidy CSV — tag each sample with a **group**
   (condition) in the Samples list and it becomes a column in the output.

**Save your work:** **💾 Save** (in Gates) writes the gate tree to JSON; the
toolbar's **🖫 Save session** stores everything — samples, gates, transforms, and
compensation — to reopen later.

### Handy extras
- **Histogram** tab: 1-D overlays of populations or samples; drag to add an
  interval gate.
- **Spillover** tab: view/edit the compensation matrix, or compute one from
  single-stain controls; the inline **⚖ Compensation** panel on the Plot tab
  adjusts the current X↔Y spillover with a live preview.
- **Keyboard:** `R`/`E`/`P`/`Q` draw tools, `G` edit, `V`/`Esc` navigate,
  `Ctrl/Cmd+Z` undo, `1`–`5` switch tabs.
- **📷 Save plot…** exports the current plot as a PNG.

## Command line

The same engine is scriptable:

```bash
flowcyto info    sample.fcs                         # metadata + channel list
flowcyto stats   sample.fcs --compensate            # per-channel summary
flowcyto popstats sample.fcs -g gates.json -o csv   # per-population stats → tidy CSV
flowcyto gate    sample.fcs -g gates.json           # gate counts
flowcyto gui     sample.fcs                          # open in the GUI
```

## Building from source

Requires the Rust toolchain.

```bash
cargo build --release          # binary at target/release/flowcyto
cargo test --release           # 88 tests
./packaging/make-macos-app.sh  # build the .app + .dmg
```
