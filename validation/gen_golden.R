#!/usr/bin/env Rscript
# Generate the frozen flowCore golden values that `flowcyto selftest` checks against.
# Run from the repo root (needs R + Bioconductor flowCore):
#     Rscript validation/gen_golden.R validation/reference.fcs validation/golden.csv
# Re-run only when the reference FCS or the validated numerics change.
suppressMessages(library(flowCore))
args <- commandArgs(trailingOnly = TRUE)
ref <- args[1]; out <- args[2]

# Transforms use explicit params matching flowcyto's defaults (NOT estimateLogicle).
COFACTOR <- 150
LG <- list(w = 0.5, t = 262144, m = 4.5, a = 0)

# Grid spanning negatives (compensated data has them), zero, and the full range.
grid <- c(-2000, -500, -100, -10, -1, 0, 1, 10, 100, 500, 1000, 5000, 26214, 131072, 262144)
asinh_g   <- arcsinhTransform(a = 0, b = 1 / COFACTOR, c = 0)(grid)
logicle_g <- logicleTransform(w = LG$w, t = LG$t, m = LG$m, a = LG$a)(grid)

rows <- list()
add <- function(kind, channel, key, val) {
  rows[[length(rows) + 1]] <<- data.frame(
    kind = kind, channel = channel, key = as.character(key),
    golden = sprintf("%.10g", val), stringsAsFactors = FALSE)
}
for (i in seq_along(grid)) {
  add("asinh", "", grid[i], asinh_g[i])
  add("logicle", "", grid[i], logicle_g[i])
}

# Reference FCS: raw values (parse) + flowCore-compensated values.
ff  <- read.FCS(ref, transformation = FALSE, truncate_max_range = FALSE)
raw <- exprs(ff)
sp  <- spillover(ff)
mat <- sp[[which(!vapply(sp, is.null, logical(1)))[1]]]   # the embedded $SPILLOVER
comp <- exprs(compensate(ff, mat))
fl <- c("FITC-A", "PE-A", "PE-Cy7-A")
K  <- 40                                                    # first 40 events
for (ch in fl) for (e in seq_len(K)) {
  add("parse", ch, e - 1, raw[e, ch])
  add("comp",  ch, e - 1, comp[e, ch])
}

# ── Gating layer: population counts + median MFI on compensated data, so flowcyto's
#    Linear-transform gates operate on identical coordinates and match flowCore's
#    rectangleGate / polygonGate exactly. The discriminating boundaries are the FSC
#    scatter threshold (uncompensated → bit-identical between the two) and the
#    fluorescence neg/pos valley at 5000 (sparse → no boundary flips from the ~1e-10
#    compensation delta). Includes one hierarchy level (Cells -> FITC+).
comp_ff <- compensate(ff, mat)
cexpr   <- exprs(comp_ff)
lgl <- function(frame, g) as(filter(frame, g), "logical")   # membership vector
# G1 "Cells": FSC-A high (2-D rectangle; SSC spans the full range)
g1 <- rectangleGate(filterId = "Cells", "FSC-A" = c(20000, 1e6), "SSC-A" = c(-1e6, 1e6))
m1 <- lgl(comp_ff, g1)
add("gate", "", "Cells_count", sum(m1))
# G2: child of G1 — FITC+ (hierarchy / within-parent count + median MFI)
g2 <- rectangleGate(filterId = "FITCpos", "FITC-A" = c(5000, 1e6))
m2 <- lgl(Subset(comp_ff, g1), g2)
add("gate", "",       "FITCpos_count",  sum(m2))
# medians over the gated population (FITC+ within Cells)
mfi <- (m1 & lgl(comp_ff, g2))
add("gate", "FITC-A", "FITCpos_median", median(cexpr[mfi, "FITC-A"]))
add("gate", "FSC-A",  "FITCpos_median", median(cexpr[mfi, "FSC-A"]))
# G3 "PEpos": polygon on (PE-A, SSC-A)
poly <- polygonGate(filterId = "PEpos", .gate = matrix(
  c(5000, -1e6,  50000, -1e6,  50000, 1e6,  5000, 1e6),
  ncol = 2, byrow = TRUE, dimnames = list(NULL, c("PE-A", "SSC-A"))))
m3 <- lgl(comp_ff, poly)
add("gate", "", "PEpos_count", sum(m3))
add("gate", "PE-A", "PEpos_median", median(cexpr[m3, "PE-A"]))
cat("gate counts: Cells =", sum(m1), " FITC+ =", sum(m2), " PE+ =", sum(m3), "\n")

golden <- do.call(rbind, rows)
write.csv(golden, out, row.names = FALSE, quote = FALSE)
cat("wrote", nrow(golden), "golden values (flowCore", as.character(packageVersion("flowCore")), ") to", out, "\n")
