#!/usr/bin/env python3
"""Independent third-oracle cross-check of flowcyto's flowCore-validated numeric layers
against **cytoflow** (Python). Mirrors `validation/gen_golden.R`: it recomputes the SAME
probes on the SAME frozen `reference.fcs` and compares them to the SAME frozen flowCore
`golden.csv`. flowcyto<->golden is already CI-gated by `flowcyto selftest`, so cytoflow
reproducing the golden closes a three-way agreement: cytoflow <-> flowCore <-> flowcyto.

Run inside a cytoflow venv (cytoflow 1.2):
    python validation/gen_golden_cytoflow.py validation/reference.fcs validation/golden.csv

Environment: cytoflow has NO wheel for macOS arm64 / Python 3.11 - it is a from-source
build. See `validation/cytoflow-requirements.txt` for the exact verified pins and the
staged install recipe (build cytoflow with --no-build-isolation, then
`pip install --no-deps -r`). Result table is saved alongside in `cytoflow_crosscheck.txt`.

What is genuinely INDEPENDENT here (honest framing, per the paper's discipline):
  * parse   - cytoflow's FCS reader (fcsparser); a different codebase from flowCore's.
  * logicle - cytoflow's Wayne Moore C++ `Logicle` extension. NOTE: this is the SAME
              Moore & Parks reference flowCore derives from, so agreement is evidence of a
              faithful port, NOT of an independent algorithm.
  * gating  - cytoflow's geometric gating ops (Range/Range2D/Polygon membership).
  * asinh   - numpy.arcsinh (trivial; independent of flowcyto's Rust, but not cytoflow-specific).
  * comp    - numpy inv(M) @ raw. cytoflow has no apply-embedded-$SPILLOVER op (only
              estimate-from-controls), so this is independent linear algebra, not cytoflow.
"""
import csv, sys, math, warnings, collections
warnings.filterwarnings("ignore")
import numpy as np
import cytoflow as flow
from cytoflow.utility.logicle_ext.Logicle import Logicle

COFACTOR = 150.0
LG = dict(T=262144.0, W=0.5, M=4.5, A=0.0)            # explicit params, matching gen_golden.R
TOL = 1e-5
ref, golden_path = sys.argv[1], sys.argv[2]


def inject_pne(src):
    """cytoflow's reader requires $PnE (amplification); the synthetic reference omits it
    (flowCore tolerates that). Write a copy with $PnE=0,0 added and DATA bytes copied
    VERBATIM (for $DATATYPE F, cytoflow only reads $PnE - the log->linear rescale is gated
    on integer data - so values are untouched)."""
    b = open(src, "rb").read()
    ts, te = int(b[10:18]), int(b[18:26])
    ds, de = int(b[26:34]), int(b[34:42])
    text, data = b[ts:te + 1], b[ds:de + 1]
    delim = text[0:1]
    d = delim.decode("latin1")
    toks = text[1:].decode("latin1").split(d)
    if toks and toks[-1] == "":
        toks = toks[:-1]
    kv = list(zip(toks[0::2], toks[1::2]))
    NP = int(dict(kv)["$PAR"])
    have = {k for k, _ in kv}
    kv = [(k, "000000000000" if k in ("$BEGINDATA", "$ENDDATA") else v) for k, v in kv]
    for i in range(1, NP + 1):
        if f"$P{i}E" not in have:
            kv.append((f"$P{i}E", "0,0"))
    body = (d + "".join(f"{k}{d}{v}{d}" for k, v in kv)).encode("latin1")
    nts, nte = 58, 58 + len(body) - 1
    nds, nde = nte + 1, nte + len(data)
    body = body.replace(b"$BEGINDATA" + delim + b"000000000000",
                        b"$BEGINDATA" + delim + f"{nds:012d}".encode())
    body = body.replace(b"$ENDDATA" + delim + b"000000000000",
                        b"$ENDDATA" + delim + f"{nde:012d}".encode())
    hdr = bytearray(b" " * 58)
    hdr[0:6] = b"FCS3.0"
    for a, bb, v in [(10, 18, nts), (18, 26, nte), (26, 34, nds), (34, 42, nde), (42, 50, 0), (50, 58, 0)]:
        hdr[a:bb] = f"{v:>8}".encode()
    out = "/tmp/_cf_ref.fcs"
    open(out, "wb").write(bytes(hdr) + body + data)
    assert open(out, "rb").read()[nds:nde + 1] == data, "DATA bytes changed!"
    return out


def read_spillover(src):
    """Embedded $SPILLOVER -> (channel list, NxN matrix)."""
    b = open(src, "rb").read()
    ts, te = int(b[10:18]), int(b[18:26])
    text = b[ts:te + 1].decode("latin1")
    d = text[0]
    toks = text[1:].split(d)
    kv = dict(zip(toks[0::2], toks[1::2]))
    s = kv["$SPILLOVER"].split(",")
    n = int(s[0])
    chans = s[1:1 + n]
    M = np.array(list(map(float, s[1 + n:1 + n + n * n]))).reshape(n, n)
    return chans, M


# --- parse via cytoflow's ImportOp ---
cfref = inject_pne(ref)
ex = flow.ImportOp(tubes=[flow.Tube(file=cfref, conditions={})], conditions={}).apply()
raw = {ch: ex.data[ch].values.astype(np.float64) for ch in ex.channels}

# --- compensation: comp = obs @ inv(M)  (reference built as obs = true @ M; flowCore: x %*% solve(M)) ---
fl, M = read_spillover(ref)
obs = np.column_stack([raw[c] for c in fl])
comp_fl = obs @ np.linalg.inv(M)
comp = {c: raw[c] for c in raw}
for j, c in enumerate(fl):
    comp[c] = comp_fl[:, j]

# --- transforms ---
lg = Logicle(LG["T"], LG["W"], LG["M"], LG["A"])
logicle_f = lambda x: lg.scale(float(x)) * LG["M"]    # Moore normalized [0,1] -> flowCore [0,M]
asinh_f = lambda x: math.asinh(float(x) / COFACTOR)

# --- gating via cytoflow ops on the compensated experiment ---
for c in fl:                                          # inject compensated fluor; scatter stays raw
    ex.data[c] = comp[c]
e = flow.Range2DOp(name="Cells", xchannel="FSC-A", xlow=20000, xhigh=1e6,
                   ychannel="SSC-A", ylow=-1e6, yhigh=1e6).apply(ex)
e = flow.RangeOp(name="FITCpos", channel="FITC-A", low=5000, high=1e6).apply(e)
e = flow.PolygonOp(name="PEpos", xchannel="PE-A", ychannel="SSC-A",
                   vertices=[(5000, -1e6), (50000, -1e6), (50000, 1e6), (5000, 1e6)]).apply(e)
cells = e.data["Cells"].values.astype(bool)
fitc = cells & e.data["FITCpos"].values.astype(bool)  # FITC+ within Cells (hierarchy)
pep = e.data["PEpos"].values.astype(bool)
med = lambda m, ch: float(np.median(e.data[ch].values[m]))
gate = {
    ("", "Cells_count"): float(cells.sum()),
    ("", "FITCpos_count"): float(fitc.sum()),
    ("FITC-A", "FITCpos_median"): med(fitc, "FITC-A"),
    ("FSC-A", "FITCpos_median"): med(fitc, "FSC-A"),
    ("", "PEpos_count"): float(pep.sum()),
    ("PE-A", "PEpos_median"): med(pep, "PE-A"),
}

# --- compare every golden probe ---
LAYER = {"parse": "Parsing", "comp": "Compensation", "asinh": "Asinh",
         "logicle": "Logicle", "gate": "Gating"}
acc = collections.OrderedDict()                       # layer -> [probes, max_rel_dev, counts_exact]
mism = []
for r in csv.DictReader(open(golden_path)):
    k, ch, key, g = r["kind"], r["channel"], r["key"], float(r["golden"])
    if k == "asinh":
        c = asinh_f(key)
    elif k == "logicle":
        c = logicle_f(key)
    elif k == "parse":
        c = raw[ch][int(key)]
    elif k == "comp":
        c = comp[ch][int(key)]
    elif k == "gate":
        c = gate[(ch, key)]
    else:
        continue
    a = acc.setdefault(LAYER[k], [0, 0.0, True])
    a[0] += 1
    is_count = (k == "gate" and key.endswith("_count"))
    if is_count and int(round(c)) != int(round(g)):
        a[2] = False
        mism.append((LAYER[k], ch, key, c, g))
    rel = abs(c - g) / (abs(g) + 1.0)
    a[1] = max(a[1], rel)
    if rel > TOL and not is_count:
        mism.append((LAYER[k], ch, key, c, g))

# --- report ---
print(f"\ncytoflow {flow.__version__} cross-check vs frozen flowCore 2.24.0 golden")
print("(same reference.fcs, same probes/params as validation/gen_golden.R)\n")
print(f"{'Layer':<14}{'Probes':>7}{'Max rel dev':>15}{'Counts':>10}{'Result':>9}")
print("-" * 55)
allok = True
for L in ("Parsing", "Compensation", "Asinh", "Logicle", "Gating"):
    if L not in acc:
        continue
    p, mx, cok = acc[L]
    ok = (mx <= TOL) and cok
    allok &= ok
    counts = "exact" if cok else "MISMATCH"
    counts = counts if L == "Gating" else "-"
    print(f"{L:<14}{p:>7}{mx:>15.2e}{counts:>10}{('PASS' if ok else 'FAIL'):>9}")
print(f"\ngate counts: Cells={int(cells.sum())} FITC+={int(fitc.sum())} "
      f"PE+={int(pep.sum())}  (golden 566 / 200 / 160)")
if mism:
    print("\nMISMATCHES (first 20):")
    for L, ch, key, c, g in mism[:20]:
        print(f"  {L:<12} {ch or '-':<10} {key:<16} cytoflow={c:.8g} golden={g:.8g}")
print("\n" + ("ALL LAYERS REPRODUCE THE flowCore GOLDEN — three-way agreement "
              "cytoflow <-> flowCore <-> flowcyto."
              if allok else "FAILED: a layer deviates from the golden (see above)."))
sys.exit(0 if allok else 1)
