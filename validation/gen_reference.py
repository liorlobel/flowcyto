#!/usr/bin/env python3
"""Generate the fixed reference FCS used by `flowcyto selftest` (with a real
$SPILLOVER so compensation is non-trivially validated). Deterministic seed so the
frozen flowCore golden values stay reproducible. Run from the repo root:
    python3 validation/gen_reference.py validation/reference.fcs
"""
import struct, random, sys

random.seed(20260608)
N = 600
chans  = ["FSC-A", "SSC-A", "FITC-A", "PE-A", "PE-Cy7-A"]
labels = ["",      "",      "CD3",    "CD19", "CD11b"]
NP = len(chans)

def blob(n, c, s):
    return [[max(0.0, random.gauss(m, sd)) for m, sd in zip(c, s)] for _ in range(n)]

# A few populations spanning scatter + the three fluorescence channels (pre-spillover
# "true" signal; the embedded matrix then mixes them, which compensation must undo).
events = []
events += blob(200, [35000, 12000, 25000,   900,   900], [5000, 2500, 6000,  400,  400])  # CD3+
events += blob(160, [35000, 12000,   900, 24000,   900], [5000, 2500,  400, 6000,  400])  # CD19+
events += blob(160, [70000, 45000,   900,   900, 26000], [6000, 6000,  400,  400, 6000])  # CD11b+
events += blob( 80, [20000,  9000,   700,   700,   700], [4000, 3000,  300,  300,  300])  # debris
random.shuffle(events); events = events[:N]

# Embedded $SPILLOVER over the three fluorescence channels (observed = true x M).
sp = {("FITC-A","PE-A"):0.08, ("FITC-A","PE-Cy7-A"):0.01,
      ("PE-A","FITC-A"):0.03, ("PE-A","PE-Cy7-A"):0.15,
      ("PE-Cy7-A","FITC-A"):0.002, ("PE-Cy7-A","PE-A"):0.05}
fl = ["FITC-A", "PE-A", "PE-Cy7-A"]
M = [[1.0 if i==j else sp.get((fl[i], fl[j]), 0.0) for j in range(3)] for i in range(3)]
fidx = {c: chans.index(c) for c in fl}
for ev in events:
    true = [ev[fidx[c]] for c in fl]
    obs = [sum(true[i]*M[i][j] for i in range(3)) for j in range(3)]  # observed = true x M
    for j, c in enumerate(fl):
        ev[fidx[c]] = obs[j]
spill_kw = "3," + ",".join(fl) + "," + ",".join(f"{M[i][j]:g}" for i in range(3) for j in range(3))

data = bytearray()
for ev in events:
    for v in ev:
        data += struct.pack("<f", float(v))

delim = "/"
PAD = "000000000000"  # fixed-width $BEGINDATA/$ENDDATA placeholders (patched below)
kv = [("$BEGINANALYSIS","0"),("$ENDANALYSIS","0"),("$BEGINSTEXT","0"),("$ENDSTEXT","0"),
      ("$NEXTDATA","0"),("$BYTEORD","1,2,3,4"),("$DATATYPE","F"),("$MODE","L"),
      ("$PAR",str(NP)),("$TOT",str(len(events))),("$CYT","flowcyto-selftest"),
      ("$BEGINDATA",PAD),("$ENDDATA",PAD),("$SPILLOVER", spill_kw)]
for i,(c,l) in enumerate(zip(chans,labels),1):
    kv += [(f"$P{i}N",c),(f"$P{i}B","32"),(f"$P{i}R","262144")]
    if l: kv.append((f"$P{i}S",l))
text = (delim + "".join(f"{k}{delim}{v}{delim}" for k,v in kv)).encode()

ts, te = 58, 58 + len(text) - 1
ds, de = te + 1, te + len(data)
# Patch $BEGINDATA/$ENDDATA in place (same width → offsets stay valid).
text = text.replace(b"$BEGINDATA/" + PAD.encode(), b"$BEGINDATA/" + f"{ds:012d}".encode())
text = text.replace(b"$ENDDATA/" + PAD.encode(),   b"$ENDDATA/"   + f"{de:012d}".encode())
hdr = bytearray(b" "*58); hdr[0:6] = b"FCS3.0"
def put(a,b,v): hdr[a:b] = f"{v:>8}".encode()
put(10,18,ts); put(18,26,te); put(26,34,ds); put(34,42,de); put(42,50,0); put(50,58,0)

with open(sys.argv[1], "wb") as f:
    f.write(hdr); f.write(text); f.write(data)
print(f"wrote {sys.argv[1]}: {len(events)} events x {NP} params; spillover over {fl}")
