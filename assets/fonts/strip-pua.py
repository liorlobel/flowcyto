#!/usr/bin/env python3
"""Strip Private-Use-Area (PUA) glyphs from the bundled Inter TTFs.

Upstream Inter maps ~745 PUA codepoints (U+E000–U+F8FF). Because Inter is the
primary proportional face in the egui font stack (index 0), those glyphs shadow
the Phosphor icon font (index 1), whose icons live in the same PUA range —
making toolbar/button icons render as blank or wrong glyphs. Dropping Inter's
PUA glyphs lets the index-1 Phosphor face resolve every icon, while all real
text glyphs (Latin, punctuation, etc.) are untouched.

Run from the repo root after updating Inter:
    python3 -m venv /tmp/ftenv && /tmp/ftenv/bin/pip install fonttools
    /tmp/ftenv/bin/python assets/fonts/strip-pua.py

Requires fonttools. Operates in place on the two Inter TTFs in this directory.
"""
import os
from fontTools.ttLib import TTFont
from fontTools.subset import Subsetter, Options

HERE = os.path.dirname(os.path.abspath(__file__))
TARGETS = ["Inter-Regular.ttf", "Inter-SemiBold.ttf"]
PUA_LO, PUA_HI = 0xE000, 0xF8FF


def strip_pua(path):
    font = TTFont(path)
    cmap = font.getBestCmap()
    keep = sorted(c for c in cmap if not (PUA_LO <= c <= PUA_HI))
    opts = Options()
    opts.layout_features = ["*"]   # keep all OpenType layout features for kept glyphs
    opts.glyph_names = True
    opts.name_IDs = ["*"]          # keep the full name table (family name, etc.)
    opts.notdef_outline = True
    opts.recalc_bounds = True
    opts.drop_tables = []
    ss = Subsetter(options=opts)
    ss.populate(unicodes=keep)
    ss.subset(font)
    font.save(path)
    after = TTFont(path).getBestCmap()
    pua = [c for c in after if PUA_LO <= c <= PUA_HI]
    assert not pua, f"{path}: {len(pua)} PUA glyphs remain"
    assert 0x61 in after and 0x41 in after and 0x20 in after, f"{path}: lost text glyphs"
    print(f"{os.path.basename(path)}: PUA stripped, text glyphs intact")


if __name__ == "__main__":
    for t in TARGETS:
        strip_pua(os.path.join(HERE, t))
