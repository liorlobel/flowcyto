#!/usr/bin/env python3
"""Generate the flowcyto app icon (packaging/icon.png).

The mark is a flow-cytometry density plot: filled viridis contour bands over a
dark navy panel, L-shaped axes, and a light gate ring around the bright
population — a nod to the app's core (density plots + gating). Fully procedural
so the icon can be regenerated and tweaked.

    python3 packaging/make-icon.py            # -> packaging/icon.png (1024px)
"""
from __future__ import annotations
import os
import numpy as np
from PIL import Image, ImageDraw, ImageFilter

S = 1024  # canvas size (px)

# ── viridis colormap (matplotlib anchor points) ────────────────────────────
_VIRIDIS = np.array([
    [0.267004, 0.004874, 0.329415], [0.282623, 0.140926, 0.457517],
    [0.253935, 0.265254, 0.529983], [0.206756, 0.371758, 0.553117],
    [0.163625, 0.471133, 0.558148], [0.127568, 0.566949, 0.550556],
    [0.134692, 0.658636, 0.517649], [0.266941, 0.748751, 0.440573],
    [0.477504, 0.821444, 0.318195], [0.741388, 0.873449, 0.149561],
    [0.993248, 0.906157, 0.143936],
])


def viridis(t: np.ndarray) -> np.ndarray:
    """Map t in [0,1] to RGB (float 0..1) via piecewise-linear viridis."""
    t = np.clip(t, 0.0, 1.0)
    x = t * (len(_VIRIDIS) - 1)
    i = np.floor(x).astype(int)
    i = np.clip(i, 0, len(_VIRIDIS) - 2)
    f = (x - i)[..., None]
    return _VIRIDIS[i] * (1 - f) + _VIRIDIS[i + 1] * f


def gaussian_blob(gx, gy, cx, cy, sx, sy, rot, amp):
    """Add a rotated 2-D gaussian to the density field."""
    ct, st = np.cos(rot), np.sin(rot)
    dx, dy = gx - cx, gy - cy
    u = dx * ct + dy * st
    v = -dx * st + dy * ct
    return amp * np.exp(-0.5 * ((u / sx) ** 2 + (v / sy) ** 2))


def rounded_mask(size, radius, ss=4):
    """Antialiased rounded-square alpha mask (0..1)."""
    m = Image.new("L", (size * ss, size * ss), 0)
    d = ImageDraw.Draw(m)
    d.rounded_rectangle([0, 0, size * ss - 1, size * ss - 1],
                        radius=radius * ss, fill=255)
    m = m.resize((size, size), Image.LANCZOS)
    return np.asarray(m, dtype=np.float32) / 255.0


def main():
    # coordinate grids
    yy, xx = np.mgrid[0:S, 0:S]
    gx = xx / S
    gy = 1.0 - yy / S  # math-y: bottom-left origin

    # ── background: vertical navy→teal-navy gradient ───────────────────────
    top = np.array([0.043, 0.094, 0.165])     # deep navy  #0b182a
    bot = np.array([0.055, 0.157, 0.204])     # teal navy  #0e2834
    g = (yy / S)[..., None]
    bg = top * (1 - g) + bot * g

    # ── density field: a main correlated population + a distinct island ────
    # The island (upper-right) is the gated population — kept separated from
    # the main cloud so the gate ring clearly "captures" it.
    main = np.zeros((S, S), dtype=np.float32)
    main += gaussian_blob(gx, gy, 0.50, 0.43, 0.165, 0.075, np.deg2rad(-27), 1.00)
    main += gaussian_blob(gx, gy, 0.37, 0.33, 0.090, 0.052, np.deg2rad(-27), 0.50)

    island = gaussian_blob(gx, gy, 0.71, 0.66, 0.072, 0.060, np.deg2rad(-18), 0.95)

    dens = np.maximum(main, island)
    dens /= dens.max()

    # crisp filled-contour banding (quantize), but keep a touch of smoothing
    LEVELS = 8
    banded = np.ceil(dens * LEVELS) / LEVELS
    color = viridis(banded)  # (S,S,3)

    # where the plot is "ink" vs background: smooth threshold for AA edges
    thr = 0.06
    soft = np.clip((dens - thr) / 0.05, 0.0, 1.0)
    soft = soft[..., None]

    img = bg * (1 - soft) + color * soft

    # ── axes: light L-shape (origin bottom-left), rounded caps ─────────────
    rgb = (np.clip(img, 0, 1) * 255).astype(np.uint8)
    im = Image.fromarray(rgb, "RGB")
    dr = ImageDraw.Draw(im, "RGBA")
    axis = (181, 199, 219, 255)   # soft slate
    w = 14
    m0, m1 = 0.22 * S, 0.80 * S   # axis extents in px
    # y axis (vertical) and x axis (horizontal)
    dr.line([(m0, 0.20 * S), (m0, m1)], fill=axis, width=w)
    dr.line([(m0, m1), (0.80 * S, m1)], fill=axis, width=w)
    r = w / 2
    for (cx, cy) in [(m0, 0.20 * S), (0.80 * S, m1)]:  # rounded end caps
        dr.ellipse([cx - r, cy - r, cx + r, cy + r], fill=axis)
    dr.ellipse([m0 - r, m1 - r, m0 + r, m1 + r], fill=axis)

    # ── gate ring around the bright satellite population ───────────────────
    # drawn on its own layer so we can soft-glow + blend
    gate = Image.new("RGBA", (S, S), (0, 0, 0, 0))
    gd = ImageDraw.Draw(gate)
    gcx, gcy = 0.71 * S, (1 - 0.66) * S
    gw, gh = 0.125 * S, 0.110 * S
    gd.ellipse([gcx - gw, gcy - gh, gcx + gw, gcy + gh],
               outline=(238, 243, 251, 240), width=10)
    im.paste(Image.alpha_composite(im.convert("RGBA"), gate).convert("RGB"))

    # ── apply rounded-square mask → transparent corners ────────────────────
    mask = rounded_mask(S, radius=int(0.225 * S))
    out = im.convert("RGBA")
    arr = np.asarray(out).astype(np.float32)
    arr[..., 3] = mask * 255.0
    out = Image.fromarray(np.clip(arr, 0, 255).astype(np.uint8), "RGBA")

    # ── macOS app-icon grid: inset the squircle to 824/1024 (~80.5%) with the
    #    standard transparent margin. The Dock renders an app icon filling its
    #    image bounds; a full-bleed icon therefore looks oversized (and a running
    #    app's tile renders slightly larger than its pinned/idle tile, which made
    #    the size jump visible). The margin makes flowcyto match other apps.
    CONTENT = 824
    off = (S - CONTENT) // 2
    canvas = Image.new("RGBA", (S, S), (0, 0, 0, 0))
    canvas.paste(out.resize((CONTENT, CONTENT), Image.LANCZOS), (off, off))
    out = canvas

    here = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(here, "icon.png")
    out.save(path)
    print(f"wrote {path}  ({out.size[0]}x{out.size[1]})")


if __name__ == "__main__":
    main()
