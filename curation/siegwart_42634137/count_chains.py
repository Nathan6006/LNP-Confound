#!/usr/bin/env python3
"""Count zig-zag carbons in a raster structure drawing, by apex detection.

Fig 1b's "commercially available cores" are drawn, not tabulated, and the only
thing that distinguishes e.g. a C12 from a C14 tail is how many apexes the
zig-zag has. Counting those by eye across TT3's six chains, OF-2's four and
Lipid 14's two is exactly the kind of thing that silently goes wrong, so it is
done here instead.

Method: inside a given box, take the ink mask, and for each x column compute the
mean y of ink. That traces the chain's centre line. Apexes are the local minima
and maxima of that trace; a bond vertex sits at every apex, so the number of
carbons in a run is (number of apexes) plus the endpoints that are not apexes.

The figure is a 600 dpi render of a 217 dpi raster, so the bonds are smooth
enough for the trace to be unambiguous. Each box is reported with its apex
x-positions so the spacing can be checked for regularity - an irregular spacing
means the box caught something other than a clean chain.
"""
import sys

import numpy as np
import PIL.Image as Image
import PIL.ImageDraw as ImageDraw

PAGE = "00_pages/p2_600.png"


def trace(img, box, ink=150):
    x0, y0, x1, y1 = box
    sub = np.array(img.convert("L")).astype(int)[y0:y1, x0:x1]
    mask = sub < ink
    xs, ys = [], []
    for c in range(mask.shape[1]):
        col = np.nonzero(mask[:, c])[0]
        if len(col):
            xs.append(c + x0)
            ys.append(col.mean() + y0)
    return np.array(xs), np.array(ys)


def apexes(xs, ys, smooth=5, min_sep=12, min_prom=3.0):
    if len(ys) < smooth * 2:
        return []
    k = np.ones(smooth) / smooth
    sm = np.convolve(ys, k, mode="same")
    out = []
    for i in range(smooth, len(sm) - smooth):
        w = sm[i - smooth:i + smooth + 1]
        if sm[i] == w.min() and (w.max() - sm[i]) >= min_prom:
            out.append((xs[i], "min"))
        elif sm[i] == w.max() and (sm[i] - w.min()) >= min_prom:
            out.append((xs[i], "max"))
    dedup = []
    for x, kind in out:
        if dedup and x - dedup[-1][0] < min_sep:
            continue
        dedup.append((x, kind))
    return dedup


def report(name, box, img, draw=None):
    xs, ys = trace(img, box)
    ap = apexes(xs, ys)
    pos = [a[0] for a in ap]
    gaps = np.diff(pos) if len(pos) > 1 else np.array([])
    print(f"{name:26} apexes {len(ap):3}  span {box[0]}-{box[2]}  "
          f"gap med {np.median(gaps) if len(gaps) else float('nan'):.1f} "
          f"cv {gaps.std()/gaps.mean() if len(gaps) else float('nan'):.3f}")
    if draw is not None:
        draw.rectangle(box, outline=(255, 0, 0), width=3)
        for x, kind in ap:
            i = int(np.searchsorted(xs, x))
            i = min(max(i, 0), len(ys) - 1)
            draw.ellipse([x - 6, ys[i] - 6, x + 6, ys[i] + 6],
                         outline=(0, 160, 0) if kind == "max" else (0, 0, 255), width=3)
    return len(ap), pos


BOXES = {
    # TT3: three identical arms, each -NH-(CH2)3-N(CnH2n+1)2
    "TT3 top-right upper": (4030, 2630, 4700, 2730),
    "TT3 top-right lower": (4080, 2700, 4720, 2790),
    "TT3 top-left upper": (1000, 2620, 1540, 2720),
    "TT3 top-left lower": (960, 2690, 1500, 2780),
    "TT3 bottom left": (1200, 3100, 1760, 3190),
    "TT3 bottom right": (1800, 3100, 2380, 3190),
    # Lipid 14
    "Lipid14 upper tail": (3700, 2880, 4320, 2960),
    "Lipid14 lower tail": (3700, 3080, 4330, 3170),
}

if __name__ == "__main__":
    img = Image.open(PAGE)
    ov = img.convert("RGB").crop((900, 2550, 4800, 3250))
    d = ImageDraw.Draw(ov)
    off = (900, 2550)
    for name, box in BOXES.items():
        n, pos = report(name, box, img)
        shifted = (box[0] - off[0], box[1] - off[1], box[2] - off[0], box[3] - off[1])
        xs, ys = trace(img, box)
        ap = apexes(xs, ys)
        d.rectangle(shifted, outline=(255, 0, 0), width=3)
        for x, kind in ap:
            i = min(max(int(np.searchsorted(xs, x)), 0), len(ys) - 1)
            cx, cy = x - off[0], ys[i] - off[1]
            d.ellipse([cx - 7, cy - 7, cx + 7, cy + 7],
                      outline=(0, 160, 0) if kind == "max" else (0, 0, 255), width=3)
    ov.save("review/fig1b_chain_counts.png")
    print("\noverlay -> review/fig1b_chain_counts.png")
