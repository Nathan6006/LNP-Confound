#!/usr/bin/env python3
"""Count carbons in a drawn zigzag alkyl chain by finding its apexes, not by eye.

Reading chain length off a structure diagram visually is the single failure mode most
likely to corrupt a lipid library (C12 vs C14 is one bond and changes every downstream
descriptor). This measures it instead: isolate the red structure pixels, take the mean
y of the stroke for each x column, and count turning points in that waveform. Each
turning point is a vertex (a carbon); the chain terminus adds one more.

Calibrated against 9A and 10A, whose identities are pinned independently by the HRMS in
Supplementary Figures 5-8 (di-n-octylamine and di-n-decylamine).

Usage: python3 -I count_vertices.py <page.png> <y0> <y1> <x0> <x1> [label]
"""
import sys
import numpy as np
from PIL import Image


def red_mask(img):
    a = np.asarray(img.convert("RGB")).astype(int)
    r, g, b = a[..., 0], a[..., 1], a[..., 2]
    return (r > 110) & (r - g > 55) & (r - b > 55)


def waveform(mask):
    """Mean y of stroke per x column; NaN where the column is empty."""
    ys = []
    for x in range(mask.shape[1]):
        col = np.flatnonzero(mask[:, x])
        ys.append(col.mean() if col.size else np.nan)
    return np.array(ys)


def turning_points(w, smooth=3, min_sep=6, min_amp=2.0):
    """Count local maxima and minima of the zigzag waveform."""
    v = w.copy()
    idx = np.flatnonzero(~np.isnan(v))
    if idx.size < 10:
        return 0, []
    v = np.interp(np.arange(len(v)), idx, v[idx])
    k = np.ones(smooth) / smooth
    v = np.convolve(v, k, mode="same")
    lo, hi = idx.min() + smooth, idx.max() - smooth
    pts = []
    for i in range(lo + 1, hi):
        seg = v[max(lo, i - min_sep):min(hi, i + min_sep) + 1]
        if v[i] == seg.max() or v[i] == seg.min():
            if not pts or i - pts[-1][0] >= min_sep:
                amp = seg.max() - seg.min()
                if amp >= min_amp:
                    pts.append((i, float(v[i])))
    return len(pts), pts


def main():
    img = Image.open(sys.argv[1])
    y0, y1, x0, x1 = (int(v) for v in sys.argv[2:6])
    label = sys.argv[6] if len(sys.argv) > 6 else ""
    sub = img.crop((x0, y0, x1, y1))
    m = red_mask(sub)
    rows = np.flatnonzero(m.any(axis=1))
    if rows.size == 0:
        print(f"{label}: no structure pixels found")
        return
    # split the band into the two chains (upper and lower) around the biggest vertical gap
    gaps = np.diff(rows)
    if gaps.size and gaps.max() > 4:
        split = rows[int(np.argmax(gaps))] + int(gaps.max()) // 2
        bands = [("upper", m[:split]), ("lower", m[split:])]
    else:
        bands = [("single", m)]
    out = []
    for name, band in bands:
        n, _p = turning_points(waveform(band))
        out.append((name, n))
    print(f"{label}: " + "  ".join(f"{n_}={c} turns -> ~{c + 1} C" for n_, c in out))


if __name__ == "__main__":
    main()
