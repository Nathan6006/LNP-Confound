#!/usr/bin/env python3
"""Digitize Figure 2b of Liu et al., Adv Healthc Mater 2026: library 1, 320 lipids.

The scale is CONTINUOUS (a colorbar running 0 to 6e4 RLU), so values are read off
exactly rather than binned, per the project rule for continuous heatmaps.

GEOMETRY is not guessed. The panel draws black gridlines between every cell, so
the grid is detected from the image: rows are the horizontal rules spanning the
body, columns the vertical ones. Detected 41 vertical lines at 348..3298 (40
columns, pitch 73.75) and 9 horizontal at 1082..2022 (8 rows, pitch 117.5).

  rows     8 amine cores, top to bottom: 70, 75, 80, 400, 113, 114, 110, T3A
  columns  8 BL groups L1..L8, each C6, C8, C10, C12, C14 (brackets under the
           axis label the groups; tick labels alternate between odd and even
           positions from group to group, which is why only three or two show)

SATURATION IS THE REAL HAZARD HERE. The colorbar stops at 6e4 but Figure 2c
plots the same lipids on a log axis reaching ~4e5, so every strong hit is
clipped to the top colour and its true value is NOT recoverable from this panel.
Those cells are emitted with `saturated=1` and a censored value, and Figure 2c
supplies the real number wherever it plots one. Treating a clipped cell as if it
read 6e4 would quietly compress the whole top of the dynamic range.
"""
import csv
import json
import os

import numpy as np
import PIL.Image as Image
import PIL.ImageDraw as ImageDraw

PAGE = "00_pages/fig_p4_600.png"
REVIEW = "review"

BODY_X = (348.0, 3298.0)        # 41 detected vertical gridlines
BODY_Y = (1082.0, 2022.0)       # 9 detected horizontal gridlines
NROW, NCOL = 8, 40

BAR_X = (3368, 3422)
# The LUT SAMPLING window: the first and last rows of the drawn colour ramp.
BAR_TOP, BAR_BOT = 1103, 1973   # outer rows of the colorbar's black border

# The VALUE anchors, which are NOT the same thing. Corrected 2026-10-30,
# review finding F-019 -- a second instance of error class 20, this time on a
# colorbar instead of an axis.
#
# The value scale used to be `VMAX + (VMIN-VMAX)*(y-BAR_TOP)/(BAR_BOT-BAR_TOP)`,
# i.e. it put 0 RLU at BAR_BOT = 1973, which is the OUTER edge of the bar's
# black border. The border runs y = 1962-1973 (verified from the pixel profile:
# dark from 1962, saturated 1964-1971), so its outer edge is about 5.5 px below
# the ramp's true zero. Stretching 0-6e4 over 870 px instead of 860 biases EVERY
# Fig 2b value HIGH -- about +1% at the 2e4 hit threshold and **+38% at the
# floor**, where it matters most because a third of library 1 sits there.
#
# These anchors come from the colorbar's own 4e4 and 2e4 tick notches,
# measured on the figure's native embedded JPEG (y = 1393.98 and 1680.60,
# 143.31 px per 1e4). Cross-checked by an independent route: the bar's bottom
# border CENTRE measures 1967.5 against the tick-derived 1967.23, agreeing to
# 0.27 px. A border-centre anchor pair (1104.5, 1967.5) gives a white-point
# value of 903.8 RLU against 913.7 here -- the two corrected routes differ by
# 9.9 RLU while the old anchor was +362 RLU out, so the correction is robust
# and the residual anchor uncertainty is ~3% of it.
#
# The top anchor is the weaker of the two and is stated as such: the bar's top
# border abuts a near-black dark-blue ramp, so its inner edge cannot be read
# from colour, which is why the ticks are the authority rather than the border.
Y_VMAX, Y_VMIN = 1107.4, 1967.6
VMAX, VMIN = 6.0e4, 0.0

AMINES = ["70", "75", "80", "400", "113", "114", "110", "T3A"]
BLS = [f"L{i}" for i in range(1, 9)]
TAILS = ["C6", "C8", "C10", "C12", "C14"]


def detect_grid(g):
    """Return the gridline positions actually present, to check the constants."""
    band = g[int(BODY_Y[0]) + 8:int(BODY_Y[1]) - 8]
    dark = (band < 120).mean(0)
    xs = [x for x in range(len(dark)) if dark[x] > 0.80]
    rowfrac = (g < 120).mean(1)
    ys = [y for y in range(g.shape[0]) if rowfrac[y] > 0.45]

    def group(v):
        out, cur = [], [v[0]]
        for x in v[1:]:
            if x - cur[-1] <= 4:
                cur.append(x)
            else:
                out.append(sum(cur) / len(cur))
                cur = [x]
        out.append(sum(cur) / len(cur))
        return out

    vx = [x for x in group(xs) if BODY_X[0] - 5 <= x <= BODY_X[1] + 5]
    vy = [y for y in group(ys) if BODY_Y[0] - 5 <= y <= BODY_Y[1] + 5]
    return vx, vy


def main():
    os.makedirs(REVIEW, exist_ok=True)
    img = np.array(Image.open(PAGE).convert("RGB")).astype(float)
    g = img.mean(2)

    vx, vy = detect_grid(g)
    assert len(vx) == NCOL + 1, f"found {len(vx)} vertical gridlines, want {NCOL+1}"
    assert len(vy) == NROW + 1, f"found {len(vy)} horizontal gridlines, want {NROW+1}"

    # Colorbar LUT. The sampling window skips 18 px at each end: the first
    # attempt sampled from +4 and picked up the border, making "top colour"
    # nearly black (0,1,46) so that no cell ever matched it and every saturated
    # cell was silently reported as a real reading.
    #
    # The value mapping uses Y_VMAX/Y_VMIN -- the TICK-derived anchors -- not
    # the border rows that bound this sampling window. Conflating the two is
    # F-019; see the constants above.
    ys = np.arange(BAR_TOP + 18, BAR_BOT - 17)
    lut = np.stack([np.median(img[y, BAR_X[0]:BAR_X[1]], axis=0) for y in ys])
    vals = VMAX + (VMIN - VMAX) * (ys - Y_VMAX) / (Y_VMIN - Y_VMAX)
    top_colour, bottom_colour = lut[0], lut[-1]

    rows, dists, nsat, nflo = [], [], 0, 0
    for r in range(NROW):
        for c in range(NCOL):
            y0, y1 = vy[r], vy[r + 1]
            x0, x1 = vx[c], vx[c + 1]
            # interior only: the gridlines are ~4 px, stay well clear of them
            a = slice(int(y0 + (y1 - y0) * 0.25), int(y0 + (y1 - y0) * 0.75))
            b = slice(int(x0 + (x1 - x0) * 0.25), int(x0 + (x1 - x0) * 0.75))
            col = np.median(img[a, b].reshape(-1, 3), axis=0)
            d = np.linalg.norm(lut - col, axis=1)
            i = int(d.argmin())
            v = float(vals[i])
            dists.append(float(d[i]))
            sat = int(np.linalg.norm(col - top_colour) < 14.0)
            flo = int(np.linalg.norm(col - bottom_colour) < 4.0)
            nsat += sat
            nflo += flo
            bl, tail = BLS[c // 5], TAILS[c % 5]
            rows.append({
                "lipid_name": f"{AMINES[r]}-{bl}-{tail}",
                "amine": AMINES[r], "bl": bl, "tail": tail,
                "panel_row": r, "panel_col": c,
                "RLU": round(v, 1),
                "saturated": sat, "floored": flo,
                "cell_rgb": f"{int(col[0])},{int(col[1])},{int(col[2])}",
                "match_distance": round(float(d[i]), 2),
                "source_figure": "Fig 2b",
                "extraction_method": "heatmap_colour_continuous_colorbar",
            })

    with open("01_fig2b_320.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    qc = {
        "paper": "Liu et al., Adv Healthc Mater 2026, 15, e02225, Fig 2b",
        "n_cells": len(rows), "expected": 320,
        "grid_detected": {"vertical_lines": len(vx), "horizontal_lines": len(vy),
                          "x": [round(vx[0], 1), round(vx[-1], 1)],
                          "y": [round(vy[0], 1), round(vy[-1], 1)],
                          "col_pitch": round((vx[-1] - vx[0]) / NCOL, 2),
                          "row_pitch": round((vy[-1] - vy[0]) / NROW, 2)},
        "colorbar": {"y_top": BAR_TOP, "y_bottom": BAR_BOT,
                     "value_top": VMAX, "value_bottom": VMIN},
        "colour_match": {"max_distance": round(max(dists), 2),
                         "median_distance": round(float(np.median(dists)), 2)},
        "saturated_cells": nsat,
        "floored_cells": nflo,
        "top_colour": [round(x, 1) for x in top_colour],
        "bottom_colour": [round(x, 1) for x in bottom_colour],
        "saturation_note": "cells at the colorbar maximum are censored at >=6e4; "
                           "Fig 2c plots the true value for the ones it shows",
        "value_range_unsaturated": [
            round(min(r["RLU"] for r in rows if not r["saturated"]), 1),
            round(max(r["RLU"] for r in rows if not r["saturated"]), 1)],
        "hit_rate_check_gt_2e4": sum(1 for r in rows if r["RLU"] > 2e4 or r["saturated"]),
    }
    json.dump(qc, open("02_fig2b_qc.json", "w"), indent=1)
    overlay(img, rows, vx, vy)

    print(f"cells {len(rows)} (expect 320)")
    print(f"  grid: {len(vx)} vertical / {len(vy)} horizontal lines detected, "
          f"pitch {qc['grid_detected']['col_pitch']} x {qc['grid_detected']['row_pitch']}")
    print(f"  colour match: median {qc['colour_match']['median_distance']}, "
          f"max {qc['colour_match']['max_distance']}")
    print(f"  saturated at >=6e4: {nsat} of {len(rows)}   floored at <=~1e3: {nflo}")
    print(f"  unsaturated range {qc['value_range_unsaturated']}")
    print(f"  overlay -> {REVIEW}/fig2b_overlay.png")


def overlay(img, rows, vx, vy):
    x0, x1, y0, y1 = vx[0], vx[-1], vy[0], vy[-1]
    pad = 60
    crop = Image.fromarray(img.astype("uint8")[int(y0) - pad:int(y1) + pad,
                                               int(x0) - pad:int(x1) + pad])
    d = ImageDraw.Draw(crop)
    ox, oy = x0 - pad, y0 - pad
    for r in range(len(vy) - 1):
        for c in range(len(vx) - 1):
            rr = rows[r * (len(vx) - 1) + c]
            a0 = vy[r] + (vy[r + 1] - vy[r]) * 0.25
            a1 = vy[r] + (vy[r + 1] - vy[r]) * 0.75
            b0 = vx[c] + (vx[c + 1] - vx[c]) * 0.25
            b1 = vx[c] + (vx[c + 1] - vx[c]) * 0.75
            colour = (220, 0, 0) if rr["saturated"] else (0, 170, 0)
            d.rectangle([b0 - ox, a0 - oy, b1 - ox, a1 - oy], outline=colour, width=3)
    d.rectangle([x0 - ox, y0 - oy, x1 - ox, y1 - oy], outline=(0, 90, 255), width=5)
    d.text((8, 8), "blue = detected body   green box = pixels sampled   "
                   "red = saturated at the colorbar maximum", fill=(0, 0, 0))
    crop.save(f"{REVIEW}/fig2b_overlay.png")


if __name__ == "__main__":
    main()
