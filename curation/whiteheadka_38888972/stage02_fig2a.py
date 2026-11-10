#!/usr/bin/env python3
"""Figure 2A -> 06_tns.csv : TNS fluorescence at pH 5 for all 580 members.

THIS IS NOT A DELIVERY MEASUREMENT AND IS NEVER EMITTED INTO THE CORPUS.
TNS fluorescence reports LNP surface ionization -- the same class as diameter,
zeta potential and uptake, which the corpus's own definition excludes (the
2026-07-23 rebuild: "same delivery-only definition, excl. diameter/zeta/
uptake"). 00_value_sources.json records it as such and DECISIONS.md D-3 queues
the question rather than answering it. What this stage IS for:

  1. it measures the grid and so turns "a library of 580" from a sentence into
     a counted 29 x 20;
  2. it supplies the TNS half of the Figure 2C cross-check, which compares the
     (TNS, in-vivo flux) pairs we built from Figures 2A and 2B against the
     point cloud the authors drew in Figure 2C -- a check that needs no
     attribution and is external to both digitisers;
  3. it tests the paper's own printed structural claim about which tails give
     the highest and lowest TNS.

WHY THIS IS WEAKER EVIDENCE THAN A DISCRETE HEATMAP, STATED UP FRONT.
CLAUDE.md section 5 ranks a CONTINUOUS colormap below a discrete-binned one,
and this is a continuous white->blue ramp. There is no legend swatch set to
nearest-neighbour against, so the lookup table has to be sampled off the
colorbar in the same raster and inverted. Every value here carries that.

THE GRID IS FITTED, NOT READ OFF THE BORDER. F-019 is a scale anchored on the
outer edge of a panel's own black border. Here the 19 internal column
boundaries and 28 internal row boundaries are DETECTED as colour
discontinuities and a uniform grid is fitted to them with axis.calibrate, whose
residual is recorded. The box border is never used as an anchor.

Run: python3 -I stage02_fig2a.py
"""
import csv
import json
import os
import sys

import numpy as np
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))

import axis                                                     # noqa: E402
import heatmap                                                  # noqa: E402

FIG = os.path.join(HERE, "source", "nihms-2006250-f0002.jpg")

BOX = {"x0": 154, "x1": 1158, "y0": 154, "y1": 1154}   # inside the rule
CBAR = {"x0": 1182, "x1": 1215, "y0": 154, "y1": 1154}
# colorbar notch rows, detected below; the values they carry, top to bottom
CBAR_TICK_VALUES = [30000.0, 20000.0, 10000.0]

ROWS = ["60", "80", "81", "110", "111", "113", "122", "123", "174", "200",
        "205", "211", "303", "304", "305", "306", "313", "347", "400", "401",
        "402", "500", "501", "503", "509", "514", "516", "601", "602"]
COLS = ["O10", "O11", "O12", "O13", "O14", "O16", "O18", "Oi10", "O8,12",
        "O6,10", "O1,1,8", "O4,8", "X1", "X2", "X3", "X4", "X5", "X6", "X7",
        "X8"]


def boundaries(rgb, lo, hi, n, along):
    """The n-1 internal cell boundaries, as colour-discontinuity peaks.

    Measured from the image's own content. The caller passes a window inset
    from the box rule by MARGIN px, because the rule's own anti-aliased inner
    edge is itself a large colour discontinuity: the first attempt detected
    x = 154.63 and x = 1156.19 -- both border edges -- and those two false
    peaks crowded out two real boundaries, which axis.calibrate then caught as
    a 26.78 px residual over 19 ticks. A border must never become an anchor
    (F-019), so it is excluded from the search rather than ranked against.
    """
    a = rgb[BOX["y0"]:BOX["y1"], lo:hi] if along == "x" else \
        rgb[lo:hi, BOX["x0"]:BOX["x1"]]
    prof = a.mean(axis=0) if along == "x" else a.mean(axis=1)   # (m, 3)
    d = np.abs(np.diff(prof, axis=0)).sum(axis=1)               # (m-1,)
    want = n - 1
    pitch = (hi - lo) / float(n - 1)
    picks = []
    order = np.argsort(d)[::-1]
    for i in order:
        if len(picks) == want:
            break
        if all(abs(i - p) > pitch * 0.5 for p in picks):
            picks.append(int(i))
    if len(picks) != want:
        raise SystemExit(f"boundaries({along}): found {len(picks)} of {want}")
    # sub-pixel: intensity-weighted centroid over the 3 px around each peak
    out = []
    for p in sorted(picks):
        s = slice(max(0, p - 1), min(len(d), p + 2))
        w = d[s]
        out.append(lo + 0.5 + float((np.arange(s.start, s.stop) * w).sum()
                                    / w.sum()))
    return out


def main():
    img = Image.open(FIG).convert("RGB")
    rgb = np.asarray(img).astype(float)
    qc = {}

    # --- 1. the grid, fitted to detected boundaries ---------------------
    M = 8   # inset from the box rule; see boundaries() for why
    xb = boundaries(rgb, BOX["x0"] + M, BOX["x1"] - M, len(COLS), "x")
    yb = boundaries(rgb, BOX["y0"] + M, BOX["y1"] - M, len(ROWS), "y")
    calx = axis.calibrate(xb, list(range(1, len(COLS))), kind="linear",
                          axis="x")
    caly = axis.calibrate(yb, list(range(1, len(ROWS))), kind="linear",
                          axis="x")   # index grows with pixel row
    qc["column_grid"] = calx.qc()
    qc["row_grid"] = caly.qc()
    qc["n_columns_from_boundaries"] = len(xb) + 1
    qc["n_rows_from_boundaries"] = len(yb) + 1
    if len(xb) + 1 != 20 or len(yb) + 1 != 29:
        raise SystemExit("fig2a: the grid is not 29 x 20")
    qc["grid_cells"] = (len(xb) + 1) * (len(yb) + 1)

    xs = [calx.to_pixel(i - 0.5) for i in range(1, len(COLS) + 1)]
    ys = [caly.to_pixel(i - 0.5) for i in range(1, len(ROWS) + 1)]

    # --- 2. the colorbar: notches, then the ramp ------------------------
    g = rgb.mean(axis=2)
    inner = g[CBAR["y0"]:CBAR["y1"], CBAR["x0"]:CBAR["x0"] + 8].mean(axis=1)
    centre = g[CBAR["y0"]:CBAR["y1"],
               CBAR["x0"] + 13:CBAR["x0"] + 21].mean(axis=1)
    diff = centre - inner
    idx = np.flatnonzero(diff > 25)
    groups, cur = [], [idx[0]]
    for i in idx[1:]:
        if i - cur[-1] <= 2:
            cur.append(i)
        else:
            groups.append(cur)
            cur = [i]
    groups.append(cur)
    notch = []
    for gp in groups:
        w = diff[gp]
        notch.append(float((np.asarray(gp) * w).sum() / w.sum()) + CBAR["y0"])
    if len(notch) != len(CBAR_TICK_VALUES):
        raise SystemExit(f"fig2a: {len(notch)} colorbar notches, expected 3: "
                         f"{np.round(notch, 2).tolist()}")
    calc = axis.calibrate(notch, CBAR_TICK_VALUES, kind="linear", axis="y")
    qc["colorbar_axis"] = calc.qc()
    qc["colorbar_notch_rows"] = [round(n, 3) for n in notch]
    qc["colorbar_notches_are_drawn_not_inferred"] = (
        "The three notches are short dark stubs the authors drew INTO the bar "
        "at 10,000 / 20,000 / 30,000. They are detected as a contrast between "
        "the bar's inner edge and its centre, so neither the bar's border nor "
        "its printed labels enters the fit.")

    # --- 3. the lookup table, sampled from the SAME raster --------------
    lut_rows = np.arange(CBAR["y0"] + 3, CBAR["y1"] - 3)
    lut_rgb = np.array([np.median(
        rgb[r - 1:r + 2, CBAR["x0"] + 6:CBAR["x1"] - 6].reshape(-1, 3), axis=0)
        for r in lut_rows])
    # Dedupe: heatmap.match refuses an exact tie, and says to dedupe the LUT
    # when the tie is between duplicate swatches. A continuous ramp sampled
    # row by row is full of duplicates.
    keyed = {}
    for r, c in zip(lut_rows, lut_rgb):
        keyed.setdefault(tuple(np.round(c, 3)), []).append(r)
    lut_c = np.array([list(k) for k in keyed])
    lut_r = np.array([float(np.mean(v)) for v in keyed.values()])
    qc["lut_rows_sampled"] = int(lut_rows.size)
    qc["lut_unique_colours"] = int(lut_c.shape[0])

    # --- 4. the cells ---------------------------------------------------
    gx, gy = np.meshgrid(xs, ys)
    cells = heatmap.sample(img, gx.ravel(), gy.ravel(), half=4)
    m = heatmap.match(cells, lut_c, max_dist=heatmap.MAX_DIST)
    vals = calc.to_value(lut_r[m["index"]])
    qc["cell_match"] = {
        "max_nearest_lut_distance": round(m["max_distance"], 3),
        "max_dist_threshold": heatmap.MAX_DIST,
        "n_within_10_of_a_second_lut_entry": m["n_ambiguous"],
        "_note": "On a CONTINUOUS ramp the 'ambiguous' count is not a defect "
                 "signal the way it is for a discrete legend -- adjacent LUT "
                 "entries are a fraction of an RGB unit apart by construction. "
                 "The number that matters is the max nearest-LUT distance: "
                 "CLAUDE.md section 5 says anything above ~5 means the "
                 "GEOMETRY is wrong, not the colours.",
    }

    rows = []
    for i, h in enumerate(ROWS):
        for j, t in enumerate(COLS):
            k = i * len(COLS) + j
            rows.append({
                "IL_name": f"{h}{t}", "head": h, "tail": t,
                "row": i + 1, "col": j + 1,
                "x_px": round(float(gx.ravel()[k]), 2),
                "y_px": round(float(gy.ravel()[k]), 2),
                "R": round(float(cells[k][0]), 1),
                "G": round(float(cells[k][1]), 1),
                "B": round(float(cells[k][2]), 1),
                "lut_distance": round(float(m["distance"][k]), 3),
                "TNS_AU": round(float(vals[k]), 1),
            })
    with open(os.path.join(HERE, "06_tns.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow(r)

    v = np.array([r["TNS_AU"] for r in rows]).reshape(len(ROWS), len(COLS))
    colmean = {t: round(float(v[:, j].mean()), 1) for j, t in enumerate(COLS)}
    qc["n_cells"] = len(rows)
    qc["TNS_range_AU"] = [round(float(v.min()), 1), round(float(v.max()), 1)]
    qc["column_mean_TNS_AU"] = colmean
    order = sorted(colmean, key=colmean.get)
    qc["columns_lowest_to_highest"] = order
    with open(os.path.join(HERE, "06_tns_qc.json"), "w") as fh:
        json.dump(qc, fh, indent=1, sort_keys=True)
    print(json.dumps({k: qc[k] for k in
                      ("column_grid", "row_grid", "grid_cells",
                       "colorbar_axis", "colorbar_notch_rows",
                       "lut_unique_colours", "cell_match", "n_cells",
                       "TNS_range_AU", "columns_lowest_to_highest")},
                     indent=1, sort_keys=True))


if __name__ == "__main__":
    main()
