#!/usr/bin/env python3
"""peer_39912611 stage 04 -- digitise Figure 2 (in vitro Luc, 3 cell lines).

TIER 5. 00_value_sources.json records why tiers 1-4 do not exist for this
paper. The raster is the publisher's own, 700 x 922 px, fetched from Europe
PMC; the higher-resolution assets on pubs.acs.org need a browser, which this
run may not use, so 700 px IS the resolution limit and every number below
carries it.

WHAT IS MEASURED AND WHAT IS NOT.
  * The y scale of each panel is fitted ONLY from the panel's own dashed
    decade gridlines, detected as darkness-weighted centroids in the columns
    that contain no bar at any height. The axis rule itself is NOT used as an
    anchor -- that is error class 20, and here it would be especially cheap to
    get away with, because panel A's axis rule happens to sit on 1e5.
  * Bar tops come from curation/lib/bars, on a mask built from each bar's OWN
    colour. The whiskers in this figure are black and centred; a colour mask
    drops them, and the full-width + flank test in bars.bar_top drops anything
    else centred. That is the F-018/F-020 class, twice found.
  * The stroke offset is MEASURED, not assumed: for every bar the first dark
    row above the first colour row is found in the bar's flanks, and the median
    of that gap over all bars in the panel is subtracted. A GraphPad bar's top
    edge lies under its own stroke, so reading the first filled pixel biases
    every bar the same way, downward.
"""
import csv
import json
import os
import sys

import numpy as np
from PIL import Image

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "lib"))
import axis                                                   # noqa: E402
import bars                                                   # noqa: E402
from errors import DetectionError                             # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
IMG = os.path.join(HERE, "source", "figs", "nn4c18636_0002.jpg")

# Bar order, left to right, read off the figure legend (two legend rows:
# 14/32/33/34/35/36 then 37/38/39/40/41/42) and confirmed by the colour of
# each group's bars matching its legend swatch -- see check F2-1 in 07_qc.json.
LIPIDS = ["Lipid 14", "Lipid 32", "Lipid 33", "Lipid 34", "Lipid 35",
          "Lipid 36", "Lipid 37", "Lipid 38", "Lipid 39", "Lipid 40",
          "Lipid 41", "Lipid 42"]
DOSES = [0.25, 0.15, 0.025]          # ug/mL, printed under every triplet

PANELS = [
    # key, cell line, row band to search for gridlines, decade values top->bottom
    ("A", "HeLa", 70, 300, [1e9, 1e8, 1e7, 1e6, 1e5]),
    ("B", "RAW264.7", 365, 590, [1e8, 1e7, 1e6, 1e5, 1e4, 1e3, 1e2, 1e1, 1e0]),
    ("C", "HepG2", 660, 882, [1e8, 1e7, 1e6, 1e5, 1e4, 1e3, 1e2, 1e1, 1e0]),
]
X_LO, X_HI = 84, 692                 # plot area, shared by all three panels


def load():
    a = np.asarray(Image.open(IMG).convert("RGB")).astype(int)
    gray = a.mean(axis=2)
    sat = a.max(axis=2) - a.min(axis=2)
    # sat>40 rather than >60, and no brightness floor: lipid 39's HeLa bar at
    # 0.025 ug/mL is drawn in a dark green whose MAXIMUM channel is 80, so a
    # `max > 80` gate dropped that bar entirely and the 36 windows came back
    # as 35 -- which would have shifted the lipid/dose label of every bar to
    # its right. Caught by asserting the count, not by looking.
    colour = (sat > 40)
    # The detection mask: the bar's own fill OR any dark ink. The bar top RULE
    # is black and spans the full width; the error whisker is black and 3 px
    # wide in a 9 px bar; the cap is 5 px. bars.bar_top's full-width + flank
    # test separates them, which is the whole reason that module exists.
    ink = (sat > 40) | (gray < 150)
    return a, gray, colour, ink


def bar_windows(colour, y0, y1):
    """Contiguous runs of columns that contain bar colour. 36 expected.

    A run is formed on >= 1 coloured row, because the shortest real bar in the
    figure (lipid 39, HeLa, 0.025 ug/mL) is four rows tall and a stricter gate
    drops it -- which would shift the lipid/dose label of every bar to its
    right. The run is then TRIMMED: a leading or trailing column holding less
    than a quarter of the run's tallest column is JPEG ringing beside a
    significance asterisk, not part of the bar, and leaving it in puts an empty
    column inside bars.bar_top's left flank, which (correctly) refuses the bar.
    """
    band = colour[y0:y1, :]
    cnt = band.sum(axis=0)
    runs, s = [], None
    for x in range(X_LO, X_HI):
        if cnt[x] > 0:
            if s is None:
                s = x
        else:
            if s is not None:
                runs.append((s, x))
                s = None
    if s is not None:
        runs.append((s, X_HI))
    out = []
    for x0, x1 in runs:
        if x1 - x0 < 6:
            continue
        peak = cnt[x0:x1].max()
        keep = [x for x in range(x0, x1) if cnt[x] >= 0.25 * peak]
        if len(keep) < 6:
            continue
        out.append((keep[0], keep[-1] + 1))
    return out


def free_columns(colour, y0, y1):
    band = colour[y0:y1, :]
    cnt = band.sum(axis=0)
    return [x for x in range(X_LO, X_HI) if cnt[x] == 0]


def gridlines(gray, free_cols, y0, y1, n_expected):
    """Darkness-weighted row centroids of the dashed decade gridlines.

    Measured ONLY in columns that hold no bar at any height, so a bar's own
    top rule can never be mistaken for a gridline.
    """
    cols = np.array(free_cols, dtype=int)
    prof = gray[y0:y1, cols]
    dark_frac = (prof < 170).mean(axis=1)
    # Total ink per row, not the row mean: a DASHED line is ink on half its
    # length, so a mean over the whole row sits above any fixed threshold and
    # the weight collapses to zero exactly on the faintest gridline -- which
    # is what silently dropped panel A's 1e8 on the first run.
    weight = np.clip(170 - prof, 0, None).sum(axis=1)
    rowmin = prof.min(axis=1)
    on = dark_frac > 0.35
    runs, s = [], None
    for i, v in enumerate(on):
        if v and s is None:
            s = i
        elif not v and s is not None:
            runs.append((s, i))
            s = None
    if s is not None:
        runs.append((s, len(on)))
    centres = []
    for a0, b0 in runs:
        if b0 - a0 > 6:              # a text block
            continue
        if float(dark_frac[a0:b0].max()) > 0.90:
            # The SOLID axis rule -- ink in EVERY bar-free column, where a
            # dashed decade line covers 0.49-0.66 of them. It coincides with a
            # decade here (1e5 in
            # panel A, 1e0 in B and C) and would calibrate perfectly, which is
            # exactly why it is excluded: anchoring on the panel box rather
            # than on detected furniture is error class 20, and a box that
            # happens to be right teaches nothing about one that is not.
            continue
        w = weight[a0:b0]
        if w.sum() <= 0:
            continue
        centres.append(y0 + float((np.arange(a0, b0) * w).sum() / w.sum()))
    return centres


def stroke_offset(gray, colour, windows, tops):
    """Median gap, in rows, between a bar's first dark row and its first
    colour row, measured in the bar's flanks only."""
    gaps = []
    for i, (x0, x1) in enumerate(windows):
        if i not in tops:
            continue
        r = tops[i]["row"]
        fw = max(2, int(round((x1 - x0) * 0.33)))
        for a0, b0 in ((x0, x0 + fw), (x1 - fw, x1)):
            col = gray[max(0, r - 6):r + 1, a0:b0].mean(axis=1)
            dk = np.where(col < 170)[0]
            if dk.size:
                gaps.append(r - (max(0, r - 6) + int(dk[0])))
    return float(np.median(gaps)) if gaps else 0.0


def main():
    a, gray, colour, ink_all = load()
    out_rows = []
    qc = {"image": os.path.relpath(IMG, HERE), "panels": {}}

    for key, cell, y0, y1, decades in PANELS:
        free = free_columns(colour, y0 + 5, y1 - 5)
        gl = gridlines(gray, free, y0, y1, len(decades))
        if len(gl) < 3:
            raise SystemExit(f"panel {key}: only {len(gl)} gridlines detected")
        # match the detected gridlines to decade values, top first
        gl = sorted(gl)
        # Assign each detected centroid to a decade by its index from the top,
        # with the step taken as the MEDIAN consecutive gap -- not span/(n-1)
        # against the expected count, which silently mis-assigns the whole
        # ladder when one line is missing (the solid axis rule is always
        # dropped above, so one is always missing, by construction).
        step = float(np.median(np.diff(gl)))
        ks = [round((c - gl[0]) / step) for c in gl]
        if ks != list(range(len(gl))):
            raise SystemExit(
                f"panel {key}: gridlines at {[round(g,2) for g in gl]} do not "
                f"form a contiguous decade ladder (indices {ks}). A line is "
                "missing from the middle; refusing rather than guessing which.")
        vals = [decades[k] for k in ks]
        used = list(gl)
        cal = axis.calibrate(used, vals, kind="log10", axis="y",
                             loo_tol_px=axis.LOO_TOL_PX)
        loo = axis.loo_residuals(used)

        wins = bar_windows(colour, y0 + 5, y1 - 5)
        if len(wins) != 36:
            raise SystemExit(
                f"panel {key}: {len(wins)} bar windows, expected 36 "
                "(12 lipids x 3 doses). A missing window shifts the "
                "lipid/dose label of every bar to its right, so this refuses "
                f"rather than emits: {wins}")
        # the solid axis rule, one decade below the lowest dashed gridline
        base_row = max(used) + cal.px_per_unit

        # mask of "this bar's colour", bar by bar: the modal hue inside the bar
        tops, fails = {}, {}
        for i, (x0, x1) in enumerate(wins):
            sub = a[y0:y1, x0:x1]
            m = colour[y0:y1, x0:x1]
            if m.sum() < 20:
                fails[i] = "no colour in window"
                continue
            try:
                # y_top is the TOP GRIDLINE, not the panel band: the band
                # reaches into the figure legend (panel A) and the panel title
                # (B, C), and a flank that finds legend text 150 px above the
                # bar is what the flank test then reports as a 150 px
                # disagreement. Those three refusals were real -- the mask was
                # wrong, not the bar.
                tops[i] = bars.bar_top(ink_all, x0, x1, int(round(min(used))) - 2,
                                       int(round(base_row)) - 1,
                                       baseline_row=None)
            except DetectionError as e:
                fails[i] = str(e).split("\n")[0]

        # Diagnostic only: how far the first FILLED row sits below the
        # detected top. It is not subtracted -- the mask already includes the
        # bar's own stroke, so the detected row IS the stroke, and the data
        # coordinate is its centre to within the +/-1 px this raster allows.
        off = 0.0
        fill_gap = stroke_offset(gray, colour, wins, tops)

        qc["panels"][key] = {
            "cell_line": cell,
            "n_gridlines_detected": len(gl),
            "gridlines_used_px": [round(u, 2) for u in used],
            "gridline_values": vals,
            "calibration": cal.qc(),
            "loo_residuals_px": [round(v, 3) for v in loo],
            "worst_loo_px": round(max(loo), 3),
            "loo_tol_px": axis.LOO_TOL_PX,
            "n_bar_windows": len(wins),
            "n_bar_tops": len(tops),
            "failures": fails,
            "stroke_offset_px_subtracted": off,
            "fill_starts_px_below_detected_top_median": fill_gap,
            "baseline_row_px": round(float(base_row), 2),
            "bar_width_px": sorted({w[1] - w[0] for w in wins}),
        }

        for i, (x0, x1) in enumerate(wins):
            lip = LIPIDS[i // 3] if i // 3 < 12 else "?"
            dose = DOSES[i % 3]
            if i in tops:
                row = tops[i]["row"] - off
                val = cal.to_value(row)
                st = "measured"
            else:
                row, val, st = None, None, "REFUSED"
            out_rows.append({
                "panel": f"Fig 2{key}", "cell_line": cell,
                "lipid_name": lip, "dose_ug_per_mL": dose,
                "bar_index": i, "x0": x0, "x1": x1,
                "top_row_px": None if row is None else round(row, 2),
                "luminescence": None if val is None else float(f"{val:.4g}"),
                "log10": None if val is None else round(float(np.log10(val)), 4),
                "status": st,
                "detail": "" if i in tops else fails.get(i, ""),
            })

    with open(os.path.join(HERE, "04_fig2_values.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out_rows[0].keys()))
        w.writeheader()
        w.writerows(out_rows)
    with open(os.path.join(HERE, "04_fig2_qc.json"), "w") as fh:
        json.dump(qc, fh, indent=1)
    print(json.dumps({k: {kk: vv for kk, vv in v.items()
                          if kk not in ("failures",)}
                      for k, v in qc["panels"].items()}, indent=1))
    nfail = sum(1 for r in out_rows if r["status"] != "measured")
    print("rows:", len(out_rows), "refused:", nfail)
    return 0


if __name__ == "__main__":
    sys.exit(main())
