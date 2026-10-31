#!/usr/bin/env python3
"""Digitize Figure 2c of Tian et al., Adv Mater 2026: the 486-LNP formulation screen.

The 486 are NOT 486 lipids. They are one lipid, C6O2B2, formulated 486 ways:

    C6O2B2 mol%   35 / 50 / 65          (3)
    DODAP  mol%   10 / 30 / 50          (3)
    helper mol%   10 / 30 / 50          (3)
    DMG-PEG mol%   1 /  2 /  4          (3)
    helper identity  DOPC DMPC DOPG DOPE EggSM DSPC   (6)

    3 * 3 * 3 * 3 * 6 = 486

shown as six 9x9 heatmaps (panels A-F, one per helper lipid). Rows are
C6O2B2 (outer, 3 groups) x DODAP (inner, 3); columns are helper mol% (outer,
3 groups) x DMG-PEG (inner, 3).

WHY THIS IS A SAFE DIGITIZATION. The legend is DISCRETE, not a gradient:
four bins, <1e4 / 1e4-1e5 / 1e5-1e6 / >1e6. So every cell is one of exactly
four printed colours, and reading it is nearest-neighbour matching against
swatches sampled from the legend IN THE SAME RASTER - no interpolation, no
inferred scale. This is the same situation as the iPhos Fig 2a heatmap, which
scored Spearman 1.000 against the curated corpus.

Geometry is taken from the cell gridlines themselves (long dark runs), not
from an assumed bounding box, then filtered to the 10 lines whose spacings are
most uniform. Cell colour is the median over the inner 50% of each cell, so
gridline pixels and antialiasing cannot swing the match.

FIVE INDEPENDENT CHECKS the paper hands us, all applied in verify():
  text  "Out of 486 LNP formulations, 240 met the threshold for acceptable
         performance, with 91 LNPs categorized as outstanding"
  2d    outstanding count by C6O2B2 level
  2e    outstanding count by DMG-PEG level
  2f    outstanding count by DODAP level
  2g    outstanding count by helper lipid
A digitization that reproduces four marginals and the total is confirmed in
aggregate, not merely asserted.
"""
import csv
import json
import os
import sys

import numpy as np
import PIL.Image as Image

PAGE_PNG = sys.argv[1] if len(sys.argv) > 1 else "00_pages/fig2_page.png"
OUT_CSV = "04b_fig2c_486.csv"
OUT_QC = "07b_fig2c_qc.json"

# legend swatch sample boxes, in page-raster pixels at 300 dpi
LEGEND_Y = (1556, 1572)
LEGEND_X = [(1640, 1740), (1780, 1880), (1915, 2020), (2040, 2155)]
BINS = [
    ("<1e4", 1, 3.16e3),          # geometric midpoint of the open bin, see note
    ("1e4-1e5", 2, 3.16e4),
    ("1e5-1e6", 3, 3.16e5),
    (">1e6", 4, 3.16e6),
]
# midpoints are geometric centres of each decade band; the open-ended outer bins
# are given the same half-decade offset. They exist only so the column is
# numerically usable - the DEFENSIBLE quantity is the ordinal bin_index.

# generous search windows, (y0, y1, x0, x1), in page-raster pixels
PANELS = [
    ("DOPC", (185, 640), (1240, 1715)),
    ("DMPC", (185, 640), (1715, 2280)),
    ("DOPG", (643, 1090), (1240, 1715)),
    ("DOPE", (643, 1090), (1715, 2280)),
    ("EggSM", (1099, 1545), (1240, 1715)),
    ("DSPC", (1099, 1545), (1715, 2280)),
]

C6O2B2_LEVELS = [35, 50, 65]     # row outer
DODAP_LEVELS = [10, 30, 50]      # row inner
HELPER_LEVELS = [10, 30, 50]     # col outer
PEG_LEVELS = [1, 2, 4]           # col inner

# values read off the printed bar charts, for verification only
REPORTED = {
    "total_acceptable": 240,
    "total_outstanding": 91,
    "by_c6o2b2": {35: 41, 50: 21, 65: 31},      # Fig 2d
    "by_peg": {1: 35, 2: 29, 4: 29},            # Fig 2e
    "by_dodap": {10: 8, 30: 31, 50: 54},        # Fig 2f
    "by_helper": {"DOPC": 3, "DMPC": 17, "DOPG": 16,
                  "DOPE": 9, "EggSM": 13, "DSPC": 35},   # Fig 2g
}


def _tot(r, c):
    """Sum of the four molar ratios for this cell."""
    return (C6O2B2_LEVELS[r // 3] + DODAP_LEVELS[r % 3]
            + HELPER_LEVELS[c // 3] + PEG_LEVELS[c % 3])


def long_runs(counts, thresh):
    """Centres of contiguous runs where counts exceed thresh."""
    out, start = [], None
    for i, v in enumerate(counts):
        if v > thresh and start is None:
            start = i
        elif v <= thresh and start is not None:
            out.append((start + i - 1) / 2.0)
            start = None
    if start is not None:
        out.append((start + len(counts) - 1) / 2.0)
    return out


def most_uniform(lines, n=10):
    """Keep the n lines whose consecutive spacings are most uniform.

    Tick brackets and axis rules sit right beside the heatmap border and survive
    any brightness threshold that also catches the faint internal gridlines of the
    pale panels. Rather than hand-tune the threshold per panel, over-detect and
    then drop whichever line's removal most reduces the spacing variance.
    """
    lines = sorted(lines)
    while len(lines) > n:
        best, best_var = None, None
        for i in range(len(lines)):
            trial = lines[:i] + lines[i + 1:]
            d = np.diff(trial)
            v = float(np.var(d))
            if best_var is None or v < best_var:
                best, best_var = i, v
        lines = lines[:best] + lines[best + 1:]
    return lines


def find_grid(gray, win):
    (y0, y1), (x0, x1) = win
    reg = gray[y0:y1, x0:x1]
    dark = reg < 185
    rs, cs = dark.sum(axis=1), dark.sum(axis=0)
    h = most_uniform([v + y0 for v in long_runs(rs, rs.max() * 0.55)])
    v = most_uniform([v + x0 for v in long_runs(cs, cs.max() * 0.55)])
    return h, v


def main():
    img = np.array(Image.open(PAGE_PNG).convert("RGB")).astype(int)
    gray = img.mean(axis=2)

    pal = np.array([np.median(img[LEGEND_Y[0]:LEGEND_Y[1], a:b].reshape(-1, 3), axis=0)
                    for a, b in LEGEND_X])

    rows, geom = [], {}
    for helper, ywin, xwin in PANELS:
        h, v = find_grid(gray, (ywin, xwin))
        dh, dv = np.diff(h), np.diff(v)
        geom[helper] = {
            "h_lines": [round(x, 1) for x in h],
            "v_lines": [round(x, 1) for x in v],
            "row_pitch_px": round(float(dh.mean()), 2),
            "col_pitch_px": round(float(dv.mean()), 2),
            "row_pitch_cv": round(float(dh.std() / dh.mean()), 4),
            "col_pitch_cv": round(float(dv.std() / dv.mean()), 4),
        }
        for r in range(9):
            for c in range(9):
                ya, yb = h[r], h[r + 1]
                xa, xb = v[c], v[c + 1]
                iy0 = int(round(ya + 0.25 * (yb - ya)))
                iy1 = int(round(ya + 0.75 * (yb - ya)))
                ix0 = int(round(xa + 0.25 * (xb - xa)))
                ix1 = int(round(xa + 0.75 * (xb - xa)))
                patch = img[iy0:iy1, ix0:ix1].reshape(-1, 3)
                col = np.median(patch, axis=0)
                d = np.linalg.norm(pal - col, axis=1)
                k = int(np.argmin(d))
                label, idx, mid = BINS[k]
                rows.append({
                    "formulation_id": f"{helper}_{C6O2B2_LEVELS[r // 3]}"
                                      f"_{DODAP_LEVELS[r % 3]}"
                                      f"_{HELPER_LEVELS[c // 3]}"
                                      f"_{PEG_LEVELS[c % 3]}",
                    "lipid_name": "C6O2B2",
                    "helper_lipid": helper,
                    # These are molar RATIOS as the paper states them, not
                    # percentages: Table S1 prints both and the ratio rows sum to
                    # 74-167, never 100. Naming them _molpct (as an earlier
                    # version did) invites a consumer to read 35 as 35 mol%,
                    # which misstates every formulation. The percentages are
                    # derived alongside so both are available unambiguously.
                    "C6O2B2_molratio": C6O2B2_LEVELS[r // 3],
                    "DODAP_molratio": DODAP_LEVELS[r % 3],
                    "helper_molratio": HELPER_LEVELS[c // 3],
                    "DMGPEG_molratio": PEG_LEVELS[c % 3],
                    "C6O2B2_molpct": round(100 * C6O2B2_LEVELS[r // 3] / _tot(r, c), 2),
                    "DODAP_molpct": round(100 * DODAP_LEVELS[r % 3] / _tot(r, c), 2),
                    "helper_molpct": round(100 * HELPER_LEVELS[c // 3] / _tot(r, c), 2),
                    "DMGPEG_molpct": round(100 * PEG_LEVELS[c % 3] / _tot(r, c), 2),
                    "panel_row": r, "panel_col": c,
                    "bin_label": label,
                    "bin_index": idx,
                    "rlu_midpoint": mid,
                    "cell_rgb": f"{int(col[0])},{int(col[1])},{int(col[2])}",
                    "match_distance": round(float(d[k]), 2),
                    "runner_up_margin": round(float(np.sort(d)[1] - d[k]), 2),
                    "assay": "luciferase, IGROV-1 cells, 20 ng mLuc per well, 24 h",
                    "readout_metric": "relative light units (RLU)",
                    "value_scale": "binned",
                    "n_replicates": 3,
                    "source_figure": "Fig 2c",
                    "extraction_method": "heatmap_color_discrete_4bin",
                    "source_pdf": "Adv Mater 2026 Tian (main)",
                })

    with open(OUT_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    qc = verify(rows, geom, pal)
    json.dump(qc, open(OUT_QC, "w"), indent=1)
    overlay(img, geom, rows, pal)
    report(rows, qc)


def overlay(img, geom, rows, pal):
    """Write review/fig2c_overlay.png: every sampled patch, outlined and tinted
    with the bin colour it was matched to, so the read can be checked by eye."""
    import PIL.ImageDraw as ImageDraw
    os.makedirs("review", exist_ok=True)
    y0, y1, x0, x1 = 140, 1600, 1150, 2300
    crop = Image.fromarray(img.astype("uint8")[y0:y1, x0:x1])
    S = 2
    crop = crop.resize((crop.width * S, crop.height * S), Image.LANCZOS)
    d = ImageDraw.Draw(crop)
    ring = [(0, 90, 200), (120, 120, 120), (230, 140, 0), (200, 0, 0)]
    by = {}
    for r in rows:
        by[(r["helper_lipid"], r["panel_row"], r["panel_col"])] = r
    for helper, _, _ in PANELS:
        g = geom[helper]
        h, v = g["h_lines"], g["v_lines"]
        d.rectangle([(v[0] - x0) * S, (h[0] - y0) * S,
                     (v[-1] - x0) * S, (h[-1] - y0) * S],
                    outline=(0, 160, 0), width=2)
        for r in range(9):
            for c in range(9):
                ya, yb = h[r], h[r + 1]
                xa, xb = v[c], v[c + 1]
                iy0 = ya + 0.25 * (yb - ya)
                iy1 = ya + 0.75 * (yb - ya)
                ix0 = xa + 0.25 * (xb - xa)
                ix1 = xa + 0.75 * (xb - xa)
                rec = by[(helper, r, c)]
                d.rectangle([(ix0 - x0) * S, (iy0 - y0) * S,
                             (ix1 - x0) * S, (iy1 - y0) * S],
                            outline=ring[rec["bin_index"] - 1], width=1)
    # legend swatch sample boxes
    for a, b in LEGEND_X:
        d.rectangle([(a - x0) * S, (LEGEND_Y[0] - y0) * S,
                     (b - x0) * S, (LEGEND_Y[1] - y0) * S],
                    outline=(255, 0, 255), width=2)
    d.text((6, 6), "green = detected panel grid   inner box = the pixels sampled "
                   "for that cell, outlined in the bin colour it matched   "
                   "magenta = legend swatches the four bin colours came from",
           fill=(0, 0, 0))
    crop.save("review/fig2c_overlay.png")


def verify(rows, geom, pal):
    out = {
        "n_cells": len(rows),
        "expected": 486,
        "legend_swatches_rgb": [[int(x) for x in p] for p in pal],
        "panel_geometry": geom,
    }
    worst = max(r["match_distance"] for r in rows)
    thin = sum(1 for r in rows if r["runner_up_margin"] < 10)
    out["colour_match"] = {
        "max_distance_to_nearest_swatch": worst,
        "cells_within_10_of_a_second_swatch": thin,
    }

    outstanding = [r for r in rows if r["bin_index"] == 4]
    acceptable = [r for r in rows if r["bin_index"] == 3]
    out["counts"] = {
        "outstanding_gt_1e6": len(outstanding),
        "acceptable_1e5_1e6": len(acceptable),
        "at_or_above_1e5": len(outstanding) + len(acceptable),
        "bin_histogram": {b[0]: sum(1 for r in rows if r["bin_index"] == b[1])
                          for b in BINS},
    }

    def marg(key):
        d = {}
        for r in outstanding:
            d[r[key]] = d.get(r[key], 0) + 1
        return d

    checks = {
        "C2_total_outstanding": (len(outstanding), REPORTED["total_outstanding"]),
        "C3_total_acceptable": (len(acceptable), REPORTED["total_acceptable"]),
    }
    out["checks_scalar"] = {k: {"digitized": a, "reported": b, "delta": a - b}
                            for k, (a, b) in checks.items()}
    out["checks_marginal"] = {}
    for name, key, rep in [("C4_fig2d_by_C6O2B2", "C6O2B2_molpct", REPORTED["by_c6o2b2"]),
                           ("C5_fig2e_by_DMGPEG", "DMGPEG_molpct", REPORTED["by_peg"]),
                           ("C6_fig2f_by_DODAP", "DODAP_molpct", REPORTED["by_dodap"]),
                           ("C7_fig2g_by_helper", "helper_lipid", REPORTED["by_helper"])]:
        got = marg(key)
        out["checks_marginal"][name] = {
            "digitized": {str(k): got.get(k, 0) for k in rep},
            "reported": {str(k): v for k, v in rep.items()},
            "delta": {str(k): got.get(k, 0) - v for k, v in rep.items()},
            "rank_order_matches": (sorted(rep, key=lambda k: -rep[k])
                                   == sorted(rep, key=lambda k: -got.get(k, 0))),
        }
    return out


def report(rows, qc):
    print(f"cells digitized: {qc['n_cells']} (expected 486)")
    print(f"bin histogram:   {qc['counts']['bin_histogram']}")
    print(f"max colour distance to nearest swatch: "
          f"{qc['colour_match']['max_distance_to_nearest_swatch']}  "
          f"(ambiguous cells: {qc['colour_match']['cells_within_10_of_a_second_swatch']})")
    print("\nscalar checks (digitized vs reported):")
    for k, v in qc["checks_scalar"].items():
        print(f"  {k:24} {v['digitized']:4}  vs {v['reported']:4}   delta {v['delta']:+d}")
    print("\nmarginal checks against the printed bar charts:")
    for k, v in qc["checks_marginal"].items():
        d = "  ".join(f"{a}:{v['digitized'][a]}/{v['reported'][a]}" for a in v["reported"])
        print(f"  {k:22} {d}   rank-order {'OK' if v['rank_order_matches'] else 'MISMATCH'}")


if __name__ == "__main__":
    main()
