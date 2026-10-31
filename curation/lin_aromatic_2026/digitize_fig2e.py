#!/usr/bin/env python3
"""Digitize Figure 2e of Lin et al., JACS 2026: the 840-lipid aromatic screen.

Library design, from the Results text and Figure 1:

    amine precursor   8A, 9A, 10A, 11A, 12A                     (5)
    aromatic analog   1B .. 7B                                  (7)
    hydrophobic tail  EP6/8/10/12/14/16/18, Ac10/12/18,
                      Ac8-1, Ac10-1                            (12)
    conjugation       nA1 = one equiv of aromatic analog,
                      nA2 = two equiv                           (2)

    5 x 7 x 12 x 2 = 840

The heatmap is laid out with 84 ROWS (aromatic x tail, "1B-EP6" .. "7B-Ac10-1")
and 10 COLUMNS (8A1, 8A2, 9A1, 9A2, 10A1, 10A2, 11A1, 11A2, 12A1, 12A2).

WHY A COLOUR READ IS EXACT HERE. The legend is DISCRETE - four bins,
<10,000 / 10,000-100,000 / 100,000-1,000,000 / >1,000,000 - so each cell is one
of exactly four printed colours and reading it is nearest-neighbour matching
against swatches sampled from the legend in the same raster. No interpolation,
no inferred scale. Same situation as the iPhos and Siegwart heatmaps, both of
which scored exactly against their own papers' marginals.

VALIDATION. Figures 2f, 2g and 2h print hit rates - the percentage of LNPs
exceeding RLU 100,000 and RLU 1,000,000 - broken down by aromatic analog (7),
by amine precursor (5) and by tail (12). That is 2 x 24 = 48 independently
printed numbers the digitization must reproduce, and they are measured from the
bars by script in digitize_fig2fgh.py rather than read by eye.

Outputs 04_fig2e_840.csv, 07_fig2e_qc.json, review/fig2e_overlay.png.
"""
import csv
import json
import os
import sys

import numpy as np
import PIL.Image as Image
import PIL.ImageDraw as ImageDraw

PAGE = sys.argv[1] if len(sys.argv) > 1 else "00_pages/p2_600.png"
REVIEW = "review"

# heatmap body, in 600 dpi page pixels (detected from the blue mask)
BODY = (3385, 733, 4627, 4142)          # x0, y0, x1, y1 (cells only: the
                                        # column-header band above y=733 has a
                                        # PINK background and must be excluded)
NROW, NCOL = 84, 10

# legend swatch sample boxes, in the same pixel space
LEGEND = [((3389, 3506), (558, 592)),   # <1e4
          ((4009, 4126), (558, 592)),   # 1e4-1e5
          ((3389, 3506), (600, 612)),   # 1e5-1e6
          ((4009, 4126), (600, 612))]   # >1e6
BINS = [("<1e4", 1, 3.16e3), ("1e4-1e5", 2, 3.16e4),
        ("1e5-1e6", 3, 3.16e5), (">1e6", 4, 3.16e6)]

AROMATICS = [f"{m}B" for m in range(1, 8)]
TAILS = ["EP6", "EP8", "EP10", "EP12", "EP14", "EP16", "EP18",
         "Ac10", "Ac12", "Ac18", "Ac8-1", "Ac10-1"]
AMINES = ["8A", "9A", "10A", "11A", "12A"]
EQUIV = ["1", "2"]

ROW_LABELS = [f"{b}-{t}" for b in AROMATICS for t in TAILS]
COL_LABELS = [f"{a}{e}" for a in AMINES for e in EQUIV]


def main():
    os.makedirs(REVIEW, exist_ok=True)
    img = np.array(Image.open(PAGE).convert("RGB")).astype(int)

    pal = np.array([np.median(img[y0:y1, x0:x1].reshape(-1, 3), axis=0)
                    for (x0, x1), (y0, y1) in LEGEND])

    x0, y0, x1, y1 = BODY
    cw, ch = (x1 - x0) / NCOL, (y1 - y0) / NROW
    rows, dists, margins = [], [], []
    for r in range(NROW):
        for c in range(NCOL):
            cy0 = int(round(y0 + ch * r + ch * 0.28))
            cy1 = int(round(y0 + ch * r + ch * 0.72))
            cx0 = int(round(x0 + cw * c + cw * 0.25))
            cx1 = int(round(x0 + cw * c + cw * 0.75))
            col = np.median(img[cy0:cy1, cx0:cx1].reshape(-1, 3), axis=0)
            d = np.linalg.norm(pal - col, axis=1)
            k = int(np.argmin(d))
            label, idx, mid = BINS[k]
            dists.append(float(d[k]))
            margins.append(float(np.sort(d)[1] - d[k]))
            arom, tail = ROW_LABELS[r].split("-", 1)
            amine_eq = COL_LABELS[c]
            amine = amine_eq[:-1]
            eq = amine_eq[-1]
            rows.append({
                "lipid_name": f"{amine}{eq}-{arom}-{tail}",
                "amine": amine, "equiv": int(eq), "aromatic": arom, "tail": tail,
                "row_label": ROW_LABELS[r], "col_label": amine_eq,
                "panel_row": r, "panel_col": c,
                "bin_label": label, "bin_index": idx, "rlu_midpoint": mid,
                "cell_rgb": f"{int(col[0])},{int(col[1])},{int(col[2])}",
                "match_distance": round(float(d[k]), 2),
                "runner_up_margin": round(float(np.sort(d)[1] - d[k]), 2),
                "assay": "firefly luciferase, IGROV1 cells",
                "readout_metric": "relative light units (RLU)",
                "value_scale": "binned",
                "source_figure": "Fig 2e",
                "extraction_method": "heatmap_color_discrete_4bin",
                "source_pdf": "JACS 2026 Lin (main)",
            })

    with open("04_fig2e_840.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    hist = {b[0]: sum(1 for r in rows if r["bin_index"] == b[1]) for b in BINS}
    qc = {
        "paper": "Lin et al., JACS 2026, 148, 6960-6969, Fig 2e",
        "n_cells": len(rows), "expected": 840,
        "grid": f"{NROW} rows (aromatic x tail) x {NCOL} cols (amine x equiv)",
        "legend_swatches_rgb": [[int(v) for v in p] for p in pal],
        "bin_histogram": hist,
        "colour_match": {"max_distance": round(max(dists), 2),
                         "cells_within_10_of_second_swatch":
                             sum(1 for m in margins if m < 10)},
        "hit_rates_digitized": hit_rates(rows),
    }
    json.dump(qc, open("07_fig2e_qc.json", "w"), indent=1)
    overlay(img, rows, pal)

    print(f"cells: {len(rows)} (expected 840)")
    print(f"bin histogram: {hist}")
    print(f"max colour distance {qc['colour_match']['max_distance']}  "
          f"ambiguous {qc['colour_match']['cells_within_10_of_second_swatch']}")
    hr = qc["hit_rates_digitized"]
    print(f"\noverall hit rate >1e5: {hr['overall']['gt_1e5']}%   "
          f">1e6: {hr['overall']['gt_1e6']}%")
    print(f"overlay -> {REVIEW}/fig2e_overlay.png")


def hit_rates(rows):
    """Percent of LNPs above each threshold, overall and by each factor.

    bin_index >= 3 means >1e5; bin_index == 4 means >1e6. These are what
    Figures 2f/2g/2h print, so they are the comparison surface.
    """
    def pct(sub, idx):
        return round(100.0 * sum(1 for r in sub if r["bin_index"] >= idx) / len(sub), 1)

    out = {"overall": {"gt_1e5": pct(rows, 3), "gt_1e6": pct(rows, 4), "n": len(rows)}}
    for name, key, order in (("by_aromatic", "aromatic", AROMATICS),
                             ("by_amine", "amine", AMINES),
                             ("by_tail", "tail", TAILS)):
        out[name] = {}
        for v in order:
            sub = [r for r in rows if r[key] == v]
            out[name][v] = {"gt_1e5": pct(sub, 3), "gt_1e6": pct(sub, 4), "n": len(sub)}
    return out


def overlay(img, rows, pal):
    x0, y0, x1, y1 = BODY
    pad = 320
    crop = Image.fromarray(img.astype("uint8")[max(0, y0 - 420):y1 + 40,
                                               max(0, x0 - pad):x1 + 40])
    d = ImageDraw.Draw(crop)
    ox, oy = max(0, x0 - pad), max(0, y0 - 420)
    ring = [(150, 150, 150), (0, 140, 255), (255, 140, 0), (200, 0, 0)]
    cw, ch = (x1 - x0) / NCOL, (y1 - y0) / NROW
    for r in rows:
        cy0 = y0 + ch * r["panel_row"] + ch * 0.28
        cy1 = y0 + ch * r["panel_row"] + ch * 0.72
        cx0 = x0 + cw * r["panel_col"] + cw * 0.25
        cx1 = x0 + cw * r["panel_col"] + cw * 0.75
        d.rectangle([cx0 - ox, cy0 - oy, cx1 - ox, cy1 - oy],
                    outline=ring[r["bin_index"] - 1], width=2)
    d.rectangle([x0 - ox, y0 - oy, x1 - ox, y1 - oy], outline=(0, 180, 0), width=4)
    for (lx0, lx1), (ly0, ly1) in LEGEND:
        d.rectangle([lx0 - ox, ly0 - oy, lx1 - ox, ly1 - oy],
                    outline=(255, 0, 255), width=4)
    d.text((8, 8), "green = detected heatmap body   inner box = pixels sampled per "
                   "cell, outlined in the bin it matched   magenta = legend swatches",
           fill=(0, 0, 0))
    crop.save(f"{REVIEW}/fig2e_overlay.png")


if __name__ == "__main__":
    main()
