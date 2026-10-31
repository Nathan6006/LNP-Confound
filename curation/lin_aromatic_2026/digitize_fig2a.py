#!/usr/bin/env python3
"""Digitize Figure 2a of Lin et al., JACS 2026: the 147-lipid initial screen.

This is a SEPARATE experiment from the 840-lipid screen in Figure 2e, and is
curated under its own Experiment_ID because it differs in three ways that
matter:

  * amines      all twelve precursors 1A-12A, where Fig 2e uses only 8A-12A.
                Seven head groups (1A-7A) appear nowhere else in the corpus.
  * tail        EP10 only, fixed. Fig 2e varies the tail twelve ways.
  * formulation lipid/DOPE/cholesterol/DMG-PEG2000 = 15/20/25/2
                (24.2/32.3/40.3/3.2 mol%) at an IL:mRNA weight ratio of 20:1,
                the INITIAL formulation. Fig 2e was run after the optimization
                in Figs 2b-2d, at 15/15/25/2 and 10:1.

Grid: 7 rows (1B-EP10 .. 7B-EP10) x 21 columns = 147. The column set is not
every amine x {1 equiv, 2 equiv}: 1A, 2A and 3A appear only as A1. The paper
explains why - "maintaining a sufficient number of reactive amine sites in the
ionizable lipid was essential to ensuring enough hydrophobic tails" - so with
two aromatic arms those small amines have too few sites left for tails. The
printed columns are taken as given rather than reconstructed.

The legend here has THREE bins (10,000-100,000 / 100,000-1,000,000 />1,000,000),
not the four of Fig 2e: no lipid in this screen fell below 10,000 RLU. Bins are
therefore indexed 2-4 on the same RLU scale as Fig 2e so the two screens remain
comparable, rather than 1-3 which would silently shift every value by a decade.
"""
import csv
import json
import os

import numpy as np
import PIL.Image as Image
import PIL.ImageDraw as ImageDraw

PAGE = "00_pages/p2_600.png"
REVIEW = "review"

# body = the CELLS only. An earlier box ran to y=1308 and sampled pure white
# on the bottom rows; the heatmap actually ends at 1072, verified by matching
# a single column against the palette and reading off the row separators.
BODY = (640, 679, 3161, 1073)
NROW, NCOL = 7, 21

# three swatches, sampled from the legend strip at y=600
# swatch interiors, inset from the borders: boxes flush to the edge picked up
# antialiasing and shifted every swatch by ~10 RGB units
LEGEND = [((1250, 1325), (568, 592)),   # 1e4-1e5
          ((1714, 1795), (568, 592)),   # 1e5-1e6
          ((2302, 2373), (568, 592))]   # >1e6
BINS = [("1e4-1e5", 2, 3.16e4), ("1e5-1e6", 3, 3.16e5), (">1e6", 4, 3.16e6)]

ROW_LABELS = [f"{m}B-EP10" for m in range(1, 8)]
COL_LABELS = ["1A1", "2A1", "3A1",
              "4A1", "4A2", "5A1", "5A2", "6A1", "6A2", "7A1", "7A2",
              "8A1", "8A2", "9A1", "9A2", "10A1", "10A2",
              "11A1", "11A2", "12A1", "12A2"]


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
            arom = ROW_LABELS[r].split("-")[0]
            ce = COL_LABELS[c]
            amine, eq = ce[:-1], ce[-1]
            rows.append({
                "lipid_name": f"{amine}{eq}-{arom}-EP10",
                "amine": amine, "equiv": int(eq), "aromatic": arom, "tail": "EP10",
                "row_label": ROW_LABELS[r], "col_label": ce,
                "panel_row": r, "panel_col": c,
                "bin_label": label, "bin_index": idx, "rlu_midpoint": mid,
                "cell_rgb": f"{int(col[0])},{int(col[1])},{int(col[2])}",
                "match_distance": round(float(d[k]), 2),
                "runner_up_margin": round(float(np.sort(d)[1] - d[k]), 2),
                "assay": "firefly luciferase, IGROV1 cells",
                "readout_metric": "relative light units (RLU)",
                "value_scale": "binned",
                "source_figure": "Fig 2a",
                "extraction_method": "heatmap_color_discrete_3bin",
                "source_pdf": "JACS 2026 Lin (main)",
            })

    with open("04b_fig2a_147.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    hist = {b[0]: sum(1 for r in rows if r["bin_index"] == b[1]) for b in BINS}
    # overlap with the 840 screen: same lipid, independent measurement
    big = {r["lipid_name"]: r for r in csv.DictReader(open("04_fig2e_840.csv"))}
    shared = [(r["lipid_name"], r["bin_index"], int(big[r["lipid_name"]]["bin_index"]))
              for r in rows if r["lipid_name"] in big]
    agree = sum(1 for _, a, b in shared if a == b)
    qc = {
        "paper": "Lin et al., JACS 2026, Fig 2a",
        "n_cells": len(rows), "expected": 147,
        "grid": f"{NROW} rows (aromatic, EP10 tail) x {NCOL} cols",
        "legend_bins": 3,
        "bin_histogram": hist,
        "colour_match": {"max_distance": round(max(dists), 2),
                         "cells_within_10_of_second_swatch":
                             sum(1 for m in margins if m < 10)},
        "overlap_with_fig2e": {
            "n_shared_lipids": len(shared),
            "same_bin": agree,
            "note": "Fig 2a and Fig 2e are different experiments (different "
                    "formulation and IL:mRNA ratio), so disagreement is a real "
                    "biological/formulation difference, not an extraction error. "
                    "Reported as information, NOT as a pass/fail check.",
            "disagreements": [s for s in shared if s[1] != s[2]][:10],
        },
    }
    json.dump(qc, open("07b_fig2a_qc.json", "w"), indent=1)
    overlay(img, rows, pal)

    print(f"cells: {len(rows)} (expected 147)")
    print(f"bin histogram: {hist}")
    print(f"max colour distance {qc['colour_match']['max_distance']}  "
          f"ambiguous {qc['colour_match']['cells_within_10_of_second_swatch']}")
    print(f"shared lipids with Fig 2e: {len(shared)}, same bin {agree} "
          f"({100*agree/len(shared):.0f}%)" if shared else "no overlap")
    print(f"overlay -> {REVIEW}/fig2a_overlay.png")


def overlay(img, rows, pal):
    x0, y0, x1, y1 = BODY
    crop = Image.fromarray(img.astype("uint8")[y0 - 120:y1 + 30, x0 - 420:x1 + 30])
    d = ImageDraw.Draw(crop)
    ox, oy = x0 - 420, y0 - 120
    ring = [(0, 140, 255), (255, 140, 0), (200, 0, 0)]
    cw, ch = (x1 - x0) / NCOL, (y1 - y0) / NROW
    for r in rows:
        cy0 = y0 + ch * r["panel_row"] + ch * 0.28
        cy1 = y0 + ch * r["panel_row"] + ch * 0.72
        cx0 = x0 + cw * r["panel_col"] + cw * 0.25
        cx1 = x0 + cw * r["panel_col"] + cw * 0.75
        d.rectangle([cx0 - ox, cy0 - oy, cx1 - ox, cy1 - oy],
                    outline=ring[r["bin_index"] - 2], width=2)
    d.rectangle([x0 - ox, y0 - oy, x1 - ox, y1 - oy], outline=(0, 180, 0), width=4)
    for (lx0, lx1), (ly0, ly1) in LEGEND:
        d.rectangle([lx0 - ox, ly0 - oy, lx1 - ox, ly1 - oy],
                    outline=(255, 0, 255), width=4)
    d.text((8, 8), "green = heatmap body   inner box = sampled pixels, outlined in "
                   "matched bin   magenta = legend swatches", fill=(0, 0, 0))
    crop.save(f"{REVIEW}/fig2a_overlay.png")


if __name__ == "__main__":
    main()
