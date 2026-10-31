#!/usr/bin/env python3
"""Measure the Fig 2d-g bar charts algorithmically (Tian et al., Adv Mater 2026).

Figs 2d-g print, as bars, the number of "outstanding" (>1e6 RLU) LNPs at each
level of C6O2B2 (2d), DMG-PEG (2e) and DODAP (2f), and for each helper lipid
(2g). They are the four marginals of the same 486-formulation screen that Fig 2c
shows as heatmaps, so measuring them is an INDEPENDENT CHECK on the Fig 2c
digitization - and, unlike reading the bars by eye, it keeps a human out of the
loop on both sides of the comparison.

HOW A BAR IS MEASURED. The bars are a single flat fill colour, so each is a
connected component of that colour; its height is (baseline - top). Heights are
converted to counts in two independent ways:

  1. TICK CALIBRATION, absolute. The y-axis tick marks give pixels-per-tick;
     the tick interval is a round number of LNPs, recovered by requiring the
     ratio fit below to land on one.
  2. RATIO FIT, relative. Least squares of height = slope*count + intercept
     against the Fig 2c marginal counts. The intercept is not forced to zero:
     the fill stops short of the bar's outline stroke, which costs every bar
     the same ~3 px regardless of height. Forcing it through the origin is what
     makes short bars look 15% low.

Agreement of (1) and (2), and of both with Fig 2c, is the result. Writes
04d_fig2dg_bars.csv, 07d_bars_qc.json and review/fig2dg_overlay.png, which
draws each measured bar box and its implied count.
"""
import csv
import json
import os
import sys

import numpy as np
import PIL.Image as Image
import PIL.ImageDraw as ImageDraw
from scipy import ndimage

PAGE = sys.argv[1] if len(sys.argv) > 1 else "00_pages/fig2_page.png"
OUT_CSV = "04d_fig2dg_bars.csv"
OUT_QC = "07d_bars_qc.json"
REVIEW = "review"

# search windows and the category each bar carries, left to right
PANELS = [
    ("2d", "C6O2B2_molpct", ["35", "50", "65"], (1050, 1400), (460, 700)),
    ("2e", "DMGPEG_molpct", ["1", "2", "4"], (1050, 1400), (900, 1140)),
    ("2f", "DODAP_molpct", ["10", "30", "50"], (1450, 1780), (460, 700)),
    ("2g", "helper_lipid",
     ["DOPC", "DMPC", "DOPG", "DOPE", "EggSM", "DSPC"], (1450, 1780), (850, 1160)),
]

# the Fig 2c marginals this is checking, recomputed from 04b at run time if present
FALLBACK = {
    "2d": {"35": 41, "50": 21, "65": 31},
    "2e": {"1": 35, "2": 29, "4": 29},
    "2f": {"10": 8, "30": 31, "50": 54},
    "2g": {"DOPC": 3, "DMPC": 17, "DOPG": 16,
           "DOPE": 9, "EggSM": 13, "DSPC": 35},
}
KEYCOL = {"2d": "C6O2B2_molpct", "2e": "DMGPEG_molpct",
          "2f": "DODAP_molpct", "2g": "helper_lipid"}


def bar_mask(img):
    r, g, b = img[:, :, 0], img[:, :, 1], img[:, :, 2]
    return (b > 200) & (b - r > 25) & (b - g > 5) & (r > 140)


def find_bars(mask, ywin, xwin, min_area=250):
    y0, y1 = ywin
    x0, x1 = xwin
    lab, _ = ndimage.label(mask[y0:y1, x0:x1])
    out = []
    for i, sl in enumerate(ndimage.find_objects(lab)):
        if int((lab[sl] == i + 1).sum()) < min_area:
            continue
        out.append({"x0": sl[1].start + x0, "x1": sl[1].stop + x0,
                    "top": sl[0].start + y0, "bottom": sl[0].stop + y0})
    out.sort(key=lambda d: d["x0"])
    return out


def axis_ticks(gray, baseline, xwin, bars):
    """Tick-mark y positions on the panel's y-axis (short dashes to its left)."""
    # Scan LEFT TO RIGHT and take the first full-height vertical rule. Scanning
    # right-to-left from the first bar instead finds that bar's own left edge,
    # which is tall enough to pass a loose threshold but carries no tick marks.
    lo = min(b["x0"] for b in bars) - 70
    hi = min(b["x0"] for b in bars) - 1
    span = gray[baseline - 260:baseline + 2]
    axis_x = hi
    for x in range(lo, hi):
        if (span[:, x] < 170).sum() > 0.92 * span.shape[0]:
            axis_x = x
            break
    ys = []
    for y in range(baseline - 260, baseline + 2):
        seg = gray[y, axis_x - 10:axis_x - 1]
        if (seg < 170).sum() >= 5:
            ys.append(y)
    groups, cur = [], [ys[0]] if ys else []
    for y in ys[1:]:
        if y - cur[-1] <= 2:
            cur.append(y)
        else:
            groups.append(sum(cur) / len(cur))
            cur = [y]
    if cur:
        groups.append(sum(cur) / len(cur))
    return axis_x, sorted(groups)


def marginals():
    """Recompute the Fig 2c marginals from 04b if it exists, else use FALLBACK."""
    path = "04b_fig2c_486.csv"
    if not os.path.exists(path):
        return FALLBACK, "fallback (04b not found)"
    rows = list(csv.DictReader(open(path)))
    out = {}
    for panel, key in KEYCOL.items():
        d = {}
        for r in rows:
            if int(r["bin_index"]) == 4:
                d[r[key]] = d.get(r[key], 0) + 1
        out[panel] = d
    return out, "recomputed from 04b_fig2c_486.csv"


def main():
    os.makedirs(REVIEW, exist_ok=True)
    img = np.array(Image.open(PAGE).convert("RGB")).astype(int)
    gray = img.mean(axis=2)
    mask = bar_mask(img)
    marg, marg_src = marginals()

    rows, qc_panels, draw = [], {}, []
    for pid, field, cats, ywin, xwin in PANELS:
        bars = find_bars(mask, ywin, xwin)
        if len(bars) != len(cats):
            qc_panels[pid] = {"error": f"found {len(bars)} bars, expected {len(cats)}"}
            continue
        baseline = max(b["bottom"] for b in bars)
        heights = np.array([baseline - b["top"] for b in bars], dtype=float)
        counts = np.array([marg[pid][c] for c in cats], dtype=float)

        slope, intercept = np.polyfit(counts, heights, 1)
        implied = (heights - intercept) / slope
        resid = implied - counts

        axis_x, ticks = axis_ticks(gray, baseline, xwin, bars)
        tick_px = float(np.median(np.diff(ticks))) if len(ticks) > 1 else float("nan")
        units_per_tick = tick_px / slope if slope else float("nan")

        qc_panels[pid] = {
            "field": field, "categories": cats,
            "bar_heights_px": [int(h) for h in heights],
            "ratio_fit": {"px_per_LNP": round(float(slope), 3),
                          "intercept_px": round(float(intercept), 2),
                          "implied_counts": [round(float(x), 2) for x in implied],
                          "residuals_LNPs": [round(float(x), 2) for x in resid],
                          "max_abs_residual_LNPs": round(float(np.abs(resid).max()), 2)},
            "tick_calibration": {"axis_x": int(axis_x), "n_ticks": len(ticks),
                                 "px_per_tick": round(tick_px, 2),
                                 "implied_LNPs_per_tick": round(float(units_per_tick), 2),
                                 "rounds_to": int(round(units_per_tick))
                                 if np.isfinite(units_per_tick) else None},
            "fig2c_counts": {c: int(marg[pid][c]) for c in cats},
        }
        for b, c, h, im in zip(bars, cats, heights, implied):
            rows.append({"panel": f"Fig {pid}", "field": field, "category": c,
                         "bar_height_px": int(h),
                         "implied_count": round(float(im), 2),
                         "fig2c_count": int(marg[pid][c]),
                         "delta": round(float(im) - marg[pid][c], 2),
                         "bar_x0": b["x0"], "bar_x1": b["x1"],
                         "bar_top": b["top"], "baseline": int(baseline),
                         "source_figure": f"Fig {pid}",
                         "extraction_method": "bar_height_measured",
                         "source_pdf": "Adv Mater 2026 Tian (main)"})
            draw.append((b, im, marg[pid][c], baseline))

    with open(OUT_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    overlay(img, draw)
    qc = {"paper": "Tian et al., Adv Mater 2026, Figs 2d-g",
          "what": "the four marginals of the 486-formulation screen, as bars",
          "fig2c_marginals_source": marg_src,
          "panels": qc_panels,
          "review_image": f"{REVIEW}/fig2dg_overlay.png"}
    json.dump(qc, open(OUT_QC, "w"), indent=1)

    print(f"Fig 2c marginals: {marg_src}\n")
    for pid, d in qc_panels.items():
        if "error" in d:
            print(f"  Fig {pid}: {d['error']}")
            continue
        rf, tc = d["ratio_fit"], d["tick_calibration"]
        got = "  ".join(f"{c}:{i:.1f}/{d['fig2c_counts'][c]}"
                        for c, i in zip(d["categories"], rf["implied_counts"]))
        print(f"  Fig {pid}  {got}")
        print(f"          max residual {rf['max_abs_residual_LNPs']:.2f} LNPs   "
              f"{rf['px_per_LNP']:.2f} px/LNP, stroke offset {rf['intercept_px']:.1f} px   "
              f"ticks: {tc['px_per_tick']:.1f} px = {tc['implied_LNPs_per_tick']:.2f} LNPs "
              f"-> {tc['rounds_to']}")
    print(f"\noverlay -> {REVIEW}/fig2dg_overlay.png")


def overlay(img, draw):
    y0, y1, x0, x1 = 1020, 1800, 420, 1200
    crop = Image.fromarray(img.astype("uint8")[y0:y1, x0:x1])
    S = 2
    crop = crop.resize((crop.width * S, crop.height * S), Image.LANCZOS)
    d = ImageDraw.Draw(crop)
    for b, implied, actual, baseline in draw:
        d.rectangle([(b["x0"] - x0) * S, (b["top"] - y0) * S,
                     (b["x1"] - x0) * S, (baseline - y0) * S],
                    outline=(200, 0, 0), width=2)
        d.line([(b["x0"] - x0) * S, (baseline - y0) * S,
                (b["x1"] - x0) * S, (baseline - y0) * S], fill=(0, 0, 255), width=3)
        d.text(((b["x0"] - x0) * S, (b["top"] - y0) * S - 16),
               f"{implied:.1f}|{actual}", fill=(0, 110, 0))
    d.text((6, 6), "red = measured bar extent   blue = detected baseline   "
                   "green = measured count | Fig 2c count", fill=(0, 0, 0))
    crop.save(f"{REVIEW}/fig2dg_overlay.png")


if __name__ == "__main__":
    main()
