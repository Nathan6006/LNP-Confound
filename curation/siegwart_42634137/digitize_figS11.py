#!/usr/bin/env python3
"""Measure the Fig S11 bars (Tian et al., Adv Mater 2026, SI page 12).

Fig S11 re-plots the SAME 486-formulation screen as Fig 2c, grouped by level:
(a) C6O2B2, (b) DMG-PEG, (c) DODAP, (d) helper lipid. Each bar is the GROUP
MEAN with every member formulation overlaid as a scatter point. It therefore
adds no per-formulation data - it is a seventh consistency check on Fig 2c, and
is kept as such rather than merged into the dataset.

TWO THINGS MAKE THIS DIFFERENT FROM FIGS 2d-g.

1. The axis is LOGARITHMIC, so a bar's HEIGHT means nothing - there is no zero
   on a log axis and the baseline is wherever the axis happens to start. The
   quantity is the bar's TOP, read against the decade ticks.

2. The scatter points are drawn OVER the bars, and where a run of grey dots
   crosses a bar's top edge it cuts the blue fill into two connected
   components. Taking a component's bounding-box top then reports the lower
   fragment: DODAP 50 came back shorter than DODAP 30, reversing the order the
   figure plainly shows. Taking the MEDIAN over per-column topmost-blue rows is
   immune to both the splitting and the dot occlusion, and restores
   50 > 30 > 10.

The bottom tick's exponent is the one human-supplied number per panel (a, b, c
start at 1e3; d starts at 1e4); the overlay draws the tick that was used so the
choice is checkable. Everything else is measured.

Writes 04f_figS11_means.csv, 07f_figS11_qc.json, review/figS11_overlay.png.
"""
import csv
import json
import os
import sys

import numpy as np
import PIL.Image as Image
import PIL.ImageDraw as ImageDraw

PAGE = sys.argv[1] if len(sys.argv) > 1 else "00_pages/figS11_400.png"
REVIEW = "review"

# panel: (field, categories, bar x-spans, axis search x-range, baseline y,
#         exponent of the BOTTOM tick)
PANELS = {
    "a": ("C6O2B2_molpct", ["35", "50", "65"],
          [(1231, 1307), (1353, 1429), (1475, 1551)], (1150, 1240), 978, 3),
    "b": ("DMGPEG_molpct", ["1", "2", "4"],
          [(1876, 1953), (1998, 2075), (2120, 2197)], (1800, 1880), 976, 3),
    "c": ("DODAP_molpct", ["10", "30", "50"],
          [(1196, 1263), (1305, 1371), (1412, 1479)], (1120, 1200), 1736, 3),
    "d": ("helper_lipid", ["DOPC", "DMPC", "DOPG", "DOPE", "EggSM", "DSPC"],
          [(1781, 1828), (1858, 1904), (1934, 1980), (2010, 2057),
           (2087, 2132), (2163, 2209)], (1700, 1790), 1736, 4),
}


def bar_mask(img):
    r, g, b = img[:, :, 0], img[:, :, 1], img[:, :, 2]
    return (b > 150) & (b - r > 40) & (b - g > 25)


def ticks_for(gray, baseline, xrange_):
    axis_x = None
    for x in range(*xrange_):
        if (gray[baseline - 460:baseline + 2, x] < 170).sum() > 0.9 * 462:
            axis_x = x
            break
    if axis_x is None:
        return None, []
    ys = [y for y in range(baseline - 470, baseline + 3)
          if (gray[y, axis_x - 12:axis_x - 1] < 170).sum() >= 6]
    groups, cur = [], []
    for y in ys:
        if cur and y - cur[-1] > 2:
            groups.append(sum(cur) / len(cur))
            cur = []
        cur.append(y)
    if cur:
        groups.append(sum(cur) / len(cur))
    return axis_x, sorted(groups)


def bar_top(mask, x0, x1, ysearch):
    """Median over columns of the topmost filled row - robust to the scatter
    points that split and occlude the bar's top edge."""
    tops = []
    for x in range(x0, x1):
        col = np.nonzero(mask[ysearch[0]:ysearch[1], x])[0]
        if len(col):
            tops.append(col.min() + ysearch[0])
    return float(np.median(tops)) if tops else None, len(tops)


def main():
    os.makedirs(REVIEW, exist_ok=True)
    img = np.array(Image.open(PAGE).convert("RGB")).astype(int)
    gray = img.mean(axis=2)
    mask = bar_mask(img)

    rows, qc, draw = [], {}, []
    for pid, (field, cats, spans, axr, baseline, bot_exp) in PANELS.items():
        axis_x, tk = ticks_for(gray, baseline, axr)
        if len(tk) < 2:
            qc[pid] = {"error": f"only {len(tk)} ticks found"}
            continue
        diffs = np.diff(tk)
        px_dec = float(np.median(diffs))
        bottom_tick = tk[-1]
        vals = []
        for cat, (x0, x1) in zip(cats, spans):
            top, ncol = bar_top(mask, x0, x1, (baseline - 460, baseline + 2))
            val = 10 ** (bot_exp + (bottom_tick - top) / px_dec)
            vals.append(val)
            rows.append({
                "panel": f"Fig S11{pid}", "field": field, "category": cat,
                "group_mean_RLU": float(f"{val:.4g}"),
                "bar_top_px": round(top, 1), "bottom_tick_px": round(bottom_tick, 1),
                "bottom_tick_exponent": bot_exp,
                "px_per_decade": round(px_dec, 2),
                "columns_measured": ncol,
                "readout_metric": "luminescence (a.u.), group mean",
                "n_formulations_in_group": 486 // len(cats),
                "source_figure": f"Fig S11{pid}",
                "extraction_method": "log_axis_bar_top_median_column",
                "source_pdf": "Adv Mater 2026 Tian (SI)",
            })
            draw.append((x0, x1, top, baseline, f"{val:.2g}"))
        qc[pid] = {
            "field": field, "categories": cats,
            "axis_x": int(axis_x), "n_ticks": len(tk),
            "px_per_decade": round(px_dec, 2),
            "tick_spacing_cv": round(float(diffs.std() / diffs.mean()), 4),
            "bottom_tick_exponent": bot_exp,
            "group_means_RLU": [float(f"{v:.4g}") for v in vals],
            "rank_order": [c for _, c in sorted(zip(vals, cats), reverse=True)],
        }

    with open("04f_figS11_means.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    cross = crosscheck(qc)
    json.dump({"paper": "Tian et al., Adv Mater 2026, Fig S11",
               "what": "group means of the same 486 formulations as Fig 2c; "
                       "a consistency check, not new per-formulation data",
               "panels": qc, "crosscheck_vs_fig2c": cross,
               "review_image": f"{REVIEW}/figS11_overlay.png"},
              open("07f_figS11_qc.json", "w"), indent=1)

    overlay(img, draw)
    for pid, d in qc.items():
        if "error" in d:
            print(f"  S11{pid}: {d['error']}")
            continue
        gm = "  ".join(f"{c}:{v:.2g}" for c, v in zip(d["categories"], d["group_means_RLU"]))
        print(f"  S11{pid} ({d['field']})  {gm}")
        print(f"        px/decade {d['px_per_decade']} (tick CV {d['tick_spacing_cv']})"
              f"   rank {' > '.join(d['rank_order'])}")
    print()
    for k, v in cross.items():
        print(f"  {k}: {v}")
    print(f"\noverlay -> {REVIEW}/figS11_overlay.png")


def crosscheck(qc):
    """Rank order of S11 group means vs mean bin index per group from Fig 2c."""
    out = {}
    if not os.path.exists("04b_fig2c_486.csv"):
        return {"status": "04b not found"}
    rows = list(csv.DictReader(open("04b_fig2c_486.csv")))
    for pid, key in (("a", "C6O2B2_molpct"), ("b", "DMGPEG_molpct"),
                     ("c", "DODAP_molpct"), ("d", "helper_lipid")):
        if pid not in qc or "error" in qc[pid]:
            continue
        acc = {}
        for r in rows:
            acc.setdefault(r[key], []).append(int(r["bin_index"]))
        order2c = sorted(acc, key=lambda k: -np.mean(acc[k]))
        orderS11 = qc[pid]["rank_order"]
        out[f"S11{pid}_rank_order"] = ("MATCH" if order2c == orderS11
                                       else f"differs: fig2c {order2c} vs S11 {orderS11}")
    return out


def overlay(img, draw):
    y0, y1, x0, x1 = 600, 1800, 1100, 2300
    crop = Image.fromarray(img.astype("uint8")[y0:y1, x0:x1])
    d = ImageDraw.Draw(crop)
    for bx0, bx1, top, baseline, label in draw:
        d.line([(bx0 - x0, top - y0), (bx1 - x0, top - y0)], fill=(255, 0, 0), width=3)
        d.rectangle([bx0 - x0, top - y0, bx1 - x0, baseline - y0],
                    outline=(255, 0, 0), width=1)
        d.text((bx0 - x0, top - y0 - 14), label, fill=(0, 110, 0))
    d.text((6, 6), "red line = measured bar top (median of per-column topmost fill)   "
                   "green = implied group mean RLU", fill=(0, 0, 0))
    crop.save(f"{REVIEW}/figS11_overlay.png")


if __name__ == "__main__":
    main()
