#!/usr/bin/env python3
"""Figure 5a: a SECOND, INDEPENDENT reading of the same 42 radiance values.

    python3 -I stage06_scatter_check.py

WHY THIS EXISTS. Figures 2b and 3b are the only place the 42 in-vivo radiances
are drawn as bars, so a geometric error there would be invisible to any check
that re-reads the bars. Figure 5a replots exactly those 42 numbers as a scatter
against pK_a -- and we hold every pK_a EXACTLY, as text, from SI Figure S3 (SI
p3, tier 4). So each marker can be assigned to a named lipid by its x
coordinate, which we did not measure, and its y read by a calibration that
shares nothing with the bar panels: different figure, different axis, different
ticks, different scale.

CLAUDE.md: "passing every printed statistic is necessary, not sufficient --
a roughly multiplicative geometric error survives ratio checks." This check is
not a recomputed marginal. A bar-panel scale error would have to be matched by
an identical error in an unrelated panel to pass it.

WHAT IT CANNOT DO. Markers overlap. A connected component whose area is a
multiple of the single-marker area is split by k-means into that many centres,
which recovers the positions but not reliably to sub-pixel accuracy; those
lipids are reported with `merged: true` and are not treated as clean
confirmations. Markers sitting ON the x axis merge with the axis rule itself
and are lost -- which is exactly what happens to 3A4V-SC5, the one bar Figure
2b does not draw. That is reported, not patched.

x IS NOT FITTED TO ANYTHING IT VALIDATES. The x calibration comes from the
pK_a axis tick marks; the pK_a values from the SI table are the ASSIGNMENT
target, and axis.calibrate is given forbid_target so that using them as the fit
input would be a type error (error class 27).
"""
import csv
import json
import os
import sys

import numpy as np
from PIL import Image
from scipy import ndimage
from scipy.optimize import linear_sum_assignment

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))

import axis   # noqa: E402

FIG = os.path.join(HERE, "source", "figures", "fig5.jpg")
PKA = os.path.join(HERE, "02_si_table_S3.csv")
BARS = os.path.join(HERE, "04_invivo_bars.json")
OUT = os.path.join(HERE, "06_scatter_check.json")

# Panel constants, read off this figure and kept next to it. Each panel's
# y-tick strip is the 6 px immediately LEFT of its axis rule: the label glyphs
# of "2.0x10^7" reach to 1 px further left than that, and including them turns
# every tick into a pair and the uniform-spacing assertion into noise.
PANELS = [
    {"name": "fig5a_3A4",
     "ytick_strip": (151, 157), "ytick_values": [8e5, 6e5, 4e5, 2e5, 0.0],
     "xtick_px": [157.5, 231.5, 306.0, 379.5, 454.5],
     "xtick_values": [5.5, 6.0, 6.5, 7.0, 7.5],
     "plot": (159, 458, 12, 215),
     "legend_box": (335, 470, 12, 100),
     "series": {"black": "3A4", "red": "3A4V", "blue": "3A4V2"}},
    {"name": "fig5a_4A3",
     "ytick_strip": (684, 690), "ytick_values": [5e7, 4e7, 3e7, 2e7, 1e7, 0.0],
     "xtick_px": [690.5, 764.5, 838.5, 913.0, 987.0],
     "xtick_values": [5.5, 6.0, 6.5, 7.0, 7.5],
     "plot": (692, 990, 12, 215),
     "legend_box": (700, 810, 12, 80),
     "series": {"black": "4A3", "blue": "4A5"}},
    {"name": "fig5a_5A4",
     "ytick_strip": (1177, 1183),
     "ytick_values": [2e7, 1.5e7, 1e7, 5e6, 0.0],
     "xtick_px": [1184.0, 1283.0, 1382.5, 1481.5],
     "xtick_values": [6.0, 6.5, 7.0, 7.5],
     "plot": (1186, 1495, 12, 215),
     "legend_box": (1400, 1500, 12, 80),
     "series": {"black": "5A4", "red": "5A5"}},
]

MIN_COMPONENT_PX = 22
# How close a marker's x must land to the SI table's pK_a for the pairing to
# be treated as established. 0.03 pK_a units is ~4.5 px on these panels.
PKA_TOL = 0.03
TICK_DARK = 170


def masks(arr):
    R, G, B = arr[:, :, 0], arr[:, :, 1], arr[:, :, 2]
    return {
        "black": (R < 90) & (G < 90) & (B < 90),
        "red": (R > 150) & (G < 110) & (B < 110),
        "blue": (B > 150) & (R < 110) & (G < 110),
    }


def centroids(mask, plot, legend, expect):
    """Marker centres, splitting merged components by area.

    A component whose area is k times the median single-marker area is k
    overlapping markers; k-means on its pixels recovers k centres. The
    alternative -- reporting one centre for two markers -- would silently drop
    a lipid and make the panel look like it confirms fewer points than it does.
    """
    x0, x1, y0, y1 = plot
    lx0, lx1, ly0, ly1 = legend
    sub = mask[y0:y1, x0:x1].copy()
    sub[max(0, ly0 - y0):ly1 - y0, max(0, lx0 - x0):lx1 - x0] = False
    lab, n = ndimage.label(sub)
    comps = []
    for i in range(1, n + 1):
        ys, xs = np.where(lab == i)
        if len(ys) < MIN_COMPONENT_PX:
            continue
        comps.append((xs, ys))
    if not comps:
        return [], None
    areas = np.array([len(c[0]) for c in comps])
    unit = float(np.median(areas))
    out = []
    for xs, ys in comps:
        k = max(1, int(round(len(xs) / unit)))
        if k == 1:
            out.append((float(xs.mean()) + x0, float(ys.mean()) + y0, False))
            continue
        pts = np.column_stack([xs, ys]).astype(float)
        cen = pts[np.linspace(0, len(pts) - 1, k).astype(int)]
        for _ in range(40):
            d = np.linalg.norm(pts[:, None, :] - cen[None, :, :], axis=2)
            a = d.argmin(axis=1)
            for j in range(k):
                if (a == j).any():
                    cen[j] = pts[a == j].mean(axis=0)
        for cx, cy in cen:
            out.append((float(cx) + x0, float(cy) + y0, True))
    # cap at `expect`: keep the largest-area-derived ones in x order
    out.sort(key=lambda t: t[0])
    return out, unit


def main():
    arr = np.asarray(Image.open(FIG).convert("RGB")).astype(float)
    gray = arr.mean(axis=2)
    mk = masks(arr)

    with open(PKA) as fh:
        si = {r["lipid"]: float(r["pKa"]) for r in csv.DictReader(fh)}
    with open(BARS) as fh:
        barjson = json.load(fh)
    barval = {}
    for panel in barjson.values():
        for b in panel["bars"]:
            barval[b["lipid"]] = b["value"]

    report = {"_what": "Figure 5a re-read, matched to lipids by their SI-table "
                       "pK_a, and compared with the Figure 2b/3b bar values.",
              "panels": {}}

    for spec in PANELS:
        lo, hi = spec["ytick_strip"]
        t = axis.detect_ticks_1d(gray[8:222, :], "y", lo, hi, dark=TICK_DARK,
                                 min_len=2, max_len=10)
        ticks = [x + 8 for x in t]
        if len(ticks) != len(spec["ytick_values"]):
            raise SystemExit(
                f"{spec['name']}: {len(ticks)} y ticks, expected "
                f"{len(spec['ytick_values'])}. REFUSING.")
        caly = axis.calibrate(ticks, spec["ytick_values"], kind="linear")

        tails = ["SC5", "SC6", "SC7", "SC8", "SC9", "SC10"]
        wanted_pka = [si[f"{c}-{t}"]
                      for c in spec["series"].values() for t in tails]
        calx = axis.calibrate(spec["xtick_px"], spec["xtick_values"],
                              kind="linear", axis="x",
                              forbid_target=wanted_pka)

        pres = {"y_calibration": caly.qc(), "x_calibration": calx.qc(),
                "y_tick_px": [round(x, 2) for x in ticks], "series": {}}

        for colour, core in spec["series"].items():
            cents, unit = centroids(mk[colour], spec["plot"],
                                    spec["legend_box"], 6)
            names = [f"{core}-{t}" for t in tails]
            pks = np.array([si[n] for n in names])
            # OPTIMAL assignment, not greedy. Greedy here mis-assigned the
            # whole 3A4V series: the lipid whose marker is LOST (it sits on
            # the x axis and merges with the axis rule) absorbed its
            # neighbour's marker and every later lipid shifted by one, giving
            # pK_a residuals up to 0.31 and ratios from 0.47 to 1.85. The
            # residual is what exposed it, which is why it is reported per row
            # and why only assignments inside PKA_TOL are called confirmed.
            if cents:
                xs = np.array([float(calx.to_value(c[0])) for c in cents])
                cost = np.abs(pks[:, None] - xs[None, :])
                ri, ci = linear_sum_assignment(cost)
                pair = dict(zip(ri.tolist(), ci.tolist()))
            else:
                pair = {}
            rows = []
            for k, name in enumerate(names):
                if k not in pair:
                    rows.append({"lipid": name, "pKa": float(pks[k]),
                                 "found": False, "confirmed": False,
                                 "why": "fewer markers than lipids in this "
                                        "series: markers overlap, or one "
                                        "sits on the axis rule and merges "
                                        "with it."})
                    continue
                cx, cy, merged = cents[pair[k]]
                resid = float(calx.to_value(cx)) - float(pks[k])
                scat = float(caly.to_value(cy))
                bar = barval.get(name)
                ok = abs(resid) <= PKA_TOL
                rows.append({
                    "lipid": name, "pKa": float(pks[k]), "found": True,
                    "x_as_pKa": round(float(calx.to_value(cx)), 3),
                    "pKa_residual": round(resid, 3),
                    "scatter_value": round(scat, 1),
                    "bar_value": bar,
                    "ratio_scatter_over_bar": (None if not bar
                                               else round(scat / bar, 4)),
                    "merged_marker": merged,
                    "confirmed": bool(ok and not merged),
                    "why": "" if ok else
                           f"|pK_a residual| {abs(resid):.3f} > {PKA_TOL}: "
                           "this marker is not the one the SI table says sits "
                           "at this pK_a, so the pairing is not established "
                           "and the value comparison is meaningless.",
                })
            pres["series"][core] = {
                "n_markers_found": len(cents),
                "single_marker_area_px": None if unit is None else round(unit, 1),
                "n_markers_unassigned": max(0, len(cents) - len(tails)),
                "rows": rows,
            }
        report["panels"][spec["name"]] = pres

    # --- the headline comparison ---------------------------------------
    rat, big, unconf = [], [], []
    for p in report["panels"].values():
        for s in p["series"].values():
            for r in s["rows"]:
                if r.get("confirmed") and r.get("ratio_scatter_over_bar"):
                    rat.append(r["ratio_scatter_over_bar"])
                    if abs(r["ratio_scatter_over_bar"] - 1) > 0.10:
                        big.append((r["lipid"], r["ratio_scatter_over_bar"]))
                elif not r.get("confirmed"):
                    unconf.append((r["lipid"], r.get("why", "")))
    report["summary"] = {
        "n_compared": len(rat),
        "median_scatter_over_bar": round(float(np.median(rat)), 4),
        "mean_scatter_over_bar": round(float(np.mean(rat)), 4),
        "sd_scatter_over_bar": round(float(np.std(rat)), 4),
        "max_abs_pKa_residual_confirmed": round(max(
            [abs(r["pKa_residual"])
             for p in report["panels"].values()
             for s in p["series"].values()
             for r in s["rows"] if r.get("confirmed")] or [0]), 3),
        "min_ratio": round(float(np.min(rat)), 4) if rat else None,
        "max_ratio": round(float(np.max(rat)), 4) if rat else None,
        "disagreements_over_10pct": big,
        "not_independently_confirmed": unconf,
        "_read": "A median ratio at 1.00 says the bar panels and the scatter "
                 "panels agree on the ABSOLUTE scale, not merely on the "
                 "ranking. A systematic multiplicative error in the bars "
                 "would show here as a median away from 1.",
    }
    with open(OUT, "w") as fh:
        json.dump(report, fh, indent=1)
    print(json.dumps(report["summary"], indent=1))
    for pn, p in report["panels"].items():
        for core, s in p["series"].items():
            print(f"--- {pn} / {core}: {s['n_markers_found']} markers, "
                  f"{s['n_markers_unassigned']} unassigned beyond the 6 lipids")
            for r in s["rows"]:
                print("   ", r)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
