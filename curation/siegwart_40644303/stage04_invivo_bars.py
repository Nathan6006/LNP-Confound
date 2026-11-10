#!/usr/bin/env python3
"""Figures 2b and 3b: the 42 in-vivo average liver radiance values.

    python3 -I stage04_invivo_bars.py

GEOMETRY, AND WHAT IT IS ANCHORED ON.

  y scale   from the y-axis TICK MARKS only, detected as dark runs in a column
            strip strictly LEFT of the axis rule, and handed to
            curation/lib/axis.calibrate, which fits the scale AND the origin
            together from those positions. There is no way to pass an origin
            separately; that is the shape F-033 took. The baseline rule is
            NOT used as the zero: it is 1-2 px below the extrapolated 0 tick
            because the rule has width, and that 1-2 px is error class 20 in
            miniature.

  x windows from the x-axis TICK MARKS below the baseline -- 18 of them in
            Figure 2b and 24 in Figure 3b, which is the bar count each panel
            must have. The bar's own left and right edges are then measured at
            a row just above the baseline, so the flank test in
            curation/lib/bars sees the FULL bar and not a trimmed window.
            F-018: 0.72 coverage of an inner-trimmed window let three
            overlapping replicate dots be read as a bar top on 15 of 60 bars.

  ink mask  from the BAR'S OWN COLOUR, sampled at its own base. The error bars,
            the significance brackets and the axis are black and drop out. A
            darkness mask would make the axis rule itself the topmost "bar".

A BAR THAT IS NOT THERE IS NOT A ZERO. Figure 2b draws no bar at all for
3A4V-SC5 and Figure 3b draws a very short one for 5A4-SC10. The first is
recorded as drawn_height_px = 0 with value_is_bound = "upper": the true mean is
somewhere in [0, one bar-line width], and calling it 0.0 would be inventing a
measurement. bars.bar_top is never asked to produce a number for it.

THE TWO PANELS ARE NOT ON ONE SCALE. Figure 2b's axis tops out at 1e6 and
Figure 3b's at 6e7. They share no lipid and no benchmark control, so nothing in
the paper ties them together. This script reports both in absolute units and
does NOT pool them; see DECISIONS.md D-1.
"""
import json
import os
import sys

import numpy as np
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))

import axis    # noqa: E402
import bars    # noqa: E402
from errors import DetectionError  # noqa: E402

OUT = os.path.join(HERE, "04_invivo_bars.json")

TAILS6 = ["SC5", "SC6", "SC7", "SC8", "SC9", "SC10"]

PANELS = {
    "fig2b": {
        "file": "fig2.jpg",
        "cores": ["3A4", "3A4V", "3A4V2"],
        "axis_col": 239,            # the y-axis rule
        "tick_strip": (231, 237),   # strictly left of it
        "tick_values": [1e6, 8e5, 6e5, 4e5, 2e5, 0.0],
        "baseline_row": 1096,
        "xtick_row": 1099,
        "plot_top": 832,
        "x_label_drop": [214, 223],  # ink belonging to the "0" label, not ticks
    },
    "fig3b": {
        "file": "fig3.jpg",
        "cores": ["4A3", "4A5", "5A4", "5A5"],
        "axis_col": 274,
        "tick_strip": (266, 272),
        "tick_values": [6e7, 4e7, 2e7, 0.0],
        "baseline_row": 1077,
        "xtick_row": 1081,
        "plot_top": 812,
        "x_label_drop": [253, 261],
    },
}

TICK_DARK = 160
COLOUR_TOL = 60.0        # RGB distance from the bar's own base colour
# "is this pixel bar ink" at the base row, as a distance from white and a floor
# on the darkest channel. Measured fills: 3A4/4A3 (231,98,84), 3A4V/4A5
# (255,208,109), 3A4V2/5A4 (115,188,213), 5A5 (160,159,164). A saturation test
# misses the GREY 5A5 series entirely (saturation 5); a "max channel < 235"
# test misses the YELLOW series, whose red channel is 255. Black ink -- the
# axis, the error bars, the significance brackets -- has min channel ~20.
WHITE_MIN_DIST = 60.0
BLACK_MIN_CHANNEL = 55


def _runs(mask):
    out, s = [], None
    for i, v in enumerate(mask):
        if v and s is None:
            s = i
        elif not v and s is not None:
            out.append((s, i))
            s = None
    if s is not None:
        out.append((s, len(mask)))
    return out


def y_ticks(gray, spec):
    lo, hi = spec["tick_strip"]
    col = gray[spec["plot_top"]:spec["baseline_row"] + 4, lo:hi].min(axis=1)
    runs = [(a, b) for a, b in _runs(col < TICK_DARK) if 2 <= b - a <= 8]
    if len(runs) != len(spec["tick_values"]):
        raise SystemExit(
            f"stage04-bars: {len(runs)} y-tick runs found in the strip "
            f"x=[{lo},{hi}), expected {len(spec['tick_values'])}. REFUSING "
            "rather than calibrating on whatever was found: a missing tick "
            "shifts the whole scale and every value with it.")
    # darkness-weighted centroid, as axis.detect_ticks_1d does
    out = []
    for a, b in runs:
        seg = col[a:b]
        w = np.clip(TICK_DARK - seg, 0, None)
        out.append(float(spec["plot_top"] + np.dot(np.arange(a, b), w) / w.sum()))
    return out


def x_ticks(gray, spec, n_expected):
    row = gray[spec["xtick_row"], :]
    runs = [(a, b) for a, b in _runs(row < TICK_DARK) if 2 <= b - a <= 12]
    cen = [0.5 * (a + b) - 0.5 for a, b in runs]
    cen = [c for c in cen
           if all(abs(c - d) > 6 for d in spec["x_label_drop"])]
    if len(cen) != n_expected:
        raise SystemExit(
            f"stage04-bars: {len(cen)} x-tick marks found at row "
            f"{spec['xtick_row']}, expected {n_expected}. REFUSING: the x "
            "ticks are what assign a bar to a named lipid.")
    return cen


def _run_at(mask, xc, off):
    for a, b in _runs(mask):
        if a + off - 1 <= xc <= b + off:
            return a + off, b + off
    return None, None


def bar_window(arr, xc, baseline, half=28):
    """A ROUGH window for the bar, measured 4 px above the baseline.

    "Not white" rather than "saturated": Figure 3b's 5A5 series is drawn in
    GREY, whose saturation is zero, and a saturation test silently returns six
    bars that are plainly on the page. Black ink (the axis, the error bars, the
    significance brackets) is excluded by the darkness floor.

    This is only the rough window. The window bars.bar_top is finally given is
    re-derived from the SAME ink mask that bar_top scans, so its flanks are the
    bar's own edges and not a JPEG-softened column outside them -- which is
    what six refusals in Figure 3b were telling us.
    """
    y = baseline - 4
    x0 = int(round(xc - half))
    x1 = int(round(xc + half)) + 1
    strip = arr[y, max(0, x0):x1]
    not_white = np.linalg.norm(strip - 255.0, axis=1) > WHITE_MIN_DIST
    not_black = strip.min(axis=1) > BLACK_MIN_CHANNEL
    return _run_at(not_white & not_black, xc, max(0, x0))


def main():
    results = {}
    for pname, spec in PANELS.items():
        path = os.path.join(HERE, "source", "figures", spec["file"])
        arr = np.asarray(Image.open(path).convert("RGB")).astype(float)
        gray = arr.mean(axis=2)

        ticks = y_ticks(gray, spec)
        cal = axis.calibrate(ticks, spec["tick_values"], kind="linear",
                             residual_tol_px=2.5)

        names = [f"{c}-{t}" for c in spec["cores"] for t in TAILS6]
        xs = x_ticks(gray, spec, len(names))

        # The baseline rule is NOT the zero. Record the discrepancy; do not
        # use it. (Error class 20: a value axis anchored on a drawn box.)
        zero_px = cal.to_pixel(0.0)
        baseline_minus_zero = float(spec["baseline_row"] - zero_px)

        rows, failures = [], {}
        for name, xc in zip(names, xs):
            x0, x1 = bar_window(arr, xc, spec["baseline_row"])
            if x0 is None:
                rows.append({
                    "lipid": name, "x_tick_px": round(xc, 2),
                    "drawn": False, "bar_top_px": None,
                    "value": None,
                    "note": "no bar is drawn at this x tick: the panel shows "
                            "no coloured pixels above the baseline here.",
                })
                continue
            base_rgb = np.median(
                arr[spec["baseline_row"] - 5:spec["baseline_row"] - 2,
                    x0 + 3:x1 - 3].reshape(-1, 3), axis=0)
            ink = (np.linalg.norm(arr - base_rgb, axis=2) < COLOUR_TOL)
            # the window bar_top scans must be the ink's own extent
            x0, x1 = _run_at(ink[spec["baseline_row"] - 4,
                                 max(0, x0 - 6):x1 + 6], xc, max(0, x0 - 6))
            if x0 is None:
                rows.append({"lipid": name, "x_tick_px": round(xc, 2),
                             "drawn": True, "bar_top_px": None, "value": None,
                             "note": "the bar's own fill colour forms no run "
                                     "at its base row."})
                continue
            try:
                top = bars.bar_top(ink, x0, x1, spec["plot_top"],
                                   spec["baseline_row"] - 1,
                                   baseline_row=spec["baseline_row"])
            except DetectionError as e:
                failures[name] = str(e).split("\n")[0]
                rows.append({"lipid": name, "x_tick_px": round(xc, 2),
                             "drawn": True, "bar_top_px": None, "value": None,
                             "note": "bars.bar_top REFUSED: " + str(e).split("\n")[0]})
                continue
            val = float(cal.to_value(top["row"]))
            rows.append({
                "lipid": name, "x_tick_px": round(xc, 2), "drawn": True,
                "bar_x0": x0, "bar_x1": x1,
                "bar_rgb": base_rgb.round(1).tolist(),
                "bar_top_px": top["row"],
                "flank_delta_px": top["flank_delta_px"],
                "coverage": top["coverage"],
                "height_px": round(float(spec["baseline_row"] - top["row"]), 1),
                "value": round(val, 1),
            })

        results[pname] = {
            "figure": spec["file"],
            "y_calibration": cal.qc(),
            "y_tick_px": [round(t, 2) for t in ticks],
            "y_tick_values": spec["tick_values"],
            "baseline_row_px": spec["baseline_row"],
            "baseline_minus_fitted_zero_px": round(baseline_minus_zero, 2),
            "baseline_note": "positive means the drawn rule sits BELOW the "
                             "fitted zero. It is reported and not used; "
                             "anchoring on it is error class 20.",
            "n_x_ticks": len(xs),
            "x_tick_spacing_px": round(float(np.median(np.diff(xs))), 2),
            "bar_width_px": sorted({r["bar_x1"] - r["bar_x0"]
                                    for r in rows if r.get("bar_x1")}),
            "n_drawn": sum(1 for r in rows if r["drawn"]),
            "n_valued": sum(1 for r in rows if r["value"] is not None),
            "refusals": failures,
            "bars": rows,
        }

    with open(OUT, "w") as fh:
        json.dump(results, fh, indent=1)
    for p, r in results.items():
        print(f"--- {p}: {r['n_valued']}/{len(r['bars'])} valued, "
              f"tick residual {r['y_calibration']['tick_fit_residual_px']} px, "
              f"baseline-zero {r['baseline_minus_fitted_zero_px']} px")
        for b in r["bars"]:
            print(f"   {b['lipid']:<12} {str(b['value']):>12} "
                  f"{b.get('note','')}")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
