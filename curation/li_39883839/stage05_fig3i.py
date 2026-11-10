#!/usr/bin/env python3
"""Digitise Figure 3I (in vivo mLuc total flux, I.M., n=3) -> 05_fig3i.json.

Tier 5. Why tier 5 is unavoidable here is argued in 00_value_sources.json.

Geometry comes only from curation/lib: the scale AND origin are fitted by
`axis.calibrate` on the six detected tick marks, and every bar top is found by
`bars.bar_top`, which requires full-width coverage plus left/right flank
agreement so a replicate marker cannot be read as a top (error class 32).

Masks are built from each bar's OWN colour, as bars.py requires: the three
series are drawn blue / red / black and their whiskers and markers share the
series colour, so a darkness mask would make the axis rule the topmost "bar".

THE EXTERNAL CHECK is in stage05_fig3h_roi.py: Fig 3H prints three ROI total
flux values straight onto the IVIS screenshot. They are single representative
mice, so each must fall inside its own bar's replicate scatter -- a number the
authors drew, never discussed, and independent of this calibration.
"""
import json, os, sys
import numpy as np
from PIL import Image

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "lib"))
import axis, bars                                            # noqa: E402
from errors import CurationLibError                          # noqa: E402

Image.MAX_IMAGE_PIXELS = None
HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "source", "fig", "pnas.2409572122fig03.jpg")

# --- panel constants, measured on this figure and recorded here, next to it --
SPINE_X = 3368              # the y-axis rule
TICK_SCAN = (3340, 3366)    # column strip OUTSIDE the rule, where ticks live
TICK_VALUES = [2.5e8, 2.0e8, 1.5e8, 1.0e8, 5.0e7, 0.0]   # top-down
PLOT_Y = (2300, 3300)       # row range that contains the whole panel
BASELINE_EXPECT = 3284

SERIES = [   # label, bar x window (full detected width), colour test
 ("H2A4B4", 3428, 3547, "blue"),
 ("A4B4-S3", 3607, 3726, "red"),
 ("A4B4-S4", 3787, 3906, "black"),
]


def masks(path):
    a = np.asarray(Image.open(path).convert("RGB")).astype(int)
    R, G, B = a[:, :, 0], a[:, :, 1], a[:, :, 2]
    return {"red": (R > 120) & (G < 95) & (B < 95),
            "blue": (B > 110) & (R < 95) & (G < 95),
            "black": (R < 110) & (G < 110) & (B < 110)}



def longest_run(row):
    """Longest contiguous True run. A significance bracket is ONE horizontal
    rule wider than a bar; a row that merely TOTALS more ink than a bar can be
    several bar edges plus the baseline, which is what a sum-based test picked
    up two pixels above the axis."""
    best = cur = 0
    for v in row:
        cur = cur + 1 if v else 0
        if cur > best:
            best = cur
    return best


def main():
    gray = np.asarray(Image.open(FIG).convert("L")).astype(float)
    mk = masks(FIG)

    ticks = axis.detect_ticks_1d(gray, "y", *TICK_SCAN,
                                 scan_lo_hi=None) if False else \
        axis.detect_ticks_1d(gray, "y", TICK_SCAN[0], TICK_SCAN[1])
    ticks = [t for t in ticks if PLOT_Y[0] <= t <= PLOT_Y[1]]
    if len(ticks) != len(TICK_VALUES):
        raise SystemExit("Fig3I: detected %d ticks, expected %d: %s"
                         % (len(ticks), len(TICK_VALUES),
                            [round(t, 1) for t in ticks]))
    cal = axis.calibrate(ticks, TICK_VALUES, kind="linear", axis="y")

    # The significance brackets are drawn in the SAME black as the A4B4-S4
    # series, span more than one bar, and sit above every bar. bars.bar_top
    # refused on them (left flank 2384 vs right flank 3208) rather than
    # returning the better-looking flank -- which is the behaviour this
    # package exists for. The scan window's top is therefore DETECTED from
    # them, not typed in: the lowest black horizontal run wider than a bar.
    BAR_W = max(x1 - x0 for _, x0, x1, _ in SERIES)
    brackets = [y for y in range(PLOT_Y[0], int(BASELINE_EXPECT) - 10)
                if longest_run(mk["black"][y, 3400:3960]) > BAR_W + 30]
    if not brackets:
        raise SystemExit("Fig3I: no significance bracket detected; the scan "
                         "top would be a literal. Refusing.")
    y_top = max(brackets) + 8

    # the baseline is the 0.0 tick itself, not the panel box.
    baseline = cal.to_pixel(0.0)
    if abs(baseline - BASELINE_EXPECT) > 3:
        raise SystemExit("Fig3I: 0.0 lands at row %.1f, expected ~%d"
                         % (baseline, BASELINE_EXPECT))

    out = {"panel": "Figure 3I", "endpoint": "whole-animal IVIS total flux (p/s)",
           "calibration": cal.qc(),
           "tick_rows_px": [round(t, 2) for t in ticks],
           "baseline_row_from_calibration": round(baseline, 2),
           "scan_top_row_px": y_top,
           "significance_bracket_rows_px": [min(brackets), max(brackets)],
           "bars": []}

    for label, x0, x1, colour in SERIES:
        try:
            top = bars.bar_top(mk[colour], x0, x1,
                               y_top=y_top, y_bot=int(baseline) - 2,
                               baseline_row=int(round(baseline)))
        except CurationLibError as e:
            out["bars"].append({"lipid": label, "REFUSED": str(e)})
            continue
        val = cal.to_value(top["row"])
        # replicate markers: the series colour inside the bar window, excluding
        # the bar's own two vertical edges and its top rule.
        inner = mk[colour][y_top:int(baseline) - 2, x0 + 12:x1 - 12]
        rows = np.nonzero(inner.sum(axis=1) > 0)[0] + y_top
        out["bars"].append({
            "lipid": label, "colour": colour,
            "bar_window_px": [x0, x1], "top_row_px": top["row"],
            "flank_delta_px": top["flank_delta_px"],
            "coverage": top["coverage"], "width_px": top["width_px"],
            "value_p_per_s": float("%.4g" % val),
            "inner_ink_row_span_px": [int(rows.min()), int(rows.max())]
            if rows.size else None,
            "inner_ink_value_span": [float("%.4g" % cal.to_value(int(rows.min()))),
                                     float("%.4g" % cal.to_value(int(rows.max())))]
            if rows.size else None})

    json.dump(out, open(os.path.join(HERE, "05_fig3i.json"), "w"), indent=1)
    print(json.dumps(out, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
