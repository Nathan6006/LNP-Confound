#!/usr/bin/env python3
"""Figure 2B -> 04_fig2b_bars.csv : total in-vivo luminescent flux, 49 LNPs.

TIER 5 (pixels). 00_value_sources.json records why tiers 1-4 do not exist for
this paper: no publisher Source Data, no values table in the SI .docx, no
repository or accession anywhere in the JATS text, and a flattened raster with
no text layer.

EVERYTHING GEOMETRIC COMES FROM curation/lib.
  * the y scale AND its origin from `axis.calibrate` on tick positions found by
    `axis.detect_ticks_1d` in the strip LEFT of the axis rule -- never from the
    panel box, a label glyph or a literal (error classes 20, 33);
  * every bar top from `bars.bar_top`, which requires >= 0.90 coverage of the
    FULL bar width over 3 contiguous rows AND left/right flank agreement, so an
    overlaid whisker or replicate marker cannot be read as a top (class 32,
    F-018, F-020).

A refusal from either is recorded as a refusal. Nothing falls through to a
literal (class 28, F-033).

THE DRAWN REFERENCE LINE. The panel carries a horizontal dashed rule labelled
"nominally high potency". The text says the authors called 1e8 the threshold
("An additional 13 LNPs were identified with total luminescence over 10^8
p/s"), so that line is a round number the authors DREW and never printed as a
coordinate. Its measured position is compared against the calibration's
prediction for 1e8 -- external to this digitiser in a way a recomputed marginal
is not (PIPELINE_V2 section 3).

Run: python3 -I stage04_fig2b.py
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
import bars                                                     # noqa: E402
from errors import CurationLibError                             # noqa: E402

FIG = os.path.join(HERE, "source", "nihms-2006250-f0002.jpg")

# --- panel constants, measured once and kept here beside the figure --------
AXIS_X = 219          # the y-axis rule (dark columns 217-221)
TICK_STRIP = (195, 216)   # columns LEFT of the rule, where the ticks live
TICK_ROWS = (1200, 1800)  # rows to search for y ticks
TICK_VALUES = [1e10, 1e9, 1e8, 1e7, 1e6, 1e5]   # top-to-bottom, as detected
BASELINE = 1744       # the x-axis rule (dark rows 1742-1746)
XTICK_STRIP = (1748, 1760)   # rows BELOW the baseline, where the x ticks live
PLOT_TOP = 1205       # above every bar, below the panel letter

# The bar colour, in RGB. Taken as the median of the saturated-blue population
# in the lower half of the figure: (102, 101, 254). The mask is built from the
# BAR'S OWN COLOUR, not from darkness -- bars.py says why: a darkness mask
# makes the black axis rule itself the topmost "bar".
BLUE_MIN_B = 150
BLUE_SEP = 50

# Transcribed left to right from the rotated x-axis labels of Figure 2B, read
# at 3x on the 300 dpi raster. The COUNT is checked against the detected tick
# count, so a dropped or doubled label fails here rather than silently shifting
# every value by one position.
LABELS = [
    "Untreated", "Naked mLuc",
    "123O10", "304O10", "205O14", "303X2", "110X7", "111O10", "122X4",
    "211O10", "303Oi10", "601O10", "503O10", "503O11", "304O13", "501O10",
    "305O10", "200O10", "174O11", "122Oi10", "205Oi10", "200O12", "347O10",
    "306O12", "509O11", "313Oi10", "509O12", "509O13", "509X3", "509O14",
    "503Oi10", "306O10", "500O13", "514O6,10", "514O8,12", "306O11",
    "516O6,10", "516O8,12", "503O8,12", "601O6,10", "509X8", "501Oi10",
    "509O4,8", "402O4,8", "200Oi10", "205O6,10", "503O6,10", "501O13",
    "500O6,10", "402O8,12", "306Oi10",
]
CONTROLS = {"Untreated", "Naked mLuc"}


def blue_mask(rgb):
    r, g, b = rgb[..., 0], rgb[..., 1], rgb[..., 2]
    return (b > BLUE_MIN_B) & (b - r > BLUE_SEP) & (b - g > BLUE_SEP)


def bar_windows(ink, centres, baseline):
    """The full left/right extent of each bar, measured, not assumed.

    Taken on the row just above the baseline, where every bar that exists is
    present. bars.bar_top explicitly forbids pre-trimming this window: an
    inner-trimmed window is what let three overlapping replicate dots be read
    as a bar top on 15 of 60 bars (F-018).
    """
    row = ink[baseline - 4]
    wins, widths = [], []
    for c in centres:
        c = int(round(c))
        if not row[c]:
            wins.append(None)
            continue
        x0 = c
        while x0 > 0 and row[x0 - 1]:
            x0 -= 1
        x1 = c
        while x1 < row.size - 1 and row[x1 + 1]:
            x1 += 1
        wins.append((x0, x1 + 1))
        widths.append(x1 + 1 - x0)
    return wins, widths


def dashed_rule_row(gray, x_lo, x_hi, y_lo, y_hi, baseline):
    """The row of the horizontal dashed 'nominally high potency' rule.

    Searched over the PANEL INTERIOR only (x_lo..x_hi is inside the axis rule,
    y_lo..y_hi stops short of the baseline), so nothing the calibration knows
    about enters this detection: the y ticks live at x < 217 and the x labels
    below the baseline, both outside the window. The bars are excluded by
    colour, not by position -- (102, 101, 254) has mean 152 and is not < 140.

    Deliberately NOT excluding rows near a y tick. An earlier version did, and
    it suppressed the only row it was looking for, because the whole point of
    this check is that the drawn rule is expected to sit ON the 1e8 tick. A
    detector that masks out the answer can only ever report a miss.

    Returns (row, coverage) or raises.
    """
    dark = gray < 140
    best, best_cov = None, 0.0
    for y in range(y_lo, y_hi):
        if abs(y - baseline) < 6:
            continue
        cov = dark[y, x_lo:x_hi].mean()
        if cov > best_cov:
            best, best_cov = y, cov
    if best is None or best_cov < 0.25:
        raise CurationLibError(
            f"dashed_rule_row: best coverage {best_cov:.3f} over "
            f"[{y_lo},{y_hi}) is below 0.25. A dashed rule spanning the panel "
            "covers far more than that; refusing rather than returning the "
            "darkest row, which would be a bar or a glyph.")
    # Sub-pixel: darkness-weighted centroid over the contiguous dark rows.
    lo = hi = best
    while lo > y_lo and dark[lo - 1, x_lo:x_hi].mean() > best_cov * 0.5:
        lo -= 1
    while hi < y_hi - 1 and dark[hi + 1, x_lo:x_hi].mean() > best_cov * 0.5:
        hi += 1
    rows = np.arange(lo, hi + 1)
    w = np.array([dark[y, x_lo:x_hi].mean() for y in rows])
    return float((rows * w).sum() / w.sum()), float(best_cov)


def main():
    img = Image.open(FIG).convert("RGB")
    rgb = np.asarray(img).astype(int)
    gray = rgb.mean(axis=2)
    ink = blue_mask(rgb)

    qc = {"figure": os.path.basename(FIG), "image_px": list(img.size)}

    # --- 1. the y axis, from its own ticks only -------------------------
    ticks = axis.detect_ticks_1d(gray[TICK_ROWS[0]:TICK_ROWS[1]], "y",
                                 TICK_STRIP[0], TICK_STRIP[1],
                                 dark=140, min_len=3, max_len=40)
    ticks = [t + TICK_ROWS[0] for t in ticks]
    if len(ticks) != len(TICK_VALUES):
        raise SystemExit(f"fig2b: {len(ticks)} y ticks, expected "
                         f"{len(TICK_VALUES)}: {np.round(ticks, 2).tolist()}")
    cal = axis.calibrate(ticks, TICK_VALUES, kind="log10", axis="y")
    qc["y_axis"] = cal.qc()
    qc["y_ticks_px"] = [round(t, 3) for t in ticks]

    # --- 2. the drawn 1e8 reference line: an EXTERNAL check -------------
    rule_row, rule_cov = dashed_rule_row(gray, 240, 2030, PLOT_TOP,
                                         BASELINE - 5, BASELINE)
    rule_value = cal.to_value(rule_row)
    qc["drawn_reference_line"] = {
        "row_px": round(rule_row, 3),
        "dark_coverage": round(rule_cov, 4),
        "value_from_calibration": float(f"{rule_value:.6g}"),
        "nominal": 1e8,
        "ratio_to_nominal": round(rule_value / 1e8, 5),
        "offset_decades": round(float(np.log10(rule_value / 1e8)), 5),
        "_why": "The authors drew this dashed rule and labelled it 'nominally "
                "high potency' without printing a coordinate; the text names "
                "10^8 p/s as the threshold. It is therefore a round number "
                "external to this digitiser, and the whole calibration is "
                "wrong by exactly the offset recorded here.",
    }

    # --- 3. the x positions, from the ticks below the baseline ----------
    g = gray.copy()
    g[:, :AXIS_X + 3] = 255
    xt = axis.detect_ticks_1d(g, "x", XTICK_STRIP[0], XTICK_STRIP[1],
                              dark=140, min_len=3, max_len=40)
    qc["n_x_ticks_detected"] = len(xt)
    qc["n_labels_transcribed"] = len(LABELS)
    if len(xt) != len(LABELS):
        raise SystemExit(
            f"fig2b: {len(xt)} x ticks but {len(LABELS)} labels transcribed. "
            "Refusing: a mismatch here shifts every value onto the wrong "
            "molecule, which is the AC_2025 failure mode.")
    gaps = np.diff(xt)
    qc["x_tick_gap_px"] = {"median": round(float(np.median(gaps)), 3),
                           "min": round(float(gaps.min()), 3),
                           "max": round(float(gaps.max()), 3)}
    if gaps.max() - gaps.min() > 1.5:
        raise SystemExit("fig2b: x ticks are not uniformly spaced; a label "
                         "may be missing.")

    # --- 4. the bars -----------------------------------------------------
    wins, widths = bar_windows(ink, xt, BASELINE)
    qc["bar_width_px"] = {"median": float(np.median(widths)),
                          "min": int(min(widths)), "max": int(max(widths))}
    results, failures, rows = {}, {}, []
    for i, (lab, w) in enumerate(zip(LABELS, wins)):
        if w is None:
            failures[lab] = ("no bar colour above the baseline at this tick "
                             "(the bar is shorter than 4 px, i.e. its value "
                             "is at the axis floor)")
            continue
        try:
            res = bars.bar_top(ink, w[0], w[1], PLOT_TOP, BASELINE - 1,
                               baseline_row=BASELINE)
        except CurationLibError as e:
            failures[lab] = str(e).split("\n")[0]
            continue
        val = cal.to_value(res["row"])
        results[lab] = res
        rows.append({
            "position": i + 1, "IL_name": lab,
            "is_control": lab in CONTROLS,
            "bar_x0": w[0], "bar_x1": w[1], "bar_width_px": w[1] - w[0],
            "top_row_px": res["row"],
            "left_row_px": res["left_row"], "right_row_px": res["right_row"],
            "flank_delta_px": res["flank_delta_px"],
            "coverage": res["coverage"],
            "total_flux_p_per_s": float(f"{val:.6g}"),
            "log10_flux": round(float(np.log10(val)), 5),
        })

    with open(os.path.join(HERE, "04_fig2b_bars.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow(r)

    qc["n_bars_read"] = len(rows)
    qc["n_refusals"] = len(failures)
    qc["refusals"] = failures
    qc["monotone_left_to_right"] = bool(
        all(a["log10_flux"] <= b["log10_flux"] + 1e-9
            for a, b in zip(rows, rows[1:])))
    with open(os.path.join(HERE, "04_fig2b_qc.json"), "w") as fh:
        json.dump(qc, fh, indent=1, sort_keys=True)
    print(json.dumps({k: qc[k] for k in
                      ("y_axis", "drawn_reference_line", "n_x_ticks_detected",
                       "x_tick_gap_px", "bar_width_px", "n_bars_read",
                       "n_refusals", "refusals", "monotone_left_to_right")},
                     indent=1, sort_keys=True))


if __name__ == "__main__":
    main()
