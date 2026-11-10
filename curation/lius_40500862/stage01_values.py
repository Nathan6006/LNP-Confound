#!/usr/bin/env python3
"""Digitise the three SI panels that carry per-member values.

Tier 5 (pixels) for all three -- 00_value_sources.json records why tiers 1-4
fail. Geometry comes from curation/lib (axis, bars); nothing here re-implements
a scale fit or a bar-top detector.

Panels
  S22C  p.34  in vitro luciferase (RLU), IGROV1, 4 bars, log 1e4..1e8
  S24   p.36  in vivo average radiance, TWO panels of 2 bars, log 1e4..1e8
  S5    p.17  spleen/liver RATIO, 11 points, log 1e0..1e3   (NOT a delivery
              magnitude -- digitised for the record, see DECISIONS.md D-1)

Every calibration is fitted from DETECTED TICK POSITIONS ONLY and is checked
with axis.loo_residuals as well as the max-residual test, because the
max-residual test is weak against a single bad tick.
"""
import json
import os
import sys

import numpy as np
from PIL import Image

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__),
                                                "..", "lib")))
import axis                                                     # noqa: E402
import bars                                                     # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))


def rgb(name):
    return np.asarray(Image.open(os.path.join(HERE, name)).convert("RGB"),
                      dtype=float)


def gray(a):
    return a.mean(2)


def colour_mask(a, lo_sat=18):
    """True where a pixel is CHROMATIC -- a bar's own colour.

    Not a darkness mask: bars.bar_top's docstring is explicit that building the
    mask from darkness makes the black axis rule the topmost 'bar'. Axis rules,
    tick labels and error whiskers here are all neutral grey/black, every bar
    is a pastel, so chroma separates them cleanly.

    SATURATION ONLY, with no ceiling on the brightest channel. The first
    version of this function also required max(R,G,B) <= 252 to exclude white,
    and that silently deleted a whole bar: S22C's third bar is (254, 224, 224),
    saturation 30, so it is unambiguously chromatic while its red channel is
    254. The run then found 3 bar windows where 4 were expected and the script
    stopped -- which is the only reason the deletion was visible. White has
    saturation 0 and the black axis rule's JPEG fringe measures 6, both far
    below the threshold, so the ceiling was doing no work except that harm.
    """
    mx, mn = a.max(2), a.min(2)
    return (mx - mn) >= lo_sat


def runs(flag, min_len):
    """Contiguous True runs of at least min_len, as [start, stop) pairs."""
    out, s = [], None
    for i, v in enumerate(flag):
        if v and s is None:
            s = i
        elif not v and s is not None:
            if i - s >= min_len:
                out.append((s, i))
            s = None
    if s is not None and len(flag) - s >= min_len:
        out.append((s, len(flag)))
    return out


def calibrate_decades(g, axis_col, top_value, n_dec, label, scan_w=10,
                      y_max=None):
    """Fit a log10 y calibration from the decade ticks left of `axis_col`."""
    cand = axis.detect_ticks_1d(g, "y", axis_col - scan_w, axis_col,
                                dark=140, min_len=1, max_len=40)
    if y_max is not None:
        cand = [c for c in cand if c <= y_max]
    cand = sorted(cand)[:n_dec]
    if len(cand) != n_dec:
        raise SystemExit(f"{label}: detected {len(cand)} ticks, need {n_dec}")
    values = [top_value / (10.0 ** i) for i in range(n_dec)]
    cal = axis.calibrate(cand, values, kind="log10", axis="y",
                         loo_tol_px=axis.LOO_TOL_PX)
    loo = axis.loo_residuals(cand)
    return cal, {"ticks_px": [round(c, 3) for c in cand],
                 "tick_values": values,
                 "tick_spacing_px": [round(b - a, 3)
                                     for a, b in zip(cand, cand[1:])],
                 "loo_residuals_px": [round(x, 4) for x in loo],
                 "loo_worst_px": round(max(loo), 4),
                 "loo_tol_px": axis.LOO_TOL_PX,
                 **cal.qc()}


# --------------------------------------------------------------------------
# S22C -- in vitro luciferase, IGROV1
# --------------------------------------------------------------------------
def do_s22c(qc):
    a = rgb("figs/p34_x281.jpeg")
    g = gray(a)
    # BASE is panel C's OWN axis rule, detected over panel C's x range only.
    # It was 440 in the first version of this file, taken from a dark-row scan
    # of the FULL page image -- row 440 is panel A's flow-cytometry box, 60 px
    # below panel C's rule. That is error class 20, an origin anchored on a
    # panel box instead of on detected ticks. The bar values did not move,
    # because every bar top (rows 197-281) sits far above both candidates, but
    # the no-height guard was aimed at the wrong line and 60 px of axis-label
    # area was inside the scan window.
    AXIS_COL, BASE = 1119, 380
    cal, ev = calibrate_decades(g, AXIS_COL, 1e8, 5, "S22C", y_max=430)
    qc["S22C"] = {"image": "figs/p34_x281.jpeg", "image_px": list(a.shape[1::-1]),
                  "axis_col": AXIS_COL, "baseline_row": BASE, "calibration": ev}

    chroma = colour_mask(a)
    plot = np.zeros_like(chroma)
    plot[90:BASE - 1, AXIS_COL + 2:1420] = chroma[90:BASE - 1,
                                                 AXIS_COL + 2:1420]
    cols = plot.sum(0) > 4
    wins = runs(cols, 12)
    qc["S22C"]["bar_windows"] = [[int(x0), int(x1)] for x0, x1 in wins]
    if len(wins) != 4:
        raise SystemExit(f"S22C: found {len(wins)} bar windows, expected 4: {wins}")

    names = ["A5-CE-C10", "A5-C10", "A5-CE-C7-6", "A5-C7-6"]
    out = []
    for nm, (x0, x1) in zip(names, wins):
        r = bars.bar_top(plot, x0, x1, 90, BASE - 1, baseline_row=BASE)
        val = cal.to_value(r["row"])
        out.append({"IL_name": nm, "value": val, **r})
    qc["S22C"]["bars"] = [{k: (round(v, 5) if isinstance(v, float) else v)
                           for k, v in b.items()} for b in out]
    return out


# --------------------------------------------------------------------------
# S24 -- in vivo average radiance, two independent panels
# --------------------------------------------------------------------------
def do_s24(qc):
    a = rgb("figs/p36_x285.jpeg")
    g = gray(a)
    panels = [("S24_left", 82, 350, 100, 360, ["A5-CE-C7-6", "A5-C7-6"]),
              ("S24_right", 424, 352, 430, 627, ["A5-CE-C10", "A5-C10"])]
    chroma = colour_mask(a)
    res = {}
    for label, acol, base, xlo, xhi, names in panels:
        cal, ev = calibrate_decades(g, acol, 1e8, 5, label, y_max=base + 5)
        plot = np.zeros_like(chroma)
        plot[40:base - 1, acol + 2:xhi] = chroma[40:base - 1, acol + 2:xhi]
        wins = runs(plot.sum(0) > 4, 10)
        if len(wins) != 2:
            raise SystemExit(f"{label}: {len(wins)} bar windows, expected 2: {wins}")
        got = []
        for nm, (x0, x1) in zip(names, wins):
            r = bars.bar_top(plot, x0, x1, 40, base - 1, baseline_row=base)
            got.append({"IL_name": nm, "value": cal.to_value(r["row"]), **r})
        qc[label] = {"image": "figs/p36_x285.jpeg", "axis_col": acol,
                     "baseline_row": base, "calibration": ev,
                     "bar_windows": [[int(x0), int(x1)] for x0, x1 in wins],
                     "bars": [{k: (round(v, 5) if isinstance(v, float) else v)
                               for k, v in b.items()} for b in got]}
        res[label] = got
    return res


# --------------------------------------------------------------------------
# S5 -- spleen/liver ratio, 11 scatter points
# --------------------------------------------------------------------------
def do_s5(qc):
    a = rgb("figs/p17_x247.jpeg")
    g = gray(a)
    AXIS_COL, BASE = 162, 435
    cal, ev = calibrate_decades(g, AXIS_COL, 1e3, 4, "S5", y_max=440)

    # The markers are a single pastel violet diamond each; no bars, so bar_top
    # does not apply. A diamond's CENTROID is its value -- a top-edge rule
    # would read the marker's half-height as signal.
    chroma = colour_mask(a)
    plot = np.zeros_like(chroma)
    plot[20:BASE - 2, AXIS_COL + 2:880] = chroma[20:BASE - 2, AXIS_COL + 2:880]
    wins = runs(plot.sum(0) > 2, 8)
    if len(wins) != 11:
        raise SystemExit(f"S5: {len(wins)} marker windows, expected 11: {wins}")

    names = ["A5-CE-C6", "A5-CE-C8", "A5-CE-C9", "A5-CE-C10", "A5-CE-C5-4",
             "A5-CE-C6-5", "A5-CE-C7-6", "A5-CE-C7-7", "A5-CE-C8-7",
             "A5-CE-C8-8", "A5-CE-C10-6"]
    out = []
    for nm, (x0, x1) in zip(names, wins):
        sub = plot[:, x0:x1]
        ys, _ = np.nonzero(sub)
        cy = float(ys.mean())
        # symmetry check: a diamond's centroid sits at its bounding-box centre
        mid = 0.5 * (ys.min() + ys.max())
        out.append({"IL_name": nm, "value": cal.to_value(cy),
                    "centroid_row": round(cy, 3),
                    "bbox_mid_row": round(mid, 3),
                    "centroid_minus_bboxmid_px": round(cy - mid, 3),
                    "marker_h_px": int(ys.max() - ys.min() + 1),
                    "marker_w_px": int(x1 - x0)})
    worst = max(abs(o["centroid_minus_bboxmid_px"]) for o in out)
    if worst > 1.5:
        raise SystemExit(f"S5: marker centroid vs bbox-centre disagrees by "
                         f"{worst} px -- the mask is catching something that "
                         "is not a symmetric diamond.")
    qc["S5"] = {"image": "figs/p17_x247.jpeg", "axis_col": AXIS_COL,
                "baseline_row": BASE, "calibration": ev,
                "marker_windows": [[int(x0), int(x1)] for x0, x1 in wins],
                "centroid_vs_bbox_worst_px": round(worst, 3),
                "points": [{k: (round(v, 5) if isinstance(v, float) else v)
                            for k, v in p.items()} for p in out]}
    return out


# --------------------------------------------------------------------------
# S23A -- spleen/liver ratio AGAIN, on a LINEAR axis, overlapping S5's lipids
# --------------------------------------------------------------------------
def do_s23a(qc):
    """The trap-4 panel. S23A reports the SAME quantity as S5 -- spleen/liver
    average radiance efficiency -- for two lipids S5 also reports. Digitised
    so the disagreement is a measurement rather than an impression."""
    a = rgb("figs/p35_x283.jpeg")
    g = gray(a)
    AXIS_COL, BASE = 157, 337
    cand = sorted(x for x in axis.detect_ticks_1d(g, "y", AXIS_COL - 8,
                                                  AXIS_COL, dark=140,
                                                  min_len=1, max_len=40)
                  if x <= BASE + 3)
    if len(cand) != 6:
        raise SystemExit(f"S23A: {len(cand)} ticks, expected 6: {cand}")
    values = [25.0, 20.0, 15.0, 10.0, 5.0, 0.0]
    cal = axis.calibrate(cand, values, kind="linear", axis="y",
                         loo_tol_px=axis.LOO_TOL_PX)
    loo = axis.loo_residuals(cand)

    chroma = colour_mask(a)
    plot = np.zeros_like(chroma)
    plot[60:BASE - 1, AXIS_COL + 2:430] = chroma[60:BASE - 1,
                                                 AXIS_COL + 2:430]
    wins = runs(plot.sum(0) > 4, 12)
    if len(wins) != 4:
        raise SystemExit(f"S23A: {len(wins)} bar windows, expected 4: {wins}")
    names = ["A5-C10", "A5-CE-C10", "A5-C7-6", "A5-CE-C7-6"]
    out = []
    for nm, (x0, x1) in zip(names, wins):
        r = bars.bar_top(plot, x0, x1, 60, BASE - 1, baseline_row=BASE)
        out.append({"IL_name": nm, "value": cal.to_value(r["row"]), **r})
    qc["S23A"] = {"image": "figs/p35_x283.jpeg", "axis_col": AXIS_COL,
                  "baseline_row": BASE,
                  "calibration": {"ticks_px": [round(c, 3) for c in cand],
                                  "tick_values": values,
                                  "tick_spacing_px": [round(b - x, 3) for x, b
                                                      in zip(cand, cand[1:])],
                                  "loo_residuals_px": [round(x, 4) for x in loo],
                                  "loo_worst_px": round(max(loo), 4),
                                  "loo_tol_px": axis.LOO_TOL_PX, **cal.qc()},
                  "bar_windows": [[int(x0), int(x1)] for x0, x1 in wins],
                  "bars": [{k: (round(v, 5) if isinstance(v, float) else v)
                            for k, v in b.items()} for b in out]}
    return out


def main():
    qc = {}
    s22 = do_s22c(qc)
    s24 = do_s24(qc)
    s5 = do_s5(qc)
    s23 = do_s23a(qc)

    # ---- trap 4: two panels, same quantity, shared lipids -> measure it ----
    by5 = {p["IL_name"]: p["value"] for p in s5}
    by23 = {b["IL_name"]: b["value"] for b in s23}
    shared = sorted(set(by5) & set(by23))
    qc["panel_disagreement_S5_vs_S23A"] = {
        "quantity": "spleen/liver average radiance efficiency ratio",
        "shared_lipids": shared,
        "per_lipid": {k: {"S5": round(by5[k], 4), "S23A": round(by23[k], 4),
                          "ratio_S23A_over_S5": round(by23[k] / by5[k], 4)}
                      for k in shared},
        "S5_orders": "A5-CE-C7-6 above A5-CE-C10",
        "S23A_orders": "A5-CE-C10 above A5-CE-C7-6",
        "rank_inverted": (by5.get("A5-CE-C7-6", 0) > by5.get("A5-CE-C10", 0)) !=
                         (by23.get("A5-CE-C7-6", 0) > by23.get("A5-CE-C10", 0)),
        "note": "NOT POOLED. Recorded in 08_contradictions.csv.",
    }
    # ---- the two S24 panels are separate axes: compare their scales --------
    l = qc["S24_left"]["calibration"]["px_per_unit"]
    r = qc["S24_right"]["calibration"]["px_per_unit"]
    qc["S24_two_panel_scale_check"] = {
        "left_px_per_decade": l, "right_px_per_decade": r,
        "relative_difference": round(abs(l - r) / ((l + r) / 2), 6),
        "left_origin_px": qc["S24_left"]["calibration"]["origin_px"],
        "right_origin_px": qc["S24_right"]["calibration"]["origin_px"],
        "note": "The two panels share NO lipid, so no per-lipid ratio between "
                "them can be measured. Their scales agree but their origins "
                "differ, so they are calibrated and reported separately and "
                "never pooled into one block.",
    }

    rows = []
    for b in s22:
        rows.append(dict(panel="S22C", page=34, IL_name=b["IL_name"],
                         value=b["value"], unit="RLU"))
    for lab, got in s24.items():
        for b in got:
            rows.append(dict(panel=lab, page=36, IL_name=b["IL_name"],
                             value=b["value"], unit="p/sec/cm2/sr"))
    for p in s5:
        rows.append(dict(panel="S5", page=17, IL_name=p["IL_name"],
                         value=p["value"], unit="ratio_spleen_over_liver"))
    for b in s23:
        rows.append(dict(panel="S23A", page=35, IL_name=b["IL_name"],
                         value=b["value"], unit="ratio_spleen_over_liver"))

    with open(os.path.join(HERE, "01_digitised.json"), "w") as fh:
        json.dump({"qc": qc, "rows": rows}, fh, indent=1)
    for r in rows:
        print(f"{r['panel']:<11} {r['IL_name']:<12} {r['value']:>12.4g} {r['unit']}")
    print(f"\n{len(rows)} digitised points written to 01_digitised.json")


if __name__ == "__main__":
    main()
