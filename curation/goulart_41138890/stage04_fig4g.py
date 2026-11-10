#!/usr/bin/env python3
"""Figure 4G -- the ONLY panel in this paper where the IONIZABLE LIPID varies.

Groups: Control | LNP4 (C12-200) | SM-102 (moderna) | ALC-0315 (Pfizer).
Bars per group: 24 h (white, black outline), 48 h (grey), 72 h (black).
Readout as printed: "Fluorescene Intensity / area (arb.units)", linear 0-5000.

Everything geometric comes from curation/lib (axis, bars). This file holds only
this panel's own constants and its own hazards, next to the figure.

TWO HAZARDS MEASURED IN THIS PANEL, both of which would have produced a
plausible wrong number:

  H1  A SIGNIFICANCE BRACKET crosses the plot at rows 3106-3109 spanning
      x=[1423,2035) -- 612 px, 9.1 bar widths -- at 0.97-1.00 coverage of the
      SM-102 and ALC-0315 bar windows. Scanning top-down without excluding it
      makes BOTH of those groups read ~3045 instead of ~1720/~1170. This is the
      "inset drawn INSIDE the plot box" shape from the brief's section 7, and
      the fix is the furniture strip it prescribes: a horizontal run longer than
      three bar widths is not a bar top. See seed S2.

  H2  Every bar carries an SEM cap. Measured cap coverage of the full bar
      window is 0.448-0.559 (12/12 bars), against bars.FRAC = 0.90. The caps
      are therefore excluded by the full-width rule with a margin of >= 0.34
      coverage, and separately by the flank test. That margin is MEASURED here
      and reported in 07_qc.json, not assumed.

Run:  python3 -I stage04_fig4g.py [--outdir DIR] [--seed NAME]
`--seed` mutates ONE thing and is used by stage09_seeds.py. It never writes to
the curation folder (F-069): --outdir is mandatory for a seeded run.
"""
import argparse
import hashlib
import json
import os
import sys

import numpy as np
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import axis                                            # noqa: E402
import bars                                            # noqa: E402
from errors import DetectionError                      # noqa: E402

FIG = os.path.join(HERE, "source", "figs", "gr4_lrg.jpg")
FIG_SHA256 = "0f597687010a420a8459ae0ddf99088d77d7f1c14b4fcf0e03a8edd797941aa0"

# --- panel constants, all MEASURED in this file's own detection, not typed ---
SPINE_SEARCH = (520, 760)      # column range holding the y-axis spine
PANEL_Y = (2820, 3400)         # row range holding the panel
TICK_VALUES = [5000.0, 4000.0, 3000.0, 2000.0, 1000.0, 0.0]
INK_LEVEL = 240                # background is pure white (255,255,255)
GROUPS = ["Control", "LNP4", "SM-102", "ALC-0315"]
TIMEPOINTS = ["24h", "48h", "72h"]


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def longest_run(mask_row):
    best = cur = start = best_start = 0
    for i, v in enumerate(mask_row):
        if v:
            if cur == 0:
                start = i
            cur += 1
            if cur > best:
                best, best_start = cur, start
        else:
            cur = 0
    return best, best_start


def detect_spine_and_baseline(ink):
    lo, hi = SPINE_SEARCH
    y0, y1 = PANEL_Y
    col = ink[y0:y1, lo:hi].sum(axis=0)
    spine = lo + int(np.argmax(col))
    row = ink[y0:y1, spine:spine + 1600].sum(axis=1)
    base = y0 + int(np.argmax(row))
    return spine, base


def detect_bar_windows(ink, baseline, probe_offsets=(2, 3, 4, 5, 6)):
    """Bar x-windows. A filled bar is one wide run; a WHITE bar is two narrow
    outline runs whose span is the window. Both shapes are handled here rather
    than by assuming a pitch, because an assumed pitch is a literal."""
    wins = None
    for off in probe_offsets:
        row = ink[baseline - off]
        runs, st = [], None
        for x in range(0, ink.shape[1]):
            if row[x] and st is None:
                st = x
            elif not row[x] and st is not None:
                if x - st >= 3:
                    runs.append((st, x))
                st = None
        runs = [r for r in runs if r[0] > 690]
        got = []
        i = 0
        while i < len(runs):
            a, b = runs[i]
            w = b - a
            if w >= 40:                       # a filled bar
                got.append((a, b))
                i += 1
            elif i + 1 < len(runs):           # a white bar: left + right outline
                c, d = runs[i + 1]
                if 40 <= d - a <= 90 and d - c < 40:
                    got.append((a, d))
                    i += 2
                else:
                    i += 1
            else:
                i += 1
        if len(got) == 12:
            wins = got
            break
    if wins is None:
        raise DetectionError(
            "detect_bar_windows: did not find 12 bar windows at any probe row. "
            "This is a REFUSAL. Do not fall back to an assumed pitch -- an "
            "assumed grid is error class 28, a detector miss falling through "
            "to a literal.")
    return wins


def furniture_rows(ink, bar_width, x_lo, x_hi, y_lo, y_hi, baseline):
    """Rows that are plot FURNITURE, not data: any row whose longest contiguous
    run exceeds three bar widths. The axis baseline is furniture too, but it is
    passed to bar_top separately, so it is kept out of this set."""
    limit = 3 * bar_width
    out = []
    for y in range(y_lo, y_hi):
        best, _ = longest_run(ink[y, x_lo:x_hi])
        if best >= limit and abs(y - baseline) > 4:
            out.append(y)
    return out, limit


def cap_coverage(ink, x0, x1, y_top, y_bot):
    """Max coverage reached by anything ABOVE the detected top rule -- i.e. the
    SEM cap. Reported so the margin against bars.FRAC is a measurement."""
    sub = ink[y_top:y_bot, x0:x1]
    return float(sub.mean(axis=1).max()) if sub.size else 0.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=HERE)
    ap.add_argument("--seed", default=None)
    a = ap.parse_args()
    if a.seed and os.path.abspath(a.outdir) == os.path.abspath(HERE):
        sys.exit("REFUSED: a seeded run must write to a temp --outdir (F-069).")
    os.makedirs(a.outdir, exist_ok=True)

    got = sha256(FIG)
    if got != FIG_SHA256:
        sys.exit(f"figure sha256 {got} != pinned {FIG_SHA256}")

    g = np.asarray(Image.open(FIG).convert("RGB")).astype(float).mean(axis=2)
    if a.seed == "roll_image":
        g = np.roll(g, 7, axis=0)
    ink = g < INK_LEVEL

    spine, baseline = detect_spine_and_baseline(ink)
    if a.seed == "erase_baseline":
        ink[baseline - 2:baseline + 3, :] = False
        spine, baseline = detect_spine_and_baseline(ink)

    # ---- 1. the scale, from detected ticks only -------------------------
    ticks = axis.detect_ticks_1d(g, "y", spine - 25, spine - 3,
                                 dark=160, min_len=3, max_len=40)
    ticks = [t for t in ticks if PANEL_Y[0] + 90 <= t <= PANEL_Y[1]]
    ticks = sorted(ticks)
    if a.seed == "tick_shift":
        ticks[2] += 3.0                 # gross: calibrate refuses (route A)
    if a.seed == "tick_shift_subtle":
        ticks[2] += 1.6                 # passes calibrate, must fail LOO
    if a.seed == "ticks_uniform_shift":
        ticks = [t + 4.0 for t in ticks]  # spacing untouched: calibrate and
        # LOO both pass perfectly. Only the UNPRINTED check -- the 0 tick
        # landing on the drawn baseline -- and the bar values can see it.
    if len(ticks) != 6:
        raise DetectionError(
            f"expected 6 y ticks (0..5000 by 1000), detected {len(ticks)}: "
            f"{[round(t,2) for t in ticks]}. REFUSAL, not a fallback.")
    cal = axis.calibrate(ticks, TICK_VALUES, kind="linear", axis="y")
    loo = axis.loo_residuals(ticks)

    # ---- 2. an UNPRINTED round number the authors DREW ------------------
    # The 0 tick must land on the drawn axis baseline. Nothing prints this;
    # it is external to our pipeline in a way a recomputed marginal is not.
    base_rows = [y for y in range(baseline - 4, baseline + 5)
                 if longest_run(ink[y, spine:spine + 1400])[0] > 1000]
    base_centroid = float(np.mean(base_rows)) if base_rows else float(baseline)
    zero_tick = max(ticks)
    zero_offset_px = abs(zero_tick - base_centroid)
    zero_offset_units = abs(cal.to_value(zero_tick) - cal.to_value(base_centroid))

    # ---- 3. bars --------------------------------------------------------
    wins = detect_bar_windows(ink, baseline)
    bar_w = int(np.median([b - a_ for a_, b in wins]))
    furn, furn_limit = furniture_rows(ink, bar_w, 690, 2100,
                                      PANEL_Y[0] - 60, baseline, baseline)
    work = ink.copy()
    if a.seed != "no_furniture_mask":
        for y in furn:
            work[y, :] = False

    # The scan's upper bound is the TOP DETECTED TICK, not a literal and not
    # the panel box. Above it live the figure legend (its 24h/48h/72h swatches
    # sit at x 1851-1991, directly over the ALC-0315 group -- the brief's
    # "legend swatches inside the plot box") and the first significance
    # bracket with its asterisks. Data cannot exceed the axis maximum, so this
    # bound is principled rather than convenient; a bar that DID reach it would
    # be censored, and that is asserted below rather than assumed.
    y_top = int(round(min(ticks)))
    censored = []
    for i, (x0, x1) in enumerate(wins):
        if work[y_top:y_top + 3, x0:x1].mean() >= bars.FRAC:
            censored.append(f"{GROUPS[i // 3]}_{TIMEPOINTS[i % 3]}")
    if censored:
        raise DetectionError(
            f"bars {censored} reach the top tick ({TICK_VALUES[0]}): they are "
            "censored by the axis and must carry value_is_bound = 'lower', "
            "not a read. REFUSAL.")
    results, caps, refusals = [], [], {}
    for i, (x0, x1) in enumerate(wins):
        grp, tp = GROUPS[i // 3], TIMEPOINTS[i % 3]
        try:
            if a.seed == "topmost_ink":
                sub = work[y_top:baseline, x0:x1]
                rows = np.where(sub.any(axis=1))[0]
                r = {"row": int(y_top + rows[0]), "coverage": None,
                     "left_row": None, "right_row": None,
                     "flank_delta_px": None, "width_px": x1 - x0,
                     "frac_used": None}
            else:
                r = bars.bar_top(work, x0, x1, y_top, baseline,
                                 baseline_row=baseline)
        except DetectionError as e:
            # A refusal is NOT a zero. The three Control bars are 3-14 px tall
            # and bars.bar_top refuses the two shortest as indistinguishable
            # from the axis rule; that refusal is recorded, never converted
            # into a number, and the Control group is not emitted anyway.
            refusals[f"{grp}_{tp}"] = str(e).split("\n")[0]
            results.append({"group": grp, "timepoint": tp, "x0": x0, "x1": x1,
                            "width_px": x1 - x0, "top_row": None,
                            "coverage": None, "flank_delta_px": None,
                            "max_cap_coverage_above_top": None,
                            "value": None, "refused": True})
            continue
        v = cal.to_value(r["row"])
        cap = cap_coverage(work, x0, x1, y_top, r["row"])
        caps.append(cap)
        results.append({"group": grp, "timepoint": tp, "x0": x0, "x1": x1,
                        "width_px": r["width_px"], "top_row": r["row"],
                        "coverage": r["coverage"],
                        "flank_delta_px": r["flank_delta_px"],
                        "max_cap_coverage_above_top": round(cap, 4),
                        "value": round(float(v), 1), "refused": False})

    out = {
        "figure": "Figure 4G",
        "figure_sha256": got,
        "readout_as_printed": "Fluorescene Intensity / area (arb.units)",
        "spine_x": spine, "baseline_y": baseline,
        "ticks_px": [round(t, 3) for t in ticks],
        "tick_values": TICK_VALUES,
        "tick_spacing_px": [round(b - a_, 3) for a_, b in zip(ticks, ticks[1:])],
        "tick_spacing_spread_px": round(float(np.ptp(np.diff(ticks))), 3),
        "calibration": cal.qc() if hasattr(cal, "qc") else str(cal),
        "px_per_1000_units": round(float(np.mean(np.diff(ticks))), 3),
        "loo_residuals_px": [round(x, 3) for x in loo],
        "loo_worst_px": round(float(max(loo)), 3),
        "loo_tol_px_strict": axis.LOO_TOL_PX,
        "loo_passes_strict": bool(max(loo) <= axis.LOO_TOL_PX),
        "zero_tick_vs_drawn_baseline_px": round(zero_offset_px, 3),
        "zero_tick_vs_drawn_baseline_units": round(zero_offset_units, 2),
        "bar_width_px": bar_w,
        "furniture_limit_px": furn_limit,
        "furniture_rows_masked": furn,
        "scan_y_top_is_top_tick": y_top,
        "bars_censored_at_axis_max": censored,
        "n_bars_read": sum(1 for r in results if not r["refused"]),
        "n_bars_refused": len(refusals),
        "refusals": refusals,
        "max_cap_coverage_over_read_bars": round(float(max(caps)), 4),
        "bars_frac_gate": bars.FRAC,
        "cap_margin_vs_gate": round(bars.FRAC - float(max(caps)), 4),
        "bars": results,
    }
    path = os.path.join(a.outdir, "04_fig4g.json")
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps({k: out[k] for k in
                      ("ticks_px", "tick_spacing_px", "loo_worst_px",
                       "loo_passes_strict", "zero_tick_vs_drawn_baseline_px",
                       "zero_tick_vs_drawn_baseline_units", "bar_width_px",
                       "furniture_rows_masked", "n_bars_read",
                       "n_bars_refused", "refusals",
                       "max_cap_coverage_over_read_bars",
                       "cap_margin_vs_gate")}, indent=1))
    for r in results:
        if r["refused"]:
            print(f"  {r['group']:<9} {r['timepoint']:<4} REFUSED")
            continue
        print(f"  {r['group']:<9} {r['timepoint']:<4} top_row={r['top_row']:<6} "
              f"w={r['width_px']:<3} flank={r['flank_delta_px']} "
              f"cap={r['max_cap_coverage_above_top']:.3f} "
              f"value={r['value']:>8.1f}")
    print("wrote", path)


if __name__ == "__main__":
    main()
