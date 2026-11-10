#!/usr/bin/env python3
"""Measure the Fig 2f/2g/2h hit-rate bars and check them against Fig 2e.

Figures 2f, 2g and 2h print the percentage of the 840 LNPs exceeding
RLU 100,000 (pale bars) and RLU 1,000,000 (dark bars), broken down by aromatic
analog (7), amine precursor (5) and hydrophobic tail (12). That is 48
independently printed numbers, and they are the marginals of the same heatmap
digitized in digitize_fig2e.py - so reproducing them validates the heatmap read
without any human sitting on both sides of the comparison.

Each bar's height is converted to a percentage by TICK CALIBRATION, and from
2026-10-30 that is actually true. It was not before, which is review finding
F-006 and error class 27.

What the code used to do: `slope`, the pixels-per-percent that converts EVERY
bar, was a through-origin least-squares fit of the measured bar heights onto
`rep_l` -- the Fig 2e hit rates this panel exists to CHECK.

    slope = sum(heights * rep) / sum(rep * rep)

The tick spacing was measured, but only to print a diagnostic, and the file
contained no `assert`. So the calibration was fitted to the answer: scaling
every digitized rate by a constant k rescales `slope` by 1/k and the measured
percentages by k, leaving the reported deltas unchanged. **A uniform
multiplicative error in the heatmap read was absorbed and invisible** -- exactly
the failure mode of class 20, which this check was meant to close. The agreement
was never worthless (a row or column permutation would still show up as a shape
mismatch over 7, 5 and 12 points), but it could not catch a scale error, and the
project cited it as if it could.

The giveaway was in its own QC output: `px_per_tick`/20 gave 5.100 / 5.150 /
5.300 while `px_per_percent` was 5.046 / 5.071 / 5.293. Had the conversion been
tick-anchored those columns would agree. They did not.

Now: a least-squares line through the detected tick centres gives pixels per
tick, the tick interval is asserted against TICK_PCT through the drawn axis
height (F-053 -- the docstring claimed this check for twelve days while no
assert mentioned TICK_PCT at all, and TICK_PCT multiplies all 48 numbers), and
`rep_l` appears
nowhere in the scale. The agreement with Fig 2e is then asserted, so a
disagreement fails the build instead of being printed.

Writes 05_fig2fgh_bars.csv, 08_fig2fgh_qc.json, review/fig2fgh_overlay.png.
"""
import csv
import json
import os

import numpy as np
import PIL.Image as Image
import PIL.ImageDraw as ImageDraw
from scipy import ndimage

PAGE = "00_pages/p2_600.png"
REVIEW = "review"

AROMATICS = [f"{m}B" for m in range(1, 8)]
AMINES = ["8A", "9A", "10A", "11A", "12A"]
TAILS = ["EP6", "EP8", "EP10", "EP12", "EP14", "EP16", "EP18",
         "Ac10", "Ac12", "Ac18", "Ac8-1", "Ac10-1"]

# panel: (field, categories, baseline y, x window)
PANELS = [
    ("2f", "aromatic", AROMATICS, 3335, (560, 1790)),
    ("2g", "amine", AMINES, 3335, (1800, 2800)),
    ("2h", "tail", TAILS, 4055, (560, 2800)),
]


def masks(img):
    r, g, b = img[:, :, 0], img[:, :, 1], img[:, :, 2]
    light = (r > 225) & (r - b > 15) & (r - g > 15) & (g > 185)
    dark = (r > 190) & (r < 235) & (r - b > 35) & (r - g > 35) & (g <= 185)
    return light, dark


# The y-axis tick interval, in percentage points. Asserted below against the
# detected tick spacing rather than recovered from a fit -- "recover the interval
# by requiring the fit to land on a round number" is how the fitted scale got in.
TICK_PCT = 20.0

# Tolerance on the Fig 2f-h vs Fig 2e agreement, in percentage points.
#
# CORRECTED 2026-11-07 (round-2 finding E2). This read 1.1, against a measured
# max of 1.00 pp and a mean of -0.52 pp with 47 of 48 deltas negative, and the
# comment attributed that one-sided bias to the BOTTOM edge -- the colour mask's
# last saturated row sitting 1-2 px above the axis rule.
#
# That diagnosis was wrong, and the bias was removable rather than inherent. The
# dominant term was the ORIGIN: F-006 re-anchored the scale on the ticks and
# left the 0% line as three hard-coded literals. On 2f and 2g the 0% tick is not
# drawn, the ladder extrapolates it 3-4 px below the panel box, and the literal
# 3335 was short by 0.61 pp on each. Deriving the origin from the tick fit, as
# the scale already was:
#
#     mean   -0.52 pp  ->  -0.042 pp
#     max |d| 1.00 pp  ->   0.40 pp
#     negative 47/48   ->   25/48      <- the one-sided signature is gone
#
# 25 of 48 is what unbiased noise looks like, which is the real evidence the
# cause was found rather than the number massaged.
#
# 0.6 is set from the MEASUREMENT FLOOR, not from the new worst value: a
# sub-pixel top edge is good to roughly half a pixel, and half a pixel is
# ~0.1 pp here, so 0.6 pp is several times the floor and still less than a
# sixth of the old tolerance. Deliberately not 0.45 -- a tolerance pinned to
# whatever the current maximum happens to be is a tolerance that cannot fail.
#
# The old comment's closing point still stands and is why the residual is
# recorded rather than driven to zero: these 48 numbers VALIDATE the Fig 2e
# heatmap read, no corpus value is computed from them (`build_rows.py` reads
# 06_merged.csv, not 05_fig2fgh_bars.csv), and tuning a residual away is exactly
# what the old fitted scale did. What changed is that this particular residual
# turned out to be an error with an address, not noise to be tolerated.
WORST_TOL = 0.6


def bars_in(mask, baseline, xwin, min_area=3000, report_rejected=False):
    """Bars whose bottom sits on this panel's baseline.

    With report_rejected, also returns the areas of components that sit on this
    baseline and were dropped by `min_area`. That list is what makes the floor
    auditable (F-051): the floor is only safe if nothing lands just underneath
    it, and a bare threshold cannot tell you that.
    """
    x0, x1 = xwin
    lab, _ = ndimage.label(mask[:, x0:x1])
    out, rejected = [], []
    for i, sl in enumerate(ndimage.find_objects(lab)):
        area = int((lab[sl] == i + 1).sum())
        if area < min_area:
            if abs(sl[0].stop - baseline) <= 6:
                rejected.append(area)
            continue
        bot, top = sl[0].stop, sl[0].start
        if abs(bot - baseline) > 6:
            continue
        out.append({"x0": sl[1].start + x0, "x1": sl[1].stop + x0,
                    "top": top, "bottom": bot, "height": bot - top,
                    "area": area})
    out.sort(key=lambda d: d["x0"])
    return (out, sorted(rejected, reverse=True)) if report_rejected else out


def subpixel_top(img, bar, baseline):
    """The bar's top edge at the HALF-INTENSITY point of its anti-aliasing ramp.

    Why this is not a refinement. The colour masks select only fully saturated
    bar pixels, so every top edge read a few pixels low -- and on ~400 px bars a
    4 px loss is about 1%. With the scale FITTED to the Fig 2e rates (F-006),
    that bias was absorbed exactly: the fitted 5.046 px/pp against the
    tick-implied 5.100 is a 1.07% difference, and the heights were ~1.1% short.
    The two cancelled, every delta came out under 1.6 pp, and the panel was
    reported as validating the heatmap read.

    Taking the scale from the ticks made the bias visible: all 48 deltas turned
    NEGATIVE, from -0.4 to -1.6 pp, which is a systematic error and not noise.
    So the fit was not merely blind to a scale error, it was CORRECTING a real
    edge-detection error with a fudge factor derived from the answer. That is
    the sharpest statement of class 27 this project has.

    Interpolates where the green channel crosses halfway between the page
    background and the bar's own saturated colour, which is where the true edge
    lies for a box filter, and is uniform across bars of any height.
    """
    g = img[:, :, 1].astype(float)
    x0, x1 = bar["x0"], bar["x1"]
    w = max(1, (x1 - x0) // 8)
    col = np.median(g[:, x0 + w:x1 - w], axis=1)      # trim anti-aliased flanks
    top, bot = bar["top"], bar["bottom"]
    bg = float(np.median(col[max(0, top - 24):max(1, top - 6)]))
    fill = float(np.median(col[top + 2:min(len(col), top + 12)])) if bot - top > 14 \
        else float(np.min(col[top:bot + 1]))
    if not np.isfinite(bg) or not np.isfinite(fill) or bg - fill < 20:
        return float(top)                             # no usable ramp
    half = (bg + fill) / 2.0
    for y in range(top, max(0, top - 20), -1):
        if col[y] <= half <= col[y - 1]:
            lo, hi = col[y], col[y - 1]
            frac = 0.0 if hi == lo else (hi - half) / (hi - lo)
            return (y - 1) + frac
    return float(top)


def axis_ticks(gray, baseline, xwin):
    """Tick y-positions on the panel's y-axis rule."""
    lo, hi = xwin[0], xwin[0] + 420
    span = gray[baseline - 560:baseline + 3]
    axis_x = None
    for x in range(lo, hi):
        if (span[:, x] < 170).sum() > 0.9 * span.shape[0]:
            axis_x = x
            break
    if axis_x is None:
        return None, []
    ys = [y for y in range(baseline - 560, baseline + 3)
          if (gray[y, axis_x - 22:axis_x - 2] < 170).sum() >= 10]
    groups, cur = [], []
    for y in ys:
        if cur and y - cur[-1] > 3:
            groups.append(sum(cur) / len(cur))
            cur = []
        cur.append(y)
    if cur:
        groups.append(sum(cur) / len(cur))
    return axis_x, sorted(groups)


def main():
    os.makedirs(REVIEW, exist_ok=True)
    img = np.array(Image.open(PAGE).convert("RGB")).astype(int)
    gray = img.mean(axis=2)
    light, dark = masks(img)
    fig2e = {r["lipid_name"]: r for r in csv.DictReader(open("04_fig2e_840.csv"))}
    allrows = list(fig2e.values())

    out_rows, qc, draw, missed = [], {}, [], []
    for pid, field, cats, baseline, xwin in PANELS:
        # F-051. The light mask keeps the 3000 floor and the dark mask uses 120,
        # and the asymmetry is REAL rather than an unrepaired oversight -- the
        # light colour also matches axis text, so at 120 every panel picks up
        # glyphs, `len(lb) != len(cats)` fires on all three, and the run emits
        # nothing. Measured: lowering it produced 0 of 48 rows.
        #
        # But the hazard the asymmetry hides is real too: a genuine light bar
        # below 3000 px^2 would be discarded in silence. So instead of moving
        # the floor, the floor is PROVEN not to cut through the bar population
        # -- `rejected_max` vs `accepted_min` below. A floor is safe when there
        # is a wide gap on either side of it, and unsafe when something sits
        # just above it, whatever its value.
        lb, lb_rejected = bars_in(light, baseline, xwin, report_rejected=True)
        # 120, not 800: EP6's >1e6 bar is ~50x9 = 450 px^2 and an 800 px^2
        # floor silently discarded it (F-007). The floor exists to reject text
        # glyphs and axis furniture, and 120 is still well above those (the next-smallest
        # component on this baseline is 3,676 px^2).
        db = bars_in(dark, baseline, xwin, min_area=120)
        axis_x, ticks = axis_ticks(gray, baseline, xwin)
        if len(lb) != len(cats) or not ticks:
            qc[pid] = {"error": f"found {len(lb)} light bars / {len(ticks)} ticks, "
                                f"expected {len(cats)} bars"}
            continue
        # F-051: the floor is auditable, not merely set. A wide gap either side
        # of it means it is separating bars from furniture; something sitting
        # just under it means it is cutting through the bars.
        accepted_min = min(b["area"] for b in lb)
        rejected_max = max(lb_rejected) if lb_rejected else 0
        assert rejected_max * 3 < accepted_min, (
            f"panel {pid}: the largest light component DROPPED by min_area is "
            f"{rejected_max} px^2 and the smallest one KEPT is {accepted_min} "
            f"px^2. They are not cleanly separated, so the floor may be "
            f"discarding a real bar in silence (F-051, class 28). Inspect the "
            f"overlay before changing the floor -- at 120 this mask picks up "
            f"axis text and every panel fails.")

        px_tick = float(np.median(np.diff(ticks)))

        # pair each dark bar to its nearest light bar on the left.
        #
        # F-052: this was `dark_for[j] = d` with no collision check, so two dark
        # components claiming the same light bar silently kept the last one. A
        # 200 px^2 dark speck on the baseline displaced a real bar and was
        # reported as a measurement, with `bars_not_detected` empty -- a
        # measurement failure reported as a measurement, which is class 28
        # verbatim. Lowering the dark floor to 120 is what raised that exposure.
        # `len(db) == len(cats)` is NOT the right check: a dark bar exists only
        # for categories with any >1e6 hits, so db is legitimately shorter.
        dark_for = {}
        collisions = []
        for d in db:
            j = int(np.argmin([abs(d["x0"] - (l["x1"] + 8)) for l in lb]))
            if j in dark_for:
                collisions.append((pid, cats[j], dark_for[j]["x0"], d["x0"]))
            dark_for[j] = d
        assert not collisions, (
            f"two dark components matched the same light bar: {collisions}. "
            "One of them is not a >1e6 bar -- most likely a speck or a glyph "
            "above the min_area floor. Do not let it overwrite a real bar: a "
            "miss is not a measurement (F-052, class 28).")

        # SCALE, FROM AXIS FURNITURE ONLY (F-006). A least-squares line through
        # the detected tick centres, which is more robust than the median of
        # successive differences and gives a residual worth asserting on. Note
        # what is NOT in these four lines: the Fig 2e rates.
        tick_y = np.array(ticks, float)
        idx = np.arange(len(tick_y), dtype=float)
        A = np.vstack([idx, np.ones_like(idx)]).T
        (px_per_tick_fit, _), res, _, _ = np.linalg.lstsq(A, tick_y, rcond=None)
        px_per_tick_fit = abs(float(px_per_tick_fit))

        # F-053: TICK_PCT, asserted against the drawn axis rather than assumed.
        # The topmost detected tick is the axis maximum (100%), so the number of
        # tick intervals between the drawn baseline and that tick must be
        # 100 / TICK_PCT. This is the check the docstring claimed and did not
        # have: nothing in the file could previously have detected a figure
        # whose ticks are every 10%, and TICK_PCT multiplies all 48 numbers.
        n_intervals = (baseline - float(min(tick_y))) / px_per_tick_fit
        assert abs(n_intervals - 100.0 / TICK_PCT) < 0.35, (
            f"panel {pid}: the drawn axis spans {n_intervals:.2f} tick "
            f"intervals from the baseline {baseline} to the top tick "
            f"{min(tick_y):.1f}, but TICK_PCT={TICK_PCT} implies "
            f"{100.0 / TICK_PCT:.2f}. Either the axis maximum is not 100% or "
            f"the percent-per-tick is wrong -- and it multiplies every value "
            f"in this panel (F-053).")

        tick_resid = float(np.max(np.abs(tick_y - A @ np.linalg.lstsq(
            A, tick_y, rcond=None)[0])))
        slope = px_per_tick_fit / TICK_PCT
        units_per_tick = TICK_PCT

        # A log or linear axis has EQUAL tick spacing; an unequal one means the
        # anchors are glyph centroids, not ticks (error class 20). 2 px at
        # 600 dpi is about 0.08 mm on the page.
        assert tick_resid <= 2.5, (
            f"panel {pid}: tick centres are not equally spaced (max residual "
            f"{tick_resid:.2f} px over {len(ticks)} ticks). Equal spacing is a "
            f"property of the axis, so unequal spacing means these are not "
            f"ticks -- see error class 20.")
        assert abs(px_tick - px_per_tick_fit) <= 2.0, (
            f"panel {pid}: median tick diff {px_tick:.2f} disagrees with the "
            f"fitted spacing {px_per_tick_fit:.2f}")

        digit = marg(allrows, field, cats)
        panel = {"field": field, "categories": cats, "axis_x": int(axis_x),
                 "n_ticks": len(ticks),
                 "px_per_tick_median_diff": round(px_tick, 2),
                 "px_per_tick_lstsq": round(px_per_tick_fit, 3),
                 "tick_fit_max_residual_px": round(tick_resid, 2),
                 "percent_per_tick": TICK_PCT,
                 "px_per_percent": round(slope, 3),
                 "axis_intervals_baseline_to_top_tick": round(n_intervals, 3),
                 "axis_intervals_implied_by_TICK_PCT": 100.0 / TICK_PCT,
                 "light_min_area_largest_rejected": rejected_max,
                 "light_min_area_smallest_kept": accepted_min,
                 "scale_source": "y-axis tick centres only; NOT fitted to the "
                                 "Fig 2e rates being validated (F-006)",
                 "measured": {}, "digitized": {}, "delta": {}}
        # THE ORIGIN, ALSO FROM AXIS FURNITURE (round-2 finding E2). F-006
        # re-anchored the SCALE on the ticks and left the ORIGIN as three
        # hard-coded literals -- error class 20 surviving the repair meant to
        # clear it, on the other end of the same axis.
        #
        # The ladder runs 100% at the top in steps of TICK_PCT, so 0% sits
        # 100/TICK_PCT tick-intervals below the first tick. On 2h the 0% tick is
        # actually drawn and the literal 4055 was nearly right (0.28 pp), which
        # is why the old comment -- which quoted only the 2h number and blamed
        # the colour mask's bottom edge -- sounded adequate. On 2f and 2g the 0%
        # tick is NOT drawn (it merges into the x-axis rule), the ladder
        # extrapolates it 3-4 px BELOW the panel box, and the axis origin, not
        # the mask, is the dominant term: 0.61 pp on each.
        #
        # Bar-count weighted that predicts (14*0.61 + 10*0.61 + 24*0.28)/48 =
        # 0.44 pp of the recorded mean -0.52 pp, with 47 of 48 deltas negative.
        # The previous comment's "recorded here instead of being tuned away" was
        # right in spirit and wrong about the cause, so it under-reported a bias
        # that was removable rather than inherent.
        #
        # `baseline` is still what the bars are DETECTED on -- they really do sit
        # on the drawn rule. Only the conversion to percent uses the tick origin.
        intercept = float(np.linalg.lstsq(A, tick_y, rcond=None)[0][1])
        zero_y = intercept + px_per_tick_fit * (100.0 / TICK_PCT)
        panel_zero_shift = zero_y - baseline
        assert abs(panel_zero_shift) < 12.0, (
            f"panel {pid}: tick-derived 0% at y={zero_y:.1f} is "
            f"{panel_zero_shift:+.1f} px from the drawn baseline {baseline}. "
            f"More than a few px means the tick ladder was mis-detected, not "
            f"that the baseline is wrong.")

        # sub-pixel top edges, uniformly for every bar (see subpixel_top)
        for bar in lb + list(dark_for.values()):
            bar["top_subpixel"] = subpixel_top(img, bar, baseline)
            bar["height"] = zero_y - bar["top_subpixel"]

        for j, c in enumerate(cats):
            m5 = lb[j]["height"] / slope
            # F-007: this used to read `else 0.0`, so a bar the detector MISSED
            # was recorded as a bar that MEASURES ZERO. EP6's >1e6 bar is about
            # 50x9 px and was dropped by the min_area floor, after which the
            # fabricated 0.0 became the whole set's largest "disagreement" --
            # a measurement failure disguised as a measurement (class 28). A
            # miss is now None, so it cannot be averaged, deltaed, or reported
            # as agreement, and it is collected and asserted on below.
            m6 = dark_for[j]["height"] / slope if j in dark_for else None
            panel["measured"][c] = [round(float(m5), 1),
                                    None if m6 is None else round(float(m6), 1)]
            panel["digitized"][c] = [digit[c]["gt_1e5"], digit[c]["gt_1e6"]]
            panel["delta"][c] = [round(float(m5) - digit[c]["gt_1e5"], 1),
                                 None if m6 is None
                                 else round(float(m6) - digit[c]["gt_1e6"], 1)]
            if m6 is None:
                missed.append((pid, c, digit[c]["gt_1e6"]))
            out_rows.append({"panel": f"Fig {pid}", "field": field, "category": c,
                             "measured_pct_gt_1e5": round(float(m5), 1),
                             "measured_pct_gt_1e6": ("" if m6 is None
                                                     else round(float(m6), 1)),
                             "fig2e_pct_gt_1e5": digit[c]["gt_1e5"],
                             "fig2e_pct_gt_1e6": digit[c]["gt_1e6"],
                             "n_lipids": digit[c]["n"],
                             "source_figure": f"Fig {pid}",
                             "extraction_method": "bar_height_measured"})
            draw.append((lb[j], dark_for.get(j), baseline, c))
        d5 = [abs(v[0]) for v in panel["delta"].values()]
        d6 = [abs(v[1]) for v in panel["delta"].values() if v[1] is not None]
        panel["max_abs_delta_pct"] = [round(max(d5), 1),
                                      round(max(d6), 1) if d6 else None]
        panel["bars_not_detected"] = [c for c, v in panel["delta"].items()
                                      if v[1] is None]
        qc[pid] = panel

    # THE ASSERTION THIS FILE DID NOT HAVE. 48 printed numbers were being
    # compared against the heatmap read and the result was printed, never
    # checked -- so "48/48 agree" was a claim about a diagnostic, not a gate.
    worst = max((abs(d) for p in qc.values() if "error" not in p
                 for v in p["delta"].values() for d in v if d is not None),
                default=float("nan"))

    # THE GATE RUNS BEFORE ANYTHING IS WRITTEN (F-050).
    #
    # It used to run after, which had two consequences. A failing run still left
    # 05_fig2fgh_bars.csv and 08_fig2fgh_qc.json on disk for the next reader to
    # pick up as if they were good. And when EVERY panel failed, `out_rows` was
    # empty and the CSV write raised IndexError at `out_rows[0]` before the gate
    # was reached -- so the thing that stopped the run was an accident of the
    # writer, not a check, and it said "list index out of range" instead of
    # saying which panels failed. Found while seeding this check to confirm it
    # fires, which is exactly what CLAUDE.md section 6 requires.
    # F-050: THE GATE DID NOT COVER A FAILED PANEL. `worst` below filters on
    # `"error" not in p`, and `missed` is only appended inside the per-category
    # loop that the error branch `continue`s past -- so a panel that failed
    # detection was silently excluded from both, and "48/48 agree" could be
    # asserted on 34. Demonstrated by shrinking one light bar below the floor:
    # the script completed, exit 0, with 17 of 24 categories compared.
    #
    # This is F-038's shape, in the file whose own comment three lines down says
    # the 48 numbers "were being compared and the result was printed, never
    # checked". Printing a count is not checking it.
    failed = {pid: p["error"] for pid, p in qc.items() if "error" in p}
    assert not failed, (
        f"{len(failed)} of {len(PANELS)} panels failed detection and would be "
        f"EXCLUDED from the delta gate rather than failing it: {failed}. A "
        f"panel that was not measured is not a panel that agreed.")
    expected_cats = sum(len(c) for _, _, c, _, _ in PANELS)
    assert len(out_rows) == expected_cats, (
        f"{len(out_rows)} categories written, expected {expected_cats}. The "
        f"delta gate below only sees the rows that exist, so a short run would "
        f"pass it on a subset (F-050).")
    n_deltas = sum(1 for p in qc.values() if "error" not in p
                   for v in p["delta"].values() for d in v if d is not None)

    assert missed == [], (
        f"bars the detector did not find, recorded as None rather than 0.0: "
        f"{missed}. A miss is not a measurement -- lower min_area or fix the "
        f"mask, do not let it fall through to a number (F-007, class 28).")
    assert worst <= WORST_TOL, (
        f"max |delta| vs Fig 2e is {worst:.2f} pp, above the {WORST_TOL} pp "
        f"tolerance ({n_deltas} deltas compared). "
        f"With the scale now taken from the ticks alone, a disagreement here is "
        f"a real disagreement between this panel and the heatmap read -- it is "
        f"no longer absorbed by a fitted scale (F-006).")

    with open("05_fig2fgh_bars.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out_rows[0]))
        w.writeheader()
        w.writerows(out_rows)
    json.dump({"paper": "Lin et al., JACS 2026, Figs 2f-2h",
               "what": "hit rates (% of LNPs above RLU 1e5 and 1e6), the marginals "
                       "of the Fig 2e heatmap",
               "panels": qc}, open("08_fig2fgh_qc.json", "w"), indent=1)

    overlay(img, draw)
    for pid, p in qc.items():
        if "error" in p:
            print(f"  Fig {pid}: {p['error']}")
            continue
        print(f"  Fig {pid} ({p['field']}): max |delta| = "
              f"{p['max_abs_delta_pct'][0]} pp (>1e5), "
              f"{p['max_abs_delta_pct'][1]} pp (>1e6)   "
              f"scale {p['px_per_percent']} px/pp from ticks "
              f"({p['px_per_tick_lstsq']} px per {p['percent_per_tick']}%, "
              f"fit residual {p['tick_fit_max_residual_px']} px)"
              + (f"   NOT DETECTED: {p['bars_not_detected']}"
                 if p["bars_not_detected"] else ""))
        for c in p["categories"]:
            print(f"      {c:8} measured {p['measured'][c]}  fig2e {p['digitized'][c]}")
    print(f"\noverlay -> {REVIEW}/fig2fgh_overlay.png")


def marg(rows, field, cats):
    out = {}
    for c in cats:
        sub = [r for r in rows if r[field] == c]
        n = len(sub)
        out[c] = {"gt_1e5": round(100 * sum(1 for r in sub
                                            if int(r["bin_index"]) >= 3) / n, 1),
                  "gt_1e6": round(100 * sum(1 for r in sub
                                            if int(r["bin_index"]) == 4) / n, 1),
                  "n": n}
    return out


def overlay(img, draw):
    y0, y1, x0, x1 = 2700, 4200, 500, 2850
    crop = Image.fromarray(img.astype("uint8")[y0:y1, x0:x1])
    d = ImageDraw.Draw(crop)
    for lb, db, baseline, cat in draw:
        d.rectangle([lb["x0"] - x0, lb["top"] - y0, lb["x1"] - x0, baseline - y0],
                    outline=(200, 0, 0), width=3)
        if db:
            d.rectangle([db["x0"] - x0, db["top"] - y0, db["x1"] - x0, baseline - y0],
                        outline=(0, 0, 220), width=3)
        d.line([lb["x0"] - x0, baseline - y0, (db or lb)["x1"] - x0, baseline - y0],
               fill=(0, 160, 0), width=3)
    d.text((8, 8), "red = measured >1e5 bar   blue = measured >1e6 bar   "
                   "green = detected baseline", fill=(0, 0, 0))
    crop.save(f"{REVIEW}/fig2fgh_overlay.png")


if __name__ == "__main__":
    main()
