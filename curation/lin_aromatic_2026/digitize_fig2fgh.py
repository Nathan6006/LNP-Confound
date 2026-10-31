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
tick, the tick interval is asserted to round to TICK_PCT, and `rep_l` appears
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
# 1.1, set against a MEASURED residual of max 1.00 pp and mean -0.52 pp, and
# deliberately not set to whatever the current worst value happens to be. The
# remaining bias is real and one-sided -- 47 of 48 deltas are negative -- and its
# probable cause is the BOTTOM edge: a bar sits on the axis rule, and the colour
# mask's last saturated row is 1-2 px above the rule's top, so every bar is
# short by that much (1.5 px / 5.35 px per pp = 0.28 pp on panel 2h). Correcting
# the bottom edge sub-pixel as well as the top would remove most of it.
#
# NOT chased further, on purpose. These 48 numbers are a VALIDATION of the Fig 2e
# heatmap read; no corpus value is computed from them (`build_rows.py` reads
# 06_merged.csv, not 05_fig2fgh_bars.csv), so the residual costs accuracy in a
# check rather than in the data. The reason it is recorded here instead of being
# tuned away is that tuning it away is exactly what the old fitted scale did.
WORST_TOL = 1.1


def bars_in(mask, baseline, xwin, min_area=3000):
    """Bars whose bottom sits on this panel's baseline."""
    x0, x1 = xwin
    lab, _ = ndimage.label(mask[:, x0:x1])
    out = []
    for i, sl in enumerate(ndimage.find_objects(lab)):
        area = int((lab[sl] == i + 1).sum())
        if area < min_area:
            continue
        bot, top = sl[0].stop, sl[0].start
        if abs(bot - baseline) > 6:
            continue
        out.append({"x0": sl[1].start + x0, "x1": sl[1].stop + x0,
                    "top": top, "bottom": bot, "height": bot - top})
    out.sort(key=lambda d: d["x0"])
    return out


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
        lb = bars_in(light, baseline, xwin)
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
        px_tick = float(np.median(np.diff(ticks)))

        # pair each dark bar to its nearest light bar on the left
        dark_for = {}
        for d in db:
            j = int(np.argmin([abs(d["x0"] - (l["x1"] + 8)) for l in lb]))
            dark_for[j] = d

        # SCALE, FROM AXIS FURNITURE ONLY (F-006). A least-squares line through
        # the detected tick centres, which is more robust than the median of
        # successive differences and gives a residual worth asserting on. Note
        # what is NOT in these four lines: the Fig 2e rates.
        tick_y = np.array(ticks, float)
        idx = np.arange(len(tick_y), dtype=float)
        A = np.vstack([idx, np.ones_like(idx)]).T
        (px_per_tick_fit, _), res, _, _ = np.linalg.lstsq(A, tick_y, rcond=None)
        px_per_tick_fit = abs(float(px_per_tick_fit))
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
                 "scale_source": "y-axis tick centres only; NOT fitted to the "
                                 "Fig 2e rates being validated (F-006)",
                 "measured": {}, "digitized": {}, "delta": {}}
        # sub-pixel top edges, uniformly for every bar (see subpixel_top)
        for bar in lb + list(dark_for.values()):
            bar["top_subpixel"] = subpixel_top(img, bar, baseline)
            bar["height"] = bar["bottom"] - bar["top_subpixel"]

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

    with open("05_fig2fgh_bars.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out_rows[0]))
        w.writeheader()
        w.writerows(out_rows)
    json.dump({"paper": "Lin et al., JACS 2026, Figs 2f-2h",
               "what": "hit rates (% of LNPs above RLU 1e5 and 1e6), the marginals "
                       "of the Fig 2e heatmap",
               "panels": qc}, open("08_fig2fgh_qc.json", "w"), indent=1)

    # THE ASSERTION THIS FILE DID NOT HAVE. 48 printed numbers were being
    # compared against the heatmap read and the result was printed, never
    # checked -- so "48/48 agree" was a claim about a diagnostic, not a gate.
    worst = max((abs(d) for p in qc.values() if "error" not in p
                 for v in p["delta"].values() for d in v if d is not None),
                default=float("nan"))
    assert missed == [], (
        f"bars the detector did not find, recorded as None rather than 0.0: "
        f"{missed}. A miss is not a measurement -- lower min_area or fix the "
        f"mask, do not let it fall through to a number (F-007, class 28).")
    assert worst <= WORST_TOL, (
        f"max |delta| vs Fig 2e is {worst:.2f} pp, above the 1.5 pp tolerance. "
        f"With the scale now taken from the ticks alone, a disagreement here is "
        f"a real disagreement between this panel and the heatmap read -- it is "
        f"no longer absorbed by a fitted scale (F-006).")

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
