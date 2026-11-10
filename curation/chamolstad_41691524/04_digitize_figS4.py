#!/usr/bin/env python3
"""Digitise Figure S4 panels A (DC2.4) and C (C2C12) of nn5c19036_si_002.pdf.

TIER 5. Justified in 00_value_sources.json: the publisher Source Data workbook
holds 3'UTR half-lives only, no SI table carries per-formulation luciferase,
and there is no author repository.

Geometry comes from curation/lib (axis, bars). Nothing here supplies an origin
separately from the scale, and nothing is fitted onto a quantity it is later
checked against.

Bar-top convention: GraphPad strokes the bar rectangle ON the data path, so the
data value is the CENTRE of the top stroke, not its first or last row. Both
edges are measured with bars.bar_top -- once on a mask of the stroke (dark) and
once on a mask of the fill -- and the stroke width that falls out is reported
per bar, so a detector tracking something other than a stroke shows up as an
absurd width rather than as a plausible number.
"""
import sys, os, json, csv
import numpy as np
from PIL import Image
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import axis, bars
from errors import DetectionError

# The raster comes out of the SI PDF in `source/`, not from a path handed in on
# the command line. A tracked output that depends on a session scratchpad is
# regenerable from nothing once the scratchpad is cleared -- error class 38,
# which has already cost this project `build_manifest.py` (F-023) and
# `lit_review/build_data.py`. The extracted PNG is written beside this script
# so a reviewer can look at exactly the pixels that were measured.
SI = os.path.join(HERE, "source", "nn5c19036_si_002.pdf")
PAGE, XREF = 11, 170
import pymupdf
_doc = pymupdf.open(SI)
_ids = [i[0] for i in _doc[PAGE - 1].get_images(full=True)]
if _ids != [XREF]:
    raise SystemExit(f"SI p.11 holds images {_ids}, expected [170]")
_px = pymupdf.Pixmap(_doc, XREF)
if _px.n > 3:
    _px = pymupdf.Pixmap(pymupdf.csRGB, _px)
IMG = os.path.join(HERE, "fig_p11_x170.png")
_px.save(IMG)
im = np.array(Image.open(IMG).convert("RGB")).astype(int)
H, W, _ = im.shape
gray = im.mean(axis=2)
mx, mn = im.max(axis=2), im.min(axis=2)
sat = mx - mn

FILL = ((sat < 25) & (gray > 70) & (gray < 190)) | (sat >= 60)   # grey or colour
# The stroke is BLACK ink. The saturation term is not cosmetic: pure blue
# (0,0,255) has mean grey 85 and pure green (0,128,0) has 43, so a darkness-only
# mask swallows the highlighted bars whole and every coloured bar's window reads
# as solid ink from the axis ceiling down.
DARK = (gray < 90) & (sat < 60)
# The BAR SILHOUETTE: its stroke, its body, and its highlight colour. This is
# what the top-rule detector runs on. A stroke-only mask puts the window's
# edge column on an anti-aliased boundary, where a single pixel is 1/7 of a
# flank and costs 14% of a 90% coverage rule -- which refuses correct bars.
# The silhouette has no such boundary: the stroke and the body it encloses are
# one object. A whisker cap is in it too, and is rejected by what is UNDER it.
SOLID = FILL | DARK
WHITE = gray > 225

LABELS = (["SM102"] + [f"1-F{i}" for i in range(1, 17)]
          + [f"2-F{i}" for i in range(1, 28)])
assert len(LABELS) == 44, len(LABELS)

PANELS = {
    "S4A": {"cell": "DC2.4", "baseline": 263},
    "S4C": {"cell": "C2C12", "baseline": 1020},
}
# 0.0 .. 5.0 with a minor tick every 0.5: 11 positions. The LABELLED values are
# 0,1,2,3,4,5; the step is passed as a scalar so calibrate asserts uniform
# spacing on the positions it was actually given.
NROW = 2
N_TICKS = 11
TICK_VALUES = [5.0 - 0.5 * i for i in range(11)]   # pixel-ascending = value-descending

report = {"image": os.path.basename(IMG), "image_px": [W, H], "panels": {}}
rows_out = []

for pid, P in PANELS.items():
    base = P["baseline"]
    scan_lo, hi = max(0, base - 300), base - 1

    # ---- the y-axis rule: the longest dark column inside the panel band
    colcount = (gray[base - 200:hi, :] < 90).sum(axis=0)
    axis_col = int(np.argmax(colcount[:400]))
    assert colcount[axis_col] > 0.8 * 199, (pid, axis_col, colcount[axis_col])

    # ---- ticks, from the strip LEFT of the rule only
    ticks = axis.detect_ticks_1d(gray, "y", axis_col - 11, axis_col - 2,
                                 dark=140, min_len=2, max_len=12)
    ticks = sorted(t for t in ticks if scan_lo <= t <= base + 4)
    if len(ticks) != N_TICKS:
        raise SystemExit(f"{pid}: detected {len(ticks)} ticks {ticks}, "
                         f"need {N_TICKS}. REFUSING.")
    loo = axis.loo_residuals(ticks)
    cal = axis.calibrate(ticks, TICK_VALUES, kind="linear", axis="y",
                         loo_tol_px=axis.LOO_TOL_PX)
    # Scale AND origin both come from this one detected vector; there is no
    # separate origin argument in axis.calibrate, by design (F-033).
    lo = int(round(min(ticks)))          # the 5.0 tick: nothing above it is data

    # ---- bar x-windows from the CATEGORY TICKS under the baseline.
    # Not from the fill runs: a bar shorter than the probe strip has no fill to
    # find, so a fill-run detector silently drops the smallest values -- the
    # ones a reader most wants to know are small. The tick ladder is uniform
    # and complete whatever the bars do, and its own LOO residual is reported.
    xt = sorted(v for v in axis.detect_ticks_1d(gray, "x", base + 3, base + 11,
                                                dark=140, min_len=2, max_len=14)
                if v > axis_col + 4)
    if len(xt) != len(LABELS):
        raise SystemExit(f"{pid}: {len(xt)} category ticks, need {len(LABELS)}. "
                         f"REFUSING. {[round(v,1) for v in xt[:4]]}")
    xloo = axis.loo_residuals(xt)
    xpitch = float(np.mean(np.diff(xt)))

    # bar width: the median width of the unambiguous fill runs just above the
    # baseline. Measured, not assumed, and asserted against the tick pitch.
    probe = (FILL | DARK)[base - 9: base - 5, :].all(axis=0)
    wruns, s0 = [], None
    for x in range(axis_col + 6, W):
        if probe[x] and s0 is None: s0 = x
        elif not probe[x] and s0 is not None:
            if 8 <= x - s0 <= int(xpitch): wruns.append(x - s0)
            s0 = None
    bar_w = int(np.median(wruns))
    if not (0.5 * xpitch < bar_w < 0.95 * xpitch):
        raise SystemExit(f"{pid}: measured bar width {bar_w} px is not a "
                         f"plausible fraction of pitch {xpitch:.2f}. REFUSING.")
    # Per-bar windows from the bar's OWN silhouette, not a constant width.
    # A constant width is off by a pixel on about half the bars, and one stray
    # column in a 7 px flank is 0.857 coverage, which fails a 0.90 rule -- that
    # is a window error presenting itself as a class-32 refusal.
    # Probe ONE row, immediately above the axis rule. Every bar reaches the
    # baseline whatever its height, so this finds the short bars too; a
    # four-row probe higher up misses them and falls back to a nominal window
    # whose rounding is a pixel off, and one stray column in a 7 px flank is
    # 0.857 coverage against a 0.90 rule -- a window error wearing a class-32
    # refusal's clothes. Gaps between bars are ~6 px, so no two merge.
    sil = SOLID[base - 2, :]
    runs, win_src = [], []
    for c in xt:
        ci = int(round(c))
        if sil[ci]:
            a = ci
            while a > 0 and sil[a - 1]: a -= 1
            b = ci
            while b < W - 1 and sil[b + 1]: b += 1
            if 0.5 * bar_w <= b + 1 - a <= 1.5 * bar_w:
                runs.append((a, b + 1)); win_src.append("silhouette"); continue
        runs.append((int(round(c - bar_w / 2)), int(round(c + bar_w / 2))))
        win_src.append("nominal")   # bar too short to cast a silhouette

    # ---- the authors' own dashed reference line: drawn, never printed.
    # Measured ONLY in the columns BETWEEN bars. A bar band has the same
    # whole-width ink density as a 50%-duty dashed rule, so a whole-width
    # density test picks the tallest bar band instead -- it did, on S4A, and
    # reported 0.2267 for a line the authors drew at 1.0. The gaps are where a
    # full-width rule and a row of bars differ completely.
    gapcols = []
    for a, b in zip(xt[:-1], xt[1:]):
        m = int(round((a + b) / 2))
        gapcols.extend(range(m - 3, m + 4))
    gapcols = [c for c in gapcols if 0 <= c < W]
    band = gray[scan_lo:base - 2, gapcols] < 110
    dens = band.sum(axis=1) / band.shape[1]
    cand = [(i + scan_lo, d) for i, d in enumerate(dens) if d > 0.30]
    cand = [c for c in cand if abs(c[0] - base) > 6]
    if not cand:
        raise SystemExit(f"{pid}: dashed reference line not found. REFUSING.")
    grp, cur = [], [cand[0]]
    for c in cand[1:]:
        if c[0] - cur[-1][0] <= 2:
            cur.append(c)
        else:
            grp.append(cur); cur = [c]
    grp.append(cur)
    grp.sort(key=lambda gg: -sum(d for _, d in gg))
    g = grp[0]
    dash_row = sum(r * d for r, d in g) / sum(d for _, d in g)
    dash_value = float(cal.to_value(dash_row))
    dash_rows = {r for r, _ in g}
    dash_duty = round(float(max(d for _, d in g)), 3)
    n_dash_groups = len(grp)

    panel = {"axis_col": axis_col, "ticks_px": [round(t, 2) for t in ticks],
             "loo_residual_px_max": round(max(loo), 3),
             "px_per_unit": round(cal.px_per_unit, 4),
             "calib_residual_px": round(cal.residual_px, 4),
             "dashed_line_row": round(dash_row, 2),
             "dashed_line_value_PREDICTED": round(dash_value, 4),
             "dashed_line_rows_masked": sorted(dash_rows),
             "dashed_duty_cycle_in_gaps": dash_duty,
             "dashed_candidate_groups": n_dash_groups,
             "bar_scan_top_row": lo,
             "n_bars": len(runs), "x_tick_loo_max_px": round(max(xloo), 3),
             "x_pitch_px": round(xpitch, 3), "bar_width_px": bar_w,
             "bars": {}, "refusals": {},
             "windows_from_silhouette": win_src.count("silhouette")}

    # ---- the detector: the TOP OF THE SOLID COLUMN.
    #
    # Measured on this figure, not assumed: (i) the error bar's LOWER cap is
    # drawn INSIDE the bar and is full width, so a fill-first detector stops
    # 15 px below the real top (2-F23, S4A) and a fill-only read of S4A's
    # SM102 gave 0.86 for a bar that is 1.000 by construction; (ii) the
    # vertical whisker runs down the middle of the bar, so fill coverage never
    # reaches 0.90 while it does; (iii) pure blue has mean grey 85 and pure
    # green 43, so a darkness mask swallows the highlighted bars whole.
    #
    # So the measured object is the bar's SILHOUETTE -- stroke, body and
    # highlight colour as one -- and the rule is: walk UP from the baseline
    # while the row is still mostly covered. A bar is solid from the baseline
    # to its top rule by construction. Everything that could be mistaken for
    # the top has a WHITE GAP beneath it: an upper whisker cap, the dashed 1.0
    # rule above a short bar, a tick label, a neighbouring panel's glyph. One
    # uncovered row stops the walk, so none of them is reachable, and nothing
    # has to be masked out of the image to hide it from the detector.
    #
    # bars.bar_top is then run top-down over the same window as an INDEPENDENT
    # cross-check with its own full-width, flank and baseline guards, and the
    # agreement count is reported in 07_qc.json. It is a check, not the source:
    # a top-down scan cannot distinguish a whisker cap from a bar top without
    # knowing what is under it, which is the quantity this rule starts from.
    COVER = 0.50
    # For the cross-check only: a top-down scan would stop on the dashed rule
    # for every bar below 1.0, so the cross-check runs with those rows cleared.
    solid_nodash = SOLID.copy()
    for rr in range(min(dash_rows) - 1, max(dash_rows) + 2):
        solid_nodash[rr, :] = False
    panel["dash_rows_cleared_for_crosscheck"] = [min(dash_rows) - 1,
                                                 max(dash_rows) + 1]

    def stroke_thickness(r, x0, x1):
        # 0.50, not 0.90: the stroke's first row is anti-aliased against white
        # and a 0.90 rule scores it 0, which then reads as "no stroke at all".
        n = 0
        while r + n < base - 1 and DARK[r + n, x0:x1].mean() >= 0.50:
            n += 1
        return n

    for name, (x0, x1), wsrc in zip(LABELS, runs, win_src):
        rec = {"x0": x0, "x1": x1, "width_px": x1 - x0, "window": wsrc}
        cov = SOLID[:, x0:x1].mean(axis=1)
        if cov[base - 2] < COVER:
            panel["refusals"][name] = ("no bar: the row above the axis rule is "
                                       f"only {cov[base-2]:.2f} covered")
            continue
        r = base - 2
        while r > lo and cov[r - 1] >= COVER:
            r -= 1
        if r <= lo + 1:
            panel["refusals"][name] = "bar reaches the axis ceiling; off scale"
            continue
        # The walk stops on the topmost row that is at least half covered, and
        # that row can be the stroke's anti-aliased edge (grey ~150, inside the
        # fill band) rather than the stroke itself. So locate the stroke within
        # the next few rows and measure it from there; refuse if there is none.
        r_dark = next((rr for rr in range(r, r + 4)
                       if DARK[rr, x0:x1].mean() >= 0.50), None)
        if r_dark is None:
            panel["refusals"][name] = ("no black top rule within 3 px of the "
                                       f"silhouette top at row {r}")
            continue
        sw = stroke_thickness(r_dark, x0, x1)
        if not (1 <= sw <= 6):
            panel["refusals"][name] = f"stroke thickness {sw}px is implausible"
            continue
        height_px = (base - 2) - r
        centre = r_dark + (sw - 1) / 2.0     # the data path is the stroke centre
        val = float(cal.to_value(centre))

        # independent cross-check, top-down, with the library's own guards
        try:
            d = bars.bar_top(solid_nodash, x0, x1, lo, base - 2, nrow=NROW,
                             baseline_row=base)
            xcheck = d["row"] - r
            xmsg = None
        except DetectionError as e:
            xcheck, xmsg = None, str(e).split("\n")[0][:90]

        near_dash = abs(centre - dash_row) <= 6
        rec.update({"row_silhouette_top": r, "row_stroke_top": r_dark,
                    "antialias_offset_px": r_dark - r, "stroke_width_px": sw,
                    "bar_height_px": height_px,
                    "row_used": round(centre, 2), "value": round(val, 4),
                    "bar_top_crosscheck_delta_px": xcheck,
                    "bar_top_crosscheck_refusal": xmsg,
                    "top_within_6px_of_dashed_1.0": bool(near_dash)})
        panel["bars"][name] = rec
        rows_out.append({"panel": pid, "cell": P["cell"], "formulation": name,
                         "value_rel_SM102": round(val, 4),
                         "row_stroke_top_px": r_dark, "stroke_width_px": sw,
                         "antialias_offset_px": r_dark - r,
                         "row_used_px": round(centre, 2),
                         "bar_height_px": height_px, "bar_width_px": x1 - x0,
                         "window_from": wsrc,
                         "bar_top_crosscheck_delta_px": xcheck,
                         "near_dashed_line": int(near_dash)})

    sws = [b["stroke_width_px"] for b in panel["bars"].values()]
    panel["stroke_width_px"] = {"min": min(sws), "max": max(sws),
                                "median": float(np.median(sws))} if sws else None
    report["panels"][pid] = panel

with open(os.path.join(HERE, "04_figS4_geometry.json"), "w") as f:
    json.dump(report, f, indent=1)
with open(os.path.join(HERE, "04_figS4_values.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows_out[0].keys())); w.writeheader()
    w.writerows(rows_out)

for pid, p in report["panels"].items():
    print(f"\n== {pid} ticks={p['ticks_px']} LOOmax={p['loo_residual_px_max']} "
          f"px/unit={p['px_per_unit']} resid={p['calib_residual_px']}")
    print(f"   DRAWN dashed line predicts {p['dashed_line_value_PREDICTED']} "
          f"(authors never printed it; a round 1.0 is the target)")
    print(f"   bars={p['n_bars']} bar_w={p['bar_width_px']}px stroke={p['stroke_width_px']}")
    print(f"   refusals={len(p['refusals'])} {list(p['refusals'])}")
    print("   " + "  ".join(f"{k}={v['value']:.2f}" for k, v in list(p["bars"].items())[:8]))
