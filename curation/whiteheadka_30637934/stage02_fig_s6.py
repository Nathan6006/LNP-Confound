#!/usr/bin/env python3
"""Figure S6 (SI p.12) -- tier 4b, the PDF's own vector coordinates.

NO PIXEL IS READ. The page is never rendered. Every number below comes from
path coordinates in the content stream, so there is no raster mask, no colour
distance and no dpi.

WHAT SEPARATES A BAR FROM A WHISKER HERE. Not a coverage threshold -- an EXACT
WIDTH. GraphPad emits, per bar, three paths in the same frame:

    bar outline     8.51 pt wide, top horizontal segment = the mean
    error-bar cap   4.19 pt wide, above the bar top      = mean + s.d.
    pattern tile    2.80 pt square                       = fill, not geometry

That is the same discriminator `curation/lib/bars.py` implements for rasters
(a top rule spans the full bar; a centred marker is in neither flank) but
exact instead of thresholded, which is why `bars.bar_top` is not called: it
takes an ink array and there is no ink. The hazard it guards -- error class 32,
a replicate marker read as a bar top (F-018, 15 of 60 bars, up to +28%) -- is
closed here by refusing any candidate whose width is not 8.51 +/- 0.05.

WHICH FRAME. The bar OUTLINE, not the fill rect. The fill rects sit 0.064 pt
inside the stroke frame (every fill baseline is 220.39 while the axis line and
the y=0 tick are both at 220.45), and the outline's legs start exactly on the
y=0 tick. Using the outline means bar tops and tick positions are in one
coordinate system, so no origin is supplied separately -- that is the shape
error class 20 / F-033 took. The 0.064 pt offset is measured and reported, not
assumed away: it is worth 1,309 RLU, which matters only for the floor bars.

CATEGORY ASSIGNMENT. From the 14 x-axis TICK MARKS, which are axis furniture,
cross-checked two independent ways (see 07_qc.json): every tick lands on a bar
centre to <=0.07 pt against a 12.81 pt pitch, and every x LABEL's right edge
sits at tick + 2.97 pt to +/-0.05 pt across all 11 O-labels. A one-slot error
would be 12.81 pt, 180x the larger residual.
"""
import csv
import json
import os
import sys

import pymupdf

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "lib"))
import axis                                                      # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
PDF = os.path.join(HERE, "source", "smll201805097-sup-0001-S1.pdf")
PAGE = 12                      # 1-based: Figure S6

BAR_W = 8.51                   # bar outline width, pt
CAP_W = 4.19                   # error-bar cap width, pt
TILE = 2.80                    # pattern tile / series marker, pt
W_TOL = 0.05

# Category order exactly as the x labels read, left to right.
PANEL_A_CATS = ["Crl", "Nkd", "Lipo MM", "O6", "O8", "O9", "Oi10", "O10",
                "O11", "O12", "O13", "O14", "O16", "O18"]
PANEL_B_CATS = ["306Oi10", "306O10", "306O11", "306O12"]
PANEL_B_SERIES = {"10 ng": (0.643, 0.333, 0.231),
                  "50 ng": (0.615, 0.228, 0.369),
                  "100 ng": (0.365, 0.479, 0.169)}


def items(page, xlo, xhi, ylo, yhi):
    """Every path item inside a box, as (kind, group, item)."""
    out = []
    for g in page.get_drawings():
        for it in g["items"]:
            if it[0] == "l":
                a, b = it[1], it[2]
                if xlo <= min(a.x, b.x) and max(a.x, b.x) <= xhi \
                        and ylo <= min(a.y, b.y) and max(a.y, b.y) <= yhi:
                    out.append(("l", g, (a.x, a.y, b.x, b.y)))
            elif it[0] == "re":
                r = it[1]
                if xlo <= r.x0 and r.x1 <= xhi and ylo <= r.y0 and r.y1 <= yhi:
                    out.append(("re", g, (r.x0, r.y0, r.x1, r.y1)))
    return out


def horiz(segs, width, tol=W_TOL):
    """Horizontal segments of a given width. Returns [(xc, y, x0, x1)]."""
    out = []
    for kind, _g, v in segs:
        if kind != "l":
            continue
        x0, y0, x1, y1 = v
        if abs(y0 - y1) > 0.05:
            continue
        w = abs(x1 - x0)
        if abs(w - width) <= tol:
            out.append(((x0 + x1) / 2.0, y0, min(x0, x1), max(x0, x1)))
    return sorted(out)


def vert_ticks(segs, axis_x, length_lo=2.0, length_hi=3.5):
    """Short vertical segments hanging off a horizontal axis line = x ticks."""
    out = []
    for kind, _g, v in segs:
        if kind != "l":
            continue
        x0, y0, x1, y1 = v
        if abs(x0 - x1) > 0.05:
            continue
        L = abs(y1 - y0)
        if length_lo <= L <= length_hi and min(y0, y1) >= axis_x - 0.1:
            out.append(round(x0, 3))
    return sorted(set(out))


def y_ticks(segs, x_lo, x_hi):
    """Short horizontal segments to the LEFT of the y spine = y ticks."""
    out = []
    for kind, _g, v in segs:
        if kind != "l":
            continue
        x0, y0, x1, y1 = v
        if abs(y0 - y1) > 0.05:
            continue
        if abs(min(x0, x1) - x_lo) < 0.1 and abs(max(x0, x1) - x_hi) < 0.1:
            out.append(round(y0, 3))
    return sorted(set(out))


def labels(page, xlo, xhi, ylo, yhi):
    return [(w[0], w[1], w[2], w[3], w[4]) for w in page.get_text("words")
            if xlo <= w[0] and w[2] <= xhi and ylo <= w[1] and w[3] <= yhi]


def main():
    doc = pymupdf.open(PDF)
    page = doc[PAGE - 1]
    qc = {"_what": "Figure S6, tier 4b. Every number is a PDF path coordinate.",
          "pdf_page": PAGE, "rendered_to_raster": False}

    # ================= PANEL A =========================================
    segA = items(page, 100, 310, 60, 262)
    ytk = y_ticks(segA, 111.71, 114.43)
    if len(ytk) != 6:
        raise SystemExit(f"panel A: {len(ytk)} y ticks, expected 6: {ytk}")
    # Axis labels read from the page text: 0.0, 5.0e5, 1.0e6 ... 2.5e6.
    calA = axis.calibrate(ytk, [2.5e6, 2.0e6, 1.5e6, 1.0e6, 5.0e5, 0.0],
                          kind="linear", axis="y")
    looA = axis.loo_residuals(ytk)
    qc["panelA_y_calibration"] = calA.qc()
    # Calibration.qc() rounds px_per_unit to 4 dp and RLU-per-pixel is 4.9e-5,
    # so it prints as 0.0. Record the same fit in readable units here rather
    # than paraphrasing it (F-017: never describe a number in prose instead of
    # recording it).
    qc["panelA_px_per_5e5_RLU"] = round(abs(calA.to_pixel(0.0)
                                            - calA.to_pixel(5.0e5)), 4)
    qc["panelA_RLU_per_px"] = round(abs(calA.to_value(0.0)
                                        - calA.to_value(1.0)), 2)
    qc["panelA_y_tick_positions_px"] = ytk
    qc["panelA_y_loo_residuals_px"] = [round(v, 4) for v in looA]
    qc["panelA_y_loo_worst_px"] = round(max(looA), 4)

    xtk = vert_ticks(segA, 220.45)
    qc["panelA_x_tick_positions_px"] = xtk
    if len(xtk) != len(PANEL_A_CATS):
        raise SystemExit(f"panel A: {len(xtk)} x ticks for "
                         f"{len(PANEL_A_CATS)} categories: {xtk}")
    pitch = [round(xtk[i + 1] - xtk[i], 3) for i in range(len(xtk) - 1)]
    qc["panelA_x_tick_pitch_px"] = pitch
    qc["panelA_x_tick_pitch_spread_px"] = round(max(pitch) - min(pitch), 4)

    tops = horiz(segA, BAR_W)
    caps = horiz(segA, CAP_W)
    zero_px = calA.to_pixel(0.0)
    qc["panelA_zero_from_tick_fit_px"] = round(zero_px, 4)
    qc["panelA_axis_line_and_zero_tick_px"] = 220.45
    fill_baselines = sorted({round(v[3], 3) for k, _g, v in segA
                             if k == "re" and abs((v[2] - v[0]) - BAR_W) <= W_TOL})
    qc["panelA_fill_rect_baselines_px"] = fill_baselines
    qc["panelA_fill_vs_stroke_frame_offset_px"] = (
        round(220.45 - fill_baselines[0], 4) if fill_baselines else None)
    qc["panelA_frame_offset_worth_RLU"] = round(
        abs(calA.to_value(220.45 - (220.45 - fill_baselines[0]))
            - calA.to_value(220.45)), 1) if fill_baselines else None

    rowsA, assign_resid, label_resid = [], [], []
    labs = {l[4]: l for l in labels(page, 100, 300, 222, 250)}
    for cat, tick in zip(PANEL_A_CATS, xtk):
        hit = [t for t in tops if abs(t[0] - tick) <= 0.30]
        cap = [c for c in caps if abs(c[0] - tick) <= 0.30]
        if len(hit) > 1:
            raise SystemExit(f"panel A {cat}: {len(hit)} bar tops at {tick}")
        if hit:
            assign_resid.append(abs(hit[0][0] - tick))
        # The label-right-edge identity is asserted only over the 11
        # HORIZONTAL O-labels. "Crl", "Nkd" and "Lipo MM" are rotated -- their
        # bounding boxes are square (w == h to 0.02 pt), so their box edge is
        # not a text edge and mixing them in would inflate the spread to
        # 1.88 pt and hide the 0.05 pt agreement that actually pins the grid.
        if cat in labs and cat.startswith("O"):
            label_resid.append(round(labs[cat][2] - tick, 3))
        mean = calA.to_value(hit[0][1]) if hit else None
        sd_top = calA.to_value(cap[0][1]) if cap else None
        rowsA.append({
            "panel": "S6A", "category": cat, "x_tick_px": tick,
            "bar_top_px": round(hit[0][1], 3) if hit else "",
            "bar_outline_width_px": round(hit[0][3] - hit[0][2], 3) if hit else "",
            "errorcap_px": round(cap[0][1], 3) if cap else "",
            "value_RLU": round(mean) if mean is not None else "",
            "mean_plus_sd_RLU": round(sd_top) if sd_top is not None else "",
            "sd_RLU": (round(sd_top - mean) if (mean is not None
                                                and sd_top is not None) else ""),
            "has_drawn_bar": bool(hit),
        })
    qc["panelA_bar_centre_vs_xtick_worst_px"] = round(max(assign_resid), 4)
    qc["panelA_label_right_edge_minus_tick_px"] = label_resid
    qc["panelA_label_offset_n"] = len(label_resid)
    qc["panelA_label_offset_mean_px"] = (
        round(sum(label_resid) / len(label_resid), 4) if label_resid else None)
    qc["panelA_label_offset_spread_px"] = (
        round(max(label_resid) - min(label_resid), 4) if label_resid else None)
    qc["panelA_bars_drawn"] = sum(r["has_drawn_bar"] for r in rowsA)
    qc["panelA_categories"] = len(PANEL_A_CATS)

    # ================= PANEL B =========================================
    segB = items(page, 320, 560, 60, 262)
    ytkB = y_ticks(segB, 360.57, 362.81)
    if len(ytkB) != 5:
        raise SystemExit(f"panel B: {len(ytkB)} y ticks: {ytkB}")
    calB = axis.calibrate(ytkB, [2.0e6, 1.5e6, 1.0e6, 5.0e5, 0.0],
                          kind="linear", axis="y")
    looB = axis.loo_residuals(ytkB)
    qc["panelB_y_calibration"] = calB.qc()
    qc["panelB_px_per_5e5_RLU"] = round(abs(calB.to_pixel(0.0)
                                            - calB.to_pixel(5.0e5)), 4)
    qc["panelB_RLU_per_px"] = round(abs(calB.to_value(0.0)
                                        - calB.to_value(1.0)), 2)
    qc["panelB_y_tick_positions_px"] = ytkB
    qc["panelB_y_loo_residuals_px"] = [round(v, 4) for v in looB]
    qc["panelB_y_loo_worst_px"] = round(max(looB), 4)

    # Series markers: the OUTER 2.80 pt tile of each symbol, by fill colour.
    # Legend swatches live at x >= 525 and are excluded by the x window.
    marks = []
    for kind, g, v in items(page, 360, 525, 60, 262):
        if kind != "re" or g.get("fill") is None:
            continue
        x0, y0, x1, y1 = v
        if abs((x1 - x0) - TILE) > 0.05 or abs((y1 - y0) - TILE) > 0.05:
            continue
        col = tuple(round(c, 3) for c in g["fill"])
        marks.append(((x0 + x1) / 2.0, (y0 + y1) / 2.0, col))
    marks = sorted(set(marks))
    qc["panelB_markers_found"] = len(marks)
    rowsB = []
    for ci, cat in enumerate(PANEL_B_CATS):
        grp = [m for m in marks if 360 + 40 * ci < m[0] < 400 + 40 * ci]
        for ser, col in PANEL_B_SERIES.items():
            hit = [m for m in grp if m[2] == col]
            if len(hit) != 1:
                raise SystemExit(f"panel B {cat} {ser}: {len(hit)} markers")
            rowsB.append({"panel": "S6B", "category": cat, "series": ser,
                          "marker_centre_px": round(hit[0][1], 3),
                          "value_RLU": round(calB.to_value(hit[0][1]))})
    # Dose monotonicity, per lipid -- also validates the legend colour order.
    mono = []
    for cat in PANEL_B_CATS:
        v = [next(r["value_RLU"] for r in rowsB
                  if r["category"] == cat and r["series"] == s)
             for s in ("10 ng", "50 ng", "100 ng")]
        mono.append({"category": cat, "values": v,
                     "monotone_increasing": v[0] < v[1] < v[2]})
    qc["panelB_dose_monotonicity"] = mono
    qc["panelB_monotone_count"] = f"{sum(m['monotone_increasing'] for m in mono)}/{len(mono)}"

    # ---- the cross-panel comparison. Measured, reported, NEVER averaged. ----
    shared = []
    for cat in PANEL_B_CATS:
        tok = cat.replace("306", "")
        a = next((r["value_RLU"] for r in rowsA if r["category"] == tok), None)
        b = next((r["value_RLU"] for r in rowsB
                  if r["category"] == cat and r["series"] == "100 ng"), None)
        shared.append({"lipid": tok, "S6A_100ng": a, "S6B_100ng": b,
                       "ratio_B_over_A": round(b / a, 4)})
    ra = {s["lipid"]: i for i, s in enumerate(
        sorted(shared, key=lambda s: -s["S6A_100ng"]))}
    rb = {s["lipid"]: i for i, s in enumerate(
        sorted(shared, key=lambda s: -s["S6B_100ng"]))}
    n = len(shared)
    d2 = sum((ra[k] - rb[k]) ** 2 for k in ra)
    qc["cross_panel_S6A_vs_S6B_100ng"] = {
        "_what": ("Both panels are HeLa / firefly-luciferase / RLU / 24 h, and "
                  "S6A's single dose IS 100 ng, so S6B's 100 ng series must "
                  "reproduce S6A on these four lipidoids. It does not: see "
                  "08_contradictions.csv. Neither panel is averaged into the "
                  "other and neither is adjusted by the other."),
        "shared": shared,
        "ratio_min": min(s["ratio_B_over_A"] for s in shared),
        "ratio_max": max(s["ratio_B_over_A"] for s in shared),
        "rank_S6A_desc": [s["lipid"] for s in sorted(shared, key=lambda s: -s["S6A_100ng"])],
        "rank_S6B_desc": [s["lipid"] for s in sorted(shared, key=lambda s: -s["S6B_100ng"])],
        "spearman_rho": round(1 - 6 * d2 / (n * (n * n - 1)), 4),
    }

    with open(os.path.join(HERE, "02_fig_s6_geometry.json"), "w") as fh:
        json.dump(qc, fh, indent=1)
    with open(os.path.join(HERE, "02_fig_s6_panelA.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rowsA[0].keys()))
        w.writeheader(); w.writerows(rowsA)
    with open(os.path.join(HERE, "02_fig_s6_panelB.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rowsB[0].keys()))
        w.writeheader(); w.writerows(rowsB)

    print(json.dumps({k: v for k, v in qc.items()
                      if k.startswith(("panelA_y", "panelA_x_tick_pitch_s",
                                       "panelA_bar_centre", "panelA_label_offset",
                                       "panelA_bars", "panelA_fill", "panelA_frame",
                                       "panelB_y_loo", "panelB_markers",
                                       "panelB_monotone"))}, indent=1))
    print("\nPANEL A")
    for r in rowsA:
        print(f"  {r['category']:8s} tick={r['x_tick_px']:7.2f} "
              f"top={str(r['bar_top_px']):8s} w={str(r['bar_outline_width_px']):6s} "
              f"RLU={str(r['value_RLU']):>9s}  sd={str(r['sd_RLU']):>9s}")
    print("\nPANEL B")
    for r in rowsB:
        print(f"  {r['category']:9s} {r['series']:7s} RLU={r['value_RLU']:>9d}")
    print("\nCROSS-PANEL")
    print(json.dumps(qc["cross_panel_S6A_vs_S6B_100ng"], indent=1))


if __name__ == "__main__":
    main()
