#!/usr/bin/env python3
"""Digitize Figure 2c: the three-tail lipids above 2e4, on a log axis.

WHY THIS MATTERS for the dataset. The Fig 2b colorbar stops at 6e4, so every
strong hit is clipped to one blue and 21 of the 320 cells carry no recoverable
value. Figure 2c re-plots those same lipids on a log axis reaching ~1e6, so it
supplies the real number for most of them. It also gives the four benchmark LNPs
(DLin-MC3-DMA, SM102, ALC-0315, Lipid 5), which are not part of the 344.

LABELS are taken from the printed x axis, read at 600 dpi, not inferred from the
heatmap: 8 lipids on amine 110, then 20 on T3A, then 4 benchmarks = 32 bars,
which matches the 32 bar runs detected. Ordering is (amine, BL, tail) ascending.

The set is itself a check on Fig 2b, and it passes for 27 of 28: every one of
the 8 lipids 2c plots for amine 110 is exactly a cell the heatmap puts above
2e4, and 19 of its 20 T3A lipids likewise. The single exception, T3A-L8-C14, is
plotted here but its heatmap cell is visibly cream (7.9e3) - an inconsistency
between two figures of the same paper, recorded rather than reconciled.

Bar top is read as the first row whose ink JOINS the two side strokes, because
the three replicate dots sit on the bar top and a coverage rule alone fires on
them. See bar_top() and provenance/review/round_2026-10-30/remeasure_QL_2026.md.
"""
import csv
import json
import os

import numpy as np
import PIL.Image as Image
import PIL.ImageDraw as ImageDraw

PAGE = "00_pages/fig_p4_600.png"
REVIEW = "review"

BASE_Y = 3280.5           # baseline rule, 3272-3289
# y-axis calibration from the TICK MARKS, re-derived 2026-10-30 (F-017).
#
# The old comment read "y-axis label centroids: 1e6 2342, 1e5 2748, 1e4 3156",
# which is an error-class-21 defect: those are not the label centroids (those
# measure 2344.3 / 2750.8 / 3158.9) and they are not the tick positions either --
# they are the tick positions ROUNDED TO INTEGERS. The rounding is what made the
# decade gaps look unequal (406 and 408), and I took that as proof the code was
# anchored on glyph centroids. It was not. **The code was right and its comment
# was wrong**, which is the opposite of what F-017 claimed, and I had marked that
# finding `confirmed` on the arithmetic of three rounded integers in a comment.
#
# Measured tick centroids, native embedded JPEG -> 600 dpi:
#   1e6  2342.78    1e5  2749.54    1e4  3156.75
# gaps 406.76 / 407.21; linear fit of log10(v) on y gives 406.99 px/decade with
# residuals +/-0.037 px (rms 2.6e-4 decades). So PX_PER_DECADE = 407.0 was right
# to 0.01 px and only Y_1E4 moves, by 0.75 px = +0.42% in value.
#
# Validated against a round number the paper never printed: Fig 2c's red dotted
# reference line has centroid y = 2872.21, which these constants put at
# 50,020 RLU -- 0.04% off a round 5e4. The stored constants give 49,804 and an
# earlier reviewer's proposed 406.0 / 3153.0 give 49,159, so this calibration is
# ~40x closer to the line the authors drew.
Y_1E4 = 3156.75
PX_PER_DECADE = 407.0
TOP_SEARCH = 2320

LAB_110 = ["L1-C10", "L2-C8", "L2-C10", "L3-C6", "L4-C8", "L7-C10",
           "L8-C8", "L8-C10"]
LAB_T3A = ["L1-C10", "L2-C8", "L2-C10", "L2-C12", "L3-C6", "L3-C10",
           "L4-C8", "L4-C10", "L4-C12", "L5-C8", "L5-C10", "L5-C14",
           "L6-C14", "L7-C10", "L7-C12", "L7-C14", "L8-C8", "L8-C10",
           "L8-C12", "L8-C14"]
BENCH = ["DLin-MC3-DMA", "SM102", "ALC-0315", "Lipid 5"]
LABELS = [f"110-{s}" for s in LAB_110] + [f"T3A-{s}" for s in LAB_T3A] + BENCH


INK_THR = 172            # "anything drawn": see bar_top() below


def _runs(idx, gap=0):
    """Split a sorted index list into runs, allowing `gap` missing entries.

    gap=0 means STRICTLY contiguous. That matters: bridging even a 2 px gap is
    enough to let a merged replicate-dot cluster pass the test below (it does
    on T3A-L2-2C10 and on ALC-0315 in Fig 3c).
    """
    out, cur = [], [idx[0]]
    for i in idx[1:]:
        if i - cur[-1] > gap + 1:
            out.append(cur)
            cur = []
        cur.append(i)
    out.append(cur)
    return out


def bar_top(ink, x0, x1, top_search, base_y, frac=0.72, nrow=2,
            band=(24, 8), min_stroke=3):
    """Topmost row of a bar's top rule, immune to the replicate dot markers.

    THE GEOMETRY, measured rather than assumed. All 60 bars in Figs 2c and 3c
    are HOLLOW: in a band just above the baseline each one shows exactly two
    fully-coloured column groups (the side strokes, 4-14 px wide) with 57-76%
    white between them. Fig 2c bars are 66-72 px wide, Fig 3c bars 46-50 px.
    None is solid. An earlier attempt at this fix read the 12-column LEFT
    STROKE of the first bar (110-L1-C10, cols 557-568) as if it were the bar
    and concluded the bars were 12 px wide and solid; the bar is cols 557-627.
    The branch below still handles a solid bar, so the two cases share one code
    path, but on this figure that branch never fires.

    THE RULE. Take the first row carrying ONE STRICTLY CONTIGUOUS run of ink
    inside [x0, x1) that spans at least `frac` of the bar width AND contains
    both side strokes' centres, sustained for `nrow` rows (2 = one native pixel
    of the 300 dpi source). A dot cluster cannot satisfy it: it does not reach
    both strokes without a gap. The bar's own top rule does, because it is the
    thing that joins them.

    `ink` must be a LUMINANCE mask ("anything drawn"), not a blue mask. Fig 2c's
    red dotted 5e4 reference line crosses the top rule of DLin-MC3-DMA and
    SM102; on a blue-only mask the rule is punched through there and DLin-MC3-DMA
    returns no top at all. The consequence is that `top_search` must now clear
    any BLACK furniture too -- see the assert in each caller.

    Returns (y_top, y_bot, frac_used), or (None, None, None).
    """
    import numpy as np
    w = x1 - x0
    # 1. side strokes, from a band just above the baseline where the hollow
    #    interior is unambiguously white.
    colfrac = ink[base_y - band[0]:base_y - band[1], x0:x1].mean(0)
    on = [i for i, v in enumerate(colfrac) if v > 0.95]
    groups = [g for g in _runs(on) if len(g) >= min_stroke] if on else []
    if len(groups) >= 2:               # hollow bar: outermost two groups
        cL = (groups[0][0] + groups[0][-1]) / 2.0
        cR = (groups[-1][0] + groups[-1][-1]) / 2.0
    elif len(groups) == 1:             # solid bar: quarter points of the fill
        a, b = groups[0][0], groups[0][-1]
        cL, cR = a + (b - a) * 0.25, a + (b - a) * 0.75
    else:                              # nothing resolvable: just inside the edges
        cL, cR = 1.0, w - 2.0

    # 2. scan down for the joined top rule. The fallback fracs never fire on
    #    this figure; `frac_used` is written to the CSV so that if one ever
    #    does, the bar is visible as needing review rather than silently rescued.
    for f in (frac, 0.66, 0.60):
        need = f * w
        ok = np.zeros(base_y - 6 - top_search, dtype=bool)
        for k in range(len(ok)):
            idx = np.flatnonzero(ink[top_search + k, x0:x1])
            if idx.size == 0:
                continue
            for sg in _runs(list(idx)):
                if (sg[-1] - sg[0] + 1) >= need and sg[0] <= cL and sg[-1] >= cR:
                    ok[k] = True
                    break
        for k in range(len(ok) - nrow + 1):
            if ok[k:k + nrow].all():
                j = k
                while j < len(ok) and ok[j]:
                    j += 1
                return top_search + k, top_search + j - 1, f
    return None, None, None


def main():
    os.makedirs(REVIEW, exist_ok=True)
    img = np.array(Image.open(PAGE).convert("RGB")).astype(int)
    sat = img.max(2) - img.min(2)
    coloured = (img[..., 2] > img[..., 0] + 40) & (sat > 60)
    # Bar RUNS are found on the blue mask (the ink mask would merge the black
    # axis into them); bar TOPS are read on the ink mask (the blue mask breaks
    # where the red 5e4 reference line crosses a top rule).
    ink = img.mean(2) <= INK_THR

    prof = np.array([coloured[TOP_SEARCH:int(BASE_Y) - 10, x].sum()
                     for x in range(530, img.shape[1])])
    on = prof >= 3
    runs, cur = [], None
    for i, v in enumerate(on):
        if v and cur is None:
            cur = i
        elif not v and cur is not None:
            if i - cur > 10:
                runs.append((530 + cur, 530 + i))
            cur = None
    if cur is not None and len(on) - cur > 10:
        runs.append((530 + cur, 530 + len(on)))

    assert len(runs) == len(LABELS), \
        f"found {len(runs)} bars, expected {len(LABELS)}"

    # TOP_SEARCH must clear every black element above the bars, because the
    # top rule is now read on a luminance mask. Guard it rather than trust it.
    for (x0, x1) in runs:
        idx = np.flatnonzero(ink[TOP_SEARCH, x0:x1])
        assert idx.size == 0 or max(r[-1] - r[0] + 1 for r in _runs(list(idx))) \
            < 0.72 * (x1 - x0), "TOP_SEARCH starts inside a full-width rule"

    rows, boxes = [], []
    for (x0, x1), lab in zip(runs, LABELS):
        # Corrected 2026-10-30, review finding F-018. The old rule took the
        # first row where >=85% of a w//8-trimmed window was coloured; three
        # overlapping replicate dots fill >=85% of THAT while covering only
        # 55-70% of the full bar, so it fired on the dots on 15 of the 60 bars
        # across Figs 2c and 3c, by up to +28% in value, and inflated the
        # cross-library offset from a true ~1.01 to 1.1116.
        ytop, ybot, fused = bar_top(ink, x0, x1, TOP_SEARCH, int(BASE_Y))
        assert ytop is not None, f"no bar top for {lab}"
        assert ytop < int(BASE_Y) - 26, f"bar top for {lab} overlaps the stroke band"
        value = 10 ** (4.0 + (Y_1E4 - ytop) / PX_PER_DECADE)
        rows.append({"lipid_name": lab, "is_benchmark": int(lab in BENCH),
                     "RLU": round(float(value), 1), "x0": x0, "x1": x1,
                     "y_top": ytop, "rule_thickness_px": ybot - ytop + 1,
                     "frac_used": fused, "source_figure": "Fig 2c",
                     "extraction_method": "bar_chart_log_axis"})
        boxes.append((x0, ytop, x1))

    with open("07_fig2c_bars.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    # how much of the heatmap's censoring does this lift?
    heat = {r["lipid_name"]: r for r in csv.DictReader(open("01_fig2b_320.csv"))}
    sat_cells = {k for k, v in heat.items() if v["saturated"] == "1"}
    rescued = sorted(sat_cells & {r["lipid_name"] for r in rows})
    still = sorted(sat_cells - set(rescued))

    qc = {
        "paper": "Liu et al., Adv Healthc Mater 2026, Fig 2c",
        "bars": len(rows), "lipids": len(rows) - len(BENCH),
        "calibration": {"y_1e4": Y_1E4, "px_per_decade": PX_PER_DECADE,
                        "anchors": "y-axis TICK centroids 2342.78/2749.54/3156.75 (F-017)"},
        "benchmark_values": {r["lipid_name"]: r["RLU"] for r in rows
                             if r["is_benchmark"]},
        "lipid_range": [min(r["RLU"] for r in rows if not r["is_benchmark"]),
                        max(r["RLU"] for r in rows if not r["is_benchmark"])],
        "saturated_cells_in_fig2b": len(sat_cells),
        "rescued_by_fig2c": len(rescued),
        "still_censored": still,
        "set_agreement_with_fig2b": {
            "plotted": len(rows) - len(BENCH),
            "also_above_2e4_in_heatmap": sum(
                1 for r in rows if not r["is_benchmark"]
                and float(heat[r["lipid_name"]]["RLU"]) > 2e4),
            "exception": "T3A-L8-C14 is plotted in 2c but its 2b cell reads "
                         "7.9e3 and is visibly cream; the two figures disagree",
        },
    }
    json.dump(qc, open("08_fig2c_qc.json", "w"), indent=1)
    overlay(img, boxes)

    print(f"bars {len(rows)}  lipids {len(rows)-len(BENCH)}  benchmarks {len(BENCH)}")
    print(f"  lipid range {qc['lipid_range'][0]:.3g} .. {qc['lipid_range'][1]:.3g} RLU")
    print("  benchmarks: " + ", ".join(f"{k} {v:.3g}"
                                       for k, v in qc["benchmark_values"].items()))
    print(f"  saturated cells in 2b: {len(sat_cells)}, rescued here: {len(rescued)}, "
          f"still censored: {len(still)} {still}")
    print(f"  set agreement with 2b: "
          f"{qc['set_agreement_with_fig2b']['also_above_2e4_in_heatmap']}"
          f"/{qc['set_agreement_with_fig2b']['plotted']}")


def overlay(img, boxes):
    y0 = TOP_SEARCH - 30
    crop = Image.fromarray(img.astype("uint8")[y0:int(BASE_Y) + 30, 450:]).copy()
    d = ImageDraw.Draw(crop)
    for x0, yt, x1 in boxes:
        d.rectangle([x0 - 450, yt - y0, x1 - 450, int(BASE_Y) - y0],
                    outline=(220, 0, 0), width=3)
    d.text((8, 8), "red = measured bar extent; top edge is the value read",
           fill=(0, 0, 0))
    crop.save(f"{REVIEW}/fig2c_overlay.png")


if __name__ == "__main__":
    main()
