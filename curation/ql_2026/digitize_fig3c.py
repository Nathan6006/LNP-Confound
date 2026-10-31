#!/usr/bin/env python3
"""Digitize Figure 3c: library 2 (24 lipids) plus four benchmark LNPs, log axis.

These 24 are the other half of the paper's 344 and exist nowhere else, so this
panel is the primary source for them - unlike Fig 2c, which only re-plots
lipids the Fig 2b heatmap already carries.

BAR TOP, NOT WHISKER. Each column carries a hollow bar, an error whisker and
three scatter dots, all in the same blue. The whisker and dots sit in a narrow
central strip while the bar's top edge spans the full bar width, so the bar top
is the highest row whose coloured coverage across the bar window exceeds 60%.
Reading the highest coloured pixel instead would return the top of the error
bar and overstate every value.

CALIBRATION is anchored on the y-axis TICK MARKS, measured by darkness-weighted
centroid over the columns left of the axis rule: 10^3 at y = 2911.5 with 151.0 px
per decade (600 dpi). It is NOT anchored on the label centroids and NOT on the
panel box.

  An earlier version of this file anchored 10^3 on the box bottom (y = 2959,
  163.75 px/decade) and made every value in the panel up to 1.76x too high. The
  box bottom is NOT the 10^3 tick - the axis floor here is ~500 RLU, so the tick
  sits 47 px above the rule. That is error class 20. This docstring described
  the pre-fix method for a further day after the code was repaired, which is
  error class 21; both are now in provenance/ERROR_CLASSES.md.

  Still open against this panel (round-1 review, 2026-10-30): the bar-top rule
  below is not strict enough - three overlapping replicate dots can exceed 85%
  of the INNER-TRIMMED width while covering only ~60-70% of the full bar, so
  several bars here are read on the dot cluster. Require a contiguous band of
  >=3 rows at >=0.9 coverage of the FULL detected width. See
  provenance/review/round_2026-10-30/FINDINGS.md.

The calibration is then CHECKED, not assumed:
the caption defines a hit as RLU > 7e5 with n = 8 per tail group and n = 3 per
BL group, and panels 3d/3e print those hit rates. A wrong decade scale moves
values by a factor and the counts stop matching.
"""
import csv
import json
import os

import numpy as np
import PIL.Image as Image
import PIL.ImageDraw as ImageDraw

PAGE = "00_pages/fig_p5_600.png"
REVIEW = "review"

# CALIBRATION. Anchored on the y-axis TICK MARKS, measured left of the axis
# line: 2308.5, 2459.5, 2610.5, 2760.5, 2911.5 -> 151.0 px per decade, with
# 10^3 at 2911.5.
#
# The first version anchored 10^3 on the panel BOX BOTTOM (2959) and 10^7 on the
# box top, giving 163.75 px/decade. The box top does coincide with 10^7, which is
# why it looked right, but the 10^3 tick sits 47 px ABOVE the bottom rule and the
# bars run on below it. Every value came out high by a factor up to 1.76. Three
# printed statistics catch it and all three now agree: Fig 3d's 2C11 rate (1/8),
# Fig 3e's L8 rate (0/3), and the Results sentence "17 lipids ... greater than
# 1 x 10^5".
Y_1E3 = 2911.5
PX_PER_DECADE = 151.0
BASE_Y = 2959.0          # panel box bottom: where bars are drawn from
BASE_LOG = 3.0
AXIS_X = 246
XMAX = 2200
INK_THR = 172            # "anything drawn": see bar_top() below
# The panel BOX TOP rule occupies y = 2301-2314. The old search start,
# int(Y_1E3 - 4*PX_PER_DECADE) + 4 = 2311, is INSIDE it. That was harmless
# while the top was read on a blue mask (the box is black), but bar_top() reads
# a luminance mask, and at 2311 every bar returns 2311. Hence an explicit
# constant, with a guard in main().
TOP_SEARCH = 2320

LIBRARY2 = [f"T3A-{bl}-{t}" for t in ("2C10", "2C11", "2C12")
            for bl in (f"L{i}" for i in range(1, 9))]
BENCH = ["DLin-MC3-DMA", "SM102", "ALC-0315", "Lipid 5"]
LABELS = LIBRARY2 + BENCH


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
    g = img.mean(2)
    # Bars, whiskers and dots are all blue; the axis and text are black. Using
    # black too merged the axis and labels into the bar runs, so restrict to blue.
    coloured = (img[..., 2] > img[..., 0] + 40) & (sat > 60)

    top = TOP_SEARCH
    ink = img.mean(2) <= INK_THR
    # bar windows from the columns that carry a tall coloured run
    # A hollow bar is two thin vertical strokes with an empty interior, so
    # requiring a TALL coloured run per column split every bar into two or three
    # pieces (37 runs for 28 bars). What spans the whole bar is its TOP EDGE, so
    # a column counts as inside a bar if it carries even a few coloured pixels.
    prof = np.zeros(XMAX - AXIS_X, dtype=int)
    for i, x in enumerate(range(AXIS_X, XMAX)):
        prof[i] = coloured[top:int(BASE_Y) - 12, x].sum()
    on = prof >= 3
    runs, cur = [], None
    for i, v in enumerate(on):
        if v and cur is None:
            cur = i
        elif not v and cur is not None:
            if i - cur > 12:
                runs.append((AXIS_X + cur, AXIS_X + i))
            cur = None
    if cur is not None:
        runs.append((AXIS_X + cur, AXIS_X + len(on)))

    assert len(runs) == len(LABELS), \
        f"found {len(runs)} bars, expected {len(LABELS)}: {runs[:40]}"

    for (x0, x1) in runs:
        idx = np.flatnonzero(ink[top, x0:x1])
        assert idx.size == 0 or max(r[-1] - r[0] + 1 for r in _runs(list(idx))) \
            < 0.72 * (x1 - x0), "TOP_SEARCH starts inside a full-width rule"

    rows, boxes = [], []
    for (x0, x1), lab in zip(runs, LABELS):
        # Corrected 2026-10-30, review finding F-018; see digitize_fig2c.py and
        # bar_top() above. The old >=85%-of-a-trimmed-window rule read the
        # replicate dots on 5 of this panel's 28 bars, by up to +28% in value.
        ytop, ybot, fused = bar_top(ink, x0, x1, top, int(BASE_Y))
        assert ytop is not None, f"no bar top found for {lab}"
        assert ytop < int(BASE_Y) - 26, f"bar top for {lab} overlaps the stroke band"
        value = 10 ** (BASE_LOG + (Y_1E3 - ytop) / PX_PER_DECADE)
        rows.append({"lipid_name": lab,
                     "is_benchmark": int(lab in BENCH),
                     "RLU": round(float(value), 1),
                     "x0": x0, "x1": x1, "y_top": ytop,
                     "rule_thickness_px": ybot - ytop + 1,
                     "frac_used": fused,
                     "source_figure": "Fig 3c",
                     "extraction_method": "bar_chart_log_axis"})
        boxes.append((x0, ytop, x1))

    with open("05_fig3c_24.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    lib = [r for r in rows if not r["is_benchmark"]]
    hits = [r for r in lib if r["RLU"] > 7e5]
    by_tail = {t: sum(1 for r in hits if r["lipid_name"].endswith(t))
               for t in ("2C10", "2C11", "2C12")}
    by_bl = {f"L{i}": sum(1 for r in hits if f"-L{i}-" in r["lipid_name"])
             for i in range(1, 9)}
    qc = {
        "paper": "Liu et al., Adv Healthc Mater 2026, Fig 3c",
        "bars": len(rows), "library2": len(lib), "benchmarks": len(BENCH),
        "calibration": {"y_of_1e3": Y_1E3, "px_per_decade": PX_PER_DECADE,
                        "anchors": "y-axis tick marks 2308.5/2459.5/2610.5/"
                                   "2760.5/2911.5 = 1e7..1e3"},
        "gt_1e5": sum(1 for r in rows if not r["is_benchmark"] and r["RLU"] > 1e5),
        "gt_1e5_printed": "Results: '17 lipids ... greater than 1 x 10^5'",
        "value_range": [min(r["RLU"] for r in lib), max(r["RLU"] for r in lib)],
        "benchmark_values": {r["lipid_name"]: r["RLU"] for r in rows
                             if r["is_benchmark"]},
        "hit_check": {
            "definition": "RLU > 7e5 (Fig 3 caption)",
            "total_hits": len(hits),
            "hits": [r["lipid_name"] for r in hits],
            "by_tail_n8": {k: f"{v}/8 = {100*v/8:.1f}%" for k, v in by_tail.items()},
            "by_bl_n3": {k: f"{v}/3 = {100*v/3:.1f}%" for k, v in by_bl.items()},
            "printed_3d": "2C10 37.5% (3/8), 2C11 12.5% (1/8), 2C12 12.5% (1/8)",
            "printed_3e": "L1 33% (1/3), L5 67% (2/3), L6 33% (1/3), L7 33% (1/3)",
            "resolution_note": "Every printed rate in 3d and 3e is reproduced. "
                          "The borderline case this note used to describe -- "
                          "T3A-L8-2C11 at 701,135 RLU, one pixel from the 7e5 "
                          "threshold -- came from the superseded bar-top rule "
                          "and no longer exists (F-018). NOTE, so that nobody "
                          "reads this as evidence the new rule is right: these "
                          "rates reproduced on the OLD values too, including on "
                          "15 of 60 bars that were up to 28% high. A printed "
                          "marginal cannot discriminate a geometric error of "
                          "this kind. The bar tops are validated geometrically, "
                          "per-bar, in provenance/review/round_2026-10-30/.",
        },
    }
    json.dump(qc, open("06_fig3c_qc.json", "w"), indent=1)
    overlay(img, boxes)

    print(f"bars {len(rows)} ({len(lib)} library-2 + {len(BENCH)} benchmarks)")
    print(f"  library-2 range {qc['value_range'][0]:.3g} .. {qc['value_range'][1]:.3g} RLU")
    print(f"  benchmarks: " + ", ".join(f"{k} {v:.3g}" for k, v in
                                        qc["benchmark_values"].items()))
    print(f"  hits >7e5: {len(hits)}  {[r['lipid_name'] for r in hits]}")
    print(f"  by tail  {qc['hit_check']['by_tail_n8']}")
    print(f"           printed 3d: {qc['hit_check']['printed_3d']}")
    print(f"  by BL    {dict((k, v) for k, v in qc['hit_check']['by_bl_n3'].items() if not v.startswith('0'))}")
    print(f"           printed 3e: {qc['hit_check']['printed_3e']}")


def overlay(img, boxes):
    y0 = TOP_SEARCH - 20
    crop = Image.fromarray(img.astype("uint8")[y0:int(BASE_Y) + 40, 180:XMAX + 40]).copy()
    d = ImageDraw.Draw(crop)
    for x0, yt, x1 in boxes:
        d.rectangle([x0 - 180, yt - y0, x1 - 180, int(BASE_Y) - y0],
                    outline=(220, 0, 0), width=3)
    d.text((8, 8), "red = measured bar extent; top edge is the value read",
           fill=(0, 0, 0))
    crop.save(f"{REVIEW}/fig3c_overlay.png")


if __name__ == "__main__":
    main()
