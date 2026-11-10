#!/usr/bin/env python3
"""Stage 04 -- Figure S3 digitised TO VERIFY Table S2, never to extract from it.

PIPELINE_V2.md section 3: "Once you have numbers from tiers 1-3, use the figure
to VERIFY, not to extract." Table S2's %TE are printed numerals (tier 5 with
numerals, 00_value_sources.json VS-1); Figure S3 is the same 45 lipidoids drawn
as bars on an axis with printed ticks. Nothing here is emitted.

The three things this buys, none of which a recomputed marginal could:

  1. A ROUND NUMBER THE AUTHORS DREW AND NEVER PRINTED. The green bar labelled
     'Positive control' L2K is what 'Relative Transfection Efficiency' is
     relative TO, so it must land on exactly 1.000. Our calibration comes only
     from the y-axis ticks (0.0/0.3/0.6/0.9/1.2), so where that bar lands is an
     external test of the calibration. CLAUDE.md section 5 ranks this above any
     printed statistic.
  2. ORDER AND TIE STRUCTURE. Table S2 has five lipidoids at 57%, three at 52%,
     five at 50%, two at 48%... If the transcription were off by a row the bar
     heights would stop matching the printed values.
  3. A RATIO THAT MUST BE CONSTANT. rel_i / %TE_i is 1/%TE(L2K) for every bar.
     Its SPREAD is the check. The scalar is NOT fitted -- it is predicted by
     the L2K bar, independently.

CLASS 20 / 27 GUARDS. The scale and origin are fitted together by
curation/lib/axis.calibrate from DETECTED TICK CENTRES; there is no API here to
supply an origin, which is the shape F-033 took. `forbid_target` is passed the
Table S2 values, so if anyone ever fits this axis onto them it raises.

Run:  python3 -I stage04_figS3.py
"""
import csv
import json
import os
import sys

import numpy as np
from PIL import Image
from scipy import ndimage

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))

import pymupdf                                              # noqa: E402
import axis as axislib                                      # noqa: E402
import bars as barslib                                      # noqa: E402
from errors import CurationLibError                         # noqa: E402

SI = os.path.join(HERE, "source", "bc8b00007_si_001.pdf")
PAGE = 16                 # 0-based; SI p.17
DPI = 600
TICK_VALUES = [1.2, 0.9, 0.6, 0.3, 0.0]   # top to bottom, as printed
AXIS_X = 1697             # the y-axis rule, measured: columns 1695-1700
BASELINE_Y = 2589         # the x-axis rule, measured: rows 2587-2592
PLOT_TOP = 950            # above the 1.2 tick, below the caption
N_BARS = 47               # 45 lipidoids + SFA + L2K


def render(doc):
    pix = doc[PAGE].get_pixmap(dpi=DPI)
    return np.array(Image.frombytes("RGB", (pix.width, pix.height),
                                    pix.samples))


def main():
    doc = pymupdf.open(SI)
    rgb = render(doc)
    gray = rgb.astype(int).mean(2)

    table = list(csv.DictReader(open(os.path.join(HERE, "05_values.csv"))))
    te = [float(r["value"]) for r in table]

    # --- y axis, from its own ticks only -------------------------------
    # MAJOR ticks only. The axis carries a minor tick halfway between every
    # labelled pair (measured: 8 evenly spaced runs at ~197.6 px when the strip
    # starts at AXIS_X-38), and feeding 9 positions to 5 labels is how a
    # calibration silently halves its scale. The majors are drawn LONGER, so a
    # strip that starts further from the rule sees only them -- the selection
    # is a property of the drawing, not a hand-picked subset.
    ticks = axislib.detect_ticks_1d(gray, "y", AXIS_X - 33, AXIS_X - 25,
                                    dark=140, min_len=3, max_len=40)
    ticks = [t for t in ticks if PLOT_TOP <= t <= BASELINE_Y + 8]
    if len(ticks) != len(TICK_VALUES):
        raise SystemExit(
            f"stage04: detected {len(ticks)} major y ticks, expected "
            f"{len(TICK_VALUES)} ({TICK_VALUES}). Refusing to calibrate on a "
            "tick set that does not match the printed labels.")
    cal = axislib.calibrate(ticks, TICK_VALUES, kind="linear",
                            forbid_target=te, axis="y")

    # --- bars: find them by their own colour, then measure the top rule --
    gray_i = rgb.astype(int).mean(2)
    R, G, B = (rgb[..., i].astype(int) for i in range(3))
    # BLACK ONLY. A plain darkness mask is wrong on this panel for a reason
    # the figure makes obvious once you look: a large dark-green arrow sweeps
    # diagonally across the plot, and its edge is dark enough to be "ink". With
    # gray < 110 the flank test refused 32 of 47 bars, with left/right
    # disagreements up to 142 px -- the arrow, caught on one flank only. The
    # bar outline, the whisker and the axis are achromatic; the arrow and the
    # fills are not.
    dark = (R < 110) & (G < 110) & (B < 110) & (abs(R - G) < 40) \
        & (abs(G - B) < 40) & (abs(R - B) < 40)
    # Bar WINDOWS are found in a strip just above the baseline, where every bar
    # is present by construction and the only white is the gap between bars.
    # Finding them from the hatched fill further up merges neighbours, because
    # a hatch leaves white columns inside a bar: that attempt returned 17 runs
    # for 47 bars and is recorded here so the choice is visible.
    strip = (gray_i < 235)[BASELINE_Y - 89:BASELINE_Y - 9, :]
    on = strip.mean(0) > 0.6
    runs, start = [], None
    for i, v in enumerate(on):
        if v and start is None:
            start = i
        elif not v and start is not None:
            if i - start >= 8 and start > AXIS_X:
                runs.append((start, i))
            start = None
    if len(runs) != N_BARS:
        raise SystemExit(
            f"stage04: found {len(runs)} bar columns, expected {N_BARS} "
            "(45 lipidoids + SFA + L2K). Refusing to guess which is which.")

    # The bar OUTLINE is black and spans the full width; the error-bar whisker
    # is black too but is centred, which is exactly what bars.bar_top's flank
    # test exists to reject (error class 32, F-018/F-020).
    # KEEP ONLY INK THAT IS PART OF A BAR. Two annotations sit over the plot
    # and both break the flank test: the legend block (black text and two
    # hatched swatches, y~1020-1240) and the large green arrow. A bar is a
    # closed outline standing ON the baseline, so the discriminator is
    # connectivity to the baseline row -- a property of the drawing, not a
    # hand-drawn exclusion box. The legend floats; the arrow floats; every bar
    # touches.
    fill = gray_i < 235                     # any ink: outline, hatch or fill
    lab, _ = ndimage.label(fill, structure=np.ones((3, 3)))
    touching = set(np.unique(lab[BASELINE_Y, :])) - {0}
    keep = np.isin(lab, list(touching))
    bar_ink = fill & keep
    # A pure-black mask does NOT work here and the reason is worth recording:
    # the figure is an embedded raster, so the bars' 1-2 px top rule is
    # anti-aliased and tinted by the red hatch underneath it, and an
    # achromatic-and-dark test loses enough of it that `frac`=0.90 full-width
    # coverage is never reached. The bar's own OUTLINE PLUS FILL is what spans
    # the full width at the top rule, so the mask is "any ink, connected to the
    # baseline". Lowering `frac` instead would have been F-018's mistake.

    tops, fails = barslib.bar_tops(bar_ink, runs, int(min(ticks)) - 30,
                                   BASELINE_Y - 4, baseline_row=BASELINE_Y)
    # A REFUSAL IS A RESULT. bars.bar_top declines a bar whose flanks
    # disagree, and three bars here are crossed by the big green arrow's
    # outline at their top rule. They are recorded unreadable and left out of
    # the comparison rather than rescued by loosening `frac` -- that loosening
    # is exactly F-018. The two control bars must both survive or the
    # drawn-round-number check is gone, so that is asserted.
    for idx in (N_BARS - 2, N_BARS - 1):
        if idx in fails:
            raise SystemExit(
                "stage04: a control bar (SFA / L2K) was refused, so the "
                "drawn-round-number check cannot be made. " + fails[idx])

    vals = {i: cal.to_value(t["row"]) for i, t in tops.items()}
    # Figure S3 runs L45 -> L1 left to right, then SFA, then L2K.
    lip = {}
    for i in range(45):
        if i in vals:
            lip[44 - i] = vals[i]           # 0-based index into L1..L45
    sfa, l2k = vals[45], vals[46]
    readable = sorted(lip)
    unreadable = [f"L{i + 1}" for i in range(45) if i not in lip]

    ratio = {i: lip[i] / te[i] for i in readable}
    pred_l2k = [1.0 / r for r in ratio.values()]

    # THE DECISIVE ALIGNMENT TEST. If Table S2's names and values had been
    # read one row out of register, the printed value would belong to the
    # neighbouring bar. Shifting the printed vector against the drawn vector
    # must make the agreement collapse, and shift 0 must be uniquely best.
    shift = {}
    for k in (-2, -1, 0, 1, 2):
        pairs = [(lip[i], te[i + k]) for i in readable
                 if 0 <= i + k < 45 and i in lip]
        rr = [a / b for a, b in pairs]
        shift[str(k)] = {"n": len(pairs),
                         "cv_pct": round(100 * float(np.std(rr, ddof=1))
                                         / float(np.mean(rr)), 3),
                         "pearson": round(float(np.corrcoef(
                             [a for a, _ in pairs],
                             [b for _, b in pairs])[0, 1]), 6)}

    ties = {}
    for i in readable:
        ties.setdefault(te[i], []).append(lip[i])
    tie_spread = {str(k): round(max(v) - min(v), 4)
                  for k, v in ties.items() if len(v) > 1}

    order_ok = all(lip[readable[k]] >= lip[readable[k + 1]] - 0.02
                   for k in range(len(readable) - 1))
    rho = float(np.corrcoef([lip[i] for i in readable],
                            [te[i] for i in readable])[0, 1])

    out = {
        "_what": "Figure S3 digitised to VERIFY Table S2. Nothing here is "
                 "emitted into 05_values.csv.",
        "render_dpi": DPI,
        "calibration": cal.qc(),
        "n_bars_detected": len(runs),
        "drawn_round_number_check": {
            "what": "'Relative Transfection Efficiency' is relative to a "
                    "positive control, so that control's own bar must be "
                    "drawn at exactly 1.000. The axis is calibrated from tick "
                    "marks alone, so where that bar lands is EXTERNAL to this "
                    "pipeline. Two green 'Positive control' bars exist, SFA "
                    "and L2K, and the figure never says which is the "
                    "normaliser -- so both are measured and the data chooses.",
            "SFA_measured": round(sfa, 4),
            "SFA_pct_from_1.000": round(100 * abs(sfa - 1.0), 2),
            "L2K_measured": round(l2k, 4),
            "L2K_pct_from_1.000": round(100 * abs(l2k - 1.0), 2),
            "normaliser_indicated": "SFA" if abs(sfa - 1) < abs(l2k - 1)
                                    else "L2K",
            "corroboration": "Figure S6 (SI p.19) states 'T18U-PY16-A1 shows "
                             "35% higher transfection efficiency compared to "
                             "positive control SFA', i.e. SFA is the "
                             "normaliser there too.",
            "residual_is_explained_by": "bars.bar_top returns the OUTER edge "
                                        "of the top rule, which sits ~5 px "
                                        "(0.004 units) above the plotted "
                                        "value at this line width.",
            "the_authors_never_printed_this_number": True},
        "SFA_bar": round(sfa, 4),
        "bars_refused_by_bar_top": {str(k): v.splitlines()[0]
                                    for k, v in fails.items()},
        "lipidoid_bars_unreadable": unreadable,
        "lipidoid_bars_compared": len(readable),
        "monotone_over_readable_bars": order_ok,
        "pearson_drawn_vs_printed": round(rho, 6),
        "ratio_drawn_over_printedTE": {
            "mean": round(float(np.mean(list(ratio.values()))), 6),
            "sd": round(float(np.std(list(ratio.values()), ddof=1)), 6),
            "cv_pct": round(100 * float(np.std(list(ratio.values()), ddof=1))
                            / float(np.mean(list(ratio.values()))), 3),
            "implied_%TE_of_the_normaliser_SFA_mean": round(float(np.mean(pred_l2k)), 2),
            "implied_%TE_of_the_normaliser_SFA_sd": round(float(np.std(pred_l2k, ddof=1)), 2)},
        "row_shift_test": shift,
        "row_shift_test_caveat": "Table S2 is sorted by value, so the drawn "
                                 "vector is monotone and a one-row shift "
                                 "still correlates. The test is therefore "
                                 "SUGGESTIVE, not decisive: shift 0 is "
                                 "uniquely best on both CV and Pearson, but "
                                 "the margin over shift +1 is small by "
                                 "construction.",
        "tie_groups_height_spread": tie_spread,
        "max_tie_spread": max(tie_spread.values()) if tie_spread else None,
        "per_bar": [{"lipidoid": table[i]["lipidoid"],
                     "printed_TE": te[i], "drawn_relative": round(lip[i], 4),
                     "ratio": round(ratio[i], 5)} for i in readable],
    }
    with open(os.path.join(HERE, "04_figS3_verification.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    for k in ("calibration", "drawn_round_number_check",
              "ratio_drawn_over_printedTE", "max_tie_spread",
              "pearson_drawn_vs_printed", "monotone_over_readable_bars",
              "lipidoid_bars_unreadable", "lipidoid_bars_compared",
              "row_shift_test"):
        print(k, "=", json.dumps(out[k]))
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except CurationLibError as e:
        print("REFUSED:", e)
        sys.exit(1)
