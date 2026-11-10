#!/usr/bin/env python3
"""ONE audited bar-top detector, immune to overlaid replicate markers.

WHY THIS FILE EXISTS. Error class 32 -- a bar-top detector reading an overlaid
replicate dot, a whisker cap or a scatter marker as the bar top -- appeared
independently in two curations written days apart:

  F-018  curation/ql_2026/digitize_fig3c.py    replicate dots read as bar tops
                                               on 15 of 60 bars, by up to +28%
                                               in value. It also contaminated
                                               three of the four benchmark
                                               anchors that set the stored
                                               cross-library offset of 1.1116,
                                               so that offset was largely
                                               measuring our own error; the
                                               corrected 1.0071 +/- 0.0247 is
                                               not distinguishable from 1.0.
  F-020  curation/xd_2026_lipidai/              the same class, 12 of 12 tops,
           stage04_fig5c.py                     mean +0.051 decades.

THE GEOMETRY THAT SEPARATES THEM. A bar's top rule spans the FULL bar width. A
replicate dot cluster, a whisker and a mean marker all sit in a narrow central
strip. So the bar top is the topmost row whose coloured coverage across the full
bar window exceeds a high threshold -- not the topmost coloured pixel, which is
the top of the error bar.

F-018's lesson is that "high" has to be measured on the FULL detected width.
digitize_fig3c.py required >= 0.72 coverage of an inner-trimmed window, and
three overlapping dots can exceed 85% of a trimmed window while covering only
60-70% of the full bar. Hence the defaults here: >= `frac` coverage of the full
window, sustained over `nrow` contiguous rows, PLUS left/right flank agreement.

THE FLANK TEST is the part no single-threshold rule can replace. Measure the
first covered row independently in the left and right thirds of the bar. A real
top rule is present in both at the same row; a centred marker is in neither.
Disagreement beyond `flank_tol_px` means the detection is not a top rule, and
this module raises instead of returning the better-looking of the two.

THIS MODULE NEVER RETURNS A FALLBACK. A miss raises DetectionError. Error
class 28 is a detector miss falling through to a literal, and F-033 is what it
costs.
"""
import numpy as np

from errors import DetectionError

FRAC = 0.90           # coverage of the FULL bar width required for a top rule
NROW = 3              # contiguous rows that must hold it
FLANK_TOL_PX = 2.0    # left/right first-covered-row agreement
FLANK_FRAC = 0.33     # width of each flank window, as a fraction of the bar


def _first_covered_row(ink, x0, x1, y_top, y_bot, frac, nrow):
    """Topmost row index with >= frac coverage over [x0,x1) for nrow rows."""
    if x1 - x0 < 2:
        return None
    sub = ink[y_top:y_bot, x0:x1]
    if sub.size == 0:
        return None
    cov = sub.mean(axis=1)
    run = 0
    for i, c in enumerate(cov):
        if c >= frac:
            run += 1
            if run >= nrow:
                return y_top + i - nrow + 1
        else:
            run = 0
    return None


def bar_top(ink, x0, x1, y_top, y_bot, frac=FRAC, nrow=NROW,
            flank_tol_px=FLANK_TOL_PX, flank_frac=FLANK_FRAC,
            baseline_row=None):
    """Topmost row of a bar's top rule over the window [x0, x1).

    ink           2-D boolean array, True where the bar's colour is present.
                  Build it from the bar's OWN colour, not from darkness: in
                  QL_2026's Fig 3c the bars, whiskers and dots are all one blue
                  while the axis and text are black, and a darkness mask makes
                  the axis rule itself the topmost "bar".
    x0, x1        the bar's full detected width. Do NOT pre-trim it -- an
                  inner-trimmed window is exactly what let F-018 through.
    y_top, y_bot  the row range to search, top-down.
    baseline_row  if given, the axis baseline row. A detection landing within
                  nrow of it means the bar has no height, which is a different
                  statement from a small value and is raised, not returned.

    Returns dict(row, coverage, left_row, right_row, flank_delta_px, frac_used).
    Raises DetectionError on a miss or on flank disagreement.
    """
    ink = np.asarray(ink, dtype=bool)
    if x1 <= x0:
        raise DetectionError(f"bar_top: empty window x0={x0} x1={x1}.")
    width = x1 - x0
    if width < 6:
        raise DetectionError(
            f"bar_top: bar window is {width} px wide. The flank test needs a "
            "bar wide enough for a left and a right third; below ~6 px the "
            "full-width coverage rule degenerates into the topmost-pixel rule "
            "that produced F-018 and F-020. Render the page at higher dpi.")

    # 1. the full-width top rule.
    row = _first_covered_row(ink, x0, x1, y_top, y_bot, frac, nrow)
    frac_used = frac
    if row is None:
        raise DetectionError(
            f"bar_top: no row in [{y_top},{y_bot}) reaches {frac:.2f} coverage "
            f"of the full {width}px bar for {nrow} contiguous rows.\n"
            "This is a REFUSAL, not a zero. A hollow bar, a bar occluded by a "
            "legend, or a mask built from the wrong colour all land here. Do "
            "not lower `frac` to make it pass: 0.72 of an inner-trimmed window "
            "is what let three overlapping replicate dots be read as a bar top "
            "on 15 of 60 bars (F-018). Fix the mask, or record the bar as "
            "unreadable in the folder's DECISIONS.md.")

    # 2. the flank test. A centred marker is in neither flank.
    fw = max(2, int(round(width * flank_frac)))
    left = _first_covered_row(ink, x0, x0 + fw, y_top, y_bot, frac, nrow)
    right = _first_covered_row(ink, x1 - fw, x1, y_top, y_bot, frac, nrow)
    if left is None or right is None:
        raise DetectionError(
            f"bar_top: full-width rule fired at row {row} but the "
            f"{'left' if left is None else 'right'} {fw}px flank has no "
            "covered row at all. A bar's top rule spans the full width, so it "
            "is present in both flanks; something present only in the middle "
            "is a replicate dot cluster, a whisker cap or a mean marker. That "
            "is error class 32 (F-018, F-020).")
    delta = abs(float(left) - float(right))
    if delta > flank_tol_px:
        raise DetectionError(
            f"bar_top: left flank first-covered row {left}, right flank "
            f"{right}, delta {delta:.1f} px > {flank_tol_px} px tolerance.\n"
            "A real top rule is at the same row in both flanks. Disagreement "
            "this large means the detection is tracking something that is not "
            "a horizontal rule -- a sloped edge, an overlapping neighbour bar, "
            "or a marker. Refusing rather than returning the better-looking "
            "flank.")

    # 3. a detection sitting on the baseline is not a measurement.
    if baseline_row is not None and abs(row - baseline_row) <= nrow:
        raise DetectionError(
            f"bar_top: detected row {row} is within {nrow} px of the baseline "
            f"{baseline_row}. The first scanned row is already the axis rule, "
            "so this is the baseline being read as a bar of zero height, not a "
            "measurement. Narrow y_top/y_bot to exclude the axis.")

    cov = float(ink[row:row + nrow, x0:x1].mean())
    return {"row": int(row), "coverage": round(cov, 4),
            "left_row": int(left), "right_row": int(right),
            "flank_delta_px": round(delta, 2), "frac_used": frac_used,
            "width_px": int(width)}


def bar_tops(ink, windows, y_top, y_bot, **kw):
    """bar_top over many windows. Collects refusals instead of stopping.

    Returns (results, failures): `results` maps the window index to its dict,
    `failures` maps the index to the refusal message. A caller that wants to
    treat a failure as a value has to write that itself.
    """
    results, failures = {}, {}
    for i, (x0, x1) in enumerate(windows):
        try:
            results[i] = bar_top(ink, x0, x1, y_top, y_bot, **kw)
        except DetectionError as e:
            failures[i] = str(e)
    return results, failures
