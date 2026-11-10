#!/usr/bin/env python3
"""peer_39912611 stage 05 -- digitise Figure 3C/D/E (in vivo avg radiance)
and Figure 3F/G/H (the log10 organ ratios), which exist only to check C/D/E.

THE REDUNDANCY THE AUTHORS DREW. Panels F, G and H are log10(lungs/liver),
log10(lungs/spleen) and log10(liver/spleen) of the very bars in C, D and E. So

    log10(C_i) - log10(D_i) = F_i      for every lipid i

is a constraint this curation did not invent, is external to our calibration,
and was never printed as a number. It is worth more than any printed statistic
because it is insensitive to everything our pipeline and the authors' pipeline
share. It also SETTLES a label this raster cannot resolve: panel C's lowest
axis label is 8 px of JPEG and reads equally as 10^5 or 10^6, while its four
major ticks are uniformly spaced -- so the ladder is 10^8..10^5 at 47.5 px per
decade, or 10^8..10^6 at 71 px with one tick unlabelled, and only the first is
geometrically possible. See 07_qc.json check F3-3.

GUARD AGAINST THE SIGNIFICANCE BRACKETS. Panels C, D and E carry horizontal
"ns"/"*"/"***" brackets that span several bars. A bracket is a full-width
horizontal rule inside a bar's window, which is exactly what bars.bar_top is
looking for. It is removed before detection by a 5-row vertical erosion: a bar
top has its own bar under it, a bracket has white.
"""
import csv
import json
import math
import os
import sys

import numpy as np
from PIL import Image

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "lib"))
import axis                                                   # noqa: E402
import bars                                                   # noqa: E402
from errors import DetectionError                             # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
IMG = os.path.join(HERE, "source", "figs", "nn4c18636_0003.jpg")

CATS = ["Untreated", "Lipid 14", "Lipid 32", "Lipid 33", "Lipid 34",
        "Lipid 35", "Lipid 36", "Lipid 37", "Lipid 38", "Lipid 39",
        "Lipid 40", "Lipid 41", "Lipid 42"]

# key, organ, axis-rule column, plot x-range, decade values top->bottom
LOG_PANELS = [
    ("C", "lungs", 33, 36, 250, [1e8, 1e7, 1e6, 1e5]),
    ("D", "liver", 267, 270, 482, [1e8, 1e7, 1e6]),
    ("E", "spleen", 497, 500, 684, [1e8, 1e7, 1e6]),
]
EROSION = 5
BASE_ROW = 761          # the x-axis rule, shared by panels C, D and E


def load():
    a = np.asarray(Image.open(IMG).convert("RGB")).astype(int)
    gray = a.mean(axis=2)
    sat = a.max(axis=2) - a.min(axis=2)
    ink = (sat > 40) | (gray < 150)
    return a, gray, ink


def erode_down(ink, n):
    """True only where ink continues for n rows. Kills a thin horizontal rule
    (a significance bracket) and keeps a bar top, which has a bar under it."""
    out = ink.copy()
    for k in range(1, n):
        out[:-k] &= ink[k:]
    return out


def support_frac(ink, x0, x1, row, base_row, min_cols=2):
    """Fraction of rows between `row` and the baseline where at least
    `min_cols` columns of the window hold ink.

    This is what separates a bar top from a significance bracket WITHOUT the
    vertical erosion, which is needed for the two cases erosion cannot serve:
    the white Untreated bar, whose interior is blank so only its two outline
    columns carry ink, and a bar four pixels tall. A bar -- filled or hollow --
    has its own verticals all the way down; a bracket drawn over a category
    with no bar has nothing under it.
    """
    sub = ink[int(row) + 2:int(base_row), x0:x1]
    if sub.size == 0:
        return 0.0
    return float((sub.sum(axis=1) >= min_cols).mean())


def major_ticks(gray, ax_col, y_lo, y_hi, max_len=9, thresh=200, min_len=5):
    """Rows whose tick reaches at least `min_len` px left of the axis rule.

    A GraphPad log axis draws eight MINOR ticks per decade. They are shorter.
    Separating them by run length is what makes the decade ladder uniform; a
    minor tick admitted here would fail calibrate's residual test, loudly.
    """
    out = []
    for y in range(y_lo, y_hi):
        n = 0
        for d in range(1, max_len + 1):
            if gray[y, ax_col - d] < thresh:
                n += 1
            else:
                break
        if n >= min_len:
            out.append((y, n))
    # merge contiguous rows into one darkness-weighted centre
    groups, cur = [], [out[0]] if out else []
    for y, n in out[1:]:
        if y - cur[-1][0] <= 1:
            cur.append((y, n))
        else:
            groups.append(cur)
            cur = [(y, n)]
    if cur:
        groups.append(cur)
    return [float(np.mean([y for y, _ in g])) for g in groups]


def category_ticks(gray, x0, x1, base_row, thresh=170):
    """Category centres, from the x-axis TICK MARKS below the baseline.

    The bars in these panels touch each other -- a run of ink just above the
    baseline is one block nine categories wide, so the bars cannot delimit
    themselves. The axis ticks can, and they are the right furniture anyway:
    they are what the drawing program placed at each category centre.
    """
    band = gray[base_row + 2:base_row + 5, :] < thresh
    cnt = band.sum(axis=0)
    runs, s = [], None
    for x in range(x0, x1):
        if cnt[x] >= 3:
            if s is None:
                s = x
        else:
            if s is not None:
                runs.append((s, x))
                s = None
    if s is not None:
        runs.append((s, x1))
    return [(a0 + b0 - 1) / 2.0 for a0, b0 in runs]


def main():
    a, gray, ink = load()
    # Pad the rows at and below the shared baseline rule before eroding. A bar
    # two pixels tall has only four ink rows before the axis runs out, so a
    # 5-row erosion deletes it -- which is how lipid 33's spleen bar turned
    # into "no bar drawn". The pad says only what the figure already says:
    # the plot ends at the axis.
    ink_pad = ink.copy()
    ink_pad[BASE_ROW:BASE_ROW + EROSION + 2, :] = True
    eroded = erode_down(ink_pad, EROSION)
    qc = {"image": os.path.relpath(IMG, HERE), "erosion_rows": EROSION,
          "panels": {}}
    out_rows = []

    for key, organ, axcol, x0, x1, decades in LOG_PANELS:
        ticks = major_ticks(gray, axcol, 605, 775)
        if len(ticks) != len(decades):
            raise SystemExit(
                f"panel {key}: {len(ticks)} major ticks {ticks} but "
                f"{len(decades)} decade labels. Refusing.")
        cal = axis.calibrate(ticks, decades, kind="log10", axis="y",
                             loo_tol_px=axis.LOO_TOL_PX)
        loo = axis.loo_residuals(ticks)
        base_row = max(ticks) if key == "C" else 761.5

        cats = category_ticks(gray, x0, x1, int(round(base_row)))
        if len(cats) != 13:
            raise SystemExit(
                f"panel {key}: {len(cats)} category ticks, expected 13 "
                f"(Untreated + 12 lipids): {cats}")
        xcal = axis.calibrate(cats, list(range(13)), kind="linear", axis="x",
                              loo_tol_px=axis.LOO_TOL_PX)
        pitch = xcal.px_per_unit          # px per category
        # The window is the FULL category pitch, not pitch-2. A hollow bar
        # carries all of its ink on its two vertical outlines, and trimming
        # one column off each side put both of them outside the window -- so
        # the Untreated control looked like an empty category with a floating
        # rule, which is exactly what support_frac is built to reject.
        # Centres from the FITTED x calibration, not from each tick's own
        # centroid: a tick is 1 or 2 px wide, so its centroid carries +/-0.5 px
        # of quantisation, and the fit over all 13 averages that away.
        ctrs = [float(xcal.to_pixel(i)) for i in range(13)]
        # The bars TOUCH, so the boundary between category i and i+1 is the
        # midpoint of their centres and the two bars share an outline column.
        # The window is inset one column on each side so a flank test can
        # never see the neighbour. The cost is the hollow Untreated control,
        # whose ink is ENTIRELY on those two shared columns: it is reported as
        # a refusal rather than measured, and it is a control, not a library
        # member, so no corpus row depends on it. See DECISIONS D-4.
        half = pitch / 2.0
        slots = [(int(round(c - half)) + 1, int(round(c + half)))
                 for c in ctrs]
        idx = list(range(13))
        tops, fails, how = {}, {}, {}
        ytop = int(round(min(ticks))) - 4
        ybot = int(round(base_row))          # rows up to the baseline rule
        for (sx, ex), i in zip(slots, idx):
            got = None
            # Route 1: the eroded mask, bars.bar_top's own 3-row rule. This is
            # the strict route and it serves every bar taller than 5 px.
            try:
                got = bars.bar_top(eroded, sx, ex, ytop, ybot,
                                   baseline_row=None)
                how[i] = "eroded_nrow3"
            except DetectionError as e:
                first = str(e).split("\n")[0]
            if got is None:
                # Routes 2-4: raw ink with a shorter run, because two kinds of
                # bar cannot satisfy a 3-row rule at all -- a HOLLOW bar, whose
                # top rule is one pixel and whose interior is white (the
                # Untreated control), and a bar 2 px tall. Each of these must
                # then pass support_frac: a bar, filled or hollow, has its own
                # verticals all the way down to the baseline; a significance
                # bracket drawn over an empty category has nothing under it.
                for nrow in (3, 2, 1):
                    try:
                        t = bars.bar_top(ink, sx, ex, ytop, ybot, nrow=nrow,
                                         baseline_row=None)
                    except DetectionError as e2:
                        last = str(e2).split("\n")[0]
                        continue
                    sf = support_frac(ink, sx, ex, t["row"], base_row)
                    if sf < 0.90:
                        last = (f"raw nrow={nrow}: top at row {t['row']} has "
                                f"support {sf:.3f} < 0.90 -- nothing under "
                                "it, so it is a significance bracket.")
                        continue
                    t["support_frac"] = round(sf, 4)
                    got = t
                    how[i] = f"raw_nrow{nrow}_support{sf:.3f}"
                    break
            if got is None:
                # "No bar drawn" is judged by SUPPORT, not by whether the
                # window holds ink: panel D's significance brackets run over
                # the empty categories 34/35/36/41/42, so an ink test says
                # they are occupied. Support says what is under the ink.
                # "Is a bar drawn here at all?" is answered at the BASELINE,
                # where every bar is anchored however short it is, and with
                # min_cols=1 so a hollow bar counts on its outline alone.
                near = ink[ybot - 8:ybot, sx:ex].sum(axis=1) >= 1
                cov = float(near.mean())
                if cov < 0.9 and i != 0:
                    # No bar was drawn in this category at all. That is a
                    # statement about the axis floor, not a missing datum:
                    # Figure 3F/G/H plot ratios for these same lipids, so the
                    # value exists and is simply below what panel D can show.
                    pass
                else:
                    fails[i] = first + " || " + last
            else:
                tops[i] = got

        qc["panels"][key] = {
            "organ": organ,
            "major_ticks_px": [round(t, 2) for t in ticks],
            "tick_values": decades,
            "calibration": cal.qc(),
            "loo_residuals_px": [round(v, 3) for v in loo],
            "worst_loo_px": round(max(loo), 3),
            "loo_tol_px": axis.LOO_TOL_PX,
            "baseline_row_px": round(float(base_row), 2),
            "floor_value_at_baseline": float(f"{cal.to_value(base_row):.4g}"),
            "category_ticks_px": [round(c, 2) for c in cats],
            "x_calibration": xcal.qc(),
            "x_worst_loo_px": round(max(axis.loo_residuals(cats)), 3),
            "bar_pitch_px": round(pitch, 3),
            "bar_window_px": [int(e - s0) for s0, e in slots],
            "bar_indices_measured": sorted(tops),
            "bar_indices_no_bar_drawn": [k for k in range(13)
                                         if k not in tops and k not in fails],
            "n_tops": len(tops),
            "detection_route": how,
            "failures": fails,
        }

        for i in range(13):
            if i in tops:
                row = tops[i]["row"]
                v = cal.to_value(row)
                st, bound = "measured", ""
            elif i in fails:
                row, v, st, bound = None, None, "REFUSED", ""
            else:
                # no bar drawn at all: the value is at or below the floor the
                # axis can show. A bound, not a measurement, not a zero.
                row, v = None, float(f"{cal.to_value(base_row):.4g}")
                st, bound = "below_axis_floor", "upper"
            out_rows.append({
                "panel": f"Fig 3{key}", "organ": organ,
                "lipid_name": CATS[i], "cat_index": i,
                "top_row_px": None if row is None else round(float(row), 2),
                "avg_radiance_p_s_cm2_sr": v,
                "log10": None if v is None else round(float(np.log10(v)), 4),
                "status": st, "value_is_bound": bound,
                "detail": fails.get(i, ""),
            })

    with open(os.path.join(HERE, "05_fig3_values.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out_rows[0].keys()))
        w.writeheader()
        w.writerows(out_rows)
    with open(os.path.join(HERE, "05_fig3_qc.json"), "w") as fh:
        json.dump(qc, fh, indent=1)
    for k, v in qc["panels"].items():
        print(k, v["organ"], "ticks", v["major_ticks_px"], "resid",
              v["calibration"]["tick_fit_residual_px"], "pxdec",
              round(v["calibration"]["px_per_unit"], 3),
              "floor", v["floor_value_at_baseline"],
              "measured", v["bar_indices_measured"],
              "no_bar", v["bar_indices_no_bar_drawn"],
              "fails", list(v["failures"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
