#!/usr/bin/env python3
"""Stage 05b -- Figure 3B: GFP mRNA transfection, 12 samples x 4 cell lines.

    python3 -I stage05b_fig3b.py

TIER 5. PIXELS. This is the only value source in the paper that is not vector,
and 00_value_sources.json records why tiers 1-4 all fail for it.

Everything geometric goes through curation/lib:
  * the y scale comes from axis.calibrate on the FOUR DETECTED TICKS, with no
    separate origin -- there is no API to supply one, which is the shape F-033
    took;
  * every tick's LEAVE-ONE-OUT residual is reported, because the max-residual
    test is weak against a single bad tick (an end tick has leverage 0.2, so
    four fifths of its error is absorbed into the fit);
  * every bar top comes from bars.bar_top, which requires >= 90% coverage of
    the FULL bar width over 3 contiguous rows AND left/right flank agreement,
    and which REFUSES rather than returning the topmost ink. The topmost ink in
    this panel is the error-bar whisker, and reading it is error class 32
    (F-018: 15 of 60 bars, up to +28%; F-020: 12 of 12 tops).

THE BAR GRID IS FITTED FROM AXIS FURNITURE, NOT FROM THE BARS. The 12 group
positions come from the x-axis tick marks below the baseline. Only the four
within-group colour offsets are measured off the bars, from the groups where
all four colours are unambiguously present, and their residuals are reported.
Doing it the other way round -- taking the grid from the bars -- would make the
short bars (87-, 90-OCholB, naked mRNA, untreated control) unlocatable exactly
where they matter.

THE MASK IS INK, NOT COLOUR. Each bar carries a thin black outline, so a mask
built from the fill colour alone starts ~2 px below the true top, a uniform
-0.4 percentage-point bias. An ink mask includes the outline. The whisker is
narrower than the bar and centred, so it cannot reach 90% full-width coverage
and cannot pass the flank test -- that is what bars.bar_top is for.
"""
import csv
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
sys.path.insert(0, HERE)

import numpy as np                                              # noqa: E402
from PIL import Image                                           # noqa: E402

import axis                                                     # noqa: E402
import bars                                                     # noqa: E402

IMG = os.path.join(HERE, "source", "pmc", "nihms-1624779-f0004.jpg")

# Panel B, located by its own long rules (see 07_qc.json fig3b.panel).
AXIS_X = 197          # the y-axis rule
BASELINE = 1223       # the x-axis rule; also the zero tick
Y_TOP = 762           # the first row below the in-panel legend band
TICK_VALUES = [90, 60, 30, 0]

SAMPLES = ["75-OCholB", "76-OCholB", "77-OCholB", "78-OCholB", "80-OCholB",
           "81-OCholB", "87-OCholB", "90-OCholB", "304-OCholB", "Lpf2k",
           "mRNA", "Control"]
SERIES = ["HeLa", "B16F10", "HEK", "NIH 3T3"]
COLOUR_CODE = ["K", "R", "B", "G"]


MIN_BAR_PX = 15      # the narrowest real bar is ~20 px at 300 dpi
MAX_WHISKER_PX = 14  # the widest black error-bar cap drawn over a short bar


def _bar_runs(s_row, x_min):
    """Bar windows from one classified row, healing whisker-split fills.

    On a bar only a few pixels tall the black error bar is drawn OVER the
    fill, so the fill reads as two short runs with a black run between them --
    the untreated control's HEK bar splits into 6 px + 5 px and the row yields
    47 runs instead of 48. Coloured runs are therefore merged across a gap
    that is entirely black and shorter than a whisker. BLACK runs are never
    merged this way: a black bar and the black outline of its neighbour are
    separated by coloured fill, and merging across that would swallow a whole
    bar.
    """
    raw = [(m.group(0)[0], m.start(), m.end())
           for m in re.finditer(r"(K+|R+|B+|G+|\.+)", s_row)
           if m.start() > x_min]
    out = []
    for c, a, b in raw:
        if c == ".":
            continue
        j = next((k for k in range(len(out) - 1, -1, -1)
                  if out[k][0] == c), None)
        # look back PAST the intervening black run: the whisker is appended to
        # `out` as its own entry, so out[-1] is that black run and not the
        # fill. Testing out[-1] alone silently failed to heal anything, and
        # the row still yielded 47 runs.
        if (j is not None and c != "K"
                and a - out[j][2] < MAX_WHISKER_PX
                and set(s_row[out[j][2]:a]) <= {"K", "."}):
            out[j] = (c, out[j][1], b)
            del out[j + 1:]
        else:
            out.append((c, a, b))
    return [r for r in out if r[2] - r[1] >= MIN_BAR_PX]



def heal(mask, xa, xb, colour, max_gap=MAX_WHISKER_PX):
    """Close gaps in ONE bar's colour mask that are bounded by its own colour.

    The black error bar is drawn OVER the fill, so inside the bar the colour
    mask has a vertical hole about 5-10 px wide. On a tall bar that costs a
    little coverage; on a short one it drops the row below bars.bar_top's 90%
    and the detector correctly refuses. Closing a gap that has the bar's OWN
    colour on both sides restores the fill without admitting anything from a
    neighbour -- the closure is computed inside [xa, xb) only, and a gap that
    runs to either edge of the window is left open.
    """
    if colour == "K":
        # NEVER heal a black bar. The error bar is black too, so there is no
        # hole to close -- and closing anyway bridges the black whisker to the
        # black top-outline of the red bar that overlaps this one's right
        # edge, filling the whole right flank five rows too high. That cost
        # three of the twelve HeLa bars a 3-5 px flank disagreement, which
        # bars.bar_top refused; it would have been a silent +0.6 to +1.1
        # percentage points had the detector returned the topmost ink.
        return mask
    sub = np.array(mask[:, xa:xb])
    w = sub.shape[1]
    for y in range(sub.shape[0]):
        row = sub[y]
        if not row.any():
            continue
        x = 0
        while x < w:
            if row[x]:
                x += 1
                continue
            j = x
            while j < w and not row[j]:
                j += 1
            if x > 0 and j < w and (j - x) <= max_gap:
                row[x:j] = True
            x = j
    out = np.array(mask)
    out[:, xa:xb] = sub
    return out


def classify(rgb):
    """Row of RGB -> one char per column. Legend-matched, not guessed."""
    r, g, b = rgb[:, 0], rgb[:, 1], rgb[:, 2]
    out = np.full(rgb.shape[0], ".", dtype="<U1")
    out[(r < 90) & (g < 90) & (b < 90)] = "K"
    out[(r > 120) & (g < 100) & (b < 100)] = "R"
    out[(b > 110) & (r < 100) & (g < 110)] = "B"
    out[(g > 90) & (r < 110) & (b < 110) & (out == ".")] = "G"
    return "".join(out)


def read_panel(y_top=None, mask_mode="colour", heal_black=False,
               tick_nudge_px=0.0):
    """The ONE implementation. main() calls it with no arguments; the seed
    harness in stage07_qc.py calls it with one thing changed, so a seeded run
    and the control differ in exactly that thing. A seed harness that
    re-implements the pipeline is error class 50 -- my first version did, and
    its 'clean' control accepted 44 of 48 bars because it had quietly dropped
    the window-refinement step."""
    y_top = Y_TOP if y_top is None else y_top
    rgb = np.asarray(Image.open(IMG).convert("RGB")).astype(int)
    gray = rgb.mean(axis=2)
    qc = {"image": os.path.relpath(IMG, HERE), "shape": list(rgb.shape)}

    # Y_TOP must sit BELOW the in-panel legend ("HeLa B16F10 HEK NIH 3T3")
    # and ABOVE every bar top. The first version of this file used 740, which
    # is inside the legend band, and bars.bar_top refused 18 of 48 bars with
    # flank deltas of 100-450 px -- the legend text filling one flank from the
    # very first scanned row. The refusals are what found it; a detector that
    # returned the topmost ink would have reported 18 bars at ~103%.
    ink_all = gray < 200
    band = [y for y in range(700, 800)
            if ink_all[y, AXIS_X + 5:1799].sum() > 200]
    if not band:
        raise SystemExit("could not find the legend band above the bars.")
    if not (max(band) < y_top <= max(band) + 4):
        qc["legend_band_warning"] = (
            f"legend band ends at row {max(band)}; y_top={y_top} is not the "
            "row just below it. main() refuses this; the seed harness sets it "
            "deliberately.")
        if y_top == Y_TOP:
            raise SystemExit(qc["legend_band_warning"])
    qc["legend_band_rows"] = [min(band), max(band)]

    # ---- y calibration, from detected ticks only ------------------------
    maj = [v + y_top for v in axis.detect_ticks_1d(
        gray[y_top:BASELINE + 20, :], "y", 178, 194, dark=140, min_len=2)]
    maj = sorted(maj)
    if tick_nudge_px:
        maj[0] += tick_nudge_px
    if len(maj) != 4:
        raise SystemExit(f"found {len(maj)} y ticks, expected 4: {maj}")
    cal = axis.calibrate(maj, TICK_VALUES, kind="linear", axis="y")
    loo = axis.loo_residuals(maj)
    qc["y_axis"] = {
        "ticks_px": [round(v, 3) for v in maj],
        "tick_values": TICK_VALUES,
        "calibration": cal.qc(),
        "leave_one_out_residuals_px": [round(float(v), 3) for v in loo],
        "max_leave_one_out_residual_px": round(float(max(loo)), 3),
        "zero_tick_vs_baseline_px": round(float(max(maj)) - BASELINE, 3),
    }

    # ---- x grid: group positions from the x-axis TICK MARKS -------------
    xt = axis.detect_ticks_1d(gray, "x", BASELINE + 4, BASELINE + 20,
                              dark=140, min_len=2)
    xt = [v for v in xt if v > AXIS_X + 40]
    # index each tick on the fitted pitch, then refit. Two of the twelve are
    # swallowed by the rotated category labels; the rest pin the ladder.
    xt = sorted(xt)
    p0 = (xt[-1] - xt[0]) / round((xt[-1] - xt[0]) / ((xt[1] - xt[0]) or 1))
    idx = [int(round((v - xt[0]) / p0)) for v in xt]
    A = np.column_stack([np.array(idx, float), np.ones(len(idx))])
    coef, *_ = np.linalg.lstsq(A, np.array(xt, float), rcond=None)
    pitch, x0 = float(coef[0]), float(coef[1])
    resid = float(np.max(np.abs(np.array(xt) - (A @ coef))))
    centres = [x0 + pitch * i for i in range(len(SAMPLES))]
    qc["x_grid"] = {
        "ticks_detected_px": [round(v, 2) for v in xt],
        "tick_index": idx, "pitch_px": round(pitch, 4),
        "first_group_px": round(x0, 3),
        "tick_fit_max_residual_px": round(resid, 3),
        "group_centres_px": [round(v, 2) for v in centres],
    }
    if resid > 2.0:
        raise SystemExit(f"x-tick ladder residual {resid:.2f} px")

    # ---- each bar's OWN x window, read at the baseline -------------------
    # A fitted grid does not survive here: the within-group bar boundaries
    # wander by up to 3 px against the group centre, and a 3 px error on a
    # 24 px bar puts a third-width flank into the inter-bar gap, where
    # bars.bar_top correctly refuses. Every bar reaches the baseline, so one
    # row above the axis rule carries all 48 fills; that row locates each bar
    # individually and needs no grid at all.
    ROW = BASELINE - 2
    s_row = classify(rgb[ROW])
    runs = _bar_runs(s_row, AXIS_X)
    if len(runs) != len(SAMPLES) * len(SERIES):
        raise SystemExit(f"row {ROW}: found {len(runs)} colour runs, expected "
                         f"{len(SAMPLES) * len(SERIES)}")
    pattern = [r[0] for r in runs]
    if pattern != COLOUR_CODE * len(SAMPLES):
        raise SystemExit(f"row {ROW}: colour pattern is not 12 x KRBG: "
                         f"{''.join(pattern)}")
    # the runs must sit on the group ladder the x TICKS define -- that is the
    # cross-check that the colour segmentation did not drift a whole bar.
    # The 48 runs must sit on the ladder the x TICKS define. Testing each
    # bar's absolute distance from its group tick is the wrong test -- the
    # four series are deliberately offset by up to ~44 px either side of it.
    # What must hold is that a given series sits at the SAME offset in every
    # group; a one-bar segmentation slip shows up there as ~24 px.
    per_series = {k: [] for k in range(len(SERIES))}
    for i, (c, a, b) in enumerate(runs):
        gi, si = i // len(SERIES), i % len(SERIES)
        per_series[si].append((a + b) / 2.0 - centres[gi])
    spread = {SERIES[k]: round(float(max(v) - min(v)), 2)
              for k, v in per_series.items()}
    qc["bar_windows"] = {
        "row_px": ROW,
        "n_runs": len(runs),
        "widths_px": [r[2] - r[1] for r in runs],
        "series_offset_from_group_tick_px":
            {SERIES[k]: round(float(np.median(v)), 2)
             for k, v in per_series.items()},
        "series_offset_spread_px": spread,
    }
    if max(spread.values()) > 6.0:
        raise SystemExit(f"a series' offset from its group tick varies by "
                         f"{max(spread.values())} px across the 12 groups; "
                         "the colour segmentation has slipped a bar.")

    # ---- read every bar -------------------------------------------------
    # THE MASK IS THE BAR'S OWN COLOUR, not ink. Adjacent bars in a group
    # OVERLAP: the black HeLa bar runs to x=393 while the red B16F10 bar that
    # is drawn over it starts at x=387, so above the red bar's top the black
    # one is still showing inside the red window. With an ink mask the red
    # bar's left flank fired ~27 px high and bars.bar_top refused -- on all
    # seven red bars whose black neighbour is taller. A per-colour mask sees
    # only this bar. The cost is the bar's own 1-2 px dark top outline, which
    # the colour mask excludes; that bias is MEASURED against the two printed
    # statistics in 07_qc.json (Lpf2k/HeLa 98%, naked mRNA/HeLa 2%) rather
    # than assumed to be zero.
    masks = {}
    r_, g_, b_ = rgb[:, :, 0], rgb[:, :, 1], rgb[:, :, 2]
    masks["K"] = (r_ < 90) & (g_ < 90) & (b_ < 90)
    masks["R"] = (r_ > 120) & (g_ < 100) & (b_ < 100)
    masks["B"] = (b_ > 110) & (r_ < 100) & (g_ < 110)
    masks["G"] = (g_ > 90) & (r_ < 110) & (b_ < 110) & ~masks["K"]
    # Refine each window to the columns that hold the bar's colour on EVERY
    # row of a six-row band just above the axis. One classified row is not
    # enough: 81-OCholB's HeLa bar reads as starting at x=891 on row 1221, but
    # x=891 is an antialiased edge column (grey 109) that is not black on any
    # other row, so the left 9 px flank topped out at 8/9 = 0.889 coverage and
    # bars.bar_top refused -- on four of the twelve black bars. Trimming the
    # flickering edge column is not pre-trimming the bar: the column is not
    # part of the bar on five rows out of six.
    BAND = 3    # rows 1219-1221; the shortest bar (untreated control, ~1.3%) is only ~6 px tall, and a six-row band reached above its top.
    refined = []
    for c, a, b in runs:
        hm = heal(masks[c], a, b, "R" if heal_black else c)
        band = hm[BASELINE - 1 - BAND:BASELINE - 1, a:b].all(axis=0)
        best, cur, start, bstart = 0, 0, 0, a
        for k, v in enumerate(band):
            if v:
                if cur == 0:
                    start = k
                cur += 1
                if cur > best:
                    best, bstart = cur, start
            else:
                cur = 0
        if best < MIN_BAR_PX:
            raise SystemExit(f"bar {c} [{a},{b}) has only {best} columns that "
                             "hold its colour on all band rows.")
        refined.append((c, a + bstart, a + bstart + best))
    qc["bar_windows"]["refined_widths_px"] = [r[2] - r[1] for r in refined]
    qc["bar_windows"]["columns_trimmed"] = [
        (o[1] - n[1], n[2] - o[2]) for o, n in zip(runs, refined)]
    runs = refined

    rows, fails = [], {}
    for i, (c, xa, xb) in enumerate(runs):
        gi, si = i // len(SERIES), i % len(SERIES)
        try:
            m = (gray < 200) if mask_mode == "ink" else heal(
                masks[c], xa, xb, "R" if heal_black else c)
            res = bars.bar_top(m, xa, xb, y_top, BASELINE - 1)
        except Exception as e:
            fails[f"{SAMPLES[gi]}/{SERIES[si]}"] = str(e).split("\n")[0]
            continue
        rows.append({
            "panel": "Fig 3B", "sample": SAMPLES[gi], "series": SERIES[si],
            "unit": "percent_GFP_positive_cells",
            "value": round(float(cal.to_value(res["row"])), 3),
            "bar_top_row_px": res["row"],
            "bar_x0_px": xa, "bar_x1_px": xb,
            "coverage": res["coverage"],
            "flank_delta_px": res["flank_delta_px"],
        })

    clipped = [r for r in rows if r["bar_top_row_px"] <= y_top + 2]
    if clipped:
        raise SystemExit(f"{len(clipped)} bar(s) detected at the top of the "
                         f"scan window; the window clips them: "
                         f"{[(r['sample'], r['series']) for r in clipped]}")
    qc["n_values"] = len(rows)
    qc["refusals"] = fails
    return rows, qc, cal


def main():
    rows, qc, cal = read_panel()
    fails = qc["refusals"]
    resid = qc["x_grid"]["tick_fit_max_residual_px"]
    pitch = qc["x_grid"]["pitch_px"]
    out = os.path.join(HERE, "05b_fig3b_values.csv")
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    with open(os.path.join(HERE, "05b_fig3b_qc.json"), "w") as fh:
        json.dump(qc, fh, indent=1)
    print(f"wrote {out}: {len(rows)}/48 bars")
    print(f"  y tick-fit residual {cal.qc()['tick_fit_residual_px']} px; "
          f"max LOO {qc['y_axis']['max_leave_one_out_residual_px']} px")
    print(f"  x ladder residual {resid:.3f} px; pitch {pitch:.3f} px")
    for k, v in fails.items():
        print("  REFUSED", k, v[:110])


if __name__ == "__main__":
    main()
