#!/usr/bin/env python3
"""Fig 4B -- liver flux per INDIVIDUAL formulation, intraduodenal, 100 ug pDNA.

Tier 5 (pixels). Why tiers 1-4 are unavailable is in 00_value_sources.json.
Geometry comes from curation/lib; this file holds only this panel's own
constants, next to the figure, where a reviewer reads them.

THE SCALE IS FITTED ON THE THREE MAJOR TICKS ONLY. The ten detected minor
ticks are held out and used as the test. They are NOT at log10(2..9): this
panel's minor ticks are EQUALLY SPACED IN LINEAR VALUE inside each decade
(2.8, 4.6, 6.4, 8.2 -- an Origin "4 minor ticks" log axis), so they are
UNEQUALLY spaced in the fitted coordinate and a wrong scale cannot reproduce
them by accident.

THREE LEVELS, ALL READ FROM THE PANEL, NONE ASSUMED: the white page, the bar's
dotted fill, and the dark border stroke. The ink mask is the midpoint of the
first two (so the detector finds where the page ends); the reported value is
the stroke CENTRE, because the calibration was fitted on tick centroids and
reporting a stroke's outer edge against a centroid scale is F-019.
"""
import json, os, sys
import numpy as np
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "lib"))
import axis, bars                                        # noqa: E402

FIG = os.environ.get("STAGE01_FIG", os.path.join(
    HERE, "..", "source", "figs_elsevier", "gr4_lrg.jpg"))
OUT = os.environ.get("STAGE01_OUT", os.path.join(HERE, "..", "01_fig4b.json"))

AXIS_X        = 214          # y-axis rule occupies columns 212-216
TICK_STRIP    = (205, 212)   # left of the rule, right of the labels
TICK_Y        = (1600, 2340)
MAJOR_MIN_LEN = 9            # measured tick lengths: 10 (major) vs 7-8 (minor).
                             # Selecting majors by SPACING fails -- the minor at
                             # 1741.70 is 100.7 px from the 1e6 tick.
BASELINE_SEARCH = (2300, 2360)   # the x-axis rule is DETECTED in here, not
                                 # hard-coded: a constant cannot be wrong in a
                                 # way the round-5,000 check can see, and a
                                 # check that cannot see its own input is the
                                 # F-038 shape.
PLOT_X        = (210, 2070)
PAGE_ROWS     = (1645, 1680) # inside the plot box, above every bar
MAJOR_VALUES  = [1e6, 1e5, 1e4]
MINOR_MULTIPLIERS = (8.2, 6.4, 4.6, 2.8)   # (10-1)/5 = 1.8 linear steps
CHROMA_MIN    = 25
MIN_LEVEL_SEPARATION = 30    # min(RGB) units, asserted per bar and reported

# --- GATES. Each is made to fire on purpose by scripts/seed_checks.py; a
# check that cannot fail is not a check (CLAUDE.md Sec 6).
GATE_DECADE_GAP_PX   = 1.0   # the two decade gaps must agree
GATE_MINOR_RESID_PX  = 1.5   # held-out minor ticks, worst case
GATE_BASELINE_PX     = 1.0   # the drawn baseline vs a round 5,000
GATE_N_BARS          = 20

LABELS = (["AI-a%d" % i for i in (1, 2, 3, 4)]
          + ["AI-b%d" % i for i in (1, 2, 3, 4)]
          + ["CII-b%d" % i for i in (1, 2, 3, 4)]
          + ["DI-a%d" % i for i in (1, 2, 3, 4)]
          + ["FIII-b%d" % i for i in (1, 2, 3, 4)])
GROUPS = ["DOTAP (AI)"] * 8 + ["DOPE (CII)"] * 4 + ["DSPC (DI)"] * 4 \
         + ["18PG (FIII)"] * 4


def detect_bar_windows(rgb, y):
    row = rgb[y, PLOT_X[0]:PLOT_X[1]]
    on = (row.max(axis=1) - row.min(axis=1)) > CHROMA_MIN
    runs, s = [], None
    for i, v in enumerate(on):
        if v and s is None:
            s = i
        elif not v and s is not None:
            runs.append((PLOT_X[0] + s, PLOT_X[0] + i)); s = None
    if s is not None:
        runs.append((PLOT_X[0] + s, PLOT_X[1]))
    return [r for r in runs if r[1] - r[0] >= 10]


def stroke_centre(prof, row_hint, bg_above, bg_below):
    """Sub-pixel row of a bar's top stroke, by two-sided half maximum.

    `prof` is darkness (255 - min(RGB)) averaged over the full bar width. The
    background differs on the two sides -- the white page above, the bar's own
    fill below -- so a single centroid is biased toward the fill. Each flank's
    half maximum is taken against ITS OWN background and the two averaged.
    """
    peak = float(max(prof[row_hint:row_hint + 6]))
    hm_up = (peak + bg_above) / 2.0
    i = int(row_hint)
    while i > 1 and prof[i] > hm_up:
        i -= 1
    if prof[i + 1] == prof[i]:
        raise SystemExit("stage01: flat upper flank at row %d" % i)
    up = i + (hm_up - prof[i]) / (prof[i + 1] - prof[i])
    hm_dn = (peak + bg_below) / 2.0
    j = int(row_hint)
    while j < len(prof) - 1 and prof[j] > hm_dn:
        j += 1
    if prof[j - 1] == prof[j]:
        raise SystemExit("stage01: flat lower flank at row %d" % j)
    dn = (j - 1) + (prof[j - 1] - hm_dn) / (prof[j - 1] - prof[j])
    return (up + dn) / 2.0, up, dn, peak


def detect_baseline(mn, lo, hi, x0, x1):
    """Darkness-weighted centroid row of the x-axis rule, in the same
    coordinate the tick centroids are measured in."""
    prof = (255.0 - mn[lo:hi, x0:x1]).mean(axis=1)
    peak = prof.max()
    on = prof > peak / 2.0
    idx = np.where(on)[0]
    if idx.size == 0:
        raise SystemExit("stage01: no x-axis rule in rows %d-%d" % (lo, hi))
    a, b = idx.min(), idx.max() + 1
    if b - a > 12:
        raise SystemExit("stage01: x-axis rule candidate is %d px thick" % (b - a))
    w = prof[a:b]
    return float(np.dot(np.arange(lo + a, lo + b), w) / w.sum())


def main():
    rgb = np.asarray(Image.open(FIG).convert("RGB")).astype(int)
    gray = rgb.mean(axis=2)
    mn = rgb.min(axis=2)
    BASELINE_Y = detect_baseline(mn, BASELINE_SEARCH[0], BASELINE_SEARCH[1],
                                 PLOT_X[0] + 400, PLOT_X[1] - 400)

    ticks = sorted(t for t in axis.detect_ticks_1d(
        gray, "y", TICK_STRIP[0], TICK_STRIP[1], dark=140, min_len=3,
        max_len=40) if TICK_Y[0] < t < TICK_Y[1])

    def tick_len(t):
        band = gray[int(round(t)) - 3:int(round(t)) + 4, :AXIS_X].min(axis=0)
        x = AXIS_X - 1
        while x > 100 and band[x] < 160:
            x -= 1
        return AXIS_X - 1 - x

    lens = {t: tick_len(t) for t in ticks}
    majors = sorted(t for t in ticks if lens[t] >= MAJOR_MIN_LEN)
    minors = sorted(t for t in ticks if lens[t] < MAJOR_MIN_LEN)
    if len(majors) != 3:
        raise SystemExit("stage01: expected 3 majors, got %r" % (majors,))
    if min(lens[t] for t in majors) - max(lens[t] for t in minors) < 1:
        raise SystemExit("stage01: tick lengths not separated: %r" % (lens,))

    cal = axis.calibrate(majors, MAJOR_VALUES, kind="log10", axis="y")

    predicted = {k * 10 ** (d - 1): cal.to_pixel(k * 10 ** (d - 1))
                 for d in (6, 5, 4) for k in MINOR_MULTIPLIERS}
    minor_res = []
    for px in minors:
        v, p = min(predicted.items(), key=lambda kv: abs(kv[1] - px))
        minor_res.append({"measured_px": round(px, 3), "value": v,
                          "predicted_px": round(p, 3),
                          "residual_px": round(px - p, 3)})
    mr = np.array([abs(m["residual_px"]) for m in minor_res])

    page = float(np.median(mn[PAGE_ROWS[0]:PAGE_ROWS[1],
                              PLOT_X[0] + 20:PLOT_X[1] - 20]))
    wins = detect_bar_windows(rgb, 2300)
    if len(wins) != 20:
        raise SystemExit("stage01: expected 20 bars, found %d" % len(wins))

    out = []
    for i, (x0, x1) in enumerate(wins):
        interior = mn[int(BASELINE_Y) - 60:int(BASELINE_Y) - 10, x0 + 8:x1 - 8]
        fill = float(np.median(interior))
        strokev = float(np.percentile(mn[1630:int(BASELINE_Y) - 2, x0:x1], 1))
        if page - fill < MIN_LEVEL_SEPARATION or \
           fill - strokev < MIN_LEVEL_SEPARATION:
            raise SystemExit(
                "stage01: bar %d (%s): levels page=%.1f fill=%.1f stroke=%.1f "
                "are not separated by %d min(RGB) units, so any threshold "
                "between them is a guess." % (i, LABELS[i], page, fill,
                                              strokev, MIN_LEVEL_SEPARATION))
        thr = (page + fill) / 2.0
        r = bars.bar_top(mn < thr, x0, x1, 1630, int(BASELINE_Y) - 2,
                         baseline_row=int(BASELINE_Y))
        prof = (255.0 - mn[:, x0:x1]).mean(axis=1)
        centre, up, dn, peak = stroke_centre(prof, r["row"], 255.0 - page,
                                             255.0 - fill)
        out.append({"label": LABELS[i], "group": GROUPS[i],
                    "x0": x0, "x1": x1,
                    "level_page": page, "level_fill": round(fill, 1),
                    "level_stroke": round(strokev, 1),
                    "threshold_minRGB": round(thr, 1),
                    "bar_top_row": r["row"],
                    "flank_delta_px": r["flank_delta_px"],
                    "coverage": r["coverage"],
                    "stroke_peak_darkness": round(peak, 1),
                    "halfmax_upper_row": round(up, 3),
                    "halfmax_lower_row": round(dn, 3),
                    "stroke_centre_row": round(centre, 3),
                    "centre_minus_bar_top_px": round(centre - r["row"], 3),
                    "value_at_bar_top_row": float(cal.to_value(r["row"])),
                    "value_photons_per_s": float(cal.to_value(centre))})

    # ------------------------- gates ----------------------------------
    gaps_d = abs((majors[1] - majors[0]) - (majors[2] - majors[1]))
    base_off = abs(BASELINE_Y - cal.to_pixel(5000.0))
    failures = []
    if gaps_d > GATE_DECADE_GAP_PX:
        failures.append("decade gaps differ by %.3f px > %.1f"
                        % (gaps_d, GATE_DECADE_GAP_PX))
    if mr.max() > GATE_MINOR_RESID_PX:
        failures.append("held-out minor tick residual %.3f px > %.1f"
                        % (mr.max(), GATE_MINOR_RESID_PX))
    if base_off > GATE_BASELINE_PX:
        failures.append("baseline is %.3f px from a round 5,000 > %.1f"
                        % (base_off, GATE_BASELINE_PX))
    print("GATE decade_gap_difference_px %.4f" % gaps_d)
    print("GATE minor_tick_residual_px_max %.4f" % mr.max())
    print("GATE baseline_offset_px %.4f" % base_off)
    print("GATE n_bars %d" % len(out))

    ppd = cal.px_per_unit
    corr = [b["centre_minus_bar_top_px"] for b in out]
    qc = {"panel": "Fig 4B",
          "calibration": cal.qc(),
          "tick_lengths_px": {str(round(t, 2)): lens[t] for t in ticks},
          "major_ticks_px": [round(t, 3) for t in majors],
          "decade_gaps_px": [round(majors[1] - majors[0], 3),
                             round(majors[2] - majors[1], 3)],
          "decade_gap_difference_px": round(
              abs((majors[1] - majors[0]) - (majors[2] - majors[1])), 4),
          "minor_ticks_held_out": minor_res,
          "minor_tick_n": len(minor_res),
          "minor_tick_residual_px_max": round(float(mr.max()), 3),
          "minor_tick_residual_px_rms": round(float(np.sqrt((mr ** 2).mean())), 3),
          "baseline_row_detected": round(BASELINE_Y, 3),
          "baseline_value_read": round(float(cal.to_value(BASELINE_Y)), 3),
          "baseline_round_number": 5000.0,
          "baseline_offset_px_from_5000": round(
              float(BASELINE_Y - cal.to_pixel(5000.0)), 4),
          "levels": {"page": page},
          "precision": {
              "px_per_decade": round(ppd, 4),
              "factor_per_px": round(10 ** (1.0 / ppd), 6),
              "worst_minor_residual_as_factor": round(
                  10 ** (mr.max() / ppd), 6)},
          "stroke_centre_correction": {
              "px_min": round(min(corr), 3), "px_max": round(max(corr), 3),
              "px_mean": round(float(np.mean(corr)), 3),
              "as_factor_at_mean": round(10 ** (float(np.mean(corr)) / ppd), 5),
              "_why": "bar_top returns the first row where the page ends, i.e. "
                      "the stroke's upper half. The calibration is fitted on "
                      "tick CENTROIDS. Reporting an edge against a centroid "
                      "scale is F-019; every value below is at the stroke "
                      "centre and this is the bias that removes."},
          "gate_failures": failures,
          "bars": out}
    with open(OUT, "w") as fh:
        json.dump(qc, fh, indent=1)
    print(json.dumps({k: v for k, v in qc.items()
                      if k not in ("bars", "minor_ticks_held_out",
                                   "tick_lengths_px")}, indent=1))
    for b in out:
        print("  %-8s %-12s top=%4d centre=%9.3f flank=%.1f  %.4e"
              % (b["label"], b["group"], b["bar_top_row"],
                 b["stroke_centre_row"], b["flank_delta_px"],
                 b["value_photons_per_s"]))
    if failures:
        raise SystemExit("stage01 GATE FAILED: " + "; ".join(failures))


if __name__ == "__main__":
    main()
