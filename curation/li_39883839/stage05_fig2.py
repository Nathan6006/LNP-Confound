#!/usr/bin/env python3
"""Digitise Figure 2H (batch 3, individual lipids) and Figure 2K -> 05_fig2.json.

Tier 5, for the reasons in 00_value_sources.json.

THE AXIS TRAP IN PANEL K, and how it is resolved WITHOUT any other panel.
Panel K's y axis is drawn with LOG-spaced minor ticks under LINEAR major
labels. Between 0 and 1e7 the minor ticks sit at pixel offsets 45, 72, 90,
105, 117, 127, 136, 143 from the 0 tick -- which is 151 px per decade, i.e.
2,3,4,...,9 x 10^6 on a log scale -- while the five LABELLED ticks 0, 1e7,
2e7, 3e7, 4e7 sit 151.0 / 150.5 / 150.5 / 150.7 px apart, which is linear.
The two readings differ by a factor of ~22 on a mid-panel bar.

It is decided on the labelled ticks alone: four equal value steps at four
equal pixel steps is a linear axis, and no log axis can carry a '0' tick at
all. The minor ticks are a drawing artifact and are recorded as contradiction
C-06. Only the five MAJOR ticks are passed to axis.calibrate, selected by
TICK LENGTH (major 23-25 px, minor 11 px) -- measured, not assumed.

WHAT IS REFUSED, and why no number is invented for it:
  Fig 2K SM-102   its bar top is in the UPPER axis segment, which carries one
                  labelled tick (1e8). axis.calibrate needs >= 3. One tick
                  cannot give a scale, so the panel does not determine this
                  value.
  Fig 2H B1, B2   drawn with no resolvable height. Reported as an UPPER BOUND
                  at the detection floor, not as a measurement of zero.
"""
import json, os, sys
import numpy as np
from PIL import Image

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "lib"))
import axis, bars                                            # noqa: E402
from errors import CurationLibError                          # noqa: E402

Image.MAX_IMAGE_PIXELS = None
HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "source", "fig", "pnas.2409572122fig02.jpg")

MAJOR_MIN_LEN = 18     # px; measured majors are 23-25, minors 11
PANELS = {
 "Fig2H": dict(spine=3470, y=(2350, 3080), baseline=3051, bar_w=101,
               tick_values=[1.5e8, 1.0e8, 5.0e7, 0.0],
               bars=[("H2A4B1", 3570), ("H2A4B2", 3721),
                     ("H2A4B3", 3872), ("H2A4B4", 4023)]),
 "Fig2K": dict(spine=3459, y=(3380, 4140), baseline=4115, bar_w=74,
               # the 4e7 tick is drawn short because it abuts the axis break,
               # so it is named here and verified against the fit, not fed in.
               tick_values=[3.0e7, 2.0e7, 1.0e7, 0.0],
               check_tick=(3512.0, 4.0e7),
               bars=[("H2A4B3", 3537), ("H2A4B4", 3653),
                     ("H2A4-Ole", 3769), ("H2A4-Butyl", 3885),
                     ("DLin-MC3-DMA", 4000.5), ("SM-102", 4116)]),
}


def tick_rows(gray, spine, ylo, yhi, maxlen=60):
    """(row, length) for every tick left of the spine. Length separates the
    labelled majors from the decorative minors; nothing else does here."""
    seen = []
    for y in range(ylo, yhi):
        x, n = spine - 4, 0
        while x > spine - maxlen and gray[y, x] < 140:
            n += 1; x -= 1
        if n:
            seen.append((y, n))
    out, cur = [], [seen[0]]
    for t in seen[1:]:
        if t[0] - cur[-1][0] <= 1:
            cur.append(t)
        else:
            out.append(cur); cur = [t]
    out.append(cur)
    return [(sum(y * n for y, n in c) / sum(n for y, n in c),
             max(n for y, n in c)) for c in out]



def longest_run(row):
    """Longest contiguous True run. A significance bracket is ONE horizontal
    rule wider than a bar; a row that merely TOTALS more ink than a bar can be
    several bar edges plus the baseline, which is what a sum-based test picked
    up two pixels above the axis."""
    best = cur = 0
    for v in row:
        cur = cur + 1 if v else 0
        if cur > best:
            best = cur
    return best


def main():
    im = Image.open(FIG)
    gray = np.asarray(im.convert("L")).astype(float)
    ink = gray < 140
    res = {}

    for pname, P in PANELS.items():
        allt = tick_rows(gray, P["spine"], *P["y"])
        major = [t for t, L in allt if L >= MAJOR_MIN_LEN]
        minor = [t for t, L in allt if L < MAJOR_MIN_LEN]
        want = len(P["tick_values"])
        major_in_segment = sorted(major)[-want:] if len(major) > want else major
        if len(major_in_segment) != want:
            raise SystemExit("%s: %d major ticks, expected %d: %s"
                             % (pname, len(major_in_segment), want,
                                [round(t, 1) for t in major]))
        cal = axis.calibrate(major_in_segment, sorted(P["tick_values"],
                                                     reverse=True),
                             kind="linear", axis="y")
        entry = {"calibration": cal.qc(),
                 "major_tick_rows_px": [round(t, 2) for t in major_in_segment],
                 "n_minor_ticks_ignored": len(minor),
                 "bars": []}

        # a named short tick is VERIFIED against the fit, never fed into it
        if "check_tick" in P:
            row, val = P["check_tick"]
            entry["short_tick_check"] = {
                "value": val, "drawn_row_px": row,
                "predicted_row_px": round(cal.to_pixel(val), 2),
                "delta_px": round(row - cal.to_pixel(val), 2)}

        # significance brackets share the bars' colour; detect, do not type in
        bw = P["bar_w"]
        x_lo = int(min(c for _, c in P["bars"]) - bw)
        x_hi = int(max(c for _, c in P["bars"]) + bw)
        br = [y for y in range(P["y"][0], P["baseline"] - 10)
              if longest_run(ink[y, x_lo:x_hi]) > bw + 30]
        y_top = (max(br) + 8) if br else P["y"][0]
        entry["significance_bracket_rows_px"] = [min(br), max(br)] if br else None
        entry["scan_top_row_px"] = y_top

        for label, centre in P["bars"]:
            x0 = int(round(centre - bw / 2))
            x1 = int(round(centre + bw / 2)) + 1
            rec = {"lipid": label, "bar_window_px": [x0, x1]}
            try:
                top = bars.bar_top(ink, x0, x1, y_top=y_top,
                                   y_bot=P["baseline"] - 2,
                                   baseline_row=P["baseline"])
                rec.update(top_row_px=top["row"],
                           flank_delta_px=top["flank_delta_px"],
                           coverage=top["coverage"],
                           value_p_per_s=float("%.4g" % cal.to_value(top["row"])))
                if top["row"] < min(major_in_segment):
                    rec["value_p_per_s"] = None
                    rec["REFUSED"] = (
                        "bar top row %d is ABOVE the highest calibrated tick "
                        "(%d). It lies in the upper segment of the broken "
                        "axis, which carries one labelled tick; a scale cannot "
                        "be fitted to one tick, so this panel does not "
                        "determine this value."
                        % (top["row"], int(min(major_in_segment))))
            except CurationLibError as e:
                # A refusal has two very different causes here. Distinguish
                # them, because one is "no height" and the other is "the panel
                # cannot express this value at all".
                top_tick = int(min(major_in_segment))
                edges = ink[P["y"][0]:top_tick - 6, [x0 + 1, x1 - 2]].sum()
                rec["REFUSED"] = str(e).split("\n")[0]
                if edges > 100:
                    rec["value_p_per_s"] = None
                    # The panel does not give a POINT value, but it does give
                    # a bound: the bar passes the axis break, and the 4e7 tick
                    # above the fit is verified by short_tick_check to 0.25 px.
                    if "check_tick" in P:
                        rec["lower_bound_p_per_s"] = P["check_tick"][1]
                        rec["bound_note"] = (
                            "The bar top is above the axis break. The top of "
                            "the calibrated lower segment is the 4e7 tick, "
                            "which this fit predicts to within 0.25 px, so "
                            "the value is a LOWER BOUND of 4e7. The upper "
                            "segment cannot be calibrated, so no point value "
                            "is available.")
                    rec["REFUSED"] = (
                        "this bar's own side edges carry %d ink px ABOVE the "
                        "highest calibrated tick (%g at row %d), so its top "
                        "lies in the UPPER segment of the broken axis. That "
                        "segment carries ONE labelled tick (1e8); "
                        "axis.calibrate needs >= 3, and two points define a "
                        "line through anything. The panel does not determine "
                        "this value and none is invented."
                        % (int(edges), max(P["tick_values"]), top_tick))
                    entry["bars"].append(rec)
                    continue
                floor = cal.to_value(P["baseline"] - bars.NROW)
                rec["upper_bound_p_per_s"] = float("%.3g" % floor)
                rec["bound_note"] = (
                    "No resolvable bar height: the top rule is within %d px of "
                    "the baseline. Reported as an UPPER BOUND at the detection "
                    "floor, which is what the drawing supports; it is not a "
                    "measurement of zero." % bars.NROW)
            entry["bars"].append(rec)
        res[pname] = entry

    json.dump(res, open(os.path.join(HERE, "05_fig2.json"), "w"), indent=1)
    for p, e in res.items():
        print("==", p, "resid %.3f px, %d majors, %d minors ignored"
              % (e["calibration"]["tick_fit_residual_px"],
                 e["calibration"]["n_ticks"], e["n_minor_ticks_ignored"]))
        if "short_tick_check" in e:
            print("   short-tick check:", e["short_tick_check"])
        for b in e["bars"]:
            print("   %-14s %-12s %s" % (b["lipid"],
                  b.get("value_p_per_s"), b.get("REFUSED", "")[:70]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
