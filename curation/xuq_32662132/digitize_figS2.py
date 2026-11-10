"""Figure S2 a/c/d - bar tops for 93-O17S, through curation/lib only.

No geometry is written here. `axis.calibrate` fits scale AND origin from
detected tick positions; `bars.bar_top` requires full-width coverage plus
left/right flank agreement so a centred s.d. cap cannot be read as a bar top
(F-018, error class 32). Both refuse rather than return a plausible number.

Panel b (CCK-8 cell viability) is measured too but is NOT a delivery endpoint
and is not emitted; it is here as an independent cross-check that the same
detector reproduces a panel whose values the SI prose constrains.
"""
import json, os, sys
import numpy as np
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "lib")))
import axis, bars                                        # noqa: E402
from errors import CalibrationError, DetectionError      # noqa: E402

RASTER = os.path.join(HERE, "source", "figS2_page5_xref101.png")

# Bar fill is a distinct mid-grey; the axis rule, the outline and the s.d.
# whiskers are near-black. Measured on the raster: fill (59,59,58), ink (0,5,3).
FILL = np.array([59, 59, 58])
FILL_TOL = 26          # L-inf; keeps (59,59,58) and rejects (0,5,3) and white

PANELS = {
    "a": dict(spine_x=(132, 137), tick_scan=(112, 131), plot_y=(78, 610),
              baseline=610, bars_x=(140, 936),
              tick_values=[8000, 6000, 4000, 2000, 0], unit="RLU",
              n_expected=14),
    "b": dict(spine_x=(1045, 1048), tick_scan=(1025, 1044), plot_y=(129, 606),
              baseline=606, bars_x=(1052, 1525),
              tick_values=[100, 50, 0], unit="percent_viability",
              n_expected=6),
    "c": dict(spine_x=(213, 216), tick_scan=(193, 212), plot_y=(954, 1432),
              baseline=1432, bars_x=(220, 693),
              tick_values=[10000, 5000, 0], unit="RLU", n_expected=8),
    "d": dict(spine_x=(976, 979), tick_scan=(956, 975), plot_y=(949, 1426),
              baseline=1426, bars_x=(983, 1456),
              tick_values=[40000, 20000, 0], unit="RLU", n_expected=6),
}


def fill_mask(rgb):
    return (np.abs(rgb.astype(int) - FILL).max(axis=2) <= FILL_TOL)


def bar_windows(mask, x0, x1, baseline, y_top=0, min_px=3, min_w=8):
    """Column windows carrying the bar fill anywhere in the plot area.

    Counting over the WHOLE column rather than a band just above the baseline:
    a band misses a bar shorter than the band is tall, and in panel a the two
    control bars are ~100 RLU, i.e. ~7 px. A band-based detector silently
    dropped them, and `n_expected` is what made that visible rather than
    producing a 12-bar panel that looked complete.
    """
    band = mask[y_top:baseline, x0:x1]
    present = band.sum(axis=0) >= min_px
    out, start = [], None
    for j, v in enumerate(present):
        if v and start is None:
            start = j
        elif not v and start is not None:
            if j - start >= min_w:
                out.append((x0 + start, x0 + j))
            start = None
    if start is not None and len(present) - start >= min_w:
        out.append((x0 + start, x0 + len(present)))
    return out


def main():
    rgb = np.array(Image.open(RASTER).convert("RGB"))
    ink_dark = (rgb.astype(int).mean(axis=2) < 150)   # fill + outline + whisker
    fill = fill_mask(rgb)
    gray = rgb.astype(float).mean(axis=2)
    report = {}

    for pid, P in PANELS.items():
        rec = {"panel": pid, "unit": P["unit"]}
        # restrict to THIS panel's rows before detecting: the tick strip column
        # range also contains the rotated x-axis labels of the panel below.
        y0, y1 = P["plot_y"][0] - 8, P["baseline"] + 10
        ticks = [t + y0 for t in axis.detect_ticks_1d(
            gray[y0:y1, :], along="y",
            scan_lo=P["tick_scan"][0], scan_hi=P["tick_scan"][1])]
        rec["ticks_px_detected"] = [round(float(t), 3) for t in ticks]
        if len(ticks) != len(P["tick_values"]):
            rec["REFUSED"] = ("tick detector found %d ticks, the axis is labelled with %d"
                              % (len(ticks), len(P["tick_values"])))
            report[pid] = rec
            continue
        cal = axis.calibrate(ticks, P["tick_values"], kind="linear", axis="y")
        rec["calibration"] = cal.qc()
        loo = axis.loo_residuals(ticks)
        rec["loo_residuals_px"] = [round(v, 4) for v in loo]
        rec["loo_worst_px"] = round(max(loo), 4)
        rec["loo_tol_px"] = axis.LOO_TOL_PX
        rec["loo_passes"] = bool(max(loo) <= axis.LOO_TOL_PX)

        zero_px = cal.to_pixel(0.0)
        rec["zero_tick_px"] = round(float(zero_px), 3)
        base = int(round(zero_px)) + 5          # bar feet sit a few px below the 0 tick
        rec["baseline_row_used"] = base
        wins = bar_windows(fill, P["bars_x"][0], P["bars_x"][1], base,
                           y_top=P["plot_y"][0])
        rec["n_bar_windows_fill"] = len(wins)
        if len(wins) != P["n_expected"]:
            # SECOND PASS, named and reported, never a silent relaxation. A bar
            # a few px tall is all outline and carries no interior fill, so the
            # fill detector cannot see it. Detect those on the outline mask,
            # above the baseline rule, and merge. The count must then reach
            # n_expected or the panel is refused.
            outline = ink_dark.copy()
            outline[int(round(zero_px)) - 2:, :] = False     # drop the axis rule
            extra = bar_windows(outline, P["bars_x"][0], P["bars_x"][1],
                                int(round(zero_px)) - 2, y_top=P["plot_y"][0],
                                min_px=2)
            merged = list(wins)
            for (a0, a1) in extra:
                if not any(not (a1 <= b0 or a0 >= b1) for (b0, b1) in wins):
                    merged.append((a0, a1))
            wins = sorted(merged)
            rec["n_bar_windows_outline_pass"] = len(extra)
            rec["windows_added_by_outline_pass"] = [w for w in wins if w not in
                                                    rec.get("_f", [])][:0] or \
                [w for w in wins if all(w != f for f in rec.get("windows_fill", []))]
        rec["n_bar_windows"] = len(wins)
        rec["n_expected"] = P["n_expected"]
        rec["windows_px"] = wins
        if len(wins) != P["n_expected"]:
            # A category can legitimately be drawn as NO bar. That is only
            # acceptable if the slot where the bar would be is measurably
            # empty, so it is measured here rather than assumed. The slot is
            # placed from the pitch of the bars that WERE found.
            miss = P["n_expected"] - len(wins)
            centres = [(a + b) / 2.0 for a, b in wins]
            ok = False
            rec["undrawn_category_test"] = {}
            if miss > 0 and len(centres) >= 2:
                pitch = (centres[-1] - centres[0]) / (len(centres) - 1)
                w = int(round(np.mean([b - a for a, b in wins])))
                slots = [centres[-1] + pitch * (k + 1) for k in range(miss)]
                counts = []
                for c in slots:
                    s0, s1 = int(round(c - w / 2)), int(round(c + w / 2))
                    s0, s1 = max(0, s0), min(ink_dark.shape[1], s1)
                    counts.append(int(ink_dark[P["plot_y"][0]:int(round(zero_px)) - 2,
                                               s0:s1].sum()))
                rec["undrawn_category_test"] = {
                    "pitch_px": round(pitch, 2), "bar_width_px": w,
                    "slot_centres_px": [round(c, 1) for c in slots],
                    "ink_pixels_in_slot": counts,
                    "verdict": ("every missing slot is empty, so the category is "
                                "drawn as no bar at all" if all(c == 0 for c in counts)
                                else "a missing slot CONTAINS INK - the detector missed "
                                     "a real bar")}
                ok = all(c == 0 for c in counts)
            if not ok:
                rec["REFUSED"] = ("found %d bar windows, expected %d" %
                                  (len(wins), P["n_expected"]))
                report[pid] = rec
                continue

        vals, fails = [], {}
        for i, (x0, x1) in enumerate(wins):
            try:
                t = bars.bar_top(ink_dark, x0, x1, P["plot_y"][0], base,
                                 baseline_row=base)
            except DetectionError as e:
                fails[i] = str(e)
                vals.append(None)
                continue
            top = float(t["row"])
            vals.append(dict(index=i, x0=int(x0), x1=int(x1), top_row=top,
                             value=round(float(cal.to_value(top)), 2), detail=t))
        rec["bars"] = vals
        rec["failures"] = fails
        report[pid] = rec

    with open(os.path.join(HERE, "05_figS2_geometry.json"), "w") as f:
        json.dump(report, f, indent=1, default=str)
    for pid, rec in report.items():
        if "REFUSED" in rec:
            print(pid, "REFUSED:", rec["REFUSED"])
        else:
            print(pid, "ticks", rec["calibration"]["n_ticks"],
                  "resid", rec["calibration"]["tick_fit_residual_px"],
                  "loo_worst", rec["loo_worst_px"],
                  "bars", rec["n_bar_windows"],
                  [b["value"] if b else None for b in rec["bars"]])

if __name__ == "__main__":
    main()
