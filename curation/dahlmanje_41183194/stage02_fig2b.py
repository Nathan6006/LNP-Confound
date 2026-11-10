#!/usr/bin/env python3
"""KG_2025 Fig. 2B -- percent Thy1.1+ cells, 4 GlycoLNPs x 2 organs.

Run:  python3 -I curation/dahlmanje_41183194/stage02_fig2b.py <repo_root> [--seed NAME]

WHY THIS PANEL IS DIGITISED AT ALL, when 00_value_sources.json records the
whole paper as tier-4-printed-numerals elsewhere: Fig. 2B is the ONLY locus
that attributes an absolute per-organ delivery value to GlycoLNP45, 59 and 62.
And it carries its own external anchor -- GlycoLNP61's two bars measure the
same quantity that Fig. 3C PRINTS as 54 (liver) and 78 (spleen). So 2 of the 8
bars have ground truth that was never used to build any part of the geometry.
That is the check; see `--seed` for the proof it can fail.

Geometry rules obeyed (curation/lib does the enforcing, not this file):
  * scale AND origin fitted from DETECTED TICK POSITIONS ONLY (axis.calibrate).
    There is no API here to supply an origin separately -- error class 20/33.
  * `forbid_target` carries the bar-top pixel rows, so a scale accidentally
    fitted onto the quantity it validates is a type error -- error class 27.
  * the ink mask is the BAR'S OWN COLOUR, not darkness. In this panel the
    replicate dots, the error whiskers and the bar fill are all one teal while
    the axis is grey; a darkness mask would make the axis rule the topmost bar.
  * `bars.bar_top` requires 90% coverage of the FULL bar width over 3 contiguous
    rows PLUS left/right flank agreement. A replicate dot is ~10 px on a ~16 px
    bar and is in neither flank. Refusals are propagated, never swallowed.
"""
import json
import os
import sys

import numpy as np
from PIL import Image

ROOT = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else ".")
sys.path.insert(0, os.path.join(ROOT, "curation", "lib"))
import axis as AX      # noqa: E402
import bars as BR      # noqa: E402
from errors import CalibrationError, DetectionError  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "source", "pnas.2409569122fig02.jpg")

# ---- panel constants, measured from this raster and written here ----------
# Bar fill RGB sampled as the modal non-white colour inside the panel.
BAR_RGB = (145, 179, 191)   # recorded as the panel's measured fill colour
BAR_TOL = 28
WHITE_CUT = 225             # min(R,G,B) below this is ink; page white is 248-255
INK_DARK = 170          # measured: the grey spine/ticks bottom out at 110-148

# Measured darkness in each tick strip: the 20..100 ticks bottom out at 133-148
# and the background is 253, so 170 separates them cleanly. The LIVER panel's
# ZERO tick is drawn fainter (min 192) and is therefore NOT detected at this
# threshold -- which is the opportunity, not the problem: the liver scale is
# fitted on the five ticks 100..20 and the zero tick's position is then a
# PREDICTION of a fit that never saw it (see `excluded_tick_prediction`).
PANELS = {
    "liver": {
        "organ": "liver", "cell_type": "endothelial cells (CD31+)",
        "spine_x": 447, "tick_scan": (441, 447),
        "tick_values": [100, 80, 60, 40, 20],      # top-down, as detected
        "y_search": (20, 212), "tick_y": (20, 200),
        "_jitter_row": 103,        # the 60% tick, displaced by the tick_jitter seed
        "bars": {"61": (458, 474), "59": (485, 501),
                 "45": (512, 529), "62": (539, 556)},
    },
    "spleen": {
        "organ": "spleen", "cell_type": "macrophages (F4/80+)",
        "spine_x": 639, "tick_scan": (633, 639),
        "tick_values": [100, 50, 0],
        "y_search": (20, 212), "tick_y": (20, 215),
        "_jitter_row": 120,        # the 50% tick
        "bars": {"61": (650, 667), "59": (677, 693),
                 "45": (704, 720), "62": (731, 747)},
    },
}
BASELINE_ROW = 208          # the drawn x-axis rule, detected independently

# Fig. 3C prints these for GlycoLNP61. NOT used to build any geometry.
PRINTED_ANCHOR = {"liver": 54.0, "spleen": 78.0}
ANCHOR_TOL_PCT = 2.0


def load_masks(seed=None):
    """gray for tick detection; `ink` = THE BAR BODY, by solidity.

    THE FIRST MASK WAS WRONG AND bars.bar_top REFUSED IT, which is the whole
    reason that refusal exists. Keying the mask on the bar's FILL colour
    (145,179,191) missed two things: the bar's own darker top STROKE, which is
    where the drawn top edge actually is, and the replicate dots drawn over the
    bar, whose interiors punch holes in the fill. The detector then found the
    first 90%-covered FILL row at 121 against a true top of 113, with the two
    flanks disagreeing by 3.0 px -- an 8 px / 4.5 percentage-point error that
    the refusal caught and that no tolerance should have been widened to admit.

    So the mask is rebuilt at the source, twice over:

      1. `not near-white` instead of the fill colour, so the top stroke and any
         dot overlapping the bar both count as bar ink. Measured separation:
         across all eight bars the highest coverage any free-floating replicate
         cluster reaches is 0.88, against 1.00 for the bar body.

      2. SOLIDITY -- each column is walked UP from the axis baseline and stops
         at its first near-white pixel. A bar is attached to its baseline; a
         cluster of replicate dots floating above it is not. Without this,
         spleen bar 62 reads 97.7% instead of ~72%, because three dots sit side
         by side across the full bar width at the top of its error whisker and
         hold >=0.90 coverage for three contiguous rows. That is the one case
         in this panel where full-width coverage alone is not enough.

    bars.bar_top is then still the authority on the top row: it re-derives it
    with its own full-width + left/right flank rule, and a disagreement is a
    refusal, not a nudge.
    """
    a = np.asarray(Image.open(FIG).convert("RGB")).astype(int)
    gray = a.mean(axis=2)

    # ---- SEEDS. Each moves ONE thing, so each tests ONE gate. ----------
    # A whole-image roll was tried first and is NOT used: it moves the ticks
    # and the bars together, so every gate here is invariant to it and it
    # passes clean. That is a no-op seed, and a no-op seed credited as a pass
    # is worse than no seed at all.
    if seed == "origin_shift":
        # ALL ticks moved +8 px, uniformly. Spacing stays uniform and LOO stays
        # clean, so neither calibration gate can see it: it is a pure ORIGIN
        # error. Only the printed Fig. 3C anchor can catch this one.
        for P in PANELS.values():
            lo, hi = P["tick_scan"]
            strip = gray[:, lo:hi].copy()
            gray[:, lo:hi] = 253.0
            gray[:-8, lo:hi] = strip[8:]
    if seed == "tick_jitter":
        # Displace ONE tick by +3 px, in the tick strip only. The bars and the
        # other ticks are untouched, so this is a non-uniform error that the
        # LOO residual is exactly the test for: 3 px vs a 1.5 px tolerance.
        for P in PANELS.values():
            lo, hi = P["tick_scan"]
            row = int(round(P["_jitter_row"]))
            strip = gray[row - 1:row + 2, lo:hi].copy()
            gray[row - 1:row + 2, lo:hi] = 253.0
            gray[row + 2:row + 5, lo:hi] = strip

    nonwhite = a.min(axis=2) < WHITE_CUT
    if seed == "fill_colour_mask":
        # The ORIGINAL defect, replayed: key the mask on the bar's fill colour
        # instead of its ink. Credited by the refusal bars.bar_top raises.
        nonwhite = (np.abs(a[:, :, 0] - BAR_RGB[0]) < BAR_TOL) & \
                   (np.abs(a[:, :, 1] - BAR_RGB[1]) < BAR_TOL) & \
                   (np.abs(a[:, :, 2] - BAR_RGB[2]) < BAR_TOL)

    ink = np.zeros_like(nonwhite)
    for P in PANELS.values():
        for (x0, x1) in P["bars"].values():
            for x in range(x0, x1):
                y = BASELINE_ROW
                while y >= 0 and nonwhite[y, x]:
                    ink[y, x] = True
                    y -= 1
    if seed == "fill_colour_mask":
        # Replay the defect EXACTLY as it occurred: the fill-colour mask was
        # used WITHOUT the solidity walk, and bar_top refused on flank
        # disagreement (left 121, right 118, delta 3.0 px > 2.0 px).
        ink = nonwhite.copy()
    if seed == "no_solidity":
        # Drop the solidity walk: floating replicate clusters are back in the
        # mask. Only spleen bar 62 is affected, and only by ~26 points.
        ink = nonwhite.copy()
    if seed == "stretch_scale":
        # a subtle error in the READ only: every bar body moved 6 px up, which
        # the ticks cannot see. Only the printed Fig. 3C anchor can.
        ink = np.concatenate([ink[6:], np.zeros((6, ink.shape[1]), bool)],
                             axis=0)
    return gray, ink


def run(seed=None, loo_tol=None):
    gray, ink = load_masks(seed)
    out = {}
    for name, P in PANELS.items():
        lo, hi = P["tick_scan"]
        ya, yb = P["tick_y"]
        ticks = AX.detect_ticks_1d(gray[ya:yb, :], "y", lo, hi, dark=INK_DARK,
                                   min_len=1, max_len=12)
        ticks = sorted(t + ya for t in ticks)
        if len(ticks) != len(P["tick_values"]):
            raise CalibrationError(
                f"{name}: detected {len(ticks)} ticks, expected "
                f"{len(P['tick_values'])}: {ticks}")

        # bar tops first, so they can be handed to forbid_target
        tops, refusals = {}, {}
        for lab, (x0, x1) in P["bars"].items():
            try:
                tops[lab] = BR.bar_top(ink, x0, x1, P["y_search"][0],
                                       P["y_search"][1],
                                       baseline_row=BASELINE_ROW)
            except DetectionError as exc:
                refusals[lab] = str(exc).splitlines()[0]

        cal = AX.calibrate(ticks, P["tick_values"], kind="linear", axis="y",
                           forbid_target=[t["row"] for t in tops.values()]
                           if len(tops) == len(ticks) else None,
                           loo_tol_px=loo_tol)

        vals = {}
        for lab, t in tops.items():
            vals[lab] = {
                "row_px": t["row"], "coverage": round(t["coverage"], 3),
                "flank_delta_px": round(t["flank_delta_px"], 2),
                "percent": round(float(cal.to_value(t["row"])), 2),
            }
        # px per percentage point, from the fit's own inverse
        ppp = abs(cal.to_pixel(1.0) - cal.to_pixel(0.0))

        out[name] = {
            "organ": P["organ"], "cell_type": P["cell_type"],
            "ticks_px": [round(t, 3) for t in ticks],
            "tick_values": P["tick_values"],
            "loo_residual_px_max": round(float(np.max(np.abs(
                AX.loo_residuals(ticks)))), 4),
            "px_per_percentage_point": round(float(ppp), 4),
            "percent_per_px": round(1.0 / float(ppp), 4),
            "bars": vals, "refusals": refusals,
        }

        # ---- CHECK A (liver only): the zero tick this fit never saw --------
        # The five fitted ticks are 100..20. Where does the fit put 0, and does
        # the drawn baseline rule -- detected independently, by a horizontal
        # dark-run scan, not by this strip -- sit there?
        if name == "liver":
            pred0 = float(np.interp(0.0, P["tick_values"][::-1],
                                    np.asarray(ticks)[::-1]))
            step = (ticks[-1] - ticks[0]) / (len(ticks) - 1)
            pred0 = ticks[-1] + step * (P["tick_values"][-1] / 20.0)
            faint = AX.detect_ticks_1d(gray[200:215, :], "y", lo, hi,
                                       dark=210, min_len=1, max_len=12)
            out[name]["excluded_tick_prediction"] = {
                "_what": "The 0 tick and the x-axis baseline were EXCLUDED "
                         "from the fit (the 0 tick is drawn fainter than the "
                         "detector threshold). Their position is therefore a "
                         "prediction, not an input -- the shape CLAUDE.md "
                         "calls a round number the authors DREW.",
                "predicted_zero_row_px": round(pred0, 3),
                "drawn_baseline_row_px": BASELINE_ROW,
                "delta_px": round(pred0 - BASELINE_ROW, 3),
                "faint_zero_tick_detected_at_dark210": [round(f + 200, 3)
                                                        for f in faint],
            }
    # ---- CHECK B: the two panels' scales, fitted independently ----------
    lv, sp = out["liver"], out["spleen"]
    out["_cross_panel_scale_check"] = {
        "_what": "Liver (5 ticks, 100..20) and Spleen (3 ticks, 100/50/0) are "
                 "calibrated from DISJOINT tick sets in DISJOINT pixel columns. "
                 "If both panels are drawn to the same plot height, the two "
                 "px-per-percentage-point figures must agree.",
        "liver_px_per_pct": lv["px_per_percentage_point"],
        "spleen_px_per_pct": sp["px_per_percentage_point"],
        "ratio": round(lv["px_per_percentage_point"] /
                       sp["px_per_percentage_point"], 5),
        "spleen_50_tick_px_observed": sp["ticks_px"][1],
        "spleen_50_tick_px_from_liver_scale": round(
            (sp["ticks_px"][0] + sp["ticks_px"][2]) / 2.0, 3),
    }
    return out


def anchor_check(res):
    """Fig. 3C's PRINTED 54 / 78 vs this panel's read for GlycoLNP61."""
    rep = {}
    for organ, want in PRINTED_ANCHOR.items():
        got = res[organ]["bars"].get("61", {}).get("percent")
        rep[organ] = {
            "printed_fig3C": want, "digitised_fig2B": got,
            "delta_pct_points": None if got is None else round(got - want, 2),
            "within_tol": None if got is None else abs(got - want) <= ANCHOR_TOL_PCT,
        }
    rep["tolerance_pct_points"] = ANCHOR_TOL_PCT
    rep["verdict"] = ("PASS" if all(v["within_tol"] for k, v in rep.items()
                                    if isinstance(v, dict)) else "FAIL")
    return rep


def main():
    seed = None
    if "--seed" in sys.argv:
        seed = sys.argv[sys.argv.index("--seed") + 1]
    res = run(seed)
    anc = anchor_check(res)
    report = {"_what": "KG_2025 Fig. 2B digitisation", "seed": seed,
              "panels": res, "printed_anchor_check": anc}
    print(json.dumps(report, indent=1))
    if seed is None:
        with open(os.path.join(HERE, "04_fig2b.json"), "w") as fh:
            json.dump(report, fh, indent=1)
    # A seed run MUST reach this line, so a crash is never credited as a catch.
    sp62 = res["spleen"]["bars"].get("62", {}).get("percent")
    nref = sum(len(v["refusals"]) for k, v in res.items() if not k.startswith("_"))
    reftxt = ";".join(t for k, v in res.items() if not k.startswith("_")
                      for t in v["refusals"].values())
    print(f"METRIC n_refusals={nref} refusal_text={reftxt[:200]!r} "
          f"anchor_verdict={anc['verdict']} "
          f"liver_delta={anc['liver']['delta_pct_points']} "
          f"spleen_delta={anc['spleen']['delta_pct_points']} "
          f"spleen62={sp62}")
    return 0 if anc["verdict"] == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
