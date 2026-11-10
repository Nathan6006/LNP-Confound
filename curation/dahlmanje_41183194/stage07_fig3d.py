#!/usr/bin/env python3
"""KG_2025 Fig. 3D -- the dose response, read ONLY to test Fig. 3A against it.

Run:  python3 -I curation/dahlmanje_41183194/stage07_fig3d.py <repo_root>

WHY THIS PANEL IS READ AT ALL. Nothing here is emitted: 00_value_sources.json
records Fig. 3D as admissible in kind but not emitted, because its 0.8 mg/kg
bars are the SAME condition as Fig. 3A/3C for the same two formulations, and
emitting both would put one formulation into one experiment twice.

That overlap is exactly what makes it worth measuring. CLAUDE.md's trap list:
"If two sources share a member, measure the ratio." Fig. 3A/3C and Fig. 3D
report the same quantity -- percent Thy1.1+ splenic macrophages, GlycoLNP61
and Control LNP61, 0.8 mg/kg, 16 h -- in two different cohorts. Fig. 3C PRINTS
78 and 66. If Fig. 3D disagrees, the disagreement is a property of the paper,
not of this curation, and it belongs in 08_contradictions.csv with a number
against it rather than as an impression.

Same discipline as stage02: ticks detected, scale and origin fitted from them
alone, mask built by solidity from the baseline, bars.bar_top the authority.
"""
import json
import os
import sys

import numpy as np
from PIL import Image

ROOT = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else ".")
sys.path.insert(0, os.path.join(ROOT, "curation", "lib"))
import axis as AX                              # noqa: E402
import bars as BR                              # noqa: E402
from errors import DetectionError              # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "source", "pnas.2409569122fig03.jpg")

TICK_SCAN = (54, 61)        # the strip left of the spine at x=61
TICK_Y = (370, 600)
TICK_VALUES = [100, 50, 0]
BASELINE_ROW = 576
Y_SEARCH = (370, 578)
WHITE_CUT = 225
BARS = {"Control LNP61 0.8 mpk": (71, 91),
        "GlycoLNP61 0.8 mpk": (100, 121),
        "GlycoLNP61 0.4 mpk": (130, 150),
        "GlycoLNP61 0.2 mpk": (159, 179)}

# Fig. 3C's printed numerals for the same two formulations, same organ, same
# dose. Not used in any geometry.
FIG3C_PRINTED = {"GlycoLNP61 0.8 mpk": 78.0, "Control LNP61 0.8 mpk": 66.0}


def main():
    a = np.asarray(Image.open(FIG).convert("RGB")).astype(int)
    gray = a.mean(axis=2)

    lo, hi = TICK_SCAN
    ya, yb = TICK_Y
    ticks = sorted(t + ya for t in AX.detect_ticks_1d(
        gray[ya:yb, :], "y", lo, hi, dark=245, min_len=1, max_len=12))
    if len(ticks) != 3:
        raise SystemExit(f"Fig3D: detected {len(ticks)} ticks: {ticks}")

    nonwhite = a.min(axis=2) < WHITE_CUT
    ink = np.zeros_like(nonwhite)
    for (x0, x1) in BARS.values():
        for x in range(x0, x1):
            y = BASELINE_ROW
            while y >= 0 and nonwhite[y, x]:
                ink[y, x] = True
                y -= 1

    tops, refusals = {}, {}
    for lab, (x0, x1) in BARS.items():
        try:
            tops[lab] = BR.bar_top(ink, x0, x1, Y_SEARCH[0], Y_SEARCH[1],
                                   baseline_row=BASELINE_ROW)
        except DetectionError as exc:
            refusals[lab] = str(exc).splitlines()[0]

    cal = AX.calibrate(ticks, TICK_VALUES, kind="linear", axis="y",
                       forbid_target=[t["row"] for t in tops.values()]
                       if len(tops) == 3 else None)

    vals = {lab: round(float(cal.to_value(t["row"])), 2)
            for lab, t in tops.items()}

    # ---- the comparison this panel was read for --------------------------
    cmp_ = {}
    for lab, printed in FIG3C_PRINTED.items():
        got = vals.get(lab)
        cmp_[lab] = {"fig3C_printed": printed, "fig3D_digitised": got,
                     "delta_points": None if got is None else round(got - printed, 2),
                     "ratio_3D_over_3C": None if got is None
                     else round(got / printed, 4)}
    g3c = FIG3C_PRINTED["GlycoLNP61 0.8 mpk"] / FIG3C_PRINTED["Control LNP61 0.8 mpk"]
    g3d = (vals["GlycoLNP61 0.8 mpk"] / vals["Control LNP61 0.8 mpk"]
           if "GlycoLNP61 0.8 mpk" in vals and "Control LNP61 0.8 mpk" in vals
           else None)

    out = {
        "_what": "Fig. 3D, read to test Fig. 3A/3C against it. NOT EMITTED.",
        "ticks_px": [round(t, 3) for t in ticks], "tick_values": TICK_VALUES,
        "loo_residual_px_max": round(float(np.max(np.abs(
            AX.loo_residuals(ticks)))), 4),
        "loo_tolerance_px_STRICT": AX.LOO_TOL_PX,
        "percent_per_px": round(abs(cal.to_value(1) - cal.to_value(0)), 4),
        "bar_tops": {lab: {"row": t["row"], "coverage": round(t["coverage"], 3),
                           "flank_delta_px": round(t["flank_delta_px"], 2),
                           "percent": vals[lab]}
                     for lab, t in tops.items()},
        "refusals": refusals,
        "vs_fig3C": cmp_,
        "glyco_over_control": {
            "_what": "The paper's own effect size, computed twice from two "
                     "panels reporting the same two formulations at the same "
                     "dose, organ, cell type and timepoint.",
            "from_fig3C_printed": round(g3c, 4),
            "from_fig3D_digitised": None if g3d is None else round(g3d, 4),
            "ratio_of_ratios": None if g3d is None else round(g3d / g3c, 4),
        },
        "monotone_dose_response": None,
    }
    if all(k in vals for k in ("GlycoLNP61 0.8 mpk", "GlycoLNP61 0.4 mpk",
                               "GlycoLNP61 0.2 mpk")):
        seq = [vals["GlycoLNP61 0.8 mpk"], vals["GlycoLNP61 0.4 mpk"],
               vals["GlycoLNP61 0.2 mpk"]]
        out["monotone_dose_response"] = {
            "doses_mg_per_kg": [0.8, 0.4, 0.2], "percent": seq,
            "strictly_decreasing": seq[0] > seq[1] > seq[2],
            "_what": "An internal consistency test the panel must pass and "
                     "that nothing in this curation imposed: delivery should "
                     "fall with dose.",
        }

    with open(os.path.join(HERE, "07_fig3d.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps(out, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
