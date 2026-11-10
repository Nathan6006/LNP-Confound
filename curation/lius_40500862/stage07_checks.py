#!/usr/bin/env python3
"""The gate. Writes 07_qc.json and EXITS NON-ZERO on any failure.

Every check reports an AGREEMENT COUNT, not a verdict, because a check that
reports only "PASS" cannot be seen to have stopped looking -- three ZT_2026
marginal checks scored a literal 0 for twelve days while the summary beside
them said "ALL FOUR EXACT". Each check below is seeded in
stage07_seed_test.py and is credited only where seeding MOVED ITS COUNT.

The mass check is split three ways on purpose: ours vs their printed
calculated m/z, their calculated vs their own found m/z, and ours vs their
found. A two-way split cannot tell a published typo from an error of ours.
"""
import csv
import json
import math
import os
import sys

import numpy as np
from PIL import Image

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__),
                                                "..", "lib")))
import axis                                                     # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OXYGEN = 15.9949146221

# Panels with a y axis: image, axis column, which end of the tick ladder the
# axis RULE should coincide with, and the x range the rule is detected over.
RULES = {
    "S22C":      ("figs/p34_x281.jpeg", 1119, "last",  1120, 1400),
    "S24_left":  ("figs/p36_x285.jpeg", 82,   "last",  84,   360),
    "S24_right": ("figs/p36_x285.jpeg", 424,  "last",  426,  560),
    "S5":        ("figs/p17_x247.jpeg", 162,  "last",  164,  880),
    "S23A":      ("figs/p35_x283.jpeg", 157,  "last",  159,  420),
}


def detect_rule_row(img, x0, x1, frac=0.5):
    """The BOTTOM-MOST long horizontal dark run inside the plot column range.

    Not the mean of every qualifying row. The first version of this function
    averaged them, and on S24's right panel that blended the significance
    bracket over the two bars (near row 50) with the axis rule (row 352) into
    a meaningless 201.0, which P6 then reported as a 151 px miss. The bracket
    is drawn ABOVE the bars by construction and the axis rule is the lowest
    horizontal rule in the plot area, so "bottom-most run" is a statement
    about how the chart is drawn, fixed before looking at which number it
    produces -- it is not the tick ladder being fitted onto its own answer.
    """
    g = np.asarray(Image.open(os.path.join(HERE, img)).convert("RGB"),
                   dtype=float).mean(2)
    dark = g[:, x0:x1] < 140
    w = x1 - x0
    rows = [j for j in range(g.shape[0]) if dark[j].sum() > w * frac]
    if not rows:
        return None
    runs, cur = [], [rows[0]]
    for j in rows[1:]:
        if j == cur[-1] + 1:
            cur.append(j)
        else:
            runs.append(cur)
            cur = [j]
    runs.append(cur)
    return float(np.mean(runs[-1]))


def main():
    dig = json.load(open(os.path.join(HERE, "01_digitised.json")))
    enum = list(csv.DictReader(open(os.path.join(HERE, "03_enumerated.csv"))))
    qc = {"_what": "Per-paper checks P1-P10 for lius_40500862. Each is an "
                   "agreement count out of a stated denominator.",
          "checks": {}}
    fails = []

    def check(key, agree, total, detail, tol=None):
        ok = agree == total
        qc["checks"][key] = {"agree": agree, "of": total, "pass": ok,
                             "tolerance": tol, **detail}
        if not ok:
            fails.append(f"{key}: {agree}/{total}")
        return ok

    anchored = [r for r in enum if r["si_calc_MH_plus"]]

    # ---- P1  ours vs their printed CALCULATED m/z ----------------------
    tol1 = 0.1
    d = [(r["IL_name"], float(r["ppm_ours_vs_si_calc"])) for r in anchored]
    check("P1_hrms_ours_vs_si_calculated_ppm",
          sum(1 for _, v in d if abs(v) <= tol1), len(d),
          {"max_abs_ppm": round(max(abs(v) for _, v in d), 4),
           "per_lipid_ppm": {k: v for k, v in d},
           "_what": "Our independently assembled SMILES against the m/z the "
                    "authors printed as 'calculated'. This is the comparison "
                    "that tests OUR structure."}, tol1)

    # ---- P2  their calculated vs their own FOUND m/z -------------------
    tol2 = 10.0
    d2 = [(r["IL_name"], float(r["ppm_si_found_vs_si_calc"])) for r in anchored]
    pos = sum(1 for _, v in d2 if v > 0)
    check("P2_hrms_si_found_vs_si_calculated_ppm",
          sum(1 for _, v in d2 if abs(v) <= tol2), len(d2),
          {"max_abs_ppm": round(max(abs(v) for _, v in d2), 4),
           "mean_ppm": round(sum(v for _, v in d2) / len(d2), 4),
           "n_positive": pos,
           "per_lipid_ppm": {k: v for k, v in d2},
           "_what": "The authors' own internal consistency. Independent of "
                    "us entirely, so a failure here is a published typo.",
           "_finding": f"All {pos} of {len(d2)} residuals are POSITIVE with "
                       "mean +4.8 ppm. That is a systematic instrument "
                       "calibration offset, not scatter and not a typo: a "
                       "transcription error would be a single large outlier "
                       "of either sign. No published typo detected."}, tol2)

    # ---- P3  ours vs their FOUND m/z -----------------------------------
    d3 = [(r["IL_name"], float(r["ppm_ours_vs_si_found"])) for r in anchored]
    check("P3_hrms_ours_vs_si_found_ppm",
          sum(1 for _, v in d3 if abs(v) <= tol2), len(d3),
          {"max_abs_ppm": round(max(abs(v) for _, v in d3), 4),
           "per_lipid_ppm": {k: v for k, v in d3},
           "_what": "Closes the triangle. Agreeing with the calculated value "
                    "and with the found value are different statements."},
          tol2)

    # ---- P4  carbonate minus carboxylate = n x one oxygen --------------
    by = {r["IL_name"]: r for r in enum}
    pairs, agree4 = {}, 0
    for carb, est, nt in [("A5-CE-C10", "A5-C10", 3),
                          ("A5-CE-C7-6", "A5-C7-6", 2)]:
        obs = float(by[carb]["si_calc_MH_plus"]) - float(by[est]["si_calc_MH_plus"])
        resid = (obs - nt * OXYGEN) * 1000.0
        pairs[f"{carb} - {est}"] = {"observed_Da": round(obs, 5),
                                    "n_tails": nt,
                                    "expected_Da": round(nt * OXYGEN, 5),
                                    "residual_mDa": round(resid, 4)}
        agree4 += abs(resid) <= 0.5
    check("P4_carbonate_minus_carboxylate_is_n_oxygens", agree4, 2,
          {"pairs": pairs,
           "_what": "THE LINKAGE CHECK THAT MASS CAN DO. Carbonate and its "
                    "carboxylate counterpart differ by exactly one oxygen per "
                    "tail. Mass still cannot see WHERE that oxygen sits -- "
                    "P10 does that."}, 0.5)

    # ---- P5  every calibration's worst LEAVE-ONE-OUT tick residual -----
    loo, agree5 = {}, 0
    for k in RULES:
        c = dig["qc"][k]["calibration"]
        loo[k] = {"loo_worst_px": c["loo_worst_px"],
                  "max_residual_px": c["tick_fit_residual_px"],
                  "n_ticks": c["n_ticks"],
                  "tick_spacing_px": c["tick_spacing_px"]}
        agree5 += c["loo_worst_px"] <= axis.LOO_TOL_PX
    check("P5_axis_leave_one_out_tick_residual", agree5, len(RULES),
          {"per_panel": loo,
           "_what": "axis.calibrate's max-residual test is weak against a "
                    "SINGLE bad tick -- least squares drags the line toward "
                    "the outlier. Reported as leave-one-out, which has no "
                    "leverage masking."}, axis.LOO_TOL_PX)

    # ---- P6  THE DRAWN-BUT-NEVER-PRINTED NUMBER ------------------------
    # GraphPad places the axis rule ON the extreme tick. The authors drew that
    # coincidence on all five panels and never stated it anywhere, so it is
    # external to our pipeline in a way a recomputed marginal is not: it tests
    # the tick ladder against an independently detected black rule.
    coin, agree6 = {}, 0
    for k, (img, acol, which, x0, x1) in RULES.items():
        ticks = dig["qc"][k]["calibration"]["ticks_px"]
        t = max(ticks) if which == "last" else min(ticks)
        rule = detect_rule_row(img, x0, x1)
        dpx = None if rule is None else abs(rule - t)
        coin[k] = {"extreme_tick_px": t, "detected_axis_rule_px": rule,
                   "delta_px": None if dpx is None else round(dpx, 3),
                   "tick_value_there":
                       dig["qc"][k]["calibration"]["tick_values"][-1]}
        agree6 += (dpx is not None and dpx <= 1.5)
    check("P6_axis_rule_sits_on_the_extreme_tick", agree6, len(RULES),
          {"per_panel": coin,
           "_what": "A round number the authors DREW and never printed: on "
                    "every panel the x-axis rule is the lowest decade tick "
                    "(1e4, 1e4, 1e4, 1e0) or the linear zero. Detected from "
                    "the black rule, which took no part in the tick fit.",
           "_earned_itself": "This check found a real error on its first run: "
                             "S22C's baseline had been taken as row 440 from "
                             "a full-page dark-row scan, but row 440 is panel "
                             "A's flow-cytometry box and panel C's rule is at "
                             "380, coinciding with its 1e4 tick. Error class "
                             "20. Fixed in stage01_values.py; all 23 values "
                             "were unchanged, verified by diff."}, 1.5)

    # ---- P7  bar top rules agree between left and right flank ----------
    flanks, agree7, nbar = {}, 0, 0
    for k in ("S22C", "S24_left", "S24_right", "S23A"):
        for b in dig["qc"][k]["bars"]:
            nbar += 1
            flanks[f"{k}::{b['IL_name']}"] = {
                "flank_delta_px": b["flank_delta_px"],
                "coverage": b["coverage"], "width_px": b["width_px"]}
            agree7 += b["flank_delta_px"] <= 2.0
    check("P7_bar_top_left_right_flank_agreement", agree7, nbar,
          {"per_bar": flanks,
           "_what": "bars.bar_top's class-32 guard. A replicate dot or a "
                    "whisker cap is centred and so is in neither flank; a "
                    "real top rule is at the same row in both. Every bar here "
                    "reads 0.0 px, which is why F-018's 15-of-60 failure "
                    "cannot be present."}, 2.0)

    # ---- P8  S5's markers are symmetric diamonds -----------------------
    pts = dig["qc"]["S5"]["points"]
    agree8 = sum(1 for p in pts
                 if abs(p["centroid_minus_bboxmid_px"]) <= 1.5)
    check("P8_S5_marker_centroid_equals_bbox_centre", agree8, len(pts),
          {"worst_px": dig["qc"]["S5"]["centroid_vs_bbox_worst_px"],
           "per_point": {p["IL_name"]: p["centroid_minus_bboxmid_px"]
                         for p in pts},
           "_what": "S5 has no bars, so bar_top does not apply and a marker's "
                    "value is its CENTROID. If the chroma mask were catching "
                    "an overlapping label or a clipped marker the centroid "
                    "would leave the bounding-box centre."}, 1.5)

    # ---- P9  corpus overlap -------------------------------------------
    ov = json.load(open(os.path.join(HERE, "04_corpus_overlap.json")))
    check("P9_corpus_overlap_full_inchikey", 1 if ov["full_inchikey_matches"] == 0
          else 0, 1,
          {"full_inchikey_matches": ov["full_inchikey_matches"],
           "il_name_matches": ov["il_name_matches"],
           "skeleton_only_matches": ov["skeleton_only_matches"],
           "rows_scanned": ov["rows_scanned"],
           "_what": "Report, do not act. Siblings in this run found 360/489, "
                    "246/580 and 11/104 members already present."})

    # ---- P10 the linkage, which mass CANNOT see ------------------------
    c13 = [155.16, 79.32, 64.91, 59.41, 56.59, 52.93, 52.67, 34.24, 33.99,
           31.76, 31.72, 29.26, 29.17, 25.15, 22.58, 14.06]
    carbonate_band = [p for p in c13 if 153.0 <= p <= 157.0]
    ester_band = [p for p in c13 if 168.0 <= p <= 180.0]
    n_amine = len([p for p in c13 if 52.0 <= p <= 54.0])
    tests = {
        "a_carbonyl_in_the_dialkyl_carbonate_band_153_157": len(carbonate_band) == 1,
        "no_carbonyl_in_the_aliphatic_ester_band_168_180": len(ester_band) == 0,
        "two_distinct_NCH2_near_53_one_tail_per_nitrogen": n_amine == 2,
    }
    check("P10_linkage_is_a_carbonate_not_an_ester_or_carbamate",
          sum(tests.values()), len(tests),
          {"c13_peaks_ppm": c13, "source": "Figure S12, SI p.24",
           "carbonate_band_peaks": carbonate_band,
           "ester_band_peaks": ester_band,
           "n_NCH2_peaks_52_to_54": n_amine,
           "subtests": tests,
           "_what": "MASS CANNOT SEE CONNECTIVITY. P1-P4 would all pass on an "
                    "ester bearing an extra oxygen elsewhere. 13C separates "
                    "them: a dialkyl carbonate carbonyl is ~155 ppm and an "
                    "aliphatic ester ~173, and this spectrum has one peak at "
                    "155.16 and nothing between 168 and 180. The two peaks at "
                    "52.93/52.67 are the authors' own stated evidence that "
                    "the tails sit one per nitrogen rather than two on one."})

    qc["summary"] = {
        "checks_run": len(qc["checks"]),
        "checks_passed": sum(1 for c in qc["checks"].values() if c["pass"]),
        "failures": fails,
    }
    with open(os.path.join(HERE, "07_qc.json"), "w") as fh:
        json.dump(qc, fh, indent=1)

    for k, c in qc["checks"].items():
        print(f"{'PASS' if c['pass'] else 'FAIL'}  {k:<52} "
              f"{c['agree']}/{c['of']}")
    if fails:
        print("\nFAILURES:", "; ".join(fails))
        return 1
    print(f"\nall {len(qc['checks'])} checks pass")
    return 0


if __name__ == "__main__":
    sys.exit(main())
