#!/usr/bin/env python3
"""Checks on the DIGITISED values. Writes 06_value_qc.json.

Every entry is a number, and each has a seeded counterpart in
stage06_seed.py that must move it.
"""
import json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))

# Figure 3H prints three ROI total-flux values straight onto the IVIS
# screenshot, left to right over the three labelled mice. The screenshot was
# embedded at roughly ONE PIXEL PER STROKE, so the leading digit and the
# exponent are legible and the mantissa's later decimals are NOT. Recorded at
# the precision actually readable, with the uncertainty stated, rather than
# transcribed to four figures that cannot be defended.
ROI = [("H2A4B4", "ROI 1", 3.0e7), ("A4B4-S3", "ROI 2", 1.7e8),
       ("A4B4-S4", "ROI 3", 4.1e7)]
ROI_LEGIBILITY = ("leading digit and exponent legible; 2nd-4th mantissa "
                  "digits are not resolvable at the embedded resolution "
                  "(~1 px per stroke). Read as 3.0e7 / 1.7e8 / 4.1e7 with an "
                  "uncertainty of about +/-0.1 in the leading decimal.")

MOLAR = {"ionizable": 50.0, "DSPC": 10.0, "cholesterol": 38.5,
         "C14-PEG2000": 1.5}


def main():
    f3 = json.load(open(os.path.join(HERE, "05_fig3i.json")))
    f2 = json.load(open(os.path.join(HERE, "05_fig2.json")))
    mut = sys.argv[sys.argv.index("--mutate") + 1] if "--mutate" in sys.argv else None

    bars3 = {b["lipid"]: b for b in f3["bars"]}
    if mut == "V1":
        # class 20: the origin anchored on the wrong tick -- every Fig 3I
        # value and every replicate span shifts by one tick interval.
        for b in bars3.values():
            if b.get("value_p_per_s"):
                b["value_p_per_s"] += 5.0e7
                b["inner_ink_value_span"] = [v + 5.0e7
                                             for v in b["inner_ink_value_span"]]
    if mut == "V2":
        # the LOG misreading of panel K's minor ticks, ~22x on a mid-panel bar
        f2["Fig2K"]["bars"] = [dict(b, value_p_per_s=(
            None if b.get("value_p_per_s") is None
            else b["value_p_per_s"] * 22.0)) for b in f2["Fig2K"]["bars"]]
    bars2K = {b["lipid"]: b for b in f2["Fig2K"]["bars"]}
    bars2H = {b["lipid"]: b for b in f2["Fig2H"]["bars"]}

    qc = {}

    # V1  the printed ROI values must fall inside their own bar's ink span ---
    rows = []
    for lip, tag, val in ROI:
        b = bars3[lip]
        hi, lo = b["inner_ink_value_span"]        # top-of-ink, bottom-of-ink
        rows.append({"lipid": lip, "roi": tag, "printed_value": val,
                     "bar_replicate_span": [lo, hi],
                     "inside": lo <= val <= hi,
                     "bar_mean_digitised": b["value_p_per_s"],
                     "ratio_to_bar_mean": round(val / b["value_p_per_s"], 3)})
    qc["V1_fig3H_printed_ROI_containment"] = {
        "n": len(rows), "n_inside": sum(r["inside"] for r in rows),
        "legibility": ROI_LEGIBILITY, "rows": rows,
        "what": "A number the authors DREW on the image and never discussed, "
                "external to our calibration. Each is one representative "
                "mouse, so the test is containment in the group's drawn "
                "spread -- NOT equality with the n=3 bar mean."}

    # V2  H2A4B4 is on Fig 2K and Fig 3I: two panels, two calibrations ------
    a, b = bars2K["H2A4B4"]["value_p_per_s"], bars3["H2A4B4"]["value_p_per_s"]
    qc["V2_H2A4B4_two_panels"] = {
        "Fig2K": a, "Fig3I": b, "ratio": round(a / b, 4),
        "pct_difference": round(100 * (a / b - 1), 2),
        "what": "Independent panels, independent tick sets, independent "
                "calibrations, one lipid at one dose and route. NOTE: this is "
                "NOT evidence for panel K's linear axis -- the linear reading "
                "was already forced by panel K's own four equally spaced "
                "labelled ticks and verified by its 4e7 short tick. It is a "
                "check on the two digitisations given that scale."}

    # V3  H2A4-Ole IS H2A4B2 -- the same molecule on two panels -------------
    ole = bars2K["H2A4-Ole"]["value_p_per_s"]
    b2 = bars2H["H2A4B2"]
    qc["V3_same_molecule_two_names"] = {
        "Fig2K_H2A4-Ole": ole,
        "Fig2H_H2A4B2_upper_bound": b2.get("upper_bound_p_per_s"),
        "Fig2H_H2A4B2_measured": b2.get("value_p_per_s"),
        "consistent": (b2.get("value_p_per_s") is None),
        "what": "H2A4-Ole and H2A4B2 have the same full InChIKey. On Fig 2H "
                "the bar has no resolvable height; on Fig 2K it is the second "
                "smallest bar. Both say 'at the floor', which is as much "
                "agreement as either panel can express."}

    # V4  the two panels that share H2A4B3 and H2A4B4 do NOT agree ----------
    sh = {}
    for lip in ("H2A4B3", "H2A4B4"):
        x, y = bars2H[lip].get("value_p_per_s"), bars2K[lip].get("value_p_per_s")
        sh[lip] = {"Fig2H": x, "Fig2K": y,
                   "ratio_H_over_K": round(x / y, 2) if x and y else None}
    qc["V4_fig2H_vs_fig2K"] = {
        "shared": sh,
        "what": "Both panels state 0.05 mg/kg mLuc, I.M., 6 h. The same two "
                "lipids differ by 2.8x and 25.8x between them, in the same "
                "direction, so the panels are different cohorts on different "
                "effective scales and must NOT share a z-score block without "
                "a decision. This is DECISIONS D-1."}

    # V5  a design constraint the authors never printed ---------------------
    tot = sum(MOLAR.values())
    qc["V5_molar_percent_sum"] = {
        "components": MOLAR, "sum": round(tot, 10),
        "exactly_100": abs(tot - 100.0) < 1e-9,
        "what": "The Methods print the ratio 50:10:38.5:1.5 and never say it "
                "is a mole-percent set. It sums to exactly 100, which fixes "
                "the formulation columns without a further assumption."}

    # V6  the batch doses: the caption and the Methods disagree. Adjudicate on
    #     a design constraint neither of them states.
    cap = {"batch1": (0.32, 16), "batch2": (0.08, 4)}
    met = {"batch1": (0.25, 16), "batch2": (0.80, 4)}
    per = lambda d: {k: round(v[0] / v[1], 6) for k, v in d.items()}
    pc, pm = per(cap), per(met)
    qc["V6_batch_dose"] = {
        "caption_mg_per_kg": {k: v[0] for k, v in cap.items()},
        "methods_mg_per_kg": {k: v[0] for k, v in met.items()},
        "caption_per_LNP": pc, "methods_per_LNP": pm,
        "caption_per_LNP_identical": abs(pc["batch1"] - pc["batch2"]) < 1e-9,
        "methods_per_LNP_ratio": round(pm["batch2"] / pm["batch1"], 2),
        "what": "Batch 1 pooled 16 LNPs per mouse and batch 2 pooled 4. The "
                "caption's two doses give EXACTLY the same per-LNP dose "
                "(0.02 mg/kg each); the Methods' two give 0.0156 and 0.200, a "
                "12.8x spread. A screen that changes the per-LNP dose 13-fold "
                "between rounds cannot compare them, so the caption is the "
                "self-consistent pair. No row is emitted from either batch "
                "(both are pooled mixtures), so this adjudicates a recorded "
                "contradiction, not a value."}

    txt = json.dumps(qc, indent=1)
    if mut:
        print(txt)
    else:
        open(os.path.join(HERE, "06_value_qc.json"), "w").write(txt)
        print("V1 ROI containment %d/%d" % (qc["V1_fig3H_printed_ROI_containment"]["n_inside"],
                                            qc["V1_fig3H_printed_ROI_containment"]["n"]))
        print("V2 H2A4B4 Fig2K/Fig3I ratio", qc["V2_H2A4B4_two_panels"]["ratio"])
        print("V4 Fig2H/Fig2K ratios", {k: v["ratio_H_over_K"] for k, v in sh.items()})
        print("V5 molar sum", qc["V5_molar_percent_sum"]["sum"])
        print("V6 caption per-LNP dose identical:",
              qc["V6_batch_dose"]["caption_per_LNP_identical"],
              "| methods per-LNP ratio", qc["V6_batch_dose"]["methods_per_LNP_ratio"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
