#!/usr/bin/env python3
"""ZG_2020 -- Figure 1c's property table and the printed yield range.

Run:  python3 -I stage04_fig1c_and_yields.py

TWO INDEPENDENT PRINTED SOURCES FOR THE SAME TEN STRUCTURES, and they were
typeset separately:

  * Appendix S1 pp.2-8 prints an HRMS line per compound, as the formula OF THE
    ADDUCT ([M+Na]+ for nine, [M+H]+ for one) with a calcd m/z.
  * Figure 1c prints a property table with the NEUTRAL formula, an average
    molecular weight and a LogP.

stage01 checks the first. This file checks the second, which matters because
an error in our build would have to survive BOTH, and because Fig 1c settles
the one disagreement stage01 found. Fig 1c is a pure BITMAP -- the JATS
declares table-count=0 -- so its cells are transcribed by eye from the 709 px
figure; they are printed glyphs, and a single misread character breaks the
10/10 match below, which is what makes the transcription checkable.

AND A FREE CHECK THE AUTHORS LEFT. Figure 1a labels the ester product box
"A10-17-2Z" and the Results sentence says the route "afforded 10 compounds,
A-10 to A-17-2Z" -- but A-10..A-17-2Z is NINE compounds, and Appendix S1 gives
a full procedure for A-9 as well. The same panel prints the yield range
"77-92%". We test that range against the ten printed yields: if A-9 does not
belong to the series, the maximum is not 92.
"""
import json, os, re, sys

from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))

# Outputs go to OUTDIR, which is HERE for a normal run. The seed harnesses
# redirect it to a temp directory: a MUTATED copy must never overwrite the
# real curation outputs. It did, and the determinism check caught it --
# reference_structures.csv was left holding seed G3's REGIOISOMER anchors.
OUTDIR = HERE

# --- Figure 1c, transcribed from the figure bitmap -------------------------
FIG1C = {
    "A-9":     ("C30H54NO8P", 587.4, 1.62),
    "A-10":    ("C31H56NO8P", 601.3, 1.97),
    "A-11":    ("C32H58NO8P", 615.5, 2.30),
    "A-12":    ("C33H60NO8P", 629.4, 2.63),
    "A-13":    ("C34H62NO8P", 642.4, 2.93),
    "A-14":    ("C35H64NO8P", 657.4, 3.20),
    "A-15":    ("C36H66NO8P", 671.5, 3.46),
    "A-16":    ("C37H68NO8P", 685.4, 3.70),
    "A-17":    ("C38H70NO8P", 699.5, 3.96),
    "A-17-2Z": ("C38H66NO8P", 695.5, 3.75),
}

# --- isolated yields, transcribed from Appendix S1 pp.2-8 (native text) ----
YIELDS = {"A-9": 92, "A-10": 78, "A-11": 85, "A-12": 77, "A-13": 84,
          "A-14": 82, "A-15": 89, "A-16": 77, "A-17": 77, "A-17-2Z": 82}
FIG1A_PRINTS_RANGE = (77, 92)
FIG1A_PRINTS_LABEL = "A10-17-2Z"
TEXT_PRINTS_COUNT = 10


def main():
    enum = json.load(open(os.path.join(HERE, "03_enumerated.json")))
    built = {p["name"]: p for p in enum["products"]}
    qc = {}

    # ---- D-a: Fig 1c neutral formula vs our independent build -------------
    rows, f_ok, mw = [], 0, []
    for name, (formula, mwt, logp) in FIG1C.items():
        ours = built[name]["neutral_formula"]
        m = Chem.MolFromSmiles(built[name]["smiles"])
        our_mw = Descriptors.MolWt(m)             # AVERAGE mass
        our_logp = Crippen.MolLogP(m)
        ok = ours == formula
        f_ok += ok
        mw.append((name, mwt, our_mw))
        rows.append({"name": name, "printed_formula": formula,
                     "our_formula": ours, "formula_match": ok,
                     "printed_MW": mwt, "our_average_MW": round(our_mw, 2),
                     "MW_delta": round(our_mw - mwt, 2),
                     "printed_LogP": logp,
                     "our_Crippen_LogP": round(our_logp, 2)})
    qc["D_a_fig1c_formulas"] = {
        "n": len(FIG1C), "exact_matches": f_ok, "rows": rows,
        "_independence": "Fig 1c prints the NEUTRAL formula; Appendix S1 "
                         "prints the ADDUCT formula. Different typesetting, "
                         "different quantity, same ten molecules.",
    }

    # ---- D-b: Fig 1c settles stage01's one disagreement -------------------
    qc["D_b_A16_adjudicated"] = {
        "SI_HRMS_prints": "calcd for C37H68NO8P [M + Na]+: 708.4575",
        "problem": "that formula carries no Na although the adduct is [M+Na]+",
        "fig1c_prints_neutral": FIG1C["A-16"][0],
        "our_neutral": built["A-16"]["neutral_formula"],
        "our_Na_adduct_calcd": built["A-16"]["our_calcd"],
        "SI_printed_calcd": built["A-16"]["printed_calcd"],
        "verdict": "The SI pasted the NEUTRAL formula into the adduct slot. "
                   "Fig 1c prints that same neutral formula for A-16, and our "
                   "C37H68NNaO8P+ reproduces the SI's own printed 708.4575. "
                   "So the FORMULA STRING is the typo and both the mass and "
                   "the structure are right -- decided from two printed "
                   "sources, not from preference (the AC_2025 lesson).",
    }

    # ---- D-c: molecular weight, offset-corrected, plus the increments ---
    # The printed column sits SYSTEMATICALLY below our computed average MW
    # (the authors' own rounding / atomic-weight convention). A gate on the
    # raw delta would therefore fail all ten for a reason that is not an
    # error. So we remove the common offset and test the RESIDUAL, and we
    # separately test the printed column against ITSELF -- each CH2 step in a
    # homologous series must be ~14.03 Da, which needs nothing from us at all.
    import statistics
    deltas = [b - a for _, a, b in mw]
    offset = statistics.median(deltas)
    resid = [{"name": n, "printed": a, "ours_average": round(b, 2),
              "delta": round(b - a, 2),
              "residual_after_offset": round((b - a) - offset, 3)}
             for n, a, b in mw]
    worst = max(resid, key=lambda r: abs(r["residual_after_offset"]))

    sat_order = ["A-9", "A-10", "A-11", "A-12", "A-13", "A-14", "A-15",
                 "A-16", "A-17"]
    pm = [FIG1C[n][1] for n in sat_order]
    steps = [round(b - a, 3) for a, b in zip(pm, pm[1:])]
    CH2 = 14.027
    bad_steps = [(sat_order[i], sat_order[i + 1], st)
                 for i, st in enumerate(steps) if abs(st - CH2) > 0.6]

    # A-13 is a KNOWN, NAMED printed typo. Naming it here (rather than
    # widening the threshold) keeps the gate able to catch a NEW one.
    KNOWN_MW_TYPO = {
        "A-13": "printed 642.4; the steps into and out of it are 13.0 and "
                "15.0 Da where every other step in the homologous series is "
                "13.9-14.2, and the two errors sum to exactly 28.0 = 2 x CH2. "
                "Its FORMULA (C34H62NO8P) is right in Fig 1c, right in the "
                "SI HRMS and right in our build, and the SI's [M+H]+ calcd "
                "644.4286 reproduces exactly, so only the Fig 1c MW NUMERAL "
                "is wrong: it should read 643.4. A one-digit typo, 3 -> 2."}
    qc["D_c_molecular_weight"] = {
        "convention": "AVERAGE, not monoisotopic. The printed column sits a "
                      f"median {offset:.2f} Da below RDKit's average MW for "
                      "all ten, which is a convention offset, not an error; "
                      "the monoisotopic masses are ~0.4 Da lower still, so "
                      "the column is certainly not monoisotopic.",
        "median_offset_Da": round(offset, 3),
        "max_abs_residual_after_offset_Da": round(
            abs(worst["residual_after_offset"]), 3),
        "worst": worst,
        "per_compound": resid,
        "printed_column_CH2_steps": dict(
            zip([f"{a}->{b}" for a, b in zip(sat_order, sat_order[1:])],
                steps)),
        "steps_off_CH2_by_more_than_0.6": bad_steps,
        "known_printed_typo": KNOWN_MW_TYPO,
        "_this_check_is_internal_to_the_paper": "The step test uses only the "
            "authors' own printed column; it would find A-13 even if our "
            "build were wrong.",
    }
    unexplained = [r["name"] for r in resid
                   if abs(r["residual_after_offset"]) > 0.25
                   and r["name"] not in KNOWN_MW_TYPO]

    # ---- D-d: LogP must be monotone in chain length ----------------------
    sat = ["A-9", "A-10", "A-11", "A-12", "A-13", "A-14", "A-15", "A-16",
           "A-17"]
    vals = [FIG1C[n][2] for n in sat]
    incs = [round(b - a, 3) for a, b in zip(vals, vals[1:])]
    qc["D_d_logp_monotone"] = {
        "saturated_series": sat, "printed_logP": vals,
        "increments": incs,
        "strictly_increasing": all(i > 0 for i in incs),
        "n_increments": len(incs),
        "unsaturated_check": {
            "A-17": FIG1C["A-17"][2], "A-17-2Z": FIG1C["A-17-2Z"][2],
            "two_cis_double_bonds_lower_logP": FIG1C["A-17-2Z"][2]
                                               < FIG1C["A-17"][2]},
        "_note": "Increments fall from 0.35 to 0.24 across the series rather "
                 "than staying constant, so the authors' LogP is not a plain "
                 "fragment sum. We report that, and do not reproduce their "
                 "estimator; Crippen values are listed alongside for scale "
                 "only and are NOT asserted to match.",
    }

    # ---- D-e: the printed yield range proves A-9 is in the series --------
    ys = YIELDS
    with_a9 = (min(ys.values()), max(ys.values()))
    without = {k: v for k, v in ys.items() if k != "A-9"}
    wo = (min(without.values()), max(without.values()))
    uniq_max = [k for k, v in ys.items() if v == max(ys.values())]
    qc["D_e_yield_range"] = {
        "fig1a_prints_range_percent": list(FIG1A_PRINTS_RANGE),
        "fig1a_prints_label": FIG1A_PRINTS_LABEL,
        "text_prints_count": TEXT_PRINTS_COUNT,
        "range_over_A9_to_A17_2Z": list(with_a9),
        "range_over_A10_to_A17_2Z_only": list(wo),
        "matches_with_A9": list(with_a9) == list(FIG1A_PRINTS_RANGE),
        "matches_without_A9": list(wo) == list(FIG1A_PRINTS_RANGE),
        "compounds_at_the_maximum": uniq_max,
        "verdict": "The printed 77-92% range is reproduced ONLY when A-9 is "
                   "included, and 92% is A-9's yield and no other compound's. "
                   "So A-9 is a member of the series and the panel's own "
                   "'A10-17-2Z' label (and the Results sentence's 'A-10 to "
                   "A-17-2Z') is an off-by-one on the range START -- which "
                   "the same sentence contradicts by saying '10 compounds', "
                   "since A-10..A-17-2Z is nine.",
        "n_in_A10_to_A17_2Z": 9, "n_in_A9_to_A17_2Z": 10,
    }

    print(f"METRIC fig1c_formula_exact={f_ok}/{len(FIG1C)}")
    print(f"METRIC mw_median_offset_Da={offset:.3f}")
    print(f"METRIC mw_unexplained_residuals={json.dumps(unexplained)}")
    print(f"METRIC mw_printed_steps_off_CH2="
          f"{json.dumps([list(b) for b in bad_steps])}")
    print(f"METRIC logp_strictly_increasing="
          f"{sum(1 for i in incs if i > 0)}/{len(incs)}")
    print(f"METRIC yield_range_matches_with_A9="
          f"{qc['D_e_yield_range']['matches_with_A9']} "
          f"without_A9={qc['D_e_yield_range']['matches_without_A9']}")

    fails = []
    if f_ok != len(FIG1C):
        fails.append(f"D-a only {f_ok}/{len(FIG1C)} Fig 1c formulas match")
    if unexplained:
        fails.append(f"D-c molecular weight residual unexplained for "
                     f"{unexplained}")
    if not all(i > 0 for i in incs):
        fails.append("D-d printed LogP is not monotone in chain length")
    if not qc["D_e_yield_range"]["matches_with_A9"]:
        fails.append("D-e the printed yield range does not reproduce")
    print("METRIC fig1c_failures=" + json.dumps(fails))

    with open(os.path.join(OUTDIR, "04_fig1c_checks.json"), "w") as fh:
        json.dump(qc, fh, indent=1)
    if fails:
        raise SystemExit("GATE FAILED: " + "; ".join(fails))
    return 0


if __name__ == "__main__":
    sys.exit(main())
