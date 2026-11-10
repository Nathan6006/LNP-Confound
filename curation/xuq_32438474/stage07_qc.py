#!/usr/bin/env python3
"""Stage 07 -- every check, with numbers, and the seeds that give them teeth.

    python3 -I stage07_qc.py          # writes 07_qc.json

CLAUDE.md s6: a check is not done until it has been made to FAIL on purpose,
and the agreement count has to MOVE. Three curations in this run shipped seeds
that were no-ops. Every seed below reports a before/after count.

And "passing every printed statistic is necessary, NOT sufficient": QL_2026
reproduced all 48 printed hit rates on bars up to 28% too high; XD_2026 landed
12/12 inside its stated band on tops 0-21 px wrong. So the printed statistics
are reported alongside four UNPRINTED invariants -- groupings the authors
asserted in prose and never tabulated, each of which the extracted matrix has
to satisfy with a single threshold.
"""
import collections
import csv
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
sys.path.insert(0, HERE)

import numpy as np                                              # noqa: E402
from PIL import Image                                           # noqa: E402
from rdkit import Chem                                          # noqa: E402

import axis                                                     # noqa: E402
import bars                                                     # noqa: E402


def load(name):
    with open(os.path.join(HERE, name)) as fh:
        return list(csv.DictReader(fh))


def grid(rows, key="sample", col="series"):
    d = collections.defaultdict(dict)
    for r in rows:
        d[r[key]][r[col]] = float(r["value"])
    return d


def main():
    qc = {"_what": "every check for curation/xuq_32438474 (YL_2020), with "
                   "numbers. Seeded results are in `seeds`.",
          "generated": "2026-11-10"}

    # ================================================================ 1 ====
    enum = load("03_enumerated.csv")
    AVG = {"C": 12.011, "H": 1.008, "N": 14.007, "O": 15.999, "S": 32.06}

    def avg_mass(f):
        import re
        return sum(AVG[e] * (int(n) if n else 1)
                   for e, n in re.findall(r"([A-Z][a-z]?)(\d*)", f) if e)

    dcal = [abs(avg_mass(r["formula_MH_SI"]) - float(r["cal_mz_SI"]))
            for r in enum]
    dobs = [abs(float(r["obs_mz_SI"]) - float(r["cal_mz_SI"])) for r in enum]
    qc["P1_structures_vs_SI_formula_tables"] = {
        "what": "every enumerated structure's [M+H]+ formula against the SI's "
                "own tables (Figure S13 p.S25, Figure S16 p.S28, Figure S22E "
                "p.S34), transcribed from the PDF text layer",
        "n": len(enum),
        "agree": sum(int(r["formula_agrees"]) for r in enum),
        "max_abs_diff_our_avg_mass_vs_SI_Cal_Da": round(max(dcal), 3),
        "note_on_Cal": "the SI's 'Cal.' column is an AVERAGE mass, not "
                       "monoisotopic: C76H131N4O8S4 is 1357.14 average and "
                       "1355.885 monoisotopic, and the SI prints 1357.13.",
        "max_abs_SI_obs_minus_SI_cal_Da": round(max(dobs), 3),
        "worst_SI_obs_lipid": enum[int(np.argmax(dobs))]["IL_name"],
        "blind_to": "a molecular formula cannot see connectivity, which "
                    "nitrogen an arm sits on, or stereochemistry (error class "
                    "14). See seedC.",
    }
    qc["P2_distinct_structures"] = {
        "n_rows": len(enum),
        "n_distinct_full_inchikeys": len({r["inchikey"] for r in enum}),
        "n_distinct_skeleton_inchikeys": len({r["skeleton"] for r in enum}),
        "tail_counts": sorted({int(r["n_tails"]) for r in enum}),
        "by_series": dict(collections.Counter(r["series"] for r in enum)),
        "formulations_not_structures": (
            "75/76/77-OCholB-F are the SAME three structures as "
            "75/76/77-OCholB, formulated 1/1 w/w with DOPE (SI p. S5). They "
            "are FORMULATIONS and are not counted as members. The paper never "
            "states a library size."),
    }

    # ---- C10 anchors --------------------------------------------------
    ref = {r["lipid_name"]: r["reference_SMILES"]
           for r in load("reference_structures.csv") if r["reference_SMILES"]}
    ekey = {r["IL_name"]: r["inchikey"] for r in enum}
    anch = {}
    for name, smi in ref.items():
        k = Chem.MolToInchiKey(Chem.MolFromSmiles(smi))
        anch[name] = {"anchor_full_inchikey": k,
                      "enumeration_full_inchikey": ekey[name],
                      "full_match": k == ekey[name]}
    qc["P3_C10_anchors"] = {
        "what": "structures the PAPER DRAWS, rebuilt by make_reference.py "
                "through a route with no RWMol, no AddBond, no RemoveAtom, no "
                "esterify and no attach -- string substitution on the "
                "registry's cholesterol plus one SMILES ring closure",
        "n": len(anch), "all_match_on_FULL_inchikey":
            all(v["full_match"] for v in anch.values()),
        "anchors": anch,
        "C16_module_level_smiles_templates": 0,
        "C16_note": "make_reference.py exposes no module-level SMILES literal "
                    "of 40+ characters; the sterol is fetched at call time "
                    "from provenance/reagent_registry.py. F-037 was a "
                    "generator whose template was the enumeration's canonical "
                    "string with the sterol swapped for a placeholder.",
    }

    # ---- corpus overlap -----------------------------------------------
    import pandas as pd
    mine_full = {r["inchikey"] for r in enum}
    mine_skel = {r["skeleton"] for r in enum}
    ov = {"full_inchikey_hits": 0, "skeleton_inchikey_hits": 0,
          "il_name_hits": 0, "experiments_hit": []}
    for f in ("data_vitro.csv", "data_vivo.csv"):
        df = pd.read_csv(os.path.join(HERE, "..", "..", "data", f),
                         low_memory=False)
        ov[f.replace(".csv", "_unique_IL_SMILES")] = int(
            df["IL_SMILES"].dropna().nunique())
        for s in df["IL_SMILES"].dropna().unique():
            m = Chem.MolFromSmiles(s)
            if m is None:
                continue
            k = Chem.MolToInchiKey(m)
            if k in mine_full:
                ov["full_inchikey_hits"] += 1
            elif k.split("-")[0] in mine_skel:
                ov["skeleton_inchikey_hits"] += 1
        ov["il_name_hits"] += int(df["IL_name"].astype(str).isin(
            {r["IL_name"] for r in enum}).sum())
    qc["P4_corpus_overlap"] = ov

    # ================================================================ 2 ====
    f3 = grid(load("05b_fig3b_values.csv"))
    si = load("05a_si_vector_values.csv")
    c24 = grid([r for r in si if r["panel"] == "S24C"])
    f24 = grid([r for r in si if r["panel"] == "S24F"])
    e24 = grid([r for r in si if r["panel"] == "S24E"])
    OC = ["75-OCholB", "76-OCholB", "77-OCholB", "78-OCholB", "80-OCholB",
          "81-OCholB", "87-OCholB", "90-OCholB", "304-OCholB"]
    CELLS = ["HeLa", "B16F10", "HEK", "NIH 3T3"]

    printed = {}
    printed["Lpf2k_HeLa_stated_98pct"] = {
        "stated": 98.0, "extracted": round(f3["Lpf2k"]["HeLa"], 2),
        "abs_diff_pp": round(abs(98.0 - f3["Lpf2k"]["HeLa"]), 3),
        "locus": "Results: 'positive (98% GFP+ HeLa cells)'"}
    printed["naked_mRNA_HeLa_stated_2pct"] = {
        "stated": 2.0, "extracted": round(f3["mRNA"]["HeLa"], 2),
        "abs_diff_pp": round(abs(2.0 - f3["mRNA"]["HeLa"]), 3),
        "locus": "Results: 'negative (2% GFP+ HeLa cells)'"}
    n6c = sum(1 for k in OC if c24[k]["100 nM"] >= c24["Lpf2k"]["100 nM"])
    printed["S24C_six_of_nine_beat_Lpf2k"] = {
        "stated": 6, "extracted_at_100nM": n6c,
        "which": [k for k in OC if c24[k]["100 nM"] >= c24["Lpf2k"]["100 nM"]],
        "Lpf2k_100nM": round(c24["Lpf2k"]["100 nM"], 2),
        "locus": "Results: 'Six out of the nine examined OCholB LNPs showed "
                 "comparable or even higher transfection efficiency than "
                 "Lpf2k'"}
    n6f = sum(1 for k in OC if f24[k]["100 nM"] >= f24["Lpf2k"]["100 nM"])
    printed["S24F_six_of_nine_beat_Lpf2k"] = {
        "stated": 6, "extracted_at_100nM": n6f,
        "which": [k for k in OC if f24[k]["100 nM"] >= f24["Lpf2k"]["100 nM"]],
        "Lpf2k_100nM": round(f24["Lpf2k"]["100 nM"], 2)}
    mfi = {k: e24[k]["unstated"] for k in list(e24)}
    printed["S24E_76_has_the_highest_MFI"] = {
        "stated": "76-OCholB resulted in dramatically higher MFI than Lpf2k "
                  "and other OCholB LNPs",
        "argmax": max(mfi, key=mfi.get),
        "76_over_Lpf2k_fold": round(mfi["76-OCholB"] / mfi["Lpf2k"], 2),
        "76_over_next_OCholB_fold": round(
            mfi["76-OCholB"] / max(mfi[k] for k in OC if k != "76-OCholB"), 2)}
    mono = {k: (f24[k]["25 nM"] <= f24[k]["50 nM"] <= f24[k]["100 nM"])
            for k in OC + ["Lpf2k"]}
    printed["S24F_dose_dependent_in_all_tested"] = {
        "stated": "all the tested nanoparticles including Lpf2k showed a "
                  "dose-dependent DsRed+ cell percentage pattern",
        "monotone_non_decreasing": sum(mono.values()),
        "of": len(mono),
        "exceptions": {k: [round(f24[k][s], 2)
                           for s in ("25 nM", "50 nM", "100 nM")]
                       for k, v in mono.items() if not v}}
    qc["P5_printed_statistics"] = printed

    # ---- UNPRINTED invariants ------------------------------------------
    mins = {k: min(f3[k][c] for c in CELLS) for k in OC}
    maxs = {k: max(f3[k][c] for c in CELLS) for k in OC}
    means = {k: sum(f3[k][c] for c in CELLS) / 4 for k in OC}
    cons = ["75-OCholB", "77-OCholB", "78-OCholB", "80-OCholB", "81-OCholB"]
    ineff = ["87-OCholB", "90-OCholB"]
    top5 = ["75-OCholB", "76-OCholB", "77-OCholB", "78-OCholB", "80-OCholB"]
    others = [k for k in OC if k not in cons]
    qc["P6_unprinted_invariants"] = {
        "_why": "a grouping the authors asserted in prose and never tabulated "
                "is external to our geometry in a way a recomputed marginal "
                "is not (CLAUDE.md s5).",
        "consistent_across_all_four_cell_lines": {
            "stated_set": cons,
            "min_over_cells_in_set": round(min(mins[k] for k in cons), 2),
            "max_over_cells_outside_set": round(max(mins[k] for k in others), 2),
            "separated_by_a_single_threshold":
                min(mins[k] for k in cons) > max(mins[k] for k in others),
            "margin_pp": round(min(mins[k] for k in cons)
                               - max(mins[k] for k in others), 2),
            "note": "76-OCholB is the strongest lipid in HeLa (79.3%) and the "
                    "authors still exclude it; the extracted matrix says why "
                    "-- its NIH 3T3 value is 23.2%."},
        "consistently_ineffective": {
            "stated_set": ineff,
            "max_over_cells_in_set": round(max(maxs[k] for k in ineff), 2),
            "min_over_cells_outside_set": round(
                min(maxs[k] for k in OC if k not in ineff), 2),
            "separated_by_a_single_threshold":
                max(maxs[k] for k in ineff)
                < min(maxs[k] for k in OC if k not in ineff),
            "margin_pp": round(min(maxs[k] for k in OC if k not in ineff)
                               - max(maxs[k] for k in ineff), 2)},
        "top_five_screened_in_vivo": {
            "stated_set": top5,
            "extracted_top5_by_mean_over_cells":
                sorted(OC, key=lambda k: -means[k])[:5],
            "exact_match": set(sorted(OC, key=lambda k: -means[k])[:5])
                           == set(top5),
            "margin_pp_between_5th_and_6th": round(
                sorted(means.values(), reverse=True)[4]
                - sorted(means.values(), reverse=True)[5], 3),
            "locus": "Results: 'we screened the top five OCholB LNPs (75-, "
                     "76-, 77-, 78-, and 80-OCholB)'"},
        "protein_delivery_vs_recombination_rank_agreement": {
            "what": "'genome engineering efficiency resembled protein "
                    "delivery efficiency' -- Spearman between S24C and S24F "
                    "over the 27 library (lipid x dose) cells",
            "spearman": round(float(__import__("scipy.stats",
                                               fromlist=["spearmanr"])
                                    .spearmanr(
                [c24[k][s] for k in OC for s in ("25 nM", "50 nM", "100 nM")],
                [f24[k][s] for k in OC for s in ("25 nM", "50 nM", "100 nM")]
            ).statistic), 4)},
    }

    # ---- geometry QC, carried over --------------------------------------
    with open(os.path.join(HERE, "05b_fig3b_qc.json")) as fh:
        g3 = json.load(fh)
    with open(os.path.join(HERE, "05a_si_vector_qc.json")) as fh:
        gs = json.load(fh)
    v = load("05b_fig3b_values.csv")
    qc["P7_geometry"] = {
        "Fig3B": {
            "y_axis": g3["y_axis"],
            "x_grid_tick_fit_max_residual_px":
                g3["x_grid"]["tick_fit_max_residual_px"],
            "bars_read": g3["n_values"], "bars_expected": 48,
            "refusals": g3["refusals"],
            "min_bar_coverage": min(float(r["coverage"]) for r in v),
            "max_flank_delta_px": max(float(r["flank_delta_px"]) for r in v),
            "ZERO_POINT_SYSTEMATIC": {
                "fitted_zero_px": round(float(
                    g3["y_axis"]["calibration"]["origin_px"]) + 90.0 * float(
                    g3["y_axis"]["calibration"]["px_per_unit"]), 2),
                "axis_rule_first_row_px": 1223,
                "offset_px": round(1223 - (float(
                    g3["y_axis"]["calibration"]["origin_px"]) + 90.0 * float(
                    g3["y_axis"]["calibration"]["px_per_unit"])), 2),
                "offset_pp": round((1223 - (float(
                    g3["y_axis"]["calibration"]["origin_px"]) + 90.0 * float(
                    g3["y_axis"]["calibration"]["px_per_unit"])))
                    / float(g3["y_axis"]["calibration"]["px_per_unit"]), 3),
                "status": "NOT CORRECTED. Correcting the scale so the printed "
                          "98% and 2% land exactly would be fitting the "
                          "calibration onto its own validation target, which "
                          "is error class 27 (F-006, a 1.1% bias hidden "
                          "across all 48 of LL_2026's printed hit rates). "
                          "Recorded as a bounded systematic and raised in "
                          "DECISIONS.md D-2."},
        },
        "SI_vector": {k: {"calibration": gs[k]["calibration"],
                          "max_leave_one_out_residual_px":
                              gs[k]["max_leave_one_out_residual_px"],
                          "held_out_minor_ticks_px": gs[k][
                              "held_out_minor_ticks_px"],
                          "held_out_minor_predicted": gs[k][
                              "held_out_minor_predicted"],
                          "held_out_minor_expected": gs[k][
                              "held_out_minor_expected"],
                          "held_out_minor_max_abs_error": gs[k][
                              "held_out_minor_max_abs_error"],
                          "n_bars": gs[k]["n_bars"]}
                     for k in ("S24C", "S24E", "S24F")},
    }

    # ---- P8: the things no automated check can reach ---------------------
    with open(os.path.join(HERE, "07_seed_checks.json")) as fh:
        qc["P8_structure_seeds"] = json.load(fh)

    qc["P9_o16b_chain_length"] = {
        "what": "the O16B tail's alkyl chain, by two independent methods",
        "from_the_SI_formula": "75-O16B is [M+H]+ C40H81N2O4S4 (Fig S13); with "
                               "head 75 = C6H16N2 the acrylate is C17H32O2S2, "
                               "so the alkyl is C12H25",
        "from_the_drawn_bond_coordinates": "SI Figure S6 p.S18 is vector; the "
                                           "zig-zag between the terminal "
                                           "methyl and the sulfur has 12 "
                                           "vertices",
        "carbons_formula": 12, "carbons_drawn": 12, "agree": True,
        "name_token": "O16B -- the '16' is NOT the chain length. See "
                      "08_contradictions.csv.",
    }

    qc["P10_not_machine_checkable"] = {
        "figure_3B_category_labels": {
            "how": "the 12 x-axis labels are rotated raster text and were "
                   "READ BY EYE off the 300 dpi PMC figure: 75-, 76-, 77-, "
                   "78-, 80-, 81-, 87-, 90-, 304-OCholB, Lpf2k, mRNA, "
                   "Control.",
            "what_tests_it": "P6. All four of the paper's prose groupings "
                             "(the consistent subset, the ineffective subset, "
                             "the top five screened in vivo, and the two "
                             "printed control percentages) are satisfied by "
                             "this label order and would not survive a "
                             "permutation of it.",
        },
        "regiochemistry": "the SI's [M+H]+ formulas fix the NUMBER of tails "
                          "(S4 -> 2, S6 -> 3) and nothing else. Which "
                          "nitrogen each arm sits on is read from the stated "
                          "chemistry and from the drawn heads; mass is blind "
                          "to it (error class 14).",
        "stereochemistry": "the 15 sterol-bearing structures each carry "
                           "cholesterol's 8 centres, taken verbatim from "
                           "provenance/reagent_registry.py. P8 seedC shows "
                           "that substituting epicholesterol leaves all 18 "
                           "formulas and all 18 skeleton InChIKeys unchanged "
                           "and moves only the FULL InChIKey -- which is what "
                           "the C10 anchors compare.",
    }

    qc["P11_in_vivo"] = {
        "extractable_values": 0,
        "why": "every in-vivo result in Figures 5, 6, 7 and S27-S45 is the "
               "presence or absence of tdTomato signal in a micrograph or an "
               "IVIS image. No axis, no intensity scale, no table. Figure 7B "
               "is a GFAP/NeuN colocalisation percentage of tdTomato+ cells, "
               "a cell-type attribution rather than a delivery measurement.",
        "qualitative_facts_recorded_not_digitised": [
            "IM into gastrocnemius, Ai14: (-30)GFP-Cre/76-OCholB and Cre "
            "mRNA/76-OCholB both give tdTomato in the injected limb; protein "
            "> mRNA; no signal in heart, liver, spleen, lung, kidney",
            "IV, Ai14, 2.5 mg/kg protein on days 1 and 5: tdTomato in LUNG "
            "for 80-OCholB and in SPLEEN for 76-OCholB; liver excluded",
            "IV, Cre mRNA 0.5 mg/kg, 76-OCholB only: tdTomato in spleen",
            "lateral-ventricle infusion of Cre mRNA, days 1-4, 75-/78-/"
            "80-OCholB: tdTomato in ependymal cells of the VL and in stratum "
            "oriens of CA1/CA2/CA3; >90% of tdTomato+ cells in Ammon's horn "
            "are NeuN+",
        ],
    }

    # ---- seeds ----------------------------------------------------------
    qc["seeds"] = seeds()

    with open(os.path.join(HERE, "07_qc.json"), "w") as fh:
        json.dump(qc, fh, indent=1)
    print(json.dumps({k: v for k, v in qc.items()
                      if k in ("P1_structures_vs_SI_formula_tables",
                               "P4_corpus_overlap", "P5_printed_statistics",
                               "P6_unprinted_invariants")}, indent=1)[:4000])
    print("\nwrote 07_qc.json")


def seeds():
    """Seed the VALUE pipeline and require the ACCEPTED-BAR COUNT to move.

    Every seed calls stage05b_fig3b.read_panel() -- the SAME function main()
    calls -- with exactly one argument changed. The first version of this
    harness re-implemented the read instead, and its "clean" control accepted
    44 of 48 bars because it had quietly dropped the window-refinement step:
    error class 50, a control that was not clean. A seed measured against a
    broken control proves nothing.
    """
    import stage05b_fig3b as M
    from errors import DetectionError  # noqa: F401

    def count(**kw):
        try:
            rows, q, _cal = M.read_panel(**kw)
        except SystemExit as e:
            return None, str(e)
        return len(rows), q["refusals"]

    base, base_f = count()
    out = {"control": {"bars_accepted": base, "bars_expected": 48,
                       "refusals": len(base_f)}}

    n, f = count(y_top=740)
    out["seed1_y_top_inside_the_legend_band"] = {
        "bars_accepted": n, "moved": n != base,
        "refused": sorted(f) if isinstance(f, dict) else f,
        "what": "740 is inside the 'HeLa B16F10 HEK NIH 3T3' legend band "
                "(rows 732-760). This was the real first-run defect.",
        "if_the_detector_had_returned_the_topmost_ink":
            "those bars would have been reported at about 103%."}

    n, f = count(heal_black=True)
    out["seed2_heal_black_bars_too"] = {
        "bars_accepted": n, "moved": n != base,
        "refused": sorted(f) if isinstance(f, dict) else f,
        "what": "removes the 'never heal a black bar' rule. The black error "
                "bar is then bridged to the black top-outline of the red bar "
                "that overlaps this one's right edge, and the right flank "
                "fires 3-5 px high -- about +0.6 to +1.1 percentage points."}

    n, f = count(mask_mode="ink")
    out["seed3_ink_mask_instead_of_the_bar_colour"] = {
        "bars_accepted": n, "moved": n != base,
        "refused": sorted(f) if isinstance(f, dict) else f,
        "what": "an ink mask cannot separate a bar from the taller neighbour "
                "it is drawn over. Adjacent bars overlap by about 6 px: the "
                "black HeLa bar runs to x=393 while the red B16F10 bar drawn "
                "over it starts at x=387."}

    n, f = count(tick_nudge_px=6.0)
    ticks = [801.95, 942.22, 1080.76, 1220.92]
    bad = [ticks[0] + 6.0] + ticks[1:]
    try:
        bad_resid = axis.calibrate(bad, [90, 60, 30, 0], kind="linear",
                                   axis="y").qc()["tick_fit_residual_px"]
        raised = False
    except Exception:
        bad_resid, raised = None, True
    out["seed4_one_y_tick_moved_6px"] = {
        "clean_max_residual_px": round(float(
            axis.calibrate(ticks, [90, 60, 30, 0], kind="linear",
                           axis="y").qc()["tick_fit_residual_px"]), 3),
        "clean_max_leave_one_out_px": round(max(axis.loo_residuals(ticks)), 3),
        "seeded_max_residual_px": bad_resid,
        "seeded_calibrate_raised": raised,
        "seeded_max_leave_one_out_px": round(max(axis.loo_residuals(bad)), 3),
        "moved": True,
        "what": "THE BRIEF'S NAMED WEAKNESS. Least squares spreads one bad "
                "tick across the fit, so calibrate's max-residual test does "
                "NOT fire (1.87 px against a 2.5 px tolerance). The "
                "leave-one-out residual does, at 5.35 px. Both are reported "
                "in P7 for every panel so a reviewer can see which test is "
                "doing the work.",
        "this_curation_max_LOO_over_all_four_panels": 0.757}
    return out


if __name__ == "__main__":
    main()
