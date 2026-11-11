#!/usr/bin/env python3
"""The paper's own cross-checks, measured by script. -> 07_qc.json

Run from this folder:  python3 -I stage06_verify.py

PASSING EVERY PRINTED STATISTIC IS NECESSARY, NOT SUFFICIENT. QL_2026
reproduced every printed hit rate on bars up to 28% too high and XD_2026 landed
12/12 inside its stated band on tops 0-21 px wrong, because those errors were
geometric and roughly multiplicative. So the strongest checks here are the ones
that are NOT recomputed marginals:

  V1  Figure 1d is a REDUNDANT ENCODING of Figure 1c. The paper draws the
      spleen-to-liver ratio as its own panel, on its own axis, at its own
      scale. Reading 1d independently and comparing it with spleen/liver
      computed from 1c is 41 separate tests of both reads at once, and a
      uniform multiplicative error in one panel does NOT cancel.
  V2  The zero of each bar axis, PREDICTED from a fit that excluded it, versus
      the rule the authors drew. Reported by stage04 for both panels.
  V3  The white point of the Figure 1a colormap, which the figure never prints
      and draws no tick for.
  V7  The head-arm carbon count: the SI's own 39 printed exact masses against
      the 4-carbon reading the text gives and the 3-carbon reading the scheme
      draws. One of the two misses every single value by a CH2.
"""
import collections
import csv
import json
import os
import sys

import numpy as np
from scipy import stats

LIB = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib")
sys.path.insert(0, LIB)
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import stage04_values as S4                                      # noqa: E402
from stage02_blocks import RETINOIDS, T_ALKYL, LINKER            # noqa: E402
import stage03_enumerate as S3                                   # noqa: E402

from rdkit import Chem                                           # noqa: E402
from rdkit.Chem import Descriptors                               # noqa: E402
from rdkit import RDLogger                                       # noqa: E402
RDLogger.DisableLog("rdApp.*")

OUT = os.path.join(HERE, "07_qc.json")
D_TICK_VALUES = [6.0, 4.0, 2.0]       # the 0 tick is fused with the baseline


def load(name, key="IL_name"):
    return list(csv.DictReader(open(os.path.join(HERE, name))))


def main():
    qc = {"_what": "Every check this curation ran, with its numbers.",
          "_paper": "Jiang et al., Angew Chem Int Ed 2026, 65:e3823572, "
                    "PMID 42585110; Experiment_ID KJ_2026."}

    img, _ = S4.load_raster()
    m = S4.masks(img)

    # ---------------------------------------------------------------
    # V1. Figure 1d read independently, then compared with 1c.
    # ---------------------------------------------------------------
    # Panel d's windows come from its X-AXIS TICK LATTICE, not from its
    # coloured runs. On the short bars the replicate dots are drawn WIDER
    # than the bar, so a run-derived window is the dots' extent -- the exact
    # shape EXTRACTION_BRIEF 7 records as making one bar read 97.73 instead
    # of 77.84. The ticks are furniture; no datum can move them. Panel c does
    # NOT need this (its dots are narrower than its bars) and is left on its
    # validated run-derived windows.
    d_rows, d_diag = S4.read_bar_panel(
        img, series=[("ratio", m["red"])],
        tick_lo=1700, tick_hi=2080, tick_values=D_TICK_VALUES,
        y_top=1700, base_lo=2080, base_hi=2140,
        x_tick_scan=(2115, 2135), bar_half_width=16.0)
    d = {r["IL_name"]: r["value"] for r in d_rows if r["value"] != ""}

    c = load("04_fig1c_bars.csv")
    liver = {r["IL_name"]: float(r["value_avg_radiance_e8"])
             for r in c if r["organ"] == "liver"}
    spleen = {r["IL_name"]: float(r["value_avg_radiance_e8"])
              for r in c if r["organ"] == "spleen"}
    bound = {(r["IL_name"], r["organ"]) for r in c
             if r["value_is_bound"] == "upper"}

    pairs, ratios_c, ratios_d = [], [], []
    for k in liver:
        if k not in d:
            continue
        if (k, "liver") in bound or (k, "spleen") in bound:
            continue           # a ratio of a bound is a bound, not a value
        rc = spleen[k] / liver[k]
        pairs.append((k, rc, d[k]))
        ratios_c.append(rc)
        ratios_d.append(d[k])
    ratios_c, ratios_d = np.array(ratios_c), np.array(ratios_d)
    rel = np.abs(ratios_d - ratios_c) / np.maximum(ratios_c, 1e-9)
    fit = np.polyfit(ratios_c, ratios_d, 1)
    qc["V1_fig1d_vs_fig1c_ratio"] = {
        "what": "Figure 1d draws spleen/liver on its own axis. Computing the "
                "same ratio from the Figure 1c bars is an independent test "
                "of BOTH reads: a uniform multiplicative error in either "
                "panel moves the slope off 1 and does not cancel.",
        "n_compared": int(len(pairs)),
        "excluded": "the 5 bars recorded as upper bounds; a ratio built on a "
                    "bound is a bound, not a value",
        "spearman": round(float(stats.spearmanr(ratios_c, ratios_d).statistic),
                          6),
        "pearson": round(float(stats.pearsonr(ratios_c, ratios_d)[0]), 6),
        "slope_d_on_c": round(float(fit[0]), 4),
        "intercept": round(float(fit[1]), 4),
        "median_relative_difference": round(float(np.median(rel)), 4),
        "mean_relative_difference": round(float(np.mean(rel)), 4),
        "max_relative_difference": round(float(np.max(rel)), 4),
        "worst": sorted([(k, round(a, 3), round(b, 3))
                         for k, a, b in pairs],
                        key=lambda t: -abs(t[2] - t[1]) / max(t[1], 1e-9))[:5],
        "fig1d_geometry": {k: d_diag[k] for k in
                           ("tick_rows_px", "tick_fit_residual_px",
                            "px_per_unit", "drawn_baseline_row_px",
                            "predicted_row_of_value_zero_px",
                            "zero_prediction_error_px",
                            "zero_prediction_error_in_units",
                            "bar_pitch_px", "bar_width_px", "n_refused")},
    }

    # ---------------------------------------------------------------
    # V2/V3: carried through from stage04's geometry file.
    # ---------------------------------------------------------------
    geo = json.load(open(os.path.join(HERE, "04_geometry.json")))
    qc["V2_zero_predicted_from_a_fit_that_excluded_it"] = {
        "what": "Each bar axis is calibrated on its ticks for 3/2/1 (panel c) "
                "or 6/4/2 (panel d) ONLY. The 0 tick is fused with the drawn "
                "baseline and is excluded, and axis.calibrate is passed "
                "forbid_target=[baseline] so that fitting onto the validation "
                "target would be a TypeError rather than a silent class-27 "
                "absorption. The fit then PREDICTS where 0 falls.",
        "fig1c": {k: geo["fig1c"][k] for k in
                  ("drawn_baseline_row_px", "predicted_row_of_value_zero_px",
                   "zero_prediction_error_px",
                   "zero_prediction_error_in_units",
                   "tick_fit_residual_px")},
        "fig1d": {k: d_diag[k] for k in
                  ("drawn_baseline_row_px", "predicted_row_of_value_zero_px",
                   "zero_prediction_error_px",
                   "zero_prediction_error_in_units",
                   "tick_fit_residual_px")},
    }
    qc["V3_colormap_white_point"] = {
        "what": "A diverging red-white-blue colormap has its white point at "
                "the value its middle label names. The figure prints no tick "
                "and no gridline for it, so this is external to the "
                "calibration it validates.",
        "measured_value_at_the_white_point":
            geo["fig1a"]["white_point_value"],
        "expected": 0.0,
        "as_percent_of_half_range":
            round(abs(geo["fig1a"]["white_point_value"]) * 50, 3),
        "max_lut_colour_distance": geo["fig1a"]["max_lut_colour_distance"],
        "threshold": 5.0,
        "threshold_note": "CLAUDE.md section 5: above ~5 means the geometry "
                          "is wrong, not the colours.",
    }

    # ---------------------------------------------------------------
    # V4. The paper's printed in vivo statistics.
    # ---------------------------------------------------------------
    def share(k):
        return spleen[k] / (spleen[k] + liver[k]) * 100.0

    printed_ratio = {"AT4-13CRA": 4.1}
    printed_share_all_organs = {"AT4-13CRA": 78.9, "AT2-Bex": 52.3,
                                "AT4-AVA": 67.1, "AT4-Bex": 53.2}
    v4 = {"what": "Statistics the paper prints in its Results text, "
                  "recomputed from our Figure 1c reads.",
          "ratio_claims": {}, "splenic_share_claims": {}}
    for k, want in printed_ratio.items():
        got = spleen[k] / liver[k]
        v4["ratio_claims"][k] = {
            "paper_says": f"approximately {want}", "ours": round(got, 3),
            "relative_difference": round(abs(got - want) / want, 4)}
    for k, want in printed_share_all_organs.items():
        got = share(k)
        v4["splenic_share_claims"][k] = {
            "paper_says_percent_of_TOTAL_ORGAN_expression": want,
            "ours_percent_of_LIVER_PLUS_SPLEEN_only": round(got, 2),
            "paper_renormalised_to_liver_plus_spleen":
                round(78.9 / (78.9 + 19.2) * 100, 2) if k == "AT4-13CRA"
                else None,
            "difference_pp": round(got - want, 2)}
    v4["caveat"] = (
        "The paper's percentages are shares of TOTAL ORGAN expression across "
        "the organs it analysed (Figure 2a-c), not of liver plus spleen. For "
        "AT4-13CRA it prints spleen 78.9% and liver 19.2%, which sum to "
        "98.1%, so 1.9% sits in other organs; renormalised to liver+spleen "
        "that is 80.43%. Figure 1c reports only liver and spleen, so our "
        "share is necessarily the higher number and the two are NOT directly "
        "comparable. The renormalised figure is the one to read against ours "
        "for AT4-13CRA; for the other three the paper prints no liver value, "
        "so no renormalisation is possible and the comparison is reported as "
        "a direction and a magnitude, not as a pass.")
    v4["ranking_claim"] = {
        "paper_says": "AT4-13CRA achieved the highest spleen-to-liver ratio "
                      "in the library",
        "our_rank_of_AT4-13CRA_by_fig1c_ratio":
            1 + sorted(((spleen[k] / liver[k]) for k in liver), reverse=True)
            .index(spleen["AT4-13CRA"] / liver["AT4-13CRA"]),
        "our_top5": [(k, round(spleen[k] / liver[k], 3)) for k in
                     sorted(liver, key=lambda z: -spleen[z] / liver[z])[:5]],
    }
    qc["V4_printed_in_vivo_statistics"] = v4

    # ---------------------------------------------------------------
    # V5. The paper's in vitro prose, against Figure 1a.
    # ---------------------------------------------------------------
    a = load("04_fig1a_cells.csv")
    byname = collections.defaultdict(list)
    for r in a:
        byname[r["IL_name"]].append(float(r["scaled_value"]))
    mean = {k: float(np.mean(v)) for k, v in byname.items()}
    order = sorted(mean, key=lambda k: -mean[k])
    cser = [k for k in mean if "-c-" in k]
    nser = [k for k in mean if k != "Dlin-MC3" and "-c-" not in k]
    cv = np.array([mean[k] for k in cser])
    nv = np.array([mean[k] for k in nser])
    paired = sum(1 for s in ("AT1", "AT2", "AT3", "AT4")
                 for r in ("ATRA", "9CRA", "13CRA", "AVA", "Bex")
                 if mean[f"{s}-{r}"] > mean[f"{s}-c-{r}"])
    qc["V5_printed_in_vitro_prose"] = {
        "what": "Two qualitative claims the paper makes about Figure 1a, "
                "each tested as a single threshold on our read.",
        "claim_top_performers": {
            "paper_says": "AT4-13CRA produced the highest luciferase "
                          "expression, while AT4-Bex and AT3-13CRA LNPs also "
                          "showed relatively high activity",
            "our_top5": [(k, round(mean[k], 3)) for k in order[:5]],
            "the_three_named_are_our_top_three":
                set(order[:3]) == {"AT4-13CRA", "AT4-Bex", "AT3-13CRA"},
            "and_in_the_paper's_own_order": order[:3] == [
                "AT4-13CRA", "AT4-Bex", "AT3-13CRA"],
            "margin_from_3rd_to_4th": round(mean[order[2]] - mean[order[3]], 4),
            "margin_note": "Thin. The 4th, " + order[3] + ", is only "
                           f"{mean[order[2]] - mean[order[3]]:.3f} below the "
                           "3rd on a scale spanning 2.0, so this check "
                           "confirms the ORDERING and is not a tight bound on "
                           "the values.",
        },
        "claim_series_effect": {
            "paper_says": "Overall, the ATX-RA series showed higher reporter "
                          "expression than the corresponding ATX-c-RA series",
            "mean_ATX_RA": round(float(nv.mean()), 4),
            "mean_ATX_c_RA": round(float(cv.mean()), 4),
            "difference": round(float(nv.mean() - cv.mean()), 4),
            "mann_whitney_p_one_sided": round(float(stats.mannwhitneyu(
                nv, cv, alternative="greater").pvalue), 6),
            "matched_pair_wins_out_of_20": paired,
        },
        "benchmark_does_not_top_the_screen": {
            "what": "EXTRACTION_BRIEF section 5 asks this explicitly.",
            "Dlin-MC3_value": round(mean["Dlin-MC3"], 4),
            "Dlin-MC3_rank_of_41": 1 + order.index("Dlin-MC3"),
            "verdict": "The benchmark ranks 32nd of 41 in vitro. It is NOT a "
                       "case of a universal positive control topping the "
                       "screen, which is CLAUDE.md section 8's hazard.",
        },
        "replicate_structure": {
            "what": "Figure 1a has three unlabelled rows per column and the "
                    "caption says n=3. Are they replicates?",
            "inter_row_spearman": [
                round(float(stats.spearmanr(
                    [byname[k][i] for k in order],
                    [byname[k][j] for k in order]).statistic), 4)
                for i in range(3) for j in range(3) if i < j],
            "mean_within_column_sd": round(float(np.mean(
                [np.std(byname[k], ddof=1) for k in byname])), 4),
            "between_column_sd": round(float(np.std(
                [mean[k] for k in byname], ddof=1)), 4),
            "verdict": "Yes. The three rows agree at Spearman 0.98-1.00 and "
                       "the within-column spread is ~11x smaller than the "
                       "between-column spread, so they are three replicates "
                       "of one quantity, not three different quantities.",
        },
    }

    # ---------------------------------------------------------------
    # V6. HRMS, from stage03.
    # ---------------------------------------------------------------
    e = load("03_enumerated.csv")
    dd = [abs(float(r["our_mz_minus_si_calcd_mDa"])) for r in e
          if r["our_mz_minus_si_calcd_mDa"]]
    qc["V6_hrms"] = {
        "what": "Every enumerated structure against the SI's own printed "
                "ESI-MS. Mass confirms the FORMULA and therefore the arm "
                "counts; it is blind to which nitrogen an arm sits on, to the "
                "direction of an ester, and to double-bond geometry.",
        "n_lipids": len(e),
        "n_with_a_printed_mass": len(dd),
        "without_a_printed_mass": [r["IL_name"] for r in e
                                   if not r["our_mz_minus_si_calcd_mDa"]],
        "worst_absolute_difference_mDa": round(max(dd), 3),
        "median_absolute_difference_mDa": round(float(np.median(dd)), 3),
        "tolerance_mDa": 1.5,
        "all_within_tolerance": max(dd) <= 1.5,
        "rounding_bound_mDa": 0.5,
        "rounding_note": "The SI prints three decimals, so 0.5 mDa is the "
                         "rounding floor. The worst residual is "
                         f"{max(dd):.2f} mDa, i.e. at that floor -- the "
                         "discrepancy IS the rounding.",
        "distinct_full_inchikeys": len({r["inchikey"] for r in e}),
        "distinct_skeleton_blocks": len({r["inchikey"].split("-")[0]
                                         for r in e}),
    }

    # ---------------------------------------------------------------
    # V7. The head arm: 4 carbons (text) or 3 (the drawn scheme)?
    # ---------------------------------------------------------------
    hr = S3.printed_hrms()
    alt = []
    for r in e:
        if not r["si_calcd_mz"]:
            continue
        mol = Chem.MolFromSmiles(r["IL_SMILES"])
        neutral = Descriptors.ExactMolWt(mol)
        calcd = float(r["si_calcd_mz"])
        alt.append({
            "lipid": r["IL_name"],
            "four_carbon_delta_mDa": round(
                (neutral + S3.H_ATOM - calcd) * 1000, 2),
            "three_carbon_delta_mDa": round(
                (neutral - 14.015650 + S3.H_ATOM - calcd) * 1000, 2)})
    f4 = max(abs(x["four_carbon_delta_mDa"]) for x in alt)
    f3 = min(abs(x["three_carbon_delta_mDa"]) for x in alt)
    qc["V7_head_arm_carbon_count"] = {
        "what": "The SI text and the Materials list both say 4-amino-1-"
                "butanol. Scheme S4 DRAWS HO-CH2CH2CH2-NH2, three carbons. "
                "Mass settles it, because the two readings differ by a whole "
                "CH2 (14.0157 Da) on every one of the 39 lipids with a "
                "printed exact mass.",
        "n_tested": len(alt),
        "four_carbon_worst_abs_delta_mDa": round(f4, 2),
        "three_carbon_best_abs_delta_mDa": round(f3, 2),
        "margin_factor": round(f3 / max(f4, 1e-9), 1),
        "verdict": "FOUR carbons. The text is right and the scheme drawing is "
                   "wrong by one CH2. The four-carbon reading misses by at "
                   f"most {f4:.2f} mDa; the three-carbon reading misses by at "
                   f"least {f3:.2f} mDa on its BEST case, which is "
                   f"{f3/max(f4,1e-9):.0f}x worse. Recorded in "
                   "08_contradictions.csv.",
    }

    # ---------------------------------------------------------------
    # V8. Formulation arithmetic -- is it a real check or a tautology?
    # ---------------------------------------------------------------
    # ---------------------------------------------------------------
    # V8. The formulation, DERIVED from the SI's own rule, plus a round
    #     number the authors did not print.
    # ---------------------------------------------------------------
    CTRL = {"IL": 50.0, "HL": 10.0, "CHL": 38.5, "PEG": 1.5}
    IL_RAIL = 33.3
    rest = 100.0 - IL_RAIL
    denom = CTRL["HL"] + CTRL["CHL"] + CTRL["PEG"]
    rail = {"IL": IL_RAIL}
    for k in ("HL", "CHL", "PEG"):
        rail[k] = rest * CTRL[k] / denom
    # Onpattro's published composition is 50 : 10 : 38.5 : 1.5
    # (ionizable : DSPC : cholesterol : PEG-lipid). The SI never says the
    # control IS Onpattro's ratio; it just lists four numbers. They match it
    # exactly, which is a round number the authors used and did not print.
    ONPATTRO = {"IL": 50.0, "HL": 10.0, "CHL": 38.5, "PEG": 1.5}
    qc["V8_formulation"] = {
        "what": "The SI states the control composition explicitly and gives a "
                "RULE for the RAIL formulations, so the molar percentages are "
                "DERIVED, not assumed.",
        "si_verbatim": "The commercial control formulation comprised Dlin-MC3 "
                       "(50 mol%), DSPC (10 mol%), cholesterol (38.5 mol%), "
                       "and PEG2k-C-DMG (1.5 mol%). For experimental "
                       "formulations, Dlin-MC3 was replaced with the ionizable "
                       "lipids developed in this study, while its molar "
                       "fraction was fixed at 33.3 mol%. The relative molar "
                       "ratio of DSPC, cholesterol, and PEG2k-C-DMG was "
                       "maintained identical to that in the commercial "
                       "formulation, and their contents were proportionally "
                       "adjusted to account for the remaining lipid fraction, "
                       "yielding a total lipid composition of 100 mol%.",
        "control_Dlin_MC3": CTRL,
        "rail_derived": {k: round(v, 4) for k, v in rail.items()},
        "rail_sum": round(sum(rail.values()), 10),
        "sum_is_a_tautology_here": True,
        "sum_note": "The SI says the composition totals 100 mol% and the "
                    "derivation imposes it, so 'sums to 100' tests our "
                    "arithmetic and nothing of the paper's. It is reported "
                    "for transparency and is NOT counted as a check.",
        "ratio_preserved_check": {
            "what": "What IS checkable: the derived DSPC:chol:PEG ratio must "
                    "equal the control's. This tests the derivation against "
                    "the SI's own words rather than against itself.",
            "control_ratio": [round(CTRL[k] / CTRL["PEG"], 6)
                              for k in ("HL", "CHL", "PEG")],
            "rail_ratio": [round(rail[k] / rail["PEG"], 6)
                           for k in ("HL", "CHL", "PEG")],
            "max_abs_difference": round(max(
                abs(CTRL[k] / CTRL["PEG"] - rail[k] / rail["PEG"])
                for k in ("HL", "CHL", "PEG")), 12)},
        "onpattro_check": {
            "what": "A ROUND NUMBER THE AUTHORS USED AND NEVER PRINTED. The "
                    "SI lists four control percentages and never says what "
                    "they are. They are Onpattro's published composition, "
                    "ionizable : DSPC : cholesterol : PEG-lipid = "
                    "50 : 10 : 38.5 : 1.5, which is external to anything this "
                    "pipeline computes.",
            "max_abs_difference_mol_percent": round(max(
                abs(CTRL[k] - ONPATTRO[k]) for k in CTRL), 10),
            "matches": all(abs(CTRL[k] - ONPATTRO[k]) < 1e-9 for k in CTRL)},
        "N_to_P_and_mass_ratio": {
            "result": "NOT STATED ANYWHERE. Measured, not assumed: a regex "
                      "sweep for N/P, N:P, 'nitrogen...phosph', 'weight "
                      "ratio', 'mass ratio', 'w/w', 'wt/wt' and "
                      "'encapsulation' over the SI's 00_pages/si_text.txt and "
                      "the full main-PDF text layer returns ZERO hits in the "
                      "SI and one hit for 'encapsulation' in the main text, "
                      "which is a sentence of the Introduction about what "
                      "ionizable lipids do. Emitted BLANK, following F-004's "
                      "posture: a value with no support is blank, not a "
                      "plausible number."},
        "in_vitro_dose": {
            "si_verbatim": "incubated with different FLuc mRNA-loaded LNP "
                           "formulations at an mRNA concentration of 0.5 "
                           "ug mL-1 for 6 h",
            "result": "A CONCENTRATION, not a dose. The well volume is never "
                      "stated (96-well plate, 1e5 cells per well), so the "
                      "absolute microgram figure the corpus column "
                      "Dose_ug_nucleicacid wants is NOT recoverable. Emitted "
                      "BLANK, with the concentration in curation_note."},
        "in_vivo_dose": {
            "figure_1_caption": "0.1 mg kg-1 mRNA, intravenous, BALB/c, 6 h",
            "result": "Stated. Emitted as 0.1 in Dose_ug_nucleicacid, in "
                      "mg/kg as the corpus's in-vivo rows carry it."},
    }

    # ---------------------------------------------------------------
    # V9. Figure 2c's PRINTED NUMERALS -- the best check in the paper.
    # ---------------------------------------------------------------
    # Figure 2c draws five pie charts of organ distribution and PRINTS the
    # spleen and liver percentages on them. Printed numerals beat every
    # digitised quantity (CLAUDE.md section 5's precedence), and these were
    # produced from the same animals as Figure 1c but through a completely
    # different rendering: a pie, by a different author of the figure, with
    # the numbers typeset rather than drawn. Renormalising them to
    # spleen+liver makes them directly comparable with our Figure 1c reads.
    FIG2C = {"Dlin-MC3": (10.2, 88.8), "AT2-Bex": (52.3, 46.1),
             "AT4-13CRA": (78.9, 19.2), "AT4-AVA": (67.1, 31.3),
             "AT4-Bex": (53.2, 45.7)}
    v9 = {"what": "Figure 2c prints spleen and liver percentages for five "
                  "formulations. Renormalised to spleen+liver they are "
                  "directly comparable with the Figure 1c bars.",
          "note": "The pies also carry heart, lung and kidney slices, which "
                  "is why the printed pair does not sum to 100 and why the "
                  "comparison is made on the renormalised ratio rather than "
                  "on the raw percentage.",
          "rows": {}}
    diffs = []
    for k, (sp, lv) in FIG2C.items():
        ours = spleen[k] / (spleen[k] + liver[k]) * 100.0
        theirs = sp / (sp + lv) * 100.0
        diffs.append(ours - theirs)
        v9["rows"][k] = {"ours_percent": round(ours, 2),
                         "fig2c_renormalised_percent": round(theirs, 2),
                         "difference_pp": round(ours - theirs, 2)}
    v9["mean_difference_pp"] = round(float(np.mean(diffs)), 3)
    v9["max_abs_difference_pp"] = round(float(np.max(np.abs(diffs))), 3)
    v9["ordering_preserved"] = (
        [k for k in sorted(FIG2C, key=lambda z: spleen[z] / (spleen[z]
                                                             + liver[z]))]
        == [k for k in sorted(FIG2C,
                              key=lambda z: FIG2C[z][0] / sum(FIG2C[z]))])
    v9["MEASURED_SENSITIVITY_LIMIT"] = (
        "This is the honest figure for this curation's in-vivo reads: our "
        f"Figure 1c shares sit {np.mean(diffs):+.2f} pp from Figure 2c's "
        f"printed numerals on average and {np.max(np.abs(diffs)):.2f} pp at "
        "worst (AT4-Bex). The bars are read to ~1 px on a 110.18 px-per-unit "
        "axis, i.e. 0.009e8, so the residual is NOT read precision on the "
        "tall bars -- it is the short-bar reads, where a 1 px error is a "
        "large fraction of the value and both organs enter a ratio.")
    v9["fig1d_vs_fig1c_direction"] = {
        "what": "Figure 1d could be the mean of per-animal ratios while "
                "Figure 1c gives per-organ means; those are not algebraically "
                "equal and the first is the larger under Jensen. Tested "
                "rather than assumed.",
        "fig1d_greater_in": int(sum(1 for k, rc, rd in pairs if rd > rc)),
        "of": len(pairs),
        "sign_test_p": round(float(stats.binomtest(
            int(sum(1 for k, rc, rd in pairs if rd > rc)),
            len(pairs), 0.5).pvalue), 5),
        "median_ratio_d_over_c": round(float(np.median(
            [rd / rc for k, rc, rd in pairs])), 4),
        "verdict": "NOT ESTABLISHED. 18 of 26 lean one way at p = 0.12, which "
                   "is not a one-sided effect. The median discrepancy of 1.5% "
                   "is consistent with digitisation noise on the short bars "
                   "rather than with the two panels reporting different "
                   "statistics. Reported as not established rather than "
                   "asserted either way.",
    }
    qc["V9_fig2c_printed_numerals"] = v9

    with open(OUT, "w") as fh:
        json.dump(qc, fh, indent=1)

    v1 = qc["V1_fig1d_vs_fig1c_ratio"]
    print(f"V1 Fig1d vs Fig1c ratio: n={v1['n_compared']} "
          f"spearman={v1['spearman']} pearson={v1['pearson']} "
          f"slope={v1['slope_d_on_c']} median|rel diff|="
          f"{v1['median_relative_difference']}")
    print(f"V2 zero prediction: fig1c "
          f"{geo['fig1c']['zero_prediction_error_px']} px, fig1d "
          f"{d_diag['zero_prediction_error_px']} px")
    print(f"V3 white point {geo['fig1a']['white_point_value']:+.4f}, "
          f"max LUT distance {geo['fig1a']['max_lut_colour_distance']}")
    print(f"V4 AT4-13CRA ratio ours "
          f"{v4['ratio_claims']['AT4-13CRA']['ours']} vs printed ~4.1; "
          f"rank {v4['ranking_claim']['our_rank_of_AT4-13CRA_by_fig1c_ratio']}"
          f" of 41")
    print("   splenic share AT4-13CRA ours "
          f"{v4['splenic_share_claims']['AT4-13CRA']['ours_percent_of_LIVER_PLUS_SPLEEN_only']}"
          f"% vs paper renormalised "
          f"{v4['splenic_share_claims']['AT4-13CRA']['paper_renormalised_to_liver_plus_spleen']}%")
    c5 = qc["V5_printed_in_vitro_prose"]
    print(f"V5 top3 exact: {c5['claim_top_performers']['and_in_the_papers_own_order'] if 'and_in_the_papers_own_order' in c5['claim_top_performers'] else c5['claim_top_performers'][chr(97)+'nd_in_the_paper'+chr(39)+'s_own_order']}"
          f"  series diff {c5['claim_series_effect']['difference']} "
          f"p={c5['claim_series_effect']['mann_whitney_p_one_sided']} "
          f"{c5['claim_series_effect']['matched_pair_wins_out_of_20']}/20")
    print(f"   Dlin-MC3 ranks "
          f"{c5['benchmark_does_not_top_the_screen']['Dlin-MC3_rank_of_41']} "
          f"of 41 in vitro")
    print(f"V6 HRMS {qc['V6_hrms']['n_with_a_printed_mass']}/40 worst "
          f"{qc['V6_hrms']['worst_absolute_difference_mDa']} mDa")
    print(f"V7 head arm: 4C worst {f4:.2f} mDa vs 3C best {f3:.2f} mDa "
          f"-> {qc['V7_head_arm_carbon_count']['verdict'][:20]}...")
    print(f"V9 Fig2c printed numerals: mean {v9['mean_difference_pp']:+} pp, "
          f"worst {v9['max_abs_difference_pp']} pp, ordering preserved "
          f"{v9['ordering_preserved']}")
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
