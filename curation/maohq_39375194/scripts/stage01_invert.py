#!/usr/bin/env python3
"""Stage 01 -- invert the one-hot, prove the inversion lossless, and check the
repository file against everything outside it.

RUN:  python3 -I scripts/stage01_invert.py <1080_LNP_OHE.csv> <Compiled_validation.xlsx> <out.json>

Nothing here digitises anything. Every value comes from the authors' repository
(00_value_sources.json tier 3). What this file does is REFUSE to believe the
file about itself: each check compares the file to something that is not the
file -- the SI's printed levels, the drawn structures, the other parameterisation
of the same quantity, or the authors' own pipeline code.

WHY THE SUM-TO-100 CHECK IS NOT HERE. The obvious "round number the authors
implied": the four molar percents sum to 100. But with

    HL  = a*h/100          IL  = a*(100-h)/100
    PEG = (100-a)*p/100    CHL = (100-a)*(100-p)/100

the sum is a + (100-a) = 100 identically, for ANY a, h, p -- including garbage.
It is a tautology of the derivation, not a property of the data, so it cannot
fail and is not credited as a check. CLAUDE.md section 6: a check that cannot
fail is not a check. It is still ASSERTED below (sum_to_100_is_tautological) so
a later reader does not mistake its absence for an oversight, and it is labelled
as what it is.

THE REAL ROUND-NUMBER CHECK is check 5: the file carries the SAME quantity in
two parameterisations -- an integer lipid:lipid ratio and a percentage -- and
the percentages are the exact images 100/(1+r) of round integers r. That is
error class 40's detector: a row-wise derivation, never a range test.
"""
import csv
import json
import math
import os
import sys

from openpyxl import load_workbook


# ---------------------------------------------------------------- loading
def load_training(path):
    with open(path, newline="", encoding="utf-8-sig") as fh:
        return list(csv.DictReader(fh))


def load_validation(path):
    wb = load_workbook(path, data_only=True)
    ws = wb.worksheets[0]
    it = ws.iter_rows(values_only=True)
    hdr = [str(h).strip() for h in next(it)]
    out = []
    for row in it:
        if all(v is None for v in row):
            continue
        out.append(dict(zip(hdr, row)))
    return out


def f(x):
    if x is None or x == "":
        return None
    return float(x)


# ------------------------------------------------------- the helper lipids
# Structures READ OFF the drawn depictions in SI Table 2 (page 2), each image
# extracted at its native embedded resolution (220-336 dpi) and paired with its
# row label by y-overlap geometry, not by assuming table order. See
# renders/p2_x6*.png and 07_qc.json check 'drawn_structure_read'.
HELPER_SMILES = {
    # 18:1/18:1 TAP, chloride.  SI p2 xref 64.
    "DOTAP": "CCCCCCCC/C=C\\CCCCCCCC(=O)OC[C@@H](COC(=O)CCCCCCC/C=C\\CCCCCCCC)C[N+](C)(C)C.[Cl-]",
    # dimethyldioctadecylammonium bromide.  SI p2 xref 65.
    "DDAB": "CCCCCCCCCCCCCCCCCC[N+](C)(C)CCCCCCCCCCCCCCCCCC.[Br-]",
    # 18:1/18:1 PE, zwitterion.  SI p2 xref 66.
    "DOPE": "CCCCCCCC/C=C\\CCCCCCCC(=O)O[C@H](COC(=O)CCCCCCC/C=C\\CCCCCCCC)COP(=O)([O-])OCC[NH3+]",
    # 18:0/18:0 PC, zwitterion.  SI p2 xref 67.
    "DSPC": "CCCCCCCCCCCCCCCCCC(=O)O[C@H](COC(=O)CCCCCCCCCCCCCCCCC)COP(=O)([O-])OCC[N+](C)(C)C",
    # DRAWN AS DIOLEOYL (18:1/18:1) PG -- see check 'reagent_18PG_conflict'.
    # SI p2 xref 68.  BOTH tails carry a cis double bond.
    "18PG": "CCCCCCCC/C=C\\CCCCCCCC(=O)O[C@H](COC(=O)CCCCCCC/C=C\\CCCCCCCC)COP(=O)([O-])OCC(O)CO",
    # 14:0/14:0 PA, sodium salt drawn.  SI p2 xref 69.
    "14PA": "CCCCCCCCCCCCCC(=O)O[C@H](COC(=O)CCCCCCCCCCCCC)COP(=O)([O-])O",
}

# SI Table 3, transcribed by hand from the page-3 text. The file must match
# these; it is the only place the paper prints the descriptors.
SI_TABLE3 = {
    #              P  N  cLogP  Hbond_A Hbond_D  Ncarb  dbonds
    "DOTAP": (1, 0, 12.5, 4, 0, 36, 2),
    "DDAB": (1, 0, 12.7, 0, 0, 36, 0),
    "DOPE": (1, 1, 9.8, 9, 2, 36, 2),
    "DSPC": (1, 1, 6.1, 8, 0, 36, 0),
    "14PA": (0, 1, 11.9, 8, 2, 28, 0),
    "18PG": (0, 1, 8.82, 10, 3, 36, 2),
}

# SI Table 1, transcribed by hand from the page-2 text.
SI_TABLE1 = {
    "HL_(IL+HL)": [0.50, 0.99, 1.96, 9.09, 50],
    "(IL+HL)": [20, 40, 60, 80],
    "PEG_(Chol+PEG)": [0.20, 0.99, 9.09],
    "NP_ratio": [4, 8, 12],
}

OHE_COLS = ["18PG", "14PA", "DOPE", "DSPC", "DOTAP", "DDAB"]
CELLS = ["HepG2", "N2a", "ARPE19", "HEK293", "B16", "PC3"]
RLU = ["RLU_" + c for c in CELLS]
DESC = ["P_charged_centers", "N_charged_centers", "cLogP", "Hbond_A",
        "Hbond_D", "Total_Carbon_Tails", "Double_bonds"]


def main(train_path, val_path, out_path):
    tr = load_training(train_path)
    va = load_validation(val_path)
    Q = {}

    # ========================================================= CHECK 1
    # THE ONE-HOT INVERSION, and its losslessness.
    #
    # Two conditions, and they are different. (a) Every row has EXACTLY one 1
    # across the six helper-lipid columns -- two or zero is a parse error, not a
    # formulation. (b) The inversion AGREES with the file's own categorical
    # Helper_lipid column, which makes the map a verified bijection rather than
    # a relabelling we invented. (b) is what makes this provable: without an
    # independent categorical column, (a) alone cannot tell a correct inversion
    # from a permuted one.
    bad_sum, disagree, inverted = [], [], []
    for i, r in enumerate(tr):
        ones = [c for c in OHE_COLS if r[c].strip() == "1"]
        zeros_ok = all(r[c].strip() in ("0", "1") for c in OHE_COLS)
        s = sum(int(r[c]) for c in OHE_COLS)
        if s != 1 or len(ones) != 1 or not zeros_ok:
            bad_sum.append((i + 2, s, ones))
            inverted.append(None)
            continue
        lvl = ones[0]
        inverted.append(lvl)
        if lvl != r["Helper_lipid"].strip():
            disagree.append((i + 2, lvl, r["Helper_lipid"]))
    Q["onehot_inversion"] = {
        "rows": len(tr),
        "group": "helper lipid, 6 columns: " + ", ".join(OHE_COLS),
        "rows_whose_group_sums_to_exactly_1": len(tr) - len(bad_sum),
        "rows_with_two_or_zero_ones": len(bad_sum),
        "offending_rows": bad_sum[:10],
        "rows_where_inversion_equals_the_file's_Helper_lipid_column":
            len(tr) - len(bad_sum) - len(disagree),
        "disagreements": disagree[:10],
        "lossless": len(bad_sum) == 0 and len(disagree) == 0,
        "levels_recovered": sorted(set(x for x in inverted if x)),
        "note": ("Lossless in BOTH directions: Helper_lipid -> one-hot is "
                 "total and one-hot -> Helper_lipid is total, so the six "
                 "columns carry no information the categorical column does "
                 "not. CLAUDE.md section 3 forbids storing one-hot columns in "
                 "the corpus; this is the proof that dropping them loses "
                 "nothing."),
    }

    # ========================================================= CHECK 2
    # Is a row a distinct ionizable lipid, or a formulation of one?
    ils = sorted(set(r["Ionizable_Lipid"].strip() for r in tr))
    hls = sorted(set(r["Helper_lipid"].strip() for r in tr))
    Q["what_a_row_is"] = {
        "distinct_Ionizable_Lipid_values": len(ils),
        "Ionizable_Lipid_values": ils,
        "distinct_Helper_lipid_values": len(hls),
        "Helper_lipid_values": hls,
        "verdict": ("FORMULATIONS. All %d rows carry the SAME single ionizable "
                    "lipid, %s. The library varies the helper lipid identity "
                    "and four compositional parameters. There is exactly ONE "
                    "distinct ionizable-lipid structure in this paper."
                    % (len(tr), ils[0] if len(ils) == 1 else "?")),
    }

    # ========================================================= CHECK 3
    # Full factorial, and the stated library size as the product of levels.
    factors = ["Helper_lipid", "HL_(IL+HL)", "(IL+HL)", "PEG_(Chol+PEG)",
               "NP_ratio"]
    levels = {c: sorted(set(r[c].strip() for r in tr)) for c in factors}
    product = 1
    for c in factors:
        product *= len(levels[c])
    cells = {}
    for r in tr:
        k = tuple(r[c].strip() for c in factors)
        cells[k] = cells.get(k, 0) + 1
    Q["factorial"] = {
        "n_levels": {c: len(levels[c]) for c in factors},
        "product_of_level_counts": product,
        "rows_in_file": len(tr),
        "product_equals_rows": product == len(tr),
        "distinct_factor_combinations": len(cells),
        "combinations_appearing_exactly_once":
            sum(1 for v in cells.values() if v == 1),
        "combinations_appearing_more_than_once":
            sum(1 for v in cells.values() if v > 1),
        "complete_and_unreplicated":
            len(cells) == product and all(v == 1 for v in cells.values()),
        "note": ("6 helper lipids x 5 HL/(IL+HL) x 4 (IL+HL) x 3 "
                 "PEG/(Chol+PEG) x 3 N/P = 1080. The stated library size "
                 "equals the product of its level counts exactly, and every "
                 "cell is occupied exactly once -- a complete, unreplicated "
                 "full factorial."),
    }

    # ========================================================= CHECK 4
    # SI Table 1's printed levels vs the file's levels.
    t1 = {}
    for col, printed in SI_TABLE1.items():
        got = sorted(set(round(float(r[col]), 6) for r in tr))
        # printed values are rounded to 2 dp in the SI
        ok = len(got) == len(printed) and all(
            abs(round(g, 2) - p) < 5e-3 for g, p in zip(got, sorted(printed)))
        t1[col] = {"printed_in_SI_Table_1": sorted(printed),
                   "in_file": got, "n_levels_in_file": len(got),
                   "n_levels_printed": len(printed), "agrees": ok}
    Q["si_table1_levels"] = t1
    Q["si_table1_all_agree"] = all(v["agrees"] for v in t1.values())

    # ========================================================= CHECK 5
    # ERROR CLASS 40 -- the same quantity in two parameterisations.
    # The file carries BOTH an integer lipid:lipid ratio and a percentage for
    # two quantities. A range check passes on a permutation, so this is tested
    # ROW BY ROW: pct == 100/(1+ratio) on every single row.
    pairs = [("Dlin-MC3_Helper lipid_ratio", "HL_(IL+HL)"),
             ("Chol_DMG-PEG_ratio", "PEG_(Chol+PEG)")]
    e40 = {}
    for ratio_col, pct_col in pairs:
        resid, worst, nbad = [], None, 0
        for i, r in enumerate(tr):
            ratio = float(r[ratio_col])
            pct = float(r[pct_col])
            pred = 100.0 / (1.0 + ratio)
            d = abs(pred - pct)
            resid.append(d)
            if worst is None or d > worst[0]:
                worst = (d, i + 2, ratio, pct, pred)
            if d > 1e-6:
                nbad += 1
        # the ratios themselves must be round integers -- that is the round
        # number the authors never printed as such
        rs = sorted(set(float(r[ratio_col]) for r in tr))
        e40[pct_col] = {
            "derived_from": ratio_col,
            "relation": "pct = 100/(1+ratio)",
            "rows_checked": len(tr),
            "rows_violating_at_1e-6": nbad,
            "max_abs_residual_pct": worst[0],
            "worst_row": {"csv_row": worst[1], "ratio": worst[2],
                          "pct_in_file": worst[3], "pct_derived": worst[4]},
            "ratio_levels": rs,
            "ratio_levels_are_exact_integers":
                all(abs(x - round(x)) < 1e-12 for x in rs),
            "note": ("ROUND NUMBER THE AUTHORS NEVER PRINTED: the SI prints "
                     "only the percentages (0.50, 0.99, 1.96, 9.09, 50 and "
                     "0.20, 0.99, 9.09), which look arbitrary. They are the "
                     "exact images 100/(1+r) of the round integer ratios "
                     "%s. Agreement to 1e-6 on every row is an external, "
                     "cheap check that the two columns were not shuffled "
                     "relative to each other -- the F-004 failure mode."
                     % rs),
        }
    Q["error_class_40_two_parameterisations"] = e40

    # ========================================================= CHECK 6
    # Helper-lipid descriptor block: each descriptor must be a FUNCTION of the
    # helper lipid (constant within level), and must match SI Table 3.
    by_hl = {}
    for r in tr:
        by_hl.setdefault(r["Helper_lipid"].strip(), []).append(r)
    nonconst, t3 = [], {}
    for hl, rows in sorted(by_hl.items()):
        vals = {}
        for c in DESC:
            u = set(rows[i][c].strip() for i in range(len(rows)))
            if len(u) != 1:
                nonconst.append((hl, c, sorted(u)))
            vals[c] = float(rows[0][c])
        exp = SI_TABLE3[hl]
        got = (vals["P_charged_centers"], vals["N_charged_centers"],
               vals["cLogP"], vals["Hbond_A"], vals["Hbond_D"],
               vals["Total_Carbon_Tails"], vals["Double_bonds"])
        # SI prints cLogP to 3 s.f.; the rest are integers
        ok = (all(abs(a - b) < 1e-9 for a, b in
                  zip([got[0], got[1], got[3], got[4], got[5], got[6]],
                      [exp[0], exp[1], exp[3], exp[4], exp[5], exp[6]]))
              and abs(got[2] - exp[2]) <= 0.005 * max(1.0, abs(exp[2])) + 0.005)
        t3[hl] = {"file": got, "si_table3": exp, "agrees": ok,
                  "n_rows": len(rows)}
    Q["helper_descriptor_block"] = {
        "constant_within_helper_lipid": len(nonconst) == 0,
        "non_constant": nonconst,
        "vs_si_table3": t3,
        "all_agree_with_si_table3": all(v["agrees"] for v in t3.values()),
        "rows_per_helper_lipid": {k: len(v) for k, v in sorted(by_hl.items())},
    }

    # ========================================================= CHECK 7
    # Measured coverage: which (formulation, cell) cells carry a number.
    by_hl_names = set(r["Helper_lipid"].strip() for r in tr)
    miss = {c: [] for c in RLU}
    for i, r in enumerate(tr):
        for c in RLU:
            if r[c].strip() == "":
                miss[c].append(r["Helper_lipid"].strip())
    cov, partial = {}, []
    for c in RLU:
        hl_missing = sorted(set(miss[c]))
        # STRUCTURAL invariant, not a frozen literal: within one cell line a
        # helper lipid is either entirely measured or entirely absent. A single
        # blanked well breaks it. Added after the seeded-mutation run: the
        # 'measurement_blanked' seed moved total_measured 5760 -> 5759 while
        # the gate stayed silent, so the gate had no coverage here at all.
        for hl in sorted(by_hl_names):
            n_hl = sum(1 for r in tr if r["Helper_lipid"].strip() == hl)
            n_bl = sum(1 for h in miss[c] if h == hl)
            if 0 < n_bl < n_hl:
                partial.append({"cell": c, "helper_lipid": hl,
                                "rows": n_hl, "blank": n_bl})
        cov[c] = {"measured": len(tr) - len(miss[c]),
                  "blank": len(miss[c]),
                  "helper_lipids_entirely_absent": hl_missing}
    Q["measurement_coverage"] = {
        "per_cell": cov,
        "total_cells_possible": len(tr) * len(RLU),
        "total_measured": sum(v["measured"] for v in cov.values()),
        "total_blank": sum(v["blank"] for v in cov.values()),
        "partially_measured_helper_lipid_in_a_cell_line": partial,
        "blank_pattern_is_all_or_nothing": len(partial) == 0,
        "note": ("The blanks are not scattered: RLU_N2a and RLU_ARPE19 are "
                 "blank on exactly the 14PA and DDAB rows. The authors' own "
                 "utilities.py::select_input_params drops the 14PA and DDAB "
                 "one-hot columns for cells 'ARPE19' and 'N2a', which is the "
                 "same fact stated in code: those two helper lipids were "
                 "never run in those two cell lines."),
    }

    # ========================================================= CHECK 8
    # What scale are the stored values on? Established by MEASUREMENT, because
    # the README states only 'natural log' and the paper's Data Preprocessing
    # describes a floor and a per-cell min-max that the file does NOT carry.
    sc = {}
    for c in RLU:
        v = [float(r[c]) for r in tr if r[c].strip() != ""]
        sc[c] = {"n": len(v), "min": min(v), "max": max(v),
                 "n_below_RLU_floor_1.5": sum(1 for x in v if x < 1.5),
                 "n_outside_unit_interval": sum(1 for x in v if x > 1.0)}
    # the ln fingerprint: exp(value) should land on half-integers for the
    # smallest values, because the stored number is ln(mean of replicates)
    small = sorted(set(round(float(r["RLU_HepG2"]), 3) for r in tr
                       if r["RLU_HepG2"].strip() != ""))[:8]
    Q["scale_is_ln_not_minmax"] = {
        "per_cell": sc,
        "cells_whose_values_all_fit_[0,1]":
            [c for c in RLU if sc[c]["n_outside_unit_interval"] == 0],
        "min_over_all_cells": min(v["min"] for v in sc.values()),
        "max_over_all_cells": max(v["max"] for v in sc.values()),
        "floor_applied_in_file": False,
        "minmax_applied_in_file": False,
        "smallest_RLU_HepG2_values": small,
        "exp_of_those": [round(math.exp(x), 4) for x in small],
        "evidence": ("Values run from %.3f to %.3f, so they are neither "
                     "floored at the paper's RLU_floor = 1.5 nor rescaled to "
                     "[0,1]. exp() of the smallest distinct values lands on "
                     "half-integers (1, 1.5, 2, 2.5, 3 ...), the fingerprint "
                     "of ln(mean of replicate raw RLU). The authors' "
                     "utilities.py::extract_training_data applies BOTH the "
                     "floor and MinMaxScaler at load time, in code, to this "
                     "same file -- so the file is upstream of both. Read from "
                     "the repo, never executed."
                     % (min(v["min"] for v in sc.values()),
                        max(v["max"] for v in sc.values()))),
    }

    # ========================================================= CHECK 9
    # The authors' own QC filters, applied to this file, so the gap between
    # 'rows in the file' and 'rows the paper modelled' is a number, not a guess.
    n_size0 = sum(1 for r in tr if float(r["Size"]) == 0)
    n_big = sum(1 for r in tr if float(r["Size"]) > 10000)
    n_pdi = sum(1 for r in tr if float(r["PDI"]) > 1)
    n_zeta0 = sum(1 for r in tr if float(r["Zeta"]) == 0)
    keep = sum(1 for r in tr
               if float(r["Size"]) != 0 and float(r["Size"]) <= 10000
               and float(r["PDI"]) <= 1 and float(r["Zeta"]) != 0)
    Q["authors_own_filters"] = {
        "source": ("utilities.py::extract_training_data, with the cutoffs in "
                   "the run directory name Final_HL_Features_PDI1_RLU1.5_"
                   "SIZE10000"),
        "Size_eq_0": n_size0, "Size_gt_10000": n_big,
        "PDI_gt_1": n_pdi, "Zeta_eq_0": n_zeta0,
        "rows_surviving_all_four": keep,
        "note": ("These filters are the AUTHORS' modelling choices, not data "
                 "quality verdicts on the transfection measurement. They are "
                 "reported, not applied: a row whose zeta reads exactly 0 "
                 "still has a real luciferase reading."),
    }

    # ========================================================= CHECK 10
    # The validation workbook: 72 rows, and its levels disjoint from training.
    vlev = {}
    for col in ["HL_(IL+HL)", "(IL+HL)", "PEG_(Chol+PEG)", "NP_ratio"]:
        tset = sorted(set(round(float(r[col]), 4) for r in tr))
        vset = sorted(set(round(float(r[col]), 4) for r in va))
        vlev[col] = {"training": tset, "validation": vset,
                     "intersection": sorted(set(tset) & set(vset))}
    vil = sorted(set(str(r["Ionizable_Lipid"]).strip() for r in va))
    vhl = sorted(set(str(r["Helper_lipid"]).strip() for r in va))
    # descriptor agreement across the two files
    dmis = []
    for hl in vhl:
        vr = [r for r in va if str(r["Helper_lipid"]).strip() == hl][0]
        tr_r = [r for r in tr if r["Helper_lipid"].strip() == hl][0]
        for c in DESC:
            if abs(float(vr[c]) - float(tr_r[c])) > 1e-9:
                dmis.append((hl, c, float(vr[c]), float(tr_r[c])))
    vcov = {}
    for c in RLU:
        v = [f(r.get(c)) for r in va]
        v = [x for x in v if x is not None]
        vcov[c] = {"measured": len(v), "min": min(v), "max": max(v)}
    Q["validation_workbook"] = {
        "rows": len(va),
        "Ionizable_Lipid_values": vil,
        "Helper_lipid_values": vhl,
        "levels_vs_training": vlev,
        "all_four_levels_disjoint_from_training":
            all(len(v["intersection"]) == 0 for v in vlev.values()),
        "helper_descriptor_mismatches_vs_training": dmis,
        "coverage": vcov,
        "note": ("The de-novo validation library. Same ONE ionizable lipid, "
                 "same six helper lipids, composition levels chosen to be "
                 "absent from the training library -- which the README states "
                 "and which is confirmed here as an empty intersection on all "
                 "four compositional factors. Values are full-precision "
                 "floats, not the training file's 3 decimals."),
    }

    # ---- shared formulations between the two files (trap 4: measure the ratio)
    def fkey(r):
        return (str(r["Helper_lipid"]).strip(),
                round(float(r["HL_(IL+HL)"]), 4),
                round(float(r["(IL+HL)"]), 4),
                round(float(r["PEG_(Chol+PEG)"]), 4),
                round(float(r["NP_ratio"]), 4))
    tkeys = set(fkey(r) for r in tr)
    vkeys = set(fkey(r) for r in va)
    Q["shared_formulations_between_files"] = {
        "n": len(tkeys & vkeys),
        "ratio_measurable": len(tkeys & vkeys) > 0,
        "note": ("No formulation appears in both files, by construction, so "
                 "there is no shared-formulation ratio to measure between "
                 "them. The two files are therefore never averaged, rescaled "
                 "or tied to each other; they are emitted as what they are."),
    }

    # ========================================================= CHECK 11
    # The tautology, labelled as a tautology.
    Q["sum_to_100_is_tautological"] = {
        "credited_as_a_check": False,
        "why": ("HL+IL+PEG+CHL = a + (100-a) = 100 identically for any "
                "(a,h,p). It cannot fail, so it is not evidence. Recorded so "
                "its absence is not read as an oversight."),
    }

    # ---------------------------------------------------- derived mol%
    rows = []
    for r in tr:
        a = float(r["(IL+HL)"])
        h = float(r["HL_(IL+HL)"])
        p = float(r["PEG_(Chol+PEG)"])
        rows.append({
            "Formula label": r["Formula label"],
            "Helper_lipid_from_onehot": [c for c in OHE_COLS
                                         if r[c].strip() == "1"][0],
            "IL_molpct": a * (100.0 - h) / 100.0,
            "HL_molpct": a * h / 100.0,
            "CHL_molpct": (100.0 - a) * (100.0 - p) / 100.0,
            "PEG_molpct": (100.0 - a) * p / 100.0,
        })
    Q["derived_molpct_examples"] = rows[:3]

    with open(out_path, "w") as fh:
        json.dump(Q, fh, indent=1)

    # ---------------------------------------------------- the gate
    fails = []
    if not Q["onehot_inversion"]["lossless"]:
        fails.append("one-hot inversion is not lossless")
    if not Q["factorial"]["complete_and_unreplicated"]:
        fails.append("factorial is not complete/unreplicated")
    if not Q["factorial"]["product_equals_rows"]:
        fails.append("product of level counts != row count")
    if not Q["si_table1_all_agree"]:
        fails.append("file levels disagree with SI Table 1")
    for k, v in Q["error_class_40_two_parameterisations"].items():
        if v["rows_violating_at_1e-6"]:
            fails.append("class-40: %s disagrees with %s on %d rows"
                         % (k, v["derived_from"], v["rows_violating_at_1e-6"]))
        if not v["ratio_levels_are_exact_integers"]:
            fails.append("class-40: %s ratio levels are not integers" % k)
    if not Q["helper_descriptor_block"]["constant_within_helper_lipid"]:
        fails.append("a helper-lipid descriptor varies within its own lipid")
    if not Q["helper_descriptor_block"]["all_agree_with_si_table3"]:
        fails.append("helper descriptors disagree with SI Table 3")
    if Q["what_a_row_is"]["distinct_Ionizable_Lipid_values"] != 1:
        fails.append("Ionizable_Lipid is no longer single-valued")
    if not Q["validation_workbook"]["all_four_levels_disjoint_from_training"]:
        fails.append("validation levels overlap training levels")
    if Q["validation_workbook"]["helper_descriptor_mismatches_vs_training"]:
        fails.append("the two files disagree on a helper-lipid descriptor")
    # --- both added because a seeded mutation moved the metric while the gate
    # --- stayed silent; see 01_seed_results.json and 07_qc.json.
    if not Q["measurement_coverage"]["blank_pattern_is_all_or_nothing"]:
        fails.append(
            "a helper lipid is PARTIALLY measured in a cell line (%s) -- "
            "missingness here is supposed to be whole-arm, so a scattered "
            "blank means a value was lost, not never run"
            % Q["measurement_coverage"][
                "partially_measured_helper_lipid_in_a_cell_line"][:3])
    flat = Q["scale_is_ln_not_minmax"]["cells_whose_values_all_fit_[0,1]"]
    if flat:
        fails.append(
            "target column(s) %s have no value above 1.0, which is what a "
            "per-cell min-max rescale looks like. The corpus needs the ln "
            "scale, upstream of the paper's MinMaxScaler; a rescaled column "
            "has had its between-formulation spread redefined." % flat)
    if Q["shared_formulations_between_files"]["n"]:
        fails.append(
            "a formulation now appears in BOTH files -- measure the ratio "
            "before emitting either; never average or rescale two sources")

    print(json.dumps({"checks_written_to": out_path, "FAILURES": fails},
                     indent=1))
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1], sys.argv[2], sys.argv[3]))
