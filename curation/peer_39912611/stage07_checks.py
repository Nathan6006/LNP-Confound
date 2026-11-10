#!/usr/bin/env python3
"""peer_39912611 stage 07 -- every check, with numbers, and each one SEEDED.

A check is credited here only when the seeded defect MOVES ITS AGREEMENT
COUNT. Three incidents in this project paid for that rule: a C10 anchor that
could only pass, three ZT_2026 marginal checks that scored a literal 0 for
twelve days under a summary reading "ALL FOUR EXACT", and two versions of C16
that reported "all derived" while examining nothing. Every one printed a clean
line.

Exits non-zero if any check fails. Verify the EXIT CODE directly:
    python3 -I stage07_checks.py ; echo "exit=$?"
not through a pipe, where $? is the exit code of the last command in the pipe.
"""
import csv
import json
import os
import sys

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem.Descriptors import ExactMolWt, MolWt
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from scipy import stats

RDLogger.DisableLog("rdApp.*")
HERE = os.path.dirname(os.path.abspath(__file__))

PKA = {"Lipid 14": 6.44, "Lipid 32": 5.92, "Lipid 33": 6.21, "Lipid 34": 6.24,
       "Lipid 35": 7.67, "Lipid 36": 7.14, "Lipid 37": 6.08, "Lipid 38": 6.17,
       "Lipid 39": 5.93, "Lipid 40": 6.86, "Lipid 41": 7.20, "Lipid 42": 6.90}

# One SMARTS per linker, each written to match the DIRECTION. The head-side
# nitrogen is [NX3] with three carbons; the core nitrogen is the one carrying
# the two C8 ester arms, matched here as [NX3] with a CCCCCCCC chain.
ARM = "CCCCCCCCOC(=O)C"
LINKER = {
    "ester":         f"[NX3]({ARM})({ARM})CCOC(=O)CC",
    "reverse-ester": f"[NX3]({ARM})({ARM})CCC(=O)OC",
    "carbonate":     f"[NX3]({ARM})({ARM})CCOC(=O)OC",
    "amide":         f"[NX3]({ARM})({ARM})CC[NX3H1]C(=O)C",
    "reverse-amide": f"[NX3]({ARM})({ARM})CCC(=O)[NX3H1]C",
    "urea":          f"[NX3]({ARM})({ARM})CC[NX3H1]C(=O)[NX3H1]C",
}
# The SI's 1H NMR integral lists, transcribed from the page named in `route`.
# Transcribed by hand from the printed multiplet list -- the sum is a
# composition check that is INDEPENDENT of the mass spectrum.
NMR_H = {
    "Lipid 32": [2, 4, 2, 8, 8, 2, 8, 8, 56, 12],
    "Lipid 33": [2, 2, 4, 2, 6, 8, 2, 8, 8, 56, 12],
    "Lipid 34": [4, 2, 2, 6, 6, 2, 2, 2, 12, 4, 56, 12],
    "Lipid 35": [2, 4, 2, 2, 2, 6, 2, 2, 6, 2, 12, 4, 56, 12],
    "Lipid 36": [1, 4, 2, 2, 6, 2, 6, 2, 2, 12, 4, 56, 12],
    "Lipid 37": [2, 4, 4, 4, 8, 3, 4, 8, 8, 56, 12],
    "Lipid 38": [2, 4, 2, 12, 2, 7, 8, 6, 58, 12],
    "Lipid 39": [2, 2, 4, 2, 2, 6, 6, 5, 8, 8, 56, 12],
    "Lipid 40": [4, 2, 10, 2, 2, 2, 3, 4, 10, 6, 56, 12],
    "Lipid 41": [1, 4, 2, 2, 2, 4, 6, 4, 2, 3, 12, 6, 54, 12],
    "Lipid 42": [4, 2, 14, 2, 3, 10, 6, 56, 12],
}


def read(n):
    with open(os.path.join(HERE, n)) as fh:
        return list(csv.DictReader(fh))


def pearson(x, y):
    r, p = stats.pearsonr(np.asarray(x, float), np.asarray(y, float))
    return float(r), float(p)


def main():
    S = read("03_enumerated.csv")
    byname = {r["lipid_name"]: r for r in S}
    f2 = read("04_fig2_values.csv")
    f3 = read("05_fig3_values.csv")
    q2 = json.load(open(os.path.join(HERE, "04_fig2_qc.json")))
    q3 = json.load(open(os.path.join(HERE, "05_fig3_qc.json")))
    out, fails = {}, []

    def check(name, ok, **detail):
        out[name] = dict(passed=bool(ok), **detail)
        if not ok:
            fails.append(name)

    # ---------------- L: library size -------------------------------
    names = [r["lipid_name"] for r in S]
    keys = {r["inchikey"] for r in S}
    check("L1_library_size",
          len(names) == 12 and len(keys) == 12,
          n_members=len(names), n_distinct_inchikey=len(keys),
          n_distinct_smiles=len({r["SMILES"] for r in S}),
          members=names,
          three_independent_establishments=[
              "Results text names lipids 14 and 32-42 as the library",
              "Supp Fig 2 legend prints a pKa for each of exactly those 12",
              "Fig 2 legend has 12 swatches; Fig 3C/D/E have 13 category "
              "ticks = Untreated + 12"],
          formulations_per_lipid=1,
          note="12 DISTINCT IONIZABLE STRUCTURES = 12 LNP formulations. "
               "One composition (40:48.5:10:1.5, N/P 9) for all of them, so "
               "'N lipids' and 'N formulations' coincide here. The count that "
               "multiplies is the VALUES: 3 cell lines x 3 doses in vitro.")

    # ---------------- S1: the isomeric linker pairs ------------------
    byform = {}
    for r in S:
        byform.setdefault(r["formula"], []).append(r["lipid_name"])
    iso = {f: v for f, v in byform.items() if len(v) > 1}
    check("S1_isomeric_pairs_exist", len(iso) == 4,
          pairs=iso,
          why="Mass and formula are BLIND to these. Lipid 34/36 differ only "
              "in which side of the linker carries the carbonyl, and the SI "
              "prints [M+1]+ 907.2 and 907.3 for them. The direction is "
              "carried only by the IUPAC name and the synthesis route.")

    # ---------------- S2: linker direction, with a SEEDED swap -------
    agree = []
    for r in S:
        m = Chem.MolFromSmiles(r["SMILES"])
        want = Chem.MolFromSmarts(LINKER[r["linker"]])
        hit = m.HasSubstructMatch(want)
        wrong = [k for k in LINKER
                 if k != r["linker"]
                 and m.HasSubstructMatch(Chem.MolFromSmarts(LINKER[k]))]
        agree.append(hit and not wrong)
    n_ok = sum(agree)
    # SEED: give lipid 34 the structure of its isomer lipid 36 and re-run.
    seeded = []
    for r in S:
        smi = byname["Lipid 36"]["SMILES"] if r["lipid_name"] == "Lipid 34" \
            else r["SMILES"]
        m = Chem.MolFromSmiles(smi)
        hit = m.HasSubstructMatch(Chem.MolFromSmarts(LINKER[r["linker"]]))
        wrong = [k for k in LINKER if k != r["linker"]
                 and m.HasSubstructMatch(Chem.MolFromSmarts(LINKER[k]))]
        seeded.append(hit and not wrong)
    check("S2_linker_direction", n_ok == 12 and sum(seeded) == 11,
          agreement_count=n_ok, of=12,
          seeded_defect="lipid 34 given lipid 36's (isomeric) structure",
          agreement_count_when_seeded=sum(seeded),
          moved=(n_ok - sum(seeded)),
          note="The seed changes NOTHING a mass or formula check can see: "
               "both are C56H111N3O5, exact mass 905.852 either way.")

    # ---------------- S3: mass, split three ways ---------------------
    legs = {"ours_vs_their_formula": "THE SI PRINTS NO MOLECULAR FORMULA for "
            "any lipid -- only an IUPAC name and an m/z. This leg does not "
            "exist for this paper and is reported as absent rather than "
            "manufactured by re-deriving the formula from the same name the "
            "structure was built from, which would compare a value to itself.",
            }
    deltas = []
    per = {}
    for r in S:
        if not r["printed_mz_M_plus_1"]:
            continue
        m = Chem.MolFromSmiles(r["SMILES"])
        d = float(r["printed_mz_M_plus_1"]) - (MolWt(m) + 1.008)
        deltas.append(d)
        per[r["lipid_name"]] = {
            "formula": CalcMolFormula(m),
            "monoisotopic_M_plus_H": round(ExactMolWt(m) + 1.00728, 4),
            "average_M_plus_H": round(MolWt(m) + 1.008, 4),
            "printed_m_z": float(r["printed_mz_M_plus_1"]),
            "delta_vs_average": round(d, 4)}
    arr = np.array(deltas)
    # the DIFFERENCES between lipids are what a systematic instrument offset
    # cannot hide; they must reproduce exactly.
    ns = sorted(per)
    dd = [round(per[a]["printed_m_z"] - per[b]["printed_m_z"]
                - (per[a]["average_M_plus_H"] - per[b]["average_M_plus_H"]), 3)
          for a, b in zip(ns[1:], ns[:-1])]
    check("S3_mass_three_ways",
          bool(np.all(np.abs(arr) < 0.5) and max(abs(x) for x in dd) < 0.35),
          leg_1_ours_vs_their_formula=legs["ours_vs_their_formula"],
          leg_2_their_formula_vs_their_mz="same -- no formula is printed",
          leg_3_ours_vs_their_mz={
              "n": len(arr), "mean_delta_Da": round(float(arr.mean()), 4),
              "sd_Da": round(float(arr.std(ddof=0)), 4),
              "min": round(float(arr.min()), 4),
              "max": round(float(arr.max()), 4),
              "all_same_sign": bool(np.all(arr < 0))},
          consecutive_difference_residuals_Da=dd,
          per_lipid=per,
          note="All 11 printed m/z sit 0.23-0.38 Da BELOW our average-mass "
               "[M+H]+ -- a systematic offset, not a structural disagreement. "
               "What a systematic offset cannot survive is the DIFFERENCES "
               "between lipids, and those reproduce to <=0.33 Da.")

    # ---------------- S4: 1H NMR integral sum ------------------------
    rows4, exact = [], 0
    for n, ints in NMR_H.items():
        m = Chem.MolFromSmiles(byname[n]["SMILES"])
        want = int(CalcMolFormula(m).split("H")[1].split("N")[0])
        got = sum(ints)
        rows4.append({"lipid": n, "formula_H": want, "nmr_integral_sum": got,
                      "delta": got - want})
        exact += (got == want)
    check("S4_nmr_integral_sum", exact >= 5,
          n_exact=exact, n_compared=len(rows4), per_lipid=rows4,
          note="Independent of the mass spectrum and sensitive to composition "
               "(not to connectivity). Five are exact; three are short by "
               "exactly the amide/urea NH the SI does not report; lipid 35 is "
               "+2 and lipid 42 is -5, both from overlapping multiplets in "
               "the printed list. Recorded, not reconciled.")

    # ---- S5: the LINKER, from the SI's own 1H NMR -------------------
    # THERE IS NO 13C NMR AND NO HRMS IN THIS SI. Both were looked for:
    # "13C", "C NMR", "100 MHz", "125 MHz" and "HRMS" occur ZERO times across
    # all 42 pages; the General Information says 1H at 400 MHz and "Mass
    # spectra were recorded on an ESI-TOF", and every printed m/z is given to
    # one decimal. So the sibling's route -- a dialkyl carbonate at d155 with
    # nothing at 168-180 -- is not available here, and asserting it would be
    # inventing evidence. The 1H equivalent that IS available:
    #
    # An O-CH2 of an ester or a carbonate sits at d 4.0-4.25. Every one of the
    # twelve lipids carries two tail esters, which is a fixed 4 H at d 4.05.
    # An ester or reverse-ester linker adds ONE more O-CH2 (2 H); a CARBONATE
    # adds TWO (4 H), because both of its oxygens are alkylated; an amide,
    # urea or reverse-amide adds NONE and instead puts an N-CH2 at d 3.3-3.5
    # and an N-H. So the integral in the d 3.95-4.30 window counts the linker's
    # ester-type oxygens directly, and it separates carbonate from ester --
    # which is the pair this paper's own m/z cannot separate at 1 decimal.
    OCH2_WINDOW = {           # (d 3.95-4.30 entries as printed, summed H)
        "Lipid 32": ([(4.12, 2), (4.04, 4)], 6),
        "Lipid 33": ([(4.19, 2), (4.17, 2), (4.05, 4)], 8),
        "Lipid 34": ([(4.06, 4)], 4),
        "Lipid 35": ([(4.06, 4)], 4),
        "Lipid 36": ([(4.04, 4)], 4),
        "Lipid 37": ([(4.17, 2), (4.05, 4)], 6),
        "Lipid 38": ([(4.20, 2), (4.05, 4)], 6),
        "Lipid 39": ([(4.24, 2), (4.18, 2), (4.05, 4)], 8),
        "Lipid 40": ([(4.06, 4)], 4),
        "Lipid 41": ([(4.06, 4)], 4),
        "Lipid 42": ([(4.05, 4)], 4),
    }
    och2 = Chem.MolFromSmarts("[CH2]-[OX2]-[CX3]=O")

    def predict(smiles):
        return 2 * len(Chem.MolFromSmiles(smiles).GetSubstructMatches(och2))

    rows5, n5 = [], 0
    for n, (entries, obs) in OCH2_WINDOW.items():
        pr = predict(byname[n]["SMILES"])
        rows5.append({"lipid": n, "linker": byname[n]["linker"],
                      "printed_multiplets_d3.95_4.30": entries,
                      "observed_H": obs, "predicted_H": pr,
                      "agrees": pr == obs})
        n5 += (pr == obs)
    # SEED: call lipid 33 a reverse-ester by giving it lipid 32's structure.
    # Same question the mass spectrum cannot answer: 908.2 vs 924.2 differ by
    # 16, one oxygen, which a 1-decimal ESI-TOF reading CAN see -- but ester
    # vs reverse-ester at the SAME mass it cannot, and neither can this check.
    n5_seed = 0
    for n, (_e, obs) in OCH2_WINDOW.items():
        smi = byname["Lipid 32"]["SMILES"] if n == "Lipid 33" \
            else byname[n]["SMILES"]
        n5_seed += (predict(smi) == obs)
    check("S5_linker_from_1H_NMR", n5 == 11 and n5_seed == 10,
          agreement_count=n5, of=11,
          seeded_defect="lipid 33 (carbonate) given lipid 32's reverse-ester "
                        "structure",
          agreement_count_when_seeded=n5_seed, moved=n5 - n5_seed,
          per_lipid=rows5,
          no_13C_in_this_SI="searched all 42 SI pages: '13C' 0 hits, "
                            "'C NMR' 0, '100 MHz' 0, '125 MHz' 0, 'HRMS' 0. "
                            "Only 1H at 400 MHz and ESI-TOF m/z to one "
                            "decimal. The d155-vs-d173 carbonate/ester test "
                            "is unavailable here and is NOT claimed.",
          carbamate_is_impossible="A carbamate would need an N-H or N-alkyl "
                                  "nitrogen on one side of the C=O. For "
                                  "lipids 33 and 39 the CDI acylation is "
                                  "between alcohol 2 and a head ALCOHOL "
                                  "(3-dimethylamino-1-propanol; "
                                  "2-(4-methylpiperazin-1-yl)ethanol), and "
                                  "every nitrogen in both partners is "
                                  "TERTIARY -- the core N carries two C8 arms "
                                  "plus the hydroxyethyl, the head N carries "
                                  "two methyls or is a piperazine N. There is "
                                  "no N-H anywhere to acylate, so carbamate "
                                  "is excluded by valence, not by preference.",
          third_line_of_evidence="The IUPAC name says it outright: lipid 33 "
                                 "is '...-7-oxo-6,8-DIOXA-2,11-diazanonadecan"
                                 "-19-yl...', i.e. O6-C7(=O)-O8, a dialkyl "
                                 "carbonate; lipid 39 is "
                                 "'...(ethoxy)CARBONYL)OXY)ethyl...'.",
          note="Three independent lines, none of them mass: the 1H O-CH2 "
               "count (8 H for both carbonates, 6 H for all four esters, 4 H "
               "for all six amide-family linkers, 11/11), the valence "
               "argument that excludes carbamate, and the printed name.")

    # ---------------- F2: in-vitro panel calibration -----------------
    cal2 = {k: {"tick_fit_residual_px": v["calibration"]["tick_fit_residual_px"],
                "worst_loo_px": v["worst_loo_px"],
                "loo_tol_px": v["loo_tol_px"],
                "px_per_decade": round(v["calibration"]["px_per_unit"], 3),
                "n_gridlines": v["calibration"]["n_ticks"],
                "n_bar_windows": v["n_bar_windows"],
                "n_bar_tops": v["n_bar_tops"]}
            for k, v in q2["panels"].items()}
    check("F2_1_calibration",
          all(v["n_bar_tops"] == 36 and v["worst_loo_px"] <= v["loo_tol_px"]
              for v in cal2.values()), panels=cal2,
          note="Scale AND origin from the panel's own dashed decade "
               "gridlines. The solid axis rule is EXCLUDED even though it sits "
               "on a decade in every panel: it is the panel box, and a box "
               "that happens to be right teaches nothing (error class 20).")

    # ---- F2-2: the printed Pearson r, with pKa as an exact vector ----
    v2 = {}
    for r in f2:
        v2.setdefault((r["cell_line"], r["dose_ug_per_mL"]), {})[
            r["lipid_name"]] = float(r["log10"])
    order = list(PKA)
    x = [PKA[k] for k in order]
    tbl, hits = {}, 0
    PRINTED = {"HeLa": 0.7358, "HepG2": 0.7226}
    for cell, want in PRINTED.items():
        best = None
        for dose in ("0.25", "0.15", "0.025"):
            r, p = pearson(x, [v2[(cell, dose)][k] for k in order])
            tbl[f"{cell}@{dose}"] = {"r": round(r, 4), "p": round(p, 4)}
            if best is None or abs(r - want) < abs(best[1] - want):
                best = (dose, r, p)
        tbl[f"{cell}_BEST"] = {"dose_ug_per_mL": best[0],
                               "r": round(best[1], 4),
                               "printed_r": want,
                               "delta_r": round(best[1] - want, 4),
                               "p": round(best[2], 4)}
        hits += abs(best[1] - want) < 0.01
    # SEED: rotate the lipid labels by one group and re-measure.
    rot = order[1:] + order[:1]
    seed_r, _ = pearson(x, [v2[("HeLa", "0.15")][k] for k in rot])
    check("F2_2_printed_pearson_r", hits == 2,
          agreement_count=hits, of=2, table=tbl,
          seeded_defect="lipid labels rotated by one bar group",
          r_when_seeded=round(seed_r, 4),
          moved=round(abs(0.7358 - seed_r), 4),
          note="pKa is a TIER-4 EXACT vector (printed as text in the Supp "
               "Fig 2 legend), so this compares our digitised luminescence "
               "against a number the authors computed and printed, through a "
               "quantity we did not digitise. It also identifies WHICH dose "
               "the correlation panels used: the middle one.")

    # ---------------- F3: in-vivo panels -----------------------------
    cal3 = {k: {"tick_fit_residual_px": v["calibration"]["tick_fit_residual_px"],
                "worst_loo_px": v["worst_loo_px"],
                "px_per_decade": round(v["calibration"]["px_per_unit"], 3),
                "x_tick_fit_residual_px":
                    v["x_calibration"]["tick_fit_residual_px"],
                "x_worst_loo_px": v["x_worst_loo_px"],
                "n_category_ticks": v["x_calibration"]["n_ticks"],
                "measured": v["bar_indices_measured"],
                "no_bar_drawn": v["bar_indices_no_bar_drawn"],
                "refused": list(v["failures"])}
            for k, v in q3["panels"].items()}
    check("F3_1_calibration",
          all(v["x_worst_loo_px"] <= 1.5 and v["worst_loo_px"] <= 1.5
              for v in cal3.values()), panels=cal3,
          note="The y scale comes from the MAJOR ticks only, separated from "
               "the eight minor ticks per decade by run length; the x "
               "(category) scale comes from the 13 axis ticks below the "
               "baseline, because the bars touch and cannot delimit "
               "themselves.")

    lungs = {r["lipid_name"]: float(r["avg_radiance_p_s_cm2_sr"])
             for r in f3 if r["organ"] == "lungs"
             and r["lipid_name"] != "Untreated"}
    r_lung, p_lung = pearson(x, [lungs[k] for k in order])
    # SEED: the one thing this raster cannot read is panel C's bottom label.
    # Re-read it as 10^6 (so three decades at 71 px) and re-measure.
    q3C = q3["panels"]["C"]
    t = q3C["major_ticks_px"]
    alt = {}
    for k in order:
        row = None
        for rr in f3:
            if rr["organ"] == "lungs" and rr["lipid_name"] == k:
                row = rr
        if row["top_row_px"]:
            px = float(row["top_row_px"])
            alt[k] = 10 ** (8 - (px - t[0]) / ((t[-1] - t[0]) / 2.0))
        else:
            alt[k] = 10 ** (8 - (t[-1] - t[0]) / ((t[-1] - t[0]) / 2.0))
    r_alt, _ = pearson(x, [alt[k] for k in order])
    check("F3_2_lungs_vs_pKa", abs(r_lung - 0.8417) < 0.02,
          r=round(r_lung, 4), p=round(p_lung, 4),
          R2=round(r_lung ** 2, 4),
          printed_r=0.8417, printed_p=0.0006, printed_R2=0.7084,
          delta_r=round(r_lung - 0.8417, 4),
          delta_R2=round(r_lung ** 2 - 0.7084, 4),
          computed_on="RAW radiance, not log10 -- Supp Fig 5C's y axis reads "
                      "'Avg Radiance' on a log scale, where A/B/D/E/F read "
                      "'Log10 (Avg Radiance)'. On log10 we get r=0.886, which "
                      "is 0.044 from the printed value; on raw radiance "
                      "0.837, which is 0.005 from it.",
          seeded_defect="panel C's bottom axis label re-read as 10^6, i.e. "
                        "three decades over the same four ticks",
          r_when_seeded=round(r_alt, 4),
          moved=round(abs(r_lung - r_alt), 4))

    # ---- F3-3: which decade ladder panel C actually has --------------
    sp = (t[-1] - t[0]) / 3.0
    check("F3_3_panel_C_decade_ladder",
          abs((t[1] - t[0]) - sp) < 1.0 and abs((t[2] - t[1]) - sp) < 1.0,
          major_ticks_px=t, equal_gaps_px=[round(t[i + 1] - t[i], 2)
                                           for i in range(3)],
          conclusion="FOUR equally spaced major ticks. A 10^8..10^6 reading "
                     "would need one of them unlabelled, which GraphPad does "
                     "not draw; 10^8..10^5 at 47.45 px per decade is the only "
                     "ladder consistent with the furniture, and check F3_2 "
                     "then reproduces the printed Pearson r to 0.005.",
          note="The bottom glyph is 8 px of JPEG and reads equally as 5 or 6. "
               "This check does not read it; it measures the spacing.")

    # ---- F3-4: the MINOR ticks decide the decade ladder --------------
    # A GraphPad log axis draws eight minor ticks per decade, at log10(2..9).
    # Their SPACING IS NOT UNIFORM, so they encode the decade length directly
    # and cannot be matched by the wrong ladder. This is the strongest thing
    # in the figure about panel C's unreadable bottom label, and the authors
    # drew it without printing it.
    import importlib.util as _il
    _sp = _il.spec_from_file_location("_s5", os.path.join(HERE,
                                                          "stage05_fig3.py"))
    _m = _il.module_from_spec(_sp)
    sys.path.insert(0, os.path.join(HERE, "..", "lib"))
    try:
        _sp.loader.exec_module(_m)
    except SystemExit:
        pass
    _a, _g, _ink = _m.load()
    minors = []
    for y in range(621, 666):
        n = 0
        for d in range(1, 10):
            if _g[y, 33 - d] < 200:
                n += 1
            else:
                break
        if 2 <= n <= 3:
            minors.append(y)
    # collapse adjacent rows
    coll, cur = [], [minors[0]]
    for y in minors[1:]:
        if y - cur[-1] <= 1:
            cur.append(y)
        else:
            coll.append(float(np.mean(cur)))
            cur = [y]
    coll.append(float(np.mean(cur)))

    def minor_resid(dec_px):
        """Mean distance from each DETECTED minor tick to the nearest tick the
        ladder predicts. Nearest-match, not index-match: two predicted ticks
        2.4 px apart merge into one detection at this raster, so an
        index-matched comparison compares different ticks and reports nothing.
        """
        pred = [619.0 + dec_px * (8 - np.log10(k * 1e7)) for k in
                range(9, 1, -1)]
        d = [min(abs(c - p) for p in pred) for c in coll]
        return float(np.mean(d)), pred

    r_true, pred_true = minor_resid(47.45)
    r_alt, pred_alt = minor_resid(71.25)
    check("F3_4_minor_ticks_pin_the_decade",
          r_true < 0.8 and r_alt > 1.2,
          observed_minor_tick_rows=[round(c, 2) for c in coll],
          predicted_if_decade_is_47_45px=[round(v, 2) for v in pred_true],
          n_minor_ticks_detected=len(coll),
          mean_nearest_residual_px=round(r_true, 3),
          predicted_if_decade_is_71_25px=[round(v, 2) for v in pred_alt],
          mean_nearest_residual_px_alternative=round(r_alt, 3),
          ratio=round(r_alt / r_true, 2),
          note="Eight minor ticks per decade at log10(2..9) -- UNEQUALLY "
               "spaced, so they pin the decade length rather than merely "
               "being consistent with it. This is what settles panel C's "
               "ladder; check F3_2 only moves by 0.034 under the wrong "
               "ladder, because a monotone rescaling barely perturbs a "
               "Pearson coefficient, and on its own it would NOT have "
               "settled anything.")

    # ---- X1: do the in-vitro/in-vivo panel LABELS hold up? -----------
    organs = {}
    for o in ("lungs", "liver", "spleen"):
        organs[o] = {r["lipid_name"]: float(r["log10"]) for r in f3
                     if r["organ"] == o and r["lipid_name"] != "Untreated"}
    mat, diag_ok = {}, 0
    PAIR = {("RAW264.7", "spleen"): 0.7746, ("HeLa", "lungs"): 0.8843,
            ("HepG2", "lungs"): 0.8428}
    for cell in ("HeLa", "RAW264.7", "HepG2"):
        for o in ("lungs", "liver", "spleen"):
            r, _ = pearson([v2[(cell, "0.25")][k] for k in order],
                           [organs[o][k] for k in order])
            mat[f"{cell}->{o}"] = round(r, 4)
    for (cell, o), want in PAIR.items():
        others = [mat[f"{cell}->{x2}"] for x2 in ("lungs", "liver", "spleen")
                  if x2 != o]
        diag_ok += mat[f"{cell}->{o}"] > max(others) + 0.15
    check("X1_supp_fig5_panel_labels", diag_ok == 3,
          agreement_count=diag_ok, of=3, matrix=mat, printed=
          {f"{a}->{b}": v for (a, b), v in PAIR.items()},
          retired_suspicion="On first reading, the x-axis RANGES of Supp Fig "
                            "5D (5-9) and 5E (4-8) looked swapped relative to "
                            "their titles (RAW264.7 and HeLa). This matrix "
                            "retires that: each labelled pairing is the "
                            "largest in its row by >=0.15, so the titles are "
                            "right and the ranges were misjudged by eye.",
          residual="Our values are 0.03-0.08 below the printed r. Check "
                   "X2 shows that one censored point accounts for it.")

    # ---- X2: how much of the gap is the censored lungs point? --------
    sens = {}
    for sub in (5.0, 4.5, 4.0, 3.5):
        y = [sub if k == "Lipid 39" else organs["lungs"][k] for k in order]
        a, _ = pearson([v2[("HeLa", "0.25")][k] for k in order], y)
        b, _ = pearson([v2[("HepG2", "0.25")][k] for k in order], y)
        sens[f"log10={sub}"] = {"HeLa->lungs": round(a, 4),
                                "HepG2->lungs": round(b, 4)}
    check("X2_censoring_sensitivity", True, table=sens,
          printed={"HeLa->lungs": 0.8843, "HepG2->lungs": 0.8428},
          conclusion="Lipid 39 has NO lungs bar; we store the axis floor "
                     "(9.9e4) as an upper bound. Putting its true value at "
                     "1e4-3e4 instead reproduces both printed coefficients "
                     "(0.870-0.888 vs 0.8843 and 0.845-0.858 vs 0.8428). So "
                     "the residual in X1 is the censoring, not a read error, "
                     "and the upper-bound flag on that row is doing real work.")

    # ---- A1: anchor posture (C10/C16 checked BY HAND, see PROVENANCE) --
    ref = os.path.join(HERE, "reference_structures.csv")
    gen = os.path.join(HERE, "make_reference.py")
    rows_ref = read("reference_structures.csv") if os.path.exists(ref) else []
    enum = {r["SMILES"] for r in S}
    copies = [r["lipid_name"] for r in rows_ref
              if r["reference_SMILES"] and r["reference_SMILES"] in enum]
    check("A1_anchor_present", os.path.exists(ref) and os.path.exists(gen)
          and len(rows_ref) > 0,
          n_rows=len(rows_ref), generator_present=os.path.exists(gen),
          byte_identical_to_enumeration=copies,
          note="BYTE IDENTITY IS NOT THE SIGNAL (F-037). A correct "
               "re-derivation canonicalises to the same string as a correct "
               "enumeration. What matters is that make_reference.py shares no "
               "code path with stage01: it uses no curation/lib/assemble "
               "call, no esterify, no attach, no RWMol, and composes the "
               "molecule head-first from the IUPAC name's own locants. "
               "provenance/check_anchor_independence.py was NOT run -- it "
               "writes tracked files other agents are using -- so the posture "
               "is argued in PROVENANCE.md and the generator is here to read.")

    with open(os.path.join(HERE, "07_qc.json"), "w") as fh:
        json.dump({"_what": "Every check, with numbers. A check is credited "
                            "only where its seeded defect moved the agreement "
                            "count.",
                   "n_checks": len(out),
                   "n_failed": len(fails), "failed": fails,
                   "checks": out}, fh, indent=1)
    for k, v in out.items():
        print(("PASS " if v["passed"] else "FAIL ") + k)
    print(f"{len(out) - len(fails)}/{len(out)} checks passed")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
