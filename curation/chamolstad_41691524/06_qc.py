#!/usr/bin/env python3
"""Every check with numbers, each SEEDED so the agreement count is shown to move.

A check that has never been made to fail has not been shown to work
(CLAUDE.md section 6). For each check below the script reports the agreement
count on the true inputs and on a deliberately corrupted copy; the pair is the
evidence, not the clean line alone.
"""
import sys, os, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
import axis

PROT = 1.00727646
enu = json.load(open(os.path.join(HERE, "02_enumerate.json")))
geo = json.load(open(os.path.join(HERE, "04_figS4_geometry.json")))
qc = {"checks": {}}

# ------------------------------------------------------------------ P1 mass
# The paper prints three ESI m/z and NO molecular formula, so only one of the
# three splits the brief asks for is available. Saying which is missing is
# part of the result.
PUB = {"4-tailed H9T6 (tR 6.7)": 1059.99,
       "tR 12.4": 1288.09, "tR 12.9": 1362.11}
def ladder(tail_acid):
    tren = Descriptors.ExactMolWt(Chem.MolFromSmiles("NCCN(CCN)CCN"))
    gly = Descriptors.ExactMolWt(Chem.MolFromSmiles("OCC1CO1"))
    tm = Descriptors.ExactMolWt(Chem.MolFromSmiles(tail_acid)) + gly - 18.0105646
    out = {f"{n}-tailed": tren + n * tm + PROT for n in range(2, 8)}
    out["5-tailed + 1 glycidol arm"] = tren + 5 * tm + gly + PROT
    return out

def mass_agreement(tail_acid, tol=0.20):
    L = ladder(tail_acid)
    hits = {}
    for label, obs in PUB.items():
        best = min(L.items(), key=lambda kv: abs(kv[1] - obs))
        if abs(best[1] - obs) <= tol:
            hits[label] = {"assignment": best[0], "computed_MH": round(best[1], 3),
                           "published_MH": obs, "delta_Da": round(obs - best[1], 3),
                           "ppm": round(1e6 * (obs - best[1]) / best[1], 1)}
    return hits

acids = {"C8": "CCCCCCCC(=O)O", "C10": "CCCCCCCCCC(=O)O",
         "C12": "CCCCCCCCCCCC(=O)O", "C14": "CCCCCCCCCCCCCC(=O)O"}
by_len = {k: mass_agreement(v) for k, v in acids.items()}
qc["checks"]["P1_LCMS_three_published_m_z"] = {
    "what": "our computed [M+H]+ for the TREN/glycidyl-ester ladder vs the three "
            "ESI m/z printed on SI p.5",
    "split_ours_vs_their_formula": "NOT AVAILABLE - the paper prints no "
                                   "molecular formula for H9T6 or any congener",
    "split_their_formula_vs_their_m_z": "NOT AVAILABLE, same reason",
    "split_ours_vs_their_m_z": by_len["C10"],
    "agreement_count_C10": len(by_len["C10"]),
    "max_abs_delta_Da": round(max(abs(h["delta_Da"]) for h in by_len["C10"].values()), 3),
    "SEEDED_with_C8": len(by_len["C8"]), "SEEDED_with_C12": len(by_len["C12"]),
    "SEEDED_with_C14": len(by_len["C14"]),
    "verdict": ("3/3 assigned at C10 and 0/3 at every other chain length the SI "
                "lists, so the tail of H9T6 is decanoate; the count moves from "
                "3 to 0 under the seed.")}

# ------------------------------------------------------- P2 NMR integrals
m = Chem.MolFromSmiles(enu["series"]["6"]["smiles"])
mh = Chem.AddHs(m)
def countH(smarts):
    patt = Chem.MolFromSmarts(smarts)
    return sum(mh.GetAtomWithIdx(i[0]).GetTotalNumHs()
               if False else sum(1 for n in mh.GetAtomWithIdx(i[0]).GetNeighbors()
                                 if n.GetSymbol() == "H")
               for i in mh.GetSubstructMatches(patt))
PRED = {
 "CH3 (0.89 t)":            countH("[CH3][CH2]"),
 "Hbeta (1.53-1.69)":       countH("[CH2]([CH2]C(=O)O)[CH2]"),
 "Halpha (2.16-2.37)":      countH("[CH2]C(=O)O"),
 "bulk CH2 (1.21-1.36)":    countH("[CH3][CH2]") // 3 * 0 + (
                             countH("[CH2;$([CH2][CH2][CH2])]")),
}
# bulk CH2 = every chain CH2 that is not alpha, not beta, counted explicitly
alpha = {i[0] for i in mh.GetSubstructMatches(Chem.MolFromSmarts("[CH2]C(=O)O"))}
beta = {i[0] for i in mh.GetSubstructMatches(
    Chem.MolFromSmarts("[CH2]([CH2]C(=O)O)[CH2]"))}
chain = {a.GetIdx() for a in mh.GetAtoms()
         if a.GetSymbol() == "C" and a.GetDegree() - sum(
             1 for n in a.GetNeighbors() if n.GetSymbol() == "H") == 2
         and all(n.GetSymbol() in ("C", "H") for n in a.GetNeighbors())
         and sum(1 for n in a.GetNeighbors() if n.GetSymbol() == "H") == 2}
bulk = chain - alpha - beta
PRED["bulk CH2 (1.21-1.36)"] = 2 * len(bulk)
PUBNMR = {"CH3 (0.89 t)": 18, "Hbeta (1.53-1.69)": 12,
          "Halpha (2.16-2.37)": 12, "bulk CH2 (1.21-1.36)": 72}
agree = sum(1 for k in PUBNMR if PRED[k] == PUBNMR[k])
# seed: the same prediction from the 5-tailed congener
m5 = Chem.AddHs(Chem.MolFromSmiles(enu["series"]["5"]["smiles"]))
def pred_for(mol):
    a = {i[0] for i in mol.GetSubstructMatches(Chem.MolFromSmarts("[CH2]C(=O)O"))}
    b = {i[0] for i in mol.GetSubstructMatches(
        Chem.MolFromSmarts("[CH2]([CH2]C(=O)O)[CH2]"))}
    ch = {x.GetIdx() for x in mol.GetAtoms()
          if x.GetSymbol() == "C"
          and all(n.GetSymbol() in ("C", "H") for n in x.GetNeighbors())
          and sum(1 for n in x.GetNeighbors() if n.GetSymbol() == "H") == 2}
    me = {x.GetIdx() for x in mol.GetAtoms()
          if x.GetSymbol() == "C"
          and sum(1 for n in x.GetNeighbors() if n.GetSymbol() == "H") == 3}
    return {"CH3 (0.89 t)": 3 * len(me), "Hbeta (1.53-1.69)": 2 * len(b),
            "Halpha (2.16-2.37)": 2 * len(a),
            "bulk CH2 (1.21-1.36)": 2 * len(ch - a - b)}
P6, P5 = pred_for(mh), pred_for(m5)
qc["checks"]["P2_1H_NMR_integrals"] = {
    "what": "integrals predicted from our H9T6 structure vs the integrals "
            "printed on SI p.5",
    "published": PUBNMR, "predicted_6_tailed": P6,
    "agreement_count": sum(1 for k in PUBNMR if P6[k] == PUBNMR[k]),
    "SEEDED_predicted_5_tailed": P5,
    "SEEDED_agreement_count": sum(1 for k in PUBNMR if P5[k] == PUBNMR[k]),
    "not_checked": ("the 2.36-5.37 multiplet, printed as 46H; our structure "
                    "predicts 48H there, and that window overlaps the "
                    "separately-printed Halpha window, so the two integrals "
                    "are not independent and the 2H gap is not evidence of "
                    "anything either way"),
    "verdict": "4/4 exact on the independent windows; 1/4 under the seed."}

# -------------------------------------------- P3 the drawn 1.0 reference line
qc["checks"]["P3_drawn_but_never_printed_reference_line"] = {
    "what": "both panels carry a dashed horizontal rule the authors drew at the "
            "SM102 benchmark and never labelled. Calibrate on the y ticks "
            "ONLY, then predict the rule's value; a round 1.0 is the target.",
    "S4A_predicted": geo["panels"]["S4A"]["dashed_line_value_PREDICTED"],
    "S4C_predicted": geo["panels"]["S4C"]["dashed_line_value_PREDICTED"],
    "worst_error_pct": round(100 * max(
        abs(geo["panels"][p]["dashed_line_value_PREDICTED"] - 1.0)
        for p in ("S4A", "S4C")), 2)}
# seed: a 2% error in the scale
for p in ("S4A", "S4C"):
    g = geo["panels"][p]
    seeded = 1.0 + (g["dashed_line_value_PREDICTED"] - 0.0) * 0.02
    qc["checks"]["P3_drawn_but_never_printed_reference_line"][
        f"SEEDED_{p}_scale_+2pct"] = round(
        g["dashed_line_value_PREDICTED"] * 1.02, 4)
qc["checks"]["P3_drawn_but_never_printed_reference_line"]["verdict"] = (
    "0.35% and 0.02% off a round 1.0 on two independently calibrated panels; "
    "a 2% scale error would put them at 1.016 and 1.020, outside that.")

# -------------------------------- P4 SM102's own bar, 1.000 by construction
qc["checks"]["P4_SM102_bar_is_unity_by_construction"] = {
    "what": "every bar is divided by SM102, so SM102's own bar is exactly "
            "1.000. It is measured by the same detector as a closed loop.",
    "S4A_measured": geo["panels"]["S4A"]["bars"].get("SM102", {}).get("value"),
    "S4C_measured": geo["panels"]["S4C"]["bars"].get("SM102", {}).get("value"),
    "worst_error_pct": round(100 * max(
        abs(geo["panels"][p]["bars"]["SM102"]["value"] - 1.0)
        for p in ("S4A", "S4C") if "SM102" in geo["panels"][p]["bars"]), 2),
    "note": ("this is an END-TO-END check: axis, window, silhouette walk, "
             "stroke centre. It is NOT used to tune anything -- no constant "
             "in 04_digitize_figS4.py was fitted to it, which is the error "
             "class 27 trap."),
    "verdict": "both panels low by ~1.1-1.4%, a systematic sub-pixel bias; "
               "reported, not corrected."}

# ----------------------------------------------- P5 axis leave-one-out
qc["checks"]["P5_axis_leave_one_out"] = {
    "what": "axis.loo_residuals on both tick ladders. calibrate's max-residual "
            "test is weak against a single bad tick; LOO is not.",
    "LOO_TOL_PX": axis.LOO_TOL_PX,
    "S4A_y": geo["panels"]["S4A"]["loo_residual_px_max"],
    "S4C_y": geo["panels"]["S4C"]["loo_residual_px_max"],
    "S4A_x": geo["panels"]["S4A"]["x_tick_loo_max_px"],
    "S4C_x": geo["panels"]["S4C"]["x_tick_loo_max_px"]}
t = [float(v) for v in geo["panels"]["S4C"]["ticks_px"]]
t_bad = list(t); t_bad[0] += 6.0
qc["checks"]["P5_axis_leave_one_out"]["SEEDED_one_tick_moved_6px"] = round(
    max(axis.loo_residuals(t_bad)), 2)
qc["checks"]["P5_axis_leave_one_out"]["verdict"] = (
    "all four ladders <= 0.24 px; moving one tick by 6 px takes the worst LOO "
    "residual to the seeded value above, far past the 1.5 px tolerance.")

# ------------------------------------- P6 independent bar-top cross-check
xs = {}
for p in ("S4A", "S4C"):
    d = [b["bar_top_crosscheck_delta_px"] for b in geo["panels"][p]["bars"].values()]
    ok = [v for v in d if v is not None]
    xs[p] = {"bars_measured": len(d), "crosscheck_returned": len(ok),
             "exact": sum(1 for v in ok if v == 0),
             "within_1px": sum(1 for v in ok if abs(v) <= 1),
             "within_3px": sum(1 for v in ok if abs(v) <= 3),
             "worst_px": max((abs(v) for v in ok), default=None)}
qc["checks"]["P6_bar_top_crosscheck_vs_library"] = {
    "what": "bars.bar_top run top-down over the same windows, with its own "
            "full-width, flank and baseline guards, against the bottom-up "
            "silhouette walk that produced the values.",
    "per_panel": xs,
    "verdict": "two detectors with different failure modes; where the library "
               "one returns at all it agrees within 1 px on most bars."}

# ------------------------------------------------- P7 Table S1 arithmetic
import csv as _csv
rows = list(_csv.DictReader(open(os.path.join(HERE, "05_values.csv"))))
sums = [round(sum(float(r[c]) for c in ("IL_molpct", "PL_molpct",
                                        "chol_molpct", "PEG_molpct")), 6)
        for r in rows]
qc["checks"]["P7_formulation_mol_percent_sums"] = {
    "what": "Table S1's four lipid mol% must sum to 100 on every row",
    "rows": len(sums), "exact_100": sum(1 for v in sums if v == 100.0),
    "SEEDED_one_row_PEG_+1": sum(1 for i, v in enumerate(sums)
                                 if (v + 1.0 if i == 0 else v) == 100.0),
    "verdict": "all rows exact; the seed moves the count by one."}

# -------------------- P8 cross-panel: Fig S3A vs Fig S4C, same formulations
from scipy.stats import spearmanr, pearsonr
s3 = json.load(open(os.path.join(HERE, "04b_figS3A_ranks.json")))["ratio_to_SM102"]
s4c = {k: v["value"] for k, v in geo["panels"]["S4C"]["bars"].items()}
s4a = {k: v["value"] for k, v in geo["panels"]["S4A"]["bars"].items()}
ks = [f"1-F{i}" for i in range(1, 17)]
pr = [(k, s3[k], s4c[k]) for k in ks if k in s4c and k in s3]
A = np.array([x for _, x, _ in pr]); B = np.array([y for _, _, y in pr])
pa = [(k, s3[k], s4a[k]) for k in ks if k in s4a and k in s3]
qc["checks"]["P8_cross_panel_same_formulations"] = {
    "what": "Fig S3A and Fig S4C both report C2C12 luciferase relative to "
            "SM102 at 24 h for formulations 1-F1..1-F16, digitised "
            "independently (S3A by height ratio, no axis fit; S4C by a "
            "calibrated axis).",
    "n_shared": len(pr),
    "spearman": round(float(spearmanr(A, B).statistic), 4),
    "spearman_p": float(f"{spearmanr(A, B).pvalue:.3g}"),
    "pearson_on_log10": round(float(pearsonr(np.log10(A), np.log10(B))[0]), 4),
    "ratio_S3A_over_S4C": {"median": round(float(np.median(A / B)), 2),
                           "min": round(float((A / B).min()), 2),
                           "max": round(float((A / B).max()), 2)},
    "different_cell_line_control_S3A_vs_S4A": {
        "n": len(pa),
        "spearman": round(float(spearmanr([x for _, x, _ in pa],
                                          [y for _, _, y in pa]).statistic), 4)},
    "verdict": ("the two panels RANK-agree at rho = 0.94 (p = 8e-06) and "
                "disagree in MAGNITUDE by a median factor of ~9. Two "
                "independent reads of two different panels tracking each "
                "other that closely is strong evidence the reads are right; "
                "the factor of 9 is the paper's, not ours, and means "
                "'relative to SM102' does not denote the same quantity in "
                "the two figures. See 08_contradictions.csv C-04.")}

qc["summary"] = {
    "rows_emitted": len(rows),
    "panels": ["Figure S4A (DC2.4)", "Figure S4C (C2C12)"],
    "bars_measured": {p: len(geo["panels"][p]["bars"]) for p in ("S4A", "S4C")},
    "bars_refused": {p: len(geo["panels"][p]["refusals"]) for p in ("S4A", "S4C")},
    "distinct_ionizable_lipid_structures": 1,
    "corpus_overlap_full_InChIKey": 0}

with open(os.path.join(HERE, "07_qc.json"), "w") as fh:
    json.dump(qc, fh, indent=1)
for k, v in qc["checks"].items():
    print("\n##", k)
    for kk, vv in v.items():
        print(f"   {kk}: {vv}")
print("\n", qc["summary"])
