#!/usr/bin/env python3
"""Every check for this curation, with numbers, plus a seeded-defect harness.

    python3 -I stage07_qc.py      ->  07_qc.json   (exit 1 if a gate fails)

THE HARNESS RULE. A seeded mutation does not count as proving a check has teeth
unless the CHECK'S OWN METRIC MOVES. A stage that prints "MISMATCH" and still
exits 0 proves nothing, and neither does a seed that merely raises somewhere
downstream. So every seed below names the metric it must change and the
direction, and `seeded_defects` records before/after for each one.

WHAT NO CHECK HERE CAN SEE. Mass and molecular formula are blind to
connectivity: they fix how many arms and of what composition, never which
nitrogen an arm sits on. P2's molar masses and P1's HRMS are both mass. The only
test of topology against something outside our own transcription is the C10
anchor (make_reference.py), and the regiochemistry argument in P8, which is a
reading of three independent namings and not a measurement.
"""
import copy
import csv
import json
import os
import statistics
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", "lib"))
sys.path.insert(0, _HERE)
import paper_library as P                                       # noqa: E402
import openpyxl                                                 # noqa: E402
from rdkit import Chem, RDLogger                                # noqa: E402
from rdkit.Chem.Descriptors import MolWt                        # noqa: E402

RDLogger.DisableLog("rdApp.*")
XLSX = os.path.join(_HERE, "source", "41467_2026_69077_MOESM6_ESM.xlsx")

# ---------------------------------------------------------------- P1 inputs
# (nothing to transcribe; stage03 already wrote the A/B/C columns)

# ------------- P2: masses the authors WEIGHED, transcribed from the SI text
# Every one is "(<mg> mg, <n> mmol)" beside a named compound, so mg/mmol is an
# implied molar mass that the authors never printed as a molar mass. This is
# external to our pipeline in the way a recomputed marginal is not -- and it is
# the ONLY independent atom-count check CBBZ0 has, since the SI gives it no
# HRMS.
#
# THE TOLERANCE IS THE PRINTED PRECISION, NOT A PERCENTAGE. "0.23 mmol" means
# the true amount is in [0.225, 0.235), so the implied molar mass is an
# INTERVAL, and the check is containment in it. A flat percentage was the first
# version and it was wrong in both directions at once: 1.5% rejected a correct
# CBB22 whose mole figure carries only two significant figures, while being far
# looser than it needed to be on "0.184 mmol", where the interval is only
# +/-0.27%. Quantities are kept as the STRINGS the SI prints so the precision
# cannot be lost on the way in -- reading 0.23 as a float and then guessing how
# many digits it had is how a tolerance stops meaning anything.
# A reagent NAME is not a structure, and a WEIGHING is of whatever salt form was
# in the bottle. `weighed_form` says which, per row -- it is not a tolerance.
# CBB22 is purchased Brilliant Blue G, and the SI's own "854 mg, 1.0 mmol"
# is 21.94 above our emitted free acid at 832.06, i.e. Na-for-H (21.9819) to
# within the authors' rounding. So the free acid is what we emit (the HRMS of
# CBB5, CBB16 and all ten CBBZ compounds is [M+H]+ of the neutral free acid and
# matched exactly) while the mono-sodium salt is what was weighed. The check
# derives the salt's mass from our own structure by one declared substitution,
# and `salt_delta_is_Na_minus_H` below asserts that substitution is what it says.
NA_MINUS_H = 22.98977 - 1.008
COCL_MINUS_H = 12.011 + 15.999 + 35.45 - 1.008
WEIGHINGS = [
 ("CBB22", "854", "1.0", "mmol",
  "SI p18, 'commercial brilliant G250 (854 mg, 1.0 mmol)'",
  "mono_sodium_salt"),
 ("CBB22", "200.0", "0.23", "mmol", "SI p28, CBB15/16 synthesis", "mono_sodium_salt"),
 ("CBB21", "200.0", "0.184", "mmol", "SI p19, CBB1 synthesis", "free_acid"),
 ("CBB21", "100.0", "0.092", "mmol", "SI p22, CBB6 synthesis", "free_acid"),
 ("CBB21", "30.0", "27.6", "umol", "SI p26, CBB12 synthesis", "free_acid"),
 ("CBB21", "100.0", "92", "umol", "SI p27, CBB13 synthesis", "free_acid"),
 ("CBB21", "40.0", "36.9", "umol", "SI p30, CBB20 synthesis", "free_acid"),
 ("CBB21", "20.0", "18.4", "umol", "SI p24, CBB10 synthesis", "free_acid"),
 ("CBB1", "20.0", "10.5", "umol", "SI p19, CBB2 synthesis", "free_acid"),
 ("CBB3", "20.0", "11.7", "umol", "SI p21, CBB4 synthesis", "free_acid"),
 ("CBB6", "16.0", "10.2", "umol", "SI p22, CBB7 synthesis", "free_acid"),
 ("CBB13", "15.0", "12.2", "umol", "SI p27, CBB14 synthesis", "free_acid"),
 ("CBBZ0", "15.0", "17.8", "umol", "SI pp43-49, CBBZ1/3/4/5/6/8/10 syntheses", "free_acid"),
 ("CBBZ0", "20.0", "23.8", "umol", "SI p44, CBBZ2 synthesis", "free_acid"),
 ("CBBZ0", "10.0", "11.9", "umol", "SI pp48-49, CBBZ7/CBBZ9 syntheses", "free_acid"),
]
# The TAIL REAGENTS, same idea. These matter more than they look: every tail
# acid here is pinned by a weighing that is independent of the HRMS formula
# strings, so a wrong tail cannot hide behind a wrong formula. Retinoic acid is
# the case that needs it -- the SI's amide name reads "2,6-dimethyl-9-(2,3,3-
# trimethylcyclohex-1-en-1-yl)" where all-trans-retinoic acid is 3,7-dimethyl-
# 9-(2,6,6-trimethylcyclohex-1-en-1-yl), and those are ISOMERS, so no mass can
# separate them; what settles it is that the SI names the REAGENT as
# all-trans-retinoic acid and weighs it at a mass only retinoic acid has.
# `delta` is a declared substitution, never a fudge factor: "sodium_salt" is
# -H +Na and "chloroformate" is -H +C(=O)Cl, and the test still has to land
# inside the printed interval afterwards.
REAGENT_WEIGHINGS = [
 ("butyryl", "sodium butyrate", "50", "0.5", "mmol", "sodium_salt", "SI p26"),
 ("lauroyl", "lauric acid", "20", "0.1", "mmol", "none", "SI p26"),
 ("palmitoyl", "palmitic acid", "256", "1", "mmol", "none", "SI p28"),
 ("triacontanoyl", "melissic acid", "25", "55", "umol", "none", "SI p24"),
 ("15-hexadecynoyl", "alkynyl palmitic acid", "4.0", "15.8", "umol", "none",
  "SI p23"),
 ("linoleoyl", "linoleic acid", "500", "1.78", "mmol", "none", "SI p25"),
 ("retinoyl", "all-trans-retinoic acid", "50", "166", "umol", "none", "SI p30"),
]
CHOL_WEIGHING = ("cholesteryl chloroformate", "400", "0.89", "mmol",
                 "chloroformate", "SI p19, the CBB1 synthesis")

# the azide partners, same idea -- one weighing per CBBZ tail
AZIDE_WEIGHINGS = [
 ("Z1", "10", "17.2"), ("Z2", "10", "15.3"), ("Z3", "10", "14.3"),
 ("Z4", "10", "13.8"), ("Z5", "10", "13.6"), ("Z6", "10", "13.4"),
 ("Z7", "10", "11.7"), ("Z8", "10", "13.0"), ("Z9", "10", "10.3"),
 ("Z10", "10", "12.0"),
]

# ------------- P3: the paper's OWN hit criterion and its OWN hit list
HIT_RULE = {
 "Fig2c": dict(control="C",
   criterion="legend: 'considered positive hits if the ratio of red-to-green "
             "(RFP/GFP) with flow cytometry was more than twice that of C'",
   stated=["CBB6", "CBB11", "CBB14", "CBB18", "CBB19", "CBB21"],
   stated_from="Results: 'Six compounds from the amide coupled library (CBB6, "
               "CBB11, CBB14, CBB18, CBB19, and CBB21) demonstrated the "
               "ability to enhance delivery of Cre'"),
 "Fig4c": dict(control="A",
   criterion="legend: 'considered positive hits if the ratio of green-to-red "
             "(eGFP/mCherry) was more than twice that of A'",
   stated=["CBB3", "CBB11"],
   stated_from="Results: 'among CBB1-22 and CBBZ1-Z7, consistent eGFP "
               "expression in reporter cells was only achieved by CBB3 "
               "(3.9 +/- 0.4%) and CBB11 (28.6 +/- 0.7%)'"),
}

# ------------- P4: mean +/- SD statements printed in the Results
PRINTED_STATS = [
 ("Fig2c", "CBB19", 22.4, 0.6, "Results: 'up to 22.4 +/- 0.6% GFP-to-RFP conversion'"),
 ("Fig4c", "CBB3", 3.9, 0.4, "Results: 'CBB3 (3.9 +/- 0.4%)'"),
 ("Fig4c", "CBB11", 28.6, 0.7, "Results: 'CBB11 (28.6 +/- 0.7%)'"),
]
# Fig 6a, in vivo: not emitted, but reading it tests the workbook layout
PRINTED_STATS_FIG6A = [
 ("TAS", 45.8, 12.9), ("DMSO", 50.6, 33.9), ("CBB11", 60.8, 14.9),
 ("CBB14", 58.5, 37.4), ("Lipo1", 227.0, 69.0), ("Lipo2", 95.3, 62.5),
]
LEGEND_N = {"Fig2c": 4, "Fig2d": 2, "Fig4c": 4, "Fig4d": 3}


# =========================================================== the checks
def p1_hrms(enum_rows):
    tot = [r for r in enum_rows if r["hrms_formula"]]
    def ok_num(v):
        try:
            return abs(float(v)) <= 0.002
        except ValueError:
            return False
    A = sum(1 for r in tot if r["check_A_formula"] == "exact")
    B = sum(1 for r in tot if ok_num(r["check_B_mz_from_si_formula"]))
    C = sum(1 for r in tot if ok_num(r["check_C_mz_from_our_formula"]))
    disagree = [dict(compound=r["compound"], A=r["check_A_formula"],
                     B=r["check_B_mz_from_si_formula"],
                     C=r["check_C_mz_from_our_formula"])
                for r in tot
                if r["check_A_formula"] != "exact"
                or not ok_num(r["check_B_mz_from_si_formula"])
                or not ok_num(r["check_C_mz_from_our_formula"])]
    return dict(n_with_hrms=len(tot), A_formula_exact=A,
                B_si_formula_matches_si_mz=B, C_our_formula_matches_si_mz=C,
                disagreements=disagree)


def _interval(s):
    """[lo, hi) implied by the number of decimals the SI actually printed."""
    half = 0.5 * 10 ** -(len(s.split(".")[1]) if "." in s else 0)
    return float(s) - half, float(s) + half


def _mw_interval(mg_s, mol_s, unit):
    """The molar-mass interval the printed quantities permit."""
    scale = 1e-3 if unit == "mmol" else 1e-6
    mg_lo, mg_hi = _interval(mg_s)
    mol_lo, mol_hi = _interval(mol_s)
    return ((mg_lo / 1000.0) / (mol_hi * scale),
            (mg_hi / 1000.0) / (mol_lo * scale))


def p2_weighings(smiles_by_name):
    rows, n_ok = [], 0
    for name, mg_s, mol_s, unit, where, form in WEIGHINGS:
        smi = smiles_by_name.get(name, "")
        if not smi:
            continue
        emitted = MolWt(Chem.MolFromSmiles(smi))
        if form == "mono_sodium_salt":
            ours = emitted + NA_MINUS_H
        elif form == "free_acid":
            ours = emitted
        else:
            raise SystemExit("unknown weighed_form: " + form)
        lo, hi = _mw_interval(mg_s, mol_s, unit)
        inside = lo <= ours <= hi
        rows.append(dict(kind="compound", name=name, weighed_form=form,
                         printed_mg=mg_s, printed_amount=mol_s + " " + unit,
                         implied_MW_interval=[round(lo, 2), round(hi, 2)],
                         our_emitted_MW=round(emitted, 2),
                         our_weighed_form_MW=round(ours, 2),
                         rel_error=round(abs(ours - (lo + hi) / 2)
                                         / ((lo + hi) / 2), 5),
                         ok=inside, where=where))
        n_ok += inside
    for zn, mg_s, mol_s in AZIDE_WEIGHINGS:
        tail = P.AZIDE_TAIL[zn][0]
        ours = MolWt(Chem.MolFromSmiles("[N-]=[N+]=N" + tail))
        lo, hi = _mw_interval(mg_s, mol_s, "umol")
        inside = lo <= ours <= hi
        rows.append(dict(kind="azide_tail", name=zn + "-azide",
                         weighed_form="neutral azide", printed_mg=mg_s,
                         printed_amount=mol_s + " umol",
                         implied_MW_interval=[round(lo, 2), round(hi, 2)],
                         our_emitted_MW=round(ours, 2),
                         our_weighed_form_MW=round(ours, 2),
                         rel_error=round(abs(ours - (lo + hi) / 2)
                                         / ((lo + hi) / 2), 5),
                         ok=inside,
                         where="SI, the CBBZ%s click step" % zn[1:]))
        n_ok += inside
    # --- the tail reagents
    import reagents as _rg
    for key, label, mg_s, mol_s, unit, delta, where in REAGENT_WEIGHINGS:
        acid = P.ACYL[key][0].replace("[CH:1](=O)", "OC(=O)")
        base = MolWt(Chem.MolFromSmiles(acid))
        ours = base + (NA_MINUS_H if delta == "sodium_salt" else 0.0)
        lo, hi = _mw_interval(mg_s, mol_s, unit)
        inside = lo <= ours <= hi
        rows.append(dict(kind="tail_reagent", name=label, tail=key,
                         weighed_form=delta, printed_mg=mg_s,
                         printed_amount=mol_s + " " + unit,
                         implied_MW_interval=[round(lo, 2), round(hi, 2)],
                         our_emitted_MW=round(base, 2),
                         our_weighed_form_MW=round(ours, 2),
                         rel_error=round(abs(ours - (lo + hi) / 2)
                                         / ((lo + hi) / 2), 5),
                         ok=inside, where=where))
        n_ok += inside
    label, mg_s, mol_s, unit, delta, where = CHOL_WEIGHING
    chol = MolWt(Chem.MolFromSmiles(_rg.smi("Cholesterol")))
    ours = chol + COCL_MINUS_H
    lo, hi = _mw_interval(mg_s, mol_s, unit)
    inside = lo <= ours <= hi
    rows.append(dict(kind="tail_reagent", name=label, tail="cholesteryl carbamate",
                     weighed_form=delta, printed_mg=mg_s,
                     printed_amount=mol_s + " " + unit,
                     implied_MW_interval=[round(lo, 2), round(hi, 2)],
                     our_emitted_MW=round(chol, 2),
                     our_weighed_form_MW=round(ours, 2),
                     rel_error=round(abs(ours - (lo + hi) / 2)
                                     / ((lo + hi) / 2), 5),
                     ok=inside, where=where))
    n_ok += inside

    free = MolWt(Chem.MolFromSmiles(smiles_by_name["CBB22"]))
    return dict(n_weighings=len(rows), n_inside_printed_interval=n_ok,
                worst_rel_error=round(max(r["rel_error"] for r in rows), 5),
                salt_delta_is_Na_minus_H=round(NA_MINUS_H, 4),
                cbb22_free_acid_MW=round(free, 2),
                cbb22_sodium_salt_MW=round(free + NA_MINUS_H, 2),
                cbb22_si_weighing_MW=854.0,
                salt_form_note="The SI weighs 854 mg as 1.0 mmol. Our free "
                               "acid is 832.06, and 854.04 is that plus "
                               "exactly Na-for-H. The salt form is therefore "
                               "read off the authors' own weighing, not "
                               "assumed, and no tolerance was widened.",
                detail=rows)


def p3_hits(vals, controls):
    out = {}
    for panel, spec in HIT_RULE.items():
        c = controls[panel][spec["control"]]["mean"]
        thr = 2.0 * c
        got = sorted(r["compound"] for r in vals
                     if r["panel"] == panel and float(r["value_mean"]) > thr)
        out[panel] = dict(
            control=spec["control"], control_mean=round(c, 4),
            threshold=round(thr, 4), criterion=spec["criterion"],
            stated_hits=spec["stated"], stated_from=spec["stated_from"],
            reproduced_hits=got,
            n_stated=len(spec["stated"]), n_reproduced=len(got),
            exact_set_match=sorted(spec["stated"]) == got,
            missing=[x for x in spec["stated"] if x not in got],
            extra=[x for x in got if x not in spec["stated"]])
    return out


def p4_printed(vals, grid6):
    rows, n_ok = [], 0
    by = {(r["panel"], r["compound"]): r for r in vals}
    for panel, cmp_, mean, sd, where in PRINTED_STATS:
        r = by[(panel, cmp_)]
        reps = [float(x) for x in r["replicates"].split(";")]
        m = statistics.fmean(reps)
        s_samp = statistics.stdev(reps)
        s_pop = statistics.pstdev(reps)
        ok = (round(m, 1) == mean and round(s_samp, 1) == sd)
        rows.append(dict(locus=panel, compound=cmp_, printed_mean=mean,
                         printed_sd=sd, our_mean=round(m, 4),
                         our_sample_sd=round(s_samp, 4),
                         our_population_sd=round(s_pop, 4),
                         sd_convention="sample (ddof=1)" if ok else "UNRESOLVED",
                         ok=ok, where=where))
        n_ok += ok
    # Fig 6a: population SD, which is a DIFFERENT convention from Figs 2c/4c
    hdr = grid6[1]
    f6, n6 = [], 0
    for name, mean, sd in PRINTED_STATS_FIG6A:
        col = list(hdr).index(name)
        reps = [float(g[col]) for g in grid6[2:12]
                if col < len(g) and g[col] not in (None, "")]
        m, ss, sp = (statistics.fmean(reps), statistics.stdev(reps),
                     statistics.pstdev(reps))
        okp = abs(sp - sd) <= 0.55 and abs(m - mean) <= 0.55
        oks = abs(ss - sd) <= 0.55
        f6.append(dict(group=name, n=len(reps), printed_mean=mean,
                       printed_sd=sd, our_mean=round(m, 3),
                       our_sample_sd=round(ss, 3),
                       our_population_sd=round(sp, 3),
                       matches_population_sd=okp, matches_sample_sd=oks))
        n6 += okp
    return dict(n_statements=len(rows), n_reproduced=n_ok, detail=rows,
                fig6a_n_statements=len(f6), fig6a_n_reproduced_pop_sd=n6,
                fig6a_detail=f6,
                fig6a_best_lipo1_uV=max(
                    float(g[list(hdr).index("Lipo1")]) for g in grid6[2:12]
                    if g[list(hdr).index("Lipo1")] not in (None, "")),
                fig6a_printed_best="353 uV, 'the best example of Lipo1-mediated "
                                   "editing achieved a b-wave amplitude of 353 uV'")


def p5_cross_panel(vals):
    by = {(r["panel"], r["compound"]): float(r["value_mean"]) for r in vals}
    pairs = []
    for a, b, c in (("Fig2c", "Fig2d", "CBB14"), ("Fig4c", "Fig4d", "CBB11")):
        if (a, c) in by and (b, c) in by:
            pairs.append(dict(compound=c, flow_panel=a, flow=by[(a, c)],
                              microscopy_panel=b, microscopy=by[(b, c)],
                              microscopy_over_flow=round(by[(b, c)] /
                                                         by[(a, c)], 3)))
    return dict(
        what="the SAME compound, cargo and cell line read by flow cytometry "
             "and by microscopy in the same figure",
        pairs=pairs,
        conclusion="the two readouts disagree by about 2x in both directions, "
                   "so Fig 2c/4c and Fig 2d/4d are four separate assay blocks "
                   "and must not be pooled into one screen")


def p6_mapping(vals, controls):
    out = {}
    for panel in LEGEND_N:
        cs = [r["compound"] for r in vals if r["panel"] == panel]
        out[panel] = dict(n_compounds=len(cs), n_distinct=len(set(cs)),
                          injective=len(cs) == len(set(cs)),
                          compounds=sorted(cs),
                          controls_excluded=sorted(controls[panel]))
    full = sorted(r["compound"] for r in vals if r["panel"] == "Fig2c")
    out["Fig2c_covers_CBB1_to_CBB22"] = (
        full == sorted("CBB%d" % i for i in range(1, 23)))
    out["Fig4c_covers_CBB1_to_CBB22"] = (
        sorted(r["compound"] for r in vals if r["panel"] == "Fig4c")
        == sorted("CBB%d" % i for i in range(1, 23)))
    return out


def p7_replicates(vals):
    rows, n_ok = [], 0
    for panel, n in LEGEND_N.items():
        got = sorted({int(r["n_replicates"]) for r in vals
                      if r["panel"] == panel})
        ok = got == [n]
        rows.append(dict(panel=panel, legend_n=n, sheet_n=got, ok=ok))
        n_ok += ok
    return dict(n_panels=len(rows), n_agreeing=n_ok, detail=rows,
                legend_source="Fig 2 legend: '(c) mean +/- SD; n = 4 "
                              "replicates and (d) mean, n = 2 replicates'; "
                              "Fig 4 legend: 'n = 4 replicates (c) and 3 "
                              "replicates (d)'")


# --------- P10: the DRAWN Markush table vs what we built from the NAMES
# Two independent statements of the same assignment. The builder read the IUPAC
# names on SI pp 18-50; this reads Supplementary Figure 1's "(Lipid1, Lipid2,
# Linker)" table and asks, by substructure search on the BUILT molecule,
# whether the blocks the figure names are actually present and in the right
# multiplicity. Anchored queries only: a bare C16 acyl chain is a substructure
# of a C30 one, so every saturated tail is matched from its carbonyl THROUGH to
# its terminal methyl, which makes the match exact-length.
def _chain(n, tail="[CH3]"):
    return "[CX3](=[OX1])" + "[CH2]" * (n - 2) + tail


LIPID_QUERY = {
 "L2":  ("smiles", P.ACYL["dha"][0].replace("[CH:1]", "C")),
 "L3":  ("smarts", _chain(16)),
 # hexadec-15-ynoate: carbonyl C1, thirteen CH2 (C2-C14), then C15#C16.
 # The first version wrote twelve and scored 0/2 on CBB8 -- the check caught
 # its own query, which is the only reason it is trustworthy on the other 20.
 "L4":  ("smarts", _chain(15, "C#[CH]")),
 "L5":  ("smarts", _chain(30)),
 "L6":  ("smiles", P.ACYL["linoleoyl"][0].replace("[CH:1]", "C")),
 "L7":  ("smarts", _chain(12)),
 "L8":  ("smarts", _chain(4)),
 "L9":  ("smarts", "[CX3](=[OX1])C(F)(F)F"),
 "L10": ("smiles", P.ACYL["retinoyl"][0].replace("[CH:1]", "C")),
 "L12": ("smarts", "[SX4](=O)(=O)[OX2H1]"),
 "H":   ("smarts", "[NX3;H2;!$(NC=O)]"),
 "-COCO-": ("smarts", "[NX3]C(=O)C(=O)[NX3]"),
 # the only carbamate in the library, and it sits on a secondary ring
 # carbon -- the sterol C3. Nothing else in the corpus of 34 has N-C(=O)-O.
 "L1":  ("smarts", "[NX3]C(=O)O[CX4H1]"),
 "L11": ("smarts", "[CH2][CH2]SC"),          # methionine side chain, 1 per CPP5
}
LINKER_QUERY = {
 "K1": ("smarts", "[CH2][NX3]([CH3])[CH2]"),
 "K2": ("smarts", "[CH2][N+]([CH3])([CH3])[CH2]"),
 "K3": ("smarts", "[SX4](=O)(=O)O[CH2]"),
}


def p10_drawn_table(smiles_by_name):
    rows, n_ok = [], 0
    for cmpd, (l1, l2, k) in sorted(P.DRAWN_TABLE.items(),
                                    key=lambda kv: int(kv[0][3:])):
        smi = smiles_by_name.get(cmpd) or smiles_by_name.get(cmpd + "_KLVPM")
        if not smi:
            rows.append(dict(compound=cmpd, drawn=(l1, l2, k),
                             ok=False, why="no structure built"))
            continue
        mol = Chem.MolFromSmiles(smi)
        want = {}
        for code in (l1, l2):
            if code is not None:
                want[code] = want.get(code, 0) + 1
        if l2 is None:                       # the "9 (-COCO-, K1)" cell
            want = {l1: 1}
        got, detail = True, {}
        for code, mult in want.items():
            kind, q = LIPID_QUERY[code]
            qm = (Chem.MolFromSmiles(q) if kind == "smiles"
                  else Chem.MolFromSmarts(q))
            if qm is None:
                raise SystemExit("P10: query for %s did not parse: %r"
                                 % (code, q))
            n = len(mol.GetSubstructMatches(qm, useChirality=False))
            detail[code] = dict(expected=mult, found=n)
            got = got and n == mult
        kind, q = LINKER_QUERY[k]
        qm = Chem.MolFromSmarts(q)
        if qm is None:
            raise SystemExit("P10: linker query %s did not parse: %r" % (k, q))
        nk = len(mol.GetSubstructMatches(qm))
        # K1/K2/K3 appear once per DERIVATISED sulfonyl arm. Only the L12
        # cells (CBB5, CBB16) leave one arm as free sulfonic acid, so only
        # they expect 1. CBB9's "(-COCO-, K1)" has l2 = None because one
        # bridge spans BOTH arms -- it still has two linkers.
        want_k = 1 if l2 == "L12" else 2
        detail[k] = dict(expected=want_k, found=nk)
        got = got and nk == want_k
        rows.append(dict(compound=cmpd, drawn=[l1, l2, k], ok=got,
                         detail=detail))
        n_ok += got
    return dict(n_cells_drawn=len(P.DRAWN_TABLE), n_agreeing=n_ok,
                source="Supplementary Figure 1, SI p2, rendered at 400 dpi and "
                       "read by eye; renders/crop_table.png and "
                       "renders/crop_linkers.png",
                what_this_adds="the builder used the IUPAC names; this uses the "
                               "drawn Markush table, so a compound-to-tail "
                               "mix-up would have to be made identically twice "
                               "from two different representations",
                detail=rows)


# --------- P11: the CBBZ drawn table (SI Fig 3) vs the azide tails from text
def _tail_from_drawing(r1, r2, n1, n2, x):
    """Compose the azide from the FIGURE's blocks, by the figure's own rule."""
    def arm(n, r):
        return "C" * (n + 1) + "C(=O)O" + P.DRAWN_R[r]
    return ("[N-]=[N+]=N" + P.DRAWN_X[x] + "N(" + arm(n1, r1) + ")"
            + arm(n2, r2))


def p11_drawn_cbbz_table():
    rows, n_ok, n_expected_mismatch = [], 0, 0
    for zn, cell in sorted(P.DRAWN_TABLE_Z.items(), key=lambda kv: int(kv[0][1:])):
        drawn = Chem.MolToSmiles(Chem.MolFromSmiles(_tail_from_drawing(*cell)))
        ours = Chem.MolToSmiles(Chem.MolFromSmiles(
            "[N-]=[N+]=N" + P.AZIDE_TAIL[zn][0]))
        same = drawn == ours
        exp = P.DRAWN_Z_EXPECTED_MISMATCH.get(zn, "")
        rows.append(dict(compound="CBBZ" + zn[1:],
                         drawn_cell=dict(R1=P.DRAWN_R_NAME[cell[0]],
                                         R2=P.DRAWN_R_NAME[cell[1]],
                                         n1=cell[2], n2=cell[3],
                                         X=P.DRAWN_X_NAME[cell[4]]),
                         agrees=same, expected_mismatch=exp,
                         resolved=bool(same or exp)))
        n_ok += same
        n_expected_mismatch += bool(exp and not same)
    return dict(n_cells_drawn=len(P.DRAWN_TABLE_Z), n_agreeing=n_ok,
                n_declared_mismatches=n_expected_mismatch,
                n_unresolved=len(P.DRAWN_TABLE_Z) - n_ok - n_expected_mismatch,
                source="Supplementary Figure 3, SI p4, rendered at 400 dpi and "
                       "read by eye; renders/crop_fig3.png",
                drawn_core_protonation="SI Fig 3 draws the dye as the "
                       "N-protonated cation with both sulfonic acids as "
                       "HO-SO2-. We emit the neutral quinoid-imine free acid, "
                       "which is what the CBBZ IUPAC names say and what the "
                       "[M+H]+ HRMS confirms. A charge-state difference, not a "
                       "connectivity one.",
                detail=rows)


# =========================================================== seeded defects
def seeded(enum_rows, vals, controls, smiles_by_name, grid6):
    """Each seed must move the named metric. Raising is not enough."""
    out = []

    def record(name, metric, before, after, direction, what):
        moved = (after < before if direction == "down" else after != before)
        out.append(dict(seed=name, metric=metric, before=before, after=after,
                        direction=direction, metric_moved=bool(moved),
                        what_was_seeded=what))

    # --- P1: break one built formula by one hydrogen (the CBB12 defect class)
    e2 = copy.deepcopy(enum_rows)
    tgt = [r for r in e2 if r["compound"] == "CBB11"][0]
    tgt["check_A_formula"] = 'MISMATCH {"H": [144, 143]}'
    tgt["check_C_mz_from_our_formula"] = "1.0078"
    record("P1_one_hydrogen_wrong_on_CBB11", "A_formula_exact",
           p1_hrms(enum_rows)["A_formula_exact"],
           p1_hrms(e2)["A_formula_exact"], "down",
           "a single H miscount on one compound, the exact class that is live "
           "in the SI's own CBB12 formula string")

    # --- P2: a tail two carbons short (the classic chain-length slip)
    s2 = dict(smiles_by_name)
    s2["CBB21"] = Chem.MolToSmiles(Chem.MolFromSmiles(
        smiles_by_name["CBB21"].replace("CCCN(C)CCCN", "CCN(C)CCN", 1)))
    record("P2_K1_linker_shortened_by_two_carbons", "n_inside_printed_interval",
           p2_weighings(smiles_by_name)["n_inside_printed_interval"],
           p2_weighings(s2)["n_inside_printed_interval"], "down",
           "CBB21's diaminopropyl linker written as diaminoethyl, so six of "
           "the weighing checks should stop agreeing")

    # --- P2b: an azide tail with the wrong ester chain
    import paper_library as PL
    keep = PL.AZIDE_TAIL["Z3"]
    before = p2_weighings(smiles_by_name)["n_inside_printed_interval"]
    PL.AZIDE_TAIL["Z3"] = (keep[0].replace("CCCCCCCCCCCCCCCC",
                                           "CCCCCCCCCCCC"), keep[1])
    after = p2_weighings(smiles_by_name)["n_inside_printed_interval"]
    PL.AZIDE_TAIL["Z3"] = keep
    record("P2b_CBBZ3_tail_hexadecyl_swapped_for_dodecyl", "n_inside_printed_interval",
           before, after, "down",
           "CBBZ3's dihexadecyl ester written as didodecyl -- a four-carbon "
           "slip on each of two arms")

    # --- P3: shift one compound below the paper's own hit threshold
    v2 = copy.deepcopy(vals)
    for r in v2:
        if r["panel"] == "Fig2c" and r["compound"] == "CBB6":
            r["value_mean"] = "0.1"
    b = p3_hits(vals, controls)["Fig2c"]["n_reproduced"]
    a = p3_hits(v2, controls)["Fig2c"]["n_reproduced"]
    record("P3_CBB6_value_zeroed", "Fig2c.n_reproduced", b, a, "down",
           "one true hit pushed under the threshold, so the reproduced hit set "
           "must stop matching the six compounds the paper names")

    # --- P3b: the control column misread, which moves the threshold
    c2 = copy.deepcopy(controls)
    c2["Fig4c"]["A"]["mean"] = 0.02
    b = p3_hits(vals, controls)["Fig4c"]["n_reproduced"]
    a = p3_hits(vals, c2)["Fig4c"]["n_reproduced"]
    record("P3b_Fig4c_control_A_misread_as_NE", "Fig4c.n_reproduced", b, a,
           "up_or_down",
           "the control mean replaced by the no-enzyme well's, which is the "
           "mistake a mis-shifted column header would make; the hit count has "
           "to change because the threshold is 2x the control")

    # --- P4: perturb one replicate so a printed mean no longer rounds right
    v3 = copy.deepcopy(vals)
    for r in v3:
        if r["panel"] == "Fig4c" and r["compound"] == "CBB11":
            reps = [float(x) for x in r["replicates"].split(";")]
            reps[0] += 2.0
            r["replicates"] = ";".join("%g" % x for x in reps)
    record("P4_one_CBB11_replicate_shifted_by_2pct", "n_reproduced",
           p4_printed(vals, grid6)["n_reproduced"],
           p4_printed(v3, grid6)["n_reproduced"], "down",
           "a 2-point change in one of four replicates, which is the size of "
           "error a one-row header offset would produce")

    # --- P7: a replicate row dropped
    v4 = [r for r in copy.deepcopy(vals)]
    for r in v4:
        if r["panel"] == "Fig2c":
            r["n_replicates"] = 3
    record("P7_Fig2c_read_as_three_replicates", "n_agreeing",
           p7_replicates(vals)["n_agreeing"], p7_replicates(v4)["n_agreeing"],
           "down", "one of Fig 2c's four replicate rows skipped, which is what "
                   "an off-by-one first_data_row does")

    # --- P10: two compounds' tails swapped, which is the mix-up the drawn
    # table exists to catch and which no mass check can see (CBB12 laurate and
    # CBB13 butyrate have different formulas, so swap two SAME-formula cells
    # instead: CBB6 and CBB15 are both palmitate, differing only in linker)
    s10 = dict(smiles_by_name)
    s10["CBB6"], s10["CBB12"] = s10["CBB12"], s10["CBB6"]
    record("P10_CBB6_and_CBB12_tails_swapped", "n_agreeing",
           p10_drawn_table(smiles_by_name)["n_agreeing"],
           p10_drawn_table(s10)["n_agreeing"], "down",
           "two compounds' structures exchanged -- the drawn table says CBB6 "
           "is palmitate and CBB12 laurate, so both cells must stop agreeing")

    # --- P11: undeclare the known Z10 mismatch; the unresolved count must rise
    keep10 = dict(P.DRAWN_Z_EXPECTED_MISMATCH)
    before = p11_drawn_cbbz_table()["n_unresolved"]
    P.DRAWN_Z_EXPECTED_MISMATCH.clear()
    after = p11_drawn_cbbz_table()["n_unresolved"]
    P.DRAWN_Z_EXPECTED_MISMATCH.update(keep10)
    record("P11_Z10_mismatch_declaration_removed", "n_unresolved",
           before, after, "up_or_down",
           "the declared Z10 exception deleted. This seeds the check's real "
           "failure mode -- a figure/text disagreement that nobody wrote down "
           "-- and the unresolved count has to go from 0 to 1")

    # --- P6: a column header shifted by one
    v5 = copy.deepcopy(vals)
    for r in v5:
        if r["panel"] == "Fig2c" and r["compound"] == "CBB22":
            r["compound"] = "CBB21"
    b = p6_mapping(vals, controls)["Fig2c"]["injective"]
    a = p6_mapping(v5, controls)["Fig2c"]["injective"]
    record("P6_two_Fig2c_columns_mapped_to_one_compound", "Fig2c.injective",
           b, a, "up_or_down",
           "a duplicated compound id, which is what a shifted header row gives")
    return out


def main():
    enum_rows = list(csv.DictReader(open(os.path.join(_HERE,
                                                      "03_enumerated.csv"))))
    vals = list(csv.DictReader(open(os.path.join(_HERE, "05_values.csv"))))
    controls = json.load(open(os.path.join(_HERE, "05_controls.json")))
    smiles_by_name = {r["compound"]: r["canonical_SMILES"]
                      for r in enum_rows if r["canonical_SMILES"]}
    wb = openpyxl.load_workbook(XLSX, read_only=True, data_only=True)
    grid6 = [list(r) for r in wb["Figure 6"].iter_rows(values_only=True)]

    qc = {
     "paper": "Zhang et al., Nat Commun 2026, PMID 41654520, PMC12992809",
     "value_tier": "1_publisher_source_data for every emitted value; "
                   "4_pdfplumber_text (SI IUPAC names + HRMS) for every "
                   "structure. No panel was digitised.",
     "cargo": "Cre recombinase protein (Fig 2) and ABE8e RNP (Fig 4). "
              "NOT mRNA -- there is no mRNA in this paper.",
     "library": dict(
       named_delivery_agents=34,
       breakdown="CBB1-CBB22 (22, of which CBB22 is commercial CBB G-250 and "
                 "CBB21 is the un-acylated bis-amine core) + CBBox11 + CBBZ0 "
                 "+ CBBZ1-CBBZ10",
       # Count COMPOUNDS, not rows. 03_enumerated.csv has 36 rows because
       # CBB18 and CBB19 each carry two sequence candidates, and counting
       # rows made "35 structures" sit next to "34 named compounds", which
       # reads as though one extra molecule appeared from somewhere.
       compounds_with_one_definite_structure=len(
           [k for k in smiles_by_name if "_" not in k]),
       compounds_with_competing_candidates=sorted(
           {k.split("_")[0] for k in smiles_by_name if "_" in k}),
       compounds_with_no_structure=["CBBox11"],
       enumeration_rows=len(enum_rows),
       structures_withheld=["CBBox11 (hydroxylation site undetermined; "
                            "DECISIONS.md D-2)",
                            "CBB18 and CBB19 (CPP5 sequence contested; two "
                            "candidates each are in 03_enumerated.csv; "
                            "DECISIONS.md D-1)"],
       no_stated_library_size="The paper states no library size and runs no "
                              "Markush enumeration: every member is named, "
                              "purified and individually characterised, so "
                              "there is no block product to reconcile against.",
       enumerated_vs_stated="34 named compounds enumerated; 34 named in the "
                            "SI. Panels carry values for 22 + 12 + 22 + 8 "
                            "compound-slots over four assay blocks."),
     "P1_hrms_formula_and_mz": p1_hrms(enum_rows),
     "P2_molar_mass_from_the_authors_own_weighings": p2_weighings(smiles_by_name),
     "P3_paper_own_hit_criterion_and_hit_list": p3_hits(vals, controls),
     "P4_printed_mean_and_sd": p4_printed(vals, grid6),
     "P5_flow_versus_microscopy": p5_cross_panel(vals),
     "P6_column_to_compound_mapping": p6_mapping(vals, controls),
     "P7_replicate_count_vs_legend": p7_replicates(vals),
     "P8_ring_methyl_regiochemistry": dict(
       what="which ring position carries the toluidine methyl",
       si_phrasings={
         "CBB21 and CBB1 (divalent)": "bis(3-methyl-4,1-phenylene)",
         "CBBZ1-CBBZ10 (divalent)": "bis(2-methyl-4,1-phenylene)",
         "CBB6/7/12/13/14/20 (monovalent, twice in each name)":
            "(4-(...amino)-2-methylphenyl) attached to the central methine, "
            "AND ...-3-methylphenyl)(ethyl)amino cited on the nitrogen"},
       resolution="The two divalent forms contradict each other under any one "
                  "locant convention. The monovalent prefixes do not: in a "
                  "monovalent prefix position 1 IS the attachment, so both of "
                  "them put the methyl ORTHO to the central methine and META "
                  "to the nitrogen. The synthesis agrees -- both series come "
                  "from I2, 3-((ethyl(m-tolyl)amino)methyl)benzenesulfonic "
                  "acid, whose N-alkyl-m-toluidine is attacked para to the "
                  "nitrogen (C4), ortho to the C3 methyl.",
       emitted="methyl ortho to the central methine in BOTH series",
       NOT_CHECKABLE_BY_MASS="every candidate is an isomer; see "
                             "08_contradictions.csv row 1"),
     "P10_drawn_markush_table_vs_iupac_names": p10_drawn_table(smiles_by_name),
     "P11_drawn_cbbz_table_vs_iupac_names": p11_drawn_cbbz_table(),
     "P9_ratio_condition": dict(
       emitted_ratio="[cargo]/[CBB(Z)] = 1:8",
       corroboration="Supplementary Figure 4 is the ratio titration (2, 4, 8, "
                     "16) on CBB2/7/11/12/13/14/18/19/Z1 and contains an 8 "
                     "row, and the Methods say 'decided to maintain the "
                     "CBB:RNP ratio of 8:1'. So the emitted panels are the "
                     "8:1 condition, read from a second place.",
       not_emitted="the titration itself -- it varies dose, not chemistry"),
    }
    qc["seeded_defects"] = seeded(enum_rows, vals, controls, smiles_by_name,
                                  grid6)

    # ------------------------------------------------------------ the gates
    p1, p2, p3, p4 = (qc["P1_hrms_formula_and_mz"],
                      qc["P2_molar_mass_from_the_authors_own_weighings"],
                      qc["P3_paper_own_hit_criterion_and_hit_list"],
                      qc["P4_printed_mean_and_sd"])
    gates = {
      "P1 our formula reproduces the SI's printed m/z for all but CBB7":
        p1["C_our_formula_matches_si_mz"] >= p1["n_with_hrms"] - 1,
      "P2 every weighing lands inside the molar-mass interval its own "
      "printed precision allows":
        p2["n_inside_printed_interval"] == p2["n_weighings"],
      "P3 both hit sets reproduce EXACTLY, no missing and no extra":
        all(v["exact_set_match"] for v in p3.values()),
      "P4 all three printed mean+/-SD statements reproduce":
        p4["n_reproduced"] == p4["n_statements"],
      "P4 all six Fig 6a statements reproduce on population SD":
        p4["fig6a_n_reproduced_pop_sd"] == p4["fig6a_n_statements"],
      "P6 every panel's column->compound map is injective":
        all(v["injective"] for k, v in qc["P6_column_to_compound_mapping"].items()
            if isinstance(v, dict)),
      "P6 Fig2c and Fig4c each cover CBB1..CBB22":
        qc["P6_column_to_compound_mapping"]["Fig2c_covers_CBB1_to_CBB22"]
        and qc["P6_column_to_compound_mapping"]["Fig4c_covers_CBB1_to_CBB22"],
      "P7 every panel's replicate count equals its legend's":
        qc["P7_replicate_count_vs_legend"]["n_agreeing"] == 4,
      "P10 all 21 drawn Markush cells agree with the structures built from "
      "the IUPAC names":
        qc["P10_drawn_markush_table_vs_iupac_names"]["n_agreeing"]
        == qc["P10_drawn_markush_table_vs_iupac_names"]["n_cells_drawn"],
      "P11 every CBBZ drawn cell either agrees or has a declared, evidenced "
      "mismatch":
        qc["P11_drawn_cbbz_table_vs_iupac_names"]["n_unresolved"] == 0,
      "every seeded defect moved its own metric":
        all(s["metric_moved"] for s in qc["seeded_defects"]),
    }
    qc["gates"] = gates
    qc["all_gates_pass"] = all(gates.values())
    json.dump(qc, open(os.path.join(_HERE, "07_qc.json"), "w"), indent=1)

    for k, v in gates.items():
        print(("  PASS  " if v else "  FAIL  ") + k)
    print("seeded defects: %d/%d moved their metric"
          % (sum(s["metric_moved"] for s in qc["seeded_defects"]),
             len(qc["seeded_defects"])))
    return 0 if qc["all_gates_pass"] else 1


if __name__ == "__main__":
    sys.exit(main())
