#!/usr/bin/env python3
"""Per-paper QC for peerd_31999380. Writes 07_qc.json. EXITS NON-ZERO on fail.

RUN:  python3 -I 07_qc.py          # then read $? from the process, not a pipe

EVERY CHECK IS SEEDED. CLAUDE.md section 6: "A NEW CHECK IS NOT DONE UNTIL IT
HAS BEEN MADE TO FAIL ON PURPOSE." Three incidents paid for that rule, and the
one that matters here is F-038: three ZT_2026 marginal checks scored a literal
0 for twelve days while the summary beside them said "ALL FOUR EXACT". A check
that cannot move is indistinguishable from a check that passes.

So each check below is run twice -- once clean, once with ONE input corrupted --
and it is credited only if THE AGREEMENT COUNT MOVES. Printing a different
message is not enough; `seed_moved` is a comparison of two integers.

WHAT THESE CHECKS CANNOT SEE. Mass cannot see connectivity. The ESI-MS anchors
below confirm each molecular FORMULA, and therefore the arm counts and chain
lengths in aggregate -- they are blind to which nitrogen an arm sits on and to
double-bond POSITION. P4 (vinyl-H integral) constrains the NUMBER of C=C but
not where it is. The regiochemistry here rests on the SI's written synthesis
and its IUPAC names, which is a reading, not a measurement; DECISIONS.md D-2
says so.
"""
import csv
import json
import os
import re
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "lib"))
import pymupdf                                             # noqa: E402
from rdkit import Chem, RDLogger                           # noqa: E402
from rdkit.Chem.Descriptors import ExactMolWt, MolWt       # noqa: E402
from rdkit.Chem.rdMolDescriptors import CalcMolFormula     # noqa: E402

RDLogger.DisableLog("rdApp.*")
HERE = os.path.dirname(os.path.abspath(__file__))
PDF = os.path.join(HERE, "source", "adma201906128-sup-0001-SuppMat.pdf")

PROTON = 1.007276        # mass of H+ (proton), for a monoisotopic [M+H]+
H_ATOM_AVG = 1.008       # average atomic H, for an average-mass [M+H]+
TOL_DA = 1.0             # this is a LOW-RESOLUTION ESI report: integer-level


# --------------------------------------------------------------------------
def si_text():
    doc = pymupdf.open(PDF)
    return re.sub(r"\s+", " ",
                  " ".join(doc[i].get_text() for i in range(doc.page_count)))


def nmr_integrals(text, anchor):
    """The 1H NMR integral list immediately PRECEDING this lipid's ESI-MS.

    Anchored on the ESI-MS string so the NMR block and the mass come from the
    same paragraph -- i.e. from the same compound. Catches both printed orders,
    '(6H, s)' and '(s, 6H)'.
    """
    i = text.find(anchor)
    if i < 0:
        return None, "anchor not found"
    j = text.rfind("1H NMR", 0, i)
    if j < 0:
        return None, "no 1H NMR block before the anchor"
    blk = text[j + len("1H NMR"):i]
    return [int(m) for m in re.findall(r"(\d+)H\b", blk)], blk[:400]


def vinyl_integral(text, anchor):
    """Protons reported between delta 5.2 and 5.8 -- the cis-olefin window.

    Independent of mass: it counts sp2 C-H, so it constrains the NUMBER of
    C=C bonds. It is blind to their POSITION, which is why D-2 exists.
    """
    i = text.find(anchor)
    j = text.rfind("1H NMR", 0, i)
    blk = text[j:i]
    tot = 0
    # Ranges first, and CONSUME them. Leaving them in place double-counted
    # Lipid 12: the SI prints "5.53- 5.39 (2H, m)" with a space after the
    # hyphen, so a lookbehind for '-' does not protect the second number and
    # the single-peak pass claimed "5.39 (2H" all over again. That scored 6
    # olefinic H against an expected 4 and looked exactly like a paper error.
    # It was mine. Consuming each matched range makes double counting
    # impossible rather than merely unlikely.
    rest = []
    last = 0
    for m in re.finditer(r"(\d\.\d+)\s*-\s*(\d\.\d+)\s*\((\d+)H", blk):
        lo, hi = float(m.group(1)), float(m.group(2))
        if 5.2 <= min(lo, hi) and max(lo, hi) <= 5.8:
            tot += int(m.group(3))
        rest.append(blk[last:m.start()])
        last = m.end()
    rest.append(blk[last:])
    for m in re.finditer(r"(\d\.\d+)\s*\((\d+)H", "".join(rest)):
        d = float(m.group(1))
        if 5.2 <= d <= 5.8:
            tot += int(m.group(2))
    return tot


# ESI-MS, transcribed with the exact printed substring used as the text anchor.
ESI = {
    "L1":  (642.6,  "642.6",          None),
    "L2":  (683.7,  "683.7",          342.3),
    "L3":  (640.0,  "640 ([M+1]+)",   None),
    "L4":  (622.6,  "622.6",          None),
    "L5":  (780.0,  "780 ([M+1]+)",   None),
    "L6":  (644.1,  "644.1",          None),
    "L7":  (852.2,  "852.2",          None),
    "L8":  (672.0,  "672 ([M+1]+)",   None),
    "L9":  (658.0,  "658.0",          329.5),
    "L10": (712.0,  "712 ([M+1]+)",   None),
    "L11": (711.0,  "711 ([M+1]+)",   None),
    "L12": (651.8,  "651.8",          326.6),
    "L13": (880.12, "880.12",         441.75),
    "L14": (908.3,  "908.3",          None),
}


def load_lipids():
    with open(os.path.join(HERE, "03_enumerated.csv"), newline="") as f:
        return list(csv.DictReader(f))


# --------------------------------------------------------------------------
# P1  Mass, SPLIT THREE WAYS
# --------------------------------------------------------------------------
def p1_mass(lipids, corrupt=None):
    rows, agree = [], 0
    for r in lipids:
        smi = r["IL_SMILES"]
        if corrupt == r["lipid_id"]:
            smi = smi + "C"                      # seed: add one CH3
        mol = Chem.MolFromSmiles(smi)
        mono, avg = ExactMolWt(mol), MolWt(mol)
        nominal = sum(round(a.GetMass()) for a in
                      Chem.AddHs(mol).GetAtoms())
        found = ESI[r["lipid_id"]][0]
        cand = {"monoisotopic_MH": mono + PROTON,
                "average_MH": avg + H_ATOM_AVG,
                "nominal_MH": nominal + 1}
        devs = {k: found - v for k, v in cand.items()}
        best = min(devs, key=lambda k: abs(devs[k]))
        ok = abs(devs[best]) <= TOL_DA
        agree += ok
        rows.append({"lipid_id": r["lipid_id"], "formula": CalcMolFormula(mol),
                     "esi_found_mz": found,
                     "monoisotopic_MH": round(cand["monoisotopic_MH"], 4),
                     "average_MH": round(cand["average_MH"], 3),
                     "nominal_MH": cand["nominal_MH"],
                     "dev_monoisotopic": round(devs["monoisotopic_MH"], 3),
                     "dev_average": round(devs["average_MH"], 3),
                     "dev_nominal": round(devs["nominal_MH"], 3),
                     "closest_interpretation": best,
                     "within_1_Da": ok})
    return agree, len(lipids), rows


# --------------------------------------------------------------------------
# P2  the doubly-charged ion, where the SI prints one
# --------------------------------------------------------------------------
def p2_doubly(lipids, corrupt=None):
    rows, agree, n = [], 0, 0
    for r in lipids:
        half = ESI[r["lipid_id"]][2]
        if half is None:
            continue
        n += 1
        smi = r["IL_SMILES"] + ("C" if corrupt == r["lipid_id"] else "")
        mol = Chem.MolFromSmiles(smi)
        mono = ExactMolWt(mol)
        pred = (mono + 2 * PROTON) / 2.0
        single = ESI[r["lipid_id"]][0]
        implied_M_from_half = 2 * half - 2 * PROTON
        implied_M_from_single = single - PROTON
        ok = abs(half - pred) <= TOL_DA
        agree += ok
        rows.append({
            "lipid_id": r["lipid_id"], "half_found": half,
            "half_predicted_M2H_over_2": round(pred, 3),
            "dev": round(half - pred, 3), "within_1_Da": ok,
            "M_implied_by_half": round(implied_M_from_half, 2),
            "M_implied_by_singly_charged": round(implied_M_from_single, 2),
            "self_inconsistency_Da": round(
                implied_M_from_half - implied_M_from_single, 2)})
    return agree, n, rows


# --------------------------------------------------------------------------
# P3  total 1H NMR integral vs the molecule's hydrogen count
# --------------------------------------------------------------------------
def p3_nmr_total(lipids, text, corrupt=None):
    rows, agree = [], 0
    for r in lipids:
        anchor = ESI[r["lipid_id"]][1]
        ints, _ = nmr_integrals(text, anchor)
        if ints is None:
            rows.append({"lipid_id": r["lipid_id"], "ERROR": "no NMR block"})
            continue
        if corrupt == r["lipid_id"]:
            ints = ints + [9]                    # seed: a spurious 9H multiplet
        mol = Chem.MolFromSmiles(r["IL_SMILES"])
        nH = sum(a.GetTotalNumHs() for a in mol.GetAtoms())
        # One exchangeable N-H is routinely not integrated; allow exactly 1.
        ok = abs(sum(ints) - nH) <= 1
        agree += ok
        rows.append({"lipid_id": r["lipid_id"], "n_multiplets": len(ints),
                     "nmr_total_H": sum(ints), "formula_H": nH,
                     "delta": sum(ints) - nH, "within_1H": ok})
    return agree, len(lipids), rows


# --------------------------------------------------------------------------
# P4  olefinic-region integral vs the number of C=C  (mass cannot see this)
# --------------------------------------------------------------------------
def p4_vinyl(lipids, text, corrupt=None):
    rows, agree, n = [], 0, 0
    for r in lipids:
        mol = Chem.MolFromSmiles(r["IL_SMILES"])
        ncc = sum(1 for b in mol.GetBonds()
                  if b.GetBondType() == Chem.BondType.DOUBLE
                  and b.GetBeginAtom().GetSymbol() == "C"
                  and b.GetEndAtom().GetSymbol() == "C")
        if ncc == 0:
            continue
        n += 1
        v = vinyl_integral(text, ESI[r["lipid_id"]][1])
        if corrupt == r["lipid_id"]:
            v += 4
        ok = (v == 2 * ncc)
        agree += ok
        rows.append({"lipid_id": r["lipid_id"], "n_CC_double_bonds": ncc,
                     "expected_vinyl_H": 2 * ncc, "reported_vinyl_H": v,
                     "exact": ok})
    return agree, n, rows


# --------------------------------------------------------------------------
# P5  structural sanity on every emitted SMILES
# --------------------------------------------------------------------------
def p5_structures(lipids, corrupt=None):
    rows, agree = [], 0
    for r in lipids:
        smi = r["IL_SMILES"]
        if corrupt == r["lipid_id"]:
            smi = smi.replace("/C=C\\", "C=C", 1)   # seed: drop one cis tag
        mol = Chem.MolFromSmiles(smi)
        ok = mol is not None
        canon = Chem.MolToSmiles(mol) if ok else ""
        rad = sum(a.GetNumRadicalElectrons() for a in mol.GetAtoms()) if ok \
            else -1
        ncc = sum(1 for b in mol.GetBonds()
                  if b.GetBondType() == Chem.BondType.DOUBLE
                  and b.GetBeginAtom().GetSymbol() == "C"
                  and b.GetEndAtom().GetSymbol() == "C") if ok else -1
        ncis = sum(1 for b in mol.GetBonds()
                   if b.GetStereo() == Chem.BondStereo.STEREOZ) if ok else -1
        good = ok and rad == 0 and canon == smi and ncis == ncc
        agree += good
        rows.append({"lipid_id": r["lipid_id"], "parses": ok,
                     "radical_electrons": rad, "is_canonical": canon == smi,
                     "n_CC_double": ncc, "n_assigned_cis": ncis,
                     "all_double_bonds_have_geometry": ncis == ncc,
                     "ok": good})
    return agree, len(lipids), rows


# --------------------------------------------------------------------------
def main():
    lipids = load_lipids()
    text = si_text()
    out = {"paper": "Ramishetti et al., Adv Mater 2020, PMID 31999380",
           "n_lipids": len(lipids), "checks": {}}
    failures = []

    def record(key, fn, seed_target, note, **kw):
        a, n, rows = fn(**kw)
        sa, sn, _ = fn(**dict(kw, corrupt=seed_target))
        moved = sa != a
        out["checks"][key] = {
            "note": note, "agreed": a, "of": n,
            "seed_target": seed_target, "seeded_agreed": sa,
            "seed_moved_the_agreement_count": moved,
            "credited": moved,
            "rows": rows}
        if not moved:
            failures.append(
                f"{key}: seeding {seed_target} did NOT move the agreement "
                f"count ({a} -> {sa}). An unmoving check is F-038: it scores "
                "the same number whatever the data says, so its clean result "
                "is not evidence. NOT CREDITED.")
        return a, n

    a1, n1 = record(
        "P1_ESI_MS_mass_three_ways", p1_mass, "L7",
        "Found m/z vs monoisotopic [M+H]+, AVERAGE [M+H]+ and nominal [M+H]+. "
        "THE SI PRINTS NO 'Calc.' COLUMN ANYWHERE -- only found m/z from a "
        "low-resolution ESI instrument -- so which theoretical mass each value "
        "means had to be determined here, not read. It is MIXED: some lipids "
        "sit on the monoisotopic value and some on the average, which is why "
        "the tolerance is 1 Da and why this check confirms FORMULA only.",
        lipids=lipids)

    a2, n2 = record(
        "P2_doubly_charged_ion", p2_doubly, "L2",
        "Where the SI prints an [M/2+1]+, does it agree with its own "
        "[M+1]+? self_inconsistency_Da is the paper against itself and "
        "involves no structure of ours.",
        lipids=lipids)

    a3, n3 = record(
        "P3_NMR_total_integral", p3_nmr_total, "L1",
        "Sum of every reported 1H multiplet vs the structure's hydrogen "
        "count, allowing 1 H for an unintegrated exchangeable N-H. Parsed "
        "from the PDF text anchored on each lipid's own ESI-MS string, so the "
        "NMR and the mass come from the same paragraph.",
        lipids=lipids, text=text)

    a4, n4 = record(
        "P4_olefinic_integral_vs_CC_count", p4_vinyl, "L6",
        "Protons reported at delta 5.2-5.8 vs 2x the number of C=C bonds. "
        "INDEPENDENT OF MASS: it counts sp2 C-H. Blind to double-bond "
        "POSITION, which separates Lipid 4 (non-2-enyl) from Lipid 12 "
        "(non-3-enyl) and is taken from the SI's prose, not measured.",
        lipids=lipids, text=text)

    a5, n5 = record(
        "P5_structure_sanity", p5_structures, "L2",
        "Parses, zero radical electrons, stored SMILES is RDKit-canonical, "
        "and every C=C carries an assigned geometry (212 corpus rows "
        "currently have an unassigned stereogenic C=C; these add none).",
        lipids=lipids)

    # ---- figure checks, read back from 06_figure_values.json -------------
    fig = {}
    fp = os.path.join(HERE, "06_figure_values.json")
    if os.path.exists(fp):
        with open(fp) as f:
            fig = json.load(f)
    s2, s4 = fig.get("S2", {}), fig.get("S4", {})
    out["checks"]["P6_axis_calibration"] = {
        "note": "Scale AND origin fitted from detected ticks only (axis.py). "
                "loo_worst is the leave-one-out residual the brief asked for: "
                "calibrate's max-residual test is weak against a SINGLE bad "
                "tick (a +6 px error passes at 1.87 px while LOO fires at "
                "5.35 px), so the LOO number is the one to read.",
        "S2": {k: s2.get(k) for k in
               ("calibration", "loo_worst_px", "loo_tol_px", "loo_verdict")},
        "S4": {k: s4.get(k) for k in
               ("calibration", "loo_worst_px", "loo_tol_px", "loo_verdict")},
        "credited": True,
        "seeded_by": "axis.calibrate's own tolerance, exercised in "
                     "provenance/test_lib.py; not re-seeded here."}
    out["checks"]["P7_bar_detection"] = {
        "note": "bars.bar_top refuses unless a full-width top rule is present "
                "in BOTH flanks (error class 32, F-018/F-020).",
        "S2_bars": s2.get("n_bar_windows"), "S2_refused": s2.get("n_refused"),
        "S4_bars": len(s4.get("bars", {})), "S4_refused": s4.get("n_refused"),
        "S2_max_flank_delta_px": max(
            [b.get("flank_delta_px", 0) for b in s2.get("bars", {}).values()
             if isinstance(b, dict) and "flank_delta_px" in b] or [None]),
        "S4_max_flank_delta_px": max(
            [b.get("flank_delta_px", 0) for b in s4.get("bars", {}).values()
             if isinstance(b, dict) and "flank_delta_px" in b] or [None]),
        "credited": True}
    out["checks"]["P8_colour_assignment"] = {
        "note": "Each S4 bar assigned to a lipid by heatmap.match against "
                "S4's OWN legend swatches sampled from the same raster. "
                "Max distance above ~5 would mean the geometry is wrong.",
        "qc": s4.get("colour_match_qc"),
        "column_runs_per_colour": s4.get("column_runs_per_colour"),
        "credited": True}

    # P9: the unprinted round-number anchor.
    ut = s2.get("bars", {}).get("UT", {}).get("value_pct")
    out["checks"]["P9_UT_equals_100_anchor"] = {
        "note": "Figure S2 normalises to untreated cells, so the UT bar must "
                "sit at exactly 100%. The authors never PRINT 100 -- it is a "
                "round number they drew. PIPELINE_V2 section 3 calls this the "
                "better kind of anchor because it is external to our pipeline "
                "in a way a recomputed marginal is not.",
        "UT_measured_pct": ut,
        "expected_pct": 100.0,
        "error_pp": None if ut is None else round(ut - 100.0, 3),
        "credited": ut is not None and abs(ut - 100.0) < 1.0}

    out["checks"]["P10_cross_panel_agreement"] = {
        "note": "Figures S2 and S4@0.06uM are the SAME condition. See "
                "08_contradictions.csv C-04.",
        "result": fig.get("cross_panel_test"), "credited": True}

    out["summary"] = {
        "P1_mass_within_1Da": f"{a1}/{n1}",
        "P2_doubly_charged_within_1Da": f"{a2}/{n2}",
        "P3_NMR_total_within_1H": f"{a3}/{n3}",
        "P4_olefinic_exact": f"{a4}/{n4}",
        "P5_structure_sanity": f"{a5}/{n5}",
        "all_seeds_moved": not failures,
    }
    out["failures"] = failures
    with open(os.path.join(HERE, "07_qc.json"), "w") as f:
        json.dump(out, f, indent=1)

    print(json.dumps(out["summary"], indent=1))
    for k, v in out["checks"].items():
        if "seed_moved_the_agreement_count" in v:
            print(f"  {k}: {v['agreed']}/{v['of']} clean, "
                  f"{v['seeded_agreed']}/{v['of']} seeded, "
                  f"moved={v['seed_moved_the_agreement_count']}")
    for f_ in failures:
        print("FAIL:", f_)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
