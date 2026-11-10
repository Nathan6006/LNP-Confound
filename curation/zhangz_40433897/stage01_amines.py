#!/usr/bin/env python3
"""Stage 01 - the 35 drawn amines of Figure S1a, and what the Figure S2 HRMS
masses can and cannot pin.

Primary sources, both from source/nn5c01378_si_001.pdf:
  * Figure S1a (p4, embedded JPEG xref 72, 220 dpi): 35 drawn amines, each
    labelled with a CAS registry number. The CAS NUMBER is the identifier used
    here - the drawing is used only to confirm it, never to trace a structure.
  * Figure S2 caption (p5-6): 8 HRMS calc/found pairs, read out of the PDF TEXT
    layer with subscript x-positions resolved against the baseline glyphs.

Nothing here writes a corpus row. It answers one question: given the drawn
reaction (Fig S1b: R-NH2 + 2 x VE-acrylate -> tertiary amine with two
tocopheryl 3-oxopropanoate arms), which of the 35 amines are consistent with
each measured exact mass, and is that consistent set a singleton?
"""
import json
import os
import sys
from collections import Counter, defaultdict

from rdkit import Chem, RDLogger
from rdkit.Chem.Descriptors import ExactMolWt
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

RDLogger.DisableLog("rdApp.*")
HERE = os.path.dirname(os.path.abspath(__file__))

# --- Figure S1a, read left-to-right, top-to-bottom at 220 dpi -----------------
# (cas, smiles, name, drawn_row, drawn_col)
# `smiles` is written from the CAS registry number; `drawn_agrees` records
# whether the structure drawn beside that CAS matches it, read by eye at 3x
# upscale. A disagreement would be a finding, not a silent correction.
AMINES = [
    ("7209-38-3",   "NCCCN1CCN(CCCN)CC1",         "1,4-bis(3-aminopropyl)piperazine", 1, 1, True),
    ("109-55-7",    "CN(C)CCCN",                   "N,N-dimethyl-1,3-propanediamine", 1, 2, True),
    ("109-85-3",    "COCCN",                       "2-methoxyethylamine", 1, 3, True),
    ("5036-48-6",   "NCCCn1ccnc1",                 "1-(3-aminopropyl)imidazole", 1, 4, True),
    ("100-36-7",    "CCN(CC)CCN",                  "N,N-diethylethylenediamine", 1, 5, True),
    ("4439-20-7",   "OCCNCCNCCO",                  "N,N'-bis(2-hydroxyethyl)ethylenediamine", 1, 6, True),

    ("10563-29-8",  "CN(C)CCCNCCCN",               "N-(3-aminopropyl)-N',N'-dimethyl-1,3-propanediamine", 2, 1, True),
    ("65604-89-9",  "CNCCN(CCNC)CCNC",             "tris(2-(methylamino)ethyl)amine", 2, 2, True),
    ("7154-73-6",   "NCCN1CCCC1",                  "1-(2-aminoethyl)pyrrolidine", 2, 3, True),
    ("51387-90-7",  "CN1CCCC1CCN",                 "2-(1-methylpyrrolidin-2-yl)ethylamine", 2, 4, True),
    ("123-00-2",    "NCCCN1CCOCC1",                "4-(3-aminopropyl)morpholine", 2, 5, True),
    ("3197-06-6",   "OCCN(CCO)CCN",                "N,N-bis(2-hydroxyethyl)ethylenediamine", 2, 6, True),

    ("24426-16-2",  "CCN(CC)CCNCCN",               "N-(2-aminoethyl)-N',N'-diethylethylenediamine", 3, 1, True),
    ("105-83-9",    "NCCCN(C)CCCN",                "bis(3-aminopropyl)methylamine", 3, 2, True),
    ("152829-22-6", "NCCN(C)CCCN(C)CCN",           "N,N'-bis(2-aminoethyl)-N,N'-dimethyl-1,3-propanediamine", 3, 3, True),
    ("115-69-5",    "OCC(C)(N)CO",                 "2-amino-2-methyl-1,3-propanediol", 3, 4, True),
    ("140-31-8",    "NCCN1CCNCC1",                 "1-(2-aminoethyl)piperazine", 3, 5, True),
    ("694-83-7",    "N[C@H]1CCCC[C@@H]1N",         "trans-1,2-diaminocyclohexane", 3, 6, True),

    ("104-78-9",    "CCN(CC)CCCN",                 "N,N-diethyl-1,3-propanediamine", 4, 1, True),
    ("107-15-3",    "NCCN",                        "ethylenediamine", 4, 2, True),
    ("57734-57-3",  "NCC1CCCN1",                   "2-(aminomethyl)pyrrolidine", 4, 3, False),
    ("4067-16-7",   "NCCNCCNCCNCCNCCN",            "pentaethylenehexamine", 4, 4, True),
    ("4246-51-9",   "NCCCOCCOCCOCCCN",             "4,7,10-trioxa-1,13-tridecanediamine", 4, 5, True),

    ("929-06-6",    "NCCOCCO",                     "2-(2-aminoethoxy)ethanol", 5, 1, True),
    ("502-32-9",    "CC(C)C[C@@H](N)CO",           "L-leucinol", 5, 2, True),
    ("31295-54-2",  "NCCN1CCN(CCNCCN)CC1",         "1-(2-aminoethyl)-4-(2-((2-aminoethyl)amino)ethyl)piperazine", 5, 3, True),
    ("1477-55-0",   "NCc1cccc(CN)c1",              "m-xylylenediamine", 5, 4, True),
    ("4097-89-6",   "NCCN(CCN)CCN",                "tris(2-aminoethyl)amine", 5, 5, True),

    ("7209-38-3",   "NCCCN1CCN(CCCN)CC1",          "1,4-bis(3-aminopropyl)piperazine (DRAWN TWICE)", 6, 1, True),
    ("120-20-7",    "COc1ccc(CCN)cc1OC",           "3,4-dimethoxyphenethylamine", 6, 2, True),
    ("7356-00-5",   "CC(N)CNCC(C)N",               "bis(2-aminopropyl)amine", 6, 3, True),
    ("111-41-1",    "NCCNCCO",                     "N-(2-aminoethyl)ethanolamine", 6, 4, True),
    ("10061-68-4",  "CCNCCCNCC",                   "N,N'-diethyl-1,3-propanediamine", 6, 5, True),

    ("2620-27-1",   "CNCCOCCNC",                   "2,2'-oxybis(N-methylethylamine)", 7, 1, True),
    ("2508-29-4",   "NCCCCCO",                     "5-amino-1-pentanol", 7, 2, True),
]

# --- Figure S2 caption, from the PDF text layer ------------------------------
# label_as_printed is EXACTLY what the caption prints; subscripts were resolved
# by matching subscript glyph x0 to the preceding baseline glyph's x-extent.
HRMS = [
    dict(panel="a", label_as_printed="VE-AC",   formula="C32H52O3Na",  adduct="[M+Na]+",
         calc=507.3809, found=507.3804),
    dict(panel="b", label_as_printed="VN3C2O",  formula="C71H123N2O8", adduct="[M+H]+",
         calc=1131.9274, found=1131.9349),
    dict(panel="c", label_as_printed="VN2C4",   formula="C70H119N2O6", adduct="[M+H]+",
         calc=1083.9063, found=1083.9053),
    dict(panel="d", label_as_printed="VN2C2",   formula="C70H121N2O6", adduct="[M+H]+",
         calc=1085.9220, found=1085.9316),
    dict(panel="e", label_as_printed="VN3C2",   formula="C71H123N2O6", adduct="[M+H]+",
         calc=1099.9376, found=1099.9362),
    dict(panel="f", label_as_printed="VN3C1",   formula="C69H119N2O6", adduct="[M+H]+",
         calc=1071.9063, found=1071.9062),
    dict(panel="g", label_as_printed="VN3C2",   formula="C71H121N2O6", adduct="[M+H]+",
         calc=1097.9220, found=1097.9215),
    dict(panel="h", label_as_printed="VN2C1",   formula="C68H117N2O6", adduct="[M+H]+",
         calc=1057.8907, found=1057.8902),
]

VE_AC = "CC1=C(C)C2=C(CCC(C)(CCC[C@H](C)CCC[C@H](C)CCCC(C)C)O2)C(C)=C1OC(=O)C=C"
ARM_ADDED = Counter(C=32, H=52, O=3)   # Michael addition: nothing is lost


def formula_counts(f):
    """Parse a Hill-ish formula string into element counts."""
    import re
    out = Counter()
    for el, n in re.findall(r"([A-Z][a-z]?)(\d*)", f):
        if el:
            out[el] += int(n) if n else 1
    return out


def exact_mass_of_formula(f):
    """Monoisotopic mass of a neutral formula string."""
    m = {"C": 12.0, "H": 1.00782503207, "N": 14.0030740048,
         "O": 15.9949146196, "Na": 22.9897692809}
    c = formula_counts(f)
    return sum(m[el] * n for el, n in c.items())


ELECTRON = 0.00054857990946


def main():
    report = {"inputs": {}, "amines": [], "hrms": [], "core_matches": {},
              "findings": []}

    # ---------- 0. the VE-AC intermediate, as an external check on the scheme
    ve = Chem.MolFromSmiles(VE_AC)
    assert ve is not None, "VE-AC SMILES does not parse"
    ve_formula = CalcMolFormula(ve)
    ve_na = ExactMolWt(ve) + 22.9897692809 - ELECTRON
    report["inputs"]["VE_AC"] = {
        "smiles": VE_AC, "formula": ve_formula,
        "calc_M_plus_Na": round(ve_na, 4),
        "paper_calc": 507.3809, "paper_found": 507.3804,
        "delta_mDa_vs_paper_calc": round((ve_na - 507.3809) * 1000, 2),
    }

    # ---------- 1. the 35 drawn amines
    by_formula = defaultdict(list)
    seen_cas = Counter()
    for cas, smi, name, row, col, drawn_ok in AMINES:
        m = Chem.MolFromSmiles(smi)
        if m is None:
            raise SystemExit("unparseable amine SMILES for CAS %s" % cas)
        f = CalcMolFormula(m)
        # how many N-H bonds does each nitrogen carry -> how many Michael arms
        # can sit on it, and is the 2-arm site unique?
        nh = []
        for a in m.GetAtoms():
            if a.GetSymbol() == "N":
                nh.append(a.GetTotalNumHs())
        rec = dict(cas=cas, smiles=smi, name=name, formula=f,
                   exact_mass=round(ExactMolWt(m), 4),
                   n_nitrogens=len(nh), nh_per_nitrogen=sorted(nh, reverse=True),
                   drawing_agrees_with_cas=drawn_ok,
                   drawn_at="row %d col %d" % (row, col))
        report["amines"].append(rec)
        by_formula[f].append(rec)
        seen_cas[cas] += 1

    dupes = {c: n for c, n in seen_cas.items() if n > 1}
    report["inputs"]["n_drawn_entries"] = len(AMINES)
    report["inputs"]["n_distinct_cas"] = len(seen_cas)
    report["inputs"]["duplicate_cas"] = dupes
    if dupes:
        report["findings"].append(
            "Figure S1a caption says '35 common commercial amino small "
            "molecules'. 35 entries ARE drawn, but CAS %s is drawn twice "
            "(row 1 col 1 and row 6 col 1), so only %d are DISTINCT."
            % (", ".join(dupes), len(seen_cas)))

    # ---------- 2. what core each HRMS formula implies
    for h in HRMS:
        if h["panel"] == "a":
            report["hrms"].append(dict(h, core_formula=None, note="intermediate"))
            continue
        ionf = formula_counts(h["formula"])
        # recomputed calc for the stated ion
        calc_recomputed = exact_mass_of_formula(h["formula"]) - ELECTRON
        neutral = Counter(ionf)
        neutral["H"] -= 1                      # [M+H]+ -> M
        core = Counter(neutral)
        for el, n in ARM_ADDED.items():        # strip two VE-AC arms
            core[el] -= 2 * n
        core_str = "".join("%s%d" % (el, core[el]) for el in ("C", "H", "N", "O")
                           if core[el])
        bad = {el: n for el, n in core.items() if n < 0}
        rec = dict(h,
                   calc_recomputed=round(calc_recomputed, 4),
                   calc_minus_paper_mDa=round((calc_recomputed - h["calc"]) * 1000, 2),
                   found_minus_calc_mDa=round((h["found"] - h["calc"]) * 1000, 2),
                   core_formula=core_str,
                   core_impossible=bool(bad))
        report["hrms"].append(rec)

        cands = by_formula.get(core_str, [])
        # a candidate is CONNECTIVITY-DETERMINATE only if exactly one nitrogen
        # can take two arms, i.e. exactly one N has >= 2 H and no other N has any
        det = []
        for c in cands:
            nh = c["nh_per_nitrogen"]
            two_arm_sites = sum(1 for x in nh if x >= 2)
            any_other_nh = sum(1 for x in nh if x >= 1) - two_arm_sites
            det.append(dict(cas=c["cas"], name=c["name"], smiles=c["smiles"],
                            two_arm_sites=two_arm_sites,
                            other_NH_sites=any_other_nh,
                            connectivity_determinate=(two_arm_sites == 1 and
                                                      any_other_nh == 0)))
        report["core_matches"][h["panel"] + ":" + h["label_as_printed"]] = dict(
            core_formula=core_str, n_candidates_in_FigS1a=len(cands),
            candidates=det,
            pinned=(len(det) == 1 and det[0]["connectivity_determinate"]))

    # ---------- 3. the duplicate label
    labels = Counter(h["label_as_printed"] for h in HRMS if h["panel"] != "a")
    for lab, n in labels.items():
        if n > 1:
            fs = [h["formula"] for h in HRMS if h["label_as_printed"] == lab]
            report["findings"].append(
                "Figure S2 prints the label '%s' on %d different panels, with "
                "DIFFERENT formulas %s. They differ by H2, so at most one of "
                "the two labels can be right." % (lab, n, fs))

    # ---------- 4. the CH2 ladder, measured not asserted
    def mass_of(lab):
        for h in HRMS:
            if h["panel"] != "a" and h["label_as_printed"] == lab:
                yield h["panel"], h["calc"]
    ladder = []
    CH2 = 14.015650
    pairs = [("VN2C1", "h", "VN3C1", "f"), ("VN2C2", "d", "VN3C2", "e"),
             ("VN2C4", "c", "VN3C2", "g")]
    for lo, lop, hi, hip in pairs:
        a = [h for h in HRMS if h["panel"] == lop][0]["calc"]
        b = [h for h in HRMS if h["panel"] == hip][0]["calc"]
        ladder.append(dict(low=lo + " (panel %s)" % lop, high=hi + " (panel %s)" % hip,
                           delta=round(b - a, 4),
                           minus_CH2_mDa=round((b - a - CH2) * 1000, 2)))
    report["ch2_ladder"] = ladder

    out = os.path.join(HERE, "01_amines.json")
    with open(out, "w") as fh:
        json.dump(report, fh, indent=1)
    print(json.dumps({k: report[k] for k in
                      ("inputs", "core_matches", "ch2_ladder", "findings")}, indent=1))
    print("\nwrote", out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
