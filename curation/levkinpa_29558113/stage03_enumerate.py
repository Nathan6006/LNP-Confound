#!/usr/bin/env python3
"""Stage 03 -- blocks and the enumerated library. Honest about 288.

MECHANISM, TESTED AGAINST THIS PAPER'S OWN SCHEME, NOT INHERITED.
The sibling curation curation/lij_42130331 (Guo 2026, same lab, 8 years later)
independently established that a thiolactone ring opens AT THE CARBONYL -- amide
to the incoming amine, free thiol released. That was carried into this paper as
a HYPOTHESIS and it is confirmed here by the paper's own drawings, which is the
only reason it is used:

  * Scheme S1c (SI p.2, vector) draws T18U + PY16 + A1 -> T18U-PY16-A1 as
    Et2N-CH2CH2-NH-C(=O)-CH(NH-C(=O)-C17H33)-CH2CH2-S-S-C16H33, plus
    2-pyridothione as the by-product. Carbonyl attack, disulfide exchange.
  * Figure S1 (SI p.15) draws the same transformation twice more, once with a
    primary amine (A5) and once with a secondary amine (A10), and its panels
    a)-f) measure the released 2-pyridothione at 370 nm -- the by-product that
    only exists if the thiol displaces pyridine from the disulfide.

So this paper's own scheme, not the sibling's finding, is what the enumerator
implements. There is NO epoxide in this chemistry (the sibling's second half
does not transfer).

ASSEMBLY uses curation/lib/assemble.attach: no atom is deleted anywhere, so the
F-013 mechanism (RWMol.RemoveAtom silently inverting an adjacent stereocentre)
cannot occur, and every build runs check_stereo_preserved and check_no_radicals.
Attachment atoms are tagged with atom MAPS, never "index 0 of the fragment".

STEREOCHEMISTRY. The core is built from DL-homocysteine thiolactone (SI p.2,
stated). The alpha carbon is therefore a RACEMIC centre and is emitted with no
stereo descriptor -- that is the correct statement, not an omission. The oleoyl
9Z is drawn and measured (stage01) and IS specified.

Run:  python3 -I stage03_enumerate.py
"""
import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))

from rdkit import Chem                                      # noqa: E402
from rdkit.Chem import Descriptors, rdMolDescriptors        # noqa: E402

import assemble                                             # noqa: E402
from errors import StructureError                           # noqa: E402

# The ring-opened homocysteine core, with three mapped attachment atoms.
# Each is written as a VALID atom carrying the hydrogen it will lose, so the
# fragment parses as a real molecule and assemble.attach never has to repair a
# valence (writing '[N:1]' instead would make a radical -- the Siegwart
# radical-carbon class, which parses and canonicalises and is invisible).
#   :1  the new amide carbonyl  -> the amine's nitrogen
#   :2  the outer disulfide S   -> the alkyl arm's C1
#   :3  the N-acyl carbonyl     -> the tail's C2
CORE = "O=[CH:1]C(CCS[SH:2])N[CH:3]=O"

THIOLACTONE = {            # block -> (acyl carbons, tail fragment = C2..Cn)
    "T14": (14, "[CH3:6]" + "C" * 12),
    "T16": (16, "[CH3:6]" + "C" * 14),
    "T18": (18, "[CH3:6]" + "C" * 16),
    # oleoyl: 17 fragment carbons, cis double bond between fragment atoms 8
    # and 9, i.e. C9=C10 of the acyl. Both the length and the POSITION are
    # measured in stage01 from the drawing's own bond coordinates.
    "T18U": (18, "[CH3:6]" + "C" * 6 + "/C=C\\" + "C" * 8),
}
DISULFIDE = {              # block -> (arm carbons, arm fragment = C1..Cn)
    "PY12": (12, "[CH3:5]" + "C" * 11),
    "PY14": (14, "[CH3:5]" + "C" * 13),
    "PY16": (16, "[CH3:5]" + "C" * 15),
    "PY18": (18, "[CH3:5]" + "C" * 17),
}
AMINES = {                 # block -> SMILES with the attacking N mapped :4
    "A1": "[NH2:4]CCN(CC)CC",
    "A5": "[NH2:4]CCCN(CC)CC",
}
N_AMINES_IN_LIBRARY = 18


def build(t, p, a):
    """Thiolactone x disulfide x amine -> one canonical lipidoid SMILES."""
    m = assemble.attach(CORE, DISULFIDE[p][1], a_map=2, b_map=5,
                        label=f"{t}-{p}-{a}:S-S-alkyl")
    m = assemble.attach(m, THIOLACTONE[t][1], a_map=3, b_map=6,
                        label=f"{t}-{p}-{a}:tail")
    m = assemble.attach(m, AMINES[a], a_map=1, b_map=4,
                        label=f"{t}-{p}-{a}:amide")
    mol = Chem.MolFromSmiles(m)
    if mol is None:
        raise StructureError(f"{t}-{p}-{a}: product does not parse")
    # The oleoyl 9Z must SURVIVE assembly. check_stereo_preserved compares
    # tetrahedral CIP codes only, so a lost double-bond configuration would be
    # invisible to it; this is the assertion that covers that gap.
    if t == "T18U":
        nst = sum(1 for b in mol.GetBonds()
                  if b.GetStereo() != Chem.BondStereo.STEREONONE)
        if nst != 1:
            raise StructureError(
                f"{t}-{p}-{a}: {nst} stereo-defined double bonds in the "
                "product, expected exactly 1 (the oleoyl 9Z). Assembly lost "
                "or duplicated the configuration.")
    return Chem.MolToSmiles(mol), mol


def main():
    rows_blocks, rows_enum = [], []
    amines_json = json.load(open(os.path.join(HERE, "01_amines.json")))
    chain_json = json.load(open(os.path.join(HERE, "02_chain_counts.json")))

    for b, (n, frag) in THIOLACTONE.items():
        cj = chain_json["blocks"][b]
        rows_blocks.append({
            "block_id": b, "role": "thiolactone (amide tail donor)",
            "n_carbons": n,
            # BUILT FROM THE SAME FRAGMENT THE ENUMERATOR USES, with the
            # atom map stripped. Writing the block's SMILES as a second
            # literal is how a block and the product built from it drift
            # apart: the HRMS check below would then pass on a correct block
            # while every product carried a wrong tail. Found by seeding
            # exactly that mutation (10_seeded.json, M4).
            "SMILES": Chem.MolToSmiles(Chem.MolFromSmiles(
                "O=C1SCCC1NC(=O)" + frag.replace("[CH3:6]", "C"))),
            "measured_by": "vector bond coordinates + HRMS carbon count",
            "geometry_carbons": cj["carbons_from_geometry"],
            "HRMS_carbons": cj["carbons_from_HRMS_formula"],
            "agree": cj["geometry_and_mass_agree"],
            "evidence": f"SI p.2 Scheme S1a; HRMS [M+H]+ "
                        f"{cj['HRMS_mz_printed']}",
            "confidence": "high", "enumerated": True})

    for b, (n, frag) in DISULFIDE.items():
        cj = chain_json["blocks"][b]
        rows_blocks.append({
            "block_id": b, "role": "pyridyl disulfide (S-S-alkyl arm donor)",
            "n_carbons": n,
            # 2-(alkyldisulfanyl)pyridine: the S-S must sit on the carbon
            # ADJACENT to the ring nitrogen. Writing it as "c1ccncc1SS..."
            # puts it meta instead -- same formula, same exact mass, invisible
            # to every mass check. Hence the explicit ring closure on n1.
            "SMILES": Chem.MolToSmiles(Chem.MolFromSmiles(
                frag.replace("[CH3:5]", "C") + "SSc1ccccn1")),
            "measured_by": "vector bond coordinates + HRMS carbon count",
            "geometry_carbons": cj["carbons_from_geometry"],
            "HRMS_carbons": cj["carbons_from_HRMS_formula"],
            "agree": cj["geometry_and_mass_agree"],
            "evidence": f"SI p.2 Scheme S1b; HRMS [M+H]+ "
                        f"{cj['HRMS_mz_printed']}",
            "confidence": "high", "enumerated": True})

    for a in ("A1", "A5", "A8", "A10"):
        aj = amines_json[a]
        rows_blocks.append({
            "block_id": a, "role": "amine head group",
            "n_carbons": "", "SMILES": aj.get("SMILES") or "",
            "measured_by": aj.get("method", ""),
            "geometry_carbons": "", "HRMS_carbons": "",
            "agree": "",
            "evidence": aj.get("name_printed") or aj.get("why")
            or aj.get("why_not_enumerated", ""),
            "confidence": aj.get("confidence") or "",
            "enumerated": bool(aj.get("enumerated"))})
    for i in range(1, N_AMINES_IN_LIBRARY + 1):
        a = f"A{i}"
        if a in ("A1", "A5", "A8", "A10"):
            continue
        rows_blocks.append({
            "block_id": a, "role": "amine head group",
            "n_carbons": "", "SMILES": "", "measured_by": "",
            "geometry_carbons": "", "HRMS_carbons": "", "agree": "",
            "evidence": "NOT IN THE SI. The amine panel is in the main "
                        "article, which this curation does not hold.",
            "confidence": "", "enumerated": False})

    # ---- the enumerated library: every one of the 4 x 4 x 18 labels ----
    for t in ("T14", "T16", "T18", "T18U"):
        for p in ("PY12", "PY14", "PY16", "PY18"):
            for i in range(1, N_AMINES_IN_LIBRARY + 1):
                a = f"A{i}"
                name = f"{t}-{p}-{a}"
                if a in AMINES:
                    smi, mol = build(t, p, a)
                    rows_enum.append({
                        "lipidoid": name, "thiolactone": t, "disulfide": p,
                        "amine": a, "IL_SMILES": smi,
                        "formula": rdMolDescriptors.CalcMolFormula(mol),
                        "exact_mass": round(Descriptors.ExactMolWt(mol), 4),
                        "MolWt": round(Descriptors.MolWt(mol), 3),
                        "structure_determined": True,
                        "confidence": amines_json[a]["confidence"],
                        "why": "" if a == "A1" else
                               "the two N-ethyl branches of A5 are read from a "
                               "drawing, not measured (01_amines.json)"})
                else:
                    rows_enum.append({
                        "lipidoid": name, "thiolactone": t, "disulfide": p,
                        "amine": a, "IL_SMILES": "", "formula": "",
                        "exact_mass": "", "MolWt": "",
                        "structure_determined": False, "confidence": "",
                        "why": ("molecular formula C7H17N3 only, no structure"
                                if a == "A8" else
                                "ring size measurement inconclusive"
                                if a == "A10" else
                                "amine not defined anywhere in the SI")})

    with open(os.path.join(HERE, "01_blocks.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows_blocks[0]))
        w.writeheader()
        w.writerows(rows_blocks)
    with open(os.path.join(HERE, "03_enumerated.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows_enum[0]))
        w.writeheader()
        w.writerows(rows_enum)

    # ---- HRMS: every block structure against its own printed [M+H]+ ----
    H = 1.00782503207
    E = 0.00054857990
    mass_checks = []
    for r in rows_blocks:
        b = r["block_id"]
        if b not in chain_json["blocks"] or not r["SMILES"]:
            continue
        m = Chem.MolFromSmiles(r["SMILES"])
        calc = Descriptors.ExactMolWt(m) + H - E
        printed = chain_json["blocks"][b]["HRMS_mz_printed"]
        mass_checks.append({
            "block": b, "formula": rdMolDescriptors.CalcMolFormula(m),
            "our_MplusH": round(calc, 4), "printed_calcd": printed,
            "delta_mDa": round(1000 * (calc - printed), 2)})
    worst = max(abs(c["delta_mDa"]) for c in mass_checks)
    if worst > 3.0:
        raise SystemExit(
            f"stage03: a block structure is {worst:.2f} mDa from the SI's own "
            "printed calcd [M+H]+. That is a composition error, not rounding. "
            "Refusing to emit.")
    with open(os.path.join(HERE, "03_block_hrms.json"), "w") as fh:
        json.dump({"_what": "Each enumerated BLOCK structure against the "
                            "SI's own printed calcd [M+H]+. Mass cannot see "
                            "connectivity -- this tests composition only.",
                   "electron_mass_subtracted": True,
                   "worst_abs_delta_mDa": round(worst, 2),
                   "checks": mass_checks}, fh, indent=1)

    det = sum(1 for r in rows_enum if r["structure_determined"])
    recon = {
        "library_size_STATED_by_the_paper": 288,
        "stated_where": "SI p.7 Scheme S2 and its caption: 'There are 16 "
                        "series, each series contains 18 lipidoids and in "
                        "total 288 lipidoids are present in the library.' "
                        "4 thiolactones x 4 disulfides = 16; 16 x 18 = 288.",
        "label_grid_ENUMERATED_here": len(rows_enum),
        "arithmetic_matches_the_claim": len(rows_enum) == 288,
        "structures_DETERMINED": det,
        "structures_not_determined": len(rows_enum) - det,
        "why": "4 x 4 = 16 backbone combinations are fully determined, but "
               f"only {len(AMINES)} of the 18 amines are. 16 x "
               f"{len(AMINES)} = {det}.",
        "unique_canonical_SMILES": len({r["IL_SMILES"] for r in rows_enum
                                        if r["IL_SMILES"]}),
        "duplicate_structures": det - len({r["IL_SMILES"] for r in rows_enum
                                           if r["IL_SMILES"]}),
        "block_HRMS_worst_abs_delta_mDa": round(worst, 2),
    }
    with open(os.path.join(HERE, "03_reconciliation.json"), "w") as fh:
        json.dump(recon, fh, indent=1)
    print(json.dumps(recon, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
