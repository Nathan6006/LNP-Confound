#!/usr/bin/env python3
"""Enumerate the lipidoids this SI actually pins, through curation/lib.

CHEMISTRY. Lipidoid = amine + 2.4 eq acrylate, aza-Michael, 80 C, 2 days
(SI p.2). Amine 1's N-H adds across the acrylate C=C, giving a
3-(alkoxy)-3-oxopropyl arm on the amine nitrogen.

ROUTE. Built as: mapped 3-propanoic acid --esterify--> the arm ester, then
attach() twice onto the amine nitrogen. No RemoveAtom near a stereocentre
(there are no stereocentres here at all), no string surgery on canonical
SMILES -- every join is an atom-mapped attach or lib.esterify.

WHAT IS PINNED AND WHAT IS NOT.
  amine 1  pinned twice over: the 1H NMR pattern (N(CH3)2 6H s; three
           2H multiplets of one propylene; no other skeleton protons) and the
           carbon balance of BOTH printed ESI-MS formulae.
  O16      pinned: C43H86N2O4 and the 56H bulk-CH2 integral both give hexadecyl.
  O16B     pinned: C39H78N2O4S4 and the 40H bulk-CH2 integral both give
           2-(dodecyldisulfanyl)ethyl -- NOT a C16 alkyl. The name token is
           counting CHAIN HEAVY ATOMS (2 C + 2 S + 12 C = 16), which is also
           why the paper can call O16B and O16 "similar chemical structures".
  O14B,
  O18B     NAMED in the SI (p.2) but given no formula, no NMR and no structure.
           Their alkyl is extrapolated from the single pinned B-member by that
           chain-atom convention: 14 -> decyl, 18 -> tetradecyl.
           structure_assumed = 1, confidence = medium. See DECISIONS.md D-2.
  O14, O18 NEVER NAMED ANYWHERE IN THIS SI. Not enumerated. Writing them would
           be inventing library members, not curating them.
  amine 2  "amine 1 or 2" (SI p.2) -- amine 2's structure appears nowhere in
           the SI. Not enumerated. See DECISIONS.md D-1.
"""
import csv
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "lib"))
import assemble  # noqa: E402
from errors import StructureError  # noqa: E402
from rdkit import Chem  # noqa: E402
from rdkit.Chem import rdMolDescriptors, Descriptors  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))

# --- source string literals, read off SI p.2-p.3 -----------------------------
# amine 1 = N,N-dimethyl-1,3-propanediamine. The primary N carries both arms.
AMINE1 = "CN(C)CCC[NH2:1]"

# the mapped Michael-donor carbon of the acrylate, in its ADDUCT form:
# the arm that ends up on nitrogen is -CH2CH2C(=O)O-.
ARM_ACID = "[CH2:2]CC(=O)O"

# alcohols. n-alkanols for the plain acrylates; 2-(alkyldisulfanyl)ethanols
# for the bioreducible ones. Written out, not generated from a template, so a
# reviewer counts the carbons on the page.
ALKANOL = {16: "OCCCCCCCCCCCCCCCC"}                      # hexadecan-1-ol, 16 C
DISULFIDE_OL = {                                          # HO-CH2CH2-S-S-CnH2n+1
    14: "OCCSSCCCCCCCCCC",      # 2-(decyldisulfanyl)ethanol      2+2+10 = 14
    16: "OCCSSCCCCCCCCCCCC",    # 2-(dodecyldisulfanyl)ethanol    2+2+12 = 16
    18: "OCCSSCCCCCCCCCCCCCC",  # 2-(tetradecyldisulfanyl)ethanol 2+2+14 = 18
}

ACRYLATES = {
    "O14B": DISULFIDE_OL[14],
    "O16B": DISULFIDE_OL[16],
    "O18B": DISULFIDE_OL[18],
    "O16": ALKANOL[16],
}

# what the SI PRINTS, for the formula anchor. Absent where the SI prints none.
PRINTED_FORMULA = {"1-O16B": "C39H78N2O4S4", "1-O16": "C43H86N2O4"}
PRINTED_CALCD = {"1-O16B": 767.3, "1-O16": 695.2}
PRINTED_FOUND = {"1-O16B": 767.7, "1-O16": 695.7}
# the bulk -CH2- envelope integral, SI p.3
PRINTED_BULK_CH2_H = {"1-O16B": 40, "1-O16": 56}
# every printed integral, in print order, for the closure test
PRINTED_INTEGRALS = {
    "1-O16B": [6, 40, 2, 6, 2, 6, 4, 4, 4, 4],
    "1-O16": [6, 56, 2, 6, 2, 6, 4, 4],
}


def _remap_amine_nitrogen(smiles, mapno=1):
    """Put the attachment map back on the ONE secondary aliphatic N.

    Done with a SMARTS match on the molecule, never by editing the canonical
    SMILES string: a sibling built acyl fragments by string surgery on
    canonical SMILES and silently mislocated every ester, with identical
    formula and exact mass to 8 dp. Refuses unless exactly one N matches.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise StructureError(f"remap: unparseable {smiles!r}")
    patt = Chem.MolFromSmarts("[NX3;H1;!$([N]C=O)]")
    hits = mol.GetSubstructMatches(patt)
    if len(hits) != 1:
        raise StructureError(
            f"remap: {len(hits)} secondary amine N in {smiles!r}, need 1. "
            "Which nitrogen takes the second arm would be a guess.")
    mol.GetAtomWithIdx(hits[0][0]).SetAtomMapNum(mapno)
    return Chem.MolToSmiles(mol)


def build(label, alcohol_smiles):
    """attach, esterify, attach, esterify -- in that order, and the order is
    forced by curation/lib, not chosen for convenience.

    The obvious route (esterify the mapped propanoic acid first, then attach
    the finished arm twice) IS REFUSED: `[CH2:2]` is valence-incomplete, so
    sanitizing the fragment assigns it an unpaired electron, and `esterify`
    does not clear that because the mapped carbon is not one of ITS
    attachment atoms -- only `attach` clears the atoms it bonds.
    check_no_radicals fired on the first run with
    "1 radical atom(s): [(0, 'C', 1)]". The refusal is correct and was
    respected rather than relaxed: doing the aza-Michael FIRST means the
    mapped carbon is always an `attach` attachment atom, so RDKit fills its
    hydrogens by valence. Each `esterify` then sees exactly one free COOH,
    which is the precondition it enforces.
    """
    m1 = assemble.attach(AMINE1, ARM_ACID, a_map=1, b_map=2,
                         label=f"{label}:aza-Michael-1")
    m2 = assemble.esterify(m1, alcohol_smiles, label=f"{label}:ester-1")
    m3 = _remap_amine_nitrogen(m2, 1)
    m4 = assemble.attach(m3, ARM_ACID, a_map=1, b_map=2,
                         label=f"{label}:aza-Michael-2")
    m5 = assemble.esterify(m4, alcohol_smiles, label=f"{label}:ester-2")
    return m5


def main():
    rows = []
    for acr, alcohol in ACRYLATES.items():
        name = f"1-{acr}"
        smi = build(name, alcohol)
        mol = Chem.MolFromSmiles(smi)
        formula = rdMolDescriptors.CalcMolFormula(mol)
        pinned = name in PRINTED_FORMULA
        rows.append({
            "IL_name": name,
            "IL_SMILES": smi,
            "amine": "amine 1 = N,N-dimethyl-1,3-propanediamine",
            "acrylate": acr,
            "acrylate_alcohol_SMILES": alcohol,
            "formula_built": formula,
            "formula_printed": PRINTED_FORMULA.get(name, ""),
            "exact_mass": round(Descriptors.ExactMolWt(mol), 4),
            "average_mass": round(Descriptors.MolWt(mol), 3),
            "n_stereocentres": len(assemble.stereocentres(mol)),
            "evidence": ("1H NMR + ESI-MS formula, SI p.3" if pinned else
                         "NAME ONLY (SI p.2). Alkyl extrapolated from O16B by "
                         "the chain-heavy-atom convention; no formula, NMR or "
                         "structure is printed for this acrylate."),
            "confidence": "high" if pinned else "medium",
            "structure_assumed": 0 if pinned else 1,
        })

    with open(os.path.join(HERE, "03_enumerated.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(json.dumps(rows, indent=1))
    return rows


if __name__ == "__main__":
    main()
