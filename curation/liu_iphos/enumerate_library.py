#!/usr/bin/env python3
"""Enumerate the 572 iPhos lipids from the confirmed amines and the nAxPm row labels.

Chemistry (all three parts independently confirmed - see 01_scheme.json):
  Pm  = 2-alkoxy-1,3,2-dioxaphospholane 2-oxide, m carbons in the O-alkyl (m = 4..16)
  nA  = amine, structures read from Fig 1d by vector bond counting
  rxn = ring-opening at an amine nitrogen:  N  ->  N-CH2CH2-O-P(=O)(OH)-O-CmH(2m+1)

`x` in nAxPm is NOT inferred - it is printed in the Fig 2a row labels, and the number of
nitrogens in each amine exactly predicts which x rows exist (19A has 2 N and only 19A2
appears; 28A has 5 N and 28A2-28A5 appear). That reproduces all 44 rows, so the
substitution rule is confirmed structurally as well as by HRMS.

ONE REAL AMBIGUITY, flagged rather than hidden: when x < (number of nitrogens), the paper
does not say WHICH nitrogens carry the phosphate - the real product is a positional
mixture. A convention is required to emit a single SMILES. We substitute primary N first,
then secondary, then tertiary (nucleophilicity/steric order), breaking ties by RDKit atom
index. Rows where x == nN are unambiguous. The `positional_ambiguity` column marks the
rest, so a mismatch against another curation on those rows is a convention difference,
not necessarily an error.
"""
import csv
import json
import re

from rdkit import Chem, RDLogger
from rdkit.Chem.Descriptors import ExactMolWt, MolWt
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

RDLogger.DisableLog("rdApp.*")

# --- amines, all confirmed in 01_scheme.json ---------------------------------
# 1A-6A   n-alkylamines              C4,6,8,10,12,14
# 7A-13A  di-n-alkylamines           C4,6,8,10,11,12,13   (irregular step confirmed twice)
# 14A-18A tri-n-alkylamines          C4,6,8,10,12
# 19A-28A polyamines, read from Fig 1d
def nalkyl(c):
    return "C" * c


AMINES = {}
for i, c in enumerate([4, 6, 8, 10, 12, 14], start=1):
    AMINES[f"{i}A"] = (f"N{nalkyl(c)}", f"n-C{c} amine")
for i, c in zip(range(7, 14), [4, 6, 8, 10, 11, 12, 13]):
    AMINES[f"{i}A"] = (f"{nalkyl(c)}N{nalkyl(c)}", f"di-n-C{c} amine")
for i, c in zip(range(14, 19), [4, 6, 8, 10, 12]):
    AMINES[f"{i}A"] = (f"{nalkyl(c)}N({nalkyl(c)}){nalkyl(c)}", f"tri-n-C{c} amine")
AMINES.update({
    "19A": ("NCCOCCOCCN", "2,2'-(ethylenedioxy)bis(ethylamine)"),
    "20A": ("NCCN1CCNCC1", "1-(2-aminoethyl)piperazine"),
    "21A": ("NCCCNCCCN(C)C", "N,N-dimethyldipropylenetriamine"),
    "22A": ("NCCCN(C)CCCN", "3,3'-diamino-N-methyldipropylamine"),
    "23A": ("NCCCN1CCN(CCCN)CC1", "1,4-bis(3-aminopropyl)piperazine"),
    "24A": ("NCCNCCNCCN", "triethylenetetramine"),
    "25A": ("NCCCNCCCCNCCCN", "spermine"),
    "26A": ("CNCCN(CCNC)CCNC", "tris(2-(methylamino)ethyl)amine"),
    "27A": ("NCCN(CCN)CCN", "tris(2-aminoethyl)amine"),
    "28A": ("NCCNCCNCCNCCN", "tetraethylenepentamine"),
})


def phosphate_arm(m, anionic=False):
    """-CH2CH2-O-P(=O)(OH)-O-alkyl, m carbons in the O-alkyl.

    `anionic` gives the deprotonated phosphodiester, used when the nitrogen being
    substituted is already tertiary: ring-opening then makes it QUATERNARY, so the
    product is a permanent zwitterion (the classic phosphocholine motif) rather than a
    protonatable amine. Charges are made explicit so RDKit's valence model accepts it
    and so downstream descriptors see the real species.
    """
    return "CCOP(=O)([O-])O" + "C" * m if anionic else "CCOP(=O)(O)O" + "C" * m


def n_rank(atom):
    """Substitution order: primary N first, then secondary, then tertiary."""
    heavy = sum(1 for nb in atom.GetNeighbors() if nb.GetAtomicNum() > 1)
    return (heavy, atom.GetIdx())


def nh_slots(mol):
    """Available N-H substitution sites per nitrogen.

    Ring-opening of Pm needs an N-H, so a PRIMARY amine offers two sites and a secondary
    one. A tertiary nitrogen offers none and can only react by quaternising, which the
    paper does use - but only where there is no N-H anywhere in the amine (14A-18A are
    trialkylamines, and those products are the quaternary zwitterions).

    This replaces an earlier, wrong model that allowed at most one arm per nitrogen and
    so quaternised tertiary nitrogens while primary ones still had a free H. That model
    reproduced every primary/secondary/tertiary amine correctly but failed the
    polyamines; the signatures of the unmatched molecules (3 arms / 0 quaternary N / 3 N,
    and so on) are only reachable if NH2 can take two arms.
    """
    slots = []
    for a in mol.GetAtoms():
        if a.GetSymbol() != "N":
            continue
        heavy = sum(1 for nb in a.GetNeighbors() if nb.GetAtomicNum() > 1)
        slots.append((a.GetIdx(), max(0, 3 - heavy)))
    return slots


def build(amine_smiles, x, m):
    mol = Chem.MolFromSmiles(amine_smiles)
    if mol is None:
        return None, None
    ns = [a for a in mol.GetAtoms() if a.GetSymbol() == "N"]
    nN = len(ns)
    slots = nh_slots(mol)
    # spend arms on N-H sites first, one per nitrogen in primary-before-secondary order,
    # then take second H's of primary nitrogens, and only quaternise if arms remain
    order = sorted(slots, key=lambda s: (-s[1], s[0]))
    plan = {i: 0 for i, _ in slots}
    left = x
    for rnd in (1, 2):
        for idx, cap in order:
            if left <= 0:
                break
            if plan[idx] < min(cap, rnd):
                plan[idx] += 1
                left -= 1
    quat = set()
    if left > 0:  # no N-H capacity remains: quaternise tertiary nitrogens
        for idx, cap in sorted(slots, key=lambda s: s[0]):
            if left <= 0:
                break
            if cap == 0:
                plan[idx] += 1
                quat.add(idx)
                left -= 1
    if left > 0:
        return None, None
    rw = Chem.RWMol(mol)
    for idx, k in sorted(plan.items()):
        for _ in range(k):
            arm = Chem.MolFromSmiles(phosphate_arm(m, anionic=(idx in quat)))
            off = rw.GetNumAtoms()
            rw.InsertMol(arm)
            rw.AddBond(idx, off, Chem.BondType.SINGLE)
        if idx in quat:
            rw.GetAtomWithIdx(idx).SetFormalCharge(1)
            rw.GetAtomWithIdx(idx).SetNoImplicit(True)
    out = rw.GetMol()
    try:
        Chem.SanitizeMol(out)
    except Exception:
        return None, None
    return Chem.MolToSmiles(out), nN


def main():
    vals = list(csv.DictReader(open("04_values.csv")))
    rows, fails = [], []
    for v in vals:
        am = v["amine"]
        x = int(v["x"])
        m = int(v["pm"][1:])
        smi, nN = build(AMINES[am][0], x, m)
        if smi is None:
            fails.append(v["lipid_name"])
            continue
        mol = Chem.MolFromSmiles(smi)
        rows.append({
            "lipid_name": v["lipid_name"], "amine": am, "x": x, "pm": v["pm"],
            "IL_SMILES": smi,
            "amine_desc": AMINES[am][1],
            "n_nitrogens": nN,
            "positional_ambiguity": int(x < nN),
            "MolWt": round(MolWt(mol), 3),
            "ExactMolWt": round(ExactMolWt(mol), 4),
            "formula": CalcMolFormula(mol),
            "bin_index": v["bin_index"], "bin_label": v["bin_label"],
            "rlu_midpoint": v["rlu_midpoint"],
            "extraction_method": v["extraction_method"],
            "source_pdf": v["source_pdf"], "page": v["page"], "figure": v["figure"],
        })
    with open("03_enumerated.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"enumerated {len(rows)} lipids ({len(fails)} failed: {fails[:5]})")
    uniq = {r["IL_SMILES"] for r in rows}
    print(f"unique canonical SMILES: {len(uniq)}")
    amb = sum(r["positional_ambiguity"] for r in rows)
    print(f"rows with positional ambiguity (x < nN): {amb}")

    # C1: HRMS anchors
    print("\nC1 HRMS anchors:")
    anchors = {"9A1P9": 491.4, "10A1P10": 561.5, "9A1P15": 575.5, "10A1P16": 645.6}
    idx = {r["lipid_name"]: r for r in rows}
    for k, rep in anchors.items():
        r = idx.get(k)
        if r:
            d = r["ExactMolWt"] - rep
            print(f"   {k:9} {r['ExactMolWt']:9.3f} vs {rep}  delta {d:+.3f}  "
                  f"{'PASS' if abs(d) < 0.15 else 'FAIL'}   {r['formula']}")


if __name__ == "__main__":
    main()
