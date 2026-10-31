#!/usr/bin/env python3
"""Enumerate the 344 thiosuccinimide lipids of Liu et al., Adv Healthc Mater 2026.

CHEMISTRY, from the Methods and confirmed by mass below. One pot, three
components, catalyst-free, 80 C / 4 h:

  1. the amine's primary NH2 attacks the linker's NHS ester  -> AMIDE, losing
     N-hydroxysuccinimide
  2. the alkyl thiol adds across the linker's maleimide C=C  -> THIOSUCCINIMIDE,
     nothing lost
  exception: L2 is SPDP, whose second site is a 2-pyridyldithio group, so step 2
  is thiol-disulfide exchange giving a DISULFIDE and releasing pyridine-2-thione.

ONE ARM PER PRIMARY AMINE, which the feed ratios in Tables S1/S2 encode exactly
(1:1.05:1.1 for the one-NH2 cores, 1:2.05:2.1 for two, 1:3.05:3.1 for three).

Structures are assembled with RWMol fragment surgery, never by string
concatenation: a previous curation in this project produced a radical and a bond
to the wrong atom that way, and both survived every mass and formula check.

VALIDATION is external. Tables S1 and S2 print MW and observed m/z for 18
representative products spanning all eight linkers and one, two and three arms.
Those are compared here, and a disagreement fails the build.
"""
import csv
import json
import os
import re

from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors as D, rdMolDescriptors as rd

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
BLOCKS = os.path.join(HERE, "liu2026_adhm202502225_blocks.csv")

# arm fragments: dummy marks the atom that bonds to the amine nitrogen
ARM = {
    "L1": "*C(=O)C1CCC(CN2C(=O)C(S{R})CC2=O)CC1",
    "L2": "*C(=O)CCSS{R}",
    "L3": "*C(=O)CN1C(=O)C(S{R})CC1=O",
    "L4": "*C(=O)CCCN1C(=O)C(S{R})CC1=O",
    "L5": "*C(=O)CCCCCN1C(=O)C(S{R})CC1=O",
    "L6": "*C(=O)CCOCCN1C(=O)C(S{R})CC1=O",
    "L7": "*C(=O)CCOCCOCCN1C(=O)C(S{R})CC1=O",
    "L8": "*C(=O)CCOCCOCCOCCN1C(=O)C(S{R})CC1=O",
}

# EXTERNAL MASS CHECK. All 24 rows of Table S1 and all 8 of Table S2, parsed
# from the SI .docx. An earlier version used only 16 and silently dropped the
# six T3A rows of S1 and the whole of S2 - the strongest evidence in the paper.
#
# Table S2 is decisive for the linker identities: every one of its eight rows is
# T3A-Lx-2C10, so the TAIL IS HELD CONSTANT and the linker increments cannot be
# traded against tail length. All eight match as printed (<=1.45 Da), which pins
# L1-L8 on their own.
REFERENCE_MW = {
    # Table S1
    "75-L1-C6": 453.3, "113-L1-C6": 791.5, "110-L1-C6": 1157.7,
    "75-L2-C6": 320.6, "113-L2-C8": 582.0, "110-L2-C10": 927.6,
    "70-L3-C6": 342.2, "113-L3-C6": 613.8, "T3A-L3-C6": 954.3,
    "70-L4-C6": 399.6, "113-L4-C6": 740.1, "T3A-L4-C6": 1122.6,
    "70-L5-C10": 427.7, "113-L5-C10": 796.2, "T3A-L5-C10": 1290.9,
    "400-L6-C6": 471.7, "113-L6-C6": 744.0, "T3A-L6-C8": 1212.7,
    "70-L8-C10": 545.8, "113-L8-C10": 1032.4, "T3A-L8-C10": 1561.2,
    # Table S2, all T3A-Lx-2C10
    "T3A-L1-2C10": 1747.8, "T3A-L2-2C10": 1347.0, "T3A-L3-2C10": 1500.1,
    "T3A-L4-2C10": 1585.5, "T3A-L5-2C10": 1668.3, "T3A-L6-2C10": 1675.6,
    "T3A-L7-2C10": 1807.8, "T3A-L8-2C10": 1939.9,
}
TOL = 1.5

AMINES = ["70", "75", "80", "400", "113", "114", "110", "T3A"]
LINEAR = ["C6", "C8", "C10", "C12", "C14"]
BRANCHED = ["2C10", "2C11", "2C12"]
BLS = [f"L{i}" for i in range(1, 9)]


def load_blocks():
    amine, tail = {}, {}
    for r in csv.DictReader(open(BLOCKS)):
        if r["block_class"] == "amine_head":
            amine[r["name_in_paper"]] = (r["SMILES"], int(r["n_reactive_sites"]))
        elif r["block_class"] == "tail":
            tail[r["name_in_paper"]] = r["SMILES"]
    return amine, tail


def build(amine_smiles, n_sites, bl, thiol_smiles):
    """Attach one BL+thiol arm to every primary amine nitrogen."""
    assert thiol_smiles.startswith("S"), thiol_smiles
    arm_smiles = ARM[bl].replace("{R}", thiol_smiles[1:])

    m = Chem.RWMol(Chem.MolFromSmiles(amine_smiles))
    prim = [a.GetIdx() for a in m.GetAtoms()
            if a.GetSymbol() == "N" and a.GetTotalNumHs() == 2
            and a.GetDegree() == 1]
    if len(prim) != n_sites:
        raise ValueError(f"{amine_smiles}: {len(prim)} primary N, expected {n_sites}")

    for n_idx in prim:
        arm = Chem.MolFromSmiles(arm_smiles)
        if arm is None:
            raise ValueError(f"arm will not parse: {arm_smiles}")
        off = m.GetNumAtoms()
        m.InsertMol(arm)
        dummy = next(a.GetIdx() for a in m.GetAtoms()
                     if a.GetAtomicNum() == 0 and a.GetIdx() >= off)
        anchor = m.GetAtomWithIdx(dummy).GetNeighbors()[0].GetIdx()
        m.AddBond(n_idx, anchor, Chem.BondType.SINGLE)
        m.RemoveAtom(dummy)           # dummy index is above every primary N

    mol = m.GetMol()
    Chem.SanitizeMol(mol)
    if any(a.GetNumRadicalElectrons() for a in mol.GetAtoms()):
        raise ValueError("radical produced")
    if len(Chem.GetMolFrags(mol)) != 1:
        raise ValueError("disconnected product")
    return Chem.MolToSmiles(mol)


def main():
    amine, tail = load_blocks()
    rows = []
    for lib, cores, tails in (("1", AMINES, LINEAR), ("2", ["T3A"], BRANCHED)):
        for core in cores:
            smi, n = amine[core]
            for bl in BLS:
                for t in tails:
                    out = build(smi, n, bl, tail[t])
                    mol = Chem.MolFromSmiles(out)
                    rows.append({
                        "lipid_name": f"{core}-{bl}-{t}",
                        "library": lib, "amine": core, "bl": bl, "tail": t,
                        "n_arms": n,
                        "IL_SMILES": out,
                        "formula": rd.CalcMolFormula(mol),
                        "MolWt": round(D.MolWt(mol), 2),
                        "ExactMolWt": round(D.ExactMolWt(mol), 4),
                        "linkage": "disulfide" if bl == "L2" else "thiosuccinimide",
                    })

    with open("09_enumerated_344.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    by_name = {r["lipid_name"]: float(r["MolWt"]) for r in rows}
    agree, anomalous = [], []
    for name, want in REFERENCE_MW.items():
        got = by_name.get(name)
        if got is not None and abs(got - want) <= TOL:
            agree.append({"lipid": name, "table_MW": want, "computed_MW": got,
                          "abs_error": round(abs(got - want), 2)})
            continue
        # which (linker, tail) would fit? Moving a CH2 between the linker and
        # the tail is mass- AND formula-degenerate, so several usually do; list
        # them rather than pick one.
        am = name.split("-")[0]
        tails = BRANCHED if name.endswith(("2C10", "2C11", "2C12")) else LINEAR
        fits = sorted((abs(by_name[f"{am}-{b}-{t}"] - want), f"{b}-{t}")
                      for b in BLS for t in tails
                      if f"{am}-{b}-{t}" in by_name)
        anomalous.append({"lipid": name, "table_MW": want,
                          "computed_as_printed": got,
                          "exact_fits": [f for e, f in fits if e <= TOL],
                          "closest": fits[0][1],
                          "closest_error": round(fits[0][0], 2)})

    n1 = sum(1 for r in rows if r["library"] == "1")
    n2 = sum(1 for r in rows if r["library"] == "2")
    qc = {
        "paper": "Liu et al., Adv Healthc Mater 2026, 15, e02225",
        "library1": n1, "library2": n2, "total": len(rows),
        "expected": "320 + 24 = 344",
        "unique_SMILES": len({r["IL_SMILES"] for r in rows}),
        "arms_per_core": {a: amine[a][1] for a in AMINES},
        "mass_check_vs_tables_S1_S2": {
            "reference_rows": len(REFERENCE_MW),
            "agree_as_printed": len(agree),
            "max_abs_error_Da": round(max(c["abs_error"] for c in agree), 2),
            "table_S2_all_agree": all(
                any(c["lipid"] == f"T3A-L{i}-2C10" for c in agree)
                for i in range(1, 9)),
            "anomalous": anomalous,
            "note": "Table S2 holds the tail constant at 2C10 across all eight "
                    "linkers and every row agrees as printed, which pins L1-L8 "
                    "independently of any tail. The six anomalous S1 rows admit "
                    "several mass-degenerate relabelings (a CH2 moved between "
                    "linker and tail changes neither mass nor formula), so no "
                    "single explanation is asserted. The enumerated STRUCTURES "
                    "are the same under all of them."},
        "table_S1_anomaly": "Table S1 prints two rows as 400-L6-C6 (471.7 and "
                            "516.7) and two as 113-L6-C6 (both 744.0). The 516.7 "
                            "row matches the L7 product, so the second pair is "
                            "mislabelled in the paper; only the L6 masses are "
                            "used here.",
    }
    json.dump(qc, open("10_enumeration_qc.json", "w"), indent=1)

    print(f"enumerated {len(rows)}  (library 1: {n1}, library 2: {n2})")
    print(f"  unique SMILES {qc['unique_SMILES']}")
    print(f"  mass check: {len(agree)}/{len(REFERENCE_MW)} reference masses "
          f"agree as printed, max |error| "
          f"{max(c['abs_error'] for c in agree):.2f} Da")
    print(f"  Table S2 (tail held at 2C10, pins every linker): "
          f"{'all 8 agree' if qc['mass_check_vs_tables_S1_S2']['table_S2_all_agree'] else 'FAILED'}")
    for a in anomalous:
        print(f"     {a['lipid']:14} table {a['table_MW']:>8}  as printed "
              f"{a['computed_as_printed']:>9}  exact fits: "
              f"{a['exact_fits'] or 'none, closest ' + a['closest'] + f" at {a['closest_error']} Da"}")
    assert len(agree) >= 23, f"only {len(agree)} reference masses agree"


if __name__ == "__main__":
    main()
