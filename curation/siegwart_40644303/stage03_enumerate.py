#!/usr/bin/env python3
"""Enumerate the library: 465 screened lipids + 24 in-vivo variants.

    python3 -I stage03_enumerate.py

THE CHEMISTRY, FROM SCHEME 1 AND THE METHODS -- NOT FROM A MASS.
  1. aza-Michael: each core N-H adds across the ACRYLATE end of
     2-(acryloyloxy)ethyl methacrylate (AEMA, "G1"), 1.1 equiv per arm, BHT,
     50 C.
  2. sulfa-Michael: the thiol adds across the METHACRYLATE end, DMPP, 50 C,
     so sulfur lands on the beta carbon and the alpha carbon keeps the methyl.
Scheme 1 draws the product 3A4-SC5 with exactly that connectivity. The swapped
assignment has the IDENTICAL molecular formula and exact mass, so no mass check
can see it; this is the Lin regiochemistry class, and the only evidence is the
written synthesis.

EVERY BOND IS MADE WITH curation/lib/assemble.attach. Nothing is deleted, so
the F-013 mechanism (RWMol.RemoveAtom dropping an atom out of a stereocentre's
neighbour list) cannot occur, and attach asserts no radicals and stereo
preservation on every one of the ~2,000 bonds this file forms. No SMILES string
is concatenated.

HOW A NITROGEN IS FOUND AGAIN FOR ITS SECOND ARM. attach clears the atom map on
the atoms it bonds, and canonicalisation renumbers everything, so a map number
cannot survive to mark the same nitrogen twice. The cores therefore carry
ISOTOPE tags (101..105), which attach does not touch and canonicalisation
preserves; the tag is copied to an atom map immediately before each bond and
the tags are stripped at the end. Bonding to "atom index 0 of the fragment" is
what the library forbids, and this is the alternative that does not need it.

PARTIAL SUBSTITUTION IS A CONVENTION, AND IS FLAGGED AS ONE. Five cores carry
fewer arms than they have N-H: 5A2 (5 of 7), 6A1 (6 of 7), 5A3 (5 of 6), 5A4
(5 of 6) and 5A5 (5 of 6). The paper says only how many arms, never which
nitrogen keeps its hydrogen. The convention used here is "primary N-H first,
then secondary from the terminus inward", stated once, applied once, and
carried onto every affected row as positional_ambiguity = 1. For 5A4, 5A5 and
5A3 the three primaries are symmetry-equivalent so that convention has a unique
product; for 5A2 and 6A1 it does not, and those rows say so.
"""
import csv
import json
import os
import sys

from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, rdMolDescriptors

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))

import assemble   # noqa: E402

RDLogger.DisableLog("rdApp.*")

BLOCKS = os.path.join(HERE, "01_blocks.csv")
OUT = os.path.join(HERE, "03_enumerated.csv")
OUT_QC = os.path.join(HERE, "03_enumerate_qc.json")

# The arm, as Scheme 1 draws it, with two attachment points:
#   map 1  the carbon that bonds the core nitrogen (the acrylate beta carbon)
#   map 2  the beta CARBON that bonds the tail's sulfur
# NOT the arm's own sulfur: the arm has none. The first version of this file
# wrote "...C(C)C[SH:2]" and bonded that S to the tail's S, producing 465
# DISULFIDES. Nothing in the assembly complained -- attach's radical and
# stereo assertions both passed and every structure parsed -- and it was the
# corpus InChIKey comparison that caught it. A fragment that carries a
# plausible but wrong attachment atom is invisible to a build-time check.
ARM = "[CH3:1]CC(=O)OCCOC(=O)C(C)[CH3:2]"

CORES_FILE_KINDS = {"core", "tail", "linker"}

VITRO_CORES = ["3A1", "3A2", "3A3", "3A4", "3A5", "4A1", "4A3", "4A4",
               "5A1", "5A2", "5A3", "5A4", "6A1", "6A2", "6A3"]
VITRO_TAILS = ([f"SC{i}" for i in range(1, 13)] + ["SC14", "SC16", "SC18"]
               + [f"SN{i}" for i in range(1, 9)]
               + [f"SO{i}" for i in range(1, 9)])
VIVO_CORES = ["3A4", "3A4V", "3A4V2", "4A3", "4A5", "5A4", "5A5"]
VIVO_TAILS = ["SC5", "SC6", "SC7", "SC8", "SC9", "SC10"]

# Imported from stage01 so the cores are stated once.
sys.path.insert(0, HERE)
from stage01_blocks import CORES  # noqa: E402


def arm_with_tail(tail_smiles):
    """Attach the tail's sulfur to the arm's beta carbon (sulfa-Michael)."""
    m = Chem.MolFromSmiles(tail_smiles)
    if m is None:
        raise SystemExit(f"stage03: unparseable tail {tail_smiles!r}")
    s = [a.GetIdx() for a in m.GetAtoms() if a.GetSymbol() == "S"]
    if len(s) != 1:
        raise SystemExit(
            f"stage03: tail {tail_smiles!r} has {len(s)} sulfur atoms, "
            "expected exactly 1. REFUSING: which sulfur the thiol reacts "
            "through is the bond this whole library is built on.")
    m.GetAtomWithIdx(s[0]).SetAtomMapNum(3)
    return assemble.attach(ARM, Chem.MolToSmiles(m), a_map=2, b_map=3,
                           label="sulfa-Michael")


def build(core_name, tail_smiles):
    smi, plan, n_nh, _src = CORES[core_name]
    arm = arm_with_tail(tail_smiles)
    cur = smi
    for tag, k in plan.items():
        for _ in range(k):
            m = Chem.MolFromSmiles(cur)
            hit = [a for a in m.GetAtoms() if a.GetIsotope() == tag]
            if len(hit) != 1:
                raise SystemExit(
                    f"stage03/{core_name}: isotope tag {tag} matches "
                    f"{len(hit)} atoms. REFUSING.")
            hit[0].SetAtomMapNum(1)
            cur = assemble.attach(Chem.MolToSmiles(m), arm, a_map=1, b_map=1,
                                  label=f"{core_name} aza-Michael")
    m = Chem.MolFromSmiles(cur)
    for a in m.GetAtoms():
        a.SetIsotope(0)
    return assemble.canonical(Chem.MolToSmiles(m))


def main():
    blocks = {r["block"]: r for r in csv.DictReader(open(BLOCKS))}
    tails = {k: v["smiles"] for k, v in blocks.items() if v["kind"] == "tail"}

    members = [(c, t, "vitro_465") for c in VITRO_CORES for t in VITRO_TAILS]
    extra = [(c, t, "vivo_42") for c in VIVO_CORES for t in VIVO_TAILS
             if (c, t, "vitro_465") not in members]
    members += extra

    rows, seen = [], {}
    for core, tail, origin in members:
        smi = build(core, tails[tail])
        mol = Chem.MolFromSmiles(smi)
        plan = CORES[core][1]
        n_arms = sum(plan.values())
        partial = n_arms < CORES[core][2]
        rows.append({
            "lipid_name": f"{core}-{tail}",
            "core": core, "tail": tail, "origin": origin,
            "n_arms": n_arms, "n_NH_in_core": CORES[core][2],
            "smiles": smi,
            "inchikey": Chem.MolToInchiKey(mol),
            "formula": rdMolDescriptors.CalcMolFormula(mol),
            "exact_mass": round(
                Descriptors.ExactMolWt(mol), 4),
            "positional_ambiguity": 1 if partial else 0,
            "confidence": "medium" if partial else "high",
            "note": ("the core carries %d arms on %d N-H and the paper does "
                     "not say which nitrogen keeps its hydrogen; the stored "
                     "structure follows the stated convention (primary N-H "
                     "first, then secondary from the terminus inward)"
                     % (n_arms, CORES[core][2])) if partial else "",
        })
        seen.setdefault(smi, []).append(f"{core}-{tail}")

    dup = {k: v for k, v in seen.items() if len(v) > 1}
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    n_vitro = sum(1 for r in rows if r["origin"] == "vitro_465")
    qc = {
        "_what": "Every library member, built bond by bond with "
                 "curation/lib/assemble.attach.",
        "n_total": len(rows),
        "n_vitro_screen": n_vitro,
        "abstract_claims": 465,
        "reconciliation":
            "15 amine cores x 31 thiol tails = 465, and the figure's own row "
            "and column labels are exactly those 15 and those 31. The "
            "abstract's 465 is reproduced exactly, with no residual.",
        "n_new_invivo_variants": len(rows) - n_vitro,
        "invivo_total_42":
            "7 cores x 6 tails = 42, matching the paper's '42 lipids in "
            "total, including the parent lipids'. 18 of the 42 (3A4, 4A3 and "
            "5A4 x SC5-SC10) are already members of the 465, so only "
            + str(len(rows) - n_vitro) + " structures are new.",
        "n_unique_smiles": len(seen),
        "duplicate_structures": dup,
        "n_positional_ambiguity": sum(r["positional_ambiguity"] for r in rows),
        "arm_fragment": ARM,
        "bonds_formed": sum(sum(CORES[r["core"]][1].values()) + 1
                            for r in rows),
    }
    with open(OUT_QC, "w") as fh:
        json.dump(qc, fh, indent=1)
    print(json.dumps(qc, indent=1))


if __name__ == "__main__":
    main()
