#!/usr/bin/env python3
"""Unit tests for enumerate.py's two assembly steps.

    python3 -I curation/siegwart_42634137/test_enumerate.py

Both functions under test have shipped a silent structural bug, and in both
cases the bug was found by a human reading figures, 537 rows after the fact.
These tests fail at the function instead.

WHY EACH TEST EXISTS

`_cholesteryl_ester` has been wrong twice:
  1. `CHOL.replace("(O)", "", 1)` plus string concatenation bonded the acyl to a
     terminal methyl of the isooctyl side chain and left C3 a radical.
  2. The RWMol fix for that bonded the acid oxygen to C3 and then deleted
     cholesterol's hydroxyl with `RemoveAtom`, which drops that oxygen out of
     C3's neighbour list WITHOUT correcting C3's chiral tag -- silently
     inverting it. All 537 corpus rows carried epicholesterol (3alpha).
     Molecular formula, exact mass and all 17 descriptor columns are
     bit-identical between the epimers, so every composition check passed, and
     the C10 anchor passed because it was a byte copy of this function's output.
  So: compare against an INDEPENDENT plain-string oracle that uses no RWMol at
  all, and compare the FULL CIP label vector, not just C3 -- a check on C3 alone
  would miss damage to the Delta-5 alkene or the side-chain centres.

`_quaternization_site` has also been wrong twice:
  1. `max(atom index)`, an artifact of SMILES writing order, which for C6
     selected the least nucleophilic nitrogen in the molecule.
  2. The repaired rule counted only beta-ESTER arms (N-C-C-C=O), so it could
     not see Lipid 14's beta-ACYLOXY arm (N-C-C-O-C=O). Both of that molecule's
     nitrogens scored 0 and the rule fell through to the atom-index tiebreak --
     reintroducing the ordering artifact it was written to remove. Writing the
     same molecule tail-first made it pick the other nitrogen.
  So: assert the choice is INVARIANT across several spellings of the same
  molecule. A test that only checks the answer on our canonical SMILES cannot
  see an ordering artifact, because our canonical SMILES is one ordering.
"""
import importlib.util
import os
import sys

from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors as rd

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
CHOLESTEROL_INCHIKEY = "HVYWMOMLDIMFJA-DPAQBDIFSA-N"   # 3beta; also in reagent_registry.json
EPICHOLESTEROL_INCHIKEY = "HVYWMOMLDIMFJA-VEIPTCAHSA-N"  # 3alpha; the shipped bug


def load():
    spec = importlib.util.spec_from_file_location(
        "sieg_enumerate", os.path.join(HERE, "enumerate.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def cip_vector(smiles):
    """Every assigned CIP label, in canonical-rank order. Orientation-free."""
    m = Chem.MolFromSmiles(smiles)
    Chem.AssignStereochemistry(m, cleanIt=True, force=True)
    ranks = list(Chem.CanonicalRankAtoms(m))
    out = [(ranks[a.GetIdx()], a.GetProp("_CIPCode"))
           for a in m.GetAtoms() if a.HasProp("_CIPCode")]
    return tuple(c for _, c in sorted(out))


def oracle(n_link, chol):
    """The same ester, built by plain string substitution. No RWMol anywhere.

    Independent by construction: it cannot share a bug with fragment assembly
    because it does no fragment assembly. This is the point -- the C10 anchor
    failed precisely because it was NOT independent of the code it checked.
    """
    acyl = "C" * (n_link - 1)
    assert chol.count("(O)") == 1, "CHOL no longer has exactly one (O) to substitute"
    return Chem.MolToSmiles(Chem.MolFromSmiles(
        chol.replace("(O)", f"(OC(=O){acyl})", 1)))


def describe_site(mol, idx):
    """The chosen nitrogen by substituent sizes -- identifies it without
    depending on an atom index, which is the thing under test."""
    sizes = []
    for nb in mol.GetAtomWithIdx(idx).GetNeighbors():
        seen, stack, n = {idx}, [nb], 0
        while stack:
            a = stack.pop()
            if a.GetIdx() in seen:
                continue
            seen.add(a.GetIdx())
            n += 1
            stack += list(a.GetNeighbors())
        sizes.append(n)
    return (mol.GetAtomWithIdx(idx).IsInRing(), tuple(sorted(sizes)))


def main():
    E = load()
    fails = []

    def check(ok, label, detail=""):
        print(f"  {'PASS' if ok else 'FAIL':4}  {label}{'  ' + detail if detail else ''}")
        if not ok:
            fails.append(label)

    print("CHOL input constant")
    k = Chem.MolToInchiKey(Chem.MolFromSmiles(E.CHOL))
    check(k == CHOLESTEROL_INCHIKEY,
          "enumerate.CHOL is cholesterol (3beta), not an epimer", k)

    print("\n_cholesteryl_ester vs an independent string oracle")
    for n in (3, 5, 7):
        got, _ = E._cholesteryl_ester(n)
        want = oracle(n, E.CHOL)
        # the tagged attachment atom is not part of the chemistry
        got_clean = Chem.MolFromSmiles(got)
        for a in got_clean.GetAtoms():
            a.SetAtomMapNum(0)
            a.SetNoImplicit(False)
            a.SetNumExplicitHs(0)
        got_clean = Chem.MolToSmiles(got_clean)
        check(got_clean == want, f"linker {n}: exact match to the oracle")
        check(cip_vector(got_clean) == cip_vector(want),
              f"linker {n}: FULL CIP vector matches (not just C3)",
              str(cip_vector(got_clean)))
        check(rd.CalcMolFormula(Chem.MolFromSmiles(got_clean))
              == rd.CalcMolFormula(Chem.MolFromSmiles(want)),
              f"linker {n}: formula matches (necessary, NOT sufficient)")
        check(E.ATTACH_MAP in [a.GetAtomMapNum() for a in Chem.MolFromSmiles(got).GetAtoms()],
              f"linker {n}: the free chain end is tagged, not found by position")

    print("\n_quaternization_site is independent of SMILES atom ordering")
    PARENTS = {
        "Lipid 14": "CCCCCCCCC(CCCCCCCC)OC(=O)CCN(CCCN(C)C)CCC(=O)OC(CCCCCCCC)CCCCCCCC",
        "MC3": "CCCCC/C=C\\C/C=C\\CCCCCCCCC(CCCCCCCC/C=C\\C/C=C\\CCCCC)OC(=O)CCCN(C)C",
        "SM-102": "CCCCCCCCC(=O)OCCCCCCCCN(CCCCCCCC)CCCCCCO",
    }
    for name, smi in PARENTS.items():
        m0 = Chem.MolFromSmiles(smi)
        spellings = [Chem.MolToSmiles(m0),
                     Chem.MolToSmiles(m0, canonical=False),
                     Chem.MolToSmiles(m0, rootedAtAtom=m0.GetNumAtoms() - 1,
                                      canonical=False),
                     Chem.MolToSmiles(m0, rootedAtAtom=0, canonical=False)]
        picks = set()
        for s in spellings:
            m = Chem.MolFromSmiles(s)
            tgt, _ = E._quaternization_site(m)
            picks.add(describe_site(m, tgt.GetIdx()) if tgt else None)
        check(len(picks) == 1,
              f"{name}: one site across {len(spellings)} spellings",
              str(sorted(picks)))

    print("\n_quaternization_site regression against the figures")
    for key in ("C5", "C6"):
        hs, _, _, narm = E.HEADS[key]
        core, _ = E.michael(hs, E.ALKYLS["O2"][0], narm)
        tgt, note = E._quaternization_site(Chem.MolFromSmiles(core))
        check(tgt is not None and tgt.IsInRing(),
              f"{key}O2 still resolves to the piperazine N drawn in Fig 2a", note)
    # Lipid 14: the paper draws no -C conjugate, so this is the chemically
    # preferred site, not a verified one. The row carries structure_assumed = 1.
    tgt, note = E._quaternization_site(Chem.MolFromSmiles(PARENTS["Lipid 14"]))
    ring, sizes = describe_site(Chem.MolFromSmiles(PARENTS["Lipid 14"]), tgt.GetIdx())
    check(sizes[:2] == (1, 1),
          "Lipid 14 resolves to the terminal N(CH3)2, not the hindered central N",
          f"substituents {sizes}")
    check("ASSUMED" in note,
          "Lipid 14's site is reported as an ASSUMPTION (the paper draws none)")

    print("\nend to end: every enumerated lipid liberates cholesterol")
    import csv
    p = os.path.join(HERE, "03_enumerated.csv")
    if os.path.exists(p):
        from rdkit.Chem import AllChem
        rx = AllChem.ReactionFromSmarts(
            "[C:1](=[O:2])[O:3][C:4]>>[C:1](=[O:2])[OH].[OH][C:4]")
        rows = list(csv.DictReader(open(p)))
        check(len(rows) == 36, f"03_enumerated.csv has 36 lipids", str(len(rows)))
        bad, nosite = [], []
        for r in rows:
            keys = set()
            for prod in rx.RunReactants((Chem.MolFromSmiles(r["IL_SMILES"]),)):
                for x in prod:
                    try:
                        Chem.SanitizeMol(x)
                    except Exception:
                        continue
                    if rd.CalcMolFormula(x) == "C27H46O":
                        keys.add(Chem.MolToInchiKey(x))
            if not keys:
                nosite.append(r["lipid_name"])     # zero sites must FAIL, not pass
            elif keys != {CHOLESTEROL_INCHIKEY}:
                bad.append((r["lipid_name"], sorted(keys)))
        check(not nosite, "every lipid has a cleavable steroid ester", str(nosite[:4]))
        check(not bad, "every liberated sterol is cholesterol (3beta)", str(bad[:3]))
    else:
        check(False, "03_enumerated.csv exists")

    print(f"\n{'ALL PASS' if not fails else f'{len(fails)} FAILURE(S): ' + '; '.join(fails)}")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
