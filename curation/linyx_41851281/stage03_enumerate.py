#!/usr/bin/env python3
"""Enumerate the 116 H-T library lipids -> 03_enumerated.csv.

CHEMISTRY. Supplementary Methods p.55: "All ionisable lipids were synthesized
through the Michael addition reaction. Taking lipid P1 as an example, H4
(3 mmol) and T1 (7.5 mmol) were stirred to react at 120 C for 24 h without
solvent ... the 3-tailed lipid P1 was purified ... from impurities and 4-tailed
byproducts". So the (n) suffix is the NUMBER OF TAIL ARMS, each arm being one
aza-Michael addition of an acrylate (ester, 'O') or acrylamide (amide, 'N')
across one N-H. The head's OH, where it has one, is not an aza-Michael site and
survives -- every H1, H3 and H11 exact mass confirms it.

ASSEMBLY. Every bond is made by curation/lib/assemble.attach, which adds a bond
and removes no atom, so the F-013 RemoveAtom mechanism cannot occur here. No
string concatenation builds a structure. An attachment atom is addressed by an
atom map, never by index, and a target nitrogen on an intermediate is re-tagged
through a SMARTS that must match exactly once OR match only symmetry-equivalent
atoms -- `tag()` raises otherwise.

WHAT IS A CONVENTION AND WHAT IS NOT.  When n equals the head's total N-H count
the arms have exactly one home and nothing is assumed. When n is less, WHICH
N-H carries which arm is not stated anywhere in the paper, is not visible to
HRMS (same formula either way -- the Lin regiochemistry lesson, CLAUDE.md s5),
and is not drawn. Those rows get positional_ambiguity = 1 and DECISIONS.md D-1
is OPEN. The convention used is MAXIMAL DISPERSION: one arm on each nitrogen in
site order before any nitrogen takes a second.
"""
import csv
import json
import os
import re
import sys

from rdkit import Chem, RDLogger

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import assemble                                            # noqa: E402
from errors import StructureError                          # noqa: E402

RDLogger.DisableLog("rdApp.*")
sys.path.insert(0, HERE)
from stage01_blocks import TAILS, names, hrms_masses  # noqa: E402

# Per head: the free head written with each reactive nitrogen carrying a
# UNIQUE ISOTOPE LABEL, plus that nitrogen's N-H capacity (2 = primary amine,
# 1 = secondary). An isotope is a property of the atom itself, so it survives
# assemble.attach (which clears only atom maps) and addresses the same nitrogen
# on every intermediate -- no index, no "the first match", no SMARTS that has
# to be rewritten for each partially substituted product. The labels are
# stripped before any mass or key is computed; `strip_isotopes` asserts none
# survive, because an isotope left in place WOULD change the exact mass and
# the InChIKey and would still parse.
SITES = {
    "H1":  ("OCCCC[15NH2]", [(15, 2)]),
    "H2":  ("CN(C)CCC[15NH2]", [(15, 2)]),
    "H3":  ("OCC[16NH]CC[15NH2]", [(15, 2), (16, 1)]),
    "H4":  ("[15NH2]CCN(C)CC[16NH2]", [(15, 2), (16, 2)]),
    "H5":  ("[15NH2]CCN(CC)CC[16NH2]", [(15, 2), (16, 2)]),
    "H6":  ("C[15NH]CCN(C)CC[16NH]C", [(15, 1), (16, 1)]),
    "H7":  ("CN(C)CC[16NH]CC[15NH2]", [(15, 2), (16, 1)]),
    "H8":  ("CN(C)CCC[16NH]CCC[15NH2]", [(15, 2), (16, 1)]),
    "H9":  ("[15NH2]CCCN(C)CCC[16NH2]", [(15, 2), (16, 2)]),
    "H10": ("[15NH2]CCN1CC[16NH]CC1", [(15, 2), (16, 1)]),
    "H11": ("[15NH2]CCN(CCCO)CC[16NH2]", [(15, 2), (16, 2)]),
    "H12": ("CC(C)(C)OC(=O)N1CCN(CC[15NH2])CC1", [(15, 2)]),
    "H13": ("[15NH2]CCN(CC[16NH2])CCCN1CCCCCC1", [(15, 2), (16, 2)]),
    "H14": ("[15NH2]CCN(CC[16NH2])CCCN1CCN(C)CC1", [(15, 2), (16, 2)]),
}


def tag(smiles, isotope, mapno, label=""):
    """Put atom map `mapno` on the atom carrying `isotope`. One, or refuse."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise StructureError("tag(%s): unparseable %r" % (label, smiles))
    hits = [a.GetIdx() for a in mol.GetAtoms() if a.GetIsotope() == isotope]
    if len(hits) != 1:
        raise StructureError(
            "tag(%s): %d atoms carry isotope %d in %s, need exactly 1."
            % (label, len(hits), isotope, smiles))
    mol.GetAtomWithIdx(hits[0]).SetAtomMapNum(mapno)
    return Chem.MolToSmiles(mol)


def strip_isotopes(smiles, label=""):
    mol = Chem.MolFromSmiles(smiles)
    for a in mol.GetAtoms():
        a.SetIsotope(0)
    out = Chem.MolToSmiles(mol)
    if any(a.GetIsotope() for a in Chem.MolFromSmiles(out).GetAtoms()):
        raise StructureError("strip_isotopes(%s) left a label behind" % label)
    return out


def build_arm(linker, tail):
    core = "[CH3:2]CC(=O)[OH:4]" if linker == "O" else "[CH3:2]CC(=O)[NH2:4]"
    return assemble.attach(core, TAILS[tail][0], a_map=4, b_map=3,
                           label="arm %s%s" % (linker, tail))


def plan(head, n):
    """Arm -> site-index assignment. MAXIMAL DISPERSION, the D-1 convention."""
    caps = [c for _, c in SITES[head][1]]
    total = sum(caps)
    if n > total:
        raise StructureError("%s has %d N-H, cannot take %d arms"
                             % (head, total, n))
    used = [0] * len(caps)
    order = []
    while len(order) < n:
        best = min(range(len(caps)), key=lambda i: (used[i], i)
                   if used[i] < caps[i] else (99, i))
        if used[best] >= caps[best]:
            raise StructureError("plan(%s,%d) stuck" % (head, n))
        used[best] += 1
        order.append(best)
    return order, total


def build_from_counts(head, linker, tail, counts):
    """Build with an explicit per-site arm count. Used by isomer_count."""
    arm = build_arm(linker, tail)
    cur, sites = SITES[head]
    for i, c in enumerate(counts):
        for _ in range(c):
            cur = assemble.attach(tag(cur, sites[i][0], 1), arm,
                                  a_map=1, b_map=2, label=head)
    return assemble.canonical(strip_isotopes(cur, head))


def isomer_count(head, linker, tail, n):
    """How many DISTINCT structures put n arms on this head?

    Not `n < total`. On H4 (two primary amines, four N-H) three arms can only
    be 2+1, and 2+1 and 1+2 are the same molecule, so a 3-tail H4 lipid is not
    ambiguous at all -- flagging it would be a false caveat. This enumerates
    every per-site assignment, builds each, and counts distinct InChIKeys.
    """
    caps = [c for _, c in SITES[head][1]]
    out = {}

    def rec(i, left, acc):
        if i == len(caps):
            if left == 0:
                k = Chem.MolToInchiKey(Chem.MolFromSmiles(
                    build_from_counts(head, linker, tail, acc)))
                out.setdefault(k, list(acc))
            return
        for c in range(min(caps[i], left) + 1):
            rec(i + 1, left - c, acc + [c])

    rec(0, n, [])
    return out


def build(head, linker, tail, n):
    arm = build_arm(linker, tail)
    order, total = plan(head, n)
    cur, sites = SITES[head]
    for i in order:
        iso = sites[i][0]
        cur_t = tag(cur, iso, 1, label="%s site %d" % (head, iso))
        cur = assemble.attach(cur_t, arm, a_map=1, b_map=2,
                              label="%s-%s%s(%d)" % (head, linker, tail, n))
        # attach cleared the MAP but not the isotope, so the same nitrogen is
        # still addressable for a second arm.
    return assemble.canonical(strip_isotopes(cur, head)), total


def main():
    nm = names()
    masses = hrms_masses()
    PNAMES = {"P1": "H4-OT1(3)", "P2": "H9-OT7(3)", "P3": "H6-OT13(2)",
              "P4": "H9-OT7(4)", "P5": "H2-OT5(2)", "P6": "H1-OT7(2)"}
    rows, nbad = [], 0
    for mlid in sorted(nm, key=int):
        name = nm[mlid]
        m = re.match(r"H(\d+)-([ON])T(\d+)\((\d)\)$", name)
        if not m:
            continue
        head, lk, tl, n = ("H" + m.group(1), m.group(2), "T" + m.group(3),
                           int(m.group(4)))
        smi, total = build(head, lk, tl, n)
        isos = isomer_count(head, lk, tl, n)
        if smi not in [None]:
            pass
        mol = Chem.MolFromSmiles(smi)
        from rdkit.Chem.Descriptors import ExactMolWt
        from rdkit.Chem.rdMolDescriptors import CalcMolFormula
        calc = ExactMolWt(mol)
        printed = masses.get(mlid)
        dev = None if printed is None else calc - printed
        if dev is not None and abs(dev) > 5e-3:
            nbad += 1
        plabel = [k for k, v in PNAMES.items() if v == name]
        rows.append(dict(
            ml_id=mlid, lipid_name=name, p_label=plabel[0] if plabel else "",
            head=head, linker={"O": "ester", "N": "amide"}[lk], tail=tl,
            n_tails=n, n_NH_sites=total,
            all_sites_used=int(n == total),
            n_distinct_isomers=len(isos),
            positional_ambiguity=int(len(isos) > 1),
            smiles=smi, formula=CalcMolFormula(mol),
            exact_mass="%.5f" % calc,
            hrms_printed="" if printed is None else "%.4f" % printed,
            hrms_dev_Da="" if dev is None else "%+.4f" % dev,
            inchikey=Chem.MolToInchiKey(mol)))
    # the six AI-designed lipids carry their own HRMS under a P label
    for p, name in PNAMES.items():
        for r in rows:
            if r["lipid_name"] == name:
                r["p_label"] = p
                pm = masses.get(p)
                if pm is not None:
                    r["hrms_printed"] = "%.4f" % pm
                    r["hrms_dev_Da"] = "%+.4f" % (float(r["exact_mass"]) - pm)
    cols = list(rows[0].keys())
    with open(os.path.join(HERE, "03_enumerated.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    uniq = len({r["inchikey"] for r in rows})
    print("enumerated %d lipids, %d unique InChIKeys, %d HRMS deviations >5 mDa"
          % (len(rows), uniq, nbad))
    amb = [r["ml_id"] for r in rows if r["positional_ambiguity"]]
    print("positional_ambiguity on %d: %s" % (len(amb), amb))


if __name__ == "__main__":
    main()
