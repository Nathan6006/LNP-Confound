#!/usr/bin/env python3
"""BUILD the C10 anchors for KJ_2026. Do not store a string -- derive it.

    python3 -I make_reference.py          -> reference_structures.csv

C10 compares each enumerated structure against the structure the paper DRAWS.
It is the only Tier-1 check that tests topology against external ground truth,
and it only tests anything if the anchor was built by a DIFFERENT CODE PATH
than the enumerator. That independence has been asserted twice in this project
and been false twice (F-013, a byte copy; F-037, where the file written to fix
F-013 re-copied the other half of the molecule).

HOW THIS ONE IS DIFFERENT FROM stage03_enumerate.py, concretely:

  stage03                              this file
  -------                              ---------
  imports curation/lib/assemble        imports NOTHING but rdkit, csv, ast
  builds with esterify() + attach(),   builds by writing one SMILES string
    i.e. RWMol, CombineMols, AddBond,    per molecule, head-first; there is
    RemoveAtom                           no RWMol anywhere
  composes from the stage02 block      takes nothing from stage02; the
    table, acid-end first                retinoid is written RING-END FIRST,
                                         the opposite traversal
  starts at the retinoid acid and      starts at the retinoid's cyclohexene
    finishes at the head arm             and finishes at the lipid tail

TWO ANCHORS, OF DIFFERENT CHAIN LENGTHS AND DIFFERENT LINKER ARCHITECTURES.
One length cannot catch an off-by-one that scales with n, and one linker
cannot catch an error in the other:

  AT1-c-ATRA  Scheme S5, drawn as a product. The 'c' linker (two esters), the
              T1 tail (C6 acid on pentadecan-8-ol), the all-trans polyene.
  AT4-Bex     Scheme S9, drawn as a product. No 'c' linker (one ester), the
              T4 tail whose ester is REVERSED relative to T1-T3, and an
              aromatic hydrophobe with no polyene at all.

BOTH GUARDS BELOW ARE MADE TO FAIL ON PURPOSE by seed_checks.py, in a child
process, with a clean control that must come back clean. A guard that has
never fired is not a guard: two earlier versions of C16 reported "all derived"
while examining nothing.
"""
import ast
import csv
import os
import sys

from rdkit import Chem
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "reference_structures.csv")
ENUM = os.path.join(HERE, "03_enumerated.csv")

# Anything that would make this file a second copy of the enumerator.
BANNED_IMPORTS = {"assemble", "stage02_blocks", "stage03_enumerate",
                  "stage04_values", "chains", "reagents"}
BANNED_CALLS = {"esterify", "attach", "RWMol", "AddBond", "RemoveAtom",
                "CombineMols", "ReplaceSubstructs", "canonical"}

# -------------------------------------------------------------------------
# THE ANCHORS, written out head-first as single SMILES strings.
#
# Read off the rendered Schemes S5 and S9 at native resolution, after
# compositing the media PNG onto white (these PNGs carry an alpha channel and
# a straight RGB conversion renders black-on-black).
#
# Traversal, deliberately opposite to the enumerator's: begin at the
# retinoid's ring, walk the polyene OUT to the carboxyl, through the linker's
# esters, into the nitrogen, out the hydroxybutyl head arm, then down the
# tail. The enumerator begins at the retinoid acid, attaches the tail, and
# finishes at the head arm.
# -------------------------------------------------------------------------
ANCHORS = [
    {
        "lipid_name": "AT1-c-ATRA",
        "source_figure": "SI Scheme S5, row 2, the drawn product",
        "smiles":
            "CC1=C(C(C)(C)CCC1)"            # the trimethylcyclohexene ring
            "/C=C/C(C)=C/C=C/C(C)=C/"       # the all-trans polyene, ring->acid
            "C(=O)"                         # C15, the retinoid carbonyl
            "OCCCCOC(=O)CCCCC"              # the 'c' linker: two esters
            "N(CCCCO)"                      # the amine and its head arm
            "CCCCCC(=O)OC(CCCCCCC)CCCCCCC",  # the T1 tail: 6-bromohexanoate
                                             # gives SIX carbons from the
                                             # nitrogen to the carbonyl
                                             # inclusive, not five -- the
                                             # first version of this line
                                             # wrote five and the InChIKey
                                             # comparison below caught it
                                             # (C54H93NO7 against the
                                             # enumeration's C55H95NO7, and
                                             # the SI's printed [M+H]+
                                             # 882.719 says the enumeration
                                             # is the correct half)
        "note":
            "Built by writing one SMILES string, ring-end first. No "
            "curation/lib import, no RWMol, no fragment assembly. The arm "
            "that is one CH2 longer here than in the AT4-Bex anchor is the "
            "linker, not the head arm.",
    },
    {
        "lipid_name": "AT4-Bex",
        "source_figure": "SI Scheme S9, bottom row, the drawn product",
        "smiles":
            "CC1=CC2=C(C=C1"                # the tetrahydronaphthalene
            "C(=C)c1ccc(cc1)C(=O)"          # the exocyclic alkene and the
                                            # para-benzoate carbonyl
            "OCCCCCC"                       # the linker, one ester only
            "N(CCCCO)"                      # the amine and its head arm
            "CCCCCCOC(=O)C(CCCCCC)CCCCCCCC)"  # the T4 tail, ester REVERSED
            "C(C)(C)CCC2(C)C",
        "note":
            "Built by writing one SMILES string, aromatic-ring first. The T4 "
            "ester points the other way than T1's and the formula is "
            "identical either way, so this anchor is the one that can see a "
            "reversed ester.",
    },
]


def _ast_guard(path):
    """Walk THIS module's AST and refuse if it imports or calls the enumerator.

    Inspects Import / ImportFrom / Call / Attribute nodes ONLY, so the ban
    list -- which is itself a set of strings in this file -- cannot match
    itself. One earlier guard elsewhere in the project failed on its own
    module for exactly that reason, and another's substring scan flagged its
    own docstring.
    """
    tree = ast.parse(open(path, encoding="utf-8").read(), filename=path)
    bad = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for al in node.names:
                if al.name.split(".")[0] in BANNED_IMPORTS:
                    bad.append(f"import {al.name}")
        elif isinstance(node, ast.ImportFrom):
            if node.module and node.module.split(".")[0] in BANNED_IMPORTS:
                bad.append(f"from {node.module} import ...")
        elif isinstance(node, ast.Call):
            f = node.func
            name = (f.attr if isinstance(f, ast.Attribute)
                    else f.id if isinstance(f, ast.Name) else None)
            if name in BANNED_CALLS:
                bad.append(f"call {name}()")
        elif isinstance(node, ast.Attribute):
            if node.attr in BANNED_CALLS:
                bad.append(f"attribute .{node.attr}")
    if bad:
        raise SystemExit(
            "make_reference: THIS FILE IS NOT INDEPENDENT OF THE ENUMERATOR. "
            f"Its own AST contains: {sorted(set(bad))}. An anchor built by "
            "the code it is supposed to check compares a value to itself, "
            "which is F-013 and F-037. Rewrite the derivation.")
    return len(tree.body)


def _template_guard(path, enum_smiles):
    """No long string literal here may be a substring of an enumerated SMILES.

    F-037's AMINE_HALF was the enumeration's canonical SMILES with the sterol
    swapped for '%' -- byte-identical over 141 characters. Byte-identity of
    the EMITTED anchor is not the signal (a correct re-derivation canonicalises
    to the same string), but a long literal lifted out of the enumerator's
    output IS.
    """
    tree = ast.parse(open(path, encoding="utf-8").read(), filename=path)
    hits = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Constant) and isinstance(node.value, str):
            lit = node.value
            if len(lit) < 25:
                continue
            for s in enum_smiles:
                if lit in s:
                    hits.append((lit[:60], len(lit)))
                    break
    if hits:
        raise SystemExit(
            "make_reference: a string literal of >=25 characters in this file "
            f"is a SUBSTRING of an enumerated SMILES: {hits}. That literal "
            "came from the enumerator's output, so the anchor is a copy of "
            "the thing it checks. That is F-037 exactly.")
    return True


def build():
    enum = {r["IL_name"]: r["IL_SMILES"] for r in csv.DictReader(open(ENUM))}
    _ast_guard(os.path.abspath(__file__))
    _template_guard(os.path.abspath(__file__), list(enum.values()))

    rows = []
    for a in ANCHORS:
        m = Chem.MolFromSmiles(a["smiles"])
        if m is None:
            raise SystemExit(f"make_reference: {a['lipid_name']} will not "
                             f"parse: {a['smiles']}")
        rows.append([a["lipid_name"], Chem.MolToSmiles(m), a["source_figure"],
                     "2026-11-10", "stage: curation/liu_42585110", "",
                     a["note"]])
    return rows


def compare():
    """Anchor vs enumeration, by FULL InChIKey. Reported, not silently passed."""
    enum = {r["IL_name"]: r["IL_SMILES"] for r in csv.DictReader(open(ENUM))}
    out = []
    for a in ANCHORS:
        ref = Chem.MolFromSmiles(a["smiles"])
        got = Chem.MolFromSmiles(enum[a["lipid_name"]])
        kr, kg = Chem.MolToInchiKey(ref), Chem.MolToInchiKey(got)
        out.append({
            "lipid": a["lipid_name"],
            "anchor_inchikey": kr, "enumeration_inchikey": kg,
            "exact_match": kr == kg,
            "skeleton_match": kr.split("-")[0] == kg.split("-")[0],
            "anchor_canonical": Chem.MolToSmiles(ref),
            "byte_identical_to_enumeration":
                Chem.MolToSmiles(ref) == Chem.MolToSmiles(got),
        })
    return out


if __name__ == "__main__":
    rows = build()
    with open(OUT, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["lipid_name", "reference_SMILES", "source_figure",
                    "read_on", "read_by", "expected_mismatch", "note"])
        w.writerows(rows)
    print(f"wrote {len(rows)} anchor(s) to {OUT}")
    bad = 0
    for c in compare():
        print(f"  {c['lipid']:<12} exact={c['exact_match']} "
              f"skeleton={c['skeleton_match']}")
        print(f"      anchor      {c['anchor_inchikey']}")
        print(f"      enumeration {c['enumeration_inchikey']}")
        print(f"      byte-identical after canonicalisation: "
              f"{c['byte_identical_to_enumeration']}  <- NOT the signal; a "
              "correct re-derivation IS byte-identical once canonicalised "
              "(F-037)")
        if not c["exact_match"]:
            bad += 1
    if bad:
        raise SystemExit(f"make_reference: {bad} anchor(s) DISAGREE with the "
                         "enumeration. That is a connectivity or "
                         "stereochemistry error in one of the two, and mass "
                         "cannot tell you which.")
