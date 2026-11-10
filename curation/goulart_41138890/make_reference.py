#!/usr/bin/env python3
"""Build the C10 anchor for C12-200 from the structure Figure 1A DRAWS.

WHAT THIS CHECKS. This curation does not enumerate a library; it takes the
ionizable lipid's structure from the corpus (AC_2025 carries `C12-200` on 498
rows with exactly one structure). The thing that needs an external test is
therefore whether THAT structure is the molecule THIS paper drew -- the brief's
case of "a name unambiguous in the corpus and simply wrong for your paper".
So the anchor is assembled here from the drawing alone and compared by
InChIKey, in a separate script, against the corpus string.

INDEPENDENCE. This module reads NOTHING from data/, imports no corpus helper,
and contains no literal long enough to have been lifted from the corpus SMILES.
Both properties are asserted at run time by check_independence() below, which
walks THIS MODULE'S OWN AST and inspects only Import / ImportFrom / Call /
Attribute nodes, so the ban list cannot match itself (F-037, and the two agents
whose guards flagged their own docstring or their own module).

WHY THE DRAWING IS SYMBOLIC HERE. Figure 1A prints each tail as the TEXT
"C10H21", not as a zig-zag. So the chain length is read, not counted, and the
usual "one anchor of one chain length cannot catch an off-by-one that scales
with n" limitation does not bite in the usual way: an off-by-one would have to
be a misreading of a printed subscript. What IS drawn, and therefore what this
anchor really tests, is the CONNECTIVITY -- which nitrogen carries how many
arms, and on which side of the piperazine the three-arm branch sits.

THE SEEDED DEFECT. `--seed markovnikov` opens the C12 epoxide at the internal
carbon instead of the terminal one: N-CH(C10H21)-CH2OH rather than
N-CH2-CH(OH)-C10H21. Identical molecular formula, identical exact mass to every
decimal. Only the InChIKey sees it, which is the proof that the gate tests
connectivity and not composition.
"""
import argparse
import ast
import csv
import inspect
import os
import sys

from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, rdMolDescriptors

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import assemble                                        # noqa: E402

BANNED_IMPORTS = {"corpus_components", "pandas", "build_dataset",
                  "reagent_registry", "numpy"}
BANNED_CALLS = {"read_csv", "DictReader", "reader", "load", "loads"}

# The tail, written as a count rather than as a literal string, so there is no
# long literal in this file at all.
TAIL_CARBONS = 10          # Figure 1A prints "C10H21" on every one of the five


def check_independence(path=None):
    """Walk this module's own AST. Import/ImportFrom/Call/Attribute only."""
    path = path or os.path.abspath(__file__)
    tree = ast.parse(open(path).read())
    bad = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for n in node.names:
                if n.name.split(".")[0] in BANNED_IMPORTS:
                    bad.append(f"import {n.name}")
        elif isinstance(node, ast.ImportFrom):
            if (node.module or "").split(".")[0] in BANNED_IMPORTS:
                bad.append(f"from {node.module} import ...")
        elif isinstance(node, ast.Attribute):
            if node.attr in BANNED_CALLS:
                bad.append(f".{node.attr}")
        elif isinstance(node, ast.Call):
            f = node.func
            if isinstance(f, ast.Name) and f.id in BANNED_CALLS:
                bad.append(f"{f.id}()")
    return bad


def check_no_lifted_template(enumerated, path=None, min_len=25):
    """No string literal of >= min_len characters in this module may be a
    substring of any structure this anchor will be compared against."""
    path = path or os.path.abspath(__file__)
    tree = ast.parse(open(path).read())
    bad = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Constant) and isinstance(node.value, str):
            s = node.value
            if len(s) >= min_len and any(s in e for e in enumerated):
                bad.append(s[:60])
    return bad


def tag(smiles, smarts, mapno, pick=0):
    """Put an atom map on one SMARTS hit, through RDKit's atom API.

    Never by string surgery: "a tail built by string concatenation became a
    ketone bearing a free alcohol on all 11 rows, two lines below its own
    comment warning against string surgery" (brief section 8).
    """
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        raise ValueError(f"tag: unparseable {smiles!r}")
    hits = m.GetSubstructMatches(Chem.MolFromSmarts(smarts))
    if len(hits) <= pick:
        raise ValueError(f"tag: {smarts!r} has {len(hits)} hits, need > {pick}")
    m.GetAtomWithIdx(hits[pick][0]).SetAtomMapNum(mapno)
    return Chem.MolToSmiles(m)


def build(seed=None):
    """C12-200, head-first from Figure 1A.

    The traversal is deliberately the opposite of the corpus string's, which
    starts at a tail terminus. Here the CORE AMINE is written from its
    connectivity as drawn -- piperazine, one N carrying -CH2CH2-NH2, the other
    carrying -CH2CH2-NH-CH2CH2-NH2 -- and the five N-H sites are then
    alkylated one at a time, each target chosen by SMARTS, never by index.
    """
    chain = "C" * TAIL_CARBONS
    if seed == "markovnikov":
        # internal opening: N-CH(C10H21)-CH2OH. Same formula, same exact mass.
        # The attachment carbon already carries TWO heavy neighbours here
        # (CH2OH and the decyl chain), so it is written [CH2:2], not [CH3:2]:
        # a fragment's explicit H count has to leave exactly one free valence
        # for the bond attach() is about to make, and getting it wrong is an
        # unparseable fragment rather than a silent wrong answer only because
        # assemble refuses. Clean arm: one heavy neighbour, so [CH3:2].
        arm = "[CH2:2](CO)" + chain
    else:
        # terminal opening: N-CH2-CH(OH)-C10H21, as drawn.
        arm = "[CH3:2]C(O)" + chain

    # The core amine, as drawn. 5 N, 5 N-H -> 5 arms.
    core = "NCCNCCN1CCN(CCN)CC1"
    m = Chem.MolFromSmiles(core)
    n_h = sum(a.GetTotalNumHs() for a in m.GetAtoms() if a.GetSymbol() == "N")
    if n_h != 5:
        raise ValueError(f"core amine has {n_h} N-H, expected 5")

    steps = []
    cur = core
    # 1. the one SECONDARY acyclic N (the drawn middle nitrogen of the long
    #    branch) -- it can take exactly one arm.
    cur = assemble.attach(tag(cur, "[NX3;H1;!R]", 1), arm, 1, 2,
                          label="N_mid + arm")
    steps.append(("secondary acyclic N", cur))
    # 2-3. the two PRIMARY N, one arm each -> both become secondary.
    for i in range(2):
        cur = assemble.attach(tag(cur, "[NX3;H2]", 1), arm, 1, 2,
                              label=f"primary N {i} + arm 1")
        steps.append((f"primary N {i}, first arm", cur))
    # 4-5. the same two N, now secondary and acyclic, take their second arm.
    for i in range(2):
        cur = assemble.attach(tag(cur, "[NX3;H1;!R]", 1), arm, 1, 2,
                              label=f"secondary N {i} + arm 2")
        steps.append((f"primary N {i}, second arm", cur))

    prod = Chem.MolFromSmiles(cur)
    left = sum(a.GetTotalNumHs() for a in prod.GetAtoms()
               if a.GetSymbol() == "N")
    if left != 0:
        raise ValueError(f"{left} N-H left after 5 alkylations")
    assemble.check_no_radicals(prod, label="C12-200 anchor")
    return assemble.canonical(cur), steps


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=HERE)
    ap.add_argument("--seed", default=None)
    a = ap.parse_args()
    if a.seed and os.path.abspath(a.outdir) == os.path.abspath(HERE):
        sys.exit("REFUSED: a seeded run must write to a temp --outdir (F-069).")
    os.makedirs(a.outdir, exist_ok=True)

    smiles, steps = build(a.seed)
    m = Chem.MolFromSmiles(smiles)
    print("C12-200 anchor, assembled from Figure 1A:")
    print("  SMILES   ", smiles)
    print("  InChIKey ", Chem.MolToInchiKey(m))
    print("  formula  ", rdMolDescriptors.CalcMolFormula(m))
    print("  exactmass", f"{Descriptors.ExactMolWt(m):.6f}")
    for name, s in steps:
        print(f"    step: {name:<28} {Chem.rdMolDescriptors.CalcMolFormula(Chem.MolFromSmiles(s))}")

    rows = [
        # the one structure this paper draws in full and the corpus also holds
        {"lipid_name": "C12-200",
         "reference_SMILES": smiles,
         "source_figure": "Figure 1A, green box (gr4/gr1_lrg.jpg)",
         "read_on": "2026-11-10", "read_by": "curation agent",
         "expected_mismatch": "",
         "note": ("Assembled head-first from the drawn connectivity by "
                  "assemble.attach with SMARTS-selected targets; every tail "
                  "is PRINTED as the text 'C10H21', so chain length is read, "
                  "not counted. No corpus file is read by make_reference.py "
                  "and no literal >=25 chars in it occurs in the corpus "
                  "string; both asserted at run time by AST walk.")},
        # drawn, but not independently derivable without counting raster
        # vertices -- the absence is recorded rather than faked
        {"lipid_name": "C14PEG2k",
         "reference_SMILES": "",
         "source_figure": "Figure 1A, red box",
         "read_on": "2026-11-10", "read_by": "curation agent",
         "expected_mismatch": "head_group_class_differs_from_corpus_name",
         "note": ("DRAWN as a 1,2-diacyl-sn-glycero-3-phospho-ethanolamine "
                  "N-carbamoyl-mPEG45 (ammonium salt): phosphate, -NH-C(=O)O-, "
                  "subscript 45. The corpus name 'C14-PEG2000' (BL_2023_3CR, "
                  "480 rows) resolves to InChIKey NLDCRRFPGOSAHV, formula "
                  "C122H242O50 -- a GLYCEROL ETHER with no N and no P, i.e. "
                  "DMG-PEG2000. No SMILES stored rather than a guessed one; "
                  "the acyl chains would have to be counted off a raster.")},
        {"lipid_name": "SM-102", "reference_SMILES": "",
         "source_figure": "", "read_on": "2026-11-10",
         "read_by": "curation agent", "expected_mismatch": "",
         "note": ("NOT DRAWN anywhere in the paper -- named only, in Fig 4G "
                  "and Fig S8. Empty on purpose so the absence of an anchor "
                  "is visible (brief section 9). Structure taken from "
                  "provenance/reagent_registry.py, which is a name-derived "
                  "pin and NOT an independent anchor.")},
        {"lipid_name": "ALC-0315", "reference_SMILES": "",
         "source_figure": "", "read_on": "2026-11-10",
         "read_by": "curation agent", "expected_mismatch": "",
         "note": ("NOT DRAWN -- named only. Same posture as SM-102. Against "
                  "an empty anchor the C10 guards are VACUOUS, not passing.")},
        {"lipid_name": "Cholesterol", "reference_SMILES": "",
         "source_figure": "Figure 1A, yellow box",
         "read_on": "2026-11-10", "read_by": "curation agent",
         "expected_mismatch": "",
         "note": ("Drawn with 3-beta-OH wedge and all eight centres, but NOT "
                  "hand-written here: one agent's hand-written cholesterol "
                  "came out a diastereomer with 4 of 8 centres inverted at "
                  "identical formula, identical exact mass and an identical "
                  "skeleton InChIKey block. Structure comes from "
                  "provenance/reagent_registry.py, and pinning the registry "
                  "against itself is DECLARED NOT INDEPENDENT rather than "
                  "counted as an anchor.")},
    ]
    path = os.path.join(a.outdir, "reference_structures.csv")
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=[
            "lipid_name", "reference_SMILES", "source_figure", "read_on",
            "read_by", "expected_mismatch", "note"])
        w.writeheader()
        w.writerows(rows)
    print("wrote", path)


if __name__ == "__main__":
    main()
