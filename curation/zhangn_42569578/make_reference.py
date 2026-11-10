#!/usr/bin/env python3
"""MY_2026 (PMID 42569578) -- the C10 anchor, and the guards that police it.

THE ANCHOR IS EMPTY, ON PURPOSE. This paper draws NO chemical structure
anywhere: Fig 1 and the graphical abstract render the four lipid classes as
coloured cartoon blobs, and Figs 2-7 are DLS scatter, Spearman/SHAP panels, bar
charts and IVIS organ images. All four ionizable lipids are catalogue reagents
bought from Avanti (SM-102, ALC-0315, DLin-MC3-DMA) and MedChemExpress
(CKK-E12). So `reference_SMILES` is stored EMPTY for all four, which is the
posture CLAUDE.md asks for: store the empty row so the ABSENCE of an anchor is
visible rather than invisible.

WHAT THAT MEANS FOR THE GUARDS. Both guards below are implemented and are
demonstrated to fire by test_seeds.py. But against an empty anchor they are
**VACUOUS, not passing**: there is no stored structure for a copy to corrupt.
Saying "C10 clean" here would be the F-013/F-037 mistake in its purest form --
comparing a value to itself and reporting a match.

WHAT IS ACTUALLY CHECKED instead, and is not vacuous: three of the four
structures come from provenance/reagent_registry.py, which writes them from the
chemical NAME, and each is then confirmed against the live corpus by FULL
InChIKey (stage03_corpus_overlap.py). CKK-E12 has no registry entry and is
taken FROM the corpus, so its corpus "match" is circular and is recorded as
NOT INDEPENDENT.

Run: python3 -I make_reference.py [--outdir DIR]
A seed harness MUST pass --outdir (F-069): a seeded run that writes into the
curation folder leaves a poisoned artifact behind, and a non-zero exit conceals
the write.
"""
import argparse
import ast
import csv
import os
import sys

from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, rdMolDescriptors

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
HEADER = ["lipid_name", "reference_SMILES", "source_figure", "read_on",
          "read_by", "expected_mismatch", "note"]

NOTHING_DRAWN = ("The paper draws no chemical structure: Fig 1 and the "
                 "graphical abstract show cartoon lipid blobs, Figs 2-7 are "
                 "DLS/SHAP/bar/IVIS panels. Catalogue reagent, no anchor "
                 "available. Structure taken from %s.")

ANCHORS = [
    ("ALC-0315", "provenance/reagent_registry.py (written from the IUPAC name)"),
    ("SM-102", "provenance/reagent_registry.py (written from the IUPAC name)"),
    ("MC-3 (DLin-MC3-DMA)",
     "provenance/reagent_registry.py (written from the IUPAC name)"),
    ("CKK-E12",
     "the live corpus via curation/corpus_components.smi -- NOT INDEPENDENT, "
     "so its corpus InChIKey agreement is circular and is not counted"),
]

# Modules/attributes this file must never reach for. Walking the AST (and only
# Import/ImportFrom/Call/Attribute nodes) means the ban list itself -- a plain
# list of strings in this module -- cannot match itself. A substring scan over
# the source WOULD flag this very docstring, which is how one sibling's guard
# failed on its own module.
BANNED_MODULES = {"stage01_stage", "corpus_components", "reagent_registry",
                  "pandas"}
BANNED_ATTRS = {"read_csv", "read_excel"}


def guard_ast(path):
    """Guard A -- independence, by walking THIS module's AST at run time."""
    tree = ast.parse(open(path, encoding="utf8").read(), filename=path)
    hits = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for al in node.names:
                if al.name.split(".")[0] in BANNED_MODULES:
                    hits.append(f"import {al.name} (line {node.lineno})")
        elif isinstance(node, ast.ImportFrom):
            if (node.module or "").split(".")[0] in BANNED_MODULES:
                hits.append(f"from {node.module} import (line {node.lineno})")
        elif isinstance(node, ast.Attribute):
            if node.attr in BANNED_ATTRS:
                hits.append(f".{node.attr} (line {node.lineno})")
        elif isinstance(node, ast.Call):
            f = node.func
            if isinstance(f, ast.Name) and f.id in BANNED_ATTRS:
                hits.append(f"{f.id}() (line {node.lineno})")
    return hits


def guard_template_provenance(path, enumerated):
    """Guard B -- no long literal here may be a substring of an emitted SMILES.

    Byte-identity of the EMITTED anchor is not the signal and testing for it is
    a trap (F-037): a correct re-derivation is byte-identical to a correct
    enumeration once canonicalised. What is diagnostic is a long literal
    sitting in the GENERATOR that came out of the enumeration.
    """
    tree = ast.parse(open(path, encoding="utf8").read(), filename=path)
    hits = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Constant) and isinstance(node.value, str):
            lit = node.value
            if len(lit) >= 25:
                for nm, smi in enumerated:
                    if lit in smi:
                        hits.append({"line": node.lineno, "len": len(lit),
                                     "inside": nm, "literal": lit[:60]})
    return hits


def load_enumerated(folder):
    p = os.path.join(folder, "03_enumerated.csv")
    if not os.path.exists(p):
        return []
    with open(p, newline="", encoding="utf8") as f:
        return [(r["IL_name_paper"], r["SMILES"]) for r in csv.DictReader(f)]


def connectivity_probe():
    """Prove the gate tests CONNECTIVITY, not composition.

    Builds a constitutional isomer of SM-102 by moving one methyl along a
    chain: identical molecular formula, identical exact mass to 6 dp, and a
    DIFFERENT InChIKey. If a check cannot separate these two it is not testing
    structure at all -- it is testing arithmetic.
    """
    a = "CCCCCCCCC(CCCCCCCC)OC(=O)CCCCCCCN(CCO)CCCCCC(=O)OCCCCCCCCCCC"
    b = "CCCCCCCCCC(CCCCCCC)OC(=O)CCCCCCCN(CCO)CCCCCC(=O)OCCCCCCCCCCC"
    ma, mb = Chem.MolFromSmiles(a), Chem.MolFromSmiles(b)
    return {
        "formula_a": rdMolDescriptors.CalcMolFormula(ma),
        "formula_b": rdMolDescriptors.CalcMolFormula(mb),
        "formula_identical":
            rdMolDescriptors.CalcMolFormula(ma) ==
            rdMolDescriptors.CalcMolFormula(mb),
        "exact_mass_a": round(Descriptors.ExactMolWt(ma), 6),
        "exact_mass_b": round(Descriptors.ExactMolWt(mb), 6),
        "exact_mass_delta_Da":
            round(abs(Descriptors.ExactMolWt(ma) -
                      Descriptors.ExactMolWt(mb)), 10),
        "inchikey_a": Chem.MolToInchiKey(ma),
        "inchikey_b": Chem.MolToInchiKey(mb),
        "inchikey_distinguishes":
            Chem.MolToInchiKey(ma) != Chem.MolToInchiKey(mb),
        "skeleton_block_distinguishes":
            Chem.MolToInchiKey(ma).split("-")[0] !=
            Chem.MolToInchiKey(mb).split("-")[0],
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=HERE)
    ap.add_argument("--self", default=os.path.abspath(__file__),
                    help="module to police (the seed harness points this at "
                         "a mutated COPY, never at the real file)")
    a = ap.parse_args()
    out = os.path.abspath(a.outdir)
    os.makedirs(out, exist_ok=True)

    enum = load_enumerated(HERE)
    ga = guard_ast(a.self)
    gb = guard_template_provenance(a.self, enum)
    probe = connectivity_probe()

    print("GUARD_A_ast_independence_violations=%d" % len(ga))
    for h in ga:
        print("   AST:", h)
    print("GUARD_B_template_provenance_violations=%d" % len(gb))
    for h in gb:
        print("   TEMPLATE:", h)
    print("CONNECTIVITY_formula_identical=%s exact_mass_delta_Da=%s "
          "inchikey_distinguishes=%s"
          % (probe["formula_identical"], probe["exact_mass_delta_Da"],
             probe["inchikey_distinguishes"]))
    print("ANCHORS_WITH_A_STRUCTURE=0 ANCHORS_EMPTY=%d" % len(ANCHORS))
    print("GUARDS_ARE_VACUOUS_AGAINST_AN_EMPTY_ANCHOR=True")

    if ga or gb:
        print("REFUSING to write reference_structures.csv: a guard fired.")
        return 1
    if not (probe["formula_identical"] and probe["inchikey_distinguishes"]
            and probe["exact_mass_delta_Da"] == 0.0):
        print("REFUSING: the connectivity probe did not behave as stated.")
        return 1

    with open(os.path.join(out, "reference_structures.csv"), "w", newline="",
              encoding="utf8") as f:
        w = csv.writer(f)
        w.writerow(HEADER)
        for nm, src in ANCHORS:
            w.writerow([nm, "", "none - no structure is drawn in this paper",
                        "2026-11-10", "curation agent (zhangn_42569578)", "",
                        NOTHING_DRAWN % src])
    print("wrote", os.path.join(out, "reference_structures.csv"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
