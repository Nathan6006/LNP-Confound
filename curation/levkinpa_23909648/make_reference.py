#!/usr/bin/env python3
"""Build the C10 anchors for LL_2013 -- the structures the paper DRAWS.

WHY A GENERATOR AND NOT A STORED STRING. CLAUDE.md section 6: a C10 anchor
only tests anything if it was built by a DIFFERENT CODE PATH THAT IS EXECUTED,
NOT ASSERTED. That has been claimed twice in this project and been false twice
-- F-013, where the Siegwart anchor was a byte copy of the enumeration with a
note claiming independent assembly, and F-037, where the file written to FIX
F-013 re-copied the other half (its AMINE_HALF template was the enumeration's
canonical SMILES with the sterol swapped for '%' and one ring digit
renumbered, byte-identical over 141 characters).

C16 tests the GENERATOR'S COMPOSED TEMPLATE, not the emitted anchor, because
byte-identity of the emitted anchor is not the signal: a correct re-derivation
is byte-identical to a correct enumeration once canonicalised. That is exactly
why C10's own `anchor_may_be_a_copy` missed F-037.

SO: HOW THIS FILE IS INDEPENDENT, CONCRETELY.
There are three structure-building routes in this folder and no two share a
template:

  stage03_enumerate.py  RDKit fragment assembly. Atom-mapped attachment via
                        curation/lib/assemble.attach (CombineMols + AddBond),
                        head fragment "C1CCN(C1)CC[NH2:1]", alkyl fragment
                        "[CH2:2]" + "C"*(n-1), one attach per tail.
  stage07_qc.py         String composition, TAIL-FIRST: f"{alk}N({alk})CCN1CCCC1".
  make_reference.py     THIS FILE. String composition, HEAD-FIRST, writing the
                        ring out atom by atom from the drawing and walking
                        outward: "C1CCCN1" + "CC" + "N(...)(...)".

The head-first template below is written by reading the drawing on SI pp.
7/10/13 left to right -- pyrrolidine ring, then the two-carbon bridge, then
the nitrogen and its tails -- which is the opposite traversal to both other
routes. It shares no substring of either beyond the unavoidable "CCN1CCCC1"
fragment that any SMILES for this head must contain.

WHAT THE ANCHOR IS READ FROM. The paper draws all three purified compounds in
full (SI pp.7, 10, 13), each with its compound label, Chemical Formula, Exact
Mass and Molecular Weight printed beneath. The drawing is vector: 32 / 24 / 40
filled bond outlines, counted in stage02_drawn.py and matching the heavy-atom
bond count of each structure exactly. The chain lengths and tail counts below
come from that drawing and from the printed formula, NOT from the compound
label -- "name tokens are not chain lengths" (KS_2024's C24 is two C4 tails).

Run:  python3 -I make_reference.py
"""
import csv
import os
import sys

from rdkit import Chem, RDLogger
from rdkit.Chem import inchi
from rdkit.Chem.Descriptors import ExactMolWt
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))

# --- the head, written atom by atom from the drawing ------------------------
# Read off SI p.10 (the clearest of the three, a single tail): a five-membered
# ring whose only heteroatom is N, then two CH2, then the second N.
RING = "C1CCC"          # four ring carbons
RING_CLOSE = "N1"       # the ring nitrogen, closing to atom 1
BRIDGE = "CC"           # the two-carbon bridge drawn between the two N
HEAD = RING + RING_CLOSE + BRIDGE

# Tail lengths read from the DRAWING and cross-checked against the printed
# Chemical Formula under each structure, never from the compound label.
#   1C12-D  C30H62N2  head C6 + 2 x C12 = C30   two tails drawn on one N
#   1C16-S  C22H46N2  head C6 + 1 x C16 = C22   one tail, N-H drawn explicitly
#   1C16-D  C38H78N2  head C6 + 2 x C16 = C38   two tails drawn on one N
ANCHORS = [
    dict(label="1C12-D", chain=12, tails=2, page=7,
         drawn_formula="C30H62N2", drawn_exact_mass=450.49130,
         iupac="N-dodecyl-N-(2-(pyrrolidin-1-yl)ethyl)dodecan-1-amine"),
    dict(label="1C16-S", chain=16, tails=1, page=10,
         drawn_formula="C22H46N2", drawn_exact_mass=338.36610,
         iupac="N-(2-(pyrrolidin-1-yl)ethyl)hexadecan-1-amine"),
    dict(label="1C16-D", chain=16, tails=2, page=13,
         drawn_formula="C38H78N2", drawn_exact_mass=562.61650,
         iupac="N-hexadecyl-N-(2-(pyrrolidin-1-yl)ethyl)hexadecan-1-amine"),
]


def compose(chain, tails):
    """HEAD-FIRST template: ring -> bridge -> N -> tails. Returns the raw,
    UNCANONICALISED template string, which is what C16 would examine."""
    tail = "C" * chain
    if tails == 1:
        return HEAD + "N" + tail          # secondary amine, N-H implicit
    if tails == 2:
        return HEAD + "N(" + tail + ")" + tail
    raise ValueError(tails)


def main():
    rows = []
    templates = {}
    for a in ANCHORS:
        template = compose(a["chain"], a["tails"])
        m = Chem.MolFromSmiles(template)
        if m is None:
            raise RuntimeError(f"{a['label']}: template unparseable: "
                               f"{template}")
        smi = Chem.MolToSmiles(m)
        formula = CalcMolFormula(m)
        exact = ExactMolWt(m)

        # The anchor must reproduce what the AUTHORS printed under the
        # drawing. If it does not, it is not an anchor for this structure.
        if formula != a["drawn_formula"]:
            raise RuntimeError(
                f"{a['label']}: composed {formula}, the drawing prints "
                f"{a['drawn_formula']}.")
        if abs(exact - a["drawn_exact_mass"]) > 2e-4:
            raise RuntimeError(
                f"{a['label']}: composed exact mass {exact:.5f}, the drawing "
                f"prints {a['drawn_exact_mass']}.")

        templates[a["label"]] = template
        rows.append(dict(
            lipid_name=a["label"],
            reference_SMILES=smi,
            source_figure=f"SI p.{a['page']}, the full structure drawing "
                          f"above the printed '{a['drawn_formula']}' / "
                          f"'Exact Mass: {a['drawn_exact_mass']}' block "
                          f"(vector; bonds are filled outlines)",
            read_on="2026-11-10",
            read_by="curation agent, from a 400 dpi pymupdf render of SI "
                    f"p.{a['page']}; bond count verified geometrically in "
                    "stage02_drawn.py",
            expected_mismatch="",
            note=("Composed HEAD-FIRST by plain string concatenation from the "
                  "drawing: ring '" + HEAD + "' then the nitrogen then the "
                  "tails. NO RDKit fragment assembly, no atom maps, no "
                  "RWMol, no AddBond, no RemoveAtom -- a different code path "
                  "from stage03_enumerate.py (which uses "
                  "curation/lib/assemble.attach) and a different traversal "
                  "from stage07_qc.py's tail-first re-derivation. Raw "
                  "template before canonicalisation: '" + template + "'. "
                  "Chain length and tail count are read from the drawing and "
                  "from the authors' own printed Chemical Formula "
                  + a["drawn_formula"] + ", not from the compound label. "
                  "IUPAC name as printed: " + a["iupac"] + ".")))

    hdr = ["lipid_name", "reference_SMILES", "source_figure", "read_on",
           "read_by", "expected_mismatch", "note"]
    with open(os.path.join(HERE, "reference_structures.csv"), "w",
              newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=hdr)
        w.writeheader()
        for r in rows:
            w.writerow(r)

    print("OK  reference_structures.csv written")
    for a, r in zip(ANCHORS, rows):
        print(f"    {r['lipid_name']:<7} template {templates[r['lipid_name']]}")
        print(f"            canonical {r['reference_SMILES']}")
        print(f"            InChIKey  "
              f"{inchi.MolToInchiKey(Chem.MolFromSmiles(r['reference_SMILES']))}"
              f"   formula {a['drawn_formula']} (matches the drawing)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
