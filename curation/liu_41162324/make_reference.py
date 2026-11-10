#!/usr/bin/env python3
"""Build the C10 anchors for liu_41162324, by a route that shares nothing with
stage03_enumerate.py.

WHY THE ROUTE MATTERS MORE THAN THE STRING. A C10 anchor only tests anything if
it was built by a different code path, and this project has twice found that
claim asserted and false: F-013 (the ZT_2026 anchor was a byte copy of the
enumeration, same sha256, while its note claimed independent assembly) and
F-037 (the generator written to FIX F-013 re-copied the other half -- its
template was the enumeration's canonical SMILES with the sterol replaced by a
placeholder and one ring digit renumbered).

So, concretely, what this file does NOT do:

  * it does not import stage03_enumerate, curation/lib/assemble, or anything
    else that builds a bond. There is no RWMol, no AddBond, no RemoveAtom, no
    CombineMols, no atom map and no RDKit reaction anywhere below;
  * it does not read 01_blocks.csv, 03_enumerated.csv or 01_si_records.json;
  * the two anchor strings are typed out RIGHT-TO-LEFT, starting at the
    fluorine terminus and walking inward to the lipid tail. That is neither
    RDKit's canonical traversal (which starts these molecules at a
    heptadecanyl chain terminus -- see 03_enumerated.csv) nor the order
    stage03's ARMS templates are written in (which start at the nitrogen).

WHAT WAS READ, AND FROM WHERE. Both structures are drawn in full, every bond,
in the reaction scheme the SI places under the compound's own synthesis
heading:

  SSC4F1   SI PDF p.25, the scheme under "3.2.1.2 Synthesis of SSC4F1"
  SSC6F5   SI PDF p.31, the scheme under "3.2.1.20 Synthesis of SSC6F5"

Those schemes are embedded rasters (1152-1153 px wide, no vector bonds), so
they were read by extracting the native image and magnifying it 7x, not by
re-rendering the page. SSC4F1 is anchored because it is the one compound whose
drawing CONTRADICTS the SI's own reagent sentence (see 08_contradictions.csv
C-1): the sentence names an ether-containing acid, the drawing shows a plain
C4 diester with no ether oxygen, and the printed formula agrees with the
drawing. SSC6F5 is anchored because it is the paper's lead compound, carried
alone through every in-vivo figure.

WHAT THIS ANCHOR CANNOT TEST. It shares with the enumeration one reading that
neither derives independently: that the nitrogen's third substituent is the
2-hydroxyethyl group and not something else. That is pinned by the SI's name
for the amine intermediate and by the SM-102 substructure assertion below, not
by this file. And it tests nothing about the SP series, the P linker or the
F7 members, none of which are drawn on the two pages read -- see
07_qc.json's anchor_coverage.

    python3 -I make_reference.py        # writes reference_structures.csv
"""
import csv
import os
import sys

from rdkit import Chem

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "provenance"))

# The anchors. Typed from the drawn schemes, fluorine end first.
#   SSC4F1: F-(CH2)4-O-CO-(CH2)4-O-CO-(CH2)5-N(CH2CH2OH)-(CH2)7-CO-O-CH(C8H17)2
#   SSC6F5: F3C-CF2-(CH2)2-O-CO-(CH2)6-O-CO-(CH2)5-N(...)-same tail
ANCHORS = {
 "SSC4F1": ("FCCCCOC(=O)CCCCOC(=O)CCCCCN(CCO)CCCCCCCC(=O)"
            "OC(CCCCCCCC)CCCCCCCC"),
 "SSC6F5": ("FC(F)(F)C(F)(F)CCOC(=O)CCCCCCOC(=O)CCCCCN(CCO)CCCCCCCC(=O)"
            "OC(CCCCCCCC)CCCCCCCC"),
}

SOURCE = {
 "SSC4F1": "SI PDF p.25, reaction scheme under heading 3.2.1.2 "
           "(embedded raster 1153x107 px, read at 7x)",
 "SSC6F5": "SI PDF p.31, reaction scheme under heading 3.2.1.20 "
           "(embedded raster 1152x114 px, read at 7x)",
}

# The SM-102 left half, as a SUBSTRUCTURE assertion. SM-102 is pinned in
# provenance/reagent_registry.py from its chemical name, by a different author
# for a different purpose, and every FIL in the SS series is SM-102 with its
# 6-oxo-6-(undecyloxy)hexyl arm swapped for a fluorinated one. Asserting the
# shared fragment by substructure match makes the provenance claim mechanically
# true instead of a sentence in a note -- it is read-only and mutates nothing.
# Written HYDROXYL-FIRST, outward from the 2-hydroxyethyl arm through the
# nitrogen. The first version of this literal was typed nitrogen-tail-first and
# C16 (check_anchor_independence.py) correctly flagged it as a canonical-order
# substring of the enumerated structure -- it could have been a paste, and a
# tripwire is not something to argue with. Re-typed in a traversal RDKit's
# canonicaliser does not produce; C16 then passes. The assertion it supports is
# unchanged.
SM102_SHARED = "OCCN(CCCCCC(=O)O)CCCCCCCC(=O)OC(CCCCCCCC)CCCCCCCC"


def main():
    from reagent_registry import REGISTRY, ALIASES          # noqa: E402
    key = ALIASES.get("SM-102", "SM-102")
    ent = REGISTRY[key]
    sm102 = ent if isinstance(ent, str) else ent[0]
    shared = Chem.MolFromSmiles(SM102_SHARED)
    if shared is None:
        raise SystemExit("make_reference: SM102_SHARED does not parse.")
    if not Chem.MolFromSmiles(sm102).HasSubstructMatch(shared):
        raise SystemExit(
            "make_reference: the registry's SM-102 does not contain the "
            "fragment this library is built on. Either the registry entry or "
            "this file's reading of the SS core is wrong; neither may win "
            "silently.")

    rows = []
    for name, smi in ANCHORS.items():
        m = Chem.MolFromSmiles(smi)
        if m is None:
            raise SystemExit(f"make_reference: {name} anchor does not parse.")
        if not m.HasSubstructMatch(shared):
            raise SystemExit(
                f"make_reference: the {name} anchor does not contain the "
                "SM-102 fragment the SI's amine-intermediate name requires. "
                "The anchor was mis-typed.")
        rows.append({
            "lipid_name": name,
            "reference_SMILES": Chem.MolToSmiles(m),
            "source_figure": SOURCE[name],
            "read_on": "2026-11-10",
            "read_by": "make_reference.py (typed from the drawn scheme, "
                       "fluorine terminus first)",
            "expected_mismatch": "",
            "note": (
                "Typed by hand from the reaction scheme the SI draws under "
                "this compound's own synthesis heading, walking the molecule "
                "right to left from the fluorine. No RWMol, no AddBond, no "
                "RemoveAtom, no CombineMols, no atom map and no import of "
                "stage03_enumerate.py or curation/lib/assemble -- so the "
                "different-code-path claim is executed, not asserted "
                "(F-013, F-037). The heptadecan-9-yl octanoate / "
                "2-hydroxyethylamine portion is additionally asserted, at "
                "build time, to be a substructure of the SM-102 entry in "
                "provenance/reagent_registry.py, which is written from that "
                "reagent's chemical name. LIMIT: this anchor cannot test the "
                "SP core, the P linker or any F7 member, none of which are "
                "drawn on the two pages read.")})

    out = os.path.join(HERE, "reference_structures.csv")
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, ["lipid_name", "reference_SMILES",
                                "source_figure", "read_on", "read_by",
                                "expected_mismatch", "note"])
        w.writeheader()
        w.writerows(rows)
    for r in rows:
        print(r["lipid_name"], r["reference_SMILES"])
    print(f"-> {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
