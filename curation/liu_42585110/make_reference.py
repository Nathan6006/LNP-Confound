#!/usr/bin/env python3
"""BUILD the C10 anchor for _scaffold_tmp_kj. Do not store a string -- derive it.

    python3 -I make_reference.py

C10 compares each enumerated structure against the structure the paper DRAWS.
It is the only Tier-1 check that tests topology against external ground truth.
It only tests anything if the anchor was built by a DIFFERENT CODE PATH than the
enumerator.

That independence has been asserted twice and false twice:

  F-013  the stored anchor was a BYTE COPY of the enumeration -- same sha256,
         223 bytes -- while its note claimed independent fragment assembly
         "through the C3 oxygen (3-beta)" and the string was 3-alpha. C10
         reported `exact_match` while comparing a value to itself, on the
         flagship curation, for the one molecule it exists to protect. 537 rows
         of epicholesterol.
  F-037  the file written in round 1 TO FIX THAT had an AMINE_HALF template
         that WAS the enumeration's canonical SMILES with the sterol swapped
         for `%` and one ring digit renumbered -- byte-identical over 141
         characters. So the repair cured the half the bug was in and re-copied
         the other half.

PREFER A ROUTE THAT SHARES AS LITTLE AS POSSIBLE with the enumerator. The
working replacement uses plain string substitution with no RWMol / AddBond /
RemoveAtom anywhere, and takes the sterol VERBATIM from reagent_registry.json
so the provenance claim is mechanically true.

If the paper draws no full structure, write a row with an empty
reference_SMILES and a note saying so, so the ABSENCE of an anchor is visible.
"""
import csv
import os
import sys

LIB = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib")
sys.path.insert(0, LIB)

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   "reference_structures.csv")


def build():
    """Return [(compound_id, reference_SMILES, source, note)].

    `source` names the FIGURE OR SCHEME the structure was read from, e.g.
    "Fig 2a, lead compound, read head-first". `note` says what route built it.
    """
    raise SystemExit(
        "make_reference.py is still the generated stub. Read the structure the "
        "paper DRAWS -- usually the lead compound or the final product of a "
        "'precise synthesis' scheme -- and build it here by a route that "
        "shares as little as possible with the enumerator.")


if __name__ == "__main__":
    rows = build()
    with open(OUT, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["lipid_name", "reference_SMILES", "source_figure",
                    "read_on", "read_by", "expected_mismatch", "note"])
        w.writerows(rows)
    print(f"wrote {len(rows)} anchor(s) to {OUT}")
