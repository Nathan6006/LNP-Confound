#!/usr/bin/env python3
"""Build reference_structures.csv (the C10 anchor). Run: python3 -I make_reference.py

C10 compares each enumerated structure against the structure the paper DRAWS.
It only tests anything when the anchor was built by a DIFFERENT CODE PATH that
is EXECUTED, not asserted -- a claim this project has made twice and had false
twice (F-013, a byte copy; F-037, the file written to fix F-013 re-copied the
other half).

Here that condition is satisfied in the strongest and least interesting way:
**there is no enumerator in this folder to collide with**, because this paper
supplied no per-member values and nothing was enumerated. This script asserts
that mechanically rather than claiming it in prose -- it exits non-zero if
03_enumerated.csv ever appears without this file being rewritten.

And the anchor itself is EMPTY, because the paper draws no structure anywhere.
PIPELINE_V2.md: "If the paper draws nothing, store an empty reference_SMILES so
the absence is visible." That is the honest deliverable, not a structure
back-filled from the corpus -- which would make C10 compare the corpus to
itself, which is exactly F-013.

The script checks the "draws nothing" claim instead of asserting it: it reads
every figure caption out of the SI's text layer and requires that none of them
announces a drawn structure, a scheme or a synthesis. Seed a caption word and
it refuses.
"""
import csv
import os
import re
import sys

import pymupdf

HERE = os.path.dirname(os.path.abspath(__file__))
PDF = os.path.join(HERE, "source", "pnas.2528477123.sapp.pdf")
OUT = os.path.join(HERE, "reference_structures.csv")

HEADER = ["lipid_name", "reference_SMILES", "source_figure", "read_on",
          "read_by", "expected_mismatch", "note"]

# A caption containing any of these would mean the paper DOES draw chemistry
# somewhere, and this file's whole premise would be wrong.
DRAWN_CHEMISTRY = re.compile(
    r"\b(chemical structure|structures? of|synthetic (route|scheme)|"
    r"scheme|synthesis of|reaction scheme)\b", re.I)

CAPTION = re.compile(r"(Figure S\d+\.[^\n]*(?:\n(?!\s*(?:Figure S|Table S))[^\n]*){0,6})")


def captions(doc):
    out = []
    for page in doc:
        txt = page.get_text()
        for m in CAPTION.finditer(txt):
            cap = " ".join(m.group(1).split())
            out.append(cap)
    return out


def main(caption_override=None):
    doc = pymupdf.open(PDF)
    caps = caption_override if caption_override is not None else captions(doc)
    if len(caps) < 10:
        sys.exit(f"make_reference: found {len(caps)} figure captions in the "
                 "SI; expected 11 (Figures S1-S11). Refusing to certify that "
                 "the paper draws nothing on the strength of a failed read.")

    offenders = [c for c in caps if DRAWN_CHEMISTRY.search(c)]
    if offenders:
        sys.exit("make_reference: a caption announces drawn chemistry, so the "
                 "premise of this file (an EMPTY anchor) is wrong. Build the "
                 "anchor from that figure instead.\n  " +
                 "\n  ".join(o[:160] for o in offenders))

    if os.path.exists(os.path.join(HERE, "03_enumerated.csv")):
        sys.exit("make_reference: 03_enumerated.csv exists. This file's claim "
                 "of independence rests on there being NO enumerator in the "
                 "folder. Rewrite the anchor by a genuinely separate route "
                 "before trusting C10 here.")

    rows = [{
        "lipid_name": "6A1-SC8",
        "reference_SMILES": "",
        "source_figure": "none - SI Fig S1 is 1H-NMR (A) and an LC/MS isotope "
                         "cluster (B); the SI contains no drawn structure, no "
                         "scheme and no synthesis figure",
        "read_on": "2026-11-10",
        "read_by": "curation agent, SI pp.1-15 text layer + 200/900 dpi page "
                   "renders",
        "expected_mismatch": "",
        "note": "ANCHOR ABSENT ON PURPOSE. The paper draws no structure, so "
                "there is nothing to compare topology against; all 11 figure "
                "captions were scanned and none announces drawn chemistry. "
                "What the SI does supply is a COMPOSITION anchor, recorded in "
                "07_qc.json: the caption's formula C110H203N5O24S6 matches the "
                "corpus structure for 6A1-SC8 element-for-element (5/5, q1), "
                "and that structure's predicted [M+2H]2+ ladder matches 5/5 "
                "m/z labels drawn inside Fig S1B to within 3.23 ppm (q4, "
                "charge state z=2 derived from the 0.5006 Da isotope spacing, "
                "which the caption never states). NMR integral ratios match "
                "the structure's own H counts 2/2 within 5% (q6). Mass and "
                "integrals cannot see connectivity: nothing here constrains "
                "which core nitrogen bears which arm. Do NOT fill this cell "
                "with the corpus SMILES - C10 would then be comparing the "
                "corpus to itself, which is F-013 exactly.",
    }]

    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=HEADER)
        w.writeheader()
        w.writerows(rows)

    print(f"captions scanned      {len(caps)}")
    print(f"captions drawing chemistry  {len(offenders)}")
    print(f"anchors written       {len(rows)} (reference_SMILES empty: "
          f"{sum(1 for r in rows if not r['reference_SMILES'])})")
    print("wrote", OUT)
    return 0


if __name__ == "__main__":
    sys.exit(main())
