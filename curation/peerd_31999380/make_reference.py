#!/usr/bin/env python3
"""Build the C10 reference anchors for peerd_31999380 INDEPENDENTLY.

RUN:  python3 -I make_reference.py      # writes reference_structures.csv

WHY THIS FILE EXISTS, AND WHAT MAKES IT INDEPENDENT
---------------------------------------------------
A C10 anchor tests nothing unless it was built by a different code path than
the enumerator. That has been claimed twice in this project and been false
twice -- F-013 (the stored anchor was a byte copy of the enumeration, same
sha256, while its note claimed independent fragment assembly) and F-037 (the
file written to FIX F-013 re-copied the other half: its AMINE_HALF was the
enumeration's canonical SMILES with the sterol swapped for '%' and a ring digit
renumbered).

So independence here is a property of PROVENANCE, not of bytes:

  1. THIS MODULE IMPORTS NO ASSEMBLER. There is no `import assemble`, no
     RDKit RWMol, no AddBond, no RemoveAtom, no attach(), no esterify().
     The only RDKit calls are MolFromSmiles / MolToSmiles / MolToInchiKey --
     parse and canonicalise, never construct.
  2. EVERY ANCHOR BELOW IS A SOURCE STRING LITERAL, typed by hand from the
     SI's own full IUPAC name for that compound and the structure the SI
     DRAWS beside it. No anchor is derived from 03_enumerated.csv; this file
     never reads that column before writing its own.
  3. THE TRAVERSAL IS REVERSED. 02_enumerate.py builds tail -> scaffold ->
     head. Every literal here is written HEAD FIRST: the dimethylamino or
     hydrazine nitrogen is the first atom typed, and the chain runs outward to
     the tails. A copy would have inherited the enumerator's traversal.

BYTE-IDENTITY IS NOT THE SIGNAL, AND TESTING FOR IT IS A TRAP. CLAUDE.md
section 6 is explicit: a CORRECT re-derivation is byte-identical to a CORRECT
enumeration once both are canonicalised, which is exactly why C10's own
`anchor_may_be_a_copy` missed F-037. If the keys below match, that is the
anchor doing its job, not evidence that it was copied. The evidence about
copying is points 1-3, which a reviewer checks by reading this file.

WHAT THE ANCHOR CAN AND CANNOT TEST. It tests CONNECTIVITY and DOUBLE-BOND
GEOMETRY against the SI's drawn structure and written name. It cannot test
double-bond POSITION independently, because the position comes from the same
IUPAC name in both files -- see DECISIONS.md D-2.
"""
import csv
import os
import sys

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")
HERE = os.path.dirname(os.path.abspath(__file__))

READ_ON = "2026-11-10"
READ_BY = "curation agent, from the rendered SI pages at 400 dpi plus the "\
          "SI's own IUPAC names"

# --------------------------------------------------------------------------
# THE ANCHORS. Hand-typed, head first. Each is accompanied by the SI's own
# name for the compound and the page that both names and draws it.
# --------------------------------------------------------------------------
ANCHORS = {
    # (2-(2-(dimethylamino)ethyl)hydrazine-1,1-diyl)bis(heptane-7,1-diyl)
    # bis(decanoate).  Me2N-CH2CH2-NH-N(< two heptyl decanoate arms).
    "Lipid 3": {
        "smiles": "CN(C)CCNN(CCCCCCCOC(=O)CCCCCCCCC)CCCCCCCOC(=O)CCCCCCCCC",
        "figure": "SI p.7, scheme 'Synthesis of Lipid 11 & 3' + full IUPAC "
                  "name in the paragraph heading",
        "note": "Head first: dimethylamino N, then the hydrazine N-N, then "
                "each heptane-7,1-diyl out to its decanoate. The enumerator "
                "builds the opposite way, decanoate inward.",
    },
    # Di((Z)-non-2-en-1-yl) 6,6'-(2-(4-(dimethylamino)butanoyl)hydrazine-1,1-
    # diyl)dihexanoate.  Carries the cis geometry, so it tests more than
    # connectivity.
    "Lipid 4": {
        "smiles": r"CN(C)CCCC(=O)NN(CCCCCC(=O)OC/C=C\CCCCCC)"
                  r"CCCCCC(=O)OC/C=C\CCCCCC",
        "figure": "SI p.12, drawn structure + full IUPAC name in the "
                  "paragraph heading",
        "note": "Chosen deliberately as one of the two anchors carrying C=C "
                "geometry. Contradiction C-05 is about this compound: its "
                "printed olefinic integral says four C=C, its two precursors "
                "and its mass say two, and this anchor encodes two.",
    },
    # di(tridecan-7-yl) 10,10'-(((4-(dimethylamino)butanoyl)oxy)azanediyl)
    # bis(decanoate).  The hydroxylamine O-ester family.
    "Lipid 7": {
        "smiles": "CN(C)CCCC(=O)ON(CCCCCCCCCC(=O)OC(CCCCCC)CCCCCC)"
                  "CCCCCCCCCC(=O)OC(CCCCCC)CCCCCC",
        "figure": "SI p.23, drawn structure + full IUPAC name in the "
                  "paragraph heading",
        "note": "The N-O bond is the feature worth anchoring: Lipid 7 and "
                "Lipid 13 differ ONLY in whether the head ester hangs off the "
                "nitrogen's own oxygen (7) or off an N-CH2CH2-O arm (13), "
                "which no mass check can distinguish -- both are C52/C54 "
                "homologues with different formulae but identical topology "
                "elsewhere.",
    },
    # ((2-((4-(dimethylamino)butanoyl)oxy)ethyl)azanediyl)bis(octane-8,1-diyl)
    # bis(2-hexyldecanoate).  The branched-tail family.
    "Lipid 14": {
        "smiles": "CN(C)CCCC(=O)OCCN(CCCCCCCCOC(=O)C(CCCCCC)CCCCCCCC)"
                  "CCCCCCCCOC(=O)C(CCCCCC)CCCCCCCC",
        "figure": "SI p.36, drawn structure + full IUPAC name in the "
                  "paragraph heading",
        "note": "Anchors the branch point: 2-hexyldecanoate puts a hexyl on "
                "C2 of a decanoate. A straight C16 acid has the same formula "
                "and the same exact mass, so ONLY a topology comparison can "
                "tell them apart.",
    },
}

# Lipids with NO anchor. Stored as empty rows on purpose so the GAP IS
# VISIBLE rather than absent -- CLAUDE.md section 5.
NO_ANCHOR = {
    "Lipid 1": "The SI draws the scheme but never writes a full IUPAC name "
               "for Lipid 1; its structure is stated only as 'N,N-dilinoleyl "
               "hydrazide' acylated with dimethylaminobutyric acid.",
    "Lipid 2": "Same as Lipid 1: named only through its two reagents.",
    "Lipid 5": "Named in full, but it shares the 2-hexyldecanoate branch with "
               "Lipid 14 and the hydrazine core with Lipid 3, both anchored; "
               "no feature of Lipid 5 is unanchored.",
    "Lipid 6": "Named only as HyAM-1 O-acylated; the N-O motif is anchored "
               "by Lipid 7.",
    "Lipid 8": "Named only as EA-1 O-acylated; the N-CH2CH2-O motif is "
               "anchored by Lipid 14.",
    "Lipid 9": "Named only through its route (EA-OMs + the alkoxide). The "
               "ETHER linkage is unique to Lipid 9 and is therefore the one "
               "connectivity in this paper with NO anchor -- see "
               "DECISIONS.md D-3.",
    "Lipid 10": "Named only as EA-1 O-acylated with the piperazine acid.",
    "Lipid 11": "Named in full; the motif is anchored by Lipid 3 (heptyl "
                "decanoate arms) and Lipid 14 (the N-CH2CH2-O head ester).",
    "Lipid 12": "Named in full; differs from the Lipid 4 anchor only in "
                "double-bond POSITION (non-3-enyl vs non-2-enyl), which no "
                "anchor can test independently -- DECISIONS.md D-2.",
    "Lipid 13": "Named in full and labelled EA-402; anchored in motif by "
                "Lipid 7 (tridecan-7-yl decanoate) and Lipid 14 (the "
                "N-CH2CH2-O head ester).",
}


def main():
    rows = []
    for name, a in ANCHORS.items():
        mol = Chem.MolFromSmiles(a["smiles"])
        if mol is None:
            print(f"FATAL: {name} anchor does not parse", file=sys.stderr)
            return 1
        rows.append({
            "lipid_name": name,
            "reference_SMILES": Chem.MolToSmiles(mol),
            "source_figure": a["figure"],
            "read_on": READ_ON,
            "read_by": READ_BY,
            "expected_mismatch": "",
            "note": a["note"],
        })
    for name, why in NO_ANCHOR.items():
        rows.append({
            "lipid_name": name, "reference_SMILES": "",
            "source_figure": "", "read_on": READ_ON, "read_by": READ_BY,
            "expected_mismatch": "",
            "note": "NO ANCHOR. " + why,
        })

    order = {f"Lipid {i}": i for i in range(1, 15)}
    rows.sort(key=lambda r: order[r["lipid_name"]])
    cols = ["lipid_name", "reference_SMILES", "source_figure", "read_on",
            "read_by", "expected_mismatch", "note"]
    with open(os.path.join(HERE, "reference_structures.csv"), "w",
              newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)

    # --- compare, AFTER writing, so the anchor is never derived from the
    # --- enumeration it is about to be compared with.
    enum = {}
    ep = os.path.join(HERE, "03_enumerated.csv")
    if os.path.exists(ep):
        with open(ep, newline="") as f:
            for r in csv.DictReader(f):
                enum[r["IL_name"]] = r["IL_SMILES"]

    print(f"{len(ANCHORS)} anchors, {len(NO_ANCHOR)} declared gaps "
          "-> reference_structures.csv")
    bad = 0
    for name, a in ANCHORS.items():
        ref = Chem.MolFromSmiles(a["smiles"])
        rk = Chem.MolToInchiKey(ref)
        got = enum.get(name)
        if got is None:
            print(f"  {name}: no enumerated structure to compare")
            continue
        gk = Chem.MolToInchiKey(Chem.MolFromSmiles(got))
        skel = rk.split("-")[0] == gk.split("-")[0]
        full = rk == gk
        verdict = ("exact_match" if full else
                   "SKELETON ok, STEREO differs" if skel else "MISMATCH")
        if not full:
            bad += 1
        print(f"  {name:>9}: {verdict}")
        print(f"             anchor {rk}")
        print(f"             enum   {gk}")
    print(f"  mismatches: {bad}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
