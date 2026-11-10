#!/usr/bin/env python3
"""BUILD the C10 anchor for peptoids 236 and 266. Writes reference_structures.csv.

Run from this folder:   python3 -I make_reference.py

WHY THIS FILE EXISTS AT ALL. C10 compares each enumerated structure against the
structure the paper DRAWS, and it only tests anything if the anchor was built by
a different code path -- one that is EXECUTED, not asserted. That claim has been
made twice in this project and been false twice: F-013, where the Siegwart
anchor had the same sha256 as the enumeration whose note claimed independent
assembly, and F-037, where the file written to FIX F-013 re-copied the other
half (an AMINE_HALF template that was the enumeration's canonical SMILES with
the sterol swapped for '%'). C16 checks the GENERATOR'S COMPOSED TEMPLATE for
exactly that, because byte-identity of the emitted anchor is not the signal.

SO HERE IS WHAT IS DIFFERENT, ITEM BY ITEM. Four independent axes:

  1. DIFFERENT INPUT. stage03_enumerate.py reads the monomer ORDER out of SI
     Table S9's text layer. This file reads it out of SI page 20's DRAWING --
     specifically the fill colours and x-coordinates of the rounded boxes that
     stand for the two side chains. Pale yellow (0.996, 0.953, 0.812) is Ehx,
     pale pink (0.929, 0.788, 0.788) is Dod; the legend at the top of the page
     says so and the legend's own two boxes are excluded by y. Nothing is read
     by eye: the order is the sort order of the boxes' x centres.

  2. DIFFERENT TRAVERSAL. stage03 builds N-terminus -> C-terminus. This builds
     C-terminus -> N-terminus, starting at the Rink amide and ending at the Apd
     nitrogen.

  3. DIFFERENT CONSTRUCTION. stage03 uses curation/lib assemble.attach, atom
     maps and RDKit RWMol bond formation. This file does plain SMILES string
     composition and never imports assemble, never builds an RWMol, never uses
     an atom map.

  4. DIFFERENTLY SPELLED FRAGMENTS. The side chains are written here from the
     far end of the chain inward; blocks.py writes them from the attachment
     atom outward. blocks.py is NOT imported.

THE HAZARD STRING CONCATENATION ACTUALLY HAS. CLAUDE.md section 6 mandated
fragment assembly after a concatenated SMILES fused piperazine to cholesterol
through a RING-CLOSURE DIGIT COLLISION. That failure needs a ring. Every
fragment here is acyclic, and the file ASSERTS that before composing -- if a
future monomer with a ring is added, this refuses rather than silently
colliding digits.

WHAT THE ANCHOR CAN AND CANNOT TEST. It tests connectivity and monomer ORDER,
which is the one thing a molecular formula cannot see and the one thing that
defines a peptoid: 236, 237, 259, 264 and 266 all have calculated exact mass
1275.09 and are five different molecules. It does NOT test the Apd headgroup's
identity independently of Figure 1B, because page 20 draws the diol explicitly
but names nothing.
"""
import csv
import os
import sys

import pymupdf
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
SI = os.path.join(HERE, "source", "nn4c05513_si_001.pdf")
PAGE = 20

# The two side-chain box fills, as the PDF stores them. From the legend at the
# top of page 20, which prints "= ... Dod (C12)" and "= ... Ehx (C8)".
FILL = {
    (0.996, 0.953, 0.812): "Ehx",
    (0.929, 0.788, 0.788): "Dod",
}

# y bands in PDF points separating the legend from the two drawn structures.
# The legend sits above y=110; peptoid 236's backbone runs near y=145 and
# peptoid 266's near y=216.
BANDS = {236: (115.0, 175.0), 266: (185.0, 245.0)}

# Side chains. A fragment is appended as N(<frag>), and SMILES branch semantics
# bond the nitrogen to the fragment's FIRST atom -- so the attachment atom must
# be written first. The first version of this file wrote them far-end-first,
# which silently produced N-O-CH2-CH(CH2OH)2 for the headgroup (an N-O bond)
# and 5-methylheptyl for Ehx. Both parsed. Both canonicalised. The anchor
# simply disagreed with the enumeration, which is the anchor doing its job --
# but it is also exactly why CLAUDE.md section 6 distrusts string composition,
# and why the comparison against the enumeration is the thing that must pass,
# never the composition on its own.
#
# Spelled differently from blocks.py where the molecule allows it: the branch
# order is reversed on Ehx and the hydroxyls are written as explicit [OH] on
# Apd. n-dodecyl has only one attachment-first spelling, and that is fine --
# the independence this file claims is in its INPUT (the drawing's geometry,
# not Table S9's text) and its TRAVERSAL, not in how a straight chain is typed.
SIDE_CHAIN = {
    # 2-ethylhexyl: -CH2-CH(C4H9)(C2H5). Branch order reversed vs blocks.py.
    "Ehx": "CC(CCCC)CC",
    # n-dodecyl
    "Dod": "CCCCCCCCCCCC",
}
# 1,3-dihydroxypropan-2-yl: attaches through the CENTRAL carbon, which is what
# page 20 draws -- the backbone nitrogen meets the carbon bearing both
# hydroxymethyls, not an oxygen.
APD = "C(C[OH])C[OH]"

C_TERMINAL_AMIDE = "N"
RESIDUE = "C(=O)CN"


def boxes():
    """[(peptoid, [monomer codes in x order])] read from the page's drawing."""
    doc = pymupdf.open(SI)
    pg = doc[PAGE - 1]
    found = {p: [] for p in BANDS}
    for path in pg.get_drawings():
        f = path.get("fill")
        if f is None:
            continue
        key = tuple(round(c, 3) for c in f)
        if key not in FILL:
            continue
        r = path["rect"]
        cy = (r.y0 + r.y1) / 2
        cx = (r.x0 + r.x1) / 2
        for pid, (lo, hi) in BANDS.items():
            if lo <= cy <= hi:
                found[pid].append((cx, FILL[key]))
    doc.close()

    out = {}
    for pid, items in found.items():
        items.sort()
        xs = [x for x, _ in items]
        if len(items) != 6:
            raise SystemExit(
                f"peptoid {pid}: found {len(items)} side-chain boxes, expected "
                "6. Refusing to compose an anchor from an incomplete drawing.")
        # The backbone is drawn on a regular pitch, so the boxes must be evenly
        # spaced in x. If they are not, the y band has caught something that is
        # not a side chain, and the ORDER -- the only thing this anchor tests --
        # would be wrong while every atom count stayed right.
        gaps = [b - a for a, b in zip(xs, xs[1:])]
        if max(gaps) - min(gaps) > 2.0:
            raise SystemExit(
                f"peptoid {pid}: side-chain box x-spacing is irregular "
                f"({[round(g, 2) for g in gaps]}). The band may have caught a "
                "box from another structure.")
        out[pid] = [m for _, m in items]
    return out


def compose(order):
    """SMILES for the peptoid, built C-terminus -> N-terminus by string."""
    frags = [SIDE_CHAIN[m] for m in order] + [APD]
    for f in frags:
        if any(ch.isdigit() for ch in f):
            raise SystemExit(
                f"fragment {f!r} contains a ring-closure digit. String "
                "composition cannot safely be used: two fragments reusing the "
                "same digit fuse to each other. That is the piperazine/"
                "cholesterol failure CLAUDE.md section 6 records. Build this "
                "monomer with curation/lib assemble instead.")
    # terminal amide N, then one -C(=O)-CH2-N(R)- unit per residue, walking
    # from the C-terminus back to the Apd nitrogen.
    smi = C_TERMINAL_AMIDE
    for m in reversed(order):
        smi += RESIDUE + "(" + SIDE_CHAIN[m] + ")"
    smi += RESIDUE + "(" + APD + ")"
    # the final nitrogen carries the Apd side chain plus the backbone CH2 and
    # therefore one implicit hydrogen -- the free N-terminal secondary amine.
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        raise SystemExit(f"composed SMILES did not parse: {smi}")
    for a in mol.GetAtoms():
        if a.GetNumRadicalElectrons():
            raise SystemExit(f"radical atom in composed anchor: {smi}")
    return smi, Chem.MolToSmiles(mol)


def main():
    order = boxes()
    rows = []
    print("side-chain order read from SI page 20's own fill-box geometry:")
    for pid in sorted(order):
        raw, canon = compose(order[pid])
        print(f"  peptoid {pid}: Apd-" + "-".join(order[pid]))
        print(f"    composed template : {raw}")
        print(f"    canonical         : {canon}")
        rows.append({
            "lipid_name": f"peptoid_{pid}",
            "reference_SMILES": canon,
            "source_figure": "SI Figure S7A (nn4c05513_si_001.pdf page 20)",
            "read_on": "2026-11-10",
            "read_by": "make_reference.py -- side-chain identity and order "
                       "taken from the fill colour and x-centre of the "
                       "drawing's own rounded boxes; backbone composed "
                       "C-terminus to N-terminus by string, no assemble, no "
                       "atom maps, blocks.py not imported",
            "expected_mismatch": "",
            "note": "Peptoid drawn in full by the authors. This anchor tests "
                    "MONOMER ORDER, which is what a formula cannot see: 236, "
                    "237, 259, 264 and 266 all have calculated exact mass "
                    "1275.09 and are five different molecules.",
        })

    # The absence of an anchor is made visible, per CLAUDE.md section 5.
    for pid in (226, 245, 263):
        rows.append({
            "lipid_name": f"peptoid_{pid}",
            "reference_SMILES": "",
            "source_figure": "",
            "read_on": "2026-11-10",
            "read_by": "make_reference.py",
            "expected_mismatch": "",
            "note": "NO ANCHOR. The paper draws this peptoid nowhere -- only "
                    "236 and 266 are drawn in full (SI Figure S7A). Its "
                    "structure rests on SI Table S9's sequence, cross-checked "
                    "against SI Table S4 (34/34 sequences agree) and against "
                    "the printed calculated exact mass (44/44 within 0.005 "
                    "Da). Mass confirms composition and not order, so this "
                    "row has no external test of its monomer order.",
        })

    out = os.path.join(HERE, "reference_structures.csv")
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=[
            "lipid_name", "reference_SMILES", "source_figure", "read_on",
            "read_by", "expected_mismatch", "note"])
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
