#!/usr/bin/env python3
"""BUILD the C10 anchor for siegwart_40644303. Do not store a string -- derive it.

    python3 -I make_reference.py

THE ANCHOR is 3A4-SC5, the one complete ionizable lipid this paper DRAWS: the
final product of Scheme 1, drawn atom by atom with every heteroatom labelled.
Figure 1b draws only building blocks, and Figures 2b/3a draw only amine cores,
so this is the single full structure available. There is no second one, and
reference_structures.csv says so rather than inventing one.

WHY THIS ROUTE, AND HOW IT DIFFERS FROM THE ENUMERATOR.

  stage03_enumerate.py   builds every lipid with curation/lib/assemble.attach,
                         i.e. RDKit RWMol + AddBond on SMILES fragments, core
                         first, one arm at a time.

  this file              writes an explicit V2000 MOL BLOCK -- a numbered atom
                         table and a numbered bond table, typed from Scheme 1 --
                         and hands it to Chem.MolFromMolBlock. There is no
                         SMILES fragment, no string concatenation of chemistry,
                         no RWMol, no AddBond and no RemoveAtom anywhere in
                         this file. The two routes share RDKit's sanitiser and
                         nothing else.

That matters because the independence of a C10 anchor has been asserted twice
and been false twice. F-013: the stored anchor was a byte copy of the
enumeration, same sha256, while its note claimed independent fragment assembly
"through the C3 oxygen (3-beta)" and the string was 3-alpha -- 537 rows of
epicholesterol, on the flagship curation, for the one molecule C10 existed to
protect. F-037: the file written in round 1 TO FIX THAT re-copied the other
half, a 141-character template that was the enumeration's canonical SMILES with
the sterol swapped for `%`.

A SMILES-based anchor would have walked straight into that here, and this is
worth recording: the enumerator's canonical SMILES for 3A4-SC5 happens to START
AT THE PENTYL TAIL --
  CCCCCSCC(C)C(=O)OCCOC(=O)CCN(CCO)CCCN(...)...
-- so the "obvious" independent route, writing the molecule tail-first from the
drawing, would have produced a byte-identical string and C16 would have been
right to call it a copy. An atom-and-bond table cannot collide that way.

WHAT THE TABLE SAYS, read off Scheme 1 head-first:
  atoms 1-8    HO-CH2-CH2-N(<)-CH2-CH2-CH2-N(<)(<)      the 3A4 core, (2-3)
  arm, 19 atoms, three times:
               N-CH2-CH2-C(=O)-O-CH2-CH2-O-C(=O)-CH(-CH3)-CH2-S-C5H11
  one arm on the former secondary nitrogen, two on the former primary one,
  which is the 3-arm pattern Scheme 1 draws and the "# of Arms = 3" bracket in
  Figure 1a states.

THE REGIOCHEMISTRY IS READ FROM THE SYNTHESIS, NOT FROM A MASS. The thiol adds
to the methacrylate's beta carbon, so sulfur sits on the CH2 and the methyl on
the alpha carbon. Swapping them gives the identical molecular formula and the
identical exact mass, so HRMS, the SI's formula and every composition check are
blind to it. Only Scheme 1's drawing and the written two-step procedure
(aza-Michael onto the acrylate, sulfa-Michael onto the methacrylate) say which.
"""
import csv
import os

from rdkit import Chem

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "reference_structures.csv")

# --- the atom table, typed from Scheme 1 ---------------------------------
CORE = ["O", "C", "C", "N", "C", "C", "C", "N"]
CORE_BONDS = [(1, 2, 1), (2, 3, 1), (3, 4, 1), (4, 5, 1),
              (5, 6, 1), (6, 7, 1), (7, 8, 1)]

# the arm, from the nitrogen outward. 19 atoms, 18 bonds, acyclic.
ARM_ATOMS = ["C", "C", "C", "O", "O", "C", "C", "O", "C", "O",
             "C", "C", "C", "S", "C", "C", "C", "C", "C"]
#            0    1    2   3=O  4-O  5    6    7-O  8   9=O
#                                                    10 CH 11 CH3 12 CH2 13 S
ARM_BONDS = [(0, 1, 1), (1, 2, 1), (2, 3, 2), (2, 4, 1), (4, 5, 1),
             (5, 6, 1), (6, 7, 1), (7, 8, 1), (8, 9, 2), (8, 10, 1),
             (10, 11, 1), (10, 12, 1), (12, 13, 1), (13, 14, 1),
             (14, 15, 1), (15, 16, 1), (16, 17, 1), (17, 18, 1)]
ARM_ON = [4, 8, 8]          # core atom numbers each arm hangs from


def molblock():
    """An explicit V2000 atom/bond table. No SMILES anywhere in this function."""
    atoms = list(CORE)
    bonds = list(CORE_BONDS)
    for anchor_atom in ARM_ON:
        off = len(atoms)
        atoms.extend(ARM_ATOMS)
        bonds.append((anchor_atom, off + 1, 1))
        for i, j, o in ARM_BONDS:
            bonds.append((off + 1 + i, off + 1 + j, o))

    if len(bonds) != len(atoms) - 1:
        raise SystemExit(
            f"make_reference: {len(atoms)} atoms and {len(bonds)} bonds. An "
            "acyclic molecule has exactly n-1 bonds; this one does not, which "
            "means the table is either disconnected or has closed a ring. "
            "3A4-SC5 has no ring, so a ring here would be the string-"
            "concatenation ring-closure collision in a different costume.")

    lines = ["3A4-SC5 reference", "  typed from Scheme 1", ""]
    lines.append(f"{len(atoms):>3}{len(bonds):>3}  0  0  0  0  0  0  0  0999 V2000")
    for el in atoms:
        lines.append(f"{0.0:>10.4f}{0.0:>10.4f}{0.0:>10.4f} {el:<3} 0  0"
                     "  0  0  0  0  0  0  0  0  0  0")
    for i, j, o in bonds:
        lines.append(f"{i:>3}{j:>3}{o:>3}  0  0  0  0")
    lines.append("M  END")
    return "\n".join(lines), atoms, bonds


def build():
    mb, atoms, bonds = molblock()
    mol = Chem.MolFromMolBlock(mb)
    if mol is None:
        raise SystemExit("make_reference: the MOL block did not parse.")
    Chem.SanitizeMol(mol)
    if mol.GetRingInfo().NumRings() != 0:
        raise SystemExit("make_reference: the product has a ring; 3A4-SC5 "
                         "has none.")
    if any(a.GetNumRadicalElectrons() for a in mol.GetAtoms()):
        raise SystemExit("make_reference: a radical atom in the product.")
    smi = Chem.MolToSmiles(mol)
    note = (
        "Built from an explicit V2000 atom/bond table typed from Scheme 1's "
        "drawn product, parsed with Chem.MolFromMolBlock. No SMILES fragment, "
        "no string concatenation of chemistry, no RWMol, no AddBond and no "
        "RemoveAtom in this file; the enumerator (stage03_enumerate.py) uses "
        "curation/lib/assemble.attach, i.e. RWMol+AddBond on SMILES fragments. "
        "Writing the anchor as a SMILES read tail-first from the drawing WOULD "
        "have collided: the enumerator's canonical string for this lipid "
        "already starts at the pentyl tail. Regiochemistry (S on the "
        "methacrylate beta carbon, methyl on alpha) is from the two-step "
        "synthesis, which is the only evidence that distinguishes it -- the "
        "swapped assignment has the same formula and the same exact mass.")
    return [("3A4-SC5", smi, "Scheme 1, drawn product", "2026-11-10",
             "curation agent, from the atom/bond table in this file", "", note)]


if __name__ == "__main__":
    rows = build()
    rows.append((
        "", "", "Figure 1b / Figures 2b, 3a", "2026-11-10",
        "curation agent", "",
        "NO SECOND ANCHOR EXISTS. Figure 1b draws the 15 amine cores and the "
        "31 thiol tails as separate building blocks and Figures 2b and 3a draw "
        "amine cores only, so Scheme 1's 3A4-SC5 is the ONLY complete lipid "
        "this paper draws. The absence is recorded here rather than left "
        "invisible: 488 of the 489 enumerated structures have no drawn "
        "counterpart to test against, and C10 therefore covers 1 of them."))
    with open(OUT, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["lipid_name", "reference_SMILES", "source_figure",
                    "read_on", "read_by", "expected_mismatch", "note"])
        w.writerows(rows)
    print(f"wrote {len(rows)} anchor row(s) to {OUT}")
    print(rows[0][1])
