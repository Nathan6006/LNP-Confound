#!/usr/bin/env python3
"""The C10 anchor: AMG1541, rebuilt from the paper's WORDS, not its picture.

    python3 -I make_reference.py      ->  reference_structures.csv

WHY THIS ROUTE. C10 only tests anything if the anchor came from a different
code path, and that has been claimed twice and been false twice: F-013 (a byte
copy of the enumeration, 537 rows of epicholesterol behind it) and F-037 (the
file written to FIX F-013 re-copied the other half -- 141 characters identical
once the sterol was swapped for a '%').

So this file shares NOTHING with the enumerator. The enumerator
(chemdraw_vector.py + stage02_blocks.py) reads filled-path coordinates out of
the PDF, builds an RDKit RWMol atom by atom, cuts it into blocks and re-joins
them with curation/lib/assemble. This file writes a SMILES string by hand from
two chemical NAMES in the paper's own prose, with no RWMol, no AddBond, no
RemoveAtom, no atom maps and no contact with any stage output:

  * "a hexacyclen headgroup (Am15)" -- Conclusions. Hexacyclen is
    1,4,7,10,13,16-hexaazacyclooctadecane: an 18-membered ring of six N-C-C
    units, written here as that repeat.
  * "monounsaturated 12 carbon tails" -- Conclusions; "Z-alkenes" -- Results.
    A 12-carbon Z-monoene acid esterified to glycidol, opened by the amine at
    the terminal CH2, gives N-CH2-CH(OH)-CH2-O-C(=O)-R. The double-bond
    POSITION is the one thing the prose does not fix, so it is taken from the
    drawing -- and that is said here rather than left implicit.

Written as a string, which CLAUDE.md warns about for good reason: concatenation
once fused piperazine to cholesterol through a ring-closure digit collision.
The ring-closure digit here is 1 and it appears exactly twice; the assertions
below check the ring count, the ring size, the formula and the absence of any
second ring closure, so a collision cannot pass silently.
"""
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from rdkit import Chem                                     # noqa: E402
from rdkit.Chem.rdMolDescriptors import CalcMolFormula     # noqa: E402

# (5Z)-dodec-5-enoyl, written from the name: C(=O) then C2 C3 C4, then the
# 5,6 Z-alkene, then C7..C12 as six more carbons ending in the methyl.
ACYL = "C(=O)CCC/C=C\\CCCCCC"

# The arm as it exists after the epoxide opening: the nitrogen's carbon, the
# secondary alcohol, then the glycidyl CH2-O carrying the acyl group.
ARM = "CC(O)CO" + ACYL

# Hexacyclen: six N-C-C units closed into one 18-membered ring. Digit 1 opens
# at the first nitrogen and closes at the last carbon; no other digit is used,
# so the ring-closure collision that once fused piperazine to cholesterol
# cannot occur -- and the assertions in main() check the ring anyway.


def main():
    # Build the 18-ring explicitly so the shape is visible in the source.
    smi = "N1(" + ARM + ")CC" + "".join(
        "N(" + ARM + ")CC" for _ in range(4)) + "N(" + ARM + ")CC1"
    mol = Chem.MolFromSmiles(smi)
    assert mol is not None, "the hand-written anchor does not parse"

    ri = mol.GetRingInfo()
    assert ri.NumRings() == 1, f"{ri.NumRings()} rings; expected exactly one"
    assert len(ri.AtomRings()[0]) == 18, \
        f"ring of {len(ri.AtomRings()[0])}; hexacyclen is 18-membered"
    ring = set(ri.AtomRings()[0])
    n_in_ring = sum(1 for i in ring if mol.GetAtomWithIdx(i).GetSymbol() == "N")
    assert n_in_ring == 6, f"{n_in_ring} ring nitrogens; hexacyclen has six"
    assert all(a.GetSymbol() != "N" or a.GetIdx() in ring
               for a in mol.GetAtoms()), "a nitrogen escaped the ring"
    dbl = [b for b in mol.GetBonds()
           if b.GetBondType() == Chem.BondType.DOUBLE
           and b.GetBeginAtom().GetSymbol() == "C"
           and b.GetEndAtom().GetSymbol() == "C"]
    assert len(dbl) == 6, f"{len(dbl)} C=C; one per arm expected"
    assert all(str(b.GetStereo()) in ("STEREOZ", "STEREOCIS") for b in dbl), \
        "an arm's alkene is not Z"
    formula = CalcMolFormula(mol)
    assert formula == "C102H186N6O18", formula
    canon = Chem.MolToSmiles(mol)

    with open(os.path.join(HERE, "reference_structures.csv"), "w",
              newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["compound_id", "reference_SMILES", "drawn_in",
                    "built_by", "note"])
        w.writerow([
            "AMG1541", canon, "SI p.109 (structure drawn beside the AMG1541 "
            "characterisation); Figure 5a and 5g in the article",
            "make_reference.py",
            "Composed as a SMILES string from two chemical names in the "
            "paper's prose -- 'hexacyclen headgroup (Am15)' and "
            "'monounsaturated 12 carbon tails' / 'Z-alkenes' -- with no "
            "RWMol, no AddBond, no RemoveAtom and no read of any stage "
            "output. The alkene POSITION (5,6) is the one feature the prose "
            "does not fix and is taken from the drawing."])
    print("anchor formula :", formula)
    print("anchor SMILES  :", canon)
    return 0


if __name__ == "__main__":
    sys.exit(main())
