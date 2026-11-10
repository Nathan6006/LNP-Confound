#!/usr/bin/env python3
"""BUILD the C10 anchors for lij_42130331. Do not store a string -- derive it.

    python3 -I make_reference.py

C10 compares each enumerated structure against the structure the paper DRAWS.
It is the only Tier-1 check that tests topology against external ground truth,
and it tests nothing unless the anchor was built by a DIFFERENT CODE PATH than
the enumerator. That independence has been asserted twice in this project and
been false twice -- F-013 (a byte copy, 537 rows of epicholesterol) and F-037
(the file written to fix F-013 re-copied the other half).

WHAT IS DIFFERENT HERE, concretely:

  * `stage03_enumerate.py` builds each molecule with RDKit `RWMol` surgery --
    `assemble.attach` (AddBond on atom-mapped atoms) and `assemble.esterify`
    (AddBond + RemoveAtom). This file performs NO molecule editing at all:
    there is no RWMol, no AddBond, no RemoveAtom and no atom map anywhere in
    it. It composes a SMILES string and parses it once.
  * The traversal is different. The enumerator starts from the amine block and
    hangs arms off it; its canonical output starts at the acyl tail. This
    template is written from the CENTRE OUT -- the tertiary amine's N-methyl
    first, then one arm, then the other -- which is how Fig. S2 reads the
    molecule left-to-right through its central nitrogen.
  * The acetamide is spelled `NC(=O)C` here and `NC(C)=O` in the enumerator,
    and the ester is reached from the carbinol carbon rather than from the
    acyl chain. Same molecule, different string.
  * Every chain is built from a REPEAT COUNT, not a copied run of `C`s, and
    each count is the number stage02 measured and 00_pages/zoom_*.png shows.

C16 checks this file's COMPOSED module-level templates against the
enumeration's canonical output, because byte-identity of the EMITTED anchor is
not the signal: a correct re-derivation is byte-identical to a correct
enumeration once canonicalised, which is exactly why C10's own
`anchor_may_be_a_copy` missed F-037. So `TEMPLATE_113` and `TEMPLATE_96` below
are module-level on purpose, for C16 to read.

WHAT THE ANCHORS DO NOT PROVE. They re-derive the drawing; they cannot test
the drawing. This paper publishes no HRMS, no NMR table, no SMILES/InChI and
no structure file, so there is no external datum against which the drawn
connectivity could be checked. That is stated on every emitted row through
`confidence = medium`, and it is why `08_contradictions.csv` and DECISIONS.md
carry the regiochemistry explicitly.
"""
import csv
import os

# --- measured counts, from stage02_blocks.py and the 2000 dpi zooms --------
N_ETHYLENE_113 = 2     # amine 113: H2N-(CH2)2-N(CH3)-(CH2)2-NH2
N_PROPYLENE_96 = 3     # amine 96:  CH3-NH-(CH2)3-NH2
R_CARBONS_C8 = 6       # 1,2-epoxyoctane leaves a hexyl on the carbinol carbon
ACYL_CARBONS_C14 = 14  # tetradecanoyl, from the printed n = 11 in Fig. S2
ACYL_CARBONS_C6 = 6    # hexanoyl, from the printed n = 3


def _aa_arm(n_link, r_carbons, acyl_carbons):
    """One beta-hydroxy-thioether arm, written carbinol-first.

    Fig. S1A fixes both bonds that a formula cannot see: the thiolactone opens
    at its CARBONYL, giving an amide to the amine and a free thiol; and the
    thiol adds to the LESS substituted epoxide carbon, so sulfur sits on the
    terminal CH2 and the hydroxyl (here acylated) on the carbon carrying R.
    """
    link = "C" * n_link
    r = "C" * r_carbons
    acyl = "C" * (acyl_carbons - 1)
    return f"{link}NC(=O)C(CCSCC({r})OC(=O){acyl})NC(=O)C"


def _n_arm(r_carbons, acyl_carbons):
    """The arm a SECONDARY amine makes by opening an epoxide directly.

    No sulfur and no amide -- Fig. S2 draws this on 96 and on 400, the only
    two amines with an N-H.
    """
    return f"CC({'C' * r_carbons})OC(=O){'C' * (acyl_carbons - 1)}"


# Composed at module level so C16 can read them.
TEMPLATE_113 = "CN({a}){a}".format(
    a=_aa_arm(N_ETHYLENE_113, R_CARBONS_C8, ACYL_CARBONS_C14))

TEMPLATE_96 = "CN({n}){a}".format(
    n=_n_arm(R_CARBONS_C8, ACYL_CARBONS_C6),
    a=_aa_arm(N_PROPYLENE_96, R_CARBONS_C8, ACYL_CARBONS_C6))

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   "reference_structures.csv")


def build():
    from rdkit import Chem
    rows = []
    for name, tpl, fig, note in (
        ("113-AA-C8C14", TEMPLATE_113, "Fig. S2, structure 113-AA-CxCy, with "
         "R and n taken from the same figure's legend (Cx = C8, Cy = C14 "
         "i.e. n = 11)",
         "Composed centre-out from the tertiary nitrogen as one SMILES "
         "string and parsed once. No RWMol, no AddBond, no RemoveAtom, no "
         "atom map in this file. Chain lengths are repeat counts of the "
         "values stage02_blocks.py measured."),
        ("96-AA-C8C6", TEMPLATE_96, "Fig. S2, structure 96-AA-CxCy, legend "
         "Cx = C8, Cy = C6 i.e. n = 3",
         "The two-arm case: amine 96's secondary N-H opens an epoxide "
         "directly, giving an arm with no sulfur and no amide. Built by the "
         "same string route, so C10 tests the arm-count rule and not only "
         "the thioether."),
    ):
        m = Chem.MolFromSmiles(tpl)
        if m is None:
            raise SystemExit(f"{name}: the composed template does not parse: "
                             f"{tpl}")
        rows.append((name, Chem.MolToSmiles(m), fig, "2026-11-10",
                     "curation agent, from the 2000 dpi renders in 00_pages/",
                     note))
    return rows


if __name__ == "__main__":
    rows = build()
    with open(OUT, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["lipid_name", "reference_SMILES", "source_figure",
                    "read_on", "read_by", "note"])
        w.writerows(rows)
    print(f"wrote {len(rows)} anchor(s) to {OUT}")
    for r in rows:
        print(" ", r[0], r[1])
