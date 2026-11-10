#!/usr/bin/env python3
"""Build the C10 anchors for li_39883839. Writes reference_structures.csv.

INDEPENDENCE, which is the only thing that makes a C10 anchor worth anything.
This has been claimed twice in this project and been false twice: F-013 (the
stored anchor was a byte copy of the enumeration, same sha256, while its note
claimed independent fragment assembly) and F-037 (the file written to FIX
F-013 re-copied the other half -- the amine template WAS the enumeration's
canonical SMILES with the sterol swapped for '%'). C16 therefore tests the
GENERATOR'S COMPOSED TEMPLATE, not the emitted anchor, because a correct
re-derivation is byte-identical to a correct enumeration once canonicalised.

So this file is built to be as unlike stage02_enumerate.py as it can be:

  stage02_enumerate.py        here
  ------------------------    ---------------------------------------------
  curation/lib/assemble       NO assemble, NO RWMol, NO AddBond, NO
  attach() + esterify()       RemoveAtom, NO atom maps. Plain "%s" string
                              substitution into one template per molecule.
  fragments are the paper's   fragments are transcribed by READING THE DRAWN
  SYNTHETIC PRECURSORS        STRUCTURE, as it is drawn, left to right
  (N-alkylformamide, primary  (SI Fig S1/S2 for the library members, Fig 3G
  alcohol)                    for A4B4-S3), never a reagent
  traversal: tail-in          traversal: HEAD-FIRST -- amine, ester carbonyl,
                              alpha-carbon, its branch, then the amide

RDKit is used only to canonicalise and to compute the InChIKey. The topology
is entirely in the template strings below.

WHAT WAS READ, AND WHERE. SI Fig S1 draws H1A1B1: a blue Me2N-CH2-C(=O)-O-
attached to a carbon that carries (a) a saturated branch running down-right
and (b) a C(=O)-N(H)- going right into a long chain. Fig 3G draws A4B4-S3 as
Me2N-(CH2)3-C(=O)-O- over an 'R =' that is the amide-to-oleyl plus the
pentamethylene-ester-to-2-hexyldecanoate.
"""
import csv, os, sys
from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")
HERE = os.path.dirname(os.path.abspath(__file__))

# --- transcribed from the drawing, head group first ------------------------
HEAD = {
 "H1": "CN(C)CC(=O)O",            # Me2N-CH2-COO-     SI Fig S1 (drawn blue)
 "H2": "CN(C)CCC(=O)O",           # Me2N-(CH2)2-COO-  SI Fig S2, Fig 3G n=2
 "S3": "CN(C)CCCC(=O)O",          # Me2N-(CH2)3-COO-  Fig 3G "Spacer = 3"
}
# the alpha carbon's own branch, as drawn, reading away from the alpha carbon
ALPHA_BRANCH = {
 "B1": "CCCCCCCCC",                               # nonyl           SI Fig S1
 "B4": "CCCCCOC(=O)C(CCCCCC)CCCCCCCC",            # (CH2)5-O-CO-CH(hexyl)octyl
}
# the amide nitrogen's chain, as drawn, reading away from the nitrogen
AMIDE_CHAIN = {
 "A1": "CCCCCCCCCCC",                             # undecyl
 "A4": r"CCCCCCCC/C=C\CCCCCCCC",                  # (Z)-octadec-9-en-1-yl
}

# ONE template. The molecule is written in the order the eye reads the figure.
TEMPLATE = "%(head)sC(%(branch)s)C(=O)N%(amide)s"

ANCHORS = [
 ("H1A1B1",  "H1", "B1", "A1", "SI Appendix Fig. S1 (top-left structure)"),
 ("H2A4B4",  "H2", "B4", "A4", "SI Appendix Fig. S2; redrawn in Fig. 3G "
                               "('Spacer = 2' head plus the shared R group)"),
 ("A4B4-S3", "S3", "B4", "A4", "Fig. 3G ('Spacer = 3' head plus the shared "
                               "R group labelled 'Top performing tail "
                               "combination by batch test')"),
]


def compose(h, b, a):
    """The composed template, BEFORE any RDKit touches it. This string is
    what C16 should examine."""
    return TEMPLATE % {"head": HEAD[h], "branch": ALPHA_BRANCH[b],
                       "amide": AMIDE_CHAIN[a]}


def main():
    enum = {r["IL_name"]: r for r in
            csv.DictReader(open(os.path.join(HERE, "03_enumerated.csv")))}
    out, bad = [], 0
    for name, h, b, a, where in ANCHORS:
        raw = compose(h, b, a)
        m = Chem.MolFromSmiles(raw)
        if m is None:
            raise SystemExit("make_reference: %s template does not parse: %s"
                             % (name, raw))
        canon = Chem.MolToSmiles(m)
        key = Chem.MolToInchiKey(m)
        e = enum[name]
        same = (key == e["InChIKey"])
        bad += (not same)
        # The independence evidence: the TEMPLATE must not be the
        # enumeration's canonical string. Byte-identity of the EMITTED anchor
        # proves nothing (F-037) -- both are canonicalised.
        copied = (raw == e["IL_SMILES"])
        out.append(dict(
            lipid_name=name, reference_SMILES=canon, source_figure=where,
            read_on="2026-11-10", read_by="curation agent li_39883839",
            expected_mismatch="",
            note=("Composed head-first by plain string substitution into "
                  "'%s' from fragments transcribed off the drawn structure; "
                  "no assemble/RWMol/atom-map code path is shared with "
                  "stage02_enumerate.py. Pre-canonicalisation template: %s . "
                  "Template differs from the enumeration string: %s. "
                  "InChIKey vs enumeration: %s."
                  % (TEMPLATE, raw, "YES" if not copied else "NO -- COPY",
                     "match" if same else "MISMATCH")))) 
        print("%-9s %s  template_is_a_copy=%s" %
              (name, "MATCH" if same else "MISMATCH", copied))
        print("   template  :", raw)
        print("   enumerator:", e["IL_SMILES"])
    with open(os.path.join(HERE, "reference_structures.csv"), "w",
              newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["lipid_name", "reference_SMILES",
                                           "source_figure", "read_on",
                                           "read_by", "expected_mismatch",
                                           "note"])
        w.writeheader(); w.writerows(out)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
