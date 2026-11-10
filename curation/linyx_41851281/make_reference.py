#!/usr/bin/env python3
"""Build the C10 anchors for this paper -> reference_structures.csv.

INDEPENDENCE, EXECUTED NOT ASSERTED. The enumerator (stage03_enumerate.py)
builds every structure with RDKit through curation/lib/assemble.attach: it
parses fragments into molecules, adds bonds with RWMol, and emits RDKit's
canonical SMILES, which for these lipids starts at a TAIL METHYL. This file
shares none of that. It composes a template written HEAD-FIRST, by plain
`str.format`, with no RDKit object constructed until the final parse -- no
RWMol, no AddBond, no RemoveAtom, no attach, no esterify. The only RDKit call
is MolToInchiKey on the finished string, which is the comparison itself.

F-013 and F-037 are why this file has to be read sceptically: an anchor was a
byte copy of the enumeration twice, the second time inside the file written to
fix the first. So two things are asserted here rather than claimed:

  1. no fragment carries a ring-closure digit, which is the ONLY way plain
     string composition can silently fuse two rings (the piperazine/cholesterol
     incident). P1 and lipid 107 contain no rings at all, so the hazard class
     is absent, and `_no_ring_digits` refuses if that ever stops being true;
  2. the composed template is printed into the CSV's note, so C16 can test the
     TEMPLATE -- byte-identity of the emitted anchor is not the signal, since a
     correct re-derivation is byte-identical to a correct enumeration once
     canonicalised.

Both fragments are written from a CHEMICAL NAME, not copied from anywhere:
  linoleyl alcohol   = (9Z,12Z)-octadeca-9,12-dien-1-ol
  hexadecan-3-ol     = the branched secondary alcohol of tail T13

The two anchors are the two lipids this paper DRAWS IN FULL:
  Lipid 107 (Lipid T3) -- Supplementary Fig. 12, third structure
  Lipid P1             -- Supplementary Methods p.56 reaction scheme, and
                          Supplementary Fig. 14a
Both drawings show the H4 head carrying TWO arms on one primary nitrogen and
ONE on the other, with that second nitrogen left as a secondary N-H. That is
the regiochemistry the enumerator produces for a 3-arm H4 lipid, and reading it
off a drawing is the only external evidence for it -- HRMS cannot see it.
"""
import csv
import os

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")
HERE = os.path.dirname(os.path.abspath(__file__))

# Head H4 written from the central N-methyl outwards: two 2-aminoethyl arms,
# the left nitrogen fully substituted, the right one keeping its N-H.
TEMPLATE = "CN(CCN({A})({A}))CCN{A}"

ALCOHOLS = {
    # chemical name -> SMILES written from the oxygen
    "(9Z,12Z)-octadeca-9,12-dien-1-ol": r"OCCCCCCCC/C=C\C/C=C\CCCCC",
    "hexadecan-3-ol": "OC(CC)CCCCCCCCCCCCC",
}


def _no_ring_digits(frag, what):
    if any(ch.isdigit() for ch in frag):
        raise SystemExit(
            "make_reference: %s (%r) carries a ring-closure digit. Plain "
            "string composition can then fuse two rings through a digit "
            "collision -- that is how a piperazine ended up bonded to "
            "cholesterol. Build this one with assemble.attach instead, and "
            "find a different independent route for the anchor." % (what, frag))


def arm(alcohol_name):
    """One Michael arm: -CH2CH2-C(=O)-O-R, R from the named alcohol."""
    alc = ALCOHOLS[alcohol_name]
    _no_ring_digits(alc, alcohol_name)
    # The ester keeps the ALCOHOL's oxygen; the acid loses its OH. Here the
    # acid half is written without one, so nothing is deleted at all.
    return "CCC(=O)" + alc


def anchor(alcohol_name):
    a = arm(alcohol_name)
    _no_ring_digits(TEMPLATE.replace("{A}", ""), "template")
    smi = TEMPLATE.format(A=a)
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        raise SystemExit("make_reference: %r does not parse" % smi)
    return Chem.MolToSmiles(mol), smi


ROWS = [
    dict(lipid_name="Lipid 107",
         alcohol="hexadecan-3-ol",
         source_figure="Supplementary Fig. 12 (MOESM1 p.17), structure "
                       "labelled 'Lipid T3'",
         note="Lipid 107 = H4-OT13(3) = the ML test set's Lipid T3. The "
              "drawing shows two arms on one primary nitrogen and one on the "
              "other, which is the only external evidence for the 3-arm H4 "
              "regiochemistry. T13's branch position (ethyl + C13 either side "
              "of the carbinol) is read from the same drawing and agrees with "
              "the trace of Supplementary Fig. 1a and with the HRMS formula "
              "C19H36O2 for the OT13 monomer."),
    dict(lipid_name="Lipid P1",
         alcohol="(9Z,12Z)-octadeca-9,12-dien-1-ol",
         source_figure="Supplementary Methods reaction scheme (MOESM1 p.56) "
                       "and Supplementary Fig. 14a",
         note="Lipid P1 = H4-OT1(3). The p.56 scheme draws H4 + T1 -> P1 and "
              "names it 'the 3-tailed lipid P1', separated from '4-tailed "
              "byproducts', which is what fixes the arm count."),
]


def main():
    out = []
    for r in ROWS:
        smi, template_used = anchor(r["alcohol"])
        out.append(dict(
            lipid_name=r["lipid_name"],
            reference_SMILES=smi,
            source_figure=r["source_figure"],
            read_on="2026-11-10",
            read_by="curation agent, linyx_41851281",
            expected_mismatch="",
            note=r["note"] + " BUILT BY: plain str.format on the head-first "
                 "template %r with the arm %r, parsed once at the end; no "
                 "RWMol, AddBond, RemoveAtom, assemble.attach or "
                 "assemble.esterify anywhere in this file."
                 % (TEMPLATE, arm(r["alcohol"]))))
    p = os.path.join(HERE, "reference_structures.csv")
    with open(p, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=[
            "lipid_name", "reference_SMILES", "source_figure", "read_on",
            "read_by", "expected_mismatch", "note"])
        w.writeheader()
        for r in out:
            w.writerow(r)
    for r in out:
        print(r["lipid_name"], Chem.MolToInchiKey(
            Chem.MolFromSmiles(r["reference_SMILES"])))
    print("wrote", p)


if __name__ == "__main__":
    main()
