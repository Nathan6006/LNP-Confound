#!/usr/bin/env python3
"""Build the C10 anchor for H9T6 by a DIFFERENT CODE PATH than the enumerator.

C10 only tests anything if the anchor was re-derived, and that claim has been
false twice in this project (F-013, a byte copy; F-037, the file written to fix
F-013 re-copied the other half). So the differences here are deliberate and
mechanical, not asserted:

  * 02_enumerate.py builds the arm by RDKit graph surgery -- RWMol.RemoveBond
    on the oxirane, then assemble.attach with atom maps, six times, each pass
    picking "the first remaining N-H" by SMARTS. This file forms NO bond with
    RDKit at all. It does plain string substitution on a template.
  * 02_enumerate.py gets the tail's stereocentre by enumerating both chiral
    tags on a flat constitution and selecting the one RDKit's CIP labels "S".
    This file writes the arm's chirality from the AMINE end, in the opposite
    traversal direction, and lets canonicalisation decide whether the two
    agree. A disagreement would show up as an InChIKey mismatch.
  * the two are compared by FULL InChIKey, so the stereo layer is in scope --
    that is the layer F-013 lived in and the skeleton block cannot see.

SOURCE: SI p.5, the named preparation of H9T6: "A solution of
Tris(2-aminoethyl)amine ... and the desired tail ... 5.7 equiv", with
0.95 x N and N = 2 x (number of primary amines) = 6 from the general procedure
on SI p.4; the tail identified as (S)-glycidyl decanoate in 07_qc.json by three
independent routes (LC-MS ladder, 1H NMR integrals, the R/S pairing of H9T2).

The paper draws NO structure of H9T6 anywhere in the SI, and the main-text
figure is not on disk. This anchor is therefore derived from the paper's
WRITTEN synthesis, not from a drawing; PROVENANCE.md records that posture.
"""
import os, csv
from rdkit import Chem

HERE = os.path.dirname(os.path.abspath(__file__))

# One arm, written from the NITROGEN outward: N-CH2-CH(OH)-CH2-O-C(=O)-C9H19.
# "%" marks where the nitrogen goes. The stereocentre is written here, in this
# traversal, and nowhere else in this file.
ARM = "%C[C@H](O)COC(=O)CCCCCCCCC"

# tris(2-aminoethyl)amine with every N-H marked. Six marks, three nitrogens,
# two each -- which is the 0.95 x N = 6 stoichiometry of the SI's procedure.
CORE = "N(CCN(%)%)(CCN(%)%)CCN(%)%"

def build():
    s = CORE
    while "%" in s:
        s = s.replace("%", ARM.replace("%", "", 1), 1) if False else \
            s.replace("%", "\x00", 1).replace("\x00", ARM.replace("%", ""), 1)
    m = Chem.MolFromSmiles(s)
    if m is None:
        raise SystemExit(f"make_reference: unparseable template product: {s}")
    return Chem.MolToSmiles(m)

if __name__ == "__main__":
    smi = build()
    m = Chem.MolFromSmiles(smi)
    Chem.AssignStereochemistry(m, cleanIt=True, force=True)
    codes = [a.GetProp("_CIPCode") for a in m.GetAtoms() if a.HasProp("_CIPCode")]
    assert len(codes) == 6 and len(set(codes)) == 1, codes
    out = os.path.join(HERE, "reference_structures.csv")
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["lipid_name", "reference_SMILES", "source_figure",
                    "read_on", "read_by", "expected_mismatch", "note"])
        w.writerow([
            "H9T6", smi, "SI p.5 (written preparation; NO structure is drawn)",
            "2026-11-10", "curation agent chamolstad_41691524", "",
            "Assembled by STRING SUBSTITUTION on a template, with no RDKit "
            "bond operation anywhere and the arm written from the nitrogen "
            "outward -- the opposite traversal to 02_enumerate.py, which uses "
            "RWMol.RemoveBond plus assemble.attach and takes the stereocentre "
            "from a CIP selection on a flat ester. The paper draws nothing, so "
            "this anchor tests the enumerator against the SI's WORDS, not "
            "against a drawing; it cannot catch an error the written "
            "procedure also contains."])
    print("reference_structures.csv written")
    print("anchor  ", smi)
    print("CIP     ", codes)
    print("InChIKey", Chem.MolToInchiKey(m))
