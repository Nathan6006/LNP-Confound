#!/usr/bin/env python3
"""Build the C10 anchor for SD_2026 -- by a route the enumerator does not share.

WHY THIS FILE EXISTS, AND WHY IT IS WRITTEN THE WAY IT IS. C10 only tests
anything if the anchor was built by a DIFFERENT CODE PATH, EXECUTED, not
asserted. That claim has been made twice in this project and was false twice:
F-013 stored a byte copy of the enumeration whose note claimed independent
assembly, and F-037 found that the file written to FIX F-013 had re-copied the
other half of the molecule. C16 exists because byte-identity of the emitted
anchor is not the signal -- a correct re-derivation is byte-identical to a
correct enumeration once canonicalised.

So what differs here is the ROUTE, and it differs in four ways:

  1. DIFFERENT DECOMPOSITION. stage03 splits LC-1 the way the FIGURE does --
     headgroup / linker / tail a / tail b. This file splits it the way the
     SYNTHESIS does, into the two amide halves that Supplementary Table 5's
     five steps actually build. No fragment here is a block of 01_blocks.csv.
  2. DIFFERENT SOURCE FOR THE STRINGS. Every fragment below is written from a
     CHEMICAL NAME printed in Supplementary Table 5 -- "2-hexylundecanoic
     acid", "6-amino-1-hexanol", "pivaldehyde", "myristic acid",
     "3-(pyrrolidin-1-yl)propan-1-amine" -- and not read from any file this
     repository generates.
  3. DIFFERENT MACHINERY. stage03 uses one mapped reaction SMARTS. This file
     uses assemble.esterify and assemble.attach, which are RWMol operations.
     There is no reaction SMARTS anywhere in this file.
  4. DIFFERENT ORDER. stage03 condenses all four components at once. This file
     forms the ester first and the tertiary amide last, which is the order the
     paper's Steps 2 and 5 run in.

WHAT THE ANCHOR IS FOR. Mass cannot see connectivity. 07_qc.json records a
seeded swap of LC-1's headgroup, H4 -> H3, that leaves the molecular formula
IDENTICAL -- the formula check is blind to it. The InChIKey SKELETON block does
move, and that is what provenance/structure_check.py compares. The anchor is
the only thing in this folder with any power over connectivity.

WHAT IT IS NOT. It does not test the Fig. 1l DRAWING against the enumeration
for the other nine named lipids, because this curation has not determined which
of the 384 enumerated members SC-1..5 and LC-2..5 are (DECISIONS.md D-2). Those
rows carry an EMPTY reference_SMILES so the gap is visible rather than absent.
"""
import csv
import os
import sys

from rdkit import Chem

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "lib"))
import assemble                                                  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))

# --- fragments, each written from a NAME in Supplementary Table 5 -----------
# Step 2's acid. "2-hexylundecanoic acid": an undecanoic (C11) acid bearing a
# hexyl at C2. Written as C3-C11 (nine carbons), then C2 with its hexyl, then
# the carboxyl -- so the C11 count is in the string and can be counted off it.
ACID_2_HEXYLUNDECANOIC = "CCCCCCCCC" "C(CCCCCC)" "C(=O)O"

# Steps 1-5 collapsed into the half of the molecule that carries the ester:
# the pivaldehyde-derived neopentyl carbon, the amide to the isocyanide-derived
# nitrogen, and the six-carbon chain from 6-amino-1-hexanol ending in its OH.
# = N-(6-hydroxyhexyl)-3,3-dimethylbutanamide. The atom map is on the carbon
# that becomes the stereocentre, so the final bond is made BY MAP, never by
# "atom index 0 of the fragment".
HALF_ISOCYANIDE_SIDE = "CC(C)(C)[CH2:2]C(=O)NCCCCCCO"

# The other half: myristic acid's amide to 3-(pyrrolidin-1-yl)propan-1-amine.
# = N-[3-(pyrrolidin-1-yl)propyl]tetradecanamide. Map on the amide nitrogen.
HALF_AMINE_SIDE = "CCCCCCCCCCCCCC(=O)[N:1]CCCN1CCCC1"

HEADER = ("lipid_name,reference_SMILES,source_figure,read_on,read_by,"
          "expected_mismatch,note")

READ_ON = "2026-11-10"
READ_BY = "curation agent (sd_2026), from Fig. 1l and Supplementary Table 5"


def build_lc1():
    """LC-1, assembled the way Supplementary Table 5 synthesises it."""
    ester = assemble.esterify(
        ACID_2_HEXYLUNDECANOIC, HALF_ISOCYANIDE_SIDE,
        label="Table 5 step 2: 2-hexylundecanoic acid + the 6-hydroxyhexyl "
              "amide")
    lc1 = assemble.attach(
        HALF_AMINE_SIDE, ester, a_map=1, b_map=2,
        label="Table 5 step 5: the Ugi tertiary amide N-C bond")
    return lc1


def main():
    lc1 = build_lc1()
    m = Chem.MolFromSmiles(lc1)
    rows = [{
        "lipid_name": "LC-1",
        "reference_SMILES": lc1,
        "source_figure": "Fig. 1l (leftmost structure) + Supplementary "
                         "Table 5 (five-step synthesis with CAS numbers)",
        "read_on": READ_ON, "read_by": READ_BY,
        "expected_mismatch": "",
        "note": "Assembled from the SYNTHESIS, not from the block figure: "
                "2-hexylundecanoic acid esterified onto "
                "N-(6-hydroxyhexyl)-3,3-dimethylbutanamide "
                "(assemble.esterify), then bonded to "
                "N-[3-(pyrrolidin-1-yl)propyl]tetradecanamide at its amide "
                "nitrogen (assemble.attach, by atom map). No reaction SMARTS "
                "and no block from 01_blocks.csv is used here, which is what "
                "makes this independent of stage03_enumerate.py. The Ugi "
                "stereocentre is left unassigned in both routes: the reaction "
                "is not stereoselective and the paper reports no resolution.",
    }]
    for nm in ("SC-1", "SC-2", "SC-3", "SC-4", "SC-5",
               "LC-2", "LC-3", "LC-4", "LC-5"):
        panel = "Fig. 1g" if nm.startswith("SC") else "Fig. 1l"
        rows.append({
            "lipid_name": nm, "reference_SMILES": "",
            "source_figure": panel, "read_on": READ_ON, "read_by": READ_BY,
            "expected_mismatch": "",
            "note": "The paper DOES draw this structure, but this curation "
                    "has not determined which of the 384 enumerated members "
                    "it is -- the SI's printed formula matches between 2 and "
                    "12 members in every case, and mass cannot choose among "
                    "them. The anchor is therefore left EMPTY on purpose so "
                    "the gap is visible; see DECISIONS.md D-2. Only LC-1 is "
                    "pinned, by Supplementary Table 5's CAS numbers.",
        })

    with open(os.path.join(HERE, "reference_structures.csv"), "w",
              newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=HEADER.split(","))
        w.writeheader()
        for r in rows:
            w.writerow(r)

    print("LC-1 anchor :", lc1)
    print("InChIKey    :", Chem.MolToInchiKey(m))
    print("skeleton    :", Chem.MolToInchiKey(m).split("-")[0])
    enum = os.path.join(HERE, "03_enumerated.csv")
    if os.path.exists(enum):
        for r in csv.DictReader(open(enum)):
            if r["lipid_id"] == "H4-L4-Ta1-Tb1":
                e = Chem.MolFromSmiles(r["canonical_SMILES"])
                ek = Chem.MolToInchiKey(e)
                ak = Chem.MolToInchiKey(m)
                print("enumeration :", r["canonical_SMILES"])
                print("byte-identical to the anchor string:",
                      r["canonical_SMILES"] == lc1,
                      "   (byte-identity is NOT the signal -- see C16)")
                print("skeleton match:", ek.split("-")[0] == ak.split("-")[0])
                print("full InChIKey match:", ek == ak)
                print("formula match:",
                      Chem.rdMolDescriptors.CalcMolFormula(e)
                      == Chem.rdMolDescriptors.CalcMolFormula(m))


if __name__ == "__main__":
    main()
