#!/usr/bin/env python3
"""Build curation/zhangz_40433897/reference_structures.csv — the C10 anchor.

THE ANCHOR POSTURE FOR THIS PAPER IS: THERE IS NO ANCHOR FOR ANY VNxCyOz LIPID.

C10 compares an enumerated structure against the structure the paper DRAWS.
This paper draws no VNxCyOz structure at all: Figure S1b's product is drawn
with a generic `R`, and Figure S2 is NMR/HRMS traces with no structure diagram.
So every VNxCyOz row carries an EMPTY `reference_SMILES` with the reason, which
is the documented way to make the ABSENCE of an anchor visible rather than
silently missing.

ONE anchor is real and is built here: **VE-AC**, the tocopheryl acrylate
intermediate, which Figure S1b DOES draw in full (the chromanol with its
`(CH2CH2CH2CH(CH3))3`-collapsed tail bracket, the phenolic oxygen, the
carbonyl, the vinyl) and for which Figure S2(a) prints an exact mass. It is
the only external, structural ground truth in the paper, and it is the right
one to pin because it is the fragment that appears TWICE in every single
library member: a wrong VE-AC is a wrong corpus, and nothing else here could
detect it.

INDEPENDENCE, BY PROVENANCE AND NOT BY BYTE-DIFFERENCE
------------------------------------------------------
F-037: byte-identity of the emitted anchor is NOT the signal, and testing for
it is a trap, because a correct re-derivation canonicalises to the same string
as a correct enumeration. Independence here is established by what this file
does and does not do:

  * it imports NO assembler — `curation/lib/assemble` is never imported, there
    is no `esterify`, no `attach`, no `RWMol`, no `AddBond`, no `RemoveAtom`;
  * its anchor is a SOURCE STRING LITERAL, typed from the drawing;
  * that literal is written in the REVERSED traversal to the enumerator's.
    `stage03b_structures.py` writes tocopherol starting at the 5-methyl of the
    aromatic ring and closes the pyran last (`Cc1c(C)c2c(c(C)c1O)CC[C@@](C)(...)O2`);
    this file starts at the far END of the phytyl tail, walks inward to the
    pyran quaternary carbon, and closes the aromatic ring last. The acrylate is
    likewise written vinyl-first here and acid-first there.
  * a SMILES fragment has a direction and mass is blind to it, so the
    traversal is also checked explicitly below: the anchor must contain the
    ester written `C(=O)O-aryl` and not `O-C(=O)` reversed, and the chromanol
    oxygen must be in the ring, not on the acyl group.

The comparison that matters is the InChIKey, computed at the bottom, plus the
exact mass against the number Figure S2(a) prints.

Run: python3 -I make_reference.py
"""
import csv
import hashlib
import os
import sys

from rdkit import Chem, RDLogger
from rdkit.Chem import inchi
from rdkit.Chem.Descriptors import ExactMolWt
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

RDLogger.DisableLog("rdApp.*")
HERE = os.path.dirname(os.path.abspath(__file__))
ELECTRON = 0.00054857990946

# --------------------------------------------------------------------------
# VE-AC, read off Figure S1b (page 4, embedded raster xref 72, 220 dpi) and
# typed here as a source literal, tail-end first.
#
#   CC(C)CCC C(C) CCC C(C) CCC  [C@] (C) (O2) ... ring ... OC(=O)C=C
#
# i.e. isopropyl terminus -> 8'-methyl -> 4'-methyl -> C2 quaternary (the one
# assigned centre, 2R) -> C3,C4 of the pyran -> the aromatic ring, closing with
# the phenolic oxygen that bears the acryloyl group.
# --------------------------------------------------------------------------
VE_AC_ANCHOR = (
    "CC(C)CCCC(C)CCCC(C)CCC[C@]1(C)CCc2c(C)c(OC(=O)C=C)c(C)c(C)c2O1"
)

ROWS = [
    dict(lipid_name="VE-AC",
         reference_SMILES=VE_AC_ANCHOR,
         source_figure="Figure S1b (SI page 4), structure of the acrylate "
                       "intermediate; exact mass in Figure S2(a) caption",
         read_on="2026-11-10",
         read_by="curation agent, 220 dpi embedded raster at 3x upscale",
         expected_mismatch="",
         note="The only structure this paper draws in full. Built in this "
              "file as a source string literal typed tail-end first, the "
              "reverse of stage03b_structures.py's traversal; this module "
              "imports no assembler. The two 4'/8' side-chain centres are "
              "left UNASSIGNED on both sides because the paper names no "
              "isomer, so C10 must compare on a stereo-tolerant basis for "
              "those two centres while C2 (2R) is asserted."),
]

# Every VNxCyOz lipid: the anchor is ABSENT, on purpose and on the record.
ABSENT = [
    ("VN2C1", "panel (h); the paper's LEAD compound"),
    ("VN2C2", "panel (d)"),
    ("VN2C3", "named only in the Figure S5 legend; no HRMS panel at all"),
    ("VN2C4", "panel (c)"),
    ("VN3C1", "panel (f)"),
    ("VN3C2", "panel (e)"),
    ("VN3C2O", "panel (b)"),
    ("PANEL-G (label disputed)", "panel (g); duplicate label, see "
                                 "08_contradictions.csv X-1"),
]
for nm, where in ABSENT:
    ROWS.append(dict(
        lipid_name=nm, reference_SMILES="",
        source_figure="NONE - no structure of this compound is drawn anywhere "
                      "in the SI (%s)" % where,
        read_on="2026-11-10",
        read_by="curation agent",
        expected_mismatch="no_anchor_paper_draws_no_structure",
        note="Figure S1b draws the product with a generic R group and Figure "
             "S2 shows only NMR/HRMS traces. The structure in 03_enumerated."
             "csv, where present, is PROPOSED from the HRMS formula plus the "
             "Figure S1a CAS list - it is not a reading of a drawing, so "
             "there is nothing independent for C10 to compare it against. "
             "See DECISIONS.md D-1."))

HEADER = ["lipid_name", "reference_SMILES", "source_figure", "read_on",
          "read_by", "expected_mismatch", "note"]


def check_direction(smi):
    """A SMILES fragment has a direction; mass is blind to it.

    Three assertions that an exact mass would pass either way:
      1. the acryloyl is an ARYL ester, written C(=O)O-aryl, not reversed;
      2. the pyran oxygen is in a RING and is not the ester oxygen;
      3. the vinyl is terminal CH2=CH- and not an internal alkene.
    """
    m = Chem.MolFromSmiles(smi)
    if m is None:
        raise SystemExit("make_reference: the anchor does not parse.")
    tests = {
        "aryl_acrylate_ester": "[CH2]=[CH]C(=O)O[c]",
        "pyran_oxygen_in_ring": "[OX2;R]([CX4;R])[c;R]",
        "terminal_vinyl": "[CH2]=[CH][CX3](=O)",
        "no_free_phenol": "[c][OX2H1]",
    }
    out = {}
    for k, sma in tests.items():
        patt = Chem.MolFromSmarts(sma)
        n = len(m.GetSubstructMatches(patt))
        out[k] = n
    if out["aryl_acrylate_ester"] != 1:
        raise SystemExit("make_reference: the acrylate is not written as an "
                         "aryl ester (%d matches). A reversed O-C(=O) has the "
                         "same formula and the same exact mass." % out["aryl_acrylate_ester"])
    if out["pyran_oxygen_in_ring"] != 1:
        raise SystemExit("make_reference: the chromanol oxygen is not in the "
                         "ring (%d matches)." % out["pyran_oxygen_in_ring"])
    if out["no_free_phenol"] != 0:
        raise SystemExit("make_reference: a free phenol survives; the acrylate "
                         "did not go onto the phenolic oxygen.")
    return out


def main():
    # Independence, checked on the AST rather than by grepping the text:
    # no module named `assemble` is imported, and no RDKit mutation call
    # appears anywhere in this file's code.
    import ast
    tree = ast.parse(open(os.path.abspath(__file__)).read())
    imported, called = set(), set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported.update(a.name.split(".")[0] for a in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imported.add(node.module.split(".")[0])
        elif isinstance(node, ast.Call):
            f = node.func
            called.add(f.attr if isinstance(f, ast.Attribute)
                       else getattr(f, "id", ""))
    banned_imports = {"assemble"} & imported
    banned_calls = {"RWMol", "AddBond", "RemoveAtom", "CombineMols",
                    "esterify", "attach", "SetChiralTag"} & called
    if banned_imports or banned_calls:
        raise SystemExit(
            "make_reference: this module must share no assembly code path "
            "with the enumerator, but it imports %s and calls %s."
            % (sorted(banned_imports), sorted(banned_calls)))
    print("independence (AST): imports=%s ; no assembler import, "
          "no RDKit mutation call" % sorted(imported))

    direction = check_direction(VE_AC_ANCHOR)
    m = Chem.MolFromSmiles(VE_AC_ANCHOR)
    na = ExactMolWt(m) + 22.9897692809 - ELECTRON

    with open(os.path.join(HERE, "reference_structures.csv"), "w",
              newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=HEADER)
        w.writeheader()
        w.writerows(ROWS)

    print("anchor            :", VE_AC_ANCHOR)
    print("canonical         :", Chem.MolToSmiles(m))
    print("formula           :", CalcMolFormula(m))
    print("InChIKey          :", inchi.MolToInchiKey(m))
    print("skeleton          :", inchi.MolToInchiKey(m).split("-")[0])
    print("[M+Na]+ calc      : %.4f   (Fig S2a prints 507.3809, found "
          "507.3804)" % na)
    print("delta vs printed  : %+.2f mDa" % ((na - 507.3809) * 1000))
    print("direction checks  :", direction)
    print("sha256 of literal :",
          hashlib.sha256(VE_AC_ANCHOR.encode()).hexdigest()[:16])
    print("rows written      : %d (1 anchor, %d explicit absences)"
          % (len(ROWS), len(ABSENT)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
