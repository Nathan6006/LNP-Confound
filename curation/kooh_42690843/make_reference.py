#!/usr/bin/env python3
"""The C10 anchors: built from the AUTHORS' NAMES, by a different code path.

WHY THIS FILE IS A GENERATOR AND NOT A STORED STRING. A C10 anchor tests
nothing unless it was produced by a route that could have disagreed with the
enumerator. That claim has been made twice in this project and was false twice:
F-013 stored a byte copy of the enumeration (537 rows of epicholesterol sailed
through), and F-037 found that the file written to FIX F-013 had re-copied the
other half of the molecule. C16 therefore inspects the GENERATOR'S COMPOSED
TEMPLATE, not the emitted string -- byte-identity of the emitted anchor is not
the signal, because a correct re-derivation is byte-identical to a correct
enumeration once canonicalised.

WHAT IS DIFFERENT HERE, concretely:

  SOURCE     stage03_enumerate.py builds from the BLOCK STRUCTURES drawn on
             SI p.S5. This file builds from the AUTHORS' OWN IUPAC NAMES,
             printed on SI p.S12 and p.S18-S19 -- a different statement of the
             molecule by different people, and the one that carries the
             regiochemistry.
  MECHANISM  stage03 uses RDKit fragment assembly (curation/lib assemble.attach:
             CombineMols + AddBond on atom-mapped atoms). This file uses PLAIN
             STRING CONCATENATION of name-derived tokens and touches no RDKit
             API until the string is complete. There is no RWMol, no AddBond,
             no RemoveAtom, no atom map anywhere below.
  TRAVERSAL  stage03 starts at the alpha carbon. This file writes each molecule
             from the far end of the acyl chain (A18B8C14) or from the
             cyclopentane (A2B1C1) inward.

STRING CONCATENATION IS THE OTHER HAZARD, so it is guarded rather than
trusted. CLAUDE.md section 6: concatenation once fused piperazine to
cholesterol through a ring-closure digit collision. `compose()` below REFUSES
to join fragments unless at most one fragment carries a ring-closure digit and
no digit appears in two fragments -- which makes that collision impossible by
construction, asserted on every run rather than asserted in a comment.

Run:  python3 -I make_reference.py            writes reference_structures.csv
      python3 -I make_reference.py --self-test
            re-runs with a deliberately collided pair of ring fragments and
            with a one-carbon-short token, and requires both to be refused.
"""
import csv
import os
import re
import sys

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "reference_structures.csv")

READ_ON = "2026-11-10"
READ_BY = "curation agent, from the SI text layer (pymupdf), SI p.S12 and "\
          "p.S18-S19, cross-read against the 600 dpi renders in 00_pages/"

# ---------------------------------------------------------------------------
# Name-derived tokens. Each line is a piece of an IUPAC name, written out as
# SMILES by hand from that piece alone. Nothing here was copied from
# 01_blocks.csv, 03_enumerated.csv or any RDKit output.
# ---------------------------------------------------------------------------
TOK = {
    # "pentane-1,5-diyl"
    "pentane_1_5_diyl":      "CCCCC",
    # "propane-1,3-diyl" (from propane-1,3-diamine)
    "propane_1_3_diyl":      "CCC",
    # "azanediyl" / "-amine" nitrogen
    "azanediyl":             "N",
    # "N,N-dimethylamino" terminus
    "dimethylamino":         "N(C)C",
    # "cyclopentyl"
    "cyclopentyl":           "C1CCCC1",
    # the C2..C9 remainder of a "nonyl" once C1 is written as the branch point
    "nonyl_C2_C9":           "CCCCCCCC",
    # "tetradecane-6,1-diyl": C1..C5 written before the C6 branch point ...
    "tetradecane_C1_C5":     "CCCCC",
    # ... and C7..C14 written as C6's branch
    "tetradecane_C7_C14":    "CCCCCCCC",
    # ester linkage "-yl ...oate"
    "ester_O":               "O",
    "carbonyl":              "C(=O)",
    # "2-hexyldecanoate": the hexyl on C2 ...
    "hexyl":                 "CCCCCC",
    # ... and C3..C10 of the decanoate
    "decanoate_C3_C10":      "CCCCCCCC",
    # a branch point: the atom that opens a parenthesis
    "branch_C_open":         "C(",
    "branch_close":          ")",
}

RING_DIGIT = re.compile(r"[0-9]")


def compose(parts, label=""):
    """Join name-derived SMILES tokens, refusing a ring-closure collision.

    The documented failure of string concatenation in this project is two
    fragments reusing the same ring-closure digit, which silently FUSES them.
    So: no digit may appear in more than one fragment, and every digit that
    appears must appear an even number of times inside its own fragment.
    """
    carriers = [i for i, p in enumerate(parts) if RING_DIGIT.search(p)]
    seen = {}
    for i in carriers:
        for d in set(RING_DIGIT.findall(parts[i])):
            if d in seen:
                raise ValueError(
                    "compose(%s): ring-closure digit %r appears in fragment "
                    "%d and fragment %d. Concatenating them would FUSE the two "
                    "rings -- that is exactly how piperazine was bonded to "
                    "cholesterol (CLAUDE.md section 6). Refusing."
                    % (label, d, seen[d], i))
            seen[d] = i
            if parts[i].count(d) % 2:
                raise ValueError(
                    "compose(%s): fragment %d has an odd number of %r, so its "
                    "ring is not closed inside itself." % (label, i, d))
    return "".join(parts)


def anchor_A18B8C14():
    """(pentane-1,5-diylbis(azanediyl))bis(tetradecane-6,1-diyl)
       bis(2-hexyldecanoate)                                  -- SI p.S18-S19.

    Read strictly left to right off the name:
      * a 2-hexyldecanoate ester ...
      * ... of the 1-position of a tetradecane-6,1-diyl ...
      * ... whose 6-position carries an azanediyl nitrogen ...
      * ... and the two nitrogens are bridged by pentane-1,5-diyl.
    Written starting at C2 of one decanoate and ending at C2 of the other.
    """
    acyl = [TOK["branch_C_open"], TOK["hexyl"], TOK["branch_close"],
            "(", TOK["decanoate_C3_C10"], ")", TOK["carbonyl"],
            TOK["ester_O"]]
    arm_in = [TOK["tetradecane_C1_C5"], TOK["branch_C_open"],
              TOK["tetradecane_C7_C14"], TOK["branch_close"]]
    arm_out = [TOK["branch_C_open"], TOK["tetradecane_C7_C14"],
               TOK["branch_close"], TOK["tetradecane_C1_C5"]]
    acyl_rev = [TOK["ester_O"], TOK["carbonyl"], TOK["branch_C_open"],
                TOK["hexyl"], TOK["branch_close"], TOK["decanoate_C3_C10"]]
    return compose(acyl + arm_in + [TOK["azanediyl"]]
                   + [TOK["pentane_1_5_diyl"]] + [TOK["azanediyl"]]
                   + arm_out + acyl_rev, "A18B8C14")


def anchor_A2B1C1():
    """N1-(1-cyclopentylnonyl)-N3,N3-dimethylpropane-1,3-diamine -- SI p.S12.

    Written starting at the cyclopentane and ending at the dimethylamino
    terminus. This anchor is here BECAUSE it has a ring: it is the one that
    exercises compose()'s ring-digit guard on real content.
    """
    return compose([TOK["cyclopentyl"], TOK["branch_C_open"],
                    TOK["nonyl_C2_C9"], TOK["branch_close"],
                    TOK["azanediyl"], TOK["propane_1_3_diyl"],
                    TOK["dimethylamino"]], "A2B1C1")


ANCHORS = [
    ("A18B8C14", anchor_A18B8C14,
     "SI p.S18-S19: the compound name '(pentane-1,5-diylbis(azanediyl))"
     "bis(tetradecane-6,1-diyl)bis(2-hexyldecanoate)' and the structure drawn "
     "beneath it", "C65H130N2O4", 1004.0103),
    ("A2B1C1", anchor_A2B1C1,
     "SI p.S12: the compound name 'N1-(1-cyclopentylnonyl)-N3,N3-"
     "dimethylpropane-1,3-diamine' and the structure drawn beneath it",
     "C19H40N2", 297.3264),
]
PROTON = 1.00727646


def self_test():
    print("-- compose() refuses a ring-closure digit collision --")
    try:
        compose(["C1CCCC1", "C1CCCCC1"], "SEEDED collision")
    except ValueError as e:
        print("   FIRED: %s" % str(e).split("\n")[0])
    else:
        raise SystemExit("compose() accepted a digit collision. The guard is "
                         "decorative, so the anchor proves nothing.")

    print("-- compose() refuses an unclosed ring --")
    try:
        compose(["C1CCCC"], "SEEDED unclosed")
    except ValueError as e:
        print("   FIRED: %s" % str(e).split("\n")[0])
    else:
        raise SystemExit("compose() accepted an unclosed ring.")

    print("-- a one-carbon-short token must break the mass check --")
    good = anchor_A2B1C1()
    bad = compose([TOK["cyclopentyl"], TOK["branch_C_open"],
                   "CCCCCCC",                      # nonyl_C2_C9 minus one C
                   TOK["branch_close"], TOK["azanediyl"],
                   TOK["propane_1_3_diyl"], TOK["dimethylamino"]], "SEEDED")
    bm = Descriptors.ExactMolWt(Chem.MolFromSmiles(bad)) + PROTON
    d = (bm - 297.3264) * 1000
    print("   good %s" % Chem.MolToSmiles(Chem.MolFromSmiles(good)))
    print("   seeded one carbon short -> [M+H]+ %.4f, %+.1f mDa from the SI's "
          "printed calcd. FIRES." % (bm, d))
    if abs(d) < 2.0:
        raise SystemExit("the mass check did not fire on a missing carbon.")
    return 0


def main():
    if "--self-test" in sys.argv:
        return self_test()

    rows = []
    for name, fn, figure, formula, calcd in ANCHORS:
        smi = fn()
        m = Chem.MolFromSmiles(smi)
        if m is None:
            raise SystemExit("anchor %s did not parse: %s" % (name, smi))
        f = rdMolDescriptors.CalcMolFormula(m)
        mh = Descriptors.ExactMolWt(m) + PROTON
        print("%-10s %s" % (name, Chem.MolToSmiles(m)))
        print("           formula %-12s (SI neutral %-12s) "
              "[M+H]+ %10.4f  SI calcd %10.4f  d %+0.3f mDa  rings %d"
              % (f, formula, mh, calcd, (mh - calcd) * 1000,
                 Chem.rdMolDescriptors.CalcNumRings(m)
                 if hasattr(Chem.rdMolDescriptors, "CalcNumRings")
                 else rdMolDescriptors.CalcNumRings(m)))
        if f != formula:
            raise SystemExit("anchor %s formula %s != SI %s" % (name, f,
                                                                formula))
        if abs(mh - calcd) * 1000 > 2.0:
            raise SystemExit("anchor %s mass disagrees with the SI's printed "
                             "calcd." % name)
        rows.append({
            "lipid_name": name,
            "reference_SMILES": Chem.MolToSmiles(m),
            "source_figure": figure,
            "read_on": READ_ON,
            "read_by": READ_BY,
            "expected_mismatch": "",
            "note": ("Derived from the compound NAME printed in the SI, not "
                     "from the block drawing the enumerator reads, and "
                     "composed by plain string concatenation of name-derived "
                     "tokens -- no RWMol, no AddBond, no RemoveAtom, no atom "
                     "map. compose() refuses a ring-closure digit shared "
                     "between fragments, which is the documented failure mode "
                     "of this route. The alpha carbon(s) are unassigned "
                     "stereocentres in both the anchor and the enumeration: "
                     "the photochemical radical addition is not "
                     "stereoselective and the SI does not resolve them, so "
                     "C10's stereo comparison has nothing to compare."),
        })
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=[
            "lipid_name", "reference_SMILES", "source_figure", "read_on",
            "read_by", "expected_mismatch", "note"])
        w.writeheader()
        w.writerows(rows)
    print("\nwrote %s (%d anchors)" % (OUT, len(rows)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
