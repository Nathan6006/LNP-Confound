#!/usr/bin/env python3
"""Build the C10 anchors, from what the paper DRAWS, by a different route.

    python3 -I make_reference.py      ->  reference_structures.csv

WHY THIS FILE IS WRITTEN THE WAY IT IS. A C10 anchor only tests anything if it
was built by a code path that is different IN FACT, not in the comment. That
claim has been false twice in this project: F-013 was a byte copy of the
enumeration with a note claiming independent fragment assembly, and F-037 was
the file written to FIX F-013, whose amine half was still the enumeration's
canonical SMILES with the sterol swapped for a placeholder. So:

  * this file does not import paper_library, curation/lib, or anything else
    the enumerator uses -- only rdkit, to canonicalise and compare;
  * there is no RWMol, AddBond, RemoveAtom or atom map anywhere in it;
  * the SMILES below is a DIFFERENT TRAVERSAL of the same molecule. The
    enumerator writes the dye head first and walks outward to the arms. This
    starts at arm A's terminal primary amine, walks inward through the linker,
    the sulfonamide, round the benzenesulfonyl ring, over the benzylic CH2 to
    the toluidine ring, out to the central methine, and only then to the
    quinoid ring and the ethoxyphenyl head, before going back down arm B.
    Every ring is entered at a different atom and the ring-closure digits are
    assigned in a different order.

Byte-identity of the EMITTED anchor is not the signal and testing for it is a
trap: a correct re-derivation canonicalises to a correct enumeration. What is
checkable is that the TEMPLATE here is not the enumerator's template, and
`python3 -I make_reference.py --show-template` prints it so a reviewer can
diff the two by eye.

READ FROM: Supplementary Figure 1 (SI p2) for CBB21 and Supplementary Figure 3
(SI p4) for CBBZ0, rendered at 400 dpi and looked at -- renders/si_p2.png,
renders/crop_core.png, renders/crop_linkers.png, renders/crop_fig3.png.
"""
import csv
import os
import sys

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")
HERE = os.path.dirname(os.path.abspath(__file__))

# --- read off the drawing, arm-first ---------------------------------------
# the K1 linker as drawn: H2N-(CH2)3-N(CH3)-(CH2)3-NH2, the far amine
# sulfonylated. Written starting from the FREE amine.
_LINKER_IN = "NCCCN(C)CCCN"
# the benzenesulfonyl ring, entered at the sulfur-bearing carbon and left at
# the meta carbon that carries the benzylic CH2
_SULFOPHENYL = "S(=O)(=O)c1cccc(c1)C"
# toluidine ring A, entered at its NITROGEN-bearing carbon; the methyl is drawn
# on the carbon three bonds round from the nitrogen and one bond from the
# central methine, so it is at position 3 counting from the nitrogen
_TOLUIDINE_FROM_N = "c2cc(C)c({centre})cc2"
# the quinoid ring and the 4-ethoxyphenyl head, reached LAST
_QUINOID_AND_HEAD = "C(=C3C=CC(=Nc4ccc(O{head})cc4)C=C3)"
# toluidine ring B, entered from the CENTRAL METHINE instead, so the methyl is
# now ortho to the entry atom and the nitrogen para to it
_TOLUIDINE_FROM_C = "c5c(C)cc(N({armB})CC)cc5"


def template(head, arm_b_tail):
    """Compose CBB21 / CBBZ0, arm A first, head last."""
    # arm B's benzenesulfonyl ring is walked the OTHER way round from arm
    # A's and from the enumerator's, so the two files share no long run.
    arm_b = "Cc6cc(S(=O)(=O)" + arm_b_tail + ")ccc6"
    tolB = _TOLUIDINE_FROM_C.format(armB=arm_b)
    centre = _QUINOID_AND_HEAD.format(head=head) + tolB
    return (_LINKER_IN + _SULFOPHENYL + "N(CC)"
            + _TOLUIDINE_FROM_N.format(centre=centre))


def template_z(head):
    """CBBZ0: both sulfonyl groups are free sulfonic acids, so there is no
    linker to start from -- walk in from one sulfonic acid instead."""
    arm_b = "Cc6cc(S(=O)(=O)O)ccc6"
    tolB = _TOLUIDINE_FROM_C.format(armB=arm_b)
    centre = _QUINOID_AND_HEAD.format(head=head) + tolB
    return ("OS(=O)(=O)c1cccc(c1)CN(CC)"
            + _TOLUIDINE_FROM_N.format(centre=centre))


ANCHORS = [
 dict(lipid_name="CBB21",
      smiles=template("CC", _LINKER_IN),
      source_figure="Supplementary Figure 1, SI p2 (the left-hand structure, "
                    "drawn in full)",
      expected_mismatch="",
      note="CBB21 is the only CBB-series compound the paper draws completely; "
           "CBB1-CBB20 are drawn as Markush cells, which stage07_qc.py P10 "
           "checks instead. Read arm-first from the drawing: the methyl sits "
           "on the ring carbon adjacent to the central methine and three bonds "
           "from the N(ethyl)(benzyl), and the sulfonamide is meta to the "
           "benzylic CH2."),
 dict(lipid_name="CBBZ0",
      smiles=template_z("CC#C"),
      source_figure="Supplementary Figure 3, SI p4 (the product structure, "
                    "drawn in full)",
      expected_mismatch="charge_state_only_figure_draws_the_cation",
      note="SI Fig 3 draws this as the N-protonated cation with both sulfonic "
           "acids as HO-SO2-. The anchor is the NEUTRAL free acid -- the same "
           "connectivity, deprotonated at the imine -- because that is what "
           "the CBBZ IUPAC names state and what the [M+H]+ HRMS of CBBZ1-"
           "CBBZ10 confirms. Connectivity is unaffected, so the skeleton must "
           "match exactly; only the protonation differs from the drawing."),
]


def main():
    if "--show-template" in sys.argv:
        print("CBB21 :", template("CC", "NCCCN(C)CCCN"))
        print("CBBZ0 :", template_z("CC#C"))
        return 0
    rows = []
    for a in ANCHORS:
        m = Chem.MolFromSmiles(a["smiles"])
        if m is None:
            raise SystemExit("anchor did not parse: " + a["lipid_name"])
        rows.append(dict(lipid_name=a["lipid_name"],
                         reference_SMILES=Chem.MolToSmiles(m),
                         source_figure=a["source_figure"],
                         read_on="2026-11-10",
                         read_by="curation agent, from a 400 dpi render of the "
                                 "SI page, read by eye",
                         expected_mismatch=a["expected_mismatch"],
                         note=a["note"]))
    hdr = ["lipid_name", "reference_SMILES", "source_figure", "read_on",
           "read_by", "expected_mismatch", "note"]
    with open(os.path.join(HERE, "reference_structures.csv"), "w",
              newline="") as f:
        w = csv.DictWriter(f, fieldnames=hdr)
        w.writeheader()
        for r in rows:
            w.writerow(r)

    # --- compare to the enumeration, and SAY SO -----------------------------
    enum = {r["compound"]: r["canonical_SMILES"]
            for r in csv.DictReader(open(os.path.join(HERE,
                                                      "03_enumerated.csv")))}
    rc = 0
    for r in rows:
        ours = enum.get(r["lipid_name"], "")
        a_mol, b_mol = (Chem.MolFromSmiles(r["reference_SMILES"]),
                        Chem.MolFromSmiles(ours) if ours else None)
        if b_mol is None:
            print("%-8s NO ENUMERATED STRUCTURE" % r["lipid_name"]); rc = 1
            continue
        ka, kb = Chem.MolToInchiKey(a_mol), Chem.MolToInchiKey(b_mol)
        same_skel = ka.split("-")[0] == kb.split("-")[0]
        print("%-8s skeleton %s   full %s   anchor_bytes==enum_bytes %s"
              % (r["lipid_name"],
                 "MATCH" if same_skel else "MISMATCH",
                 "match" if ka == kb else "differ",
                 r["reference_SMILES"] == ours))
        if not same_skel:
            rc = 1
    return rc


if __name__ == "__main__":
    sys.exit(main())
