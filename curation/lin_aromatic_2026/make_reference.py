#!/usr/bin/env python3
"""Build the C10 anchor for `10A2-3B-EP10` from the paper, by a different route.

WHY THIS EXISTS. `reference_structures.csv` held this structure as a typed
literal with no generator, so C10 reported `anchor_may_be_a_copy`: a reference
copied out of the enumeration is byte-identical to an honestly re-derived one
after canonicalisation, and can therefore only pass. The anchor's VALUE turned
out to be correct -- a round-2 reviewer derived it from scratch and got the same
InChIKey, `JUAPTPKSYXADRC-UHFFFAOYSA-N` -- so this is not F-013 repeating. What
was missing was that the provenance claim had never been *executed*.

"A different code path" has to be run, not asserted (CLAUDE.md §6). So:

  * no import of `enumerate_library.py` or anything else in this folder;
  * no `RWMol`, `AddBond`, `RemoveAtom` or `InsertMol` anywhere -- the molecule
    is built by plain string concatenation of three fragments written from the
    reagents' chemical names;
  * RDKit is used ONLY to canonicalise and to weigh, never to edit.

Fragment assembly is mandated elsewhere in this project precisely because string
concatenation once fused a piperazine to cholesterol through a ring-closure
collision. That hazard does not apply here: there is exactly one ring system in
each fragment (a benzene, written `c1ccccc1`), the fragments are concatenated
and never nested, and the assert against the stored string would catch a
collision immediately. The enumerator uses templates and explicit ring-closure
renumbering, which is the code path this file deliberately does not share.

WHAT PINS EACH FEATURE, from the SI's precise synthesis (Fig S17 and pp. 29-32):

  spermine core, 3-4-3        traced off the Fig S17 final product
  one 3B arm per PRIMARY N    the ROUTE, not an assumption: P1 is bis-TFA
                              spermine on the PRIMARIES ([M+H]+ 395.19, and the
                              1H NMR integrates 4H next to the amide N -- it
                              would be 8H had the TFA gone on the secondaries);
                              P2/P3 put Boc on the secondaries (403.33,
                              18H Boc singlet); 3B is then added to P3, whose
                              only free N-H are the two primaries; P5 (503.36)
                              confirms exactly two 3B arms
  4 EP10 tails on the rest    final step P5 + EP10 (6 eq)
  3B opens at the terminal    Fig S17: exactly 3 carbons between the ether O and
  CH2                         the N, OH on the middle one
  EP10 = N-CH2-CH(OH)-C8H17   Fig S17 vertex trace of one tail
  no stereocentres asserted   Fig S17 draws no wedge on the final product
                              (unlike Fig 1b's 5B), and 3B is achiral

THE ASSERT THAT MATTERS is the HRMS one. Checking the built string against the
stored string only proves we are self-consistent; checking its [M+H]+ against
the 1127.9648 the SI prints makes the anchor testable against the PAPER. That is
the whole point of a C10 anchor, and it is what the old typed literal could not
do.

    python3 -I curation/lin_aromatic_2026/make_reference.py
"""
import csv
import os
import sys

from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors as D, rdMolDescriptors as rd

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "reference_structures.csv")

LIPID = "10A2-3B-EP10"

# --- fragments, written from the reagents' chemical names -------------------
# 1,2-epoxydecane opened by an amine at the terminal (unsubstituted) C1, giving
# a secondary carbinol at C2 and an octyl remainder:  N-CH2-CH(OH)-C8H17
EP10 = "CC(O)CCCCCCCC"
# glycidyl phenyl ether (3B), likewise opened at C1:  N-CH2-CH(OH)-CH2-O-Ph
B3 = "CC(O)COc1ccccc1"

# Spermine (10A), N1-(3-aminopropyl)butane-1,4-diamine: a 3-4-3 methylene
# backbone over four nitrogens. One 3B on each PRIMARY nitrogen (N1 and N4);
# the four remaining N-H -- two on the primaries, one on each secondary -- each
# take an EP10 tail. Written starting from one 3B arm so that no fragment is
# ever nested inside another.
PROTON = 1.007276466


def build():
    return ("c1ccccc1OCC(O)C"                 # 3B arm on primary N1
            "N(" + EP10 + ")CCC"              # N1 + its EP10 tail, 3 CH2 to N2
            "N(" + EP10 + ")CCCC"             # N2 (secondary) + tail, 4 CH2
            "N(" + EP10 + ")CCC"              # N3 (secondary) + tail, 3 CH2
            "N(" + EP10 + ")" + B3)           # N4 + its tail and 3B arm


def main():
    smi = build()
    m = Chem.MolFromSmiles(smi)
    if m is None:
        sys.exit("make_reference: the assembled SMILES does not parse")
    canon = Chem.MolToSmiles(m)

    formula = rd.CalcMolFormula(m)
    mh = D.ExactMolWt(m) + PROTON
    key = Chem.MolToInchiKey(m)

    # 1. against the PAPER: SI Fig S21 HRMS, [M+H]+ calc 1127.9648
    assert abs(mh - 1127.9648) < 0.0015, f"[M+H]+ {mh:.4f} != 1127.9648 (SI)"
    assert formula == "C68H126N4O8", formula
    # 2. the structural facts the route pins, checked on the built molecule
    assert sum(1 for a in m.GetAtoms() if a.GetSymbol() == "N") == 4
    assert rd.CalcNumAromaticRings(m) == 2, "two 3B arms, so two phenyl rings"
    assert sum(1 for a in m.GetAtoms()
               if a.GetSymbol() == "O" and a.GetTotalNumHs() == 1) == 6, \
        "six secondary carbinols: 4 EP10 + 2 3B"
    assert not m.GetSubstructMatches(Chem.MolFromSmarts("[NX3;H1,H2]")), \
        "every N-H must be capped"
    # 3. no stereochemistry is asserted anywhere
    assert "@" not in canon, "Fig S17 draws no wedge on the final product"

    note = (
        "Final product of the PRECISE protection/deprotection synthesis, so the "
        "regiochemistry is pinned by the route rather than assumed: one 3B arm "
        "on EACH primary nitrogen, four EP10 tails on the remaining N-H. "
        "Rebuilt 2026-11-07 by make_reference.py, which concatenates three "
        "fragments written from the reagents' chemical names and uses no "
        "RWMol/AddBond/RemoveAtom and no import from this folder -- so the "
        "different-code-path claim is executed, not asserted (F-008 round). The "
        "anchor is checked against the SI's own HRMS ([M+H]+ 1127.9648), which "
        "makes it testable against the paper and not merely against us.")

    rows = []
    if os.path.exists(OUT):
        rows = list(csv.DictReader(open(OUT)))
    hdr = ["lipid_name", "reference_SMILES", "source_figure", "read_on",
           "read_by", "note"]
    rec = {"lipid_name": LIPID, "reference_SMILES": canon,
           "source_figure": "Fig S17 (+ Fig S21 HRMS)",
           "read_on": "2026-11-07",
           "read_by": "make_reference.py (fragment strings from chemical names)",
           "note": note}
    prev = next((r for r in rows if r["lipid_name"] == LIPID), None)
    if prev is not None:
        # The anchor's value must not move. If it does, say so loudly rather
        # than overwrite -- a silently changing C10 reference is worse than none.
        if Chem.MolToSmiles(Chem.MolFromSmiles(prev["reference_SMILES"])) != canon:
            sys.exit("make_reference: REBUILT ANCHOR DIFFERS FROM THE STORED ONE\n"
                     f"  stored : {prev['reference_SMILES']}\n"
                     f"  rebuilt: {canon}\n"
                     "  Resolve this before writing anything.")
        rows = [rec if r["lipid_name"] == LIPID else r for r in rows]
    else:
        rows.append(rec)

    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=hdr, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    print(f"{LIPID}")
    print(f"  SMILES   {canon}")
    print(f"  formula  {formula}   [M+H]+ {mh:.4f}  (SI prints 1127.9648)")
    print(f"  InChIKey {key}")
    print(f"  wrote    {OUT}")


if __name__ == "__main__":
    main()
