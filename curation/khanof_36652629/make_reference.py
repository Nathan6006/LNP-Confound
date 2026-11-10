"""Build the C10 reference anchors for PMID 36652629, by a DIFFERENT code path.

WHY THIS FILE EXISTS AT ALL. A C10 anchor only tests something if it was built
by a route that does not share its reasoning with the enumerator. That claim
has been made twice in this project and been false twice: F-013 (a byte copy of
the enumeration, annotated as independent fragment assembly) and F-037 (the
file written to FIX F-013, whose amine half was the enumeration's canonical
SMILES with the sterol swapped for a '%'). So the difference is spelled out
here and it is mechanical, not asserted:

  the enumerator (stage01_enumerate.py)   this file
  ------------------------------------    ---------------------------------
  RDKit RWMol + CombineMols + AddBond      plain Python string concatenation
  atom-mapped attachment atoms             no atom maps at all
  built CORE first, then arm, then tails   written TAIL first, outside in
  tails added one nitrogen at a time       whole arm written as one literal
  imports curation/lib assemble            imports no project module

NO RDKIT MOLECULE IS EDITED HERE. RDKit is used only to parse the finished
string and print an InChIKey, which is the comparison, not the construction.

WHY STRING ASSEMBLY IS SAFE HERE SPECIFICALLY. String concatenation was banned
after it fused piperazine to cholesterol through a ring-closure digit
collision. That hazard needs rings. Every structure in this paper is ACYCLIC --
zero ring-closure digits in any of the eight SMILES, which the script asserts
below before writing anything. Where the hazard cannot occur, the route that
shares least with the enumerator is the better anchor.

WHAT WAS READ, AND FROM WHERE. The drawn structures, at 600 dpi:
  Fig S1A (SI p2)  C3-K2-E14, C3-K2-A14
  Fig S2A (SI p3)  C3-K2-E10, C3-K2-E15
  Fig S3A (SI p3)  C3-K3-E10
  Fig S4A (SI p4)  OC2-K3-E10, OC3-K3-E10
  Fig S7B (SI p7)  OC2-K2-A10
Each drawing shows the core amine, two propanamide arms, the bridge, and four
tails drawn with a bracket-and-subscript repeat. The subscript is what carries
the chain length; a bond count cannot, because the chain is abbreviated.

STEREOCHEMISTRY IS DELIBERATELY ABSENT. The carbinols come from a racemic
1,2-epoxyalkane, so the compound is a stereoisomer mixture and no single
configuration is right. The anchors are flat, like the enumeration, and C10's
stereo comparison is therefore not informative here -- said plainly rather than
left to look like agreement.
"""
import os, sys, csv

HERE = os.path.dirname(os.path.abspath(__file__))

# --- fragments, written out by hand from the drawings ----------------------
# A tail, written OUTWARD from its terminal methyl, the way the eye reads the
# drawing: the terminal CH3, the (CH2)n of the bracket, the carbinol, the CH2
# that meets the nitrogen.
def tail_E(total_carbons):
    """Epoxide-opened tail CH3-(CH2)n-CH(OH)-CH2- , n = total - 3.

    Returned as a PAIR, because a SMILES fragment has a direction and this
    tail is not symmetric: `out` is written ending at the attachment carbon
    (for the branch that runs INTO a nitrogen) and `into` is written starting
    from it (for a branch hanging OFF a nitrogen). Writing one where the other
    belongs puts the nitrogen on the far end of the chain, past the carbinol.
    """
    n = total_carbons - 3
    return ("C" * (n + 1) + "C(O)C", "CC(O)" + "C" * (n + 1))


def tail_A(total_carbons):
    """Plain CH3-(CH2)n- chain. Direction-symmetric, so both forms are equal
    -- which is exactly why the first version of this file passed on
    C3-K2-A14 and failed on all six epoxide lipids."""
    s = "C" * total_carbons
    return (s, s)


# A whole arm, written as one literal. The molecule is written as a single
# chain through the core nitrogen, so the two arms run in OPPOSITE directions:
# the first tail -> core, the second core -> tail.
#
# Direction matters TWICE here, and both were got wrong on the way to this
# version. Writing both arms outward gave a string with the SAME MOLECULAR
# FORMULA that parsed, canonicalised and carried no radical -- and bonded the
# core nitrogen to a TAIL instead of to the propanoyl. Fixing that left the
# tails themselves reversed inside the branches, which put each nitrogen past
# the carbinol at the far end of the chain; the formula was again identical,
# and C3-K2-A14 still passed because a plain alkyl tail reads the same in
# either direction. Nothing but the InChIKey comparison saw either one. That
# is the case for C10 in two lines: mass cannot see connectivity, and a
# formula check cannot either.
def arm_to_core(tail, bridge_carbons):
    """tail, tail, bridge, amide N, carbonyl, CH2CH2 -> meets the core N."""
    out, into = tail
    return out + "N(" + into + ")" + "C" * bridge_carbons + "NC(=O)CC"


def arm_from_core(tail, bridge_carbons):
    """core N -> CH2CH2, carbonyl, amide N, bridge, tail, tail."""
    _out, into = tail
    return "CCC(=O)N" + "C" * bridge_carbons + "N(" + into + ")" + into


# The core amine's substituent, written ATTACHMENT-FIRST: the nitrogen bonds
# to the first atom. Third instance of the direction class in this one file --
# "OCC" would hang the core nitrogen off the OXYGEN, giving N-O-CH2CH3 in
# place of N-CH2CH2-OH. Same atoms, same formula, same exact mass; a different
# molecule, and an amine-N-oxide ether at that.
CORE = {"C3": "CCC", "OC2": "CCO", "OC3": "CCCO"}


LIPIDS = {
    #  name          core   bridge  tail
    "C3-K2-E14":  ("C3",  2, tail_E(14), "Fig S1A"),
    "C3-K2-A14":  ("C3",  2, tail_A(14), "Fig S1A"),
    "C3-K2-E15":  ("C3",  2, tail_E(15), "SI p11 scheme"),
    "C3-K2-E10":  ("C3",  2, tail_E(10), "Fig S2A"),
    "C3-K3-E10":  ("C3",  3, tail_E(10), "Fig S3A"),
    "OC2-K3-E10": ("OC2", 3, tail_E(10), "Fig S4A"),
    "OC3-K3-E10": ("OC3", 3, tail_E(10), "Fig S4A"),
    "OC2-K2-A10": ("OC2", 2, tail_A(10), "Fig S7B"),
}

# Drawn in this paper but NOT synthesised in it. Their structures come from
# provenance/reagent_registry.py (check C11) and were confirmed against the
# corpus on all three keys in 04_corpus_overlap.json. No anchor was read off
# the drawing, so the gap is recorded rather than hidden.
NO_ANCHOR = {
    "SM-102": ("Fig S1A", "Drawn in Fig S1A but not made in this paper. "
                          "Structure taken from provenance/reagent_registry "
                          "(C11) and matched to the corpus on full InChIKey, "
                          "skeleton and name in AC_2025/WT_2026 (vitro) and "
                          "SB_2024/SS_2018/SY_2025/XD_2026 (vivo). The "
                          "drawing was NOT independently traced."),
    "DLin-MC3-DMA": ("none", "Named on the Fig S7A axis as 'MC3' and in "
                             "Table S2 rows 25-27; never drawn in the SI. "
                             "Structure from the registry (C11); corpus "
                             "match on all three keys across 5 vitro and 6 "
                             "vivo experiments."),
    "ALC-0315": ("none", "Appears ONLY in the Fig S11 pKa row (6.27). Never "
                         "drawn, and absent from Table S2 entirely, so this "
                         "paper formulated it without listing it. No value "
                         "row is emitted for it."),
}


def build(name):
    core, bridge, tl, _fig = LIPIDS[name]
    return (arm_to_core(tl, bridge) + "N(" + CORE[core] + ")"
            + arm_from_core(tl, bridge))


def main():
    from rdkit import Chem, RDLogger
    from rdkit.Chem import rdMolDescriptors
    RDLogger.DisableLog("rdApp.*")

    enum = {r["IL_name"]: r for r in
            csv.DictReader(open(os.path.join(HERE, "03_enumerated.csv")))}

    rows, agree = [], 0
    for name, (_c, _b, _t, fig) in LIPIDS.items():
        smi = build(name)
        assert not any(ch.isdigit() for ch in smi), (
            f"{name}: the hand-written anchor contains a ring-closure digit. "
            "Every structure in this paper is acyclic, so a digit here means "
            "the string is wrong -- and a ring-closure collision is the exact "
            "failure that banned string assembly in the first place.")
        m = Chem.MolFromSmiles(smi)
        assert m is not None, f"{name}: anchor does not parse"
        ours = Chem.MolFromSmiles(enum[name]["IL_SMILES"])
        same = Chem.MolToInchiKey(m) == Chem.MolToInchiKey(ours)
        agree += same
        f_anchor = rdMolDescriptors.CalcMolFormula(m)
        print(f"  {name:12s} {f_anchor:14s} SI {enum[name]['formula_SI']:14s} "
              f"anchor==enumeration: {same}   "
              f"byte-identical: {smi == enum[name]['IL_SMILES']}")
        rows.append(dict(
            lipid_name=name,
            reference_SMILES=smi,
            source_figure=fig,
            read_on="2026-11-10",
            read_by="curation agent, 600 dpi render of the SI page",
            expected_mismatch="",
            note=("Hand-written from the drawn structure, tail-first, by "
                  "plain string concatenation; no RDKit editing, no atom "
                  "maps, no curation/lib import -- a different code path from "
                  "stage01_enumerate.py, which builds core-first with "
                  "assemble.attach. Flat: the four carbinols are racemic, so "
                  "C10's stereo comparison carries no information here.")))

    for name, (fig, why) in NO_ANCHOR.items():
        rows.append(dict(lipid_name=name, reference_SMILES="",
                         source_figure=fig, read_on="2026-11-10",
                         read_by="not read",
                         expected_mismatch="no_anchor_read", note=why))

    out = os.path.join(HERE, "reference_structures.csv")
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=[
            "lipid_name", "reference_SMILES", "source_figure", "read_on",
            "read_by", "expected_mismatch", "note"])
        w.writeheader(); w.writerows(rows)
    print(f"\nanchors written: {len(LIPIDS)}; "
          f"agree with enumeration: {agree}/{len(LIPIDS)}; "
          f"declared gaps: {len(NO_ANCHOR)}")
    print(f"wrote {out}")
    return 0 if agree == len(LIPIDS) else 1


if __name__ == "__main__":
    sys.exit(main())
