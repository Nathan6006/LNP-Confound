#!/usr/bin/env python3
"""BUILD the C10 anchors for namjh_39264000. Do not store a string -- derive it.

    python3 -I make_reference.py      ->  reference_structures.csv

C10 compares each enumerated structure against the structure the paper
specifies. It is the only Tier-1 check that tests TOPOLOGY against external
ground truth, and it only tests anything if the anchor was built by a DIFFERENT
CODE PATH than the enumerator. That independence has been asserted twice in
this project and was false twice -- F-013 (a byte copy with the same sha256,
while its note claimed independent fragment assembly) and F-037 (the file
written to fix F-013 re-copied the other half of the molecule).

HOW THIS ONE IS DIFFERENT, IN BOTH EVIDENCE AND MECHANISM.

  EVIDENCE.  The enumerator's numbers are MEASURED FROM A DRAWING. stage01
             traces Fig 3a's skeleton, counts Douglas-Peucker vertices, and
             reads ring locants off a hexagon's own vertex angles. Nothing
             below looks at a pixel. Every count here is read out of a
             REAGENT NAME or an IUPAC NAME printed as prose in the SI:

               "3,5-dihydroxybenzoic acid"            (SI 1.2.1)
               "glutaric anhydride"                   (SI 1.2.2)
               "2-hexyl-1-decanol"                    (SI 1.2.2)
               "2-octyl-1-dodecanol"                  (SI 1.2.2)
               "1-(2-((tert-butyldimethylsilyl)oxy)ethyl)piperazine" (1.2.3)
               "3-iodopropane-1,2-diol"               (SI 1.2.4)
               "O,O'-(3-(4-(2-((tert-butyldimethylsilyl)oxy)ethyl)piperazin-
                1-yl)propane-1,2-diyl)bis(2-hexyldecyl)diglutarate"  (1.2.5.1)
               "O,O'-(3-(4-(2-((3,5-bis((tert-butyldimethylsilyl)oxy)benzoyl)
                oxy)ethyl)piperazin-1-yl)propane-1,2-diyl)bis(2-hexyldecyl)
                diglutarate"                          (SI 1.2.5.1)

             The last two names state the whole connectivity: which oxygen of
             propane-1,2-diyl carries which ester, that the piperazine is
             attached at its 1-position with the hydroxyethyl arm at its
             4-position, and that the benzoyl is 3,5-disubstituted.

  MECHANISM. Plain string substitution into a template, traversed HEAD FIRST.
             There is no RWMol, no AddBond, no RemoveAtom, no CombineMols, no
             call into curation/lib `assemble`, and no import of stage01 or
             stage03. The enumerator traverses tail first and bonds atom-
             mapped fragments with assemble.attach.

THE TEMPLATE IS NOT A COPY OF THE ENUMERATION. The enumerator's canonical
output for FH-5-EB-5T-1 begins 'CCCCCCCCC(CCCCCC)COC(=O)...' -- a tail. The
template below begins 'O=C(c1...' -- the functional head. They are required to
canonicalise to the same InChIKey, which is the whole test; what must NOT be
shared is the route, and `check_anchor_independence.py` (C16) inspects THIS
file's composed template, because byte-identity of the EMITTED anchor is not
the signal.
"""
import csv
import os

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "reference_structures.csv")


# --------------------------------------------------------------------------
# Every integer below is read out of a NAME, with the reading written beside
# it. None of them comes from a measurement of a figure.
# --------------------------------------------------------------------------
def chain(n):
    """n carbons, written out. No molecule object is constructed."""
    return "C" * n


# "glutaric anhydride" -> glutaric acid is PENTANEdioic acid: 5 carbons, two of
# which are the carbonyls, so 3 methylenes between them.
GLUTARYL = "C(=O)" + chain(5 - 2) + "C(=O)"

# "2-hexyl-1-decanol" -> a DECAN-1-ol (10 carbons, OH on C1) carrying a HEXYL
# (6 carbons) on C2. As an ester alkyl, written from C1: CH2, then C2 bearing
# the hexyl branch, then C3..C10 = 8 more carbons.
HEXYLDECYL = "C" + "C(" + chain(6) + ")" + chain(10 - 2)

# "2-octyl-1-dodecanol" -> DODECAN-1-ol (12 carbons) carrying an OCTYL (8) on
# C2. From C1: CH2, C2 with the octyl branch, then C3..C12 = 10 more.
OCTYLDODECYL = "C" + "C(" + chain(8) + ")" + chain(12 - 2)

# "3,5-dihydroxybenzoic acid" -> benzoyl with OH at ring positions 3 and 5,
# i.e. both meta to the carbonyl. Written as an explicit six-membered aromatic
# ring: C1(acyl), C2, C3(OH), C4, C5(OH), C6.
DIHYDROXYBENZOYL = "O=C(c1cc(O)cc(O)c1)"

# "1-(2-((TBS)oxy)ethyl)piperazine" + "3-iodopropane-1,2-diol", then TBAF.
# The piperazine carries the 2-hydroxyethyl arm on one nitrogen and the
# propane-1,2-diol's C3 on the other. Ring digit 2 is used so that the aryl
# ring's digit 1 is never reused in the same string -- reusing a ring number
# across two fragments is how a piperazine once ended up fused into
# cholesterol (CLAUDE.md section 6).
TEMPLATE = "@FH@OCCN2CCN(CC2)CC(O@EBT@)CO@EBT@"


def anchor(fh_acyl, alkyl):
    """Substitute, and nothing else. Two str.replace calls."""
    ebt = GLUTARYL + "O" + alkyl
    return TEMPLATE.replace("@FH@", fh_acyl).replace("@EBT@", ebt)


def build():
    """Return [(lipid_name, SMILES, source_figure, read_on, read_by,
                expected_mismatch, note)]."""
    rows = []
    for name, alkyl, alkyl_name in (
            ("FH-5-EB-5T-1", HEXYLDECYL, "2-hexyl-1-decanol"),
            ("FH-5-EB-5T-2", OCTYLDODECYL, "2-octyl-1-dodecanol")):
        rows.append((
            name,
            anchor(DIHYDROXYBENZOYL, alkyl),
            "SI section 1.2.5.1 / 1.2.5.2 (IUPAC names of compounds 4-7) "
            "and SI Figure S9 (HRMS)",
            "2026-11-10",
            "derived from the SI's PROSE NAMES, not from a figure: "
            f"'3,5-dihydroxybenzoic acid' + 'glutaric anhydride' + "
            f"'{alkyl_name}' + the piperazine/propane-1,2-diol core, "
            "assembled head-first by plain string substitution into a "
            "template with no RWMol, no AddBond and no import of the "
            "enumerator",
            "",
            "The paper specifies this member completely: an IUPAC name for "
            "the protected precursor, a 1H and a 13C NMR spectrum (Figures "
            "S5-S8) and an HRMS peak (Figure S9). It is the only member with "
            "an external anchor of any kind."))
    rows.append((
        "FH-1 / FH-2 / FH-4 / FH-6 and EB-6 / EB-7 members",
        "",
        "none",
        "2026-11-10",
        "no anchor exists",
        "",
        "DELIBERATELY EMPTY so the ABSENCE is visible. 24 of the 36 enumerated "
        "members contain at least one block that the paper names in no "
        "synthetic procedure, no reagent list, no NMR assignment and no mass "
        "spectrum. Their only evidence is the Fig 3a drawing. C10 therefore "
        "tests 2 of 36 members against external ground truth, and the other "
        "34 rest on stage01's tracer, which is itself calibrated on the five "
        "blocks the SI does name."))
    return rows


if __name__ == "__main__":
    rows = build()
    with open(OUT, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["lipid_name", "reference_SMILES", "source_figure",
                    "read_on", "read_by", "expected_mismatch", "note"])
        w.writerows(rows)
    print(f"wrote {len(rows)} anchor row(s) to {OUT}")
    for r in rows:
        print(f"  {r[0]}: {r[1] or '(empty -- no anchor exists)'}")
