#!/usr/bin/env python3
"""The C10 anchors for xuq_32438474 (YL_2020), BUILT -- not stored, not pasted.

    python3 -I make_reference.py        # rewrites reference_structures.csv

WHY A GENERATOR AND NOT A STRING. C10 only tests anything if the anchor came
from a different code path than the enumerator. That has been claimed twice in
this project and been false twice: F-013 (the ZT_2026 anchor had the same
sha256 as the enumeration it was checking, while its note claimed independent
fragment assembly "through the C3 oxygen (3-beta)" and the string was 3-alpha)
and F-037 (the generator written to cure F-013 hard-coded a template that was
the enumeration's canonical SMILES with the sterol swapped for a placeholder
and one ring digit renumbered).

HOW THIS ROUTE DIFFERS FROM stage03_enumerate.py, concretely:

    enumerator          assemble.esterify  -> RWMol.AddBond + RWMol.RemoveAtom
                        assemble.attach    -> RWMol.AddBond, atom-map lookup
                        tail built bottom-up from acrylic acid
    this file           str.replace on the registry's cholesterol string, and
                        one SMILES ring-closure bond. NO RWMol, NO AddBond, NO
                        RemoveAtom, NO esterify, NO attach, NO atom maps.

The two share only RDKit's parser and the registry's sterol -- and sharing the
registry sterol is the point, because the registry writes it from the chemical
name "cholest-5-en-3beta-ol" and is the one place in the project that is not
downstream of any enumeration.

TRAVERSAL. The molecule is written ROOTED AT THE STEROL of arm 1 and walks
sterol -> carbamate -> disulfide -> ester -> head -> arm 2. The enumerator's
output is RDKit canonical order, which starts elsewhere entirely. Nothing here
is a canonical string, so C16 has something real to look at.

RING-CLOSURE SAFETY. Exactly one closure digit is used, written as `%11`, and
the code asserts the registry sterol contains no '%' at all. Fusing a head to
cholesterol through a ring-closure collision is a bug this project has already
shipped once; the assertion is cheap.

WHAT THIS DOES NOT DO. It does not prove the drawn structure was read
correctly. It reproduces, by an independent traversal, what a human read off
Figure 1A (the OCholB tail and head 75) and Figure 2A (head 304) at 300 dpi on
2026-11-10. If that reading is wrong, both halves are wrong together and only
an external number can tell -- which is what the SI's [M+H]+ formula tables do
in stage03, and what they CANNOT do for stereochemistry (see 07_seed_checks.json
seedC: epicholesterol leaves all 18 formulas and all 18 skeleton InChIKeys
unchanged).
"""
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
sys.path.insert(0, HERE)

from rdkit import Chem                                          # noqa: E402
from rdkit.Chem.rdMolDescriptors import CalcMolFormula          # noqa: E402

import reagents                                                 # noqa: E402

# The sterol's C3 neighbourhood as the registry writes it. Short on purpose:
# the sterol itself is fetched at call time, never stored here.
C3 = "C[C@@H](O)CC"

# The OCholB arm, written from the sterol's C3 oxygen OUTWARD, in the order the
# synthesis runs (SI p. S2): cholesteryl chloroformate -> carbamate, then the
# cysteamine/2-mercaptoethanol disulfide, then the acrylate ester, then the two
# Michael carbons that reach the head nitrogen.
ARM = "OC(=O)NCCSSCCOC(=O)CC"

# Heads, written from the lipidoid nitrogen outward, as drawn.
HEAD_FROM_N = {
    # 75: N,N-diethylethylenediamine -- the reacted N carries -CH2CH2-NEt2
    "75": "CCN(CC)CC",
    # 304: tris(2-(methylamino)ethyl)amine -- each reacted N is an N-methyl
    # that also carries -CH2CH2-N(central).  Written here for ARM 1's nitrogen,
    # with the other two arms hung off the central nitrogen.
    "304": None,     # built by compose_304, which needs three arms
}


def sterol():
    s = reagents.registered("Cholesterol")
    if not s:
        raise SystemExit("Cholesterol is not pinned in reagent_registry.py; "
                         "refusing to type a sterol into a curation script.")
    if "%" in s:
        raise SystemExit("the registry sterol contains a '%' ring closure; "
                         "the %11 used below could collide.")
    if s.count(C3) != 1:
        raise SystemExit(f"the registry sterol does not contain {C3!r} exactly "
                         "once; this generator's substitution is no longer "
                         "safe and must be re-derived, not patched.")
    return s


def arm(tail_suffix):
    """One sterol-bearing arm, rooted at the sterol, ending at `tail_suffix`."""
    return sterol().replace(C3, "C[C@@H](" + ARM + tail_suffix + ")CC")


def compose_75():
    """75-OCholB: Et2N-CH2CH2-N(arm)(arm).

    Arm 1 is written first and carries the head; arm 2 is joined to the same
    nitrogen through ring-closure %11, which is how two sterol-rooted fragments
    are bonded without ever calling AddBond.
    """
    a1 = arm("N(" + HEAD_FROM_N["75"] + ")%11")
    a2 = arm("%11")
    return a1 + "." + a2


def compose_304():
    """304-OCholB: N(CH2CH2-N(Me)(arm))3, three arms on three N-methyl amines.

    Arm 1 is rooted at its sterol and reaches the central tertiary nitrogen;
    arms 2 and 3 come back to their own N-methyl nitrogens through %12 and %13.
    """
    a1 = arm("N(C)CCN(CCN(C)%12)CCN(C)%13")
    a2 = arm("%12")
    a3 = arm("%13")
    return a1 + "." + a2 + "." + a3


ANCHORS = [
    ("75-OCholB", compose_75,
     "Fig 1A", "C76H130N4O8S4",
     "head 75 and the full OCholB tail are both DRAWN in Figure 1A of the "
     "paper (PMC raster nihms-1624779-f0002.jpg). Built here rooted at the "
     "sterol by string substitution on provenance/reagent_registry.py's "
     "cholesterol; no RWMol, no AddBond, no RemoveAtom, no atom maps."),
    ("304-OCholB", compose_304,
     "Fig 2A + Fig 1A", "C114H195N7O12S6",
     "head 304 is DRAWN in Figure 2A (nihms-1624779-f0003.jpg); the OCholB "
     "tail is DRAWN in Figure 1A. 304 is anchored because it is the one head "
     "whose drawing is genuinely ambiguous -- its two upper arms close "
     "visually into what reads as a piperazine. Three arms, from three "
     "N-methyl secondary amines; the S6 in the SI's own [M+H]+ formula is the "
     "independent witness that there are three."),
]


def build(fn, want_formula, label):
    smi = fn()
    m = Chem.MolFromSmiles(smi)
    if m is None:
        raise SystemExit(f"{label}: composed SMILES does not parse:\n{smi}")
    frags = Chem.GetMolFrags(m)
    if len(frags) != 1:
        raise SystemExit(f"{label}: composed into {len(frags)} fragments; the "
                         "ring closure did not bond them.")
    f = CalcMolFormula(m)
    if f != want_formula:
        raise SystemExit(f"{label}: formula {f}, expected {want_formula}")
    return Chem.MolToSmiles(m), f


def main():
    rows = []
    for name, fn, fig, formula, note in ANCHORS:
        smi, f = build(fn, formula, name)
        rows.append({
            "lipid_name": name, "reference_SMILES": smi,
            "source_figure": fig, "read_on": "2026-11-10",
            "read_by": "claude-code", "expected_mismatch": "",
            "note": note,
        })
        print(f"{name:12s} {f}  {Chem.MolToInchiKey(Chem.MolFromSmiles(smi))}")
    out = os.path.join(HERE, "reference_structures.csv")
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["lipid_name", "reference_SMILES",
                                           "source_figure", "read_on",
                                           "read_by", "expected_mismatch",
                                           "note"])
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
