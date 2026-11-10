#!/usr/bin/env python3
"""Build the C10 anchor for this folder: 306Oi10, as Figure 1A DRAWS it.

WHAT AN ANCHOR IS FOR. C10 compares each enumerated structure against the
structure the paper draws. Its entire value rests on the anchor having been
built by a DIFFERENT CODE PATH than the enumerator -- copy the enumeration in
and C10 compares a value to itself and can only pass. That has happened twice
in this project (F-013, a byte copy with a note claiming otherwise; F-037, the
file written to fix F-013 re-copying the other half), which is why C16 now
inspects the GENERATOR'S template rather than the emitted string.

WHY 306Oi10 IS THE RIGHT ANCHOR HERE. Figure 1A does not merely draw a product;
it draws the whole reaction -- the amine `306`, "+ 4" of the acrylate `Oi10`,
"3 days / 90 C", and the full structure of `306Oi10` with every bond. It is the
only place in the paper where the substitution COUNT is drawn rather than
implied, and the substitution count is the one structural claim that molecular
formula and exact mass cannot test (they are consequences of it). So this
anchor tests exactly what nothing else can.

HOW THIS ROUTE DIFFERS FROM stage03_enumerate.py, concretely:

  stage03            RDKit. Starts from the head SMILES `NCCCN(C)CCCN`, finds
                     N-H by valence, tags the nitrogen with an atom map, and
                     forms each bond with assemble.attach (CombineMols +
                     AddBond + Sanitize). The tail arrives as an ester built by
                     assemble.esterify from a map-tagged propanoic acid.
  this file          NO RDKit ASSEMBLY AT ALL. No RWMol, no AddBond, no
                     RemoveAtom, no CombineMols, no reaction. Four fragment
                     STRINGS, written out from the drawing, substituted into a
                     template by str.replace, and parsed once at the end only
                     to canonicalise and to compute the InChIKey.

  The traversal is also deliberately different: the template is written
  TAIL-FIRST, starting at the isopropyl terminus of one 8-methylnonyl chain and
  walking inward to the central N-methyl -- the order a person reading
  Figure 1A left to right would write. RDKit's canonical order is not that.

WHAT WAS READ OFF THE FIGURE, ONCE:
  * head 306 is drawn H2N-CH2CH2CH2-N(CH3)-CH2CH2CH2-NH2: two primary amines,
    one tertiary N-methyl, three-carbon arms;
  * FOUR tails are drawn, one on each N-H, two per primary nitrogen;
  * each tail is drawn -CH2CH2-C(=O)-O-CH2(CH2)6CH(CH3)2, i.e. the Michael
    adduct of an acrylate ester of 8-methylnonan-1-ol (the terminal isopropyl
    is drawn explicitly in both Figure 1A and Figure 1C);
  * the central nitrogen carries NO tail, as it must, having no N-H.

LIMIT, STATED: `Oi10` is sold as "isodecyl acrylate", a branched-ISOMER
MIXTURE. The 8-methylnonyl branch is what the authors drew, not what a vendor
guarantees, so this anchor pins the CONNECTIVITY of the lipidoid and the
substitution count, and inherits the drawing's claim about the branch. Every
Oi10 row in 03_enumerated.csv carries structure_assumed = 1 for that reason.

Run: python3 -I make_reference.py
"""
import csv
import hashlib
import os

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

HERE = os.path.dirname(os.path.abspath(__file__))

# --- fragments, written from Figure 1A ------------------------------------
# 8-methylnonyl, written from its FAR end inward: (CH3)2CH- then six CH2 then
# the CH2 that bears the ester oxygen.
ISODECYL_OUT = "CC(C)CCCCCCC"
# the same chain written the other way, O-first, for the three inward arms
ISODECYL_IN = "CCCCCCCC(C)C"
# the Michael-adduct linker: -O-C(=O)-CH2-CH2- toward the nitrogen
ARM_OUT = "OC(=O)CC"          # from the chain's CH2-O toward N
ARM_IN = "CCC(=O)O"           # from N toward the chain's O

# The molecule, tail-first. {A} is the outgoing chain + arm; {B} is an inward
# arm + chain. The three-carbon propylene arms of head 306 and its N-methyl are
# written literally, as drawn.
TEMPLATE = "{A}N({B})CCCN(C)CCCN({B}){B}"


def build():
    a = ISODECYL_OUT + ARM_OUT
    b = ARM_IN + ISODECYL_IN
    return TEMPLATE.replace("{A}", a).replace("{B}", b)


def main():
    raw = build()
    m = Chem.MolFromSmiles(raw)
    if m is None:
        raise SystemExit(f"make_reference: unparseable template: {raw}")

    # --- independent assertions on the anchor, from the drawing ----------
    n_ester = len(m.GetSubstructMatches(
        Chem.MolFromSmarts("[CX3](=O)[OX2][CX4]")))
    n_n = sum(1 for at in m.GetAtoms() if at.GetSymbol() == "N")
    n_nh = sum(at.GetTotalNumHs() for at in m.GetAtoms()
               if at.GetSymbol() == "N")
    n_iso = len(m.GetSubstructMatches(Chem.MolFromSmarts("[CH3][CH]([CH3])C")))
    rad = sum(at.GetNumRadicalElectrons() for at in m.GetAtoms())
    if n_ester != 4:
        raise SystemExit(f"anchor: {n_ester} esters, Figure 1A draws 4")
    if n_n != 3:
        raise SystemExit(f"anchor: {n_n} nitrogens, head 306 has 3")
    if n_nh != 0:
        raise SystemExit(f"anchor: {n_nh} N-H left; Figure 1A draws none")
    if n_iso != 4:
        raise SystemExit(f"anchor: {n_iso} terminal isopropyls, need 4")
    if rad:
        raise SystemExit(f"anchor: {rad} radical electrons")

    can = Chem.MolToSmiles(m)
    key = Chem.MolToInchiKey(m)
    formula = rdMolDescriptors.CalcMolFormula(m)
    mass = Descriptors.ExactMolWt(m)

    note = (
        "Derived from Figure 1A, which draws the complete reaction 306 + 4 "
        "Oi10 -> 306Oi10 with every bond. Built by STRING SUBSTITUTION of four "
        "short fragment literals into a tail-first template (make_reference.py)"
        ": no RWMol, no AddBond, no RemoveAtom, no CombineMols and no RDKit "
        "reaction anywhere, so it shares no assembly code with "
        "stage03_enumerate.py, which builds every bond through "
        "curation/lib assemble.attach from an atom-mapped head. The traversal "
        "order is also different on purpose -- the template starts at an "
        "isopropyl terminus and walks inward, which RDKit's canonical order "
        "does not. WHAT IT TESTS THAT NOTHING ELSE CAN: the SUBSTITUTION COUNT "
        "(four tails, two on each primary nitrogen, none on the tertiary one). "
        "Formula and exact mass are consequences of that count, not evidence "
        "for it. LIMIT: Oi10 is 'isodecyl acrylate', a commercial branched-"
        "isomer mixture; the 8-methylnonyl branch is the one the authors DREW, "
        "so this anchor inherits that claim rather than testing it, and every "
        "Oi10 row in 03_enumerated.csv carries structure_assumed = 1.")

    path = os.path.join(HERE, "reference_structures.csv")
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["lipid_name", "reference_SMILES", "source_figure",
                    "read_on", "read_by", "expected_mismatch", "note"])
        w.writerow(["306Oi10", can, "Fig 1A", "2026-11-10",
                    "main session; read at 300 dpi, 3x", "", note])

    print("template      :", raw)
    print("canonical     :", can)
    print("InChIKey      :", key)
    print("formula       :", formula, " exact mass", round(mass, 4))
    print("sha256(canon) :", hashlib.sha256(can.encode()).hexdigest()[:16])
    print("esters 4, N 3, N-H 0, terminal isopropyls 4  -- all asserted")

    enum = os.path.join(HERE, "03_enumerated.csv")
    if os.path.exists(enum):
        for r in csv.DictReader(open(enum)):
            if r["IL_name"] == "306Oi10":
                e = Chem.MolFromSmiles(r["IL_SMILES"])
                print("enumeration   :", r["IL_SMILES"])
                print("MATCH (InChIKey):", Chem.MolToInchiKey(e) == key)
                print("template is a substring of the enumeration's "
                      "canonical SMILES:",
                      TEMPLATE.replace("{A}", ISODECYL_OUT + ARM_OUT)
                      .replace("{B}", ARM_IN + ISODECYL_IN) in r["IL_SMILES"])
                break


if __name__ == "__main__":
    main()
