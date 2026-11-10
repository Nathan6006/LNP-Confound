#!/usr/bin/env python3
"""The building blocks of the 40-member RAIL library, read from the SI schemes.

Run from this folder:   python3 -I stage02_blocks.py      -> 01_blocks.csv

WHERE EACH BLOCK COMES FROM. Every structure below was read off a scheme the
authors DRAW in full, at native resolution, after compositing the PNG onto
white (the media PNGs carry an alpha channel; converting straight to RGB turns
the black bond lines into black-on-black and shows an empty panel):

  Scheme S1  word/media/image1.png  the five retinoid acids, drawn in full,
                                    and the VI1 linker Br(CH2)5CO2(CH2)4OH
  Scheme S2  word/media/image2.png  the VI2 linker Br(CH2)6OH
  Scheme S3  word/media/image3.png  T1-T4
  Scheme S4  word/media/image4.png  AT1-AT4
  Scheme S5  word/media/image5.png  the AT1-c-RA coupling, drawn as a product

THE DECOMPOSITION IS NOT THE SYNTHESIS, AND SAYS SO. The real route alkylates a
secondary amine with an alkyl BROMIDE. `assemble` has no function that expels a
leaving group across a C-N bond, and this folder is forbidden to add one, so the
two halves that carry a bromide in the paper are declared here with the bromide
already replaced by the bond that will be formed:

  * a T arm is declared as its alkyl fragment with an atom map on the carbon
    that bore the bromine (`[CH2:2]...`);
  * a VI linker is declared as the corresponding AMINE, `[NH2:1]...`, not as a
    mapped carbon radical -- a nitrogen with two explicit hydrogens and one bond
    is valence-complete, so `esterify` can put the retinoid on it without
    tripping its own `check_no_radicals`. A mapped `[CH2:2]` would be a radical
    and esterify would refuse, correctly.

The product connectivity is identical either way; only the bookkeeping differs.
Nothing here is built by string concatenation and nothing calls RWMol.

NAME TOKENS ARE NOT CHAIN LENGTHS, and here they are not even indices. `AT1`,
`AT2`, `AT3`, `AT4` are four DIFFERENT linker/tail motifs, not a chain walking
by one carbon: AT1 is a C6 acid on pentadecan-8-ol, AT2 a C8 acid on
heptadecan-9-ol, AT3 a C10 acid on heptadecan-9-ol, and AT4 is the ester
REVERSED -- a 6-bromohexyl ester of 2-hexyldecanoic acid, i.e. the carbonyl
points the other way. A SMILES fragment has a direction that mass cannot see,
and AT4 is where that would bite: see stage03's per-lipid InChIKey table.
"""
import csv
import json
import os
import sys

LIB = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib")
sys.path.insert(0, LIB)
import assemble                                                  # noqa: E402

from rdkit import Chem                                           # noqa: E402
from rdkit.Chem import Descriptors, rdMolDescriptors             # noqa: E402
from rdkit import RDLogger                                       # noqa: E402
RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "01_blocks.csv")

# ---------------------------------------------------------------------------
# 1. THE FIVE RETINOID ACIDS -- the varying hydrophobic domain.
#
# Read from Scheme S1, where all five are drawn in full beside the product.
# The three retinoic acids are CONSTITUTIONAL TWINS: identical connectivity,
# differing only in the geometry of one polyene double bond. Mass is blind to
# that, which is exactly why stage03 asserts the E/Z pattern explicitly and
# pins each one by FULL InChIKey rather than by formula.
# ---------------------------------------------------------------------------
# Where the Z bond must sit, as the topological distance in bonds from the
# CARBOXYL carbon to the nearer of the two atoms of the Z double bond. Counting
# E and Z is not enough: 13-cis and 9-cis both have 3E/1Z, so a count alone
# cannot tell them apart and would pass with the two swapped.
#   13CRA: Z is C13=C14, and C14 is ONE bond from the acid carbon C15.
#    9CRA: Z is  C9=C10, and C10 is FIVE bonds from C15.
Z_BOND_DISTANCE_FROM_ACID = {"13CRA": 1, "9CRA": 5}

RETINOIDS = {
    # all-trans: 7E,9E,11E,13E
    "ATRA":  ("OC(=O)/C=C(C)/C=C/C=C(C)/C=C/C1=C(C)CCCC1(C)C",
              "all-trans-retinoic acid", "C20H28O2", (4, 0)),
    # 13-cis (isotretinoin): the Z bond is the one NEXT TO the carboxyl
    # Every directional marker DOWNSTREAM of a flipped bond has to flip too,
    # or the flip propagates into the next double bond. Writing only the first
    # two as "\\" gave 2E/2Z and the geometry assertion below caught it.
    "13CRA": ("OC(=O)/C=C(C)\\C=C\\C=C(C)\\C=C\\C1=C(C)CCCC1(C)C",
              "13-cis-retinoic acid (isotretinoin)", "C20H28O2", (3, 1)),
    # 9-cis (alitretinoin): the Z bond is the 9,10-ene, further from the acid
    "9CRA":  ("OC(=O)/C=C(C)/C=C/C=C(C)\\C=C\\C1=C(C)CCCC1(C)C",
              "9-cis-retinoic acid (alitretinoin)", "C20H28O2", (3, 1)),
    # acitretin: an aromatic ring replaces the trimethylcyclohexene
    "AVA":   ("OC(=O)/C=C(C)/C=C/C=C(C)/C=C/c1c(C)c(C)c(OC)cc1C",
              "acitretin", "C21H26O3", (4, 0)),
    # bexarotene: no polyene at all; an exocyclic 1,1-disubstituted alkene,
    # which has no E/Z to assign
    "Bex":   ("OC(=O)c1ccc(cc1)C(=C)c1cc2c(cc1C)C(C)(C)CCC2(C)C",
              "bexarotene", "C24H28O2", (0, 0)),
}

# ---------------------------------------------------------------------------
# 2. THE FOUR ATX TAIL/LINKER MOTIFS -- declared as alkyl fragments, with the
#    atom map on the carbon that carries the bromine in Scheme S3.
# ---------------------------------------------------------------------------
T_ALKYL = {
    "AT1": ("[CH2:2]CCCCC(=O)OC(CCCCCCC)CCCCCCC",
            "from T1 = 8-pentadecyl 6-bromohexanoate (Scheme S3 row 1); "
            "C6 acid, pentadecan-8-ol, ester points acid->alcohol"),
    "AT2": ("[CH2:2]CCCCCCC(=O)OC(CCCCCCCC)CCCCCCCC",
            "from T2 = 9-heptadecyl 8-bromooctanoate (Scheme S3 row 2); "
            "C8 acid, heptadecan-9-ol"),
    "AT3": ("[CH2:2]CCCCCCCCC(=O)OC(CCCCCCCC)CCCCCCCC",
            "from T3 = 9-heptadecyl 10-bromodecanoate (Scheme S3 row 3); "
            "C10 acid, heptadecan-9-ol"),
    "AT4": ("[CH2:2]CCCCCOC(=O)C(CCCCCC)CCCCCCCC",
            "from T4 = 6-bromohexyl 2-hexyldecanoate (Scheme S3 row 4); "
            "ESTER REVERSED relative to AT1-AT3 -- the bromoalkyl is the "
            "ALCOHOL half here and 2-hexyldecanoic acid is the acyl half"),
}

# ---------------------------------------------------------------------------
# 3. THE TWO LINKER ARCHITECTURES -- the 'c' in AT1-c-ATRA.
#    Declared as the amine, standing in for the bromide (see the docstring).
# ---------------------------------------------------------------------------
LINKER = {
    "c":    ("[NH2:1]CCCCCC(=O)OCCCCO",
             "VI1 (Scheme S1): Br(CH2)5CO2(CH2)4OH. Two ester bonds between "
             "the amine and the retinoid -- the longer, doubly-ester linker. "
             "The library name writes this as '-c-'."),
    "":     ("[NH2:1]CCCCCCO",
             "VI2 (Scheme S2): Br(CH2)6OH. One ester bond to the retinoid. "
             "The library name writes this as nothing at all, e.g. AT1-ATRA."),
}

# ---------------------------------------------------------------------------
# 4. THE HEAD ARM -- 4-amino-1-butanol's hydroxybutyl, mapped at the carbon
#    that bonds to nitrogen.
#
# THE PAPER'S OWN SCHEME DISAGREES WITH ITS OWN TEXT HERE. The SI text and the
# Materials list both say 4-amino-1-butanol; Scheme S4 DRAWS HO-CH2CH2CH2-NH2,
# three carbons, i.e. 3-amino-1-propanol. The 39 printed HRMS values settle it
# for the four-carbon reading and against the three-carbon one by a full CH2 --
# see stage03 and 08_contradictions.csv. Four carbons is what is built.
# ---------------------------------------------------------------------------
HEAD_ARM = ("[CH2:2]CCCO",
            "4-hydroxybutyl, from 4-amino-1-butanol (SI Materials and the "
            "Synthesis of AT1 paragraph). Scheme S4 draws three carbons; the "
            "HRMS says four. See 08_contradictions.csv row 1.")


def _mol(smiles, what):
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        raise SystemExit(f"stage02: {what} does not parse: {smiles}")
    return m


def _stereo_counts(mol):
    """(#E, #Z) over assigned double bonds."""
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    e = z = 0
    for b in mol.GetBonds():
        s = b.GetStereo()
        if s == Chem.BondStereo.STEREOE:
            e += 1
        elif s == Chem.BondStereo.STEREOZ:
            z += 1
    return e, z


def _z_bond_distance_from_acid(mol):
    """Bonds from the carboxyl carbon to the nearer atom of the Z double bond.

    Returns None where there is no Z bond. Raises if there is more than one,
    because every molecule here has at most one by construction.
    """
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    zs = [b for b in mol.GetBonds() if b.GetStereo() == Chem.BondStereo.STEREOZ]
    if not zs:
        return None
    if len(zs) > 1:
        raise SystemExit("stage02: more than one Z bond; cannot locate it.")
    acid = mol.GetSubstructMatches(Chem.MolFromSmarts("[CX3](=O)[OX2H1]"))
    if len(acid) != 1:
        raise SystemExit(f"stage02: {len(acid)} carboxyl groups, need 1.")
    c_acid = acid[0][0]
    dm = Chem.GetDistanceMatrix(mol)
    b = zs[0]
    return int(min(dm[c_acid][b.GetBeginAtomIdx()],
                   dm[c_acid][b.GetEndAtomIdx()]))


def main():
    rows = []
    keys = {}

    # --- the retinoids, each checked against its drawn formula and geometry --
    for name, (smi, longname, formula, want_ez) in RETINOIDS.items():
        m = _mol(smi, name)
        got_formula = rdMolDescriptors.CalcMolFormula(m)
        if got_formula != formula:
            raise SystemExit(
                f"stage02: {name} formula {got_formula}, expected {formula}. "
                "The expected formula is not a guess -- it is what the SI's "
                "own ESI-MS for VI1-{name} requires.")
        got_ez = _stereo_counts(m)
        if got_ez != want_ez:
            raise SystemExit(
                f"stage02: {name} has {got_ez[0]}E/{got_ez[1]}Z double bonds, "
                f"expected {want_ez[0]}E/{want_ez[1]}Z. CIP stereocentre "
                "counting is blind to double-bond geometry, so this is "
                "asserted separately and on purpose: ATRA, 13CRA and 9CRA are "
                "constitutional twins and nothing else in this pipeline can "
                "tell them apart.")
        want_d = Z_BOND_DISTANCE_FROM_ACID.get(name)
        if want_d is not None:
            got_d = _z_bond_distance_from_acid(m)
            if got_d != want_d:
                raise SystemExit(
                    f"stage02: {name}'s Z double bond is {got_d} bonds from "
                    f"the carboxyl carbon, expected {want_d}. 13-cis and "
                    "9-cis retinoic acid BOTH have 3E/1Z, so the count alone "
                    "passes with the two isomers swapped; this is what "
                    "distinguishes them.")
        elif _z_bond_distance_from_acid(m) is not None:
            raise SystemExit(f"stage02: {name} has a Z bond and should not.")

        key = Chem.MolToInchiKey(m)
        if key in keys:
            raise SystemExit(
                f"stage02: {name} and {keys[key]} have the SAME InChIKey "
                f"{key}. Two of the three retinoic acid isomers collapsed "
                "into one structure -- the E/Z markers are wrong.")
        keys[key] = name
        rows.append(["retinoid", name, longname, smi, got_formula,
                     f"{Descriptors.ExactMolWt(m):.6f}", key,
                     f"{got_ez[0]}E/{got_ez[1]}Z",
                     "Scheme S1, drawn in full beside the product"])

    # --- the three retinoic acids must share one skeleton -------------------
    ra = [k for k, v in keys.items() if v in ("ATRA", "13CRA", "9CRA")]
    blocks = {k.split("-")[0] for k in ra}
    if len(blocks) != 1:
        raise SystemExit(
            f"stage02: the three retinoic acid isomers have {len(blocks)} "
            "different InChIKey skeleton blocks. They are geometric isomers "
            "of ONE connectivity; more than one block means a connectivity "
            "error, not a geometry one.")

    # --- the tail motifs ----------------------------------------------------
    for name, (smi, note) in T_ALKYL.items():
        m = _mol(smi, name)
        rows.append(["atx_tail_motif", name, "", smi,
                     rdMolDescriptors.CalcMolFormula(m), "", "", "", note])

    for name, (smi, note) in LINKER.items():
        m = _mol(smi, name or "(no linker token)")
        rows.append(["linker", name or "(none)", "", smi,
                     rdMolDescriptors.CalcMolFormula(m), "", "", "", note])

    m = _mol(HEAD_ARM[0], "head arm")
    rows.append(["head_arm", "4-hydroxybutyl", "", HEAD_ARM[0],
                 rdMolDescriptors.CalcMolFormula(m), "", "", "", HEAD_ARM[1]])

    with open(OUT, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["block_class", "block_name", "long_name", "SMILES",
                    "formula", "exact_mass", "inchikey", "double_bond_geometry",
                    "source"])
        w.writerows(rows)

    # the factorial the paper claims, proved rather than restated
    n = len(T_ALKYL) * len(LINKER) * len(RETINOIDS)
    print(f"wrote {len(rows)} blocks to {OUT}")
    print(f"factorial: {len(T_ALKYL)} ATX motifs x {len(LINKER)} linkers x "
          f"{len(RETINOIDS)} retinoids = {n}")
    if n != 40:
        raise SystemExit(
            f"stage02: the blocks multiply out to {n}, but the abstract, the "
            "Results section and Scheme 1 all say 40. EXTRACTION_BRIEF "
            "section 2: if a factorial is claimed, prove it.")
    print("            = 40, which is the library size the abstract, "
          "Results 2.1 and Scheme 1(a) all state.")
    print("NOTE: the paper calls these 'two molecular scaffolds (ATX-c-RA and "
          "ATX-RA), four linker motifs, and five ... hydrophobic domains'.")
    print("      Its 'scaffold' is our `linker` (2) and its 'linker motif' is "
          "our `atx_tail_motif` (4). Same 2 x 4 x 5; the words are swapped.")


if __name__ == "__main__":
    main()
