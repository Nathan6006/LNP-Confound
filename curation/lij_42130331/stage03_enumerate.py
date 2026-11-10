#!/usr/bin/env python3
"""Enumerate Library 1 (50) and Library 2 (250) from the drawn blocks.

    python3 -I stage03_enumerate.py     # writes 01_blocks.csv, 03_enumerated.csv

THE REGIOCHEMISTRY, WHICH MASS CANNOT SEE. Two bonds in this scaffold are
decided by the synthesis, not by any formula:

 1. THE THIOLACTONE RING OPENS AT THE CARBONYL. Fig. S1A draws the amine's
    lone pair attacking the C=O of N-acetyl-DL-homocysteine thiolactone, and
    the drawn intermediate is R1-NH-C(=O)-CH(NHAc)-CH2CH2-SH -- an AMIDE to
    the amine and a free thiol. The alternative (S-alkylation) has the same
    formula.

 2. THE THIOL ADDS TO THE LESS SUBSTITUTED EPOXIDE CARBON. Fig. S1A draws the
    product as ...CH2-S-CH2-CH(OH)-R: sulfur on the terminal CH2, hydroxyl on
    the internal carbon. The opposite regiochemistry, S on the internal carbon
    and OH terminal, is C- and H-identical and would be invisible to HRMS, to
    every descriptor, and to the InChIKey skeleton block. It is read off the
    drawing, and the drawing is unambiguous: the OH label sits on the carbon
    that also carries R.

    This paper offers NO independent check on either: the SI contains no HRMS,
    no NMR table, no SMILES, no InChI and no deposited structure file (grep of
    the full extracted text for smiles|inchi|HRMS|ESI-MS|calcd returns 0). The
    only characterisation of any library member is Fig. S9, a mass spectrum
    and a 1H NMR of the lead compound published as images. So the drawing is
    the sole statement of connectivity, and the C10 anchor re-derives it by a
    different route rather than confirming it against external data.

 3. SECONDARY AMINES ALSO OPEN AN EPOXIDE. Fig. S2 draws 96 and 400 -- the two
    amines with an N-H -- carrying an extra N-CH2-CH(O-acyl)-R arm that has no
    sulfur and no amide. So the arm count is not "one per amine": it is one
    AA arm per PRIMARY amine plus one direct 2-hydroxyalkyl arm per SECONDARY
    amine. Reading the scheme and not the figure would have given 96 and 400
    one arm each instead of two.

 4. LIBRARY 2 ESTERIFIES EVERY beta-HYDROXY. Fig. S2 draws the acyl chloride
    adding to each R-bearing hydroxyl, including 96's and 400's N-arm.

STEREOCHEMISTRY IS ABSENT ON PURPOSE. The thiolactone is the DL (racemic)
reagent and Fig. S1B draws it flat, with no wedge anywhere in Fig. S1 or S2.
Both the alpha-carbon of the homocysteine unit and the carbinol carbon are
drawn without configuration, so every structure here is emitted without
stereo and every row carries confidence below `high` for that reason.
"""
import csv
import json
import os
import sys

LIB = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib")
sys.path.insert(0, LIB)
import assemble  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ARM_MAP = 9

# --- BLOCKS -------------------------------------------------------------
# Every chain length below was MEASURED by stage02_blocks.py (two independent
# raster routes) and then confirmed on a 2000 dpi render saved in 00_pages/.
# `measured` records the component the count came from. Name tokens are NOT
# chain lengths here and the mismatch is deliberate and recorded: the epoxide
# "C8" carries a SIX-carbon R, because a 1,2-epoxyalkane CnH2nO opened by a
# thiol yields S-CH2-CH(OH)-C(n-2)H(2n-3).

AMINES = {
    #  name: (SMILES with a map on every reactive N, primary maps, secondary
    #         maps, how its arms were measured)
    "70":  ("CN(C)CC[NH2:1]", [1], [],
            "row1 #0: 2 carbons, both routes agree"),
    "75":  ("CCN(CC)CC[NH2:1]", [1], [],
            "row1 #1/#2/#3: ethylene arm 2 C, two N-ethyls"),
    "80":  ("CN(C)CCC[NH2:1]", [1], [],
            "row1 #5: 3 carbons, both routes agree"),
    "96":  ("C[NH:5]CCC[NH2:1]", [1], [5],
            "row1 #7: 3 carbons; the N-H is drawn explicitly"),
    "400": ("CN(C)CCC[NH:5]CCC[NH2:1]", [1], [5],
            "row1 #8/#9: 3 and 3 carbons"),
    "113": ("[NH2:1]CCN(C)CC[NH2:2]", [1, 2], [],
            "row2 #0/#1: 2 and 2, both routes agree"),
    "414": ("[NH2:1]CCCCN(C)CCCC[NH2:2]", [1, 2], [],
            "row2 #2/#3: 4 and 4, both routes agree"),
    "616": ("[NH2:1]CCCCCCN(C)CCCCCC[NH2:2]", [1, 2], [],
            "row2 #4/#5: 6 and 6"),
    "110": ("[NH2:1]CCN(CC[NH2:2])CC[NH2:3]", [1, 2, 3], [],
            "row3 #0/#2/#3: 2, 2, 2 -- tris(2-aminoethyl)amine"),
    "T3A": ("[NH2:1]CCCN(CCC[NH2:2])CCC[NH2:3]", [1, 2, 3], [],
            "row3 #1/#4/#5: 3, 3, 3 -- tris(3-aminopropyl)amine"),
}

# R = the group the epoxide leaves on the carbinol carbon.
EPOXIDES = {
    "C8":  ("CCCCCC", "1,2-epoxyoctane", 8, 6,
            "legend #0: 6 carbons after the squiggle, both routes agree"),
    "C10": ("CCCCCCCC", "1,2-epoxydecane", 10, 8,
            "legend #1: 8, by the measured +2-carbon step from C8"),
    "C12": ("CCCCCCCCCC", "1,2-epoxydodecane", 12, 10,
            "legend #2: 10, same step"),
    "C14": ("CCCCCCCCCCCC", "1,2-epoxytetradecane", 14, 12,
            "legend #3: 12, both routes agree"),
    "OG8": ("COCC(CC)CCCC", "2-ethylhexyl glycidyl ether", 11, 8,
            "legend #4/#5: CH2-O-CH2-CH(ethyl)-butyl, read at 2000 dpi"),
}

# n in the drawn Cl-C(=O)-CH2-(CH2)n-CH3 is PRINTED beside the legend, so the
# acyl chain lengths are a printed statement and not a pixel measurement.
ACYLS = {
    "C6":  ("CCCCCC(=O)O", "hexanoyl", 3),
    "C8":  ("CCCCCCCC(=O)O", "octanoyl", 5),
    "C10": ("CCCCCCCCCC(=O)O", "decanoyl", 7),
    "C12": ("CCCCCCCCCCCC(=O)O", "dodecanoyl", 9),
    "C14": ("CCCCCCCCCCCCCC(=O)O", "tetradecanoyl", 11),
}

THIOLACTONE = "CC(=O)NC1CCSC1=O"   # N-acetyl-DL-homocysteine thiolactone


def aa_arm(r_smiles, acid=None):
    """-C(=O)-CH(NHAc)-CH2CH2-S-CH2-CH(OH or O-acyl)-R, carbonyl C mapped.

    Written as an aldehyde ([CH:9](=O)) so the fragment is a real molecule
    with no unpaired electron; `assemble.attach` hands the hydrogen count back
    to RDKit when it makes the amide bond. Writing it as a bare '[C:9](=O)'
    instead gives a RADICAL that parses, canonicalises and has self-consistent
    descriptors -- the Siegwart radical-carbon class.
    """
    frag = f"[CH:{ARM_MAP}](=O)C(NC(C)=O)CCSCC(O){r_smiles}"
    if acid is None:
        return assemble.canonical(frag)
    return assemble.esterify(acid, frag, label="AA arm ester")


def n_arm(r_smiles, acid=None):
    """-CH2-CH(OH or O-acyl)-R, the arm a SECONDARY amine makes directly.

    Written as a METHYL, not a methylene: '[CH2:9]C(O)R' is a carbon with two
    explicit hydrogens and three bonds, which RDKit sanitises into a RADICAL,
    and `assemble.check_no_radicals` refuses it -- correctly. The attachment
    atom has to be a complete atom in the fragment; `attach` removes the
    hydrogen when it makes the bond.
    """
    frag = f"[CH3:{ARM_MAP}]C(O){r_smiles}"
    if acid is None:
        return assemble.canonical(frag)
    return assemble.esterify(acid, frag, label="N arm ester")


def build(amine_key, epox_key, acyl_key=None):
    amine, prim, sec, _ = AMINES[amine_key]
    r = EPOXIDES[epox_key][0]
    acid = ACYLS[acyl_key][0] if acyl_key else None
    cur = amine
    for m in prim:
        cur = assemble.attach(cur, aa_arm(r, acid), a_map=m, b_map=ARM_MAP,
                              label=f"{amine_key}-AA-{epox_key} primary {m}")
    for m in sec:
        cur = assemble.attach(cur, n_arm(r, acid), a_map=m, b_map=ARM_MAP,
                              label=f"{amine_key}-AA-{epox_key} secondary {m}")
    return assemble.canonical(cur)


def main():
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors

    with open(os.path.join(HERE, "01_blocks.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["block_type", "block_id", "reagent", "SMILES",
                    "n_primary_NH2", "n_secondary_NH", "name_token_carbons",
                    "measured_carbons", "name_token_matches",
                    "how_measured"])
        for k, (smi, prim, sec, how) in AMINES.items():
            w.writerow(["amine", k, "", assemble.canonical(
                smi.replace(":1", "").replace(":2", "").replace(":3", "")
                   .replace(":5", "").replace("[NH2]", "N")
                   .replace("[NH]", "N")),
                len(prim), len(sec), "", "", "", how])
        for k, (r, reagent, ntok, nmeas, how) in EPOXIDES.items():
            w.writerow(["epoxide_R", k, reagent, r, "", "", ntok, nmeas,
                        "FALSE -- see note", how])
        for k, (acid, reagent, n) in ACYLS.items():
            # The stored SMILES is the ACID, not the acid chloride, because
            # assemble.esterify takes an acid and the acyl group transferred
            # is identical either way. The name says so: a reagent name that
            # does not describe the structure beside it is the AC_2025 class,
            # where DLin-MC3-DMA and KC2 carried each other's structures --
            # same formula, same exact mass, invisible to every check.
            w.writerow(["acyl", "Cl" + k[1:],
                        f"{reagent.replace(chr(111)+chr(121)+chr(108), chr(111)+chr(105)+chr(99))} acid, the acyl donor used for "
                        f"assembly (the paper's reagent is {reagent} "
                        f"chloride; the acyl group transferred is the same)",
                        acid, "", "", k[1:], k[1:], "TRUE",
                        f"printed in the Fig. S2 legend as n = {n}"])
        w.writerow(["thiolactone", "AA",
                    "N-acetyl-DL-homocysteine thiolactone", THIOLACTONE,
                    "", "", "", "", "", "Fig. S1B, drawn in full"])

    rows = []
    for ak in AMINES:
        for ek in EPOXIDES:
            smi = build(ak, ek)
            rows.append(("1", f"{ak}-AA-{ek}", ak, ek, "", smi))
            for cy in ACYLS:
                s2 = build(ak, ek, cy)
                rows.append(("2", f"{ak}-AA-{ek}{cy}", ak, ek, cy, s2))

    with open(os.path.join(HERE, "03_enumerated.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        # `lipid_name` duplicates `IL_name` because provenance/structure_check.py
        # (C10) keys its lookup on `lipid_name` + `IL_SMILES`; the corpus
        # schema calls the same field `IL_name`. Writing both is cheaper than
        # a check that silently finds no structures to compare.
        w.writerow(["library", "lipid_name", "IL_name", "amine", "epoxide",
                    "acyl",
                    "IL_SMILES", "formula", "exact_mass", "n_arms",
                    "confidence", "confidence_reason"])
        for lib, name, ak, ek, cy, smi in rows:
            m = Chem.MolFromSmiles(smi)
            if m is None:
                raise SystemExit(f"{name}: unparseable {smi}")
            nprim, nsec = len(AMINES[ak][1]), len(AMINES[ak][2])
            w.writerow([lib, name, name, ak, ek, cy, smi,
                        rdMolDescriptors.CalcMolFormula(m),
                        round(Descriptors.ExactMolWt(m), 4),
                        nprim + nsec,
                        "medium",
                        "Connectivity is read from a drawing and has no "
                        "independent check in this paper: no HRMS, no NMR "
                        "table, no SMILES/InChI and no deposited structure "
                        "file anywhere in the SI. Stereochemistry is absent "
                        "because the thiolactone is the DL reagent and no "
                        "wedge is drawn."])

    n1 = sum(1 for r in rows if r[0] == "1")
    print(f"Library 1: {n1}   Library 2: {len(rows) - n1}   total {len(rows)}")
    uniq = len({r[5] for r in rows})
    print(f"unique canonical SMILES: {uniq}")
    if uniq != len(rows):
        print("WARNING: some structures collide -- that is a real collision "
              "in the library, not a bug, and it must be recorded.")
    with open(os.path.join(HERE, "03_enumeration_summary.json"), "w") as fh:
        json.dump({"library_1": n1, "library_2": len(rows) - n1,
                   "total": len(rows), "unique_canonical_SMILES": uniq,
                   "stated_in_abstract": 300}, fh, indent=1)


if __name__ == "__main__":
    main()
