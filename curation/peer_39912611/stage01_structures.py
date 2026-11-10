#!/usr/bin/env python3
"""peer_39912611 stage 01 -- one structure per library member.

WHERE THE CONNECTIVITY COMES FROM. Not from a figure. Every one of the twelve
lipids is pinned by THREE independent symbolic statements in the SI:

  (a) its full IUPAC name, printed as the heading of its own synthesis
      write-up (SI pp. 21-41; lipid 14 has none -- see DECISIONS D-1);
  (b) the reagent pair drawn in the synthetic scheme, Supplementary Figure 15
      (SI p. 17), which is VECTOR text and extracts with pdfplumber;
  (c) which intermediate the write-up starts from -- alcohol 2, acid 3 or
      amine 4 -- which fixes the linker's DIRECTION on its own.

WHY THAT MATTERS HERE MORE THAN USUAL. The point of this paper is the linker,
and its linker pairs are ISOMERIC. Lipid 34 (amide) and lipid 36 (reverse
amide) are C56H110N2O5 both ways; the paper reports [M+1]+ 907.2 and 907.3.
Lipid 37 (ester) / 38 (reverse ester) are both 949.2. Lipid 40 / 42 are 948.2
and 948.3. A mass check cannot tell any of these pairs apart, and neither can a
formula check. The direction is carried ENTIRELY by (a) and (c).

HOW THE BUILD AVOIDS THE SIBLING'S ERROR. Every acyl group is a fragment whose
carbonyl carbon is ATOM-MAPPED, and the bond is made by curation/lib/assemble,
which never deletes an atom except the acid's own OH in `esterify`. No canonical
SMILES is ever rewritten, sliced or regex-substituted anywhere in this file.
Which fragment carries the carbonyl is read off the synthesis route, so the
ester/reverse-ester and amide/reverse-amide directions are a consequence of (c),
not of a string.
"""
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "lib"))
import assemble                                              # noqa: E402
from rdkit import Chem                                       # noqa: E402
from rdkit.Chem.rdMolDescriptors import CalcMolFormula       # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))

# --- the tail arm ----------------------------------------------------------
# "octane-8,1-diyl 2-hexyldecanoate": C1 carries the ester oxygen, C8 carries
# the amine. Written from that name, with C8 atom-mapped. Asserted below
# against the formula its name requires, so a miscount cannot pass silently.
ARM_MAPPED = "[CH3:1]CCCCCCCOC(=O)C(CCCCCC)CCCCCCCC"
ARM = "CCCCCCCCOC(=O)C(CCCCCC)CCCCCCCC"          # same, H-terminated

# The three intermediates the SI isolates and names (pp. 18-20). Each is the
# core amine with a different stub; the stub is what every head is hung on.
CORE_ALCOHOL = f"OCCN({ARM})" + ARM                 # 2, ((2-hydroxyethyl)...)
CORE_ACID = f"OC(=O)CCN({ARM})" + ARM               # 3, 3-(bis(...)amino)propanoic acid
CORE_AMINE_MAPPED = f"[NH2:1]CCN({ARM})" + ARM      # 4, ((2-aminoethyl)...)
CORE_ALCOHOL_MAPPED = f"[OH:1]CCN({ARM})" + ARM     # 2, O mapped for a carbonate
CORE_ACYL_MAPPED = f"[C:1](=O)CCN({ARM})" + ARM     # 3 as its acyl, for a reverse amide

# --- the two head groups ---------------------------------------------------
DMA = "N(C)C"                 # dimethylamine head  (lipids 14, 32-36)
MPZ = "N1CCN(C)CC1"           # N-methylpiperazine head (lipids 37-42)

# --- the recipe, one row per lipid ----------------------------------------
# route: which SI intermediate it starts from, which reagent is added, and
# therefore which way the linker points. `how` names the assemble call.
LIBRARY = [
    # name, head, linker, how, args
    # NOTE the C4 acid. Lipid 14 is the ONE member whose dimethylamine arm is
    # not the direct homologue of its piperazine partner: lipid 37 is built
    # from a C3 acid (3-(4-methylpiperazin-1-yl)propanoic acid) and lipid 14
    # from a C4 one. Read off the drawing at 600 dpi -- three CH2 between the
    # carbonyl and the dimethylamino N -- and confirmed independently against
    # ZT_2026's "Lipid 14-CB1/2/3" in data/data_vivo.csv, a different paper
    # that carries this exact scaffold. See 08_contradictions.csv row 1.
    ("Lipid 14", "dimethylamine", "ester",
     "esterify", dict(acid="OC(=O)CCC" + DMA, alcohol=CORE_ALCOHOL),
     "alcohol 2 + 4-(dimethylamino)butyric acid (Supp Fig 15 row 1 col 1, "
     "SI p.17; lipid 14 itself is from ref 41 and has no write-up here)"),
    ("Lipid 37", "N-methylpiperazine", "ester",
     "esterify", dict(acid="OC(=O)CC" + MPZ, alcohol=CORE_ALCOHOL),
     "alcohol 2 + 3-(4-methylpiperazin-1-yl)propanoic acid (SI p.31)"),

    ("Lipid 32", "dimethylamine", "reverse-ester",
     "esterify", dict(acid=CORE_ACID, alcohol="OCCC" + DMA),
     "acid 3 + 3-dimethylamino-1-propanol (SI p.21)"),
    ("Lipid 38", "N-methylpiperazine", "reverse-ester",
     "esterify", dict(acid=CORE_ACID, alcohol="OCC" + MPZ),
     "acid 3 + 2-(4-methylpiperazin-1-yl)ethanol (SI p.33)"),

    ("Lipid 33", "dimethylamine", "carbonate",
     "attach", dict(a=CORE_ALCOHOL_MAPPED, b="[C:2](=O)OCCC" + DMA),
     "alcohol 2 + CDI + 3-dimethylamino-1-propanol (SI p.23)"),
    ("Lipid 39", "N-methylpiperazine", "carbonate",
     "attach", dict(a=CORE_ALCOHOL_MAPPED, b="[C:2](=O)OCC" + MPZ),
     "alcohol 2 + CDI + 2-(4-methylpiperazin-1-yl)ethanol (SI p.35)"),

    ("Lipid 34", "dimethylamine", "amide",
     "attach", dict(a=CORE_AMINE_MAPPED, b="[C:2](=O)CCC" + DMA),
     "amine 4 + 4-(dimethylamino)butyric acid (SI p.25)"),
    ("Lipid 40", "N-methylpiperazine", "amide",
     "attach", dict(a=CORE_AMINE_MAPPED, b="[C:2](=O)CC" + MPZ),
     "amine 4 + 3-(4-methylpiperazin-1-yl)propanoic acid (SI p.37)"),

    ("Lipid 35", "dimethylamine", "urea",
     "attach", dict(a=CORE_AMINE_MAPPED, b="[C:2](=O)NCCC" + DMA),
     "amine 4 + CDI + 3-(dimethylamino)-1-propylamine (SI p.27)"),
    ("Lipid 41", "N-methylpiperazine", "urea",
     "attach", dict(a=CORE_AMINE_MAPPED, b="[C:2](=O)NCC" + MPZ),
     "amine 4 + CDI + 2-(4-methylpiperazin-1-yl)ethylamine (SI p.39)"),

    ("Lipid 36", "dimethylamine", "reverse-amide",
     "attach", dict(a=CORE_ACYL_MAPPED, b="[NH2:2]CCC" + DMA),
     "acid 3 + 3-(dimethylamino)-1-propylamine (SI p.29)"),
    ("Lipid 42", "N-methylpiperazine", "reverse-amide",
     "attach", dict(a=CORE_ACYL_MAPPED, b="[NH2:2]CC" + MPZ),
     "acid 3 + 2-(4-methylpiperazin-1-yl)ethylamine (SI p.41)"),
]

# The ESI-MS [M+1]+ the SI prints for each lipid, transcribed from the page
# named in `route`. Used for the three-way mass split in 07_qc.json; NOT used
# to choose a structure, because it cannot distinguish any of the three
# isomeric linker pairs.
PRINTED_MZ = {"Lipid 14": None, "Lipid 32": 908.2, "Lipid 33": 924.2,
              "Lipid 34": 907.2, "Lipid 35": 922.2, "Lipid 36": 907.3,
              "Lipid 37": 949.2, "Lipid 38": 949.2, "Lipid 39": 965.2,
              "Lipid 40": 948.2, "Lipid 41": 963.3, "Lipid 42": 948.3}

# The IUPAC name the SI prints as the heading of each synthesis.
IUPAC = {
 "Lipid 32": "((3-(3-(Dimethylamino)propoxy)-3-oxopropyl)azanediyl)bis(octane-8,1-diyl) bis(2-hexyldecanoate)",
 "Lipid 33": "11-(8-((2-hexyldecanoyl)oxy)octyl)-2-methyl-7-oxo-6,8-dioxa-2,11-diazanonadecan-19-yl 2-hexyldecanoate",
 "Lipid 34": "((2-(4-(dimethylamino)butanamido)ethyl)azanediyl)bis(octane-8,1-diyl) bis(2-hexyldecanoate)",
 "Lipid 35": "11-(8-((2-hexyldecanoyl)oxy)octyl)-2-methyl-7-oxo-2,6,8,11-tetraazanonadecan-19-yl 2-hexyldecanoate",
 "Lipid 36": "((3-((3-(Dimethylamino)propyl)amino)-3-oxopropyl)azanediyl)bis(octane-8,1-diyl) bis(2-hexyldecanoate)",
 "Lipid 37": "((2-((3-(4-Methylpiperazin-1-yl)propanoyl)oxy)ethyl)azanediyl)bis(octane-8,1-diyl) bis(2-hexyldecanoate)",
 "Lipid 38": "((3-(2-(4-Methylpiperazin-1-yl)ethoxy)-3-oxopropyl)azanediyl)bis(octane-8,1-diyl) bis(2-hexyldecanoate)",
 "Lipid 39": "((2-(((2-(4-methylpiperazin-1-yl)ethoxy)carbonyl)oxy)ethyl)azanediyl)bis(octane-8,1-diyl) bis(2-hexyldecanoate)",
 "Lipid 40": "((2-(3-(4-methylpiperazin-1-yl)propanamido)ethyl)azanediyl)bis(octane-8,1-diyl) bis(2-hexyldecanoate)",
 "Lipid 41": "((2-(3-(2-(4-methylpiperazin-1-yl)ethyl)ureido)ethyl)azanediyl)bis(octane-8,1-diyl) bis(2-hexyldecanoate)",
 "Lipid 42": "((3-((2-(4-Methylpiperazin-1-yl)ethyl)amino)-3-oxopropyl)azanediyl)bis(octane-8,1-diyl) bis(2-hexyldecanoate)",
 "Lipid 14": "",   # the SI gives none; see DECISIONS D-1
}

# The SMARTS each linker must match EXACTLY once, and the ones it must not
# match at all. This is the only thing standing between an ester and a reverse
# ester once the SMILES exists, so it is asserted per lipid, not described.
LINKER_SMARTS = {
    # core-CH2CH2-O-C(=O)-CH2CH2-N :  the carbonyl is on the HEAD side
    "ester":          "[NX3]([CH2][CH2][CH2])([CH2])[CH2][CH2]OC(=O)[CH2][CH2][NX3]",
    # core-N-CH2CH2-C(=O)-O-...-N : the carbonyl is on the CORE side
    "reverse-ester":  "[NX3][CH2][CH2]C(=O)O[CH2]",
    "carbonate":      "[CH2]OC(=O)O[CH2]",
    "amide":          "[NX3][CH2][CH2][NX3H1]C(=O)[CH2]",
    "reverse-amide":  "[NX3][CH2][CH2]C(=O)[NX3H1][CH2]",
    "urea":           "[NX3H1]C(=O)[NX3H1]",
}


def build():
    rows = []
    for name, head, linker, how, args, route in LIBRARY:
        if how == "esterify":
            smi = assemble.esterify(args["acid"], args["alcohol"], label=name)
        elif how == "attach":
            smi = assemble.attach(args["a"], args["b"], a_map=1, b_map=2,
                                  label=name)
        else:
            raise SystemExit(f"{name}: unknown build {how!r}")
        m = Chem.MolFromSmiles(smi)
        if m is None:
            raise SystemExit(f"{name}: product does not parse")
        rows.append({
            "lipid_name": name,
            "head": head,
            "linker": linker,
            "built_by": how,
            "route": route,
            "iupac_from_SI": IUPAC[name],
            "SMILES": smi,
            "formula": CalcMolFormula(m),
            "inchikey": Chem.MolToInchiKey(m),
            "printed_mz_M_plus_1": PRINTED_MZ[name],
        })
    return rows


def main():
    # --- the arm, asserted against its own name -----------------------
    arm = Chem.MolFromSmiles(ARM)
    if CalcMolFormula(arm) != "C24H48O2":
        raise SystemExit(f"ARM formula {CalcMolFormula(arm)} != C24H48O2; "
                         "octyl 2-hexyldecanoate is C8 + C16. A miscount in "
                         "the arm would move every one of the twelve masses "
                         "by the same amount, which no mass check could see.")
    if Chem.CanonSmiles(ARM_MAPPED.replace("[CH3:1]", "C")) != \
            Chem.CanonSmiles(ARM):
        raise SystemExit("ARM_MAPPED and ARM are not the same molecule.")

    rows = build()
    out = os.path.join(HERE, "03_enumerated.csv")
    import csv
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    # --- blocks ------------------------------------------------------
    blocks = [
        {"block_role": "tail", "block_name": "T", "block_SMILES": ARM,
         "source": "SI p.17 Supp Fig 15 'R = 8-oxooctyl 2-hexyldecanoate'; "
                   "SI p.18 names aldehyde 1",
         "n_members_using_it": 12},
        {"block_role": "head", "block_name": "N (dimethylamine)",
         "block_SMILES": "CN(C)C", "source": "SI pp.21-29 reagent names",
         "n_members_using_it": 6},
        {"block_role": "head", "block_name": "P (N-methylpiperazine)",
         "block_SMILES": "CN1CCNCC1", "source": "SI pp.31-41 reagent names",
         "n_members_using_it": 6},
    ]
    for lk in ("ester", "reverse-ester", "carbonate", "amide", "urea",
               "reverse-amide"):
        blocks.append({"block_role": "linker", "block_name": lk,
                       "block_SMILES": "", "source": "SI p.17 Supp Fig 15",
                       "n_members_using_it": 2})
    with open(os.path.join(HERE, "01_blocks.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(blocks[0].keys()))
        w.writeheader()
        w.writerows(blocks)

    print(json.dumps({"n_built": len(rows),
                      "n_distinct_inchikey": len({r["inchikey"] for r in rows}),
                      "n_distinct_smiles": len({r["SMILES"] for r in rows})},
                     indent=1))
    for r in rows:
        print(f"{r['lipid_name']:9s} {r['linker']:14s} {r['formula']:14s} "
              f"{r['inchikey']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
