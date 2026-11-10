#!/usr/bin/env python3
"""Enumerate the 148 Passerini lipids and check every one against Table S1.

THE REACTION, read off the scheme and NOT inferred from mass. The Passerini
3CR is atom-economical: acid + aldehyde + isocyanide -> product, no atoms lost.
Which component contributes which fragment is therefore invisible to the
molecular formula, and the scheme (Fig 1A) plus the drawn library (SI Figs
S1-S9) settle it:

    R_head-COOH  +  R_B-CHO  +  R_A-NC
        -> R_A-NH-C(=O)-CH(R_B)-O-C(=O)-R_head

  * the ACID becomes the ester acyl group (the ionizable head sits here),
  * the ALDEHYDE's carbonyl carbon becomes the alpha-CH, keeping its chain
    and acquiring the ester oxygen,
  * the ISOCYANIDE's N becomes the amide N-H and its terminal C the amide C=O.

Read off SI Fig S1 (H1A1B1 etc.): the blue Me2N-CH2-C(=O)-O- is the head acid
as an ester; the alpha-carbon carries a nonyl branch (decanal minus CHO); the
amide N carries the undecyl chain (the isocyanide's alkyl).

ASSEMBLY ROUTE. Deliberately NOT a string template. Every bond is made by
curation/lib/assemble, which asserts no-radicals and stereo-preservation on
every build. The two fragments used are the paper's own synthetic precursors:

  isocyanide R_A-NC   <- its N-alkylformamide R_A-NH-CHO   (SI p.3: the
                         isocyanides are made by dehydrating this formamide)
  aldehyde   R_B-CHO  <- its primary alcohol R_B-CH2-OH    (SI p.3: the
                         aldehydes are made by DMP oxidation of this alcohol)

so `attach` joins the formamide's formyl carbon to the alcohol's carbinol
carbon (both keep every atom; RDKit drops the two hydrogens by valence), and
`esterify` then condenses the head acid onto the resulting secondary alcohol.
Net: acid + aldehyde + isocyanide, exactly.
"""
import csv, json, os, sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "lib"))
import assemble                                              # noqa: E402
from rdkit import Chem, RDLogger                             # noqa: E402
from rdkit.Chem import rdMolDescriptors, Descriptors         # noqa: E402

RDLogger.DisableLog("rdApp.*")
HERE = os.path.dirname(os.path.abspath(__file__))

# ---------------------------------------------------------------- blocks ---
# Every SMILES below is written FROM THE CHEMICAL NAME printed in SI Section A
# (or, for H3, from the structure drawn in SI Fig S3). None is copied from the
# corpus or from another curation folder.

HEAD = {            # name in SI, SMILES, where the name is printed
 "H1": ("N,N-dimethylglycine",                  "CN(C)CC(=O)O",
        "SI Fig S1 caption + drawn head group"),
 "H2": ("3-(dimethylamino)propanoic acid",      "CN(C)CCC(=O)O",
        "SI Fig S2 caption"),
 "H3": ("N,N-diethylglycine",                   "CCN(CC)CC(=O)O",
        "SI Fig S3 DRAWN head group (its caption misreads 'dimethylglycine', "
        "which is H1's name -- see 08_contradictions.csv C-01)"),
 "H4": ("1-methylpiperidine-4-carboxylic acid", "CN1CCC(CC1)C(=O)O",
        "SI Fig S4 caption"),
 "H5": ("3-(piperidin-1-yl)propanoic acid",     "OC(=O)CCN1CCCCC1",
        "SI Fig S5 caption"),
 "H6": ("3-(pyrrolidin-1-yl)propanoic acid",    "OC(=O)CCN1CCCC1",
        "SI Fig S6 caption"),
 "H7": ("4-(dimethylamino)benzoic acid",        "CN(C)c1ccc(cc1)C(=O)O",
        "SI Fig S7 caption"),
 "H8": ("4-(diethylamino)benzoic acid",         "CCN(CC)c1ccc(cc1)C(=O)O",
        "SI Fig S8 caption"),
 "H9": ("4-(4-methylpiperazin-1-yl)benzoic acid",
        "CN1CCN(CC1)c1ccc(cc1)C(=O)O", "SI Fig S9 caption"),
 "S3": ("4-(dimethylamino)butanoic acid",       "CN(C)CCCC(=O)O",
        "SI p.5, A4B4-S3 synthesis"),
 "S4": ("5-(dimethylamino)pentanoic acid",      "CN(C)CCCCC(=O)O",
        "SI p.5, A4B4-S4 synthesis"),
}

# formamide precursors; map 1 = the formyl carbon that becomes the amide C=O
TAILA = {
 "A1": ("1-isocyanoundecane",            "CCCCCCCCCCCN[CH:1]=O",      "SI p.3"),
 "A2": ("1-isocyanododecane",            "CCCCCCCCCCCCN[CH:1]=O",     "SI p.3"),
 "A3": ("1-isocyanohexadecane",          "CCCCCCCCCCCCCCCCN[CH:1]=O", "SI p.3"),
 "A4": ("(Z)-1-isocyanooctadec-9-ene",
        r"CCCCCCCC/C=C\CCCCCCCCN[CH:1]=O", "SI p.3"),
}

# primary-alcohol precursors; map 2 = the carbinol carbon that becomes alpha-CH
TAILB = {
 "B1": ("decanal",                       "CCCCCCCCC[CH2:2]O",         "SI p.3"),
 "B2": ("olealdehyde ((Z)-octadec-9-enal)",
        r"CCCCCCCC/C=C\CCCCCCC[CH2:2]O", "SI p.3"),
 "B3": ("6-oxohexyl dodecanoate",
        "CCCCCCCCCCCC(=O)OCCCCC[CH2:2]O", "SI pp.3-4"),
 "B4": ("6-oxohexyl 2-hexyldecanoate",
        "CCCCCCCCC(CCCCCC)C(=O)OCCCCC[CH2:2]O", "SI p.4"),
 "Bbu": ("6-oxohexyl 2-butyloctanoate",
         "CCCCCCC(CCCC)C(=O)OCCCCC[CH2:2]O", "SI p.4"),
}

# the four lipids made outside the 9x4x4 grid
EXTRA = [("H2A4-Ole",  "H2", "A4", "B2"),
         ("H2A4-Butyl", "H2", "A4", "Bbu"),
         ("A4B4-S3",   "S3", "A4", "B4"),
         ("A4B4-S4",   "S4", "A4", "B4")]


def formula(smi):
    return rdMolDescriptors.CalcMolFormula(Chem.MolFromSmiles(smi))


def build(h, a, b, label):
    alpha = assemble.attach(TAILA[a][1], TAILB[b][1], a_map=1, b_map=2,
                            label=label + ":alpha")
    return assemble.esterify(HEAD[h][1], alpha, label=label + ":ester")


def main():
    # ---- block sanity: every precursor must give back the SI's own reagent --
    blocks = []
    for k, (nm, smi, where) in HEAD.items():
        blocks.append(dict(block=k, role="head_acid", si_name=nm,
                           smiles=assemble.canonical(smi),
                           formula=formula(smi), source=where))
    for k, (nm, smi, where) in TAILA.items():
        f = Chem.MolFromSmiles(smi.replace(":1", ""))
        blocks.append(dict(block=k, role="tailA_isocyanide",
                           si_name=nm, smiles=smi,
                           formula=rdMolDescriptors.CalcMolFormula(f) +
                           " (as the N-alkylformamide precursor)",
                           source=where))
    for k, (nm, smi, where) in TAILB.items():
        f = Chem.MolFromSmiles(smi.replace(":2", ""))
        blocks.append(dict(block=k, role="tailB_aldehyde",
                           si_name=nm, smiles=smi,
                           formula=rdMolDescriptors.CalcMolFormula(f) +
                           " (as the primary-alcohol precursor)",
                           source=where))
    with open(os.path.join(HERE, "01_blocks.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(blocks[0]))
        w.writeheader(); w.writerows(blocks)

    # ---- enumerate -------------------------------------------------------
    rows = []
    for h in ["H%d" % i for i in range(1, 10)]:
        for a in ["A%d" % i for i in range(1, 5)]:
            for b in ["B%d" % i for i in range(1, 5)]:
                nm = h + a + b
                smi = build(h, a, b, nm)
                rows.append(dict(IL_name=nm, head=h, tailA=a, tailB=b,
                                 in_library=1, IL_SMILES=smi))
    for nm, h, a, b in EXTRA:
        rows.append(dict(IL_name=nm, head=h, tailA=a, tailB=b,
                         in_library=0, IL_SMILES=build(h, a, b, nm)))

    for r in rows:
        m = Chem.MolFromSmiles(r["IL_SMILES"])
        r["Formula"] = rdMolDescriptors.CalcMolFormula(m)
        r["MW_calc"] = round(Descriptors.MolWt(m), 2)
        r["ExactMolWt"] = round(Descriptors.ExactMolWt(m), 4)
        r["InChIKey"] = Chem.MolToInchiKey(m)

    with open(os.path.join(HERE, "03_enumerated.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["IL_name", "head", "tailA", "tailB",
                                           "in_library", "Formula", "MW_calc",
                                           "ExactMolWt", "InChIKey",
                                           "IL_SMILES"])
        w.writeheader(); w.writerows(rows)
    print("03_enumerated.csv: %d structures (%d library + %d extra), "
          "%d distinct canonical SMILES, %d distinct InChIKeys"
          % (len(rows), sum(r["in_library"] for r in rows),
             len(rows) - sum(r["in_library"] for r in rows),
             len({r["IL_SMILES"] for r in rows}),
             len({r["InChIKey"] for r in rows})))
    return 0


if __name__ == "__main__":
    sys.exit(main())
