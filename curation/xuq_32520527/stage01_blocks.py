#!/usr/bin/env python3
"""XQ_2020 stage 01 -- the building blocks, written from measured evidence.

    python3 -I stage01_blocks.py

Two kinds of block:

  HEAD   the eleven commercial amines drawn in full in Figure 1C.
  TAIL   O16CBA, the acrylate-terminated cyclic-benzylidene-acetal tail,
         drawn in Figure 1B and in Figure S1, and built here from the
         Experimental Section's NAMED reagents rather than from the drawing.

WHY THE TAIL IS BUILT FROM NAMES, NOT FROM PIXELS. The Experimental Section
names every reactant and gives a gram/millimole pair for each. Dividing one by
the other recovers each reagent's molecular weight to about 1%, which is an
external, arithmetic check on the identity we assign -- and it is independent of
any drawing. The hexyl chain in particular is pinned by "1-bromohexane (6.8 g,
41.2 mmol)" and never has to be counted off a raster. CLAUDE.md section 5:
"name tokens are not chain lengths" -- here the chain comes from a named
reagent with a mass, which is a different and stronger thing than a name token.

Writes 01_blocks.csv and 01_blocks_qc.json.
"""
import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))

from rdkit import Chem, RDLogger                                    # noqa: E402
from rdkit.Chem.Descriptors import ExactMolWt                       # noqa: E402
from rdkit.Chem.rdMolDescriptors import CalcMolFormula              # noqa: E402
import assemble                                                     # noqa: E402

RDLogger.DisableLog("rdApp.*")

# ----------------------------------------------------------------------------
# HEADS -- read off Figure 1C (nihms-1799810-f0001.jpg, panel C) at 3-4x.
# `arms` is NOT read off the drawing: it is derived in stage03 from the
# molecular formula Figure 2F prints for each lipidoid. It is absent here on
# purpose.
# ----------------------------------------------------------------------------
HEADS = [
    # id,  SMILES,                    chemical name,                               what Figure 1C draws
    ("75",  "CCN(CC)CCN",
     "N,N-diethylethylenediamine",
     "diethylamino-N, two-carbon linker, terminal NH2"),
    ("76",  "NCCN1CCCC1",
     "1-(2-aminoethyl)pyrrolidine",
     "pyrrolidine N, two-carbon linker, terminal NH2"),
    ("77",  "CN1CCCC1CCN",
     "2-(2-aminoethyl)-1-methylpyrrolidine",
     "N-methylpyrrolidine substituted at ring C2, two-carbon linker, terminal NH2"),
    ("78",  "CN(C)CCN",
     "N,N-dimethylethylenediamine",
     "dimethylamino-N, two-carbon linker, terminal NH2"),
    ("80",  "CN(C)CCCN",
     "N,N-dimethyl-1,3-propanediamine",
     "dimethylamino-N, three-carbon linker, terminal NH2"),
    ("81",  "CCN(CC)CCCN",
     "N,N-diethyl-1,3-propanediamine",
     "diethylamino-N, three-carbon linker, terminal NH2"),
    ("82",  "NCCN1CCCCC1",
     "1-(2-aminoethyl)piperidine",
     "piperidine N, two-carbon linker, terminal NH2"),
    ("93",  "NCCCn1ccnc1",
     "1-(3-aminopropyl)imidazole",
     "imidazole N1, three-carbon linker, terminal NH2"),
    ("113", "CN(CCN)CCN",
     "2,2'-(methylimino)bis(ethylamine)  [N-methyl-N,N-bis(2-aminoethyl)amine]",
     "central N bearing a methyl stub drawn downward and two 2-aminoethyl arms"),
    ("306", "CN(CCCN)CCCN",
     "3,3'-diamino-N-methyldipropylamine",
     "central N bearing a methyl stub drawn downward and two 3-aminopropyl arms"),
    ("400", "CN(C)CCCNCCCN",
     "N,N-dimethyldipropylenetriamine",
     "dimethylamino-N, propyl, an explicitly drawn N-H, propyl, terminal NH2"),
]

# ----------------------------------------------------------------------------
# TAIL -- O16CBA.
#
# Built as a sequence of three NAMED reagents, each with its own printed
# gram/millimole pair, exactly as the Experimental Section describes it:
#
#   (a) 4-hydroxy-2,6-dimethoxybenzaldehyde + 1-bromohexane / K2CO3
#           -> HexDMBA          (Williamson ether on the 4-OH)
#   (b) HexDMBA + 1,1,1-tris(hydroxymethyl)ethane / p-TsOH, 5A sieves
#           -> HexDMBAH         (cyclic acetal; a 1,3-dioxane closes onto TWO
#                                of the three CH2OH groups, leaving the third)
#   (c) HexDMBAH + acryloyl chloride / Et3N
#           -> O16CBA           (acrylate ester on the remaining CH2OH)
#
# WHICH TWO OXYGENS CLOSE THE RING is the question the run brief flags as
# invisible to molecular formula. It is settled here by the reagent, not by a
# guess: 1,1,1-tris(hydroxymethyl)ethane is C(CH3)(CH2OH)3, so every pair of
# its three CH2OH groups is related by the molecule's own C3 symmetry and all
# three pairings give ONE constitution -- a 2-aryl-5-methyl-5-(hydroxymethyl)-
# 1,3-dioxane. There is no regiochemical choice to make. That symmetry claim is
# not asserted; stage03 proves it by building all three pairings and comparing.
# ----------------------------------------------------------------------------
SYNTHESIS = [
    # reagent,                               SMILES,                       grams, mmol
    ("4-hydroxy-2,6-dimethoxybenzaldehyde",  "O=Cc1c(OC)cc(O)cc1OC",       5.1,  28.0),
    ("potassium carbonate",                  None,                         4.6,  33.3),
    ("1-bromohexane",                        "CCCCCCBr",                   6.8,  41.2),
    ("HexDMBA (isolated)",                   "O=Cc1c(OC)cc(OCCCCCC)cc1OC", 5.2,  19.5),
    ("1,1,1-tris(hydroxymethyl)ethane",      "OCC(C)(CO)CO",               6.7,  55.9),
    ("p-toluenesulfonic acid",               "Cc1ccc(cc1)S(=O)(=O)O",      0.44, 2.56),
    ("HexDMBAH (isolated)",
     "OCC1(C)COC(-c2c(OC)cc(OCCCCCC)cc2OC)OC1",                            5.4,  14.7),
    ("triethylamine",                        "CCN(CC)CC",                  2.2,  22.1),
    ("acryloyl chloride",                    "C=CC(=O)Cl",                 2.00, 22.1),
    ("O16CBA (isolated)",
     "C=CC(=O)OCC1(C)COC(-c2c(OC)cc(OCCCCCC)cc2OC)OC1",                    4.5,  None),
]

# potassium carbonate is inorganic; RDKit's organic MW is not the right
# comparison object, so its expected MW is stated here as a literal from the
# periodic table and the check is done against that.
K2CO3_MW = 138.205

O16CBA = "C=CC(=O)OCC1(C)COC(-c2c(OC)cc(OCCCCCC)cc2OC)OC1"


def avg_mw(smiles):
    from rdkit.Chem.Descriptors import MolWt
    return MolWt(Chem.MolFromSmiles(smiles))


def main():
    qc = {"_what": "stage 01 QC for XQ_2020: the eleven amine heads and the "
                   "O16CBA tail, each checked against evidence outside the "
                   "drawing it was read from."}

    # ---- every head parses, is canonical, and has the formula its name implies
    rows, head_qc = [], []
    for hid, smi, name, drawn in HEADS:
        m = Chem.MolFromSmiles(smi)
        if m is None:
            raise SystemExit(f"head {hid}: unparseable")
        canon = Chem.MolToSmiles(m)
        nh = sum(1 for a in m.GetAtoms()
                 if a.GetSymbol() == "N" and a.GetTotalNumHs() > 0)
        slots = sum(a.GetTotalNumHs() for a in m.GetAtoms()
                    if a.GetSymbol() == "N" and a.GetTotalNumHs() > 0)
        head_qc.append({
            "id": hid, "name": name, "formula": CalcMolFormula(m),
            "exact_mass": round(ExactMolWt(m), 4),
            "n_NH_nitrogens": nh,
            "n_michael_slots": slots,
            "stereocentres_legacy": assemble.stereocentres(m),
            "stereo_new_perception": [
                {"atom": int(e.centeredOn), "type": str(e.type),
                 "specified": str(e.specified)}
                for e in Chem.FindPotentialStereo(m)],
        })
        rows.append({
            "block_id": hid, "block_type": "amine_head", "SMILES": canon,
            "formula": CalcMolFormula(m),
            "exact_mass": f"{ExactMolWt(m):.4f}",
            "chemical_name": name,
            "n_michael_slots": slots,
            "source": "Figure 1C (nihms-1799810-f0001.jpg panel C)",
            "what_is_drawn": drawn,
        })
    qc["heads"] = head_qc

    # ---- the synthesis: every named reagent's MW vs its own printed g/mmol
    syn = []
    worst = 0.0
    for name, smi, g, mmol in SYNTHESIS:
        if mmol is None:
            syn.append({"reagent": name, "note": "yield only, no mmol printed"})
            continue
        implied = g / (mmol / 1000.0)
        ours = K2CO3_MW if smi is None else avg_mw(smi)
        dev = abs(implied - ours) / ours * 100.0
        worst = max(worst, dev)
        syn.append({"reagent": name, "grams": g, "mmol": mmol,
                    "MW_implied_by_their_own_g_per_mmol": round(implied, 2),
                    "MW_of_the_structure_we_assign": round(ours, 2),
                    "deviation_pct": round(dev, 2)})
    qc["named_reagent_mw_crosscheck"] = {
        "_what": "The Experimental Section prints a gram and a millimole figure "
                 "for every reagent. g/mmol recovers the MW the AUTHORS were "
                 "working with. Comparing it to the MW of the structure WE "
                 "assign tests the identity of each block without using any "
                 "drawing at all. This is the check that pins the hexyl chain "
                 "length without counting a single pixel.",
        "n_checked": sum(1 for s in syn if "deviation_pct" in s),
        "worst_deviation_pct": round(worst, 2),
        "tolerance_pct": 2.5,
        "why_that_tolerance": "the grams are printed to 2 significant figures "
                              "('2.2 g', '6.8 g'), so a 1-2% rounding spread is "
                              "expected and is not evidence of a wrong structure.",
        "pass": worst <= 2.5,
        "rows": syn,
    }
    if worst > 2.5:
        raise SystemExit(f"stage01: reagent MW cross-check FAILED at {worst:.2f}%")

    # ---- the tail itself
    t = Chem.MolFromSmiles(O16CBA)
    rows.append({
        "block_id": "O16CBA", "block_type": "acrylate_tail",
        "SMILES": Chem.MolToSmiles(t), "formula": CalcMolFormula(t),
        "exact_mass": f"{ExactMolWt(t):.4f}",
        "chemical_name": "2-(4-(hexyloxy)-2,6-dimethoxyphenyl)-5-methyl-1,3-"
                         "dioxan-5-yl)methyl prop-2-enoate",
        "n_michael_slots": 0,
        "source": "Figure 1B + Figure S1 + the Experimental Section's named reagents",
        "what_is_drawn": "acrylate ester on a CH2 of a 5-methyl-1,3-dioxane "
                         "whose C2 bears 4-hexyloxy-2,6-dimethoxyphenyl",
    })
    qc["tail"] = {
        "SMILES": Chem.MolToSmiles(t),
        "formula": CalcMolFormula(t),
        "exact_mass": round(ExactMolWt(t), 4),
        "stereocentres_seen_by_assemble_stereocentres": assemble.stereocentres(t),
        "stereocentres_seen_by_FindPotentialStereo": [
            {"atom": int(e.centeredOn), "type": str(e.type),
             "specified": str(e.specified)}
            for e in Chem.FindPotentialStereo(t)],
        "THE_STEREO_FINDING": (
            "The two perceptions DISAGREE, and the disagreement is the point. "
            "assemble.stereocentres() -- which is what assemble.attach's "
            "check_stereo_preserved relies on -- returns {} for this tail. "
            "Chem.FindPotentialStereo returns TWO unspecified tetrahedral "
            "centres: the acetal carbon C2 and the quaternary ring carbon C5. "
            "They are para-stereocentres: neither is a CIP centre alone, "
            "because each sees two constitutionally identical ring branches, "
            "but TOGETHER they define the cis/trans relationship between the "
            "C2-aryl group and the C5 substituents. The legacy "
            "flagPossibleStereoCenters perception does not see that pair. "
            "So curation/lib's stereo guard is BLIND to the one stereochemical "
            "feature this paper's chemistry turns on, and the guard passing is "
            "not evidence here. Recorded rather than worked around; see "
            "DECISIONS.md D-2."),
    }

    with open(os.path.join(HERE, "01_blocks.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    with open(os.path.join(HERE, "01_blocks_qc.json"), "w") as f:
        json.dump(qc, f, indent=1)
    print(f"stage01: {len(rows)} blocks "
          f"({len(HEADS)} heads + 1 tail); worst reagent MW deviation "
          f"{worst:.2f}%")
    return 0


if __name__ == "__main__":
    sys.exit(main())
