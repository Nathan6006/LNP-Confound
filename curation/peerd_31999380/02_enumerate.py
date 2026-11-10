#!/usr/bin/env python3
"""Enumerate the 14 ionizable lipids of Ramishetti et al., Adv Mater 2020
(PMID 31999380) from the written synthesis in the Supporting Information.

RUN:  python3 -I 02_enumerate.py

WHERE THE STRUCTURES COME FROM. Every one is built from the SI's own prose:
the named product of each step (e.g. "(2-(2-(dimethylamino)ethyl)hydrazine-
1,1-diyl)bis(heptane-7,1-diyl) bis(decanoate) (Lipid 3)", SI p.7), the reagents
named in that paragraph, and the 1H NMR integrals that pin the substituent
counts. Nothing here is read off a figure: pages 3-36 carry native-vector
reaction schemes, but the SI SPELLS OUT each product's IUPAC name, which is a
symbolic source and therefore strictly better than tracing the drawing
(PIPELINE_V2 section 3: prefer symbolic over pixel).

HOW THEY ARE BUILT. Through curation/lib/assemble only -- `attach` for every
C-N and C-O bond and `esterify` for every ester, so the acid loses its OH and
the alcohol keeps its oxygen, and so every build runs
check_stereo_preserved + check_no_radicals. There is no RWMol, no RemoveAtom
and no string surgery on a canonical SMILES anywhere in this file. That matters
because the lipids here carry 2-4 cis double bonds whose geometry is exactly
what a careless assembly drops (error class 48 / F-013's mechanism).

ASSEMBLY ORDER IS NOT SYNTHETIC ORDER. For Lipid 3 and Lipid 5 the arms are
attached before the dimethylaminoethyl group, which is the reverse of the
paper's sequence. The product is the same molecule; the order was chosen so
that every intermediate has exactly ONE nitrogen matching the SMARTS used to
place the next bond, and `remap` below refuses when it does not. Mirroring the
synthetic order instead would have required picking between two chemically
equivalent nitrogens by match order, which is the argmin-by-row-order shape
that heatmap.match refuses for colours and that this file refuses for atoms.
"""
import csv
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "lib"))

import assemble                                            # noqa: E402
from errors import StructureError                          # noqa: E402
from rdkit import Chem, RDLogger                           # noqa: E402
from rdkit.Chem.Descriptors import ExactMolWt, MolWt       # noqa: E402
from rdkit.Chem.rdMolDescriptors import CalcMolFormula     # noqa: E402

RDLogger.DisableLog("rdApp.*")
HERE = os.path.dirname(os.path.abspath(__file__))


# --------------------------------------------------------------------------
# Caller-side glue. NOT an assembler: it forms no bond and deletes no atom.
# It only marks WHICH atom the next lib call should bond to, and it refuses
# when that atom is not unique -- so an ambiguous site becomes a crash, never
# a silent choice by match order.
# --------------------------------------------------------------------------
def remap(smiles, smarts, mapno):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise StructureError(f"remap: unparseable input {smiles!r}")
    patt = Chem.MolFromSmarts(smarts)
    hits = [h[0] for h in mol.GetSubstructMatches(patt)]
    if len(hits) != 1:
        raise StructureError(
            f"remap: SMARTS {smarts!r} matched {len(hits)} atoms in "
            f"{smiles!r}, need exactly 1. Two matches means the next bond "
            "would be placed by match order, which is not a chemical "
            "statement. Re-order the assembly so the site is unique.")
    mol.GetAtomWithIdx(hits[0]).SetAtomMapNum(mapno)
    return Chem.MolToSmiles(mol)


# N with exactly one H and no H, used to walk a nitrogen through alkylation.
N_1H = "[NX3;H1]"
N_2H = "[NX3;H2]"
OH_ONLY = "[OX2H1]"


# --------------------------------------------------------------------------
# BLOCKS. Each is a source literal typed from the SI's own compound name.
# `[...:1]` marks the atom that bonds to the core; `[...:2]` the core atom.
# --------------------------------------------------------------------------
TAILS = {
    # (9Z,12Z)-octadeca-9,12-dien-1-yl, from "linoleic alcohol" -> "linoleic
    # aldehyde" -> reductive amination (SI p.16). 18 C, cis,cis at 9 and 12.
    "linoleyl": {
        "smiles": r"CCCCC/C=C\C/C=C\CCCCCCC[CH2:1]",
        "from": "linoleyl alcohol via linoleic aldehyde, SI p.16",
        "formula_of_group": "C18H33",
    },
    # 7-(decanoyloxy)heptyl, from "7-oxoheptyl decanoate" (SI p.6-7):
    # 1,7-heptanediol -> PCC -> decanoic acid/EDC.
    "heptyl_decanoate": {
        "smiles": "CCCCCCCCCC(=O)OCCCCCC[CH2:1]",
        "from": "7-oxoheptyl decanoate, SI p.6-7",
        "formula_of_group": "C17H33O2",
    },
    # 6-(((Z)-non-2-en-1-yl)oxy-6-oxo)hexyl, from "di((Z)-non-2-en-1-yl)
    # 6,6'-...dihexanoate" (SI p.10-12). cis-2-nonen-1-ol ester of a hexanoate
    # whose C6 carries the nitrogen.
    "hexanoate_non2enyl": {
        "smiles": r"CCCCCC/C=C\COC(=O)CCCC[CH2:1]",
        "from": "cis-2-nonen-1-ol hexanoate arm, SI p.10-12",
        "formula_of_group": "C15H27O2",
    },
    # Same hexanoate arm but esterified with cis-3-nonen-1-ol (SI p.29-31).
    # The DOUBLE BOND POSITION IS THE ONLY DIFFERENCE from the Lipid 4 arm.
    "hexanoate_non3enyl": {
        "smiles": r"CCCCC/C=C\CCOC(=O)CCCC[CH2:1]",
        "from": "cis-3-nonen-1-ol hexanoate arm, SI p.29-31",
        "formula_of_group": "C15H27O2",
    },
    # 10-((tridecan-7-yl)oxy)-10-oxodecyl, from "tridecan-7-yl 10-oxodecanoate"
    # (SI p.19-20). The branch point is C7 of a C13 chain: a symmetric
    # di-hexyl carbinol from tridecan-7-one/NaBH4.
    "decanoate_tridecan7yl": {
        "smiles": "CCCCCCC(CCCCCC)OC(=O)CCCCCCCC[CH2:1]",
        "from": "tridecan-7-yl 10-oxodecanoate, SI p.19-20",
        "formula_of_group": "C23H45O2",
    },
    # 8-((2-hexyldecanoyl)oxy)octyl, from "8-oxooctyl 2-hexyldecanoate"
    # (SI p.33-35). 1,8-octanediol backbone, branched C16 acid.
    "octyl_2hexyldecanoate": {
        "smiles": "CCCCCCCCC(CCCCCC)C(=O)OCCCCCCC[CH2:1]",
        "from": "8-oxooctyl 2-hexyldecanoate, SI p.33-35",
        "formula_of_group": "C24H47O2",
    },
    # 6-((2-hexyldecanoyl)oxy)hexyl, from the Lipid 5 diol + 2-hexyldecanoic
    # acid (SI p.13-15). Same branched acid, C6 instead of C8 spacer.
    "hexyl_2hexyldecanoate": {
        "smiles": "CCCCCCCCC(CCCCCC)C(=O)OCCCCC[CH2:1]",
        "from": "hexane-6,1-diyl 2-hexyldecanoate arm, SI p.13-15",
        "formula_of_group": "C22H43O2",
    },
}

HEADS = {
    # 4-(dimethylamino)butanoyl. Written as the ALDEHYDE form so `attach` can
    # bond the carbonyl carbon and hand the H count back to RDKit; the lib
    # clears noImplicit on the mapped atom, so the product is the amide/ester
    # carbonyl, not an aldehyde. Used for an AMIDE (lipids 1, 4).
    "DMBA_acyl": {"smiles": "O=[C:1]CCCN(C)C",
                  "from": "N,N-dimethylaminobutyric acid, used throughout"},
    # The same head as a free acid, for `esterify` (lipids 6,7,8,11,12,13,14).
    "DMBA_acid": {"smiles": "CN(C)CCCC(=O)O",
                  "from": "4-(dimethylamino)butanoic acid"},
    # 3-(4-methylpiperazin-1-yl)propanoyl / acid (lipids 2, 10).
    "MPzP_acyl": {"smiles": "O=[C:1]CCN1CCN(C)CC1",
                  "from": "N-methyl piperazine propionic acid, SI p.5"},
    "MPzP_acid": {"smiles": "CN1CCN(CCC(=O)O)CC1",
                  "from": "N-methyl piperazine propionic acid, SI p.5"},
    # 2-(dimethylamino)ethyl, placed on a hydrazine N (lipids 3, 5).
    "DMAE": {"smiles": "CN(C)C[CH2:1]",
             "from": "dimethylaminoethyl hydrazine, SI p.7 / p.14"},
    # 4-(dimethylamino)butyl, placed on an ETHER oxygen (lipid 9).
    "DMAB": {"smiles": "CN(C)CCC[CH2:1]",
             "from": "4-(dimethylamino)butan-1-ol, SI p.27"},
}

SCAFFOLDS = {
    # hydrazine; map 2 on ONE nitrogen so the first alkylation is unambiguous
    # even though the two nitrogens are equivalent in the free base.
    "hydrazine": "[NH2:2]N",
    # hydroxylamine, N-alkylated twice then O-acylated (lipids 6, 7).
    "hydroxylamine": "[NH2:2]O",
    # ethanolamine, N-alkylated twice then O-acylated (lipids 8,10,11,12,13,14)
    "ethanolamine": "[NH2:2]CCO",
}


# --------------------------------------------------------------------------
# Assembly routes
# --------------------------------------------------------------------------
def dialkyl_scaffold(scaffold, tail, label):
    """Put two identical tails on the scaffold's mapped nitrogen."""
    mono = assemble.attach(TAILS[tail]["smiles"], SCAFFOLDS[scaffold],
                           a_map=1, b_map=2, label=f"{label} arm1")
    di = assemble.attach(TAILS[tail]["smiles"], remap(mono, N_1H, 2),
                         a_map=1, b_map=2, label=f"{label} arm2")
    return mono, di


def build(lipid):
    """Return (smiles, [(step, product)]) for one lipid. Raises on anything
    the lib refuses; nothing here falls through to a literal."""
    trail = []

    def step(name, smi):
        trail.append((name, smi))
        return smi

    # ---- hydrazide family: N,N-dialkylhydrazine, N'-acylated -------------
    if lipid in ("L1", "L2", "L4"):
        tail = "linoleyl" if lipid in ("L1", "L2") else "hexanoate_non2enyl"
        mono, di = dialkyl_scaffold("hydrazine", tail, lipid)
        step("N-monoalkylhydrazine", mono)
        step("N,N-dialkylhydrazine", di)
        head = "MPzP_acyl" if lipid == "L2" else "DMBA_acyl"
        return step("acylhydrazide", assemble.attach(
            HEADS[head]["smiles"], remap(di, N_2H, 2),
            a_map=1, b_map=2, label=f"{lipid} amide")), trail

    # ---- hydrazine family with a dimethylaminoethyl on nitrogen ----------
    if lipid == "L3":
        # N1 carries both ester arms (the "hydrazine-1,1-diyl" of the name);
        # N2 keeps one H and takes the dimethylaminoethyl group.
        mono, di = dialkyl_scaffold("hydrazine", "heptyl_decanoate", lipid)
        step("N-monoalkylhydrazine", mono)
        step("N,N-dialkylhydrazine (1,1-diyl)", di)
        return step("N'-(dimethylaminoethyl)", assemble.attach(
            HEADS["DMAE"]["smiles"], remap(di, N_2H, 2),
            a_map=1, b_map=2, label="L3 DMAE")), trail

    if lipid == "L5":
        # "hydrazine-1,2-diyl": ONE arm on each nitrogen. The DMAE goes on
        # first precisely so each subsequent nitrogen is unique to SMARTS.
        dmae = assemble.attach(HEADS["DMAE"]["smiles"],
                               SCAFFOLDS["hydrazine"],
                               a_map=1, b_map=2, label="L5 DMAE")
        step("(dimethylaminoethyl)hydrazine", dmae)
        a1 = assemble.attach(TAILS["hexyl_2hexyldecanoate"]["smiles"],
                             remap(dmae, N_1H, 2),
                             a_map=1, b_map=2, label="L5 arm on N1")
        step("arm on N1 (bearing DMAE)", a1)
        return step("arm on N2", assemble.attach(
            TAILS["hexyl_2hexyldecanoate"]["smiles"], remap(a1, N_2H, 2),
            a_map=1, b_map=2, label="L5 arm on N2")), trail

    # ---- hydroxylamine family: N,N-dialkyl, O-acylated -------------------
    if lipid in ("L6", "L7"):
        tail = "linoleyl" if lipid == "L6" else "decanoate_tridecan7yl"
        mono, di = dialkyl_scaffold("hydroxylamine", tail, lipid)
        step("N-monoalkylhydroxylamine", mono)
        step("N,N-dialkylhydroxylamine", di)
        return step("O-acyl (ester)", assemble.esterify(
            HEADS["DMBA_acid"]["smiles"], di, label=f"{lipid} O-acyl")), trail

    # ---- aminoethanol family: N,N-dialkyl, O-acylated --------------------
    if lipid in ("L8", "L10", "L11", "L12", "L13", "L14"):
        tail = {"L8": "linoleyl", "L10": "linoleyl",
                "L11": "heptyl_decanoate", "L12": "hexanoate_non3enyl",
                "L13": "decanoate_tridecan7yl",
                "L14": "octyl_2hexyldecanoate"}[lipid]
        mono, di = dialkyl_scaffold("ethanolamine", tail, lipid)
        step("N-monoalkylethanolamine", mono)
        step("N,N-dialkylethanolamine", di)
        head = "MPzP_acid" if lipid == "L10" else "DMBA_acid"
        return step("O-acyl (ester)", assemble.esterify(
            HEADS[head]["smiles"], di, label=f"{lipid} O-acyl")), trail

    # ---- the one ether ---------------------------------------------------
    if lipid == "L9":
        mono, di = dialkyl_scaffold("ethanolamine", "linoleyl", lipid)
        step("N-monoalkylethanolamine", mono)
        step("EA-1 (N,N-dilinoleyl aminoethanol)", di)
        return step("ether (via mesylate)", assemble.attach(
            HEADS["DMAB"]["smiles"], remap(di, OH_ONLY, 2),
            a_map=1, b_map=2, label="L9 ether")), trail

    raise StructureError(f"build: no route defined for {lipid!r}")


# --------------------------------------------------------------------------
# The SI's reported ESI-MS, transcribed verbatim with its page.
# NOTE the SI prints NO calculated ("Calc.") column anywhere -- only FOUND
# m/z from a low-resolution ESI instrument. Which theoretical mass each value
# corresponds to is therefore something this curation has to DETERMINE, not
# something it can read. 07_qc.json does that three ways.
# --------------------------------------------------------------------------
ESI_MS = {
    "L1":  {"mz": 642.6,  "label": "(M+1).+", "page": 4},
    "L2":  {"mz": 683.7,  "label": "(M+1)+",  "page": 5,
            "second": {"mz": 342.3, "label": "(M/2+1)+"}},
    "L3":  {"mz": 640.0,  "label": "([M+1]+)", "page": 7},
    "L4":  {"mz": 622.6,  "label": "([M+1]+)", "page": 12},
    "L5":  {"mz": 780.0,  "label": "([M+1]+)", "page": 15},
    "L6":  {"mz": 644.1,  "label": "[M+1]+",   "page": 17,
            "second": {"mz": 643.1, "label": "[M]+"}},
    "L7":  {"mz": 852.2,  "label": "([M+1]+)", "page": 23},
    "L8":  {"mz": 672.0,  "label": "([M+1]+)", "page": 25,
            "second": {"mz": 671.0, "label": "([M]+)"}},
    "L9":  {"mz": 658.0,  "label": "(M+1)+",   "page": 27,
            "second": {"mz": 329.5, "label": "(M/2+1)+"}},
    "L10": {"mz": 712.0,  "label": "([M+1]+)", "page": 28},
    "L11": {"mz": 711.0,  "label": "([M+1]+)", "page": 9},
    "L12": {"mz": 651.8,  "label": "([M+1]+)", "page": 32,
            "second": {"mz": 326.6, "label": "([M/2+1]+)"}},
    "L13": {"mz": 880.12, "label": "([M+1]+)", "page": 22,
            "second": {"mz": 441.75, "label": "([M/2+1]+)"}},
    "L14": {"mz": 908.3,  "label": "([M+1]+)", "page": 36},
}

# The SI's own IUPAC name for each lipid, transcribed from the page given.
# Where the SI names only the precursor ("Lipid 1", "Lipid 2"), that is said so.
IUPAC = {
    "L1": ("not named in the SI; described as N,N-dilinoleyl hydrazide "
           "acylated with N,N-dimethyl aminobutyric acid", 3),
    "L2": ("not named in the SI; N,N-dilinoleyl hydrazide + N-methyl "
           "piperazine propionic acid", 4),
    "L3": ("(2-(2-(dimethylamino)ethyl)hydrazine-1,1-diyl)bis(heptane-7,1-"
           "diyl) bis(decanoate)", 7),
    "L4": ("Di((Z)-non-2-en-1-yl) 6,6'-(2-(4-(dimethylamino)butanoyl)"
           "hydrazine-1,1-diyl)dihexanoate", 12),
    "L5": ("(1-(2-(dimethylamino)ethyl)hydrazine-1,2-diyl)bis(hexane-6,1-"
           "diyl) bis(2-hexyldecanoate)", 15),
    "L6": ("not named in the SI; HyAM-1 (N,N-dilinoleyl-hydroxylamine) "
           "O-acylated with 4-(dimethylamino)butanoic acid", 17),
    "L7": ("di(tridecan-7-yl) 10,10'-(((4-(dimethylamino)butanoyl)oxy)"
           "azanediyl)bis(decanoate)", 23),
    "L8": ("not named in the SI; EA-1 (N,N-dilinoleyl-aminoethanol) "
           "O-acylated with 4-(dimethylamino)butanoic acid", 25),
    "L9": ("not named in the SI; EA-1 mesylate displaced by "
           "4-(dimethylamino)butan-1-olate", 27),
    "L10": ("not named in the SI; EA-1 O-acylated with N-methyl "
            "piperazinepropanoic acid", 28),
    "L11": ("((2-((4-(dimethylamino)butanoyl)oxy)ethyl)azanediyl)bis(heptane-"
            "7,1-diyl) bis(decanoate)", 9),
    "L12": ("di((Z)-non-3-en-1-yl) 6,6'-((2-((4-dimethylamino)butanoyl)oxy)"
            "ethyl)azanediyl)dihexanoate", 31),
    "L13": ("di(tridecan-7-yl) 10,10'-((2-((4-(dimethylamino)butanoyl)oxy)"
            "ethyl)azanediyl)bis(decanoate)", 21),
    "L14": ("((2-((4-(dimethylamino)butanoyl)oxy)ethyl)azanediyl)bis(octane-"
            "8,1-diyl) bis(2-hexyldecanoate)", 36),
}

ALIASES = {"L13": "EA-402", "L8": "from EA-1", "L6": "from HyAM-1",
           "L9": "from EA-1 via EA-OMs", "L10": "from EA-1"}

ORDER = [f"L{i}" for i in range(1, 15)]


def main():
    rows, blocks, trails = [], [], {}
    for lid in ORDER:
        smi, trail = build(lid)
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            raise StructureError(f"{lid}: product does not re-parse")
        trails[lid] = trail
        name, page = IUPAC[lid]
        rows.append({
            "lipid_id": lid,
            "IL_name": f"Lipid {int(lid[1:])}",
            "IL_SMILES": Chem.MolToSmiles(mol),
            "si_iupac_name": name,
            "si_name_page": page,
            "alias": ALIASES.get(lid, ""),
            "formula": CalcMolFormula(mol),
            "exact_mass_monoisotopic": round(ExactMolWt(mol), 4),
            "average_mw": round(MolWt(mol), 2),
            "n_stereocentres": len(assemble.stereocentres(mol)),
            "n_cis_double_bonds": sum(
                1 for b in mol.GetBonds()
                if b.GetStereo() == Chem.BondStereo.STEREOZ),
            "esi_ms_found_mz": ESI_MS[lid]["mz"],
            "esi_ms_label_as_printed": ESI_MS[lid]["label"],
            "esi_ms_page": ESI_MS[lid]["page"],
            "inchikey": Chem.MolToInchiKey(mol),
            "inchikey_skeleton": Chem.MolToInchiKey(mol).split("-")[0],
            "assembly_steps": len(trail),
        })

    for group, d in (("tail", TAILS), ("head", HEADS)):
        for k, v in d.items():
            blocks.append({"block_role": group, "block_name": k,
                           "block_SMILES": v["smiles"],
                           "source_in_SI": v.get("from", ""),
                           "group_formula": v.get("formula_of_group", "")})
    for k, v in SCAFFOLDS.items():
        blocks.append({"block_role": "scaffold", "block_name": k,
                       "block_SMILES": v,
                       "source_in_SI": "named starting material",
                       "group_formula": ""})

    with open(os.path.join(HERE, "01_blocks.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(blocks[0].keys()))
        w.writeheader()
        w.writerows(blocks)
    with open(os.path.join(HERE, "03_enumerated.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    with open(os.path.join(HERE, "03_assembly_trail.json"), "w") as f:
        json.dump(trails, f, indent=1)

    print(f"{len(blocks)} blocks -> 01_blocks.csv")
    print(f"{len(rows)} lipids  -> 03_enumerated.csv")
    for r in rows:
        print(f"  {r['lipid_id']:>3} {r['formula']:>14} "
              f"mono {r['exact_mass_monoisotopic']:>9} "
              f"avg {r['average_mw']:>8}  ESI {r['esi_ms_found_mz']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
