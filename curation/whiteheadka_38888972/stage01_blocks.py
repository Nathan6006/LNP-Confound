#!/usr/bin/env python3
"""The 29 amine heads and the 20 acrylate tails -> 01_blocks.csv.

RULE THIS FILE OBEYS. Every structure here is derived FROM THE CHEMICAL NAME
printed in Table S1 / Table S2 of the Supporting Information, never from the
code token and never from the drawing. `Oi10` is not "isomer 10"; `O1,1,8` is
not a chain-length triple (it is 3,7-dimethyloctyl, a C10); `O8,12` is
2-octyldodecyl, a C20. CLAUDE.md section 5: name tokens are not chain lengths.

The DRAWING (Figure 1B, 1C) is used only as a cross-check, and where it
disagrees with the name that disagreement is recorded in 08_contradictions.csv
and the NAME wins -- because these are catalogue reagents bought under those
names from Sigma / TCI / Enamine / MP Bio, and a vendor ships what the name
says. Three head drawings disagree; see the contradictions file.

The acrylate tails are built with curation/lib assemble.esterify, from acrylic
acid plus the ALCOHOL the Table S2 name names. The acid loses its OH and the
alcohol keeps its oxygen, which is what the chemistry does and what keeps the
F-013 mechanism out of reach.

Run: python3 -I stage01_blocks.py
"""
import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))

from rdkit import Chem                                          # noqa: E402
from rdkit.Chem import Descriptors, rdMolDescriptors            # noqa: E402

import assemble                                                 # noqa: E402

ACRYLIC_ACID = "C=CC(=O)O"
# Propanoic acid with an atom map on C3. Esterified with the same alcohol this
# gives the MICHAEL ADDUCT fragment -- the tail as it appears once it has added
# across an N-H -- with the attachment carbon tagged. assemble.attach resolves
# by map number, never by atom index.
MAPPED_PROPANOIC_ACID = "[CH3:2]CC(=O)O"

# ---------------------------------------------------------------------------
# TABLE S1 -- 29 amine heads. (code, Table S1 name verbatim, vendor, SMILES
# written from that name, one-line derivation).
# ---------------------------------------------------------------------------
HEADS = [
    ("60", "N,N'-Dimethylethylenediamine", "Sigma", "CNCCNC",
     "ethane-1,2-diamine with one methyl on each N"),
    ("80", "3-(Dimethylamino)-1-propylamine", "Sigma", "CN(C)CCCN",
     "propan-1-amine with NMe2 on C3"),
    ("81", "3-(Diethylamino)propylamine", "Sigma", "CCN(CC)CCCN",
     "propan-1-amine with NEt2 on C3"),
    ("110", "Tris(2-aminoethyl)amine", "Sigma", "NCCN(CCN)CCN",
     "one central N bearing three 2-aminoethyl arms"),
    ("111", "Tetraethylenepentamine", "Sigma", "NCCNCCNCCNCCN",
     "five N in a chain joined by four ethylene units"),
    ("113", "N'-methyl-2,2'-diaminodiethylamine", "TCI", "CN(CCN)CCN",
     "bis(2-aminoethyl)amine METHYLATED on the central N (CAS 4097-88-5). "
     "Figure 1B draws the central N with an H and no methyl -- see "
     "08_contradictions.csv C-1."),
    ("122", "1,2-Bis(3-aminopropylamino)ethane", "Sigma", "NCCCNCCNCCCN",
     "ethane-1,2-diyl bridging two 3-aminopropylamino groups"),
    ("123", "N,N-Diethyldiethylenetriamine", "Sigma", "CCN(CC)CCNCCN",
     "diethylenetriamine with both ethyls on one terminal N"),
    ("174", "Bis(3-aminopropyl)amine", "Sigma", "NCCCNCCCN",
     "one secondary N bearing two 3-aminopropyl arms"),
    ("200", "N1-(2-(4-(2-aminoethyl)piperazin-1-yl)ethyl)ethane-1,2-diamine",
     "Enamine", "NCCNCCN1CCN(CCN)CC1",
     "ethane-1,2-diamine whose N1 carries a 2-(piperazin-1-yl)ETHYL arm; the "
     "piperazine's 4-position carries 2-aminoethyl. Figure 1B draws ONE "
     "methylene between the piperazine N and the NH, which would be an "
     "aminal -- see 08_contradictions.csv C-2."),
    ("205", "1-(2-Aminoethyl)piperazine", "Sigma", "NCCN1CCNCC1",
     "piperazine with a 2-aminoethyl on N1"),
    ("211", "Hexamethylenediamine", "Sigma", "NCCCCCCN", "hexane-1,6-diamine"),
    ("303", "N,N-Dimethyldipropylenetriamine", "Sigma", "CN(C)CCCNCCCN",
     "dipropylenetriamine (= bis(3-aminopropyl)amine) with both methyls on "
     "one terminal N (CAS 10563-29-8). Figure 1B draws ONE methylene between "
     "the NH and the NMe2 -- see 08_contradictions.csv C-3."),
    ("304", "Tris[2-(methylamino)ethyl]amine", "Sigma", "CNCCN(CCNC)CCNC",
     "one central N bearing three 2-(methylamino)ethyl arms"),
    ("305", "1,2-Diaminopropane, 99%", "Sigma", "CC(N)CN",
     "propane-1,2-diamine"),
    ("306", "3,3'-Diamino-N-methyldipropylamine", "Sigma", "NCCCN(C)CCCN",
     "one central N-methyl bearing two 3-aminopropyl arms"),
    ("313", "1,4-Bis(3-aminopropyl)piperazine", "TCI", "NCCCN1CCN(CCCN)CC1",
     "piperazine with a 3-aminopropyl on each ring N"),
    ("347", "m-Xylylenediamine", "Sigma", "NCc1cccc(CN)c1",
     "benzene-1,3-diyl bearing two aminomethyl groups. Table S1 and the "
     "Figure 2A row label both read 347; Figure 1B labels this structure "
     "374 -- see 08_contradictions.csv C-4."),
    ("400", "2-methylpentane-1,5-diamine", "Sigma", "CC(CN)CCCN", "as named"),
    ("401", "1,3-diaminopropan-2-ol", "Sigma", "NCC(O)CN",
     "as named; the secondary OH is NOT a Michael acceptor site -- the "
     "AA_2008 head 87 precedent in the corpus keeps its hydroxyls free"),
    ("402", "4-(2-aminopropan-2-yl)-1-methylcyclohexan-1-amine", "Sigma",
     "CC1(N)CCC(C(C)(C)N)CC1", "as named"),
    ("601", "Tris[2-(isopropylamino)ethyl]amine", "Sigma",
     "CC(C)NCCN(CCNC(C)C)CCNC(C)C",
     "one central N bearing three 2-(isopropylamino)ethyl arms"),
    ("602", "1-Amino-4-methylpiperazine", "Sigma", "CN1CCN(N)CC1",
     "4-methylpiperazin-1-amine; the NH2 sits ON a ring nitrogen"),
    ("500", "N1,N1'-(ethane-1,2-diyl)bis(N1-methylethane-1,2-diamine)",
     "In house", "NCCN(C)CCN(C)CCN",
     "two N1-methylethane-1,2-diamine units joined through their N1 by an "
     "ethane-1,2-diyl bridge"),
    ("501", "N1,N1'-(ethane-1,2-diyl)bis(N1-ethylethane-1,2-diamine)",
     "In house", "NCCN(CC)CCN(CC)CCN", "as 500 with ethyl for methyl"),
    ("503", "N1,N1'-(propane-1,3-diyl)bis(N1-methylethane-1,2-diamine)",
     "In house", "NCCN(C)CCCN(C)CCN", "as 500 with a propane-1,3-diyl bridge"),
    ("509", "N1,N1'-(ethane-1,2-diyl)bis(N1-butylethane-1,2-diamine)",
     "In house", "NCCN(CCCC)CCN(CCCC)CCN", "as 500 with butyl for methyl"),
    ("514", "N1,N1'-(propane-1,3-diyl)bis(N1-ethylethane-1,2-diamine)",
     "In house", "NCCN(CC)CCCN(CC)CCN",
     "as 503 with ethyl for methyl. This is the paper's lead lipidoid head."),
    ("516", "N1,N1'-(but-2-ene-1,4-diyl)bis(N1-ethylethane-1,2-diamine)",
     "In house", "NCCN(CC)CC=CCN(CC)CCN",
     "as 501 with a but-2-ene-1,4-diyl bridge. The name gives no E/Z AND "
     "Figure 1B draws that double bond with a WAVY bond, which is an explicit "
     "statement that the geometry is undefined. Stored with no geometry, "
     "which is what both sources say."),
]

# ---------------------------------------------------------------------------
# TABLE S2 -- 20 acrylate tails. (code, Table S2 name verbatim, vendor,
# ALCOHOL SMILES written from that name, derivation).
# The acrylate itself is BUILT, not typed.
# ---------------------------------------------------------------------------
TAILS = [
    ("O10", "n-Decyl Acrylate", "Sigma", "CCCCCCCCCCO", "decan-1-ol, C10"),
    ("O11", "Undecyl Acrylate", "MP Bio", "CCCCCCCCCCCO", "undecan-1-ol, C11"),
    ("O12", "n-Lauryl Acrylate", "Sigma", "CCCCCCCCCCCCO",
     "lauryl = dodecan-1-ol, C12"),
    ("O13", "Tridecyl Acrylate", "Sigma", "CCCCCCCCCCCCCO",
     "tridecan-1-ol, C13"),
    ("O14", "Tetradecyl Acrylate", "TCI Chemicals", "CCCCCCCCCCCCCCO",
     "tetradecan-1-ol, C14"),
    ("O16", "Hexadecyl Acrylate", "TCI Chemicals", "CCCCCCCCCCCCCCCCO",
     "hexadecan-1-ol, C16"),
    ("O18", "Stearyl Acrylate", "Sigma", "CCCCCCCCCCCCCCCCCCO",
     "stearyl = octadecan-1-ol, C18"),
    ("Oi10", "Isodecyl Acrylate", "Sigma", "CC(C)CCCCCCCO",
     "8-methylnonan-1-ol. 'Isodecyl acrylate' (CAS 1330-61-6) is sold as a "
     "BRANCHED-ISOMER MIXTURE; Figure 1A and 1C both draw it with a terminal "
     "isopropyl, i.e. the 8-methylnonyl isomer, so that is what is stored and "
     "every row carrying it is flagged structure_assumed."),
    ("O8,12", "2-Octyl-Dodecyl Acrylate", "Sigma", "OCC(CCCCCCCC)CCCCCCCCCC",
     "2-octyldodecan-1-ol, C20; the code is a BRANCH DESCRIPTOR, not a chain "
     "length pair"),
    ("O6,10", "2-Hexyl-Decyl Acrylate", "In house", "OCC(CCCCCC)CCCCCCCC",
     "2-hexyldecan-1-ol, C16. This is the tail of the paper's lead lipidoid."),
    ("O1,1,8", "3,7-Dimethyl-Octanyl Acrylate", "Sigma", "OCCC(C)CCCC(C)C",
     "3,7-dimethyloctan-1-ol, C10 in total. The code is a methyl-position "
     "descriptor and contains no chain length at all."),
    ("O4,8", "2-Butyl-Octanyl Acrylate", "Sigma", "OCC(CCCC)CCCCCC",
     "2-butyloctan-1-ol, C12"),
    ("X1", "cis-4-Decenyl Acrylate", "Sigma", "OCCC/C=C\\CCCCC",
     "(Z)-dec-4-en-1-ol"),
    ("X2", "trans-2-Tridecenyl Acrylate", "In house", "OC/C=C/CCCCCCCCCC",
     "(E)-tridec-2-en-1-ol"),
    ("X3", "trans-2-Undecenyl Acrylate", "In house", "OC/C=C/CCCCCCCC",
     "(E)-undec-2-en-1-ol"),
    ("X4", "10-Undecenyl Acrylate", "Sigma", "OCCCCCCCCCC=C",
     "undec-10-en-1-ol; the alkene is TERMINAL so it has no E/Z"),
    ("X5", "trans-3,7-Dimethyl-2,6-Octaneyl Acrylate", "Sigma",
     "OC/C=C(C)/CCC=C(C)C",
     "geraniol, (E)-3,7-dimethylocta-2,6-dien-1-ol. The 6-ene bears two "
     "methyls on C7 and therefore has NO E/Z to assign; only the 2-ene does, "
     "and the name says trans. ('Octaneyl' is the paper's typo for "
     "octadienyl.)"),
    ("X6", "cis-9-Octadecenyl Acrylate", "In house",
     "OCCCCCCCC/C=C\\CCCCCCCC", "oleyl alcohol, (Z)-octadec-9-en-1-ol"),
    ("X7", "cis-9,12-Octadecenyl Acrylate", "In house",
     "OCCCCCCCC/C=C\\C/C=C\\CCCCC",
     "linoleyl alcohol, (9Z,12Z)-octadeca-9,12-dien-1-ol"),
    ("X8", "3,7,11,15-Tetramethyl-2-Hexadecenyl Acrylate", "In house",
     "OC/C=C(C)/CCCC(C)CCCC(C)CCCC(C)C",
     "3,7,11,15-tetramethylhexadec-2-en-1-ol (phytol). The name gives NO E/Z "
     "for the 2-ene while every other X name states cis or trans; Figure 1C "
     "draws it trans, and E is stored on the drawing's authority -- flagged "
     "structure_assumed. C7 and C11 are stereocentres that neither the name "
     "nor the drawing assigns, and they are left unassigned."),
]


def nh_sites(smiles):
    """N atoms carrying at least one H, and the total N-H count.

    Counted by RDKit, never by reading the drawing. The O-H of head 401 is
    deliberately not a site: the AA_2008 head 87 precedent in the corpus
    (NCCCN(CCO)CCO -> two ester tails, both hydroxyls free) shows the acrylate
    adds across N-H only.
    """
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        raise SystemExit(f"unparseable: {smiles}")
    m = Chem.AddHs(m)
    sites, total = 0, 0
    for a in m.GetAtoms():
        if a.GetSymbol() != "N":
            continue
        h = sum(1 for n in a.GetNeighbors() if n.GetSymbol() == "H")
        if h:
            sites += 1
            total += h
    return sites, total


def main():
    rows = []
    qc = {"heads": {}, "tails": {}}

    for code, name, vendor, smi, deriv in HEADS:
        m = Chem.MolFromSmiles(smi)
        if m is None:
            raise SystemExit(f"head {code}: unparseable {smi!r}")
        sites, nh = nh_sites(smi)
        can = Chem.MolToSmiles(m)
        rows.append({
            "block_type": "amine_head", "code": code, "table": "S1",
            "name_verbatim": name, "vendor": vendor,
            "SMILES": can, "reagent_SMILES": can,
            "formula": rdMolDescriptors.CalcMolFormula(m),
            "exact_mass": round(Descriptors.ExactMolWt(m), 4),
            "n_nitrogen": sum(1 for a in m.GetAtoms() if a.GetSymbol() == "N"),
            "n_NH_bonds": nh, "n_NH_sites": sites,
            "inchikey": Chem.MolToInchiKey(m),
            "derivation": deriv,
        })
        qc["heads"][code] = {"n_NH_bonds": nh, "n_NH_sites": sites,
                             "formula": rdMolDescriptors.CalcMolFormula(m)}

    for code, name, vendor, alcohol, deriv in TAILS:
        alc = Chem.MolFromSmiles(alcohol)
        if alc is None:
            raise SystemExit(f"tail {code}: unparseable alcohol {alcohol!r}")
        acrylate = assemble.esterify(ACRYLIC_ACID, alcohol,
                                     label=f"{code} acrylate")
        adduct = assemble.esterify(MAPPED_PROPANOIC_ACID, alcohol,
                                   label=f"{code} Michael adduct")
        # The adduct must be the acrylate with its C=C saturated and the
        # terminal carbon tagged. Check that by DERIVATION, not by eye: strip
        # the map, and the two must differ by exactly one degree of
        # unsaturation on the acyl side.
        am = Chem.MolFromSmiles(acrylate)
        dm = Chem.MolFromSmiles(adduct)
        if dm is None or am is None:
            raise SystemExit(f"tail {code}: build produced unparseable SMILES")
        if (rdMolDescriptors.CalcNumHBD(am) != 0
                or rdMolDescriptors.CalcNumHBD(dm) != 0):
            raise SystemExit(f"tail {code}: a free OH survived esterification")
        n_vinyl = len(am.GetSubstructMatches(
            Chem.MolFromSmarts("[CH2]=[CH][CX3](=O)[OX2]")))
        if n_vinyl != 1:
            raise SystemExit(
                f"tail {code}: built acrylate has {n_vinyl} acrylate groups, "
                "need exactly 1")
        if len(dm.GetSubstructMatches(
                Chem.MolFromSmarts("[CH2]=[CH][CX3](=O)[OX2]"))) != 0:
            raise SystemExit(f"tail {code}: adduct still carries a vinyl")
        n_c_alc = sum(1 for a in alc.GetAtoms() if a.GetSymbol() == "C")
        rows.append({
            "block_type": "acrylate_tail", "code": code, "table": "S2",
            "name_verbatim": name, "vendor": vendor,
            "SMILES": acrylate,
            "reagent_SMILES": acrylate,
            "alcohol_SMILES": Chem.MolToSmiles(alc),
            "formula": rdMolDescriptors.CalcMolFormula(am),
            "exact_mass": round(Descriptors.ExactMolWt(am), 4),
            "adduct_fragment_SMILES": adduct,
            "alcohol_carbons": n_c_alc,
            "n_CC_double_bonds_in_alcohol": len(alc.GetSubstructMatches(
                Chem.MolFromSmarts("[CX3]=[CX3]"))),
            "inchikey": Chem.MolToInchiKey(am),
            "derivation": deriv,
        })
        qc["tails"][code] = {"alcohol_carbons": n_c_alc,
                             "formula": rdMolDescriptors.CalcMolFormula(am),
                             "adduct": adduct}

    cols = ["block_type", "code", "table", "name_verbatim", "vendor",
            "SMILES", "reagent_SMILES", "alcohol_SMILES",
            "adduct_fragment_SMILES", "formula", "exact_mass",
            "alcohol_carbons", "n_CC_double_bonds_in_alcohol",
            "n_nitrogen", "n_NH_bonds", "n_NH_sites", "inchikey",
            "derivation"]
    out = os.path.join(HERE, "01_blocks.csv")
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)

    n_heads = sum(1 for r in rows if r["block_type"] == "amine_head")
    n_tails = sum(1 for r in rows if r["block_type"] == "acrylate_tail")
    # The library size claim, derived here and asserted.
    if n_heads != 29 or n_tails != 20:
        raise SystemExit(f"blocks: {n_heads} heads x {n_tails} tails")
    if n_heads * n_tails != 580:
        raise SystemExit("blocks: product is not 580")
    # Every code must be unique, and every structure distinct.
    codes = [r["code"] for r in rows]
    if len(set(codes)) != len(codes):
        raise SystemExit("blocks: duplicate code")
    for bt in ("amine_head", "acrylate_tail"):
        keys = [r["inchikey"] for r in rows if r["block_type"] == bt]
        if len(set(keys)) != len(keys):
            raise SystemExit(f"blocks: two {bt} blocks are the same molecule")

    qc["n_heads"] = n_heads
    qc["n_tails"] = n_tails
    qc["product"] = n_heads * n_tails
    with open(os.path.join(HERE, "01_blocks_qc.json"), "w") as fh:
        json.dump(qc, fh, indent=1, sort_keys=True)
    print(f"01_blocks.csv: {n_heads} heads x {n_tails} tails = "
          f"{n_heads * n_tails}")
    print("N-H bond counts:", {r["code"]: r["n_NH_bonds"] for r in rows
                               if r["block_type"] == "amine_head"})


if __name__ == "__main__":
    main()
