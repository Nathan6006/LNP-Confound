#!/usr/bin/env python3
"""Enumerate the 15 ionizable lipids of Han et al., ACS Nano 2026 (PMID 41555650).

ROUTE. Every structure is built by AMINE CORE + 5 (or 4) IDENTICAL 2-HYDROXYALKYL
ARMS, one bond at a time, through curation/lib/assemble.attach. No string
concatenation anywhere (a sibling's string-built tail produced a ketone bearing a
free alcohol instead of an acrylate on 11 rows). No RemoveAtom: the epoxide
ring-opening loses no atom -- the oxirane O becomes the arm's OH -- so `attach`
is the correct primitive and the F-013 mechanism cannot arise.

ARM COUNT IS NOT ASSERTED, IT IS COUNTED. The loop attaches one arm per
remaining N-H until no N-H is left. So the arm count is a PROPERTY of the
derived core, and it is then compared against the count stated in the paper's
own IUPAC name (check P3). The two agreeing is evidence; neither is assumed.

AMINE CORES are derived from the SI's systematic names (pp. 2-4) and then
INDEPENDENTLY pinned by the paper's own stoichiometry (check P1), which is why
they are not a guess. See 07_qc.json.
"""
import csv, json, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib"))
import assemble  # curation/lib/assemble.py -- the ONE audited assembler
from rdkit import Chem
from rdkit.Chem import inchi
from rdkit.Chem.Descriptors import MolWt, ExactMolWt
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

HERE = os.path.dirname(os.path.abspath(__file__))

# ---------------------------------------------------------------- amine cores
# Each is the FREE BASE. Written from the SI systematic name, atom by atom.
# `arms_in_name` is what the IUPAC name says, transcribed separately and used
# ONLY as a check -- the builder never reads it.
CORES = {
    "200": dict(
        smiles="NCCN1CCN(CCNCCN)CC1",
        iupac_fragment="2-(4-(2-((2-aminoethyl)(H)amino)ethyl)piperazin-1-yl)ethan-1-amine",
        arms_in_name=5,          # bis + 1 + bis
        si_amine_mg_mmol=[(26.0, 0.121), (64.6, 0.300), (26.0, 0.121)],
        si_epoxide_mmol_equiv=[(0.724, 6.0), (1.80, 6.0), (0.724, 6.0)],
        note="1,4-disubstituted piperazine carrying a 2-aminoethyl arm on one N "
             "and a 2-((2-aminoethyl)amino)ethyl arm on the other. Five N-H.",
    ),
    "494": dict(
        smiles="NCCOCCN1CCN(CCNCCOCCN)CC1",
        iupac_fragment="the 200 skeleton with an ethoxy (-O-) inserted in each terminal arm",
        arms_in_name=5,
        si_amine_mg_mmol=[(40.0, 0.132), (40.0, 0.132), (40.0, 0.132)],
        si_epoxide_mmol_equiv=[(0.791, 6.0), (0.791, 6.0), (0.791, 6.0)],
        note="Same five-nitrogen skeleton as 200 with two ether oxygens.",
    ),
    "306": dict(
        smiles="CN(CCCN)CCCN",
        iupac_fragment="N1-(3-aminopropyl)-N1-methylpropane-1,3-diamine",
        arms_in_name=4,          # tetrakis
        si_amine_mg_mmol=[(43.6, 0.300), (43.6, 0.300), (43.6, 0.300)],
        si_epoxide_mmol_equiv=[(2.10, 7.0), (2.10, 7.0), (2.10, 7.0)],
        note="methyl-bis(3-aminopropyl)amine. Three nitrogens, FOUR N-H, so four arms.",
    ),
    "114": dict(
        smiles="NCCNCCN",
        iupac_fragment="diethylenetriamine (N1-(2-aminoethyl)ethane-1,2-diamine)",
        arms_in_name=5,          # tetrakis + 1
        si_amine_mg_mmol=[(30.9, 0.300), (30.9, 0.300), (30.9, 0.300)],
        si_epoxide_mmol_equiv=[(2.10, 7.0), (2.10, 7.0), (2.10, 7.0)],
        note="Diethylenetriamine. Five N-H, so five arms: the central secondary "
             "N takes one, which is why the name reads 'tetrakis' PLUS a "
             "(2-hydroxyalkyl)azanediyl.",
    ),
    "490": dict(
        smiles="NCC(OCC)CN1CCN(CCNCC(OCC)CN)CC1",
        iupac_fragment="the 200 skeleton with each terminal 2-aminoethyl replaced by "
                       "3-amino-2-ethoxypropyl",
        arms_in_name=5,
        si_amine_mg_mmol=[(25.0, 0.0750), (25.0, 0.0750), (25.0, 0.0750)],
        si_epoxide_mmol_equiv=[(0.528, 7.0), (0.528, 7.0), (0.528, 7.0)],
        note="Two 2-ethoxypropyl bridges in place of 200's ethylene bridges. "
             "Carries two unassigned stereocentres of its OWN (the ethoxy-bearing "
             "C), independent of the arm stereocentres.",
    ),
}

# ------------------------------------------------------------------ tails
# 2-alkyloxirane + N-H -> N-CH2-CH(OH)-alkyl. The arm is the ring-opened
# product, mapped at the CH2 that bonds to nitrogen.
TAILS = {
    "C10": dict(arm="[CH2:2]C(O)CCCCCCCC",      epoxide="CCCCCCCCC1CO1",
                epoxide_name="2-octyloxirane",  n_carbon_arm=10),
    "C12": dict(arm="[CH2:2]C(O)CCCCCCCCCC",    epoxide="CCCCCCCCCCC1CO1",
                epoxide_name="2-decyloxirane",  n_carbon_arm=12),
    "C14": dict(arm="[CH2:2]C(O)CCCCCCCCCCCC",  epoxide="CCCCCCCCCCCCC1CO1",
                epoxide_name="2-dodecyloxirane", n_carbon_arm=14),
}

# SI-printed epoxide mass (mg) per (tail, amine) -- transcribed from pp. 2-4.
EPOXIDE_MG = {
    ("C10", "200"): 113.0, ("C12", "200"): 332.0, ("C14", "200"): 154.0,
    ("C10", "494"): 124.0, ("C12", "494"): 146.0, ("C14", "494"): 168.0,
    ("C10", "306"): 328.0, ("C12", "306"): 387.0, ("C14", "306"): 446.0,
    ("C10", "114"): 328.0, ("C12", "114"): 387.0, ("C14", "114"): 446.0,
    ("C10", "490"): 82.5,  ("C12", "490"): 97.3,  ("C14", "490"): 112.0,
}
YIELD_PCT = {
    ("C10","200"):46.5, ("C12","200"):58.0, ("C14","200"):50.0,
    ("C10","494"):52.0, ("C12","494"):38.5, ("C14","494"):43.2,
    ("C10","306"):33.1, ("C12","306"):33.8, ("C14","306"):61.0,
    ("C10","114"):14.9, ("C12","114"):6.02, ("C14","114"):51.6,
    ("C10","490"):22.6, ("C12","490"):78.52, ("C14","490"):3.33,
}


def n_NH(mol):
    """Total N-H count. This is what decides the arm count."""
    return sum(a.GetTotalNumHs() for a in mol.GetAtoms() if a.GetSymbol() == "N")


def build(core_smiles, arm_smiles, label):
    """Attach one arm per remaining N-H, through assemble.attach only.

    The target N is chosen structurally (first N still bearing an H, in RDKit
    canonical order), never by a literal atom index -- 'atom index 0 of the
    fragment' is the construct CLAUDE.md section 6 forbids.
    """
    cur = Chem.MolFromSmiles(core_smiles)
    if cur is None:
        raise SystemExit(f"{label}: core SMILES does not parse")
    n_arms = 0
    while True:
        m = Chem.MolFromSmiles(Chem.MolToSmiles(cur))
        target = None
        for a in m.GetAtoms():
            if a.GetSymbol() == "N" and a.GetTotalNumHs() > 0:
                target = a.GetIdx()
                break
        if target is None:
            break
        for a in m.GetAtoms():
            a.SetAtomMapNum(0)
        m.GetAtomWithIdx(target).SetAtomMapNum(1)
        cur = Chem.MolFromSmiles(
            assemble.attach(Chem.MolToSmiles(m), arm_smiles,
                            a_map=1, b_map=2, label=f"{label}:arm{n_arms+1}"))
        n_arms += 1
        if n_arms > 12:
            raise SystemExit(f"{label}: runaway arm loop")
    return Chem.MolToSmiles(cur), n_arms


def main():
    rows, blocks = [], []
    for amine, cd in CORES.items():
        cm = Chem.MolFromSmiles(cd["smiles"])
        if cm is None:
            raise SystemExit(f"amine {amine}: core does not parse")
        blocks.append(dict(
            block_type="amine_core", block_id=amine, SMILES=Chem.MolToSmiles(cm),
            formula=CalcMolFormula(cm), MW=round(MolWt(cm), 3),
            n_NH=n_NH(cm), arms_in_name=cd["arms_in_name"],
            iupac_fragment=cd["iupac_fragment"], note=cd["note"]))
    for tail, td in TAILS.items():
        em = Chem.MolFromSmiles(td["epoxide"])
        blocks.append(dict(
            block_type="epoxide_reagent", block_id=tail, SMILES=Chem.MolToSmiles(em),
            formula=CalcMolFormula(em), MW=round(MolWt(em), 3), n_NH="",
            arms_in_name="", iupac_fragment=td["epoxide_name"],
            note=f"arm after ring opening = {td['arm']}, {td['n_carbon_arm']} carbons"))

    for amine, cd in CORES.items():
        for tail, td in TAILS.items():
            name = f"{tail}-{amine}"
            smi, n_arms = build(cd["smiles"], td["arm"], name)
            m = Chem.MolFromSmiles(smi)
            assemble.check_no_radicals(m, label=name)
            sc = assemble.stereocentres(m)
            rows.append(dict(
                IL_name=name, IL_SMILES=smi,
                amine_core=amine, tail=tail,
                arm_carbons=td["n_carbon_arm"],
                n_arms_built=n_arms,
                n_arms_in_iupac_name=cd["arms_in_name"],
                core_n_NH=n_NH(Chem.MolFromSmiles(cd["smiles"])),
                formula=CalcMolFormula(m),
                MW=round(MolWt(m), 4), ExactMolWt=round(ExactMolWt(m), 4),
                InChIKey=inchi.MolToInchiKey(m),
                InChIKey_skeleton=inchi.MolToInchiKey(m).split("-")[0],
                n_stereocentres=len(sc),
                n_stereocentres_assigned=sum(1 for v in sc.values() if v != "?"),
                epoxide_name=td["epoxide_name"],
                si_yield_pct=YIELD_PCT[(tail, amine)],
                si_amine_mg=cd["si_amine_mg_mmol"][list(TAILS).index(tail)][0],
                si_amine_mmol_printed=cd["si_amine_mg_mmol"][list(TAILS).index(tail)][1],
                si_epoxide_mg=EPOXIDE_MG[(tail, amine)],
                si_epoxide_mmol_printed=cd["si_epoxide_mmol_equiv"][list(TAILS).index(tail)][0],
                si_epoxide_equiv=cd["si_epoxide_mmol_equiv"][list(TAILS).index(tail)][1],
                source_locus="SI 'Ionizable Lipid Synthesis', pp. 2-4, systematic name",
            ))

    with open(os.path.join(HERE, "01_blocks.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(blocks[0]))
        w.writeheader(); w.writerows(blocks)
    with open(os.path.join(HERE, "03_enumerated.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)
    print(f"blocks: {len(blocks)}   lipids: {len(rows)}")
    for r in rows:
        print(f"  {r['IL_name']:9s} {r['formula']:14s} MW {r['MW']:9.3f}  "
              f"arms built {r['n_arms_built']} (name says {r['n_arms_in_iupac_name']})  "
              f"stereo {r['n_stereocentres_assigned']}/{r['n_stereocentres']}  {r['InChIKey']}")


if __name__ == "__main__":
    main()
