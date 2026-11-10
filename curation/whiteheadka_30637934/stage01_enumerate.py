#!/usr/bin/env python3
"""306-family lipidoid enumeration for Hajj et al. 2019 (PMID 30637934).

CHEMISTRY, read from the SI Methods (p.2), not inferred:

    "Lipidoids were synthesized as previously described [16] by reacting the
     amines 2[4-2(2-aminoethyl)amino)ethylpiperazine-1-YL)ethan-1-amine (200,
     Enamine), Tris[2-(methylamino)ethyl]amine (304, Sigma Aldrich), and
     3,3-Diamino-N-methyldipropylamine (306, Sigma Aldrich) with alkyl-acrylate
     tails of varying lengths (Sartomer, TCI Chemicals) ... the amine was
     combined with acrylate at a molar ratio of 1:4 ... In vitro experiments
     were conducted with crude lipidoids. For in vivo experiments, the
     lipidoids were purified ... to isolate the fully substituted lipidoid
     product."

Only the 306 family carries a value anywhere in the SI, so only the 306 family
is enumerated here. Amine 306 = 3,3'-diamino-N-methyldipropylamine, two primary
amines = FOUR N-H, and 1:4 amine:acrylate is exactly stoichiometric for the
tetra-substituted product. Every arm is a Michael adduct of a primary/secondary
amine onto an alkyl acrylate: N-CH2CH2-C(=O)O-R.

HOW THIS IS BUILT. Through curation/lib/assemble only:
  * the ESTER by `esterify`, so the acid loses its OH and the alcohol keeps its
    oxygen (the F-013 asymmetry);
  * every N-C bond by `attach` on ATOM-MAPPED attachment atoms, never on an
    atom index -- a sibling built acyl fragments by string surgery on canonical
    SMILES and mislocated every ester with the exact mass identical to 8 dp.
`tag()` below only SETS a map number on a uniquely-matched atom; it forms and
breaks no bonds, and it asserts the match is unique so a silent second site is
impossible.

NAME TOKEN -> CHAIN LENGTH. "O<n>" is an n-carbon n-alkyl ester; "Oi10" is
ISODECYL. That reading is not assumed from the token (CLAUDE.md: a name token
is not a chain length) -- it is pinned by two prior, independent curations
already in the corpus, which store these exact names with these exact
structures: NC_2024 (A306_O10..O14, A306_Oi10) and KW_2014 (A306_C10..C14).
stage03_corpus_overlap.py does that comparison and it is the paper's only
external structure anchor, since the SI prints no HRMS, no formula and no MW.
"""
import json
import os
import sys

from rdkit import Chem
from rdkit import RDLogger

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "lib"))
import assemble                                                  # noqa: E402

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))

# 3,3'-diamino-N-methyldipropylamine (306), one primary N pre-tagged.
AMINE_306 = "[NH2:2]CCCN(C)CCCN"
AMINE_306_NAME = "3,3'-diamino-N-methyldipropylamine (306)"

# Propanoic acid with map 1 on C3 -- the carbon that bonds to nitrogen. It is
# written as a VALENCE-COMPLETE methyl on purpose: '[CH2:1]CC(=O)O' is a
# radical the moment it is parsed, and `attach` hands the hydrogen count back
# to RDKit on the bonded atom, so CH3 becomes CH2 when the N-C bond forms.
ACID = "[CH3:1]CC(=O)O"

# The alcohol half of each alkyl acrylate, written from the reagent's chemical
# NAME. Numbers are carbons in the alkyl group, matching the "O<n>" token.
TAILS = {
    "O6":   ("hexan-1-ol",            "CCCCCCO",          6, "n-hexyl"),
    "O8":   ("octan-1-ol",            "CCCCCCCCO",        8, "n-octyl"),
    "O9":   ("nonan-1-ol",            "CCCCCCCCCO",       9, "n-nonyl"),
    "Oi10": ("8-methylnonan-1-ol",    "CC(C)CCCCCCCO",   10, "isodecyl (8-methylnonyl)"),
    "O10":  ("decan-1-ol",            "CCCCCCCCCCO",     10, "n-decyl"),
    "O11":  ("undecan-1-ol",          "CCCCCCCCCCCO",    11, "n-undecyl"),
    "O12":  ("dodecan-1-ol",          "CCCCCCCCCCCCO",   12, "n-dodecyl"),
    "O13":  ("tridecan-1-ol",         "CCCCCCCCCCCCCO",  13, "n-tridecyl"),
    "O14":  ("tetradecan-1-ol",       "CCCCCCCCCCCCCCO", 14, "n-tetradecyl"),
    "O16":  ("hexadecan-1-ol",        "C" * 16 + "O",    16, "n-hexadecyl"),
    "O18":  ("octadecan-1-ol",        "C" * 18 + "O",    18, "n-octadecyl"),
}

# Order as the paper plots and tabulates them (Table S1, Figure S6A).
ORDER = ["O6", "O8", "O9", "Oi10", "O10", "O11", "O12", "O13", "O14",
         "O16", "O18"]


def tag(smiles, smarts, mapno):
    """Set map `mapno` on the one atom matching `smarts`. Raises if not unique.

    This only labels an atom. It adds no bond, removes no atom and touches no
    hydrogen count, so it cannot be the F-013 mechanism. The uniqueness
    assertion is the point: `assemble` refuses an ambiguous attachment atom and
    this keeps that refusal reachable one step earlier.
    """
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        raise ValueError(f"tag: unparseable {smiles!r}")
    hits = m.GetSubstructMatches(Chem.MolFromSmarts(smarts))
    if len(hits) != 1:
        raise AssertionError(
            f"tag: {smarts!r} matched {len(hits)} atoms in {smiles!r}, need 1. "
            "A second match means the next arm could go to either nitrogen, "
            "which is a regiochemistry decision, not a default.")
    m.GetAtomWithIdx(hits[0][0]).SetAtomMapNum(mapno)
    return Chem.MolToSmiles(m)


def build(tail_token):
    """Assemble one tetra-substituted 306 lipidoid. Returns (smiles, trace)."""
    alc_name, alc_smiles, ncarb, alkyl = TAILS[tail_token]
    arm = assemble.esterify(ACID, alc_smiles, label=f"{tail_token} arm ester")
    trace = {"arm_ester": arm, "alcohol": alc_smiles,
             "alcohol_name": alc_name, "alkyl": alkyl, "tail_carbons": ncarb}

    # Arm 1 onto the pre-tagged primary nitrogen.
    p = assemble.attach(arm, AMINE_306, a_map=1, b_map=2,
                        label=f"306{tail_token} arm1")
    steps = [p]
    # Arm 2 onto the SAME nitrogen, now the molecule's only secondary amine.
    p = assemble.attach(arm, tag(p, "[NX3;H1]", 2), 1, 2,
                        label=f"306{tail_token} arm2")
    steps.append(p)
    # Arm 3 onto the remaining primary amine.
    p = assemble.attach(arm, tag(p, "[NX3;H2]", 2), 1, 2,
                        label=f"306{tail_token} arm3")
    steps.append(p)
    # Arm 4 onto that nitrogen, again the only secondary amine left.
    p = assemble.attach(arm, tag(p, "[NX3;H1]", 2), 1, 2,
                        label=f"306{tail_token} arm4")
    steps.append(p)
    trace["intermediates"] = steps[:-1]
    return p, trace


def main():
    from rdkit.Chem import Descriptors, rdMolDescriptors

    blocks, rows = [], []
    for tok in ORDER:
        alc_name, alc_smiles, ncarb, alkyl = TAILS[tok]
        blocks.append({"block_type": "tail", "block_name": tok,
                       "reagent": f"{alkyl} acrylate",
                       "alcohol_name": alc_name, "alcohol_SMILES": alc_smiles,
                       "tail_carbons": ncarb,
                       "arm_after_michael_addition":
                           f"N-CH2CH2-C(=O)O-{alkyl}"})
    blocks.append({"block_type": "amine_core", "block_name": "306",
                   "reagent": AMINE_306_NAME,
                   "alcohol_name": "", "alcohol_SMILES": "",
                   "tail_carbons": "",
                   "arm_after_michael_addition":
                       "4 N-H, tetra-substituted at 1:4 amine:acrylate"})

    for tok in ORDER:
        smi, trace = build(tok)
        m = Chem.MolFromSmiles(smi)
        if m is None:
            raise SystemExit(f"{tok}: product does not parse")
        assemble.check_no_radicals(m, label=f"306{tok} final")
        n_arms = len(m.GetSubstructMatches(
            Chem.MolFromSmarts("[NX3]CCC(=O)O[CX4]")))
        if n_arms != 4:
            raise SystemExit(f"306{tok}: {n_arms} Michael arms, expected 4")
        rows.append({
            "IL_name": f"A306_{tok}",
            "paper_name": f"306{tok}",
            "IL_SMILES": assemble.canonical(smi),
            "formula": rdMolDescriptors.CalcMolFormula(m),
            "exact_mass": round(Descriptors.ExactMolWt(m), 6),
            "avg_mw": round(Descriptors.MolWt(m), 3),
            "n_michael_arms": n_arms,
            "amine_core": "306",
            "tail_token": tok,
            "tail_carbons": trace["tail_carbons"],
            "tail_alkyl": trace["alkyl"],
            "arm_ester_SMILES": trace["arm_ester"],
            "substitution_state_stored": "tetra (fully substituted)",
            "material_measured_in_vitro": "crude (SI Methods p.2)",
            "material_measured_in_vivo": "purified, fully substituted (SI Methods p.2)",
        })

    import csv
    with open(os.path.join(HERE, "01_blocks.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(blocks[0].keys()))
        w.writeheader()
        w.writerows(blocks)
    with open(os.path.join(HERE, "03_enumerated.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    # The authors' own unprinted identity: 306O10 and 306Oi10 "have ten-carbon
    # tails and identical molecular weights" (abstract). Nothing in the SI
    # prints a mass, so this is the only mass statement the paper makes.
    o10 = next(r for r in rows if r["tail_token"] == "O10")
    oi10 = next(r for r in rows if r["tail_token"] == "Oi10")
    identity = {
        "claim": "abstract: 306O10 and 306Oi10 have ten-carbon tails and identical molecular weights",
        "O10_formula": o10["formula"], "Oi10_formula": oi10["formula"],
        "formulae_identical": o10["formula"] == oi10["formula"],
        "O10_exact_mass": o10["exact_mass"], "Oi10_exact_mass": oi10["exact_mass"],
        "exact_mass_delta_Da": round(abs(o10["exact_mass"] - oi10["exact_mass"]), 9),
        "O10_tail_carbons": o10["tail_carbons"],
        "Oi10_tail_carbons": oi10["tail_carbons"],
        "both_ten_carbon": o10["tail_carbons"] == oi10["tail_carbons"] == 10,
        "isomers_not_identical": o10["IL_SMILES"] != oi10["IL_SMILES"],
    }
    with open(os.path.join(HERE, "01_mass_identity.json"), "w") as fh:
        json.dump(identity, fh, indent=1)
    print(json.dumps(identity, indent=1))
    print(f"\n{len(rows)} structures built, {len(blocks)} blocks.")
    for r in rows:
        print(f"  {r['IL_name']:12s} {r['formula']:14s} "
              f"{r['exact_mass']:12.4f}  arms={r['n_michael_arms']}")


if __name__ == "__main__":
    main()
