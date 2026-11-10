#!/usr/bin/env python3
"""Stage 03 -- enumerate the 36 selected cyclic-diamine ionizable lipids.

    python3 -I stage03_enumerate.py

THE LIBRARY. 6 functional heads x 3 ester bridges x 2 tails = 36, exactly the
count the paper states for Figure 3a ("Schematic diagram of the 36 ionizable
lipids") and the count of cells in Fig S2c, Fig S2e, Fig 4b and Table S3.

THE ASSEMBLY ORDER IS THE SYNTHESIS ORDER (SI 1.2.5):

  1. anhydride + branched alcohol            -> EB-NT-M acid   (mono-ester)
  2. EB-NT-M acid x2 + the core's 1,2-diol   -> compound 4/6
  3. TBAF                                    -> FH-3 product (free 1-degree OH)
  4. FH acid + that primary OH               -> FH-N-EB-NT-M

Every bond is made with `assemble.attach`, which bonds two atom-MAPPED atoms
and REMOVES NOTHING. There is no RWMol.RemoveAtom anywhere in this file, so the
F-013 mechanism -- deleting an atom out of a stereocentre's neighbour list and
inverting its chiral tag on 537 rows -- cannot occur. `attach` runs
check_no_radicals and check_stereo_preserved on every build.

`esterify` is deliberately NOT used: it requires exactly one free COOH and one
free OH, and it would have to be handed a diacid (two COOH) and a triol (three
OH). Which group reacts would then be decided by a substructure match order
rather than by the chemistry. Mapping the atoms states it instead.

STEREOCHEMISTRY IS LEFT UNSPECIFIED, ON PURPOSE. Three sp3 centres in every
member are unassigned in the paper: the glycerol CH (built from racemic
3-iodopropane-1,2-diol, SI 1.2.4), and the branch carbon of each tail
(2-hexyl-1-decanol and 2-octyl-1-dodecanol are sold as racemates). FH-1 adds a
fourth. The paper assigns none of them, draws none of them with a wedge, and
reports no optical rotation. An unspecified centre is the honest encoding; a
specified one would be a fabrication.

WHAT THIS STAGE DOES NOT CLAIM. The paper's screen is 213 LNPs built from 14
functional heads and 14 EB-T chains (main text, ML section). Those 14x14 blocks
are drawn NOWHERE: Fig S1 shows only the generic scaffold with R1 and R2, and
Fig 3a draws only the six heads, three bridges and two tails used for the
36-member follow-up. So 36 is all that can be enumerated, and 177 of the 213
LNPs have no recoverable structure. See 08_contradictions.csv and DECISIONS.md
D-1.
"""
import csv
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import assemble                                             # noqa: E402

HEADS = ["FH-1", "FH-2", "FH-3", "FH-4", "FH-5", "FH-6"]
BRIDGES = ["EB-5", "EB-6", "EB-7"]
TAILS = ["T-1", "T-2"]

# The core, with its three hydroxyl oxygens atom-mapped so that WHICH oxygen
# receives WHICH acyl group is stated here, not inferred from a match order.
#   map 2 = the hydroxyethyl arm's primary OH  -> the FH acid (or nothing)
#   map 3 = the glycerol secondary OH          -> an EB-T acid
#   map 4 = the glycerol primary OH            -> an EB-T acid
CORE_MAPPED = "[OH:2]CCN1CCN(CC1)CC([OH:3])C[OH:4]"


def blocks():
    with open(os.path.join(HERE, "01_blocks.csv")) as fh:
        return {r["block_id"]: r for r in csv.DictReader(fh)}


def ebt_acyl(bl, eb, t):
    """EB-N + T-M -> the atom-mapped acyl fragment of the EB-T mono-ester.

    Both halves are strings stage01 BUILT from measured counts, so joining
    them at the EB fragment's terminal ester oxygen is a statement about
    atoms, not a rewrite of a canonical SMILES.

    THE FIRST VERSION OF THIS FUNCTION WAS WRONG AND NOTHING COMPOSITIONAL
    CAUGHT IT. It built the fragment by stripping 'C(=O)O' off the END of the
    canonical EB-T acid and prefixing '[C:1](=O)'. In SMILES the atom after a
    prefix bonds to the FIRST atom of the remainder, and canonicalisation had
    moved the tail to the front -- so the glycerol ester landed on the far end
    of the tail and the diacid became a butyryl stub. Formula unchanged, both
    HRMS anchors still matching at 0.3 and 0.5 ppm.

    The guards below are what a string join has to earn:
      - the EB fragment must end in the ester oxygen it is being joined at;
      - neither string may carry a ring-closure digit, because fusing two
        fragments that reuse ring number 1 is how piperazine once ended up
        bonded into cholesterol (CLAUDE.md section 6).
    """
    acyl, alkyl = bl[eb]["fragment"], bl[t]["fragment"]
    if not acyl.endswith("C(=O)O"):
        raise SystemExit(
            f"ebt_acyl: {eb}'s fragment {acyl!r} does not end in a free "
            "carboxylic acid; refusing to guess where the tail attaches.")
    bare = [re.sub(r"\[[^\]]*\]", "", x) for x in (acyl, alkyl)]
    if any(ch.isdigit() for x in bare for ch in x):
        raise SystemExit(
            f"ebt_acyl: a ring-closure digit in {acyl!r} or {alkyl!r}. "
            "Joining two fragments that both use ring number 1 fuses them at "
            "the wrong atoms.")
    return acyl + alkyl


def member(bl, fh, eb, t):
    acyl = ebt_acyl(bl, eb, t)
    # the glycerol 1,2-diol, both oxygens, in synthesis order
    mol = assemble.attach(CORE_MAPPED, acyl, a_map=3, b_map=1,
                          label=f"core-2-O-{eb}{t}")
    mol = assemble.attach(mol, acyl, a_map=4, b_map=1,
                          label=f"core-1-O-{eb}{t}")
    if fh == "FH-3":
        # "None FH (1 degree alcohol)": TBAF leaves the arm as a free primary
        # alcohol and nothing is coupled to it. The map must still be cleared.
        return assemble.canonical(mol.replace("[OH:2]", "O"))
    return assemble.canonical(
        assemble.attach(mol, bl[fh]["fragment"], a_map=2, b_map=1,
                        label=f"{fh}-O-core"))


# --------------------------------------------------------------------------
# CONNECTIVITY ASSERTIONS -- the ones mass cannot make
# --------------------------------------------------------------------------
def connectivity(smiles, fh, eb, t, bl):
    """Counts that change when a bond moves but a formula does not.

    `tertiary_C` is the paper's own FP3. Figure 2b prints FP3 = 3 for exactly
    36 of the 213 LNPs and FP3 = 2 for the other 177, and 36 is the number of
    FH-1 LNPs -- FH-1 is the only functional head with a branch point. So the
    required value here is read off a printed statistic, not chosen.

    `esters` is 5 for a member with a functional head (two tail esters, two
    glycerol esters, one head ester) and 4 for FH-3, whose arm stays a free
    primary alcohol.
    """
    from rdkit import Chem
    m = Chem.MolFromSmiles(smiles)
    q = {
        "tertiary_C": len(m.GetSubstructMatches(
            Chem.MolFromSmarts("[CX4H1]([#6])([#6])[#6]"))),
        "esters": len(m.GetSubstructMatches(
            Chem.MolFromSmarts("[CX3](=O)[OX2H0][#6]"))),
        "free_acids": len(m.GetSubstructMatches(
            Chem.MolFromSmarts("[CX3](=O)[OX2H1]"))),
        "free_OH": len(m.GetSubstructMatches(
            Chem.MolFromSmarts("[OX2H1][CX4]"))),
        "phenols": len(m.GetSubstructMatches(
            Chem.MolFromSmarts("[OX2H1][c]"))),
        "ring_N": sum(1 for a in m.GetAtoms()
                      if a.GetSymbol() == "N" and a.IsInRing()),
        "longest_acyl_chain_from_glycerol_ester": _acyl_span(m),
    }
    want = {
        "tertiary_C": 2 + (1 if fh == "FH-1" else 0),
        "esters": 4 if fh == "FH-3" else 5,
        "free_acids": 0,
        "free_OH": 1 if fh == "FH-3" else 0,
        "phenols": {"FH-4": 2, "FH-5": 2, "FH-6": 2}.get(fh, 0),
        "ring_N": 2,
        "longest_acyl_chain_from_glycerol_ester": int(eb.split("-")[1]),
    }
    bad = {k: (q[k], want[k]) for k in want if q[k] != want[k]}
    if bad:
        raise SystemExit(
            f"stage03: {fh}/{eb}/{t} FAILS a connectivity assertion "
            f"(measured, required): {bad}\n{smiles}\n"
            "These are exactly the counts a molecular formula cannot see. "
            "tertiary_C is the paper's own FP3, which Figure 2b prints as 3 "
            "for 36 of 213 LNPs and 2 for the other 177.")
    return q


def _acyl_span(m):
    """Carbons in the acyl group hanging off a glycerol ester oxygen.

    Walks from a glycerol ester's carbonyl carbon along carbons only, never
    crossing an oxygen, and returns the size of that component. For a correct
    member this is the ester bridge's own backbone length: carbonyl, the CH2
    run, and the second carbonyl, i.e. 5 / 6 / 7. When the first version of
    `ebt_acyl` put the glycerol ester on the far end of the tail, this came
    back as 18 instead of 5.
    """
    from rdkit import Chem
    patt = Chem.MolFromSmarts("[CX4;!R][OX2H0][CX3](=O)")
    hits = m.GetSubstructMatches(patt)
    if not hits:
        return -1
    spans = []
    for _c_gly, _o, c_carbonyl, _o2 in hits:
        seen, stack = set(), [c_carbonyl]
        while stack:
            i = stack.pop()
            if i in seen:
                continue
            a = m.GetAtomWithIdx(i)
            if a.GetSymbol() != "C":
                continue
            seen.add(i)
            for nb in a.GetNeighbors():
                if nb.GetSymbol() == "C" and nb.GetIdx() not in seen:
                    stack.append(nb.GetIdx())
        spans.append(len(seen))
    return min(spans)


def main():
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Descriptors, rdMolDescriptors
    RDLogger.DisableLog("rdApp.*")

    bl = blocks()
    rows, seen = [], {}
    for fh in HEADS:
        for t in TAILS:
            for eb in BRIDGES:
                smi = member(bl, fh, eb, t)
                conn = connectivity(smi, fh, eb, t, bl)
                m = Chem.MolFromSmiles(smi)
                if m is None:
                    raise SystemExit(f"stage03: {fh}/{eb}/{t} does not parse.")
                if "[" in smi and ":" in smi:
                    raise SystemExit(
                        f"stage03: {fh}/{eb}/{t} still carries an atom map: "
                        f"{smi}. An unconsumed map means a hydroxyl that "
                        "should have been acylated was not.")
                name = f"{fh}-{eb}{t}"
                if smi in seen:
                    raise SystemExit(
                        f"stage03: {name} is the same structure as "
                        f"{seen[smi]}. 36 distinct members were expected.")
                seen[smi] = name
                unspec = Chem.FindMolChiralCenters(
                    m, includeUnassigned=True, useLegacyImplementation=False)
                rows.append({
                    "IL_name": name,
                    "FH": fh, "EB": eb, "T": t,
                    "IL_SMILES": smi,
                    "formula": rdMolDescriptors.CalcMolFormula(m),
                    "exact_mass": round(Descriptors.ExactMolWt(m), 4),
                    "MolWt": round(Descriptors.MolWt(m), 3),
                    "n_unassigned_stereocentres":
                        sum(1 for _, c in unspec if c == "?"),
                    "tertiary_C_FP3": conn["tertiary_C"],
                    "n_esters": conn["esters"],
                    "n_phenols": conn["phenols"],
                    "acyl_span_from_glycerol_ester":
                        conn["longest_acyl_chain_from_glycerol_ester"],
                    "head_confidence": bl[fh]["confidence"],
                    "bridge_confidence": bl[eb]["confidence"],
                    "tail_confidence": bl[t]["confidence"],
                })

    # the member confidence is the WEAKEST of its own blocks -- never a literal
    order = {"high": 2, "medium": 1, "low": 0}
    inv = {v: k for k, v in order.items()}
    for r in rows:
        r["confidence"] = inv[min(order[r["head_confidence"]],
                                  order[r["bridge_confidence"]],
                                  order[r["tail_confidence"]])]

    with open(os.path.join(HERE, "03_enumerated.csv"), "w", newline="") as fh_:
        w = csv.DictWriter(fh_, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"03_enumerated.csv: {len(rows)} members, "
          f"{len(set(r['IL_SMILES'] for r in rows))} distinct SMILES")
    for r in rows:
        if r["IL_name"].startswith("FH-5-EB-5"):
            print(f"  {r['IL_name']:>18}  {r['formula']:>14}  "
                  f"exact {r['exact_mass']}")


if __name__ == "__main__":
    main()
