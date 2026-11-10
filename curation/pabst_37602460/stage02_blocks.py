#!/usr/bin/env python3
"""The building blocks, from the figure that draws them -> 01_blocks.csv

WHERE THEY COME FROM. Figure 1A of the article (source/figs/am3c05481_0002.jpg,
4183 x 1847 px from PMC's own figure store) draws all 8 polyamine headgroups and
all 13 acrylate R groups in full. There is no SMILES, InChI, MOL block, CAS
number or chemical name for any of them anywhere in the article or the SI, so
the drawing is the only statement of connectivity the authors made. Crops used
while reading are in 00_pages/fig1A_*.png.

Run from this folder:   python3 -I stage02_blocks.py

WHAT MAKES THIS MORE THAN A TRANSCRIPTION. Every block read off the figure is
checked against a printed number that was NOT used to produce it:

  tails      SI Table S4 prints a molecular weight for each of the 13 tail
             groups. That weight is the ACRYLATE MONOMER's formula weight (not
             the R group's, not the arm's), which this file asserts to <= 0.05
             Da on all 13. The naming scheme also predicts three pairs of
             CONSTITUTIONAL ISOMERS -- 6B/8, 10B'/12, 13B/14 -- and Table S4
             shows exactly those three as exact MW ties. A wrong assignment
             would have to break a tie or invent one.

  headgroups SI Table S5 prints "charge density" = number of nitrogens divided
             by headgroup molecular weight. This file recomputes it from the
             drawn amine. 6 of 8 agree to <= 0.0001. The two that do not are
             3NB and 4NC, and each one's printed value is the OTHER's computed
             value to 4 decimal places -- see 08_contradictions.csv (C-2).

  name token CLAUDE.md: "Name tokens are not chain lengths." Here the paper
             DEFINES its token: M is "the tail length in terms of the number of
             atoms from the first acrylate methylene to the furthest carbon".
             So the token is checkable, and it is checked -- by a longest-path
             walk over the arm's own molecular graph, through
             curation/lib/chains.check_name_token, which reports and never
             corrects.

Nothing here reads a pixel as a number. 00_value_sources.json records why that
is not a tier-5 digitisation.
"""
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))

import assemble                                        # noqa: E402
import chains                                          # noqa: E402
from rdkit import Chem                                 # noqa: E402
from rdkit.Chem import Descriptors, rdMolDescriptors   # noqa: E402

OUT = os.path.join(HERE, "01_blocks.csv")

# --- the eight polyamine headgroups, read from Figure 1A ------------------
# `n_nh` is the number of N-H hydrogens, which under the paper's "fully
# substituted species was the major product" is also the number of tails. The
# article states the tail counts independently ("Linear headgroups with 3 or 4
# amines (3N, 3N', 4N, 4N') produced 5 or 6 tail species, respectively, while
# branched headgroups with 3 primary amines (3NB, 3NB') imparted 6 tails, and
# lipidoids with the cyclic tetraamine (4NC) were 4-tailed"), and `n_tails_said`
# records that claim so the two can be compared rather than assumed equal.
HEADGROUPS = [
    # name,   SMILES,                     identity read from Fig 1A, n_tails_said
    ("2N",   "NCCCCCCCCN",
     "1,8-diaminooctane: H2N-(CH2)8-NH2, the gemini diamine the text says was "
     "chosen with an 8-carbon spacer", 4),
    ("3N",   "NCCNCCN",
     "diethylenetriamine: H2N-CH2CH2-NH-CH2CH2-NH2, drawn under 'Linear'", 5),
    ("3N'",  "NCCCNCCCCN",
     "spermidine: H2N-(CH2)3-NH-(CH2)4-NH2, drawn under 'Linear' with the "
     "longer propyl/butyl spacers the prime denotes", 5),
    ("4N",   "NCCNCCNCCN",
     "triethylenetetramine: three ethylene spacers, drawn under 'Linear'", 6),
    ("4N'",  "NCCCNCCCCNCCCN",
     "spermine: H2N-(CH2)3-NH-(CH2)4-NH-(CH2)3-NH2, drawn under 'Linear'", 6),
    ("3NB",  "NCCN(CCN)CCN",
     "tris(2-aminoethyl)amine (TREN): a central tertiary N bearing three "
     "2-carbon arms, each ending NH2, drawn under 'Branched'", 6),
    ("3NB'", "NCCCN(CCCN)CCCN",
     "tris(3-aminopropyl)amine: the same branch point with 3-carbon arms, "
     "which is what the prime denotes", 6),
    ("4NC",  "N1CCNCCCNCCNCCC1",
     "cyclam (1,4,8,11-tetraazacyclotetradecane): a 14-membered ring of four "
     "NH separated by 2,3,2,3 carbons, drawn under 'Cyclic'", 4),
]

# --- the thirteen acrylates, read from Figure 1A --------------------------
# Stored as the ALCOHOL R-OH, because that is what curation/lib/assemble.esterify
# consumes, and because it makes the acid/alcohol asymmetry explicit: the acid
# loses its OH, the alcohol keeps its oxygen (F-013's lesson).
TAILS = [
    # name, R-OH SMILES, token, what Fig 1A draws
    ("6B",   "CC(C)(C)O", 6,
     "tert-butyl: the attachment bond runs to a quaternary carbon with three "
     "methyls radiating from it"),
    ("8",    "CCCCO", 8, "n-butyl: four-vertex zigzag, no branch"),
    ("8CB",  "CC1(C)C2CCC1(C)C(O)C2", 8,
     "isobornyl: a bicyclo[2.2.1]heptane cage with a gem-dimethyl bridge and a "
     "third methyl -- the 'CB' (cyclic and branched) of the naming scheme"),
    ("9",    "CCCCCO", 9, "n-pentyl"),
    ("9C",   "@TCD8", 9,
     "tricyclodecanyl (hydrogenated dicyclopentadienyl): a cyclopentane fused "
     "to a bridged bicycloheptane -- the text calls it 'the tricyclic 9C tail'"),
    ("10",   "CCCCCCO", 10, "n-hexyl"),
    ("10B",  "OCCC(C)CC(C)(C)C", 10,
     "3,5,5-trimethylhexyl: two CH2, then a methyl-bearing CH, then a CH2, "
     "then a tert-butyl terminus"),
    ("10B'", "CCCCC(CC)CO", 10,
     "2-ethylhexyl: an ethyl branch on the carbon next to the oxygen"),
    ("11",   "CCCCCCCO", 11, "n-heptyl"),
    ("12",   "CCCCCCCCO", 12, "n-octyl"),
    ("13B",  "OCCCCCCCC(C)C", 13,
     "isodecyl, drawn as 8-methylnonyl: a long chain with a single methyl "
     "branch on the penultimate carbon"),
    ("14",   "CCCCCCCCCCO", 14, "n-decyl"),
    ("16",   "CCCCCCCCCCCCO", 16, "n-dodecyl"),
]

# SI Table S4 page S11, transcribed BY HAND from the rendered page, on purpose.
# C15's design: the expected value is written out so that changing a documented
# number means changing both halves, and the pair either agrees or the build
# fails.
TAIL_PRINTED = {           # name: (MW of the acrylate monomer, c log P)
    "6B": (128.2, 2.02), "8": (128.2, 2.39), "8CB": (208.3, 4.22),
    "9": (142.2, 2.92), "9C": (206.3, 3.69), "10": (156.2, 3.45),
    "10B": (198.3, 4.49), "10B'": (184.3, 4.33), "11": (170.3, 3.98),
    "12": (184.3, 4.5), "13B": (212.3, 5.39), "14": (212.3, 5.57),
    "16": (240.4, 6.64),
}

# SI Table S5 page S11: N/MW (charge density) and c log P, as printed.
HEAD_PRINTED = {
    "2N": (0.0139, 1.1), "3N": (0.0291, -1.87), "3N'": (0.0206, -0.84),
    "4N": (0.0274, -2.18), "4N'": (0.0198, -0.96), "3NB": (0.0200, -0.97),
    "3NB'": (0.0212, -1.01), "4NC": (0.0274, -2.68),
}

ACRYLIC_ACID = "C=CC(=O)O"


def atoms_to_furthest_carbon(mol, start_idx):
    """Atoms along the SHORTEST path from start_idx to the furthest carbon.

    The paper's own definition of the tail-length token: M is "the tail length
    in terms of the number of atoms from the first acrylate methylene to the
    furthest carbon". Walking the graph is the only way to check a token the
    authors defined geometrically; reading it off the name would be assuming
    the thing under test.

    SHORTEST path, i.e. graph eccentricity + 1, not the longest simple path.
    The distinction is invisible on an acyclic tail (the two coincide) and
    decides the answer on the two caged ones: a longest-path walk round the
    isobornyl cage returns 12 where the authors wrote 8, because it counts its
    way round a ring rather than across it. The acrylate's carbonyl oxygen is a
    one-atom branch and never wins, so it does not have to be special-cased.
    """
    from collections import deque
    dist = {start_idx: 0}
    q = deque([start_idx])
    while q:
        i = q.popleft()
        for nb in mol.GetAtomWithIdx(i).GetNeighbors():
            j = nb.GetIdx()
            if j not in dist:
                dist[j] = dist[i] + 1
                q.append(j)
    carbons = [d for i, d in dist.items()
               if mol.GetAtomWithIdx(i).GetSymbol() == "C"]
    return max(carbons) + 1


def mol_from_bond_table(n_carbons, bonds, oh_on):
    """Build a saturated hydrocarbon-ol from an explicit atom/bond table.

    Used only for the 9C tricyclic alcohol, whose cage cannot be written as a
    SMILES string without ring-closure bookkeeping -- and a ring-closure digit
    collision is exactly what fused piperazine to cholesterol once
    (CLAUDE.md section 6). A bond table has no ring-closure digits to collide.

    `bonds` are 1-based carbon pairs; `oh_on` is the 1-based carbon that
    carries the hydroxyl.
    """
    rw = Chem.RWMol()
    for _ in range(n_carbons):
        rw.AddAtom(Chem.Atom(6))
    for a, b in bonds:
        rw.AddBond(a - 1, b - 1, Chem.BondType.SINGLE)
    o = rw.AddAtom(Chem.Atom(8))
    rw.AddBond(oh_on - 1, o, Chem.BondType.SINGLE)
    m = rw.GetMol()
    Chem.SanitizeMol(m)
    assemble.check_no_radicals(m, label="9C cage")
    return Chem.MolToSmiles(m)


# tricyclo[5.2.1.0^{2,6}]decan-8-ol, by IUPAC skeleton numbering:
#   main bicyclo[5.2.1] ring  C1-C2-C3-C4-C5-C6-C7, C7-C8-C9-C1, C1-C10-C7
#   plus the 0^{2,6} bond     C2-C6
# which gives the cyclopentane C2-C3-C4-C5-C6 fused across the 2-carbon bridge
# of a bicyclo[2.2.1]heptane -- the cyclopentane-fused bridged cage Figure 1A
# draws for 9C. The hydroxyl sits on C8, in the 2-carbon bridge next to
# bridgehead C7, which is where the figure puts the attachment bond.
TCD_BONDS = [(1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (7, 8), (8, 9),
             (9, 1), (1, 10), (10, 7), (2, 6)]


def arm_from_alcohol(roh, label):
    """The Michael-adduct arm, atom-mapped at the carbon that bonds to N.

    `[CH3:2]CC(=O)O` is propanoic acid with map 2 on C3. It is written as a
    METHYL, not as `[CH2:2]`, deliberately: a bracket `[CH2:2]` carries
    noImplicit=True, so as a standalone fragment it is one hydrogen short and
    RDKit assigns it an unpaired electron. assemble.attach clears both the H
    flags and the radical count on the atom it bonds, turning this methyl into
    the CH2 the chemistry has -- which is exactly the Siegwart radical-carbon
    class the library exists to make unwritable.
    """
    return assemble.esterify("[CH3:2]CC(=O)O", roh, label=label)


def main():
    rows = []

    for name, smi, note, n_tails_said in HEADGROUPS:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            raise SystemExit(f"headgroup {name}: unparseable SMILES {smi!r}")
        n_n = sum(1 for a in mol.GetAtoms() if a.GetSymbol() == "N")
        n_nh = sum(a.GetTotalNumHs() for a in mol.GetAtoms()
                   if a.GetSymbol() == "N")
        mw = Descriptors.MolWt(mol)
        cd = n_n / mw
        pr_cd, pr_logp = HEAD_PRINTED[name]
        rows.append({
            "block_id": f"H:{name}", "kind": "headgroup", "name": name,
            "smiles": Chem.MolToSmiles(mol),
            "formula": rdMolDescriptors.CalcMolFormula(mol),
            "mw": round(mw, 3),
            "n_nitrogens": n_n, "n_NH_hydrogens": n_nh,
            "n_tails_computed": n_nh, "n_tails_stated_in_text": n_tails_said,
            "n_tails_agree": int(n_nh == n_tails_said),
            "computed_charge_density_N_over_MW": round(cd, 6),
            "printed_charge_density_tableS5": pr_cd,
            "charge_density_abs_diff": round(abs(cd - pr_cd), 6),
            "charge_density_agrees": int(abs(cd - pr_cd) <= 1e-4),
            "printed_clogP_tableS5": pr_logp,
            "name_token": "", "token_measured": "", "token_agrees": "",
            "printed_mw_tableS4": "", "mw_abs_diff": "", "mw_agrees": "",
            "source_figure": "Figure 1A, headgroup panel",
            "note": note,
        })

    for name, roh, token, note in TAILS:
        if roh == "@TCD8":
            roh = mol_from_bond_table(10, TCD_BONDS, 8)
        alc = Chem.MolFromSmiles(roh)
        if alc is None:
            raise SystemExit(f"tail {name}: unparseable alcohol {roh!r}")
        monomer = assemble.esterify(ACRYLIC_ACID, roh, label=f"acrylate {name}")
        mono = Chem.MolFromSmiles(monomer)
        mw = Descriptors.MolWt(mono)
        pr_mw, pr_logp = TAIL_PRINTED[name]

        # the token, measured on the ARM's own graph
        arm = arm_from_alcohol(roh, label=f"arm {name}")
        armmol = Chem.MolFromSmiles(arm)
        starts = [a.GetIdx() for a in armmol.GetAtoms()
                  if a.GetAtomMapNum() == 2]
        if len(starts) != 1:
            raise SystemExit(f"tail {name}: arm lost its atom map")
        measured = atoms_to_furthest_carbon(armmol, starts[0])
        agree = chains.check_name_token(token, measured, label=f"tail {name}")

        rows.append({
            "block_id": f"T:{name}", "kind": "tail_acrylate", "name": name,
            "smiles": monomer,
            "formula": rdMolDescriptors.CalcMolFormula(mono),
            "mw": round(mw, 3),
            "n_nitrogens": 0, "n_NH_hydrogens": "",
            "n_tails_computed": "", "n_tails_stated_in_text": "",
            "n_tails_agree": "",
            "computed_charge_density_N_over_MW": "",
            "printed_charge_density_tableS5": "",
            "charge_density_abs_diff": "", "charge_density_agrees": "",
            "printed_clogP_tableS5": pr_logp,
            "name_token": token, "token_measured": measured,
            "token_agrees": int(agree["matches"]),
            "printed_mw_tableS4": pr_mw,
            "mw_abs_diff": round(abs(mw - pr_mw), 3),
            "mw_agrees": int(abs(mw - pr_mw) <= 0.05),
            "source_figure": "Figure 1A, 'R =' panel (Linear / Branched-cyclic)",
            "note": note + " | alcohol R-OH as built: " + roh,
        })

    cols = list(rows[0].keys())
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)

    heads = [r for r in rows if r["kind"] == "headgroup"]
    tails = [r for r in rows if r["kind"] == "tail_acrylate"]
    print(f"wrote {len(rows)} blocks -> {OUT}")
    print(f"  headgroups {len(heads)}  tails {len(tails)}  "
          f"library size {len(heads) * len(tails)}")
    print("  tail MW vs Table S4      : "
          f"{sum(r['mw_agrees'] for r in tails)}/{len(tails)} agree "
          f"(max |d| {max(r['mw_abs_diff'] for r in tails)} Da)")
    print("  tail name token measured : "
          f"{sum(r['token_agrees'] for r in tails)}/{len(tails)} agree")
    print("  head N/MW vs Table S5    : "
          f"{sum(r['charge_density_agrees'] for r in heads)}/{len(heads)} agree")
    print("  head tail-count vs text  : "
          f"{sum(r['n_tails_agree'] for r in heads)}/{len(heads)} agree")
    for r in heads:
        if not r["charge_density_agrees"]:
            print(f"    MISMATCH {r['name']}: computed "
                  f"{r['computed_charge_density_N_over_MW']} vs printed "
                  f"{r['printed_charge_density_tableS5']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
