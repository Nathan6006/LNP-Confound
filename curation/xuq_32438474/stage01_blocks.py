#!/usr/bin/env python3
"""Stage 01 -- the building blocks of the xuq_32438474 (YL_2020) library.

    python3 -I stage01_blocks.py

Two block families:

  AMINE HEADS   9 of them, DRAWN in Figure 1A (75, 76, 77) and Figure 2A
                (78, 80, 81, 87, 90, 304). Read off the 300 dpi PMC rasters.
  ACRYLATE TAILS  3 of them, DRAWN in Figure 1A and detailed in Figures S1
                (OCholB), S6 (O16B) and S7 (OChol), with the SI's synthesis
                text naming every reagent.

EVERY BLOCK IS CHECKED AGAINST AN EXTERNAL NUMBER, not against itself. The SI
gives the [M+H]+ MOLECULAR FORMULA of every fully substituted lipidoid in two
text tables (Figure S13 p. S25, Figure S16 p. S28) and of the three single-
tailed ones (Figure S22E p. S34). Subtracting n x tail from the lipidoid
formula gives the head formula WITHOUT ever consulting our drawing. That is the
check stage03 runs. It is necessary and NOT sufficient: a formula cannot see
connectivity, cannot see which nitrogen an arm sits on, and cannot see
stereochemistry (CLAUDE.md s6; error class 14).

THE STEROL IS NOT TYPED HERE. It is taken verbatim from
provenance/reagent_registry.py via curation/lib/reagents.registered(). This
project built epicholesterol into 537 rows (F-013) and the epimers agree on
formula, exact mass to 4 dp, the InChIKey SKELETON block and all 17 descriptor
columns -- so the structure has to come from the one place that writes it from
its chemical name, and the product has to be checked on its FULL InChIKey.
"""
import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))

from rdkit import Chem                                          # noqa: E402
from rdkit.Chem.rdMolDescriptors import CalcMolFormula          # noqa: E402

import assemble                                                 # noqa: E402
import reagents                                                 # noqa: E402

# ---------------------------------------------------------------- heads ----
# Read off the drawn structures. source_figure names the panel and the file
# the structure was read from, so a reviewer can re-read the same pixels.
HEADS = [
    ("75",  "CCN(CC)CCN",        "Fig 1A (f0002)", "N,N-diethylethylenediamine"),
    ("76",  "NCCN1CCCC1",        "Fig 1A (f0002)", "1-(2-aminoethyl)pyrrolidine"),
    ("77",  "CN1CCCC1CCN",       "Fig 1A (f0002)", "2-(2-aminoethyl)-1-methylpyrrolidine; C2 of the ring is a stereocentre and the paper draws NO wedge"),
    ("78",  "CN(C)CCN",          "Fig 2A (f0003)", "N,N-dimethylethylenediamine"),
    ("80",  "CN(C)CCCN",         "Fig 2A (f0003)", "N,N-dimethyl-1,3-propanediamine"),
    ("81",  "CCN(CC)CCCN",       "Fig 2A (f0003)", "N,N-diethyl-1,3-propanediamine"),
    ("87",  "OCCN(CCO)CCCN",     "Fig 2A (f0003)", "N,N-bis(2-hydroxyethyl)-1,3-propanediamine"),
    ("90",  "NCCN1CCOCC1",       "Fig 2A (f0003)", "4-(2-aminoethyl)morpholine"),
    ("304", "CNCCN(CCNC)CCNC",   "Fig 2A (f0003)", "tris(2-(methylamino)ethyl)amine; three arms, NOT a piperazine -- see 08_contradictions.csv C-1"),
]

# ------------------------------------------------- tail precursor acids ----
# Each is a carbamic / carboxylic acid that esterifies ONTO an alcohol. The
# acid loses its OH and the alcohol keeps its oxygen (assemble.esterify), which
# is both what the chemistry does and the only way to touch a sterol's C3
# without the F-013 mechanism.
ACRYLIC = "C=CC(=O)O"
# propanoic acid carrying an atom map on C3: the SATURATED form of acrylic
# acid, used to build the Michael-adduct handle. Formula differs from ACRYLIC
# by exactly H2, which stage03 asserts.
PROPANOIC_MAPPED = "[CH3:2]CC(=O)O"

CARBAMIC = {
    # OCholB: cholesteryl chloroformate + Py-SS-NH2, then 2-mercaptoethanol
    # displaces the pyridyl (SI p. S2). Net alcohol = HO-CH2CH2-S-S-CH2CH2-NH-
    # carbamate-cholesteryl.
    "OCholB": ("OC(=O)NCCSSCCO",
               "SI p. S2 'Synthesis of Py-SS-Chol / OH-SS-Chol / OCholB'; "
               "Figure S1"),
    # OChol: cholesteryl chloroformate + 6-amino-1-hexanol, named and weighed
    # in the SI text (1.29 g, 11.0 mmol).
    "OChol": ("OC(=O)NCCCCCCO",
              "SI p. S3 'Synthesis of OChol'; Figure S7"),
}
# O16B carries no sterol: 1-dodecanethiol -> Py-SS-C12 -> 2-mercaptoethanol
# displacement -> HO-CH2CH2-S-S-C12H25. The alkyl is DODECYL, established two
# independent ways (see 07_qc.json check 'o16b_chain'); the NAME TOKEN '16' is
# not the chain length (CLAUDE.md s5).
O16B_ALCOHOL = "OCCSSCCCCCCCCCCCC"


def build_tails():
    """The three Michael acceptors, and their mapped saturated handles."""
    sterol = reagents.registered("Cholesterol")
    if not sterol:
        raise SystemExit("Cholesterol is not pinned in reagent_registry.py.")
    out = {}
    alcohols = {}
    for name, (acid, prov) in CARBAMIC.items():
        alcohols[name] = assemble.esterify(acid, sterol,
                                           label=f"{name}-carbamate")
    alcohols["O16B"] = assemble.canonical(O16B_ALCOHOL)
    for name in ("OCholB", "O16B", "OChol"):
        alc = alcohols[name]
        acrylate = assemble.esterify(ACRYLIC, alc, label=f"{name}-acrylate")
        handle = assemble.esterify(PROPANOIC_MAPPED, alc,
                                   label=f"{name}-handle")
        out[name] = {"alcohol": alc, "acrylate": acrylate, "handle": handle,
                     "sterol_from_registry": name != "O16B",
                     "provenance": CARBAMIC.get(name, (None, "SI p. S3 "
                                 "'Synthesis of O16B'; Figure S6"))[1]}
    return out, sterol


def main():
    tails, sterol = build_tails()
    rows = []
    for num, smi, fig, note in HEADS:
        m = Chem.MolFromSmiles(smi)
        if m is None:
            raise SystemExit(f"head {num} unparseable")
        rows.append({
            "block_type": "amine_head", "block_id": num,
            "SMILES": Chem.MolToSmiles(m),
            "formula": CalcMolFormula(m),
            "source_figure": fig, "note": note,
        })
    for name, d in tails.items():
        am = Chem.MolFromSmiles(d["acrylate"])
        rows.append({
            "block_type": "acrylate_tail", "block_id": name,
            "SMILES": d["acrylate"],
            "formula": CalcMolFormula(am),
            "source_figure": d["provenance"],
            "note": ("sterol taken VERBATIM from provenance/reagent_registry.py"
                     if d["sterol_from_registry"] else
                     "no sterol; dodecyl alkyl, see 07_qc.json o16b_chain"),
        })
    out = os.path.join(HERE, "01_blocks.csv")
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    with open(os.path.join(HERE, "01_tails.json"), "w") as fh:
        json.dump({"sterol_from_registry": sterol, "tails": tails}, fh, indent=1)
    print(f"wrote {out}: {len(rows)} blocks "
          f"({sum(1 for r in rows if r['block_type']=='amine_head')} heads, "
          f"{sum(1 for r in rows if r['block_type']=='acrylate_tail')} tails)")
    for name, d in tails.items():
        print(f"  {name:7s} acrylate {CalcMolFormula(Chem.MolFromSmiles(d['acrylate']))}")


if __name__ == "__main__":
    main()
