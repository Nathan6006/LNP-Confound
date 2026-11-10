#!/usr/bin/env python3
"""Build the C10 reference anchors for curation/xuq_24574196, and gate on them.

HOW INDEPENDENCE IS ESTABLISHED HERE -- BY PROVENANCE, NOT BY BYTES.
A correct hand-derivation is byte-identical to a correct enumeration once
canonicalised, so testing for byte-difference is a trap: that is exactly how
C10's own `anchor_may_be_a_copy` missed F-037, where `make_reference.py` --
the file written to FIX a byte-copied anchor -- carried a template that was the
enumeration's canonical SMILES with the sterol swapped for `%`. What is
claimed here instead, and is mechanically checkable by reading this file:

  * this module imports NOTHING from curation/lib and does not import
    stage03_enumerate. `assemble`, `attach` and `esterify` appear nowhere.
  * every fragment below is a SOURCE STRING LITERAL transcribed from the 1H
    NMR assignment list on SI p.3, with the shift and integral that fixes it
    written on the same line.
  * the traversal is REVERSED relative to the enumerator. stage03 starts at
    the amine and adds arms; this starts at the dimethylamino terminus and
    writes outward to the tail methyls, so the atom ordering in the template
    is the opposite of the one the enumerator produces.
  * there are NO RINGS in this family, so plain concatenation cannot produce
    the ring-closure collision that made fragment assembly mandatory
    elsewhere; `assert` below confirms zero rings in every product.

WHAT THE ANCHOR IS ANCHORED TO. This SI DRAWS NO STRUCTURE -- not one. Its
figures are two NMR spectra, two TEM micrographs and a viability bar chart.
So the external ground truth is the printed NMR ASSIGNMENT LIST, which states
connectivity directly (e.g. "4.35 (t, 4H, COOCH2CH2SS)" puts the ester oxygen
on a CH2 that is two bonds from a disulfide). That is weaker than a drawn
structure and is recorded as such in reference_structures.csv, rather than
being dressed up as a figure read.

  python3 -I make_reference.py          build, verify, write, exit 0 on success
  python3 -I make_reference.py --seed   seed a defect; MUST exit 1
"""
import csv
import os
import sys

from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors, inchi

RDLogger.DisableLog("rdApp.*")
HERE = os.path.dirname(os.path.abspath(__file__))

# --- fragments, read off the 1H NMR assignment list, SI p.3 ------------------
# 2.21 (s, 6H, N(CH3)2) and 2.26 / 1.37 / 2.49, three 2H multiplets, one
# propylene bridging the two nitrogens. Written from the dimethylamino end.
HEAD = "CN(C)CCCN"
# 2.92 (t, 4H, NCH2CH2CO) and 2.72 (t, 4H, NCH2CH2CO), then the ester carbonyl
# and its oxygen -- the Michael adduct of an acrylate on the second nitrogen.
ARM = "CCC(=O)O"
# the alcohol half of each ester, written outward to the terminal methyl.
TAILS = {
    # 4.05 (m, 4H, COOCH2CH2) + 56H bulk / 2 tails / 2 H = 14 CH2 + 0.87 CH3
    "1-O16": "CCCCCCCCCCCCCCCC",
    # 4.35 (t, 4H, COOCH2CH2SS) + 2.81 (t, 4H, COOCH2CH2SS) + the disulfide
    # + 2.49 (SSCH2, inside the 6H dd) + 40H bulk / 2 / 2 = 10 CH2 + 0.91 CH3
    "1-O16B": "CCSSCCCCCCCCCCCC",
}
PRINTED_FORMULA = {"1-O16B": "C39H78N2O4S4", "1-O16": "C43H86N2O4"}

NO_ANCHOR = {
    "1-O14B": "named in SI p.2 but no structure, formula, NMR or drawing "
              "exists for it anywhere in this SI",
    "1-O18B": "named in SI p.2 but no structure, formula, NMR or drawing "
              "exists for it anywhere in this SI",
}


def compose(tail):
    """Head-first template. Both arms hang off the SECOND nitrogen -- the one
    the NMR shows as fully substituted (no N-H, and the 2.92/2.72 pair
    integrates 4H each, i.e. two equivalent arms), never off the
    dimethylamino nitrogen, which the 6H singlet shows is still N(CH3)2."""
    return f"{HEAD}({ARM}{tail}){ARM}{tail}"


def main(seed=False):
    rows, failures = [], []
    for name, tail in TAILS.items():
        if seed and name == "1-O16B":
            # SEEDED DEFECT: read the name token "16" as a C16 S-alkyl, which
            # is the exact mistake this paper's nomenclature invites.
            tail = "CCSSCCCCCCCCCCCCCCCC"
        smi = compose(tail)
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            failures.append(f"{name}: template does not parse: {smi}")
            continue
        assert mol.GetRingInfo().NumRings() == 0, f"{name}: unexpected ring"
        assert all(a.GetNumRadicalElectrons() == 0 for a in mol.GetAtoms()), \
            f"{name}: radical in anchor"
        canon = Chem.MolToSmiles(mol)
        f = rdMolDescriptors.CalcMolFormula(mol)
        if f != PRINTED_FORMULA[name]:
            failures.append(
                f"{name}: anchor formula {f} != printed {PRINTED_FORMULA[name]}")
        rows.append({
            "lipid_name": name,
            "reference_SMILES": canon,
            "source_figure": "SI p.3, 1H NMR assignment list (this SI draws no "
                             "structure anywhere; see note)",
            "read_on": "2026-11-10",
            "read_by": "curation agent, curation/xuq_24574196",
            "expected_mismatch": "",
            "note": "Head-first (reversed) traversal from source string "
                    "literals; this module imports no assembler and no stage "
                    "module. Byte-identity to the enumeration after "
                    "canonicalisation is EXPECTED and is not evidence of "
                    "copying (F-037). Anchored to the printed NMR assignment "
                    "list, not to a drawn structure -- this SI contains no "
                    "structural drawing at all, which is weaker ground truth "
                    "than C10 normally rests on.",
            "anchor_formula": f,
            "anchor_inchikey": inchi.MolToInchiKey(mol),
        })
    for name, why in NO_ANCHOR.items():
        rows.append({
            "lipid_name": name, "reference_SMILES": "",
            "source_figure": "", "read_on": "2026-11-10",
            "read_by": "curation agent, curation/xuq_24574196",
            "expected_mismatch": "",
            "note": f"NO ANCHOR: {why}. Stored empty so the ABSENCE is "
                    "visible rather than silent. These two structures are "
                    "extrapolations and carry confidence = medium, "
                    "structure_assumed = 1.",
            "anchor_formula": "", "anchor_inchikey": ""})

    if failures:
        for f in failures:
            print("ANCHOR GATE FAILED:", f, file=sys.stderr)
        return 1

    header = ["lipid_name", "reference_SMILES", "source_figure", "read_on",
              "read_by", "expected_mismatch", "note"]
    with open(os.path.join(HERE, "reference_structures.csv"), "w",
              newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=header, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    with open(os.path.join(HERE, "reference_structures_detail.csv"), "w",
              newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=header + ["anchor_formula",
                                                    "anchor_inchikey"])
        w.writeheader()
        w.writerows(rows)
    for r in rows:
        print(f"{r['lipid_name']:10s} {r['anchor_formula']:16s} "
              f"{r['reference_SMILES'] or '(no anchor)'}")
    return 0


if __name__ == "__main__":
    sys.exit(main(seed="--seed" in sys.argv))
