#!/usr/bin/env python3
"""Blocks (01_blocks.csv) and the enumerated library (03_enumerated.csv).

WHAT CAN BE WRITTEN AND WHAT CANNOT.

The BLOCKS can be written: every amine and every tail lipid in this paper is a
catalogue reagent named in full in the SI design tables, so each SMILES here is
written from the CHEMICAL NAME and then checked against the formula that name
implies.  Two of my own first drafts failed that check and were corrected by it
(an epoxide one carbon too long; a macrocycle whose nitrogens did not fall on
1,4,8,12), which is the whole point of checking a name-derived structure rather
than trusting the typing.

The PRODUCTS cannot.  A lipomer is an oligoamine or polyamine carrying an
unstated number of tails:

  * 7 of the 15 lipomers are built on PEI 600 or PEI 1200.  Branched
    poly(ethylenimine) of nominal MW 600 or 1200 is a POLYDISPERSE MIXTURE of
    branched oligomers -- there is no single molecule to write, exactly as
    CLAUDE.md section 3 says of PEG-lipids ("the block length is a
    distribution, so no single SMILES is right").  The SI's own scheme admits
    this: it draws PEI_600 as diethylamine with the label underneath.
  * For the other 8, the amine IS a defined molecule, but the paper never
    states the degree of substitution.  Spermine has 6 N-H; a spermine/lauryl
    acrylate lipomer may carry 1 to 6 tails, i.e. 442 to 1,653 Da.  That is not
    a feature to flag on a row, it is the molecule being unknown.

So `product_SMILES` is EMPTY on all 15 rows and `blocked_reason` says which of
the two reasons applies.  Guessing full substitution would put a specific
invented molecule into a corpus keyed on structure.
"""
import csv
import json
import os
import sys

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors as rdmd

HERE = os.path.dirname(os.path.abspath(__file__))

# name -> (SMILES written from the chemical NAME, expected formula, role, note)
AMINES = {
    "PEI 600": (
        "", "",
        "amine",
        "Branched poly(ethylenimine), nominal MW 600 Da. A POLYDISPERSE "
        "MIXTURE of branched oligomers with primary, secondary and tertiary "
        "nitrogens in a distribution; no single SMILES is correct. SI Fig "
        "S2a/S3d/S3f/S5a draw it as diethylamine with 'PEI 600' written "
        "underneath, which is an abstraction, not a structure."),
    "PEI 1200": (
        "", "",
        "amine",
        "Branched poly(ethylenimine), nominal MW 1200 Da. Same reason as "
        "PEI 600."),
    "Spermine Free Base": (
        "NCCCNCCCCNCCCN", "C10H26N4",
        "amine",
        "N,N'-bis(3-aminopropyl)butane-1,4-diamine. 6 N-H available for "
        "alkylation."),
    "Spermidine": (
        "NCCCNCCCCN", "C7H19N3",
        "amine",
        "N-(3-aminopropyl)butane-1,4-diamine. 5 N-H available."),
    "Triethylenetetramine": (
        "NCCNCCNCCN", "C6H18N4",
        "amine",
        "N,N'-bis(2-aminoethyl)ethane-1,2-diamine. 6 N-H available."),
    "1,4,8,12-Tetraazacyclopentadecane": (
        "N1CCNCCCNCCCNCCC1", "C11H26N4",
        "amine",
        "15-membered macrocycle; the four ring nitrogens are asserted to sit "
        "at positions 1, 4, 8 and 12 and stage05 WALKS THE RING to confirm it "
        "rather than trusting the SMILES as typed. 4 N-H available."),
}

TAILS = {
    "C12 Epoxy": (
        "CCCCCCCCCCC1CO1", "C12H24O",
        "tail",
        "1,2-epoxydodecane. The name token is the TOTAL chain length: the "
        "alkyl substituent on the oxirane is C10, which is what SI Fig S5a "
        "labels 'R1=C10' for this same reagent."),
    "C15 Epoxy": (
        "CCCCCCCCCCCCCC1CO1", "C15H30O",
        "tail",
        "1,2-epoxypentadecane; alkyl substituent C13, which is SI Fig "
        "S2a/S3f's 'R2=C13'."),
    "Lauryl Acrylate": (
        "C=CC(=O)OCCCCCCCCCCCC", "C15H28O2",
        "tail", "Dodecyl prop-2-enoate; Michael acceptor."),
    "Stearyl Methacrylate": (
        "CC(=C)C(=O)OCCCCCCCCCCCCCCCCCC", "C22H42O2",
        "tail", "Octadecyl 2-methylprop-2-enoate; Michael acceptor."),
}

EXCIPIENTS = {
    # taken verbatim from provenance/reagent_registry.py, not retyped here
    "Cholesterol": "reagent_registry.py",
    "DOPE": "reagent_registry.py",
    "DOPC": "reagent_registry.py",
}

# lipomer id -> (amine, tail).  Read off the design tables, not assumed.
LIPOMERS = {
    "7C1": ("PEI 600", "C15 Epoxy"),
    "104": ("PEI 600", "C12 Epoxy"),
    "68":  ("PEI 600", "Lauryl Acrylate"),
    "69":  ("PEI 600", "Stearyl Methacrylate"),
    "106": ("PEI 1200", "C12 Epoxy"),
    "78":  ("PEI 1200", "Lauryl Acrylate"),
    "79":  ("PEI 1200", "Stearyl Methacrylate"),
    "100": ("Spermine Free Base", "C12 Epoxy"),
    "48":  ("Spermine Free Base", "Lauryl Acrylate"),
    "102": ("Spermidine", "C12 Epoxy"),
    "58":  ("Spermidine", "Lauryl Acrylate"),
    "108": ("Triethylenetetramine", "C12 Epoxy"),
    "88":  ("Triethylenetetramine", "Lauryl Acrylate"),
    "92":  ("1,4,8,12-Tetraazacyclopentadecane", "C12 Epoxy"),
    "8":   ("1,4,8,12-Tetraazacyclopentadecane", "Lauryl Acrylate"),
}

NH_SITES = {
    "Spermine Free Base": 6, "Spermidine": 5,
    "Triethylenetetramine": 6, "1,4,8,12-Tetraazacyclopentadecane": 4,
}


def ring_nitrogen_positions(smiles):
    """Walk the single ring and return the 1-based N positions."""
    m = Chem.MolFromSmiles(smiles)
    ring = m.GetRingInfo().AtomRings()[0]
    order, cur, prev = [], ring[0], None
    for _ in range(len(ring)):
        order.append(m.GetAtomWithIdx(cur).GetSymbol())
        nb = [a.GetIdx() for a in m.GetAtomWithIdx(cur).GetNeighbors()
              if a.GetIdx() in ring and a.GetIdx() != prev]
        prev, cur = cur, nb[0]
    return [i + 1 for i, s in enumerate(order) if s == "N"], len(ring)


def main():
    rows, report = [], {"formula_checks": [], "failures": []}
    for name, (smi, formula, role, note) in list(AMINES.items()) + list(TAILS.items()):
        rec = {"block": name, "role": role, "SMILES": "",
               "canonical_SMILES": "", "formula": "", "exact_mass": "",
               "formula_expected_from_name": formula,
               "formula_matches_name": "", "note": note}
        if smi:
            m = Chem.MolFromSmiles(smi)
            if m is None:
                report["failures"].append({"block": name, "why": "unparseable"})
                rows.append(rec)
                continue
            f = rdmd.CalcMolFormula(m)
            ok = (f == formula)
            rec.update({"SMILES": smi,
                        "canonical_SMILES": Chem.MolToSmiles(m),
                        "formula": f,
                        "exact_mass": round(Descriptors.ExactMolWt(m), 4),
                        "formula_matches_name": "yes" if ok else "NO"})
            report["formula_checks"].append(
                {"block": name, "expected": formula, "got": f, "match": ok})
            if not ok:
                report["failures"].append(
                    {"block": name, "why": "formula %s != %s" % (f, formula)})
        rows.append(rec)

    # the macrocycle's nitrogen positions, walked not trusted
    pos, size = ring_nitrogen_positions(
        AMINES["1,4,8,12-Tetraazacyclopentadecane"][0])
    report["macrocycle"] = {"ring_size": size, "n_positions_measured": pos,
                            "n_positions_required_by_name": [1, 4, 8, 12],
                            "match": pos == [1, 4, 8, 12] and size == 15}
    if not report["macrocycle"]["match"]:
        report["failures"].append({"block": "1,4,8,12-Tetraazacyclopentadecane",
                                   "why": "ring N positions %s" % pos})

    with open(os.path.join(HERE, "01_blocks.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        for r in rows:
            w.writerow(r)

    # ---- the enumerated library ----------------------------------------
    enum = []
    for lid, (am, tl) in sorted(LIPOMERS.items()):
        polymeric = am.startswith("PEI")
        if polymeric:
            reason = ("amine is a polydisperse polymer (%s): no single SMILES "
                      "exists for the amine, let alone for the product" % am)
        else:
            reason = ("degree of substitution not stated: %s offers %d N-H and "
                      "the paper never says how many carry a tail, so the "
                      "product's own molecular formula is open" %
                      (am, NH_SITES[am]))
        enum.append({
            "lipomer_id": lid,
            "amine": am,
            "amine_SMILES": AMINES[am][0],
            "tail": tl,
            "tail_SMILES": TAILS[tl][0],
            "reaction": ("epoxide ring-opening" if "Epoxy" in tl
                         else "aza-Michael addition"),
            "conditions_as_printed": "90 C, 72 hrs, EtOH (SI Fig S2a/S3d/S3f/S5a)",
            "product_SMILES": "",
            "degree_of_substitution_stated": "no",
            "n_NH_sites": "" if polymeric else NH_SITES[am],
            "blocked_reason": reason,
            "libraries": "",
        })
    # which libraries each lipomer appears in, read from the design tables
    with open(os.path.join(HERE, "01_design_tables.csv"), newline="") as fh:
        des = list(csv.DictReader(fh))
    for e in enum:
        libs = sorted({r["library"] for r in des if r["Lipomer"] == e["lipomer_id"]})
        e["libraries"] = "|".join(libs)
        e["n_formulations"] = sum(1 for r in des
                                  if r["Lipomer"] == e["lipomer_id"])

    with open(os.path.join(HERE, "03_enumerated.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(enum[0]))
        w.writeheader()
        for e in enum:
            w.writerow(e)

    report["n_blocks"] = len(rows)
    report["n_blocks_with_a_structure"] = sum(1 for r in rows if r["SMILES"])
    report["n_blocks_without_a_structure"] = sum(1 for r in rows if not r["SMILES"])
    report["n_lipomers"] = len(enum)
    report["n_lipomers_with_a_product_SMILES"] = 0
    report["n_lipomers_blocked_by_polymer_amine"] = sum(
        1 for e in enum if e["amine"].startswith("PEI"))
    report["n_lipomers_blocked_by_unstated_substitution"] = sum(
        1 for e in enum if not e["amine"].startswith("PEI"))
    report["n_formulations_covered"] = sum(e["n_formulations"] for e in enum)
    with open(os.path.join(HERE, "05_blocks_qc.json"), "w") as fh:
        json.dump(report, fh, indent=1)
    print(json.dumps(report, indent=1))
    return 1 if report["failures"] else 0


if __name__ == "__main__":
    sys.exit(main())
