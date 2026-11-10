#!/usr/bin/env python3
"""THE C10 ANCHOR -- 19T-O12, built independently of the enumerator.

WHAT AN ANCHOR IS FOR. C10 compares an enumerated structure against the
structure the paper DRAWS. It only tests anything if the anchor came from a
different route. That claim has been made twice in this project and was false
twice: F-013 (the stored anchor was a byte copy of the enumeration, same
sha256, while its note claimed independent fragment assembly) and F-037 (the
file written to FIX F-013 re-copied the other half -- its AMINE_HALF template
was the enumeration's canonical SMILES with the sterol swapped for `%` and a
ring digit renumbered). C16 exists because of that second one, and it tests the
GENERATOR'S TEMPLATE, not the emitted string.

SO WHAT IS ACTUALLY INDEPENDENT HERE, AND WHAT IS NOT. Stated plainly, because
over-claiming is the failure mode:

  INDEPENDENT  the FIGURE          Supplementary Fig. 1a (SI p.3, 400 dpi)
                                   vs. the enumerator's article Fig. 1a
               the MEASUREMENTS    02_si_linker_measurements.json, taken on
                                   that page with its own calibration, vs.
                                   01_chain_measurements.json
               the ASSEMBLY        plain string composition below. No
                                   `assemble`, no `attach`, no atom maps, no
                                   RWMol, no AddBond, no RemoveAtom anywhere in
                                   this file.
               the TRAVERSAL       written outward from the tertiary nitrogen.
                                   The enumerator's output is written from a
                                   tail inward, because it grew the molecule
                                   by bonding an arm onto an amine.

  NOT INDEPENDENT  the chemistry of the reaction. Both routes encode the same
                   reading of how 2-iminothiolane opens and how a thiol adds to
                   an acrylate. An anchor cannot test that; the two printed
                   masses do (07_qc.json P1 and P1b).

A correct re-derivation canonicalises to the same string as a correct
enumeration. Byte-identity of the EMITTED anchor is therefore not evidence of
copying and not evidence of independence -- which is exactly why C10's own
`anchor_may_be_a_copy` missed F-037. What this file offers a reviewer instead
is a template composed out of chemically named pieces that does not contain the
enumerator's output as a substring, and a `--verify` mode that re-derives it.

Run: python3 -I make_reference.py            # writes reference_structures.csv
     python3 -I make_reference.py --verify   # compares against 03_enumerated
"""
import csv
import json
import os
import sys

from rdkit import Chem

HERE = os.path.dirname(os.path.abspath(__file__))
SI = os.path.join(HERE, "02_si_linker_measurements.json")
ENUM = os.path.join(HERE, "03_enumerated.csv")
OUT = os.path.join(HERE, "reference_structures.csv")


def measured():
    with open(SI) as fh:
        d = json.load(fh)
    return d["core_from_this_drawing"], d["linker"], d["measurements"]


def compose():
    """19T-O12 as a SMILES string, written outward from the tertiary N.

    Every piece below is named for what it is in the drawing, and every repeat
    count comes from the SI-page measurement. Nothing is copied from
    03_enumerated.csv -- `--verify` reads that file, this function does not.
    """
    core, link, ms = measured()

    methyl = "C" * core["methyl_carbons"]                    # N-CH3
    spacer = "C" * core["N_to_amidine_NH_carbons"]           # N-CH2CH2-NH
    amidine = "NC(=N)"                                       # -NH-C(=NH)-
    thio_run = "C" * link["CH2_between_amidine_C_and_S"]     # -(CH2)3-
    sulfur = "S"
    michael = "C" * link["CH2_between_S_and_carbonyl"]       # -S-(CH2)2-
    ester = "C(=O)O"
    tail = "C" * ms["upper_arm_tail"]["heavy_atoms"]         # O-C12H25

    arm = spacer + amidine + thio_run + sulfur + michael + ester + tail
    # Outward from the tertiary nitrogen: methyl first, then the two arms.
    return methyl + "N(" + arm + ")" + arm


def main():
    core, link, ms = measured()
    smi = compose()
    m = Chem.MolFromSmiles(smi)
    if m is None:
        raise SystemExit(f"make_reference: composed SMILES does not parse: "
                         f"{smi}")
    canon = Chem.MolToSmiles(m)

    rows = [{
        "paper": "Wang et al., Nat Chem 2024 (PMID 38982196)",
        "compound": "19T-O12",
        "drawn_in": "Supplementary Fig. 1a, SI p.3 (the model T-MCR product)",
        "reference_SMILES": canon,
        "InChIKey": Chem.MolToInchiKey(m),
        "built_by": "make_reference.py compose(): plain string composition, "
                    "outward from the tertiary nitrogen. No assemble/attach, "
                    "no atom maps, no RWMol, no AddBond, no RemoveAtom.",
        "measurements_used": "02_si_linker_measurements.json only -- the "
                             "400 dpi SI page, its own calibration. The "
                             "enumerator uses 01_chain_measurements.json, the "
                             "article figure, a different scale.",
        "what_it_cannot_test": "The reaction chemistry itself: both routes "
                               "encode the same reading of how 2-iminothiolane "
                               "opens and how a thiol adds to an acrylate. "
                               "That is tested by the two printed masses "
                               "(07_qc.json P1, P1b), not by this anchor.",
        "note": "The one structure in this paper drawn in full. The article's "
                "Fig. 1a draws the 25 amines and 4 acrylates but shows the "
                "PRODUCT only as a scheme with R1/R2, so no library member "
                "other than this one is drawn.",
    }]
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {OUT}")
    print("composed template :", compose())
    print("canonical         :", canon)
    print("InChIKey          :", rows[0]["InChIKey"])


def verify():
    smi = Chem.MolToSmiles(Chem.MolFromSmiles(compose()))
    with open(ENUM) as fh:
        hit = [r for r in csv.DictReader(fh) if r["lipid_name"] == "19T-O12"]
    if not hit:
        raise SystemExit("make_reference --verify: 19T-O12 not enumerated.")
    other = hit[0]["IL_SMILES"]
    ok = Chem.MolToSmiles(Chem.MolFromSmiles(other)) == smi
    tmpl = compose()
    print(json.dumps({
        "anchor_canonical": smi,
        "enumerated_canonical": Chem.MolToSmiles(Chem.MolFromSmiles(other)),
        "match": ok,
        "anchor_InChIKey": Chem.MolToInchiKey(Chem.MolFromSmiles(smi)),
        "enumerated_InChIKey": hit[0]["InChIKey"],
        "template_contains_the_enumerated_string": other in tmpl,
        "template": tmpl,
        "C16_note": "A correct re-derivation is byte-identical to a correct "
                    "enumeration once canonicalised, so `match` being True is "
                    "the expected result and is NOT evidence of independence. "
                    "The checkable statement is the one above it: the "
                    "generator's composed template does not contain the "
                    "enumerator's output.",
    }, indent=1))
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    if "--verify" in sys.argv:
        verify()
    else:
        main()
