#!/usr/bin/env python3
"""The C10 anchor for PL_2018, BUILT -- not copied from the enumeration.

WHY THIS FILE EXISTS AT ALL. C10 compares each enumerated structure against the
structure the paper DRAWS. It tests nothing if the stored anchor came out of
the enumerator. That has been claimed twice in this project and been false
twice: F-013 (a byte copy, sha256-identical, 537 rows of epicholesterol behind
it) and F-037 -- inside the file written to FIX F-013, whose amine half was the
enumeration's canonical SMILES with the sterol replaced by '%' and one ring
digit renumbered. C16 therefore inspects THE GENERATOR'S COMPOSED TEMPLATE, not
the emitted string, because a correct re-derivation is byte-identical to a
correct enumeration once canonicalised.

HOW THIS ROUTE DIFFERS FROM THE ENUMERATOR, concretely:

  stage03_enumerate.py        this file
  --------------------------  ----------------------------------------
  curation/lib/assemble       no assemble, no RWMol, no AddBond, no
  .attach x 3, atom maps,     atom maps. Plain string concatenation of
  CombineMols on a shared     fragments that are each written from a
  ring-opened CORE template   CHEMICAL NAME, with no shared template.
  traversal: tail -> core     traversal: HEAD FIRST, in the order the
  -> arm -> head              eye reads Scheme S1c left to right.

The fragments below are written from the names the SI prints, read off the
drawing on SI p.2:

  "N1,N1-diethyl ethylene diamine"  -> CCN(CC)CC + N
  the new amide                     -> C(=O)
  the homocysteine alpha carbon     -> C(...)(N...)
  its two-carbon thioether arm      -> CCSS
  "2-(hexadecyldisulfanyl) pyridine" contributes hexadecyl -> 16 C
  "N-(2-oxotetrahydrothiophen-3-yl) oleamide" contributes oleoyl ->
       C(=O) + 7 C + cis C=C + 8 C

Chain lengths are NOT typed from the names: the 16 and the 18-with-9Z come
from stage01_chains.py, which measured them from the PDF's own bond
coordinates, and this file asserts against that file rather than hard-coding.

Run:  python3 -I make_reference.py
"""
import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

from rdkit import Chem                                      # noqa: E402
from rdkit import RDLogger                                  # noqa: E402

RDLogger.DisableLog("rdApp.*")

READ_ON = "2026-11-10"
READ_BY = "curation agent, from a 600 dpi render of SI p.2 Scheme S1c"


def fragments(n_arm, n_tail, dbl):
    """Every piece written from its chemical name. No shared template."""
    if dbl != [9, 10]:
        raise SystemExit(f"make_reference: oleoyl unsaturation measured at "
                         f"C{dbl[0]}=C{dbl[1]}; this builder writes 9Z.")
    before = dbl[0] - 2          # carbons between the carbonyl and C9
    after = n_tail - dbl[1]      # carbons after C10
    return {
        "diethylamino_head": "CCN(CC)",      # N,N-diethyl
        "ethylene_linker": "CC",             # ...ethylene diamine
        "new_amide": "NC(=O)",               # the ring-opened carbonyl
        "alpha_carbon_open": "C(",
        "thioether_arm": "CCSS",             # -CH2CH2-S-S-
        "hexadecyl": "C" * n_arm,
        "alpha_carbon_close": ")",
        "acylamino": "NC(=O)",               # the original oleamide N-H
        "oleoyl_pre": "C" * before,
        "oleoyl_cis": "/C=C\\",
        "oleoyl_post": "C" * after,
    }


def main():
    ch = json.load(open(os.path.join(HERE, "02_chain_counts.json")))["blocks"]
    n_arm = ch["PY16"]["carbons_from_geometry"]
    n_tail = ch["T18U"]["carbons_from_geometry"]
    dbl = ch["T18U"]["double_bond_between_carbons"]

    f = fragments(n_arm, n_tail, dbl)
    # HEAD FIRST, the order Scheme S1c reads: the diethylamine on the left,
    # through the new amide, out along the disulfide arm, then back down the
    # oleamide tail.
    built = (f["diethylamino_head"] + f["ethylene_linker"] + f["new_amide"]
             + f["alpha_carbon_open"] + f["thioether_arm"] + f["hexadecyl"]
             + f["alpha_carbon_close"] + f["acylamino"] + f["oleoyl_pre"]
             + f["oleoyl_cis"] + f["oleoyl_post"])

    mol = Chem.MolFromSmiles(built)
    if mol is None:
        raise SystemExit(f"make_reference: built string does not parse: {built}")
    canon = Chem.MolToSmiles(mol)
    key = Chem.MolToInchiKey(mol)

    enum = {r["lipidoid"]: r for r in csv.DictReader(
        open(os.path.join(HERE, "03_enumerated.csv")))}["T18U-PY16-A1"]
    emol = Chem.MolFromSmiles(enum["IL_SMILES"])
    ekey = Chem.MolToInchiKey(emol)

    report = {
        "built_string_before_canonicalisation": built,
        "canonical": canon,
        "inchikey": key,
        "enumerated_smiles": enum["IL_SMILES"],
        "enumerated_inchikey": ekey,
        "inchikey_match": key == ekey,
        "skeleton_match": key.split("-")[0] == ekey.split("-")[0],
        "built_string_is_byte_identical_to_enumeration":
            built == enum["IL_SMILES"],
        "chain_lengths_taken_from": "02_chain_counts.json (measured), not "
                                    "typed into this file",
        "n_arm": n_arm, "n_tail": n_tail, "double_bond": dbl,
    }
    if not report["inchikey_match"]:
        raise SystemExit(
            "make_reference: the anchor and the enumeration are DIFFERENT "
            f"molecules.\n  anchor      {canon}\n  enumeration "
            f"{enum['IL_SMILES']}\nThat is what C10 exists to find. Do not "
            "edit the anchor to agree -- find out which half is wrong.")
    if report["built_string_is_byte_identical_to_enumeration"]:
        raise SystemExit(
            "make_reference: the string this file BUILDS is byte-identical to "
            "the enumeration's canonical SMILES. That is the F-037 signature "
            "-- the anchor would be comparing a value to itself.")

    with open(os.path.join(HERE, "09_anchor_build.json"), "w") as fh:
        json.dump(report, fh, indent=1)

    rows = [
        {"lipid_name": "T18U-PY16-A1",
         "reference_SMILES": canon,
         "source_figure": "SI p.2, Scheme S1c (vector; the only full lipidoid "
                          "structure the SI draws)",
         "read_on": READ_ON, "read_by": READ_BY,
         "expected_mismatch": "",
         "note": "Built head-first by string concatenation of fragments each "
                 "written from the chemical name the SI prints; no RDKit "
                 "assembly, no atom maps, no shared template with "
                 "stage03_enumerate.py. Chain lengths and the 9Z position "
                 "come from 02_chain_counts.json, measured from the PDF's own "
                 "bond coordinates."},
        {"lipid_name": "T18-PY18-A5",
         "reference_SMILES": "",
         "source_figure": "SI p.15, Figure S1 top scheme (raster)",
         "read_on": READ_ON, "read_by": READ_BY,
         "expected_mismatch": "",
         "note": "NO ANCHOR STORED, deliberately, so the absence is visible. "
                 "Figure S1's top scheme does draw a complete primary-amine "
                 "product, but it is a JPEG raster whose chain lengths could "
                 "not be measured (the same noise floor that defeated A5's "
                 "N-ethyl branches, 01_amines.json), and an anchor whose "
                 "chain lengths were counted by eye is error class 48 with "
                 "extra steps."},
        {"lipid_name": "T18U-PY12-A8",
         "reference_SMILES": "",
         "source_figure": "none -- the SI draws no structure for A8",
         "read_on": READ_ON, "read_by": READ_BY,
         "expected_mismatch": "",
         "note": "NO ANCHOR POSSIBLE. A8 is known only as the molecular "
                 "formula C7H17N3 (01_amines.json). This is the paper's "
                 "SECOND-best lipidoid, so the gap is recorded here rather "
                 "than left to be noticed later. DECISIONS.md D-2."},
    ]
    with open(os.path.join(HERE, "reference_structures.csv"), "w",
              newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=[
            "lipid_name", "reference_SMILES", "source_figure", "read_on",
            "read_by", "expected_mismatch", "note"])
        w.writeheader()
        w.writerows(rows)

    print(json.dumps({k: report[k] for k in (
        "inchikey", "enumerated_inchikey", "inchikey_match",
        "built_string_is_byte_identical_to_enumeration",
        "n_arm", "n_tail", "double_bond")}, indent=1))
    print("built :", built)
    print("enum  :", enum["IL_SMILES"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
