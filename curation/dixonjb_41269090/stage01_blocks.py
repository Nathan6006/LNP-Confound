#!/usr/bin/env python3
"""Stage 01 -- the library's component blocks, and the Screening-1 grid.

Two outputs:

  01_blocks.csv       one row per named lipid component in the library, with
                      the structure we would use for it and WHERE that
                      structure came from (registry / guarded corpus lookup /
                      none), so a reader can see which components have a
                      pinned structure and which do not.

  01_screening1.csv   the 30-member Screening-1 grid, lifted from the SI page-2
                      NATIVE table (tier 4, text layer, no pixels).

No structure is written from memory here. Benchmarks come from
provenance/reagent_registry.py; anything else comes from
curation/corpus_components.py's guarded lookup, which RAISES on an ambiguous
name (error class 30).
"""
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "curation"))
sys.path.insert(0, os.path.join(ROOT, "provenance"))

import pdfplumber                                        # noqa: E402
from corpus_components import load, smi                  # noqa: E402
import reagent_registry as rr                            # noqa: E402

SI = os.path.join(HERE, "source", "nn5c12080_si_001.pdf")

# Every component named in the paper's library description (Results, "Generation
# of LNP Library") or in the Screening-1 table. role is as the paper uses it.
COMPONENTS = [
    ("cKK-E12", "ionizable", "CKK-E12",
     "BroadPharm BP-29590, named in Methods 'LNP Formulation'"),
    ("DLin-MC3-DMA", "ionizable", None,
     "MedKoo 555308, named in Methods as 'MC3'; benchmark"),
    ("C12S", "ionizable", None,
     "stereopure lipomer, 'provided from the Dahlman lab'; NO structure disclosed"),
    ("C12R", "ionizable", None,
     "stereopure lipomer, 'provided from the Dahlman lab'; NO structure disclosed"),
    ("C13S", "ionizable", None,
     "stereopure lipomer, 'provided from the Dahlman lab'; NO structure disclosed"),
    ("C13R", "ionizable", None,
     "stereopure lipomer, 'provided from the Dahlman lab'; NO structure disclosed"),
    ("Cholesterol", "sterol", "Cholesterol", "Avanti"),
    ("20alpha-hydroxycholesterol", "sterol", None,
     "named in Results; not in the registry and not in the corpus"),
    ("C14PEG2K", "PEG", "C14-PEG2000",
     "Avanti; a PEG block length is a DISTRIBUTION, so no single SMILES is right"),
    ("C18PEG2K", "PEG", None, "named in Results; Screening 2 only"),
    ("DSPC", "helper", "DSPC", "Avanti"),
    ("DOPE", "helper", "DOPE", "Avanti"),
    ("18:1 Cap PE", "helper", None,
     "named in Results and the Screening-1 table as 18:1CAPPE; not registered"),
    ("18:0 DDAB", "helper", "DDAB",
     "Screening-1 table token '18:0DDAB'; DDAB in the corpus is the C18 "
     "dioctadecyl quaternary ammonium, which is what 18:0 DDAB means"),
    ("DOTAP", "helper", "DOTAP", "Avanti"),
    ("DDAP", "helper", None, "named in Results; not in the Screening-1 table"),
    ("DOTMA", "helper", None, "named in Results as a Screening-2 component"),
]


def main():
    rows = []
    for f in ("data_vitro.csv", "data_vivo.csv"):
        with open(os.path.join(ROOT, "data", f)) as fh:
            rows += list(csv.DictReader(fh))
    resolved, ambiguous = load(rows)

    out = []
    for name, role, corpus_name, why in COMPONENTS:
        if name in rr.REGISTRY:
            s, authority = rr.REGISTRY[name]
            src = "reagent_registry"
        elif corpus_name:
            s = smi(resolved, ambiguous, corpus_name)   # RAISES if ambiguous
            src = f"corpus_components[{corpus_name}]"
            authority = why
        else:
            s, src, authority = "", "none", why
        out.append({"component": name, "role": role, "SMILES": s,
                    "structure_source": src, "note": authority})

    with open(os.path.join(HERE, "01_blocks.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["component", "role", "SMILES",
                                           "structure_source", "note"])
        w.writeheader()
        w.writerows(out)

    # --- Screening 1 grid, tier 4 ----------------------------------------
    with pdfplumber.open(SI) as doc:
        tabs = doc.pages[1].find_tables()
        if len(tabs) != 1:
            raise SystemExit(f"stage01: expected 1 native table on SI page 2, "
                             f"found {len(tabs)}. Refusing to guess.")
        grid = tabs[0].extract()
    hdr, body = grid[0], grid[1:]
    if len(body) != 30:
        raise SystemExit(f"stage01: Screening 1 should hold 30 LNPs, got "
                         f"{len(body)}.")
    with open(os.path.join(HERE, "01_screening1.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["LNP", "lipomer", "cholesterol", "PEG", "helper_lipid",
                    "lipomer_molpct", "cholesterol_molpct", "PEG_molpct",
                    "helper_molpct"])
        for r in body:
            w.writerow([c.strip() for c in r])

    n_struct = sum(1 for r in out if r["SMILES"])
    print(f"01_blocks.csv       {len(out)} components, {n_struct} with a "
          f"pinned structure, {len(out) - n_struct} without")
    print(f"01_screening1.csv   {len(body)} LNPs from the SI page-2 native "
          f"table (header: {hdr})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
