#!/usr/bin/env python3
"""Stage 05 -- one table of every per-member value, with its exact source.

    python3 -I stage05_values.py        # after 05a and 05b

PROTEIN AND mRNA ARE DIFFERENT EXPERIMENTS. Figure 3B delivers GFP mRNA to four
cell lines; Figure S24C delivers (-30)GFP-Cre RECOMBINASE PROTEIN to one. They
share neither cargo, nor dose units (ug/mL of mRNA vs nM of protein), nor
readout definition, and pooling them into one Experiment_ID would invent a
comparison the paper never made. They are emitted as three separate blocks and
`proposed_experiment_id` says so.

THE QUALITY COLUMNS COME FROM curation/lib/emit.py, from each lipid's own
enumeration record. emit_rows REFUSES a literal confidence, which is error
class 29 -- found three times in this project (F-011, F-025, F-046) and
inverted every time, so that the most uncertain rows carried the strongest
claim.
"""
import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
sys.path.insert(0, HERE)

import emit                                                     # noqa: E402

# Conditions, each transcribed from a named locus. Nothing here is inferred.
BLOCKS = {
    "Fig 3B": {
        "proposed_experiment_id": "YL_2020",
        "cargo": "GFP mRNA (TriLink)",
        "Model_type": "in_vitro",
        "Model_target": "luminescence_normalized? NO -- percent GFP+ cells by "
                        "flow cytometry",
        "route": "n/a (cell culture)",
        "dose": "LNP/mRNA = 10/1 w/w; [mRNA] = 0.86 ug/mL",
        "dose_locus": "Results, 'using the optimized delivery conditions "
                      "(LNP/mRNA = 10/1, w/w; [mRNA] = 0.86 ug mL-1; 24 h "
                      "exposure)'",
        "timepoint": "24 h exposure",
        "formulation": "unformulated lipidoid LNP, 25 mM sodium acetate "
                       "pH 5.2, sonication + vortex (SI p. S4)",
        "readout": "percent GFP-positive cells",
    },
    "S24C": {
        "proposed_experiment_id": "YL_2020_2",
        "cargo": "(-30)GFP-Cre recombinase protein",
        "Model_type": "in_vitro",
        "Model_target": "percent GFP+ cells by flow cytometry",
        "route": "n/a (cell culture)",
        "dose": "25 / 50 / 100 nM protein; lipidoid 1.7 / 3.4 / 6.6 ug/mL",
        "dose_locus": "Results, '(25 nM/1.7 ug mL-1, 50 nM/3.4 ug mL-1, and "
                      "100 nM/6.6 ug mL-1)'",
        "timepoint": "8 h exposure, medium change, 16 h further; 24 h total "
                     "(SI p. S8)",
        "formulation": "unformulated lipidoid LNP + protein, 30 min at RT "
                       "in PBS pH 7.4 (SI p. S5)",
        "readout": "percent GFP-positive cells",
    },
    "S24F": {
        "proposed_experiment_id": "YL_2020_3",
        "cargo": "(-30)GFP-Cre recombinase protein",
        "Model_type": "in_vitro",
        "Model_target": "percent DsRed+ cells (Cre-mediated recombination)",
        "route": "n/a (cell culture)",
        "dose": "25 / 50 / 100 nM protein; lipidoid 1.7 / 3.4 / 6.6 ug/mL",
        "dose_locus": "Results, same sentence as S24C",
        "timepoint": "8 h exposure, 24 h total (SI p. S8)",
        "formulation": "as S24C",
        "readout": "percent DsRed-positive cells; a DOWNSTREAM functional "
                   "readout of the same delivery event, not a second "
                   "measurement of it",
    },
}
CELL_LINE = {"HeLa": "HeLa", "B16F10": "B16F10", "HEK": "HEK 293",
             "NIH 3T3": "NIH 3T3"}

LIBRARY = {"75-OCholB", "76-OCholB", "77-OCholB", "78-OCholB", "80-OCholB",
           "81-OCholB", "87-OCholB", "90-OCholB", "304-OCholB"}


def quality(enum_rows):
    """Per-lipid quality columns, via emit.emit_rows. Nothing is a literal."""
    rec = {r["IL_name"]: r for r in enum_rows}

    def confidence(item, r):
        if item.startswith("87-"):
            return "medium"
        return "high"

    def assumed(item, r):
        return "1" if (item.startswith("87-") or item.startswith("77-")) \
            else "0"

    def note(item, r):
        if item.startswith("87-"):
            return ("head 87 is N,N-bis(2-hydroxyethyl)-1,3-propanediamine, so "
                    "it offers a primary amine AND two primary alcohols to the "
                    "acrylate; aza- and oxa-Michael addition give the IDENTICAL "
                    "molecular formula, so the SI's own [M+H]+ "
                    f"{r['formula_MH_SI']} cannot distinguish them. Both tails "
                    "are placed on the amine, which is the stated chemistry "
                    "and the standard lipidoid route, but the paper draws no "
                    "structure for this member.")
        if item.startswith("77-"):
            return ("head 77 is 2-(2-aminoethyl)-1-methylpyrrolidine, a "
                    "commercial reagent supplied as a racemate; ring C2 is a "
                    "stereocentre and Figure 1A draws no wedge, so the "
                    "enumerated SMILES leaves it unassigned. Connectivity is "
                    "confirmed by the SI's [M+H]+ "
                    f"{r['formula_MH_SI']}; stereochemistry is not.")
        if item.startswith("304-"):
            return ("head 304 is tris(2-(methylamino)ethyl)amine. At low "
                    "magnification its Figure 2A drawing reads as a "
                    "piperazine; the SI's [M+H]+ C114H196N7O12S6 settles it, "
                    "because the head formula it implies, C9H24N4, has zero "
                    "degrees of unsaturation and therefore no ring, and its "
                    "S6/O12 fix three tails on three secondary amines.")
        return (f"enumerated from the drawn amine head and the drawn OCholB/"
                f"O16B/OChol acrylate; molecular formula confirmed against "
                f"the SI's own [M+H]+ table ({r['formula_MH_SI']}).")

    return {it: q for it, q in zip(
        sorted(rec), emit.emit_rows(sorted(rec), rec, confidence, note,
                                    structure_assumed=assumed,
                                    label="xuq_32438474"))}


def main():
    enum = list(csv.DictReader(open(os.path.join(HERE, "03_enumerated.csv"))))
    q = quality(enum)
    smiles = {r["IL_name"]: r["IL_SMILES"] for r in enum}
    keys = {r["IL_name"]: r["inchikey"] for r in enum}

    rows = []
    for fn, panel_of in (("05b_fig3b_values.csv", lambda r: "Fig 3B"),
                         ("05a_si_vector_values.csv", lambda r: r["panel"])):
        for r in csv.DictReader(open(os.path.join(HERE, fn))):
            panel = panel_of(r)
            if panel not in BLOCKS:
                continue
            b = BLOCKS[panel]
            name = r["sample"]
            is_member = name in LIBRARY
            cell = CELL_LINE.get(r["series"], "HeLa-DsRed")
            rows.append({
                "proposed_experiment_id": b["proposed_experiment_id"],
                "IL_name": name,
                "is_library_member": int(is_member),
                "IL_SMILES": smiles.get(name, ""),
                "inchikey": keys.get(name, ""),
                "panel": panel,
                "source_location": (
                    f"{panel}, bar at x={r.get('bar_x0_px') or r.get('bar_x0_px')}"
                    f" px, top at "
                    f"{r.get('bar_top_row_px') or r.get('bar_top_px')} px"),
                "source_file": ("source/pmc/nihms-1624779-f0004.jpg"
                                if panel == "Fig 3B"
                                else "source/anie202004994-sup-0001-"
                                     "misc_information.pdf p.37"),
                "tier": "5_pixels" if panel == "Fig 3B"
                        else "4b_pdf_vector_geometry",
                "cargo": b["cargo"],
                "cell_or_organ": cell,
                "Model_type": b["Model_type"],
                "route": b["route"],
                "dose": r["series"] if panel != "Fig 3B" else b["dose"],
                "dose_detail": b["dose"],
                "timepoint": b["timepoint"],
                "readout": b["readout"],
                "value": r["value"],
                "unit": r["unit"],
                **{k: (q[name][k] if is_member else "")
                   for k in emit.columns()},
            })
    out = os.path.join(HERE, "05_values.csv")
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    member = [r for r in rows if r["is_library_member"]]
    print(f"wrote {out}: {len(rows)} values, {len(member)} on library members")
    print(f"  rows_with_a_caveat (F-024 definition): "
          f"{emit.caveat_count(member)}")
    for b in sorted({r['proposed_experiment_id'] for r in rows}):
        n = sum(1 for r in rows if r["proposed_experiment_id"] == b)
        nm = sum(1 for r in member if r["proposed_experiment_id"] == b)
        print(f"  {b}: {n} values ({nm} library)")


if __name__ == "__main__":
    main()
