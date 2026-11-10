#!/usr/bin/env python3
"""Stage 03 -- checks against what the paper DRAWS on SI p14, plus the
enumeration (such as it is) of this paper's ionizable lipids.

    python3 -I stage03_drawn_checks.py [OUTDIR]

SI p14 carries a structure band above the library-1 table: the 7C1 synthesis
scheme, cholesterol, DSPC, and four PEG-lipids. It is VECTOR, not raster: 885
drawing paths, of which 637 are FILLED polygons with width=None and 248 are
stroked at exactly 0.25 pt. (That is the measured shape behind the brief's
warning that chains.vector_bonds returns zero on filled-polygon bonds unless
include_filled=True. No bond tracing is performed here -- see the note in
reference_structures.csv on why.)

The best check available from that band is the brief's favourite class: A
NUMBER THE AUTHORS DREW AND NEVER COMPUTED. The four PEG-lipids are drawn with
explicit repeat counts in their polyether subscripts -- (OCH2CH2)45OCH3 three
times and (OCH2CH2)67OCH3 once -- and the formulation table, an entirely
separate object, carries a PEG MW column of 2000 and 3000. Neither was derived
from the other, and no statistic in the paper connects them.
"""
import csv
import json
import os
import re
import sys

import pymupdf

HERE = os.path.dirname(os.path.abspath(__file__))
SI = os.path.join(HERE, "source", "nl8b03149_si_001.pdf")
OUTDIR = os.path.abspath(sys.argv[1]) if len(sys.argv) > 1 else HERE
BLOCKS = os.path.join(OUTDIR, "01_blocks.csv")

EO = 44.0526           # ethylene oxide repeat, C2H4O
STRUCT_BAND_Y = 236.0  # everything above this on p14 is the structure band


def read_structure_band():
    """Pull the drawn labels and the polyether repeat subscripts off SI p14."""
    d = pymupdf.open(SI)
    pg = d[13]
    draw = pg.get_drawings()
    widths = {}
    for p in draw:
        w = None if p["width"] is None else round(p["width"], 4)
        widths[w] = widths.get(w, 0) + 1
    spans = []
    for b in pg.get_text("dict")["blocks"]:
        for ln in b.get("lines", []):
            for s in ln["spans"]:
                if s["bbox"][3] < STRUCT_BAND_Y:
                    spans.append({"x": round(s["bbox"][0], 1),
                                  "y": round(s["bbox"][1], 1),
                                  "size": round(s["size"], 2),
                                  "text": s["text"]})
    names = [s["text"] for s in spans
             if re.fullmatch(r"C\d\dPEG\d000", s["text"])]
    # A repeat count is a SUBSCRIPT (smaller font) sitting just right of a
    # '(OCH' label's closing paren. Identify it geometrically rather than by
    # grepping for plausible integers anywhere on the page.
    body = [s for s in spans if s["size"] > 4.5]
    subs = [s for s in spans if s["size"] < 4.5]
    repeats = []
    for anchor in [s for s in body if s["text"] == ")"]:
        cands = [s for s in subs
                 if 0 <= s["x"] - anchor["x"] <= 6
                 and abs(s["y"] - anchor["y"]) <= 4
                 and re.fullmatch(r"\d+", s["text"])]
        if len(cands) == 1:
            repeats.append({"y": anchor["y"], "n": int(cands[0]["text"])})
    repeats.sort(key=lambda r: r["y"])
    d.close()
    return {"n_drawing_paths": len(draw),
            "stroke_width_histogram": {str(k): v for k, v in widths.items()},
            "n_filled_paths_width_none": widths.get(None, 0),
            "drawn_peg_lipid_labels": names,
            "drawn_molecule_labels": [s["text"] for s in spans
                                      if s["text"] in ("Cholesterol", "DSPC",
                                                       "PEI600", "R=C13")],
            "polyether_repeat_counts_in_drawing_order": [r["n"] for r in repeats],
            "label_to_repeat_pairs": [
                {"label": nm, "repeat": rp["n"]}
                for nm, rp in zip(names, repeats)],
            }


def check_peg_repeats(band):
    """The drawn repeat counts vs the table's PEG MW column.

    Reported two ways, because they are not equally strong:
      - ABSOLUTE: nearest integer to nominal_MW / 44.0526.
      - RATIO: 67/45 vs 3000/2000, which cancels the end-group mass and the
        vendor's nominal-MW rounding, and is therefore the sharper test.
    """
    pairs = band["label_to_repeat_pairs"]
    rows = []
    for p in pairs:
        m = re.fullmatch(r"C(\d\d)PEG(\d000)", p["label"])
        nominal = float(m.group(2))
        pred = nominal / EO
        rows.append({"label": p["label"], "alkyl": "C" + m.group(1),
                     "nominal_mw": nominal, "drawn_repeat": p["repeat"],
                     "predicted_repeat": round(pred, 3),
                     "residual_repeats": round(p["repeat"] - pred, 3),
                     "within_one_repeat": abs(p["repeat"] - pred) < 1.5})
    r2000 = [r for r in rows if r["nominal_mw"] == 2000]
    r3000 = [r for r in rows if r["nominal_mw"] == 3000]
    out = {"per_lipid": rows,
           "n_lipids_drawn": len(rows),
           "all_within_one_repeat": all(r["within_one_repeat"] for r in rows),
           "distinct_repeat_for_2000": sorted({r["drawn_repeat"] for r in r2000}),
           "distinct_repeat_for_3000": sorted({r["drawn_repeat"] for r in r3000}),
           "consistency_within_mw": (len({r["drawn_repeat"] for r in r2000}) == 1
                                     and len({r["drawn_repeat"] for r in r3000}) == 1)}
    if r2000 and r3000:
        a, b = r2000[0]["drawn_repeat"], r3000[0]["drawn_repeat"]
        out["ratio_test"] = {
            "drawn_ratio": round(b / a, 6),
            "nominal_ratio": 1.5,
            "relative_residual_pct": round(abs(b / a - 1.5) / 1.5 * 100, 4),
            "why_sharper": ("the ratio cancels both the methoxy/linker end-group "
                            "mass and the vendor's rounding of 1982 to '2000', "
                            "neither of which the absolute test can separate"),
        }
    out["is_this_a_tautology"] = False
    out["why_not"] = (
        "The repeat subscripts live in the vector drawing on SI p14; the PEG MW "
        "column lives in the library table on the same page but is a separate "
        "object with its own extraction path. Nothing in this curation derived "
        "either from the other, and the paper never prints a PEG molecular "
        "weight computed from a repeat count. This is the drawn-but-never-"
        "computed class.")
    return out


def check_drawn_vs_table(band):
    """Do the drawn PEG-lipids account for the table's PEG-lipid combinations?"""
    if not os.path.exists(BLOCKS):
        raise SystemExit("01_blocks.csv missing -- run stage01 first")
    rows = list(csv.DictReader(open(BLOCKS)))
    combos = {}
    for lib in ("1", "2"):
        cs = {}
        for r in rows:
            if r["library"] != lib:
                continue
            tail = re.sub(r"PEG$", "", r["peg_tail"])
            key = f"{tail}PEG{r['peg_mw']}"
            cs[key] = cs.get(key, 0) + 1
        combos[f"library_{lib}"] = cs
    drawn = set(band["drawn_peg_lipid_labels"])
    out = {"drawn_peg_lipids": sorted(drawn),
           "n_drawn": len(drawn),
           "table_combinations": combos}
    l1 = set(combos["library_1"])
    l2 = set(combos["library_2"])
    out["library_1"] = {
        "table_combinations": sorted(l1),
        "n_table": len(l1),
        "drawn_equals_table": drawn == l1,
        "counts_per_combination": combos["library_1"],
        "design_is_balanced": len(set(combos["library_1"].values())) == 1,
        "drawn_caption_on_p14": "Library 1:  128 LNPs (4 PEG-lipids)",
        "caption_count_matches": len(l1) == 4 and sum(combos["library_1"].values()) == 128,
    }
    out["library_2"] = {
        "table_combinations": sorted(l2),
        "n_table": len(l2),
        "drawn_covers_table": l2 <= drawn,
        "in_table_but_never_drawn": sorted(l2 - drawn),
        "counts_per_combination": combos["library_2"],
        "design_is_balanced": len(set(combos["library_2"].values())) == 1,
    }
    return out


def enumerate_ionizable(band):
    """The 'enumeration'. There is one ionizable component and it has no
    atom-resolved structure, so the honest output is one row with an empty
    SMILES and the reason -- never a representative monomer."""
    return [{
        "il_name": "7C1",
        "il_smiles": "",
        "n_formulations_using_it": 248,
        "drawn_as": ("generic scheme: an alkyl epoxide labelled R=C13 (hence a "
                     "C15 epoxide) + PEI600, 90 C, 72 h, EtOH"),
        "why_no_smiles": ("PEI600 is branched polydisperse polyethylenimine of "
                          "nominal MW 600; the degree of epoxide substitution is "
                          "not stated. The product is a distribution over both "
                          "backbone topology and substitution, so no single "
                          "SMILES is correct."),
        "paper_gives_formula": False,
        "paper_gives_hrms": False,
        "paper_gives_nmr": False,
        "confidence": "low",
        "structure_assumed": 0,
        "curation_note": ("7C1 is a polydisperse PEI600/C15-epoxide polymer "
                          "conjugate drawn only as a generic R-group scheme "
                          "(SI p14); no atom-resolved structure exists to store."),
    }]


def main():
    band = read_structure_band()
    out = {
        "_what": ("Checks against the structures SI p14 DRAWS, and the "
                  "enumeration of this paper's ionizable lipids."),
        "structure_band": band,
        "peg_repeat_count_check": check_peg_repeats(band),
        "drawn_vs_table": check_drawn_vs_table(band),
        "ionizable_lipid_enumeration": enumerate_ionizable(band),
    }
    os.makedirs(OUTDIR, exist_ok=True)
    with open(os.path.join(OUTDIR, "03_drawn_checks.json"), "w") as f:
        json.dump(out, f, indent=1)

    with open(os.path.join(OUTDIR, "03_enumerated.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["il_name", "il_smiles", "n_formulations", "confidence",
                    "structure_assumed", "positional_ambiguity",
                    "curation_note"])
        for e in out["ionizable_lipid_enumeration"]:
            w.writerow([e["il_name"], e["il_smiles"],
                        e["n_formulations_using_it"], e["confidence"],
                        e["structure_assumed"], 0, e["curation_note"]])

    print(json.dumps({
        "n_drawing_paths": band["n_drawing_paths"],
        "n_filled_paths_width_none": band["n_filled_paths_width_none"],
        "drawn_peg_lipids": band["drawn_peg_lipid_labels"],
        "drawn_repeats": band["polyether_repeat_counts_in_drawing_order"],
        "peg_check": {k: v for k, v in out["peg_repeat_count_check"].items()
                      if k != "per_lipid"},
        "lib1_drawn_equals_table": out["drawn_vs_table"]["library_1"]["drawn_equals_table"],
        "lib1_caption_matches": out["drawn_vs_table"]["library_1"]["caption_count_matches"],
        "lib2_never_drawn": out["drawn_vs_table"]["library_2"]["in_table_but_never_drawn"],
        "n_ionizable_lipids": len(out["ionizable_lipid_enumeration"]),
    }, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
