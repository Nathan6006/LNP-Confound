#!/usr/bin/env python3
"""Formulation + metadata, and the full value ledger. -> 05_values.csv

Run from this folder:  python3 -I stage05_formulation.py

EVERY VALUE THIS PAPER OFFERS GETS A ROW, admissible or not, with
`corpus_admissible` and a reason, and the counts sum. A value that is left out
silently is indistinguishable from a value nobody looked for.

THE FORMULATION IS DERIVED, NOT ASSUMED. The SI states the Dlin-MC3 control as
50 / 10 / 38.5 / 1.5 mol% (ionizable / DSPC / cholesterol / PEG2k-C-DMG) and
gives an explicit rule for the RAIL formulations: the ionizable lipid is fixed
at 33.3 mol% and the other three keep the control's relative ratio, scaled to
the remaining 66.7 mol%. stage06 checks the derivation against that rule and
also records that the control's four numbers are Onpattro's published
composition exactly -- a round number the authors used and never printed.

WHAT IS NOT RECOVERABLE, stated rather than invented:
  * N:P and any mass or weight ratio. ZERO occurrences in the SI or the main
    text; emitted blank.
  * the in-vitro dose in micrograms. The SI gives a CONCENTRATION, 0.5 ug/mL,
    and never states the well volume, so the absolute figure the corpus column
    wants is not recoverable. Emitted blank, concentration in the note.
"""
import collections
import csv
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "05_values.csv")

# --- from the SI's RAIL-LNP formation paragraph and Figure 1's caption -----
CONTROL_MOLPCT = {"IL": 50.0, "HL": 10.0, "CHL": 38.5, "PEG": 1.5}
RAIL_IL_MOLPCT = 33.3

META = {
    "vitro": {
        "Model": "in_vitro", "Model_type": "DC2.4", "Model_target": "in_vitro",
        "Route_of_administration": "in_vitro",
        "Cargo": "mRNA", "Cargo_type": "FLuc",
        "Experiment_method": "luminescence_normalized",
        "timepoint": "6 h incubation",
        "assay": "Varioskan Flash microplate reader luminescence",
        "dose_as_printed": "0.5 ug mL-1 mRNA, 96-well plate, 1e5 cells/well",
        "n": 3,
    },
    "vivo": {
        "Model": "in_vivo", "Model_type": "BALB/c",
        "Route_of_administration": "intravenous",
        "Cargo": "mRNA", "Cargo_type": "FLuc",
        "Experiment_method": "luminescence_normalized",
        "timepoint": "6 h after IV administration",
        "assay": "IVIS Spectrum ex vivo organ bioluminescence after IP "
                 "D-luciferin 150 mg/kg; signal normalised to organ area",
        "dose_as_printed": "0.1 mg kg-1 mRNA",
        "n": 3,
    },
    "Aqueous_buffer": "citrate",          # 10 mM citrate, pH 4.0
    "Dialysis_buffer": "PBS",             # 1x PBS, 4 C, 2 h, 20 kDa MWCO
    "Mixing_method": "microfluidic",      # aqueous:ethanol 3:1 by volume
}


def rail_molpct():
    rest = 100.0 - RAIL_IL_MOLPCT
    denom = sum(CONTROL_MOLPCT[k] for k in ("HL", "CHL", "PEG"))
    out = {"IL": RAIL_IL_MOLPCT}
    for k in ("HL", "CHL", "PEG"):
        out[k] = rest * CONTROL_MOLPCT[k] / denom
    return out


def main():
    rows = []
    a = list(csv.DictReader(open(os.path.join(HERE, "04_fig1a_cells.csv"))))
    c = list(csv.DictReader(open(os.path.join(HERE, "04_fig1c_bars.csv"))))

    for r in a:
        rows.append({
            "value_id": f"fig1a|{r['IL_name']}|rep{r['replicate_row']}",
            "panel": "Figure 1a", "tier": "5_pixels",
            "source_location":
                f"main PDF page 4, image xref 81 (2051x2352 px, 300 dpi "
                f"native), heat map column {r['panel_column']} of 41, row "
                f"{r['replicate_row']} of 3",
            "IL_name": r["IL_name"],
            "readout_as_printed": "Luciferase expression in DC2.4 cells "
                                  "(colour scale -1.0 to 1.0, UNITLESS and "
                                  "unlabelled; the figure prints no axis "
                                  "title for the colorbar)",
            "units": "none (min-max rescaled onto [-1, +1])",
            "value": r["scaled_value"],
            "cargo": META["vitro"]["Cargo_type"],
            "cell_or_organ": META["vitro"]["Model_type"],
            "route": "in vitro",
            "dose": META["vitro"]["dose_as_printed"],
            "timepoint": META["vitro"]["timepoint"],
            "assay": META["vitro"]["assay"],
            "pooled": "no", "pool_size": "",
            "corpus_admissible": "yes",
            "corpus_admissible_reason":
                "Luciferase reporter expression in a cell line requires "
                "translation of the delivered mRNA, which is the corpus's "
                "delivery endpoint (EXTRACTION_BRIEF 3). The CAVEAT is the "
                "scale, not the endpoint: see DECISIONS.md D-2.",
        })

    for r in c:
        bound = r["value_is_bound"] == "upper"
        rows.append({
            "value_id": f"fig1c|{r['IL_name']}|{r['organ']}",
            "panel": "Figure 1c", "tier": "5_pixels",
            "source_location":
                f"main PDF page 4, image xref 81, grouped bar chart, group "
                f"{r['panel_column']} of 41, {r['organ']} series",
            "IL_name": r["IL_name"],
            "readout_as_printed": "Avg Radiance (x 10^8 p/s/cm2/sr)",
            "units": "1e8 photons s-1 cm-2 sr-1",
            "value": r["value_avg_radiance_e8"],
            "cargo": META["vivo"]["Cargo_type"],
            "cell_or_organ": r["organ"],
            "route": "intravenous, BALB/c",
            "dose": META["vivo"]["dose_as_printed"],
            "timepoint": META["vivo"]["timepoint"],
            "assay": META["vivo"]["assay"],
            "pooled": "no", "pool_size": "",
            "corpus_admissible": "yes",
            "corpus_admissible_reason":
                ("Whole-organ luciferase bioluminescence requires translation "
                 "of the delivered mRNA. This is delivery, not "
                 "biodistribution: the readout is light from the expressed "
                 "enzyme, not recovered cargo.")
                + (" RECORDED AS AN UPPER BOUND: the bar's whole body lies "
                   "within 4 px of the drawn baseline and is not separable "
                   "from the axis rule's own thickness." if bound else ""),
        })

    # Figure 1d, read in stage06 as a cross-check and deliberately NOT emitted
    qc = json.load(open(os.path.join(HERE, "07_qc.json")))
    for k, rc, rd in qc["V1_fig1d_vs_fig1c_ratio"]["worst"]:
        pass   # the full table lives in 07_qc.json
    rows.append({
        "value_id": "fig1d|ALL 41 spleen-to-liver ratios",
        "panel": "Figure 1d", "tier": "5_pixels",
        "source_location": "main PDF page 4, image xref 81, lower panel",
        "IL_name": "(all 41)",
        "readout_as_printed": "Ratio of Spleen/Liver",
        "units": "dimensionless ratio",
        "value": "read, not emitted",
        "cargo": META["vivo"]["Cargo_type"],
        "cell_or_organ": "spleen / liver",
        "route": "intravenous, BALB/c",
        "dose": META["vivo"]["dose_as_printed"],
        "timepoint": META["vivo"]["timepoint"],
        "assay": META["vivo"]["assay"],
        "pooled": "no", "pool_size": "",
        "corpus_admissible": "no",
        "corpus_admissible_reason":
            "A SELECTIVITY CONTRAST, not a poolable readout. CLAUDE.md "
            "section 8 rules this out by the molea.csv 'mTP_Difference' "
            "precedent, and two earlier curations held back spleen/liver "
            "ratio panels for the same reason. It is also arithmetically "
            "derived from the two Figure 1c quantities we DO emit, so "
            "emitting it would double-count them. Read anyway, as an "
            "independent 26-way check on the Figure 1c reads (stage06 V1).",
    })

    for panel, what, why in [
        ("Figure S30", "RAIL molar-fraction optimisation in vivo",
         "A ONE-LIPID FORMULATION SCREEN: the varying factor is the RAIL mol "
         "fraction, not the ionizable-lipid chemistry. EXTRACTION_BRIEF 2 -- "
         "a corpus keyed on ionizable-lipid chemistry gains nothing from it "
         "however many rows exist."),
        ("Figures S31-S33", "particle size, PDI, zeta potential",
         "Not delivery (EXTRACTION_BRIEF 3)."),
        ("Figure S34", "cell viability",
         "Not delivery (EXTRACTION_BRIEF 3)."),
        ("Figure S35", "endosomal-escape quantification",
         "Not delivery; a trafficking readout upstream of translation."),
        ("Figure 2d,e", "time-dependent cellular uptake",
         "CELLULAR UPTAKE is explicitly not delivery, and the paper's own "
         "abstract says expression 'was not explained by total cellular "
         "uptake alone'."),
        ("Figure 2g", "relative FLuc expression, 4 cell lines, 2 LNPs",
         "Delivery in kind, but only TWO formulations (Dlin-MC3 and "
         "AT4-13CRA) and a relative scale with no stated reference. It adds "
         "no ionizable-lipid chemistry the Figure 1a rows do not already "
         "carry, and pooling a 2-point relative panel with the 41-point "
         "screen would need a tie neither panel provides."),
        ("Figure 2i-k", "transcriptomics, volcano, GSEA",
         "Host response, not delivery."),
        ("Figure 3 onward", "SORT comparison, Ai9, OVA vaccination, tumour",
         "Figure 3's splenic quantification is a THREE-formulation SORT "
         "comparison (Dlin-MC3, Dlin-MC3/18PA, AT4-13CRA) at a different "
         "formulation; the rest is immunology and tumour growth. Neither "
         "adds ionizable-lipid chemistry beyond the 40 already emitted."),
        ("Figures S36-S37", "corona proteomics, serum biochemistry",
         "Not delivery."),
    ]:
        rows.append({
            "value_id": f"{panel}|NOT EMITTED", "panel": panel,
            "tier": "n/a", "source_location": "main PDF / SI",
            "IL_name": "", "readout_as_printed": what, "units": "",
            "value": "not extracted", "cargo": "", "cell_or_organ": "",
            "route": "", "dose": "", "timepoint": "", "assay": "",
            "pooled": "", "pool_size": "",
            "corpus_admissible": "no", "corpus_admissible_reason": why})

    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    cnt = collections.Counter(r["corpus_admissible"] for r in rows)
    print(f"wrote {len(rows)} ledger rows to {OUT}")
    print(f"  corpus_admissible: {dict(cnt)}  (sums to {sum(cnt.values())})")
    print(f"  admissible individual values: "
          f"{sum(1 for r in rows if r['corpus_admissible'] == 'yes')} "
          f"= 123 Figure 1a cells + 82 Figure 1c bars")
    f = rail_molpct()
    print(f"  RAIL formulation derived: IL {f['IL']:.1f} / DSPC {f['HL']:.2f} "
          f"/ chol {f['CHL']:.3f} / PEG {f['PEG']:.3f} = "
          f"{sum(f.values()):.4f} mol%")
    print(f"  Dlin-MC3 control: {CONTROL_MOLPCT} (= Onpattro's published "
          "composition, which the SI never says)")


if __name__ == "__main__":
    main()
