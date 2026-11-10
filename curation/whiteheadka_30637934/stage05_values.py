#!/usr/bin/env python3
"""05_values.csv -- every delivery value this paper resolves, emitted or not.

26 rows: the 14 categories of Figure S6A and the 12 points of Figure S6B.
Eleven are emitted into 06_rows.csv. The other 15 carry a per-row reason.

Panels located but NOT digitised are in 00_value_sources.json, with the reason
there; the one that matters is Figure S4, which this stage DID measure
structurally: it is a STACKED bar chart (4 x-categories = Crl + three doses,
six organ segments summed within each bar), not the grouped-per-organ chart its
legend reads like. A naive bar-top read of it returns the six-organ SUM while
believing it is the liver. Recorded in 08_contradictions.csv.
"""
import csv
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))

COMMON = {
    "cargo": "firefly luciferase mRNA, UNMODIFIED (TriLink)",
    "cargo_type": "mRNA",
    "cell_or_organ": "HeLa (ATCC), 15,000 cells/well, white 96-well plate",
    "route": "in vitro, added to wells 6-8 h after seeding",
    "timepoint": "24 h post-transfection",
    "assay": "Bright-Glo Luciferase Assay System (Promega)",
    "readout_as_printed": "Relative Light Units",
    "units_as_printed": "RLU (arbitrary; plate-reader gain dependent -- the SI "
                        "warns of exactly this for its TNS readout, p.3)",
    "n": "4 (caption: 'Error bars represent s.d. (n = 4)')",
}

CONTROL_REASON = {
    "Crl": ("Untreated control, not an LNP. No bar is drawn at its x tick "
            "(122.94 pt) at all, so there is no value to read."),
    "Nkd": ("Naked mRNA control, not an LNP. No bar is drawn at its x tick "
            "(135.76 pt) at all, so there is no value to read."),
    "Lipo MM": ("Lipofectamine MessengerMax, a commercial transfection "
                "reagent, not a lipid nanoparticle: it has no ionizable "
                "lipid, no excipients and no formulation in the corpus "
                "schema's sense. Drawn in a different colour from the eleven "
                "LNP bars. Value measured (346,411 RLU) and kept here as the "
                "panel's internal scale reference."),
}
PANEL_B_REASON = (
    "Figure S6B is a SECOND PLATE. Its 100 ng series shares four lipidoids "
    "and the identical nominal dose with Figure S6A and disagrees with it by "
    "0.50x to 1.28x, rank-inverting the four (Spearman rho = -0.20). Emitting "
    "both under one Experiment_ID would pool two assays that demonstrably "
    "disagree, which is error class 22. Whether S6B becomes a separate "
    "Experiment_ID is DECISIONS.md D-3; neither panel is averaged into the "
    "other and neither is rescaled by the other.")


def main():
    with open(os.path.join(HERE, "02_fig_s6_panelA.csv"), newline="") as fh:
        pa = list(csv.DictReader(fh))
    with open(os.path.join(HERE, "02_fig_s6_panelB.csv"), newline="") as fh:
        pb = list(csv.DictReader(fh))
    geo = json.load(open(os.path.join(HERE, "02_fig_s6_geometry.json")))

    rows = []
    for r in pa:
        cat = r["category"]
        is_lnp = cat not in CONTROL_REASON
        rows.append({
            "value_id": f"S6A_{cat.replace(' ', '')}",
            "panel": "Figure S6A",
            "source_location": (
                "SI p.12 (smll201805097-sup-0001-S1.pdf), Figure S6 panel A; "
                f"bar outline top at y = {r['bar_top_px'] or 'no bar drawn'} pt, "
                f"x tick at {r['x_tick_px']} pt, width "
                f"{r['bar_outline_width_px'] or '-'} pt"),
            "extraction_tier": "4b_pdf_vector_geometry",
            "category": cat,
            "lipid": f"306{cat}" if is_lnp else "",
            "IL_name": f"A306_{cat}" if is_lnp else "",
            "value": r["value_RLU"],
            "sd": r["sd_RLU"],
            "mean_plus_sd": r["mean_plus_sd_RLU"],
            "dose": "100 ng mRNA per well (0.81 nM); 0.1 ug",
            "emitted_to_06_rows": "yes" if (is_lnp and r["value_RLU"]) else "no",
            "reason_if_not_emitted": CONTROL_REASON.get(cat, ""),
            **COMMON,
        })
    for r in pb:
        rows.append({
            "value_id": f"S6B_{r['category']}_{r['series'].replace(' ', '')}",
            "panel": "Figure S6B",
            "source_location": (
                "SI p.12, Figure S6 panel B; 2.80 pt series marker centred at "
                f"y = {r['marker_centre_px']} pt"),
            "extraction_tier": "4b_pdf_vector_geometry",
            "category": r["category"], "lipid": r["category"],
            "IL_name": "A306_" + r["category"].replace("306", ""),
            "value": r["value_RLU"], "sd": "", "mean_plus_sd": "",
            "dose": f"{r['series']} mRNA per well",
            "emitted_to_06_rows": "no",
            "reason_if_not_emitted": PANEL_B_REASON,
            **COMMON,
        })

    cols = ["value_id", "panel", "source_location", "extraction_tier",
            "category", "lipid", "IL_name", "readout_as_printed",
            "units_as_printed", "value", "sd", "mean_plus_sd", "cargo",
            "cargo_type", "cell_or_organ", "route", "dose", "timepoint",
            "assay", "n", "emitted_to_06_rows", "reason_if_not_emitted"]
    with open(os.path.join(HERE, "05_values.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader(); w.writerows(rows)

    qc = {
        "values_located": len(rows),
        "panelA_categories": len(pa),
        "panelA_with_a_drawn_bar": sum(1 for r in pa if r["value_RLU"]),
        "panelB_points": len(pb),
        "emitted": sum(1 for r in rows if r["emitted_to_06_rows"] == "yes"),
        "not_emitted": sum(1 for r in rows if r["emitted_to_06_rows"] == "no"),
        "counts_sum": (sum(1 for r in rows if r["emitted_to_06_rows"] == "yes")
                       + sum(1 for r in rows
                             if r["emitted_to_06_rows"] == "no")) == len(rows),
        "every_not_emitted_row_has_a_reason": all(
            r["reason_if_not_emitted"].strip()
            for r in rows if r["emitted_to_06_rows"] == "no"),
        "three_counts_that_differ": {
            "_what": ("Distinct structures, formulations, and individually "
                      "measured values are three different numbers and this "
                      "paper separates them cleanly."),
            "distinct_structures_in_the_306_family": 11,
            "distinct_structures_in_the_whole_library":
                ("NOT DETERMINABLE from what we hold. The SI Methods name "
                 "THREE amines (200, 304, 306) and the abstract names 'a "
                 "small library' plus 'a secondary lipidoid library', but "
                 "every value-bearing table and figure in the SI is "
                 "306-family only and the main text is behind a "
                 "subscription. So the library is larger than 11 by an "
                 "unknown factor and no value for a 200- or 304-amine "
                 "lipidoid is recoverable."),
            "formulations": 1,
            "formulation_detail": ("One formulation for the whole screen: "
                                   "35 : 16 : 46.5 : 2.5 mol%, 10:1 w/w "
                                   "lipidoid:mRNA. No formulation variable is "
                                   "screened anywhere in this paper."),
            "individually_measured_in_vitro_delivery_values": 11,
            "individually_measured_counting_the_second_plate": 11 + 12,
            "individually_measured_physicochemical_only": {
                "Table_S1_LNP_characterisation": 11,
                "Table_S2_lipidoid_pKa_and_buffering": 10,
                "Figure_S7_printed_pKa": 8},
            "in_vivo_values_resolved_per_lipidoid":
                ("Figure S3D: 3 lipidoids x a time course. Figure S4: 2 "
                 "lipidoids x 3 doses, as a six-organ STACKED total. Figure "
                 "S5: 2 lipidoids x 3 timepoints, EPO protein not luciferase. "
                 "None digitised; see DECISIONS.md D-2."),
        },
        "cross_panel_disagreement": geo["cross_panel_S6A_vs_S6B_100ng"],
    }
    with open(os.path.join(HERE, "05_values_qc.json"), "w") as fh:
        json.dump(qc, fh, indent=1)
    print(json.dumps({k: v for k, v in qc.items()
                      if k != "cross_panel_disagreement"}, indent=1))


if __name__ == "__main__":
    main()
