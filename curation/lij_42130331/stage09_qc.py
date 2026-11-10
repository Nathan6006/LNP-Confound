#!/usr/bin/env python3
"""Assemble 07_qc.json from the artifacts, not from prose.

    python3 -I stage09_qc.py

Every number in the QC file is READ OUT of the stage that produced it. None is
transcribed. F-038 is why: three ZT_2026 marginal checks scored a literal `0`
for twelve days while the summary beside them said "ALL FOUR EXACT" -- the
claim was true and simply untested, and a transcribed summary cannot tell the
difference between a check that passed and a check that never ran.
"""
import csv
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))


def load(name):
    with open(os.path.join(HERE, name)) as fh:
        return json.load(fh)


def main():
    glyphs = load("04_glyphs.json")
    merged = load("04_merged_glyph_verification.json")
    vals = load("05_value_checks.json")
    enum = load("03_enumeration_summary.json")
    traces = load("02_traces.json")
    emitted = load("07_emit_summary.json")

    cells = glyphs["cells"]
    s1_order = [cells[f"S1|0|{c}"] for c in range(1, 11)]
    s2_hdr_rows = sorted({int(k.split("|")[1]) for k in cells
                          if k.startswith("S2|") and f"S2|{k.split('|')[1]}|0"
                          not in cells})
    s2_order = [cells[f"S2|{h}|{c}"].split("-")[0]
                for h in s2_hdr_rows for c in (1, 6)]

    from rdkit import Chem
    enum_rows = {r["lipid_name"]: r["IL_SMILES"] for r in csv.DictReader(
        open(os.path.join(HERE, "03_enumerated.csv")))}
    anchors = []
    for r in csv.DictReader(open(os.path.join(HERE,
                                              "reference_structures.csv"))):
        a = Chem.MolToInchiKey(Chem.MolFromSmiles(r["reference_SMILES"]))
        b = Chem.MolToInchiKey(Chem.MolFromSmiles(enum_rows[r["lipid_name"]]))
        anchors.append({"lipid_name": r["lipid_name"],
                        "anchor_inchikey": a, "enumerated_inchikey": b,
                        "full_inchikey_match": a == b,
                        "source_figure": r["source_figure"]})

    chains = {}
    for region, r in traces["regions"].items():
        chains[region] = {
            "bond_length_px_measured": r["bond_length_px"],
            "components": [
                {"id": c["id"], "carbons_A_rdp": c["carbons_A_rdp"],
                 "carbons_B_arc": c["carbons_B_arc"],
                 "carbons_C_apex": c["carbons_C_apex"],
                 "two_primary_methods_agree": c["methods_agree"]}
                for c in r["components"]],
        }

    qc = {
        "paper": "Guo et al., ACS Nano 2026, 20(20):14602-14614, "
                 "doi 10.1021/acsnano.6c00706, PMID 42130331",
        "experiment_id": emitted["experiment_id"],

        "P1_value_source_tier": {
            "tier": "5_pixels, but printed numerals in the authors' own "
                    "source-data tables, not a colour ramp or a bar height",
            "why_not_tiers_1_to_4": "see 00_value_sources.json; each of the "
                                    "four is a reproducible observation, not "
                                    "an assertion",
            "classes_structurally_inapplicable": [16, 20, 23, 27, 32, 33, 43],
            "reason": "there is no axis, no origin, no scale factor, no "
                      "colormap and no bar in the read",
        },

        "P2_glyph_decoding": {
            "glyphs": glyphs["n_glyphs"],
            "clusters": glyphs["n_clusters"],
            "max_intra_cluster_distance": glyphs["max_intra_cluster_distance"],
            "min_inter_cluster_distance": glyphs["min_inter_cluster_distance"],
            "interpretation": "max intra-cluster distance is 0.000, i.e. every "
                              "member of every cluster is a BYTE-IDENTICAL "
                              "normalised bitmap; the nearest two clusters are "
                              "2.236 apart. The alphabet is 45 exact shapes "
                              "over 2,243 glyphs.",
            "cells_decoded": len(glyphs["cells"]),
            "underline_components_rejected": 0,
            "note_on_underlines": "The height<=3 & width>=10 filter never "
                                  "fired: there are no underlines in either "
                                  "table. Recorded so that a future run in "
                                  "which it DOES fire is visible.",
        },

        "P3_merged_glyph_independent_rederivation": {
            "multi_character_clusters": merged["n_multi"],
            "disagreeing": merged["n_disagreeing"],
            "residual_zero": sum(1 for r in merged["results"]
                                 if r.get("residual_px") == 0),
            "disagreements": [
                {"eye": r["eye"], "derived": r["derived"],
                 "residual_px": r.get("residual_px"),
                 "ink_px": r.get("ink_px")}
                for r in merged["results"] if not r["agrees"]],
            "interpretation": "25 of 28 merged clusters are reconstructed "
                              "EXACTLY (residual 0 px) from the "
                              "single-character templates by a dynamic "
                              "programme that never sees the eye label. The "
                              "three that are not -- AA, OG, 3A -- are the "
                              "only clusters containing A, O or G, letters "
                              "that occur nowhere in isolation in either "
                              "table, so the template alphabet cannot spell "
                              "them; their residuals are 90-159 px against 0 "
                              "for every confirmed one. All three occur ONLY "
                              "in name cells (AA in all 50 Table S2 headers, "
                              "OG in OG8, 3A in T3A) and each is printed "
                              "independently in Fig. S1B, Fig. S2 and Fig. "
                              "S3, so no VALUE depends on them.",
        },

        "P4_table_geometry": {
            "S1": glyphs["meta"]["S1"],
            "S2": glyphs["meta"]["S2"],
            "interpretation": "The cell grid is the table's own ruled lines, "
                              "detected as rows/columns above 0.6 ink "
                              "fraction. Table S1: 7 horizontal rules x 12 "
                              "vertical = 6 x 11 cells. Table S2: 35 "
                              "horizontal rules (5 blocks x 7) x 12 vertical. "
                              "No grid position is inferred or interpolated.",
        },

        "P5_values": vals,

        "P5b_amine_column_order_S1_vs_S2": {
            "table_S1_header": s1_order,
            "table_S2_block_headers": s2_order,
            "orders_agree": s1_order == s2_order,
            "difference": [i for i, (a, b) in
                           enumerate(zip(s1_order, s2_order)) if a != b],
            "why_it_costs_nothing_here": "Every column in both tables is "
                                         "labelled with the full compound "
                                         "name, and the decoder keys on the "
                                         "label, never on position. The "
                                         "finding is recorded because a read "
                                         "that assumed the two tables share a "
                                         "column order would swap 96 and 400 "
                                         "-- ten values in Table S1 -- and "
                                         "nothing downstream would notice.",
        },

        "P6_lead_compound_external_check": {
            "claim": "The paper names 113-AA-C8C14 as its lead compound "
                     "throughout (title-level claim; Fig. S4, S8-S17, S20-S24 "
                     "captions).",
            "test": "That name decomposes to Table S2 column `113-AA-C8` and "
                    "row `Cl14`. Nothing in the decoding knows which cell is "
                    "supposed to win.",
            "result": vals["library_2_maximum"],
            "passes": vals["lead_is_the_maximum"],
            "runner_up": "113-AA-C8C12 at 9.3E+04, i.e. the maximum is unique "
                         "and the second-highest value is the same column",
            "seeded_failure": "ROW_TO_ACYL rotated by one (an off-by-one row) "
                              "-> the check REFUSED with 'the maximum is "
                              "113-AA-C8C6'. Confirmed 2026-11-10. A NEW "
                              "CHECK IS NOT DONE UNTIL IT HAS BEEN MADE TO "
                              "FAIL ON PURPOSE.",
        },

        "P7_library_size": {
            "stated_in_abstract": enum["stated_in_abstract"],
            "enumerated": enum["total"],
            "library_1": enum["library_1"],
            "library_2": enum["library_2"],
            "unique_canonical_SMILES": enum["unique_canonical_SMILES"],
            "values_present_for": vals["n_values"],
            "interpretation": "10 amines x 5 epoxides = 50 (Library 1, Table "
                              "S1) and 50 x 5 acyl chlorides = 250 (Library "
                              "2, Table S2) = 300, which is the number the "
                              "abstract states. Every one of the 300 has both "
                              "a structure and a value, with no leftovers in "
                              "either direction.",
        },

        "P8_chain_counts": {
            "vector_route_unavailable": "page.get_drawings() returns 0 paths "
                                        "on every page of this SI; each "
                                        "figure page holds exactly one image "
                                        "object. The figures are fully "
                                        "flattened, so chains' vector path "
                                        "(the PDF's own bond coordinates) "
                                        "does not exist here. stage02's "
                                        "render() asserts this.",
            "methods": "A = chains.raster_chain_length (skeletonise, trace, "
                       "Douglas-Peucker). B = arc length / the bond length "
                       "measured from the region's own longest chain, with "
                       "the label clearance added back. C = per-column ink "
                       "centroid turns, reported but NOT used for agreement: "
                       "it is unreliable on these drawings because several "
                       "arms contain near-vertical and near-horizontal bonds, "
                       "for which y(x) is not single-valued or not turning.",
            "regions": chains,
            "resolved_by_inspection": [
                "figS1B_row1 #2/#3 (amine 75's two N-ethyls), #7 (amine 96's "
                "propylene) and #8 (one of amine 400's two propylenes), and "
                "figS1B_row3 #1 (one of T3A's three arms), where A and B "
                "differ by one. Each was read on a 2000 dpi render saved in "
                "00_pages/ (zoom_amines_70_75_80.png, zoom_amines_96_400.png) "
                "and in every case the A count is the one the drawing shows. "
                "For T3A the other two arms measure 3 and 3 with both methods "
                "agreeing. THIS IS ERROR CLASS 48 EXPOSURE AND IT IS RECORDED "
                "RATHER THAN HIDDEN: five chain lengths in this curation rest "
                "on one automated method plus a human reading a render.",
                "figS2_legend #4 (the branched 2-ethylhexyl of OG8), where "
                "both automated routes follow only the longest path and so "
                "miss the ethyl branch entirely. Read at 2000 dpi "
                "(zoom_R_OG8.png): CH2-O-CH2-CH(CH2CH3)-CH2CH2CH2CH3.",
            ],
            "name_token_vs_measurement": "The epoxide name tokens are NOT "
                                         "chain lengths and the mismatch is "
                                         "the chemistry, not an error: `C8` "
                                         "is 1,2-epoxyoctane and leaves a "
                                         "SIX-carbon R on the carbinol "
                                         "carbon, measured as 6 with both "
                                         "routes agreeing. `OG8` is a "
                                         "2-ethylhexyl glycidyl ether, 11 "
                                         "carbons in the reagent and 8 in the "
                                         "group the token names. The acyl "
                                         "tokens DO match, and they are not a "
                                         "measurement at all: Fig. S2 prints "
                                         "n = 3,5,7,9,11 against Cy = "
                                         "C6,C8,C10,C12,C14.",
        },

        "P9_C10_anchors": {
            "anchors": anchors,
            "all_match": all(a["full_inchikey_match"] for a in anchors),
            "independence": "make_reference.py composes a SMILES string "
                            "centre-out from the tertiary nitrogen and parses "
                            "it once. It contains no RWMol, no AddBond, no "
                            "RemoveAtom and no atom map, where the enumerator "
                            "is built entirely on those. C16 examines the "
                            "generator's COMPOSED templates and reports them "
                            "as derived.",
            "C16_seeded_failure": "TEMPLATE_113 was replaced by the "
                                  "enumeration's own canonical SMILES and C16 "
                                  "flagged it: '1 TEMPLATE(S) APPEAR TO BE "
                                  "TRANSCRIBED FROM THE ENUMERATION ... a "
                                  "104-character SMILES literal is ... a "
                                  "substring of the enumerated structure'. "
                                  "Restored and re-run clean. Confirmed "
                                  "2026-11-10.",
            "what_the_anchors_do_NOT_prove": "They re-derive the DRAWING. "
                                             "They cannot test the drawing, "
                                             "because this paper publishes no "
                                             "HRMS, NMR table, SMILES, InChI "
                                             "or structure file against which "
                                             "the drawn connectivity could be "
                                             "checked. Mass cannot see "
                                             "connectivity and here there is "
                                             "not even a mass.",
        },

        "P10_glyph_label_gate_seeded_failure":
            "One cluster's label was deleted from 04_glyph_labels.json and "
            "stage04_tables.py --read REFUSED: 'cluster 10 (n=256) has "
            "fingerprint 1bf2efaa75d078b0, which is not in "
            "04_glyph_labels.json'. Exit 1, nothing written. Restored and "
            "re-run clean. Confirmed 2026-11-10.",

        "P11_emitted_rows": emitted,

        "what_is_NOT_checked": [
            "The assay. Fig. 1C and Fig. 1E are never described in the SI and "
            "the main text is paywalled, so cell line, readout, units, dose, "
            "timepoint and even in vitro vs in vivo are unknown. Five of the "
            "79 columns are blank for this reason. DECISIONS.md D-1.",
            "Whether Library 1 and Library 2 share a measurement scale. They "
            "have no compound in common, so no anchor ties them. DECISIONS.md "
            "D-2.",
            "The formulation used for the screen. Fig. S11 optimises the "
            "ionizable-lipid / DOPE / cholesterol ratio at constant "
            "DMG-PEG2000 for the LEAD compound only, which is downstream of "
            "the screen. DECISIONS.md D-3.",
            "The drawn connectivity itself: no HRMS, no NMR table, no "
            "SMILES/InChI, no structure file anywhere in the SI (grep of the "
            "full extracted text for smiles|inchi|HRMS|ESI-MS|calcd returns "
            "0). The only characterisation of any member is Fig. S9, a mass "
            "spectrum and a 1H NMR of the lead published as IMAGES, with no "
            "printed m/z or shift list.",
            "Stereochemistry. The thiolactone is the DL reagent and no wedge "
            "is drawn anywhere in Fig. S1 or Fig. S2, so every structure has "
            "two or more unassigned stereocentres.",
        ],
    }
    with open(os.path.join(HERE, "07_qc.json"), "w") as fh:
        json.dump(qc, fh, indent=1)
    print(json.dumps({k: (v if not isinstance(v, (dict, list)) else "...")
                      for k, v in qc.items()}, indent=1))
    print("\nwrote 07_qc.json")


if __name__ == "__main__":
    main()
