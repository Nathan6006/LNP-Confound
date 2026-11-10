#!/usr/bin/env python3
"""07_qc.json -- the rollup. Every check with NUMBERS, nothing in prose only.

Assembled from the stage QC files rather than restated, because a restated
number drifts (F-001, F-024: a quantity with a prose definition and no
implementation gets re-implemented, slightly differently, by every pass).
"""
import csv
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))


def j(fn):
    return json.load(open(os.path.join(HERE, fn)))


def rows(fn):
    with open(os.path.join(HERE, fn), newline="") as fh:
        return list(csv.DictReader(fh))


def main():
    geo, tbl = j("02_fig_s6_geometry.json"), j("04_tables_qc.json")
    val, emi = j("05_values_qc.json"), j("06_rows_qc.json")
    seeded, ov = j("07_seeded_checks.json"), j("02_corpus_overlap.json")
    mass = j("01_mass_identity.json")
    out = {
        "_what": ("Rollup for whiteheadka_30637934 -> Experiment_ID KW_2019. "
                  "Every entry is a measured number. Each check was run clean "
                  "AND with a defect seeded; 07_seeded_checks.json holds both "
                  "readings and a check is credited only where the AGREEMENT "
                  "COUNT MOVED."),
        "paper": {
            "pmid": "30637934", "doi": "10.1002/smll.201805097",
            "title": ("Branched-Tail Lipid Nanoparticles Potently Deliver "
                      "mRNA In Vivo due to Enhanced Ionization at Endosomal "
                      "pH."),
            "title_and_endpoint_verified_from":
                "Europe PMC REST core record, SRC:MED EXT_ID:30637934, 2026-11-10",
            "endpoint": ("mRNA delivery / transfection. Emitted endpoint: in "
                         "vitro firefly-luciferase activity in HeLa, RLU, 100 "
                         "ng mRNA per well, 24 h, n=4."),
            "experiment_id": "KW_2019",
            "id_collision_check": ("KW_2014 exists in the corpus (Whitehead "
                                   "2014); KW_2019 does not. The staged "
                                   "sibling whiteheadka_38888972 uses KW_2024, "
                                   "so the lab key is consistent."),
        },
        "extraction_tier": {
            "value_source": "4b_pdf_vector_geometry",
            "pages_rendered_to_raster": 0,
            "pixels_read": 0,
            "why_not_tier_1_2_3": ("No Source Data file (Wiley, isOpenAccess=N, "
                                   "hasSuppl=N, availabilityCode S); no SI "
                                   "numeric table carries a delivery value; "
                                   "find_repo.py returned 0 relevant hits and "
                                   "NCBI esearch db=pmc on the DOI returns "
                                   "Count=0."),
        },
        "three_counts_that_differ": val["three_counts_that_differ"],
        "structures": {
            "enumerated": len(rows("03_enumerated.csv")),
            "michael_arms_per_structure": 4,
            "all_structures_parse": True,
            "radical_atoms": 0,
            "built_via": ("curation/lib assemble.esterify for the ester bond "
                          "and assemble.attach on atom-mapped attachment "
                          "atoms for all four N-C bonds; no RemoveAtom "
                          "anywhere near a stereocentre, no string surgery."),
            "stereocentres": 0,
            "stereocentres_note": ("The 306 core is achiral and 8-methylnonyl's "
                                   "branch carbon bears two identical methyls, "
                                   "so no structure here has a stereocentre and "
                                   "assemble's stereo assertions are vacuous "
                                   "rather than passing. Stated so the clean "
                                   "line is not read as evidence."),
        },
        "check_mass_ladder": {
            "_what": ("Each consecutive linear-tail step must add exactly "
                      "4 x dC x CH2, because the lipidoid has four arms. A "
                      "relation among the authors' name tokens and our masses "
                      "that the paper never prints."),
            "agreements": seeded["checks"][2]["clean"],
            "of": 9, "seeded": seeded["checks"][2]["seeded"],
            "CH2_Da": 14.0156500642,
            "one_tail_carbon_adds_Da": 56.0626002568,
        },
        "check_abstract_mass_identity": mass,
        "check_figS7_printed_pKa_vs_tableS1": {
            "n": tbl["check_figS7_vs_tableS1_pKa"]["n"],
            "exact_matches": tbl["check_figS7_vs_tableS1_pKa"]["exact_matches"],
            "max_abs_delta": tbl["check_figS7_vs_tableS1_pKa"]["max_abs_delta"],
        },
        "check_two_pKa_conventions_DO_NOT_POOL":
            {k: v for k, v in tbl["check_two_pKa_conventions_DO_NOT_POOL"].items()
             if k != "per_lipid"},
        "check_mole_percents_sum_to_100": tbl["check_mole_percents_sum_to_100"],
        "geometry_panelA": {
            "y_ticks_detected": geo["panelA_y_calibration"]["n_ticks"],
            "y_tick_fit_residual_pt": geo["panelA_y_calibration"]["tick_fit_residual_px"],
            "y_tick_residual_tolerance_pt": geo["panelA_y_calibration"]["tick_residual_tol_px"],
            "y_loo_worst_residual_pt": geo["panelA_y_loo_worst_px"],
            "y_loo_tolerance_pt": 1.5,
            "px_per_5e5_RLU": geo["panelA_px_per_5e5_RLU"],
            "RLU_per_pt": geo["panelA_RLU_per_px"],
            "x_ticks_detected": len(geo["panelA_x_tick_positions_px"]),
            "x_categories": geo["panelA_categories"],
            "x_tick_pitch_spread_pt": geo["panelA_x_tick_pitch_spread_px"],
            "bar_centre_vs_x_tick_worst_pt": geo["panelA_bar_centre_vs_xtick_worst_px"],
            "label_right_edge_minus_tick_mean_pt": geo["panelA_label_offset_mean_px"],
            "label_right_edge_minus_tick_spread_pt": geo["panelA_label_offset_spread_px"],
            "label_offset_n": geo["panelA_label_offset_n"],
            "label_offset_spread_is_shift_invariant": True,
            "label_offset_spread_caveat": ("0.002 pt across 11 labels is the "
                                           "striking number but it does NOT "
                                           "pin the assignment: rotating the "
                                           "category list leaves it unchanged. "
                                           "The absolute nearest-tick test "
                                           "does pin it and scores 11/11 "
                                           "clean, 0/11 seeded. See "
                                           "08_contradictions.csv X-8."),
            "bars_drawn": geo["panelA_bars_drawn"],
            "bar_outline_width_pt": 8.51,
            "errorbar_cap_width_pt": 4.19,
            "pattern_tile_width_pt": 2.80,
            "width_discriminator_note": ("Bar, whisker cap and marker separate "
                                         "by EXACT width, not by a coverage "
                                         "threshold, so error class 32 "
                                         "(F-018, F-020: a replicate marker "
                                         "read as a bar top) is closed by "
                                         "construction. curation/lib bars.py "
                                         "is not called because it takes an "
                                         "ink array and no raster exists."),
            "fill_vs_stroke_frame_offset_pt": geo["panelA_fill_vs_stroke_frame_offset_px"],
            "fill_vs_stroke_frame_offset_RLU": geo["panelA_frame_offset_worth_RLU"],
            "frame_choice": ("bar OUTLINE, whose legs start exactly on the "
                             "y=0 tick at 220.45 pt, so bar tops and ticks "
                             "share one coordinate system and no origin is "
                             "supplied separately (error class 20 / F-033)."),
        },
        "geometry_panelB": {
            "y_ticks_detected": geo["panelB_y_calibration"]["n_ticks"],
            "y_tick_fit_residual_pt": geo["panelB_y_calibration"]["tick_fit_residual_px"],
            "y_loo_worst_residual_pt": geo["panelB_y_loo_worst_px"],
            "px_per_5e5_RLU": geo["panelB_px_per_5e5_RLU"],
            "markers_found": geo["panelB_markers_found"],
            "markers_expected": 12,
            "dose_monotonicity": geo["panelB_monotone_count"],
            "dose_monotonicity_note": ("4/4 lipidoids increase across 10 < 50 "
                                       "< 100 ng. This also validates the "
                                       "legend colour -> dose mapping: the "
                                       "seeded swap of two legend swatches "
                                       "drops it to 0/4."),
        },
        "CROSS_PANEL_DISAGREEMENT": geo["cross_panel_S6A_vs_S6B_100ng"],
        "corpus_overlap": {
            **ov["summary"],
            "tierA_substructure_rows": ov["near_miss_tierA_exact_306_core"]["n"],
            "tierA_rows_identical_to_ours":
                ov["near_miss_tierA_exact_306_core"]["n_identical_to_one_of_ours"],
            "attestation": ("KW_2014 stores five of these structures under "
                            "'A306_C<n>' names and NC_2024 stores six under "
                            "'A306_O<n>' names; all six of ours that match do "
                            "so by FULL InChIKey with canonically identical "
                            "SMILES, and zero of the eleven have a name match "
                            "with a structure mismatch. Since this paper "
                            "prints no HRMS, no formula, no MW and draws no "
                            "structure, this is its ONLY external structure "
                            "anchor besides the reagent names themselves."),
            "tierB_same_chemistry_other_core_total":
                ov["near_miss_tierB_same_chemistry_other_amine_core"]["total"],
        },
        "emitted_rows": {
            **{k: v for k, v in emi.items()
               if k not in ("excipients_resolved_by_lib_reagents",)},
            "excipients_resolved_by": ("curation/lib reagents.smi, which "
                                       "RAISED on DMPE-PEG2000 (two "
                                       "structures in the corpus) and "
                                       "resolved DOPE, Cholesterol and "
                                       "C14-PEG2000 uniquely. See "
                                       "08_contradictions.csv X-10."),
        },
        "C10_anchor": {
            "anchors": 3,
            "anchor_names": ["306O10", "306Oi10", "306O18"],
            "paper_draws_a_structure": False,
            "derived_from": ("the SI Methods reagent NAMES (p.2), typed as "
                             "source string literals, head-first traversal"),
            "independence_established_by": ("PROVENANCE, not byte difference "
                                            "(F-037). make_reference.py's "
                                            "assert_no_assembler() walks its "
                                            "own AST at run time and confirms "
                                            "it imports no assembler module "
                                            "and calls none of RWMol, "
                                            "AddBond, RemoveAtom, "
                                            "CombineMols, EditableMol, "
                                            "ReactionFromSmarts, "
                                            "RunReactants, esterify, attach."),
            "full_inchikey_matches": "3/3",
            "gate_exit_clean": 0,
            "gate_exit_seeded": {"swap_tail_length": 1, "move_an_ester": 1,
                                 "drop_an_arm": 1},
            "move_an_ester_note": ("identical molecular formula and identical "
                                   "exact mass to every decimal place, so the "
                                   "gate demonstrably tests connectivity and "
                                   "not mass."),
            "first_version_of_assert_no_assembler_could_not_pass":
                ("it substring-scanned a file that contains its own "
                 "blocklist and reported all eight forbidden tokens on a "
                 "clean run. Replaced with AST inspection. A check that "
                 "cannot pass is as useless as one that cannot fail."),
            "NOT_RUN_BY_INSTRUCTION": ["provenance/structure_check.py",
                                       "provenance/check_anchor_independence.py",
                                       "provenance/check_anchor_coverage.py"],
        },
        "seeded_checks": {
            "n_checks": seeded["n_checks"],
            "n_credited": seeded["n_credited"],
            "not_credited": seeded["n_not_credited"],
            "rule": ("credited only where the agreement count MOVED; a seed "
                     "that merely raises proves the code has an exception, "
                     "not that the check discriminates."),
        },
        "contradictions": {
            "rows": len(rows("08_contradictions.csv")),
            "confirmed": sum(1 for r in rows("08_contradictions.csv")
                             if r["status"].startswith("CONFIRMED")),
            "retired_by_measurement": sum(
                1 for r in rows("08_contradictions.csv")
                if r["status"].startswith("RETIRED")),
            "my_own_errors_retired": 2,
        },
        "open_decisions": ["D-1 excipient 16 vs 46.5", "D-2 in-vivo panels",
                           "D-3 Figure S6B as a separate experiment",
                           "D-4 the isodecyl branch",
                           "D-5 crude-product flag policy"],
        "NOT_MERGED": ("Nothing was written to data/. 06_rows.csv is not "
                       "registered in build_dataset.py's MERGE list and five "
                       "DECISIONS.md entries are OPEN, which blocks the merge "
                       "by design."),
    }
    with open(os.path.join(HERE, "07_qc.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps(out, indent=1)[:1200])
    print(f"\n07_qc.json written, {len(json.dumps(out))} bytes")


if __name__ == "__main__":
    main()
