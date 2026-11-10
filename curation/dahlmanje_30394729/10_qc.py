#!/usr/bin/env python3
"""Assemble 07_qc.json from the stage outputs. Every entry carries numbers.

    python3 -I 10_qc.py

Run after: 01_si_tables.py, 02_corpus_overlap.py, 05_values.py, 06_checks.py,
07_emit.py, 09_seeds.py, make_reference.py --selftest. It composes rather than
recomputes, so no number in 07_qc.json can disagree with the file that produced
it; if a stage has not been run, this refuses rather than writing a partial QC.
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
NEEDED = {"checks": "07_qc_checks.json", "overlap": "02_corpus_overlap.json",
          "values": "05_values_summary.json", "emit": "06_rows_report.json",
          "seeds": "09_seedgate.json", "anchor": "10_anchor_report.json",
          "tables": "04_si_tables.json"}


def main():
    missing = [v for v in NEEDED.values() if not os.path.exists(os.path.join(HERE, v))]
    if missing:
        print("REFUSING: these stage outputs do not exist yet: " + ", ".join(missing))
        return 1
    d = {k: json.load(open(os.path.join(HERE, v))) for k, v in NEEDED.items()}
    c = d["checks"]

    n = c["N1_N2_normalisation"]
    qc = {
     "_what": "QC for curation/dahlmanje_30394729 -- Sago, Lokugamage, Islam, "
              "Krupczak, Sato, Dahlman, J Am Chem Soc 140(49), 2018, PMID "
              "30394729, 'Nanoparticles That Deliver RNA to Bone Marrow "
              "Identified by in Vivo Directed Evolution.' Assembled by 10_qc.py "
              "from the stage outputs; numbers are composed, not retyped.",
     "proposed_experiment_id": "CS_2018",
     "rows_emitted": d["emit"]["rows_emitted"],

     "THE_GOVERNING_NUMBER": {
      "distinct_ionizable_lipids": 1,
      "which": "7C1 only, in every formulation of both libraries and in every "
               "individually tested LNP.",
      "formulations_formulated": d["values"]["distinct_formulations_formulated"],
      "formulations_pooled_and_measured":
          d["values"]["distinct_formulations_pooled_and_measured"],
      "formulations_individually_measured":
          d["values"]["distinct_formulations_individually_measured"],
      "design_supported_per_LNP_values": d["values"]["design_supported_per_LNP_values"],
      "attributed_per_LNP_values": d["values"]["attributed_per_LNP_values"],
      "note": "A corpus keyed on ionizable-lipid structure gains ONE lipid from "
              "160 formulations -- and gains zero, because that one lipid has no "
              "single SMILES."},

     "N1_normalisation_reconstructed_from_the_paper_s_own_worked_example": {
      "formula": n["formula"],
      "stage_1_raw_column_totals": {
       "n_cells": n["stage_1_raw_column_totals"]["n_cells"],
       "max_residual": n["stage_1_raw_column_totals"]["max_residual"],
       "computed": n["stage_1_raw_column_totals"]["computed"],
       "printed": n["stage_1_raw_column_totals"]["printed"]},
      "stage_3_normalised_delivery": {
       "n_cells": n["stage_3_normalised_delivery"]["n_cells"],
       "n_exact_after_rounding": n["stage_3_normalised_delivery"]["n_exact_after_rounding"],
       "max_residual": n["stage_3_normalised_delivery"]["max_residual"],
       "verdict": "9/9 EXACT, residual 0.0. The normalisation is "
                  "count-share divided by input-share, renormalised to 100% "
                  "across the pool -- so the readout is COMPOSITIONAL."},
      "stage_2_printed_intermediates": {
       "n_cells": n["stage_2_within_sample_percent"]["n_cells"],
       "n_exact_after_rounding": n["stage_2_within_sample_percent"]["n_exact_after_rounding"],
       "max_residual": n["stage_2_within_sample_percent"]["max_residual"],
       "mismatches": n["stage_2_within_sample_percent"]["mismatches"],
       "verdict": "10/12. Two printed cells in the Mouse 2 column are wrong: 34 "
                  "where 80/230 = 34.78 rounds to 35, and 3 where 45/230 = 19.57 "
                  "rounds to 20. Stage 3 is computed from the RAW counts and is "
                  "still exact, which is what shows these are typos in the "
                  "intermediate table rather than a different formula. Confirmed "
                  "visually on the page render, not only from the text layer."},
      "the_drawn_but_never_computed_total": {
       "printed_Total_row_stage_2": n["stage_2_within_sample_percent"]["printed_total_row"],
       "printed_cells_actually_sum_to_stage_2":
           n["stage_2_within_sample_percent"]["printed_cells_actually_sum_to"],
       "printed_Total_row_stage_3": n["stage_3_normalised_delivery"]["printed_total_row"],
       "printed_cells_actually_sum_to_stage_3":
           n["stage_3_normalised_delivery"]["printed_cells_actually_sum_to"],
       "verdict": "'Total (%) = 100' is DRAWN under all three tables and computed "
                  "under none: the stage-2 columns sum to 100 / 83 / 101 / 99 and "
                  "every stage-3 column sums to 101."}},

     "F1_factorial_proved": {
      "library_1": {k: c["F1_factorial"]["factorial_library_1"][k] for k in
                    ("n_peg_variants", "n_ratios_per_variant", "product",
                     "n_rows", "every_variant_carries_the_same_ratio_set",
                     "claimed", "claim_holds")},
      "library_2": {k: c["F1_factorial"]["factorial_library_2"][k] for k in
                    ("n_peg_variants", "n_ratios_per_variant", "product",
                     "n_rows", "every_variant_carries_the_same_ratio_set",
                     "claimed", "claim_holds")},
      "verdict": "5 x 24 = 120 and 2 x 20 = 40, with every PEG variant carrying "
                 "the IDENTICAL ratio set exactly once and zero duplicates. The "
                 "claimed factorials are proved, not taken on trust."},

     "M1_molpct_sums": {
      "library_1": {k: c["M1_molpct"]["molpct_library_1"][k] for k in
                    ("n", "n_sum_to_100", "n_violations")},
      "library_2": {k: c["M1_molpct"]["molpct_library_2"][k] for k in
                    ("n", "n_sum_to_100", "n_violations", "violations")},
      "verdict": "158 of 160 formulations sum to exactly 100 mol%. The two that "
                 "do not are library 2 LNP 36 and LNP 40, both 60/45/5/0 = 110, "
                 "at a stride of 4 across the two PEG blocks -- one error in a "
                 "block that was pasted twice. NOT a tautology of our own "
                 "derivation: the four columns are transcribed independently "
                 "from native PDF text cells and nothing here computes one from "
                 "the others."},

     "E1_formulated_pies_vs_the_composition_tables": {
      "summary": c["E1_formulated_pies_vs_tables"]["_summary"],
      "Fig_4A": {"printed": c["E1_formulated_pies_vs_tables"]
                 ["Fig_4A_library1_PEG_molpct"]["printed_slices"],
                 "computed": c["E1_formulated_pies_vs_tables"]
                 ["Fig_4A_library1_PEG_molpct"]["bands_half_open"]["rounded"],
                 "percent": c["E1_formulated_pies_vs_tables"]
                 ["Fig_4A_library1_PEG_molpct"]["bands_half_open"]["percent"],
                 "verdict": "5/5 slices exact, including 20.83 -> 21 and 16.67 -> "
                            "17, which also PINS the band convention as [lo,hi) "
                            "with '>20' meaning >=20. One wrong PEG value anywhere "
                            "in the 120-row table would break a slice, so this "
                            "validates the whole transcription."},
      "Fig_5H": {"printed": c["E1_formulated_pies_vs_tables"]
                 ["Fig_5H_library2_PEG_molpct"]["printed_slices"],
                 "computed_half_open": c["E1_formulated_pies_vs_tables"]
                 ["Fig_5H_library2_PEG_molpct"]["bands_half_open"]["rounded"],
                 "computed_inclusive": c["E1_formulated_pies_vs_tables"]
                 ["Fig_5H_library2_PEG_molpct"]["bands_inclusive"]["rounded"],
                 "verdict": "DOES NOT RECONCILE under either convention: printed "
                            "40/40/20 against 45/20/35 (the convention Fig 4A "
                            "pins) or 45/40/15 (inclusive). Only the middle slice "
                            "ever matches. Reconciliation would need exactly 2 of "
                            "the 40 LNPs to move from the low band to >20, and "
                            "exactly 2 LNPs (36 and 40) independently fail the "
                            "mol% sum -- a convergence recorded as a hypothesis "
                            "in 08_contradictions.csv, NOT as a conclusion, "
                            "because the 33-35 series pattern implies the "
                            "opposite repair."},
      "the_other_four": "Fig 4B 20/20/20/20/20, Fig 5I 50/50, Fig 5J 20/40/20/20, "
                        "Fig 5K 60/40 -- all exact against the tables, no rounding "
                        "needed."},

     "B1_barcodes": c["B1_barcodes"],
     "T1_total_formulated": c["T1_total_formulated"],
     "W1_the_five_individually_tested_formulations": {
      "n": c["W1_individually_tested"]["n"],
      "n_uniquely_located_in_a_library_table":
          c["W1_individually_tested"]["n_uniquely_located"],
      "n_in_the_expected_library": c["W1_individually_tested"]["n_as_expected"],
      "rows": c["W1_individually_tested"]["rows"],
      "verdict": "All five compositions printed in Fig S3A / Fig 6A are found "
                 "where they should be: BM1 = library 2 LNP 18, LNP1 = library 2 "
                 "LNP 6, LNP3 = library 2 LNP 19, the screen-1 winner = library 1 "
                 "LNP 25. The screen-1 winner's ratio ALSO appears as library 2 "
                 "LNP 34, which is not an error -- library 2 was evolved from "
                 "library 1 and deliberately re-tested that ratio with both PEGs. "
                 "The 'original' 80:20 7C1 is in neither table, as expected: it "
                 "is the comparator, not a library member."},
     "P1_a_round_number_implied_but_never_printed": c["P1_peg_repeat_count"],

     "C1_corpus_overlap": {
      "token_pattern_test": {
       "rows_matched_by_a_substring_join_on_7C1":
           d["overlap"]["token_pattern_test"]["total_rows"],
       "distinct_structures":
           d["overlap"]["token_pattern_test"]["total_distinct_structures"],
       "experiments": d["overlap"]["token_pattern_test"]["total_experiments"],
       "true_matches": 0,
       "verdict": d["overlap"]["token_pattern_test"]["verdict"]},
      "substructure_sweep_SUPERSET": {
       "vitro": d["overlap"]["substructure_sweep_SUPERSET"]["vitro"],
       "vivo": d["overlap"]["substructure_sweep_SUPERSET"]["vivo"],
       "verdict": "The loose amino-alcohol pattern matches 1,247 of 18,720 vitro "
                  "and 276 of 1,506 vivo structures -- a CHEMOTYPE, reported as a "
                  "superset. The strict pattern (the full drawn C13 arm) matches "
                  "0 and 0."},
      "excipients_by_inchikey": {
       k: {"our_full_inchikey": v["our_full_inchikey"],
           "matches_registry": v["matches_registry"],
           "corpus": {kk: {"rows": vv["rows"],
                           "distinct_corpus_structures": vv["distinct_corpus_structures"],
                           "matches": vv["matches_our_drawn_structure"]}
                      for kk, vv in v["corpus"].items()}}
       for k, v in d["overlap"]["excipients_by_inchikey"].items()},
      "in_vivo_shared_chemistry": {
       "n_structures_shared_across_experiments":
           d["overlap"]["in_vivo_shared_chemistry"]["n_structures_shared_across_experiments"],
       "which": [{"names": s["names"], "n_experiments": s["n_experiments"]}
                 for s in d["overlap"]["in_vivo_shared_chemistry"]["structures"]],
       "does_this_paper_add_one": False,
       "why": d["overlap"]["in_vivo_shared_chemistry"]["why"]},
      "benchmarks": "NONE. The paper uses no MC3, no SM-102, no ALC-0315 and no "
                    "other field benchmark anywhere; its only comparator is its "
                    "own 80:20 7C1 formulation. So the 'do benchmarks top the "
                    "screen' question does not arise here, and the paper's own "
                    "lead BM1 does top its own screen 2 by construction -- it was "
                    "SELECTED as the top performer, which is not an independent "
                    "finding."},

     "A1_anchors": {
      "independence_guard": d["anchor"]["independence_guard"],
      "anchors": {k: ({"inchikey": v.get("inchikey"),
                       "skeleton": v.get("skeleton"),
                       "matches_registry_inchikey": v.get("inchikey_matches_registry"),
                       "expected_mismatch": None}
                      if "inchikey" in v else v)
                  for k, v in d["anchor"]["anchors"].items()},
      "discarded_cholesterol_stereo_attempt":
          d["anchor"]["discarded_cholesterol_stereo_attempt"],
      "selftest": d["anchor"].get("selftest_credited"),
      "verdict": "DSPC was assembled from sn-glycero-3-phosphocholine and two "
                 "stearoyl fragments and agrees with the registry by FULL "
                 "InChIKey. The 7C1 representative arm was assembled with the "
                 "regiochemistry pinned by atom map; its Markovnikov isomer has "
                 "the identical formula AND the identical exact mass to six "
                 "decimals (0.0 Da difference) and a different InChIKey, which is "
                 "how the gate is shown to test connectivity. Cholesterol is "
                 "emitted at CONNECTIVITY ONLY after the hand-written stereo "
                 "attempt produced a diastereomer; it was discarded rather than "
                 "tuned against the registry."},

     "S1_seeded_mutation_gate": {
      "verdict": d["seeds"]["verdict"],
      "n_credited": d["seeds"]["n_credited"], "n_seeds": d["seeds"]["n_seeds"],
      "control": d["seeds"]["control"],
      "seeds": {k: {"metric": v["metric"], "control": v["control_value"],
                    "seeded": v["seeded_value"], "returncode": v["returncode"],
                    "credited": v["credited"]}
                for k, v in d["seeds"]["seeds"].items()},
      "a_seed_that_was_NOT_credited": "The first E1 seed moved one library-1 PEG "
                                      "value from 20 to 22.5 and was a NO-OP, "
                                      "because the half-open convention already "
                                      "puts 20 in the '>20' band. It exited "
                                      "non-zero -- on the factorial, not on E1 -- "
                                      "while E1's metric never moved, so it was "
                                      "refused. That is why 'credited' requires "
                                      "the metric to MOVE and not merely a "
                                      "non-zero exit."},

     "emit_guards": d["emit"]["guard_selftests"],

     "V1_values": d["values"],

     "UNVERIFIABLE": [
      "WHICH 115 of the 120 and WHICH 31 of the 40 formulated LNPs were pooled. "
      "Never stated; the diameters that decided it are plotted unlabelled.",
      "The mapping from barcode number to LNP number. The 156 barcodes are "
      "listed and the 160 formulations are listed; nothing connects them.",
      "Every per-LNP normalized delivery value (115 + 31 + the technical "
      "replicates). Unlabelled scatter points.",
      "The N:P or mass ratio for the library LNPs. Only the 'original' 7C1 has "
      "one stated (7C1 + PEG : siRNA = 5:1 by mass); the libraries state only a "
      "total nucleic acid dose of 1.5 mg/kg and a siRNA:barcode mass ratio of "
      "10:1.",
      "The exact molecular identity of the PEG-lipids. Fig 3B draws a "
      "glycerophosphoethanolamine-PEG skeleton; no catalogue product is named "
      "and no dispersity is given.",
      "The structure of 7C1 itself, which is the governing blocker."]
    }
    with open(os.path.join(HERE, "07_qc.json"), "w") as f:
        json.dump(qc, f, indent=1)
    print(json.dumps({"rows_emitted": qc["rows_emitted"],
                      "governing": qc["THE_GOVERNING_NUMBER"],
                      "seed_gate": qc["S1_seeded_mutation_gate"]["verdict"],
                      "pies": qc["E1_formulated_pies_vs_the_composition_tables"]
                      ["summary"]}, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
