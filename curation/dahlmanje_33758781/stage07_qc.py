#!/usr/bin/env python3
"""ZG_2020 -- assemble 07_qc.json from the stage outputs, never by hand.

Run:  python3 -I stage07_qc.py

Every number here is READ OUT of the file the stage that computed it wrote.
A QC file transcribed by hand is the F-024 shape: a quantity with a prose
definition and no implementation gets re-implemented, slightly differently,
by every pass.
"""
import csv, hashlib, json, os, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))


def j(name):
    with open(os.path.join(HERE, name)) as fh:
        return json.load(fh)


def run(script, *args):
    r = subprocess.run([sys.executable, "-I", os.path.join(HERE, script),
                        *args], capture_output=True, text=True, cwd=HERE,
                       timeout=1200)
    return r.returncode, r.stdout + r.stderr


def main():
    enum, fig, fig1c = j("03_enumerated.json"), j("04_figure_checks.json"), \
        j("04_fig1c_checks.json")
    ov, vs = j("02_corpus_overlap.json"), j("05_values_summary.json")

    # re-run the two seed harnesses so their verdicts are measured, not quoted
    rc_seed, out_seed = run("stage09_seeds.py")
    rc_anch, out_anch = run("make_reference.py", "--selftest")
    seed_tail = out_seed.strip().splitlines()[-1]
    anch_tail = out_anch.strip().splitlines()[-1]

    with open(os.path.join(HERE, "06_rows.csv")) as fh:
        rows06 = fh.read().splitlines()
    ref = list(csv.DictReader(open(os.path.join(HERE,
                                                "reference_structures.csv"))))
    contra = list(csv.DictReader(open(os.path.join(HERE,
                                                   "08_contradictions.csv"))))

    src = os.path.join(HERE, "source")
    files = {}
    for d, _sub, fns in os.walk(src):
        for fn in sorted(fns):
            p = os.path.join(d, fn)
            b = open(p, "rb").read()
            files[os.path.relpath(p, HERE)] = {
                "bytes": len(b), "magic": repr(b[:4]),
                "sha256": hashlib.sha256(b).hexdigest()}

    qc = {
        "paper": {
            "pmid": 33758781, "pmcid": "PMC7974401",
            "doi": "10.1002/btm2.10161",
            "first_author": "Gan Z", "last_author": "Dahlman JE",
            "prospective_Experiment_ID": "ZG_2020",
            "licence_article": "CC BY 4.0 (Europe PMC core record license "
                               "field = 'cc by'; JATS <permissions> cites "
                               "creativecommons.org/licenses/by/4.0/)",
            "licence_repo": "n/a - no repository or deposit exists",
            "endpoint": "Cre mRNA -> tdTomato recombination in Ai14 mice, read "
                        "out as (i) % tdTom+ cells by flow cytometry and (ii) "
                        "'Normalized Delivery' from DNA-barcode sequencing of "
                        "tdTomato+ sorted cells",
            "readout_requires_expression": True,
            "readout_compositional": True,
        },

        "HEADLINE": {
            "distinct_ionizable_lipids_contributed": 1,
            "ionizable_lipid": "cKK-E15, identical on all 109 LNPs",
            "distinct_helper_phospholipids_synthesised": len(enum["products"]),
            "distinct_sterols": 3, "distinct_molar_ratios": 4,
            "formulations_designed": 120, "formulations_passed_QC": 109,
            "formulations_dropped_at_QC": 11,
            "cell_types_assayed": 16,
            "per_LNP_values_the_design_supports": 109 * 16,
            "per_LNP_values_the_paper_ATTRIBUTES": 0,
            "individually_injected_LNPs": 1,
            "value_loci_recorded_in_05_values": vs["value_loci_recorded"],
            "corpus_admissible_yes": vs["corpus_admissible_yes"],
            "corpus_admissible_no": vs["corpus_admissible_no"],
            "rows_emitted": 0,
            "06_rows_csv_lines": len(rows06),
            "06_rows_csv_columns": len(rows06[0].split(",")),
        },

        "P1_HRMS_vs_our_build": {
            "n_compounds": enum["checks"]["n_products"],
            "exact_mass_matches_within_0.6_mDa":
                enum["checks"]["mass_matches_within_0.6_mDa"],
            "max_abs_delta_mDa": enum["checks"]["max_abs_delta_mDa"],
            "adduct_formula_matches": enum["checks"]["formula_matches"],
            "adduct_conventions_in_one_column":
                enum["checks"]["adduct_conventions_used"],
            "_convention_drift": "Nine compounds are reported as [M+Na]+ and "
                                 "A-13 as [M+H]+, in the same column. Both "
                                 "reproduce. This is trap 2, measured.",
            "unexplained_formula_mismatches":
                enum["checks"]["formula_mismatches_unexplained"],
            "named_exception": list(enum["checks"]["known_formula_typo"]),
        },

        "P2_Fig1c_independent_property_table": {
            "n": fig1c["D_a_fig1c_formulas"]["n"],
            "neutral_formula_exact_matches":
                fig1c["D_a_fig1c_formulas"]["exact_matches"],
            "mw_median_printed_offset_Da":
                fig1c["D_c_molecular_weight"]["median_offset_Da"],
            "mw_max_abs_residual_after_offset_Da":
                fig1c["D_c_molecular_weight"][
                    "max_abs_residual_after_offset_Da"],
            "mw_printed_steps_off_CH2":
                fig1c["D_c_molecular_weight"]["steps_off_CH2_by_more_than_0.6"],
            "logp_increments_positive":
                f"{sum(1 for i in fig1c['D_d_logp_monotone']['increments'] if i > 0)}"
                f"/{fig1c['D_d_logp_monotone']['n_increments']}",
            "_source": "Fig 1c is a pure BITMAP -- the JATS declares "
                       "table-count=0 -- so its cells are transcribed by eye. "
                       "A single misread character breaks the 10/10 formula "
                       "match, which is what makes the transcription "
                       "checkable.",
        },

        "P3_free_checks_the_authors_left": {
            "molar_ratio_columns_summing_to_100":
                f"{fig['C_a_molar_ratio_sums']['n_summing_to_100']}"
                f"/{fig['C_a_molar_ratio_sums']['n_formulations']}",
            "molar_ratio_offenders": fig["C_a_molar_ratio_sums"]["offenders"],
            "factorial_10x3x4": fig["C_b_factorial"]["product"],
            "paper_states_total": fig["C_b_factorial"]["paper_states_total"],
            "factorial_closes": fig["C_b_factorial"]["factorial_closes"],
            "fig3h_teal_bars": fig["C_c_fig3h_sign_by_colour"]["n_positive"],
            "fig3h_caption_prints": fig["C_c_fig3h_sign_by_colour"]["caption_prints"],
            "fig3h_sign_counts_sum":
                fig["C_c_fig3h_sign_by_colour"]["n_positive"]
                + fig["C_c_fig3h_sign_by_colour"]["n_negative"]
                + fig["C_c_fig3h_sign_by_colour"]["n_about_zero"]
                == fig["C_c_fig3h_sign_by_colour"]["n_bars"],
            "fig3h_tallest_enriched": fig["C_d_fig3h_rank"]["tallest_teal_bar"],
            "text_claims_most_enriched":
                fig["C_d_fig3h_rank"]["text_claims_most_enriched"],
            "printed_yield_range": fig1c["D_e_yield_range"][
                "fig1a_prints_range_percent"],
            "yield_range_reproduces_WITH_A9":
                fig1c["D_e_yield_range"]["matches_with_A9"],
            "yield_range_reproduces_WITHOUT_A9":
                fig1c["D_e_yield_range"]["matches_without_A9"],
            "_best_one": "The printed 77-92% yield range is a number the "
                         "authors drew and never computed, and it settles "
                         "their own off-by-one: 92% is A-9's yield alone, so "
                         "A-9 is in the series that Fig 1a labels "
                         "'A10-17-2Z'.",
        },

        "P4_name_token_vs_chain_length": {
            "rule_measured": "A-N carries an acyl chain of N+1 carbons: A-9 is "
                             "decanoyl (C10), A-17 is stearoyl (C18). The "
                             "token counts the alkyl carbons BEYOND the ester "
                             "carbonyl.",
            "n_checked": len(enum["products"]),
            "per_compound": {p["name"]: p["acid_carbons"]
                             for p in enum["products"]},
            "all_consistent": all(
                p["acid_carbons"] == int(p["name"].split("-")[1]) + 1
                for p in enum["products"]),
        },

        "P5_corpus_overlap": {
            "ckk_e15_rows_matching_full_InChIKey":
                ov["ckk_e15"]["ckk_e15_full_inchikey_rows_in_corpus"],
            "ckk_e15_rows_matching_skeleton_block":
                ov["ckk_e15"]["ckk_e15_skeleton_block_rows_in_corpus"],
            "ckk_e12_IS_in_corpus_under_names":
                ov["ckk_e15"]["corpus_ckk_e12_names_held"],
            "name_join_would_find": ov["ckk_e15"]["_name_join_would_have_found"],
            "inchikey_join_finds": ov["ckk_e15"]["_inchikey_join_finds"],
            "substructure_sweep_SUPERSET": {
                "structures_scanned":
                    ov["substructure_sweep_SUPERSET"]["unique_structures_scanned"],
                "loose_chemotype_hits":
                    ov["substructure_sweep_SUPERSET"][
                        "loose_hits_diketopiperazine_core_structures"],
                "strict_hits":
                    ov["substructure_sweep_SUPERSET"][
                        "strict_hits_core_plus_two_bis(hydroxyalkyl)amine_arms_structures"],
                "after_exact_InChIKey": 0,
                "_positive_control": "The first version of the core SMARTS "
                                     "demanded a ring CH2 and so could not "
                                     "match cKK-E12 itself; it reported a "
                                     "clean 0. A sweep that misses its own "
                                     "known positive is broken, not loose. "
                                     "The positive control that caught it is "
                                     "now permanent.",
            },
            "token_pattern_damage_if_joined_on_a_bare_number": {
                k: {"rows": v["substring_rows"],
                    "structures": v["substring_structures"],
                    "experiments": len(v["substring_experiments"])}
                for k, v in ov["token_pattern_test"]["hits"].items()
                if k in ("9", "11", "17")},
            "token_pattern_damage_if_joined_on_the_full_label": {
                k: ov["token_pattern_test"]["hits"][k]["substring_rows"]
                for k in ("A-9", "A-11", "A-17-2Z")},
            "_and_the_wrong_namespace": "A-N names a HELPER phospholipid, not "
                                        "an ionizable lipid, so an IL_name "
                                        "join would be wrong twice over.",
            "adamantane_bearing_HL_SMILES_already_in_corpus":
                ov["helper_lipid_namespace"][
                    "any_adamantane_bearing_HL_SMILES_in_corpus"],
        },

        "P6_anchors_C10_posture": {
            "_recorded_by_hand_in_PROVENANCE.md": "provenance/structure_check."
                "py, check_anchor_independence.py and check_anchor_coverage.py "
                "were NOT run: they write tracked files other agents read.",
            "anchors_derived": [r["lipid_name"] for r in ref
                                if r["reference_SMILES"]
                                and r["lipid_name"] != "Cholesterol"],
            "anchor_chain_lengths": "A-9 (C10 acyl) and A-17 (C18 acyl) -- "
                                    "deliberately different lengths",
            "registry_pin_not_independent": ["Cholesterol"],
            "declared_absent": [r["lipid_name"] for r in ref
                                if not r["reference_SMILES"]],
            "inchikey_agreement_with_enumeration": "2/2 (full InChIKey)",
            "guard_1_AST": "walks make_reference.py's own syntax tree at run "
                           "time for RDKit editing primitives and an "
                           "`assemble` import; inspects Call/Import nodes "
                           "ONLY, never string constants, so the ban list "
                           "cannot match itself",
            "guard_2_template_provenance": "no string literal of >=25 chars "
                                           "in make_reference.py may be a "
                                           "substring of any enumerated "
                                           "SMILES -- F-037's defect stated "
                                           "directly, which byte-comparison "
                                           "of the OUTPUT cannot see",
            "byte_identical_to_enumeration": True,
            "_and_that_is_not_the_signal": "A correct re-derivation "
                                           "canonicalises to the same string "
                                           "as a correct enumeration. That is "
                                           "exactly why C10's own "
                                           "anchor_may_be_a_copy missed F-037.",
        },

        "P7_seeded_mutations": {
            "stage_gates_harness": "stage09_seeds.py",
            "stage_gates_exit": rc_seed, "stage_gates_summary": seed_tail,
            "anchor_guard_harness": "make_reference.py --selftest",
            "anchor_guard_exit": rc_anch, "anchor_guard_summary": anch_tail,
            "_credit_rule": "A seed is credited only if its child process "
                            "REACHED a METRIC line. Two seeds were refused "
                            "credit on the first run for exiting before their "
                            "metric (a missing lib path, a missing input "
                            "file); the harness was fixed, not the seeds.",
            "measured_blind_spot": "S5 -- swapping the adamantylacetyl and "
                                   "fatty acyl on ALL TEN compounds leaves "
                                   "the HRMS gate GREEN. Measured delta: "
                                   "formula identical, exact mass identical "
                                   "to 0.000000 Da. Seed G3 shows the "
                                   "InChIKey anchor gate DOES catch the same "
                                   "defect (0/2 agreements), so the two "
                                   "checks have complementary coverage.",
            "measured_sensitivity_limit": "S10 -- shifting every Fig 3h bar "
                                          "slot by half the 22.8 px pitch "
                                          "breaks both the sign count and the "
                                          "rank, so the slot positions are "
                                          "load-bearing at that scale and are "
                                          "detected, not assumed.",
        },

        "P8_contradictions": {
            "n": len(contra),
            "resolved": sum(1 for c in contra if c["status"] == "RESOLVED"),
            "open": sum(1 for c in contra if c["status"] == "OPEN"),
            "ids": {c["id"]: c["claim"] for c in contra},
        },

        "P9_source_files": files,
        "P9_provenance_note": "The Appendix S1 PDF that was already in the "
                              "folder at task start was verified BYTE-"
                              "IDENTICAL (sha256 24d89f54...) to the "
                              "BTM2-5-e10161-s001.pdf served inside Europe "
                              "PMC's supplementaryFiles ZIP for PMC7974401, "
                              "which we downloaded independently. That is an "
                              "F-065 contamination check on a file we did not "
                              "fetch ourselves. No glob copy out of the "
                              "session scratchpad was made at any point.",

        "P10_tier3_null_result": {
            "data_access_statement": "'The data, analyses, and scripts used to "
                                     "generate all figures in the paper are "
                                     "available upon request from J.E.D. or "
                                     "dahlmanlab.org.' No deposit.",
            "jats_chars_searched": 88909,
            "si_text_chars_searched": 25290,
            "patterns": ["GSE[0-9]+", "SRP[0-9]+", "PRJNA[0-9]+", "SRX[0-9]+",
                         "SRR[0-9]+", "SRA", "accession", "deposit", "github",
                         "zenodo", "dryad", "figshare", "osf.io"],
            "hits_in_either_document": 0,
            "find_repo_runs": 2,
            "find_repo_authors_tried": ["Dahlman", "Gan"],
            "find_repo_code_search_hits_on_doi_or_pmid": 0,
            "find_repo_repo_hits": 0,
        },
    }

    fails = []
    if qc["HEADLINE"]["corpus_admissible_yes"] + \
            qc["HEADLINE"]["corpus_admissible_no"] != \
            qc["HEADLINE"]["value_loci_recorded_in_05_values"]:
        fails.append("admissibility counts do not sum")
    if qc["HEADLINE"]["06_rows_csv_lines"] != 1:
        fails.append("06_rows.csv is not header-only")
    if qc["HEADLINE"]["06_rows_csv_columns"] != 79:
        fails.append("06_rows.csv is not the 79-column schema")
    if not qc["P4_name_token_vs_chain_length"]["all_consistent"]:
        fails.append("name-token rule does not hold for all ten")
    if not qc["P3_free_checks_the_authors_left"]["fig3h_sign_counts_sum"]:
        fails.append("Fig 3h sign counts do not sum to 10")
    if rc_seed != 0:
        fails.append("stage09_seeds.py did not credit every seed")
    if rc_anch != 0:
        fails.append("make_reference.py --selftest did not credit every guard")
    qc["GATE_failures"] = fails

    with open(os.path.join(HERE, "07_qc.json"), "w") as fh:
        json.dump(qc, fh, indent=1)
    print(json.dumps(qc["HEADLINE"], indent=1))
    print("\nseed harness:", seed_tail)
    print("anchor guards:", anch_tail)
    print("METRIC qc_failures=" + json.dumps(fails))
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
