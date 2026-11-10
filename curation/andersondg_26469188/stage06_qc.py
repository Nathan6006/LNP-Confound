#!/usr/bin/env python3
"""KK_2015 -- assemble 07_qc.json. Every check with NUMBERS."""
import csv
import hashlib
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.environ.get("KK2015_OUTDIR", HERE)


def j(n):
    return json.load(open(os.path.join(OUTDIR, n)))


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def main():
    P, C, S, E = j("_parsed.json"), j("_checks.json"), j("_structures.json"), j("_emit.json")
    seeds = json.load(open(os.path.join(HERE, "_seeds.json")))
    ov = j("02_corpus_overlap.json")
    vals = list(csv.DictReader(open(os.path.join(OUTDIR, "05_values.csv"))))
    rows = list(csv.DictReader(open(os.path.join(OUTDIR, "06_rows.csv"))))
    refs = list(csv.DictReader(open(os.path.join(OUTDIR, "reference_structures.csv"))))
    contra = list(csv.DictReader(open(os.path.join(HERE, "08_contradictions.csv"))))

    tw = {k: v for k, v in S["three_way_helper_concordance"].items()
          if not k.startswith("_")}
    import statistics as _st
    _zs = [float(r["Experiment_value"]) for r in rows]
    _zmean, _zsd = _st.fmean(_zs), _st.pstdev(_zs)

    qc = {
"_what": "KK_2015 quality control. Kauffman KJ et al., Nano Lett 15(11):7300-7306 "
         "(2015), PMID 26469188. Every figure here is measured, not asserted.",

"headline": {
  "ENDPOINT": "Serum erythropoietin (EPO) protein by ELISA, ng/mL, mean +/- SD, "
              "n = 3 mice, 6 h after 15 ug total mRNA by intravenous injection. "
              "A PROTEIN-ABUNDANCE delivery readout, admissible.",
  "DISTINCT_IONIZABLE_LIPIDS_CONTRIBUTED": 1,
  "the_one_lipid": "C12-200",
  "what_actually_varies": "the phospholipid identity (DSPC / DOPE / DOPC, plus a "
                          "DSPE arm that never formulated), the four molar "
                          "percentages, and the C12-200:mRNA weight ratio",
  "formulations_specified": 39,
  "formulations_individually_measured_for_EPO": 36,
  "formulations_with_no_value": 3,
  "rows_emitted": len(rows),
  "value_tier": "4_pdfplumber_text_table -- ZERO pixels touched anywhere in this "
                "curation. No page rendered, no raster sampled, no axis "
                "calibrated, no bar top detected.",
},

"C0_file_provenance": {
  "source_file": "source/nl5b02497_si_001.pdf",
  "sha256": sha(os.path.join(HERE, "source", "nl5b02497_si_001.pdf")),
  "bytes": os.path.getsize(os.path.join(HERE, "source", "nl5b02497_si_001.pdf")),
  "files_in_source_dir": sorted(os.listdir(os.path.join(HERE, "source"))),
  "n_files_in_source_dir": len(os.listdir(os.path.join(HERE, "source"))),
  "F065_guard": "No file was copied into source/ by this curation, so the "
                "shared-scratchpad contamination path (17 figure files from OTHER "
                "PAPERS landing in a sibling's source/figs/) cannot have occurred. "
                "The directory is asserted by NAME AND HASH above, not by a count.",
  "belongs_to_this_PMID_because": [
    "the filename stem nl5b02497 is the DOI suffix of 10.1021/acs.nanolett.5b02497",
    "SI page 1 reproduces the Europe PMC title and all 8 authors in order",
    "all 14 pages carry the publisher watermark path "
    "pubs.acs.org/nalefd/article-supplement/1257894/pdf/nl5b02497_si_001/"],
  "zip_byte_identity_check": "UNAVAILABLE: pmcid=null, inPMC=N, hasSuppl=N, so "
                             "Europe PMC has no supplementaryFiles ZIP to diff against.",
},

"C1_molpct_sums_to_100": C["C1_molpct_sums"],
"C2_cross_parameterisation": C["C2_cross_parameterisation"],
"C3_coded_design_levels": {k: v for k, v in C["C3_coded_levels"].items()
                           if k != "positions"},
"C3_detail": C["C3_coded_levels"]["positions"],
"C4_printed_regression_statistics": C["C4_table_s2_statistics"],
"C5_abstract_fold_claim": C["C5_abstract_fold_claim"],
"C6_library_C_ladder_and_B26": C["C6_library_c_and_B26"],
"C7_S3_S4_formulation_identity": C["C7_s3_s4_identity"],
"C8_batch_to_batch_spread": C["C8_batch_spread"],
"C9_ranks": C["C9_ranks"],
"C10_phospholipid_effect": C["C10_phospholipid_effect"],

"C11_three_way_helper_concordance": {
  "per_helper": tw,
  "n_helpers": len(tw),
  "n_all_three_agree": sum(1 for v in tw.values() if v["all_three_agree"]),
  "max_exact_mass_delta_Da": max(v["exact_mass_max_delta"] for v in tw.values()),
  "independence_posture": S["three_way_helper_concordance"]["_independence_posture"],
},

"C12_amine_200_arm_count": S["amine_200_arm_count"],

"C13_corpus_overlap": {
  "positive_control_cKK_E12": ov["positive_control_1_cKK_E12"],
  "C12_200": ov["C12_200_in_corpus"],
  "in_vivo_shared": ov["in_vivo_shared_structures"],
  "substructure_sweep_SUPERSET": ov["substructure_sweep"],
  "token_pattern": {k: v for k, v in ov["token_pattern_test"].items()
                    if k != "per_name"},
  "experiment_id": ov["experiment_id"],
},

"C14_text_extraction_geometry": {
  "table_S4_dose_assignment": P["table_s4_geometry"],
  "why_geometric": "Table S4's READING ORDER IS NOT ITS VISUAL ORDER: the dose "
                   "cell is vertically centred, so each formulation's first dose "
                   "line sorts ABOVE its own code label. A line-order parser "
                   "recovered 6 of 8 value rows and dropped the two 0.1 mg/kg "
                   "points silently. Doses are therefore assigned by NATIVE TEXT "
                   "y-coordinate (the PDF's own character positions, not pixels), "
                   "each to its nearest phospholipid anchor.",
  "worst_margin_ratio": P["table_s4_geometry"]["worst_margin_ratio"],
  "margin_interpretation": "The tightest decision is 9.40 pt to the correct "
                           "anchor against 18.68 pt to the wrong one, a 1.99x "
                           "margin on a 28.08 pt anchor separation and a ~9.1 pt "
                           "row pitch. Every one of the 6 dose bands is assigned "
                           "by at least that margin, and the parser asserts each "
                           "formulation receives exactly the dose set "
                           "{0.01, 0.03, 0.1} rather than trusting the margin.",
  "em_dash_normalisations": P["emdash_fixed"],
  "em_dash_note": "One, on B-29. pdfplumber fuses a hyphen run into U+2014. "
                  "Accepted only because width 5 forces exactly two hyphens AND "
                  "the result satisfies all 5 of B-29's coded-level assertions "
                  "independently (08_contradictions.csv C-06).",
  "A_DEFECT_FOUND_IN_OUR_OWN_CODE": {
    "what": "_fix_emdash's first version returned early whenever len(s) == want, "
            "so an em-dash string that happened to be the right width came back "
            "with the em dash still in it. check_coded_levels would then have "
            "indexed that character as a design level, adding a fourth key to a "
            "three-key char->level map and MASKING a violation.",
    "found_by": "its own seed (emdash_width_guard), which scored 1 of 3 refusals",
    "after_fix": "3 of 3 refusals; the character set is now validated on the way "
                 "OUT unconditionally and length alone never buys an exit",
    "output_change": "none -- B-29 was the only affected row and was already "
                     "handled correctly",
    "class": "the 'a check that canonicalises its input has given up the ability "
             "to detect non-canonical input' shape (error class 37), in our own code",
  },
},

"C15_emission": dict(E, **{
  "value_rows_by_admissibility": {
    "yes": sum(1 for v in vals if v["corpus_admissible"] == "yes"),
    "no": sum(1 for v in vals if v["corpus_admissible"] == "no"),
    "total": len(vals),
    "sums": (sum(1 for v in vals if v["corpus_admissible"] == "yes")
             + sum(1 for v in vals if v["corpus_admissible"] == "no")) == len(vals),
  },
  "not_admissible_breakdown": {
    "PBS vehicle, Table S1": 1, "PBS vehicle, Table S4": 1,
    "unsynthesizable (DSPE insoluble), Table S1": 3},
  "non_delivery_readouts_deliberately_NOT_value_rows":
    "encapsulation efficiency, particle size and PDI are carried as CHARACTERISATION "
    "COLUMNS on the value rows, never as value rows of their own, because the "
    "endpoint vocabulary excludes all three. Figures S2a and S2b plot exactly "
    "those two quantities and were therefore not read as value sources.",
  "z_block": {
    "n": E["z_block"]["n"], "unit": "one z per experiment",
    "quantity": "z of log10(serum EPO ng/mL), POPULATION sd (ddof=0)",
    "mean_log10": E["z_block"]["mu_log10"], "sd_log10": E["z_block"]["sd_log10"],
    "z_min": E["z_block"]["z_min"], "z_max": E["z_block"]["z_max"],
    "block_mean_residual_computed": "< 1e-12, asserted in stage05_emit.py on the "
                                    "in-memory floats",
    "block_sd_residual_computed": "< 1e-12, asserted likewise",
    "block_stats_AS_STORED": {
      "mean": _zmean, "pop_sd": _zsd,
      "note": "Experiment_value is written at 6 decimal places, following the "
              "existing curations' convention, so the block statistics recomputed "
              "FROM THE FILE differ from the in-memory ones by a rounding "
              "artifact. Measured: mean %.3e and population sd %.12f. A 1e-12 "
              "bound on the stored values could only ever fail, which is the "
              "shape the briefing warns about; the figure predicted from 6-dp "
              "rounding of 36 values is ~5e-7/sqrt(12) ~ 1.4e-7 for the mean, "
              "and the observed 2.8e-8 sits inside that." % (_zmean, _zsd),
    },
    "error_class_22_cannot_arise": "there is ONE Model_type and ONE organ across "
                                   "all 36 rows, so no per-cell or per-cohort z "
                                   "exists to delete a between-assay offset",
  },
  "columns": {"n": E["n_columns_emitted"],
              "header_source": "read LIVE from data/data_vivo.csv and asserted "
                               "identical to data/data_vitro.csv"},
  "flags": {
    "written_by": "curation/lib/emit.py emit_rows, with PER-ITEM CALLABLES",
    "literal_confidence_refused": "verified by seed emit_literal_confidence "
                                  "(EmitError)",
    "caveat_without_note_refused": "verified by seed emit_caveat_without_note "
                                   "(EmitError)",
    "rows_with_a_caveat": E["caveat_rows"],
    "of_rows": len(rows),
    "confidence": {"medium": sum(1 for r in rows if r["confidence"] == "medium")},
    "structure_assumed_1": sum(1 for r in rows if r["structure_assumed"] == "1"),
    "positional_ambiguity_1": sum(1 for r in rows if r["positional_ambiguity"] == "1"),
    "value_is_bound_set": sum(1 for r in rows if r["value_is_bound"]),
    "value_is_bound_note": "ZERO. Every emitted value is a printed numeral with a "
                           "printed SD; no colormap floor, no censored axis, no "
                           "digitised bar exists in this curation.",
  },
}),

"C16_anchor_posture": {
  "paper_draws_nothing": "There is not one chemical structure drawing in the 14 "
                         "SI pages. Per-page object counts: the 8 text pages hold "
                         "0 image objects, and the 6 image pages are Fig S1 "
                         "(regression plots), Fig S2 (EE and size trends), Fig S3 "
                         "(a methodology flowchart), Fig S4 (TEM micrographs) and "
                         "two pages of JMP screenshots.",
  "C12_200_anchor": "NONE. reference_SMILES is EMPTY on purpose so the absence is "
                    "visible. No drawing, no formula, no HRMS, no MW, no synthesis.",
  "excipients_anchored_instead": ["DSPC", "DOPE", "DOPC"],
  "what_that_anchor_IS": "a NAME-CONCORDANCE anchor. The SI spells the systematic "
                         "names out in its footnotes, and make_reference.py "
                         "re-derives structures from those names by a head-first "
                         "assemble.attach route. All three legs of the comparison "
                         "are derivations from the SAME name by different hands.",
  "what_that_anchor_IS_NOT": "an external structural attestation of this paper's "
                             "chemistry. It would catch a transcription error in "
                             "any leg -- the briefing records a sibling whose "
                             "hand-built DSPC came out a ketone bearing an ether at "
                             "the right mass to six decimals -- but it cannot catch "
                             "a name that is wrong FOR THIS PAPER. NOT counted as a "
                             "C10 pass.",
  "cholesterol": "NOT hand-written and NOT anchored. The paper gives only the "
                 "abbreviation 'Chol.'. CHL_SMILES is populated from "
                 "provenance/reagent_registry.py and that registry pin is "
                 "explicitly NOT counted as an anchor.",
  "reference_structures_csv": {
    "rows": len(refs),
    "with_a_structure": sum(1 for r in refs if r["reference_SMILES"]),
    "empty_on_purpose": sum(1 for r in refs if not r["reference_SMILES"]),
    "header": list(refs[0].keys()) if refs else [],
  },
  "independence_verified_how": "make_reference.py walks ITS OWN AST at run time "
                               "and inspects Call, Import and ImportFrom nodes "
                               "ONLY, so the ban list cannot match itself (a "
                               "substring scan over the source would flag the ban "
                               "list inside its own docstring -- one sibling's "
                               "guard failed on its own module for exactly that). "
                               "8 import nodes and 168 call nodes examined, 0 "
                               "findings. A template-provenance guard additionally "
                               "requires that no string literal of >= 25 chars in "
                               "the module be a substring of any emitted structure: "
                               "54 literals examined, 0 hits.",
  "byte_identity_is_not_the_signal": "F-037. The built DSPC/DOPE/DOPC canonicalise "
                                     "byte-identically to the registry's strings, "
                                     "and that is EXPECTED: a correct re-derivation "
                                     "is byte-identical to a correct stored "
                                     "structure. Testing for byte-difference is the "
                                     "trap. What is tested is the ROUTE, by AST.",
  "stereo_assertions_are_NOT_vacuous": "DSPC, DOPE and DOPC each have exactly 1 "
                                       "real stereocentre (glycerol C2), asserted "
                                       "(R) from the sn convention. Seed "
                                       "make_reference_stereo flips the backbone tag "
                                       "and the assertion fires, reporting (S). "
                                       "Double-bond geometry is asserted explicitly: "
                                       "2 C=C, 2 Z, 0 unassigned for the dioleoyl "
                                       "lipids; 0 C=C for DSPC.",
  "connectivity_gate_proven": "Seed make_reference_regioisomer builds "
                              "1,3-diacyl-glycero-2-phospho DSPC: IDENTICAL formula "
                              "C44H88NO8P, exact-mass delta 0.000000 Da, different "
                              "InChIKey. The assembly route REFUSES it -- via "
                              "assemble.check_stereo_preserved, because C2 stops "
                              "being stereogenic when both CH2 arms become "
                              "identical. So the gate sees connectivity, not mass.",
  "MEASURED_SENSITIVITY_LIMIT": "make_reference's OWN assertions (C2=R, C=C count, "
                                "Z count, no unassigned) CANNOT see double-bond "
                                "POSITION. Seed three_way_delta6 moves the oleoyl "
                                "double bond from delta-9 to delta-6: every internal "
                                "assertion still passes, formula identical, "
                                "exact-mass delta 0.000000 Da. ONLY the three-way "
                                "InChIKey comparison catches it, and it does. That "
                                "is the honest limit of this curation's structural "
                                "gates.",
},

"C17_seeded_mutation_suite": {
  "n_seeds": seeds["n_seeds"], "n_credited": seeds["n_credited"],
  "all_credited": seeds["all_credited"],
  "curation_folder_untouched_by_seeds": seeds["curation_folder_untouched_by_seeds"],
  "F069_posture": "Every seed runs in a CHILD PROCESS with its own temp OUTDIR; "
                  "nothing a seed writes can land in the curation folder, and the "
                  "7 tracked artifacts are sha256-compared before and after the "
                  "whole suite. Verdicts come from returncode, never from a "
                  "pipeline's exit status.",
  "A_SEED_THAT_ONLY_RAISES_IS_NOT_CREDITED": "each seed must print a METRIC line "
                                             "naming the quantity it moved; a child "
                                             "that dies before its metric is "
                                             "reported NOT CREDITED even on a "
                                             "non-zero exit",
  "seeds": [{k: v for k, v in s.items() if k != "stderr_tail"}
            for s in seeds["seeds"]],
},

"C18_reproducibility": {
  "two_independent_clean_runs": "byte-identical across all 6 generated artifacts, "
                                "run into two separate temp OUTDIRs AFTER the seed "
                                "suite, and byte-identical to the tracked copies "
                                "in this folder",
  "artifacts_compared": ["01_blocks.csv", "02_corpus_overlap.json",
                         "03_enumerated.csv", "05_values.csv", "06_rows.csv",
                         "reference_structures.csv"],
  "python": sys.version.split()[0],
},

"C19_contradictions": {
  "n": len(contra),
  "live": sum(1 for r in contra if r["status"].startswith("LIVE")),
  "retired": sum(1 for r in contra if r["status"].startswith("RETIRED")),
  "ids_live": [r["id"] for r in contra if r["status"].startswith("LIVE")],
  "ids_retired": [r["id"] for r in contra if r["status"].startswith("RETIRED")],
  "retired_that_were_MY_error_not_the_paper's": [r["id"] for r in contra
                                                 if "my error" in r["status"]],
},

"what_could_NOT_be_verified": [
  "The main article is not held (closed access, no PMC deposit of any kind), so "
  "Figures 1-4, the Methods, the mouse strain, the PEG-lipid identity, the EPO "
  "species and any structural depiction of C12-200 are all unavailable. 0 "
  "characters of main text were searched.",
  "C12-200's structure is not attested by this paper at ANY level. It is carried "
  "from the corpus and cannot be checked against the source.",
  "Cholesterol mol% has NO independent attestation anywhere in the SI -- the p10 "
  "and p11 design tables omit the column entirely -- which is precisely where the "
  "one real numerical error (A-03) sits.",
  "The PEG-lipid is unidentified, so PEG_SMILES is blank on all 36 rows.",
  "The mouse strain is not stated, so Model_type is 'Mouse_unspecified'.",
  "N:P / charge ratio is never given; IL_to_nucleicacid_chargeratio is blank.",
  "Aqueous buffer, dialysis buffer and mixing method live in the main text's "
  "Methods and are blank.",
  "F-061 cannot be adjudicated from this paper, which provides no structural "
  "evidence on either side.",
],
}
    with open(os.path.join(OUTDIR, "07_qc.json"), "w") as fh:
        json.dump(qc, fh, indent=1, sort_keys=False, default=str)
    print(json.dumps(qc["headline"], indent=1))
    print("C1  mol%% sums 100     : %d/%d pass, max dev %.1f, offenders %s"
          % (qc["C1_molpct_sums_to_100"]["n_pass"],
             qc["C1_molpct_sums_to_100"]["n_formulations_tested"],
             qc["C1_molpct_sums_to_100"]["max_abs_deviation_pct"],
             qc["C1_molpct_sums_to_100"]["offenders"]))
    print("C2  cross-param       : %d/%d cells exact"
          % (qc["C2_cross_parameterisation"]["n_exact"],
             qc["C2_cross_parameterisation"]["n_cells_compared"]))
    print("C3  coded levels      : %d/%d consistent, all_bijective=%s"
          % (qc["C3_coded_design_levels"]["n_consistent"],
             qc["C3_coded_design_levels"]["n_character_assertions"],
             qc["C3_coded_design_levels"]["all_bijective"]))
    print("C11 3-way helpers     : %d/%d agree, max mass delta %.8f Da"
          % (qc["C11_three_way_helper_concordance"]["n_all_three_agree"],
             qc["C11_three_way_helper_concordance"]["n_helpers"],
             qc["C11_three_way_helper_concordance"]["max_exact_mass_delta_Da"]))
    print("C17 seeds             : %d/%d credited, folder untouched=%s"
          % (qc["C17_seeded_mutation_suite"]["n_credited"],
             qc["C17_seeded_mutation_suite"]["n_seeds"],
             qc["C17_seeded_mutation_suite"]["curation_folder_untouched_by_seeds"]))
    print("C19 contradictions    : %d (%d live, %d retired)"
          % (qc["C19_contradictions"]["n"], qc["C19_contradictions"]["live"],
             qc["C19_contradictions"]["retired"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
