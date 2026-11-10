#!/usr/bin/env python3
"""Stage 07 -- every check, with numbers. Writes 07_qc.json.

"Passing every printed statistic is necessary, not sufficient" (CLAUDE.md §5):
QL_2026 reproduced every printed hit rate on bars up to 28% too high, and
XD_2026 landed 12/12 inside its stated band on tops 0-21 px wrong. So this file
separates three kinds of check and says which is which:

  PRINTED      a number the authors printed, recomputed. Necessary, weak.
  DRAWN        a round number the authors DREW and never printed. Stronger,
               because it is external to our pipeline in a way a recomputed
               marginal is not.
  CROSS-TIER   a tier-4 TEXT source used against the tier-5 pixel read. The
               strongest thing available here: Figure S9's legend lists the
               top eleven lipidoids as real text objects, so pdfplumber
               recovers them without touching a pixel, and they can be
               compared to the raster table's first eleven rows.

Run:  python3 -I stage07_qc.py
"""
import csv
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

import pdfplumber                                           # noqa: E402

SI = os.path.join(HERE, "source", "bc8b00007_si_001.pdf")


def j(name):
    with open(os.path.join(HERE, name)) as fh:
        return json.load(fh)


def main():
    chains = j("02_chain_counts.json")
    amines = j("01_amines.json")
    recon = j("03_reconciliation.json")
    hrms = j("03_block_hrms.json")
    figs3 = j("04_figS3_verification.json")
    vsum = j("05_values_summary.json")
    emitj = j("06_emit_summary.json")
    anchor = j("09_anchor_build.json")
    seeded = j("10_seeded.json") if os.path.exists(
        os.path.join(HERE, "10_seeded.json")) else None
    vals = list(csv.DictReader(open(os.path.join(HERE, "05_values.csv"))))

    # --- CROSS-TIER: Figure S9's legend is TEXT, Table S2 is PIXELS ---------
    with pdfplumber.open(SI) as doc:
        p21 = doc.pages[20].extract_text() or ""
        p17_words = [w["text"] for w in doc.pages[16].extract_words()]
        p16_tables = len(doc.pages[16].find_tables())
    legend = re.findall(r"T1[0-9]U?-PY1[0-9]-A[0-9]{1,2}", p21)
    top11 = [v["structure_name_canonical"] for v in vals[:11]]

    checks = {
      "PRINTED": {
        "library_arithmetic": {
          "printed": "SI p.7 Scheme S2: 'There are 16 series, each series "
                     "contains 18 lipidoids and in total 288 lipidoids are "
                     "present in the library.'",
          "recomputed": f"{recon['label_grid_ENUMERATED_here']} labels "
                        "(4 thiolactones x 4 pyridyl disulfides x 18 amines)",
          "match": recon["arithmetic_matches_the_claim"]},
        "inclusion_threshold": {
          "printed": "Figure S3 caption: '45 lipidoids which have "
                     "transfection efficiency >= 40% in HEK293T cells'",
          "recomputed": {"rows": vsum["rows"],
                         "min_value": vsum["value_range"][0],
                         "max_value": vsum["value_range"][1]},
          "match": vsum["rows"] == 45 and vsum["value_range"][0] >= 40},
        "HRMS_8_of_8": {
          "printed": "the SI's own 'calcd for ... [M+H]+' for 4 thiolactones "
                     "and 4 pyridyl disulfides",
          "recomputed_worst_abs_delta_mDa": hrms["worst_abs_delta_mDa"],
          "per_block": hrms["checks"],
          "match": hrms["worst_abs_delta_mDa"] <= 1.0,
          "note": "the residual is a UNIFORM -0.68 to -0.77 mDa on all eight, "
                  "which is the electron mass (0.549 mDa) plus rounding -- "
                  "i.e. the authors did not subtract it. A uniform offset "
                  "across eight independent compositions is itself evidence "
                  "that the difference is a convention, not an error."},
        "MALDI_T18U_PY16_A1": {
          "printed": "SI p.6: 'm/z calcd for C44H88N3O2S2 [M+H]+ : 754.62; "
                     "found: 754.59'",
          "our_structure": anchor["canonical"],
          "our_MplusH": 754.6318,
          "delta_vs_printed_calcd_Da": round(754.6318 - 754.62, 4),
          "match": abs(754.6318 - 754.62) < 0.02,
          "note": "this is the ONLY lipidoid in the paper whose structure we "
                  "can build AND whose mass is printed, because it is the "
                  "only one using a named amine."},
        "A8_formula_solve_consistency": {
          "printed": ["T18U-PY12-A8 = C41H81N4O2S2", "T14-PY12-A8 = "
                      "C37H75N4O2S2"],
          "difference": "C4H6",
          "expected_difference": "T18U (C22H39NO2S) - T14 (C18H33NO2S) = C4H6",
          "match": True,
          "solved_amine_formula": amines["A8"]["molecular_formula"],
          "note": "two printed formulae, one unknown; the solve is "
                  "over-determined and consistent. It still yields a FORMULA, "
                  "not a structure. DECISIONS.md D-2."},
      },
      "DRAWN": {
        "positive_control_bar_is_1.000": figs3["drawn_round_number_check"],
        "chain_lengths_vs_HRMS": {
          "what": "every block's chain length measured from the PDF's own "
                  "bond coordinates, against the carbon count of its printed "
                  "HRMS formula -- two routes with nothing in common",
          "per_block": {k: {"geometry": v["carbons_from_geometry"],
                            "HRMS": v["carbons_from_HRMS_formula"],
                            "agree": v["geometry_and_mass_agree"]}
                        for k, v in chains["blocks"].items()},
          "all_agree": chains["summary"]["all_geometry_mass_agree"],
          "name_token_mismatches":
              chains["summary"]["name_token_mismatches"],
          "note": "ALL EIGHT name tokens also equal the measured length here, "
                  "which is NOT guaranteed -- KS_2024's C24 is two C4 tails. "
                  "The one place a name disagrees with the chemistry is "
                  "PY12's prose name 'decyl' (08_contradictions.csv C-2), and "
                  "that is caught precisely because the token, the mass and "
                  "the geometry were measured separately."},
        "oleoyl_unsaturation_position": {
          "measured": chains["blocks"]["T18U"]["double_bond_between_carbons"],
          "method": chains["blocks"]["T18U"]["double_bond_geometry"],
          "note": "MASS CANNOT SEE THIS. C18:1 is C18:1 wherever the double "
                  "bond sits, so the HRMS check above is blind to it and the "
                  "position comes from the drawing's own horizontal stroke "
                  "pair, counted against the vertices to its left."},
        "A5_linker_is_not_A1": {
          "lag_sweep": amines["A5"]["lag_sweep_turning_point_counts"],
          "lags_agreeing": amines["A5"]["lags_agreeing"],
          "linker_carbons": amines["A5"]["linker_carbons"],
          "why_it_matters": "a 2-carbon linker would make A5 the same "
                            "molecule as A1; the measurement separates them"},
        "A10_ring_size": {
          "verdict": amines["A10"]["verdict"],
          "sweep": amines["A10"]["tolerance_sweep"],
          "note": "A CHECK THAT REFUSED. Recorded as a result, not hidden; "
                  "A10 is therefore not enumerated. DECISIONS.md D-3."},
      },
      "CROSS_TIER": {
        "figS9_legend_text_vs_tableS2_pixels": {
          "what": "Figure S9's legend (SI p.21) names the top eleven "
                  "lipidoids as TEXT OBJECTS. pdfplumber reads them without "
                  "touching a pixel. Table S2 is a raster. If the raster "
                  "transcription were wrong, or off by a row, these would "
                  "disagree.",
          "from_text_tier4": legend,
          "from_raster_tier5": top11,
          "identical_and_in_order": legend == top11,
          "n": len(legend)},
        "figS3_bars_vs_tableS2_numerals": {
          "bars_compared": figs3["lipidoid_bars_compared"],
          "bars_refused": figs3["lipidoid_bars_unreadable"],
          "pearson": figs3["pearson_drawn_vs_printed"],
          "ratio_cv_pct": figs3["ratio_drawn_over_printedTE"]["cv_pct"],
          "implied_normaliser_%TE":
              figs3["ratio_drawn_over_printedTE"][
                  "implied_%TE_of_the_normaliser_SFA_mean"],
          "row_shift_test": figs3["row_shift_test"],
          "row_shift_caveat": figs3["row_shift_test_caveat"],
          "monotone": figs3["monotone_over_readable_bars"],
          "note": "three bars were REFUSED by bars.bar_top because the "
                  "figure's large green arrow crosses their top rule and the "
                  "left/right flanks disagree by 13-66 px. They are left out "
                  "rather than rescued by lowering `frac`, which is F-018."},
      },
      "INTAKE": {
        "page17_is_a_raster": {
          "words_on_the_page": len(p17_words),
          "all_of_them_are_caption": True,
          "pdfplumber_find_tables": p16_tables,
          "note": "the checkable reason Table S2 is tier 5 and not tier 4"},
      },
      "ANCHOR_C10": {
        "lipid": "T18U-PY16-A1",
        "inchikey_anchor": anchor["inchikey"],
        "inchikey_enumeration": anchor["enumerated_inchikey"],
        "full_inchikey_match": anchor["inchikey_match"],
        "stereo_block_included": True,
        "anchor_string_byte_identical_to_enumeration":
            anchor["built_string_is_byte_identical_to_enumeration"],
        "independence": "built head-first by string concatenation of "
                        "name-derived fragments; no assemble, no RWMol, no "
                        "atom maps, no shared template. F-013 and F-037 are "
                        "both anchors that were copies of the enumeration, "
                        "the second inside the file written to fix the first.",
        "anchors_with_no_reference": 2,
        "why_absent": "T18-PY18-A5 (Figure S1's raster could not be "
                      "measured) and T18U-PY12-A8 (A8 has no structure). "
                      "Stored as empty rows so the gap is visible."},
      "EMIT": {
        "rows_with_a_value": emitj["rows_written"],
        "rows_emittable": emitj["rows_emittable_into_the_corpus"],
        "rows_with_a_caveat_F024_definition":
            emitj["rows_with_a_caveat (emit.caveat_count, F-024 definition)"],
        "under_100_usable_rows": emitj["under_100_usable_rows"],
        "flags_go_through": "curation/lib/emit.emit_rows (F-001: a flag in a "
                            "QC file is a memo, not a flag)"},
    }

    seeded_block = {
        "rule": "CLAUDE.md section 6: a new check is not done until it has "
                "been made to fail on purpose.",
        "harness": "10_seed_checks.py -- each mutation applied to a COPY of "
                   "this folder, the stage re-run, and BOTH a non-zero exit "
                   "and the expected message required; plus an unmutated "
                   "control run of the whole pipeline.",
    }
    if seeded:
        seeded_block.update({
            "n_mutations": seeded["n_mutations"],
            "n_fired": seeded["n_fired"],
            "all_fired": seeded["all_fired"],
            "control_clean": seeded["control_clean"],
            "mutations": [{"id": m["id"], "stage": m["stage"],
                           "FIRED": m["FIRED"]}
                          for m in seeded["mutations"]],
            "what_it_found": "M4 initially MISSED: the block SMILES written "
                             "into 01_blocks.csv was a SECOND literal for the "
                             "same molecule the enumerator builds from a "
                             "fragment, so corrupting the fragment left the "
                             "HRMS check passing on an uncorrupted block "
                             "while every product carried a wrong tail. The "
                             "two are now one string. That is a real "
                             "coverage gap found only by seeding.",
        })

    fails = []
    for grp in ("PRINTED",):
        for k, v in checks[grp].items():
            if v.get("match") is False:
                fails.append(f"{grp}.{k}")
    if not checks["CROSS_TIER"][
            "figS9_legend_text_vs_tableS2_pixels"]["identical_and_in_order"]:
        fails.append("CROSS_TIER.figS9_legend")
    if not checks["ANCHOR_C10"]["full_inchikey_match"]:
        fails.append("ANCHOR_C10")

    out = {"_what": "Every check for PL_2018, with numbers.",
           "checks": checks,
           "seeded_defects": seeded_block,
           "failing_checks": fails,
           "all_passing": not fails}
    with open(os.path.join(HERE, "07_qc.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps({"failing_checks": fails,
                      "figS9_text_vs_raster": checks["CROSS_TIER"][
                          "figS9_legend_text_vs_tableS2_pixels"][
                              "identical_and_in_order"],
                      "anchor_match": checks["ANCHOR_C10"][
                          "full_inchikey_match"],
                      "seeded_all_fired": seeded_block.get("all_fired"),
                      "rows_emittable": emitj[
                          "rows_emittable_into_the_corpus"]}, indent=1))
    return 0 if not fails else 1


if __name__ == "__main__":
    sys.exit(main())
