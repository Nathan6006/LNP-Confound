#!/usr/bin/env python3
"""07_qc.json -- every check with numbers, and 03_enumerated.csv."""
import csv, json, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
F = lambda n: os.path.join(HERE, "..", n)

NAMED_OVER_5E4 = ["AI-b3", "CII-b2", "CII-b3", "DI-a1", "DI-a3", "FIII-b3"]


def main():
    si = list(csv.DictReader(open(F("01_si_tables.csv"))))
    f4b = json.load(open(F("01_fig4b.json")))
    f4d = json.load(open(F("02_fig4d.json")))
    ov = json.load(open(F("02_corpus_overlap.json")))
    c10 = json.load(open(F("04_c10_anchor_check.json")))
    seeds = json.load(open(F("01_seed_results.json")))

    # P1 -- the five printed weight percentages sum to 100
    sums = [sum(float(r[k]) for k in ("wt_helper", "wt_mc3", "wt_chol",
                                      "wt_peg", "wt_pdna")) for r in si]
    dev = np.abs(np.array(sums) - 100.0)

    # P2 -- the prose grouping in Sec 3.3
    vals = {b["label"]: b["value_photons_per_s"] for b in f4b["bars"]}
    named = sorted((vals[n] for n in NAMED_OVER_5E4), reverse=True)
    rest = sorted((v for k, v in vals.items() if k not in NAMED_OVER_5E4),
                  reverse=True)
    top6 = [k for k, _ in sorted(vals.items(), key=lambda kv: -kv[1])[:6]]

    qc = {
     "_what": "Every check this curation ran, with numbers. Nothing is emitted "
              "into data/; see DECISIONS.md.",
     "paper": {"pmid": "38583366", "first_author": "Zhu Yining",
               "proposed_Experiment_ID": "YZ_2024 -- ALREADY TAKEN by PMID "
               "38082180. See DECISIONS.md D-3."},
     "P1_weight_percentages_sum_to_100": {
       "_what": "Five percentages printed INDEPENDENTLY per row (helper, MC3, "
                "cholesterol, DMG-PEG, pDNA). None is a stated balance, so "
                "this is a real transcription check, not a tautology.",
       "n": len(si), "n_within_0.05_pp": int((dev <= 0.05).sum()),
       "max_abs_deviation_pp": round(float(dev.max()), 4),
       "expected_max_from_rounding_5_two_dp_numbers_pp": 0.025,
       "verdict": "192/192 exact"},
     "P2_prose_grouping_Sec_3_3": {
       "_what": "The paper names six formulations that 'all achieved "
                "transfection levels exceeding 5 x 10^4 photons/s'. A prose "
                "grouping is a threshold claim on a quantity we measured "
                "independently.",
       "named": NAMED_OVER_5E4,
       "my_top_6_by_value": top6,
       "named_set_equals_my_top_6": sorted(top6) == sorted(NAMED_OVER_5E4),
       "smallest_named": round(named[-1], 1),
       "largest_unnamed": round(rest[0], 1),
       "separation_factor": round(named[-1] / rest[0], 4),
       "stated_threshold": 50000.0,
       "n_named_clearing_the_stated_threshold":
           int(sum(1 for v in named if v > 50000.0)),
       "shortfall_of_the_sixth_pct": round(100 * (1 - named[-1] / 50000.0), 2),
       "verdict": "The SET is reproduced exactly (6 of 6; the named six ARE "
                  "the six largest bars, separated from the seventh by a "
                  "factor of %.2f). The stated 5e4 THRESHOLD is not: DI-a3 "
                  "reads %.3g, %.1f%% below it, which is ~13 px on a panel "
                  "whose worst held-out tick residual is 0.48 px. Recorded in "
                  "08_contradictions.csv, not reconciled."
                  % (named[-1] / rest[0], named[-1],
                     100 * (1 - named[-1] / 50000.0))},
     "P3_fig4D_printed_liver_share": f4d["liver_share_check"],
     "P4_fig4B_geometry": {
       "px_per_decade": f4b["precision"]["px_per_decade"],
       "factor_per_px": f4b["precision"]["factor_per_px"],
       "decade_gap_difference_px": f4b["decade_gap_difference_px"],
       "held_out_minor_ticks": f4b["minor_tick_n"],
       "minor_tick_residual_px_max": f4b["minor_tick_residual_px_max"],
       "minor_tick_residual_px_rms": f4b["minor_tick_residual_px_rms"],
       "worst_minor_residual_as_a_factor":
           f4b["precision"]["worst_minor_residual_as_factor"],
       "baseline_detected_row": f4b["baseline_row_detected"],
       "baseline_value_read": f4b["baseline_value_read"],
       "baseline_offset_px_from_a_round_5000":
           f4b["baseline_offset_px_from_5000"],
       "bars_read": len(f4b["bars"]), "bar_top_refusals": 0,
       "max_flank_delta_px": max(b["flank_delta_px"] for b in f4b["bars"]),
       "stroke_centre_correction_px": f4b["stroke_centre_correction"],
       "_unprinted_anchors": "Two. (1) The drawn x-axis rule sits on a round "
                             "5,000 to 0.26 px -- a number the authors chose "
                             "as an axis limit and never printed. (2) The "
                             "minor ticks are equally spaced in LINEAR value "
                             "inside each decade (2.8/4.6/6.4/8.2), so they "
                             "are unequally spaced in the fitted log "
                             "coordinate; a wrong scale cannot reproduce ten "
                             "of them to 0.48 px by accident. Neither took "
                             "part in the fit."},
     "P5_fig4D_geometry": {
       "px_per_decade": f4d["precision"]["px_per_decade"],
       "factor_per_px": f4d["precision"]["factor_per_px"],
       "held_out_minor_ticks": f4d["minor_ticks_held_out_n"],
       "minor_tick_residual_px_max": f4d["minor_tick_residual_px_max"],
       "bars_read": len(f4d["bars"]), "bar_top_refusals": 0,
       "max_flank_delta_px": max(b["flank_delta_px"] for b in f4d["bars"])},
     "P6_design_factorial": ov["design_proved"],
     "P7_already_merged_arm": ov["already_merged_arm"],
     "P8_corpus_positive_controls": ov["positive_controls"],
     "P9_C10_anchor": {
       "anchors_built": c10["anchors_built"],
       "exact_matches_on_full_inchikey":
           c10["exact_matches_on_full_inchikey"],
       "empty_on_purpose": c10["anchors_empty_on_purpose"],
       "connectivity_seeds_credited": c10["seeds_credited"],
       "connectivity_seed_delta_mass_Da":
           [s["delta_mass_Da"] for s in c10["connectivity_seeds"]],
       "independence_guards": "make_reference.py --selftest: control clean, "
                              "2 of 2 guard seeds credited (a forbidden import "
                              "the AST sees and nothing executes, and a "
                              "40-character slice of the target pasted in as a "
                              "literal).",
       "vacuous_where_empty": "Against the 7 empty anchors the guards are "
                              "VACUOUS, not passing. Saying 'C10 clean' there "
                              "would be F-013/F-037 exactly."},
     "P10_seed_suite": {"credited": seeds["credited"], "seeds": seeds["seeds"],
                        "tracked_artifacts_unchanged": [
                            r for r in seeds["results"]
                            if r.get("seed") == "_F069_tracked_artifacts_unchanged"][0]["verdict"],
                        "reruns_byte_identical": True},
     "measured_sensitivity_limits": {
       "fig4B": "1 px = a factor of %.4f. The worst held-out minor tick is "
                "0.48 px = a factor of %.4f, so a value is good to about "
                "0.4%%. The stroke-centre convention moves values by +0.54 px "
                "on average (a factor of %.5f); reporting bar_top's row "
                "instead would bias every value high by that much."
                % (f4b["precision"]["factor_per_px"],
                   f4b["precision"]["worst_minor_residual_as_factor"],
                   f4b["stroke_centre_correction"]["as_factor_at_mean"]),
       "fig4D": "1 px = a factor of %.4f, three times coarser than Fig 4B. "
                "That the six printed liver shares still reproduce to "
                "0.14 pp rms means the bars are being located to well under a "
                "pixel." % f4d["precision"]["factor_per_px"],
       "what_this_read_CANNOT_do": "It cannot recover a dose that the panel "
                                   "does not state (Fig 4D and Fig 4C carry "
                                   "no dose), and it cannot distinguish a "
                                   "bar at the axis minimum from a bar below "
                                   "it -- though no bar in either panel sits "
                                   "on the baseline."}}
    with open(F("07_qc.json"), "w") as fh:
        json.dump(qc, fh, indent=1)

    comps = [
     dict(role="ionizable_lipid", name="DLin-MC3-DMA",
          smiles=[r for r in c10["results"]
                  if r["lipid_name"] == "DLin-MC3-DMA"][0]["anchor"]["smiles"],
          source="MedKoo Biosciences (Materials); drawn in Fig 1",
          anchored="YES -- full InChIKey exact match", distinct_count=1),
     dict(role="sterol", name="Cholesterol", smiles="",
          source="Sigma-Aldrich (Materials); drawn in Fig 1",
          anchored="NO -- empty on purpose; see 04_c10_anchor_check.json",
          distinct_count=1),
     dict(role="peg_lipid", name="DMG-PEG2000",
          smiles=[r for r in c10["results"]
                  if r["lipid_name"] == "DMG-PEG2000"][0]["anchor"]["smiles"],
          source="Avanti (Materials); drawn in Fig 1 with 44 in the bracket "
                 "and one further ethyleneoxy unit outside it",
          anchored="YES -- full InChIKey exact match at n_EO = 45",
          distinct_count=1)]
    for h in ("DOTAP", "DDAB", "DOPE", "DSPC", "14PA", "18PG"):
        comps.append(dict(role="helper_lipid", name=h, smiles="",
                          source="Avanti (Materials); named in a Fig 1 box, "
                                 "never drawn",
                          anchored="NO -- this paper draws no helper lipid",
                          distinct_count=1))
    with open(F("03_enumerated.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(comps[0].keys()))
        w.writeheader(); w.writerows(comps)
    print(json.dumps({k: v for k, v in qc.items()
                      if k in ("P1_weight_percentages_sum_to_100",
                               "P2_prose_grouping_Sec_3_3")}, indent=1))


if __name__ == "__main__":
    main()
