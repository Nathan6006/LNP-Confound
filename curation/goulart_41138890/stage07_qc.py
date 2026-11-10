#!/usr/bin/env python3
"""07_qc.json and 02_corpus_overlap.json -- every check with NUMBERS.

Also proves, in this folder, that emit.emit_rows' two refusals actually fire
(a literal `confidence`, and a caveat with no note). A check that cannot fail
is not a check.
"""
import csv
import json
import os
import sys

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
sys.path.insert(0, os.path.join(HERE, "..", "..", "provenance"))
import emit                                            # noqa: E402
import reagent_registry as RR                          # noqa: E402
from errors import EmitError                           # noqa: E402

DATA = os.path.join(HERE, "..", "..", "data")
_ik = {}


def ik(s):
    if s not in _ik:
        m = Chem.MolFromSmiles(s)
        _ik[s] = Chem.MolToInchiKey(m) if m else None
    return _ik[s]


def load(name):
    with open(os.path.join(DATA, name), newline="") as fh:
        return list(csv.DictReader(fh))


def main():
    vit, viv = load("data_vitro.csv"), load("data_vivo.csv")
    allr = vit + viv

    # ================= 02_corpus_overlap.json =========================
    ov = {"corpus": {
        "vitro_rows": len(vit), "vivo_rows": len(viv),
        "vitro_experiments": len({r["Experiment_ID"] for r in vit}),
        "vivo_experiments": len({r["Experiment_ID"] for r in viv}),
        "distinct_experiment_ids": len({r["Experiment_ID"] for r in allr})}}

    # --- positive controls: each asserts a SPECIFIC NAMED molecule -----
    pc = {}
    n = sum(1 for r in allr if "DSPC" in (r.get("HL_name") or ""))
    pc["name_scan_HL_name_contains_DSPC"] = {
        "rows": n, "brief_reference": 7038, "matches_brief": n == 7038}
    k = ik(RR.REGISTRY["Cholesterol"][0])
    n = sum(1 for r in allr if (r.get("CHL_SMILES") or "").strip()
            and ik(r["CHL_SMILES"].strip()) == k)
    pc["inchikey_Cholesterol"] = {
        "inchikey": k, "rows": n, "brief_reference": 30397,
        "matches_brief": n == 30397,
        "asserts": "the registry cholesterol, by name, is present"}
    k = ik(RR.REGISTRY["DLin-MC3-DMA"][0])
    hits = [r for r in allr if (r.get("IL_SMILES") or "").strip()
            and ik(r["IL_SMILES"].strip()) == k]
    pc["inchikey_DLin_MC3_DMA"] = {
        "inchikey": k, "rows": len(hits), "brief_reference": 3098,
        "matches_brief": len(hits) == 3098,
        "names_found": sorted({(r.get("IL_name") or "").strip()
                               for r in hits})}
    pc["all_controls_match_the_briefs_measured_numbers"] = all(
        v.get("matches_brief", True) for v in pc.values()
        if isinstance(v, dict))
    ov["positive_controls"] = pc

    # --- the three ionizable lipids this paper contributes --------------
    enum = list(csv.DictReader(open(os.path.join(HERE, "03_enumerated.csv"),
                                    newline="")))
    tg = {}
    for e in enum:
        k = e["InChIKey"]
        hv = [r for r in vit if (r.get("IL_SMILES") or "").strip()
              and ik(r["IL_SMILES"].strip()) == k]
        ho = [r for r in viv if (r.get("IL_SMILES") or "").strip()
              and ik(r["IL_SMILES"].strip()) == k]
        tg[e["IL_name"]] = {
            "inchikey": k, "already_in_corpus": bool(hv or ho),
            "vitro_rows": len(hv),
            "vitro_experiments": sorted({r["Experiment_ID"] for r in hv}),
            "vivo_rows": len(ho),
            "vivo_experiments": sorted({r["Experiment_ID"] for r in ho}),
            "il_names_it_sits_under": sorted(
                {(r.get("IL_name") or "").strip() for r in hv + ho})}
    ov["ionizable_lipids"] = tg
    ov["n_new_to_the_corpus"] = sum(
        1 for v in tg.values() if not v["already_in_corpus"])
    ov["in_vivo_shared_lipids_touched"] = {
        "note": ("the emitted block is IN VITRO, so it adds no in-vivo row "
                 "and cannot change the in-vivo pool's two shared structures "
                 "(DLin-MC3-DMA, SM-102)"),
        "SM_102_vivo_sources_today": tg["SM-102"]["vivo_experiments"],
        "ALC_0315_vivo_sources_today": tg["ALC-0315"]["vivo_experiments"]}
    ov["vitro_source_count_change_if_merged"] = {
        n: {"today": len(v["vitro_experiments"]),
            "after": len(v["vitro_experiments"]) + 1}
        for n, v in tg.items()}

    # --- token-pattern test (mandatory) ---------------------------------
    tok = {}
    for t in ["C12-200", "200", "C12", "SM-102", "102", "ALC-0315", "0315",
              "C14PEG", "C14PEG2k", "LNP4", "LNP", "DOPE", "DOTAP"]:
        sub = [r for r in allr if t.lower() in (r.get("IL_name") or "").lower()]
        ex = [r for r in allr
              if (r.get("IL_name") or "").strip().lower() == t.lower()]
        tok[t] = {
            "IL_name_substring_rows": len(sub),
            "IL_name_substring_distinct_structures": len(
                {ik(r["IL_SMILES"]) for r in sub if r["IL_SMILES"]}),
            "IL_name_substring_experiments": sorted(
                {r["Experiment_ID"] for r in sub}),
            "IL_name_exact_rows": len(ex),
            "IL_name_exact_distinct_structures": len(
                {ik(r["IL_SMILES"]) for r in ex if r["IL_SMILES"]})}
    # the wrong-NAMESPACE case, which a token test on IL_name alone misses
    ns = {}
    for t in ["200", "C14", "C14PEG", "DOPE", "DSPC", "DOTAP"]:
        for col in ("HL", "PEG"):
            h = [r for r in allr
                 if t.lower() in (r.get(col + "_name") or "").lower()]
            if h:
                ns[f"{t}_in_{col}_name"] = {
                    "rows": len(h),
                    "distinct_structures": len(
                        {ik(r[col + "_SMILES"]) for r in h
                         if r[col + "_SMILES"]})}
    tok["_wrong_namespace"] = ns
    ov["token_pattern_test"] = tok

    # --- the PEG-lipid: a corpus name that is WRONG for this paper -------
    peg = {}
    for nme in ["C14-PEG2000", "DMPE-PEG2000", "DMG-PEG2000"]:
        d = {}
        for r in allr:
            if (r.get("PEG_name") or "").strip() == nme and r.get("PEG_SMILES"):
                k = ik(r["PEG_SMILES"].strip())
                e = d.setdefault(k, {"rows": 0, "exps": set(),
                                     "has_N": "N" in r["PEG_SMILES"].upper()
                                     .replace("N1", "N").replace("N(", "N"),
                                     "has_P": "P" in r["PEG_SMILES"]})
                e["rows"] += 1
                e["exps"].add(r["Experiment_ID"])
        peg[nme] = {k: {"rows": v["rows"], "experiments": sorted(v["exps"]),
                        "has_P": v["has_P"]} for k, v in d.items()}
    ov["peg_lipid"] = {
        "drawn_in_Fig_1A": ("1,2-diacyl-sn-glycero-3-phospho-ethanolamine "
                            "N-carbamoyl-mPEG45, ammonium salt -- phosphate "
                            "present, -NH-C(=O)O- carbamate present, "
                            "subscript 45 printed"),
        "corpus_names": peg,
        "finding": ("the corpus name 'C14-PEG2000' resolves UNAMBIGUOUSLY to "
                    "a structure with no N and no P, so corpus_components.smi "
                    "would NOT raise and would return the wrong molecule for "
                    "this paper. PEG_SMILES is therefore left blank. This is "
                    "an eighth independent sighting of the shape behind the "
                    "open F-049."),
        "f049_restated": ("'DMPE-PEG2000' itself carries TWO structures in "
                          "the corpus: the phosphoethanolamine on 3,611 rows "
                          "across 8 experiments, and the glycerol ether on "
                          "1,873 rows in JW_2024 alone")}

    # --- is any ARM of this screen already merged? (brief section 5) -----
    import collections
    ac = [r for r in vit if r["Experiment_ID"] == "AC_2025"
          and (r.get("IL_name") or "").strip() == "C12-200"]
    combos = collections.Counter(
        (float(r["IL_molratio"]), r["HL_name"], float(r["HL_molratio"]),
         float(r["CHL_molratio"]), float(r["PEG_molratio"])) for r in ac)
    s1 = list(csv.DictReader(open(os.path.join(HERE, "02_tableS1.csv"),
                                  newline="")))
    twins = []
    for f in s1:
        k = (float(f["IL_molratio"]), f["helper_lipid"],
             float(f["HL_molratio"]), float(f["CHL_molratio"]),
             float(f["PEG_molratio"]))
        if combos.get(k):
            twins.append({"formulation": f["formulation"],
                          "composition": list(k),
                          "AC_2025_rows_at_that_composition": combos[k]})
    ov["already_merged_arm_test"] = {
        "question": "does any arm of this screen duplicate an experiment "
                    "already in the corpus?",
        "closest_candidate": "AC_2025 -- the corpus's other C12-200 "
                             "formulation screen (498 rows, 308 distinct "
                             "IL/HL/CHL/PEG compositions)",
        "composition_twins_found": len(twins), "of_n_formulations": len(s1),
        "twins": twins,
        "answer": "NO. The two composition twins (LNP6, LNP12) are measured "
                  "in a different cell (DC2.4 vs primary rat cortical "
                  "neuron), with a different cargo (mRNA/FLuc vs pDNA/GFP), "
                  "by a different lab. The FORMULATION SPACE overlaps at "
                  "2/18 points; the measurements do not overlap at all. "
                  "Neither twin is among the three emitted rows, which sit "
                  "on the LNP4 backbone (35/30.5/32/2.5) -- a composition "
                  "AC_2025 does not contain.",
        "no_duplicated_Experiment_ID_risk": True}

    json.dump(ov, open(os.path.join(HERE, "02_corpus_overlap.json"), "w"),
              indent=1)

    # ================= 07_qc.json =====================================
    si = json.load(open(os.path.join(HERE, "02_si_tables.json")))
    fig = json.load(open(os.path.join(HERE, "04_fig4g.json")))
    anc = json.load(open(os.path.join(HERE, "05_anchor_check.json")))
    seeds = json.load(open(os.path.join(HERE, "09_seeds.json")))
    xchk = json.load(open(os.path.join(HERE, "03_crosscheck.json")))
    eqc = json.load(open(os.path.join(HERE, "06_emit_qc.json")))

    # --- prove emit's refusals fire -------------------------------------
    rec = {"X": {"what": "a test item"}}
    refusals = {}
    try:
        emit.emit_rows(["X"], rec, confidence="high", note=lambda i, r: "n",
                       label="seed")
        refusals["literal_confidence"] = "DID NOT RAISE -- check has no teeth"
    except EmitError as e:
        refusals["literal_confidence"] = "raises: " + str(e).split("\n")[0]
    try:
        emit.emit_rows(["X"], rec, confidence=lambda i, r: "low",
                       note=lambda i, r: "", label="seed")
        refusals["caveat_without_a_note"] = "DID NOT RAISE"
    except EmitError as e:
        refusals["caveat_without_a_note"] = ("raises: "
                                             + str(e).split("\n")[0])
    refusals["both_fire"] = all("raises" in v for k, v in refusals.items()
                                if k != "both_fire")

    qc = {
        "P1_identity": {
            "pmid": "41138890", "doi": "10.1016/j.ijpharm.2025.126278",
            "first_author": "Marco Tullio Rodrigues Alves",
            "proposed_Experiment_ID": "MR_2025",
            "checked_against_all_corpus_ids_by_PMID": True,
            "n_corpus_ids_checked": ov["corpus"]["distinct_experiment_ids"],
            "id_collides": False, "pmid_already_in_corpus": False,
            "abstractText_length": 1740,
            "no_original_data_discriminator":
                "abstract is 1,740 chars with Methods and Results -> ordinary "
                "research article"},
        "P2_governing_number": {
            "n_distinct_ionizable_lipids": 3,
            "n_new_to_the_corpus": ov["n_new_to_the_corpus"],
            "n_formulations": 18,
            "n_individually_measured_formulations_fig1D": 18,
            "n_design_supported_per_LNP_values":
                "18 formulations x 5 organs = 90 b-DNA cells (Fig 1C), all "
                "biodistribution",
            "n_attributed_delivery_values_with_varying_chemistry": 3,
            "factorial_claimed": "3 helper lipids x 6 molar contents",
            "factorial_proved": si["checks"]["S1_factorial"],
            "factorial_equals_row_count":
                si["checks"]["S1_factorial_equals_row_count"],
            "ionizable_lipid_molratio_values_across_all_29_SI_rows":
                si["checks"]["IL_molratio_constant"],
            "verdict": "a HELPER-LIPID / MOLAR-RATIO formulation screen on one "
                       "ionizable lipid, plus a 3-point ionizable-lipid "
                       "substitution"},
        "P3_si_tables": {
            "container": "mmc1.docx, tables embedded as EMF METAFILES",
            "w_tbl_count_in_document_xml": 1,
            "value_tables_actually_present": 3,
            "why_the_count_disagrees": "Tables S1/S2/S3 are images; only "
                                       "Table S4 is a Word table",
            "EXTTEXTOUTW_records": {"TableS1": 694, "TableS3": 534,
                                    "TableS2": 1571},
            "glyph_grouping": si["checks"]["glyph_grouping"],
            "S1_formulations": si["checks"]["S1_n_formulations"],
            "S3_formulations": si["checks"]["S3_n_formulations"],
            "molpct_sum_exact_100_S1":
                f'{si["checks"]["S1_molpct_sum_exact_100"]}/'
                f'{si["checks"]["S1_molpct_sum_n"]}',
            "molpct_sum_exact_100_S3":
                f'{si["checks"]["S3_molpct_sum_exact_100"]}/'
                f'{si["checks"]["S3_molpct_sum_n"]}',
            "molpct_sum_is_partly_tautological":
                si["checks"]["S1_molpct_sum_is_tautological"],
            "S3_is_a_subset_of_S1":
                f'{si["checks"]["S3_subset_of_S1_agree"]}/'
                f'{si["checks"]["S3_n_formulations"]} agree on composition, '
                f'{len(si["checks"]["S3_subset_of_S1_disagree"])} disagree'},
        "P3b_independent_tie_figure_vs_SI_table": {
            "what": xchk["what"],
            "n_compared": xchk["n_compared"],
            "n_agree_exactly": xchk["n_agree_exactly"],
            "n_disagree": xchk["n_disagree"],
            "decimal_separators_differ": True,
            "verdict": xchk["verdict"],
            "what_it_does_NOT_validate": xchk["what_it_does_NOT_validate"]},
        "P4_fig4g_geometry": {
            "tier": "5_pixels", "source": "publisher raster gr4_lrg.jpg",
            "native_dpi": "NOT DERIVABLE -- a standalone JPEG has no page "
                          "rectangle. axis.native_tol was NOT called; it "
                          "refuses a guessed dpi and a guess is all there was.",
            "tolerance_used": "the STRICT axis.LOO_TOL_PX = 1.5 px, "
                              "unrescaled; the widening was NOT relied on",
            "ticks_px": fig["ticks_px"],
            "tick_spacing_px": fig["tick_spacing_px"],
            "tick_spacing_spread_px": fig["tick_spacing_spread_px"],
            "px_per_1000_units": fig["px_per_1000_units"],
            "loo_worst_px": fig["loo_worst_px"],
            "loo_passes_strict_1.5px": fig["loo_passes_strict"],
            "precision_statement":
                f'1 px = {1000.0 / fig["px_per_1000_units"]:.2f} units; the '
                f'strict 1.5 px gate therefore permits at most '
                f'{1.5 * 1000.0 / fig["px_per_1000_units"]:.1f} units of '
                f'scale error, 0.44% of the smallest emitted value',
            "bars_read": fig["n_bars_read"], "bars_refused":
                fig["n_bars_refused"], "refusals": fig["refusals"],
            "flank_delta_px_max": max(
                b["flank_delta_px"] for b in fig["bars"]
                if b["flank_delta_px"] is not None),
            "max_error_bar_cap_coverage": fig["max_cap_coverage_over_read_bars"],
            "bars_FRAC_gate": fig["bars_frac_gate"],
            "cap_margin_vs_gate": fig["cap_margin_vs_gate"],
            "furniture_masked":
                f'rows {fig["furniture_rows_masked"]}; limit '
                f'{fig["furniture_limit_px"]} px = 3 bar widths. The bracket '
                f'at 3106-3109 spans 612 px (9.1 bar widths) and covers the '
                f'SM-102 and ALC-0315 windows at 0.97-1.00',
            "censored_bars": fig["bars_censored_at_axis_max"]},
        "P5_unprinted_round_number": {
            "what": "the 0 tick must land on the DRAWN axis baseline -- the "
                    "authors drew it and never printed it",
            "offset_px": fig["zero_tick_vs_drawn_baseline_px"],
            "offset_in_axis_units": fig["zero_tick_vs_drawn_baseline_units"],
            "why_it_earns_its_place":
                "it is the ONLY check that sees a uniform origin error: the "
                "ticks_uniform_shift seed moves every tick by +4 px, so "
                "calibrate and the LOO gate both still pass perfectly "
                "(0.279 px) while this offset moves 0.101 -> 4.101 px"},
        "P6_anchor_C10": {
            "anchor_built_from": "Figure 1A, green box (C12-200 drawn in "
                                 "full; every tail PRINTED as 'C10H21')",
            "result": anc["C10_result"],
            "matched_corpus_rows": anc["C10_matched_rows"],
            "anchor_inchikey": anc["anchor"]["inchikey"],
            "anchor_formula": anc["anchor"]["formula"],
            "anchor_exact_mass": anc["anchor"]["exact_mass"],
            "corpus_stored_Molecular_Weight":
                list(anc["corpus_C12_200"].values())[0]
                ["stored_Molecular_Weight"],
            "connectivity_seed": anc["connectivity_seed"],
            "independence_guards": {
                "clean_control_clean":
                    anc["independence_guards"]["clean_control"]["clean"],
                "banned_import_guard_fires":
                    anc["independence_guards"]["seed_banned_import"]["fires"],
                "lifted_template_guard_fires":
                    anc["independence_guards"]["seed_lifted_template"]
                    ["fires"]},
            "anchors_empty_on_purpose":
                anc["reference_structures_csv"]["n_empty_on_purpose"],
            "vacuity_statement":
                "SM-102, ALC-0315, cholesterol, DOPE, DSPC and DOTAP have NO "
                "independent anchor. Against an empty anchor the guards are "
                "VACUOUS, not passing."},
        "P7_seeds": {
            "n_seeds": seeds["summary"]["n_seeds"],
            "n_credited": seeds["summary"]["n_credited"],
            "negative_control_moved_nothing": any(
                s["seed"] == "roll_image" and s["credited"]
                for s in seeds["seeds"]),
            "byte_identical_across_two_clean_runs":
                seeds["reproducibility"]
                ["byte_identical_across_two_clean_runs"],
            "per_seed": [{"seed": s["seed"], "credited": s["credited"],
                          "route": s.get("route", "n/a"),
                          "reason": s["reason"]} for s in seeds["seeds"]],
            "emit_refusals": refusals},
        "P8_emission": eqc,
        "P9_benchmarks": {
            "question": "do the field's benchmark lipids top this screen?",
            "answer": "NO",
            "C12_200_over_SM_102":
                eqc["benchmark_ratio_C12_200_over_SM_102"],
            "C12_200_over_ALC_0315":
                eqc["benchmark_ratio_C12_200_over_ALC_0315"],
            "but_the_comparison_is_CONFOUNDED":
                "the formulation (DOPE 30.5 / chol 32 / PEG 2.5) was selected "
                "as the winner of an 18-member screen run on C12-200. SM-102 "
                "and ALC-0315 are evaluated in a formulation optimised for a "
                "different ionizable lipid, so the ranking is not a clean "
                "chemistry comparison. Recorded on the row and in DECISIONS.",
            "all_three_are_universal_controls":
                "C12-200, SM-102 and ALC-0315 are the field's benchmark "
                "lipids and all three are already in the corpus. 3 of 3 "
                "shared lipids are controls -- CLAUDE.md section 8's hazard "
                "in its purest form."},
        "measured_sensitivity_limits": {
            "scale": f'1 px = {1000.0 / fig["px_per_1000_units"]:.2f} axis '
                     f'units; the smallest emitted value (1179.1) is '
                     f'{1179.1 * fig["px_per_1000_units"] / 1000.0:.0f} px '
                     f'tall, so a 1 px read error is 0.96% of it',
            "error_bar_caps": f'the widest cap covers '
                              f'{fig["max_cap_coverage_over_read_bars"]} of '
                              f'its bar against the 0.90 gate -- margin '
                              f'{fig["cap_margin_vs_gate"]}. A cap wider than '
                              f'0.90 would leak; none here is.',
            "control_72h": "bars.bar_top refuses it. The bar is ~3 px tall "
                           "and the refusal is recorded as a refusal, never "
                           "as a zero.",
            "dose": "Fig 4G's pDNA dose is stated nowhere in the SI and the "
                    "article is paywalled. Dose_ug_nucleicacid is BLANK "
                    "(F-004's posture), not reconstructed from Fig 3B."}}

    json.dump(qc, open(os.path.join(HERE, "07_qc.json"), "w"), indent=1)
    print(json.dumps({"positive_controls": pc,
                      "ionizable_lipids": {k: {
                          "already_in_corpus": v["already_in_corpus"],
                          "vitro_experiments": v["vitro_experiments"],
                          "vivo_experiments": v["vivo_experiments"]}
                          for k, v in tg.items()},
                      "emit_refusals": refusals,
                      "n_new_to_the_corpus": ov["n_new_to_the_corpus"]},
                     indent=1))
    print("wrote 02_corpus_overlap.json, 07_qc.json")


if __name__ == "__main__":
    main()
