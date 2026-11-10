#!/usr/bin/env python3
"""Consolidate every check in this curation into 07_qc.json, WITH NUMBERS."""
import csv, json, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
L = lambda n: json.load(open(os.path.join(HERE, n)))
dig, vq, ov = L("04_figS1S2_measurements.json"), L("05_values_qc.json"), L("02_corpus_overlap.json")
ref, seeds = L("reference_report.json"), L("07_seed_results.json")
enum = list(csv.DictReader(open(os.path.join(HERE, "03_enumerated.csv"))))
rows = list(csv.DictReader(open(os.path.join(HERE, "06_rows.csv"))))
bars = {(b["lnp"], b["organ"]): b["value"] for b in dig["bars"] if b["panel"] == "S2A"}
L8 = ["13", "14", "70", "96", "97", "98", "45", "102"]
ORG = ["Brain", "Heart", "Lungs", "Liver", "Kidneys", "Spleen"]

qc = {"_paper": "Han et al., ACS Nano 2026, PMID 41555650, DOI 10.1021/acsnano.5c19138",
      "_lab": "Michael J. Mitchell (UPenn). The extraction brief misattributed "
              "this to the Dahlman lab; the folder was corrected to "
              "curation/mitchellmj_41555650.",
      "_endpoint": "mRNA delivery. Pooled mRNA-BARCODE screen (main article) + a "
                   "separate SINGLE-LNP IVIS luciferase validation (Fig 3, SI "
                   "Figs S1/S2/S5). Only the latter is corpus-admissible."}

# ---- P1 amine-core stoichiometry
from rdkit import Chem
from rdkit.Chem.Descriptors import MolWt
CORES = {"200": "NCCN1CCN(CCNCCN)CC1", "494": "NCCOCCN1CCN(CCNCCOCCN)CC1",
         "306": "CN(CCCN)CCCN", "114": "NCCNCCN",
         "490": "NCC(OCC)CN1CCN(CCNCC(OCC)CN)CC1"}
MM = {"200": (26.0, 0.121, 0.724, 6.0), "494": (40.0, 0.132, 0.791, 6.0),
      "306": (43.6, 0.300, 2.10, 7.0), "114": (30.9, 0.300, 2.10, 7.0),
      "490": (25.0, 0.0750, 0.528, 7.0)}
p1 = {}
for k, s in CORES.items():
    mw = MolWt(Chem.MolFromSmiles(s))
    mg, mmol_pr, emmol, eq = MM[k]
    p1[k] = dict(derived_formula_MW=round(mw, 3),
                 MW_from_printed_amine_mmol=round(mg / mmol_pr, 2),
                 dev_pct_printed=round(100 * ((mg / mmol_pr) - mw) / mw, 3),
                 MW_from_epoxide_equivalents=round(mg / (emmol / eq), 2),
                 dev_pct_via_equivalents=round(100 * ((mg / (emmol / eq)) - mw) / mw, 3))
qc["P1_amine_core_pinned_by_the_papers_own_stoichiometry"] = dict(
    per_core=p1,
    worst_dev_pct_printed=round(max(abs(v["dev_pct_printed"]) for v in p1.values()), 3),
    worst_dev_pct_via_equivalents=round(max(abs(v["dev_pct_via_equivalents"]) for v in p1.values()), 3),
    what=("Each amine core is DERIVED from the SI systematic name, then its "
          "formula is tested against the mass and moles the SI prints for that "
          "amine. Amine 490's printed mmol (0.0750) gives +0.55%, an outlier; "
          "the SI's own epoxide line gives 0.528 mmol at 7.0 equivalents, i.e. "
          "a true amine charge of 0.0754 mmol, and the deviation collapses to "
          "-0.02%. The equivalents column is a SECOND, independent statement "
          "of the same quantity and it is the one that is self-consistent."),
    SEEDED=[s for s in seeds if s["check"].startswith("P1")][0])

# ---- P2 epoxide reagent MW
EPMW = {"C10": 156.26, "C12": 184.32, "C14": 212.37}
EPMG = {("C10", "200"): 113.0, ("C12", "200"): 332.0, ("C14", "200"): 154.0,
        ("C10", "494"): 124.0, ("C12", "494"): 146.0, ("C14", "494"): 168.0,
        ("C10", "306"): 328.0, ("C12", "306"): 387.0, ("C14", "306"): 446.0,
        ("C10", "114"): 328.0, ("C12", "114"): 387.0, ("C14", "114"): 446.0,
        ("C10", "490"): 82.5, ("C12", "490"): 97.3, ("C14", "490"): 112.0}
EPMMOL = {"200": {"C10": 0.724, "C12": 1.80, "C14": 0.724},
          "494": {c: 0.791 for c in EPMW}, "306": {c: 2.10 for c in EPMW},
          "114": {c: 2.10 for c in EPMW}, "490": {c: 0.528 for c in EPMW}}
devs = []
for (t, a), mg in EPMG.items():
    devs.append(abs(100 * ((mg / EPMMOL[a][t]) - EPMW[t]) / EPMW[t]))
qc["P2_epoxide_reagent_MW_from_printed_mass_and_moles"] = dict(
    n=len(devs), worst_dev_pct=round(max(devs), 3),
    median_dev_pct=round(float(np.median(devs)), 3),
    what=("All 15 epoxide charges reproduce the 2-alkyloxirane molecular "
          "weight, confirming the TAIL identity independently of the lipid name."))

# ---- P3 arm count
agree = sum(1 for r in enum if r["n_arms_built"] == r["n_arms_in_iupac_name"])
qc["P3_arms_built_vs_arms_in_IUPAC_name"] = dict(
    n=len(enum), n_agree=agree,
    per_core={r["amine_core"]: dict(core_n_NH=r["core_n_NH"],
                                    built=r["n_arms_built"],
                                    name_says=r["n_arms_in_iupac_name"])
              for r in enum},
    what=("The builder attaches one arm per remaining N-H and COUNTS the "
          "result; the IUPAC name's 'tetrakis'/'bis' is transcribed separately "
          "and never read by the builder. Agreement is evidence."),
    SEEDED=[s for s in seeds if s["check"].startswith("P3")][0])

# ---- P4-P7 tables
for k in ("P4_block_map_vs_printed_labels",
          "P5_design_grid_matches_printed_molar_ratios",
          "P6_helper_and_sterol_constant_across_library",
          "P7_molar_percents_sum_to_100"):
    qc[k] = vq[k]
qc["P5_design_grid_matches_printed_molar_ratios"]["SEEDED"] = \
    [s for s in seeds if s["check"].startswith("P5")][0]

# ---- P8-P12 geometry
qc["P8_axis_calibration"] = {k: dig["panels"][k] for k in dig["panels"]}
qc["P8_axis_calibration"]["SEEDED"] = [s for s in seeds if s["check"].startswith("P8")][0]
qc["P10_P11_P12_cross_panel"] = dig["cross_panel"]

# ---- benchmark question
bt = {}
for o in ORG:
    rank = sorted(L8, key=lambda l: -bars[(l, o)])
    bt[o] = dict(order=rank, top=rank[0],
                 benchmark_SM102_rank=rank.index("102") + 1,
                 SM102_over_best_novel=round(bars[("102", o)] / bars[(rank[1] if rank[0] == "102" else rank[0], o)], 4))
qc["benchmark_tops_the_screen"] = dict(
    per_organ=bt,
    answer=("YES, on absolute delivery. SM-102 (b-LNP 102) is ranked 1st of 8 "
            "in ALL SIX organs, including the brain, where it is 2.091x the "
            "best novel formulation (b-LNP 70, C12-200)."),
    but=("The paper's claim is about SELECTIVITY, not absolute brain signal. "
         "Brain/liver ratio: b-LNP 14 (C14-306) = 0.0271 vs SM-102 = 0.0020, a "
         "13.6x advantage. So the lead wins on the ratio while losing on every "
         "raw organ value."),
    brain_over_liver={l: round(bars[(l, "Brain")] / bars[(l, "Liver")], 5) for l in L8},
    why_this_matters=("CLAUDE.md section 8 has an open decision on whether "
                      "benchmark-only pairs belong in Stage 5's pooled "
                      "concordance. This is a case where the universal positive "
                      "control tops the screen outright in every organ, which "
                      "is the pattern that can push cross-lab concordance "
                      "toward +1 for reasons unrelated to transferable SAR."))

# ---- corpus overlap
qc["corpus_overlap"] = dict(
    summary=ov["exact_structure_match_summary"],
    exact_hits=[{k: e[k] for k in ("IL_name", "InChIKey", "corpus_rows_skeleton",
                                   "corpus_experiments", "corpus_names")}
                for e in ov["exact_structure_match"] if e["full_inchikey_match"]],
    near_miss_same_amine_core={k: {kk: v[kk] for kk in
                                   ("near_miss_unique_structures_excluding_ours",
                                    "near_miss_rows", "near_miss_experiments")}
                               for k, v in ov["amine_core_near_miss_sweep"].items()},
    superset_sweep=ov["superset_sweep_2hydroxyalkyl_amine"])

# ---- anchor
qc["P13_C10_anchor"] = ref
qc["P13_C10_anchor"]["SEEDED"] = [s for s in seeds if s["check"].startswith("C10")][0]

# ---- emitted rows
qc["emitted_rows"] = dict(
    n_rows=len(rows),
    n_distinct_IL_structures=len({r["IL_InChIKey"] for r in rows}),
    n_distinct_formulations=len({r["b_LNP"] for r in rows}),
    n_organs=len({r["Model_target"] for r in rows}),
    confidence_counts={c: sum(1 for r in rows if r["confidence"] == c)
                       for c in ("high", "medium", "low", "")},
    rows_with_a_caveat_F024=sum(
        1 for r in rows if r["confidence"] in ("medium", "low")
        or r["positional_ambiguity"] == "1" or r["structure_assumed"] == "1"
        or r["value_is_bound"] in ("upper", "lower")),
    z_mean=round(float(np.mean([float(r["Experiment_value"]) for r in rows])), 12),
    z_sd=round(float(np.std([float(r["Experiment_value"]) for r in rows], ddof=0)), 9),
    SEEDED_emit=[s for s in seeds if s["check"].startswith("emit")
                 or s["check"].startswith("conf_rule")])

qc["P14_pool_membership_from_text_COLOUR"] = {
    k: v for k, v in L("09_pool_membership.json").items() if k != "per_bLNP"}
qc["pooled_screen_admissibility"] = dict(
    corrected=True,
    first_verdict="inadmissible by design (a pooled value is competitive)",
    corrected_verdict=("ADMISSIBLE IN KIND. data/data_vivo.csv has an "
                       "Experiment_batching column, vocabulary "
                       "{barcoded, individual}; LR_2023 contributes 407 "
                       "barcoded in-vivo rows over 383 unique ILs "
                       "(protein_adundance_normalized, liver). Pooling alone "
                       "is not disqualifying -- the corpus declares it on the "
                       "row. The binding constraint is ACQUISITION."),
    pooled_values_emitted=0,
    pooled_values_emitted_because="main-article Figure 2 was not acquired",
    normalisation_residual_computed=False,
    normalisation_residual_not_computed_because=(
        "the formula, the worked example and the values are all in the "
        "unacquired main article"))
qc["sibling_same_lab_cross_check"] = dict(
    sibling="curation/mitchellmj_38982196 (Wang et al., Nat Chem 2024, MW_2024)",
    same_lab=True, n_sibling_structures=100,
    full_inchikey_overlap=0, skeleton_inchikey_overlap=0,
    sibling_chemistry="amidine-incorporated degradable lipids; thioether-ester "
                      "tails (e.g. CCCCCCCCCCNC(=N)CCCSCCC(=O)OCCCCCCCCC)",
    sibling_naming_grammar="<amine index>T-O<tail index>, e.g. 1T-O9 .. 25T-O16",
    finding=("ZERO shared chemistry despite a shared lab, and the two amine "
             "catalogues are DIFFERENT AND LOCALLY SCOPED: MW_2024's amines are "
             "numbered 1-25 from its own building-block set, while this paper's "
             "200/306/114/490/494 are the Anderson/Akinc lipidoid catalogue "
             "shared with KL_2010, KW_2014 and AC_2025. Same lab is not same "
             "numbering -- a leading token is not an identity."))
qc["value_counts"] = vq["value_counts"]
qc["seed_summary"] = dict(
    n_checks_seeded=len(seeds),
    n_credited=sum(1 for s in seeds if s["CREDITED"]),
    rule=("A check is credited only when the AGREEMENT COUNT MOVES under the "
          "seed, or a clean run exits 0 where the seeded run exits non-zero. "
          "One seed (P1) initially exited non-zero from an unrelated "
          "TypeError and was briefly credited; 07_seed_checks.py now refuses "
          "to credit a seed that did not reach its own metric."))

qc["licence"] = dict(
    article=dict(journal="ACS Nano (American Chemical Society)",
                 isOpenAccess="N", europepmc_license=None,
                 fullTextUrl_availability="Subscription required",
                 verdict="NOT CC0/CC-BY. Closed, subscription. CLAUDE.md "
                         "section 8 requires CC0/CC-BY compatibility before "
                         "corpus RELEASE; this source does not have it, and "
                         "that is a release-time question recorded in "
                         "DECISIONS.md D-6."),
    si_pdf=dict(note="Downloaded through an entitled institutional session "
                     "(watermark: PHILLIPS EXETER ACADEMY, 10 November 2026) "
                     "before this run began. No separate licence is stated on "
                     "the file."),
    repository=dict(exists=False, searched="provenance/tools/find_repo.py, "
                                           "DOI + PMID code search + 3 repo-name "
                                           "queries, 0 hits"))

json.dump(qc, open(os.path.join(HERE, "07_qc.json"), "w"), indent=1)
print("P1 worst dev (printed / via equivalents): "
      f"{qc['P1_amine_core_pinned_by_the_papers_own_stoichiometry']['worst_dev_pct_printed']}% / "
      f"{qc['P1_amine_core_pinned_by_the_papers_own_stoichiometry']['worst_dev_pct_via_equivalents']}%")
print("P2 epoxide worst dev:", qc["P2_epoxide_reagent_MW_from_printed_mass_and_moles"]["worst_dev_pct"], "%")
print("P3:", agree, "/", len(enum))
print("P7 sums to 100:", qc["P7_molar_percents_sum_to_100"]["n_exactly_100"], "/ 93")
print("values:", qc["value_counts"]["total"], "admissible", qc["value_counts"]["admissible"])
print("seeds credited:", qc["seed_summary"]["n_credited"], "/", qc["seed_summary"]["n_checks_seeded"])
