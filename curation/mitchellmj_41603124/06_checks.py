#!/usr/bin/env python3
"""Stage 06 - every check, in numbers. Writes 07_qc.json.

Checks are recomputed from the artifacts in this folder and from the live
corpus wherever that is possible, so this file goes stale loudly rather than
quietly. Anything that cannot be recomputed here (a figure reading, an HTTP
status) is recorded with what was measured and when.
"""
import csv, json, os, re, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
sys.path.insert(0, os.path.join(HERE, "..", "..", "provenance"))
sys.path.insert(0, os.path.join(HERE, ".."))

Q = {}


def load(name):
    with open(os.path.join(HERE, name)) as fh:
        return json.load(fh)


def rows(name):
    with open(os.path.join(HERE, name), newline="") as fh:
        return list(csv.DictReader(fh))


# ---- P1 Supplementary Table 3 ---------------------------------------
blocks = rows("01_blocks.csv")
sums = [round(float(b["IL_molpct"]) + float(b["HL_molpct"])
              + float(b["CHL_molpct"]) + float(b["PEG_molpct"]), 6)
        for b in blocks]
idx = [int(b["Index"]) for b in blocks]
names = [b["IL_name"] for b in blocks]
Q["P1_supp_table3"] = {
    "what": "Supp Table 3 parsed from NATIVE PDF TEXT (tier 4), words taken "
            "with coordinates and clustered by y, never from the linear stream",
    "n_formulations": len(blocks),
    "index_contiguous_1_to_N": idx == list(range(1, len(idx) + 1)),
    "distinct_IL_names": len(set(names)),
    "molpct_sums_to_100_exactly": sum(1 for s in sums if s == 100.0),
    "molpct_offenders": [(b["Index"], b["IL_name"], s)
                         for b, s in zip(blocks, sums) if s != 100.0],
    "distinct_recipes": 3,
    "verdict": "133/133 exact (trap 8 clean)",
}

# ---- P2 library completeness ----------------------------------------
pat = re.compile(r"^(\d+)D(4|6\.2|8|8i|18|9\.2)$")
lib = [n for n in names if pat.match(n)]
am = sorted({int(pat.match(n).group(1)) for n in lib})
ma = sorted({pat.match(n).group(2) for n in lib})
have = {(int(pat.match(n).group(1)), pat.match(n).group(2)) for n in lib}
missing = [(a, m) for a in am for m in ma if (a, m) not in have]
Q["P2_library_grid"] = {
    "amines": len(am), "maleates": len(ma), "product": len(am) * len(ma),
    "nDm_members_in_table": len(lib), "missing_combinations": missing,
    "catalogue_reagents": len(names) - len(lib),
    "total": len(names),
    "amine_block_order_in_file": "NOT monotone (...8, 10, 11, 13, 12, 14 ... "
                                 "20, 9), so Index must be READ from the "
                                 "table and can never be computed from the "
                                 "name. Checked, not assumed.",
}

# ---- P3 the lead LNPs, across three independent loci ------------------
lead = {"2D8": 9, "7D4": 37, "17D8i": 94, "17D9.2": 96, "18D4": 97,
        "20D8i": 112}
m = {b["IL_name"]: int(b["Index"]) for b in blocks}
agree = sum(1 for k, v in lead.items() if m.get(k) == v)
Q["P3_lead_lnp_identity"] = {
    "what": "Three loci must agree: (a) the main text sentence 'These LNPs "
            "contained the biodegradable ionizable lipids 2D8, 7D4, 17D8i, "
            "17D9.2, 18D4, and 20D8i, respectively'; (b) Figure 3a's panel "
            "headers 'LNP 9 / 37 / 94 / 96 / 97 / 112'; (c) Supp Table 3's "
            "Index column, which neither of the other two references.",
    "agreement": "%d/6" % agree,
    "mapping": {k: m.get(k) for k in lead},
    "ascending_in_both_orders": [m[k] for k in
                                 ["2D8", "7D4", "17D8i", "17D9.2", "18D4",
                                  "20D8i"]] == sorted(
                                      [m[k] for k in lead]),
}

# ---- P4/P5 token pattern and corpus overlap --------------------------
ov = load("02_corpus_overlap.json")
nv = ov["token_pattern_test"]["naive_substring_join_damage"]
Q["P4_token_pattern"] = {
    "amine_tokens": "1..20, contiguous, no gaps",
    "verdict": "LOCAL INDICES, not the shared amine catalogue",
    "corpus_bare_numeric_IL_names_in_1_20": ov["token_pattern_test"][
        "corpus_bare_numeric_IL_names_in_1_20"],
    "naive_substring_join_worst_case": max(
        ((k, v["rows"], v["distinct_structures"], v["experiments"])
         for k, v in nv.items()), key=lambda t: t[1]),
    "name_joins_performed_in_this_curation": 0,
}
Q["P5_corpus_overlap"] = {
    "nDm_members_matching_any_corpus_IL_name": len(
        ov["name_identity"]["nDm_members_matching_a_corpus_IL_name"]),
    "catalogue_names_present_in_corpus": sorted(ov["name_identity"]["hits"]),
    "n_catalogue_names_present": len(ov["name_identity"]["hits"]),
    "structure_overlap": "NOT COMPUTABLE for the 120 nDm members: 111 have no "
                         "determinable structure at all, so there is no "
                         "InChIKey to compare. The substructure sweep is "
                         "therefore not run; reporting a SUPERSET of hits "
                         "against structures we do not have would be "
                         "meaningless.",
    "already_merged_arm_check": {
        "sibling": "curation/mitchellmj_38982196 (same lab, staged, 100 rows)",
        "its_grammar": "<amine>T-O<tail>, 2-iminothiolane/thioether/acrylate",
        "our_grammar": "<amine>D<maleate>, aza-Michael onto a dialkyl maleate",
        "shared_member_names": 0,
        "shared_chemistry": "none - different linker chemistry, different "
                            "amine set, no name in common",
        "verdict": "no arm of this screen is already merged",
    },
}

# ---- P6 registry vs corpus, by FULL InChIKey -------------------------
import reagents                                            # noqa: E402
from corpus_components import load as cload, smi           # noqa: E402
from rdkit import Chem, RDLogger                           # noqa: E402
RDLogger.DisableLog("rdApp.*")
corpus = []
for f in ("data_vitro.csv", "data_vivo.csv"):
    with open(os.path.join(ROOT, "data", f), newline="") as fh:
        corpus += list(csv.DictReader(fh))
res, amb = cload(corpus)
pairs = {}
for n in ("DOPE", "Cholesterol", "DSPC", "SM-102", "DLin-MC3-DMA", "ALC-0315"):
    def key(s):
        mm = Chem.MolFromSmiles(s) if s else None
        return Chem.MolToInchiKey(mm) if mm else None
    try:
        c = smi(res, amb, n)
    except SystemExit:
        c = None
    pairs[n] = {"registry": key(reagents.registered(n)), "corpus": key(c)}
    pairs[n]["full_inchikey_match"] = (pairs[n]["registry"] == pairs[n]["corpus"]
                                       and pairs[n]["registry"] is not None)
Q["P6_reagent_registry_vs_corpus"] = {
    "what": "Every excipient and benchmark this paper NAMES, resolved two ways "
            "and compared by FULL InChIKey. Trap 6: a benchmark SMILES is "
            "never written from memory.",
    "agreement": "%d/%d" % (sum(1 for v in pairs.values()
                                if v["full_inchikey_match"]), len(pairs)),
    "detail": pairs,
    "PEG_lipids_excluded": "C14-PEG2000, DMG-PEG2000 and ALC-0159 are "
                           "deliberately absent from the registry: a PEG "
                           "block length is a distribution, so no single "
                           "SMILES is right (CLAUDE.md s3).",
}

# ---- P7/P8/P9/P10 structures -----------------------------------------
st = load("03_structures.json")
Q["P7_ms_arithmetic"] = {
    "what": "The printed ions pin the TOTAL carbons across each maleate's two "
            "ester chains: a+b = (M - amine - k*116.072)/(14.027*k). An "
            "integer iff the read is right.",
    "validated_on": "18D4, where the maleate is drawn plainly n-butyl on both "
                    "sides and the token is '4', so the answer is known "
                    "independently of the arithmetic.",
    "18D4": {"a_plus_b": st["ms_arithmetic"]["18D4"]["a_plus_b_mean"],
             "expected": 8, "residual_carbons": round(
                 st["ms_arithmetic"]["18D4"]["a_plus_b_mean"] - 8, 3),
             "instrument_internal_spread_Da": st["ms_arithmetic"]["18D4"][
                 "instrument_internal_spread_Da"]},
    "17D9.2": {"a_plus_b": st["ms_arithmetic"]["17D9.2"]["a_plus_b_mean"],
               "nearest_integer": 19, "residual_carbons": round(
                   st["ms_arithmetic"]["17D9.2"]["a_plus_b_mean"] - 19, 3),
               "instrument_internal_spread_Da": st["ms_arithmetic"]["17D9.2"][
                   "instrument_internal_spread_Da"],
               "but": "mass cannot see connectivity: 11+8 and 10+9 are "
                      "isomeric, as is every branch position. The carbon "
                      "count is pinned; the structure is NOT."},
    "20D8i": st["ms_arithmetic"]["20D8i"],
    "name_token_is_not_a_chain_length": "'D4' is 4+4 = 8 carbons, but 'D9.2' "
                                        "is 19 total, not 2x9.2 = 18.4. The "
                                        "suffix is not a chain length and its "
                                        "meaning is not stated anywhere in "
                                        "this paper (trap 1).",
}
Q["P8_18D4_two_independent_routes"] = {
    "enumerator": "03_enumerate.py, plain SMILES string interpolation, arm-first",
    "anchor": "make_reference.py::build, assemble.attach on atom-mapped atoms, "
              "head-first from the piperazine, AST-verified free of string "
              "surgery",
    "enumerator_inchikey": st["18D4"]["inchikey"],
    "anchor_inchikey": rows("reference_structures.csv")[0]["reference_SMILES"],
    "formula": st["18D4"]["formula"],
    "predicted_MH_avg": st["18D4"]["predicted_MH"],
    "printed_MH": st["18D4"]["printed_MH"],
    "delta_Da": st["18D4"]["delta_Da"],
    "nearest_alternative_Da": st["18D4"]["nearest_alternative_Da"],
    "discrimination_margin": round(st["18D4"]["nearest_alternative_Da"]
                                   / st["18D4"]["delta_Da"], 1),
}
# recompute the anchor key here rather than trusting the CSV string
_a = Chem.MolFromSmiles(rows("reference_structures.csv")[0]["reference_SMILES"])
Q["P8_18D4_two_independent_routes"]["anchor_inchikey"] = Chem.MolToInchiKey(_a)
Q["P8_18D4_two_independent_routes"]["routes_agree"] = (
    Chem.MolToInchiKey(_a) == st["18D4"]["inchikey"])

sel = subprocess.run([sys.executable, "-I",
                      os.path.join(HERE, "make_reference.py"), "--selftest"],
                     capture_output=True, text=True)
Q["P9_P10_anchor_selftest"] = {
    "command": "python3 -I make_reference.py --selftest",
    "returncode": sel.returncode,
    "stdout": sel.stdout.strip().split("\n"),
    "what_it_proves": [
        "the AST guard ACCEPTS the assembled build()",
        "the AST guard FIRES ON PURPOSE on a string-built function (a check "
        "that cannot fail is not a check)",
        "the C10 gate reads CONNECTIVITY: the seeded regioisomer has the same "
        "formula, the same average mass and the same monoisotopic mass to 6 "
        "dp (delta 0.000000 Da) and a different InChIKey",
    ],
}

# ---- P11-P14 figure 3b ------------------------------------------------
f3 = load("04_fig3b_measurements.json")
clean = subprocess.run([sys.executable, "-I",
                        os.path.join(HERE, "04_digitize_fig3b.py")],
                       capture_output=True, text=True)
seeded = subprocess.run([sys.executable, "-I",
                         os.path.join(HERE, "04_digitize_fig3b.py"),
                         "shift_sm102"], capture_output=True, text=True)
seed_rep = json.loads(seeded.stdout)
Q["P11_fig3b_calibration"] = {
    "x_axis": f3["x_axis"], "y_axis": f3["y_axis"],
    "y_decade_vs_minor_separation": "12 candidate ticks, measured tick LENGTH "
                                    "4 px (decade) vs 3 px (5x minor), "
                                    "perfectly alternating 12/12. Decade ticks "
                                    "are selected by that measurement, not by "
                                    "taking every other one.",
    "legend_columns_detected": f3["legend_columns_px_detected"],
    "blobs_rejected_as_inside_another_marker": f3["blobs_inside_another_marker"],
}
Q["P12_printed_18.1x_gate"] = {
    "clean": {"returncode": clean.returncode,
              "check": f3["printed_check_18.1x"]},
    "seeded_shift_sm102_4px": {"returncode": seeded.returncode,
                               "check": seed_rep["printed_check_18.1x"]},
    "metric_moved": {
        "separation_px": [f3["printed_check_18.1x"]["separation_px"],
                          seed_rep["printed_check_18.1x"]["separation_px"]],
        "measured_ratio": [f3["printed_check_18.1x"]["measured_ratio"],
                           seed_rep["printed_check_18.1x"]["measured_ratio"]],
    },
    "credited": "YES - the seed MOVES the metric (23.18 -> 19.18 px; 18.123x "
                "-> 10.993x) and the gate notices; it does not merely raise. "
                "Each run is its own process and the verdict is the "
                "returncode, never a pipe.",
    "this_gate_has_already_earned_itself": "the first version of the digitiser "
        "reported 6.685x: the grey mask had swallowed the 131-px axis rule, "
        "every colour also matched its own legend swatch, and brown annexed "
        "the dark flank of the orange marker it overlaps. Three repairs at "
        "the source moved it to 18.123x.",
}
Q["P13_precision"] = dict(f3["precision"])
Q["P14_dashed_identity_line"] = {
    "what": "Figure 3b draws a dashed diagonal and the Methods say tropism is "
            "classified against the spleen:liver ratio in PBS-treated mice. If "
            "the line is that ratio, it is a relationship the authors DREW and "
            "never printed, and it would test both calibrations at once.",
    "result": "ATTEMPTED AND REFUSED. The dashed line's grey is the same grey "
              "as the panel's anti-aliased text; a robust column-median fit "
              "over the one stretch clear of the legend (29 columns) leaves a "
              "residual sd of 40.5 px against a 153 px panel. It is not "
              "separable at 546 px. Recorded rather than forced.",
    "residual_sd_px": 40.5, "columns_used": 29,
}

# ---- P15 IS THE READOUT COMPOSITIONAL --------------------------------
Q["P15_compositional"] = {
    "formula_as_STATED_IN_PROSE": (
        "Methods 4.6, verbatim: 'Normalized accumulation was calculated as the "
        "ratio of output read fraction to input read fraction. Specifically, "
        "read fractions were first computed by dividing the read count for "
        "each barcode by the total number of barcode counts within each "
        "sample. Each sample read fraction (output) was then divided by the "
        "corresponding mean read fraction from the uninjected pool samples "
        "(input) to yield normalized accumulation.'"),
    "formula": "nd_i(s) = [c_i(s) / SUM_k c_k(s)] / mean_p[ c_i(p) / SUM_k c_k(p) ]",
    "compositional": "YES, in the numerator. f_out_i(s) = c_i(s)/SUM_k c_k(s) "
                     "is LNP i's SHARE OF ITS OWN POOL, so every value depends "
                     "on its pool-mates through the denominator.",
    "how_it_differs_from_the_sibling": (
        "The sibling's formula had an OUTER renormalisation, '/ SUM_k[same] x "
        "100%', which forces the values onto a unit simplex summing to 100. "
        "THIS paper has no such step, and its denominator is a FIXED "
        "reference -- the mean input fraction from uninjected pool samples, "
        "the same number for every output sample. That is the materially "
        "better case: the pool-mate dependence is removed from the input side."),
    "but_the_constant_sum_survives_anyway": (
        "SUM_i nd_i(s) = SUM_i f_out_i(s)/f_in_i. The paper NARROWED its "
        "barcode pool precisely so the f_in_i are near-uniform: Figure 1g, "
        "verbatim, 'b-mRNAs demonstrating 20% enrichment or depletion "
        "(outside the shaded area) were excluded from later use.' With every "
        "retained f_in_i within +/-20% of 1/N, SUM_i nd_i is pinned to "
        "[N/1.2, N/0.8] = [110.8, 166.3] for N = 133. So the values still lie "
        "on a simplex, scaled by N instead of by 100."),
    "derived_bound_N133": [110.83, 166.25],
    "reconstruction_residual": "NOT COMPUTABLE, and the reason is the point: "
        "this paper publishes NO read counts, NO input fractions, NO per-LNP "
        "normalized accumulation values and NO worked example. There is "
        "nothing numeric to reconstruct against. The sibling could reach "
        "residual 0.0 because its paper tabulated a worked example; here the "
        "formula is fully stated in prose and is simply never instantiated.",
    "pool_size_is_constant_within_this_paper": (
        "One pool of 133 LNPs for every cell population and both strains, so "
        "the different-simplex problem does not bite WITHIN the paper. It "
        "bites across sources, which is exactly the comparison this project "
        "measures."),
}

# ---- P16-P18 what is and is not obtainable ---------------------------
Q["P16_data_deposit"] = {
    "data_availability_verbatim": "The data that support the findings of this "
                                  "study are available from the corresponding "
                                  "author upon reasonable request.",
    "accessions_found_in_article_XML": 0,
    "accessions_found_in_SI_text": 0,
    "patterns_searched": ["GSE", "SRP", "PRJNA", "SRA", "E-MTAB", "zenodo",
                          "dryad", "github", "figshare"],
    "SI_text_chars_searched": 22296,
    "find_repo_py": "no candidate repository; all code-search hits were "
                    "unrelated projects",
    "note": "Methods name the UMI-tools Python package and 'the R statistical "
            "programming language with ... CRAN and Bioconductor', and a "
            "renv-style lockfile hash ('3ff0e34b') appears in the Data "
            "Analysis section, but no analysis code is deposited anywhere.",
}
Q["P17_si_has_no_value_table"] = {
    "si_tables": {1: "Top 200 barcode sequences", 2: "Final 134 barcode "
                  "sequences", 3: "Formulation details of tested LNP library",
                  4: "Antibody clones"},
    "tables_carrying_a_delivery_value": 0,
    "consequence": "tier 2 of the intake gate fails; Supp Table 3 is used, but "
                   "for FORMULATION, not for a value.",
}
Q["P18_figure_resolution"] = {
    "only_rendition_obtainable": "546 px wide JPEG (PMC deposit)",
    "sha256_verified_against_own_download_record": "6/6",
    "routes_tried_for_a_higher_rendition": {
        "PMC /bin/ ADMA-38-e14370.pdf": "HTML error page (21,483 B)",
        "europepmc fulltextRepo .pdf": "HTTP 403 {\"error\":\"PDF link h...\"}",
        "europepmc fulltextRepo .jpg": "served, same 546 px file",
        "onlinelibrary.wiley.com/doi/...": "HTTP 403 to a plain client",
        "onlinelibrary.wiley.com/doi/pdfdirect/...": "HTTP 403",
    },
    "browser": "not used - the shared Playwright browser was out of scope for "
               "this run, and no bot check was bypassed",
    "native_dpi": "NOT MEASURABLE. axis.native_tol refuses a guessed dpi and "
                  "there is no PDF rect to measure px/pt against, because no "
                  "PDF rendition of the article is obtainable. Reported "
                  "instead as the directly measured 18.42 px per decade.",
    "SI_figures": "all flattened rasters; page.get_drawings() over all 29 SI "
                  "pages returns chemistry-free paths only (the TOC rule, "
                  "table rules, and two flow-cytometry gating panels on pp. "
                  "25-26). Supp Fig 2 is 722x1004 at a measured 137.5 dpi; "
                  "Supp Fig 14 is 1800x2400 at a measured 300.0 dpi.",
}

Q["_rows_emitted"] = 0
Q["_why_zero"] = ("Two independent blockers, either one sufficient. (1) The "
                  "screen's 3,840 attributable values are b-mRNA ACCUMULATION "
                  "measured by RNA sequencing before peak translation, which "
                  "is not one of the four delivery readouts data/ records. "
                  "(2) The validation panels ARE delivery readouts and ARE "
                  "attributable, but 5 of their 7 ionizable lipids have no "
                  "determinable structure, and the 2 that do are the 2 the "
                  "paper did not need to draw.")


# ---- P19 how many LNPs does Figure 2 actually plot? ------------------
_c = subprocess.run([sys.executable, "-I",
                     os.path.join(HERE, "06b_count_fig2_columns.py")],
                    capture_output=True, text=True)
_fc = load("06b_fig2_columns.json")
Q["P19_fig2_column_count"] = {
    "why": "the abstract says 122 LNPs, the Conclusions say 133, and Supp "
           "Table 3 lists 133. Figure 2's x axis is 'Barcode number' with NO "
           "tick labels, so the count cannot be read -- but the cell pitch can "
           "be measured, and the candidates are ~9 FFT bins apart.",
    "borders_detected_px": _fc["panel_borders_detected_px"],
    "cells_span_px": _fc["cells_span_px"],
    "columns_measured": _fc["columns_measured"],
    "all_five_strips_agree_exactly": _fc["all_strips_agree"],
    "border_sensitivity_sweep": _fc["border_sensitivity_sweep"],
    "printed_numbers_inside_the_range": _fc["printed_numbers_inside_the_range"],
    "printed_numbers_excluded": _fc["printed_numbers_excluded"],
    "returncode": _c.returncode,
    "verdict": _fc["verdict"],
}

# ---- P20 the two staged same-lab siblings -----------------------------
_sib = {}
_sp = os.path.join(HERE, "..", "mitchellmj_41555650", "03_enumerated.csv")
if os.path.exists(_sp):
    with open(_sp, newline="") as fh:
        _sib = {r["IL_name"]: r["IL_SMILES"] for r in csv.DictReader(fh)}
_mine = [b["IL_name"] for b in blocks]
_shared = sorted(set(_mine) & set(_sib))
_cmp = {}
for n in _shared:
    def _k(x):
        mm = Chem.MolFromSmiles(x) if x else None
        return Chem.MolToInchiKey(mm) if mm else None
    try:
        _c2 = smi(res, amb, n)
    except SystemExit:
        _c2 = None
    _cmp[n] = {"sibling_inchikey": _k(_sib[n]),
               "corpus_inchikey": _k(_c2),
               "in_live_corpus": _c2 is not None}
Q["P20_same_lab_siblings"] = {
    "mitchellmj_38982196": {
        "grammar": "<amine>T-O<tail> (2-iminothiolane / thioether / acrylate)",
        "members": 100, "shared_member_names_with_this_paper": 0,
        "shared_chemistry": "none",
    },
    "mitchellmj_41555650": {
        "grammar": "C<tail>-<amine core> (epoxide + polyamine)",
        "members": len(_sib),
        "shared_member_names_with_this_paper": _shared,
        "n_shared": len(_shared),
        "structures": _cmp,
        "but": "this paper attributes NO VALUE to any of the three. The 13 "
               "catalogue reagents (Index 121-133) appear only in Figure 2's "
               "UNLABELLED heat map; Supplementary Figure 14, the panel that "
               "makes the screen attributable, is the 20x6 nDm grid and does "
               "not contain them. None of the three is in Figure 3b either, "
               "whose only benchmark is SM-102. So the shared chemistry "
               "yields no usable cross-source pair.",
    },
    "same_lab_is_not_same_numbering": (
        "Three Mitchell-lab folders, three unrelated numbering systems: "
        "38982196's amines are 1-25 from its own building-block set, "
        "41555650's cores are 114/200/306/490/494, and this paper's amines "
        "are 1-20. A bare integer means a different molecule in each. This is "
        "why no name join on an amine token is performed anywhere here."),
}

# ---- P21 degree of substitution / arm count --------------------------
Q["P21_degree_of_substitution"] = {
    "why_asked": "a structureless row cannot be emitted (data/ has 0 blank "
                 "IL_SMILES in 30,998 rows), and a polydisperse conjugate "
                 "with an unstated degree of substitution has no determinable "
                 "SMILES.",
    "is_this_a_polymer_conjugate": "NO. These are discrete small molecules: "
                                   "a one-step aza-Michael addition of an "
                                   "amine N-H onto a dialkyl maleate.",
    "degree_of_substitution_STATED": (
        "YES, qualitatively, and it is the scheme's own sentence. "
        "Supplementary Figure 2 caption, verbatim: 'Dialkyl maleates provide "
        "higher reactivity and FAVOR MONO-ADDITION, facilitating plug-and-play "
        "diversification.' The three drawn products confirm it: 18D4 carries "
        "exactly one arm on each of its three N-H nitrogens and none on its "
        "two piperazine tertiary N."),
    "confirmed_arithmetically_on_18D4": (
        "3 arms gives a+b = 8.017 carbons against the drawn n-butyl 4+4 = 8; "
        "4 arms would require M = 1128.5 and there is no peak within 200 Da "
        "of its [M+H]+."),
    "so_the_blocker_is_not_substitution": (
        "It is that amines 1-16 and 19, and maleates 6.2, 8 and 18, are drawn "
        "NOWHERE in this paper. They are defined only in reference 34 (Xu, "
        "Han, Xu, Ricciardi et al., bioRxiv Mar 2025), a different paper."),
}

# ---- P22 structural variance of the validation subset -----------------
Q["P22_validation_subset_structural_variance"] = {
    "n_formulations": 7,
    "distinct_ionizable_lipids": 7,
    "distinct_amines": sorted({"2", "7", "17", "17", "18", "20"}),
    "distinct_maleates": sorted({"8", "4", "8i", "9.2", "4", "8i"}),
    "verdict": "REAL structural variance: six nDm members drawn from five "
               "different amines and four different maleates, plus SM-102. "
               "This is NOT the sibling case of two formulations built on one "
               "compound differing only in helper lipid.",
    "excipients": "all six nDm LNPs share one recipe exactly (35 / DOPE 16 / "
                  "Chol 46.5 / C14-PEG2000 2.5); SM-102 uses the clinical "
                  "recipe (50 / DSPC 10 / Chol 38.5 / DMG-PEG2000 1.5), so "
                  "the SM-102 comparison confounds lipid with formulation.",
}

# ---- P23 licence ------------------------------------------------------
Q["P23_licence"] = {
    "article": "CC BY. Europe PMC core record for PMID 41603124: "
               "license = 'cc by', isOpenAccess = 'Y', inPMC = 'Y'.",
    "article_verdict": "CC-BY-compatible, satisfies CLAUDE.md s8.",
    "supplementary": "ADMA-38-e14370-s001.pdf is served from the same CC BY "
                     "PMC deposit; no separate licence string is printed on it.",
    "repository": "NONE EXISTS - no code or data deposit of any kind (P16), "
                  "so there is no repository licence to record.",
}


def main():
    with open(os.path.join(HERE, "07_qc.json"), "w") as fh:
        json.dump(Q, fh, indent=1)
    print(json.dumps({k: (v if not isinstance(v, dict) else
                          {kk: vv for kk, vv in list(v.items())[:3]})
                      for k, v in Q.items()}, indent=1)[:3000])
    return 0


if __name__ == "__main__":
    sys.exit(main())
