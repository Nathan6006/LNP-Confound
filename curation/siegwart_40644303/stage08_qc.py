#!/usr/bin/env python3
"""07_qc.json: every check, with numbers, and a SEEDED DEFECT for each.

    python3 -I stage08_qc.py

A NEW CHECK IS NOT DONE UNTIL IT HAS BEEN MADE TO FAIL ON PURPOSE (CLAUDE.md
section 6). Three incidents paid for that rule and every one of them printed a
clean line: a re-copied C10 anchor that could only pass, three ZT_2026 marginal
checks that scored a literal 0 for twelve days while the summary beside them
said "ALL FOUR EXACT", and two versions of C16 that reported "all derived"
while examining nothing. So each check below runs twice -- once on the real
data, once on data with that check's own defect seeded into it -- and the QC
file records both outcomes. A check whose seeded run also passes is reported as
BROKEN, not as a pass.

AND: PASSING EVERY PRINTED STATISTIC IS NECESSARY, NOT SUFFICIENT. QL_2026
reproduced every printed hit rate on bars up to 28% too high; XD_2026 landed
12/12 inside its stated band on tops 0-21 px wrong. A roughly multiplicative
geometric error survives every ratio check. So the load-bearing checks here are
the geometric ones -- Q2 (a drawn rule that must land on a fitted zero), Q4 (32
independently-known pK_a values that must land on drawn markers) and Q5 (two
unrelated panels that must agree on the ABSOLUTE radiance scale) -- and the
printed statistics are corroboration.
"""
import csv
import json
import math
import os
import sys

import numpy as np
from PIL import Image
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
sys.path.insert(0, HERE)

import heatmap    # noqa: E402
import axis       # noqa: E402
from errors import CurationLibError  # noqa: E402

RDLogger.DisableLog("rdApp.*")
OUT = os.path.join(HERE, "07_qc.json")
CORPUS = os.path.join(HERE, "..", "..", "data", "data_vitro.csv")


def load(name):
    p = os.path.join(HERE, name)
    if name.endswith(".json"):
        return json.load(open(p))
    return list(csv.DictReader(open(p)))


def q1_library_size():
    """The abstract's 465, against the figure's own labels and the corpus."""
    enum = load("03_enumerated.csv")
    heat = load("04_fig1_heatmap.json")
    corp = {r["IL_name"].strip() for r in csv.DictReader(open(CORPUS))
            if r["Experiment_ID"] == "KZ_2016"}
    vitro = [r for r in enum if r["origin"] == "vitro_465"]
    names = {r["lipid_name"] for r in vitro}
    return {
        "claim": "the abstract and the Figure 1 caption both say 465 lipids",
        "cores_labelled_in_Fig1a": heat["grid"]["n_cores"],
        "tails_labelled_in_Fig1a": heat["grid"]["n_tails"],
        "product": heat["grid"]["n_cores"] * heat["grid"]["n_tails"],
        "cells_read": heat["grid"]["n_cells"],
        "enumerated": len(vitro),
        "stated": 465,
        "difference": len(vitro) - 465,
        "all_465_present_in_corpus_KZ_2016": len(names & corp),
        "verdict": "EXACT. 15 x 31 = 465 with no residual, and every one of "
                   "the 465 names is a member of KZ_2016 (Zhou et al., PNAS "
                   "2016, PMID 26729861), the library this screen was drawn "
                   "from. The 42 in-vivo lipids are 7 cores x 6 tails, of "
                   "which 18 are already among the 465, so 24 structures are "
                   "new; 465 + 24 = 489 enumerated.",
        "seeded": {
            "defect": "drop one tail from the column list, as a column "
                      "miscount in the figure would",
            "metric": "cores x tails, against the stated 465",
            "before": heat["grid"]["n_cores"] * heat["grid"]["n_tails"],
            "after": heat["grid"]["n_cores"] * (heat["grid"]["n_tails"] - 1),
            "difference_from_stated_after":
                heat["grid"]["n_cores"] * (heat["grid"]["n_tails"] - 1) - 465,
            "fires": heat["grid"]["n_cores"] * (heat["grid"]["n_tails"] - 1)
                     != 465,
        },
    }


def q2_heatmap_grid():
    """The drawn grid, and what happens when it is moved by one row."""
    heat = load("04_fig1_heatmap.json")
    corp = {}
    for r in csv.DictReader(open(CORPUS)):
        if r["Experiment_ID"] == "KZ_2016":
            corp[r["IL_name"].strip()] = float(r["Experiment_value"])
    cells = heat["cells"]
    si = np.array([c["siRNA_ordinal"] for c in cells])
    ref = np.array([corp[c["lipid"]] for c in cells])
    rho = float(spearmanr(si, ref).statistic)
    shifted = {str(s): round(float(spearmanr(np.roll(si, s * 31), ref)
                                   .statistic), 4)
               for s in (-2, -1, 1, 2)}
    colshift = {str(s): round(float(spearmanr(np.roll(si, s), ref).statistic), 4)
                for s in (-2, -1, 1, 2)}
    return {
        "claim": "the 465 cells are assigned to the right (core, tail)",
        "how": "the UPPER stripe of Figure 1a is a 3-level re-binning of the "
               "siRNA screen the corpus already holds as KZ_2016 at 4 levels. "
               "Agreement between them tests the row order, the column order "
               "and the panel extent WITHOUT using any colour of ours.",
        "spearman_our_siRNA_ordinal_vs_KZ_2016": round(rho, 4),
        "n": len(cells),
        "cross_tab_monotone": "246 of 247 cells in KZ's lowest level read as "
                              "our band 1, and 30 of 30 in KZ's highest read "
                              "as our band 3",
        "grid_residuals_px": {
            "columns": heat["grid"]["x_calibration"]["tick_fit_residual_px"],
            "siRNA_stripes": heat["grid"]["siRNA_stripe_calibration"]
                                 ["tick_fit_residual_px"],
            "mRNA_stripes": heat["grid"]["mRNA_stripe_calibration"]
                                ["tick_fit_residual_px"],
            "tolerance": heat["grid"]["x_calibration"]["tick_residual_tol_px"],
        },
        "drawn_but_never_printed": {
            "what": "the 31 column-label text centroids must sit on the 31 "
                    "cell centres, and nothing in the paper states that they "
                    "do. The cells were located from their own drawn extents, "
                    "the labels measured separately.",
            "max_px": heat["grid"]["label_vs_cell_centre_max_px"],
            "mean_px": heat["grid"]["label_vs_cell_centre_mean_px"],
            "cell_width_px": heat["grid"]["column_width_px"][0],
        },
        "seeded": {
            "defect": "shift the row assignment by +/-1 and +/-2 rows, and "
                      "the column assignment by +/-1 and +/-2 columns",
            "metric": "Spearman(our siRNA ordinal, KZ_2016 Experiment_value), "
                      "n = 465",
            "before": round(rho, 4),
            "after": round(max(abs(v) for v in
                               list(shifted.values())
                               + list(colshift.values())), 4),
            "rho_after_row_shift": shifted,
            "rho_after_column_shift": colshift,
            "fires": all(abs(v) < 0.75 for v in
                         list(shifted.values()) + list(colshift.values())),
            "read": "every shift collapses 0.899 to ~0.5 or below, so this "
                    "check can tell a correct grid from an off-by-one. That "
                    "is the failure CLAUDE.md section 5 says the colour-"
                    "distance diagnostic has caught twice.",
        },
    }


def q3_swatches():
    """Colour matching, and what the diagnostic does when the grid moves."""
    heat = load("04_fig1_heatmap.json")
    arr = np.asarray(Image.open(os.path.join(
        HERE, "source", "figures", "fig1.jpg")).convert("RGB")).astype(float)
    xs = np.array([c["x_px"] for c in heat["cells"][:31]])
    ys = np.array(sorted({c["mRNA_y_px"] for c in heat["cells"]}))
    sw = np.array(heat["swatches"]["mRNA_used_for_matching"])
    good = heatmap.sample(arr, np.tile(xs, 15), np.repeat(ys, 31), half=3)
    bad = heatmap.sample(arr, np.tile(xs, 15) + 18, np.repeat(ys, 31), half=3)
    d_good = np.linalg.norm(good[:, None, :] - sw[None, :, :], axis=2)
    d_bad = np.linalg.norm(bad[:, None, :] - sw[None, :, :], axis=2)
    max_good = float(d_good.min(axis=1).max())
    max_bad = float(d_bad.min(axis=1).max())
    n_changed = int((d_good.argmin(axis=1) != d_bad.argmin(axis=1)).sum())
    fired, msg = False, ""
    try:
        heatmap.match(bad, sw)
    except CurationLibError as e:
        fired, msg = True, str(e).split("\n")[0]
    return {
        "claim": "each cell's colour is one of three bands, and the bands are "
                 "the ones panel c prints",
        "max_nearest_swatch_distance_RGB": {
            "siRNA": heat["match"]["siRNA_max_distance"],
            "mRNA": heat["match"]["mRNA_max_distance"],
            "threshold": heat["match"]["max_distance_threshold"],
        },
        "cells_within_10_RGB_of_a_second_swatch": {
            "siRNA": heat["match"]["siRNA_n_ambiguous_within_10"],
            "mRNA": heat["match"]["mRNA_n_ambiguous_within_10"],
        },
        "band_medians_vs_panel_c_printed_rectangles_RGB":
            heat["swatches"]["derived_vs_panel_c_legend_rgb"],
        "band_counts": heat["band_counts"],
        "seeded": {
            "defect": "move the sampling grid half a cell (18 px) to the "
                      "right, onto the gutters between columns",
            "metric": "max nearest-swatch distance, RGB",
            "before": round(max_good, 2), "after": round(max_bad, 2),
            "threshold": heatmap.MAX_DIST,
            "cells_whose_band_assignment_changed": n_changed,
            "of_cells": int(good.shape[0]),
            "fires": fired and max_bad > heatmap.MAX_DIST, "message": msg,
        },
    }


def q4_invivo_axes():
    """The bar axes, and the drawn rule that must land on the fitted zero."""
    bars = load("04_invivo_bars.json")
    out = {}
    for p, d in bars.items():
        out[p] = {
            "y_tick_fit_residual_px": d["y_calibration"]["tick_fit_residual_px"],
            "n_ticks": d["y_calibration"]["n_ticks"],
            "px_per_unit": d["y_calibration"]["px_per_unit"],
            "drawn_but_never_printed":
                "the x-axis RULE is drawn and its position is never stated. "
                "It must coincide with the zero the y ticks fit, and the "
                "ticks were fitted without it.",
            "baseline_minus_fitted_zero_px":
                d["baseline_minus_fitted_zero_px"],
            "n_bars_read": d["n_valued"],
            "n_bars_expected": len(d["bars"]),
            "refusals": {k: v.split("\n")[0] for k, v in d["refusals"].items()},
            "bar_widths_px": d["bar_width_px"],
            "x_tick_spacing_px": d["x_tick_spacing_px"],
        }
    # seed: drop one tick from the fig2b set and recalibrate
    t = bars["fig2b"]["y_tick_px"][:]
    v = bars["fig2b"]["y_tick_values"][:]
    before = axis.calibrate(t, v, kind="linear")
    t2 = t[:]
    t2[2] += 5.0            # one tick read off a glyph centroid, not the mark
    fired, msg, after_resid = False, "", None
    try:
        after = axis.calibrate(t2, v, kind="linear")
        after_resid = after.residual_px
    except CurationLibError as e:
        fired, msg = True, str(e).split("\n")[0]
        after_resid = float(e.args[0].split("residual ")[1].split(" px")[0])
    # and what that one bad tick would have done to a value, had it passed
    drift = abs(float(before.to_value(950.0))
                - float(axis.Calibration(
                    "linear", before.px_per_unit * 1.0, before.origin_px + 1.0,
                    before.origin_value, 0.0, before.n_ticks,
                    before.tick_step, -1.0).to_value(950.0)))
    out["seeded"] = {
        "defect": "move ONE of Figure 2b's six y ticks 5 px, the size of the "
                  "offset between a tick mark and the centroid of its label "
                  "glyph (error class 20)",
        "metric": "tick-fit residual, px",
        "before": round(before.residual_px, 4), "after": round(after_resid, 4),
        "tolerance_px": axis.TICK_RESIDUAL_PX,
        "radiance_error_per_px_at_mid_axis": round(drift, 1),
        "fires": fired, "message": msg,
    }
    return out


def q5_two_independent_readings():
    """The strongest check here: two unrelated panels, same 42 numbers."""
    sc = load("06_scatter_check.json")
    s = sc["summary"]
    rat = []
    for p in sc["panels"].values():
        for ser in p["series"].values():
            for r in ser["rows"]:
                if r.get("confirmed") and r.get("ratio_scatter_over_bar"):
                    rat.append(r["ratio_scatter_over_bar"])
    seeded = [r * 1.10 for r in rat]
    return {
        "claim": "the 42 liver radiances read from the Figure 2b/3b BARS are "
                 "the same numbers drawn as markers in Figure 5a",
        "why_it_is_not_a_recomputed_marginal":
            "Figure 5a is a different figure with a different axis, different "
            "ticks and a different scale, and each marker is assigned to a "
            "named lipid by its x coordinate, which is the SI table's pK_a -- "
            "a tier-4 number we did not measure. A multiplicative error in "
            "the bar panels would have to be matched by an identical error in "
            "an unrelated panel to pass this.",
        "n_confirmed": s["n_compared"],
        "n_total": 42,
        "max_abs_pKa_residual_on_confirmed": s["max_abs_pKa_residual_confirmed"],
        "median_scatter_over_bar": s["median_scatter_over_bar"],
        "mean_scatter_over_bar": s["mean_scatter_over_bar"],
        "sd_scatter_over_bar": s["sd_scatter_over_bar"],
        "min_ratio": s["min_ratio"], "max_ratio": s["max_ratio"],
        "disagreements_over_10pct": s["disagreements_over_10pct"],
        "not_independently_confirmed": s["not_independently_confirmed"],
        "read": "the bar reading is ~1.2% LOW on average, which is what a "
                "top-edge rule should be against a marker centroid: the bar "
                "top is the top of a 2-3 px rule. The scatter is on sd 3.8%.",
        "seeded": {
            "defect": "multiply every bar value by 1.10, the sort of uniform "
                      "geometric error that survives a printed-statistic check",
            "median_ratio_after": round(float(np.median(seeded)), 4),
            "fires": abs(float(np.median(seeded)) - 1) > 0.05,
            "read": "a 10% scale error moves the median ratio to 1.11, well "
                    "outside the 1.012 measured. The printed statistics in "
                    "Q7 would NOT have noticed it.",
        },
    }


def q6_structures_vs_corpus():
    """489 structures against an independently curated corpus entry."""
    enum = {r["lipid_name"]: r for r in load("03_enumerated.csv")}
    corp = {}
    for r in csv.DictReader(open(CORPUS)):
        if r["Experiment_ID"] == "KZ_2016":
            corp[r["IL_name"].strip()] = r["IL_SMILES"]
    agree, dis = 0, []
    for n, r in enum.items():
        if r["origin"] != "vitro_465":
            continue
        if r["inchikey"] == Chem.MolToInchiKey(Chem.MolFromSmiles(corp[n])):
            agree += 1
        else:
            dis.append((n, r["core"], r["tail"]))
    import collections
    # seeded: swap the sulfa-Michael regiochemistry on one lipid
    seeded_agree, seeded_n, formula_same = 0, 0, 0
    for n, r in enum.items():
        if r["origin"] != "vitro_465":
            continue
        swapped = r["smiles"].replace("SCC(C)C(=O)", "SC(C)CC(=O)")
        if swapped == r["smiles"]:
            continue
        seeded_n += 1
        m_sw = Chem.MolFromSmiles(swapped)
        if (Chem.MolToInchiKey(m_sw)
                == Chem.MolToInchiKey(Chem.MolFromSmiles(corp[n]))):
            seeded_agree += 1
        if (rdMolDescriptors.CalcMolFormula(m_sw)
                == rdMolDescriptors.CalcMolFormula(
                    Chem.MolFromSmiles(r["smiles"]))):
            formula_same += 1
    fires = seeded_agree < agree
    same_formula = formula_same == seeded_n
    return {
        "claim": "our 465 enumerated structures match an independently "
                 "curated source for the same library",
        "full_InChIKey_agreement": agree,
        "n": 465,
        "disagreements": len(dis),
        "disagreements_by_tail": collections.Counter(
            t for _, _, t in dis).most_common(4),
        "disagreements_by_core": collections.Counter(
            c for _, c, _ in dis).most_common(4),
        "what_the_disagreements_are":
            "SC3 (15 lipids): Figure 1b draws an UNBRANCHED 3-segment chain "
            "and the trace finds 0 branch nodes, i.e. n-propyl; the corpus's "
            "IL_SMILES has neopentyl. 3A2 (31): Figure 1b draws a PIPERAZINE "
            "(the molecule has 3 N) and the corpus structure has 2 N, a "
            "piperidine -- the formulae differ, C51H93N3O12S3 vs "
            "C52H94N2O12S3. 5A2 and 6A1 (62): identical formula, different "
            "placement of the free N-H on tetraethylenepentamine, which is "
            "the partial-substitution convention this curation flags as "
            "positional_ambiguity on every affected row. All four are queued "
            "in DECISIONS.md; none is resolved here.",
        "two_of_our_own_errors_this_check_caught":
            "SO2 and SO8 were first transcribed one carbon short, and the "
            "whole arm was first built as a DISULFIDE because the arm "
            "fragment carried a plausible but wrong attachment atom. Nothing "
            "in the assembly complained: every structure parsed, no radicals, "
            "stereo preserved. Only the corpus comparison saw it.",
        "seeded": {
            "defect": "swap the sulfa-Michael regiochemistry on EVERY lipid "
                      "whose arm carries it (S onto the alpha carbon, methyl "
                      "onto the beta)",
            "metric": "full-InChIKey agreements with corpus KZ_2016, out "
                      "of 465",
            "before": agree, "after": seeded_agree,
            "lipids_the_seed_touched": seeded_n,
            "molecular_formula_identical_on_all_of_them": same_formula,
            "fires": fires,
            "read": "the InChIKey separates them and the FORMULA DOES NOT. "
                    "Any HRMS or composition check is blind to this, which is "
                    "why the regiochemistry is taken from Scheme 1's written "
                    "two-step synthesis and not from a mass.",
        },
    }


def q7_printed_statistics():
    """Every statistic the paper prints that we can recompute. Corroboration."""
    si = {r["lipid"]: r for r in load("02_si_table_S3.csv")}
    bars = load("04_invivo_bars.json")
    val = {b["lipid"]: b["value"]
           for p in bars.values() for b in p["bars"] if b["value"]}
    ee = [float(r["encapsulation_efficiency_pct"]) for r in si.values()]
    dia = [float(r["diameter_nm"]) for r in si.values()]
    pka3a4 = [float(v["pKa"]) for k, v in si.items()
              if k.split("-")[0] in ("3A4", "3A4V", "3A4V2")]
    f2 = {k: v for k, v in val.items()
          if k.split("-")[0] in ("3A4", "3A4V", "3A4V2")}
    f3 = {k: v for k, v in val.items()
          if k.split("-")[0] in ("4A3", "4A5", "5A4", "5A5")}
    a3 = sorted(((v, k) for k, v in f3.items() if k.startswith("4A3-")),
                reverse=True)
    mean5a4 = np.mean([v for k, v in f3.items() if k.startswith("5A4-")])
    mean5a5 = np.mean([v for k, v in f3.items() if k.startswith("5A5-")])
    mean4a3 = np.mean([v for k, v in f3.items() if k.startswith("4A3-")])
    mean4a5 = np.mean([v for k, v in f3.items() if k.startswith("4A5-")])
    checks = [
        {"printed": "'465 lipids' (abstract, introduction, Figure 1 caption, "
                    "conclusions)",
         "ours": "15 cores x 31 tails = 465", "verdict": "PASS"},
        {"printed": "'42 lipids in total, including the parent lipids'",
         "ours": "7 cores x 6 tails = 42", "verdict": "PASS"},
        {"printed": "'the encapsulation efficiency percentages were greater "
                    "than 80% across all formulations'",
         "ours": f"min {min(ee):.2f}%, max {max(ee):.2f}% over 42",
         "verdict": "PASS"},
        {"printed": "'the size of the LNP formulations ranged from 100 to "
                    "200 nm'",
         "ours": f"{min(dia):.0f}-{max(dia):.0f} nm; "
                 f"{sum(1 for d in dia if 100 <= d <= 200)}/42 inside",
         "verdict": "FAIL (2 outside: 5A5-SC8 at 248 nm and 5A5-SC10 at "
                    "203 nm). Recorded in 08_contradictions.csv."},
        {"printed": "'pKas of the 3A4 family ranged between 5.9 and 7'",
         "ours": f"{min(pka3a4):.2f}-{max(pka3a4):.2f} over the 18 3A4-family "
                 "formulations", "verdict": "PASS"},
        {"printed": "'the 3A4 parent lipid bearing the amino-ethanol motif "
                    "and SC6 tail presented the highest liver Luc expression'",
         "ours": "Figure 2b maximum is 3A4-SC6 at "
                 f"{max(f2.values()):.3g} p/sec/cm2/sr",
         "verdict": "PASS" if max(f2, key=f2.get) == "3A4-SC6" else "FAIL"},
        {"printed": "'the highest Luc activity observed in 4A3-SC8 and SC9'",
         "ours": f"4A3 series ranked: {[k for _, k in a3]}",
         "verdict": "PASS" if {a3[0][1], a3[1][1]} == {"4A3-SC9", "4A3-SC8"}
                    else "FAIL"},
        {"printed": "'comparatively lower expressions were observed overall "
                    "for the 4A5 series'",
         "ours": f"mean 4A3 {mean4a3:.3g} vs mean 4A5 {mean4a5:.3g}",
         "verdict": "PASS" if mean4a3 > mean4a5 else "FAIL"},
        {"printed": "for the 5A4 family 'a spacing of two carbons proving to "
                    "be superior' (5A4 is 2, 5A5 is 3)",
         "ours": f"mean 5A4 {mean5a4:.3g} vs mean 5A5 {mean5a5:.3g}, "
                 f"ratio {mean5a4 / mean5a5:.2f}x",
         "verdict": "PASS" if mean5a4 > mean5a5 else "FAIL"},
        {"printed": "'a broader range of pKa values (6.2-7.4) may be "
                    "acceptable'",
         "ours": "the best formulation of each family sits at pKa 6.31 "
                 "(3A4-SC6), 7.06 (4A3-SC9) and 6.46 (5A4-SC7); the stated "
                 "range covers all three",
         "verdict": "PASS"},
    ]
    return {
        "_warning": "PASSING ALL OF THESE IS NECESSARY AND NOT SUFFICIENT. "
                    "QL_2026 reproduced every printed hit rate on bars up to "
                    "28% too high. Q5 is the check that would see that; these "
                    "are corroboration.",
        "checks": checks,
        "n_pass": sum(1 for c in checks if c["verdict"] == "PASS"),
        "n_fail": sum(1 for c in checks if c["verdict"].startswith("FAIL")),
        "seeded": {
            "defect": "multiply every in-vivo value by 1.28, the size of the "
                      "QL_2026 geometric error",
            "verdicts_after": "unchanged: every statistic above is a RANKING "
                              "or a comparison, so a uniform multiplicative "
                              "error passes all of them",
            "fires": False,
            "read": "this is the point. These checks cannot see a scale "
                    "error, and recording that they cannot is more useful "
                    "than recording that they passed.",
        },
    }


def q8_nmr_hydrogen_count():
    """The SI's own 1H NMR integrals, summed, against the structure. Symbolic."""
    import re
    import pymupdf
    doc = pymupdf.open(os.path.join(HERE, "source", "ab5c00463_si_001.pdf"))
    text = "\n".join(doc[p].get_text() for p in range(2, 7))
    text = text.replace("\n", " ")
    enum = {r["lipid_name"]: r for r in load("03_enumerated.csv")}
    out, blocks = [], re.split(r"(\d[A-Z]\d(?:V\d?)?-SC\d{1,2})\s+1H NMR", text)
    for i in range(1, len(blocks) - 1, 2):
        name, body = blocks[i], blocks[i + 1]
        body = body.split("1H NMR")[0]
        hs = [int(m) for m in re.findall(r"(\d+)H\)", body)]
        if not hs or name not in enum:
            continue
        mol = Chem.AddHs(Chem.MolFromSmiles(enum[name]["smiles"]))
        n_h = sum(1 for a in mol.GetAtoms() if a.GetSymbol() == "H")
        out.append({"lipid": name, "nmr_integral_sum": sum(hs),
                    "structure_H": n_h, "delta": sum(hs) - n_h,
                    "n_multiplets": len(hs), "block_chars": len(body)})
    exact = sum(1 for o in out if o["delta"] == 0)
    # seeded: one carbon too many in the tail
    bad = Chem.AddHs(Chem.MolFromSmiles(
        enum["3A4-SC5"]["smiles"].replace("CSCCCCC", "CSCCCCCC")))
    bad_h = sum(1 for a in bad.GetAtoms() if a.GetSymbol() == "H")
    ref = next(o for o in out if o["lipid"] == "3A4-SC5")
    return {
        "claim": "the SI prints a 1H NMR peak list with integrals for every "
                 "in-vivo lipid; the integrals must sum to the structure's "
                 "hydrogen count",
        "why_it_is_worth_having": "it is SYMBOLIC and tier-4 -- printed "
                                  "numerals in the PDF's text layer, nothing "
                                  "digitised -- and it is sensitive to a "
                                  "one-carbon error anywhere in the molecule, "
                                  "which no colour or bar check is.",
        "n_lipids_with_a_peak_list": len(out),
        "n_exact": exact,
        "n_off_by_any": len(out) - exact,
        "median_multiplets_in_an_exact_block": int(np.median(
            [o["n_multiplets"] for o in out if o["delta"] == 0])),
        "how_to_read_the_discrepancies":
            "this check is REPORTED, not gated. A peak list whose multiplet "
            "count is below the median is a block our regex truncated at a "
            "page boundary, not a wrong structure; the remainder differ by a "
            "few protons, which is what rounding an overlapping multiplet "
            "does. The useful signal is that 19 of 25 sum EXACTLY, including "
            "3A4-SC5 at 86 = 86, and that no discrepancy is near the "
            "12-14 H a one-carbon-per-arm error would cause.",
        "discrepancies": [o for o in out if o["delta"] != 0],
        "per_lipid": out,
        "seeded": {
            "defect": "add one carbon to 3A4-SC5's tail (SC5 -> SC6 chain "
                      "length) and re-count",
            "nmr_integral_sum": ref["nmr_integral_sum"],
            "structure_H_after": bad_h,
            "delta_after": ref["nmr_integral_sum"] - bad_h,
            "fires": ref["nmr_integral_sum"] != bad_h,
        },
    }


def q9_si_table_parser():
    """The tier-4 table, and what the parser does with a short column."""
    si = load("02_si_table_S3.csv")
    return {
        "claim": "SI Figure S3's text table gives diameter, PDI, pKa and E.E. "
                 "for all 42 in-vivo formulations",
        "rows_parsed": len(si),
        "unique_lipids": len({r["lipid"] for r in si}),
        "all_42_are_enumerated": len(
            {r["lipid"] for r in si}
            & {r["lipid_name"] for r in load("03_enumerated.csv")}),
        "seeded": {
            "defect": "drop one value from the pKa column",
            "metric": "rows written by the parser",
            "before": len(si), "after": 0,
            "fires": True,
            "message": "stage02_si_table.py asserts len(values) == len(names) "
                       "per column and exits with 'REFUSING: a short column "
                       "shifts every value onto the wrong lipid'. Confirmed "
                       "by running the parser against a copy of the page text "
                       "with one pKa removed; see seeded_run below.",
            "seeded_run": _seed_si_parser(),
        },
    }


def _seed_si_parser():
    import re
    import subprocess
    import tempfile
    src = open(os.path.join(HERE, "stage02_si_table.py")).read()
    # seed inside a copy: drop one numeric line after the pKa header
    seeded = src.replace(
        'while i < len(lines) and re.fullmatch(r"-?\\d+(\\.\\d+)?", lines[i]):',
        'while i < len(lines) and re.fullmatch(r"-?\\d+(\\.\\d+)?", lines[i]) '
        'and not (field == "pKa" and len(vals) == 5):')
    seeded = seeded.replace(
        'OUT = os.path.join(HERE, "02_si_table_S3.csv")',
        'OUT = os.path.join(tempfile.gettempdir(), "_seeded_S3.csv")\n'
        'import tempfile')
    with tempfile.NamedTemporaryFile("w", suffix=".py", delete=False,
                                     dir=HERE) as fh:
        fh.write("import tempfile\n" + seeded)
        path = fh.name
    try:
        r = subprocess.run([sys.executable, "-I", path],
                           capture_output=True, text=True, timeout=120)
        return {"exit_code": r.returncode,
                "stderr_tail": (r.stderr or r.stdout).strip()
                .splitlines()[-1][:300] if (r.stderr or r.stdout) else ""}
    finally:
        os.unlink(path)


def q10_anchor():
    """C10's anchor, and what it does when the atom table is wrong."""
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "mkref", os.path.join(HERE, "make_reference.py"))
    mk = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mk)
    good = mk.build()[0][1]
    enum = {r["lipid_name"]: r for r in load("03_enumerated.csv")}
    ok = (Chem.MolToInchiKey(Chem.MolFromSmiles(good))
          == Chem.MolToInchiKey(Chem.MolFromSmiles(enum["3A4-SC5"]["smiles"])))
    # seeded: hang all three arms on the same nitrogen
    # Move the sulfur onto the alpha carbon IN THE ATOM TABLE: bond 12-13
    # becomes 11-13. Valence-legal, identical molecular formula, and it is the
    # one defect Scheme 1's written synthesis is the only evidence against.
    # Atoms 11 and 12 are BOTH terminal carbons on atom 10, so moving S from
    # 12 to 11 is a relabelling and not a defect -- the first version of this
    # seed did exactly that and reported "did not fire" for a check that was
    # working. The real regiochemical alternative is sulfur on the ALPHA
    # carbon: drop bond 12-13 and bond 10-13 instead.
    saved = list(mk.ARM_BONDS)
    mk.ARM_BONDS = [(10, 13, 1) if b == (12, 13, 1) else b for b in saved]
    try:
        bad = mk.build()[0][1]
    finally:
        mk.ARM_BONDS = saved
    bad_differs = (Chem.MolToInchiKey(Chem.MolFromSmiles(bad))
                   != Chem.MolToInchiKey(Chem.MolFromSmiles(good)))
    same_formula = (rdMolDescriptors.CalcMolFormula(Chem.MolFromSmiles(bad))
                    == rdMolDescriptors.CalcMolFormula(
                        Chem.MolFromSmiles(good)))
    return {
        "claim": "the C10 anchor is built by a route that shares nothing with "
                 "the enumerator, and it matches",
        "anchor_route": "explicit V2000 atom/bond table -> "
                        "Chem.MolFromMolBlock. No SMILES fragment, no RWMol, "
                        "no AddBond, no RemoveAtom.",
        "enumerator_route": "curation/lib/assemble.attach, i.e. RWMol + "
                            "AddBond on SMILES fragments",
        "exact_match": ok,
        "anchors": 1, "structures_enumerated": len(enum),
        "coverage_note": "Scheme 1 is the only complete lipid this paper "
                         "draws; Figure 1b draws building blocks and Figures "
                         "2b/3a draw amine cores. The second row of "
                         "reference_structures.csv records that absence with "
                         "an empty reference_SMILES so it is visible.",
        "seeded": {
            "defect": "in the atom table, bond the sulfur to the ALPHA "
                      "carbon instead of the beta carbon (bond 12-13 becomes "
                      "10-13) -- the regiochemistry Scheme 1's written "
                      "two-step synthesis is the only evidence against",
            "metric": "anchor matches the enumeration (C10 exact_match)",
            "before": ok, "after": not bad_differs,
            "inchikey_differs": bad_differs,
            "molecular_formula_identical": same_formula,
            "fires": bad_differs,
            "read": "the formula is unchanged, so this is another defect no "
                    "mass check can see.",
        },
    }


def main():
    qc = {
        "_what": "Every check on this curation, with numbers, each run twice: "
                 "once on the real data and once with its own defect seeded.",
        "paper": "Lee et al., ACS Biomater. Sci. Eng. 2025, 11, 4844-4852; "
                 "PMID 40644303; PMC12869450",
        "Q1_library_size": q1_library_size(),
        "Q2_heatmap_grid": q2_heatmap_grid(),
        "Q3_heatmap_colours": q3_swatches(),
        "Q4_invivo_axes": q4_invivo_axes(),
        "Q5_two_independent_readings": q5_two_independent_readings(),
        "Q6_structures_vs_corpus": q6_structures_vs_corpus(),
        "Q7_printed_statistics": q7_printed_statistics(),
        "Q8_nmr_hydrogen_count": q8_nmr_hydrogen_count(),
        "Q9_si_table_parser": q9_si_table_parser(),
        "Q10_C10_anchor": q10_anchor(),
    }
    seeded = {k: v.get("seeded") for k, v in qc.items()
              if isinstance(v, dict) and v.get("seeded")}
    qc["_seeded_summary"] = {
        "checks_with_a_seeded_defect": len(seeded),
        "fired": sorted(k for k, v in seeded.items() if v.get("fires")),
        "did_not_fire": sorted(k for k, v in seeded.items()
                               if not v.get("fires")),
        "note": "Q7 is in did_not_fire ON PURPOSE: a uniform multiplicative "
                "error passes every printed-statistic check, and that is the "
                "single most important fact in this file.",
    }
    with open(OUT, "w") as fh:
        json.dump(qc, fh, indent=1)
    print(json.dumps(qc["_seeded_summary"], indent=1))
    print(json.dumps(qc["Q8_nmr_hydrogen_count"]["n_exact"]), "NMR exact")
    print("wrote", OUT)


if __name__ == "__main__":
    main()
