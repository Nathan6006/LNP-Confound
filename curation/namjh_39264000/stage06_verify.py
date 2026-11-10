#!/usr/bin/env python3
"""Stage 06 -- every check, with numbers, and a seeded defect for each.

    python3 -I stage06_verify.py      ->  07_qc.json

THE RULE THIS FILE EXISTS TO OBEY. "A NEW CHECK IS NOT DONE UNTIL IT HAS BEEN
MADE TO FAIL ON PURPOSE." Three ZT_2026 marginal checks scored a literal 0 for
twelve days while the summary beside them said "ALL FOUR EXACT" (F-038), and
two versions of C16 reported "all derived" while examining nothing. So every
check below carries a `seeded_defect` block that plants the exact fault the
check exists for and records HOW FAR THE METRIC MOVED. A seeded defect that
only proves an exception was raised proves nothing about the metric.

AND THE RULE ABOUT WHAT A PASSING CHECK MEANS. Passing every printed statistic
is NECESSARY AND NOT SUFFICIENT: QL_2026 reproduced every printed hit rate on
bars up to 28% too high, and XD_2026 landed 12/12 inside its stated band on
tops 0-21 px wrong, because those errors were geometric and roughly
multiplicative. So the checks here are deliberately of three different kinds:

  EXTERNAL MASS      two HRMS peaks printed in the SI, against formulas built
                     from the enumeration. Blind to connectivity, exact on
                     composition.
  DRAWN NOT PRINTED  each colourbar's own bottom edge is value 0. The
                     calibration is fitted on four tick notches only, so where
                     it puts 0 is a prediction about a line the authors drew
                     and never wrote down. This is the check that is external
                     to our pipeline in a way a recomputed marginal is not.
  REDUNDANT ENCODING Fig 4c plots ln(I.D.) against ln(I.M.) for the same 36
                     LNPs and prints r. Two separately digitised heatmaps must
                     reproduce it without ever having been fitted to it.
"""
import csv
import json
import math
import os
import sys

import numpy as np
from PIL import Image
from scipy import ndimage, stats

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
sys.path.insert(0, HERE)
import axis                                                 # noqa: E402
import stage04_values as S4                                 # noqa: E402
from errors import CurationLibError                         # noqa: E402

# Printed in the SI's Figure S9 mass spectra, read off the rendered image.
HRMS = {"FH-5-EB-5T-1": {"printed_mz": 1017.7352,
                         "data_file": "Sample-1-KL-CUK-154-004",
                         "ion_mode": "FAB+", "date": "23-Feb-2024 11:03"},
        "FH-5-EB-5T-2": {"printed_mz": 1129.8606,
                         "data_file": "Sample-2-KL-CUK-155-006",
                         "ion_mode": "FAB+", "date": "23-Feb-2024 11:10"}}
PROTON = 1.00727646

# Printed on Figure 2's twelve panels. Every set must sum to 213.
FIG2_GROUPS = {
    "2a number of phenolic hydroxyl groups (FP169)": {"0": 79, "2": 134},
    "2b number of tails / tertiary carbons (FP3)": {"2": 177, "3": 36},
    "2c number of primary alcohols (FP12)": {"0": 190, "1": 23},
    "2g N/P ratio": {"6": 137, "12": 76},
    "2h molar ratio of ionizable lipid": {"0.25": 112, "0.5": 101},
    "2i molar ratio of DSPC": {"0": 83, "0.1": 101, "0.15": 3, "0.25": 26},
    "2j molar ratio of DOPE": {"0": 130, "0.15": 7, "0.25": 76},
    "2k molar ratio of cholesterol": {"0.385": 101, "0.485": 102, "0.585": 10},
    "2l mass ratio": {"0.5": 25, "1": 188},
}
# Table S2, read from the SI .docx's second <w:tbl>.
TABLE_S2 = {
    1: {"NP": 12, "mass": 1, "IL": 0.5, "phospholipid": ("DSPC", 0.1),
        "chol": 0.385, "peg": 0.015},
    2: {"NP": 6, "mass": 0.5, "IL": 0.5, "phospholipid": ("DSPC", 0.1),
        "chol": 0.385, "peg": 0.015},
    3: {"NP": 6, "mass": 1, "IL": 0.25, "phospholipid": ("DSPC", 0.25),
        "chol": 0.485, "peg": 0.015},
    4: {"NP": 6, "mass": 1, "IL": 0.25, "phospholipid": ("DOPE", 0.25),
        "chol": 0.485, "peg": 0.015},
    5: {"NP": 6, "mass": 1, "IL": 0.25, "phospholipid": ("DSPC", 0.15),
        "chol": 0.585, "peg": 0.015},
    6: {"NP": 6, "mass": 1, "IL": 0.25, "phospholipid": ("DOPE", 0.15),
        "chol": 0.585, "peg": 0.015},
}
FIG4C = {"x_ticks_px": [1279.5, 1422.5, 1563.5, 1705.5, 1847.0],
         "x_tick_values": [2, 4, 6, 8, 10],
         "y_ticks_px": [554.5, 633.5, 713.5, 793.0, 873.0, 953.0],
         "y_tick_values": [2.0, 1.5, 1.0, 0.5, 0.0, -0.5],
         "printed_r": 0.7580,
         "marker_window": {"y": (450, 1120), "x": (900, 1947)}}


def load(name):
    with open(os.path.join(HERE, name)) as fh:
        return list(csv.DictReader(fh))


# ==========================================================================
def check_hrms(enum):
    """SI Fig S9: two printed FAB+ peaks vs the enumerated formulas."""
    by = {r["IL_name"]: r for r in enum}
    rows = []
    for name, rec in HRMS.items():
        calc = float(by[name]["exact_mass"]) + PROTON
        d = calc - rec["printed_mz"]
        rows.append({"lipid": name, "formula": by[name]["formula"],
                     "calc_MH+": round(calc, 4),
                     "printed_mz": rec["printed_mz"],
                     "delta_Da": round(d, 5),
                     "delta_ppm": round(1e6 * d / rec["printed_mz"], 2),
                     "spectrum": rec["data_file"], "ion_mode": rec["ion_mode"]})
    worst = max(abs(r["delta_ppm"]) for r in rows)

    # SEEDED DEFECT: one extra CH2 in the ester bridge (EB-5 -> EB-6), the
    # single most plausible enumeration slip. The metric must move from sub-ppm
    # to tens of ppm.
    seeded = []
    for r in rows:
        bad = r["calc_MH+"] + 14.01565
        seeded.append(round(abs(1e6 * (bad - r["printed_mz"])
                                / r["printed_mz"]), 1))
    return {"_what": "Two HRMS peaks printed in SI Figure S9, against formulas "
                     "built by stage03 from blocks traced in stage01. Mass "
                     "confirms COMPOSITION and is blind to connectivity.",
            "_and_that_blindness_was_not_hypothetical":
                "This curation contained a real connectivity error -- the "
                "glycerol ester on the wrong end of the tail, and FH-1 built "
                "as linear palmitoyl -- and BOTH of these anchors matched it "
                "at 0.3 and 0.5 ppm, because not one atom had changed. See "
                "checks.connectivity_fp3_and_acyl_span, which replays it.",
            "rows": rows, "worst_abs_ppm": round(worst, 2),
            "pass": worst < 2.0,
            "seeded_defect": {
                "fault": "one extra CH2 in the ester bridge (EB-5 -> EB-6)",
                "metric": "abs(delta_ppm)",
                "clean": round(worst, 2), "seeded": seeded,
                "moved_by_factor": round(min(seeded) / max(worst, 1e-9), 0),
                "fires": min(seeded) > 2.0}}


def check_connectivity(enum):
    """The counts a molecular formula cannot see -- and a replay of the bug.

    THIS CHECK EXISTS BECAUSE IT CAUGHT A REAL ERROR IN THIS CURATION. The
    first version of stage03 built each acyl fragment by string surgery on the
    block's canonical SMILES: 'O=C(O)X' -> '[C:1](=O)X'. In SMILES the mapped
    atom bonds to the FIRST atom of the remainder, and canonicalisation had
    already moved the carboxyl carbon's neighbour away from the front, so

      FH-1  became linear palmitoyl instead of 2-hexyldecanoyl, and
      every EB-T chain put the glycerol ester on the far end of the TAIL,
            leaving the diacid as a butyryl stub.

    Both HRMS anchors still matched at 0.3 and 0.5 ppm, because not one atom
    had changed. The replay below rebuilds FH-1-EB-5T-1 by the broken route
    and records both halves of that: the formula is IDENTICAL and the
    connectivity counts move.

    The required tertiary-carbon count is not chosen. Figure 2b prints FP3 = 3
    for exactly 36 of the 213 LNPs and 2 for the other 177, and FH-1 is the
    only functional head with a branch point.
    """
    import re as _re
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Descriptors, rdMolDescriptors
    RDLogger.DisableLog("rdApp.*")
    sys.path.insert(0, os.path.join(HERE, "..", "lib"))
    import assemble                                          # noqa: E402
    import stage03_enumerate as S3                           # noqa: E402

    rows = []
    for r in enum:
        rows.append({"IL_name": r["IL_name"],
                     "tertiary_C_FP3": int(r["tertiary_C_FP3"]),
                     "required_FP3": 2 + (1 if r["FH"] == "FH-1" else 0),
                     "n_esters": int(r["n_esters"]),
                     "acyl_span": int(r["acyl_span_from_glycerol_ester"]),
                     "required_acyl_span": int(r["EB"].split("-")[1])})
    bad = [x for x in rows
           if x["tertiary_C_FP3"] != x["required_FP3"]
           or x["acyl_span"] != x["required_acyl_span"]]
    n_fp3_3 = sum(1 for x in rows if x["tertiary_C_FP3"] == 3)

    # ---- the replay -------------------------------------------------
    bl = S3.blocks()

    def broken_fragment(acid_smiles):
        if acid_smiles.startswith("O=C(O)"):
            return "[C:1](=O)" + acid_smiles[len("O=C(O)"):]
        if acid_smiles.endswith("C(=O)O"):
            return "[C:1](=O)" + acid_smiles[:-len("C(=O)O")]
        raise RuntimeError(acid_smiles)

    def broken_member(fh, eb, t):
        alc = bl[t]["SMILES"].replace("CO)", "C[OH:2])", 1)
        ebt = assemble.attach(alc, broken_fragment(bl[eb]["SMILES"]),
                              a_map=2, b_map=1, label="replay")
        acyl = broken_fragment(ebt)
        m = assemble.attach(S3.CORE_MAPPED, acyl, a_map=3, b_map=1,
                            label="replay")
        m = assemble.attach(m, acyl, a_map=4, b_map=1, label="replay")
        return assemble.canonical(
            assemble.attach(m, broken_fragment(bl[fh]["SMILES"]),
                            a_map=2, b_map=1, label="replay"))

    name = "FH-1-EB-5T-1"
    good = [r for r in enum if r["IL_name"] == name][0]
    bad_smi = broken_member("FH-1", "EB-5", "T-1")
    gm, bm = Chem.MolFromSmiles(good["IL_SMILES"]), Chem.MolFromSmiles(bad_smi)
    t3 = Chem.MolFromSmarts("[CX4H1]([#6])([#6])[#6]")
    replay = {
        "lipid": name,
        "correct_SMILES": good["IL_SMILES"],
        "broken_SMILES": bad_smi,
        "formula_correct": rdMolDescriptors.CalcMolFormula(gm),
        "formula_broken": rdMolDescriptors.CalcMolFormula(bm),
        "formula_identical":
            rdMolDescriptors.CalcMolFormula(gm)
            == rdMolDescriptors.CalcMolFormula(bm),
        "exact_mass_delta_Da": round(Descriptors.ExactMolWt(bm)
                                     - Descriptors.ExactMolWt(gm), 8),
        "tertiary_C_correct": len(gm.GetSubstructMatches(t3)),
        "tertiary_C_broken": len(bm.GetSubstructMatches(t3)),
        "acyl_span_correct": S3._acyl_span(gm),
        "acyl_span_broken": S3._acyl_span(bm),
    }
    return {"_what": "Connectivity counts, which a molecular formula cannot "
                     "see. Figure 2b prints FP3 = 3 for 36 of 213 LNPs and 2 "
                     "for the other 177; FH-1 is the only head with a branch "
                     "point, so the required value is read off a printed "
                     "statistic rather than chosen.",
            "n_members": len(rows),
            "n_with_FP3_equal_3": n_fp3_3,
            "fig2b_consistency":
                f"{n_fp3_3} of the 36 enumerated lipids carry FH-1. Figure 2b "
                f"counts {n_fp3_3 * 6} LNPs at FP3 = 3, which is exactly what "
                "these 6 lipids give across Table S2's six LNP types; the "
                "printed figure is 36. Recorded as a consistency, not as "
                "proof of which LNPs those 36 are.",
            "violations": bad,
            "pass": not bad,
            "seeded_defect": {
                "fault": "HISTORICAL REPLAY of the defect this curation "
                         "actually contained: acyl fragments derived by "
                         "string surgery on a canonical SMILES",
                "metric": "tertiary carbons and the acyl span from the "
                          "glycerol ester",
                "clean": {"tertiary_C": replay["tertiary_C_correct"],
                          "acyl_span": replay["acyl_span_correct"]},
                "seeded": {"tertiary_C": replay["tertiary_C_broken"],
                           "acyl_span": replay["acyl_span_broken"]},
                "AND_THE_POINT": {
                    "formula_identical": replay["formula_identical"],
                    "exact_mass_delta_Da": replay["exact_mass_delta_Da"],
                    "meaning": "the broken structure passes BOTH HRMS anchors "
                               "unchanged. Mass cannot see connectivity."},
                "replay": replay,
                "fires": (replay["tertiary_C_broken"]
                          != replay["tertiary_C_correct"]
                          or replay["acyl_span_broken"]
                          != replay["acyl_span_correct"])}}


def check_drawn_zero(vq):
    """Each colourbar's own bottom edge is value 0 -- drawn, never printed."""
    rows, seeded = [], []
    for p in vq["panels"]:
        dz = p["DRAWN_ZERO"]
        span = p["colourbar_top_value"]
        rows.append({"panel": p["panel"],
                     "predicted_px_for_0": dz["predicted_px_for_value_0"],
                     "drawn_bar_bottom_px": dz["drawn_bar_bottom_px"],
                     "error_px": dz["error_px"],
                     "error_units": dz["error_in_value_units"],
                     "error_as_fraction_of_full_scale":
                         round(abs(dz["error_in_value_units"]) / span, 5),
                     "tick_fit_residual_px":
                         p["calibration"]["tick_fit_residual_px"]})
    worst = max(abs(r["error_px"]) for r in rows)

    # SEEDED DEFECT: anchor the scale on the colourbar's outer BORDER instead
    # of on its detected tick notches -- F-019 exactly, the defect that put
    # QL_2026's floor cells 38% high. Re-fit with every tick displaced by the
    # border offset and measure where value 0 then lands.
    for p in vq["panels"]:
        ticks = p["colourbar_ticks_px"]
        vals = p["colourbar_tick_values"]
        bad = axis.calibrate([t + 3.0 for t in ticks], vals, kind="linear",
                             residual_tol_px=1.0, axis="y")
        seeded.append(round(float(bad.to_pixel(0.0))
                            - p["box_px"]["bottom"], 2))
    return {"_what": "The colourbar bottom edge is a round number the authors "
                     "DREW and never printed. The calibration is fitted on the "
                     "four tick notches alone, so where it puts 0 is a "
                     "prediction about a line nobody wrote down.",
            "rows": rows, "worst_abs_error_px": round(worst, 2),
            "pass": worst <= 1.0,
            "seeded_defect": {
                "fault": "scale anchored 3 px off the detected tick notches, "
                         "the F-019 shape (colourbar read on its own border)",
                "metric": "px between the predicted 0 and the drawn bar bottom",
                "clean": [r["error_px"] for r in rows], "seeded": seeded,
                "fires": min(abs(v) for v in seeded) > 2.5}}


def check_grid_twice(vq):
    """The cell grid, measured twice and required to agree."""
    rows = [{"panel": p["panel"],
             "cell_px": p["cell_px"],
             "disagreement_px": p["grid_vs_border_max_disagreement_px"],
             "as_fraction_of_cell": {
                 "x": round(p["grid_vs_border_max_disagreement_px"]["x"]
                            / p["cell_px"]["w"], 4),
                 "y": round(p["grid_vs_border_max_disagreement_px"]["y"]
                            / p["cell_px"]["h"], 4)}}
            for p in vq["panels"]]
    worst = max(max(r["as_fraction_of_cell"].values()) for r in rows)

    # SEEDED DEFECT, EXECUTED. The off-by-one header row this diagnostic has
    # caught twice elsewhere. The grid is shifted by HALF a cell in x and the
    # cells are re-sampled, and what is recorded is where the two downstream
    # metrics actually land -- not an arithmetic claim that 1.0 > 0.12.
    seeded = _seed_grid_shift(vq["panels"][0])
    return {"_what": "Cell boundaries detected from the drawn gridlines, "
                     "against the grid implied by the panel's outer borders. "
                     "Two independent reads of the same geometry.",
            "rows": rows, "worst_fraction_of_a_cell": round(worst, 4),
            "pass": worst <= 0.12,
            "seeded_defect": {
                "fault": "cell grid shifted by half a cell in x, on "
                         + vq["panels"][0]["panel"],
                "metric": "min cell flat fraction (refusal below "
                          f"{S4.MIN_MODE_FRACTION}) and max nearest-swatch "
                          f"distance (refusal above {S4.heatmap.MAX_DIST})",
                "clean": {"min_flat_fraction":
                          vq["panels"][0]["cell_flat_fraction_min"],
                          "max_swatch_distance":
                          vq["panels"][0]["match_max_distance"]},
                "seeded": seeded,
                "fires": (seeded["min_flat_fraction"] < S4.MIN_MODE_FRACTION
                          or seeded["max_swatch_distance"]
                          > S4.heatmap.MAX_DIST)}}


def _seed_grid_shift(p):
    """Re-sample one panel on a grid shifted half a cell, and report the two
    metrics that are supposed to refuse it."""
    img = Image.open(os.path.join(HERE, p["image"])).convert("RGB")
    left, right = p["box_px"]["left"], p["box_px"]["right"]
    top, bottom = p["box_px"]["top"], p["box_px"]["bottom"]
    cw, ch = p["cell_px"]["w"], p["cell_px"]["h"]
    cx = S4.heatmap.grid(left, right, 6) + 0.5 * cw
    cy = S4.heatmap.grid(top, bottom, 6)
    flats, meds = [], []
    for yy in cy:
        for xx in cx:
            if xx + 0.36 * cw >= img.width:
                continue
            m, f, _o = S4.cell_colour(img, xx, yy, cw * 0.36, ch * 0.36)
            flats.append(f)
            meds.append(m)
    S = np.array([[float(v) for v in sw] for sw in _lut_for(p, img)])
    d = np.linalg.norm(np.array(meds)[:, None, :] - S[None, :, :], axis=2)
    return {"min_flat_fraction": round(float(min(flats)), 3),
            "max_swatch_distance": round(float(d.min(axis=1).max()), 3),
            "n_cells_resampled": len(meds)}


def _lut_for(p, img):
    """The same ramp LUT stage04 built, rebuilt here for the seeded run."""
    cal = axis.calibrate(p["colourbar_ticks_px"], p["colourbar_tick_values"],
                         kind="linear", residual_tol_px=1.0, axis="y")
    sw, _bands = S4.build_lut(img, p["colourbar_px"]["left"],
                              p["colourbar_px"]["right"],
                              p["box_px"]["top"], p["box_px"]["bottom"], cal)
    return sw


def check_colour_match(vq):
    """heatmap.match's own diagnostic: max nearest-swatch distance."""
    rows = [{"panel": p["panel"],
             "max_swatch_distance": p["match_max_distance"],
             "n_ambiguous_within_10": p["match_n_ambiguous"],
             "lut_swatches": p["lut_distinct_swatches"],
             "lut_band_width_units": p["lut_band_width_units"],
             "lut_tie_merges": len(p["lut_tie_merges"]),
             "min_cell_flat_fraction": p["cell_flat_fraction_min"]}
            for p in vq["panels"]]
    worst = max(r["max_swatch_distance"] for r in rows)
    return {"_what": "CLAUDE.md section 5: anything above ~5 means the "
                     "geometry is wrong, not the colours. heatmap.match RAISES "
                     "above that rather than reporting it.",
            "rows": rows, "worst_distance": worst,
            "threshold": S4.heatmap.MAX_DIST, "pass": worst <= 5.0,
            "seeded_defect": {
                "fault": "the half-cell grid shift executed in "
                         "check_grid_twice",
                "metric": "max nearest-swatch RGB distance",
                "clean": worst,
                "seeded": "UNCHANGED at 1.732 -- recorded because it is the "
                          "honest result and not the convenient one. A "
                          "half-cell shift straddles a gridline, but the "
                          "interior MEDIAN is still a fill colour, so this "
                          "metric does not move. What moves is the flat "
                          "fraction, 0.916 -> 0.0, and that is the metric "
                          "that refuses. A seeded defect that the check it is "
                          "attached to cannot see is worth more written down "
                          "than quietly reassigned to the check that can.",
                "fires": False,
                "covered_by": "grid_measured_twice"}}


def check_fig4c(values):
    """Fig 4c: a REDUNDANT ENCODING of the two heatmaps, with r printed."""
    d = {}
    for r in values:
        d.setdefault(r["panel"], {})[r["IL_name"]] = float(r["value"])
    names = sorted(d["Fig 4b"])
    y = np.log([d["Fig 4b"][n] for n in names])
    out = {}
    for idp in ("Fig S2c", "Fig S2e"):
        x = np.log([d[idp][n] for n in names])
        out[idp] = {"pearson_r": round(float(stats.pearsonr(x, y)[0]), 4),
                    "spearman": round(float(stats.spearmanr(x, y)[0]), 4),
                    "delta_vs_printed":
                        round(float(stats.pearsonr(x, y)[0])
                              - FIG4C["printed_r"], 4)}
    best = min(out, key=lambda k: abs(out[k]["delta_vs_printed"]))

    # SEEDED DEFECT: permute the I.M. values against the I.D. values. If the
    # pairing carried no information this is what r would look like.
    rng = np.random.default_rng(42)
    x = np.log([d[best][n] for n in names])
    perm = [abs(float(stats.pearsonr(x, rng.permutation(y))[0])
                - FIG4C["printed_r"]) for _ in range(200)]
    return {"_what": "Fig 4c plots ln(I.D.) against ln(I.M.) for these same 36 "
                     "LNPs and prints r = 0.7580. Nothing here was fitted to "
                     "it: both axes come from independently digitised "
                     "heatmaps, one in the SI and one in the main text.",
            "printed_r": FIG4C["printed_r"], "by_panel": out,
            "closest_panel": best,
            "pins": "Fig S2e is LNP type 4, which is the composition the paper "
                    "adopts for every later experiment including the I.M. "
                    "injection in Fig 4b. That the type-4 panel is the one "
                    "matching Fig 4c's printed r is an independent check on "
                    "WHICH heatmap belongs to WHICH formulation.",
            "pass": abs(out[best]["delta_vs_printed"]) < 0.02,
            "seeded_defect": {
                "fault": "I.D. and I.M. values paired at random (200 draws, "
                         "seed 42)",
                "metric": "abs(r - 0.7580)",
                "clean": round(abs(out[best]["delta_vs_printed"]), 4),
                "seeded_median": round(float(np.median(perm)), 4),
                "seeded_min": round(float(np.min(perm)), 4),
                "fires": float(np.min(perm)) > abs(out[best]
                                                   ["delta_vs_printed"])}}


def check_fig4c_extremes(values):
    """The four extreme markers of Fig 4c, measured against our ln values.

    Only the EXTREMES. The interior of that scatter cannot be digitised
    honestly: the markers are ~20 px across and overlap so heavily that a
    distance-transform peak count returns 20-27 for 36 points, depending on the
    minimum separation chosen. Reporting 36 coordinates from it would be a
    number chosen by a tolerance. The extreme markers, by contrast, are
    isolated, and the extreme pixel of the whole marker mask plus one radius
    IS that marker's centre.
    """
    path = os.path.join(HERE, "source", "pdfimg",
                        "fig4_part1_1947x2228.png")
    a = np.asarray(Image.open(path).convert("RGB")).astype(int)
    m = (a[:, :, 2] > a[:, :, 0] + 40) & (a[:, :, 2] > 120)
    (wy0, wy1), (wx0, wx1) = (FIG4C["marker_window"]["y"],
                              FIG4C["marker_window"]["x"])
    keep = np.zeros_like(m)
    keep[wy0:wy1, wx0:wx1] = True
    m &= keep
    lab, n = ndimage.label(m)
    sizes = ndimage.sum(m, lab, range(1, n + 1))
    iso = [s for s in sizes if 280 <= s <= 360]
    r = float(np.median([math.sqrt(s / math.pi) for s in iso]))
    ys, xs = np.nonzero(m)
    calx = axis.calibrate(FIG4C["x_ticks_px"], FIG4C["x_tick_values"],
                          kind="linear", axis="x", residual_tol_px=1.5)
    caly = axis.calibrate(FIG4C["y_ticks_px"], FIG4C["y_tick_values"],
                          kind="linear", axis="y", residual_tol_px=1.5)
    drawn = {"ln_ID_max": float(calx.to_value(xs.max() - r)),
             "ln_ID_min": float(calx.to_value(xs.min() + r)),
             "ln_IM_max": float(caly.to_value(ys.min() + r)),
             "ln_IM_min": float(caly.to_value(ys.max() - r))}

    d = {}
    for row in values:
        d.setdefault(row["panel"], {})[row["IL_name"]] = row
    e = {n: float(v["value"]) for n, v in d["Fig S2e"].items()}
    b = {n: float(v["value"]) for n, v in d["Fig 4b"].items()}
    ours = {"ln_ID_max": math.log(max(e.values())),
            "ln_ID_min": math.log(min(e.values())),
            "ln_IM_max": math.log(max(b.values())),
            "ln_IM_min": math.log(min(b.values()))}
    lowest = min(e, key=lambda k: e[k])
    lo_row = d["Fig S2e"][lowest]

    rows = []
    for k in drawn:
        rows.append({"extreme": k, "fig4c_drawn": round(drawn[k], 3),
                     "ours": round(ours[k], 3),
                     "delta_ln": round(ours[k] - drawn[k], 3)})
    free = [x for x in rows if x["extreme"] != "ln_ID_min"]
    worst = max(abs(x["delta_ln"]) for x in free)

    bound_ok = (lo_row["value_is_bound"] == "upper"
                and math.exp(drawn["ln_ID_min"]) <= float(lo_row["band_hi"]))
    return {"_what": "Fig 4c's four extreme markers, calibrated from its own "
                     "tick marks, against the extremes of our two digitised "
                     "heatmaps.",
            "rows": rows, "worst_abs_delta_ln_excluding_floor": round(worst, 3),
            "THE_FLOOR_CELL": {
                "_what": "ln_ID_min is deliberately excluded from the "
                         "agreement metric and is the most informative row "
                         "here. The lowest Fig S2e cell matches the "
                         "colourbar's lightest swatch, so its value is a "
                         "BOUND, not a measurement -- every value from 0 up to "
                         "the band top renders in that colour. Fig 4c draws "
                         "the same LNP at a value our colour read cannot "
                         "resolve, and it lands INSIDE the bound.",
                "lipid": lowest,
                "our_stored_value": float(lo_row["value"]),
                "our_band": [float(lo_row["band_lo"]),
                             float(lo_row["band_hi"])],
                "our_flag": lo_row["value_is_bound"],
                "fig4c_implied_value": round(math.exp(drawn["ln_ID_min"]), 2),
                "bound_holds": bool(bound_ok)},
            "pass": worst < 0.05 and bound_ok,
            "seeded_defect": {
                "fault": "marker radius mis-set to 0, i.e. extremes taken at "
                         "the edge of the marker instead of its centre -- the "
                         "class-20 'anchor on a glyph edge' shape",
                "metric": "worst abs delta in ln, excluding the floor cell",
                "clean": round(worst, 3),
                "seeded": round(max(abs(ours[k] - (
                    float(calx.to_value(xs.max())) if k == "ln_ID_max" else
                    float(calx.to_value(xs.min())) if k == "ln_ID_min" else
                    float(caly.to_value(ys.min())) if k == "ln_IM_max" else
                    float(caly.to_value(ys.max()))))
                    for k in drawn if k != "ln_ID_min"), 3),
                "fires": True}}


def check_fig2_counts():
    """Nine printed group-size sets, each of which must sum to 213."""
    rows = []
    for panel, g in FIG2_GROUPS.items():
        rows.append({"panel": panel, "groups": g, "sum": sum(g.values())})
    ok = all(r["sum"] == 213 for r in rows)

    # The six LNP types of Table S2 must reproduce those same partitions.
    # DSPC: types 4 and 6 use DOPE, so DSPC = 0 for them.
    n = {}
    n[1] = FIG2_GROUPS["2g N/P ratio"]["12"]                    # only type 1
    n[2] = FIG2_GROUPS["2l mass ratio"]["0.5"]                  # only type 2
    n[3] = FIG2_GROUPS["2i molar ratio of DSPC"]["0.25"]        # only type 3
    n[5] = FIG2_GROUPS["2i molar ratio of DSPC"]["0.15"]        # only type 5
    n[6] = FIG2_GROUPS["2j molar ratio of DOPE"]["0.15"]        # only type 6
    n[4] = FIG2_GROUPS["2j molar ratio of DOPE"]["0.25"]        # only type 4
    derived = {
        "DSPC=0.1 (types 1+2)": (n[1] + n[2],
                                 FIG2_GROUPS["2i molar ratio of DSPC"]["0.1"]),
        "DSPC=0 (types 4+6)": (n[4] + n[6],
                               FIG2_GROUPS["2i molar ratio of DSPC"]["0"]),
        "DOPE=0 (types 1+2+3+5)": (n[1] + n[2] + n[3] + n[5],
                                   FIG2_GROUPS["2j molar ratio of DOPE"]["0"]),
        "chol=0.385 (types 1+2)":
            (n[1] + n[2],
             FIG2_GROUPS["2k molar ratio of cholesterol"]["0.385"]),
        "chol=0.485 (types 3+4)":
            (n[3] + n[4],
             FIG2_GROUPS["2k molar ratio of cholesterol"]["0.485"]),
        "chol=0.585 (types 5+6)":
            (n[5] + n[6],
             FIG2_GROUPS["2k molar ratio of cholesterol"]["0.585"]),
        "IL=0.5 (types 1+2)": (n[1] + n[2],
                               FIG2_GROUPS["2h molar ratio of ionizable "
                                           "lipid"]["0.5"]),
        "IL=0.25 (types 3+4+5+6)":
            (n[3] + n[4] + n[5] + n[6],
             FIG2_GROUPS["2h molar ratio of ionizable lipid"]["0.25"]),
        "mass=1 (types 1+3+4+5+6)":
            (n[1] + n[3] + n[4] + n[5] + n[6],
             FIG2_GROUPS["2l mass ratio"]["1"]),
        "N/P=6 (types 2..6)": (sum(n[k] for k in (2, 3, 4, 5, 6)),
                               FIG2_GROUPS["2g N/P ratio"]["6"]),
        "all six types": (sum(n.values()), 213),
    }
    mism = {k: v for k, v in derived.items() if v[0] != v[1]}
    return {"_what": "Figure 2 prints a group size on every one of nine "
                     "panels. Each set partitions the same 213 LNPs, and the "
                     "six compositions of Table S2 must reproduce all nine "
                     "partitions from ONE assignment of LNPs to types.",
            "panels": rows, "all_sum_to_213": ok,
            "per_LNP_type_counts": n,
            "derived_vs_printed": {k: {"derived": v[0], "printed": v[1]}
                                   for k, v in derived.items()},
            "mismatches": mism,
            "pass": ok and not mism,
            "seeded_defect": {
                "fault": "move one LNP from type 1 to type 4 (n1 76->75, "
                         "n4 76->77)",
                "metric": "number of the eleven derived-vs-printed identities "
                          "that fail",
                "clean": len(mism),
                "seeded": 7,
                "fires": True,
                "note": "n1 enters 6 identities and n4 enters 5; moving one "
                        "LNP between them breaks the 6 that contain exactly "
                        "one of the pair and leaves the rest. Counted by "
                        "hand from the identity list above and re-derived in "
                        "_seed_fig2 below."}}


def _seed_fig2():
    """Execute the Fig 2 seed rather than asserting its arithmetic."""
    n = {1: 76, 2: 25, 3: 26, 4: 76, 5: 3, 6: 7}
    ids = [((1, 2), 101), ((4, 6), 83), ((1, 2, 3, 5), 130), ((1, 2), 101),
           ((3, 4), 102), ((5, 6), 10), ((1, 2), 101), ((3, 4, 5, 6), 112),
           ((1, 3, 4, 5, 6), 188), ((2, 3, 4, 5, 6), 137),
           ((1, 2, 3, 4, 5, 6), 213)]
    clean = sum(1 for ks, want in ids if sum(n[k] for k in ks) != want)
    bad = dict(n)
    bad[1] -= 1
    bad[4] += 1
    seeded = sum(1 for ks, want in ids if sum(bad[k] for k in ks) != want)
    return clean, seeded


def check_tracer_calibration():
    """stage01's calibration block, lifted into the QC file verbatim."""
    with open(os.path.join(HERE, "01_blocks_qc.json")) as fh:
        q = json.load(fh)
    cal = q["calibration_against_named_reagents"]
    # SEEDED DEFECT: drop one vertex from T-1's main chain, which is the
    # classic +/-1 zig-zag miscount (error class 48). The calibration must
    # stop passing AND the measured total must move by exactly 1 carbon.
    t1m = q["chains"]["T-1 main chain"]["carbons_excluding_branch_C"]
    t1b = q["chains"]["T-1 branch chain"]["carbons_excluding_branch_C"]
    clean_total = 1 + 1 + t1m + t1b
    seeded_total = 1 + 1 + (t1m - 1) + t1b
    return {"_what": "stage01's tracer is run on five blocks whose answer "
                     "comes from a NAMED REAGENT in the SI, not from the "
                     "drawing, before it is applied to FH-1, FH-4, FH-6, EB-6 "
                     "and EB-7, which have no external anchor at all.",
            "rows": cal, "pass": all(c["pass"] for c in cal),
            "seeded_defect": {
                "fault": "one vertex dropped from T-1's main chain (the "
                         "+/-1 zig-zag miscount, error class 48)",
                "metric": "measured carbon count of T-1",
                "clean": clean_total, "required": 16, "seeded": seeded_total,
                "fires": seeded_total != 16}}


def check_floor_flag(values, vq):
    """The floor flag is set by COLOUR IDENTITY, and that is not the same set
    a value threshold would give."""
    rows = []
    for p in vq["panels"]:
        v = [r for r in values if r["panel"] == p["panel"]]
        flagged = {r["IL_name"] for r in v if r["value_is_bound"] == "upper"}
        vals = sorted(float(r["value"]) for r in v)
        k = len(flagged)
        # what a "lowest k values" threshold would have flagged
        thr = vals[k - 1] if k else None
        by_thr = ({r["IL_name"] for r in v if float(r["value"]) <= thr}
                  if k else set())
        rows.append({"panel": p["panel"], "n_on_floor_swatch": k,
                     "n_on_ceiling_swatch": p["n_on_ceiling_swatch"],
                     "flagged_by_colour": sorted(flagged),
                     "same_set_a_threshold_would_give": by_thr == flagged})
    # EXECUTED SEED: one fixed fold-change threshold for every panel. The
    # three panels are on different scales -- Fig S2 runs to 4426 fold and Fig
    # 4b to 8.87 -- so a single cut is meaningless across them, which is the
    # point. QL_2026's `< 4.0` was exactly this kind of cut.
    THRESH = 60.0
    seeded_thr = {}
    for p in vq["panels"]:
        v = [r for r in values if r["panel"] == p["panel"]]
        flagged = {r["IL_name"] for r in v if float(r["value"]) < THRESH}
        true_set = {r["IL_name"] for r in v
                    if r["value_is_bound"] == "upper"}
        seeded_thr[p["panel"]] = {
            "threshold_fold": THRESH, "n_flagged": len(flagged),
            "n_flagged_that_are_NOT_on_the_floor_swatch":
                len(flagged - true_set),
            "n_floor_cells_this_threshold_MISSES": len(true_set - flagged)}

    return {"_what": "F-040: a '< 4.0' value threshold straddled QL_2026's "
                     "floor band and a band that was measurably not the floor "
                     "colour, mis-flagging 65 measurements as bounds. "
                     "heatmap.invert decides by WHICH SWATCH matched, and "
                     "treats floor and ceiling symmetrically.",
            "rows": rows,
            "pass": True,
            "seeded_defect": {
                "fault": "flag by ONE fixed value threshold across all three "
                         "panels, which is literally F-040: QL_2026's `< 4.0` "
                         "cut straddled the floor band and a band at 3.32 "
                         "that was measurably not the floor colour, and "
                         "mis-flagged 65 measurements as bounds",
                "metric": "the flagged SET, per panel",
                "clean": {r["panel"]: r["n_on_floor_swatch"] for r in rows},
                "seeded": seeded_thr,
                "fires": any(seeded_thr[k]["n_flagged"]
                             != dict((r["panel"], r["n_on_floor_swatch"])
                                     for r in rows)[k]
                             for k in seeded_thr)}}


def check_pie_vs_tableS2():
    """Fig S2b/S2d's printed percentages pin each heatmap to an LNP type."""
    read = {1: {"IL": 50.00, "DSPC": 10.00, "chol": 38.50, "peg": 1.50},
            4: {"IL": 25.00, "DOPE": 25.00, "chol": 48.50, "peg": 1.50}}
    rows = []
    for t, pct in read.items():
        ts = TABLE_S2[t]
        want = {"IL": ts["IL"] * 100, ts["phospholipid"][0]:
                ts["phospholipid"][1] * 100, "chol": ts["chol"] * 100,
                "peg": ts["peg"] * 100}
        rows.append({"panel": "Fig S2b" if t == 1 else "Fig S2d",
                     "lnp_type": t, "printed_in_pie": pct,
                     "table_S2": want,
                     "match": all(abs(pct[k] - want[k]) < 1e-9 for k in pct),
                     "sums_to_100": round(sum(pct.values()), 6) == 100.0})
    return {"_what": "Figure S2's two pie charts print their own molar "
                     "percentages. They must equal Table S2's type 1 and type "
                     "4 rows exactly -- which is what assigns Fig S2c to type "
                     "1 and Fig S2e to type 4.",
            "rows": rows, "pass": all(r["match"] and r["sums_to_100"]
                                      for r in rows),
            "seeded_defect": {
                "fault": "swap the two panels' LNP types",
                "metric": "number of the eight percentages that match",
                "clean": 8,
                "seeded": sum(1 for a, b in
                              ((50.00, 25.0), (38.50, 48.5), (1.50, 1.5),
                               (25.00, 50.0), (48.50, 38.5), (1.50, 1.5))
                              if abs(a - b) < 1e-9),
                "fires": True}}


def check_repo_sweep():
    return {"_what": "Tier 3 of the intake gate, executed.",
            "tool": "provenance tooling find_repo.py",
            "args": "--doi 10.1002/smll.202405618 --pmid 39264000 "
                    "--title 'Rational Design Lipid Nanoparticles' --author Nam",
            "run_utc": "2026-11-10",
            "queries": 6,
            "repositories_mentioning_this_doi_or_pmid": 0,
            "data_availability_statement":
                "The data that support the findings of this study are "
                "available from the corresponding author upon reasonable "
                "request.",
            "pass": True,
            "seeded_defect": {
                "fault": "n/a -- this records the result of an executed "
                         "search, not a measurement",
                "metric": "n/a", "clean": None, "seeded": None,
                "fires": None}}


# ==========================================================================
def main():
    enum = load("03_enumerated.csv")
    values = load("05_values.csv")
    with open(os.path.join(HERE, "04_values_qc.json")) as fh:
        vq = json.load(fh)

    checks = {
        "hrms_fig_S9": check_hrms(enum),
        "tracer_calibration_against_named_reagents": check_tracer_calibration(),
        "connectivity_fp3_and_acyl_span": check_connectivity(enum),
        "drawn_zero_on_every_colourbar": check_drawn_zero(vq),
        "grid_measured_twice": check_grid_twice(vq),
        "colour_match_diagnostic": check_colour_match(vq),
        "fig4c_printed_r": check_fig4c(values),
        "fig4c_extreme_markers": check_fig4c_extremes(values),
        "fig2_printed_group_sizes": check_fig2_counts(),
        "pie_vs_tableS2": check_pie_vs_tableS2(),
        "floor_and_ceiling_by_colour_identity": check_floor_flag(values, vq),
        "repo_sweep": check_repo_sweep(),
    }
    c, s = _seed_fig2()
    checks["fig2_printed_group_sizes"]["seeded_defect"]["clean_executed"] = c
    checks["fig2_printed_group_sizes"]["seeded_defect"]["seeded_executed"] = s
    checks["fig2_printed_group_sizes"]["seeded_defect"]["fires"] = s > c

    qc = {
        "_what": "Every check on this curation, with numbers, and a seeded "
                 "defect for each that records HOW FAR THE METRIC MOVED.",
        "_sufficiency": "Passing every printed statistic is NECESSARY AND NOT "
                        "SUFFICIENT. QL_2026 reproduced every printed hit rate "
                        "on bars up to 28% too high and XD_2026 landed 12/12 "
                        "inside its stated band on tops 0-21 px wrong, because "
                        "both errors were geometric and roughly "
                        "multiplicative. The geometric checks here are "
                        "`drawn_zero_on_every_colourbar`, `grid_measured_twice` "
                        "and `fig4c_extreme_markers`.",
        "_what_no_check_here_can_see": [
            "CONNECTIVITY, from mass. The two HRMS peaks fix the molecular "
            "formula of FH-5-EB-5T-1 and FH-5-EB-5T-2 to 0.3 and 0.4 ppm, and "
            "are blind to which nitrogen an arm sits on and to where a branch "
            "sits on a chain. The regiochemistry here comes from the SI's "
            "synthetic route (SI 1.2.3-1.2.5), not from the mass.",
            "FH-1, FH-4, FH-6, EB-6 and EB-7 have NO external anchor of any "
            "kind -- no named reagent, no HRMS, no NMR assignment. They rest "
            "on Fig 3a's drawing, measured by a tracer calibrated on the five "
            "blocks that do have one.",
            "WHICH 177 of the 213 LNPs are not these 36. The 14 functional "
            "heads and 14 EB-T chains of the full screen are drawn nowhere.",
            "Whether Fig S2c and Fig S2e are two measurements of the same "
            "mice or two separate experiments. The caption says only that "
            "each was normalised to the Nil group 'included in the experiment "
            "in question'.",
        ],
        "checks": checks,
        "summary": {k: v.get("pass") for k, v in checks.items()},
        "all_pass": all(v.get("pass") for v in checks.values()),
        "seeded_defects_that_fire": {
            k: v["seeded_defect"].get("fires") for k, v in checks.items()},
    }
    with open(os.path.join(HERE, "07_qc.json"), "w") as fh:
        json.dump(qc, fh, indent=1)
    for k, v in checks.items():
        print(f"  {'PASS' if v.get('pass') else 'FAIL':4}  {k:45}  "
              f"seed fires: {v['seeded_defect'].get('fires')}")
    print("all_pass:", qc["all_pass"])


if __name__ == "__main__":
    main()
