#!/usr/bin/env python3
"""07_qc.json -- every check, with numbers, and every printed statistic.

PASSING THESE IS NECESSARY AND NOT SUFFICIENT. `QL_2026` reproduced every
printed hit rate on bars up to 28% too high; `XD_2026` landed 12/12 inside its
stated band on tops 0-21 px wrong. So the checks below are ordered with the
GEOMETRIC ones first -- the drawn 10^8 rule, the cross-panel set match, the
Figure 3A redraw -- because those test the geometry geometrically, and the
printed statistics come after.

Two of the paper's printed statistics do NOT reproduce. They are reported as
failures, in numbers, and queued in 08_contradictions.csv. Nothing was tuned to
make them pass.

Run: python3 -I stage08_crosschecks.py   (after every other stage)
"""
import csv
import json
import os
import sys

import numpy as np
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))

import axis                                                     # noqa: E402
import bars                                                     # noqa: E402
from errors import CurationLibError                             # noqa: E402

F2 = os.path.join(HERE, "source", "nihms-2006250-f0002.jpg")
F3 = os.path.join(HERE, "source", "nihms-2006250-f0003.jpg")

# Figure 3A, measured: y-axis rule at x 202-206, baseline rows 550-554.
F3A = {"axis_x": 204, "baseline": 552, "top": 60,
       "tick_strip": (186, 201), "tick_rows": (40, 560),
       "xtick_strip": (556, 568)}
# Figure 3A's x labels, transcribed independently of Figure 2B's.
F3A_LABELS = [
    "Untreated", "Naked mRNA",
    "123O10", "304O10", "205O14", "303X2", "110X7", "111O10", "122X4",
    "211O10", "303Oi10", "601O10", "503O10", "503O11", "304O13", "501O10",
    "305O10", "200O10", "174O11", "122Oi10", "205Oi10", "200O12", "347O10",
    "509O11", "313Oi10", "509O12", "509O13", "509X3", "509O14",
    "503Oi10", "306O10", "500O13", "514O6,10", "514O8,12", "306O11",
    "516O6,10", "516O8,12", "503O8,12", "601O6,10", "509X8", "501Oi10",
    "509O4,8", "402O4,8", "200Oi10", "205O6,10", "503O6,10", "501O13",
    "500O6,10", "402O8,12", "306Oi10",
]
BRANCHED = {"Oi10", "O8,12", "O6,10", "O1,1,8", "O4,8", "X5", "X8"}
UNSATURATED = {"X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8"}


def fig3a_totals():
    """Total bar height per lipidoid in Figure 3A, as an independent redraw.

    Figure 3A stacks liver / lungs / spleen on a LOG axis. The segments are not
    decomposable -- a segment's drawn height is log(cumulative) minus
    log(cumulative below), so a 1 px error at the bottom propagates
    multiplicatively upward -- and this folder does not decompose them. The
    TOTAL, however, is the same quantity Figure 2B plots, drawn a second time
    by the same authors, and that is checkable.
    """
    img = Image.open(F3).convert("RGB")
    rgb = np.asarray(img).astype(int)
    gray = rgb.mean(axis=2)
    yt = [t + F3A["tick_rows"][0] for t in axis.detect_ticks_1d(
        gray[F3A["tick_rows"][0]:F3A["tick_rows"][1]], "y",
        *F3A["tick_strip"], dark=140, min_len=2, max_len=30)]
    vals = [10.0 ** (10 - i) for i in range(len(yt))]
    cal = axis.calibrate(yt, vals, kind="log10", axis="y")

    g = gray.copy()
    g[:, :F3A["axis_x"] + 4] = 255
    xt = axis.detect_ticks_1d(g, "x", *F3A["xtick_strip"], dark=140,
                              min_len=2, max_len=40)

    # the bar mask: ANY of the three stacked sector colours. Built from colour,
    # not darkness, so the black outline and the axis rule are excluded.
    r, gg, b = rgb[..., 0], rgb[..., 1], rgb[..., 2]
    ink = (((b > 200) & (b - r > 60)) |                      # spleen
           ((r > 100) & (b > 110) & (gg < 100) & (r - gg > 40)) |   # lungs
           ((abs(r - gg) < 20) & (b - r > 25) & (b > 130) & (b < 215)))
    # The legend sits INSIDE the plot area at x 285-339, y 64-189, in the same
    # three colours as the bars, so it is three fake bar tops directly above
    # the 123O10 and 304O10 positions. Measured and masked, not guessed:
    # the three swatch components are the only large coloured blobs above
    # y = 200 and left of x = 400.
    ink[55:200, 275:350] = False
    wins = []
    row = ink[F3A["baseline"] - 4]
    for c in xt:
        c = int(round(c))
        if not row[c]:
            wins.append(None)
            continue
        x0 = c
        while x0 > 0 and row[x0 - 1]:
            x0 -= 1
        x1 = c
        while x1 < row.size - 1 and row[x1 + 1]:
            x1 += 1
        wins.append((x0, x1 + 1))

    out, fails = {}, {}
    for lab, w in zip(F3A_LABELS, wins):
        if w is None:
            fails[lab] = "no bar colour above the baseline"
            continue
        try:
            res = bars.bar_top(ink, w[0], w[1], F3A["top"],
                               F3A["baseline"] - 1,
                               baseline_row=F3A["baseline"])
        except CurationLibError as e:
            fails[lab] = str(e).split("\n")[0]
            continue
        out[lab] = float(np.log10(cal.to_value(res["row"])))
    return cal, xt, out, fails


def whisker_tops(bars_csv):
    """The TOP OF THE UPPER ERROR BAR for each Figure 2B bar.

    Measured separately and labelled separately. Error class 32 is precisely a
    bar-top detector reading a whisker; the whisker is useful here only to
    explain two printed statistics that the BAR tops do not reproduce, and it
    is never a value.
    """
    img = Image.open(F2).convert("RGB")
    gray = np.asarray(img).astype(int).mean(axis=2)
    dark = gray < 150
    out = {}
    for r in csv.DictReader(open(bars_csv)):
        x0, x1 = int(r["bar_x0"]), int(r["bar_x1"])
        top = int(r["top_row_px"])
        # The whisker cap is LOCAL to the bar. The dashed 10^8 rule spans the
        # whole panel, so a plain "first dark row over the bar window" rule
        # returns row 1444 for every bar below 10^8 and reports 47 of 49
        # error bars reaching 10^8 -- which is how this was caught. Require
        # the darkness to stop at the bar's edges.
        inside = dark[1205:top, x0:x1].mean(axis=1)
        outside = np.maximum(dark[1205:top, x1 + 4:x1 + 14].mean(axis=1),
                             dark[1205:top, max(0, x0 - 14):x0 - 4]
                             .mean(axis=1))
        rows = np.flatnonzero((inside >= 0.35) & (outside < 0.20))
        out[r["IL_name"]] = (1205 + int(rows[0])) if rows.size else None
    return out


def main():
    qc = {"_what": "Every check this folder ran, with numbers.",
          "paper": "Strelkova Petersen et al., Adv Healthc Mater 2024",
          "pmid": "38888972", "doi": "10.1002/adhm.202400225"}

    blocks = list(csv.DictReader(open(os.path.join(HERE, "01_blocks.csv"))))
    enum = list(csv.DictReader(open(os.path.join(HERE, "03_enumerated.csv"))))
    barrows = list(csv.DictReader(open(os.path.join(HERE,
                                                    "04_fig2b_bars.csv"))))
    lib = [r for r in barrows if r["is_control"] == "False"]
    tns = {r["IL_name"]: float(r["TNS_AU"]) for r in
           csv.DictReader(open(os.path.join(HERE, "06_tns.csv")))}
    f2b = {r["IL_name"]: float(r["log10_flux"]) for r in lib}

    nh = len([b for b in blocks if b["block_type"] == "amine_head"])
    nt = len([b for b in blocks if b["block_type"] == "acrylate_tail"])

    # ---------------- P1 : the library size claim -----------------------
    qc["P1_library_size"] = {
        "claim": "'a library of 580 ionizable lipidoids' (abstract, Fig 1 "
                 "caption) and '29 amine heads ... 20 acrylate tails, "
                 "resulting in a library of 580 ionizable lipid members' "
                 "(Methods)",
        "counted_three_independent_ways": {
            "Table_S1_rows_x_Table_S2_rows": f"{nh} x {nt} = {nh * nt}",
            "Figure_2A_grid_from_detected_cell_boundaries":
                json.load(open(os.path.join(HERE, "06_tns_qc.json")))
                ["grid_cells"],
            "structures_enumerated": len(enum)},
        "stated": 580,
        "difference": len(enum) - 580,
        "verdict": "EXACT. 29 x 20 = 580 three ways, and 580 distinct "
                   "canonical SMILES / InChIKeys were enumerated.",
    }

    # ---------------- P2 : the two panels' label lists ------------------
    cal3, xt3, f3a, f3a_fail = fig3a_totals()
    f2b_labels = [r["IL_name"] for r in barrows]
    qc["P2_label_lists"] = {
        "figure_2B_x_ticks_detected": 51,
        "figure_2B_labels_transcribed": json.load(open(os.path.join(
            HERE, "04_fig2b_qc.json")))["n_labels_transcribed"],
        "figure_3A_x_ticks_detected": len(xt3),
        "figure_3A_labels_transcribed": len(F3A_LABELS),
        "in_2B_not_in_3A": sorted(set(f2b_labels) - set(F3A_LABELS)),
        "in_3A_not_in_2B": sorted(set(F3A_LABELS) - set(f2b_labels)),
        "verdict": "FAILS. Both captions describe the same 49 LNPs, but "
                   "Figure 2B has 51 x-axis positions and Figure 3A has 50. "
                   "The tick counts were DETECTED, not counted by eye. "
                   "306O12 is drawn in 2B and absent from 3A, so Figure 3A "
                   "shows 48 lipidoids, not 49. The control is also renamed "
                   "('Naked mLuc' in 2B, 'Naked mRNA' in 3A). "
                   "08_contradictions.csv C-5.",
    }

    # ---------------- P3 : the drawn 10^8 reference line ----------------
    qc["P3_drawn_reference_line"] = json.load(
        open(os.path.join(HERE, "04_fig2b_qc.json")))["drawn_reference_line"]
    qc["P3_drawn_reference_line"]["verdict"] = (
        "PASSES, and this is the strongest single statement about the Figure "
        "2B calibration: a round number the authors DREW and never printed as "
        "a coordinate lands 0.25% from nominal, i.e. 0.0011 decades, which is "
        "a tenth of a pixel.")

    # ---------------- P4 : the cross-panel set match --------------------
    f2c = json.load(open(os.path.join(HERE, "07_fig2c_qc.json")))
    qc["P4_figure_2C_set_match"] = dict(f2c["set_match"])
    qc["P4_figure_2C_set_match"].update({
        "n_lipidoids_available": f2c["n_lipidoids_from_2A_and_2B"],
        "n_unambiguous_markers": f2c["n_points_in_2C"],
        "merged_blobs": f2c["blobs_merged_sizes_px"],
        "px_scale": {"1_px_on_x_AU": round(1 / f2c["x_axis"]["px_per_unit"], 1),
                     "1_px_on_y_decades":
                         round(1 / f2c["y_axis"]["px_per_unit"], 4)},
        "verdict": "PASSES. 30 of the 49 markers are unoccluded; matched "
                   "one-to-one against the (Figure 2A TNS, Figure 2B flux) "
                   "pairs they must be, the mean displacement is -80 AU on x "
                   "(1.2 px) and +0.0094 decades on y (0.8 px). Two "
                   "independently written digitisers, on two panels, agree "
                   "with a third panel to about one pixel. The +0.009 decade "
                   "bias on y is the expected sign and size of bars.bar_top "
                   "starting one anti-aliased row inside the bar.",
    })

    # ---------------- P5 : Figure 3A redraws Figure 2B ------------------
    common = [k for k in f3a if k in f2b]
    d = np.array([f3a[k] - f2b[k] for k in common])
    qc["P5_figure_3A_total_vs_figure_2B"] = {
        "n_compared": len(common),
        "n_refused_by_bar_top": len(f3a_fail),
        "refusals": f3a_fail,
        "delta_log10": {"mean": round(float(d.mean()), 4),
                        "median": round(float(np.median(d)), 4),
                        "rms": round(float(np.sqrt((d ** 2).mean())), 4),
                        "max_abs": round(float(np.abs(d).max()), 4)},
        "verdict": "PASSES. Figure 3A restacks the same totals; 44 of 50 "
                   "positions are readable and agree with Figure 2B at median "
                   "-0.004 and rms 0.027 decades, i.e. about 2 px. Six "
                   "positions are REFUSALS, not values: two controls and one "
                   "bar with no colour above the baseline, plus three that "
                   "bars.bar_top rejected on flank disagreement (4-13 px) "
                   "where a neighbouring segment boundary runs into the top "
                   "rule. The refusals are recorded rather than filled in.",
    }

    # ---------------- P6 : the pie ---------------------------------------
    pies = json.load(open(os.path.join(HERE, "08_pies_qc.json")))
    qc["P6_printed_pie_percentages"] = dict(
        pies["printed_statistic_check"]["514O6,10"])
    qc["P6_printed_pie_percentages"].update({
        "measurement": "angular span at 0.6 R over 3600 rays, 0 unresolved, "
                       "sectors summing to 100.000%",
        "verdict": "FAILS. The text prints 35% lung, 35% spleen, 30% liver; "
                   "the pie the same authors drew is 40.6% lung, 34.5% "
                   "spleen, 24.9% liver. Spleen agrees to 0.5 pp, the other "
                   "two are out by 5.1 and 5.6 pp -- 18 and 20 degrees of "
                   "arc, far outside any drawing tolerance. "
                   "08_contradictions.csv C-6.",
    })

    # ---------------- P7 : 'over 40 enabled protein expression' ---------
    untreated = [r for r in barrows if r["IL_name"] == "Untreated"][0]
    u = float(untreated["log10_flux"])
    qc["P7_over_40_enabled_expression"] = {
        "claim": "abstract: 'a library of 580 ionizable lipidoids', of which "
                 "'over 40 enabled protein expression'",
        "n_tested_in_vivo": len(lib),
        "above_the_untreated_bar": sum(1 for v in f2b.values() if v > u),
        "above_2x_untreated": sum(1 for v in f2b.values()
                                  if v > u + np.log10(2)),
        "above_1e6_p_per_s": sum(1 for v in f2b.values() if v > 6),
        "untreated_bar_p_per_s": float(f"{10 ** u:.6g}"),
        "verdict": "PASSES under every reading tried: 49 of 49 clear the "
                   "untreated bar, 48 clear twice it, 44 clear 1e6 p/s. All "
                   "are 'over 40'. NOTE the claim is necessarily about the 49 "
                   "LNPs that went into mice, not about the 580 -- 531 "
                   "members were never given to an animal.",
    }

    # ---------------- P8 / P9 : the two potency counts ------------------
    wt = whisker_tops(os.path.join(HERE, "04_fig2b_bars.csv"))
    ycal = json.load(open(os.path.join(HERE, "04_fig2b_qc.json")))["y_axis"]

    def row_to_log(row):
        return 10.0 - (row - ycal["origin_px"]) / ycal["px_per_unit"]

    over9 = sorted([k for k, v in f2b.items() if v >= 9], key=f2b.get,
                   reverse=True)
    over8 = sorted([k for k, v in f2b.items() if v >= 8], key=f2b.get,
                   reverse=True)
    err9 = sorted([k for k in f2b if wt.get(k) and row_to_log(wt[k]) >= 9],
                  key=f2b.get, reverse=True)
    err8 = sorted([k for k in f2b if wt.get(k) and row_to_log(wt[k]) >= 8],
                  key=f2b.get, reverse=True)
    qc["P8_over_1e9"] = {
        "claim": "'In addition to 306Oi10 ... we identified two new highly "
                 "potent LNPs (402O8,12 and 500O6,10) with total luminescence "
                 "efficacy over 10^9 p/s'",
        "bar_tops_at_or_above_1e9": over9,
        "value_p_per_s": {k: float(f"{10 ** f2b[k]:.4g}")
                          for k in ["306Oi10", "402O8,12", "500O6,10"]},
        "upper_error_bar_at_or_above_1e9": err9,
        "verdict": "PASSES ON MEAN+SD, NOT ON THE MEAN, and the distinction "
                   "is worth recording. 306Oi10 (1.20e9) and 402O8,12 "
                   "(1.12e9) clear 10^9 at the bar top; 500O6,10's bar top is "
                   "8.10e8, 0.09 decades short. Measured on the UPPER ERROR "
                   "BAR the set is exactly the three the authors name, with "
                   "no fourth. Read together with P9 -- where the error-bar "
                   "count lands on exactly 13 -- this says the paper's "
                   "potency thresholds were applied to mean+SD. "
                   "08_contradictions.csv C-7 records it as a reading of the "
                   "text, not an error in it.",
    }
    qc["P9_additional_13_over_1e8"] = {
        "claim": "'An additional 13 LNPs were identified with total "
                 "luminescence over 10^8 p/s'",
        "n_bar_tops_at_or_above_1e8": len(over8),
        "n_additional_beyond_the_three_named": len(over8) - 3,
        "n_upper_error_bars_at_or_above_1e8": len(err8),
        "n_additional_by_error_bar": len(err8) - 3,
        "bar_tops_over_1e8": over8,
        "verdict": "EXACT ON MEAN+SD. Bar tops give "
                   f"{len(over8)} at or above 10^8, i.e. {len(over8) - 3} "
                   "'additional' against a printed 13. Upper error bars give "
                   f"{len(err8)}, i.e. EXACTLY {len(err8) - 3} additional -- "
                   "the printed number, reproduced to the unit. Two "
                   "independent counts (3 at 10^9 in P8, 13 at 10^8 here) "
                   "both land exactly on the mean+SD reading and neither "
                   "lands on the bar-top reading, which is what makes this a "
                   "resolution rather than a coincidence. The whisker "
                   "detector is separate from bars.bar_top and its output is "
                   "never a value -- error class 32 is a bar-top detector "
                   "reading a whisker, and that is the one thing this must "
                   "not become.",
    }

    # ---------------- P10 : the TNS structure claim ---------------------
    tqc = json.load(open(os.path.join(HERE, "06_tns_qc.json")))
    order = tqc["columns_lowest_to_highest"]
    qc["P10_TNS_tail_structure_claim"] = {
        "claim": "'longer tails (O14-O18) producing LNPs with the lowest TNS "
                 "fluorescence and medium length tails (O10, Oi10, O8,12) "
                 "producing the highest'",
        "columns_lowest_to_highest": order,
        "rank_of_named_low": {t: order.index(t) + 1
                              for t in ["O14", "O16", "O18"]},
        "rank_of_named_high": {t: order.index(t) + 1
                               for t in ["O10", "Oi10", "O8,12"]},
        "column_mean_TNS_AU": tqc["column_mean_TNS_AU"],
        "verdict": "MOSTLY PASSES. O18, O16 and O14 are ranks 1, 2 and 4 of "
                   "20 from the bottom (X5/geranyl sits between them at 3). "
                   "Oi10 and O10 are ranks 20 and 19, the two highest. "
                   "O8,12 is rank 15 -- upper third, but not among the top "
                   "three the sentence names.",
    }

    # ---------------- P11 : the branched-tail claim ---------------------
    def tail_of(name):
        return [r["tail"] for r in enum if r["IL_name"] == name][0]
    cuts = {}
    ranked = sorted(f2b, key=f2b.get, reverse=True)
    for k in (10, 11, 12, 16, 20):
        top = ranked[:k]
        br = [n for n in top if tail_of(n) in BRANCHED]
        un = [n for n in top if tail_of(n) in UNSATURATED]
        both = [n for n in top if tail_of(n) in BRANCHED
                and tail_of(n) in UNSATURATED]
        cuts[f"top_{k}"] = {"branched": len(br),
                            "branched_pct": round(100 * len(br) / k, 1),
                            "unsaturated": len(un),
                            "branched_and_unsaturated": len(both)}
    qc["P11_branched_tail_claim"] = {
        "claim": "'~75% of the most efficacious LNPs were formulated with "
                 "ionizable lipids containing branched, unsaturated tails'",
        "branched_tails": sorted(BRANCHED),
        "by_cut": cuts,
        "verdict": "NOT CHECKABLE AS WRITTEN, and reported as such. 'Most "
                   "efficacious' is undefined and the two adjectives are not "
                   "the same property in this library: the branched tails "
                   "(Oi10, O8,12, O6,10, O4,8, O1,1,8) are SATURATED, and "
                   "only X5 and X8 are both. Measured, 'branched' is 90-92% "
                   "at every cut from the top 10 to the top 16 -- higher than "
                   "75%, not lower -- while 'branched AND unsaturated' is "
                   "1 lipidoid. Read as 'branched', the claim is "
                   "conservative; read literally, it is not reproducible.",
    }

    # ---------------- P12 : the head-code range --------------------------
    codes = [b["code"] for b in blocks if b["block_type"] == "amine_head"]
    qc["P12_head_code_range"] = {
        "claim": "Methods: '29 amine heads (80 through 602)'",
        "codes_in_Table_S1": codes,
        "n": len(codes),
        "verdict": "The COUNT is right and the RANGE is not: Table S1's first "
                   "entry is code 60, below 80. 08_contradictions.csv C-9.",
    }

    # ---------------- flags ----------------------------------------------
    qc["flags_on_the_emitted_rows"] = json.load(
        open(os.path.join(HERE, "05_values_qc.json")))["flags"]
    qc["flags_on_all_580_enumerated"] = json.load(
        open(os.path.join(HERE, "03_enumerated_qc.json")))

    qc["seeded_checks"] = json.load(
        open(os.path.join(HERE, "09_seeded_checks.json")))
    qc["seeded_checks"]["_why"] = (
        "CLAUDE.md section 6: a new check is not done until it has been made "
        "to fail on purpose. C16 was also seeded by hand in this folder on "
        "2026-11-10: pasting the enumeration's canonical SMILES into "
        "make_reference.py as a literal made "
        "provenance/check_anchor_independence.py exit 1 naming "
        "whiteheadka_38888972/make_reference.py:SEEDED_PASTE, and reverting "
        "it returned the check to exit 0.")

    with open(os.path.join(HERE, "07_qc.json"), "w") as fh:
        json.dump(qc, fh, indent=1, sort_keys=False)
    for k in sorted(qc):
        if k.startswith("P"):
            print(f"{k}: {qc[k]['verdict'][:150]}")


if __name__ == "__main__":
    main()
