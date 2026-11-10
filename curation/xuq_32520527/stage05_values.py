#!/usr/bin/env python3
"""XQ_2020 stage 05 -- the delivery values, from Figure 5.

    python3 -I stage05_values.py

Figure 5 is a grouped bar chart: eleven amine groups x three formulations,
plus three controls behind an x-axis break. The y axis is 'Delivery efficiency
/ %', linear, with four long tick marks.

EVERYTHING GEOMETRIC HERE COMES FROM curation/lib. The scale and the origin are
fitted together by axis.calibrate from DETECTED tick positions -- there is no
API to supply an origin separately, which is the shape error class 33 took. The
bar tops come from bars.bar_top, whose left/right flank agreement is the defence
against reading an error-bar cap: on this figure the error bars are drawn in the
BAR'S OWN COLOUR, so a colour mask alone does not separate them and the flank
test is doing real work.

HOW THE BAR WINDOWS ARE FOUND, and why not by a literal. Three rules, in order,
each derived from something measured:

  1. a column is a candidate if the colour is present in the single row
     immediately above the detected baseline -- i.e. something of that colour
     is standing on the axis;
  2. candidates are grouped into contiguous runs;
  3. within a run, a column is kept only if its total colour count is at least
     half the run's maximum. This removes JPEG colour ringing at a bar's edge.
     It is not cosmetic: without it the red window of group 5 absorbed one
     column of ringing from the neighbouring black bar and bars.bar_top
     REFUSED it with a 241 px flank disagreement. The refusal was right and
     the window was wrong.

The resulting windows are then checked for uniform width and uniform pitch,
which neither rule 1 nor rule 3 was given.

Writes 05_values.csv and 05_values_qc.json.
"""
import csv
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))

import numpy as np                                                  # noqa: E402
from PIL import Image                                               # noqa: E402
import axis                                                         # noqa: E402
import bars                                                         # noqa: E402
import emit                                                         # noqa: E402

FIG = os.path.join(HERE, "source", "nihms-1799810-f0005.jpg")

# The eleven groups in drawn left-to-right order, read off the x-axis labels of
# Figure 5. The ORDER is not taken on trust: stage 07 re-derives it from five
# printed point values that identify groups 6, 8, 9, 10 and 11 by name.
GROUPS = ["75", "76", "77", "78", "80", "81", "82", "93", "113", "306", "400"]

# Table 1 of the paper, as JATS markup -- not digitised.
SERIES = {
    "black": {"suffix": "",    "lipidoid_chol_dope_weight_ratio": "4/0/0"},
    "red":   {"suffix": "-F1", "lipidoid_chol_dope_weight_ratio": "4/1/1"},
    "blue":  {"suffix": "-F2", "lipidoid_chol_dope_weight_ratio": "4/1/4"},
}

CONDITIONS = {
    "cargo": "GFP mRNA (TriLink)",
    "cell_line": "HeLa",
    "plate": "48-well, 20 000 cells/well in 250 uL DMEM, seeded 24 h before dosing",
    "route": "in vitro, added to culture medium (20 uL of mRNA-loaded LNP per well)",
    "dose_lipidoid": "7.4 ug/mL",
    "dose_mRNA": "0.74 ug/mL",
    "lipidoid_to_mRNA_weight_ratio": "10/1",
    "timepoint": "24 h exposure",
    "assay": "flow cytometry (Invitrogen Attune NxT), percentage of GFP-positive cells",
    "replication": "'Experiments were repeated at least three times'; error bars are s.d.",
    "source_of_all_of_the_above": "RESULTS 'R-O16CBA LNPs for mRNA Delivery', "
                                  "the Figure 5 caption, and EXPERIMENTAL "
                                  "SECTION 'Intracellular Delivery of mRNA' -- "
                                  "all read as JATS text, not from pixels",
}


def colour_masks(rgb):
    R, G, B = rgb[:, :, 0], rgb[:, :, 1], rgb[:, :, 2]
    return {
        # Black ink and the grey control bar are split at 85 rather than at
        # 120. The Lpf2k bar is filled with exactly (128,128,128) but its
        # bottom row blends with the baseline rule to (94,94,94); a grey band
        # starting at 120 missed that row, the baseline seed found no grey
        # column standing on the axis, and the one control bar that carries a
        # printed anchor went unmeasured. Moving the split was verified not to
        # move any of the 33 library values -- see 07_qc.json ->
        # threshold_sensitivity.
        "black": (R < 85) & (G < 85) & (B < 85),
        "red":   (R > 150) & (G < 90) & (B < 90),
        "blue":  (B > 130) & (R < 90) & (G < 90),
        "grey":  (abs(R - G) < 25) & (abs(G - B) < 25) & (abs(R - B) < 25)
                 & (R >= 85) & (R < 215),
    }


def windows(mask, top, base, xlo, xhi, min_width=10, keep_frac=0.5):
    m = np.zeros(mask.shape, bool)
    m[top:base] = mask[top:base]
    m[:, :xlo] = False
    m[:, xhi:] = False
    cnt = m.sum(0)
    seed = [x for x in range(xlo, xhi) if m[base - 1, x]]
    groups, s, prev = [], None, None
    for x in seed:
        if s is None:
            s = x
        elif x != prev + 1:
            groups.append((s, prev + 1))
            s = x
        prev = x
    if s is not None:
        groups.append((s, prev + 1))
    out = []
    for a, b in groups:
        mx = cnt[a:b].max()
        keep = [x for x in range(a, b) if cnt[x] >= keep_frac * mx]
        if keep and keep[-1] + 1 - keep[0] >= min_width:
            out.append((keep[0], keep[-1] + 1))
    return out, m


def main():
    rgb = np.asarray(Image.open(FIG).convert("RGB")).astype(int)
    grey_img = np.asarray(Image.open(FIG).convert("L")).astype(float)
    qc = {"_what": "stage 05 QC for XQ_2020: Figure 5 geometry and the 36 bar "
                   "measurements it yields.",
          "image": {"file": os.path.relpath(FIG, HERE),
                    "shape_hw": list(rgb.shape[:2])}}

    # ---- the axis spine and the baseline, detected ------------------------
    # The spine and the baseline are found on a BLACK-ONLY mask, not on
    # grayscale darkness. Pure red (254,0,0) has luminance ~76 and pure blue
    # ~29, so every coloured bar is "dark"; a row-sum over grayscale darkness
    # makes the bar region compete with the axis rule and the walk that finds
    # the rule's top edge then runs 100 px up into the bars. Measured here:
    # grayscale gives rows 344-348 a count of ~700 against the baseline's 949
    # and the walk terminates at row 245; the black-only mask gives those rows
    # ~240 against 945 and the walk terminates at 349, the real top edge.
    masks = colour_masks(rgb)
    dark = masks["black"]
    spine_x = int(np.argmax(dark.sum(0)))
    baseline_y = int(np.argmax(dark.sum(1)[spine_x:].size and dark.sum(1)))
    # the baseline is the longest horizontal dark run; take the row with the
    # most dark pixels at or below the vertical mid-point of the spine.
    spine_rows = np.where(dark[:, spine_x])[0]
    lower = dark.sum(1).copy()
    lower[: (spine_rows.min() + spine_rows.max()) // 2] = 0
    thick = int(np.argmax(lower))
    # The baseline rule is 2-3 px thick. Take its TOP row: a bar's last row of
    # fill is the row above that, and seeding one row below the top of the rule
    # would seed every column in the image, because the rule itself is ink.
    baseline_y = thick
    while baseline_y > 0 and lower[baseline_y - 1] >= 0.5 * lower[thick]:
        baseline_y -= 1
    baseline_thickness = 1
    while (baseline_y + baseline_thickness < lower.size
           and lower[baseline_y + baseline_thickness] >= 0.5 * lower[thick]):
        baseline_thickness += 1

    # ---- y calibration, from the LONG tick marks only ---------------------
    # The y axis carries eight ticks: four long (labelled 0/30/60/90) and four
    # short minor ticks halfway between them. The long ticks reach further into
    # the plot, so a strip BEYOND the minor ticks' reach selects exactly the
    # four labelled ones. The strip bound is measured, not assumed: it is the
    # right-hand end of a long tick, found by walking right from the spine.
    probe_rows = axis.detect_ticks_1d(grey_img, "y", spine_x + 2, spine_x + 4,
                                      dark=140, min_len=2, max_len=10)
    lengths = []
    for pr in probe_rows:
        r = int(round(pr))
        x = spine_x
        while dark[r, x + 1]:
            x += 1
        # The 0 tick sits ON the drawn baseline, so walking right from it
        # traverses the whole baseline rule and measures ~800 px. That is the
        # baseline's length, not a tick's, and letting it set `long` would
        # push the strip off the end of every real tick. It is excluded by
        # the DETECTED baseline row, not by a length cut-off.
        if abs(r - baseline_y) <= 2:
            continue
        lengths.append(x)
    long_tick_end = int(max(lengths))
    short_tick_end = int(min(lengths))
    if long_tick_end - short_tick_end < 2:
        raise SystemExit(
            f"stage05: all {len(probe_rows)} y ticks reach the same depth "
            f"({short_tick_end}..{long_tick_end}); the long/short split that "
            "picks out the four LABELLED ticks does not exist. Refusing.")
    qc_tick_lengths = {"probe_rows": [round(p, 2) for p in probe_rows],
                       "right_edge_x_excluding_the_baseline_tick": lengths,
                       "long": long_tick_end, "short": short_tick_end,
                       "strip_px": [short_tick_end + 1, long_tick_end]}
    # a strip strictly beyond the short ticks' reach and inside the long ones'
    ticks = axis.detect_ticks_1d(grey_img, "y", short_tick_end + 1,
                                 long_tick_end, dark=140,
                                 min_len=2, max_len=10)
    if len(ticks) != 4:
        raise SystemExit(f"stage05: found {len(ticks)} long y ticks, need 4 "
                         "(the labelled 0/30/60/90). Refusing to calibrate.")
    loo = axis.loo_residuals(ticks)
    cal = axis.calibrate(ticks, [90.0, 60.0, 30.0, 0.0], kind="linear",
                         axis="y", loo_tol_px=axis.LOO_TOL_PX)
    qc["y_axis"] = {
        "spine_x": spine_x, "baseline_y_top_of_rule": baseline_y,
        "baseline_rule_thickness_px": baseline_thickness,
        "tick_length_probe": qc_tick_lengths,
        "tick_rows_px": [round(t, 3) for t in ticks],
        "tick_values": [90, 60, 30, 0],
        "px_per_percent": round(cal.px_per_unit, 4),
        "origin": f"{cal.origin_value} @ {cal.origin_px:.2f} px",
        "max_residual_px": round(cal.residual_px, 4),
        "loo_residuals_px": [round(v, 3) for v in loo],
        "worst_loo_residual_px": round(max(loo), 3),
        "loo_tol_px": axis.LOO_TOL_PX,
        "why_loo_is_reported": "axis.calibrate's max-residual test is weak "
                               "against a SINGLE bad tick -- least squares "
                               "drags the line toward the outlier. A "
                               "leave-one-out residual has no leverage "
                               "masking.",
        "baseline_vs_zero_tick_px": round(
            abs(baseline_y - cal.to_pixel(0.0)), 3),
        "zero_tick_lands_on_the_drawn_baseline": bool(
            abs(baseline_y - cal.to_pixel(0.0)) < 2.0),
    }

    # ---- where the library ends and the controls begin ---------------------
    # The x axis carries a '//' break; right of it are Untreated, naked mRNA
    # and Lpf2k, which are not library members.
    #
    # The widest-gap rule does NOT work here and the failure is worth keeping:
    # the gap between the naked-mRNA bar and the Lpf2k bar (17 px) is wider
    # than the gap between the last blue bar and the Untreated bar (12 px), so
    # "widest gap" put the boundary at x = 954 and swept TWO control bars into
    # the library, giving 12 black windows where there are 11. A rule that is
    # merely plausible is worth no more than the data it happens to fit.
    #
    # The rule used instead is a property of the figure's own encoding: RED and
    # BLUE are the two helper-lipid formulations and appear ONLY on library
    # members -- every control bar is black or grey. So the library ends where
    # the coloured series end. That boundary is then TESTED, not trusted: each
    # of the three series must give exactly 11 windows of uniform width and
    # uniform pitch inside it.
    top = 30                     # below the in-plot legend swatches (rows 7-26)
    base = baseline_y
    foot = np.zeros(rgb.shape[1], bool)
    for nm in ("red", "blue"):
        foot |= masks[nm][base - 1]
    break_x = int(np.where(foot)[0].max()) + 1
    qc["library_region_px"] = [spine_x + 2, break_x]
    qc["how_the_library_region_was_bounded"] = (
        "the right-hand end of the red/blue footprint on the baseline row; "
        "red and blue mark the two helper-lipid formulations and never appear "
        "on a control bar. The widest-gap rule was tried first and was wrong "
        "-- see the comment in stage05_values.py.")
    qc["legend_band_excluded_rows"] = [0, top]

    # ---- the bars ---------------------------------------------------------
    measured, refused = {}, {}
    for nm in ("black", "red", "blue"):
        w, m = windows(masks[nm], top, base, spine_x + 2, break_x)
        if len(w) != len(GROUPS):
            raise SystemExit(f"stage05: {nm} gave {len(w)} bar windows, "
                             f"expected {len(GROUPS)}.")
        widths = [b - a for a, b in w]
        pitch = [w[i + 1][0] - w[i][0] for i in range(len(w) - 1)]
        if max(widths) - min(widths) > 2 or max(pitch) - min(pitch) > 2:
            raise SystemExit(f"stage05: {nm} windows are not uniform: "
                             f"widths {widths}, pitch {pitch}.")
        res, fail = bars.bar_tops(m, w, top, base, baseline_row=base)
        qc.setdefault("series_geometry", {})[nm] = {
            "n_windows": len(w), "widths_px": widths, "pitch_px": pitch,
            "windows": [list(x) for x in w],
            "flank_deltas_px": [res[i]["flank_delta_px"] for i in sorted(res)],
            "worst_flank_delta_px": max(
                (res[i]["flank_delta_px"] for i in res), default=None),
            "n_refused": len(fail),
            "refusals": {str(k): v.splitlines()[0] for k, v in fail.items()},
        }
        for i, (x0, x1) in enumerate(w):
            key = (nm, GROUPS[i])
            if i in fail:
                refused[key] = fail[i]
            else:
                measured[key] = (res[i], x0, x1)

    if refused:
        raise SystemExit(f"stage05: {len(refused)} library bars refused: "
                         f"{list(refused)}")

    # ---- the three controls, measured but NOT emitted ---------------------
    controls = {}
    for nm in ("black", "red", "blue", "grey"):
        w, m = windows(masks[nm], top, base, break_x, rgb.shape[1],
                       min_width=10)
        res, fail = bars.bar_tops(m, w, top, base, baseline_row=base)
        for i, (x0, x1) in enumerate(w):
            if i in res:
                controls[f"{nm}@{x0}-{x1}"] = {
                    "value_pct": round(cal.to_value(res[i]["row"]), 3),
                    "row": res[i]["row"],
                    "flank_delta_px": res[i]["flank_delta_px"]}
            else:
                controls[f"{nm}@{x0}-{x1}"] = {
                    "refused": fail[i].splitlines()[0]}
    qc["controls_right_of_the_axis_break"] = {
        "_what": "Untreated, naked mRNA and Lpf2k. These are NOT corpus rows: "
                 "Lpf2k is a commercial transfection reagent, not a library "
                 "member, and the other two are blanks. They are measured "
                 "only so the Lpf2k bar can test the calibration against the "
                 "printed '~91% GFP+ cells'.",
        "measured": controls,
        "untreated_and_mRNA": "the paper states both are <1% GFP+. At 3.07 px "
                              "per percentage point those bars are ~3 px tall "
                              "and bars.bar_top refuses them as "
                              "baseline-adjacent, which is the correct answer: "
                              "a refusal, not a zero.",
    }

    # ---- assemble the value table ----------------------------------------
    rows, records = [], {}
    for colour, meta in SERIES.items():
        for gid in GROUPS:
            det, x0, x1 = measured[(colour, gid)]
            pct = cal.to_value(det["row"])
            item = f"{gid}-O16CBA{meta['suffix']}"
            records[item] = {
                "item": item, "colour": colour, "group": gid,
                "window_px": [x0, x1], "row_px": det["row"],
                "coverage": det["coverage"],
                "flank_delta_px": det["flank_delta_px"],
                "value_pct": round(pct, 3),
            }
            rows.append({
                "formulation": item,
                "lipidoid": f"{gid}-O16CBA",
                "amine_head": gid,
                "formulation_suffix": meta["suffix"] or "(none)",
                "lipidoid_chol_DOPE_weight_ratio":
                    meta["lipidoid_chol_dope_weight_ratio"],
                "value_pct_GFP_positive": f"{pct:.2f}",
                "log10_value": f"{math.log10(pct):.6f}",
                "source_location":
                    f"Figure 5 (source/nihms-1799810-f0005.jpg), "
                    f"{colour} bar of group '{gid}', "
                    f"x[{x0},{x1}) top row {det['row']}",
                "cargo": CONDITIONS["cargo"],
                "cell_line": CONDITIONS["cell_line"],
                "route": CONDITIONS["route"],
                "dose_lipidoid": CONDITIONS["dose_lipidoid"],
                "dose_mRNA": CONDITIONS["dose_mRNA"],
                "timepoint": CONDITIONS["timepoint"],
                "assay": CONDITIONS["assay"],
            })

    # ---- the within-experiment z of log10, the corpus convention ----------
    lg = np.array([math.log10(float(r["value_pct_GFP_positive"]))
                   for r in rows])
    z = (lg - lg.mean()) / lg.std(ddof=0)
    for r, v in zip(rows, z):
        r["Experiment_value_z_of_log10"] = f"{v:.6f}"

    # ---- the quality columns, through curation/lib/emit -------------------
    enum = {}
    with open(os.path.join(HERE, "03_enumerated.csv")) as f:
        for e in csv.DictReader(f):
            enum[e["lipidoid"]] = e

    def confidence(item, rec):
        """Derived from THIS lipidoid's own enumeration record, never a literal.

        emit.emit_rows refuses a literal `confidence` because a hard-coded
        "high" once relabelled 24 rows whose record read "ASSUMED: the paper
        draws no structure" as the project's highest confidence (F-011). The
        rule here is one line and it reads off `03_enumerated.csv`:

          medium  the head carries a configurational feature the paper draws
                  flat and never states  (structure_assumed = 1)
          high    the head is drawn in full, the arm count is forced by the
                  printed molecular formula, and the arm PLACEMENT is the
                  unique solution of an exhaustive enumeration
        """
        e = enum[f"{rec['group']}-O16CBA"]
        if e["positional_ambiguity"] == "1":
            return "low"
        return "medium" if e["structure_assumed"] == "1" else "high"

    def note(item, rec):
        e = enum[f"{rec['group']}-O16CBA"]
        base_note = (
            "Value digitised from the Figure 5 bar chart, the only place these "
            "numbers are published; the acetal ring's cis/trans configuration "
            "is not drawn anywhere in the paper and the stored SMILES carries "
            "no stereochemistry.")
        if e["structure_assumed"] == "1":
            return ("ASSUMED: amine head 77 is 2-(2-aminoethyl)-1-"
                    "methylpyrrolidine, whose pyrrolidine C2 is a stereocentre "
                    "that Figure 1C draws flat and the paper never specifies; "
                    "the commercial reagent is not stated to be resolved. "
                    + base_note)
        return base_note

    quality = emit.emit_rows(
        items=[r["formulation"] for r in rows],
        records=records,
        confidence=confidence,
        note=note,
        positional_ambiguity=lambda i, rec: int(
            enum[f"{rec['group']}-O16CBA"]["positional_ambiguity"]),
        structure_assumed=lambda i, rec: int(
            enum[f"{rec['group']}-O16CBA"]["structure_assumed"]),
        value_is_bound=lambda i, rec: "",
        value_normalization_unit=lambda i, rec: "experiment",
        label="XQ_2020 Figure 5",
    )
    for r, q in zip(rows, quality):
        r.update(q)

    qc["n_values"] = len(rows)
    qc["rows_with_a_caveat"] = emit.caveat_count(quality)
    qc["measurements"] = records

    with open(os.path.join(HERE, "05_values.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    with open(os.path.join(HERE, "05_values_qc.json"), "w") as f:
        json.dump(qc, f, indent=1)
    print(f"stage05: {len(rows)} values; px/%% {cal.px_per_unit:.4f}; "
          f"max tick residual {cal.residual_px:.3f} px; worst LOO "
          f"{max(loo):.3f} px; {qc['rows_with_a_caveat']} rows with a caveat")
    return 0


if __name__ == "__main__":
    sys.exit(main())
