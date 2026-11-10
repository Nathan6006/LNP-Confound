#!/usr/bin/env python3
"""Digitise Figure S2 and Figure S4 (U266 cell viability, XTT) to test ONE
thing: the two panels report the SAME experimental condition, so they must
agree.

RUN:  python3 -I 06_digitize_S2_S4.py

THE CROSS-PANEL TEST. Figure S2 is "LNPs containing NC-siRNA at 0.06 uM for
72hrs ... XTT". Figure S4's 0.06 uM group is "LNPs containing siNC at different
doses for 72hrs ... XTT". siNC and NC-siRNA are the same negative-control
siRNA, the dose, the cell line, the incubation and the assay are identical, so
the twelve lipids these panels share should carry the same number twice. This
is an UNPRINTED check: the authors never state the comparison, so it cannot be
satisfied by construction the way a recomputed marginal can.

NEITHER VALUE IS CORPUS-ADMISSIBLE -- see 00_value_sources.json. Cell viability
under a negative-control siRNA is a safety readout, not delivery. This panel is
digitised to CHECK THE PAPER, not to emit rows, which is why the output goes to
06_*.json and never to an emit call.

TIER. 5_pixels, and 00_value_sources.json says why tiers 1-4 fail: no Source
Data exists for this DOI, no SI table carries these numbers, no repo/PMC, and
page.get_drawings() returns 0 items on both pages with one DCTDecode JPEG each,
so there is no vector geometry to read.

NO RE-RASTERISATION. The embedded JPEGs are read at their native 1367x874 and
1400x870 rather than re-rendering the page, so no resampling sits between the
authors' pixels and the measurement.
"""
import io
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "lib"))
import numpy as np                                         # noqa: E402
import pymupdf                                             # noqa: E402
from PIL import Image                                      # noqa: E402

import axis                                                # noqa: E402
import bars                                                # noqa: E402
import heatmap                                             # noqa: E402
from errors import CurationLibError                        # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
PDF = os.path.join(HERE, "source", "adma201906128-sup-0001-SuppMat.pdf")

# Panel constants. These live HERE, next to the figure, not in curation/lib --
# PIPELINE_V2 section 2: a paper's own panel bounds stay where a reviewer reads
# them beside the page.
S2 = {"xref": 258, "page": 38,
      "axis_rule_x": 197, "tick_scan": (180, 196),
      "tick_values": [150, 125, 100, 75, 50, 25, 0],
      "bar_rgb": (88, 87, 250), "plot_y": (40, 700)}
S4 = {"xref": 266, "page": 40,
      "axis_rule_x": 138, "tick_scan": (120, 136),
      "tick_values": [150, 125, 100, 75, 50, 25, 0],
      "plot_y": (270, 755),
      # The legend is a 6x2 grid in the top right. These are the SWATCH
      # columns inside that crop; the label text sits at 193-206 / 422-436
      # and had to be excluded by position, because a glyph is non-white too
      # and a threshold alone put text pixels in the swatch median.
      "legend_crop": (0, 226, 840, None),
      "legend_bands": [(42, 64), (72, 93), (101, 122),
                       (131, 152), (160, 181), (189, 210)],
      "legend_swatch_x": ((100, 145), (328, 374)),
      # S4's OWN legend order, read off S4's own legend. It is NOT S3's order:
      # positions 4 and 5 of the right column are swapped between the two
      # panels. See 08_contradictions.csv C-01.
      "legend_order": [(1, 7), (2, 9), (3, 11), (4, 13), (5, 12), (6, 14)],
      "dose_groups": [0.0, 0.012, 0.03, 0.06, 0.12],
      "matched_group_index": 3}          # the 0.06 uM group, matching Fig S2

# Figure S2's x order, read off the printed category labels (p.38).
S2_LABELS = ["UT", "MC3"] + [f"Lipid {i}" for i in range(1, 15)]


def native(xref):
    doc = pymupdf.open(PDF)
    blob = doc.extract_image(xref)
    img = Image.open(io.BytesIO(blob["image"])).convert("RGB")
    return np.array(img).astype(int), blob["ext"]


def colour_mask(rgb_img, target, tol=60):
    """Boolean ink mask from the BAR'S OWN COLOUR, never from darkness.

    bars.bar_top's docstring is explicit about why: a darkness mask makes the
    axis rule itself the topmost 'bar'. Here the bars are blue and every piece
    of furniture -- axis, ticks, labels, error whiskers' caps -- is black.
    """
    d = np.linalg.norm(rgb_img - np.asarray(target, dtype=float), axis=2)
    return d < tol


def bar_windows(mask, y_lo, y_hi, min_width=8):
    """Contiguous column runs where the bar colour is present."""
    cols = mask[y_lo:y_hi, :].sum(axis=0)
    on = cols > 3
    runs, start = [], None
    for i, v in enumerate(on):
        if v and start is None:
            start = i
        elif not v and start is not None:
            if i - start >= min_width:
                runs.append((start, i))
            start = None
    if start is not None and len(on) - start >= min_width:
        runs.append((start, len(on)))
    return runs


def calibrate_panel(gray, scan, values, label):
    lo, hi = scan
    ticks = axis.detect_ticks_1d(gray, "y", lo, hi, dark=140,
                                 min_len=2, max_len=40)
    ticks = sorted(ticks)
    if len(ticks) != len(values):
        raise CurationLibError(
            f"{label}: detected {len(ticks)} ticks, expected {len(values)} "
            f"at {[round(t, 2) for t in ticks]}")
    # values run top-down in pixel order: 150 at the smallest row index.
    cal = axis.calibrate(ticks, values, kind="linear",
                         forbid_target=None, axis="y")
    # The brief's instruction, and the lib's own warning: the max-residual
    # test is weak against a SINGLE bad tick. Report the leave-one-out worst.
    loo = axis.loo_residuals(ticks)
    return cal, ticks, loo


def do_S2():
    img, ext = native(S2["xref"])
    gray = img.mean(axis=2)
    cal, ticks, loo = calibrate_panel(gray, S2["tick_scan"],
                                      S2["tick_values"], "S2")
    mask = colour_mask(img, S2["bar_rgb"])
    y_lo, y_hi = S2["plot_y"]
    wins = [w for w in bar_windows(mask, y_lo, y_hi)
            if w[0] > S2["axis_rule_x"]]
    res, fail = bars.bar_tops(mask, wins, y_lo, y_hi,
                              baseline_row=int(round(cal.to_pixel(0))))
    out = {}
    for i, w in enumerate(wins):
        if i in fail:
            out[S2_LABELS[i] if i < len(S2_LABELS) else f"win{i}"] = {
                "REFUSED": fail[i]}
            continue
        out[S2_LABELS[i] if i < len(S2_LABELS) else f"win{i}"] = {
            "value_pct": round(float(cal.to_value(res[i]["row"])), 2),
            "row": res[i]["row"], "x0": w[0], "x1": w[1],
            "width_px": res[i]["width_px"],
            "flank_delta_px": res[i]["flank_delta_px"]}
    return {"native_px": list(img.shape[:2]), "image_format": ext,
            "calibration": cal.qc(),
            "tick_rows": [round(t, 3) for t in ticks],
            "loo_residuals_px": [round(v, 3) for v in loo],
            "loo_worst_px": round(max(loo), 3),
            "loo_tol_px": axis.LOO_TOL_PX,
            "loo_verdict": ("PASS" if max(loo) <= axis.LOO_TOL_PX
                            else "FAIL -- one tick is mis-detected"),
            "n_bar_windows": len(wins), "n_refused": len(fail),
            "bars": out}


def s4_swatches(img):
    """Figure S4's twelve legend swatches, sampled FROM THIS SAME RASTER.

    heatmap.match requires that and says why: a colour named in a caption is
    not a swatch. Sampled at fixed swatch columns inside each legend band,
    because the label glyphs are non-white too and a pure threshold put text
    pixels into the swatch median on two of the six bands.
    """
    y0, y1, x0, _ = S4["legend_crop"]
    sub = img[y0:y1, x0:]
    out = {}
    for (b0, b1), (lo_n, hi_n) in zip(S4["legend_bands"],
                                      S4["legend_order"]):
        for (sx0, sx1), lipid in zip(S4["legend_swatch_x"], (lo_n, hi_n)):
            px = sub[b0 + 5:b1 - 4, sx0:sx1].reshape(-1, 3)
            out[lipid] = tuple(int(v) for v in np.median(px, axis=0))
    return out


def do_S4():
    img, ext = native(S4["xref"])
    gray = img.mean(axis=2)
    cal, ticks, loo = calibrate_panel(gray, S4["tick_scan"],
                                      S4["tick_values"], "S4")
    sw = s4_swatches(img)
    y_lo, y_hi = S4["plot_y"]
    base = int(round(cal.to_pixel(0)))

    # One mask PER COLOUR. The bars inside a dose group touch, so a single
    # "is a bar" mask would merge all twelve into one window; each lipid's own
    # colour separates them, and it is also the mask bars.bar_top asks for.
    per_lipid, runcounts, refused = {}, {}, {}
    for lipid, rgb in sorted(sw.items()):
        m = colour_mask(img, rgb, tol=28)
        wins = [w for w in bar_windows(m, y_lo, y_hi)
                if w[0] > S4["axis_rule_x"]]
        runcounts[lipid] = len(wins)
        if len(wins) != len(S4["dose_groups"]):
            refused[lipid] = (
                f"colour {rgb} gives {len(wins)} column runs, expected "
                f"{len(S4['dose_groups'])} (one per dose group). Refusing "
                "rather than guessing which run is the 0.06 uM bar.")
            continue
        w = sorted(wins)[S4["matched_group_index"]]
        try:
            r = bars.bar_top(m, w[0], w[1], y_lo, y_hi, baseline_row=base)
        except CurationLibError as e:
            refused[lipid] = str(e)
            continue
        per_lipid[lipid] = {
            "value_pct": round(float(cal.to_value(r["row"])), 2),
            "row": r["row"], "x0": w[0], "x1": w[1],
            "width_px": r["width_px"],
            "flank_delta_px": r["flank_delta_px"],
            "swatch_rgb": list(rgb)}

    # Independent confirmation that each measured bar really is the colour the
    # legend says. heatmap.match refuses on an exact tie and above max_dist.
    match_qc = None
    if per_lipid:
        lips = sorted(per_lipid)
        cells = []
        for lp in lips:
            d = per_lipid[lp]
            xm = (d["x0"] + d["x1"]) // 2
            ym = d["row"] + 12
            cells.append(img[ym, xm])
        order = sorted(sw)
        res = heatmap.match(np.array(cells),
                            np.array([sw[k] for k in order]))
        assigned = [order[i] for i in res["index"]]
        match_qc = {
            "n_cells": len(cells),
            "max_nearest_swatch_distance": round(res["max_distance"], 3),
            "n_ambiguous": res["n_ambiguous"],
            "assignment_agrees_with_mask": assigned == lips,
            "assigned": assigned, "expected": lips}

    return {"native_px": list(img.shape[:2]), "image_format": ext,
            "calibration": cal.qc(),
            "tick_rows": [round(t, 3) for t in ticks],
            "loo_residuals_px": [round(v, 3) for v in loo],
            "loo_worst_px": round(max(loo), 3),
            "loo_tol_px": axis.LOO_TOL_PX,
            "loo_verdict": ("PASS" if max(loo) <= axis.LOO_TOL_PX
                            else "FAIL -- one tick is mis-detected"),
            "legend_swatches": {str(k): list(v) for k, v in sorted(sw.items())},
            "column_runs_per_colour": runcounts,
            "colour_match_qc": match_qc,
            "n_refused": len(refused), "refusals": refused,
            "dose_uM": S4["dose_groups"][S4["matched_group_index"]],
            "bars": {f"Lipid {k}": v for k, v in sorted(per_lipid.items())}}


def s4_legend_swatches_UNUSED(img):
    """Sample Figure S4's TWELVE legend swatches from THIS SAME raster.

    heatmap.match requires that; a colour named in a caption is not a swatch.
    The legend is a 6x2 block in the top right. Rows are found as saturated
    ink bands, and each band's two swatch runs are taken left and right.
    """
    h, w, _ = img.shape
    sub = img[:int(h * 0.26), int(w * 0.60):]
    mx, mn = sub.max(2), sub.min(2)
    m = (mx > 60) & (mx < 255) & ((mx - mn) > 12)
    rows = np.where(m.sum(1) > 8)[0]
    bands, cur = [], [rows[0]]
    for r in rows[1:]:
        if r - cur[-1] <= 2:
            cur.append(r)
        else:
            bands.append(cur)
            cur = [r]
    bands.append(cur)
    swatches = []
    for g in bands:
        band = m[g[0]:g[-1] + 1]
        cols = np.where(band.sum(0) > len(g) * 0.5)[0]
        segs, cc = [], [cols[0]]
        for c in cols[1:]:
            if c - cc[-1] <= 3:
                cc.append(c)
            else:
                segs.append(cc)
                cc = [c]
        segs.append(cc)
        segs = [s for s in segs if len(s) >= 10]
        rgbs = [tuple(int(v) for v in np.median(
            sub[g[0]:g[-1] + 1, s[0]:s[-1] + 1].reshape(-1, 3), 0))
            for s in segs]
        swatches.append({"band_rows": [int(g[0]), int(g[-1])], "rgb": rgbs})
    return swatches


def main():
    out = {"what": "Figures S2 and S4, U266 cell viability (%), XTT",
           "tier": "5_pixels",
           "corpus_admissible": False,
           "purpose": "cross-panel agreement test, not emission"}
    try:
        out["S2"] = do_S2()
    except CurationLibError as e:
        out["S2"] = {"REFUSED": str(e)}

    try:
        out["S4"] = do_S4()
    except CurationLibError as e:
        out["S4"] = {"REFUSED": str(e)}

    # --- the cross-panel test -------------------------------------------
    s2 = out.get("S2", {}).get("bars", {})
    s4 = out.get("S4", {}).get("bars", {})
    shared = sorted(set(s2) & set(s4),
                    key=lambda s: int(s.split()[-1]))
    pairs = [(k, s2[k]["value_pct"], s4[k]["value_pct"]) for k in shared
             if "value_pct" in s2[k] and "value_pct" in s4[k]]
    if pairs:
        a = np.array([p[1] for p in pairs], dtype=float)
        b = np.array([p[2] for p in pairs], dtype=float)
        ratio = b / a
        ra = a.argsort().argsort()
        rb = b.argsort().argsort()
        rho = float(np.corrcoef(ra, rb)[0, 1])
        out["cross_panel_test"] = {
            "claim": ("Figure S2 (NC-siRNA, 0.06 uM, 72 h, U266, XTT) and "
                      "Figure S4's 0.06 uM group (siNC, 72 h, U266, XTT) are "
                      "the same condition and must agree."),
            "n_shared_lipids": len(pairs),
            "shared": [p[0] for p in pairs],
            "S2_values": [round(float(x), 2) for x in a],
            "S4_values": [round(float(x), 2) for x in b],
            "ratio_S4_over_S2": [round(float(x), 4) for x in ratio],
            "ratio_mean": round(float(ratio.mean()), 4),
            "ratio_sd": round(float(ratio.std(ddof=0)), 4),
            "ratio_min": round(float(ratio.min()), 4),
            "ratio_max": round(float(ratio.max()), 4),
            "max_abs_difference_pp": round(float(np.abs(b - a).max()), 2),
            "mean_abs_difference_pp": round(float(np.abs(b - a).mean()), 2),
            "spearman_rank_corr": round(rho, 4),
        }
    with open(os.path.join(HERE, "06_figure_values.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps(out.get("S2", {}), indent=1)[:4000])
    return 0


if __name__ == "__main__":
    sys.exit(main())
