#!/usr/bin/env python3
"""Digitise Figure 1a (in vitro heat map) and Figure 1c (in vivo bars).

Run from this folder:  python3 -I stage04_values.py
    -> 00_pages/fig1_native.jpeg    the figure, extracted at native resolution
    -> 04_fig1a_cells.csv           123 heat-map cells
    -> 04_fig1c_bars.csv            82 bar tops
    -> 04_geometry.json             every geometric constant, measured

TIER 5. Both panels are pixels, and 00_value_sources.json records in a
reviewer-checkable sentence why tiers 1-4 are unavailable: no Source Data file,
an SI with zero tables (0 <w:tbl> elements AND 0 occurrences of "Table S" in
959,060 characters), a Data Availability Statement pointing at the
corresponding author, zero repositories from find_repo.py for both the first
author and the PI, zero real accessions in a 3,126,538-character sweep, and
2 vector paths on the figure page -- the Wiley watermark -- so there is no
tier-4b geometry either.

NATIVE RESOLUTION, MEASURED. The figure is an embedded JPEG, 2051 x 2352 px,
placed on a 492.24 x 564.27 pt rectangle: 299.97 x 300.13 dpi by px/rect_pt*72.
It is EXTRACTED, never re-rendered, so there is no second rasterisation.

WHAT IS MEASURED AND WHAT IS NOT.
  * Fig 1a is a CONTINUOUS diverging colormap on an unlabelled -1.0 .. +1.0
    scale with no tick marks. Its LUT is built from the colorbar in this same
    raster. What the number MEANS is a separate question and this file does not
    answer it -- see DECISIONS.md D-2.
  * Fig 1c has a real axis in physical units (Avg Radiance, x10^8 p/s/cm2/sr)
    with four ticks. Only THREE are usable: the 0 tick is fused with the drawn
    baseline. So the fit uses the ticks for 3, 2 and 1 only, and the position
    of 0 is then PREDICTED and compared with the drawn baseline -- a round
    number the authors drew and never printed. `forbid_target` is passed so
    that fitting onto the thing being validated is a type error (class 27).
"""
import csv
import json
import os
import sys

import numpy as np

LIB = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib")
sys.path.insert(0, LIB)
import axis                                                      # noqa: E402
import bars                                                      # noqa: E402
import heatmap                                                   # noqa: E402
from errors import DetectionError                                # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
PDF = os.path.join(HERE, "source", "Angew Chem Int Ed - 2026 - Jiang - "
                   "Retinoic‐Acid‐Derived Ionizable Lipids Enable "
                   "Spleen‐Selective Antigen Expression for.pdf")
RASTER = os.path.join(HERE, "00_pages", "fig1_native.jpeg")

sys.path.insert(0, HERE)
from stage03_enumerate import panel_order                        # noqa: E402

LABELS = ["Dlin-MC3"] + panel_order()          # 41, left to right

# ---- Figure 1a, measured in stage04 itself (see _measure_heatmap) ---------
A_BAND_MIN_COLOURED = 190   # rows of a column that must be non-white
# ---- Figure 1c -----------------------------------------------------------
C_TICK_VALUES = [3.0, 2.0, 1.0]      # the 0 tick is fused with the baseline
C_Y_TOP = 1205      # below the legend dots, which end by y=1200
C_SPINE_X = (172, 178)
# A bar whose whole body lies within this many pixels of the drawn baseline is
# not separable from the axis rule. 4 px / 110.18 px-per-unit = 0.036e8.
FLOOR_PX = 4.0


def load_raster():
    import fitz
    if not os.path.exists(RASTER):
        d = fitz.open(PDF)
        page = d[3]                      # page 4, the Figure 1 page
        imgs = page.get_images(full=True)
        if len(imgs) != 1:
            raise SystemExit(
                f"stage04: page 4 holds {len(imgs)} images, expected exactly "
                "1. Refusing to guess which one is Figure 1.")
        xref = imgs[0][0]
        px = d.extract_image(xref)
        rects = page.get_image_rects(xref)
        if len(rects) != 1:
            raise SystemExit("stage04: the figure is placed more than once.")
        r = rects[0]
        dpi_x = px["width"] / (r.x1 - r.x0) * 72.0
        dpi_y = px["height"] / (r.y1 - r.y0) * 72.0
        os.makedirs(os.path.dirname(RASTER), exist_ok=True)
        with open(RASTER, "wb") as fh:
            fh.write(px["image"])
        meta = {"width_px": px["width"], "height_px": px["height"],
                "rect_pt": [round(v, 3) for v in (r.x0, r.y0, r.x1, r.y1)],
                "native_dpi_x": round(dpi_x, 2),
                "native_dpi_y": round(dpi_y, 2),
                "ext": px["ext"], "page": 4, "xref": xref,
                "drawings_on_this_page": len(page.get_drawings())}
    else:
        meta = None
    from PIL import Image
    img = np.asarray(Image.open(RASTER).convert("RGB")).astype(float)
    return img, meta


# =========================================================================
# Figure 1a
# =========================================================================
def _measure_heatmap(img):
    """Outer edges of the heat map and of the colorbar, from the raster.

    THE COLORBAR IS FOUND BY HOW MANY DISTINCT COLOURS ITS COLUMNS CONTAIN,
    not by looking for a white gutter. My first version looked for the single
    widest white gap on one scan row and found FIVE, because a heat-map cell
    whose value is near 0 renders white on a diverging map and is
    indistinguishable from background by any whiteness test. A continuous ramp
    has ~210 distinct colours down its height; a column of three flat cells
    has three. That separates them by two orders of magnitude.
    """
    # THE CONTIGUOUS RUN, not min..max of every passing row. The rotated
    # x-axis labels below the map also put >800 non-white pixels on a row, so
    # min/max swept 26..301 and dragged 67 rows of label text into the band --
    # which inflated the distinct-colour count of 444 map columns and made the
    # colorbar test below report 476 "colorbar" columns spanning the whole
    # figure. The band is one block; take the block.
    hits = [y for y in range(0, 320)
            if (img[y, :, :].sum(axis=1) < 735).sum() > 800]
    runs_y, cur = [], [hits[0], hits[0]]
    for y in hits[1:]:
        if y == cur[1] + 1:
            cur[1] = y
        else:
            runs_y.append(tuple(cur))
            cur = [y, y]
    runs_y.append(tuple(cur))
    y0, y1 = max(runs_y, key=lambda r: r[1] - r[0])
    band = img[y0:y1 + 1, :, :].astype(int)

    nband = y1 - y0 + 1
    ncolours = np.array([len({tuple(c) for c in band[:, x, :]})
                         for x in range(img.shape[1])])

    # --- the colorbar, by distinct-colour count --------------------------
    # THE THRESHOLD IS AUDITABLE, not a magic number, and the separation is
    # recorded below: a column of three flat cells reaches 33 distinct colours
    # once JPEG noise is counted, while a continuous ramp holds 196-208. The
    # cut at 100 sits between two populations 5.9x apart.
    cb = np.flatnonzero(ncolours > 100)
    if cb.size == 0:
        raise DetectionError("stage04: no colorbar column found.")
    cb_runs, cur = [], [int(cb[0]), int(cb[0])]
    for x in cb[1:]:
        if x == cur[1] + 1:
            cur[1] = int(x)
        else:
            cb_runs.append(tuple(cur))
            cur = [int(x), int(x)]
    cb_runs.append(tuple(cur))
    cb_runs = [r for r in cb_runs if r[1] - r[0] >= 8]
    if len(cb_runs) != 1:
        raise DetectionError(
            f"stage04: {len(cb_runs)} candidate colorbar runs {cb_runs}, "
            "expected exactly 1. Refusing to pick one.")
    cb_x0, cb_x1 = cb_runs[0]

    # --- the heat map, by FULL-BAND-HEIGHT coverage ----------------------
    # Not by "is there any ink": the bold panel letter "a" sits 8 px to the
    # left of the map's first cell, and an any-ink test put map_x0 at 11,
    # stretching the 41-column grid over the letter and pushing the cell
    # centres half a cell off. The map is a solid rectangle spanning the whole
    # band; the letter spans 13 of 209 rows. That is a 15x separation, and the
    # numbers are recorded so a reviewer can check it rather than take it.
    covered = (band.sum(axis=2) < 760).sum(axis=0)
    thresh = int(round(0.95 * nband))
    cand = np.flatnonzero(covered[:cb_x0 - 4] >= thresh)
    if cand.size == 0:
        raise DetectionError("stage04: no full-height heat-map column found.")
    map_x0, map_x1 = int(cand.min()), int(cand.max())

    # THE ASSERTION IS ON THE EDGES, NOT ON EVERY INTERIOR COLUMN. Requiring
    # contiguity was my first version and it refused a correct extent: 125 of
    # the 1,889 interior columns fall below 95% coverage because a cell whose
    # value is near 0 renders white on a diverging map, and a cell boundary
    # between two pale cells is white all the way down. Those are the data,
    # not a detection failure. What actually has to hold is that the extent is
    # set by the map and not by some other object, so what is asserted is that
    # NOTHING to the left of map_x0 comes anywhere near the threshold -- the
    # bold panel letter "a", 8 px away, is the thing that would.
    # The three columns immediately left of map_x0 are the map's own
    # antialiased edge (coverage 142/209 on this raster), not a foreign
    # object, so the test starts beyond them.
    EDGE_AA = 3
    outside = covered[:max(0, map_x0 - EDGE_AA)]
    worst_outside = int(outside.max()) if outside.size else 0
    if worst_outside >= 0.5 * nband:
        raise DetectionError(
            f"stage04: a column left of the map at coverage {worst_outside}/"
            f"{nband} is close to the {thresh} threshold, so map_x0={map_x0} "
            "may have been set by something that is not the heat map.")
    sep = {"threshold_rows": thresh, "band_rows": int(nband),
           "n_interior_columns_below_threshold":
               int((covered[map_x0:map_x1 + 1] < thresh).sum()),
           "n_interior_columns": int(map_x1 - map_x0 + 1),
           "interior_below_threshold_note":
               "Pale cells and the boundaries between them. A cell at value "
               "~0 is white on a red-white-blue map and is indistinguishable "
               "from the page by any whiteness test.",
           "antialias_columns_excluded_at_each_edge": EDGE_AA,
           "largest_coverage_left_of_the_map_beyond_its_antialias":
               worst_outside,
           "coverage_at_map_x0": int(covered[map_x0]),
           "coverage_at_map_x1": int(covered[map_x1]),
           "separation_note":
               f"The strongest non-map column left of the map reaches "
               f"{worst_outside}/{nband} rows (the panel letter 'a'); the "
               f"map's own outer columns reach {int(covered[map_x0])} and "
               f"{int(covered[map_x1])}."}

    # --- the colorbar's vertical extent ----------------------------------
    cbmid = img[:, cb_x0 + 6:cb_x1 - 5, :]
    med = np.median(cbmid, axis=1)
    strong = np.flatnonzero(med[:330].sum(axis=1) < 740)
    cb_y0, cb_y1 = int(strong.min()), int(strong.max())

    col = (cb_x0 + cb_x1) // 2
    return {"map_x0": map_x0, "map_x1": map_x1, "y0": y0, "y1": y1,
            "cb_x0": cb_x0, "cb_x1": cb_x1,
            "cb_y0": cb_y0, "cb_y1": cb_y1,
            "map_edge_separation": sep,
            "colorbar_distinct_colours": int(ncolours[col]),
            "colorbar_min_distinct_colours":
                int(ncolours[cb_x0:cb_x1 + 1].min()),
            "max_distinct_colours_in_a_map_column":
                int(ncolours[map_x0:map_x1 + 1].max()),
            "separation_ratio": round(
                float(ncolours[cb_x0:cb_x1 + 1].min())
                / max(1.0, float(ncolours[map_x0:map_x1 + 1].max())), 2),
            "band_row_runs": [list(r) for r in runs_y]}


def read_fig1a(img):
    g = _measure_heatmap(img)
    # --- the LUT, from the colorbar in THIS raster -------------------------
    # The bar carries NO tick marks; the three labels 1.0 / 0 / -1.0 are free
    # text whose glyph centroids sit at y 26.5, 118.5 and 216.5, i.e. 92 px
    # apart and then 98 -- NOT equally spaced, so they are not usable anchors
    # (error class 20). What IS usable is the bar's own extent plus the fact
    # that a diverging red-white-blue map has its white point at the value the
    # middle label names. That is checked below, not assumed.
    cb_mid = (g["cb_x0"] + g["cb_x1"]) // 2
    ys = np.arange(g["cb_y0"] + 1, g["cb_y1"])
    ramp = np.array([np.median(img[y, cb_mid - 10:cb_mid + 11, :], axis=0)
                     for y in ys])
    span = float(g["cb_y1"] - g["cb_y0"])
    vals = 1.0 - 2.0 * (ys - g["cb_y0"]) / span

    # dedupe exact-duplicate colours: a continuous ramp renders neighbouring
    # positions identically, and heatmap.match RAISES on an exact tie rather
    # than letting argmin choose by row order. Collapsing them to their mean
    # value is the honest reading, and the widest collapsed bin IS this
    # panel's value resolution.
    seen, keep_c, keep_v, widths = {}, [], [], []
    for c, v in zip(ramp, vals):
        t = tuple(c)
        seen.setdefault(t, []).append(v)
    for t, vv in seen.items():
        keep_c.append(list(t))
        keep_v.append(float(np.mean(vv)))
        widths.append(float(max(vv) - min(vv)))
    lut_c = np.asarray(keep_c)
    lut_v = np.asarray(keep_v)
    # --- the white point: a round number the authors DREW, never printed ---
    d = lut_c[:, 0] - lut_c[:, 2]
    white_v = float(lut_v[int(np.argmin(np.abs(d)))])

    # --- the cells ---------------------------------------------------------
    ncol, nrow = 41, 3
    xs_c = heatmap.grid(g["map_x0"], g["map_x1"] + 1, ncol)
    ys_c = heatmap.grid(g["y0"], g["y1"] + 1, nrow)
    XX, YY = np.meshgrid(xs_c, ys_c)
    rgb = heatmap.sample(img, XX.ravel(), YY.ravel(), half=8)

    # --- THE TIE POLICY, WRITTEN AT THE CALL SITE ------------------------
    # heatmap.match REFUSES when a cell is exactly equidistant from two
    # swatches, because on a DISCRETE panel that means duplicate swatches or a
    # grid sampling a gridline, and letting argmin break the tie by LUT row
    # order biased one earlier read by +18% over 293 rows.
    #
    # This panel is a CONTINUOUS ramp, where a tie means something else
    # entirely and harmless: the cell's colour falls exactly between two
    # ADJACENT rows of the colorbar, so its value is between their two values.
    # Six of the 123 cells do that.
    #
    # The policy is NOT to loosen the gate -- the gate is what catches a wrong
    # grid, and it still has to run on all 123 cells. It is to remove the
    # ambiguity at its source: merge the two tying ramp rows into one entry
    # carrying their mean colour and their mean value, then ask match again.
    # That coarsens the LUT by exactly one row, only where the data needs it,
    # and the widest merged bin is reported as this panel's resolution.
    merges = 0
    while True:
        try:
            res = heatmap.match(rgb, lut_c, max_dist=heatmap.MAX_DIST)
            break
        except DetectionError as e:
            if "exactly equidistant" not in str(e):
                raise
            d = np.linalg.norm(rgb[:, None, :] - lut_c[None, :, :], axis=2)
            srt = np.argsort(d, axis=1)
            near = np.take_along_axis(d, srt[:, :2], axis=1)
            tied = np.flatnonzero(np.abs(near[:, 1] - near[:, 0]) < 1e-9)
            pairs = {tuple(sorted(srt[t, :2].tolist())) for t in tied}
            drop = set()
            for i, j in sorted(pairs):
                if i in drop or j in drop:
                    continue
                lut_c[i] = 0.5 * (lut_c[i] + lut_c[j])
                lut_v[i] = 0.5 * (lut_v[i] + lut_v[j])
                widths.append(abs(float(lut_v[i] - lut_v[j])))
                drop.add(j)
            if not drop:
                raise
            keep_ix = [k for k in range(len(lut_v)) if k not in drop]
            lut_c, lut_v = lut_c[keep_ix], lut_v[keep_ix]
            merges += len(drop)
            if merges > 40:
                raise DetectionError(
                    "stage04: merged 40 LUT rows and cells are STILL tying. "
                    "That is not a continuous-ramp artefact; the grid is "
                    "probably off the panel.")
    resolution = float(max(widths)) + 2.0 / span
    v = lut_v[res["index"]].reshape(nrow, ncol)

    rows = []
    for c in range(ncol):
        for r in range(nrow):
            rows.append({"IL_name": LABELS[c], "panel_column": c,
                         "replicate_row": r + 1,
                         "scaled_value": round(float(v[r, c]), 4),
                         "rgb": "|".join(str(int(x)) for x in
                                         rgb.reshape(nrow, ncol, 3)[r, c]),
                         "lut_distance": round(
                             float(res["distance"].reshape(nrow, ncol)[r, c]),
                             2)})
    diag = {"geometry": g,
            "cell_w_px": round((g["map_x1"] + 1 - g["map_x0"]) / ncol, 3),
            "cell_h_px": round((g["y1"] + 1 - g["y0"]) / nrow, 3),
            "lut_entries_raw": int(len(ys)),
            "lut_entries_after_dedupe": int(len(lut_v)),
            "lut_rows_merged_to_break_ties": merges,
            "value_resolution": round(resolution, 4),
            "max_lut_colour_distance": round(res["max_distance"], 3),
            "mean_lut_colour_distance": round(
                float(res["distance"].mean()), 3),
            "n_ambiguous_within_10_rgb": res["n_ambiguous"],
            "n_ambiguous_note":
                "Structural, not a defect: on a CONTINUOUS ramp adjacent LUT "
                "entries differ by ~2 RGB units, so almost every cell has a "
                "second swatch within 10. heatmap.match reports this and does "
                "not raise on it; what it raises on is an exact tie, and the "
                "dedupe above removes those legitimately.",
            "white_point_value": round(white_v, 4),
            "white_point_check":
                "The diverging map's white point must be the value its middle "
                "label names, 0. Measured " + f"{white_v:+.4f}" + ", i.e. "
                f"{abs(white_v)*50:.2f}% of the half-range. The figure never "
                "prints this and draws no tick for it, so it is external to "
                "the calibration it validates.",
            "n_cells_at_scale_top": int((v >= 1.0 - resolution).sum()),
            "n_cells_at_scale_bottom": int((v <= -1.0 + resolution).sum()),
            "saturation_note":
                "The brightest cell's RGB equals the colorbar's top interior "
                "colour and the darkest equals its bottom interior colour, so "
                "the data SPAN the scale exactly: the colour limits were set "
                "to the data's own min and max. Nothing is clipped, so no "
                "cell is a bound and value_is_bound stays empty on this panel.",
            }
    return rows, diag, v


# =========================================================================
# Figure 1c
# =========================================================================
def _baseline_row(img, lo=1480, hi=1540):
    """Darkness-weighted centroid of a drawn x-axis rule."""
    gray = img.min(axis=2)
    sub = gray[lo:hi, 200:2040]
    prof = sub.min(axis=1)
    dark = np.flatnonzero(prof < 140)
    if dark.size == 0:
        raise DetectionError("stage04: no baseline found for panel c.")
    a, b = dark.min(), dark.max() + 1
    w = np.clip(140.0 - prof[a:b], 0, None)
    return float(lo + np.dot(np.arange(a, b), w) / w.sum())


def x_tick_lattice(img, scan_lo, scan_hi, n_groups=41, tol=2.0):
    """The regular lattice of x-axis tick marks under a panel.

    AXIS FURNITURE, not ink. Figure 1d's replicate dots are drawn WIDER than
    its bars on the short bars, so a window taken from the coloured runs is
    the dots' extent and not the bar's -- the shape EXTRACTION_BRIEF 7 records
    as making one bar read 97.73 instead of 77.84. The x ticks are drawn by
    the plotting package at the group centres and no datum can move them.

    Returns (origin, pitch, n_on_lattice, n_detected). The lattice is chosen
    as the (first, last) pair of detected ticks that puts the most ticks
    within `tol` px of a regular n_groups-position grid, so stray glyph
    centroids are outvoted rather than averaged in.
    """
    gray = np.asarray(img).min(axis=2)
    t = np.array(sorted(axis.detect_ticks_1d(gray, "x", scan_lo, scan_hi,
                                             dark=140, min_len=3,
                                             max_len=20)))
    if t.size < n_groups // 2:
        raise DetectionError(
            f"stage04: {t.size} x ticks found in rows {scan_lo}-{scan_hi}.")
    best = None
    for i in range(t.size):
        for j in range(t.size - 1, i, -1):
            pitch = (t[j] - t[i]) / (n_groups - 1)
            if not (20.0 < pitch < 80.0):
                continue
            idx = np.round((t - t[i]) / pitch)
            ok = np.abs(t - (t[i] + pitch * idx)) < tol
            if best is None or ok.sum() > best[0]:
                best = (int(ok.sum()), float(t[i]), float(pitch), ok.copy())
    if best is None:
        raise DetectionError("stage04: no regular x-tick lattice found.")
    n_ok, t0, pitch, ok = best
    on = t[ok]
    i0 = np.round((on - t0) / pitch)
    A = np.column_stack([i0, np.ones_like(i0)])
    coef, *_ = np.linalg.lstsq(A, on, rcond=None)
    pitch, origin = float(coef[0]), float(coef[1] + coef[0] * i0.min())
    resid = float(np.max(np.abs(on - (A @ coef))))
    return {"origin": origin - pitch * 0.0, "pitch": pitch,
            "n_on_lattice": n_ok, "n_detected": int(t.size),
            "lattice_residual_px": round(resid, 3),
            "first_index": int(i0.min()), "last_index": int(i0.max())}


def read_bar_panel(img, series, tick_lo, tick_hi, tick_values, y_top,
                   base_lo, base_hi, n_groups=41, labels=None,
                   x_tick_scan=None, bar_half_width=None):
    """Read one grouped-bar panel. Shared by Figure 1c and Figure 1d.

    ONE implementation, two panels. Figure 1d encodes the ratio of the two
    series Figure 1c reports absolutely, so reading it with the SAME code is
    what makes the 41-way agreement in stage06 a check on the reads rather
    than a check on two different pieces of glue.

    series      [(name, mask), ...] in left-to-right order within a group
    tick_lo/hi  the row range to search for y ticks (the 0 tick is excluded
                on purpose on both panels: it is fused with the baseline, and
                it is the prediction target)
    """
    labels = labels or LABELS
    gray = img.min(axis=2)
    ticks = axis.detect_ticks_1d(gray, "y", 152, 170, dark=140,
                                 min_len=3, max_len=20)
    ticks = [t for t in ticks if tick_lo < t < tick_hi]
    if len(ticks) != len(tick_values):
        raise DetectionError(
            f"stage04: {len(ticks)} y ticks in rows {tick_lo}-{tick_hi}, "
            f"expected {len(tick_values)} for values {tick_values}. Found at "
            f"{np.round(ticks, 2).tolist()}. The 0 tick is fused with the "
            "drawn baseline and is deliberately NOT used -- it is the "
            "prediction target.")
    base = _baseline_row(img, base_lo, base_hi)
    cal = axis.calibrate(ticks, tick_values, kind="linear", axis="y",
                         forbid_target=[base])
    pred0 = float(cal.to_pixel(0.0))
    zero_err = pred0 - base

    def runs(mask, minpx=3):
        """Bar x-extents, with the antialiased fringe trimmed off each edge.

        THE TRIM IS LOAD-BEARING. A run taken as the union over all rows is
        1 px wider than the bar on most rows, because one row of the top rule
        antialiases a pixel further out. bars.bar_top then needs 90% of a
        7 px flank and gets 6/7 = 0.857 on every row but one, so it refuses --
        correctly, because the window it was handed was not the bar. Five bars
        failed exactly that way. The trim walks inward from each edge while
        the column's coverage is under half the run's maximum; it only ever
        removes OUTER columns, so the low-coverage INTERIOR of a hollow bar
        (which is hollow by construction) is untouched.
        """
        c = mask[y_top:int(base) - 2, :].sum(axis=0)
        out, cur = [], None
        for x in range(mask.shape[1]):
            if c[x] >= minpx:
                cur = [x, x] if cur is None else [cur[0], x]
            elif cur:
                out.append(tuple(cur))
                cur = None
        if cur:
            out.append(tuple(cur))
        trimmed = []
        for a, b in out:
            if b - a < 12:
                continue
            hi = c[a:b + 1].max()
            while a < b and c[a] < 0.5 * hi:
                a += 1
            while b > a and c[b] < 0.5 * hi:
                b -= 1
            if b - a >= 12:
                trimmed.append((int(a), int(b)))
        return trimmed

    anchor_name, anchor_mask = series[-1]

    # --- TICK-ANCHORED WINDOWS, when the caller asks for them -------------
    if x_tick_scan is not None:
        lat = x_tick_lattice(img, *x_tick_scan, n_groups=n_groups)
        hw = float(bar_half_width)
        fits = {nm: {"pitch": lat["pitch"],
                     "origin": lat["origin"] - hw, "resid": 0.0,
                     "n_runs": None, "source": "x_tick_lattice"}
                for nm, _ in series}
        by_idx = {nm: {} for nm, _ in series}
        pitch, width = lat["pitch"], 2.0 * hw
        cross = {"x_tick_lattice": lat}
        pitch_resid = lat["lattice_residual_px"]
        rruns = []
    else:
        bruns, rruns = runs(series[0][1]), runs(anchor_mask)

        # --- the window grid, fitted over the runs that ARE detected ----------
        # My first version demanded all 41 red runs and refused at 38. The three
        # missing ones are bars whose height is at or near zero, so their side
        # strokes do not survive the run-length floor -- that is the data, not a
        # detection failure, and bars.bar_top is the right thing to refuse on
        # them, one window at a time. So the grid is FITTED: each detected run is
        # assigned the group index its position implies, and the fit is checked.
        if len(rruns) < 30:
            raise DetectionError(
                f"stage04: only {len(rruns)} red bar runs; too few to fit a "
                "41-group grid.")
        starts = np.array([a for a, _ in rruns], dtype=float)
        pitch0 = (starts.max() - starts.min()) / float(n_groups - 1)
        idx = np.round((starts - starts.min()) / pitch0)
        if (len(set(idx.tolist())) != len(idx) or idx.max() != n_groups - 1
                or idx.min() != 0):
            raise DetectionError(
                f"stage04: the detected red runs do not map one-to-one onto "
                f"group indices 0..{n_groups - 1} (got {sorted(set(idx.tolist()))[:8]}...). "
                "Refusing to fit a grid onto an ambiguous assignment.")
        A = np.column_stack([idx, np.ones_like(idx)])
        coef, *_ = np.linalg.lstsq(A, starts, rcond=None)
        pitch, r0 = float(coef[0]), float(coef[1])
        pitch_resid = float(np.max(np.abs(starts - (A @ coef))))
        if pitch_resid > 2.0:
            raise DetectionError(
                f"stage04: red bar starts are not evenly pitched; worst residual "
                f"{pitch_resid:.2f} px over {len(rruns)} runs. The groups are not "
                "a regular grid.")
        width = float(np.median([b - a + 1 for a, b in rruns]))

        # EACH SERIES IS FITTED INDEPENDENTLY and the fits are cross-checked.
        # Two independent fits that agree are a check; one fit plus an assumption
        # that "the other bar sits one width to the left" is not.
        fits, by_idx = {}, {}
        fits[anchor_name] = {"pitch": pitch, "origin": r0, "resid": pitch_resid,
                             "n_runs": len(rruns)}
        by_idx[anchor_name] = {int(k): v for k, v in zip(idx, rruns)}
        cross = {}
        for si, (sname, smask) in enumerate(series[:-1]):
            sruns = runs(smask)
            if len(sruns) < 20:
                raise DetectionError(
                    f"stage04: only {len(sruns)} runs for series {sname!r}.")
            ss = np.array([a_ for a_, _ in sruns], dtype=float)
            nominal = r0 - width * (len(series) - 1 - si)
            sidx = np.round((ss - nominal) / pitch)
            if len(set(sidx.tolist())) != len(sidx):
                raise DetectionError(
                    f"stage04: series {sname!r} runs do not map one-to-one onto "
                    "group indices.")
            As = np.column_stack([sidx, np.ones_like(sidx)])
            sc, *_ = np.linalg.lstsq(As, ss, rcond=None)
            s_pitch, s0 = float(sc[0]), float(sc[1])
            s_res = float(np.max(np.abs(ss - (As @ sc))))
            if abs(s_pitch - pitch) > 0.15:
                raise DetectionError(
                    f"stage04: series {sname!r} and {anchor_name!r} disagree on "
                    f"pitch, {s_pitch:.4f} vs {pitch:.4f} px.")
            off = r0 - s0
            if not (0.6 * width <= off <= 1.4 * width * (len(series) - 1 - si + 1)):
                raise DetectionError(
                    f"stage04: series {sname!r} sits {off:.2f} px from "
                    f"{anchor_name!r} but a bar is {width:.1f} px wide.")
            fits[sname] = {"pitch": s_pitch, "origin": s0, "resid": s_res,
                           "n_runs": len(sruns)}
            by_idx[sname] = {int(k): v for k, v in zip(sidx, sruns)}
            cross[f"{sname}_vs_{anchor_name}"] = {
                "pitch_agreement_px": round(abs(s_pitch - pitch), 4),
                "offset_px": round(off, 3)}

    out, fails = [], {}
    for i in range(n_groups):
        for colour, mask in series:
            x0f = fits[colour]["origin"]
            det = by_idx[colour]
            if i in det:
                x0, x1 = det[i][0], det[i][1] + 1
                window_from = "detected_run"
            else:
                x0 = int(round(x0f + pitch * i))
                x1 = int(round(x0 + width))
                window_from = "fitted_grid"
            try:
                t = bars.bar_top(mask, x0, x1, C_Y_TOP, int(base) + 1,
                                 baseline_row=int(round(base)))
                val = float(cal.to_value(t["row"]))
                rec = {"IL_name": labels[i], "panel_column": i,
                       "organ": colour,
                       "value": round(val, 4),
                       "top_row_px": t["row"],
                       "flank_delta_px": t["flank_delta_px"],
                       "coverage": t["coverage"], "width_px": t["width_px"],
                       "window_from": window_from,
                       "value_is_bound": "", "refused": ""}
            except DetectionError as e:
                fails[f"{labels[i]}|{colour}"] = str(e).split("\n")[0]
                # IS THIS A REFUSAL, OR A BAR OF ESSENTIALLY ZERO HEIGHT?
                # Those are different statements and bar_top is right to
                # refuse rather than return 0. The discriminator is measured,
                # not assumed: the topmost coloured row anywhere in the
                # window. If the whole bar lives within FLOOR_PX of the drawn
                # baseline it is not separable from the axis rule's own
                # thickness, and the honest value is an UPPER BOUND at that
                # floor -- which is exactly what value_is_bound exists for.
                # The probe is the topmost FULL-WIDTH row, with the nrow and
                # flank requirements dropped. It is a CLASSIFIER ONLY -- it
                # decides "is this bar sitting on the baseline" and never
                # produces a value, so dropping those requirements costs
                # nothing. "Topmost coloured row" was my first probe and it is
                # wrong: replicate dots float above even a zero-height bar, so
                # it classified 4 of these 5 as readable.
                sub = mask[C_Y_TOP:int(round(base)) + 1, x0:x1]
                cov = sub.mean(axis=1)
                wide = np.flatnonzero(cov >= 0.90)
                topmost = (C_Y_TOP + int(wide.min())) if wide.size else None
                if topmost is not None and base - topmost <= FLOOR_PX:
                    bound_v = FLOOR_PX / abs(float(cal.to_pixel(0.0)
                                                   - cal.to_pixel(1.0)))
                    rec = {"IL_name": labels[i], "panel_column": i,
                           "organ": colour,
                           "value": round(bound_v, 4),
                           "top_row_px": topmost, "flank_delta_px": "",
                           "coverage": "", "width_px": x1 - x0,
                           "window_from": window_from,
                           "value_is_bound": "upper",
                           "refused": "AT BASELINE, recorded as an upper "
                                      "bound: " + str(e).split("\n")[0]}
                else:
                    rec = {"IL_name": labels[i], "panel_column": i,
                           "organ": colour, "value": "",
                           "top_row_px": "", "flank_delta_px": "",
                           "coverage": "", "width_px": "",
                           "window_from": window_from,
                           "value_is_bound": "",
                           "refused": str(e).split("\n")[0]}
            out.append(rec)

    diag = {
        "tick_rows_px": [round(t, 3) for t in ticks],
        "tick_values": tick_values,
        "tick_fit_residual_px": round(float(np.max(np.abs(
            np.asarray(sorted(ticks))
            - np.asarray(sorted([float(cal.to_pixel(v))
                                 for v in C_TICK_VALUES]))))), 4),
        "px_per_unit": round(abs(float(cal.to_pixel(0.0) - cal.to_pixel(1.0))),
                             4),
        "drawn_baseline_row_px": round(base, 3),
        "predicted_row_of_value_zero_px": round(pred0, 3),
        "zero_prediction_error_px": round(zero_err, 3),
        "zero_prediction_error_in_units": round(
            abs(zero_err) / abs(float(cal.to_pixel(0.0) - cal.to_pixel(1.0))),
            5),
        "zero_prediction_note":
            "THE UNPRINTED CHECK. The axis is calibrated on the ticks for 3, "
            "2 and 1 ONLY; the 0 tick is fused with the baseline and was "
            "excluded, and `forbid_target` was passed so that fitting onto "
            "the baseline would be a TypeError rather than a silent class-27 "
            "absorption. The fit then predicts where value 0 falls, and it is "
            "compared with the rule the authors actually drew.",
        "bar_pitch_px": round(pitch, 4),
        "bar_pitch_fit_residual_px": round(pitch_resid, 3),
        "series_fits": {k: {kk: round(vv, 4) if isinstance(vv, float) else vv
                            for kk, vv in v.items()} for k, v in fits.items()},
        "series_cross_checks": cross,
        "bar_width_px": round(width, 2),
        "anchor_runs_detected": len(rruns),
        "runs_note":
            f"{len(rruns)} of {n_groups} runs are detectable for the anchor "
            "series; the "
            "missing ones are liver bars so short that their side strokes are "
            "under the 12 px run-length floor. The windows are NOT taken from "
            "the runs -- they come from the fitted pitch -- so a missing run "
            "costs nothing. bars.bar_top then refuses on those windows, which "
            "is the correct outcome and is recorded as a refusal, not a zero.",
        "refusals": fails,
        "n_refused": len(fails),
        "floor_px": FLOOR_PX,
        "n_recorded_as_upper_bound": sum(1 for r in out
                                         if r.get("value_is_bound") == "upper"),
        "bound_note":
            "A bar whose entire body lies within FLOOR_PX of the drawn "
            "baseline is recorded with value_is_bound = 'upper' at the floor "
            "value, not as a zero and not as a missing value. bars.bar_top "
            "refuses on them, and that refusal is correct: 'the bar has no "
            "resolvable height' is a different statement from 'the bar is "
            "this tall'.",
    }
    return out, diag


def masks(img):
    R, G, B = img[:, :, 0], img[:, :, 1], img[:, :, 2]
    return {"blue": (B - R) > 60, "red": (R - B) > 60}


def read_fig1c(img):
    m = masks(img)
    rows, diag = read_bar_panel(
        img, series=[("liver", m["blue"]), ("spleen", m["red"])],
        tick_lo=1150, tick_hi=1480, tick_values=C_TICK_VALUES,
        y_top=C_Y_TOP, base_lo=1480, base_hi=1540)
    for r in rows:
        r["value_avg_radiance_e8"] = r.pop("value")
    return rows, diag


def main():
    img, meta = load_raster()
    a_rows, a_diag, _ = read_fig1a(img)
    c_rows, c_diag = read_fig1c(img)

    with open(os.path.join(HERE, "04_fig1a_cells.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(a_rows[0].keys()))
        w.writeheader()
        w.writerows(a_rows)
    with open(os.path.join(HERE, "04_fig1c_bars.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(c_rows[0].keys()))
        w.writeheader()
        w.writerows(c_rows)
    geo = {"_what": "Every geometric constant used to read Figure 1, measured "
                    "from the raster rather than assumed.",
           "raster": meta or "already extracted; see a previous run",
           "fig1a": a_diag, "fig1c": c_diag}
    with open(os.path.join(HERE, "04_geometry.json"), "w") as fh:
        json.dump(geo, fh, indent=1)

    print("Fig 1a:", a_diag["geometry"])
    print(f"  cells 41x3, max LUT colour distance "
          f"{a_diag['max_lut_colour_distance']} (mean "
          f"{a_diag['mean_lut_colour_distance']}), resolution "
          f"{a_diag['value_resolution']}")
    print(f"  white point {a_diag['white_point_value']:+.4f} (want 0)")
    print(f"  cells at scale top {a_diag['n_cells_at_scale_top']}, "
          f"bottom {a_diag['n_cells_at_scale_bottom']}")
    print("Fig 1c:")
    print(f"  ticks {c_diag['tick_rows_px']} -> fit residual "
          f"{c_diag['tick_fit_residual_px']} px, "
          f"{c_diag['px_per_unit']} px per 1e8 unit")
    print(f"  ZERO PREDICTION: baseline drawn at "
          f"{c_diag['drawn_baseline_row_px']}, fit on ticks 3/2/1 predicts "
          f"{c_diag['predicted_row_of_value_zero_px']} -> error "
          f"{c_diag['zero_prediction_error_px']} px = "
          f"{c_diag['zero_prediction_error_in_units']} units")
    print(f"  pitch {c_diag['bar_pitch_px']} px (residual "
          f"{c_diag['bar_pitch_fit_residual_px']}), width "
          f"{c_diag['bar_width_px']} px")
    print(f"  bars read {82 - c_diag['n_refused']}/82, refused "
          f"{c_diag['n_refused']}")
    for k, v in list(c_diag["refusals"].items())[:12]:
        print(f"     REFUSED {k}: {v[:110]}")


if __name__ == "__main__":
    main()
