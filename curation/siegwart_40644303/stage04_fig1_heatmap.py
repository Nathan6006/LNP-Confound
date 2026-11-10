#!/usr/bin/env python3
"""Figure 1a: the 465 x 2 discrete delivery heatmap.

    python3 -I stage04_fig1_heatmap.py

WHAT THE PANEL IS. 15 amine cores (rows 3A1..6A3) x 31 thiol tails (columns
SC1..SC18, SN1..SN8, SO1..SO8). 15 x 31 = 465, which is the number the abstract
claims. Each cell is split into two stripes: the UPPER stripe is siRNA knockdown
(red, 3 bands) and the LOWER stripe is mRNA luminescence (blue, 3 bands).

THE GRID COMES FROM THE DRAWN CELLS, NOT FROM THE PANEL BOX. The cells are
separated by white gutters in both directions, so the stripe extents are
measurable directly. F-019 is what anchoring on a panel border costs. The
detected column centres are then CHECKED against the 31 column-label text
centroids, and the residual is recorded rather than described.

THE SWATCHES ARE CHECKED AGAINST PANEL c, NOT NAMED IN A CAPTION.
Panel c prints six legend rectangles ("High Delivery" / "Intermediate
Delivery" x three hues). Two of them are the siRNA high and intermediate bands
and two are the mRNA high and intermediate bands, sampled from THIS SAME
raster. The third (lowest) band of each cargo has no legend rectangle.

The swatch actually used for matching is the MEDIAN OF THE CELLS IN ITS OWN
BAND, found by a 3-median clustering seeded from the two legend rectangles plus
the leftovers' median. That is deliberate and it is a measurement choice, not a
relaxation: the figure is a lossy JPEG with chroma subsampling, so a flat fill
is reproduced slightly differently in panel c than in panel a, and the right
nearest-neighbour problem is "which band is this cell in", asked in the space
the cells are measured in. The legend is then used as what it is -- an EXTERNAL
ANCHOR -- and each derived swatch is asserted to lie within 5.0 RGB of its
printed rectangle. Measured: siRNA high 2.45, siRNA intermediate 2.00, mRNA
high 1.00, mRNA intermediate 0.50. The band ASSIGNMENT is identical either way
and identical at both patch sizes tried (counts 35/153/277 and 68/72/325), so
nothing about the result rests on this choice; only the distance diagnostic
does. Matching the cells straight against the legend rectangles is also
recorded, and its max distance is 5.00 for both cargoes -- which passes, but
with no headroom, and a diagnostic with no headroom cannot tell a geometry
error from a compression artefact.

The derived lowest band additionally has to survive three tests: it must be a
single tight cluster, it must be lighter than the intermediate swatch in the
same hue, and it must sit near the pale legend TEXT colour printed beside panel
a. All three are asserted, not assumed.

PATCH SIZE. The stripes are 7-8 px tall and the columns 36 px wide, so half=3
(a 7x7 median) is the largest square that stays inside a stripe. At half=2 the
max nearest-swatch distance is 7.3-7.9 and at half=3 it is 4.2-4.6, entirely
from JPEG chroma noise: the assignments do not change. The threshold was not
touched -- the measurement was made better.

WHY A DISCRETE PANEL IS THE BEST CASE FOR A PIXEL READ. There is no colour ramp
to invert, so there is no scale to get wrong -- only a grid. CLAUDE.md section
5: a discrete legend means nearest-neighbour matching against swatches sampled
from the same raster, which is exact rather than estimated, and anything above
~5 colour distance means the geometry is wrong, not the colours.
heatmap.match() raises there; it is not a number this file prints and ignores.
"""
import json
import os
import sys

import numpy as np
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))

import axis      # noqa: E402
import heatmap   # noqa: E402

FIG = os.path.join(HERE, "source", "figures", "fig1.jpg")
OUT = os.path.join(HERE, "04_fig1_heatmap.json")

# --- the library's own order, read off the figure's labels -----------------
CORES = ["3A1", "3A2", "3A3", "3A4", "3A5",
         "4A1", "4A3", "4A4",
         "5A1", "5A2", "5A3", "5A4",
         "6A1", "6A2", "6A3"]
TAILS = ([f"SC{i}" for i in range(1, 13)] + ["SC14", "SC16", "SC18"]
         + [f"SN{i}" for i in range(1, 9)]
         + [f"SO{i}" for i in range(1, 9)])

# Panel a's interior, generous: the run detectors find the real extents.
PANEL_Y0, PANEL_Y1 = 131, 516
PANEL_X0, PANEL_X1 = 160, 1377
# The column-label band, above the panel. Used only to CHECK the grid.
LABEL_Y0, LABEL_Y1 = 40, 127
LABEL_MIN_W = 15          # "# of Arms" is 64 px; a bracket tick is 4 px

# Panel c legend rectangles. Located by colour, inside these two bands.
LEGEND_BANDS = {"high": (840, 872), "intermediate": (885, 920)}
LEGEND_X0, LEGEND_X1 = 1000, 1500
# The pale legend TEXT beside panel a, used as the third-band cross-check.
TEXT_LEGEND = {"siRNA_low": (260, 272), "mRNA_low": (365, 378)}
TEXT_X0, TEXT_X1 = 1385, 1500

# The bands the legend names, in the order [high, intermediate, low].
SIRNA_BANDS = [(70.0, 100.0), (25.0, 70.0), (0.0, 25.0)]      # % knockdown
MRNA_BANDS = [(2.0, None), (1.0, 2.0), (0.0, 1.0)]            # log10 luminescence


def _runs(mask):
    out, s = [], None
    for i, v in enumerate(mask):
        if v and s is None:
            s = i
        elif not v and s is not None:
            out.append((s, i))
            s = None
    if s is not None:
        out.append((s, len(mask)))
    return out


def _centres(runs, offset=0):
    return [offset + 0.5 * (a + b) - 0.5 for a, b in runs]


def detect_grid(arr):
    """Column and stripe centres, from the drawn cells' own extents."""
    R, B = arr[:, :, 0], arr[:, :, 2]
    red = (R - B) > 12
    blue = (B - R) > 12

    ysl = slice(PANEL_Y0, PANEL_Y1)
    xsl = slice(PANEL_X0, PANEL_X1)

    red_rows = _runs(red[ysl, xsl].mean(axis=1) > 0.5)
    blue_rows = _runs(blue[ysl, xsl].mean(axis=1) > 0.5)
    cols = _runs((red | blue)[ysl, :].mean(axis=0) > 0.5)

    if len(red_rows) != len(CORES):
        raise SystemExit(
            f"stage04: {len(red_rows)} red stripes detected, expected "
            f"{len(CORES)}. REFUSING. The row grid is the thing that assigns a "
            "value to a lipid; a miscount here silently relabels every row "
            "below it, which is the off-by-one header row CLAUDE.md section 5 "
            "says the colour-distance diagnostic has caught twice.")
    if len(blue_rows) != len(CORES):
        raise SystemExit(
            f"stage04: {len(blue_rows)} blue stripes detected, expected "
            f"{len(CORES)}. REFUSING, same reason.")
    if len(cols) != len(TAILS):
        raise SystemExit(
            f"stage04: {len(cols)} columns detected, expected {len(TAILS)}. "
            "REFUSING.")

    return {
        "red_y": _centres(red_rows, PANEL_Y0),
        "blue_y": _centres(blue_rows, PANEL_Y0),
        "x": _centres(cols, 0),
        "col_runs": cols,
        "red_runs": [(a + PANEL_Y0, b + PANEL_Y0) for a, b in red_rows],
        "blue_runs": [(a + PANEL_Y0, b + PANEL_Y0) for a, b in blue_rows],
    }


def label_centres(arr):
    """x centroids of the 31 rotated column labels. A CHECK, never the grid."""
    g = arr.mean(axis=2)
    band = g[LABEL_Y0:LABEL_Y1, :] < 128
    runs = [r for r in _runs(band.sum(axis=0) > 0)
            if LABEL_MIN_W <= r[1] - r[0] <= 30]
    return _centres(runs, 0)


def legend_swatches(arr):
    """The six panel-c legend rectangles, by colour, from this same raster."""
    out = {}
    for level, (y0, y1) in LEGEND_BANDS.items():
        row = arr[y0:y1, LEGEND_X0:LEGEND_X1]
        R, G, B = row[:, :, 0], row[:, :, 1], row[:, :, 2]
        col = (np.abs(R - B) > 25) | (np.abs(R - G) > 25)
        runs = [r for r in _runs(col.sum(axis=0) > 5) if r[1] - r[0] > 10]
        if len(runs) != 3:
            raise SystemExit(
                f"stage04: panel c '{level}' band has {len(runs)} legend "
                "rectangles, expected 3 (red, orange, blue). REFUSING rather "
                "than matching cells against a swatch set we cannot locate.")
        vals = []
        for s, e in runs:
            patch = row[5:-5, s + 4:e - 4].reshape(-1, 3)
            vals.append(np.median(patch, axis=0))
        out[level] = vals        # [red, orange, blue] in drawn order
    return out


def text_legend_colour(arr, y0, y1):
    """Median colour of the pale legend TEXT beside panel a (a cross-check)."""
    reg = arr[y0:y1, TEXT_X0:TEXT_X1]
    R, G, B = reg[:, :, 0], reg[:, :, 1], reg[:, :, 2]
    m = (np.abs(R - B) > 25) | (np.abs(R - G) > 25)
    if m.sum() < 30:
        raise SystemExit("stage04: legend text band carries no coloured ink.")
    return np.median(reg[m], axis=0)


PATCH_HALF = 3            # 7x7 median; the stripes are 7-8 px tall
LEGEND_TOL_RGB = 5.0      # derived swatch vs its printed legend rectangle


def cluster_swatches(cells, hi, mid, lo_seed, label):
    """3-median clustering seeded from the two legend rectangles + leftovers.

    Returns (centres, counts). Refuses if a cluster empties -- an empty band
    means the seeds do not describe this panel.
    """
    cen = np.vstack([hi, mid, lo_seed]).astype(float)
    idx = None
    for _ in range(50):
        d = np.linalg.norm(cells[:, None, :] - cen[None, :, :], axis=2)
        new = d.argmin(axis=1)
        if idx is not None and np.array_equal(new, idx):
            break
        idx = new
        for k in range(3):
            if not (idx == k).any():
                raise SystemExit(
                    f"stage04/{label}: cluster {k} emptied during the median "
                    "refinement. The three legend bands do not describe this "
                    "panel. REFUSING.")
            cen[k] = np.median(cells[idx == k], axis=0)
    return cen, [int((idx == k).sum()) for k in range(3)]


def third_band(cells, hi, mid, label):
    """Derive the lowest band's swatch from the cells it is the only home of.

    Panel c prints rectangles for 'High' and 'Intermediate' only. Rather than
    invent the third colour, take every cell that is closer to NEITHER printed
    swatch than to the cluster of leftovers, and require the leftovers to be a
    single tight cluster. A second cluster here would mean the panel has four
    bands, not three, and that is a refusal, not a rounding choice.
    """
    d_hi = np.linalg.norm(cells - hi, axis=1)
    d_mid = np.linalg.norm(cells - mid, axis=1)
    rest = cells[(d_hi > 40) & (d_mid > 40)]
    if rest.shape[0] < 10:
        raise SystemExit(
            f"stage04/{label}: only {rest.shape[0]} cells are far from both "
            "printed swatches, which is too few to define the third band.")
    centre = np.median(rest, axis=0)
    spread = float(np.max(np.linalg.norm(rest - centre, axis=1)))
    if spread > 12.0:
        raise SystemExit(
            f"stage04/{label}: the leftover cells span {spread:.1f} in RGB "
            "around their own median. A discrete band is a flat fill; this "
            "much spread means there is more than one band in the leftovers, "
            "so the panel does not have three levels and this script's band "
            "list is wrong. REFUSING.")
    return centre, spread, int(rest.shape[0])


def main():
    arr = np.asarray(Image.open(FIG).convert("RGB")).astype(float)
    grid = detect_grid(arr)

    # --- the grid is uniform, and that is asserted on measured floats ----
    cal_x = axis.calibrate(grid["x"], list(range(len(TAILS))), kind="linear",
                           residual_tol_px=2.5)
    cal_red = axis.calibrate(grid["red_y"], list(range(len(CORES))),
                             kind="linear", residual_tol_px=2.5)
    cal_blue = axis.calibrate(grid["blue_y"], list(range(len(CORES))),
                              kind="linear", residual_tol_px=2.5)

    # --- the independent check on the column grid ------------------------
    labs = label_centres(arr)
    if len(labs) != len(TAILS):
        raise SystemExit(
            f"stage04: {len(labs)} column labels found, expected "
            f"{len(TAILS)}. REFUSING: without them the column grid has no "
            "independent check.")
    dx = np.abs(np.asarray(labs) - np.asarray(grid["x"]))
    if dx.max() > 6.0:
        raise SystemExit(
            f"stage04: a column label centroid is {dx.max():.1f} px from its "
            "cell centre. The labels and the cells disagree about where the "
            "columns are, which is an off-by-one in one of them. REFUSING.")

    # --- swatches ---------------------------------------------------------
    leg = legend_swatches(arr)
    si_hi, _orange_hi, mr_hi = leg["high"]
    si_mid, _orange_mid, mr_mid = leg["intermediate"]

    xs = np.asarray(grid["x"])
    cells_si = heatmap.sample(
        arr, np.tile(xs, len(CORES)),
        np.repeat(np.asarray(grid["red_y"]), len(TAILS)), half=PATCH_HALF)
    cells_mr = heatmap.sample(
        arr, np.tile(xs, len(CORES)),
        np.repeat(np.asarray(grid["blue_y"]), len(TAILS)), half=PATCH_HALF)

    si_lo, si_spread, si_n = third_band(cells_si, si_hi, si_mid, "siRNA")
    mr_lo, mr_spread, mr_n = third_band(cells_mr, mr_hi, mr_mid, "mRNA")

    # the pale legend TEXT is anti-aliased, so it is lighter than the fill;
    # it is a direction check, not an equality check.
    si_txt = text_legend_colour(arr, *TEXT_LEGEND["siRNA_low"])
    mr_txt = text_legend_colour(arr, *TEXT_LEGEND["mRNA_low"])
    for name, lo, mid, txt in (("siRNA", si_lo, si_mid, si_txt),
                               ("mRNA", mr_lo, mr_mid, mr_txt)):
        if not float(np.sum(lo)) > float(np.sum(mid)):
            raise SystemExit(
                f"stage04/{name}: the derived lowest band is not lighter than "
                "the printed intermediate swatch. The band order is wrong.")
        if float(np.linalg.norm(lo - txt)) > 45.0:
            raise SystemExit(
                f"stage04/{name}: derived low swatch {lo.round(1).tolist()} is "
                f"{np.linalg.norm(lo - txt):.1f} from the pale legend text "
                f"colour {txt.round(1).tolist()}. REFUSING.")

    # --- the swatches used for matching: the cells' own band medians ------
    si_sw, si_counts = cluster_swatches(cells_si, si_hi, si_mid, si_lo, "siRNA")
    mr_sw, mr_counts = cluster_swatches(cells_mr, mr_hi, mr_mid, mr_lo, "mRNA")

    # --- and the external anchor they must agree with ---------------------
    legend_check = {}
    for cargo, sw, printed in (("siRNA", si_sw, [si_hi, si_mid]),
                               ("mRNA", mr_sw, [mr_hi, mr_mid])):
        for lvl, (derived, rect) in enumerate(zip(sw[:2], printed)):
            dist = float(np.linalg.norm(np.asarray(derived) - np.asarray(rect)))
            legend_check[f"{cargo}_{'high' if lvl == 0 else 'intermediate'}"] \
                = round(dist, 3)
            if dist > LEGEND_TOL_RGB:
                raise SystemExit(
                    f"stage04/{cargo}: the band median derived from the cells "
                    f"{np.asarray(derived).round(1).tolist()} is {dist:.2f} "
                    f"RGB from the panel-c legend rectangle "
                    f"{np.asarray(rect).round(1).tolist()}, over the "
                    f"{LEGEND_TOL_RGB} tolerance. The legend is the only "
                    "EXTERNAL statement of what these colours mean; if the "
                    "cells do not agree with it, the band identification is "
                    "not established and no value may be emitted. REFUSING.")

    m_si = heatmap.match(cells_si, si_sw)
    m_mr = heatmap.match(cells_mr, mr_sw)

    # Recorded, not acted on: matching the cells straight against the printed
    # rectangles. It passes at exactly 5.00 for both cargoes, i.e. with no
    # headroom at all, which is why the derived medians are what `match` is
    # given. A diagnostic sitting on its own threshold cannot distinguish a
    # geometry error from a compression artefact.
    legend_direct = {}
    for cargo, cells, sw in (("siRNA", cells_si,
                              np.vstack([si_hi, si_mid, si_sw[2]])),
                             ("mRNA", cells_mr,
                              np.vstack([mr_hi, mr_mid, mr_sw[2]]))):
        d = np.linalg.norm(cells[:, None, :] - sw[None, :, :], axis=2)
        legend_direct[cargo] = {
            "max_distance": round(float(d.min(axis=1).max()), 3),
            "assignment_identical_to_derived":
                bool(np.array_equal(d.argmin(axis=1),
                                    (m_si if cargo == "siRNA" else m_mr)
                                    ["index"])),
        }

    # --- bands -> ordinal values -----------------------------------------
    # The mRNA high band is OPEN above (">log2"), so invert() returns its
    # finite edge and marks it `lower`. The lowest band of each cargo is the
    # colormap floor: every value inside it renders the same colour, so it is
    # an upper bound on the truth, flagged by COLOUR IDENTITY (floor_index),
    # never by comparing a value to a threshold -- that asymmetry is F-040.
    # NO floor_index. The lowest band of each cargo is CLOSED -- the legend
    # prints "0-25%" and "log 0-1", not "<25%" and "<log1" -- so it is a band,
    # not a colormap floor, and flagging it as `value_is_bound = upper` would
    # repeat F-040, where 65 rows that were measurements were marked bounds.
    # Only the mRNA top band is open ("> log2"), and invert marks that "lower"
    # from the band list itself. Whether the authors' lowest bin is really
    # closed at 0 is queued as D-5.
    inv_si = heatmap.invert(m_si["index"], SIRNA_BANDS)
    inv_mr = heatmap.invert(m_mr["index"], MRNA_BANDS)

    rows = []
    for ri, core in enumerate(CORES):
        for ci, tail in enumerate(TAILS):
            k = ri * len(TAILS) + ci
            rows.append({
                "lipid": f"{core}-{tail}",
                "core": core, "tail": tail, "row": ri, "col": ci,
                "x_px": round(float(xs[ci]), 2),
                "siRNA_y_px": round(float(grid["red_y"][ri]), 2),
                "mRNA_y_px": round(float(grid["blue_y"][ri]), 2),
                "siRNA_band": int(m_si["index"][k]),
                "siRNA_ordinal": 3 - int(m_si["index"][k]),
                "siRNA_dist": round(float(m_si["distance"][k]), 2),
                "mRNA_band": int(m_mr["index"][k]),
                "mRNA_ordinal": 3 - int(m_mr["index"][k]),
                "mRNA_dist": round(float(m_mr["distance"][k]), 2),
                "mRNA_lo": inv_mr[k]["lo"], "mRNA_hi": inv_mr[k]["hi"],
                "mRNA_bound": inv_mr[k]["bound"],
                "siRNA_lo": inv_si[k]["lo"], "siRNA_hi": inv_si[k]["hi"],
                "siRNA_bound": inv_si[k]["bound"],
            })

    out = {
        "_what": "Figure 1a read cell by cell. ordinal 3 = highest band.",
        "figure": os.path.relpath(FIG, HERE),
        "grid": {
            "n_cores": len(CORES), "n_tails": len(TAILS),
            "n_cells": len(rows),
            "x_calibration": cal_x.qc(),
            "siRNA_stripe_calibration": cal_red.qc(),
            "mRNA_stripe_calibration": cal_blue.qc(),
            "label_vs_cell_centre_max_px": round(float(dx.max()), 2),
            "label_vs_cell_centre_mean_px": round(float(dx.mean()), 2),
            "stripe_height_px": [int(b - a) for a, b in grid["red_runs"][:3]],
            "column_width_px": [int(b - a) for a, b in grid["col_runs"][:3]],
        },
        "swatches": {
            "patch_half_px": PATCH_HALF,
            "siRNA_used_for_matching": si_sw.round(1).tolist(),
            "mRNA_used_for_matching": mr_sw.round(1).tolist(),
            "cluster_counts_siRNA": si_counts,
            "cluster_counts_mRNA": mr_counts,
            "derived_vs_panel_c_legend_rgb": legend_check,
            "derived_vs_panel_c_tolerance": LEGEND_TOL_RGB,
            "matching_straight_against_the_printed_rectangles": legend_direct,
            "siRNA_high_from_panel_c": si_hi.round(1).tolist(),
            "siRNA_intermediate_from_panel_c": si_mid.round(1).tolist(),
            "siRNA_low_seed_derived": si_lo.round(1).tolist(),
            "siRNA_low_cluster_spread_rgb": round(si_spread, 2),
            "siRNA_low_n_cells_used": si_n,
            "siRNA_low_legend_text_colour": si_txt.round(1).tolist(),
            "mRNA_high_from_panel_c": mr_hi.round(1).tolist(),
            "mRNA_intermediate_from_panel_c": mr_mid.round(1).tolist(),
            "mRNA_low_seed_derived": mr_lo.round(1).tolist(),
            "mRNA_low_cluster_spread_rgb": round(mr_spread, 2),
            "mRNA_low_n_cells_used": mr_n,
            "mRNA_low_legend_text_colour": mr_txt.round(1).tolist(),
        },
        "match": {
            "siRNA_max_distance": round(m_si["max_distance"], 3),
            "siRNA_n_ambiguous_within_10": m_si["n_ambiguous"],
            "mRNA_max_distance": round(m_mr["max_distance"], 3),
            "mRNA_n_ambiguous_within_10": m_mr["n_ambiguous"],
            "max_distance_threshold": heatmap.MAX_DIST,
        },
        "band_counts": {
            "siRNA": {str(b): int((m_si["index"] == b).sum())
                      for b in (0, 1, 2)},
            "mRNA": {str(b): int((m_mr["index"] == b).sum())
                     for b in (0, 1, 2)},
        },
        "cells": rows,
    }
    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps({k: v for k, v in out.items() if k != "cells"}, indent=1))
    print(f"wrote {len(rows)} cells to {OUT}")


if __name__ == "__main__":
    main()
