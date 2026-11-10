#!/usr/bin/env python3
"""Stage 04 -- the three heatmaps, read through curation/lib.

    python3 -I stage04_values.py

WHAT IS BEING READ.

  Fig S2c   I.D. Renilla-luciferase relative expression, LNP type 1, 36 cells
  Fig S2e   I.D. Renilla-luciferase relative expression, LNP type 4, 36 cells
  Fig 4b    I.M. serum human-EPO relative expression,     LNP type 4, 36 cells

All three are fold change versus the Nil group (Fig S2 caption: "The results of
all experiments were compared by substituting fold change based on the values
of the Nil group included in the experiment in question"). The units are FOLD
CHANGE and nothing is converted.

A CONTINUOUS COLOURBAR IS STILL A LOOK-UP TABLE, AND ITS RESOLUTION IS FINITE.

The ramp is sampled row by row from the bar itself -- the same raster as the
cells, never a colour named in a caption -- and then COLLAPSED ON EXACT RGB
EQUALITY. Two bar rows with identical 8-bit RGB are one swatch covering the
value interval they span, and that interval is the band carried into
`heatmap.invert`. So a cell's value is the midpoint of the band of values that
render in its colour, with `lo`/`hi` travelling beside it (F-032), rather than
a single number pretending to a precision the colour ramp does not have.

Collapsing on equality is also what makes `heatmap.match` usable: it RAISES on
an exact LUT tie, and an un-deduped ramp is nothing but ties in its flat
regions.

THE FLOOR AND THE CEILING ARE DECIDED BY COLOUR IDENTITY, SYMMETRICALLY. A cell
matching the bar's lightest swatch gets `value_is_bound = upper`; a cell
matching the darkest gets `lower`. Never by comparing a value to a threshold --
F-040 is 65 rows mis-flagged by a `< 4.0` threshold that straddled the floor
band and a band that was measurably not the floor colour.

GEOMETRY. The scale and origin of every colourbar come from its own four
detected tick notches via `axis.calibrate`; there is no way to supply an origin
separately (F-033). The cell grid is taken from `heatmap.grid` between the
panel's two outer border lines, and then CHECKED against the internal
gridlines detected independently in the raster, so the grid is measured twice.

THE ASTERISKS ARE ON THE CELLS. Significance markers are drawn in the middle of
many cells, so a centre sample would read a glyph. Each cell's colour is taken
as the modal RGB of its interior, `heatmap.sample` is called at an offset point
chosen to miss the glyph, and the two are required to be equal. A cell whose
modal colour covers less than MIN_MODE_FRACTION of its interior is refused
rather than averaged.
"""
import csv
import json
import os
import sys
from collections import Counter

import numpy as np
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import axis                                                 # noqa: E402
import heatmap                                              # noqa: E402
from errors import DetectionError                           # noqa: E402

MIN_MODE_FRACTION = 0.55   # of the cell interior that must sit on the flat fill
FLAT_TOL = 3.0             # RGB euclidean; the PDF figure streams carry +/-1 dither
SAMPLE_TOL = 2.0           # max |interior median - offset sample|, RGB euclidean

# Row order is the same in all three panels: T-1 over EB-5/6/7, then T-2.
ROWS = [("T-1", "EB-5"), ("T-1", "EB-6"), ("T-1", "EB-7"),
        ("T-2", "EB-5"), ("T-2", "EB-6"), ("T-2", "EB-7")]
COLS = ["FH-1", "FH-2", "FH-3", "FH-4", "FH-5", "FH-6"]

PANELS = [
    {"panel": "Fig S2c", "image": "source/si_media/figS2_ID_RLuc.tif",
     "search": {"x": (780, 1180), "y": (560, 950)},
     "bar_x": (1170, 1230), "tick_values": [4000, 3000, 2000, 1000],
     "readout": "Renilla luciferase relative expression (I.D.)",
     "unit": "fold change vs Nil group",
     "lnp_type": 1, "route": "intradermal"},
    {"panel": "Fig S2e", "image": "source/si_media/figS2_ID_RLuc.tif",
     "search": {"x": (776, 1176), "y": (1100, 1490)},
     "bar_x": (1165, 1225), "tick_values": [4000, 3000, 2000, 1000],
     "readout": "Renilla luciferase relative expression (I.D.)",
     "unit": "fold change vs Nil group",
     "lnp_type": 4, "route": "intradermal"},
    {"panel": "Fig 4b", "image": "source/pdfimg/fig4_part1_1947x2228.png",
     "search": {"x": (160, 740), "y": (490, 1070)},
     "bar_x": (735, 800), "tick_values": [8, 6, 4, 2],
     "readout": "human EPO relative expression (I.M., serum, 6 h)",
     "unit": "fold change vs Nil group",
     "lnp_type": 4, "route": "intramuscular"},
]


# --------------------------------------------------------------------------
def runs(idx, gap=2):
    """Group consecutive indices; return each group's centre."""
    out, cur = [], [idx[0]]
    for v in idx[1:]:
        if v - cur[-1] <= gap:
            cur.append(v)
        else:
            out.append(cur)
            cur = [v]
    out.append(cur)
    return [float(np.mean(c)) for c in out]


def find_box(dark, search):
    x0, x1 = search["x"]
    y0, y1 = search["y"]
    sub = dark[y0:y1, x0:x1]
    w, h = sub.shape[1], sub.shape[0]
    rows = np.flatnonzero(sub.sum(axis=1) > 0.75 * w)
    cols = np.flatnonzero(sub.sum(axis=0) > 0.75 * h)
    if len(rows) < 2 or len(cols) < 2:
        raise DetectionError("find_box: no full-width/height border lines.")
    rc, cc = runs(rows), runs(cols)
    return (cc[0] + x0, cc[-1] + x0, rc[0] + y0, rc[-1] + y0)


def internal_gridlines(img, left, right, top, bottom, n, axis_name):
    """Cell edges detected from the drawn gridlines, independent of the box.

    A boundary shows itself in one of two ways and BOTH are needed:

      dip   a thin grey rule, darker than the quietest pixel a short way off
            on both sides. This is the only signal where two neighbouring
            cells carry the same fill.
      step  the two plateaus a short way either side are different colours.
            This is the only signal where the rule is LIGHTER than at least
            one of its neighbours -- measured in Fig 4b, where the FH-5/FH-6
            rule is (161,160,161) between fills of (130,132,192) and
            (117,120,185), i.e. brighter than both, so the dip test alone
            found 3 of the 5 internal boundaries.

    Two earlier detectors were wrong in instructive ways and the comments stay
    because the failures are the evidence that this one works:
      - "darker than the strip median" marked a whole dark COLUMN (FH-5 in
        Fig S2c) and returned its centre as a gridline -- a detector miss
        dressed as a detection, error class 28;
      - a MEAN along the perpendicular axis let the significance glyphs make
        six phantom boundaries at the cell centres of Fig S2e. The median
        along that axis sees a rule, which spans the panel, and not a glyph,
        which covers a tenth of it.
    """
    a = img.astype(float)
    # Sample only the leading 22% of each cell along the perpendicular axis.
    # The significance glyphs are centred in their cells, so this strip is
    # glyph-free by construction. Taking the whole cell instead let white
    # asterisks drag the column median far enough to fire the step test four
    # times INSIDE the FH-5 and FH-6 cells of Fig 4b.
    cw, ch = (right - left) / n, (bottom - top) / n
    if axis_name == "x":
        rows = np.concatenate([np.arange(int(top + k * ch) + 4,
                                         int(top + k * ch + 0.22 * ch))
                               for k in range(n)])
        band = a[rows, int(left):int(right) + 1, :3]
        rgb = np.median(band, axis=0)
    else:
        cols = np.concatenate([np.arange(int(left + k * cw) + 4,
                                         int(left + k * cw + 0.22 * cw))
                               for k in range(n)])
        band = a[int(top):int(bottom) + 1, :, :3][:, cols, :]
        rgb = np.median(band, axis=1)
    strip = rgb.mean(axis=1)
    look, margin, step_tol = 9, 8.0, 10.0
    dips = []
    # The outer borders are excluded by construction: both tests need `look`
    # px of fill on BOTH sides, which the panel's first and last boundary do
    # not have.
    edge = look + 5          # clear of the panel's own border rule
    for i in range(edge, len(strip) - edge):
        lo, hi = strip[i - look:i - 2], strip[i + 3:i + look + 1]
        is_dip = (lo.size and hi.size and strip[i] < lo.min() - margin
                  and strip[i] < hi.min() - margin)
        pl = np.median(rgb[i - look:i - 2], axis=0)
        pr = np.median(rgb[i + 3:i + look + 1], axis=0)
        is_step = float(np.linalg.norm(pl - pr)) > step_tol
        if is_dip or is_step:
            dips.append(i)
    dips = np.array(dips, dtype=int)
    if dips.size == 0:
        raise DetectionError(f"internal_gridlines({axis_name}): none found.")
    span = (right - left) if axis_name == "x" else (bottom - top)
    off = left if axis_name == "x" else top
    # Cluster with a tolerance that scales with the cell, so the THICK rule
    # between the T-1 and T-2 blocks -- three dark rows plus a colour step on
    # each side -- comes back as one boundary and not two.
    cell = span / n
    centres = [c + off for c in runs(dips, gap=max(3, int(0.12 * cell)))]
    if len(centres) != n - 1:
        raise DetectionError(
            f"internal_gridlines({axis_name}): {len(centres)} lines detected "
            f"at {np.round(centres, 1).tolist()}, expected {n - 1} internal "
            f"cell boundaries. span={span:.1f}")
    return centres


def bar_ticks(img, bar_left, bar_right):
    """Tick notches drawn INSIDE the colourbar, as departures from the ramp.

    The ramp is horizontally uniform, so a tick is a row where the pixels near
    one edge differ from that row's own median colour. Detected on both edges
    and required to agree.
    """
    a = img.astype(float)
    w = bar_right - bar_left
    inner = (int(bar_left + 0.30 * w), int(bar_left + 0.70 * w))
    sides = {"left": (int(bar_left + 0.07 * w), int(bar_left + 0.25 * w)),
             "right": (int(bar_left + 0.75 * w), int(bar_left + 0.93 * w))}
    found = {}
    inset = 5          # skip the bar's own border rows at top and bottom
    for name, (sx0, sx1) in sides.items():
        dev = []
        for y in range(inset, a.shape[0] - inset):
            ramp = np.median(a[y, inner[0]:inner[1], :3], axis=0)
            patch = a[y, sx0:sx1, :3]
            dev.append(float(np.linalg.norm(patch - ramp, axis=1).max()))
        dev = np.array(dev)
        hit = np.flatnonzero(dev > 60)
        found[name] = [v + inset for v in runs(hit, gap=3)] if hit.size else []
    return found


def cell_colour(img, cx, cy, half_w, half_h):
    """(interior median RGB, flat fraction, heatmap.sample at an offset point).

    MEDIAN, NOT MODE. The modal 8-bit triple looked like the right statistic
    and is not: the article PDF's figure streams carry +/-1 dither, so Fig 4b's
    FH-3/T-1/EB-5 fill is spread over (202,200,228), (203,199,228) and
    (202,199,228) and the strict mode covered only 0.46 of the interior. A
    mode that splits is a flatness test that fails on a flat cell.

    `flat` is the fraction of interior pixels within FLAT_TOL of that median,
    which is what "the fill is uniform" actually means. The significance
    glyphs are white and are exactly the pixels this fraction excludes.

    The offset point sits left of and below centre, where no significance
    glyph is drawn in any of the three panels, and `heatmap.sample`'s own
    patch median must agree with the interior median.
    """
    a = np.asarray(img)
    y0, y1 = int(cy - half_h), int(cy + half_h) + 1
    x0, x1 = int(cx - half_w), int(cx + half_w) + 1
    patch = a[y0:y1, x0:x1, :3].reshape(-1, 3).astype(float)
    med = np.median(patch, axis=0)
    flat = float(np.mean(np.linalg.norm(patch - med, axis=1) <= FLAT_TOL))
    off = heatmap.sample(img, [cx - 0.62 * half_w], [cy + 0.60 * half_h],
                         half=3)[0]
    return med, flat, off


def build_lut(img, bar_left, bar_right, top, bottom, cal):
    """Distinct ramp colours, each with the value band that renders in it."""
    a = np.asarray(img).astype(int)
    w = bar_right - bar_left
    x0, x1 = int(bar_left + 0.32 * w), int(bar_left + 0.68 * w)
    seen = {}
    for y in range(int(top) + 3, int(bottom) - 2):
        rgb = tuple(int(v) for v in np.median(a[y, x0:x1, :3], axis=0))
        v = float(cal.to_value(y))
        seen.setdefault(rgb, []).append(v)
    items = sorted(((rgb, min(vs), max(vs)) for rgb, vs in seen.items()),
                   key=lambda t: t[1])
    swatches = [list(rgb) for rgb, _, _ in items]
    bands = [(lo, hi) for _, lo, hi in items]
    return swatches, bands


# --------------------------------------------------------------------------
def resolve_ties(cells, swatches, bands, max_rounds=8):
    """Merge LUT entries that a cell sits exactly between, and say so.

    `heatmap.match` RAISES on an exact tie, because argmin would otherwise
    resolve it by LUT row order -- an arbitrary choice presented as a
    measurement. That refusal is correct and it fires here for a real reason:
    the article PDF's figure streams carry +/-1 dither, so a cell's interior
    median can land exactly halfway between two neighbouring ramp colours.

    The honest reading of that is not "pick one". It is that the cell's value
    lies somewhere in the UNION of the two colours' bands. So the two LUT
    entries are merged into one, whose band spans both and whose midpoint is
    the value. Every merge is recorded and reported in the QC file; the cost
    is a wider band on the affected cells, which is what `lo`/`hi` are for.
    """
    merged = []
    for _ in range(max_rounds):
        S = np.asarray(swatches, dtype=float)
        d = np.linalg.norm(np.asarray(cells, float)[:, None, :]
                           - S[None, :, :], axis=2)
        srt = np.sort(d, axis=1)
        tied = np.flatnonzero(np.abs(srt[:, 1] - srt[:, 0]) < 1e-9)
        if tied.size == 0:
            return swatches, bands, merged
        i = int(tied[0])
        a, b = sorted(np.argsort(d[i])[:2].tolist())
        merged.append({"cell_index": i,
                       "cell_rgb": [round(float(v), 1) for v in cells[i]],
                       "swatch_a": swatches[a], "swatch_b": swatches[b],
                       "band_a": [round(bands[a][0], 2), round(bands[a][1], 2)],
                       "band_b": [round(bands[b][0], 2), round(bands[b][1], 2)]})
        newrgb = [0.5 * (swatches[a][k] + swatches[b][k]) for k in range(3)]
        newband = (min(bands[a][0], bands[b][0]), max(bands[a][1], bands[b][1]))
        swatches = ([s for j, s in enumerate(swatches) if j not in (a, b)]
                    + [newrgb])
        bands = ([bd for j, bd in enumerate(bands) if j not in (a, b)]
                 + [newband])
        order = sorted(range(len(bands)), key=lambda j: bands[j][0])
        swatches = [swatches[j] for j in order]
        bands = [bands[j] for j in order]
    raise DetectionError(
        "resolve_ties: still tied after %d merges. The grid is probably "
        "sampling a gridline rather than a fill." % max_rounds)


def read_panel(spec):
    path = os.path.join(HERE, spec["image"])
    img = Image.open(path).convert("RGB")
    arr = np.asarray(img).astype(int)
    dark = arr.sum(axis=2) < 200

    left, right, top, bottom = find_box(dark, spec["search"])
    xe = internal_gridlines(arr, left, right, top, bottom, 6, "x")
    ye = internal_gridlines(arr, left, right, top, bottom, 6, "y")
    gx = heatmap.grid(left, right, 6, kind="edges")
    gy = heatmap.grid(top, bottom, 6, kind="edges")
    dx = float(np.max(np.abs(np.array(xe) - gx[1:-1])))
    dy = float(np.max(np.abs(np.array(ye) - gy[1:-1])))
    cellw = (right - left) / 6.0
    cellh = (bottom - top) / 6.0
    if dx > 0.12 * cellw or dy > 0.12 * cellh:
        raise DetectionError(
            f"{spec['panel']}: the detected gridlines and the grid implied by "
            f"the outer borders disagree by {dx:.1f}/{dy:.1f} px "
            f"(cell {cellw:.1f}x{cellh:.1f}). Two independent reads of the "
            "same geometry must agree before any colour is sampled.")

    # --- colourbar ------------------------------------------------------
    bl, br = spec["bar_x"]
    # Restricted to THIS panel's own rows. Over the whole image height the
    # window catches the other panel's colourbar too, and bcent[0]/bcent[-1]
    # then straddle two different bars -- a geometry error that would have
    # produced a plausible calibration from the wrong pixels.
    bsub = dark[int(top):int(bottom) + 1, bl:br + 1]
    bcols = np.flatnonzero(bsub.sum(axis=0) > 0.6 * (bottom - top))
    if bcols.size < 2:
        raise DetectionError(f"{spec['panel']}: colourbar borders not found.")
    bcent = runs(bcols + bl)
    bar_left, bar_right = bcent[0], bcent[-1]
    tk = bar_ticks(arr[int(top):int(bottom) + 1, :, :], bar_left, bar_right)
    tk = {k: [v + top for v in vs] for k, vs in tk.items()}
    if len(tk["left"]) != 4 or len(tk["right"]) != 4:
        raise DetectionError(
            f"{spec['panel']}: {len(tk['left'])} left / {len(tk['right'])} "
            "right colourbar tick notches; expected 4 of each.")
    if max(abs(a - b) for a, b in zip(tk["left"], tk["right"])) > 1.5:
        raise DetectionError(
            f"{spec['panel']}: the colourbar's left and right tick notches "
            f"disagree: {tk['left']} vs {tk['right']}.")
    ticks = [0.5 * (a + b) for a, b in zip(tk["left"], tk["right"])]
    cal = axis.calibrate(ticks, spec["tick_values"], kind="linear",
                         residual_tol_px=1.0, axis="y")

    # THE ROUND NUMBER THE AUTHORS DREW AND NEVER PRINTED: the bar's own
    # bottom edge. If the calibration is right it lands on exactly 0.
    zero_px = float(cal.to_pixel(0.0))
    zero_err_px = zero_px - bottom
    top_value = float(cal.to_value(top))

    swatches, bands = build_lut(img, bar_left, bar_right, top, bottom, cal)

    # --- cells ----------------------------------------------------------
    cx = heatmap.grid(left, right, 6)
    cy = heatmap.grid(top, bottom, 6)
    modes, offs, fracs, cells = [], [], [], []
    for ri, (t, eb) in enumerate(ROWS):
        for ci, fh in enumerate(COLS):
            m, f, o = cell_colour(img, cx[ci], cy[ri],
                                  cellw * 0.36, cellh * 0.36)
            if f < MIN_MODE_FRACTION:
                raise DetectionError(
                    f"{spec['panel']} {fh}/{t}/{eb}: only {f:.2f} of the cell "
                    "interior sits within FLAT_TOL of its median. The fill is not flat "
                    "-- the grid is probably off by a row or a column.")
            if float(np.linalg.norm(m - o)) > SAMPLE_TOL:
                raise DetectionError(
                    f"{spec['panel']} {fh}/{t}/{eb}: the interior median {m} and "
                    f"heatmap.sample at the offset point {o} disagree. One of "
                    "them is sitting on a significance glyph.")
            modes.append(m)
            offs.append(o)
            fracs.append(f)
            cells.append((fh, t, eb))

    swatches, bands, merged = resolve_ties(np.array(modes), swatches, bands)
    res = heatmap.match(np.array(modes), np.array(swatches),
                        max_dist=heatmap.MAX_DIST)
    inv = heatmap.invert(res["index"], bands,
                         floor_index=0, ceiling_index=len(bands) - 1)

    # The floor and ceiling swatches are OPEN BANDS. The LUT is sampled a few
    # rows inside the bar's own border, so values between 0 and the lowest
    # sampled row render in a colour the LUT does not contain and fall onto
    # the floor swatch. Their true value is therefore anywhere in
    # [0, band_hi] -- which is exactly what `value_is_bound = upper` asserts,
    # and leaving band_lo at the swatch's own lo would contradict the flag.
    # Fig 4c draws this curation's floor cell at 12.1 against a stored 53.1,
    # inside the widened band and outside the narrow one.
    top_v = float(cal.to_value(top))
    for v in inv:
        if v["bound"] == "upper":
            v["lo"] = 0.0
        elif v["bound"] == "lower":
            v["hi"] = top_v

    out = []
    for (fh, t, eb), v, d in zip(cells, inv, res["distance"]):
        out.append({"panel": spec["panel"], "IL_name": f"{fh}-{eb}{t}",
                    "FH": fh, "EB": eb, "T": t,
                    "lnp_type": spec["lnp_type"], "route": spec["route"],
                    "readout": spec["readout"], "unit": spec["unit"],
                    "value": round(v["value"], 2),
                    "band_lo": round(v["lo"], 2), "band_hi": round(v["hi"], 2),
                    "value_is_bound": v["bound"],
                    "swatch_distance": round(float(d), 3)})
    qc = {"panel": spec["panel"], "image": spec["image"],
          "box_px": {"left": left, "right": right, "top": top,
                     "bottom": bottom},
          "cell_px": {"w": round(cellw, 3), "h": round(cellh, 3)},
          "gridlines_detected_x": [round(v, 1) for v in xe],
          "gridlines_detected_y": [round(v, 1) for v in ye],
          "grid_vs_border_max_disagreement_px": {"x": round(dx, 2),
                                                 "y": round(dy, 2)},
          "colourbar_px": {"left": bar_left, "right": bar_right},
          "colourbar_ticks_px": [round(v, 2) for v in ticks],
          "colourbar_tick_values": spec["tick_values"],
          "calibration": cal.qc(),
          "DRAWN_ZERO": {
              "_what": "The colourbar's own bottom edge is a round number the "
                       "authors DREW and never printed. The calibration is "
                       "fitted ONLY on the four tick notches, so where it puts "
                       "value 0 is a prediction, not a fit.",
              "predicted_px_for_value_0": round(zero_px, 2),
              "drawn_bar_bottom_px": round(bottom, 2),
              "error_px": round(zero_err_px, 2),
              "error_in_value_units": round(abs(zero_err_px)
                                            / cal.px_per_unit, 3)},
          "colourbar_top_value": round(top_value, 1),
          "lut_distinct_swatches": len(swatches),
          "lut_tie_merges": merged,
          "lut_band_width_units": {
              "median": round(float(np.median([h - l for l, h in bands])), 3),
              "max": round(float(max(h - l for l, h in bands)), 3)},
          "match_max_distance": round(res["max_distance"], 3),
          "match_n_ambiguous": res["n_ambiguous"],
          "cell_flat_fraction_min": round(float(min(fracs)), 3),
          "n_on_floor_swatch": sum(1 for r in out
                                   if r["value_is_bound"] == "upper"),
          "n_on_ceiling_swatch": sum(1 for r in out
                                     if r["value_is_bound"] == "lower")}
    return out, qc


def main():
    rows, qcs = [], []
    for spec in PANELS:
        r, q = read_panel(spec)
        rows.extend(r)
        qcs.append(q)
        print(f"{q['panel']}: {len(r)} cells, max swatch distance "
              f"{q['match_max_distance']}, LUT {q['lut_distinct_swatches']} "
              f"swatches, drawn-zero error {q['DRAWN_ZERO']['error_px']} px "
              f"({q['DRAWN_ZERO']['error_in_value_units']} units), "
              f"floor {q['n_on_floor_swatch']} ceiling "
              f"{q['n_on_ceiling_swatch']}")
    with open(os.path.join(HERE, "05_values.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    with open(os.path.join(HERE, "04_values_qc.json"), "w") as fh:
        json.dump({"panels": qcs}, fh, indent=1)
    print(f"05_values.csv: {len(rows)} values")


if __name__ == "__main__":
    main()
