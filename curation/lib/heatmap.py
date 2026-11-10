#!/usr/bin/env python3
"""ONE audited heatmap reader: cell grid, swatch matching, LUT inversion.

WHY THIS FILE EXISTS. Three separate findings, one class each:

  F-019  QL_2026's colorbar was anchored on the OUTER EDGE of its own black
         border. Every Fig 2b value came out high -- ~1% at the hit threshold
         and +38% at the floor -- and the 48 white cells were stored at 1,310
         RLU when their whole feasible band is [0, ~982], i.e. ABOVE the
         `value_is_bound = upper` flag the same review round had just given
         them. The geometry of a colorbar is the geometry of an axis, so the
         calibration lives in axis.py and this module only consumes it.

  F-032  band midpoints. A discrete legend (`<1e4 / 1e4-1e5 / ...`) gives a
         cell a BAND, not a value. Storing an edge biases every cell one way;
         the midpoint is the only choice that does not, and the bounds have to
         travel with it so a consumer can tell a band from a measurement.

  F-040  65 QL_2026 rows flagged `value_is_bound = upper` are measurements. A
         `< 4.0` colour threshold straddled the floor band and a band at 3.32
         that is measurably not white (distance 1.73), and those 65 cells are
         bounded BELOW as well as above -- so the flag was both wrong and
         asymmetric. Hence: `invert()` reports floor and ceiling the same way,
         by COLOUR IDENTITY to the extreme swatch, never by a value threshold.

THE DIAGNOSTIC THAT HAS EARNED ITSELF TWICE. `match()` returns the maximum
colour distance over all cells and the count of cells within `ambiguous_within`
of a second swatch. CLAUDE.md section 5: "anything above ~5 means the geometry
is wrong, not the colours." It has caught an off-by-one header row twice, so
`match()` RAISES above the threshold rather than reporting it.
"""
import numpy as np

from errors import CalibrationError, DetectionError

MAX_DIST = 5.0          # max permitted nearest-swatch distance, RGB euclidean
AMBIGUOUS_WITHIN = 10.0  # a cell this close to a 2nd swatch is reported


def grid(lo, hi, n, kind="centres"):
    """n evenly spaced cell centres (or n+1 edges) between two axis anchors.

    `lo` and `hi` are the OUTER edges of the first and last cell, measured from
    the axis labels -- not from the panel border. F-019 is what reading an outer
    border costs. Returns a float array; the caller rounds.
    """
    if n < 1:
        raise CalibrationError(f"grid: n={n}.")
    if hi == lo:
        raise CalibrationError("grid: lo == hi, zero-extent axis.")
    edges = np.linspace(float(lo), float(hi), n + 1)
    if kind == "edges":
        return edges
    if kind != "centres":
        raise CalibrationError(f"grid: kind={kind!r}, want centres|edges.")
    return 0.5 * (edges[:-1] + edges[1:])


def sample(img, xs, ys, half=2):
    """Median RGB over a (2*half+1)^2 patch at each (x, y). Shape (n, 3).

    Median, not mean: a mean over a patch straddling a gridline returns a
    colour that is in the LUT nowhere, and then `match` reports a large
    distance for a correct grid. Median survives one dark row.
    """
    a = np.asarray(img)
    if a.ndim != 3 or a.shape[2] < 3:
        raise DetectionError(f"sample: expected an RGB image, got {a.shape}.")
    out = []
    for x, y in zip(np.asarray(xs).ravel(), np.asarray(ys).ravel()):
        xi, yi = int(round(float(x))), int(round(float(y)))
        y0, y1 = max(0, yi - half), min(a.shape[0], yi + half + 1)
        x0, x1 = max(0, xi - half), min(a.shape[1], xi + half + 1)
        if y1 <= y0 or x1 <= x0:
            raise DetectionError(
                f"sample: patch at ({xi},{yi}) is outside the {a.shape[1]}x"
                f"{a.shape[0]} image. The grid is off the panel.")
        out.append(np.median(a[y0:y1, x0:x1, :3].reshape(-1, 3), axis=0))
    return np.asarray(out, dtype=float)


def match(cells_rgb, swatch_rgb, max_dist=MAX_DIST,
          ambiguous_within=AMBIGUOUS_WITHIN):
    """Nearest-swatch index per cell, with the geometry diagnostic enforced.

    `swatch_rgb` must be sampled FROM THE SAME RASTER as the cells -- a legend
    swatch read off the same rendered page, not a colour named in the caption.
    Returns dict(index, distance, max_distance, n_ambiguous, ambiguous_rows).
    Raises DetectionError above `max_dist`.
    """
    C = np.asarray(cells_rgb, dtype=float)
    S = np.asarray(swatch_rgb, dtype=float)
    if C.ndim != 2 or C.shape[1] != 3:
        raise DetectionError(f"match: cells shape {C.shape}, want (n,3).")
    if S.ndim != 2 or S.shape[1] != 3 or S.shape[0] < 2:
        raise DetectionError(f"match: swatches shape {S.shape}, want (k>=2,3).")

    d = np.linalg.norm(C[:, None, :] - S[None, :, :], axis=2)   # (n, k)
    idx = d.argmin(axis=1)
    near = d.min(axis=1)

    # A tie must not be resolved by argmin's "first smallest wins". That is the
    # open QL E-1 shape: d.argmin() returning the topmost tying LUT row biased
    # a read by up to +18% over 293 rows.
    srt = np.sort(d, axis=1)
    second = srt[:, 1]
    tied = np.flatnonzero(np.abs(second - near) < 1e-9)
    if tied.size:
        raise DetectionError(
            f"match: {tied.size} cell(s) are exactly equidistant from two "
            f"swatches (rows {tied[:10].tolist()}). argmin would resolve the "
            "tie by LUT row order, which is an arbitrary choice presented as a "
            "measurement. Either the two swatches are duplicates -- dedupe the "
            "LUT -- or the grid is sampling a gridline.")

    amb = np.flatnonzero(second - near < ambiguous_within)
    mx = float(near.max())
    if mx > max_dist:
        worst = int(near.argmax())
        raise DetectionError(
            f"match: max nearest-swatch distance {mx:.2f} > {max_dist} "
            f"(cell {worst}, rgb {C[worst].round(1).tolist()}, nearest swatch "
            f"{S[idx[worst]].round(1).tolist()}).\n"
            "CLAUDE.md section 5: anything above ~5 means THE GEOMETRY IS "
            "WRONG, NOT THE COLOURS. This diagnostic has caught an off-by-one "
            "header row twice. Check the grid's outer edges against the axis "
            "LABELS, not the panel border (F-019), and confirm the swatches "
            "were sampled from this same raster.")
    return {"index": idx, "distance": near, "max_distance": mx,
            "n_ambiguous": int(amb.size),
            "ambiguous_rows": amb[:50].tolist()}


def invert(index, bands, floor_index=None, ceiling_index=None):
    """Discrete band index -> (value, lo, hi, bound) per cell.

    `bands` is a list of (lo, hi) in data units, one per swatch, in swatch
    order. Use None for an open end: (None, 1e4) for "<1e4", (1e7, None) for
    ">1e7".

    The returned value is the BAND MIDPOINT (geometric for a positive band,
    arithmetic otherwise) -- F-032. An open band has no midpoint, so its value
    is its finite edge and `bound` says which way it is open.

    `bound` is "upper", "lower" or "" and is decided by WHICH SWATCH the cell
    matched, never by comparing a value to a threshold. F-040: a `< 4.0`
    threshold straddled the floor band and a band at 3.32 that is measurably
    not the floor colour, mis-flagging 65 measurements as bounds. Colour
    identity is the only thing that actually says "this cell is on the floor".

    `floor_index` / `ceiling_index` let the caller name the extreme swatches
    explicitly when the band list does not make them obvious. They are treated
    SYMMETRICALLY: F-040's 65 cells were bounded below as well as above and
    only one direction was recorded.
    """
    idx = np.asarray(index).ravel().astype(int)
    if len(bands) == 0:
        raise DetectionError("invert: empty band list.")
    if idx.size and (idx.max() >= len(bands) or idx.min() < 0):
        raise DetectionError(
            f"invert: swatch index range [{idx.min()},{idx.max()}] does not fit "
            f"{len(bands)} bands. The LUT and the band list disagree.")

    out = []
    for i in idx:
        lo, hi = bands[i]
        if lo is None and hi is None:
            raise DetectionError(f"invert: band {i} is open at both ends.")
        if lo is None:
            val, bound = float(hi), "upper"
        elif hi is None:
            val, bound = float(lo), "lower"
        elif lo > 0 and hi > 0:
            val, bound = float(np.sqrt(lo * hi)), ""
        else:
            val, bound = 0.5 * (float(lo) + float(hi)), ""
        if floor_index is not None and i == floor_index:
            bound = "upper"
        if ceiling_index is not None and i == ceiling_index:
            bound = "lower"
        out.append({"value": val, "lo": lo, "hi": hi, "bound": bound,
                    "swatch": int(i)})
    return out

def cross_panel_collisions(panel_a, panel_b, max_dist=None):
    """Label pairs whose swatches are closer than MAX_DIST ACROSS two panels.

    WHY. `MAX_DIST` is calibrated for matching cells against swatches sampled
    from THE SAME raster, where a legend swatch and its cells differ only by
    JPEG noise. It is not a statement about two different figures, and across
    panels it becomes actively dangerous.

    Measured 2026-11-10 (`peerd_31999380`): Figure S3 and Figure S4 give
    different colours AND different legend orders to the same two lipids.

        S3 "Lipid 13" = (127, 4, 127)
        S4 "Lipid 12" = (128, 0, 129)
        euclidean distance 4.583  <  MAX_DIST 5.0

    So a cross-panel colour match would have confidently assigned S4's Lipid 12
    values to Lipid 13 and **never tripped the guard that exists to catch
    exactly this**. The curation agent avoided it by matching only within each
    panel; this function makes the hazard visible instead of relying on that.

    Pass `{label: (r, g, b)}` for each panel. Returns the DIFFERENT-label pairs
    that collide -- each one a place where a cross-panel match is unsafe.
    Same-label pairs are not collisions: those are the figures agreeing.
    """
    lim = MAX_DIST if max_dist is None else max_dist
    out = []
    for la, ca in panel_a.items():
        for lb, cb in panel_b.items():
            d = sum((x - y) ** 2 for x, y in zip(ca, cb)) ** 0.5
            if d < lim and str(la).strip() != str(lb).strip():
                out.append({"a": la, "a_rgb": tuple(ca), "b": lb,
                            "b_rgb": tuple(cb), "distance": round(d, 3)})
    out.sort(key=lambda r: r["distance"])
    return out
