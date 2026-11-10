#!/usr/bin/env python3
"""ONE audited chain-length counter: vector first, raster second, agree or raise.

WHY THIS FILE EXISTS. Counting a zig-zag by eye is this project's documented
+/-1 failure mode (error class 48, formerly 19). It cost the iPhos apex count,
and it is the mechanism behind the `OF-2` tail-length error -- all four tails
read as C17 when the reagent is C18, while `LL_2025` held the same reagent
correctly as `OF-02` (F-031, and see reagents.cross_check).

Two counters exist in the repo and each is used by exactly one paper. This
module promotes both and makes them check each other:

  VECTOR   curation/ng_2025/stage01_vector.py reads the PDF's own bond
           coordinates. 20 of 20 structures exact. Strictly better than any
           vision-based read, because the drawing's geometry IS the answer --
           there is no estimation step at all.
  RASTER   curation/xd_2026_lipidai/stage02_trace.py skeletonises, traces, and
           counts Douglas-Peucker vertices. The fallback for a flattened
           figure, which is most of these PDFs.

THE CONTRACT. `count_carbons` requires TWO INDEPENDENT COUNTS THAT AGREE.
Given a vector count and a raster count it returns the number only if they
match; given one it says which, and the caller must record that the second was
unavailable. It does not average them, prefer one, or return the "more
confident" one -- a chain is an integer, and two methods differing by one means
one of them is wrong, not that the answer is 17.5.

CLAUDE.md section 5, and this is the trap that is NOT a geometry problem:
"Name tokens are not chain lengths." KS_2024's `C24` means two C4 tails,
KZ_2016's `SC1/SC2/SC3` are isobutyl/isopentyl/neopentyl, YX_2024's `C1/C2/C3`
are component indices. `check_name_token` makes the comparison explicit instead
of letting a token stand in for a measurement; audit.py C2 does the corpus-wide
version over 634 series.
"""
import math

from errors import DetectionError

VECTOR_TOL_PT = 1.2    # endpoints closer than this are one vertex (ng_2025)


def vector_bonds(page, bbox, min_width=0.01, include_filled=False):
    """Stroked line segments inside `bbox` of a pymupdf page. [(x0,y0,x1,y1)].

    A STROKED path is a bond. A FILLED path with no width is part of a letter
    outline -- heteroatom labels are drawn as glyph outlines, not text -- so it
    is label evidence, never a bond. Promoted from ng_2025/stage01_vector.py.

    `include_filled` EXISTS BECAUSE THAT RULE IS NOT UNIVERSAL. Measured in
    andersondg_41203968 (Nat Nanotechnol 2025 SI, 2026-11-10): **all 13,271
    drawing paths are FILLED with `width = None`**, so the stroked-only filter
    returned zero segments on every page and this function was useless there.
    Some ChemDraw export paths emit every bond as a filled polygon.

    Turning it on is not free: with filled paths included, glyph outlines come
    back as "bonds" too, so the caller MUST separate them -- by line thickness
    (that SI is cleanly bimodal, 12,723 paths at 0.3-0.6 pt and 402 at 1.0-2.3
    with nothing between), by aspect ratio, or by proximity to a label box.
    `vector_chain_length` will happily count glyph strokes as vertices if you
    do not. Report which mode you used in the folder's QC file.
    """
    X0, Y0, X1, Y1 = bbox
    out = []
    for p in page.get_drawings():
        if not include_filled and (not p.get("width")
                                   or p["width"] < min_width):
            continue
        for it in p["items"]:
            if it[0] != "l":
                continue
            a, b = it[1], it[2]
            if (X0 <= a.x <= X1 and X0 <= b.x <= X1
                    and Y0 <= a.y <= Y1 and Y0 <= b.y <= Y1):
                out.append((a.x, a.y, b.x, b.y))
    return out


def _merge_vertices(bonds, tol):
    pts = []
    for s in bonds:
        pts.extend([(s[0], s[1]), (s[2], s[3])])
    reps = []
    for p in pts:
        for i, r in enumerate(reps):
            if math.hypot(p[0] - r[0], p[1] - r[1]) <= tol:
                break
        else:
            reps.append(p)
    return reps


def vector_chain_length(bonds, tol=VECTOR_TOL_PT, drop_double=True,
                        double_offset_pt=2.2):
    """Vertices in a chain of stroked bonds, after merging coincident endpoints.

    `drop_double` removes the inner line of a double bond -- a near-parallel
    segment within `double_offset_pt` of another -- which would otherwise add a
    phantom vertex. Returns the vertex count; the caller maps vertices to
    carbons for its own drawing convention (a terminal methyl is a vertex).
    """
    segs = list(bonds)
    if not segs:
        raise DetectionError(
            "vector_chain_length: no stroked segments in the bbox. Either the "
            "figure is raster -- use the raster path and say so -- or the bbox "
            "is wrong. Do not fall back to counting by eye without recording "
            "it: that is error class 48, the +/-1 failure mode.")
    if drop_double:
        keep, dropped = [], set()
        for i, a in enumerate(segs):
            for j in range(i + 1, len(segs)):
                if j in dropped:
                    continue
                b = segs[j]
                la = math.hypot(a[2] - a[0], a[3] - a[1])
                lb = math.hypot(b[2] - b[0], b[3] - b[1])
                if la == 0 or lb == 0:
                    continue
                cross = abs((a[2] - a[0]) * (b[3] - b[1])
                            - (a[3] - a[1]) * (b[2] - b[0])) / (la * lb)
                if cross >= 0.08:
                    continue       # not parallel
                # Decompose the midpoint displacement in a's frame. A double
                # bond's inner line is offset PERPENDICULAR to the bond and
                # overlaps it along its length. Testing only the straight-line
                # midpoint distance also deletes the next-but-one bond of a
                # zig-zag, which is parallel to this one and roughly two bond
                # lengths away -- it collapsed a 7-vertex chain to 3 in this
                # module's own smoke test. Longitudinal displacement is the
                # discriminator, so it is tested separately.
                ux, uy = (a[2] - a[0]) / la, (a[3] - a[1]) / la
                dmx = (b[0] + b[2]) / 2 - (a[0] + a[2]) / 2
                dmy = (b[1] + b[3]) / 2 - (a[1] + a[3]) / 2
                along = abs(dmx * ux + dmy * uy)
                perp = abs(-dmx * uy + dmy * ux)
                if perp <= double_offset_pt and along <= 0.5 * la:
                    dropped.add(j)
        segs = [s for i, s in enumerate(segs) if i not in dropped]
    return len(_merge_vertices(segs, tol))


def raster_chain_length(path_xy, eps_px=3.0):
    """Turning points of a traced skeleton path, via Douglas-Peucker.

    `path_xy` is an ordered pixel path from a skeletonised component. Promoted
    from xd_2026_lipidai/stage02_trace.py. Returns the simplified vertex count.
    """
    pts = [tuple(map(float, p)) for p in path_xy]
    if len(pts) < 2:
        raise DetectionError(
            f"raster_chain_length: path has {len(pts)} point(s).")

    def rdp(p, eps):
        if len(p) < 3:
            return p
        x0, y0 = p[0]
        x1, y1 = p[-1]
        dx, dy = x1 - x0, y1 - y0
        nrm = math.hypot(dx, dy)
        best_i, best_d = 0, -1.0
        for i in range(1, len(p) - 1):
            x, y = p[i]
            d = abs(dy * (x - x0) - dx * (y - y0)) / nrm if nrm else \
                math.hypot(x - x0, y - y0)
            if d > best_d:
                best_i, best_d = i, d
        if best_d <= eps:
            return [p[0], p[-1]]
        return rdp(p[:best_i + 1], eps)[:-1] + rdp(p[best_i:], eps)

    return len(rdp(pts, eps_px))


def count_carbons(vector=None, raster=None, name_token=None, label=""):
    """The chain length, only if two independent methods agree.

    vector / raster  integer counts from the two routes above. At least one.
    name_token       the length implied by the compound's NAME, if any. It is
                     compared and reported; it is never accepted as a count.

    Returns dict(n, methods, agreed, name_token_matches).
    Raises DetectionError when two supplied counts disagree.
    """
    have = {k: v for k, v in (("vector", vector), ("raster", raster))
            if v is not None}
    if not have:
        raise DetectionError(
            f"count_carbons({label!r}): no count supplied. A chain length read "
            "by eye is error class 48 and is not an acceptable input here. Use "
            "vector_chain_length on the PDF's own coordinates, or "
            "raster_chain_length on a skeletonised trace.")
    vals = sorted(set(int(v) for v in have.values()))
    if len(vals) > 1:
        raise DetectionError(
            f"count_carbons({label!r}): methods disagree -- "
            + ", ".join(f"{k}={int(v)}" for k, v in have.items()) + ".\n"
            "A chain length is an integer and a one-carbon disagreement means "
            "one method is wrong, not that the answer is between them. "
            "Refusing to average or to prefer a route. Render higher, check "
            "the bbox, and look at the annotated overlay.")
    n = vals[0]
    res = {"n": n, "methods": sorted(have), "agreed": len(have) > 1,
           "name_token_matches": None}
    if name_token is not None:
        res["name_token_matches"] = (int(name_token) == n)
    return res


def check_name_token(token_length, measured, label=""):
    """Compare a name's implied length to the measured one. Reports, never sets.

    CLAUDE.md section 5: name tokens are not chain lengths. KS_2024's `C24` is
    two C4 tails. This returns a dict for the QC file so a mismatch is a
    recorded fact a reviewer can adjudicate, rather than either an assumption
    or a silent correction.
    """
    return {"label": label, "name_token_length": int(token_length),
            "measured_length": int(measured),
            "matches": int(token_length) == int(measured),
            "note": "" if int(token_length) == int(measured) else
                    "Name token does not equal the measured chain length. This "
                    "is NOT automatically an error -- it is how KS_2024's C24 "
                    "(two C4 tails) and KZ_2016's SC1/SC2/SC3 (isobutyl/"
                    "isopentyl/neopentyl) are written. Decide from the paper "
                    "which the token means and record it in DECISIONS.md."}
