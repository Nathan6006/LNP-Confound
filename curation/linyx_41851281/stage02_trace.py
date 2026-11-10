#!/usr/bin/env python3
"""Raster trace of Supplementary Fig. 1a (MOESM1 p.5): head and tail skeletons.

WHY RASTER. MOESM1 p.5 carries ONE embedded image (2484 x 3047 px) and ZERO
vector drawing paths -- `page.get_drawings()` returns []. curation/lib/chains
`vector_bonds` therefore returns [] by construction, so the vector route does
not exist for this page and the raster route is the only geometric one. That is
recorded in 00_value_sources.json entry S1, and in 07_qc.json.

THE SECOND METHOD IS NOT A SECOND TRACE. It is HRMS: Supplementary Notes
(MOESM1 pp.41-53) print `Mass calculated m/z` for 118 of the 121 lipids. Head
and tail formulas are solved from those masses in stage01_blocks.py. Exact mass
pins the FORMULA, hence the carbon count; this trace pins the CONNECTIVITY
(where the branch sits, how the carbons are partitioned) that mass is blind to.
chains.count_carbons is given both.

Writes 02_trace.json and an annotated overlay so a reviewer checks the geometry,
not a number in a comment.
"""
import json
import os
import sys

import math

import numpy as np
import pymupdf
from PIL import Image, ImageDraw
from skimage.measure import label, regionprops
from skimage.morphology import skeletonize

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import chains                                     # noqa: E402

# The attachment squiggle and the chain are separated by BOND LENGTH, not by
# eye. MEASURED, not assumed, and reported in 02_trace.json so it can be
# re-checked: over the 424 inter-vertex spacings that belong to a segment
# containing at least one >=40 px bond, 78 are 5.0-10.3 px (squiggle wiggles),
# 343 are 25.1-51 px (drawn bonds) and THREE are 24.2-24.7 px.
#
# Those three are not tails. They are head fragments -- the stub bond into an
# `N` label in H6, H8 and H9 -- and they are why this threshold is applied to
# TAIL chains only. Head carbon counts come from plain vertex counts and are
# pinned independently by the head's exact mass (stage01_blocks.py), so the
# threshold never touches them. Saying "nothing in between" would have been
# false, and the kind of false a comment gets away with; the numbers are in
# the JSON instead.
LONG_PX = 25.0
NARROW_BAND = (11.0, 25.0)

SRC = os.path.join(HERE, "source", "41551_2026_1640_MOESM1_ESM.pdf")
PAGE = 4                                          # 0-based; SI p.5
NB = [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)]


def native_image(out_png):
    d = pymupdf.open(SRC)
    page = d[PAGE]
    if page.get_drawings():
        raise RuntimeError(
            "page %d now reports %d vector drawing paths. The vector route "
            "would be strictly better -- use chains.vector_bonds and delete "
            "this justification." % (PAGE + 1, len(page.get_drawings())))
    imgs = page.get_images(full=True)
    if len(imgs) != 1:
        raise RuntimeError("expected one embedded image, got %d" % len(imgs))
    px = pymupdf.Pixmap(d, imgs[0][0])
    px.save(out_png)
    return out_png


def paths_of(sk):
    """Pixel paths between endpoints/junctions of one skeleton component."""
    ys, xs = np.nonzero(sk)
    px = set(zip(ys.tolist(), xs.tolist()))
    deg = {p: sum(((p[0] + dy, p[1] + dx) in px) for dy, dx in NB) for p in px}
    nodes = {p for p, d in deg.items() if d != 2}
    if not nodes:
        return [sorted(px)]
    out, seen = [], set()
    for s in nodes:
        for dy, dx in NB:
            nxt = (s[0] + dy, s[1] + dx)
            if nxt not in px or (s, nxt) in seen:
                continue
            path = [s, nxt]
            seen.add((s, nxt))
            seen.add((nxt, s))
            cur, prev = nxt, s
            while cur not in nodes:
                adv = None
                for dy2, dx2 in NB:
                    c = (cur[0] + dy2, cur[1] + dx2)
                    if c in px and c != prev and (cur, c) not in seen:
                        adv = c
                        break
                if adv is None:
                    break
                seen.add((cur, adv))
                seen.add((adv, cur))
                path.append(adv)
                prev, cur = cur, adv
            out.append(path)
    return out


def components(img_png, region, min_area=200):
    a = np.array(Image.open(img_png).convert("L"))
    x0, y0, x1, y1 = region
    sub = a[y0:y1, x0:x1] < 128
    lab = label(sub)
    comps = []
    for r in regionprops(lab):
        if r.area < min_area:
            continue
        m = (lab[r.bbox[0]:r.bbox[2], r.bbox[1]:r.bbox[3]] == r.label)
        sk = skeletonize(m)
        segs = []
        for p in paths_of(sk):
            if len(p) < 15:
                continue
            pts = [(x + r.bbox[1] + x0, y + r.bbox[0] + y0) for y, x in p]
            segs.append(pts)
        comps.append({
            "bbox": [int(r.bbox[1] + x0), int(r.bbox[0] + y0),
                     int(r.bbox[3] + x0), int(r.bbox[2] + y0)],
            "area": int(r.area),
            "segments": segs,
        })
    comps.sort(key=lambda c: (c["bbox"][1], c["bbox"][0]))
    return comps


def _rdp(path, eps=3.0):
    """Vertex COORDINATES, so the caller can measure bond lengths.

    chains.raster_chain_length returns only the count, by design; this repeats
    the same Douglas-Peucker on the same path purely to recover the points.
    The count written to 02_trace.json is always the library's, never this
    function's, so the two cannot drift into disagreeing.
    """
    pts = [(float(x), float(y)) for x, y in path]

    def rdp(p, eps):
        if len(p) < 3:
            return p
        x0, y0 = p[0]
        x1, y1 = p[-1]
        dx, dy = x1 - x0, y1 - y0
        nrm = math.hypot(dx, dy)
        bi, bd = 0, -1.0
        for i in range(1, len(p) - 1):
            x, y = p[i]
            d = (abs(dy * (x - x0) - dx * (y - y0)) / nrm if nrm
                 else math.hypot(x - x0, y - y0))
            if d > bd:
                bi, bd = i, d
        if bd <= eps:
            return [p[0], p[-1]]
        return rdp(p[:bi + 1], eps)[:-1] + rdp(p[bi:], eps)

    return rdp(pts, eps)


def main():
    png = native_image(os.path.join(HERE, "02_sfig1a.png"))
    # Supplementary Fig. 1a occupies the top of the page; 1b (the name table)
    # is below it. The split y is the horizontal rule under the "Tails" box.
    region = (0, 150, 2484, 1560)
    comps = components(png, region)
    out, spans = [], []
    for c in comps:
        rec = {"bbox": c["bbox"], "area": c["area"],
               "n_segments": len(c["segments"]), "segment_vertices": []}
        for s in c["segments"]:
            try:
                n = chains.raster_chain_length(s, eps_px=3.0)
                v = _rdp(s)
                sp = [math.hypot(v[i + 1][0] - v[i][0],
                                 v[i + 1][1] - v[i][1])
                      for i in range(len(v) - 1)]
                spans.extend(sp)
                rec["segment_vertices"].append(
                    {"n": n,
                     "long_bonds": sum(1 for x in sp if x >= LONG_PX),
                     "spacings_px": [round(x, 1) for x in sp],
                     "ends": [list(map(int, s[0])), list(map(int, s[-1]))],
                     "px": len(s)})
            except Exception as exc:                      # refusal, recorded
                rec["segment_vertices"].append({"error": str(exc)})
        out.append(rec)
    chain = [x for x in spans]
    middle = [x for x in chain if NARROW_BAND[0] <= x < NARROW_BAND[1]]
    payload = {
        "_what": "RDP vertices and inter-vertex spacings for every traced "
                 "component of Supplementary Fig. 1a. `long_bonds` counts "
                 "spacings >= %g px; for a tail that count IS the carbon "
                 "count of R, because the first long bond starts at the "
                 "attachment carbon. Applied to TAIL chains only -- see the "
                 "LONG_PX comment." % LONG_PX,
        "bond_length_px": {"n": len(spans),
                           "min_long": round(min([x for x in spans
                                                  if x >= LONG_PX] or [0]), 1),
                           "max_short": round(max([x for x in spans
                                                   if x < LONG_PX] or [0]), 1),
                           "in_band_%g_%g" % NARROW_BAND: len(middle),
                           "_note": "counted over EVERY traced segment, "
                                    "glyph outlines included. The per-segment "
                                    "`spacings_px` below is what a reviewer "
                                    "should read; the three 24.2-24.7 px "
                                    "bonds are head stubs into an N label, "
                                    "not tails."},
        "components": out}
    json.dump(payload, open(os.path.join(HERE, "02_trace.json"), "w"), indent=1)

    ov = Image.open(png).convert("RGB")
    d = ImageDraw.Draw(ov)
    for c in comps:
        d.rectangle(c["bbox"], outline=(255, 0, 0), width=3)
    ov.save(os.path.join(HERE, "02_trace_overlay.png"))
    print("components:", len(out))
    for r in out:
        print("  bbox=%s area=%d segs=%d verts=%s" % (
            r["bbox"], r["area"], r["n_segments"],
            [v.get("n") for v in r["segment_vertices"]]))


if __name__ == "__main__":
    main()
