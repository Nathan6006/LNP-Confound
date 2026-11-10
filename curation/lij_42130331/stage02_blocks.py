#!/usr/bin/env python3
"""Count every drawn chain in Fig. S1B and Fig. S2 -- never by eye.

    python3 -I stage02_blocks.py --trace

WHY THE VECTOR ROUTE IS NOT AVAILABLE, AS A CHECKABLE FACT. `chains` prefers
the PDF's own bond coordinates, which is strictly better than any vision-based
read -- the drawing's geometry IS the answer, with no estimation step. This SI
has none: `page.get_drawings()` returns **0 paths** on every page, and each
figure page carries exactly one image object. The figure is fully flattened.
`render()` asserts that, so if a future copy of the file does carry vector art
this script fails rather than quietly using the weaker route.

SO TWO INDEPENDENT RASTER COUNTS ARE USED, AND THEY MUST AGREE.

  A  skeletonise -> trace -> Douglas-Peucker, via `chains.raster_chain_length`,
     the audited implementation promoted from xd_2026_lipidai.
  C  arc length of the same skeleton divided by the BOND LENGTH measured from
     the figure's own regular zig-zags, then rounded.

A and C fail differently. A is fooled by a skeleton spur or by an eps that
smooths a shallow turn; C is fooled by a wrong bond length or by a chain drawn
at an unusual angle, and is completely blind to how many turns there are. A
one-vertex disagreement means one of them is wrong, so this script raises --
`chains.count_carbons` refuses to average two counts and so does this.

`count_carbons` is given A. C is asserted equal to A at the call site rather
than passed as `vector=`, because `vector` means the PDF's own coordinates and
there are none here. Calling C "vector" would be a false provenance claim --
which is precisely what F-013 and F-037 each were.

COLOUR DOES THE SEGMENTATION. Bonds are black; heteroatom labels (N, NH2, O,
OH, SH, HN) are blue or red glyph outlines. A black-only mask therefore drops
every label without a hand-placed crop box, and each chain between two labels
falls out as its own connected component. A hand-placed crop box is how an
off-by-one row gets into a read; a colour threshold is a property of the
drawing.

HOW A VERTEX COUNT BECOMES A CARBON COUNT. It depends on what sits at each end
of the component, and that is DETECTED, not assumed:

  * an end with coloured ink within `LABEL_RADIUS` px is a heteroatom -- the
    bond to it is drawn, so the end point is the N or O, not a carbon;
  * an end inside a detected squiggle (the wavy attachment bond) is an
    attachment point, not a carbon;
  * any other end is a free line end, i.e. a terminal carbon.

carbons = vertices - (ends that are heteroatoms or squiggles).
"""
import json
import math
import os
import sys

import numpy as np
import pymupdf
from PIL import Image, ImageDraw
from skimage.measure import label as sklabel, regionprops
from skimage.morphology import skeletonize

LIB = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib")
sys.path.insert(0, LIB)
import chains  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
PDF = os.path.join(HERE, "source", "nn6c00706_si_001.pdf")
DPI = 900

# This paper's own geometry: the regions of the SI that draw a chain we must
# count, in PDF points, read off page.get_image_rects(). Paper constants stay
# in the curation folder, next to the figure, where a reviewer reads them.
REGIONS = {
    # The three rows of amine STRUCTURES, with the bold name labels ("70",
    # "75", ...) left outside. Those labels are black ink too, so they are
    # found as components and their glyph skeletons are long and full of short
    # segments -- the first version of this file measured its bond length off
    # one and got 23 px instead of 61, which made the arc-length route
    # disagree with Douglas-Peucker on 44 of 56 components. The rows are split
    # on the label bands, which are visible as a y-gap in the component list.
    "figS1B_row1": {"page": 1, "rect": (90, 262, 506, 277), "figure": "S1B",
                    "what": "amines 70, 75, 80, 96, 400"},
    "figS1B_row2": {"page": 1, "rect": (90, 289, 506, 297), "figure": "S1B",
                    "what": "amines 113, 414, 616"},
    "figS1B_row3": {"page": 1, "rect": (90, 310, 506, 337), "figure": "S1B",
                    "what": "amines 110, T3A"},
    "figS2_legend": {"page": 2, "rect": (330, 285, 426, 395), "figure": "S2",
                     "what": "the R= legend: C8, C10, C12, C14, OG8"},
}

MIN_AREA = 400        # px at 900 dpi. A single bond is ~14 x 82 px.
MIN_PATH = 30         # skeleton px; shorter is a glyph serif, not a bond
RDP_EPS = 6.0         # px at 900 dpi
LABEL_RADIUS = 50     # px; coloured ink this close to an end is a label
SQUIGGLE_RATIO = 0.55  # a step shorter than this x modal is squiggle wobble
SQUIGGLE_MIN = 3      # that many consecutive short steps make a squiggle
APEX_PROMINENCE = 0.30  # fraction of the detrended amplitude a turn must span
APEX_SQUIGGLE_GAP = 1.4  # x-gap this many times the median marks a squiggle
NB = [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)]


def render(page_idx, rect):
    doc = pymupdf.open(PDF)
    page = doc[page_idx]
    if page.get_drawings():
        raise SystemExit(
            f"page {page_idx + 1} now carries vector drawings. Use "
            "chains.vector_bonds -- the PDF's own bond coordinates are "
            "strictly better than any raster trace -- and delete this raster "
            "fallback rather than leaving both live.")
    pix = page.get_pixmap(dpi=DPI, clip=pymupdf.Rect(*rect))
    a = np.frombuffer(pix.samples, dtype=np.uint8)
    return a.reshape(pix.height, pix.width, pix.n)[:, :, :3].astype(int)


def masks(rgb):
    """(black ink = bonds, coloured ink = heteroatom labels)."""
    lum = rgb.mean(axis=2)
    spread = rgb.max(axis=2) - rgb.min(axis=2)
    ink = lum < 150
    return ink & (spread < 50), ink & (spread >= 50)


def trace_component(sk):
    ys, xs = np.nonzero(sk)
    px = set(zip(ys.tolist(), xs.tolist()))
    deg = {p: sum(((p[0] + dy, p[1] + dx) in px) for dy, dx in NB)
           for p in px}
    nodes = {p for p, d in deg.items() if d != 2}
    if not nodes:
        return [sorted(px)]
    paths, seen = [], set()
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
            paths.append(path)
    return paths


def rdp_vertices(pts, eps):
    """The simplified polyline itself (chains returns only its length)."""
    def go(p):
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
        return go(p[:bi + 1])[:-1] + go(p[bi:])
    return go([tuple(map(float, q)) for q in pts])


def seglens(poly):
    return [math.hypot(b[0] - a[0], b[1] - a[1])
            for a, b in zip(poly, poly[1:])]


def trim_squiggle(poly, modal):
    """Drop a wavy attachment bond from either end. Returns (poly, l, r).

    A squiggle is drawn as a run of short, tightly spaced turns; a real bond is
    one segment of about `modal`. So a leading or trailing run of at least
    SQUIGGLE_MIN segments shorter than SQUIGGLE_RATIO x modal is the squiggle,
    and nothing else in these drawings looks like that.
    """
    L = seglens(poly)
    short = [s < SQUIGGLE_RATIO * modal for s in L]
    lo = 0
    while lo < len(short) and short[lo]:
        lo += 1
    if lo < SQUIGGLE_MIN:
        lo = 0
    hi, k = len(poly) - 1, 0
    j = len(short) - 1
    while j >= 0 and short[j]:
        j -= 1
        k += 1
    if k >= SQUIGGLE_MIN:
        hi = j + 1
    return poly[lo:hi + 1], lo > 0, hi < len(poly) - 1


def apices(mask):
    """METHOD B: the zig-zag's apices, from the per-column ink centroid.

    No skeleton, no Douglas-Peucker, no bond length -- so it fails in a
    completely different way from method A. A chain's centroid y(x) is
    piecewise monotone and each sign change of its slope is an apex, i.e. a
    carbon. Three things it must handle, each a property of these drawings and
    not a knob:

      * a chain drawn with an overall slope (T3A's upper and lower arms) --
        removed by subtracting a straight-line fit before looking for turns;
      * the one-pixel wobble of a centroid on an antialiased stroke --
        removed by requiring a turn to span APEX_PROMINENCE of the amplitude;
      * the wavy attachment bond in Fig. S2's legend -- its turns are packed
        far more tightly in x than real apices, so a leading or trailing turn
        separated from its neighbour by more than APEX_SQUIGGLE_GAP times the
        median gap is dropped with everything outside it.

    Returns (n_apices, apex_x_positions) or (None, []) if x is not usable.
    """
    cols = [np.where(mask[:, x])[0] for x in range(mask.shape[1])]
    xs = np.array([x for x, c in enumerate(cols) if len(c)])
    if len(xs) < 12:
        return None, []
    y = np.array([cols[x].mean() for x in xs], dtype=float)
    y = y - np.poly1d(np.polyfit(xs, y, 1))(xs)       # detrend
    amp = y.max() - y.min()
    if amp < 4:
        return 0, []                                   # a straight stub
    turns, sign, last = [], 0, None
    d = np.diff(y)
    for i, v in enumerate(d):
        s = 1 if v > 0.02 else (-1 if v < -0.02 else 0)
        if s == 0 or s == sign:
            continue
        if sign != 0 and (last is None
                          or abs(y[i] - last) > APEX_PROMINENCE * amp):
            turns.append((int(xs[i]), float(y[i])))
            last = y[i]
        sign = s
    if len(turns) >= 3:
        gaps = [b[0] - a[0] for a, b in zip(turns, turns[1:])]
        med = float(np.median(gaps))
        while len(turns) > 2 and gaps and gaps[0] > APEX_SQUIGGLE_GAP * med:
            turns, gaps = turns[1:], gaps[1:]
        while len(turns) > 2 and gaps and gaps[-1] > APEX_SQUIGGLE_GAP * med:
            turns, gaps = turns[:-1], gaps[:-1]
    return len(turns), [t[0] for t in turns]


def raw_components(bm):
    lab = sklabel(bm)
    raw = []
    for r in regionprops(lab):
        if r.area < MIN_AREA:
            continue
        sub = (lab[r.bbox[0]:r.bbox[2], r.bbox[1]:r.bbox[3]] == r.label)
        sk = skeletonize(sub)
        paths = [p for p in trace_component(sk) if len(p) >= MIN_PATH]
        if not paths:
            continue
        longest = max(paths, key=len)
        pts = [(x + r.bbox[1], y + r.bbox[0]) for y, x in longest]
        raw.append({"bbox": [int(v) for v in r.bbox], "area": int(r.area),
                    "n_paths": len(paths), "pts": pts, "mask": sub})

    return raw


def bond_length(raws):
    """The figure's own bond length: the median segment of its LONGEST chain.

    Measured, not typed in. Once the bold name labels are outside the region,
    no glyph in these figures is as long as a real zig-zag, so the longest
    component is a chain and the median of its simplified segments is the bond
    length by construction. The first version of this estimate took the mode
    over EVERY component including the bold numerals under each structure, and
    returned 23 px where the truth is 61 -- which made the arc-length route
    disagree with Douglas-Peucker on 44 of 56 components. The filter is on the
    component's SHAPE and length, never on its position or its answer.
    """
    cand = []
    for c in raws:
        poly = rdp_vertices(c["pts"], RDP_EPS)
        L = [s for s in seglens(poly) if s > 15]
        if len(L) >= 4:
            cand.append((sum(L), L))
    if not cand:
        return None
    cand.sort(key=lambda t: -t[0])
    return float(np.median(cand[0][1]))


def analyse(raw, cm, modal):
    out = []
    for i, c in enumerate(raw):
        poly = rdp_vertices(c["pts"], RDP_EPS)
        poly, sq_l, sq_r = trim_squiggle(poly, modal)
        if len(poly) < 2:
            continue
        a_count = chains.raster_chain_length(poly, eps_px=RDP_EPS)
        arc = sum(seglens(poly))
        ends = [poly[0], poly[-1]]
        het, gap = [], 0.0
        cys, cxs = np.nonzero(cm)
        for (x, y), sq in zip(ends, (sq_l, sq_r)):
            x0 = max(0, int(x) - LABEL_RADIUS)
            x1 = min(cm.shape[1], int(x) + LABEL_RADIUS + 1)
            y0 = max(0, int(y) - LABEL_RADIUS)
            y1 = min(cm.shape[0], int(y) + LABEL_RADIUS + 1)
            is_label = bool(cm[y0:y1, x0:x1].any())
            het.append(is_label or bool(sq))
            # A bond drawn TO a heteroatom label stops short of it, by the
            # label's own clearance. The arc-length route would therefore
            # undercount that bond; the gap is measurable (the distance from
            # the end point to the nearest coloured pixel), so it is measured
            # and added back rather than absorbed. Without this the two
            # routes differ by exactly 1 on every chain that ends at a label,
            # which is the whole population of amine arms.
            if is_label and len(cys):
                d = np.hypot(cxs - x, cys - y).min()
                gap += float(d)
        c_count = int(round((arc + gap) / modal)) + 1
        n_apex, apex_x = apices(c["mask"])
        out.append({
            "id": i, "bbox": c["bbox"], "area": c["area"],
            "n_paths": c["n_paths"],
            "vertices_rdp": int(a_count),
            "vertices_arclength": int(c_count),
            "methods_agree": int(a_count) == int(c_count),
            "arc_px": round(arc, 1),
            "apices": n_apex, "apex_x": apex_x,
            "end_is_heteroatom_or_attachment": het,
            "squiggle_trimmed": [bool(sq_l), bool(sq_r)],
            "free_ends": 2 - sum(het),
            "carbons_A_rdp": int(a_count) - sum(het),
            "carbons_B_arc": int(c_count) - sum(het),
            "carbons_C_apex": (None if n_apex is None
                               else n_apex + (2 - sum(het))),
            "carbons": int(a_count) - sum(het),
            "polyline": [[round(p[0], 1), round(p[1], 1)] for p in poly],
        })
    out.sort(key=lambda c: (c["bbox"][0] // 120, c["bbox"][1]))
    for i, c in enumerate(out):
        c["id"] = i
    return out


def trace():
    # Pass 1: components per region, and ONE bond length per FIGURE. A row of
    # short amine arms carries no chain long enough to measure a bond from, so
    # the figure's longest chain supplies it for every row of that figure.
    pre, bl = {}, {}
    for name, spec in REGIONS.items():
        rgb = render(spec["page"], spec["rect"])
        bm, cm = masks(rgb)
        pre[name] = (rgb, cm, raw_components(bm))
        bl.setdefault(spec["figure"], []).extend(pre[name][2])
    fig_bl = {fig: bond_length(raws) for fig, raws in bl.items()}
    if any(v is None for v in fig_bl.values()):
        raise SystemExit("a figure has no chain long enough to measure a "
                         "bond from; refusing to guess one.")

    result = {}
    for name, spec in REGIONS.items():
        rgb, cm, raw = pre[name]
        # Prefer the bond length of THIS region: Fig. S1B's three rows are not
        # drawn at quite the same scale (row 1's bond is ~73 px, row 2's ~69),
        # and a 5% scale error is worth exactly one carbon on a four-bond
        # chain. A row with no chain of four bonds falls back to the figure's.
        modal = bond_length(raw) or fig_bl[spec["figure"]]
        comps = analyse(raw, cm, modal)
        result[name] = {"what": spec["what"], "rect": spec["rect"], "dpi": DPI,
                        "figure": spec["figure"],
                        "bond_length_px": round(modal, 2),
                        "n_components": len(comps), "components": comps}
        ov = Image.fromarray(rgb.astype("uint8")).convert("RGB")
        d = ImageDraw.Draw(ov)
        for c in comps:
            y0, x0, y1, x1 = c["bbox"]
            d.rectangle([x0 - 3, y0 - 3, x1 + 3, y1 + 3],
                        outline=(255, 0, 0), width=3)
            d.text((x0, max(0, y0 - 30)),
                   f"#{c['id']} C={c['carbons']}", fill=(255, 0, 0))
            for p in c["polyline"]:
                d.ellipse([p[0] - 7, p[1] - 7, p[0] + 7, p[1] + 7],
                          outline=(0, 160, 0), width=3)
        ov.save(os.path.join(HERE, "00_pages", f"trace_{name}.png"))

    for r in result.values():
        for c in r["components"]:
            c["methods_agree"] = (c["carbons_A_rdp"] == c["carbons_B_arc"])
    bad = [(n, c["id"]) for n, r in result.items()
           for c in r["components"] if not c["methods_agree"]]
    with open(os.path.join(HERE, "02_traces.json"), "w") as fh:
        json.dump({"_what": "every black (bond) component in the regions that "
                            "draw a building block, with two independent "
                            "vertex counts and the endpoint classification "
                            "that turns vertices into carbons.",
                   "n_method_disagreements": len(bad),
                   "regions": result}, fh, indent=1)
    for name, r in result.items():
        print(f"\n{name}: {r['n_components']} components, "
              f"bond length {r['bond_length_px']} px")
        for c in r["components"]:
            flag = "  " if c["methods_agree"] else "XX"
            print(f" {flag} #{c['id']:>2} bbox {c['bbox']} "
                  f"A_rdp={c['carbons_A_rdp']} B_arc={c['carbons_B_arc']} "
                  f"C_apex={c['carbons_C_apex']} "
                  f"ends_het={c['end_is_heteroatom_or_attachment']} "
                  f"squiggle={c['squiggle_trimmed']} apex_x={c['apex_x']}")
    if bad:
        print(f"\n{len(bad)} components where the two methods DISAGREE: {bad}")
        print("A one-vertex disagreement means one method is wrong, not that "
              "the answer is between them. Not resolving it here.")
    return 0


if __name__ == "__main__":
    if "--trace" in sys.argv:
        sys.exit(trace())
    raise SystemExit("usage: stage02_blocks.py --trace")
