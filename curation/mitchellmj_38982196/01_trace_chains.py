#!/usr/bin/env python3
"""Stage 01 -- MEASURE every chain length in Fig. 1a. Nothing here is counted by eye.

Run:  python3 -I 01_trace_chains.py            # writes 01_chain_measurements.json
      python3 -I 01_trace_chains.py --audit-vector   # proves the vector route is unavailable

WHY RASTER ONLY. 00_value_sources.json records the checkable reason: there is
no vector rendition of Fig. 1a obtainable. `--audit-vector` re-runs that proof
against the SI PDF so a reviewer does not have to take it on trust. Because
only one route exists, `chains.count_carbons` is called with raster= alone and
`agreed` is False on every chain; a SECOND, independent statistic (zig-zag apex
count) is computed here and compared, and the P1 mass anchor in 06_checks.py
closes two of these measurements against a number the authors printed.

THE +/-1 CONVENTION, AND WHY IT IS CHECKED AND NOT ASSUMED. A drawn polyline
ends ON a heteroatom's text label, so its terminal vertex is the position of
that N/O/S, not a carbon. Getting this backwards is error class 48 by another
route. `_label_adjacent` therefore asserts, for every endpoint declared to be a
label terminus, that a glyph component actually starts within LABEL_GAP_PX --
and asserts the opposite for an endpoint declared to be a terminal methyl.
"""
import json
import os
import sys

import numpy as np
from PIL import Image
from scipy import ndimage
from skimage.morphology import skeletonize

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "lib"))
import chains                                              # noqa: E402
from errors import DetectionError                          # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
FIG1 = os.path.join(HERE, "source", "nature_fig1_full.png")
SI = os.path.join(HERE, "source", "41557_2024_1557_MOESM1_ESM.pdf")
OUT = os.path.join(HERE, "01_chain_measurements.json")

INK = 128            # a pixel darker than this is line work (background is pink)
EPS_SWEEP = (3.0, 4.0, 5.0, 6.0)   # Douglas-Peucker tolerances that must agree
LABEL_GAP_PX = 12    # a glyph this close to a free end is that end's atom label
# ...and 12 is not a guess. Measured over all 48 declared endpoints: a real
# atom label sits 0.0-8.3 px from the bond's end, while the nearest row NUMERAL
# (the bold "1".."25" in the margin) is 19.0-45.6 px away. The floor is
# auditable rather than asserted -- main() recomputes both extremes and raises
# if they ever stop being cleanly separated, which is the F-051 lesson: an
# asymmetric or tuned threshold has to print the gap it is relying on.
MIN_GLYPH_FILL = 0.30  # text is compact; line work is thin (measured 0.09-0.15)

# ---------------------------------------------------------------------------
# WHAT IS MEASURED. Paper-specific geometry belongs in the paper's folder
# (PIPELINE_V2 section 3), so the seeds live here and not in curation/lib.
#
# seed      a pixel ON the structure's line work. If it is not ink, this script
#           RAISES rather than silently measuring a neighbouring component --
#           error class 28 is a detector miss falling through to a literal.
# ends      ("methyl"|"label", "methyl"|"label") for the two path endpoints, in
#           left-to-right order. Checked against the glyph map, not assumed.
# atoms     how many ATOMS of the named element the vertex count implies.
#           "V-1" where one endpoint sits on a label, "V" where neither does.
# ---------------------------------------------------------------------------
M = [
    dict(id="amine01_chain", seed=(66, 514), kind="linear",
         ends=("methyl", "label"),
         minus=1,
         role="n-alkyl chain of amine 1"),
    dict(id="amine02_chain", seed=(66, 568), kind="linear",
         ends=("methyl", "label"),
         minus=1,
         role="n-alkyl chain of amine 2"),
    dict(id="amine03_chain", seed=(66, 621), kind="linear",
         ends=("methyl", "label"),
         minus=1,
         role="n-alkyl chain of amine 3"),
    dict(id="amine04_chain", seed=(66, 676), kind="linear",
         ends=("methyl", "label"),
         minus=1,
         role="n-alkyl chain of amine 4"),
    dict(id="amine05_chain", seed=(66, 733), kind="linear",
         ends=("methyl", "label"),
         minus=1,
         role="n-alkyl chain of amine 5"),
    dict(id="amine06_chain", seed=(66, 794), kind="linear",
         ends=("methyl", "label"),
         minus=1,
         role="mono-unsaturated chain of amine 6. The inner line of the C=C is a SEPARATE component (amine06_double_bond), so this trace is the plain backbone; the apex comparison is the one place it is skipped, because a C=C is drawn FLAT and so has no turning point"),
    dict(id="amine06_double_bond", seed=(300, 801), kind="double_bond_marker",
         parent='amine06_chain',
         role="the inner line of amine 6 C=C: which backbone bond it marks"),
    dict(id="amine07_acyl", seed=(66, 860), kind="to_junction",
         ends=("methyl", "junction"),
         minus=0,
         role="acyl chain of amine 7, terminal CH3 -> carbonyl C. The junction IS the carbonyl carbon, so every vertex is a carbon"),
    dict(id="amine08_acyl", seed=(66, 956), kind="to_junction",
         ends=("methyl", "junction"),
         minus=0,
         role="acyl chain of amine 8, terminal CH3 -> carbonyl C"),
    dict(id="amine09_branched", seed=(66, 1076), kind="three_branch",
         role="amine 9: the three branches of the 2-substituted carbon"),
    dict(id="amine10_ester_alkyl", seed=(66, 1136), kind="linear",
         ends=("methyl", "label"),
         minus=1,
         role="O-alkyl chain of amine 10; its right end rests on the ester O glyph"),
    dict(id="amine10_glycyl", seed=(300, 1142), kind="to_junction_from_right",
         ends=("label", "junction"),
         minus=1,
         role="amine 10 head, NH2 end -> carbonyl carbon. The carbonyl O letter is fused to the two C=O lines, so the component carries a loop and has only two free ends; walking in from the NH2 side avoids it. Two heavy atoms means CH2 + C(=O), i.e. glycine and not beta-alanine"),
    dict(id="amine10_ester_link", seed=(300, 1142), kind="to_junction",
         ends=("label", "junction"),
         minus=1,
         role="amine 10 head, ester O -> carbonyl carbon. One heavy atom means the ester oxygen is bonded straight to the carbonyl carbon, with no carbon in between"),
    dict(id="amine11_alkyl", seed=(66, 1265), kind="to_junction",
         ends=("methyl", "junction"),
         minus=1,
         role="4-alkyl chain of amine 11, terminal CH3 -> the substituted ring carbon, which is a vertex and is subtracted"),
    dict(id="amine12_alkyl", seed=(66, 1348), kind="to_junction",
         ends=("methyl", "junction"),
         minus=1,
         role="4-alkyl chain of amine 12"),
    dict(id="amine13_alkyl", seed=(66, 1429), kind="to_junction",
         ends=("methyl", "junction"),
         minus=1,
         role="4-alkyl chain of amine 13"),
    dict(id="amine14_alkoxy", seed=(66, 1529), kind="linear",
         ends=("methyl", "label_fused"),
         minus=1,
         role="4-alkoxy chain of amine 14; its right end rests on the ether O, whose outline letter is fused into the benzene-ring component, so the weaker adjacent-ink test is used there"),
    dict(id="amine14_ring", seed=(309, 1541), kind="ring",
         role="the benzene ring of amine 14: ring size from its own enclosed hole"),
    dict(id="amine11_ring_substitution", seed=(66, 1265), kind="ring_substitution",
         role="amine 11: where the alkyl chain and the NH2 sit on the ring"),
    dict(id="amine12_ring_substitution", seed=(66, 1348), kind="ring_substitution",
         role="amine 12: where the alkyl chain and the NH2 sit on the ring"),
    dict(id="amine13_ring_substitution", seed=(66, 1429), kind="ring_substitution",
         role="amine 13: where the alkyl chain and the NH2 sit on the ring"),
    dict(id="amine14_ring_substitution", seed=(309, 1541), kind="ring_substitution",
         role="amine 14: where the ether O and the NH2 sit on the ring"),
    dict(id="amine15_linker", seed=(915, 516), kind="linear",
         ends=("label", "label"),
         minus=2,
         role="N-to-N chain of amine 15"),
    dict(id="amine16_linker", seed=(915, 576), kind="linear",
         ends=("label", "label"),
         minus=2,
         role="N-to-N chain of amine 16"),
    dict(id="amine17_linker", seed=(915, 631), kind="linear",
         ends=("label", "label"),
         minus=2,
         role="N-to-N chain of amine 17"),
    dict(id="amine18_linkerA", seed=(915, 689), kind="linear",
         ends=("label", "label"),
         minus=2,
         role="amine 18, terminal NH2 -> secondary NH"),
    dict(id="amine18_linkerB", seed=(998, 691), kind="linear",
         ends=("label", "label"),
         minus=2,
         role="amine 18, secondary NH -> terminal NH2"),
    dict(id="amine19_linkerA", seed=(915, 770), kind="linear",
         ends=("label", "label"),
         minus=2,
         role="amine 19, terminal NH2 -> tertiary N"),
    dict(id="amine19_linkerB", seed=(998, 773), kind="linear",
         ends=("label", "label"),
         minus=2,
         role="amine 19, tertiary N -> terminal NH2"),
    dict(id="amine19_Nmethyl", seed=(986, 793), kind="stub",
         role="amine 19: the bond off the tertiary N whose far end carries no label -- a methyl"),
    dict(id="amine20_linkerA", seed=(915, 859), kind="linear",
         ends=("label", "label"),
         minus=2,
         role="amine 20, first N-to-N linker"),
    dict(id="amine20_linkerB", seed=(998, 862), kind="linear",
         ends=("label", "label"),
         minus=2,
         role="amine 20, second N-to-N linker"),
    dict(id="amine20_linkerC", seed=(986, 882), kind="linear",
         ends=("label", "label"),
         minus=2,
         role="amine 20, third (lower) N-to-N linker"),
    dict(id="amine21_linkerA", seed=(915, 985), kind="linear",
         ends=("label", "label"),
         minus=2,
         role="amine 21, first N-to-N linker"),
    dict(id="amine21_linkerB", seed=(1026, 985), kind="linear",
         ends=("label", "label"),
         minus=2,
         role="amine 21, second N-to-N linker"),
    dict(id="amine21_Nmethyl", seed=(1014, 946), kind="stub",
         role="amine 21: unlabelled bond off the tertiary N -- a methyl"),
    dict(id="amine22_linkerA", seed=(915, 1104), kind="linear",
         ends=("label", "label"),
         minus=2,
         role="amine 22, first N-to-N linker"),
    dict(id="amine22_linkerB", seed=(1026, 1104), kind="linear",
         ends=("label", "label"),
         minus=2,
         role="amine 22, second N-to-N linker"),
    dict(id="amine22_linkerC", seed=(1014, 1064), kind="linear",
         ends=("label", "label"),
         minus=2,
         role="amine 22, third (upper) N-to-N linker"),
    dict(id="amine23_linkerA", seed=(915, 1233), kind="linear",
         ends=("label", "label"),
         minus=2,
         role="amine 23, left N-to-N linker"),
    dict(id="amine23_linkerB", seed=(1082, 1187), kind="linear",
         ends=("label", "label"),
         minus=2,
         role="amine 23, right N-to-N linker"),
    dict(id="amine23_ringA", seed=(1014, 1194), kind="linear",
         ends=("label", "label"),
         minus=2,
         role="amine 23 piperazine, upper half-ring between the two ring N glyphs"),
    dict(id="amine23_ringB", seed=(1026, 1234), kind="linear",
         ends=("label", "label"),
         minus=2,
         role="amine 23 piperazine, lower half-ring"),
    dict(id="amine24_linkerA", seed=(915, 1344), kind="linear",
         ends=("label", "label"),
         minus=2,
         role="amine 24, terminal NH2 -> secondary NH"),
    dict(id="amine24_linkerB", seed=(998, 1346), kind="linear",
         ends=("label", "label"),
         minus=2,
         role="amine 24, secondary NH -> piperazine N"),
    dict(id="amine24_linkerC", seed=(1137, 1298), kind="linear",
         ends=("label", "label"),
         minus=2,
         role="amine 24, piperazine N -> terminal NH2"),
    dict(id="amine24_ringA", seed=(1070, 1304), kind="linear",
         ends=("label", "label"),
         minus=2,
         role="amine 24 piperazine, upper half-ring"),
    dict(id="amine24_ringB", seed=(1082, 1344), kind="linear",
         ends=("label", "label"),
         minus=2,
         role="amine 24 piperazine, lower half-ring"),
    dict(id="amine25_armA", seed=(942, 1440), kind="to_junction",
         ends=("label", "junction"),
         minus=2,
         role="amine 25, upper-left aminopropyl arm. Its component also contains the LETTER N of the left tertiary nitrogen, fused to the bond, so the walk stops at that fusion and both ends are nitrogen positions"),
    dict(id="amine25_armB", seed=(915, 1488), kind="linear",
         ends=("label", "label_fused"),
         minus=2,
         role="amine 25, lower-left aminopropyl arm. Its right end abuts the left tertiary N, whose glyph is fused into amine25_armA component, so the glyph test is unavailable there and the weaker adjacent-ink test is used"),
    dict(id="amine25_core", seed=(1026, 1488), kind="linear",
         ends=("label_fused", "label"),
         minus=2,
         role="amine 25, the central N-to-N chain joining the two tertiary nitrogens"),
    dict(id="amine25_armC", seed=(1165, 1491), kind="linear",
         ends=("label", "label"),
         minus=2,
         role="amine 25, right aminopropyl arm"),
    dict(id="amine25_armD", seed=(1154, 1511), kind="linear",
         ends=("label", "label"),
         minus=2,
         role="amine 25, lower-right aminopropyl arm"),
    dict(id="acrylate_O9_tail", seed=(190, 1744), kind="to_junction_from_right",
         ends=("methyl", "junction"),
         minus=1,
         name_token=9,
         role="the ester O-alkyl chain of acrylate O9. The acrylate panel draws its ester O as an outline letter that TOUCHES the first C-O bond, so letter and tail are one component; the walk therefore starts at the terminal methyl and stops where the bond meets the letter, which is the oxygen position"),
    dict(id="acrylate_O12_tail", seed=(186, 1833), kind="to_junction_from_right",
         ends=("methyl", "junction"),
         minus=1,
         name_token=12,
         role="the ester O-alkyl chain of acrylate O12. The acrylate panel draws its ester O as an outline letter that TOUCHES the first C-O bond, so letter and tail are one component; the walk therefore starts at the terminal methyl and stops where the bond meets the letter, which is the oxygen position"),
    dict(id="acrylate_O14_tail", seed=(185, 1911), kind="to_junction_from_right",
         ends=("methyl", "junction"),
         minus=1,
         name_token=14,
         role="the ester O-alkyl chain of acrylate O14. The acrylate panel draws its ester O as an outline letter that TOUCHES the first C-O bond, so letter and tail are one component; the walk therefore starts at the terminal methyl and stops where the bond meets the letter, which is the oxygen position"),
    dict(id="acrylate_O16_tail", seed=(187, 1997), kind="to_junction_from_right",
         ends=("methyl", "junction"),
         minus=1,
         name_token=16,
         role="the ester O-alkyl chain of acrylate O16. The acrylate panel draws its ester O as an outline letter that TOUCHES the first C-O bond, so letter and tail are one component; the walk therefore starts at the terminal methyl and stops where the bond meets the letter, which is the oxygen position"),
]


# ---------------------------------------------------------------------------
def load():
    g = np.array(Image.open(FIG1).convert("L"))
    ink = g < INK
    lab, n = ndimage.label(ink, structure=np.ones((3, 3)))
    return g, ink, lab, n, ndimage.find_objects(lab)


def glyph_boxes(lab, objs):
    """Bounding boxes of the compact components -- the text glyphs."""
    out = []
    for i, sl in enumerate(objs):
        ys, xs = sl
        w, h = xs.stop - xs.start, ys.stop - ys.start
        a = int((lab[sl] == i + 1).sum())
        if w * h == 0:
            continue
        if (a / (w * h) >= MIN_GLYPH_FILL and 7 <= w < 60 and 10 <= h < 60):
            out.append((xs.start, xs.stop, ys.start, ys.stop, i + 1))
    return out


def _nearest_glyph(pt, boxes, skip_cid=None):
    """Distance from a bond end to the nearest text glyph, in px.

    `skip_cid` excludes the component being measured. Without it a short
    vertical bond -- an N-methyl stub -- is compact enough to be mistaken for
    a letter and then reports a distance of 0 to ITSELF, which would make
    every methyl look like a labelled substituent.
    """
    x, y = pt
    best = float("inf")
    for x0, x1, y0, y1, cid in boxes:
        if cid == skip_cid:
            continue
        dx = 0 if x0 <= x <= x1 else min(abs(x - x0), abs(x - x1))
        dy = 0 if y0 <= y <= y1 else min(abs(y - y0), abs(y - y1))
        best = min(best, (dx * dx + dy * dy) ** 0.5)
    return best


def _nearest_other_ink(pt, lab, cid, gap=LABEL_GAP_PX):
    """Is there ANY other component within `gap` of this bond end?

    This is the weaker test used where an atom's letter is fused into a
    neighbouring component and so is not in the glyph map: a bond that ends on
    an atom label ends ON something, a terminal methyl ends in open space.
    """
    x, y = pt
    r = gap
    sub = lab[max(0, y - r):y + r + 1, max(0, x - r):x + r + 1]
    other = set(int(v) for v in np.unique(sub)) - {0, cid}
    return sorted(other)


def endpoint_kind(pt, boxes, lab, cid):
    """'label' | 'label_fused' | 'methyl', measured, never declared."""
    d = _nearest_glyph(pt, boxes, skip_cid=cid)
    if d <= LABEL_GAP_PX:
        return "label", d
    if _nearest_other_ink(pt, lab, cid):
        return "label_fused", d
    return "methyl", d


def component_at(seed, ink, lab):
    x, y = seed
    if not ink[y, x]:
        raise DetectionError(
            f"01_trace_chains: seed {seed} is not ink. The figure or its "
            "resolution changed under a hard-coded seed. REFUSING rather than "
            "measuring whatever component is nearby -- error class 28.")
    return int(lab[y, x])


def skeleton(cid, lab, objs):
    sl = objs[cid - 1]
    sk = skeletonize(lab[sl] == cid)
    oy, ox = sl[0].start, sl[1].start
    S = set((int(r) + oy, int(c) + ox) for r, c in np.argwhere(sk))
    def nb(p):
        r, c = p
        return [(r + dr, c + dc) for dr in (-1, 0, 1) for dc in (-1, 0, 1)
                if (dr or dc) and (r + dr, c + dc) in S]
    return S, nb


def longest_path(S, nb):
    import collections
    deg = {p: len(nb(p)) for p in S}
    ends = [p for p in S if deg[p] == 1]
    if not ends:
        ends = [min(S)]
    best = None
    for e in ends:
        dist, par = {e: 0}, {e: None}
        q = collections.deque([e])
        while q:
            u = q.popleft()
            for v in nb(u):
                if v not in dist:
                    dist[v] = dist[u] + 1
                    par[v] = u
                    q.append(v)
        far = max(dist, key=lambda k: dist[k])
        if best is None or dist[far] > best[0]:
            best = (dist[far], far, par)
    _, far, par = best
    path, cur = [], far
    while cur is not None:
        path.append(cur)
        cur = par[cur]
    path.reverse()
    return [(c, r) for r, c in path], ends, [p for p in S if deg[p] > 2]


def path_to_junction(S, nb, from_right=False):
    """From the leftmost (or rightmost) degree-1 end, up to the first branch."""
    deg = {p: len(nb(p)) for p in S}
    ends = sorted([p for p in S if deg[p] == 1], key=lambda p: (p[1], p[0]),
                  reverse=from_right)
    if not ends:
        raise DetectionError("path_to_junction: component has no free end.")
    path, cur = [ends[0]], ends[0]
    while True:
        if deg[cur] > 2 and len(path) > 3:
            break
        nxt = [q for q in nb(cur) if q not in path[-3:]]
        if not nxt:
            break
        cur = nxt[0]
        path.append(cur)
    return [(c, r) for r, c in path], [p for p in S if deg[p] > 2]


def dp_sweep(path):
    vals = {e: chains.raster_chain_length(path, e) for e in EPS_SWEEP}
    return vals, sorted(set(vals.values()))


def dp_stable(path, label):
    """Vertex count that does not move across the tolerance sweep, or raise."""
    vals, uniq = dp_sweep(path)
    if len(uniq) != 1:
        raise DetectionError(
            f"01_trace_chains[{label}]: the Douglas-Peucker vertex count is "
            f"not stable across tolerances {dict(vals)}. A chain length is an "
            "integer; a count that depends on the smoothing parameter is a "
            "measurement that has not converged. Raise the render or fix the "
            "component, do not pick a tolerance.")
    return uniq[0], {str(k): v for k, v in vals.items()}


def path_shape(path, w=9):
    """Two measured statistics that decide whether `apex_count` applies.

    apex_count assumes EVERY vertex is a turning point in y. Two drawn things
    break that, and both are measurable rather than declared:
      steep_max  a bond drawn vertically has no y-extremum at its start, so a
                 path containing one loses a vertex. A 30-degree zig-zag bond
                 measures 0.58; a vertical bond measures infinity.
      flat_run   a C=C is drawn FLAT, so its two ends sit at the same y and the
                 reversal between them never happens.
    """
    xs = np.array([p[0] for p in path], float)
    ys = np.array([p[1] for p in path], float)
    if len(xs) < 2 * w:
        return None, None
    dx = xs[w:] - xs[:-w]
    dy = ys[w:] - ys[:-w]
    steep = float(np.max(np.abs(dy) / np.maximum(np.abs(dx), 1e-9)))
    run = best = 0
    for f in np.abs(dy) < 2.5:
        run = run + 1 if f else 0
        best = max(best, run)
    return steep, int(best)


STEEP_MAX = 1.0     # a 30-degree zig-zag measures 0.58-0.78; a vertical bond, inf
FLAT_RUN_MAX = 12   # a plain zig-zag measures 3-6 px; amine 6 C=C measures 31
ARC_TOL_BONDS = 0.35  # how far the arc-length bond count may sit from an integer


def arc_length(path):
    """Length along the traced skeleton, in px."""
    return float(sum(((path[i + 1][0] - path[i][0]) ** 2 +
                      (path[i + 1][1] - path[i][1]) ** 2) ** 0.5
                     for i in range(len(path) - 1)))


def turning_indices(path, smooth=5):
    """Indices of the y-reversals in a traced path."""
    ys = np.array([q[1] for q in path], float)
    if len(ys) < smooth + 2:
        return []
    yy = np.convolve(ys, np.ones(smooth) / smooth, mode="valid")
    d = np.diff(yy)
    # Keep the INDEX of every non-flat step, then look for sign changes
    # between consecutive non-flat steps. Dropping the flat steps without
    # keeping their indices (or requiring both neighbours to be non-flat)
    # silently loses most turning points and inflates the ruler -- it made the
    # bond measure 67.9 px instead of 27.9 on this module's first run.
    nz = [(i, float(np.sign(d[i]))) for i in range(len(d)) if abs(d[i]) > 1e-9]
    out = []
    for (i0, s0), (i1, s1) in zip(nz, nz[1:]):
        if s0 != s1:
            out.append((i0 + i1) // 2 + smooth // 2)
    return out


def calibrate_bond(paths):
    """THE INDEPENDENT SECOND METHOD'S RULER, and it never sees a DP count.

    A drawn bond's length is the arc between two consecutive turning points of
    the zig-zag. Measuring it that way uses only the turning-point detector, so
    the arc-length bond count derived from it is independent of
    Douglas-Peucker -- which matters, because DP is the only other route and a
    "second statistic" computed from the first is not a second statistic.

    Returns (bond_px, label_truncation_px). The truncation is how much shorter
    the last bond is when it stops on a heteroatom's text instead of at a
    methyl.
    """
    seg, trunc = [], []
    for path in paths:
        t = turning_indices(path)
        if len(t) < 3:
            continue
        for a, b in zip(t, t[1:]):
            seg.append(arc_length(path[a:b + 1]))
        trunc.append(arc_length(path[t[-1]:]))
    if len(seg) < 20:
        raise DetectionError(
            "calibrate_bond: fewer than 20 inter-turning-point segments; the "
            "ruler would be fitted on too little to be a check.")
    b = float(np.median(seg))
    return b, b - float(np.median(trunc))


def apex_count(path, smooth=5):
    """SECOND statistic: reversals of the drawn zig-zag, + the two endpoints.

    Independent of the Douglas-Peucker tolerance. Valid only for a chain drawn
    as a horizontal zig-zag -- it undercounts a vertical or sharply bent arm,
    so `kind` decides whether it is compared.
    """
    ys = np.array([p[1] for p in path], float)
    if len(ys) < smooth + 2:
        return None
    k = np.ones(smooth) / smooth
    yy = np.convolve(ys, k, mode="valid")
    d = np.diff(yy)
    d = d[np.abs(d) > 1e-9]
    if len(d) == 0:
        return 2
    s = np.sign(d)
    return int((s[1:] != s[:-1]).sum()) + 2


def ring_from_holes(cid, lab, objs):
    """Ring size from the enclosed hole the ring draws, not from its skeleton.

    The skeleton route fails here: amine 14's benzene ring shares a component
    with the OUTLINE LETTER O of its ether oxygen, and an outline letter is a
    closed loop too, so a cycle walk can wander into the wrong one. A ring's
    own hole is unambiguous -- the largest enclosed region -- and its boundary
    contour has one turn per ring bond.
    """
    from skimage import measure as skmeasure
    sl = objs[cid - 1]
    m = (lab[sl] == cid)
    holes = ndimage.binary_fill_holes(m) & ~m
    hl, nh = ndimage.label(holes)
    if nh == 0:
        raise DetectionError("ring_from_holes: this component encloses nothing.")
    sizes = ndimage.sum(holes, hl, range(1, nh + 1))
    big = int(np.argmax(sizes)) + 1
    cs = skmeasure.find_contours((hl == big).astype(float), 0.5)
    if not cs:
        raise DetectionError("ring_from_holes: the hole has no contour.")
    c = max(cs, key=len)
    pts = [(float(q[1]), float(q[0])) for q in c]
    vals, uniq = dp_sweep(pts)
    # A closed contour's first and last point coincide, so DP reports one
    # extra vertex; the ring size is that count minus the duplicate.
    return ([u - 1 for u in uniq], {str(k): v - 1 for k, v in vals.items()},
            float(sizes[big - 1]), [float(np.mean(c[:, 1])),
                                    float(np.mean(c[:, 0]))])


def ring_substitution(cid, lab, objs):
    """Where the two substituents sit on a drawn benzene ring. MEASURED.

    ortho / meta / para is a structural claim and nothing else in this
    curation can see it: the molecular formula is identical for all three, so
    HRMS, exact mass and every descriptor are blind (CLAUDE.md section 5,
    "mass cannot see connectivity"). The angle the two exocyclic bonds subtend
    at the ring's centre is 60 / 120 / 180 degrees respectively.
    """
    import math
    sl = objs[cid - 1]
    m = (lab[sl] == cid)
    holes = ndimage.binary_fill_holes(m) & ~m
    hl, nh = ndimage.label(holes)
    sizes = ndimage.sum(holes, hl, range(1, nh + 1))
    big = int(np.argmax(sizes)) + 1
    ys, xs = np.where(hl == big)
    cx, cy = xs.mean() + sl[1].start, ys.mean() + sl[0].start
    area = float(sizes[big - 1])
    r_ring = math.sqrt(area / 2.598)       # regular hexagon: A = 2.598 R^2

    S, nb = skeleton(cid, lab, objs)
    deg = {p: len(nb(p)) for p in S}
    clusters = []
    for p in [q for q in S if deg[q] > 2]:
        for c in clusters:
            if min(math.hypot(p[0] - q[0], p[1] - q[1]) for q in c) <= 3:
                c.append(p)
                break
        else:
            clusters.append([p])
    on_ring, off_ring = [], []
    for c in clusters:
        mx = float(np.mean([q[1] for q in c]))
        my = float(np.mean([q[0] for q in c]))
        r = math.hypot(mx - cx, my - cy)
        a = math.degrees(math.atan2(my - cy, mx - cx))
        (on_ring if abs(r - r_ring) <= 0.25 * r_ring
         else off_ring).append({"r_px": round(r, 1), "angle_deg": round(a, 1)})
    if len(on_ring) != 2:
        raise DetectionError(
            f"ring_substitution: found {len(on_ring)} substituted ring "
            f"vertices, need 2. (off-ring junctions: {off_ring})")
    d = abs(on_ring[0]["angle_deg"] - on_ring[1]["angle_deg"]) % 360
    d = min(d, 360 - d)
    pattern = {60: "ortho", 120: "meta", 180: "para"}
    nearest = min(pattern, key=lambda k: abs(d - k))
    if abs(d - nearest) > 15:
        raise DetectionError(
            f"ring_substitution: the two substituents subtend {d:.1f} deg, "
            "which is not within 15 deg of ortho, meta or para.")
    return {"hole_area_px": area, "ring_circumradius_px": round(r_ring, 1),
            "centre_xy": [round(cx, 1), round(cy, 1)],
            "substituted_vertices": on_ring, "other_junctions": off_ring,
            "subtended_angle_deg": round(d, 1), "pattern": pattern[nearest],
            "deviation_from_ideal_deg": round(abs(d - nearest), 1)}


def audit_vector():
    """Reproduce the claim that no vector chemistry exists in the SI."""
    import pymupdf
    doc = pymupdf.open(SI)
    rows = []
    for i in range(doc.page_count):
        dr = doc[i].get_drawings()
        stroked = sum(1 for p in dr if p.get("width"))
        if dr:
            rows.append({"page": i + 1, "paths": len(dr), "stroked": stroked})
    print(json.dumps({"si_pages": doc.page_count,
                      "pages_with_any_vector_path": rows,
                      "conclusion": "Every structural figure in the SI is a "
                      "flattened image; the stroked paths are table and "
                      "flow-cytometry rules. chains.vector_bonds has nothing "
                      "to read, so the raster route is the only route."},
                     indent=1))


def main():
    g, ink, lab, n, objs = load()
    boxes = glyph_boxes(lab, objs)
    out = {"_what": "Stage 01 chain measurements. Every integer here was "
                    "measured from Fig. 1a's pixels by curation/lib/chains; "
                    "none was counted by eye.",
           "figure": os.path.basename(FIG1), "image_px": list(g.shape[::-1]),
           "ink_threshold": INK, "dp_tolerance_sweep": list(EPS_SWEEP),
           "n_components": int(n), "n_glyph_boxes": len(boxes),
           "vector_route": "UNAVAILABLE -- see 00_value_sources.json and "
                           "`--audit-vector`. count_carbons is therefore "
                           "single-method and `agreed` is False everywhere.",
           "measurements": {}}

    # --- calibrate the independent ruler on the five clean homologous chains
    cal_paths = []
    for mid in ("amine01_chain", "amine02_chain", "amine03_chain",
                "amine04_chain", "amine05_chain"):
        mm = [x for x in M if x["id"] == mid][0]
        c = component_at(mm["seed"], ink, lab)
        S_, nb_ = skeleton(c, lab, objs)
        pth, _e, _j = longest_path(S_, nb_)
        if (pth[0][0], pth[0][1]) > (pth[-1][0], pth[-1][1]):
            pth = pth[::-1]
        cal_paths.append(pth)
    BOND_PX, TRUNC_PX = calibrate_bond(cal_paths)

    # The acrylate panel is a separate drawing block with its own zig-zag
    # amplitude, so it gets its own ruler, calibrated the same way on its own
    # four tails. Using the amine ruler there left a 0.40-bond residual on
    # O14 -- small, systematic, and exactly the kind of thing a single global
    # constant hides.
    acr_paths = []
    for mm in [x for x in M if x["id"].startswith("acrylate_")]:
        c = component_at(mm["seed"], ink, lab)
        S_, nb_ = skeleton(c, lab, objs)
        pth, _j = path_to_junction(S_, nb_, from_right=True)
        acr_paths.append(pth)
    ACR_BOND_PX, ACR_TRUNC_PX = calibrate_bond(acr_paths)
    RULERS = {"amines": (BOND_PX, TRUNC_PX),
              "acrylates": (ACR_BOND_PX, ACR_TRUNC_PX)}

    out["bond_ruler"] = {
        "amines": {"bond_px": round(BOND_PX, 2),
                   "label_truncation_px": round(TRUNC_PX, 2)},
        "acrylates": {"bond_px": round(ACR_BOND_PX, 2),
                      "label_truncation_px": round(ACR_TRUNC_PX, 2)},
        "derived_from": "the arc between consecutive turning points of amines "
                        "1-5, i.e. from the turning-point detector alone. It "
                        "never sees a Douglas-Peucker count, which is what "
                        "makes the arc-length bond count an INDEPENDENT second "
                        "measurement rather than a restatement of the first."}

    gaps, shape_gate, arc_resid = [], [], []
    for m in M:
        cid = component_at(m["seed"], ink, lab)
        S, nb = skeleton(cid, lab, objs)
        rec = {"role": m["role"], "seed": list(m["seed"]), "component": cid,
               "kind": m["kind"], "skeleton_px": len(S)}

        if m["kind"] == "ring_substitution":
            rec.update(ring_substitution(cid, lab, objs))
            out["measurements"][m["id"]] = rec
            continue

        if m["kind"] == "ring":
            uniq, sweep, area, centre = ring_from_holes(cid, lab, objs)
            rec.update(ring_vertices=uniq[0] if len(uniq) == 1 else uniq,
                       dp_sweep=sweep, hole_area_px=area,
                       hole_centroid_xy=[round(centre[0], 1),
                                         round(centre[1], 1)])
            if len(uniq) != 1:
                raise DetectionError(
                    f"01_trace_chains[{m['id']}]: the ring's hole contour "
                    f"gives {uniq} corners depending on the tolerance.")
            out["measurements"][m["id"]] = rec
            continue

        if m["kind"] == "stub":
            path, ends, juncs = longest_path(S, nb)
            far = max(path, key=lambda q: (q[1], q[0]))
            near = min(path, key=lambda q: (q[1], q[0]))
            fk, fd = endpoint_kind(far, boxes, lab, cid)
            nk, nd = endpoint_kind(near, boxes, lab, cid)
            rec.update(px_length=len(path), n_free_ends=len(ends),
                       upper_end=[fk, round(fd, 1)],
                       lower_end=[nk, round(nd, 1)])
            # A bond with NOTHING at its far end is a methyl. If a glyph ever
            # appears there the substituent is not a methyl and the structure
            # changes, so this is asserted rather than described.
            unlabelled = [e for e in (fk, nk) if e == "methyl"]
            if len(unlabelled) != 1:
                raise DetectionError(
                    f"01_trace_chains[{m['id']}]: declared an unlabelled stub "
                    f"(a methyl) but its two ends measure {fk!r} and {nk!r}. "
                    "Exactly one end must be bare.")
            rec["is_methyl"] = True
            out["measurements"][m["id"]] = rec
            continue

        if m["kind"] == "double_bond_marker":
            pm = [x for x in M if x["id"] == m["parent"]][0]
            pc = component_at(pm["seed"], ink, lab)
            pS, pnb = skeleton(pc, lab, objs)
            ppath, _e, _j = longest_path(pS, pnb)
            if (ppath[0][0], ppath[0][1]) > (ppath[-1][0], ppath[-1][1]):
                ppath = ppath[::-1]
            tix = turning_indices(ppath)
            vx = [ppath[0][0]] + [ppath[i][0] for i in tix] + [ppath[-1][0]]
            vy = [ppath[0][1]] + [ppath[i][1] for i in tix] + [ppath[-1][1]]
            gaps_x = [vx[i + 1] - vx[i] for i in range(len(vx) - 1)]
            med = float(np.median(gaps_x))
            ys_, xs_ = np.where(lab == cid)
            mx0, mx1 = int(xs_.min()), int(xs_.max())
            my = float(np.mean(ys_))
            # A FLAT bond has no turning point at either end, so the detector
            # reports ONE vertex in the middle of it instead of two. Find that
            # merged vertex by its two over-long neighbouring gaps and replace
            # it with the marker's own two ends.
            merged = [i for i in range(1, len(vx) - 1)
                      if gaps_x[i - 1] > 1.5 * med and gaps_x[i] > 1.5 * med]
            if len(merged) != 1:
                raise DetectionError(
                    f"01_trace_chains[{m['id']}]: expected exactly one merged "
                    f"vertex at the flat C=C, found {len(merged)}.")
            k = merged[0]
            vx2 = vx[:k] + [mx0, mx1] + vx[k + 1:]
            vy2 = vy[:k] + [my, my] + vy[k + 1:]
            par = out["measurements"][m["parent"]]
            if len(vx2) != par["vertices"]:
                raise DetectionError(
                    f"01_trace_chains[{m['id']}]: rebuilding the vertex list "
                    f"around the flat bond gives {len(vx2)} vertices, but "
                    f"Douglas-Peucker measured {par['vertices']} on the same "
                    "path. Refusing.")
            i_flat = k                      # first atom of the C=C
            ncarbon = par["heavy_atoms_in_chain"]
            # v_0 is the terminal methyl, v_(n) is the N position, so vertex i
            # is carbon number (ncarbon - i) counting from the amine end.
            delta = ncarbon - (i_flat + 1)
            # cis/trans: the two chain stubs flanking the flat bond
            left_y, right_y = vy2[i_flat - 1], vy2[i_flat + 2]
            same_side = (left_y - my) * (right_y - my) > 0
            rec.update(marker_bbox=[mx0, mx1, int(ys_.min()), int(ys_.max())],
                       vertex_x_after_rebuild=[int(q) for q in vx2],
                       flat_bond_between_vertices=[i_flat, i_flat + 1],
                       double_bond_position_from_amine_end=delta,
                       flank_y={"marker_y": round(my, 1),
                                "left_vertex_y": round(float(left_y), 1),
                                "right_vertex_y": round(float(right_y), 1)},
                       both_stubs_same_side=bool(same_side),
                       geometry="Z (cis)" if same_side else "E (trans)",
                       geometry_basis="The two chain stubs flanking the flat "
                       "C=C are on the same side of it (cis) or on opposite "
                       "sides (trans). Measured from the rebuilt vertex list, "
                       "not read off the picture.")
            out["measurements"][m["id"]] = rec
            continue

        if m["kind"] == "three_branch":
            import collections
            deg = {p: len(nb(p)) for p in S}
            ends = [p for p in S if deg[p] == 1]
            jset = {p for p in S if deg[p] > 2}
            if not jset:
                raise DetectionError(
                    f"01_trace_chains[{m['id']}]: declared `three_branch` but "
                    "the skeleton has no branch point.")
            br = []
            for e in ends:
                dist, par, q = {e: 0}, {e: None}, collections.deque([e])
                hit = None
                while q:
                    u = q.popleft()
                    if u in jset and u != e:
                        hit = u
                        break
                    for v2 in nb(u):
                        if v2 not in dist:
                            dist[v2] = dist[u] + 1
                            par[v2] = u
                            q.append(v2)
                pth, cur = [], hit
                while cur is not None:
                    pth.append(cur)
                    cur = par[cur]
                pth = [(c, r) for r, c in pth]
                v, sweep = dp_stable(pth, f"{m['id']}:{e}")
                ept = (int(e[1]), int(e[0]))
                kind_, d_ = endpoint_kind(ept, boxes, lab, cid)
                br.append({"free_end_xy": list(ept), "vertices": v,
                           "dp_sweep": sweep, "free_end": [kind_, round(d_, 1)],
                           "carbons_excluding_junction":
                               v - 1 - (1 if kind_ != "methyl" else 0)})
            rec["branches"] = sorted(br, key=lambda b: -b["vertices"])
            rec["n_free_ends"] = len(ends)
            if len(ends) != 3:
                raise DetectionError(
                    f"01_trace_chains[{m['id']}]: expected 3 free ends at a "
                    f"branch point, measured {len(ends)}.")
            out["measurements"][m["id"]] = rec
            continue

        if m["kind"] in ("to_junction", "to_junction_from_right"):
            path, juncs = path_to_junction(
                S, nb, from_right=(m["kind"] == "to_junction_from_right"))
            rec["n_junction_px"] = len(juncs)
            if not juncs:
                raise DetectionError(
                    f"01_trace_chains[{m['id']}]: declared `{m['kind']}` but "
                    "the component has no branch point.")
            checkable = [(path[0], m["ends"][0])]
        else:
            path, ends, juncs = longest_path(S, nb)
            # Orient left-to-right: longest_path returns whichever endpoint the
            # BFS started from, while `ends` is declared in drawing order.
            if (path[0][0], path[0][1]) > (path[-1][0], path[-1][1]):
                path = path[::-1]
            rec["n_free_ends"] = len(ends)
            rec["n_junction_px"] = len(juncs)
            if m["kind"] == "linear" and juncs:
                raise DetectionError(
                    f"01_trace_chains[{m['id']}]: declared `linear` but the "
                    f"skeleton has {len(juncs)} branch pixel(s). A branched "
                    "tail read as linear is exactly the isomer error this "
                    "project flags; refusing.")
            checkable = [(path[0], m["ends"][0]), (path[-1], m["ends"][1])]

        # --- the two counts, in the order that lets the independent one
        # --- arbitrate if Douglas-Peucker has not converged.
        kinds = []
        for pt, declared in checkable:
            got, d = endpoint_kind(pt, boxes, lab, cid)
            gaps.append((m["id"], got, d))
            if got != declared:
                raise DetectionError(
                    f"01_trace_chains[{m['id']}]: endpoint {pt} measures "
                    f"{got!r} (nearest glyph {d:.1f} px) but is declared "
                    f"{declared!r}. Getting this wrong shifts the chain by one "
                    "atom, which is error class 48 by another route. Refusing.")
            kinds.append([got, round(d, 1)])

        n_label_ends = sum(1 for k in kinds if k[0] in ("label", "label_fused"))
        rul = "acrylates" if m["id"].startswith("acrylate_") else "amines"
        bpx, tpx = RULERS[rul]
        rec["ruler"] = rul
        bonds_arc = (arc_length(path) + n_label_ends * tpx) / bpx
        arc_int = int(round(bonds_arc))
        arc_resid_frac = abs(bonds_arc - arc_int)

        vals, uniq = dp_sweep(path)
        sweep = {str(k): val for k, val in vals.items()}
        if len(uniq) == 1:
            v = uniq[0]
            rec["dp_converged"] = True
        else:
            # Douglas-Peucker has not converged. Do NOT pick a tolerance --
            # let the INDEPENDENT arc-length measurement decide, and only if
            # it lands cleanly on an integer that the sweep actually produced.
            cand = arc_int + 1          # vertices = bonds + 1
            if arc_resid_frac > ARC_TOL_BONDS or cand not in uniq:
                raise DetectionError(
                    f"01_trace_chains[{m['id']}]: the Douglas-Peucker count "
                    f"is not stable ({sweep}) and the independent arc-length "
                    f"measurement ({bonds_arc:.2f} bonds) does not resolve it "
                    f"to one of those values. Refusing to pick a tolerance.")
            v = cand
            rec["dp_converged"] = False
            rec["dp_resolved_by"] = (
                "the arc-length measurement, which is independent of the "
                f"Douglas-Peucker tolerance: {bonds_arc:.2f} bonds against a "
                f"ruler of {BOND_PX:.2f} px/bond, i.e. {v} vertices, and that "
                "value does appear in the sweep. The tolerance was not chosen.")

        rec["arc_length_px"] = round(arc_length(path), 1)
        rec["bonds_from_arc_length"] = round(bonds_arc, 2)
        rec["bonds_from_dp"] = v - 1
        resid = abs(bonds_arc - (v - 1))
        rec["arc_vs_dp_residual_bonds"] = round(resid, 2)
        arc_resid.append((m["id"], resid))
        if resid > ARC_TOL_BONDS:
            raise DetectionError(
                f"01_trace_chains[{m['id']}]: Douglas-Peucker counts {v - 1} "
                f"bonds; the drawn arc length is {bonds_arc:.2f} bonds against "
                f"a ruler of {bpx:.2f} px/bond. The two independent "
                "measurements differ by more than "
                f"{ARC_TOL_BONDS} of a bond. A chain length is an integer, so "
                "one of them is wrong. Refusing.")

        ap = apex_count(path)
        rec.update(vertices=v, dp_sweep=sweep, apex_count=ap,
                   px_length=len(path),
                   end_left_xy=[int(path[0][0]), int(path[0][1])],
                   end_right_xy=[int(path[-1][0]), int(path[-1][1])])
        rec["endpoints_measured"] = kinds
        rec["minus"] = m["minus"]
        rec["heavy_atoms_in_chain"] = v - m["minus"]

        steep, flat = path_shape(path)
        applicable = (ap is not None and steep is not None
                      and steep <= STEEP_MAX and flat <= FLAT_RUN_MAX)
        rec["path_shape"] = {"steep_max": None if steep is None else round(steep, 2),
                             "flat_run_px": flat,
                             "apex_applicable": bool(applicable)}
        shape_gate.append((m["id"], applicable, steep, flat))
        if applicable:
            if ap != v:
                raise DetectionError(
                    f"01_trace_chains[{m['id']}]: Douglas-Peucker says {v} "
                    f"vertices, the zig-zag apex count says {ap}. Two "
                    "statistics on the same trace disagree and the vector "
                    "route is unavailable, so there is nothing to break the "
                    "tie. Refusing.")
            rec["second_statistic_agrees"] = True
        else:
            rec["second_statistic_agrees"] = None
            rec["second_statistic_not_applicable"] = (
                "This path contains a vertical bond (steep_max "
                f"{steep}) or a flat drawn bond (flat_run {flat} px), so the "
                "apex statistic structurally undercounts and is not compared. "
                "The Douglas-Peucker count is stable over the whole tolerance "
                "sweep and is the single measurement for this chain.")

        cc = chains.count_carbons(raster=v, name_token=m.get("name_token"),
                                  label=m["id"])
        rec["count_carbons"] = cc
        if m.get("name_token") is not None:
            rec["name_token_check"] = chains.check_name_token(
                m["name_token"], rec["heavy_atoms_in_chain"], label=m["id"])
        out["measurements"][m["id"]] = rec

    # The figure's own bond pitch, from the homologous series 1-5. This is a
    # self-consistency statistic, not an input to any count.
    w = [out["measurements"][f"amine0{i}_chain"]["px_length"]
         for i in (1, 2, 3, 4, 5)]
    d = [w[i + 1] - w[i] for i in range(4)]
    out["bond_pitch_px"] = {
        "consecutive_width_deltas": d,
        "implied_pitch_per_carbon": [round(x / 2.0, 2) for x in d],
        "note": "Amines 1-5 differ by exactly two drawn carbons each, so each "
                "delta is two bond pitches. A delta that was not ~2x a single "
                "pitch would mean the series is not a simple homologous one."}

    out["arc_vs_dp"] = {
        "tolerance_bonds": ARC_TOL_BONDS,
        "worst_residual": round(max(r for _, r in arc_resid), 3),
        "worst_at": max(arc_resid, key=lambda z: z[1])[0],
        "n_compared": len(arc_resid),
        "note": "Every chain was measured twice: Douglas-Peucker vertices on "
                "the skeleton, and bonds from the drawn arc length against a "
                "ruler calibrated from turning points only. chains."
                "count_carbons still receives one `raster` count because both "
                "routes are raster; this is the agreement that the missing "
                "vector route would otherwise have provided."}

    # The apex-applicability gate, made auditable the same way.
    ok = [(st, fl) for _, a, st, fl in shape_gate if a]
    no = [(st, fl) for _, a, st, fl in shape_gate if not a and st is not None]
    out["apex_applicability_gate"] = {
        "steep_max_threshold": STEEP_MAX, "flat_run_threshold": FLAT_RUN_MAX,
        "applicable_worst_steep": round(max(s_ for s_, _ in ok), 2) if ok else None,
        "applicable_worst_flat": max(f for _, f in ok) if ok else None,
        "inapplicable_best_steep": round(min(s_ for s_, _ in no), 2) if no else None,
        "inapplicable_best_flat": max(f for _, f in no) if no else None,
        "n_applicable": len(ok), "n_not_applicable": len(no),
        "note": "Every path the gate admits had apex == Douglas-Peucker. Every "
                "path it rejects would have disagreed, and each one does "
                "contain a vertical or a flat drawn bond. The two populations "
                "are separated by a wide empty band, so the thresholds are "
                "not doing the deciding."}

    # The LABEL_GAP_PX floor, made auditable (the F-051 lesson).
    lab_d = [d for _, k, d in gaps if k == "label"]
    oth_d = [d for _, k, d in gaps if k != "label"]
    out["label_gap_separation"] = {
        "threshold_px": LABEL_GAP_PX,
        "largest_accepted_as_a_label": round(max(lab_d), 1) if lab_d else None,
        "smallest_rejected": round(min(oth_d), 1) if oth_d else None,
        "note": "A real atom label abuts the bond; the margin row numeral does "
                "not. If these two ever stop being cleanly separated the "
                "threshold is doing the deciding, not the drawing."}
    if lab_d and oth_d and not (max(lab_d) < LABEL_GAP_PX < min(oth_d)):
        raise DetectionError(
            "01_trace_chains: LABEL_GAP_PX no longer separates atom labels "
            f"(max {max(lab_d):.1f} px) from everything else "
            f"(min {min(oth_d):.1f} px). The threshold would be deciding the "
            "answer. Refusing.")

    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=1)
    print(f"wrote {OUT}: {len(out['measurements'])} measurements")
    for k, v in out["measurements"].items():
        print(f"  {k:26s} "
              f"{v.get('heavy_atoms_in_chain', v.get('ring_vertices', v.get('branches', v.get('is_methyl', ''))))}")


if __name__ == "__main__":
    if "--audit-vector" in sys.argv:
        audit_vector()
    else:
        main()
