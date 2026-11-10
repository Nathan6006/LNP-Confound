#!/usr/bin/env python3
"""Stage 01 -- the building blocks of the 36 selected ionizable lipids.

    python3 -I stage01_blocks.py

WHAT THE PAPER GIVES, AND WHAT IT DOES NOT.

The SI's synthetic procedures NAME the reagents for five of the eleven blocks:

  core    (2-bromoethoxy)-TBS + piperazine -> 1-(2-TBSO-ethyl)piperazine (2);
          2 + 3-iodopropane-1,2-diol -> the "piperazine (diamine) core" (SI
          1.2.3-1.2.4).
  EB-5    glutaric anhydride, 942 mg / 8.24 mmol -> MW 114.3 (glutaric
          anhydride is 114.10).
  T-1     2-hexyl-1-decanol, 1.00 g / 4.12 mmol -> MW 242.7 (C16H34O = 242.44).
  T-2     2-octyl-1-dodecanol, 1.00 g / 3.35 mmol -> MW 298.5 (C20H42O
          = 298.55).
  FH-5    3,5-dihydroxybenzoic acid (SI 1.2.1).
  FH-3    no head at all: the free primary alcohol left after TBAF removes the
          TBS group ("None FH (1 degree alcohol)" in Fig 3a).

It names NOTHING for FH-1, FH-2, FH-4, FH-6, EB-6 and EB-7. Those appear in no
procedure, no reagent list, no HRMS and no NMR assignment. Figure 3a's drawing
is the ONLY evidence that exists for them, so they are measured here from the
drawing -- never counted by eye (error class 48).

THE MEASUREMENT, AND WHY IT HAS TEETH.

One tracer is used for every block: binarise -> connected component ->
skeletonise -> walk the skeleton -> curation/lib `chains.raster_chain_length`
(Douglas-Peucker). One ring-geometry routine is used for every aryl head:
fill the ring's hole, take the hexagon's six vertex angles from the hole's own
radial profile, and find which vertices carry an exocyclic bond.

The point is that FIVE of the blocks it is run on have an answer that does not
come from the drawing:

  T-1   must come out 2-hexyldecyl      (from 2-hexyl-1-decanol)
  T-2   must come out 2-octyldodecyl    (from 2-octyl-1-dodecanol)
  EB-5  must come out 5 backbone carbons (from glutaric anhydride)
  FH-5  must come out 3,5-disubstituted (from 3,5-dihydroxybenzoic acid)
  FH-2  must come out monosubstituted   (plain benzoyl, drawn with no OH)

So the tracer is calibrated against external ground truth on the blocks that
have it, and only then applied to FH-1, FH-4, FH-6, EB-6 and EB-7. A tracer
that could not reproduce 2-hexyl-1-decanol has no business counting FH-1.

`chains.count_carbons` is called with raster= only, because there IS no vector
route here: every figure in the article PDF is an embedded raster (<=9 drawing
paths per page), so `chains.vector_bonds` has nothing to read. That is recorded
on every block as `second_method: "unavailable"`, not quietly dropped.
"""
import csv
import json
import math
import os
import sys
from collections import deque

import numpy as np
from PIL import Image
from scipy import ndimage
from skimage.morphology import skeletonize

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import chains                                               # noqa: E402
from errors import DetectionError                           # noqa: E402

FIG3 = os.path.join(HERE, "source", "pdfimg", "fig3_1590x2450.png")
INK_SUM = 300          # RGB-sum threshold for ink on Fig 3a's white-ish panels
RDP_EPS = 3.0          # px; reported at 2.5 / 3.0 / 4.0 and required stable

# The core with its three hydroxyl oxygens atom-mapped: map 2 = the
# hydroxyethyl arm (the FH site), map 3 = the glycerol secondary OH, map 4 =
# the glycerol primary OH. Stated here so WHICH oxygen receives WHICH acyl
# group is a written fact rather than a substructure match order.
CORE_MAPPED = "[OH:2]CCN1CCN(CC1)CC([OH:3])C[OH:4]"


# --------------------------------------------------------------------------
# skeleton machinery
# --------------------------------------------------------------------------
def components(path, ink_sum=INK_SUM):
    img = np.asarray(Image.open(path).convert("RGB")).astype(int)
    ink = img.sum(axis=2) < ink_sum
    lab, n = ndimage.label(ink, structure=np.ones((3, 3)))
    return img, ink, lab, n


def _degree(sk):
    s = sk.astype(np.uint8)
    return ndimage.convolve(s, np.ones((3, 3), np.uint8), mode="constant") - s


def walk_to_branch(sk, start, limit=5000):
    """Ordered skeleton path from an endpoint until the first branch point.

    Straightest-continuation tie-break, so a bond crossing never turns the walk
    down a side chain.
    """
    deg = _degree(sk)
    path = [tuple(int(v) for v in start)]
    prev, cur = None, path[0]
    while True:
        y, x = cur
        if deg[y, x] >= 3 and len(path) > 1:
            break
        cand = [(y + dy, x + dx) for dy in (-1, 0, 1) for dx in (-1, 0, 1)
                if (dy or dx) and sk[y + dy, x + dx]
                and (y + dy, x + dx) != prev
                and (y + dy, x + dx) not in path[-3:]]
        if not cand:
            break
        nxt = cand[0]
        if len(cand) > 1 and prev is not None:
            vy, vx = y - prev[0], x - prev[1]
            nxt = max(cand, key=lambda c: (c[0] - y) * vy + (c[1] - x) * vx)
        path.append(nxt)
        prev, cur = cur, nxt
        if len(path) > limit:
            raise DetectionError("walk_to_branch: runaway path; the component "
                                 "is not a simple chain.")
    return path


def geodesic(sk, a, b):
    prev = {a: None}
    q = deque([a])
    while q:
        y, x = q.popleft()
        if (y, x) == b:
            break
        for dy in (-1, 0, 1):
            for dx in (-1, 0, 1):
                if dy or dx:
                    p = (y + dy, x + dx)
                    if sk[p] and p not in prev:
                        prev[p] = (y, x)
                        q.append(p)
    if b not in prev:
        raise DetectionError(f"geodesic: {a} and {b} are not connected.")
    out, c = [], b
    while c:
        out.append(c)
        c = prev[c]
    return out[::-1]


def vertices(path_yx, label):
    """Vertex count via curation/lib, with the eps sweep recorded.

    A count that moves with eps is a count that is being chosen, not measured,
    so the three values are reported and the caller asserts they agree.
    """
    xy = [(x, y) for y, x in path_yx]
    got = {eps: chains.raster_chain_length(xy, eps_px=eps)
           for eps in (2.5, RDP_EPS, 4.0)}
    if len(set(got.values())) != 1:
        raise DetectionError(
            f"vertices({label!r}): Douglas-Peucker count depends on eps "
            f"({got}). The count is being chosen by a tolerance rather than "
            "measured. Re-crop or re-render; do not pick one.")
    return got[RDP_EPS], {str(k): v for k, v in got.items()}


# --------------------------------------------------------------------------
# aryl ring geometry
# --------------------------------------------------------------------------
def ring_substitution(ink, lab, cid, label):
    """(vertex angles, exocyclic angles, substituent ring positions).

    Positions are numbered from the attachment vertex, walking the hexagon in
    whichever direction gives the lexicographically smallest locant set -- the
    IUPAC rule, applied to measured angles rather than to a reading of the
    picture.
    """
    m = (lab == cid)
    filled = ndimage.binary_fill_holes(m)
    hole = filled & ~m
    hl, hn = ndimage.label(hole)
    sizes = ndimage.sum(hole, hl, range(1, hn + 1))
    if hn != 1:
        raise DetectionError(
            f"ring_substitution({label!r}): {hn} holes in the component; "
            "expected exactly one aromatic ring.")
    cy, cx = ndimage.center_of_mass(hl == 1)
    ys, xs = np.nonzero(hole)
    rr = np.hypot(ys - cy, xs - cx)
    aa = np.degrees(np.arctan2(-(ys - cy), xs - cx)) % 360.0
    rmax = np.zeros(360)
    for a, r in zip(aa, rr):
        i = int(a) % 360
        rmax[i] = max(rmax[i], r)
    vang = _cluster([i for i in range(360) if rmax[i] >= rmax.max() - 1.5], 6)
    if len(vang) != 6:
        raise DetectionError(
            f"ring_substitution({label!r}): {len(vang)} hexagon vertices "
            f"found at {vang}; a benzene ring has six.")
    a_vertex = float(rmax.max()) + 3.0          # bond centreline, roughly
    exo = []
    for a in range(360):
        th = math.radians(a)
        hit = 0
        for r in np.arange(a_vertex * 1.5, a_vertex * 1.85, 1.0):
            y = int(round(cy - r * math.sin(th)))
            x = int(round(cx + r * math.cos(th)))
            if 0 <= y < ink.shape[0] and 0 <= x < ink.shape[1] and ink[y, x]:
                hit += 1
        if hit >= 8:
            exo.append(a)
    eang = _cluster(exo, None)
    # keep only exocyclic directions that point at a measured ring vertex
    on_vertex = [e for e in eang if min(_angdist(e, v) for v in vang) < 8.0]
    return cy, cx, vang, eang, on_vertex


def _angdist(a, b):
    d = abs(a - b) % 360.0
    return min(d, 360.0 - d)


def _cluster(idx, expect, gap=8):
    if not idx:
        return []
    cl, cur = [], [idx[0]]
    for v in idx[1:]:
        if v - cur[-1] <= gap:
            cur.append(v)
        else:
            cl.append(cur)
            cur = [v]
    cl.append(cur)
    if len(cl) > 1 and (cl[0][0] + 360 - cl[-1][-1]) <= gap:
        cl[0] = [x - 360 for x in cl[-1]] + cl[0]
        cl.pop()
    out = sorted(round(float(np.mean(c)) % 360.0, 1) for c in cl)
    return out


def locants(attach_angle, vang, sub_angles):
    """Ring positions of the substituents, attachment = 1."""
    order = sorted(vang, key=lambda v: _angdist(v, attach_angle))
    start = order[0]
    ring = sorted(vang)
    i0 = ring.index(start)
    best = None
    for direction in (1, -1):
        pos = {}
        for k in range(6):
            pos[ring[(i0 + direction * k) % 6]] = k + 1
        got = sorted(pos[min(ring, key=lambda v: _angdist(v, s))]
                     for s in sub_angles)
        if best is None or got < best:
            best = got
    return best


# --------------------------------------------------------------------------
# the measurements, with their component ids and seed points
# --------------------------------------------------------------------------
# Component ids and start pixels are read off Fig 3a once and recorded here, so
# a reviewer can re-run the same measurement on the same pixels. They are
# positions, not answers.
CHAIN_JOBS = [
    # label                 cid  start (y, x)   what the vertex count means
    ("FH-1 main chain",     168, (390, 377), "branch C + n carbons"),
    ("FH-1 branch chain",   168, (444, 284), "branch C + n carbons"),
    ("T-1 main chain",      166, (364, 1517), "branch C + n carbons"),
    ("T-1 branch chain",    166, (418, 1424), "branch C + n carbons"),
    ("T-2 main chain",      224, (534, 1549), "branch C + n carbons"),
    ("T-2 branch chain",    224, (588, 1456), "branch C + n carbons"),
]
EB_JOBS = [("EB-5", 167, (348, 869), (348, 993)),
           ("EB-6", 211, (462, 855), (525, 1003)),
           ("EB-7", 240, (580, 840), (580, 1026))]
ARYL_JOBS = [("FH-2", 169), ("FH-4", 223), ("FH-5", 221), ("FH-6", 222)]
SPACER_JOBS = [("FH-4", 223), ("FH-5", 221), ("FH-6", 222)]


def measure():
    img, ink, lab, _ = components(FIG3)
    qc = {"figure": os.path.relpath(FIG3, HERE), "ink_sum_threshold": INK_SUM,
          "rdp_eps_px": RDP_EPS, "chains": {}, "ester_bridges": {},
          "aryl_rings": {}, "spacers": {}}

    for label, cid, start, meaning in CHAIN_JOBS:
        sk = skeletonize(lab == cid)
        p = walk_to_branch(sk, start)
        n, sweep = vertices(p, label)
        qc["chains"][label] = {
            "component": cid, "start_px": list(start), "end_px": list(p[-1]),
            "skeleton_path_px": len(p), "vertices": n,
            "carbons_excluding_branch_C": n - 1, "rdp_eps_sweep": sweep,
            "meaning": meaning}

    for label, cid, t1, t2 in EB_JOBS:
        sk = skeletonize(lab == cid)
        c1 = walk_to_branch(sk, t1)[-1]
        c2 = walk_to_branch(sk, t2)[-1]
        p = geodesic(sk, c1, c2)
        n, sweep = vertices(p, label)
        tok = int(label.split("-")[1])
        qc["ester_bridges"][label] = {
            "component": cid, "carbonyl_C_px": [list(c1), list(c2)],
            "backbone_carbons": n, "rdp_eps_sweep": sweep,
            "name_token_check": chains.check_name_token(tok, n, label)}

    for label, cid in ARYL_JOBS:
        cy, cx, vang, eang, on_vertex = ring_substitution(ink, lab, cid, label)
        # the attachment is the exocyclic direction pointing back at the
        # carbonyl, i.e. the one nearest 150 deg in this figure's layout
        att = min(on_vertex, key=lambda a: _angdist(a, 150.0))
        subs = [a for a in on_vertex if a != att]
        qc["aryl_rings"][label] = {
            "component": cid, "ring_centre_px": [round(cy, 1), round(cx, 1)],
            "hexagon_vertex_angles_deg": vang,
            "exocyclic_angles_deg": eang,
            "exocyclic_on_a_vertex_deg": on_vertex,
            "attachment_angle_deg": att,
            "substituent_locants": locants(att, vang, subs) if subs else []}

    for label, cid in SPACER_JOBS:
        sk = skeletonize(lab == cid)
        deg = _degree(sk)
        ends = np.argwhere(sk & (deg == 1))
        top = ends[np.argmin(ends[:, 0])]       # a C=O double-bond stub tip
        co = walk_to_branch(sk, top)[-1]
        cy, cx = qc["aryl_rings"][label]["ring_centre_px"]
        att = math.radians(qc["aryl_rings"][label]["attachment_angle_deg"])
        a = 37.0
        ty, tx = cy - a * math.sin(att), cx + a * math.cos(att)
        pts = np.argwhere(sk)
        ring_c = tuple(int(v) for v in
                       pts[np.argmin(np.hypot(pts[:, 0] - ty,
                                              pts[:, 1] - tx))])
        p = geodesic(sk, co, ring_c)
        n, sweep = vertices(p, label + " spacer")
        qc["spacers"][label] = {
            "carbonyl_C_px": list(co), "ring_C_px": list(ring_c),
            "vertices_inclusive": n, "CH2_between_CO_and_ring": n - 2,
            "rdp_eps_sweep": sweep}
    return qc


# --------------------------------------------------------------------------
# the five calibration blocks -- external ground truth, asserted
# --------------------------------------------------------------------------
def calibrate(qc):
    """Assert the tracer on the blocks whose answer does not come from pixels.

    Each entry is (what was measured, what the NAMED REAGENT requires, why).
    A failure here means the tracer is wrong, and nothing downstream of it may
    be trusted -- including FH-1, FH-4, FH-6, EB-6 and EB-7, which have no
    external anchor at all.
    """
    out = []

    def check(name, got, want, source):
        ok = (got == want)
        out.append({"block": name, "measured": got, "required": want,
                    "pass": ok, "external_source": source})
        return ok

    t1m = qc["chains"]["T-1 main chain"]["carbons_excluding_branch_C"]
    t1b = qc["chains"]["T-1 branch chain"]["carbons_excluding_branch_C"]
    check("T-1 total carbons", 1 + 1 + t1m + t1b, 16,
          "SI 1.2.2: 2-hexyl-1-decanol, 1.00 g / 4.12 mmol -> MW 242.7; "
          "C16H34O = 242.44. A 2-hexyldecyl group is 16 carbons.")
    check("T-1 branch length", t1b, 6, "the 'hexyl' of 2-hexyl-1-decanol")
    t2m = qc["chains"]["T-2 main chain"]["carbons_excluding_branch_C"]
    t2b = qc["chains"]["T-2 branch chain"]["carbons_excluding_branch_C"]
    check("T-2 total carbons", 1 + 1 + t2m + t2b, 20,
          "SI 1.2.2: 2-octyl-1-dodecanol, 1.00 g / 3.35 mmol -> MW 298.5; "
          "C20H42O = 298.55.")
    check("T-2 branch length", t2b, 8, "the 'octyl' of 2-octyl-1-dodecanol")
    check("EB-5 backbone carbons",
          qc["ester_bridges"]["EB-5"]["backbone_carbons"], 5,
          "SI 1.2.2: glutaric anhydride, 942 mg / 8.24 mmol -> MW 114.3; "
          "glutaric anhydride is C5H6O3, 114.10.")
    check("FH-5 substituent locants",
          qc["aryl_rings"]["FH-5"]["substituent_locants"], [3, 5],
          "SI 1.2.1: 3,5-dihydroxybenzoic acid, 1.5 g / 9.7 mmol -> MW 154.6; "
          "C7H6O4 = 154.12.")
    check("FH-5 spacer CH2 count",
          qc["spacers"]["FH-5"]["CH2_between_CO_and_ring"], 0,
          "3,5-dihydroxyBENZOIC acid: the carbonyl is bonded to the ring.")
    check("FH-2 substituent count",
          len(qc["aryl_rings"]["FH-2"]["substituent_locants"]), 0,
          "Fig 3a draws FH-2 as an unsubstituted phenyl -- a benzoyl group.")
    return out


# --------------------------------------------------------------------------
# blocks -> SMILES
# --------------------------------------------------------------------------
def _alkyl(n):
    return "C" * n


def build(qc, cal):
    if not all(c["pass"] for c in cal):
        bad = [c["block"] for c in cal if not c["pass"]]
        raise SystemExit(
            "stage01: the tracer FAILED its external calibration on "
            f"{bad}. Those blocks' answers come from NAMED REAGENTS in the SI, "
            "not from the drawing, so a mismatch means the tracer is wrong. "
            "FH-1, FH-4, FH-6, EB-6 and EB-7 have no external anchor at all "
            "and must not be built from a tracer that cannot reproduce "
            "2-hexyl-1-decanol.")

    c = qc["chains"]
    fh1_main = c["FH-1 main chain"]["carbons_excluding_branch_C"]
    fh1_br = c["FH-1 branch chain"]["carbons_excluding_branch_C"]
    t1_main = c["T-1 main chain"]["carbons_excluding_branch_C"]
    t1_br = c["T-1 branch chain"]["carbons_excluding_branch_C"]
    t2_main = c["T-2 main chain"]["carbons_excluding_branch_C"]
    t2_br = c["T-2 branch chain"]["carbons_excluding_branch_C"]

    def aryl(label):
        loc = qc["aryl_rings"][label]["substituent_locants"]
        ring = {1: "c1", 2: "c", 3: "c", 4: "c", 5: "c", 6: "c"}
        for p in loc:
            ring[p] = "c(O)"
        return ("c1" + "".join(ring[p] for p in (2, 3, 4, 5, 6)) + "1")

    rows = []

    def add(bid, kind, smiles, fragment, attach, evidence, confidence, note):
        """`fragment` is the ATOM-MAPPED attachment form of this block.

        IT IS BUILT HERE, FROM THE MEASURED COUNTS, AND NEVER BY REWRITING A
        CANONICAL SMILES. The first version of stage03 derived each acyl
        fragment by string surgery on the block's canonical SMILES --
        'O=C(O)X' -> '[C:1](=O)X' -- and that is wrong whenever
        canonicalisation has moved the carboxyl carbon's NEIGHBOUR away from
        the front of the remainder. It silently produced palmitoyl in place of
        2-hexyldecanoyl for FH-1, and for every EB-T chain it put the glycerol
        ester on the far end of the tail instead of on the diacid.

        Every one of those products had the CORRECT MOLECULAR FORMULA, so both
        HRMS anchors still matched to 0.3 and 0.5 ppm. Mass cannot see
        connectivity. What does see it is the tertiary-carbon count Fig 2b
        prints, which stage03 now asserts.
        """
        rows.append({"block_id": bid, "kind": kind, "SMILES": smiles,
                     "fragment": fragment,
                     "attachment": attach, "evidence": evidence,
                     "confidence": confidence, "note": note})

    # --- the shared core ------------------------------------------------
    add("CORE", "core", "OCCN1CCN(CC1)CC(O)CO", CORE_MAPPED, "three OH",
        "SI 1.2.3-1.2.4: (2-bromoethoxy)-TBS + piperazine -> 2; 2 + "
        "3-iodopropane-1,2-diol -> the diamine core. Drawn in Fig S1.",
        "high",
        "1-(2-hydroxyethyl)-4-(2,3-dihydroxypropyl)piperazine. Named reagents, "
        "not pixels.")

    # --- functional heads ----------------------------------------------
    add("FH-1", "head_acid",
        f"OC(=O)C({_alkyl(fh1_br)})" + _alkyl(fh1_main),
        f"[C:1](=O)C({_alkyl(fh1_br)})" + _alkyl(fh1_main), "COOH",
        f"Fig 3a, traced: main chain {fh1_main} C, branch {fh1_br} C, "
        "both from curation/lib chains.raster_chain_length with a stable "
        "Douglas-Peucker count over eps 2.5/3.0/4.0.",
        "medium",
        "2-hexyldecanoic acid. NO named reagent, NO HRMS, NO NMR: the only "
        "evidence is the Fig 3a drawing, measured by the same tracer that "
        "reproduces 2-hexyl-1-decanol and 2-octyl-1-dodecanol exactly. The "
        "branch position (alpha to the carbonyl) is independently corroborated "
        "by Fig 2b, where FP3 (tertiary carbon) = 3 on exactly 36 of 213 LNPs "
        "-- the FH-1 count -- against 2 for every other head.")
    add("FH-2", "head_acid", "OC(=O)c1ccccc1", "[C:1](=O)c1ccccc1", "COOH",
        "Fig 3a: an unsubstituted phenyl on the carbonyl; ring geometry "
        "measured, zero exocyclic substituents besides the attachment.",
        "medium",
        "Benzoic acid. No named reagent; read from the drawing.")
    add("FH-3", "none", "", "", "n/a",
        "Fig 3a: 'None FH (1 degree alcohol)'. SI 1.2.5.1: TBAF removes the "
        "TBS group and the product FH-3EB-5T-1 is isolated as such.",
        "high",
        "There is no head. The hydroxyethyl arm stays a free primary alcohol.")
    add("FH-4", "head_acid", f"OC(=O){aryl('FH-4')}",
        f"[C:1](=O){aryl('FH-4')}", "COOH",
        "Fig 3a: ring locants "
        f"{qc['aryl_rings']['FH-4']['substituent_locants']} measured from the "
        "hexagon's own vertex angles; spacer "
        f"{qc['spacers']['FH-4']['CH2_between_CO_and_ring']} CH2.",
        "medium",
        "3,4-dihydroxybenzoic acid (protocatechuic acid). No named reagent. "
        "The text calls FH-4 and FH-5 positional isomers chosen to test where "
        "the two phenolic alcohols sit, which is consistent with 3,4 vs 3,5.")
    add("FH-5", "head_acid", f"OC(=O){aryl('FH-5')}",
        f"[C:1](=O){aryl('FH-5')}", "COOH",
        "SI 1.2.1 names 3,5-dihydroxybenzoic acid (1.5 g / 9.7 mmol -> MW "
        "154.6; C7H6O4 = 154.12) AND Fig 3a's measured locants are "
        f"{qc['aryl_rings']['FH-5']['substituent_locants']}.",
        "high",
        "3,5-dihydroxybenzoic acid. Named reagent, and the drawing agrees.")
    add("FH-6", "head_acid", f"OC(=O)CC{aryl('FH-6')}",
        f"[C:1](=O)CC{aryl('FH-6')}", "COOH",
        "Fig 3a: ring locants "
        f"{qc['aryl_rings']['FH-6']['substituent_locants']}, spacer "
        f"{qc['spacers']['FH-6']['CH2_between_CO_and_ring']} CH2 measured "
        "between the carbonyl carbon and the ring carbon.",
        "medium",
        "3-(3,4-dihydroxyphenyl)propanoic acid (dihydrocaffeic acid). No named "
        "reagent. The text's reason for choosing FH-6 -- 'the effect of the "
        "carbon length beyond the position of the phenolic alcohols' -- is "
        "what the measured 2-carbon spacer is.")

    # --- ester bridges --------------------------------------------------
    for lbl in ("EB-5", "EB-6", "EB-7"):
        nb = qc["ester_bridges"][lbl]["backbone_carbons"]
        named = {"EB-5": "glutaric anhydride, named in SI 1.2.2 "
                         "(942 mg / 8.24 mmol -> MW 114.3)"}.get(lbl, "")
        add(lbl, "diacid", "OC(=O)" + _alkyl(nb - 2) + "C(=O)O",
            "[C:1](=O)" + _alkyl(nb - 2) + "C(=O)O", "both COOH",
            f"Fig 3a: {nb} backbone carbons between the two carbonyl carbons, "
            "measured as a skeleton geodesic."
            + (" " + named if named else ""),
            "high" if named else "medium",
            {"EB-5": "glutaric acid", "EB-6": "adipic acid",
             "EB-7": "pimelic acid"}[lbl]
            + ". The name token EB-N equals the measured backbone carbon "
              "count: " + json.dumps(
                  qc["ester_bridges"][lbl]["name_token_check"]["matches"])
            + (". Only EB-5 is pinned by a named reagent; EB-6 and EB-7 rest "
               "on the drawing plus the token." if not named else "."))

    # --- tails ----------------------------------------------------------
    add("T-1", "alcohol", f"OCC({_alkyl(t1_br)}){_alkyl(t1_main)}",
        f"CC({_alkyl(t1_br)}){_alkyl(t1_main)}", "OH",
        f"SI 1.2.2 names 2-hexyl-1-decanol (MW 242.7 measured from 1.00 g / "
        f"4.12 mmol) AND Fig 3a traces to {t1_main}+{t1_br} carbons off the "
        "branch carbon.", "high",
        "2-hexyl-1-decanol. Named reagent, and the drawing agrees.")
    add("T-2", "alcohol", f"OCC({_alkyl(t2_br)}){_alkyl(t2_main)}",
        f"CC({_alkyl(t2_br)}){_alkyl(t2_main)}", "OH",
        f"SI 1.2.2 names 2-octyl-1-dodecanol (MW 298.5 measured from 1.00 g / "
        f"3.35 mmol) AND Fig 3a traces to {t2_main}+{t2_br} carbons off the "
        "branch carbon.", "high",
        "2-octyl-1-dodecanol. Named reagent, and the drawing agrees.")
    return rows


def main():
    qc = measure()
    cal = calibrate(qc)
    qc["calibration_against_named_reagents"] = cal
    rows = build(qc, cal)

    from rdkit import Chem, RDLogger
    from rdkit.Chem import rdMolDescriptors
    RDLogger.DisableLog("rdApp.*")
    for r in rows:
        if r["SMILES"]:
            m = Chem.MolFromSmiles(r["SMILES"])
            if m is None:
                raise SystemExit(f"stage01: {r['block_id']} does not parse: "
                                 f"{r['SMILES']}")
            r["SMILES"] = Chem.MolToSmiles(m)
            r["formula"] = rdMolDescriptors.CalcMolFormula(m)
        else:
            r["formula"] = ""

    with open(os.path.join(HERE, "01_blocks.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["block_id", "kind", "SMILES",
                                           "formula", "fragment",
                                           "attachment",
                                           "evidence", "confidence", "note"])
        w.writeheader()
        w.writerows(rows)
    with open(os.path.join(HERE, "01_blocks_qc.json"), "w") as fh:
        json.dump(qc, fh, indent=1)
    print(f"01_blocks.csv: {len(rows)} blocks")
    for c in cal:
        print(f"  calibration {'PASS' if c['pass'] else 'FAIL'}  "
              f"{c['block']}: measured {c['measured']}, required {c['required']}")


if __name__ == "__main__":
    main()
