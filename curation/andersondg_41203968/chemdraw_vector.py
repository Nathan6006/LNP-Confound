#!/usr/bin/env python3
"""This paper's drawings: ChemDraw vector art whose bonds are FILLED paths.

WHY THIS IS NOT curation/lib/chains.py. `chains.vector_bonds` selects stroked
paths and documents filled paths as "part of a letter outline -- label
evidence, never a bond". That is true of the PDF it was promoted from
(ng_2025). It is FALSE of this one: every bond in 41565_2025_2044_MOESM1 is a
filled polygon and `width` is None on all 13,271 of them, so
`chains.vector_bonds` returns zero segments on every page here. Measured, not
assumed. This module is the filled-path counterpart; `chains.count_carbons`
and `chains.check_name_token` are still used, on the segments this produces.

WHAT THE GEOMETRY ACTUALLY LOOKS LIKE, measured over SI pages 45-121:
  * 13,271 filled paths. Thickness is cleanly bimodal -- 12,723 at 0.3-0.6 pt
    (plain bonds) and 402 at 1.0-2.3 pt (stereo wedges, tapered quads). NOTHING
    between 0.6 and 1.0, so one threshold separates them with no judgement.
  * Bond spans cluster at 5, 7 and 10 pt: a bond drawn full length is ~10 pt
    and one TRIMMED at a heteroatom label is shorter. The shortest real bond
    found is 2.82 pt -- a wedge trimmed at a nitrogen -- which is why BOND_MIN
    is low and why the m/z gate below, not a span floor, is what protects the
    output.
  * Above 20 pt there is nothing until 34: those are page rules and arrows.

THE SNAPPING RULE, AND THE BUG THAT SET IT. ChemDraw trims a bond at the
label's box, so the drawn endpoint is short of the atom -- but the bond's own
line still POINTS AT the label centre. Snapping to the nearest label box broke
a 1,4,7-triazonane into an aziridine: two ring nitrogens 6 pt apart have
overlapping padded boxes and the wrong one was 0.5 pt nearer. Ray geometry is
not a tie-break on that, it is a different measurement. The box fallback is
kept for bonds too short for their direction to be reliable, but it must then
produce a plausible bond length -- without that it invented a 17.1 pt bond to
a nitrogen 3 pt from a wedge's narrow end, gave that N a valence of four, and
split AMG61 in two.

WHAT CHECKS IT. Nothing in this module. Every structure it produces is gated
in stage01 against the SI's own HRMS m/z for that compound (139 of them).
Mass cannot see connectivity -- but it does see the HYDROGEN COUNT, and that
is what distinguishes a monocyclic headgroup from a spiro one. It earned that
on its first run: with wedges excluded, AMG510's headgroup reconstructed as
4-(aminomethyl)piperidine and the formula came out C64H122N2O10; with the
wedges in it is the spiro diamine and C66H124N2O10, which is what the paper
reports.
"""
import math, re
import pymupdf

# Lengths that scale with the DRAWING are given as multiples of the page's
# own median bond span. The Am15 series is set at 5 pt type, ~0.7x the rest of
# the SI, and absolute tolerances there merged distinct oxygens into one atom.
TOL_F     = 0.17    # x scale: vertices closer than this are one atom
LBL_PAD_F = 0.30    # x scale: pad a glyph box when snapping a short bond
RAY_PERP_F= 0.22    # x scale: how far off the bond's ray a label centre may sit
DBL_PERP  = 3.6     # pt: parallel offset of a double-bond inner line
BOND_MIN  = 1.5     # pt: the shortest real bond measured here is a WEDGE
                    # trimmed at a nitrogen label, span 2.82 pt (AMG61). A
                    # floor of 3.5 dropped it and split that piperazine in two.
BOND_MAX  = 20.0    # pt  (nothing between 16 and 34; above that are rules and arrows)
THICK_MAX = 3.0     # pt  (plain bond ~0.43; a stereo WEDGE is a tapered quad up to ~2.5)
TRIM_SLACK_LO = 0.22  # x scale: how much SHORT of the predicted trim a bond may stop
TRIM_SLACK_HI = 0.48  # x scale: and how much past it. Relative because the
                      # Am15 series is drawn at 5 pt type, ~0.7x the rest.

MAX_BONDS = {"O": 2, "N": 3, "S": 2, "P": 5, "F": 1, "Cl": 1, "Br": 1,
             "I": 1, "Si": 4, "B": 3, "C": 4}

ELEM = {"O": "O", "N": "N", "S": "S", "P": "P", "F": "F", "Cl": "Cl",
        "Br": "Br", "I": "I", "Si": "Si", "B": "B"}


def _spine(p):
    V = []
    for it in p["items"]:
        if it[0] == "l":
            V += [(it[1].x, it[1].y), (it[2].x, it[2].y)]
        elif it[0] == "re":
            r = it[1]; V += [(r.x0, r.y0), (r.x1, r.y0), (r.x1, r.y1), (r.x0, r.y1)]
        elif it[0] == "c":
            V += [(it[1].x, it[1].y), (it[4].x, it[4].y)]
    if len(V) < 2: return None
    n = len(V)
    cx = sum(v[0] for v in V)/n; cy = sum(v[1] for v in V)/n
    sxx = sum((v[0]-cx)**2 for v in V); syy = sum((v[1]-cy)**2 for v in V)
    sxy = sum((v[0]-cx)*(v[1]-cy) for v in V)
    th = 0.5*math.atan2(2*sxy, sxx-syy)
    ux, uy = math.cos(th), math.sin(th)
    pr = [((v[0]-cx)*ux + (v[1]-cy)*uy, v) for v in V]
    lo, hi = min(pr)[0], max(pr)[0]
    span = hi-lo
    if span <= 0: return None
    loG = [v for t, v in pr if t <= lo+0.25*span]
    hiG = [v for t, v in pr if t >= hi-0.25*span]
    a = (sum(v[0] for v in loG)/len(loG), sum(v[1] for v in loG)/len(loG))
    b = (sum(v[0] for v in hiG)/len(hiG), sum(v[1] for v in hiG)/len(hiG))
    thick = 2*max(abs(-(v[0]-cx)*uy + (v[1]-cy)*ux) for v in V)
    return a, b, span, thick


def bond_segments(page, bbox=None):
    """Bond-like filled paths, as (a, b, span, thickness). Refuses nothing;
    the caller sees what was rejected."""
    keep, rej = [], []
    for p in page.get_drawings():
        if p.get("type") not in ("f", "fs"):
            rej.append(("not-filled", p["rect"])); continue
        s = _spine(p)
        if s is None:
            rej.append(("no-spine", p["rect"])); continue
        a, b, span, thick = s
        if bbox is not None:
            X0, Y0, X1, Y1 = bbox
            if not (X0 <= a[0] <= X1 and X0 <= b[0] <= X1
                    and Y0 <= a[1] <= Y1 and Y0 <= b[1] <= Y1):
                continue
        if not (BOND_MIN <= span <= BOND_MAX) or thick > THICK_MAX:
            rej.append((f"span={span:.2f} thick={thick:.2f}", p["rect"]))
            continue
        keep.append(s)
    return keep, rej


def glyph_boxes(page, bbox=None):
    """Per-character boxes of the Arial heteroatom labels, plus compound tags.

    An 'OH' label's attachment point is its O, not the span centre; splitting
    the span into per-character boxes is what makes 'OH' and 'HO' different.
    """
    atoms, tags = [], []
    for blk in page.get_text("dict")["blocks"]:
        if blk["type"] != 0: continue
        for ln in blk["lines"]:
            for sp in ln["spans"]:
                # Element labels are ArialMT; compound tags are ArialMT on
                # most pages and Helvetica on SI pages 110-116 (the Am15
                # series). Keying on "Arial" alone silently lost all 17 of
                # those compounds -- they reconstructed fine and had no name.
                if not (sp["font"].startswith("Arial")
                        or sp["font"].startswith("Helvetica")): continue
                t = sp["text"]
                if not t.strip(): continue
                x0, y0, x1, y1 = sp["bbox"]
                if bbox is not None:
                    X0, Y0, X1, Y1 = bbox
                    if not (X0 <= (x0+x1)/2 <= X1 and Y0 <= (y0+y1)/2 <= Y1):
                        continue
                st = t.strip()
                if re.fullmatch(r"(AMG|Am|G)\s?[0-9][0-9A-Za-z\-]*", st):
                    tags.append((st, sp["bbox"])); continue
                w = (x1-x0)/max(len(t), 1)
                i = 0
                while i < len(t):
                    two = t[i:i+2]
                    if two in ("Cl", "Br", "Si"):
                        el, k = two, 2
                    elif t[i] in ELEM:
                        el, k = t[i], 1
                    else:
                        i += 1; continue
                    atoms.append((el, (x0+i*w, y0, x0+(i+k)*w, y1)))
                    i += k
    return atoms, tags


def _dedupe_doubles(segs):
    n = len(segs); out = []; used = set()
    for i in range(n):
        if i in used: continue
        a0, a1, la, _ = segs[i]
        best = None
        for j in range(i+1, n):
            if j in used: continue
            b0, b1, lb, _ = segs[j]
            ux, uy = (a1[0]-a0[0])/la, (a1[1]-a0[1])/la
            bx, by = (b1[0]-b0[0])/lb, (b1[1]-b0[1])/lb
            if abs(ux*by - uy*bx) >= 0.09: continue
            dmx = (b0[0]+b1[0])/2 - (a0[0]+a1[0])/2
            dmy = (b0[1]+b1[1])/2 - (a0[1]+a1[1])/2
            along = abs(dmx*ux + dmy*uy); perp = abs(-dmx*uy + dmy*ux)
            if perp < 0.7: continue
            if perp <= DBL_PERP and along <= 0.55*max(la, lb):
                if best is None or perp < best[1]:
                    best = (j, perp)
        if best is None:
            out.append((segs[i][0], segs[i][1], 1))
        else:
            j = best[0]; used.add(j)
            k = i if segs[i][2] >= segs[j][2] else j
            out.append((segs[k][0], segs[k][1], 2))
    return out


def graph(segs, atom_labels):
    """(nodes, edges). nodes = [(element, x, y)]; edges = {(i,j): order}."""
    spans = sorted(s[2] for s in segs if 4.0 <= s[2] <= 20.0)
    scale = spans[len(spans)//2] if spans else 10.0
    tol, pad, ray_perp = TOL_F*scale, LBL_PAD_F*scale, RAY_PERP_F*scale
    nodes = []
    for el, bx in atom_labels:
        nodes.append([el, (bx[0]+bx[2])/2, (bx[1]+bx[3])/2, bx])
    nlabel = len(nodes)

    def _boxdist(p, bx):
        x0, y0, x1, y1 = bx
        dx = max(x0-p[0], 0.0, p[0]-x1)
        dy = max(y0-p[1], 0.0, p[1]-y1)
        return math.hypot(dx, dy)

    def candidates(p, q):
        """Ranked (residual, label_index) for endpoint p of the bond p<-q.

        ChemDraw TRIMS a bond at the heteroatom label's box, so p stops short
        of the atom -- but the bond's own line still POINTS AT the label
        centre, and the distance it stops short by is where that ray LEAVES
        the box. Both are predictions, so the fit has a residual and the
        candidates can be RANKED.

        Nearest-box snapping, which has neither, broke a 1,4,7-triazonane into
        an aziridine: two ring nitrogens 6 pt apart have overlapping padded
        boxes and the wrong one was 0.5 pt nearer.
        """
        out = []
        L = math.hypot(p[0]-q[0], p[1]-q[1])
        ux = uy = 0.0
        if L > 0.2:
            ux, uy = (p[0]-q[0])/L, (p[1]-q[1])/L
            for k in range(nlabel):
                cx, cy = nodes[k][1], nodes[k][2]
                t = (cx-q[0])*ux + (cy-q[1])*uy
                if t <= 0:
                    continue
                d = abs(-(cx-q[0])*uy + (cy-q[1])*ux)
                if d > ray_perp:
                    continue
                bw = (nodes[k][3][2]-nodes[k][3][0]) / 2.0
                bh = (nodes[k][3][3]-nodes[k][3][1]) / 2.0
                dux = abs(cx-q[0]) / t
                duy = abs(cy-q[1]) / t
                exit_d = min(bw/dux if dux > 1e-6 else 1e9,
                             bh/duy if duy > 1e-6 else 1e9)
                gap = t - L
                if (gap < exit_d - TRIM_SLACK_LO*scale
                        or gap > exit_d + TRIM_SLACK_HI*scale):
                    continue
                out.append((math.hypot(d, gap-exit_d), k))
        # Fallback for a bond too short for its direction to be reliable: the
        # label box, but only if the bond would come out a plausible length
        # and the label sits roughly along it. Without those two conditions
        # this fallback invented a 17.1 pt bond to a nitrogen 3 pt from a
        # wedge's narrow end, gave that N a valence of four, and split AMG61
        # in two.
        for k in range(nlabel):
            d = _boxdist(p, nodes[k][3])
            if d > pad:
                continue
            cx, cy = nodes[k][1], nodes[k][2]
            if math.hypot(cx-q[0], cy-q[1]) > 1.35*scale:
                continue
            if L > 0.2 and abs(-(cx-q[0])*uy + (cy-q[1])*ux) > 0.45*scale:
                continue
            out.append((d + 10.0, k))        # ranked below any ray fit
        out.sort()
        return out

    def plain_vertex(p):
        for k in range(nlabel, len(nodes)):
            if math.hypot(p[0]-nodes[k][1], p[1]-nodes[k][2]) <= tol:
                return k
        nodes.append(["C", p[0], p[1], None])
        return len(nodes)-1

    # --- assignment, best fit first, respecting each element's valence -----
    # A label is assigned in order of how well the bond fits it, and a label
    # already carrying its maximum number of bonds is skipped so the next-best
    # candidate gets the endpoint. Without this, six hydroxyls drawn 5 pt apart
    # on the Am15 hexacyclen head collapsed onto one oxygen with five bonds and
    # every Am15 lipid failed to parse.
    bonds = _dedupe_doubles(segs)
    ends = []
    for bi, (a, b, o) in enumerate(bonds):
        ends.append((candidates(a, b), bi, 0, o))
        ends.append((candidates(b, a), bi, 1, o))
    order = sorted(range(len(ends)),
                   key=lambda i: ends[i][0][0][0] if ends[i][0] else 1e9)
    load = [0]*nlabel
    assign = {}
    for i in order:
        cand, bi, side, order_ = ends[i]
        for res, k in cand:
            # the budget is VALENCE, not bond count: a carbonyl oxygen drawn
            # with a double bond spends two
            if load[k] + order_ <= MAX_BONDS.get(nodes[k][0], 4):
                load[k] += order_
                assign[(bi, side)] = k
                break
    edges = {}
    for bi, (a, b, o) in enumerate(bonds):
        i = assign.get((bi, 0))
        j = assign.get((bi, 1))
        if i is None: i = plain_vertex(a)
        if j is None: j = plain_vertex(b)
        if i == j: continue
        key = (min(i, j), max(i, j))
        edges[key] = max(o, edges.get(key, 0))
    return [(n[0], n[1], n[2]) for n in nodes], edges


def components(nodes, edges):
    import collections
    adj = collections.defaultdict(set)
    for (i, j) in edges: adj[i].add(j); adj[j].add(i)
    seen, out = set(), []
    for v in range(len(nodes)):
        if v in seen or v not in adj: continue
        st, c = [v], set()
        while st:
            u = st.pop()
            if u in c: continue
            c.add(u); st += list(adj[u])
        seen |= c; out.append(sorted(c))
    return out


def to_mol(nodes, edges, keep):
    from rdkit import Chem
    keep = sorted(keep)
    idx = {v: k for k, v in enumerate(keep)}
    rw = Chem.RWMol()
    for v in keep:
        rw.AddAtom(Chem.Atom(nodes[v][0]))
    for (i, j), o in edges.items():
        if i in idx and j in idx:
            rw.AddBond(idx[i], idx[j],
                       Chem.BondType.DOUBLE if o == 2 else Chem.BondType.SINGLE)
    m = rw.GetMol()
    Chem.SanitizeMol(m)

    # --- C=C geometry, from the drawing's own coordinates ----------------
    # The paper states the tails carry "Z-alkenes"; that is a claim to CHECK,
    # not to copy in. A planar conformer built from the bond coordinates lets
    # RDKit read cis/trans off the picture. Atom chirality is then cleared
    # explicitly: with every z = 0 there is no chiral volume to read, and the
    # SI says so too -- "all AMG lipids thus generated were obtained as
    # inseparable diastereomeric mixtures" (p.5) -- so an assigned centre here
    # would be an invention.
    conf = Chem.Conformer(m.GetNumAtoms())
    for k, v in enumerate(keep):
        conf.SetAtomPosition(k, (float(nodes[v][1]), -float(nodes[v][2]), 0.0))
    m.AddConformer(conf, assignId=True)
    Chem.AssignStereochemistryFrom3D(m)
    for a in m.GetAtoms():
        a.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    m.RemoveAllConformers()
    Chem.AssignStereochemistry(m, cleanIt=True, force=True)
    return m
