"""Extract molecular graphs from the VECTOR line art of main-text Fig. 1.

Fig. 1b (phenolic tails A-D) and Fig. 1c (amine cores a1-a16) are vector, not
raster, so the bonds are exact coordinates: chain length is counted from the
drawing's own geometry instead of by eye, which is this project's documented
+/-1 failure mode (ERROR_CLASSES row 19).

Stroked paths (width ~0.35) are bonds. Heteroatom LABELS are drawn as filled
glyph outlines, not text, so they come back as closed fill paths; this script
reports each glyph cluster's bbox and the vertices nearest to it, and the atom
letters themselves are read off a rendered crop by eye - the one step that
needs eyes, and the one that is not length-sensitive.
"""
import sys, json, math, collections
import pymupdf

PDF = sys.argv[1]
PAGE = int(sys.argv[2])
X0, Y0, X1, Y1 = (float(v) for v in sys.argv[3:7])
TOL = 1.2            # pt: endpoints closer than this are the same vertex
DBL = 2.2            # pt: parallel offset below this = inner line of a double bond


def collect(page):
    bonds, glyph = [], []
    for p in page.get_drawings():
        for it in p['items']:
            if it[0] == 'c' and not p.get('width'):
                # a filled CURVE is part of a letter outline (O, C, S, G...).
                # Bonds are never curved here, so these are label evidence only.
                xs = [q.x for q in it[1:]]; ys = [q.y for q in it[1:]]
                if X0 <= min(xs) and max(xs) <= X1 and Y0 <= min(ys) and max(ys) <= Y1:
                    glyph.append((min(xs), min(ys), max(xs), max(ys), 'c'))
                continue
            if it[0] != 'l':
                continue
            a, b = it[1], it[2]
            seg = (a.x, a.y, b.x, b.y)
            if not all(X0 <= x <= X1 for x in (seg[0], seg[2])):
                continue
            if not all(Y0 <= y <= Y1 for y in (seg[1], seg[3])):
                continue
            if p.get('width'):
                bonds.append(seg)
            else:
                glyph.append((min(seg[0], seg[2]), min(seg[1], seg[3]),
                              max(seg[0], seg[2]), max(seg[1], seg[3]), 'l'))
    return bonds, glyph


def length(s):
    return math.hypot(s[2] - s[0], s[3] - s[1])


def midpoint(s):
    return ((s[0] + s[2]) / 2, (s[1] + s[3]) / 2)


def parallel_pairs(segs):
    """Mark the shorter partner of each near-parallel close pair (double bond)."""
    drop = set()
    for i, a in enumerate(segs):
        if i in drop:
            continue
        va = (a[2] - a[0], a[3] - a[1])
        for j in range(i + 1, len(segs)):
            if j in drop:
                continue
            b = segs[j]
            vb = (b[2] - b[0], b[3] - b[1])
            cross = va[0] * vb[1] - va[1] * vb[0]
            if abs(cross) > 0.08 * length(a) * length(b):
                continue                      # not parallel
            ma, mb = midpoint(a), midpoint(b)
            d = math.hypot(ma[0] - mb[0], ma[1] - mb[1])
            if d < DBL:
                drop.add(j if length(b) <= length(a) else i)
                break
    return drop


def main():
    doc = pymupdf.open(PDF)
    page = doc[PAGE]
    bonds, glyph = collect(page)
    dbl = parallel_pairs(bonds)
    single = [s for i, s in enumerate(bonds) if i not in dbl]

    # cluster endpoints into vertices
    pts = []
    for s in single:
        pts += [(s[0], s[1]), (s[2], s[3])]
    verts = []
    vid = []
    for p in pts:
        for k, v in enumerate(verts):
            if math.hypot(p[0] - v[0], p[1] - v[1]) < TOL:
                vid.append(k)
                break
        else:
            verts.append(p)
            vid.append(len(verts) - 1)

    adj = collections.defaultdict(set)
    order = {}
    for i, s in enumerate(single):
        u, w = vid[2 * i], vid[2 * i + 1]
        if u == w:
            continue
        adj[u].add(w)
        adj[w].add(u)
        order[tuple(sorted((u, w)))] = 1
    # a dropped parallel partner upgrades its host bond to order 2
    for j in sorted(dbl):
        b = bonds[j]
        mb = midpoint(b)
        best, bd = None, 9e9
        for key in order:
            u, w = key
            mu = ((verts[u][0] + verts[w][0]) / 2, (verts[u][1] + verts[w][1]) / 2)
            d = math.hypot(mu[0] - mb[0], mu[1] - mb[1])
            if d < bd:
                best, bd = key, d
        if best and bd < DBL + 1.0:
            order[best] = 2

    # connected components
    seen, comps = set(), []
    for v in range(len(verts)):
        if v in seen or v not in adj:
            continue
        stack, comp = [v], []
        seen.add(v)
        while stack:
            u = stack.pop()
            comp.append(u)
            for w in adj[u]:
                if w not in seen:
                    seen.add(w)
                    stack.append(w)
        comps.append(sorted(comp))

    # glyph clusters (heteroatom labels), merged by proximity
    clusters = []
    for b in glyph:
        for c in clusters:
            if b[0] < c[2] + 1.0 and c[0] < b[2] + 1.0 and b[1] < c[3] + 1.0 and c[1] < b[3] + 1.0:
                c[0] = min(c[0], b[0]); c[1] = min(c[1], b[1])
                c[2] = max(c[2], b[2]); c[3] = max(c[3], b[3])
                c[4] += (b[4] == 'l'); c[5] += (b[4] == 'c')
                break
        else:
            clusters.append([b[0], b[1], b[2], b[3], int(b[4] == 'l'), int(b[4] == 'c')])
    # merge again (transitive)
    changed = True
    while changed:
        changed = False
        for i in range(len(clusters)):
            for j in range(i + 1, len(clusters)):
                a, c = clusters[i], clusters[j]
                if a[0] < c[2] + 1.0 and c[0] < a[2] + 1.0 and a[1] < c[3] + 1.0 and c[1] < a[3] + 1.0:
                    a[0] = min(a[0], c[0]); a[1] = min(a[1], c[1])
                    a[2] = max(a[2], c[2]); a[3] = max(a[3], c[3])
                    a[4] += c[4]; a[5] += c[5]
                    clusters.pop(j); changed = True
                    break
            if changed:
                break

    # Heteroatom labels break the drawn skeleton: a bond toward a labelled atom
    # is drawn SHORT, stopping at the label box, so its endpoint is not a carbon -
    # it IS that atom. Endpoints inside GAP of a glyph box are therefore MERGED
    # into the glyph node, not linked to it. Linking instead would add one phantom
    # carbon per heteroatom bond, which is the same +/-1 chain error this whole
    # geometric route exists to remove.
    GAP = float(sys.argv[7]) if len(sys.argv) > 7 else 3.2
    gnode = {gi: len(verts) + gi for gi in range(len(clusters))}
    merge = {}
    for vi, v in enumerate(verts):
        best, bd = None, 9e9
        for gi, c in enumerate(clusters):
            dx = max(c[0] - v[0], 0, v[0] - c[2])
            dy = max(c[1] - v[1], 0, v[1] - c[3])
            d = math.hypot(dx, dy)
            if d < bd:
                best, bd = gi, d
        if best is not None and bd < GAP:
            merge[vi] = gnode[best]

    def rep(v):
        return merge.get(v, v)

    adj2 = collections.defaultdict(set)
    order2 = {}
    for (u, w), o in order.items():
        a, b = rep(u), rep(w)
        if a == b:
            continue
        adj2[a].add(b); adj2[b].add(a)
        key = tuple(sorted((a, b)))
        order2[key] = max(order2.get(key, 1), o)
    adj, order = adj2, order2

    seen, comps = set(), []
    for v in list(adj):
        if v in seen:
            continue
        stack, comp = [v], []
        seen.add(v)
        while stack:
            u = stack.pop()
            comp.append(u)
            for w in adj[u]:
                if w not in seen:
                    seen.add(w); stack.append(w)
        comps.append(sorted(comp))

    out = {
        "pdf": PDF, "page": PAGE, "clip": [X0, Y0, X1, Y1],
        "n_bond_segments": len(bonds), "n_double_partners": len(dbl),
        "vertices": [[round(v[0], 2), round(v[1], 2)] for v in verts]
                     + [[round((c[0] + c[2]) / 2, 2), round((c[1] + c[3]) / 2, 2)] for c in clusters],
        "n_carbon_vertices": len(verts),
        "bonds": [[u, w, o] for (u, w), o in sorted(order.items())],
        "merged_endpoints": len(merge),
        "components": comps,
        "glyphs": [[round(x, 2) for x in c] for c in clusters],
    }
    json.dump(out, sys.stdout, indent=1)


if __name__ == "__main__":
    main()
