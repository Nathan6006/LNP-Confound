#!/usr/bin/env python3
"""Count alkyl-chain carbons from the PDF's VECTOR bond objects, exactly.

The figure is vector art, so every bond is a stroked path with real coordinates. Instead
of rasterising and estimating apexes (which proved unreliable - it disagreed with both
HRMS anchors), this reads the bond segments directly, snaps their endpoints into shared
vertices, and counts the atoms of each connected component. That is an exact graph, so
chain length is counted rather than inferred.

Validated against 9A and 10A, whose identities are fixed independently by the HRMS in
Supplementary Figs 5-8 (di-n-octylamine, di-n-decylamine).
"""
import json
import re
import sys
from collections import defaultdict

import pdfplumber

PDF = "s41563-020-00886-0.pdf"
PAGE = 1
TOL = 1.2          # pt: endpoints closer than this are the same vertex
MIN_BOND = 2.0     # pt: ignore hairline artefacts


def segments(page, bbox):
    """Stroked segments inside bbox, as ((x0,y0),(x1,y1))."""
    x0, y0, x1, y1 = bbox
    segs = []
    for c in list(page.curves or []) + list(page.lines or []):
        pts = c.get("pts") or []
        if len(pts) < 2:
            continue
        for a, b in zip(pts, pts[1:]):
            if not all(x0 <= p[0] <= x1 and y0 <= p[1] <= y1 for p in (a, b)):
                continue
            d = ((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2) ** 0.5
            if d >= MIN_BOND:
                segs.append((tuple(a), tuple(b)))
    return segs


def graph(segs):
    """Snap endpoints to shared vertices; return adjacency."""
    verts = []

    def vid(p):
        for i, q in enumerate(verts):
            if abs(p[0] - q[0]) <= TOL and abs(p[1] - q[1]) <= TOL:
                return i
        verts.append(p)
        return len(verts) - 1

    adj = defaultdict(set)
    for a, b in segs:
        ia, ib = vid(a), vid(b)
        if ia != ib:
            adj[ia].add(ib)
            adj[ib].add(ia)
    return verts, adj


def components(adj, n):
    seen, comps = set(), []
    for s in range(n):
        if s in seen:
            continue
        stack, comp = [s], []
        seen.add(s)
        while stack:
            u = stack.pop()
            comp.append(u)
            for v in adj.get(u, ()):
                if v not in seen:
                    seen.add(v)
                    stack.append(v)
        comps.append(comp)
    return comps


def main():
    with pdfplumber.open(PDF) as pdf:
        page = pdf.pages[PAGE]
        labels = {}
        for w in page.extract_words():
            if re.fullmatch(r"\d{1,2}A", w["text"]):
                labels.setdefault(w["text"], w)
        # y in pdfplumber 'pts' space is bottom-up; convert label tops
        H = page.height
        order = [f"{i}A" for i in range(1, 14)]
        rows = sorted(((labels[k]["top"], k) for k in order))
        res = {}
        for i, (top, k) in enumerate(rows):
            nxt = rows[i + 1][0] if i + 1 < len(rows) else top + 22
            # band in bottom-up coords
            yb, yt = H - (nxt - 2.0), H - (top - 5.0)
            bbox = (labels[k]["x0"] + 10, yb, 342.0, yt)
            segs = segments(page, bbox)
            verts, adj = graph(segs)
            comps = components(adj, len(verts))
            comps = [c for c in comps if len(c) >= 3]
            comps.sort(key=len, reverse=True)
            res[k] = {"n_components": len(comps),
                      "sizes": [len(c) for c in comps],
                      "n_segments": len(segs)}
        anchors = {"9A": 8, "10A": 10}
        print(f"{'amine':>6} {'segs':>5} {'components (vertex counts)':<30} anchor")
        for k in order:
            d = res[k]
            note = ""
            if k in anchors:
                note = f"expect two chains of {anchors[k]} C"
            print(f"{k:>6} {d['n_segments']:>5} {str(d['sizes']):<30} {note}")
        json.dump(res, open("02_vector_chains.json", "w"), indent=1)


if __name__ == "__main__":
    main()
