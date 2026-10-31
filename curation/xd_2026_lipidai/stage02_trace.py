"""Trace zig-zag skeletal structures in Fig 5A and report polyline vertices.

Rationale (CLAUDE.md error ledger): manual carbon counting is the documented
+-1 failure mode. This derives vertices geometrically instead, and writes an
annotated overlay so the result can be checked by eye against the figure.
"""
import sys, json
import numpy as np
from PIL import Image, ImageDraw
from skimage.morphology import skeletonize
from skimage.measure import label, regionprops

D = '/Users/nathanliu/LNP/confound/curation/xd_2026_lipidai/00_pages/png'

def rdp(pts, eps):
    pts = np.asarray(pts, float)
    if len(pts) < 3:
        return list(range(len(pts)))
    a, b = pts[0], pts[-1]
    ab = b - a
    n = np.hypot(*ab)
    if n == 0:
        d = np.hypot(*(pts - a).T)
    else:
        d = np.abs(np.cross(ab, pts - a)) / n
    i = int(np.argmax(d))
    if d[i] > eps:
        l = rdp(pts[:i+1], eps); r = rdp(pts[i:], eps)
        return l[:-1] + [i + j for j in r]
    return [0, len(pts) - 1]

NB = [(-1,-1),(-1,0),(-1,1),(0,-1),(0,1),(1,-1),(1,0),(1,1)]

def trace_component(sk):
    """Return list of pixel paths between endpoints/junctions of one skeleton."""
    ys, xs = np.nonzero(sk)
    px = set(zip(ys.tolist(), xs.tolist()))
    deg = {p: sum(((p[0]+dy, p[1]+dx) in px) for dy, dx in NB) for p in px}
    nodes = {p for p, d in deg.items() if d != 2}
    paths = []
    if not nodes:                      # pure loop
        return [sorted(px)]
    seen = set()
    for s in nodes:
        for dy, dx in NB:
            nxt = (s[0]+dy, s[1]+dx)
            if nxt not in px or (s, nxt) in seen:
                continue
            path = [s, nxt]; seen.add((s, nxt)); seen.add((nxt, s))
            cur, prev = nxt, s
            while cur not in nodes:
                adv = None
                for dy2, dx2 in NB:
                    c = (cur[0]+dy2, cur[1]+dx2)
                    if c in px and c != prev and (cur, c) not in seen:
                        adv = c; break
                if adv is None:
                    break
                seen.add((cur, adv)); seen.add((adv, cur))
                path.append(adv); prev, cur = cur, adv
            paths.append(path)
    return paths

def main(img_path, out_json, out_png, min_area=60, eps=3.0, min_path=18):
    im = Image.open(img_path).convert('L')
    a = np.array(im)
    bw = a < 128
    lab = label(bw)
    comps = []
    for r in regionprops(lab):
        if r.area < min_area:
            continue
        sub = (lab[r.bbox[0]:r.bbox[2], r.bbox[1]:r.bbox[3]] == r.label)
        sk = skeletonize(sub)
        paths = [p for p in trace_component(sk) if len(p) >= min_path]
        verts = []
        for p in paths:
            pts = [(x + r.bbox[1], y + r.bbox[0]) for y, x in p]
            idx = rdp(pts, eps)
            verts.append([pts[i] for i in idx])
        comps.append(dict(bbox=[int(v) for v in r.bbox], area=int(r.area),
                          n_paths=len(paths),
                          polylines=[[[int(x), int(y)] for x, y in v] for v in verts]))
    comps.sort(key=lambda c: -c['area'])
    json.dump(comps, open(out_json, 'w'), indent=1)
    ov = Image.open(img_path).convert('RGB')
    d = ImageDraw.Draw(ov)
    for c in comps:
        for pl in c['polylines']:
            for (x, y) in pl:
                d.ellipse([x-5, y-5, x+5, y+5], outline=(255, 0, 0), width=3)
    ov.save(out_png)
    print(f'{len(comps)} components >= {min_area}px')
    for c in comps[:12]:
        print(' area', c['area'], 'bbox', c['bbox'], 'paths', c['n_paths'],
              'verts', [len(p) for p in c['polylines']])

if __name__ == '__main__':
    main(*sys.argv[1:4])
