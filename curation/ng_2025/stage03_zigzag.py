"""Count the vertices of an alkyl zig-zag by its own geometry.

Skeletonising these Supplementary rasters fails: the strokes are 2-3 px wide
and anti-aliased at ~106 dpi, so the skeleton comes back branched everywhere.
A zig-zag does not need a skeleton. Take the ink centroid of each pixel column,
smooth it, and count the turning points: for a chain drawn left-to-right every
carbon is one turn, and the two ends are the remaining vertices. The count is
therefore turns + 2, derived from the drawing rather than from counting by eye.

Usage: stage03_zigzag.py IMG x0,y0,x1,y1 [min_area]
prints one line per connected component in the crop.
"""
import sys, collections
import numpy as np
from PIL import Image
from skimage.measure import label, regionprops

img, box = sys.argv[1], [int(v) for v in sys.argv[2].split(',')]
min_area = int(sys.argv[3]) if len(sys.argv) > 3 else 60

im = Image.open(img).convert('RGB').split()[0].crop(box)
a = np.array(im)
bw = a > int(a.min()) + 6
lab = label(bw)


def turns(mask, axis_horizontal=True):
    m = mask if axis_horizontal else mask.T
    cols = []
    for x in range(m.shape[1]):
        ys = np.nonzero(m[:, x])[0]
        cols.append(ys.mean() if len(ys) else np.nan)
    y = np.array(cols, float)
    good = ~np.isnan(y)
    if good.sum() < 5:
        return 0, []
    xs = np.nonzero(good)[0]
    y = np.interp(np.arange(len(y)), xs, y[good])
    k = max(3, len(y) // 40) | 1
    ker = np.ones(k) / k
    ys = np.convolve(y, ker, mode='same')
    ys[:k] = y[:k]; ys[-k:] = y[-k:]
    d = np.sign(np.diff(ys))
    d = d[d != 0]
    sw = [i for i in range(1, len(d)) if d[i] != d[i - 1]]
    return len(sw), sw


rows = []
for r in regionprops(lab):
    if r.area < min_area:
        continue
    sub = (lab[r.bbox[0]:r.bbox[2], r.bbox[1]:r.bbox[3]] == r.label)
    h, w = sub.shape
    horiz = w >= h
    n, _ = turns(sub, horiz)
    rows.append((r.bbox, r.area, w, h, 'H' if horiz else 'V', n, n + 2))
rows.sort(key=lambda t: -t[1])
for bbox, area, w, h, o, n, v in rows:
    print(f'bbox(y0,x0,y1,x1)={tuple(int(b) for b in bbox)} area={area:6.0f} {w:4d}x{h:<4d} {o} turns={n:3d} vertices={v:3d}')
