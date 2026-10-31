#!/usr/bin/env python3
"""Read the 51 printed values out of Figure 1d (Tian et al., Adv Mater 2026).

Figure 1d is a 3 x 17 grid = 51 cells. Rows are the cholesterol tails B1/B2/B3;
columns are the 17 lipid cores (5 commercial: MC3-C, SM-102-C, Lipid 14-C,
OF-2-C, TT3-C; plus 12 synthesized: C1-C6 x O1/O2). Every cell PRINTS its value
as "D.D x 10^E", so this is symbolic reading, not pixel digitization of a colour.

WHY NOT OCR. Tesseract misreads the raised exponent (1.2x10^5 came back as
"109", and a decimal point was dropped). Instead each cell is segmented into
glyph connected-components, which the rigid layout makes unambiguous:

      d1  .  d2  x  1  0  ^E        <- always exactly 7 glyphs

All 51 cells segment to exactly 7, which is itself a check on the geometry.
Only d1, d2 and E carry information; '.', 'x', '1' and '0' are fixed furniture
whose shape is asserted. Glyph bitmaps are binarized (killing the antialiasing
that varies with cell fill), normalized with ASPECT RATIO PRESERVED, and grouped
by complete-linkage hierarchical clustering on Jaccard distance. A human labels
10 mantissa + 4 exponent cluster representatives instead of 153 glyphs; the
montages written to review/ are that labelling surface.

Aspect preservation and binarization both matter: stretching glyphs to a square
collapsed 0/9, 6/8 and 5/6 into shared clusters, and 5/6 in an exponent is a
ten-fold error. Complete linkage matters too - greedy centroid growth chained
6 into 8.

INDEPENDENT CHECK. Fig 1d also carries a CONTINUOUS white-to-blue colour bar
(0 to 3e6). Cell fill colour is sampled and mapped through that bar, giving a
value estimate that never touches the glyphs. The bar saturates below ~1e5, so
it checks the large values and the rank order, not the small ones.

Outputs 04c_fig1d_51.csv, 07c_fig1d_qc.json, review/fig1d_overlay.png.
"""
import collections
import csv
import json
import os
import sys

import numpy as np
import PIL.Image as Image
import PIL.ImageDraw as ImageDraw
from scipy import ndimage
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import pdist, squareform
from scipy.stats import spearmanr

PAGE = sys.argv[1] if len(sys.argv) > 1 else "00_pages/p2_600.png"
OUT_CSV = "04c_fig1d_51.csv"
OUT_QC = "07c_fig1d_qc.json"
REVIEW = "review"

CELL_X0, CELL_W, NCOL = 860, 160.0, 17
ROW_BANDS = [(4712, 4800), (4820, 4908), (4932, 5020)]
TAILS = ["B1", "B2", "B3"]
CORES = ["MC3-C", "SM-102-C", "Lipid 14-C", "OF-2-C", "TT3-C",
         "C1O1", "C1O2", "C2O1", "C2O2", "C3O1", "C3O2",
         "C4O1", "C4O2", "C5O1", "C5O2", "C6O1", "C6O2"]
COMMERCIAL = set(CORES[:5])

CBAR_Y, CBAR_X, CBAR_MAX = (5112, 5140), (868, 3580), 3.0e6
CROSS = np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]])

N_MANTISSA_CLUSTERS = 14
N_EXPONENT_CLUSTERS = 4

# Assigned once by eye from review/fig1d_glyphs_*.png, then fixed.
# Empty = first pass: write the montages and stop.
MANTISSA_LABELS = {0: 1, 1: 6, 2: 2, 3: 1, 4: 7, 5: 3, 6: 3, 7: 4,
                   8: 0, 9: 2, 10: 4, 11: 9, 12: 5, 13: 8}
EXPONENT_LABELS = {0: 4, 1: 5, 2: 6, 3: 3}


def segment_cell(gray, y0, y1, ci):
    """Glyph boxes for one cell, left to right, with touching pairs split."""
    x0 = int(CELL_X0 + CELL_W * ci) + 6
    x1 = int(CELL_X0 + CELL_W * (ci + 1)) - 6
    sub = gray[y0:y1, x0:x1]
    # Cell fills run from near-white to deep blue, so a fixed offset from the
    # background either misses ink in dark cells or floods pale ones. Threshold
    # at a fixed fraction of each cell's own ink-to-fill contrast, with the ink
    # level taken as a low percentile rather than the minimum, which is noise.
    bg, ink = np.median(sub), np.percentile(sub, 2)
    mask = sub < bg - 0.45 * (bg - ink)
    lab, _ = ndimage.label(mask, structure=CROSS)
    boxes = []
    for i, sl in enumerate(ndimage.find_objects(lab)):
        if int((lab[sl] == i + 1).sum()) < 12:
            continue
        boxes.append([sl[1].start, sl[0].start, sl[1].stop, sl[0].stop])
    boxes.sort()

    # '0' and the raised exponent touch diagonally in most cells. A merged
    # component is far wider than any single glyph; cut it at the column with
    # the fewest ink pixels, searched only in the interior so the cut cannot
    # shave an edge off either glyph.
    out = []
    for bx, by, ex, ey in boxes:
        if ex - bx > 26:
            colsum = mask[by:ey, bx:ex].sum(axis=0)
            lo, hi = 8, (ex - bx) - 8
            cut = bx + lo + int(np.argmin(colsum[lo:hi]))
            for a, b in ((bx, cut), (cut, ex)):
                ys, xs = np.nonzero(mask[by:ey, a:b])
                if len(ys):
                    out.append([a + xs.min(), by + ys.min(),
                                a + xs.max() + 1, by + ys.max() + 1])
        else:
            out.append([bx, by, ex, ey])
    out.sort()
    return out, mask, x0, y0


def norm_bitmap(mask, box, size=40):
    """Binarized, aspect-preserved, centred glyph bitmap."""
    bx, by, ex, ey = box
    patch = mask[by:ey, bx:ex].astype("uint8") * 255
    h, w = patch.shape
    s = (size - 8) / max(h, w)
    nh, nw = max(1, int(round(h * s))), max(1, int(round(w * s)))
    im = Image.fromarray(patch).resize((nw, nh), Image.LANCZOS)
    canvas = Image.new("L", (size, size), 0)
    canvas.paste(im, ((size - nw) // 2, (size - nh) // 2))
    return (np.array(canvas, dtype=float) > 110).astype(float).ravel()


def group(vecs, k):
    """Complete-linkage clustering into exactly k groups, with a tightness
    score per group (max within-group Jaccard distance)."""
    d = pdist(vecs, "jaccard")
    labels = fcluster(linkage(d, "complete"), k, "maxclust")
    sq = squareform(d)
    groups, tight = [], []
    for g in range(1, k + 1):
        idx = [i for i, l in enumerate(labels) if l == g]
        groups.append(idx)
        tight.append(float(sq[np.ix_(idx, idx)].max()) if len(idx) > 1 else 0.0)
    order = np.argsort([-len(g) for g in groups])
    return [groups[i] for i in order], [tight[i] for i in order]


def montage(path, masks, boxes, groups, tight, cols=14):
    th = 34
    im = Image.new("RGB", (th * cols + 150, th * len(groups) + 10), "white")
    d = ImageDraw.Draw(im)
    for gi, mem in enumerate(groups):
        d.text((4, th * gi + 6), f"c{gi}  n={len(mem)}", fill="black")
        d.text((4, th * gi + 18), f"tight={tight[gi]:.2f}", fill="gray")
        for j, k in enumerate(mem[:cols]):
            bx = boxes[k]
            patch = masks[k][bx[1]:bx[3], bx[0]:bx[2]]
            g = Image.fromarray(((1 - patch) * 255).astype("uint8")).resize((28, 28))
            im.paste(g, (144 + th * j, th * gi + 4))
    im.save(path)


def main():
    os.makedirs(REVIEW, exist_ok=True)
    img = np.array(Image.open(PAGE).convert("RGB")).astype(int)
    gray = img.mean(axis=2)

    cells = []
    mv, mb, mref, ev, eb, eref = [], [], [], [], [], []
    shape_issues, furniture_fail = [], []

    for ri, (y0, y1) in enumerate(ROW_BANDS):
        for ci in range(NCOL):
            boxes, mask, ox, oy = segment_cell(gray, y0, y1, ci)
            rec = {"row": ri, "col": ci, "n_glyphs": len(boxes),
                   "origin": (ox, oy), "boxes": boxes, "mask": mask}
            cells.append(rec)
            if len(boxes) != 7:
                shape_issues.append([TAILS[ri], CORES[ci], len(boxes)])
                continue
            hs = [b[3] - b[1] for b in boxes]
            ws = [b[2] - b[0] for b in boxes]
            ok = (hs[1] <= 8 and ws[1] <= 10             # the decimal point
                  and 14 <= ws[3] <= 22 and hs[3] <= 22  # the multiplication x
                  and ws[4] <= 16                        # the '1' of 10
                  and boxes[6][1] < boxes[5][1])         # exponent sits raised
            rec["furniture_ok"] = bool(ok)
            if not ok:
                furniture_fail.append([TAILS[ri], CORES[ci]])
            for slot in (0, 2):
                mv.append(norm_bitmap(mask, boxes[slot]))
                mb.append(boxes[slot])
                mref.append((mask, ri, ci, slot))
            ev.append(norm_bitmap(mask, boxes[6]))
            eb.append(boxes[6])
            eref.append((mask, ri, ci, 6))

    mg, mt = group(np.array(mv), N_MANTISSA_CLUSTERS)
    eg, et = group(np.array(ev), N_EXPONENT_CLUSTERS)
    montage(f"{REVIEW}/fig1d_glyphs_mantissa.png", [r[0] for r in mref], mb, mg, mt)
    montage(f"{REVIEW}/fig1d_glyphs_exponent.png", [r[0] for r in eref], eb, eg, et)

    print(f"cells {len(cells)}  7-glyph {sum(1 for c in cells if c['n_glyphs'] == 7)}/51"
          f"  shape issues {shape_issues}  furniture fails {furniture_fail}")
    print(f"mantissa {len(mv)} glyphs -> sizes {[len(g) for g in mg]}  "
          f"tightness {[round(t, 2) for t in mt]}")
    print(f"exponent {len(ev)} glyphs -> sizes {[len(g) for g in eg]}  "
          f"tightness {[round(t, 2) for t in et]}")

    if not MANTISSA_LABELS or not EXPONENT_LABELS:
        print(f"\nfirst pass. Label clusters in {REVIEW}/fig1d_glyphs_*.png, "
              "fill MANTISSA_LABELS / EXPONENT_LABELS, re-run.")
        return

    dig = {}
    for gi, mem in enumerate(mg):
        for k in mem:
            dig[("m", k)] = MANTISSA_LABELS[gi]
    for gi, mem in enumerate(eg):
        for k in mem:
            dig[("e", k)] = EXPONENT_LABELS[gi]

    parsed = {}
    for k, (_, ri, ci, slot) in enumerate(mref):
        parsed.setdefault((ri, ci), {})[slot] = dig[("m", k)]
    for k, (_, ri, ci, _slot) in enumerate(eref):
        parsed.setdefault((ri, ci), {})[6] = dig[("e", k)]

    strip = np.median(img[CBAR_Y[0]:CBAR_Y[1], CBAR_X[0]:CBAR_X[1]].astype(float), axis=0)
    rows = []
    for rec in cells:
        p = parsed.get((rec["row"], rec["col"]))
        if not p:
            continue
        mant = float(f"{p[0]}.{p[2]}")
        expo = int(p[6])
        x0 = int(CELL_X0 + CELL_W * rec["col"]) + 8
        x1 = int(CELL_X0 + CELL_W * (rec["col"] + 1)) - 8
        y0, y1 = ROW_BANDS[rec["row"]]
        patch = img[y0:y1, x0:x1].reshape(-1, 3)
        light = patch[patch.mean(axis=1) > 150]
        fill = np.median(light, axis=0) if len(light) else np.array([255., 255., 255.])
        dd = np.linalg.norm(strip - fill, axis=1)
        i = int(np.argmin(dd))
        core = CORES[rec["col"]]
        rows.append({
            "lipid_name": f"{core}{TAILS[rec['row']]}",
            "core": core, "chol_tail": TAILS[rec["row"]],
            "core_origin": "commercial" if core in COMMERCIAL else "synthesized",
            "value_printed": mant * (10 ** expo),
            "mantissa": mant, "exponent": expo,
            "value_from_colourbar": round(CBAR_MAX * i / (len(strip) - 1), 1),
            "colourbar_match_distance": round(float(dd[i]), 1),
            "furniture_ok": rec.get("furniture_ok"),
            "cell_fill_rgb": f"{int(fill[0])},{int(fill[1])},{int(fill[2])}",
            "assay": "whole-brain luciferase, C57BL/6, I.V. 0.1 mg/kg mLuc, 6 h",
            "readout_metric": "brain total flux (photons/s)",
            "value_scale": "raw", "n_replicates": 3,
            "source_figure": "Fig 1d",
            "extraction_method": "printed_numeral_glyph_cluster",
            "source_pdf": "Adv Mater 2026 Tian (main)",
        })

    with open(OUT_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    overlay(img, cells, rows)
    report(rows, shape_issues, furniture_fail, mg, mt, eg, et)


def overlay(img, cells, rows):
    """Draw every sampled glyph box, the cell fill patch, and the parsed value."""
    y0 = ROW_BANDS[0][0] - 140
    y1 = CBAR_Y[1] + 50
    xoff = CELL_X0 - 110
    S = 1
    crop = Image.fromarray(img.astype("uint8")[y0:y1, xoff:3660])
    crop = crop.resize((crop.width * S, crop.height * S), Image.LANCZOS)
    d = ImageDraw.Draw(crop)
    by = {(r["chol_tail"], r["core"]): r for r in rows}
    for rec in cells:
        ox, oy = rec["origin"]
        for j, b in enumerate(rec["boxes"]):
            col = (200, 0, 200) if j in (0, 2, 6) else (255, 150, 0)
            d.rectangle([(ox + b[0] - xoff) * S, (oy + b[1] - y0) * S,
                         (ox + b[2] - xoff) * S, (oy + b[3] - y0) * S],
                        outline=col, width=1)
        r = by.get((TAILS[rec["row"]], CORES[rec["col"]]))
        if r:
            d.text(((CELL_X0 + CELL_W * rec["col"] + 8 - xoff) * S,
                    (ROW_BANDS[rec["row"]][1] - y0) * S - 8),
                   f"{r['value_printed']:.3g}", fill=(200, 0, 0))
    d.rectangle([(CBAR_X[0] - xoff) * S, (CBAR_Y[0] - y0) * S,
                 (CBAR_X[1] - xoff) * S, (CBAR_Y[1] - y0) * S],
                outline=(0, 0, 255), width=2)
    d.text((6, 6), "magenta = information glyphs (d1, d2, exponent)   "
                   "orange = asserted furniture (. x 1 0)   "
                   "green = parsed value   blue = colour-bar sample", fill=(0, 0, 0))
    crop.save(f"{REVIEW}/fig1d_overlay.png")


def report(rows, shape_issues, furniture_fail, mg, mt, eg, et):
    v = np.array([r["value_printed"] for r in rows])
    c = np.array([r["value_from_colourbar"] for r in rows])
    big = v > 2e5
    out = {
        "paper": "Tian et al., Adv Mater 2026, Fig 1d",
        "n_cells": len(rows), "expected": 51,
        "grid": "3 cholesterol tails (B1,B2,B3) x 17 cores "
                "(5 commercial + 12 synthesized C1-C6 x O1/O2)",
        "method": "printed numerals; glyph connected components, binarized and "
                  "aspect-preserved, complete-linkage clustered on Jaccard distance",
        "cells_segmenting_to_7_glyphs": f"{51 - len(shape_issues)}/51",
        "shape_issues": shape_issues,
        "furniture_assertions_failed": furniture_fail,
        "mantissa_cluster_sizes": [len(g) for g in mg],
        "mantissa_cluster_tightness": [round(t, 3) for t in mt],
        "exponent_cluster_sizes": [len(g) for g in eg],
        "exponent_cluster_tightness": [round(t, 3) for t in et],
        "exponent_distribution": dict(collections.Counter(
            int(r["exponent"]) for r in rows)),
        "independent_colourbar_check": {
            "what": "cell fill mapped through the continuous white-to-blue bar "
                    "(0-3e6), independent of the glyphs. Saturates below ~1e5, so "
                    "it checks the large values and rank order, not the small ones.",
            "spearman_all": round(float(spearmanr(v, c).statistic), 4),
            "spearman_above_2e5": round(float(spearmanr(v[big], c[big]).statistic), 4)
            if big.sum() > 2 else None,
            "n_above_2e5": int(big.sum()),
        },
        "max_value": {"lipid": rows[int(np.argmax(v))]["lipid_name"],
                      "value": float(v.max())},
        "review_image": f"{REVIEW}/fig1d_overlay.png",
    }
    json.dump(out, open(OUT_QC, "w"), indent=1)
    print(f"\n51 values read. max = {out['max_value']['lipid']} at {v.max():.3g}")
    print(f"exponent distribution: {out['exponent_distribution']}")
    print(f"colour-bar cross-check: Spearman all "
          f"{out['independent_colourbar_check']['spearman_all']}, "
          f">2e5 (n={int(big.sum())}) "
          f"{out['independent_colourbar_check']['spearman_above_2e5']}")
    print(f"overlay -> {REVIEW}/fig1d_overlay.png")


if __name__ == "__main__":
    main()
