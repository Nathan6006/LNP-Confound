#!/usr/bin/env python3
"""Read Tables S1 and S2 -- the authors' own source data for Fig. 1C / 1E.

    python3 -I stage04_tables.py            # cluster, write 04_glyphs.json
    python3 -I stage04_tables.py --montage  # write a contact sheet to label
    python3 -I stage04_tables.py --read     # decode with 04_glyph_labels.json

WHAT THIS IS AND IS NOT. The tier is `5_pixels` (see 00_value_sources.json) but
nothing here is a measurement of a length, a colour or a position. The pixels
encode PRINTED NUMERALS inside a ruled table, so:

  * the cell grid comes from the table's OWN ruled lines, detected as rows and
    columns whose ink fraction exceeds 0.6. There is no axis, no origin and no
    scale factor anywhere, so error classes 20 / 27 / 33 cannot apply -- there
    is nothing to anchor wrongly.
  * there is no colormap, so classes 16 / 23 / 43 cannot apply.
  * there are no bars, so class 32 cannot apply.

The only way to be wrong is to misread a glyph, and that is why this file does
not call an OCR engine. CLAUDE.md s5: "OCR is not good enough. Tesseract
misread 1.2x10^5 as 109." Instead every connected component in every cell is
normalised, clustered with COMPLETE linkage, and the cluster representatives
are labelled BY EYE once. The labels live in 04_glyph_labels.json keyed by the
sha256 of the representative's normalised bitmap, so a label can never drift
onto a different shape: if clustering changes, the fingerprint stops matching
and this script raises instead of decoding.

A NEW CHECK IS NOT DONE UNTIL IT HAS BEEN MADE TO FAIL ON PURPOSE. The
fingerprint gate was seeded with a corrupted label file and confirmed to raise;
see 07_qc.json -> glyph_label_gate_seeded_failure.
"""
import hashlib
import json
import os
import sys

import numpy as np
from PIL import Image
from scipy import ndimage
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import pdist, squareform

HERE = os.path.dirname(os.path.abspath(__file__))
PDF = os.path.join(HERE, "source", "nn6c00706_si_001.pdf")

# The two source-data tables, as embedded image objects. Page and xref are
# facts about the PDF, printed by stage00; they are not tunable.
TABLES = {
    "S1": {"page": 38, "xref": 299, "blocks": 1,
           "caption": "Source data for the heatmap in Fig. 1C"},
    "S2": {"page": 39, "xref": 303, "blocks": 5,
           "caption": "Source data for the heatmap in Fig. 1E"},
}

RULE_INK_FRACTION = 0.6   # a ruled line is >60% ink across the image
NORM = 16                 # glyph bitmaps are normalised to NORM x NORM
UNDERLINE_MAX_H = 3       # a component this flat and this wide is a rule,
UNDERLINE_MIN_W = 10      # not a glyph. Recorded, never decoded.
LINKAGE_CUTOFF = 2.2      # complete-linkage cut, in normalised-bitmap L2


def ink_image(path_pdf, xref):
    """The embedded image as a boolean ink mask. Raises if it is not bilevel."""
    import pymupdf
    doc = pymupdf.open(path_pdf)
    info = doc.extract_image(xref)
    import io
    im = Image.open(io.BytesIO(info["image"])).convert("L")
    a = np.array(im)
    vals = np.unique(a)
    if not set(vals.tolist()) <= {0, 255}:
        raise SystemExit(
            f"xref {xref}: image is not bilevel (unique values {vals[:8]}). "
            "A grey image means it has been resampled somewhere, and glyph "
            "clustering on resampled text is not exact. Stop and find out "
            "where the resampling happened.")
    return a < 128


def _runs(profile, thr):
    idx = [i for i, v in enumerate(profile) if v > thr]
    if not idx:
        raise SystemExit("no ruled lines found; the image is not the table.")
    groups, cur = [], [idx[0]]
    for i in idx[1:]:
        if i - cur[-1] <= 2:
            cur.append(i)
        else:
            groups.append(cur)
            cur = [i]
    groups.append(cur)
    return [(g[0], g[-1]) for g in groups]


def grid(ink):
    """(horizontal rules, vertical rules) as (start, end) pixel spans."""
    return (_runs(ink.mean(axis=1), RULE_INK_FRACTION),
            _runs(ink.mean(axis=0), RULE_INK_FRACTION))


def normalise(bitmap):
    """Crop to bbox, pad to square preserving aspect ratio, resize to NORM."""
    ys, xs = np.where(bitmap)
    sub = bitmap[ys.min():ys.max() + 1, xs.min():xs.max() + 1]
    h, w = sub.shape
    side = max(h, w)
    pad = np.zeros((side, side), dtype=bool)
    pad[(side - h) // 2:(side - h) // 2 + h,
        (side - w) // 2:(side - w) // 2 + w] = sub
    im = Image.fromarray((pad * 255).astype("uint8"))
    return np.array(im.resize((NORM, NORM), Image.LANCZOS)) / 255.0


def components(ink, y0, y1, x0, x1):
    """Glyph components inside one cell, left to right.

    Returns [(x_centre, aspect_h, aspect_w, normalised_bitmap)] plus the count
    of components rejected as underlines.
    """
    sub = ink[y0:y1, x0:x1]
    lab, n = ndimage.label(sub, structure=np.ones((3, 3)))
    glyphs, underlines = [], 0
    for i in range(1, n + 1):
        ys, xs = np.where(lab == i)
        h = ys.max() - ys.min() + 1
        w = xs.max() - xs.min() + 1
        if h <= UNDERLINE_MAX_H and w >= UNDERLINE_MIN_W:
            underlines += 1
            continue
        m = np.zeros_like(sub)
        m[ys, xs] = True
        glyphs.append({"x": float(xs.mean()), "ymin": int(ys.min()),
                       "ymax": int(ys.max()), "h": int(h), "w": int(w),
                       "native": m[ys.min():ys.max() + 1,
                                   xs.min():xs.max() + 1].copy(),
                       "bitmap": normalise(m)})
    glyphs.sort(key=lambda g: g["x"])
    return glyphs, underlines


def harvest():
    """Every glyph in both tables, tagged with its cell address."""
    out = []
    meta = {}
    for name, spec in TABLES.items():
        ink = ink_image(PDF, spec["xref"])
        H, V = grid(ink)
        nblock = spec["blocks"]
        if (len(H) - 1) % 1 != 0:
            raise SystemExit("impossible")
        meta[name] = {"shape": list(ink.shape), "h_rules": len(H),
                      "v_rules": len(V), "blocks": nblock}
        for r in range(len(H) - 1):
            for c in range(len(V) - 1):
                y0, y1 = H[r][1] + 1, H[r + 1][0]
                x0, x1 = V[c][1] + 1, V[c + 1][0]
                if y1 - y0 < 4:
                    continue
                gl, ul = components(ink, y0, y1, x0, x1)
                for k, g in enumerate(gl):
                    out.append({"table": name, "row": r, "col": c, "pos": k,
                                "h": g["h"], "w": g["w"], "ymin": g["ymin"],
                                "ymax": g["ymax"], "bitmap": g["bitmap"],
                                "native": g["native"], "underlines": ul})
    return out, meta


def cluster(glyphs):
    X = np.array([g["bitmap"].ravel() for g in glyphs])
    D = pdist(X, metric="euclidean")
    Z = linkage(D, method="complete")
    lab = fcluster(Z, t=LINKAGE_CUTOFF, criterion="distance")
    S = squareform(D)
    stats = {}
    for cid in sorted(set(lab.tolist())):
        idx = np.where(lab == cid)[0]
        intra = S[np.ix_(idx, idx)].max() if len(idx) > 1 else 0.0
        other = np.setdiff1d(np.arange(len(glyphs)), idx)
        inter = S[np.ix_(idx, other)].min() if len(other) else float("inf")
        stats[int(cid)] = {"n": int(len(idx)), "max_intra": float(intra),
                           "min_inter": float(inter),
                           "rep": int(idx[0])}
    return lab, stats


def fingerprint(bitmap):
    q = (np.round(bitmap * 255)).astype("uint8")
    return hashlib.sha256(q.tobytes()).hexdigest()[:16]


def verify_merged(glyphs, lab, stats):
    """Re-derive every MULTI-CHARACTER cluster's label from single-char parts.

    At this resolution some characters TOUCH and so form one connected
    component -- 28 of the 45 clusters are two or three characters. Those
    labels were read off a montage by eye, which puts one person on both sides
    of the comparison. This re-derives each of them by a route that uses only
    the SINGLE-character templates (each of which came from an isolated,
    unambiguous component) and pixel arithmetic: a left-to-right dynamic
    programme places templates, allowing 0-3 px of overlap where the ink of
    adjacent glyphs merged, and keeps the string with the fewest mismatched
    pixels. It never sees the eye label.

    Reports best string, its residual, and the runner-up's residual, so a
    reader can see the margin rather than a bare 'agrees'.
    """
    labels = json.load(open(os.path.join(HERE, "04_glyph_labels.json")))
    by_fp = labels["labels"]
    singles = []
    for cid, s in stats.items():
        fp = fingerprint(glyphs[s["rep"]]["bitmap"])
        ch = by_fp.get(fp)
        if ch is not None and len(ch) == 1:
            singles.append((ch, glyphs[s["rep"]]["native"]))
    print(f"single-character templates: {''.join(sorted(c for c, _ in singles))}")

    out = []
    for cid, s in sorted(stats.items()):
        fp = fingerprint(glyphs[s["rep"]]["bitmap"])
        eye = by_fp[fp]
        if len(eye) < 2:
            continue
        comp = glyphs[s["rep"]]["native"]
        H_, W_ = comp.shape
        best = {}          # x -> (cost, string)
        best[0] = (0, "")
        for x in range(W_ + 1):
            if x not in best:
                continue
            cost, sofar = best[x]
            for ch, tpl in singles:
                th, tw = tpl.shape
                if th != H_ or x + tw > W_ + 3:
                    continue
                win = np.zeros((H_, tw), dtype=bool)
                xs = min(tw, W_ - x)
                if xs <= 0:
                    continue
                win[:, :xs] = comp[:, x:x + xs]
                c = int((win ^ tpl).sum())
                for ov in (0, 1, 2, 3):
                    nx = x + tw - ov
                    if nx > W_ + 1 or nx <= x:
                        continue
                    cand = (cost + c + ov, sofar + ch)
                    if nx not in best or cand[0] < best[nx][0]:
                        best[nx] = cand
        ends = [best[x] for x in (W_ - 1, W_, W_ + 1) if x in best]
        ends.sort()
        if not ends:
            out.append({"cluster": cid, "eye": eye, "derived": None,
                        "agrees": False, "note": "no decomposition reached "
                                                 "the component's right edge"})
            continue
        derived, resid = ends[0][1], ends[0][0]
        runner = next((e for e in ends[1:] if e[1] != derived), None)
        out.append({"cluster": cid, "n": s["n"], "eye": eye,
                    "derived": derived, "residual_px": resid,
                    "runner_up": runner[1] if runner else None,
                    "runner_up_residual": runner[0] if runner else None,
                    "ink_px": int(comp.sum()),
                    "agrees": derived == eye})
    bad = [o for o in out if not o["agrees"]]
    for o in out:
        flag = "OK " if o["agrees"] else "XX "
        print(f"{flag}cluster {o['cluster']:>2} n={o.get('n'):>3} "
              f"eye={o['eye']:<4} derived={o['derived']!s:<4} "
              f"residual={o.get('residual_px')}/{o.get('ink_px')}")
    path = os.path.join(HERE, "04_merged_glyph_verification.json")
    with open(path, "w") as fh:
        json.dump({"_what": "Every multi-character glyph cluster, re-derived "
                            "from the single-character templates alone by a "
                            "dynamic programme that never sees the eye label.",
                   "n_multi": len(out), "n_disagreeing": len(bad),
                   "results": out}, fh, indent=1)
    print(f"\n{len(out)} multi-character clusters, {len(bad)} disagree "
          f"-> {os.path.basename(path)}")
    return 1 if bad else 0


def main(argv):
    glyphs, meta = harvest()
    lab, stats = cluster(glyphs)
    print(f"{len(glyphs)} glyphs, {len(stats)} clusters")
    worst = max(s["max_intra"] for s in stats.values())
    tight = min(s["min_inter"] for s in stats.values())
    print(f"max intra-cluster distance {worst:.3f}; "
          f"min inter-cluster distance {tight:.3f}")

    # Within a cluster every NORMALISED bitmap is identical. That alone does
    # not prove the members are the same character: two different shapes could
    # in principle coincide after an aspect-preserving resize. So the native
    # (h, w) must also agree, and a disagreement raises.
    for cid, s in stats.items():
        idx = [i for i, l in enumerate(lab) if l == cid]
        hw = {(glyphs[i]["h"], glyphs[i]["w"]) for i in idx}
        if len(hw) > 1:
            raise SystemExit(
                f"cluster {cid}: members share a normalised bitmap but have "
                f"different native sizes {sorted(hw)}. Normalisation has "
                "merged two different shapes; raise NORM or tighten the cut.")

    if "--montage" in argv:
        reps = sorted(stats.items(), key=lambda kv: -kv[1]["n"])
        cols = 6
        scale = 10
        pad = 10
        cellw = max(glyphs[s["rep"]]["w"] for _, s in reps) * scale + 3 * pad
        cellh = max(glyphs[s["rep"]]["h"] for _, s in reps) * scale + 3 * pad
        rows = (len(reps) + cols - 1) // cols
        sheet = Image.new("L", (cols * cellw, rows * cellh), 255)
        for i, (cid, s) in enumerate(reps):
            nb = glyphs[s["rep"]]["native"]
            img = Image.fromarray((255 - nb * 255).astype("uint8")).resize(
                (nb.shape[1] * scale, nb.shape[0] * scale), Image.NEAREST)
            sheet.paste(img, ((i % cols) * cellw + pad,
                              (i // cols) * cellh + pad))
        out = os.path.join(HERE, "00_pages", "glyph_montage.png")
        sheet.save(out)
        order = [{"slot": i, "cluster": cid, "n": s["n"],
                  "native_h": glyphs[s["rep"]]["h"],
                  "native_w": glyphs[s["rep"]]["w"],
                  "fingerprint": fingerprint(glyphs[s["rep"]]["bitmap"])}
                 for i, (cid, s) in enumerate(reps)]
        with open(os.path.join(HERE, "04_glyph_clusters.json"), "w") as fh:
            json.dump({"_what": "cluster -> montage slot, with the sha256 of "
                                "the representative's normalised bitmap. The "
                                "label file is keyed on this fingerprint so a "
                                "label can never drift onto a different shape.",
                       "cols": cols, "cell_px": [cellw, cellh],
                       "max_intra": worst, "min_inter": tight,
                       "clusters": order}, fh, indent=1)
        print(f"wrote {out} ({rows} rows x {cols} cols)")
        return 0

    if "--verify" in argv:
        return verify_merged(glyphs, lab, stats)

    if "--read" not in argv:
        print("run with --montage to produce the sheet to label, then "
              "--read once 04_glyph_labels.json exists.")
        return 0

    lf = os.path.join(HERE, "04_glyph_labels.json")
    if not os.path.exists(lf):
        raise SystemExit("04_glyph_labels.json missing: label the montage.")
    labels = json.load(open(lf))["labels"]   # fingerprint -> character

    # --- the gate: every cluster's representative must be labelled, BY
    # --- FINGERPRINT. An unknown fingerprint raises; it never defaults.
    decoded = {}
    for cid, s in stats.items():
        fp = fingerprint(glyphs[s["rep"]]["bitmap"])
        if fp not in labels:
            raise SystemExit(
                f"cluster {cid} (n={s['n']}) has fingerprint {fp}, which is "
                "not in 04_glyph_labels.json. The clustering has changed, or "
                "a cluster was never labelled. Refusing to decode: a glyph "
                "that falls through to a default is how 1.2x10^5 becomes 109.")
        decoded[cid] = labels[fp]

    cells = {}
    for g, cid in zip(glyphs, lab):
        key = (g["table"], g["row"], g["col"])
        cells.setdefault(key, []).append(decoded[int(cid)])
    text = {f"{t}|{r}|{c}": "".join(v) for (t, r, c), v in cells.items()}

    with open(os.path.join(HERE, "04_glyphs.json"), "w") as fh:
        json.dump({"_what": "every table cell, decoded from clustered glyph "
                            "bitmaps. Key is table|row|col in the image's own "
                            "ruled-line grid.",
                   "meta": meta,
                   "n_glyphs": len(glyphs), "n_clusters": len(stats),
                   "max_intra_cluster_distance": worst,
                   "min_inter_cluster_distance": tight,
                   "cells": text}, fh, indent=1)
    print(f"wrote 04_glyphs.json: {len(text)} non-empty cells")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
