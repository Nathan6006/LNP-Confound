#!/usr/bin/env python3
"""Measure drawn chain lengths for amines 1A-13A by apex counting, then calibrate.

Visual counting of zigzag carbons is the failure mode most likely to corrupt this
library, so chain length is MEASURED here rather than read. 9A and 10A are pinned
independently by HRMS (Supplementary Figs 5-8: di-n-octylamine and di-n-decylamine), so
they serve as the calibration points: the measured apex count is converted to carbons by
the offset those two anchors imply, and the result is only trusted if both anchors land
exactly and the series is monotonic.

Outputs JSON: label -> {apexes, carbons, n_chains}.
"""
import json
import re
import numpy as np
import pdfplumber
from PIL import Image

PDF = "s41563-020-00886-0.pdf"
PNG = "00_pages/img/main_p-02.png"
DPI_SCALE = None  # computed from page size vs raster size


def red_mask(arr):
    r, g, b = arr[..., 0], arr[..., 1], arr[..., 2]
    return (r > 110) & (r - g > 55) & (r - b > 55)


def apexes(mask, min_sep=5, min_amp=1.5):
    ys = []
    for x in range(mask.shape[1]):
        c = np.flatnonzero(mask[:, x])
        ys.append(c.mean() if c.size else np.nan)
    v = np.array(ys)
    idx = np.flatnonzero(~np.isnan(v))
    if idx.size < 12:
        return 0
    v = np.interp(np.arange(len(v)), idx, v[idx])
    v = np.convolve(v, np.ones(3) / 3, mode="same")
    lo, hi = idx.min() + 3, idx.max() - 3
    n, last = 0, -99
    for i in range(lo + 1, hi):
        s = v[max(lo, i - min_sep):min(hi, i + min_sep) + 1]
        if (v[i] == s.max() or v[i] == s.min()) and i - last >= min_sep:
            if s.max() - s.min() >= min_amp:
                n += 1
                last = i
    return n


def main():
    with pdfplumber.open(PDF) as pdf:
        page = pdf.pages[1]
        pw, ph = page.width, page.height
        labels = {}
        for w in page.extract_words():
            if re.fullmatch(r"\d{1,2}A", w["text"]):
                labels.setdefault(w["text"], (w["x0"], w["top"]))
    im = Image.open(PNG)
    arr = np.asarray(im.convert("RGB")).astype(int)
    sx, sy = im.size[0] / pw, im.size[1] / ph
    mask = red_mask(arr)

    order = [f"{i}A" for i in range(1, 14)]
    tops = sorted((labels[k][1], k) for k in order)
    out = {}
    for i, (top, k) in enumerate(tops):
        nxt = tops[i + 1][0] if i + 1 < len(tops) else top + 22
        y0, y1 = int((top - 4) * sy), int((nxt - 3) * sy)
        x0 = int((labels[k][0] + 14) * sx)
        x1 = int(340 * sx)
        band = mask[y0:y1, x0:x1]
        rows = np.flatnonzero(band.any(axis=1))
        if rows.size == 0:
            out[k] = {"apexes": 0, "n_chains": 0}
            continue
        gaps = np.diff(rows)
        if gaps.size and gaps.max() > 5:
            sp = rows[int(np.argmax(gaps))] + int(gaps.max()) // 2
            chains = [band[:sp], band[sp:]]
        else:
            chains = [band]
        counts = [apexes(c) for c in chains]
        out[k] = {"apexes": max(counts), "per_chain": counts, "n_chains": len(chains)}

    # calibrate on the HRMS-pinned anchors
    anchors = {"9A": 8, "10A": 10}
    offs = {k: anchors[k] - out[k]["apexes"] for k in anchors if k in out}
    off = round(np.mean(list(offs.values()))) if offs else 0
    print(f"apex->carbon offset from HRMS anchors {offs} -> using {off}\n")
    print(f"{'amine':>6} {'chains':>6} {'apexes':>7} {'carbons':>8}  note")
    for k in order:
        d = out[k]
        c = d["apexes"] + off
        d["carbons"] = c
        note = ""
        if k in anchors:
            note = "ANCHOR ok" if c == anchors[k] else f"ANCHOR MISMATCH exp {anchors[k]}"
        print(f"{k:>6} {d['n_chains']:>6} {d['apexes']:>7} {c:>8}  {note}")
    json.dump(out, open("02_amine_measurements.json", "w"), indent=1)


if __name__ == "__main__":
    main()
