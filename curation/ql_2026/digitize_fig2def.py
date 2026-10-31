#!/usr/bin/env python3
"""Digitize Figure 2d/2e/2f, the printed hit-rate bars, as an INDEPENDENT check
on the Figure 2b heatmap reading.

These panels are not new measurements - they summarise the same 320 lipids - so
they are not curated as rows. They are worth digitizing because they are an
arithmetic consequence of the heatmap, which makes them a hard test of it:
a hit is RLU > 2e4, and every bar must land on an exact k/40 or k/64.

Decoding 2e is what makes the check sharp. Its three bars read 12.5 / 15.0 /
62.5 %, which look nothing like the group hit rates (160, 80 and 80 lipids).
They resolve exactly when the denominator is a flat 40 - the caption's "n = 40":
62.5% = 25/40, and 25 is precisely T3A's 17 hits plus 110's 8 from panel 2d.
So 2d and 2e are two views of the same counts, and both must be reproduced.

Bars are hollow and colour-outlined, so a bar top is found as the highest row
where that bar's x-window holds a run of saturated non-grey pixels.
"""
import csv
import json
import os

import numpy as np
import PIL.Image as Image
import PIL.ImageDraw as ImageDraw

PAGE = "00_pages/fig_p4_600.png"
REVIEW = "review"

# panel -> (x range to search, y range, axis max, bar labels)
PANELS = {
    "2d": dict(x=(386, 1780), y=(3600, 4700), ymax=50.0,
               labels=["70", "75", "80", "400", "113", "114", "110", "T3A"]),
    "2e": dict(x=(2074, 2720), y=(3600, 4700), ymax=80.0,
               labels=["1", "2", "3"]),
    "2f": dict(x=(3024, 3900), y=(3600, 4700), ymax=25.0,
               labels=["C6", "C8", "C10", "C12", "C14"]),
}


def find_axes(g, xr, yr):
    """Baseline = the long horizontal dark rule; y-axis = the long vertical one."""
    sub = g[yr[0]:yr[1], xr[0]:xr[1]]
    dark = sub < 140
    rows = dark.mean(1)
    base = max((y for y in range(len(rows)) if rows[y] > 0.45),
               key=lambda y: (rows[y], y))
    cols = dark[:base].mean(0)
    axis = max((x for x in range(len(cols)) if cols[x] > 0.45), default=None,
               key=lambda x: cols[x])
    return yr[0] + base, xr[0] + axis


def main():
    os.makedirs(REVIEW, exist_ok=True)
    img = np.array(Image.open(PAGE).convert("RGB")).astype(int)
    g = img.mean(2)
    sat = img.max(2) - img.min(2)
    coloured = (sat > 60) | (g < 110)

    out, boxes = [], []
    for pname, P in PANELS.items():
        base, axis = find_axes(g, P["x"], P["y"])
        # tick values on the y axis give the scale: find the topmost tick label
        # rule by locating the axis top (highest dark pixel on the axis column)
        colmask = g[P["y"][0]:base, axis - 2:axis + 3].min(1) < 140
        top = P["y"][0] + int(np.argmax(colmask))
        px_per_unit = (base - top) / P["ymax"]

        # Bars are HOLLOW colour-outlined rectangles, so probing a band just
        # above the baseline finds outline strokes, not solid fills, and an
        # earlier version merged all eight into one run. Instead build the
        # TOP-EDGE PROFILE: for each column, the highest coloured pixel above
        # the baseline. A bar shows up as a plateau, the gaps between bars as
        # columns with nothing at all. This also handles the zero-height bar
        # (amine 70), which simply has no plateau.
        lo, hi = top - 40, base - 15
        prof = np.full(P["x"][1] - P["x"][0], -1)
        for i, x in enumerate(range(P["x"][0], P["x"][1])):
            colm = coloured[lo:hi, x]
            w = np.nonzero(colm)[0]
            if len(w) >= 3:
                prof[i] = lo + w[0]
        runs, cur = [], None
        for i, v in enumerate(prof):
            if v >= 0 and cur is None:
                cur = i
            elif v < 0 and cur is not None:
                if i - cur > 20:
                    runs.append((P["x"][0] + cur, P["x"][0] + i))
                cur = None
        if cur is not None and len(prof) - cur > 20:
            runs.append((P["x"][0] + cur, P["x"][0] + len(prof)))
        runs = [r for r in runs if r[0] > axis + 6]

        # a bar of height zero leaves no plateau; the panels label it anyway
        if len(runs) == len(P["labels"]) - 1 and P["labels"][0] == "70":
            runs = [(axis + 10, runs[0][0] - 10)] + runs

        assert len(runs) == len(P["labels"]), \
            f"{pname}: found {len(runs)} bars, expected {len(P['labels'])}: {runs}"

        for (x0, x1), lab in zip(runs, P["labels"]):
            seg = prof[x0 - P["x"][0]:x1 - P["x"][0]]
            seg = seg[seg >= 0]
            ytop = int(np.median(seg)) if len(seg) else base
            value = max((base - ytop) / px_per_unit, 0.0)
            out.append({"panel": pname, "label": lab,
                        "value_pct": round(value, 3),
                        "x0": x0, "x1": x1, "y_top": ytop, "y_base": base})
            boxes.append((x0, ytop, x1, base))

    with open("03_fig2def_bars.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0]))
        w.writeheader()
        w.writerows(out)

    # snap to the exact fractions the panels must be made of
    def snap(v, n):
        k = round(v * n / 100)
        return k, round(100 * k / n, 4), abs(v - 100 * k / n)

    snapped = {}
    for r in out:
        n = 64 if r["panel"] == "2f" else 40
        k, exact, err = snap(r["value_pct"], n)
        snapped[f"{r['panel']}:{r['label']}"] = {
            "read_pct": r["value_pct"], "n": n, "implied_count": k,
            "exact_pct": exact, "abs_error_pct": round(err, 3)}

    qc = {"panels": "Fig 2d (by amine, n=40), 2e (by tail number, n=40 flat), "
                    "2f (by tail length, n=64)",
          "hit_definition": "RLU > 2e4",
          "bars": snapped,
          "max_snap_error_pct": round(max(v["abs_error_pct"]
                                          for v in snapped.values()), 3),
          "implied_total_hits_from_2d": sum(v["implied_count"]
                                            for k, v in snapped.items()
                                            if k.startswith("2d")),
          "implied_total_hits_from_2e": sum(v["implied_count"]
                                            for k, v in snapped.items()
                                            if k.startswith("2e")),
          "implied_total_hits_from_2f": sum(v["implied_count"]
                                            for k, v in snapped.items()
                                            if k.startswith("2f"))}
    json.dump(qc, open("04_fig2def_qc.json", "w"), indent=1)

    crop_overlay(img, boxes)
    for k, v in snapped.items():
        print(f"  {k:10} read {v['read_pct']:6.2f}%  -> {v['implied_count']:2}/{v['n']}"
              f" = {v['exact_pct']:6.2f}%  (err {v['abs_error_pct']:.2f})")
    print(f"\nmax snap error {qc['max_snap_error_pct']}%")
    print(f"implied total hits: 2d={qc['implied_total_hits_from_2d']}  "
          f"2e={qc['implied_total_hits_from_2e']}  2f={qc['implied_total_hits_from_2f']}")
    print(f"overlay -> {REVIEW}/fig2def_overlay.png")


def crop_overlay(img, boxes):
    y0 = min(b[1] for b in boxes) - 120
    y1 = max(b[3] for b in boxes) + 60
    crop = Image.fromarray(img.astype("uint8")[y0:y1, 200:3900]).copy()
    d = ImageDraw.Draw(crop)
    for x0, yt, x1, yb in boxes:
        d.rectangle([x0 - 200, yt - y0, x1 - 200, yb - y0],
                    outline=(220, 0, 0), width=4)
    d.text((8, 8), "red = measured bar extent (top edge = value)", fill=(0, 0, 0))
    crop.save(f"{REVIEW}/fig2def_overlay.png")


if __name__ == "__main__":
    main()
