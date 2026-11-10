#!/usr/bin/env python3
"""Figure 3C -> 08_pies.csv : organ fractions, as a PRINTED-STATISTIC check.

The Results section prints one number this panel can be held to:

    "514O6,10 produced 35% of the signal in the lung and spleen each, with only
     30% in the liver."

Three printed percentages for one drawn pie. Nothing here is emitted into the
corpus -- Figure 3A's stacked log bars are not decomposable (see
00_value_sources.json) and these pies are fractions, not fluxes. This stage
exists only so that a statistic the authors printed is reproduced by
measurement rather than asserted.

HOW THE FRACTION IS MEASURED. Each pie is a filled disc in three flat colours.
The fraction of a disc a sector covers is its ANGLE over 360 degrees, and for a
sector of a disc that equals its AREA fraction exactly -- so counting classified
pixels inside the disc is the measurement, with no angle-finding and no centre
estimate needed. The disc is found as a connected component of non-white pixels;
its black outline and the gaps between sectors are excluded by classifying only
pixels within 40 RGB units of one of the three legend colours, and the
unclassified fraction is reported so a bad classification cannot hide.

THE LEGEND COLOURS ARE SAMPLED FROM THIS SAME RASTER, from the legend swatches
beside panel C, never named from the caption.

Run: python3 -I stage07_fig3c_pie.py
"""
import csv
import json
import os

import numpy as np
from PIL import Image
from scipy import ndimage

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "source", "nihms-2006250-f0003.jpg")

BAND = (1350, 1976)      # rows containing panel C's three pies
ORDER = ["514O6,10", "503Oi10", "503O8,12"]
PRINTED = {"514O6,10": {"Lungs": 35.0, "Spleen": 35.0, "Liver": 30.0}}
TOL_RGB = 40.0


def legend_swatches(rgb, band):
    """The three legend swatch colours, sampled from this raster.

    The legend sits to the RIGHT of the pies; its swatches are the three small
    filled squares. They are found as the three non-white components in the
    right margin whose bounding box is roughly square and small.
    """
    sub = rgb[band[0]:band[1], :]
    nonwhite = sub.mean(axis=2) < 235
    lab, n = ndimage.label(nonwhite)
    objs = ndimage.find_objects(lab)
    out = []
    for i, sl in enumerate(objs, start=1):
        h = sl[0].stop - sl[0].start
        w = sl[1].stop - sl[1].start
        # Measured: the swatches are 37 x 62 RECTANGLES, not squares. The first
        # attempt required squareness and found none.
        if not (25 <= h <= 60 and 40 <= w <= 90):
            continue
        if sl[1].start < sub.shape[1] * 0.6:
            continue
        m = (lab[sl] == i)
        px = sub[sl][m]
        med = np.median(px, axis=0)
        # A legend swatch is COLOURED. Black text glyphs pass the size and
        # squareness tests -- the first attempt returned (18,18,18),
        # (53,53,53), (54,54,54), i.e. three letters, and classified 96% of
        # every pie as unclassified. A grey channel spread of 0 is the tell.
        if float(med.max() - med.min()) < 20.0:
            continue
        out.append((sl[0].start, med))
    out.sort(key=lambda t: t[0])
    if len(out) < 3:
        raise SystemExit(f"legend_swatches: found {len(out)}, need 3.")
    return [c for _, c in out[:3]]


def main():
    img = Image.open(FIG).convert("RGB")
    rgb = np.asarray(img).astype(float)
    sw = legend_swatches(rgb, BAND)
    names = ["Spleen", "Lungs", "Liver"]      # legend order, top to bottom

    sub = rgb[BAND[0]:BAND[1], :]

    # THE SECTOR COLOURS ARE NOT THE LEGEND COLOURS, and that had to be
    # measured rather than assumed. The spleen swatch is (124, 172, 212) and
    # the spleen SECTOR is (124, 175, 254) -- 42 RGB units apart, so a 40-unit
    # match against the swatch classified 0.78% of the first pie as spleen
    # against a printed 35%, and the three fractions summed to 62.6% instead
    # of 100. Raising the tolerance would have been the wrong repair: it would
    # have hidden the fact that the two are different colours.
    #
    # So the three SECTOR colours are taken from the discs themselves, as the
    # three dominant quantised colours, and are then LABELLED by nearest legend
    # swatch. The labelling must be injective and unambiguous, and that is
    # asserted below rather than hoped for.
    pie_band = sub.mean(axis=2) < 235
    q = (sub // 16 * 16).astype(int)
    keys, counts = np.unique(q[pie_band].reshape(-1, 3), axis=0,
                             return_counts=True)
    top = keys[np.argsort(counts)[::-1][:3]]
    sector = []
    for k in top:
        m = pie_band & np.all(q == k[None, None, :], axis=2)
        sector.append(np.median(sub[m], axis=0))
    DD = np.array([[float(np.linalg.norm(c - w)) for w in sw] for c in sector])
    from scipy.optimize import linear_sum_assignment as _lsa
    ri, ci = _lsa(DD)
    for r, c in zip(ri, ci):
        others = sorted(DD[r][j] for j in range(3) if j != c)
        if DD[r][c] >= others[0]:
            raise SystemExit(
                f"fig3c: sector colour {sector[r].round(0).tolist()} is not "
                f"unambiguously nearest one legend swatch "
                f"({DD[r].round(1).tolist()}). Refusing to label it.")
    sw_ordered = [None] * 3
    for r, c in zip(ri, ci):
        sw_ordered[c] = sector[r]
    sector_vs_swatch = {k: round(float(DD[r][c]), 1)
                        for (r, c), k in zip(sorted(zip(ri, ci),
                                                    key=lambda t: t[1]),
                                             ["Spleen", "Lungs", "Liver"])}

    D = np.stack([np.linalg.norm(sub - c[None, None, :], axis=2)
                  for c in sw_ordered])
    near = D.argmin(axis=0)
    ok = D.min(axis=0) < TOL_RGB
    classified = np.where(ok, near, -1)

    # The discs: connected components of NON-WHITE, not of sector colour. The
    # sectors are separated by black rules, so a mask built from sector colour
    # alone splits one pie into three pieces -- the first attempt found 0 discs
    # for exactly that reason.
    nonwhite = sub.mean(axis=2) < 235
    nonwhite = ndimage.binary_closing(nonwhite, np.ones((5, 5)))
    lab, n = ndimage.label(nonwhite)
    sizes = ndimage.sum(nonwhite, lab, range(1, n + 1))
    discs = []
    for i, s in enumerate(sizes, start=1):
        if s < 20000:
            continue
        sl = ndimage.find_objects(lab == i)[0]
        h = sl[0].stop - sl[0].start
        w = sl[1].stop - sl[1].start
        if abs(h - w) > 0.12 * max(h, w):
            continue
        discs.append((sl[1].start, i, sl, int(s)))
    discs.sort()
    if len(discs) != len(ORDER):
        raise SystemExit(f"fig3c: {len(discs)} discs, expected {len(ORDER)}: "
                         f"{[(d[0], d[3]) for d in discs]}")

    rows, qc = [], {"legend_swatch_rgb": {names[k]: [round(float(v), 1)
                                                     for v in sw[k]]
                                          for k in range(3)},
                    "sector_rgb_measured_in_the_pies": {
                        names[k]: [round(float(v), 1) for v in sw_ordered[k]]
                        for k in range(3)},
                    "sector_to_swatch_distance_rgb": sector_vs_swatch,
                    "tol_rgb": TOL_RGB}
    for (x0, i, sl, size), lipid in zip(discs, ORDER):
        m = (lab == i)
        tot = int(m.sum())
        px_frac = {nm: 100.0 * float(((classified == k) & m).sum()) / tot
                   for k, nm in enumerate(names)}

        # THE MEASUREMENT: angular span, not pixel count.
        # A pie sector's share IS its angle, and an angular measurement is not
        # biased by the black outline and the sector dividers, which cover
        # ~6% of the disc and sit exactly where two sectors meet. Rays are cast
        # at 3600 angles from the disc centroid at 0.6 R; a ray that lands on a
        # divider is unclassified and is resolved by its two angular
        # neighbours, which must agree -- if they do not, the ray is dropped
        # and counted, so an unreadable pie cannot pass as a clean one.
        cy, cx = ndimage.center_of_mass(m)
        R = 0.5 * (sl[0].stop - sl[0].start)
        ang = np.linspace(0, 2 * np.pi, 3600, endpoint=False)
        yy = np.rint(cy + 0.6 * R * np.sin(ang)).astype(int)
        xx = np.rint(cx + 0.6 * R * np.cos(ang)).astype(int)
        ray = classified[yy, xx].astype(int)
        unres = 0
        fixed = ray.copy()
        for j in np.flatnonzero(ray < 0):
            a, b = ray[(j - 1) % 3600], ray[(j + 1) % 3600]
            if a >= 0 and a == b:
                fixed[j] = a
            else:
                unres += 1
        good = fixed >= 0
        frac = {nm: 100.0 * float((fixed[good] == k).sum()) / int(good.sum())
                for k, nm in enumerate(names)}
        rows.append({"lipidoid": lipid, "disc_px": tot,
                     "height_px": sl[0].stop - sl[0].start,
                     "width_px": sl[1].stop - sl[1].start,
                     "rays": 3600, "rays_unresolved": int(unres),
                     **{f"{k}_pct": round(v, 2) for k, v in frac.items()},
                     "sum_pct": round(sum(frac.values()), 3),
                     **{f"{k}_pct_by_pixel_count": round(v, 2)
                        for k, v in px_frac.items()}})
    with open(os.path.join(HERE, "08_pies.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow(r)

    chk = {}
    for lipid, want in PRINTED.items():
        got = [r for r in rows if r["lipidoid"] == lipid][0]
        chk[lipid] = {k: {"printed_pct": v,
                          "measured_pct": got[f"{k}_pct"],
                          "delta_pp": round(got[f"{k}_pct"] - v, 2)}
                      for k, v in want.items()}
    qc["printed_statistic_check"] = chk
    qc["rows"] = rows
    qc["_what_this_does_not_prove"] = (
        "A pie reproduces a FRACTION. It says nothing about the absolute "
        "fluxes in Figure 3A, which are stacked on a log axis and are not "
        "decomposed anywhere in this folder.")
    with open(os.path.join(HERE, "08_pies_qc.json"), "w") as fh:
        json.dump(qc, fh, indent=1, sort_keys=True)
    print(json.dumps(qc, indent=1, sort_keys=True))


if __name__ == "__main__":
    main()
