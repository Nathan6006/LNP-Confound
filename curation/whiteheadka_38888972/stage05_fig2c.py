#!/usr/bin/env python3
"""Figure 2C (top) -> 07_fig2c_points.csv, and THE CROSS-PANEL SET CHECK.

WHY THIS STAGE EXISTS. Figures 2A and 2B were digitised independently, by
different code, off different panels. Figure 2C (top) is a THIRD panel in which
the authors plotted the same 49 LNPs as (TNS fluorescence at pH 5, in-vivo
luminescent flux) pairs. So the point cloud here must be the same SET of pairs
we assembled from 2A and 2B.

That is the strongest check available in this paper, for three reasons:
  * it needs no attribution -- the points are unlabelled, and a SET comparison
    does not care which point is which;
  * it is external to both digitisers, so a uniform scale error in either one
    shows up as a systematic displacement rather than being absorbed;
  * it is a comparison against something the AUTHORS drew, not against a
    statistic we recomputed from our own numbers. PIPELINE_V2: "passing every
    printed statistic is necessary, not sufficient."

The matching is a minimum-cost assignment in (TNS AU, log10 flux) space with
scipy.optimize.linear_sum_assignment, so no point is used twice and no greedy
order effect can flatter the result.

This panel's values are NEVER emitted. 00_value_sources.json records why: the
points carry no lipidoid identity, and inventing one is the AC_2025 failure.

Run: python3 -I stage05_fig2c.py   (after stage02 and stage04)
"""
import csv
import json
import os
import sys

import numpy as np
from PIL import Image
from scipy import ndimage
from scipy.optimize import linear_sum_assignment

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))

import axis                                                     # noqa: E402

FIG = os.path.join(HERE, "source", "nihms-2006250-f0002.jpg")

# `top` is the 10^10 tick row + 5. The first attempt used 70 and swept in the
# panel TITLE as 24 extra "markers", all of them above 10^10 -- which is how
# it was caught: a point above the axis maximum is not a point.
PANEL = {"x_axis_px": 1461, "baseline": 523, "top": 125, "right": 2090}
Y_TICK_STRIP = (1448, 1459)
Y_TICK_ROWS = (60, 528)
Y_TICK_VALUES = [1e10, 1e9, 1e8, 1e7, 1e6, 1e5]
X_TICK_STRIP = (527, 538)
X_TICK_VALUES = [0.0, 15000.0, 30000.0]

# A single marker measures ~210 px. A blob materially larger is a marker fused
# with a dashed-rule segment or with a neighbouring marker, and this panel has
# both: the markers OVERLAP, so the 49 points cannot all be recovered as
# distinct centroids. That is reported, not worked around -- a merged blob is
# excluded from the comparison rather than split by an area quotient, which
# would be a guess presented as a measurement.
MIN_BLOB_PX = 140
SINGLE_MAX_PX = 260
MAX_BLOB_PX = 2000


def main():
    img = Image.open(FIG).convert("RGB")
    rgb = np.asarray(img).astype(float)
    gray = rgb.mean(axis=2)
    qc = {}

    # --- axes, from detected ticks only ---------------------------------
    yt = [t + Y_TICK_ROWS[0] for t in axis.detect_ticks_1d(
        gray[Y_TICK_ROWS[0]:Y_TICK_ROWS[1]], "y", *Y_TICK_STRIP,
        dark=140, min_len=2, max_len=25)]
    if len(yt) != len(Y_TICK_VALUES):
        raise SystemExit(f"fig2c: {len(yt)} y ticks: {np.round(yt, 2).tolist()}")
    caly = axis.calibrate(yt, Y_TICK_VALUES, kind="log10", axis="y")

    g = gray.copy()
    g[:, :1452] = 255
    g[:, PANEL["right"]:] = 255
    xt = axis.detect_ticks_1d(g, "x", *X_TICK_STRIP, dark=140, min_len=2,
                              max_len=30)
    if len(xt) != len(X_TICK_VALUES):
        raise SystemExit(f"fig2c: {len(xt)} x ticks: {np.round(xt, 2).tolist()}")
    calx = axis.calibrate(xt, X_TICK_VALUES, kind="linear", axis="x")
    qc["y_axis"] = caly.qc()
    qc["x_axis"] = calx.qc()
    qc["y_ticks_px"] = [round(t, 3) for t in yt]
    qc["x_ticks_px"] = [round(t, 3) for t in xt]

    # --- the markers ----------------------------------------------------
    sub = rgb[PANEL["top"]:PANEL["baseline"] - 1,
              PANEL["x_axis_px"] + 2:PANEL["right"]]
    sg = sub.mean(axis=2)
    ink = sg < 190
    lab, n = ndimage.label(ink)
    sizes = ndimage.sum(ink, lab, range(1, n + 1))
    keep = [i + 1 for i, s in enumerate(sizes)
            if MIN_BLOB_PX <= s <= SINGLE_MAX_PX]
    merged = [int(s) for s in sizes if SINGLE_MAX_PX < s <= MAX_BLOB_PX]
    cents = ndimage.center_of_mass(ink, lab, keep)
    qc["blobs_total"] = int(n)
    qc["blobs_kept"] = len(keep)
    qc["blob_size_px"] = {"kept_min": int(min(sizes[np.array(keep) - 1])),
                          "kept_max": int(max(sizes[np.array(keep) - 1]))}
    qc["blobs_merged_sizes_px"] = sorted(merged)
    qc["estimated_markers_in_merged_blobs"] = int(
        round(sum(merged) / 210.0)) if merged else 0
    qc["_marker_recovery"] = (
        "The panel's markers overlap, so not all 49 appear as separate "
        "connected components. Only unambiguous single markers enter the "
        "comparison; the merged blobs are counted and left out.")

    pts = []
    for (cy, cx), k in zip(cents, keep):
        y = cy + PANEL["top"]
        x = cx + PANEL["x_axis_px"] + 2
        pts.append({
            "x_px": round(float(x), 2), "y_px": round(float(y), 2),
            "area_px": int(sizes[k - 1]),
            "TNS_AU": round(float(calx.to_value(x)), 1),
            "log10_flux": round(float(np.log10(caly.to_value(y))), 5),
        })
    pts.sort(key=lambda p: (p["TNS_AU"], p["log10_flux"]))
    with open(os.path.join(HERE, "07_fig2c_points.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(pts[0].keys()))
        w.writeheader()
        for p in pts:
            w.writerow(p)

    # --- the set comparison ---------------------------------------------
    tns = {r["IL_name"]: float(r["TNS_AU"])
           for r in csv.DictReader(open(os.path.join(HERE, "06_tns.csv")))}
    bars = [r for r in csv.DictReader(
        open(os.path.join(HERE, "04_fig2b_bars.csv")))
        if r["is_control"] == "False"]
    ours = [(r["IL_name"], tns[r["IL_name"]], float(r["log10_flux"]))
            for r in bars if r["IL_name"] in tns]
    qc["n_lipidoids_from_2A_and_2B"] = len(ours)
    qc["n_points_in_2C"] = len(pts)

    if len(pts) and len(ours):
        A = np.array([[t, f] for _, t, f in ours])
        B = np.array([[p["TNS_AU"], p["log10_flux"]] for p in pts])
        # Scale each axis to its own range so neither dominates the cost.
        sc = np.array([max(float(np.ptp(A[:, 0])), 1.0), max(float(np.ptp(A[:, 1])), 1e-6)])
        C = np.linalg.norm(A[:, None, :] / sc - B[None, :, :] / sc, axis=2)
        # Rectangular: every detected marker is matched to a DISTINCT one of
        # our 49, and the unmatched remainder of our 49 are the points the
        # panel draws underneath another marker.
        ri, ci = linear_sum_assignment(C.T if C.shape[0] > C.shape[1] else C)
        if C.shape[0] > C.shape[1]:
            ri, ci = ci, ri
        d_tns = B[ci, 0] - A[ri, 0]
        d_log = B[ci, 1] - A[ri, 1]
        qc["set_match"] = {
            "n_matched": int(len(ri)),
            "TNS_AU_delta": {
                "mean": round(float(d_tns.mean()), 1),
                "median": round(float(np.median(d_tns)), 1),
                "rms": round(float(np.sqrt((d_tns ** 2).mean())), 1),
                "max_abs": round(float(np.abs(d_tns).max()), 1)},
            "log10_flux_delta": {
                "mean": round(float(d_log.mean()), 4),
                "median": round(float(np.median(d_log)), 4),
                "rms": round(float(np.sqrt((d_log ** 2).mean())), 4),
                "max_abs": round(float(np.abs(d_log).max()), 4)},
            "_reading": "A mean delta near zero on BOTH axes means neither "
                        "digitiser carries a uniform scale or origin error. "
                        "The rms is the per-point agreement and is bounded "
                        "below by the panels' own drawing precision: a 1 px "
                        "marker-centroid error is ~68 AU on x and ~0.012 "
                        "decades on y.",
        }
    with open(os.path.join(HERE, "07_fig2c_qc.json"), "w") as fh:
        json.dump(qc, fh, indent=1, sort_keys=True)
    print(json.dumps(qc, indent=1, sort_keys=True))


if __name__ == "__main__":
    main()
