#!/usr/bin/env python3
"""Stage 04 - Figure 3b, the only individually-measured DELIVERY panel.

Total spleen luminescence (p/s) vs total liver luminescence (p/s). NanoLuc
mRNA, C57BL/6, i.v., 0.1 mg/kg encapsulated mRNA, 6 h. Eight groups:
PBS + LNP 9 / 37 / 94 / 96 / 97 / 112 + SM-102.

TIER 5 (pixels). Why tiers 1-4 fail, each a checkable sentence:
  1 publisher Source Data - Wiley published none. The JATS declares exactly one
    supplementary item, ADMA-38-e14370-s001.pdf = source/adma...SuppMat.pdf.
  2 SI tables - the SI has four tables (barcodes x2, formulations, antibody
    clones); none carries a luminescence value (06_checks.py::si_value_tables).
  3 author repository - no GEO/SRA/PRJNA/Zenodo/Dryad/GitHub/figshare string
    occurs in the article XML or the 22,296-char SI text; Data Availability
    reads "available from the corresponding author upon reasonable request";
    provenance/tools/find_repo.py returned no candidate.
  4 pdf text table - Figure 3 exists only as a flattened JPEG. No PDF rendition
    of the ARTICLE is obtainable: PMC /bin/ serves an HTML error page,
    europepmc fulltextRepo answers HTTP 403 {"error":"PDF link h..."}, and
    onlinelibrary.wiley.com returns HTTP 403 to a plain client.

EVERY POSITION COMES FROM DETECTED FURNITURE. Decade ticks are detected in a
strip OUTSIDE each rule and, on y, separated from the 5x minor ticks by MEASURED
TICK LENGTH (decade 4 px, minor 3 px, 12/12 alternating) rather than by taking
every other one. The printed decade labels say only WHICH decade a detected tick
is, never where one sits. That is the class-20 boundary.

THE PANEL'S OWN PRINTED "18.1x" IS THE GATE, AND IT HAS ALREADY EARNED ITSELF.
The first version of this file took the mean of each colour mask. The grey mask
swallowed the 131-px-wide axis rule, every colour also matched its own LEGEND
swatch, and brown additionally matched the dark flank of the orange marker it
overlaps. That version reported 6.685x against a printed 18.1x. Three repairs at
the source -- connected components instead of a mask mean, winner-take-all
nearest-colour assignment so two series cannot claim one pixel, and a saturation
gate that drops anti-aliased black text -- move it to 18.17x, a 0.03 px error.
"""
import json, math, os, sys
import numpy as np
from PIL import Image
from scipy import ndimage

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import axis                                              # noqa: E402

FIG = os.path.join(HERE, "source", "figs", "ADMA-38-e14370-g002.jpg")
OUT = os.path.join(HERE, "04_fig3b_measurements.json")

# Marker fills, as printed in the panel's own legend.
SERIES = {
    "LNP 9": (227, 30, 36), "LNP 37": (50, 110, 185), "LNP 94": (55, 165, 75),
    "LNP 96": (125, 60, 152), "LNP 97": (245, 145, 30),
    "LNP 112": (255, 242, 0), "SM-102": (165, 75, 40),
}
# PBS is deliberately absent: its fill is neutral grey, which at this
# resolution is the same colour as the panel's anti-aliased black text. All 30
# of its candidate blobs have saturation 0. It is a vehicle control carrying no
# ionizable lipid, so it is never a corpus row; the refusal costs nothing and
# is recorded rather than worked around.
COLOUR_TOL = 70          # max distance from a series fill for a pixel to vote
MIN_SAT = 40             # max(RGB)-min(RGB); anti-aliased black text is ~0
MARKER_PX = (6, 40)
MARKER_SIDE = (3, 8)
MIN_MARKER_SEP = 3.0     # px between two DATA markers' centres
DECADE_TICK_LEN = 4      # measured; the 5x minor ticks are 3 px


def rules(gray):
    dark = gray < 140
    sub = dark[0:200, 330:546]
    cols = [i + 330 for i, c in enumerate(sub.sum(0)) if c > 100]
    rows = [i for i, r in enumerate(sub.sum(1)) if r > 80]
    if len(cols) != 1 or len(rows) != 1:
        raise SystemExit("panel rules not unique: cols=%s rows=%s" % (cols, rows))
    return cols[0], rows[0]


def tick_len(gray, row, x_rule, dark=210, span=16):
    n = 0
    for x in range(x_rule - 1, x_rule - span, -1):
        rows = [row] + ([row - 1] if row > 0 else []) + \
               ([row + 1] if row + 1 < gray.shape[0] else [])
        if any(gray[r, x] < dark for r in rows):
            n += 1
        else:
            break
    return n


def fill_missed(ticks):
    gaps = list(np.diff(ticks))
    med = float(np.median(gaps))
    out, added = [], 0
    for i, v in enumerate(ticks):
        out.append(v)
        if i + 1 < len(ticks):
            n = int(round(gaps[i] / med))
            for k in range(1, n):
                out.append(v + gaps[i] * k / n)
                added += 1
    return sorted(out), added


def assign(img, x_rule, y_rule):
    """Winner-take-all: each pixel votes for its NEAREST series fill only.

    Two series claiming one pixel is what let brown annex the dark flank of the
    orange marker it overlaps, producing a second spurious SM-102 blob 24 px
    from the real one -- on the very axis the 18.1x bracket is drawn along.
    """
    names = list(SERIES)
    d = np.stack([np.sqrt(((img - np.array(SERIES[n])) ** 2).sum(2))
                  for n in names])
    best = d.argmin(0)
    ok = d.min(0) < COLOUR_TOL
    sat = img.max(2) - img.min(2)
    ok &= sat >= MIN_SAT
    ok[:, :x_rule + 1] = False
    ok[y_rule + 1:, :] = False
    return {n: (ok & (best == i)) for i, n in enumerate(names)}


def blobs(mask):
    lab, k = ndimage.label(mask)
    out = []
    for i in range(1, k + 1):
        ys, xs = np.nonzero(lab == i)
        w, h = int(xs.max() - xs.min() + 1), int(ys.max() - ys.min() + 1)
        if not (MARKER_PX[0] <= len(xs) <= MARKER_PX[1]):
            continue
        if not (MARKER_SIDE[0] <= w <= MARKER_SIDE[1]
                and MARKER_SIDE[0] <= h <= MARKER_SIDE[1]):
            continue
        out.append({"n_px": int(len(xs)), "x": float(xs.mean()),
                    "y": float(ys.mean()), "w": w, "h": h})
    return out


def main(seed=None):
    img = np.array(Image.open(FIG).convert("RGB")).astype(int)
    gray = np.array(Image.open(FIG).convert("L")).astype(float)
    x_rule, y_rule = rules(gray)
    rep = {"_panel": "Figure 3b", "figure_px": [img.shape[1], img.shape[0]],
           "rules": {"x_rule_px": x_rule, "y_rule_px": y_rule},
           "seed": seed}

    # ---- x axis (liver) ---------------------------------------------
    xt = [v for v in axis.detect_ticks_1d(gray, along="x", scan_lo=y_rule + 3,
                                          scan_hi=y_rule + 5, dark=150,
                                          min_len=1, max_len=99)
          if x_rule < v < 546]
    xt, x_added = fill_missed(xt)
    calx = axis.calibrate(xt, [10.0 ** (5 + i) for i in range(len(xt))],
                          kind="log10", axis="x")
    rep["x_axis"] = {"decade_ticks_px": [round(v, 2) for v in xt],
                     "recovered_from_a_double_gap": x_added,
                     "px_per_decade": round(calx.px_per_unit, 4),
                     "residual_px": round(calx.residual_px, 4),
                     "loo_residuals_px": [round(r, 3) for r in axis.loo_residuals(xt)]}

    # ---- y axis (spleen): decade ticks separated by MEASURED LENGTH --
    cand = [v for v in axis.detect_ticks_1d(gray, along="y", scan_lo=402,
                                            scan_hi=405, dark=200, min_len=1,
                                            max_len=99) if v < y_rule]
    lens = [tick_len(gray, int(round(v)), x_rule) for v in cand]
    yt = [v for v, L in zip(cand, lens) if L >= DECADE_TICK_LEN]
    rep["y_axis_tick_lengths"] = [[round(v, 1), L] for v, L in zip(cand, lens)]
    caly = axis.calibrate(yt, [10.0 ** (9 - i) for i in range(len(yt))],
                          kind="log10", axis="y")
    rep["y_axis"] = {"decade_ticks_px": [round(v, 2) for v in yt],
                     "px_per_decade": round(abs(caly.px_per_unit), 4),
                     "residual_px": round(caly.residual_px, 4),
                     "loo_residuals_px": [round(r, 3) for r in axis.loo_residuals(yt)]}

    ppd = calx.px_per_unit
    rep["precision"] = {
        "px_per_decade_liver_x": round(ppd, 4),
        "px_per_decade_spleen_y": round(abs(caly.px_per_unit), 4),
        "value_factor_per_px_x": round(10 ** (1.0 / ppd), 4),
        "library_tick_gate_px": axis.TICK_RESIDUAL_PX,
        "value_factor_that_gate_permits_x": round(
            10 ** (axis.TICK_RESIDUAL_PX / ppd), 3),
        "why_this_matters": (
            "axis.TICK_RESIDUAL_PX = 2.5 px was set for 600 dpi renders, where "
            "it is ~0.1 pt. A decade on this x axis is %.2f px, so the SAME "
            "gate would pass a factor-%.2f error in a reported value. The gate "
            "is not wrong; it is being asked to protect a panel it was not "
            "scaled for, which is why the printed 18.1x is the real gate here."
            % (ppd, 10 ** (axis.TICK_RESIDUAL_PX / ppd))),
    }

    # ---- markers ------------------------------------------------------
    masks = assign(img, x_rule, y_rule)
    allb = {n: blobs(m) for n, m in masks.items()}
    from collections import Counter
    c = Counter(round(b["x"]) for bl in allb.values() for b in bl)
    legx = {x for x, n in c.items() if n >= 3}
    rep["legend_columns_px_detected"] = sorted(legx)

    # Two DATA markers cannot have centres closer than a marker radius. The
    # dark anti-aliased flank of the orange LNP 97 marker is nearer brown than
    # orange, so it survives winner-take-all as a 6 px SM-102 blob whose centre
    # is 1.63 px from LNP 97's -- inside LNP 97's own 4 px footprint. The
    # smaller of any such pair is dropped, and the drop is recorded.
    flat = [(n, b) for n, bl in allb.items() for b in bl]
    swallowed = []
    for n, b in flat:
        for m, o in flat:
            if m == n and o is b:
                continue
            d = math.hypot(b["x"] - o["x"], b["y"] - o["y"])
            if d < MIN_MARKER_SEP and o["n_px"] > b["n_px"]:
                swallowed.append({"series": n, "x": round(b["x"], 2),
                                  "y": round(b["y"], 2), "n_px": b["n_px"],
                                  "inside": m, "centre_sep_px": round(d, 2)})
                break
    drop = {(s_["series"], s_["x"], s_["y"]) for s_ in swallowed}
    rep["blobs_inside_another_marker"] = swallowed

    pts, refused = {}, {}
    for n, bl in allb.items():
        keep = [b for b in bl if round(b["x"]) not in legx
                and (n, round(b["x"], 2), round(b["y"], 2)) not in drop]
        if seed == "shift_sm102" and n == "SM-102":
            for b in keep:
                b["x"] -= 4.0
        if len(keep) != 1:
            refused[n] = {"n_blobs": len(keep), "blobs": keep}
            continue
        b = keep[0]
        pts[n] = {"x_px": round(b["x"], 2), "y_px": round(b["y"], 2),
                  "n_px": b["n_px"], "w": b["w"], "h": b["h"],
                  "liver_p_per_s": float(calx.to_value(b["x"])),
                  "spleen_p_per_s": float(caly.to_value(b["y"]))}
    rep["markers"] = pts
    rep["markers_refused"] = refused
    rep["PBS"] = ("not detected: neutral grey is the colour of the panel's "
                  "anti-aliased black text at 546 px; all 30 candidate blobs "
                  "have saturation 0. PBS carries no ionizable lipid and is "
                  "never a corpus row.")

    # ---- THE GATE: the printed "18.1x" --------------------------------
    chk = {"printed": 18.1}
    if "LNP 112" in pts and "SM-102" in pts:
        dpx = pts["SM-102"]["x_px"] - pts["LNP 112"]["x_px"]
        meas = 10 ** (dpx / ppd)
        chk.update({"separation_px": round(dpx, 3),
                    "px_predicted_by_18.1x": round(ppd * math.log10(18.1), 3),
                    "px_error": round(dpx - ppd * math.log10(18.1), 3),
                    "measured_ratio": round(float(meas), 3),
                    "ratio_of_ratios": round(float(meas / 18.1), 4),
                    "PASS": bool(abs(dpx - ppd * math.log10(18.1)) <= 1.0)})
    else:
        chk["NOT_COMPUTABLE"] = "a marker was refused"
        chk["PASS"] = False
    rep["printed_check_18.1x"] = chk

    if seed is None:
        with open(OUT, "w") as fh:
            json.dump(rep, fh, indent=1)
    print(json.dumps(rep, indent=1))
    return 0 if chk.get("PASS") else 1


if __name__ == "__main__":
    sys.exit(main(seed=(sys.argv[1] if len(sys.argv) > 1 else None)))
