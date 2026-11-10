#!/usr/bin/env python3
"""Stage 05a -- Figures S24C, S24E and S24F, from the PDF's OWN vector geometry.

    python3 -I stage05a_si_vector.py

NO PIXEL IS SAMPLED. The SI's bar charts are native PDF drawings: every bar is
a filled rectangle and every tick is a stroked line, both carried as exact
coordinates in PDF user space. So there is no rasterisation, no JPEG, no
antialiasing and no colour matching -- and in particular none of the hazards
that produced F-018 (a bar-top detector reading an overlaid replicate dot),
F-019 (a colorbar anchored on the outer edge of its own border) or F-020.

The scale still comes from axis furniture only (curation/lib/axis.calibrate),
fitted on the FOUR MAJOR TICKS and nothing else. The panel's four MINOR ticks
are held out of the fit and then predicted: they must land on 15/45/75/105.
That is an invariant the authors DREW and never printed, external to this
pipeline in the way CLAUDE.md s5 asks for, and it is what this stage trusts
more than any recomputed marginal.
"""
import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
sys.path.insert(0, HERE)

import pymupdf                                                  # noqa: E402
import numpy as np                                              # noqa: E402

import axis                                                     # noqa: E402

PDF = os.path.join(HERE, "source",
                   "anie202004994-sup-0001-misc_information.pdf")
PAGE = 36                                      # 0-based; printed page S36

# Panels, identified by the y of their own baseline rule (which is also the
# zero tick). Nothing here is a literal scale or origin -- see F-033.
# `x_range` is required, not optional: panels E and F share a baseline to
# within 0.137 pt, so a y-only selector silently merged them (48 bars where
# there are 12) on this file's first run.
PANELS = {
    "S24C": {"baseline": 340.231, "axis_x": 140.154, "n_series": 3,
             "x_range": (140.0, 292.0),
             "series": ["25 nM", "50 nM", "100 nM"],
             "unit": "percent_GFP_positive_cells",
             "what": "(-30)GFP-Cre protein delivery to HeLa-DsRed"},
    "S24E": {"baseline": 510.621, "axis_x": 139.558, "n_series": 1,
             "x_range": (139.0, 292.0),
             "series": ["unstated"],
             "unit": "mean_green_fluorescence_au",
             "what": "mean green fluorescence intensity, HeLa-DsRed"},
    "S24F": {"baseline": 510.758, "axis_x": 336.813, "n_series": 3,
             "x_range": (336.0, 488.0),
             "series": ["25 nM", "50 nM", "100 nM"],
             "unit": "percent_DsRed_positive_cells",
             "what": "Cre-mediated gene recombination, HeLa-DsRed"},
}
MAJOR_LEN = 2.18        # pt, measured; minor ticks are 1.09
TOL = 0.05    # pt. E and F differ by 0.137, so this must stay well under it.


def page():
    return pymupdf.open(PDF)[PAGE]


def ticks(p, axis_x, y_lo, y_hi):
    """(major, minor) tick y-positions on the vertical rule at axis_x."""
    maj, mnr = [], []
    for g in p.get_drawings():
        if g["type"] != "s":
            continue
        for it in g["items"]:
            if it[0] != "l":
                continue
            a, b = it[1], it[2]
            if abs(a.y - b.y) > 0.05 or abs(a.x - axis_x) > 0.05:
                continue
            ln = abs(b.x - a.x)
            if ln < 0.3:
                continue
            (maj if ln > (MAJOR_LEN + 1.09) / 2 else mnr).append(round(a.y, 3))
    return sorted(set(maj)), sorted(set(mnr))


def bars(p, baseline, x_range):
    """Filled rectangles sitting on `baseline` inside `x_range`, left to right."""
    lo, hi = x_range
    out = []
    for g in p.get_drawings():
        if g["type"] != "fs" or not any(i[0] == "re" for i in g["items"]):
            continue
        r = g["rect"]
        if abs(r.y1 - baseline) > TOL or r.x0 < lo or r.x1 > hi:
            continue
        out.append({"x0": r.x0, "x1": r.x1, "top": r.y0,
                    "fill": tuple(round(c, 3) for c in (g.get("fill") or ()))})
    out.sort(key=lambda b: b["x0"])
    return out


def labels(p, y_centre, x_lo, x_hi):
    """Rotated category labels under a panel, left to right."""
    out = []
    for blk in p.get_text("dict")["blocks"]:
        for line in blk.get("lines", []):
            for s in line["spans"]:
                t = s["text"].strip()
                if not t:
                    continue
                x, y = s["bbox"][0], s["bbox"][1]
                if abs(y - y_centre) < 1.5 and x_lo <= x <= x_hi:
                    out.append((x, t))
    out.sort()
    return [t for _x, t in out]


def main():
    p = page()
    qc = {}
    rows = []
    for pid, cfg in PANELS.items():
        base = cfg["baseline"]
        maj, mnr = ticks(p, cfg["axis_x"], base - 120, base + 1)
        bb = bars(p, base, cfg["x_range"])
        # drop legend swatches: they do not sit on the baseline, but guard
        # anyway by requiring the count to factor into n_series.
        if len(bb) % cfg["n_series"]:
            raise SystemExit(f"{pid}: {len(bb)} bars is not a multiple of "
                             f"{cfg['n_series']} series.")
        ngroup = len(bb) // cfg["n_series"]

        # --- calibration: MAJOR ticks only --------------------------------
        vals = [800, 600, 400, 200, 0] if pid == "S24E" else [90, 60, 30, 0]
        if len(maj) != len(vals):
            raise SystemExit(f"{pid}: found {len(maj)} major ticks, "
                             f"expected {len(vals)}: {maj}")
        cal = axis.calibrate(maj, vals, kind="linear", axis="y",
                             forbid_target=[b["top"] for b in bb][:len(maj)])
        loo = axis.loo_residuals(maj)
        # --- the held-out invariant: the MINOR ticks ----------------------
        # The minors took no part in the fit. Each must land on an ODD multiple
        # of half the major step -- i.e. exactly half way between two majors.
        # The expected value is DERIVED from the fitted step, not listed by
        # hand: a hand-written list is how this check first reported a 30.0
        # error, because panel F has three minors and panel C has four and the
        # two lists were zipped out of alignment.
        half = abs(vals[0] - vals[1]) / 2.0
        pred = [float(cal.to_value(m)) for m in sorted(mnr)]
        minor_expect = [round(v / half) * half for v in pred]
        odd = [int(round(v / half)) % 2 == 1 for v in minor_expect]
        if not all(odd):
            raise SystemExit(f"{pid}: a held-out minor tick did not land "
                             f"half way between majors: {pred} -> "
                             f"{minor_expect}")
        minor_err = [abs(a - b) for a, b in zip(pred, minor_expect)]
        qc[pid] = {
            "what": cfg["what"],
            "n_bars": len(bb), "n_groups": ngroup,
            "major_ticks_px": maj, "major_tick_values": vals,
            "calibration": cal.qc(),
            "leave_one_out_residuals_px": [round(float(v), 4) for v in loo],
            "max_leave_one_out_residual_px": round(float(max(abs(v) for v in loo)), 4),
            "held_out_minor_ticks_px": sorted(mnr),
            "held_out_minor_predicted": [round(v, 3) for v in pred],
            "held_out_minor_expected": minor_expect,
            "held_out_minor_are_odd_half_steps": all(odd),
            "held_out_minor_max_abs_error": round(max(minor_err), 4)
            if minor_err else None,
        }
        # --- labels --------------------------------------------------------
        y_lab = {"S24C": 341.3, "S24E": 510.8, "S24F": 510.9}[pid]
        x_lo, x_hi = cfg["x_range"][0] - 20, cfg["x_range"][1]
        labs = labels(p, y_lab, x_lo, x_hi)
        if len(labs) != ngroup:
            raise SystemExit(f"{pid}: {len(labs)} labels for {ngroup} groups: "
                             f"{labs}")
        qc[pid]["labels"] = labs
        # --- values ---------------------------------------------------------
        for gi in range(ngroup):
            for si in range(cfg["n_series"]):
                b = bb[gi * cfg["n_series"] + si]
                rows.append({
                    "panel": pid, "sample": labs[gi],
                    "series": cfg["series"][si],
                    "unit": cfg["unit"],
                    "value": round(float(cal.to_value(b["top"])), 3),
                    "bar_top_px": round(b["top"], 3),
                    "bar_x0_px": round(b["x0"], 3),
                    "fill_rgb": str(b["fill"]),
                })
    out = os.path.join(HERE, "05a_si_vector_values.csv")
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    with open(os.path.join(HERE, "05a_si_vector_qc.json"), "w") as fh:
        json.dump(qc, fh, indent=1)
    print(f"wrote {out}: {len(rows)} values")
    for pid in PANELS:
        q = qc[pid]
        print(f"  {pid}: {q['n_bars']} bars / {q['n_groups']} groups; "
              f"tick-fit residual {q["calibration"]["tick_fit_residual_px"]:.4f} px; "
              f"max LOO {q['max_leave_one_out_residual_px']:.4f} px; "
              f"held-out minor ticks max err {q['held_out_minor_max_abs_error']}")


if __name__ == "__main__":
    main()
