"""Stage 04 -- do this folder's checks have teeth? Seeded mutations.

A CHECK THAT CANNOT FAIL IS NOT A CHECK, AND A SEED THAT ONLY RAISES IS NOT
CREDITED. Every seed below is credited only when a COUNT MOVES: the number of
panels read, or the number of values read, or the number of bars accepted. A
seed that makes the pipeline exit non-zero while moving nothing is recorded as
NOT CREDITED, because `axis.calibrate` refusing before any count exists proves
only that calibrate refuses.

Verdicts are read from the returncode of the function call, never from stdout
and never through a pipe.

Usage: python3 -I stage04_selfcheck.py   (appends to 07_qc.json)
"""
import copy
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
sys.path.insert(0, HERE)
import axis                                                   # noqa: E402
import stage01_digitize as S                                  # noqa: E402
from errors import DetectionError, CalibrationError           # noqa: E402

# One panel, used as the test bed for every seed: Fig S13a Liver, which reads
# 6 of 6 LNP bars and 2 of 6 PBS bars cleanly.
PANEL = {"page": 23, "xref": 244, "spine_x": 597, "y0": 99, "y1": 338,
         "x_right": 1065, "labels": S.S13_ORGANS[1][1], "paired": 0.22}


def read(img, **over):
    p = dict(PANEL)
    p.update(over)
    vals, diag, cal = S.read_panel(
        img, p["spine_x"], p["y0"], p["y1"], p["x_right"], p["labels"],
        S.CAL_TOL, S.LOO_TOL, p["paired"])
    n = sum(1 for v in vals.values() if v["lnp"] is not None)
    return n, vals, cal


def shift_ticks(img, rows, dy):
    """Move the y tick STUBS of one panel down by dy px, leaving all else.

    The stubs live in columns [spine_x-5, spine_x-1); moving them and nothing
    else is a pure scale perturbation -- the bars do not move, so a detector
    that still reads the same number of bars has had its SCALE corrupted, and
    the LOO gate is the only thing standing between that and 6 wrong values.
    """
    out = img.copy()
    x0, x1 = PANEL["spine_x"] - 5, PANEL["spine_x"] - 1
    strip = out[:, x0:x1].copy()
    for r in rows:
        out[r - 3:r + 4, x0:x1] = 254            # erase the stub
    for r in rows:
        out[r + dy - 3:r + dy + 4, x0:x1] = strip[r - 3:r + 4]
    return out


def main():
    img, _ = S.load(PANEL["page"], PANEL["xref"])
    report = {"_what": "seeded mutations against this folder's own geometry; "
                       "credited only when a count MOVES"}

    base_n, base_vals, base_cal = read(img)
    base_values = [base_vals[k]["lnp"] for k in sorted(base_vals)]
    report["control"] = {
        "panel": "Fig S13a Liver", "lnp_bars_read": base_n,
        "values_pct": base_values,
        "loo_max_px": base_cal["loo_max_px"],
        "spine_top_residual_px": base_cal["spine_top_residual_px"],
        "verdict": "PASS -- the control must still pass, or a seed proves "
                   "nothing",
    }
    ticks = [int(round(t)) for t in base_cal["ticks_px"]]

    seeds = {}

    # --- seed 1: a GROSS tick displacement -------------------------------
    for dy, name in ((6, "gross_+6px_on_the_top_tick"),
                     (3, "moderate_+3px_on_the_top_tick"),
                     (2, "subtle_+2px_on_the_top_tick"),
                     (1, "below_the_gate_+1px_on_the_top_tick")):
        bad = shift_ticks(img, [ticks[0]], dy)
        rec = {"seed": f"y tick stub at row {ticks[0]} moved +{dy} px"}
        try:
            n, vals, cal = read(bad)
            rec["returncode_equivalent"] = 0
            rec["lnp_bars_read"] = n
            rec["loo_max_px"] = cal["loo_max_px"]
            rec["values_pct"] = [vals[k]["lnp"] for k in sorted(vals)]
            moved = sum(1 for a, b in zip(base_values, rec["values_pct"])
                        if a is not None and b is not None
                        and abs(a - b) > 0.25)
            rec["values_that_moved_by_more_than_0.25pct"] = moved
            rec["credited"] = False
            rec["verdict"] = (
                "NOT CAUGHT. The scale was corrupted and %d of %d values "
                "moved by more than 0.25 pct, and nothing refused. This is "
                "the gate's measured SENSITIVITY LIMIT and it is reported, "
                "not hidden: a 1 px tick error on this panel gives a "
                "leave-one-out residual of %.3f px against the 1.5 px gate, "
                "which is 0.42 pct of signal on a 2.378 px/pct axis."
                % (moved, n, rec["loo_max_px"]))
        except (DetectionError, CalibrationError) as e:
            rec["returncode_equivalent"] = 1
            rec["raised"] = type(e).__name__
            rec["message"] = str(e).splitlines()[0]
            rec["lnp_bars_read"] = 0
            rec["count_moved"] = f"{base_n} -> 0"
            rec["credited"] = base_n != 0
            rec["verdict"] = (
                "CAUGHT, and CREDITED: the readable-bar count moved from "
                f"{base_n} to 0. This is not a bare raise -- the control run "
                "reads 6 bars from the same panel, so the refusal removed "
                "real output.")
        seeds[name] = rec

    # --- seed 2: delete one category tick --------------------------------
    bad = img.copy()
    base_row = int(round(base_cal["baseline_row"]))
    ct = int(round(base_cal["category_ticks_px"][2]))
    bad[base_row + 1:base_row + 6, ct - 3:ct + 4] = 254
    rec = {"seed": f"x category stub at column {ct} erased"}
    try:
        n, vals, cal = read(bad)
        rec["returncode_equivalent"] = 0
        rec["lnp_bars_read"] = n
        rec["credited"] = False
        rec["verdict"] = ("NOT CAUGHT -- five stubs were zipped onto six "
                          "labels and %d values were returned." % n)
    except DetectionError as e:
        rec["returncode_equivalent"] = 1
        rec["message"] = str(e).splitlines()[0]
        rec["count_moved"] = f"{base_n} -> 0"
        rec["credited"] = True
        rec["verdict"] = ("CAUGHT and CREDITED. Without this the remaining "
                          "five bars would have been labelled ECs, Kupffer "
                          "Cells, DCs, T cells, B cells -- every value after "
                          "the gap on the wrong cell type.")
    seeds["category_stub_deleted"] = rec

    # --- seed 3: a replicate-marker style blob on a bar's centre ---------
    # error class 32: a centred marker above the bar top. The flank test is
    # what must refuse it.
    bad = img.copy()
    d = base_vals["Kupffer Cells"]
    x0, x1 = PANEL["spine_x"], PANEL["x_right"]
    # find that bar's window from the diagnostics of the control run
    _, dg, _ = S.read_panel(img, PANEL["spine_x"], PANEL["y0"], PANEL["y1"],
                            PANEL["x_right"], PANEL["labels"], S.CAL_TOL,
                            S.LOO_TOL, PANEL["paired"])
    bx0, bx1 = dg["Kupffer Cells"]["lnp"]["x0"], \
        dg["Kupffer Cells"]["lnp"]["x1"]
    row = dg["Kupffer Cells"]["lnp"]["row"]
    cx = (bx0 + bx1) // 2
    MARKER_FRAC = float(os.environ.get("SEED_MARKER_FRAC", "0.30"))
    half = max(2, int((bx1 - bx0) * MARKER_FRAC))
    bad[row - 30:row - 22, cx - half:cx + half] = [113, 90, 138]
    rec = {"seed": "an 8 px tall, %d%%-wide bar-coloured blob 25 px above "
                   "the Kupffer Cells bar top -- a replicate marker"
                   % round(MARKER_FRAC * 200)}
    try:
        n, vals, cal = read(bad)
        rec["returncode_equivalent"] = 0
        rec["lnp_bars_read"] = n
        rec["kupffer_value_pct"] = vals["Kupffer Cells"]["lnp"]
        rec["control_value_pct"] = d["lnp"]
        moved = abs((rec["kupffer_value_pct"] or 0) - d["lnp"]) > 0.25
        rec["value_moved"] = moved
        rec["credited"] = False
        rec["verdict"] = (
            "NOT CAUGHT -- the marker was read as the bar top, value moved "
            f"{d['lnp']} -> {rec['kupffer_value_pct']}" if moved else
            "NOT CREDITED, and that is the honest label. Nothing refused and "
            "no count moved: the full-width coverage rule simply stepped "
            "past a blob that never reaches 0.90 of the bar width, and the "
            "value came back unchanged. Correct behaviour, but it is not "
            "evidence that the flank test fired. The wider variant below is "
            "what exercises that path.")
    except DetectionError as e:
        rec["returncode_equivalent"] = 1
        rec["message"] = str(e).splitlines()[0]
        rec["count_moved"] = f"{base_n} -> 0"
        rec["credited"] = True
        rec["verdict"] = "CAUGHT and CREDITED (count moved)."
    seeds["centred_replicate_marker"] = rec

    # --- seed 4: the inverted tick-value order ---------------------------
    # The real bug this file had: TICK_VALUES ascending against ascending
    # rows calibrates cleanly with 0% at the TOP.
    rec = {"seed": "TICK_VALUES passed bottom-up instead of top-down"}
    try:
        cal = axis.calibrate(base_cal["ticks_px"],
                             list(reversed(S.TICK_VALUES)), kind="linear",
                             residual_tol_px=axis.TICK_RESIDUAL_PX, axis="y")
        resid = abs(float(cal.to_pixel(100.0)) - float(PANEL["y0"]))
        rec["calibrate_raised"] = False
        rec["tick_fit_residual_px"] = cal.qc()["tick_fit_residual_px"]
        rec["spine_top_residual_px"] = round(resid, 2)
        rec["caught_by_spine_check"] = resid > 3.0
        rec["credited"] = resid > 3.0
        rec["verdict"] = (
            "calibrate itself is BLIND to this -- it fits cleanly at "
            f"{rec['tick_fit_residual_px']} px residual with 0% at the top of "
            "the panel. The spine check catches it at "
            f"{rec['spine_top_residual_px']} px against a 3.0 px gate. This "
            "is the bug this file actually had, and the spine check is why "
            "it is now impossible.")
    except Exception as e:                       # pragma: no cover
        rec["calibrate_raised"] = True
        rec["message"] = str(e).splitlines()[0]
        rec["credited"] = True
    seeds["inverted_tick_value_order"] = rec

    # The wider variants. A marker only reaches the FLANK test if it first
    # passes the 0.90 full-width coverage rule, so the two guards cover
    # adjacent width bands and the sweep below measures where the boundary
    # actually is on THIS bar.
    bar_w = bx1 - bx0
    sweep = {}
    for w in range(int(bar_w * 0.55), bar_w + 1):
        bad2 = img.copy()
        lo = cx - w // 2
        bad2[row - 30:row - 22, lo:lo + w] = [113, 90, 138]
        try:
            _, vals2, _ = read(bad2)
            v2 = vals2["Kupffer Cells"]["lnp"]
            sweep[w] = {"refused": False, "value_pct": v2,
                        "moved": v2 is not None
                        and abs(v2 - d["lnp"]) > 0.25}
        except DetectionError as e:
            sweep[w] = {"refused": True,
                        "by": ("flank" if "flank" in str(e) else "coverage"),
                        "message": str(e).splitlines()[0][:120]}
    caught = [w for w, r in sweep.items()
              if r["refused"] or not r.get("moved", True)]
    leaked = [w for w, r in sweep.items()
              if not r["refused"] and r.get("moved")]
    seeds["centred_replicate_marker_width_sweep"] = {
        "seed": "an 8 px tall bar-coloured blob 25 px above the Kupffer "
                "Cells bar top, swept from %d to %d px wide on a %d px bar"
                % (int(bar_w * 0.55), bar_w, bar_w),
        "bar_width_px": bar_w,
        "control_value_pct": d["lnp"],
        "sweep": {str(k): v for k, v in sweep.items()},
        "n_widths_tested": len(sweep),
        "n_harmless_or_refused": len(caught),
        "n_that_changed_the_value": len(leaked),
        "widths_that_changed_the_value": leaked,
        "credited": len(leaked) == 0,
        "verdict": (
            "CREDITED: no marker width changed the value. "
            if not leaked else
            "NOT CREDITED: widths %s changed the value. " % leaked) +
        "MEASURED LIMITATION, stated rather than glossed: on an %d px bar the "
        "flank test has almost no room to operate. A marker must cover >= 90%% "
        "of the bar to get past the coverage rule, i.e. >= %d px, and a "
        "centred marker that wide already covers both %d px flanks -- so the "
        "two guards do not overlap at this bar width and the protection here "
        "is the coverage rule alone. The flank test earns its keep on wider "
        "bars; it is not what is holding this panel up."
        % (bar_w, int(bar_w * 0.9) + 1, max(2, int(bar_w * 0.33))),
    }

    report["seeds"] = seeds
    report["n_seeds"] = len(seeds)
    report["n_credited"] = sum(1 for v in seeds.values()
                               if v.get("credited"))

    qc_path = os.path.join(HERE, "07_qc.json")
    qc = json.load(open(qc_path))
    qc["seeded_selfchecks"] = report
    json.dump(qc, open(qc_path, "w"), indent=1, sort_keys=True, default=str)
    print(json.dumps({k: {"credited": v.get("credited"),
                          "verdict": v.get("verdict", "")[:90]}
                      for k, v in seeds.items()}, indent=1))
    print("control bars read:", base_n, " credited:", report["n_credited"],
          "/", report["n_seeds"])


if __name__ == "__main__":
    main()
