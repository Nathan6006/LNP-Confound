#!/usr/bin/env python3
"""Digitise SI Figure S2A -- RAW luminescence total flux (p/s) per organ for the
eight single-LNP validation formulations (13, 14, 70, 96, 97, 98, 45, 102).

TIER 5 (pixels). Tiers 1-4 fail for a checkable reason: there is no publisher
Source Data and no author repository (00_value_sources.json), and SI p.14 carries
page.get_drawings() == 0 with exactly ONE embedded raster, so no vector geometry
and no text layer exist inside the figure.

NATIVE RESOLUTION IS MEASURED, NOT GUESSED: 2134 px / (432.0 pt) * 72 = 355.7 dpi.
curation/lib/axis's tolerances are written for a 600 dpi render, so they are
rescaled with axis.native_tol and BOTH figures are reported. The widened
tolerance is SEEDED in this panel's own geometry (check P8), because the
module's own docstring says the margin numbers do not transfer between panels.

SCALE COMES FROM AXIS FURNITURE ONLY -- the seven detected decade tick marks
protruding left of the spine. `forbid_target` is passed so that fitting the
scale onto the quantity being validated is a type error (error class 27).

PANEL B IS NOT DIGITISED. It is DiR fluorescence = particle biodistribution.
Uptake is not delivery.
"""
import json, os, sys
import numpy as np
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import axis, bars                                   # noqa: E402
from errors import CalibrationError, DetectionError  # noqa: E402

# Two panels are read. They carry THE SAME 48 measurements in two
# parameterisations -- S2A is raw total flux, S1A is the same values as a
# percentage of each animal's organ total -- which is what makes the
# cross-panel identity in check P10/P11 a real external test rather than a
# restatement.
#
# THE COLOUR LUT IS SAMPLED PER PANEL, NEVER SHARED. S1 and S2 are separate
# rasters with separate JPEG histories and their fills differ by up to 14 per
# channel (e.g. Brain 249,150,157 on S1 vs 235,155,159 on S2). Reusing one
# panel's LUT on the other is the cross-panel colour hazard the library warns
# about; here it would simply mis-segment.
PANELS = {
    "S2A": dict(
        png="p14_figS2.png", native_dpi=355.7, axis_x=185,
        y_range=(90, 640), x_right=1891,
        tick_values=[10.0 ** e for e in (10, 9, 8, 7, 6, 5, 4)],
        quantity="Total flux (p/s)", units="p/s",
        fill={"Brain": (235, 155, 159), "Heart": (241, 205, 165),
              "Lungs": (235, 238, 186), "Liver": (178, 204, 216),
              "Kidneys": (182, 179, 231), "Spleen": (222, 179, 227)},
        locus="SI Figure S2A (p. 14)"),
    "S1A": dict(
        png="p13_figS1.png", native_dpi=220.2, axis_x=110,
        y_range=(25, 372), x_right=1160,
        tick_values=[100.0, 10.0, 1.0, 0.1, 0.01],
        quantity="% of total luminescence", units="percent", min_bar_w=5,
        fill={"Brain": (249, 150, 157), "Heart": (248, 203, 159),
              "Lungs": (234, 238, 180), "Liver": (171, 205, 218),
              "Kidneys": (183, 179, 235), "Spleen": (230, 177, 230)},
        locus="SI Figure S1A (p. 13)"),
}
PNG = os.path.join(HERE, "work", "p14_figS2.png")
NATIVE_DPI = 355.7           # measured: see check P7
AXIS_X = 185                 # y-axis spine, measured
ORGANS = ["Brain", "Heart", "Lungs", "Liver", "Kidneys", "Spleen"]
LNPS = ["13", "14", "70", "96", "97", "98", "45", "102"]
# Legend LUT, sampled from this same raster. Organ identity is confirmed a
# SECOND way (left-to-right order inside each group); the two must agree.
FILL = {
    "Brain":   (235, 155, 159),
    "Heart":   (241, 205, 165),
    "Lungs":   (235, 238, 186),
    "Liver":   (178, 204, 216),
    "Kidneys": (182, 179, 231),
    "Spleen":  (222, 179, 227),
}
PANEL_A = (90, 640)          # y range of panel A's plot box
# The plot box's right edge is TAKEN FROM AXIS FURNITURE, not guessed: the drawn
# baseline rule runs to x = 1889 and the first legend swatch starts at x = 1917.
# A guessed 1760 silently clipped the 8th group for four of six organs and the
# script reported "found 7 windows, expected 8" -- a refusal, which is why it
# was visible at all.
PLOT_X = (AXIS_X + 2, 1891)


def load():
    a = np.asarray(Image.open(PNG).convert("RGB")).astype(int)
    return a


def detect_ticks(dark, lo, hi, ax=AXIS_X):
    """Tick centroids from marks protruding LEFT of the spine. Detected, never
    supplied -- axis.calibrate has no API to hand it an origin, by design."""
    seg = dark[lo:hi, ax - 18:ax - 2]
    rows = np.where(seg.sum(1) >= 10)[0] + lo
    if rows.size == 0:
        raise DetectionError("no y-axis tick marks found")
    grp, cur = [], [rows[0]]
    for r in rows[1:]:
        if r - cur[-1] <= 2:
            cur.append(r)
        else:
            grp.append(cur)
            cur = [r]          # a NEW list. cur.clear() here aliases the list
            # already appended to grp, collapsing every group to the last one --
            # and the collapsed ladder yields LOO residual 0.0, i.e. a PERFECT
            # looking calibration from a detector that saw nothing. Found by
            # this script's own tick_gaps_px print, which read all zeros.
    grp.append(cur)
    return [float(np.mean(g)) for g in grp]


def colour_mask(a, rgb, tol=26):
    d = np.abs(a - np.array(rgb)).sum(2)
    return d <= tol


def close_interior(mask, x0, x1):
    """Fill interior holes in a bar's colour mask, row by row, inside [x0,x1).

    WHY, measured on this panel. GraphPad draws the error-bar STEM as a black
    vertical rule down the CENTRE of the bar, over the fill. The colour mask
    therefore reads '####....######' -- 10 of 14 px -- for every row of the
    bar's upper section, which is below bars.FRAC, so the full-width rule
    skipped the true top and fired lower down. Five bars refused this way.

    A bar is a SOLID rectangle. An interior gap bounded by bar colour on both
    sides in the same row is an overlay, not an edge, so it is filled. Anything
    NOT bounded on both sides is left alone, which is what keeps this safe:

      - an error-bar CAP above the bar top has no bar colour either side of it,
        so nothing is filled and the cap can never become a top (error class 32);
      - a stem ABOVE the bar top likewise has no flanking colour;
      - replicate dots sitting on the top corners are not interior and are
        untouched.

    The fix is applied to the MASK, which is what bars.bar_top's own docstring
    instructs ('Fix the mask, or record the bar as unreadable'). bars.py is not
    modified. Check P9 confirms every bar that already read cleanly reads the
    SAME row afterwards.
    """
    out = mask.copy()
    sub = mask[:, x0:x1]
    for y in range(sub.shape[0]):
        idx = np.where(sub[y])[0]
        if idx.size >= 2:
            out[y, x0 + idx[0]:x0 + idx[-1] + 1] = True
    return out


def bar_windows(mask, y0, y1, x0, x1, min_w=8):
    """Contiguous x runs where this colour appears inside the plot box."""
    col = mask[y0:y1, x0:x1].sum(0)
    on = col > 3
    runs, s = [], None
    for i, v in enumerate(on):
        if v and s is None:
            s = i
        elif not v and s is not None:
            if i - s >= min_w:
                runs.append((s + x0, i + x0))
            s = None
    if s is not None and len(on) - s >= min_w:
        runs.append((s + x0, len(on) + x0))
    return runs


def whisker_top(dark, x0, x1, y_lo, y_hi):
    """Topmost row of the error-bar CAP inside the bar's x window.

    The cap is a short horizontal black rule NARROWER than the bar. It is found
    by requiring a dark run that spans >= 35% of the window, scanning down from
    the plot top. Returns None rather than a number when nothing qualifies.
    """
    w = x1 - x0
    for y in range(y_lo, y_hi):
        run = dark[y, x0:x1]
        if run.sum() >= 0.35 * w:
            return y
    return None


def read_panel(key, cfg):
    """Digitise one panel. Returns (qc dict, list of bar records)."""
    a = np.asarray(Image.open(os.path.join(HERE, "work", cfg["png"])
                              ).convert("RGB")).astype(int)
    dark = (a.sum(2) < 250)
    y0, y1 = cfg["y_range"]
    ax = cfg["axis_x"]
    dpi = cfg["native_dpi"]

    ticks = detect_ticks(dark, 0, y1 + 10, ax=ax)
    tol_loo = axis.native_tol(axis.LOO_TOL_PX, native_dpi=dpi)
    tol_res = axis.native_tol(axis.TICK_RESIDUAL_PX, native_dpi=dpi)
    loo = axis.loo_residuals(ticks)
    qc = {
        "locus": cfg["locus"], "quantity": cfg["quantity"],
        "native_dpi_measured": dpi,
        "native_dpi_method": "image px width / image bbox pt width * 72",
        "LOO_TOL_PX_at_600dpi": axis.LOO_TOL_PX,
        "LOO_TOL_PX_rescaled_to_native": round(tol_loo, 4),
        "TICK_RESIDUAL_PX_at_600dpi": axis.TICK_RESIDUAL_PX,
        "TICK_RESIDUAL_PX_rescaled_to_native": round(tol_res, 4),
        "n_ticks_detected": len(ticks),
        "n_ticks_expected": len(cfg["tick_values"]),
        "tick_px": [round(t, 2) for t in ticks],
        "tick_gaps_px": [round(b - x, 2) for x, b in zip(ticks, ticks[1:])],
        "loo_worst_px": round(float(np.max(np.abs(loo))), 4),
        "loo_passes_at_STRICT_600dpi_tol": bool(np.max(np.abs(loo)) <= axis.LOO_TOL_PX),
        "loo_passes_at_native_tol": bool(np.max(np.abs(loo)) <= tol_loo),
    }
    if len(ticks) != len(cfg["tick_values"]):
        raise DetectionError(f"{key}: detected {len(ticks)} ticks, "
                             f"expected {len(cfg['tick_values'])}")
    cal = axis.calibrate(ticks, cfg["tick_values"], kind="log10", axis="y",
                         residual_tol_px=tol_res, loo_tol_px=tol_loo,
                         forbid_target=None)
    qc["calibration_px_per_decade"] = round(abs(float(cal.px_per_unit)), 4)

    base = np.where(dark[y0:y1, ax + 40:cfg["x_right"] - 40].sum(1) > 600)[0] + y0
    if base.size:
        bc = float(np.mean([base.min(), base.max()]))
        qc["drawn_baseline_rows"] = [int(base.min()), int(base.max())]
        qc["lowest_tick_px"] = round(max(ticks), 2)
        qc["baseline_centre_px"] = round(bc, 2)
        qc["lowest_tick_minus_drawn_baseline_px"] = round(max(ticks) - bc, 3)

    rows, failures, p9 = [], {}, []
    xs_by_organ = {}
    for organ in ORGANS:
        mask = colour_mask(a, cfg["fill"][organ])
        wins = bar_windows(mask, y0, y1, ax + 2, cfg["x_right"],
                           min_w=cfg.get("min_bar_w", 8))
        if len(wins) != len(LNPS):
            failures[organ] = f"found {len(wins)} windows, expected {len(LNPS)}"
            continue
        xs_by_organ[organ] = [w[0] for w in wins]
        closed = mask
        for w in wins:
            closed = close_interior(closed, w[0], w[1])
        res, fails = bars.bar_tops(closed, wins, y0, y1)
        res_raw, _ = bars.bar_tops(mask, wins, y0, y1)
        for i, w in enumerate(wins):
            if i in fails:
                failures[f"{organ}/{LNPS[i]}"] = fails[i]
                continue
            top = res[i]["row"]
            if i in res_raw:
                p9.append([f"{organ}/{LNPS[i]}", int(res_raw[i]["row"]), int(top)])
            val = float(cal.to_value(top))
            rows.append(dict(
                panel=key, lnp=LNPS[i], organ=organ,
                x0=int(w[0]), x1=int(w[1]), top_px=float(top),
                value=val, log10_value=float(np.log10(val)),
                value_raw_mask=(float(cal.to_value(res_raw[i]["row"]))
                                if i in res_raw else None),
                coverage=res[i]["coverage"],
                flank_delta_px=res[i]["flank_delta_px"],
                width_px=res[i]["width_px"]))
    qc["n_bars_read"] = len(rows)
    qc["n_refusals"] = len(failures)
    qc["refusals"] = failures

    moved = [t for t in p9 if t[1] != t[2]]
    qc["P9_effect_of_mask_closing"] = dict(
        n_readable_from_raw_mask=len(p9), n_moved=len(moved),
        n_moved_upward=sum(1 for t in moved if t[2] < t[1]),
        n_moved_downward=sum(1 for t in moved if t[2] > t[1]),
        max_move_px=max([abs(t[2] - t[1]) for t in moved], default=0),
        moved=moved,
        what=("The raw colour mask is occluded by the error-bar STEM, which "
              "GraphPad draws down the centre of the bar. EVERY move must be "
              "UPWARD: an occlusion can only make a detector under-read, so a "
              "single downward move would mean the closing is filling "
              "something that is not a bar. The uniform sign is the evidence "
              "that this is a fix and not a nudge."))

    if not failures:
        ok = all(all(xs_by_organ[ORGANS[k]][g] < xs_by_organ[ORGANS[k + 1]][g]
                     for k in range(len(ORGANS) - 1))
                 for g in range(len(LNPS)))
        qc["organ_order_matches_colour_LUT"] = bool(ok)
    return qc, rows


def main():
    out = {"panels": {}, "bars": []}
    for key, cfg in PANELS.items():
        qc, rows = read_panel(key, cfg)
        out["panels"][key] = qc
        out["bars"].extend(rows)
        print(f"{key}: {qc['n_bars_read']} bars, {qc['n_refusals']} refusals, "
              f"LOO {qc['loo_worst_px']} px (strict pass="
              f"{qc['loo_passes_at_STRICT_600dpi_tol']}), "
              f"closing moved {qc['P9_effect_of_mask_closing']['n_moved']} "
              f"(up {qc['P9_effect_of_mask_closing']['n_moved_upward']}, "
              f"down {qc['P9_effect_of_mask_closing']['n_moved_downward']})")

    # ---- P10/P11: the two panels are the SAME measurements, two parameterisations
    v = {}
    for b in out["bars"]:
        v.setdefault(b["panel"], {})[(b["lnp"], b["organ"])] = b["value"]
    cross = {"_what": (
        "SI Fig S1A states '% of total signal = (organ signal)/(total organ "
        "signal)*100'. So S1A is a CLOSED-FORM function of S2A. Neither panel "
        "was used to calibrate the other; the two were digitised from separate "
        "rasters at different dpi with separately sampled colour LUTs.")}
    if "S1A" in v and "S2A" in v:
        per_lnp = {}
        missing = [k for l in LNPS for o in ORGANS
                   if (k := (l, o)) not in v.get("S1A", {})]
        cross["S1A_bars_refused"] = ["/".join(k) for k in missing]
        for l in LNPS:
            tot = sum(v["S2A"][(l, o)] for o in ORGANS)
            pred = {o: 100.0 * v["S2A"][(l, o)] / tot for o in ORGANS}
            obs = {o: v["S1A"][(l, o)] for o in ORGANS if (l, o) in v["S1A"]}
            if len(obs) < len(ORGANS):
                per_lnp[l] = dict(INCOMPLETE=True,
                                  missing=[o for o in ORGANS if o not in obs],
                                  obs_over_pred={o: round(obs[o] / pred[o], 4)
                                                 for o in obs})
                continue
            ratios = [obs[o] / pred[o] for o in ORGANS]
            per_lnp[l] = dict(
                S1A_sum_over_organs=round(sum(obs.values()), 3),
                predicted=  {o: round(pred[o], 4) for o in ORGANS},
                observed=   {o: round(obs[o], 4) for o in ORGANS},
                obs_over_pred={o: round(obs[o] / pred[o], 4) for o in ORGANS},
                worst_abs_log10_dev=round(
                    float(max(abs(np.log10(r)) for r in ratios)), 4))
        cross["P11_S1A_vs_100xS2A_over_sum"] = per_lnp
        allr = [r for l in LNPS for r in per_lnp[l]["obs_over_pred"].values()]
        cross["P11_summary"] = dict(
            n=len(allr),
            median_obs_over_pred=round(float(np.median(allr)), 4),
            worst_obs_over_pred=round(float(max(allr, key=lambda x: abs(np.log10(x)))), 4),
            n_within_1pct=sum(1 for r in allr if abs(r - 1) <= 0.01),
            n_within_5pct=sum(1 for r in allr if abs(r - 1) <= 0.05),
            n_within_10pct=sum(1 for r in allr if abs(r - 1) <= 0.10),
            n_within_20pct=sum(1 for r in allr if abs(r - 1) <= 0.20))
        sums = [per_lnp[l]["S1A_sum_over_organs"] for l in LNPS
                if "S1A_sum_over_organs" in per_lnp[l]]
        cross["P10_S1A_percentages_sum_to_100"] = dict(
            per_lnp={l: per_lnp[l]["S1A_sum_over_organs"] for l in LNPS
                     if "S1A_sum_over_organs" in per_lnp[l]},
            n_lnps_complete=len(sums),
            worst_abs_dev_from_100=round(float(max(abs(s - 100) for s in sums)), 3),
            what=("The authors never print this. It is implied by their own "
                  "definition of the quantity, and it tests the S1A read "
                  "WITHOUT reference to S2A."))
        print("P10 worst |sum-100| =",
              cross["P10_S1A_percentages_sum_to_100"]["worst_abs_dev_from_100"])
        print("P11", cross["P11_summary"])
    # ---- P12: does the mask closing IMPROVE agreement, or merely change it?
    # A fix must move the AGREEMENT COUNT, not just the numbers. Recomputing
    # P11 on the RAW (un-closed) S2A read is the seed: if closing were wrong,
    # the raw read would agree with S1A at least as well.
    vraw = {}
    for b in out["bars"]:
        if b["panel"] == "S2A" and b.get("value_raw_mask") is not None:
            vraw[(b["lnp"], b["organ"])] = b["value_raw_mask"]
    # only LNPs whose SIX organs all read from the raw mask can be compared,
    # because the percentage needs a complete organ total.
    LNP_OK = [l for l in LNPS if all((l, o) in vraw for o in ORGANS)]
    if "S1A" in v and LNP_OK:
        def agree(src):
            rs = []
            for l in LNP_OK:
                tot = sum(src[(l, o)] for o in ORGANS)
                for o in ORGANS:
                    if (l, o) in v["S1A"]:
                        rs.append(v["S1A"][(l, o)] / (100.0 * src[(l, o)] / tot))
            return rs
        a_closed, a_raw = agree(v["S2A"]), agree(vraw)
        cross["P12_closing_improves_agreement_with_S1A"] = dict(
            n=len(a_closed), lnps_compared=LNP_OK,
            lnps_excluded=[l for l in LNPS if l not in LNP_OK],
            why_excluded=("the raw, un-closed mask REFUSED at least one of "
                          "that LNP's six organs, so no organ total can be "
                          "formed for it from the raw read"),
            closed=dict(
                median_abs_log10_dev=round(float(np.median([abs(np.log10(r)) for r in a_closed])), 5),
                n_within_10pct=sum(1 for r in a_closed if abs(r - 1) <= 0.10),
                n_within_20pct=sum(1 for r in a_closed if abs(r - 1) <= 0.20)),
            raw_unclosed=dict(
                median_abs_log10_dev=round(float(np.median([abs(np.log10(r)) for r in a_raw])), 5),
                n_within_10pct=sum(1 for r in a_raw if abs(r - 1) <= 0.10),
                n_within_20pct=sum(1 for r in a_raw if abs(r - 1) <= 0.20)),
            what=("THE AGREEMENT COUNT MUST MOVE. This is the seed for the "
                  "mask-closing repair: the un-closed read is the defect, and "
                  "S1A -- a different raster, different dpi, different colour "
                  "LUT, different calibration -- is the external judge. If the "
                  "two columns were equal, the closing would be unproven."))
        print("P12 closed:", cross["P12_closing_improves_agreement_with_S1A"]["closed"])
        print("P12 raw   :", cross["P12_closing_improves_agreement_with_S1A"]["raw_unclosed"])
    out["cross_panel"] = cross
    json.dump(out, open(os.path.join(HERE, "04_figS1S2_measurements.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
