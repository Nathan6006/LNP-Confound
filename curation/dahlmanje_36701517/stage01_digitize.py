"""Digitise the tdTomato+ bar panels of Radmand 2023 (PMID 36701517).

TIER 5 (pixels). 00_value_sources.json states why tiers 1-4 do not exist for
this paper: the publisher deposited exactly one PDF, it has no table object,
and no repository or accession is named anywhere.

GEOMETRY RULES OBEYED (curation/lib):
  * every y scale is fitted by axis.calibrate from DETECTED tick positions
    only; there is no API here to supply an origin, and none is supplied.
  * every bar top goes through bars.bar_top, which refuses on a flank
    disagreement rather than returning the better-looking flank.
  * the panels' native resolution is MEASURED (220 dpi, see 07_qc.json), and
    axis.native_tol rescales the library's 600-dpi tolerances. Both numbers and
    a seeded displacement are recorded; nothing is widened silently.

Usage:  python3 -I stage01_digitize.py <out.json>
"""
import json
import os
import sys

import numpy as np
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import axis                                              # noqa: E402
import bars                                              # noqa: E402
from errors import DetectionError                        # noqa: E402

SI = os.path.join(HERE, "source", "nl2c04479_si_001.pdf")

# --- measured, not assumed. page -> (xref, native dpi) ----------------------
PANELS_PAGE = {"S13ab": (23, 244, 220.15), "S13c": (24, 249, 220.10),
               "S15": (26, 262, 219.80)}

NATIVE_DPI = 220.0
LOO_TOL = axis.native_tol(axis.LOO_TOL_PX, native_dpi=NATIVE_DPI)   # 4.09 px
CAL_TOL = axis.native_tol(axis.TICK_RESIDUAL_PX, native_dpi=NATIVE_DPI)

# In TOP-DOWN order, because detect_ticks_1d returns ascending ROW indices and
# axis.calibrate pairs the two lists element-wise after a common sort. Passing
# [0,20,...,100] against ascending rows silently fits an inverted axis: it
# calibrated cleanly (residual 0.150 px) with 0% at the TOP of the panel.
TICK_VALUES = [100.0, 80.0, 60.0, 40.0, 20.0, 0.0]

# Cell types per organ, IN DRAWN LEFT-TO-RIGHT ORDER, read off the rendered
# panel. The order DIFFERS between Fig S13 and Fig S15 (S13 lung is
# ECs/Monocytes/DCs/T/B, Fig 3b and S15 lung is ECs/DCs/T/B/Monocytes), which
# is exactly why every cross-check in stage03 joins on the LABEL, never on the
# column index.
S13_ORGANS = [
    ("lung",   ["ECs", "Monocytes", "DCs", "T cells", "B cells"]),
    ("liver",  ["ECs", "Kupffer Cells", "DCs", "T cells", "B cells",
                "Hepatocytes"]),
    ("kidney", ["ECs", "Monocytes", "DCs", "T cells", "B cells"]),
    ("spleen", ["DCs", "B cells", "T cells"]),
    ("heart",  ["ECs", "Immune cells"]),
]
S15_ORGANS = [
    ("lung",  ["ECs", "DCs", "T cells", "B cells", "Monocytes"]),
    ("liver", ["ECs", "DCs", "T cells", "B cells", "Kupffer Cells",
               "Hepatocytes"]),
    ("heart", ["ECs", "Immune cells"]),
]


def load(page, xref):
    import pymupdf
    doc = pymupdf.open(SI)
    info = doc.extract_image(xref)
    rect = doc[page].get_image_rects(xref)[0]
    dpi = info["width"] / (rect.width / 72.0)
    import io
    img = Image.open(io.BytesIO(info["image"])).convert("RGB")
    return np.asarray(img).astype(int), dpi


def masks(a):
    """Colour masks. Built from the bars' OWN colours, never from darkness."""
    r, g, b = a[..., 0], a[..., 1], a[..., 2]
    mx, mn = a.max(axis=2), a.min(axis=2)
    nearwhite = mn > 246
    chroma = mx - mn
    # purple / blue family: blue channel leads, green trails.
    coloured = (~nearwhite) & (chroma >= 10) & (b > g + 6) & (b >= r - 8)
    # grey family: achromatic but not the near-black of axis rules and text.
    grey = (~nearwhite) & (chroma < 10) & (mn >= 45)
    return coloured, grey


def find_spines(a):
    """Achromatic near-black vertical rules >= 120 px long -> y-axis spines."""
    mx, mn = a.max(axis=2), a.min(axis=2)
    # 140, not 110: MEASURED, not chosen for convenience. The same figure
    # draws its y-axis rules at two different ink levels -- the Fig S13a Lung
    # spine is grey 83 and the Fig S13a Spleen spine is grey 123-129 -- so a
    # 110 threshold finds three of five sub-panels and silently loses two.
    # Coloured bar outlines are excluded by the achromatic test (the Cat-LNP
    # outline is (113,90,138), chroma 48), not by this level.
    black = (mx < 140) & ((mx - mn) < 28)
    H, W = black.shape
    # Short runs first: a spine can JOG one pixel sideways partway down. The
    # Fig S13c Liver spine is column 596 for rows 100-265 and column 597 for
    # rows 265-343, so a per-column length test reports a 165 px rule and the
    # tick search then clips the ladder to four stubs.
    runs = []
    for x in range(W):
        col = black[:, x]
        s = None
        for y in range(H + 1):
            v = bool(col[y]) if y < H else False
            if v and s is None:
                s = y
            elif not v and s is not None:
                if y - s >= 20:
                    runs.append([x, s, y])
                s = None
    runs.sort(key=lambda t: (t[1], t[0]))
    groups = []
    for x, y0, y1 in runs:
        hit = None
        for g in groups:
            if (abs(x - g["x"]) <= 2 and y0 <= g["y1"] + 3
                    and y1 >= g["y0"] - 3):
                hit = g
                break
        if hit is None:
            groups.append({"x": x, "y0": y0, "y1": y1})
        else:
            hit["x"] = min(hit["x"], x)
            hit["y0"] = min(hit["y0"], y0)
            hit["y1"] = max(hit["y1"], y1)
    return [(g["x"], g["y0"], g["y1"]) for g in groups
            if g["y1"] - g["y0"] >= 120]


def panel_ticks(a, spine_x, y0, y1):
    gray = a.mean(axis=2)
    # The tick stubs are 1-3 px tall and sit in the 5 px gap between the
    # labels' right edge and the spine. MEASURED on the Fig S13a Lung panel:
    # the six stubs occupy rows 99-100, 146-148, 194-195, 241-243, 289-290,
    # 336-338 with minima 77, 59, 123, 56, 95, 60. The library default
    # (dark=140, min_len=2) finds ONE of the six, because five of them have
    # only a single row below 140. min_len=1 is therefore required here, and
    # the darkness weighting is what keeps a 2-row stub's centre sub-pixel.
    lo, hi = max(0, spine_x - 5), spine_x - 1
    centres = axis.detect_ticks_1d(gray, "y", lo, hi, dark=240,
                                   min_len=1, max_len=6)
    return [c for c in centres if y0 - 8 <= c <= y1 + 8]


def category_ticks(gray, baseline_row, lo=2, hi=4, dark=200):
    """x positions of the category stubs below the x-axis rule.

    A column counts only if it is dark in EVERY row of the window, not in the
    minimum over the window. MEASURED, twice: the rotated category labels
    start one or two rows below the stubs, and a min-over-rows detector picked
    up the top of "Monocytes" as a sixth stub in the Fig S15d Lung panel
    (389.1 beside the real 392.8) and the Liver panel's own "0" label in the
    Fig S13b Lung panel (576.7, 584.9). A stub spans the whole window; a glyph
    top does not.
    """
    win = gray[baseline_row + lo:baseline_row + hi, :]
    col = (win < dark).all(axis=0)
    weight = np.clip(dark - win.mean(axis=0), 0, None)
    out, s = [], None
    for i in range(len(col) + 1):
        v = bool(col[i]) if i < len(col) else False
        if v and s is None:
            s = i
        elif not v and s is not None:
            w = weight[s:i]
            if w.sum() > 0:
                out.append(float(np.dot(np.arange(s, i), w) / w.sum()))
            s = None
    return out


def bar_windows(mask, x_lo, x_hi, y_top, y_bot, min_w=5):
    """Contiguous column runs of `mask` inside the plot box."""
    sub = mask[y_top:y_bot, x_lo:x_hi]
    present = sub.any(axis=0)
    runs, s = [], None
    for i, v in enumerate(present):
        if v and s is None:
            s = i
        elif not v and s is not None:
            if i - s >= min_w:
                runs.append((x_lo + s, x_lo + i))
            s = None
    if s is not None and len(present) - s >= min_w:
        runs.append((x_lo + s, x_lo + len(present)))
    return runs


def band_runs(mask, x_lo, x_hi, baseline_row, lo=3, hi=5, min_w=5):
    """Column runs inked in EVERY row of [baseline-hi, baseline-lo)."""
    sub = mask[baseline_row - hi:baseline_row - lo + 1, x_lo:x_hi]
    if sub.size == 0:
        return []
    col = sub.all(axis=0)
    runs, s = [], None
    for i in range(len(col) + 1):
        v = bool(col[i]) if i < len(col) else False
        if v and s is None:
            s = i
        elif not v and s is not None:
            if i - s >= min_w:
                runs.append((x_lo + s, x_lo + i))
            s = None
    return runs


def read_panel(a, spine_x, y0, y1, x_right, labels, cal_tol,
               loo_tol, paired_offset):
    """Return {label: {series: value}} plus the calibration record."""
    ticks = panel_ticks(a, spine_x, y0, y1)
    # Five is allowed, six is the norm. MEASURED: the Fig S13b Spleen panel
    # has no 100% tick stub drawn at all (rows 1360-1374 of the raster carry
    # the "100" label glyphs and nothing between them and the spine), while
    # every other panel draws six. Five uniformly spaced stubs still define
    # the scale; what makes that safe is the spine check below, which pins
    # which VALUE the top stub carries instead of assuming the missing one is
    # the top one.
    if not 5 <= len(ticks) <= len(TICK_VALUES):
        raise DetectionError(
            f"panel at x={spine_x} y=[{y0},{y1}): detected {len(ticks)} y "
            f"ticks, expected 5 or 6. Refusing rather than calibrating on a "
            "partial ladder.")
    tick_values = TICK_VALUES[-len(ticks):]          # top-down, ending at 0
    # Gate on the UNRESCALED library tolerance. The rescaled one (loo_tol) is
    # carried only so the record says what the weaker fallback would have
    # been; it is never the operative gate unless a panel genuinely needs it,
    # and none does. Reporting a pass against 1.5 px when 4.09 px was
    # available is strictly the stronger statement.
    cal = axis.calibrate(ticks, tick_values, kind="linear",
                         residual_tol_px=axis.TICK_RESIDUAL_PX, axis="y")
    loo = axis.loo_residuals(ticks)
    gate = axis.LOO_TOL_PX
    if max(loo) > gate:
        if max(loo) > loo_tol:
            raise DetectionError(
                f"panel at x={spine_x}: worst leave-one-out tick residual "
                f"{max(loo):.2f} px > {loo_tol:.2f} px, the 1.5 px library "
                f"tolerance rescaled for a measured {NATIVE_DPI:.0f} dpi "
                "panel. Both gates fail.")
        raise DetectionError(
            f"panel at x={spine_x}: worst leave-one-out tick residual "
            f"{max(loo):.2f} px > the unrescaled {gate} px gate. It would "
            f"pass the {loo_tol:.2f} px rescaled gate, but widening is a "
            "decision to record in DECISIONS.md, not one to take here.")
    baseline = cal.to_pixel(0.0)

    # THE ROUND NUMBER THE AUTHORS DREW AND NEVER PRINTED. Every y-axis spine
    # in these figures is stroked from the 100% level to the 0% level. So the
    # row the fitted scale puts 100% at must coincide with the TOP of the
    # detected spine -- a prediction made from the tick stubs alone, tested
    # against a piece of furniture that took no part in the fit. It is also
    # what makes a five-stub ladder safe: it pins which value the top stub
    # carries rather than assuming the missing stub was the top one.
    spine_top_resid = abs(float(cal.to_pixel(100.0)) - float(y0))
    if spine_top_resid > 3.0:
        raise DetectionError(
            f"panel at x={spine_x}: the fitted scale puts 100% at row "
            f"{cal.to_pixel(100.0):.1f} but the y-axis spine starts at row "
            f"{y0} ({spine_top_resid:.1f} px apart). The spine is drawn from "
            "100% to 0%, so this means the tick ladder was mis-valued -- "
            "most likely a missing stub that is not the top one.")

    coloured, grey = masks(a)
    plot_top, plot_bot = int(min(ticks)) - 2, int(round(baseline)) - 1

    # Category slots come from the X-AXIS TICKS, not from the colour mask.
    # A slot whose bar is zero-height has no coloured pixels at all, so
    # deriving the slots from the mask loses the slot and then zips the
    # remaining labels onto the wrong bars -- measured here: the Fig S13a
    # Kidney panel yields 3 coloured runs for 5 cell types.
    gray = a.mean(axis=2)
    cats = [c for c in category_ticks(gray, int(round(baseline)))
            if spine_x < c < x_right]
    if len(cats) != len(labels):
        raise DetectionError(
            f"panel at x={spine_x}: {len(cats)} x-axis category ticks but "
            f"{len(labels)} cell-type labels. Refusing: silently zipping a "
            "short list onto a long one is how values land on the wrong "
            "category.")
    pitch = float(np.median(np.diff(cats))) if len(cats) > 1 else 79.0

    # Bar WIDTHS come from a two-row band just above the baseline, where a bar
    # of any appreciable height is solid, and a column counts only if it is
    # inked in BOTH rows. Not a trim toward the middle -- the window is the
    # bar's full width, measured where it is unambiguous. Taking the width
    # from "any row in the plot has this colour" instead picks up a column
    # that is inside the bar near the baseline and white at the top: measured
    # on the Fig S13c Liver ECs bar, that gave a 19 px window whose right
    # flank first fires 226 rows below its left flank, and bars.bar_top
    # correctly refused all six bars of the panel.
    col_runs = band_runs(coloured, spine_x + 1, x_right, int(round(baseline)))
    grey_runs = band_runs(grey, spine_x + 1, x_right, int(round(baseline)))

    def pick(runs, want_centre):
        best, bd = None, None
        for r in runs:
            c = 0.5 * (r[0] + r[1])
            d = abs(c - want_centre)
            if bd is None or d < bd:
                best, bd = r, d
        if best is None or bd > pitch * 0.30:
            return None
        return best

    vals, diag = {}, {}
    for lab, ct in zip(labels, cats):
        rec = {"category_tick_px": round(ct, 2), "pitch_px": round(pitch, 2)}
        series_spec = (
            (("lnp", coloured, col_runs, ct + paired_offset * pitch),
             ("pbs", grey, grey_runs, ct - paired_offset * pitch))
            if paired_offset else
            (("lnp", coloured, col_runs, ct),))
        for series, m, runs, want in series_spec:
            win = pick(runs, want)
            if win is None:
                rec[series] = {
                    "refused": "no bar at least 4 px tall centred within "
                               f"{pitch * 0.30:.0f} px of {want:.0f} -- the "
                               "bar is below the width-measurement floor",
                    "want_centre_px": round(want, 1)}
                continue
            wx0, wx1 = win
            try:
                t = bars.bar_top(m, wx0, wx1, plot_top, plot_bot + 1,
                                 baseline_row=int(round(baseline)))
                v = float(cal.to_value(t["row"]))
                rec[series] = {"row": t["row"], "value_pct": round(v, 3),
                               "flank_delta_px": t["flank_delta_px"],
                               "width_px": t["width_px"],
                               "coverage": t["coverage"],
                               "x0": wx0, "x1": wx1}
            except DetectionError as e:
                rec[series] = {"refused": str(e).split("\n")[0],
                               "x0": wx0, "x1": wx1}
        diag[lab] = rec
        vals[lab] = {k: rec.get(k, {}).get("value_pct")
                     for k in ("lnp", "pbs")}
    return vals, diag, {"spine_top_residual_px": round(spine_top_resid, 3),
                        "n_y_ticks": len(ticks),
                        "category_ticks_px": [round(c, 2) for c in cats],
                        "category_pitch_px": round(pitch, 2),
                        "ticks_px": [round(t, 2) for t in ticks],
                        "loo_max_px": round(max(loo), 3),
                        "loo_tol_px": round(loo_tol, 3),
                        "loo_tol_px_unrescaled": axis.LOO_TOL_PX,
                        "px_per_percent": round(cal.px_per_unit, 4),
                        "baseline_row": round(baseline, 2),
                        "calibration": cal.qc()}


def main(out_path):
    result = {"native_dpi_measured": {}, "panels": {}}
    a13, dpi13 = load(23, 244)
    a13c, dpi13c = load(24, 249)
    a15, dpi15 = load(26, 262)
    result["native_dpi_measured"] = {"S13ab": round(dpi13, 2),
                                     "S13c": round(dpi13c, 2),
                                     "S15": round(dpi15, 2)}
    for nm, d in result["native_dpi_measured"].items():
        if abs(d - NATIVE_DPI) > 1.5:
            raise SystemExit(f"{nm}: measured dpi {d} is not the {NATIVE_DPI} "
                             "this module's tolerances were rescaled for.")

    jobs = []
    # Fig S13a (Cat-LNP) and S13b (Neu-LNP) share one raster.
    # paired_offset: Fig S13 draws a PBS/LNP PAIR per category with the stub
    # between them, so each bar's centre is ~0.22 pitch either side of the
    # stub. Fig S15 draws ONE bar per category, centred on the stub. MEASURED:
    # Fig S15d Lung stubs at 203.3/250.7/298.1/345.5/392.8 against bar windows
    # (187,221) and (377,411), i.e. centres 204.0 and 394.0.
    jobs.append(("S13a", "Cat-LNP", a13, 0, 740, S13_ORGANS, 0.22))
    jobs.append(("S13b", "Neu-LNP", a13, 850, 1749, S13_ORGANS, 0.22))
    jobs.append(("S13c", "An-LNP", a13c, 0, 1500, S13_ORGANS, 0.22))
    # The y windows for the Fig S15 panels exclude the composition boxes of
    # S15a and S15c, whose "Diam. | 87" table rules are 193 px vertical
    # achromatic lines at x=489 and would be counted as a fourth axis spine.
    jobs.append(("S15b", "Neu-LNP", a15, 300, 700, S15_ORGANS, 0.0))
    jobs.append(("S15d", "An-LNP", a15, 1050, 1450, S15_ORGANS, 0.0))

    for panel, lnp, img, ylo, yhi, organs, paired_offset in jobs:
        spines = [s for s in find_spines(img) if ylo <= s[1] < yhi]
        spines.sort(key=lambda s: (round(s[1] / 60.0), s[0]))
        if len(spines) != len(organs):
            raise SystemExit(
                f"{panel}: found {len(spines)} axis spines in y=[{ylo},{yhi}) "
                f"but the panel has {len(organs)} organ sub-panels: {spines}")
        W = img.shape[1]
        out = {}
        for (organ, labels), (sx, sy0, sy1) in zip(organs, spines):
            others = [s[0] for s in spines
                      if abs(s[1] - sy0) < 40 and s[0] > sx]
            x_right = (min(others) - 30) if others else (W - 5)
            vals, diag, cal = read_panel(img, sx, sy0, sy1, x_right, labels,
                                         CAL_TOL, LOO_TOL, paired_offset)
            out[organ] = {"values_pct": vals, "diagnostics": diag,
                          "axis": cal, "spine_x": sx,
                          "spine_y": [sy0, sy1]}
        result["panels"][panel] = {"lnp": lnp, "organs": out}

    with open(out_path, "w") as fh:
        json.dump(result, fh, indent=1, sort_keys=True)
    print("wrote", out_path)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else
         os.path.join(HERE, "04_digitized.json"))
