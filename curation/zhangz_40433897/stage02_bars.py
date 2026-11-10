#!/usr/bin/env python3
"""Stage 02 - bar digitisation of every SI panel that carries a DELIVERY value.

Tier 5 (pixels). 00_value_sources.json records why tiers 1-4b are unavailable:
no Source Data, no SI numeric table of values, no author repo, no text numerals,
and 0 vector paths on every figure page (one flattened DCTDecode JPEG per
figure).

NATIVE DPI. These are embedded rasters at 173.8-220.5 dpi, not a 600 dpi page
render. Every px tolerance is therefore rescaled with axis.native_tol(), and
BOTH the 600-dpi-equivalent figure and the native-px figure are recorded per
panel. 07_qc.json seeds the widened tolerance to show it still bites.

GEOMETRY. Bars are hollow: a coloured outline with a white interior. The top of
such a bar is a full-width horizontal coloured rule, which is exactly what
bars.bar_top looks for; the replicate dots, which sit on and above the top rule
(including on both top corners), are not full-width and are rejected by the
flank test. Bar windows are taken from the two wall clusters each bar leaves in
a band just above the baseline, never from an assumed pitch.

COLOUR IS A CHECK, NOT AN IDENTIFIER. Bars are identified by their x-order
inside a group, which the legend fixes. The measured wall colour is then matched
against legend swatches sampled FROM THE SAME PANEL and the assignment must come
out as the legend order, or this script refuses. Colours are never matched
across panels (the trap that gave a sibling a 4.583 swatch distance, under
heatmap.MAX_DIST, on two lipids that had swapped colours).
"""
import json
import os
import sys

import numpy as np
from PIL import Image

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "lib"))
import axis            # noqa: E402
import bars            # noqa: E402
from errors import DetectionError, CalibrationError  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
SI = os.path.join(HERE, "source", "nn5c01378_si_001.pdf")

INK_THRESHOLD = 60      # 255 - min(R,G,B) > this  => not paper white


# --------------------------------------------------------------------------
# Panel descriptions. Every number here was measured, not assumed; the
# measurement that produced it is named.
# --------------------------------------------------------------------------
PANELS = [
    dict(
        pid="FigS5a", page=7, xref=79, native_dpi=219.5,
        figure="Figure S5a", readout="Luciferase expression (%)",
        cell="4T1-Luc", cargo="siLuc",
        spine_x=86, baseline_row=259,
        tick_values=[120.0, 80.0, 40.0, 0.0],
        plot=dict(x0=89, x1=640, y_top=55),
        max_bar_width=18,
        groups=["PEG-0.5", "PEG-1", "PEG-2"],
        series=["VN3C2O", "VN2C4", "VN2C2", "VN3C2", "VN3C1", "VN2C3"],
        legend_x=(514, 552), legend_y=(70, 230),
    ),
    dict(
        pid="FigS5b", page=7, xref=79, native_dpi=219.5,
        figure="Figure S5b", readout="Relative luciferase expression (%)",
        cell="MC38-Luc", cargo="siLuc",
        spine_x=697, baseline_row=254,
        tick_values=[50.0, 40.0, 30.0, 20.0, 10.0, 0.0],
        plot=dict(x0=700, x1=1267, y_top=58),
        max_bar_width=22,
        groups=["PEG-0.5", "PEG-1", "PEG-2"],
        series=["VN3C2O", "VN2C4", "VN2C2", "VN3C2", "VN3C1", "VN2C3"],
        legend_x=(1124, 1162), legend_y=(60, 215),
    ),
    dict(
        pid="FigS4b", page=6, xref=77, native_dpi=219.9,
        figure="Figure S4b", readout="Luciferase expression (%)",
        cell="4T1-Luc", cargo="siLuc",
        spine_x=87, baseline_row=505,
        tick_values=[150.0, 100.0, 50.0, 0.0],
        plot=dict(x0=95, x1=1130, y_top=292),
        max_bar_width=16,
        groups=[""],
        series=["LNP%d" % i for i in range(1, 61)] + ["Lip3000", "Naked siRNA"],
        legend_x=None, legend_y=None,
    ),
    dict(
        pid="FigS7", page=7, xref=81, native_dpi=219.5,
        figure="Figure S7", readout="Percentage of GFP+ cells (%)",
        cell=None, cargo="pGFP",
        spine_x=77, baseline_row=307,
        tick_values=[60.0, 40.0, 20.0, 0.0],
        plot=dict(x0=80, x1=436, y_top=58),
        max_bar_width=38,
        groups=[""],
        series=["5:1", "10:1", "15:1", "20:1", "30:1", "Lip3000", "Naked pGFP"],
        legend_x=None, legend_y=None,
    ),
    dict(
        pid="FigS8a", page=8, xref=84, native_dpi=216.2,
        figure="Figure S8a", readout="Luciferase intensity (RLU)",
        cell=None, cargo="mFLuc",
        spine_x=121, baseline_row=299,
        tick_values=[800000.0, 600000.0, 400000.0, 200000.0, 0.0],
        plot=dict(x0=124, x1=380, y_top=35),
        max_bar_width=42,
        groups=[""],
        series=["PBS", "VE", "TPGS", "VN2C1", "Lip3000"],
        legend_x=None, legend_y=None,
    ),
    dict(
        pid="FigS8b", page=8, xref=84, native_dpi=216.2,
        figure="Figure S8b", readout="Luciferase intensity (RLU)",
        cell=None, cargo="mFLuc",
        spine_x=529, baseline_row=286,
        tick_values=[800000.0, 600000.0, 400000.0, 200000.0, 0.0],
        plot=dict(x0=532, x1=790, y_top=22),
        max_bar_width=42,
        groups=[""],
        series=["PBS", "VE+SM102", "TPGS+SM102", "SM102", "VN2C1"],
        legend_x=None, legend_y=None,
    ),
    dict(
        pid="FigS13b", page=9, xref=89, native_dpi=189.95,
        figure="Figure S13b", readout="GFP+ cells (%)",
        cell="RAW264.7", cargo="mGFP",
        spine_x=489, baseline_row=274,
        tick_values=[60.0, 40.0, 20.0, 0.0],
        plot=dict(x0=492, x1=648, y_top=65),
        max_bar_width=38,
        groups=[""],
        series=["PBS", "VN2C1 (day 0)", "VN2C1 (day 7)"],
        legend_x=None, legend_y=None,
    ),
]


def load_raster(page, xref):
    import pymupdf
    doc = pymupdf.open(SI)
    info = doc.extract_image(xref)
    pg = doc[page - 1]
    rect = pg.get_image_rects(xref)[0]
    dpi_x = info["width"] / rect.width * 72.0
    dpi_y = info["height"] / rect.height * 72.0
    import io
    im = Image.open(io.BytesIO(info["image"])).convert("RGB")
    return np.asarray(im).astype(int), dpi_x, dpi_y


def ink_mask(rgb):
    return (255 - rgb.min(axis=2)) > INK_THRESHOLD


def detect_y_ticks(gray, spine_x, lo=0, hi=None):
    """Tick marks are short horizontal rules OUTSIDE the spine, to its left.

    Measured on the spine's own left side so the tick centroid cannot be a
    glyph centroid (error class 20): the band is 6 px wide ending 2 px before
    the spine, which is narrower than any digit in this figure.
    """
    hi = gray.shape[0] if hi is None else hi
    band = gray[lo:hi, spine_x - 8:spine_x - 2] < 150
    rows = band.sum(axis=1)
    on = [i + lo for i in range(len(rows)) if rows[i] >= 4]
    if not on:
        raise DetectionError("detect_y_ticks: no tick rows found left of the "
                             "spine at x=%d." % spine_x)
    groups, cur = [], [on[0]]
    for v in on[1:]:
        if v - cur[-1] <= 1:
            cur.append(v)
        else:
            groups.append(cur)
            cur = [v]
    groups.append(cur)
    return [float(np.mean(g)) for g in groups], [(g[0], g[-1]) for g in groups]


def bar_windows(ink, x0, x1, baseline_row, max_width):
    """Bar windows from the wall clusters each bar leaves above the baseline.

    A hollow bar meets the baseline with two vertical walls. Consecutive wall
    clusters are merged while the merged span still fits inside max_width, so a
    bar whose wall has fused with an overlaid replicate dot still yields one
    window. The merge never crosses max_width, which is why it cannot silently
    swallow a neighbouring bar.
    """
    band = ink[baseline_row - 6:baseline_row - 1, x0:x1]
    cols = band.sum(axis=0)
    on = [i + x0 for i in range(len(cols)) if cols[i] >= 4]
    if not on:
        raise DetectionError("bar_windows: no ink above the baseline.")
    clusters, cur = [], [on[0]]
    for v in on[1:]:
        if v - cur[-1] <= 2:
            cur.append(v)
        else:
            clusters.append((cur[0], cur[-1]))
            cur = [v]
    clusters.append((cur[0], cur[-1]))

    windows, i = [], 0
    while i < len(clusters):
        start, end = clusters[i]
        j = i + 1
        while j < len(clusters) and clusters[j][1] - start + 1 <= max_width:
            end = clusters[j][1]
            j += 1
        windows.append((start, end + 1))
        i = j
    return windows, clusters


def trim_halo(intensity, baseline_row, windows, floor=0.5):
    """Drop window edge columns that are only a JPEG halo of the wall.

    WHY. A boolean ink mask makes a column that is 23% covered by the pen
    indistinguishable from one that is 100% covered, so a window taken from it
    can be one column too wide. That column then caps the full-width coverage
    at 13/14 = 0.929 and the 5 px flank at 4/5 = 0.80, and bars.bar_top
    correctly refuses - on 21 of FigS4b's 62 bars.

    The fix belongs in the WINDOW, not in the threshold, and not in a trim of
    the bar's interior (an inner-trimmed window is what let F-018 through).
    A column is part of the drawn wall if, over a band where the wall is solid,
    its mean pen intensity is at least `floor` of the strongest wall column's -
    i.e. the pixel is at least half covered by the pen. Only columns at the two
    OUTER edges are ever removed; nothing inside the window is touched, and the
    trim stops at the first column that qualifies.

    `intensity` is 255 - min(R,G,B), i.e. distance from paper white, NOT a
    boolean mask - the whole point is that the partial coverage is visible.
    """
    band = intensity[baseline_row - 9:baseline_row - 2]
    out = []
    for (x0, x1) in windows:
        prof = band[:, x0:x1].mean(axis=0)
        peak = prof.max()
        lo, hi = 0, len(prof) - 1
        while lo < hi and prof[lo] < floor * peak:
            lo += 1
        while hi > lo and prof[hi] < floor * peak:
            hi -= 1
        out.append((x0 + lo, x0 + hi + 1))
    return out


def strip_long_rules(ink, max_bar_width):
    """Remove horizontal runs longer than any bar can be.

    WHY. Figure S4b carries an INSET TABLE inside the plot box (the four
    printed Pearson r). Its three horizontal rules are ~300 px long and cross
    the search window of every bar beneath them. On LNP43 the detector read the
    table's top rule at row 300 and returned 149.38 on a 0-150 axis - a value
    no bar in the panel comes near.

    A bar's top rule is at most one bar wide. A continuous horizontal run three
    times wider than the widest bar is therefore furniture - a table rule, an
    axis, a significance bracket - and never a bar top. Only the pixels of such
    runs are cleared; nothing else is touched, so a bar whose top happens to sit
    on the same ROW as a table rule keeps its own ink.
    """
    limit = 3 * max_bar_width
    out = ink.copy()
    n_cleared = 0
    for y in range(ink.shape[0]):
        row = ink[y]
        x = 0
        while x < row.size:
            if not row[x]:
                x += 1
                continue
            j = x
            while j < row.size and row[j]:
                j += 1
            if j - x > limit:
                out[y, x:j] = False
                n_cleared += j - x
            x = j
    return out, n_cleared


def rule_thickness(ink, windows, y_top, y_bot, frac=bars.FRAC):
    """Measured thickness, in rows, of each bar's drawn top rule.

    Independent of the `nrow` that bars.bar_top will be given: it scans with a
    one-row requirement, then counts how many contiguous rows actually hold
    `frac` coverage. That is a property of the PEN the figure was drawn with,
    at the resolution the raster actually has, and it is the only honest basis
    for choosing `nrow`. A panel whose rule is 2 rows thick cannot satisfy a
    3-row contiguity requirement no matter how correct the digitisation is.
    """
    out = []
    for (x0, x1) in windows:
        sub = ink[y_top:y_bot, x0:x1]
        cov = sub.mean(axis=1)
        hit = np.where(cov >= frac)[0]
        if hit.size == 0:
            out.append(0)
            continue
        i = int(hit[0])
        n = 0
        while i + n < len(cov) and cov[i + n] >= frac:
            n += 1
        out.append(n)
    return out


def sample_legend(rgb, lx, ly):
    """Median RGB of each legend swatch ring, sampled in THIS panel only."""
    sub = rgb[ly[0]:ly[1], lx[0]:lx[1]]
    mask = (255 - sub.min(axis=2)) > INK_THRESHOLD
    rows = mask.sum(axis=1)
    on = [i for i in range(len(rows)) if rows[i] > 0]
    groups, cur = [], [on[0]]
    for v in on[1:]:
        if v - cur[-1] <= 3:
            cur.append(v)
        else:
            groups.append(cur)
            cur = [v]
    groups.append(cur)
    out = []
    for g in groups:
        m = mask[g[0]:g[-1] + 1]
        px = sub[g[0]:g[-1] + 1][m]
        # the ring is antialiased toward white; take the third of pixels
        # furthest from white, which is the pen colour itself
        d = 255 * 3 - px.sum(axis=1)
        keep = px[d >= np.percentile(d, 67)]
        out.append(np.median(keep, axis=0))
    return [np.asarray(c, dtype=float) for c in out]


def bar_colour(rgb, ink, x0, x1, top_row, baseline_row):
    """Median pen colour of a bar's walls, over its own height."""
    y0 = max(top_row + 2, 0)
    y1 = baseline_row - 1
    if y1 - y0 < 2:
        y0, y1 = top_row, baseline_row
    sub = rgb[y0:y1, x0:x1]
    m = ink[y0:y1, x0:x1]
    if m.sum() < 5:
        return None
    px = sub[m]
    d = 255 * 3 - px.sum(axis=1)
    keep = px[d >= np.percentile(d, 67)]
    return np.median(keep, axis=0)


def main():
    out = {"panels": {}}
    rows = []

    for P in PANELS:
        rgb, dpi_x, dpi_y = load_raster(P["page"], P["xref"])
        gray = rgb.mean(axis=2)
        ink = ink_mask(rgb)
        native = (dpi_x + dpi_y) / 2.0

        rec = {"measured_dpi_x": round(dpi_x, 1),
               "measured_dpi_y": round(dpi_y, 1),
               "declared_native_dpi": P["native_dpi"]}

        # ---- 1. ticks -> calibration ---------------------------------
        ticks, spans = detect_y_ticks(gray, P["spine_x"],
                                      lo=0, hi=P["baseline_row"] + 4)
        rec["tick_rows"] = [round(t, 2) for t in ticks]
        rec["tick_spans"] = spans
        if len(ticks) != len(P["tick_values"]):
            raise DetectionError(
                "%s: detected %d ticks, panel declares %d values %s."
                % (P["pid"], len(ticks), len(P["tick_values"]),
                   P["tick_values"]))

        # tolerances: the library's numbers are for a 600 dpi render.
        cal_tol_600 = axis.TICK_RESIDUAL_PX
        cal_tol_native = axis.native_tol(cal_tol_600, native_dpi=native)
        loo_tol_600 = 1.5
        loo_tol_native = axis.native_tol(loo_tol_600, native_dpi=native)
        rec["tolerances"] = {
            "calibrate_residual_px_at_600dpi": cal_tol_600,
            "calibrate_residual_px_native": round(cal_tol_native, 3),
            "loo_px_at_600dpi": loo_tol_600,
            "loo_px_native": round(loo_tol_native, 3),
            "scale_factor": round(600.0 / native, 3),
        }

        cal = axis.calibrate(ticks, P["tick_values"], kind="linear",
                             residual_tol_px=cal_tol_native, axis="y")
        loo = axis.loo_residuals(ticks)
        rec["calibration"] = {
            "px_per_unit": round(cal.px_per_unit, 5),
            "residual_px": round(cal.residual_px, 4),
            "loo_residuals_px": [round(v, 4) for v in loo],
            "worst_loo_px_native": round(max(loo), 4),
            "worst_loo_px_600dpi_equivalent": round(max(loo) * 600.0 / native, 4),
            "loo_passes_native_tol": bool(max(loo) <= loo_tol_native),
            "loo_passes_unrescaled_tol": bool(max(loo) <= loo_tol_600),
        }

        # the zero tick must BE the baseline - a round number the authors drew
        zero_idx = int(np.argmin(np.abs(np.asarray(P["tick_values"]))))
        rec["zero_tick_vs_baseline_px"] = round(
            float(ticks[zero_idx] - P["baseline_row"]), 3)

        # ---- 2. bar windows ------------------------------------------
        wins_raw, clusters = bar_windows(ink, P["plot"]["x0"], P["plot"]["x1"],
                                         P["baseline_row"], P["max_bar_width"])
        ink, n_cleared = strip_long_rules(ink, P["max_bar_width"])
        rec["furniture_pixels_cleared"] = int(n_cleared)
        intensity = 255 - rgb.min(axis=2)
        wins = trim_halo(intensity, P["baseline_row"], wins_raw)
        rec["halo_columns_trimmed"] = sum(
            (b[0] - a[0]) + (a[1] - b[1]) for a, b in zip(wins_raw, wins))
        rec["n_wall_clusters"] = len(clusters)
        rec["n_bar_windows"] = len(wins)
        expect = len(P["groups"]) * len(P["series"])
        rec["n_bars_expected"] = expect
        if len(wins) != expect:
            rec["windows"] = wins
            out["panels"][P["pid"]] = rec
            raise DetectionError(
                "%s: segmented %d bar windows, expected %d (%d groups x %d "
                "series). Windows: %s"
                % (P["pid"], len(wins), expect, len(P["groups"]),
                   len(P["series"]), wins))
        rec["bar_widths_px"] = [w[1] - w[0] for w in wins]

        # ---- 3. legend colours, this panel only ----------------------
        if P["legend_x"] is None:
            legend = []
            rec["legend_swatches_rgb"] = None
            rec["colour_check_note"] = (
                "This panel has no colour legend; bars are identified by their "
                "x-order against the printed category labels, so there is no "
                "colour assignment to check.")
        else:
            legend = sample_legend(rgb, P["legend_x"], P["legend_y"])
            rec["legend_swatches_rgb"] = [[int(v) for v in c] for c in legend]
            if len(legend) != len(P["series"]):
                raise DetectionError(
                    "%s: %d legend swatches sampled, %d series declared."
                    % (P["pid"], len(legend), len(P["series"])))

        # ---- 4. bar tops ---------------------------------------------
        thick = rule_thickness(ink, wins, P["plot"]["y_top"],
                               P["baseline_row"] - 1)
        nz = [t for t in thick if t > 0]
        med = int(np.median(nz)) if nz else 0
        rec["rule_thickness_rows"] = {
            "per_bar": thick, "median_nonzero": med,
            "min_nonzero": min(nz) if nz else None,
            "max": max(thick) if thick else None}
        nrow = int(P.get("nrow", bars.NROW))
        rec["nrow_used"] = nrow
        rec["nrow_default"] = bars.NROW
        results, failures = bars.bar_tops(
            ink, wins, P["plot"]["y_top"], P["baseline_row"] - 1,
            nrow=nrow, baseline_row=P["baseline_row"])

        # Per-bar retry at the thickness THAT BAR's rule actually has.
        #
        # Across these panels the drawn rule is 2-10 rows thick (per-bar figures
        # recorded above). A bar whose rule is 2 rows cannot satisfy a 3-row
        # contiguity requirement however correct the read is, and bars.bar_top
        # refuses - correctly, on its own terms. The retry lowers `nrow` only
        # to that bar's MEASURED thickness, never below 2, and only for bars
        # that failed the coverage clause; a bar that failed the FLANK clause is
        # never retried, because the flank test is the one that separates a top
        # rule from a replicate marker (error class 32) and relaxing it is
        # exactly the move F-018 made. Every retried bar is marked
        # `nrow_relaxed` and carries confidence = medium downstream.
        retried = {}
        for i in list(failures):
            if "no row in" not in failures[i]:
                continue              # a flank failure is never retried
            t = thick[i]
            if t < 2:
                continue
            try:
                r = bars.bar_top(ink, wins[i][0], wins[i][1],
                                 P["plot"]["y_top"], P["baseline_row"] - 1,
                                 nrow=t, baseline_row=P["baseline_row"])
            except DetectionError:
                continue
            results[i] = r
            retried[i] = t
            del failures[i]
        rec["nrow_relaxed_bars"] = {str(k): v for k, v in retried.items()}
        rec["n_refused"] = len(failures)
        rec["refusals"] = {str(k): v.splitlines()[0] for k, v in failures.items()}

        vals, colour_check = [], []
        for i, (x0, x1) in enumerate(wins):
            gi, si = divmod(i, len(P["series"]))
            name = P["series"][si]
            if not legend:
                if i not in results:
                    vals.append(dict(index=i, group=P["groups"][gi],
                                     series=name, x0=x0, x1=x1, value=None,
                                     status="REFUSED_BY_bars.bar_top"))
                    continue
                r = results[i]
                vals.append(dict(index=i, group=P["groups"][gi], series=name,
                                 x0=x0, x1=x1, top_row=r["row"],
                                 width_px=r["width_px"],
                                 coverage=r["coverage"],
                                 flank_delta_px=r["flank_delta_px"],
                                 value=round(float(cal.to_value(r["row"])), 3),
                                 nrow_relaxed=retried.get(i),
                                 status="ok"))
                continue
            if i not in results:
                vals.append(dict(index=i, group=P["groups"][gi], series=name,
                                 x0=x0, x1=x1, value=None,
                                 status="REFUSED_BY_bars.bar_top"))
                continue
            r = results[i]
            v = cal.to_value(r["row"])
            c = bar_colour(rgb, ink, x0, x1, r["row"], P["baseline_row"])
            dists = [float(np.linalg.norm(c - L)) for L in legend] if c is not None else []
            nearest = int(np.argmin(dists)) if dists else -1
            colour_check.append(dict(index=i, expected=si, nearest=nearest,
                                     dist=round(min(dists), 2) if dists else None,
                                     runner_up=round(sorted(dists)[1], 2)
                                     if len(dists) > 1 else None))
            vals.append(dict(index=i, group=P["groups"][gi], series=name,
                             x0=x0, x1=x1,
                             top_row=r["row"], width_px=r["width_px"],
                             coverage=r["coverage"],
                             flank_delta_px=r["flank_delta_px"],
                             value=round(float(v), 3),
                             nrow_relaxed=retried.get(i),
                             status="ok"))
        rec["colour_check"] = colour_check
        n_colour_ok = sum(1 for c in colour_check if c["expected"] == c["nearest"])
        rec["colour_order_agreement"] = "%d/%d" % (n_colour_ok, len(colour_check))
        rec["values"] = vals
        out["panels"][P["pid"]] = rec

        for v in vals:
            rows.append(dict(panel=P["pid"], figure=P["figure"],
                             readout=P["readout"], cell=P["cell"],
                             cargo=P["cargo"], **v))

    with open(os.path.join(HERE, "02_bars.json"), "w") as fh:
        json.dump(out, fh, indent=1)

    print(json.dumps({p: {k: v for k, v in d.items()
                          if k not in ("values", "colour_check")}
                      for p, d in out["panels"].items()}, indent=1))
    for p, d in out["panels"].items():
        print("\n== %s ==" % p)
        for v in d["values"]:
            print("  %-8s %-8s %s" % (v["group"], v["series"],
                                      v.get("value", v["status"])))
    return 0


if __name__ == "__main__":
    sys.exit(main())
