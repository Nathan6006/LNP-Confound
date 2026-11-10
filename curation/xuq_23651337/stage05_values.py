#!/usr/bin/env python3
"""Stage 05 -- every value located in the only file we have.

THE FILE WE HAVE is the 3-page SI. It carries three figures:

  Figure S1  TEM micrographs of 16II/DNA and 16II/mRNA   -- no numbers at all
  Figure S2  EB exclusion assay vs heparin concentration -- numbers, not delivery
  Figure S3  Cell Viability (%) of HeLa vs lipidoid 1-16 -- numbers, not delivery

So the count of corpus-admissible values is ZERO, and the one panel that
carries per-lipidoid numbers is digitised anyway -- because "we looked" has to
be a measurement, not an assertion, and because the panel is what pins the
library size at 16.

TIER 4b. Nothing is rendered and no pixel is read. The panel is native vector:
page.get_drawings() returns 16 filled rectangles (the bars), 16 stroked
rectangles at stroke width 0.4596 pt (their outlines), error bars at 0.2757 pt
and axis furniture at 1.1949 pt. BAR BODY, OUTLINE, ERROR BAR AND AXIS ARE
SEPARATED BY EXACT STROKE WIDTH, which closes error class 32 by construction.
Because it is vector there is no native dpi and no pixel tolerance applies;
axis.native_tol() is deliberately NOT called (it refuses a guessed dpi, and a
guess is the only thing we could give it).

UNITS. Every coordinate is converted from PDF points to the 600-dpi pixels
that curation/lib's tolerances are written for (x 600/72 = x 8.3333). That
conversion is exact and lossless on vector coordinates, and it is what makes
axis.LOO_TOL_PX mean what it says. Tolerances this file sets itself are in
POINTS and are stated with their 600-dpi-pixel equivalent, because a vector
read has no detector noise: its only slack is PDF coordinate rounding, so a
1.5 px tolerance would be ~75x too loose to catch anything.

Run:  python3 -I stage05_values.py            (clean)
      python3 -I stage05_values.py --seed X   (seeded; must exit non-zero)
"""
import argparse
import csv
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))

import pymupdf                                             # noqa: E402
import axis                                                # noqa: E402

PDF = os.path.join(HERE, "source", "sb300023h_si_001.pdf")
PAGE = 2                                      # 0-based; SI page S3
PT_TO_PX600 = 600.0 / 72.0

# stroke widths, measured -- used as the SELECTOR, which is the point
SW_OUTLINE = 0.4596
SW_ERRORBAR = 0.2757
SW_AXIS = 1.1949
SW_TOL = 1e-3

# y-axis labels, read off the text layer with their glyph-centre y
Y_LABELS = [("100", 119.495), ("80", 149.975), ("60", 180.455),
            ("40", 211.055), ("20", 241.415), ("0", 272.015)]

# Tolerances THIS FILE sets, in points. See the UNITS note above.
TOL_PT_COINCIDENT = 0.01      # two drawn objects meant to be at the same y
# (an exact-coincidence tolerance was tried and the file does not satisfy it;
# see CF. Removed rather than loosened.)
TOL_PT_LABEL_MARGIN = 5.0     # a label must beat the runner-up tick by this


def _d(x):
    return x * PT_TO_PX600


class Geometry:
    """Everything read out of the vector page, and nothing else."""

    def __init__(self, seed=None):
        self.seed = seed
        doc = pymupdf.open(PDF)
        dr = doc[PAGE].get_drawings()

        self.fills, self.outlines = [], []
        segs = {}
        for p in dr:
            w = p.get("width") or 0.0
            if p["type"] == "f" and p.get("fill") == (0.0, 0.0, 0.0) \
                    and all(i[0] == "re" for i in p["items"]):
                self.fills.append(p["items"][0][1])
            elif p["type"] == "s" and all(i[0] == "re" for i in p["items"]) \
                    and abs(w - SW_OUTLINE) < SW_TOL:
                self.outlines.append(p["items"][0][1])
            elif p["type"] == "s":
                for it in p["items"]:
                    if it[0] == "l":
                        segs.setdefault(round(w, 4), []).append((it[1], it[2]))
        self.segs = segs

        if len(self.fills) != 16 or len(self.outlines) != 16:
            raise RuntimeError(
                f"expected 16 bar fills and 16 outlines, found "
                f"{len(self.fills)} and {len(self.outlines)}. Refusing rather "
                "than digitising a panel whose object count is not what the "
                "figure shows.")

        self.fills.sort(key=lambda r: r.x0)
        self.outlines.sort(key=lambda r: r.x0)

        # --- axis furniture, at SW_AXIS -----------------------------------
        ax = self.segs.get(SW_AXIS, [])
        hz = [s for s in ax if abs(s[0].y - s[1].y) < 1e-6
              and abs(s[0].x - s[1].x) > 1e-6]
        vt = [s for s in ax if abs(s[0].x - s[1].x) < 1e-6
              and abs(s[0].y - s[1].y) > 1e-6]
        # y ticks point inward from the axis; majors are ~2.88 pt, minors ~1.44
        self.y_major = sorted({round(s[0].y, 4) for s in hz
                               if 2.5 < abs(s[1].x - s[0].x) < 3.3})
        self.y_minor = sorted({round(s[0].y, 4) for s in hz
                               if 1.0 < abs(s[1].x - s[0].x) < 2.0})
        # x ticks hang below the baseline; majors ~3.36 pt, minors ~1.68
        self.x_major = sorted({round(s[0].x, 4) for s in vt
                               if 3.0 < abs(s[1].y - s[0].y) < 3.8})
        self.x_minor = sorted({round(s[0].x, 4) for s in vt
                               if 1.3 < abs(s[1].y - s[0].y) < 2.0})
        # the plot frame: the two full-width horizontal rules
        self.frame_y = sorted({round(s[0].y, 4) for s in hz
                               if abs(s[1].x - s[0].x) > 150})

        # --- error bars, at SW_ERRORBAR -----------------------------------
        eb = self.segs.get(SW_ERRORBAR, [])
        self.eb_stems = [s for s in eb if abs(s[0].x - s[1].x) < 1e-6]
        self.eb_caps = [s for s in eb if abs(s[0].y - s[1].y) < 1e-6]

        self._apply_seed()

    # -- seeding ----------------------------------------------------------
    def _apply_seed(self):
        s = self.seed
        if s is None:
            return
        if s == "seed_major_tick":
            # +0.25 pt (2.08 px at 600 dpi) on one INTERIOR major y tick.
            # Deliberately SUBTLE: least squares has leverage ~0.19 on an
            # interior tick, so calibrate's max-residual test absorbs most of
            # it and PASSES, while the leave-one-out residual -- which has no
            # such masking -- does not. This is the seed that exercises CB
            # itself and moves its agreement count.
            self.y_major[2] += 0.25
        elif s == "seed_major_tick_gross":
            # +3 pt. calibrate REFUSES outright (residual 20.92 px vs its 2.5
            # px limit) before any agreement count can be computed, so this
            # seed is credited by the REFUSAL, not by a count moving. Kept
            # separate from seed_major_tick so the two detection routes are
            # never conflated.
            self.y_major[2] += 3.0
        elif s == "seed_minor_tick":
            self.y_minor[3] += 0.5
        elif s == "seed_bartop_fill":
            r = self.fills[5]
            self.fills[5] = pymupdf.Rect(r.x0, r.y0 + 1.0, r.x1, r.y1)
        elif s == "seed_errbar_origin":
            # move bar 8's stem origin away from its top
            tgt = round(self.outlines[7].y0, 4)
            for i, (a, b) in enumerate(self.eb_stems):
                if abs(round(b.y, 4) - tgt) < 1e-3:
                    self.eb_stems[i] = (a, pymupdf.Point(b.x, b.y + 1.0))
                elif abs(round(a.y, 4) - tgt) < 1e-3:
                    self.eb_stems[i] = (pymupdf.Point(a.x, a.y + 1.0), b)
        elif s == "seed_baseline":
            r = self.fills[11]
            self.fills[11] = pymupdf.Rect(r.x0, r.y0, r.x1, r.y1 + 1.0)
        elif s == "seed_xtick":
            # more than half the 12.48 pt pitch, so the ASSIGNMENT flips.
            # A 1 pt displacement would NOT flip it and would NOT move CF's
            # agreement count -- a seed that leaves the count where it was is
            # a no-op seed and credits nothing.
            self.x_major[9] += 7.0
        elif s == "seed_grid":
            r = self.outlines[6]
            self.outlines[6] = pymupdf.Rect(r.x0, r.y0, r.x1 + 0.5, r.y1)
        elif s == "seed_label":
            # swap two labels' glyph centres: the ASSIGNMENT becomes wrong
            global Y_LABELS
            Y_LABELS = list(Y_LABELS)
            Y_LABELS[1], Y_LABELS[2] = Y_LABELS[2], Y_LABELS[1]
        elif s == "seed_frame_top":
            self.frame_y[0] += 2.0
        else:
            raise SystemExit(f"unknown seed {s!r}")


def run(seed=None):
    g = Geometry(seed=seed)
    checks, fails = {}, []

    # ---------------- C-A  label -> tick assignment, ABSOLUTE -------------
    # Nearest-tick criterion with an explicit margin. NOT a spread about the
    # mean of its own set: that is shift-invariant and would score clean on a
    # uniformly displaced read (the defect a sibling's category check had).
    agree, margins, worst_lab = 0, [], None
    for i, (txt, cy) in enumerate(Y_LABELS):
        d = sorted((abs(cy - t), j) for j, t in enumerate(g.y_major))
        nearest, runner = d[0], d[1]
        if nearest[1] == i:
            agree += 1
        margin = runner[0] - nearest[0]
        margins.append(margin)
        if worst_lab is None or margin < worst_lab[1]:
            worst_lab = (txt, margin)
    checks["CA_y_label_to_tick"] = {
        "criterion": "each y-axis label's glyph-centre y must be nearest to "
                     "the tick it is assigned to, beating the runner-up tick "
                     f"by >= {TOL_PT_LABEL_MARGIN} pt",
        "agreement": f"{agree}/6", "n_agree": agree, "n_total": 6,
        "worst_margin_pt": round(worst_lab[1], 4),
        "worst_label": worst_lab[0],
        "tol_margin_pt": TOL_PT_LABEL_MARGIN,
        "mean_label_minus_tick_offset_pt": round(
            sum(cy - g.y_major[i] for i, (_, cy) in enumerate(Y_LABELS)) / 6, 4),
        "note": "the ~0.33 pt systematic label-minus-tick offset is font "
                "metrics (glyph bbox centre vs baseline), not a read error. "
                "It is reported, and it is NOT used for the scale: the scale "
                "comes from tick positions only.",
    }
    if agree != 6 or worst_lab[1] < TOL_PT_LABEL_MARGIN:
        fails.append(f"CA: {agree}/6 labels on their own tick, worst margin "
                     f"{worst_lab[1]:.4f} pt")

    # ---------------- C-B  y calibration, from TICKS ONLY ----------------
    if len(g.y_major) != 6:
        fails.append(f"CB: {len(g.y_major)} major y ticks, expected 6")
        cal = None
    else:
        ticks_px = [_d(t) for t in g.y_major]
        minor_px = [_d(t) for t in g.y_minor]
        cal = axis.calibrate(
            ticks_px, [100.0, 80.0, 60.0, 40.0, 20.0, 0.0],
            kind="linear", axis="y",
            # the minor ticks are the hold-out this calibration is checked
            # against; naming them makes error class 27 a type error
            forbid_target=minor_px)
        loo = axis.loo_residuals(ticks_px)
        checks["CB_y_calibration"] = dict(
            cal.qc(),
            units="600-dpi pixels (PDF points x 600/72); the conversion is "
                  "exact on vector coordinates and is what makes "
                  "axis.LOO_TOL_PX mean what it was written to mean",
            tick_positions_pt=[round(t, 4) for t in g.y_major],
            tick_values=[100, 80, 60, 40, 20, 0],
            loo_worst_px=round(max(loo), 4),
            loo_worst_pt=round(max(loo) / PT_TO_PX600, 4),
            loo_tol_px=axis.LOO_TOL_PX,
            n_agree=sum(1 for r in loo if r <= axis.LOO_TOL_PX),
            n_total=len(loo),
            agreement=f"{sum(1 for r in loo if r <= axis.LOO_TOL_PX)}/{len(loo)}",
            native_tol_called=False,
            native_tol_reason="the panel is VECTOR. There is no native dpi to "
                              "measure, and axis.native_tol() refuses a "
                              "guessed one. No tolerance was widened.",
        )
        if max(loo) > axis.LOO_TOL_PX:
            fails.append(f"CB: LOO worst {max(loo):.4f} px > "
                         f"{axis.LOO_TOL_PX} px")

    # ---------------- C-C  minor-tick HOLD-OUT ---------------------------
    # The authors drew six minor ticks and printed none of them. They are
    # external to this read in a way a recomputed marginal is not.
    if cal is not None:
        want = [110.0, 90.0, 70.0, 50.0, 30.0, 10.0]
        agree_m, errs = 0, []
        for v, obs in zip(want, g.y_minor):
            pred = cal.to_pixel(v) / PT_TO_PX600
            e = abs(pred - obs)
            errs.append(e)
            if e * PT_TO_PX600 <= axis.LOO_TOL_PX:
                agree_m += 1
        checks["CC_minor_tick_holdout"] = {
            "criterion": "the 6 minor y ticks -- HELD OUT of the fit -- must "
                         "fall at 110/90/70/50/30/10 %, within "
                         f"{axis.LOO_TOL_PX} px (600 dpi)",
            "agreement": f"{agree_m}/6", "n_agree": agree_m, "n_total": 6,
            "worst_err_pt": round(max(errs), 4),
            "worst_err_px600": round(max(errs) * PT_TO_PX600, 4),
            "tol_px600": axis.LOO_TOL_PX,
            "predicted_values_pct": want,
            "observed_minor_tick_y_pt": [round(t, 4) for t in g.y_minor],
            "why_this_is_an_unprinted_check": "no minor tick carries a "
                "printed numeral anywhere on the page; the authors drew them "
                "and left them unlabelled. Nothing in our pipeline could have "
                "fitted to them, because they were passed to calibrate() as "
                "forbid_target.",
        }
        if agree_m != 6:
            fails.append(f"CC: {agree_m}/6 minor ticks predicted, worst "
                         f"{max(errs) * PT_TO_PX600:.4f} px")

    # ---------------- C-D  bar top, THREE independent objects -------------
    # fill rect y0, stroked outline y0, and the error bar's stem origin. Three
    # separately drawn objects at three different stroke widths. Error class
    # 32 cannot arise: the replicate/whisker geometry is a DIFFERENT object,
    # selected out by its own stroke width, not thresholded away.
    stem_origin = {}
    for a, b in g.eb_stems:
        # the two stems of one error bar meet at the bar top; the shared
        # endpoint is whichever y both stems have
        for p in (a, b):
            stem_origin.setdefault(round(p.x, 2), []).append(round(p.y, 4))
    agree3, spreads, detail = 0, [], []
    for i in range(16):
        f_y, o_y = g.fills[i].y0, g.outlines[i].y0
        cx = round((g.outlines[i].x0 + g.outlines[i].x1) / 2, 2)
        cand = [v for k, v in stem_origin.items() if abs(k - cx) < 1.0]
        eb_y = None
        if cand:
            ys = cand[0]
            shared = [y for y in set(ys) if ys.count(y) >= 2]
            eb_y = shared[0] if shared else None
        vals = [f_y, o_y] + ([eb_y] if eb_y is not None else [])
        spread = max(vals) - min(vals)
        spreads.append(spread)
        ok = spread <= TOL_PT_COINCIDENT and eb_y is not None
        agree3 += bool(ok)
        detail.append({"lipidoid": i + 1, "fill_y0_pt": round(f_y, 4),
                       "outline_y0_pt": round(o_y, 4),
                       "errorbar_stem_origin_y_pt": eb_y,
                       "spread_pt": round(spread, 5)})
    checks["CD_bar_top_three_way"] = {
        "criterion": "for every bar, the FILL rectangle's top, the STROKED "
                     "OUTLINE's top and the ERROR BAR's stem origin must "
                     f"agree within {TOL_PT_COINCIDENT} pt "
                     f"({TOL_PT_COINCIDENT * PT_TO_PX600:.3f} px at 600 dpi)",
        "agreement": f"{agree3}/16", "n_agree": agree3, "n_total": 16,
        "worst_spread_pt": round(max(spreads), 5),
        "worst_spread_px600": round(max(spreads) * PT_TO_PX600, 4),
        "tol_pt": TOL_PT_COINCIDENT,
        "why_error_class_32_cannot_arise": "the error bar is a DIFFERENT "
            "drawing object at a DIFFERENT stroke width (0.2757 vs 0.4596 pt) "
            "and is selected out by that width. Nothing is thresholded, so "
            "there is no threshold to be wrong. F-018 and F-020 were raster "
            "bar-top detectors reading an overlaid marker; that failure mode "
            "has no analogue here.",
        "per_bar": detail,
    }
    if agree3 != 16:
        fails.append(f"CD: {agree3}/16 bars agree three ways, worst spread "
                     f"{max(spreads):.5f} pt")

    # ---------------- C-E  baseline is the 0 % tick ----------------------
    zero_tick = g.y_major[-1]
    agree_b, errs_b = 0, []
    for i in range(16):
        e = abs(g.fills[i].y1 - zero_tick)
        errs_b.append(e)
        if e <= TOL_PT_COINCIDENT:
            agree_b += 1
    checks["CE_baseline_is_zero_tick"] = {
        "criterion": "every bar's BOTTOM must sit on the 0 % major tick "
                     f"within {TOL_PT_COINCIDENT} pt",
        "agreement": f"{agree_b}/16", "n_agree": agree_b, "n_total": 16,
        "worst_err_pt": round(max(errs_b), 5),
        "zero_tick_y_pt": round(zero_tick, 4),
        "why_it_matters": "a bar chart whose bars do not start at the axis "
            "zero is not a bar chart of the quantity it looks like. This is "
            "the cheapest possible guard against a shifted origin, which is "
            "error class 20's signature.",
    }
    if agree_b != 16:
        fails.append(f"CE: {agree_b}/16 bars on the zero tick, worst "
                     f"{max(errs_b):.5f} pt")

    # ---------------- C-F  bar -> label assignment, ABSOLUTE --------------
    # NEAREST-TICK with a margin, not exact coincidence.
    #
    # WHY NOT COINCIDENCE. The first version of this check required the bar
    # centre to sit on its x tick within 0.02 pt and REFUSED at 0.0607 pt on
    # 8 of 16 bars. The threshold was not loosened; the source was found.
    # The panel's x geometry is QUANTISED: fit a uniform grid to the 16 bar
    # RIGHT edges and the max residual is 0.00111 pt, while the x ticks, the
    # bar LEFT edges and the minor ticks each alternate by +/-0.06 pt about
    # the same grid -- and 0.06 pt is half of 0.1191 pt, exactly the
    # difference between the two bar widths the file uses (6.2391 and 6.1200).
    # The generator snapped some series to a coarser device grid; the drawn
    # centre of a quantised rectangle is therefore NOT on its tick, and a
    # coincidence criterion was testing a claim the file does not make.
    # For a CATEGORICAL axis the claim that matters is the ASSIGNMENT, and the
    # right criterion for that is nearest tick with an absolute margin.
    # CH below measures the quantisation instead of hiding it.
    agree_x, margins_x, offs = 0, [], []
    if len(g.x_major) != 16:
        fails.append(f"CF: {len(g.x_major)} major x ticks, expected 16")
    else:
        for i in range(16):
            c = (g.outlines[i].x0 + g.outlines[i].x1) / 2
            d = sorted((abs(c - t), j) for j, t in enumerate(g.x_major))
            if d[0][1] == i:
                agree_x += 1
            margins_x.append(d[1][0] - d[0][0])
            offs.append(c - g.x_major[i])
        checks["CF_bar_to_label_assignment"] = {
            "criterion": "each bar's centre x must be NEAREST to the x major "
                         "tick carrying its own label, beating the runner-up "
                         f"tick by >= {TOL_PT_LABEL_MARGIN} pt",
            "agreement": f"{agree_x}/16", "n_agree": agree_x, "n_total": 16,
            "worst_margin_pt": round(min(margins_x), 4),
            "tol_margin_pt": TOL_PT_LABEL_MARGIN,
            "measured_centre_minus_tick_offset_pt": {
                "max_abs": round(max(abs(o) for o in offs), 5),
                "values": [round(o, 5) for o in offs],
                "explanation": "alternating +/-0.06 pt. This is DEVICE "
                    "QUANTISATION, measured in CH, not a read error: the "
                    "bar right edges lie on a uniform grid to 0.0011 pt "
                    "while ticks and left edges alternate about it, and the "
                    "alternation equals half the 0.1191 pt difference between "
                    "the file's two bar widths.",
            },
            "why_it_matters": "this is what makes the bar -> lipidoid-number "
                "assignment MEASURED rather than assumed from left-to-right "
                "order. A dropped or duplicated bar, or a tick displaced by "
                "more than half the 12.48 pt pitch, would break it.",
        }
        if agree_x != 16 or min(margins_x) < TOL_PT_LABEL_MARGIN:
            fails.append(f"CF: {agree_x}/16 bars on their own tick, worst "
                         f"margin {min(margins_x):.4f} pt")

    # ---------------- C-G  axis top is the 110 % minor tick --------------
    # A round number the authors DREW and never printed: the plot frame's top
    # rule lands on the unlabelled 110 % minor tick.
    if cal is not None and g.frame_y:
        top = g.frame_y[0]
        pred110 = cal.to_pixel(110.0) / PT_TO_PX600
        e = abs(top - pred110)
        checks["CG_frame_top_is_110pct"] = {
            "criterion": "the plot frame's top rule must land on the 110 % "
                         "level predicted by the major-tick fit",
            "frame_top_y_pt": round(top, 4),
            "predicted_110pct_y_pt": round(pred110, 4),
            "err_pt": round(e, 5), "err_px600": round(e * PT_TO_PX600, 4),
            "tol_px600": axis.LOO_TOL_PX,
            "agreement": f"{1 if e * PT_TO_PX600 <= axis.LOO_TOL_PX else 0}/1",
            "n_agree": 1 if e * PT_TO_PX600 <= axis.LOO_TOL_PX else 0,
            "n_total": 1,
            "why_this_is_external": "the axis upper limit is a CHOICE the "
                "authors made and never printed. 110 % is not a value any "
                "measurement would produce and nothing in our read could have "
                "been fitted to it.",
        }
        if e * PT_TO_PX600 > axis.LOO_TOL_PX:
            fails.append(f"CG: frame top off 110 % by "
                         f"{e * PT_TO_PX600:.4f} px")

    # ---------------- C-H  ONE x grid, four independent series ------------
    # An unprinted relation the authors drew: the panel's x geometry is a
    # single uniform grid, and the plot box is EXACTLY 16 bins wide.
    def grid(v):
        n = len(v)
        idx = [float(i) for i in range(n)]
        mi = sum(idx) / n
        mv = sum(v) / n
        num = sum((idx[i] - mi) * (v[i] - mv) for i in range(n))
        den = sum((idx[i] - mi) ** 2 for i in range(n))
        p = num / den
        c = mv - p * mi
        r = max(abs(v[i] - (p * i + c)) for i in range(n))
        return p, c, r

    x1s = [g.outlines[i].x1 for i in range(16)]
    x0s = [g.outlines[i].x0 for i in range(16)]
    p1, c1, r1 = grid(x1s)                      # the EXACT series
    p0, _, r0 = grid(x0s)
    pm, _, rm = grid(g.x_minor) if len(g.x_minor) > 2 else (None, None, None)
    pM, _, rM = grid(g.x_major) if len(g.x_major) > 2 else (None, None, None)
    TOL_GRID_EXACT = 0.01       # the right-edge series must be this exact
    TOL_PITCH_AGREE = 0.01      # independent series must agree on the pitch
    box_w = (g.frame_y and (max(s[1].x for s in g.segs[SW_AXIS]
                                if abs(s[0].y - s[1].y) < 1e-6
                                and abs(s[1].x - s[0].x) > 150)
                            - min(s[0].x for s in g.segs[SW_AXIS]
                                 if abs(s[0].y - s[1].y) < 1e-6
                                 and abs(s[1].x - s[0].x) > 150))) or None
    sub = []
    sub.append(("right_edge_series_is_an_exact_grid", r1 <= TOL_GRID_EXACT))
    sub.append(("minor_tick_pitch_agrees",
                pm is not None and abs(pm - p1) <= TOL_PITCH_AGREE))
    sub.append(("left_edge_pitch_agrees", abs(p0 - p1) <= TOL_PITCH_AGREE))
    sub.append(("major_tick_pitch_agrees",
                pM is not None and abs(pM - p1) <= TOL_PITCH_AGREE))
    sub.append(("plot_box_is_exactly_16_bins_wide",
                box_w is not None and abs(box_w - 16 * p1) <= TOL_GRID_EXACT))
    n_ok = sum(1 for _, ok in sub if ok)
    checks["CH_one_x_grid"] = {
        "criterion": "four independently drawn series (bar right edges, bar "
                     "left edges, minor x ticks, major x ticks) must share "
                     f"one pitch within {TOL_PITCH_AGREE} pt, the right-edge "
                     f"series must be uniform to {TOL_GRID_EXACT} pt, and the "
                     "plot box must be exactly 16 bins wide",
        "agreement": f"{n_ok}/5", "n_agree": n_ok, "n_total": 5,
        "subchecks": {k: bool(v) for k, v in sub},
        "pitch_pt": {"bar_right_edges": round(p1, 6),
                     "bar_left_edges": round(p0, 6),
                     "minor_x_ticks": None if pm is None else round(pm, 6),
                     "major_x_ticks": None if pM is None else round(pM, 6)},
        "max_grid_residual_pt": {
            "bar_right_edges": round(r1, 5),
            "bar_left_edges": round(r0, 5),
            "minor_x_ticks": None if rm is None else round(rm, 5),
            "major_x_ticks": None if rM is None else round(rM, 5)},
        "plot_box_width_pt": None if box_w is None else round(box_w, 4),
        "sixteen_bins_pt": round(16 * p1, 4),
        "box_minus_16_bins_pt": None if box_w is None
        else round(box_w - 16 * p1, 5),
        "the_unprinted_number": "the plot box is 199.6795 pt wide and 16 bins "
            "at the measured 12.4800 pt pitch is 199.6799 pt -- agreeing to "
            "0.0005 pt. Nobody printed either figure; the authors chose an "
            "axis range of exactly 16 bins and the drawing records it. This "
            "is external to our read in the way a recomputed marginal is not.",
        "the_quantisation_this_measures": "bar right edges are uniform to "
            "0.0011 pt; bar left edges, minor ticks and major ticks each "
            "alternate about the same grid by up to 0.090 pt. The alternation "
            "is half the 0.1191 pt difference between the file's two bar "
            "widths (6.2391 / 6.1200), i.e. the generator snapped some series "
            "to a coarser device grid. This is why CF uses nearest-tick "
            "assignment and not coincidence.",
    }
    if n_ok != 5:
        fails.append(f"CH: {n_ok}/5 x-grid subchecks, residual "
                     f"{r1:.5f} pt on the right-edge series")

    # ---------------- the values ------------------------------------------
    rows = []
    if cal is not None:
        for i in range(16):
            top = g.outlines[i].y0
            v = cal.to_value(_d(top))
            # error bar, from the caps
            cx = (g.outlines[i].x0 + g.outlines[i].x1) / 2
            caps = sorted(round(a.y, 4) for a, b in g.eb_caps
                          if abs((a.x + b.x) / 2 - cx) < 2.0)
            hi = lo = None
            if len(caps) == 2:
                hi = cal.to_value(_d(caps[0]))
                lo = cal.to_value(_d(caps[1]))
            rows.append({
                "lipidoid_label": i + 1,
                "value": round(v, 2),
                "err_low": None if lo is None else round(lo, 2),
                "err_high": None if hi is None else round(hi, 2),
            })

    out = {
        "panel": "Figure S3 (its own caption; the SI's table of contents "
                 "calls it Figure S1 -- see 08_contradictions.csv)",
        "tier": "4b_pdf_vector_geometry",
        "readout_as_printed": "Cell Viability (%)",
        "n_values": len(rows),
        "checks": checks,
        "checks_failed": fails,
        "values": rows,
    }
    return out, fails


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", default=None)
    ap.add_argument("--quiet", action="store_true")
    ap.add_argument("--json", action="store_true",
                    help="print the checks dict as JSON, for the seeded-check "
                         "harness to compare AGREEMENT COUNTS against the "
                         "clean run. A seed that leaves every count where it "
                         "was is a no-op seed and credits nothing.")
    a = ap.parse_args()

    out, fails = run(seed=a.seed)

    if a.seed is None:
        # -------- 05_values.csv: EVERY located value, with admissibility ----
        cols = ["value_id", "panel", "panel_per_caption", "source_location",
                "item_label", "readout_as_printed", "units_as_printed",
                "value", "err_low", "err_high", "cargo", "cell_or_organ",
                "route", "dose", "timepoint", "n_replicates", "assay",
                "tier", "corpus_admissible", "admissible_reason"]
        common = dict(
            panel="S3", panel_per_caption="Figure S3",
            source_location="SI p. S3, bar chart, vector drawing objects",
            readout_as_printed="Cell Viability", units_as_printed="%",
            cargo="NOT STATED on this panel",
            cell_or_organ="HeLa (stated in the caption)",
            route="in vitro", dose="NOT STATED anywhere in the SI",
            timepoint="NOT STATED anywhere in the SI",
            n_replicates="NOT STATED; two error-bar caps are drawn per bar, "
                         "which implies replication but states no n",
            assay="cell viability (the SI prints no assay name -- the table of "
                  "contents says 'Cytotoxicity assay', the caption says 'Cell "
                  "variability' [sic]; neither names MTT, MTS or similar)",
            tier="4b_pdf_vector_geometry",
            corpus_admissible="no",
            admissible_reason="VIABILITY / CYTOTOXICITY. data/'s "
                "Experiment_method vocabulary is delivery-only (luminescence "
                "/ protein-abundance / editing-efficiency / knockdown). "
                "Viability is explicitly not a delivery endpoint, so this row "
                "cannot be emitted no matter how exactly it is measured.")
        rows = []
        for r in out["values"]:
            rows.append(dict(common, value_id=f"S3_lip{r['lipidoid_label']}",
                             item_label=str(r["lipidoid_label"]),
                             value=r["value"], err_low=r["err_low"],
                             err_high=r["err_high"]))
        # the two panels that supply no per-item number, recorded so the
        # accounting is complete rather than silently short
        rows.append({**{c: "" for c in cols},
                     "value_id": "S2_curve", "panel": "S2",
                     "panel_per_caption": "Figure S2",
                     "source_location": "SI p. S2, two-series line plot",
                     "item_label": "16I and 16II",
                     "readout_as_printed": "Relative Fluorescence Intensity",
                     "units_as_printed": "%",
                     "value": "NOT DIGITISED",
                     "cargo": "plasmid DNA", "cell_or_organ": "none (cuvette)",
                     "route": "n/a", "dose": "heparin 0-16 ug/mL (the x axis)",
                     "timepoint": "NOT STATED", "n_replicates": "NOT STATED",
                     "assay": "ethidium bromide exclusion / heparin displacement",
                     "tier": "4b_pdf_vector_geometry",
                     "corpus_admissible": "no",
                     "admissible_reason": "NUCLEIC-ACID BINDING, not delivery. "
                        "A cuvette polyanion-competition assay: no cells, no "
                        "cargo expression, no knockdown. Its x axis (heparin "
                        "concentration) has no corpus column. Not digitised "
                        "because no admissible row could come of it."},
                    )
        rows.append({**{c: "" for c in cols},
                     "value_id": "S1_TEM", "panel": "S1",
                     "panel_per_caption": "Figure S1",
                     "source_location": "SI p. S2, 4 embedded Indexed-colour "
                                        "rasters, 1681 x ~150 px each",
                     "item_label": "16II/DNA and 16II/mRNA",
                     "readout_as_printed": "none -- micrographs",
                     "value": "NO NUMBER EXISTS",
                     "assay": "transmission electron microscopy",
                     "tier": "5_pixels",
                     "corpus_admissible": "no",
                     "admissible_reason": "PARTICLE MORPHOLOGY. Micrographs "
                        "carry no numeric value; no scale-bar value and no "
                        "stated diameter appears in the page text layer. Size "
                        "is explicitly not a delivery endpoint."})
        with open(os.path.join(HERE, "05_values.csv"), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols)
            w.writeheader()
            for r in rows:
                w.writerow({c: r.get(c, "") for c in cols})

        out["accounting"] = {
            "rows_in_05_values_csv": len(rows),
            "per_item_numeric_values_measured": len(out["values"]),
            "panels_recorded_without_a_per_item_number": 2,
            "sum_check": f"{len(out['values'])} + 2 = {len(rows)}",
            "corpus_admissible_yes": 0,
            "corpus_admissible_no": len(rows),
            "admissible_plus_inadmissible_equals_total":
                0 + len(rows) == len(rows),
        }
        with open(os.path.join(HERE, "05_values_qc.json"), "w") as fh:
            json.dump(out, fh, indent=2)

    if a.json:
        print(json.dumps({k: {"n_agree": c.get("n_agree"),
                              "n_total": c.get("n_total")}
                          for k, c in out["checks"].items()}))
        return 1 if fails else 0

    if not a.quiet:
        for k, c in out["checks"].items():
            print(f"{k}: {c.get('agreement', c.get('skipped', 'n/a'))}"
                  + (f"  worst={c['worst_err_pt']} pt"
                     if "worst_err_pt" in c else
                     (f"  worst_spread={c['worst_spread_pt']} pt"
                      if "worst_spread_pt" in c else
                      (f"  loo={c['loo_worst_px']} px"
                       if "loo_worst_px" in c else ""))))
        if out["values"]:
            print("values (Cell Viability %, NOT corpus-admissible):")
            print("  " + "  ".join(f"{r['lipidoid_label']}:{r['value']:.1f}"
                                   for r in out["values"]))
    if fails:
        print("FAILED:", *fails, sep="\n  - ")
        return 1
    print("all checks pass" + (f" (seed={a.seed})" if a.seed else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())
