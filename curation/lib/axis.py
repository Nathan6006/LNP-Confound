#!/usr/bin/env python3
"""ONE audited pixel -> units calibration, anchored on axis furniture only.

WHY THIS FILE EXISTS. Error class 20 -- a scale or an origin anchored on the
panel box, a border, or a glyph centroid instead of on detected tick marks --
has FOUR independent instances in four separately written files:

  curation/ql_2026/digitize_fig3c.py        10^3 on the box bottom, 47 px below
                                            the tick. Values up to 1.76x high.
  curation/ql_2026/digitize_fig2b.py        colorbar anchored on the OUTER EDGE
                                            of its own black border (F-019).
                                            +1% at the hit threshold, +38% at
                                            the floor.
  curation/lin_aromatic_2026/              scale fitted onto the hit rates the
    digitize_fig2fgh.py                     panel validated (F-006, class 27).
  the same file, after that repair          origin left as three hard-coded
                                            literals (F-033). Mean delta
                                            -0.52 pp; 47 of 48 deltas negative.

F-033 is the one that justifies a library rather than a convention: it landed in
a file in a folder where the class had already been found and fixed. A fix in
one folder is a fix in one folder.

THE CONTRACT

  cal = calibrate(ticks_px, tick_values, kind="log10"|"linear")
  value = cal.to_value(y_px)
  y_px  = cal.to_pixel(value)

1. Scale AND origin both come from `ticks_px`. There is no way to pass an origin
   separately, because that is exactly the shape F-033 took.
2. `calibrate` asserts the tick spacing is uniform ON THE POSITIONS IT WAS GIVEN.
   A log axis's decade gaps must be equal; an unequal one means the anchors are
   glyph centroids, not ticks.
3. It never returns a fallback. Too few ticks, non-monotone ticks, or residuals
   over tolerance raise CalibrationError.
4. `forbid_target=` makes class 27 a type error: pass the vector you intend to
   validate against and calibrate refuses if it is also the fit input.

F-017 is why rule 2 says "on the positions it was given": I once marked a
finding confirmed because a COMMENT listed anchors whose decade gaps were 406
and 408. They were the tick positions rounded to integers, and the code was
right while its comment was wrong. Verifying a claim about code by reading a
comment is the same error the comment was an instance of. So the assertion here
runs on floats the caller measured, and `Calibration.residual_px` is returned
for the caller to record -- not described in prose.
"""
import math

import numpy as np

from errors import CalibrationError, ValidationTargetError

# Default tolerance on the tick-fit residual, in pixels. 2.5 px at 600 dpi is
# ~0.1 pt -- tighter than any tick mark's own width, loose enough for the
# half-pixel quantisation of a darkness-weighted centroid. digitize_fig2fgh.py
# uses 2.5 for exactly this reason; it is the value that made F-006's bias
# visible instead of absorbing it.
TICK_RESIDUAL_PX = 2.5


class Calibration:
    """An immutable pixel <-> unit mapping, plus the evidence for it."""

    __slots__ = ("kind", "px_per_unit", "origin_px", "origin_value",
                 "residual_px", "n_ticks", "tick_step", "_sign")

    def __init__(self, kind, px_per_unit, origin_px, origin_value,
                 residual_px, n_ticks, tick_step, sign):
        self.kind = kind
        self.px_per_unit = px_per_unit
        self.origin_px = origin_px
        self.origin_value = origin_value
        self.residual_px = residual_px
        self.n_ticks = n_ticks
        self.tick_step = tick_step
        self._sign = sign

    # -- conversion ------------------------------------------------------
    def to_value(self, px):
        """Pixel coordinate -> data value. Vectorised over array input."""
        px = np.asarray(px, dtype=float)
        units = self._sign * (px - self.origin_px) / self.px_per_unit
        out = self.origin_value + units if self.kind == "linear" else \
            10.0 ** (math.log10(self.origin_value) + units)
        return float(out) if out.ndim == 0 else out

    def to_pixel(self, value):
        """Data value -> pixel coordinate. The inverse of to_value."""
        value = np.asarray(value, dtype=float)
        units = value - self.origin_value if self.kind == "linear" else \
            np.log10(value) - math.log10(self.origin_value)
        out = self.origin_px + self._sign * units * self.px_per_unit
        return float(out) if out.ndim == 0 else out

    # -- evidence --------------------------------------------------------
    def qc(self):
        """The dict to write into the folder's QC JSON. Record it, do not
        paraphrase it -- F-017 was a comment that disagreed with its code."""
        return {
            "kind": self.kind,
            "px_per_unit": round(self.px_per_unit, 4),
            "origin_px": round(self.origin_px, 3),
            "origin_value": self.origin_value,
            "tick_step_units": self.tick_step,
            "n_ticks": self.n_ticks,
            "tick_fit_residual_px": round(self.residual_px, 4),
            "tick_residual_tol_px": TICK_RESIDUAL_PX,
        }

    def __repr__(self):
        return (f"<Calibration {self.kind} {self.px_per_unit:.3f} px/unit "
                f"origin {self.origin_value}@{self.origin_px:.1f} "
                f"resid {self.residual_px:.3f}px n={self.n_ticks}>")


# Leave-one-out tolerance, in pixels. Tighter than TICK_RESIDUAL_PX because a
# LOO residual has no leverage masking: the tick being tested took no part in
# the fit that predicts it.
LOO_TOL_PX = 1.5          # written for a 600 dpi RENDER -- see native_tol()
RENDER_DPI = 600          # what every tolerance in this module assumes


def native_tol(tol_px, native_dpi, render_dpi=RENDER_DPI):
    """`tol_px` rescaled for a panel whose real resolution is `native_dpi`.

    WHY. Every pixel tolerance here was written for a page WE render at 600 dpi,
    where 1.5 px is a sub-pixel quantity. A flattened raster embedded in an SI
    has its own resolution, and rendering the page at 600 dpi does not create
    information that is not in the JPEG -- it upsamples. Measured 2026-11-10
    (`khanof_36652629`): six SI panels are embedded rasters at **188-221 dpi**,
    so this module's 1.5 px is **0.55 of a real pixel** and three panels failed
    a test no correct digitisation could pass.

    Rescaling is therefore not a loosening, it is the same physical tolerance
    expressed in the units the panel actually has. But A WIDENED TOLERANCE MUST
    BE SEEDED -- that curation seeded a **-3.25 px** displacement and confirmed
    it still fails LOO after rescaling, which is the only thing that makes the
    widening honest rather than convenient. Record BOTH numbers and the measured
    dpi beside the result; never rescale silently.

        tol = axis.native_tol(axis.LOO_TOL_PX, native_dpi=204)   # -> 4.41

    TWO THINGS MEASURED HERE THAT A CALLER MUST KNOW, 2026-11-10:

    1. **`calibrate`'s own hard 2.5 px residual limit is the real blocker, and
       it is NOT rescaled by this function.** Seeding a 4-tick log ladder at
       ~47 px/decade, a +5.0 px displacement on one tick makes `calibrate`
       RAISE at 3.50 px before `loo_residuals` is ever reached. So a low-dpi
       panel can fail at the calibrate step for the same units reason; rescale
       that threshold too, by passing an explicit tolerance, and say so.

    2. **The rescaled LOO tolerance is WEAK at the margin and that is not
       hidden.** On the same ladder at a native 204 dpi (tolerance 4.412 px):

           seeded tick error   LOO worst   vs 1.5   vs 4.412
                    +1.00 px       1.333     pass      pass
                    +2.00 px       2.667     FAIL      pass
                    +3.25 px       4.333     FAIL      **pass**

       A 3.25 px tick error slips through by 0.08 px. The rescale is correct in
       units and still costs real sensitivity, so it is not a free move: seed
       the widened tolerance **in the panel's own geometry**, because the
       numbers above depend on decade spacing and do not transfer.
    """
    if not native_dpi or native_dpi <= 0:
        raise ValueError(
            "axis.native_tol: native_dpi must be measured from the embedded "
            "image (width_px / width_pt * 72), not assumed. A tolerance "
            "rescaled by a guessed dpi is worse than one that is simply wrong, "
            "because it looks derived.")
    return tol_px * float(render_dpi) / float(native_dpi)


def loo_residuals(ticks_px):
    """Each tick's residual against a fit that EXCLUDES it.

    WHY THIS EXISTS. The max-residual test in `calibrate` is weak against a
    SINGLE bad tick, because least squares drags the line toward the outlier
    and hides most of the error. Measured 2026-11-10 on a clean 5-tick,
    100 px-pitch ladder: moving ONE tick by +6 px gives a max residual of only
    2.40 px and PASSES a 2.5 px tolerance. An end tick has leverage 0.2, so
    four fifths of its error is absorbed into the fit.

    A leave-one-out residual has no such masking: refit on the other n-1 ticks
    and predict the held-out one. The same +6 px error shows up as ~7.5 px.
    """
    import numpy as _np
    t = _np.asarray(ticks_px, dtype=float).ravel()
    t = t[_np.argsort(t)]
    n = t.size
    out = []
    for i in range(n):
        idx = _np.delete(_np.arange(n, dtype=float), i)
        y = _np.delete(t, i)
        A = _np.column_stack([idx, _np.ones_like(idx)])
        coef, *_ = _np.linalg.lstsq(A, y, rcond=None)
        out.append(float(abs(t[i] - (coef[0] * i + coef[1]))))
    return out


def calibrate(ticks_px, tick_values, kind="log10",
              residual_tol_px=TICK_RESIDUAL_PX, forbid_target=None,
              axis="y", loo_tol_px=None):
    """Fit a pixel <-> unit mapping from DETECTED TICK POSITIONS ONLY.

    ticks_px      measured tick coordinates, in pixels, any order
    tick_values   the data value at each tick, same order. Either an explicit
                  per-tick list, or a single number read as a uniform step
                  (decades for kind="log10", units for "linear").
    kind          "log10" or "linear"
    forbid_target the quantity this calibration will be checked against, if
                  any. Passing it makes error class 27 impossible: if the fit
                  input is that same vector, this raises.
    axis          "y" (pixels increase downward, the image convention) or "x".

    Returns a Calibration. Raises CalibrationError rather than falling back.
    """
    ticks = np.asarray(ticks_px, dtype=float).ravel()
    if ticks.size < 3:
        raise CalibrationError(
            f"calibrate: {ticks.size} tick position(s) given, need >= 3. Two "
            "points define a line through anything, so a 2-point fit cannot "
            "detect that an anchor is a glyph centroid rather than a tick -- "
            "which is error class 20. Detect more ticks, or state in the "
            "folder's DECISIONS.md that this panel cannot be calibrated.")

    # --- class 27 guard -------------------------------------------------
    if forbid_target is not None:
        tgt = np.asarray(forbid_target, dtype=float).ravel()
        if tgt.size == ticks.size and np.allclose(np.sort(tgt), np.sort(ticks)):
            raise ValidationTargetError(
                "calibrate: the fit input IS the validation target. A scale "
                "fitted onto the quantity it will be checked against absorbs "
                "any uniform error in the read and reports a clean comparison "
                "-- that is error class 27, and it hid a ~1.1% bias across all "
                "48 of LL_2026's printed hit rates (F-006). Take the scale "
                "from tick marks, labelled gridlines or printed axis anchors.")

    order = np.argsort(ticks)
    ticks = ticks[order]

    # --- tick values ----------------------------------------------------
    if np.isscalar(tick_values) or np.asarray(tick_values).size == 1:
        step = float(np.asarray(tick_values).ravel()[0])
        if step == 0:
            raise CalibrationError("calibrate: tick step is 0.")
        vals = None
        tick_step = step
    else:
        vals = np.asarray(tick_values, dtype=float).ravel()[order]
        if vals.size != ticks.size:
            raise CalibrationError(
                f"calibrate: {ticks.size} ticks but {vals.size} tick values.")
        units = np.log10(vals) if kind == "log10" else vals
        gaps = np.diff(units)
        if np.any(gaps == 0):
            raise CalibrationError("calibrate: two ticks carry the same value.")
        tick_step = float(np.median(gaps))
        # The ticks must be EQUALLY spaced in data units for a single linear
        # fit to be the right model at all.
        if np.max(np.abs(gaps - tick_step)) > 1e-9 * max(1.0, abs(tick_step)):
            raise CalibrationError(
                f"calibrate: tick VALUES are not equally spaced in "
                f"{'log10 ' if kind == 'log10' else ''}units: gaps {gaps}. "
                "A single scale cannot describe this axis. If it is a broken "
                "axis, calibrate each segment separately.")

    idx = np.arange(ticks.size, dtype=float)

    # --- the fit: scale AND origin, together, from the ticks -------------
    # There is deliberately no way to supply the origin. F-033 is what happens
    # when the origin is maintained separately from the scale: the scale was
    # re-derived from ticks and the origin stayed as three literals.
    A = np.column_stack([idx, np.ones_like(idx)])
    coef, *_ = np.linalg.lstsq(A, ticks, rcond=None)
    fitted = A @ coef
    residual = float(np.max(np.abs(ticks - fitted)))

    if residual > residual_tol_px:
        raise CalibrationError(
            f"calibrate: tick-fit residual {residual:.2f} px exceeds "
            f"{residual_tol_px} px over {ticks.size} ticks at "
            f"{np.round(ticks, 2).tolist()}.\n"
            "Equally-valued ticks on a log or linear axis are equally spaced "
            "in PIXELS. Unequal spacing means these positions are not tick "
            "marks -- most often glyph centroids of the axis labels, whose "
            "spacing tracks the label text, not the axis. That is error "
            "class 20. Re-detect the ticks from the darkness-weighted centroid "
            "of the columns OUTSIDE the axis rule.")

    # Leave-one-out outlier test. Opt-in via loo_tol_px until the running
    # extraction batch drains, because tightening a shared check mid-batch is
    # the hazard that bit this project on 2026-11-10; pass LOO_TOL_PX to arm it.
    if loo_tol_px is not None and ticks.size >= 4:
        loo = loo_residuals(ticks)
        worst = max(loo)
        if worst > loo_tol_px:
            bad = int(max(range(len(loo)), key=lambda i: loo[i]))
            raise CalibrationError(
                f"calibrate: tick {bad} at {ticks[bad]:.2f} px is "
                f"{worst:.2f} px from where the OTHER {ticks.size - 1} ticks "
                f"predict it (tolerance {loo_tol_px} px), even though the "
                f"overall fit residual is only {residual:.2f} px.\n"
                "Least squares drags the line toward a single outlier and "
                "hides most of its error -- on a 5-tick ladder one tick can be "
                "6 px out and still pass the max-residual test. One tick is "
                "mis-detected: re-measure it, do not widen the tolerance.")

    px_per_tick = abs(float(coef[0]))
    if px_per_tick <= 0:
        raise CalibrationError("calibrate: degenerate fit, 0 px per tick.")

    px_per_unit = px_per_tick / abs(tick_step)

    # Pixel direction. On a y axis, pixel row increases downward while value
    # increases upward, so the sign is negative; on x, positive.
    if vals is None:
        # Uniform step: the first tick in pixel order is the anchor, and the
        # caller's step sign says which way values run.
        sign = -1.0 if axis == "y" else 1.0
        origin_px = float(ticks[0])
        origin_value = None
        raise CalibrationError(
            "calibrate: a uniform step alone does not pin the origin -- it "
            "gives the scale and leaves the anchor unknown, which is the F-033 "
            "shape. Pass explicit tick_values (e.g. [1e3, 1e4, 1e5, 1e6]) so "
            "the origin is fitted from the same data as the scale.")

    # Values ascending in pixel order means pixels run with values.
    sign = 1.0 if (vals[-1] > vals[0]) else -1.0
    # Anchor on the FITTED position of the first tick, not its measured one:
    # the whole point of the least-squares fit is that no single tick carries
    # the origin alone.
    origin_px = float(fitted[0])
    origin_value = float(vals[0])

    cal = Calibration(kind, px_per_unit, origin_px, origin_value,
                      residual, int(ticks.size), float(tick_step), sign)

    # --- round trip. A mapping that does not invert is not a mapping. ----
    back = cal.to_pixel(vals)
    rt = float(np.max(np.abs(np.asarray(back) - fitted)))
    if rt > 1e-6:
        raise CalibrationError(
            f"calibrate: round trip failed, {rt:.3g} px. Internal error.")
    return cal


def detect_ticks_1d(gray, along, scan_lo, scan_hi, dark=140, min_len=3,
                    max_len=40):
    """Tick centres from a dark-ink strip OUTSIDE the axis rule.

    `gray` is a 2-D array, dark ink low. For a y axis pass along="y" with
    [scan_lo, scan_hi) the COLUMN range left of the axis line; the returned
    values are row coordinates, as darkness-weighted centroids of each
    contiguous run.

    Darkness-weighted, not run-midpoint, because an anti-aliased 600 dpi tick
    is 2-4 px tall with unequal edges and the midpoint quantises to half a
    pixel. The weighting is what makes the uniform-spacing assertion in
    `calibrate` tight enough to catch a glyph centroid.
    """
    g = np.asarray(gray, dtype=float)
    strip = g[:, scan_lo:scan_hi] if along == "y" else g[scan_lo:scan_hi, :]
    prof = strip.min(axis=1) if along == "y" else strip.min(axis=0)
    ink = prof < dark
    weight = np.clip(dark - prof, 0, None)

    runs, start = [], None
    for i, on in enumerate(ink):
        if on and start is None:
            start = i
        elif not on and start is not None:
            runs.append((start, i))
            start = None
    if start is not None:
        runs.append((start, len(ink)))

    centres = []
    for a, b in runs:
        if not (min_len <= b - a <= max_len):
            continue
        w = weight[a:b]
        if w.sum() <= 0:
            continue
        centres.append(float(np.dot(np.arange(a, b), w) / w.sum()))
    return centres
