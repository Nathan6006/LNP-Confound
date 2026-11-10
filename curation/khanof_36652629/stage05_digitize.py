"""Stage 05 -- read the six in-vivo bar panels of PMID 36652629.

TIER 5 (pixels). 00_value_sources.json records why tiers 1-4 and 4b all fail:
no flux value is typeset anywhere in the 43,853-character SI text layer, and
the figure pages carry 0 filled paths and 0 rectangle items, i.e. the panels
were flattened to raster before deposit.

GEOMETRY, and why it is shaped this way.

Every bar is a PALE FILL inside a SATURATED STROKE of the same hue, and the
stroke is THICK -- measured at 16 to 24 px at 600 dpi, i.e. 2 to 3 pt, which is
13 to 19 per cent of a bar's width on each side. Overlaid on the bar are
replicate dots and an error-bar cap in that same stroke colour, and in this
paper the dots frequently sit ON THE BAR'S TOP CORNERS (the orange C3-K2-A14
bar of Fig S1b carries one at each corner, ~20 px above the top rule).

Three things follow, each learned by a refusal rather than chosen up front.

1. THE MASK IS NOT A COLOUR MATCH. A first version matched the bar's modal
   fill RGB within a tolerance. It refused on Fig S3b because a thick stroke
   takes 33 per cent of the bar width, so the FILL alone can never reach the
   90 per cent full-width coverage bars.bar_top requires. The mask used here is
   instead "ink that belongs to a bar": not white, and either coloured or
   bright. It needs no per-bar colour and no tolerance, so it cannot delete a
   bar the way a chroma ceiling can.

2. A SATURATION TEST CANNOT FIND THESE BARS. Fig S1b's teal fill is
   (215, 219, 218): a channel spread of FOUR. Any saturation threshold able to
   see it also sees grey antialiasing. The first version used spread > 28 and
   deleted every lavender and teal bar; lowering it to 8 still split the teal
   bar of Fig S1b into two runs at the one column where the spread dipped to 4.
   The brightness branch of the mask is what fixes that.

3. THE FLANK TEST IS NOT ENOUGH HERE, so there is a second check.
   bars.bar_top rejects a CENTRED marker by requiring the flanks to agree --
   but two dots placed symmetrically on both top corners agree with each other
   perfectly, so they would pass. What separates them from a real top rule is
   that a bar is SOLID from its top rule all the way down to the baseline,
   while a dot cluster has white underneath it. So after bar_top returns, this
   script requires the mask to stay above MIN_INSIDE_COVER for EVERY row from
   the detected top to the baseline, and refuses otherwise. That check is
   seeded in 07_qc.json: a synthetic dot pair is pasted above a known bar and
   the refusal is confirmed to fire.

THE VALUE IS THE STROKE'S CENTRELINE. bar_top on this mask returns the stroke's
OUTER edge. GraphPad centres the stroke on the data path, so

    value_row = outer_row + t/2

with `t` measured from the bar's OWN SIDE WALLS -- the gap between the mask's
outer edge and the fill colour's edge -- in a row band in the bottom 30 per
cent of the panel, where no dot, cap or whisker reaches. Left and right walls
are measured separately and must agree; an asymmetric pen is F-031.

Nothing here lowers a threshold to make a refusal pass.
"""
import os, sys, json, math
import numpy as np
import pymupdf

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import axis, bars                                            # noqa: E402
from errors import DetectionError                            # noqa: E402

PDF = os.path.join(HERE, "source", "ja2c10670_si_001.pdf")
DPI = 600

DARK = 110            # a pixel this dark is axis/text ink
WHITE_MIN = 242       # a pixel whose darkest channel is at least this is paper
SPREAD = 8            # channel spread above this is "coloured"
PALE_MIN = 150        # a FILL pixel's darkest channel is at least this
BAND_FRAC = 0.30      # bottom fraction of the panel used for colour + wall work
MIN_INSIDE_COVER = 0.80   # mask coverage required at EVERY row inside a bar
WALL_ASYM_TOL_PX = 5.0    # left vs right stroke width
RULE_THICK_FRAC = 0.70    # a top rule must be this much of the side-wall width


# ---------------------------------------------------------------------------
# Panel register. The clip is in PDF POINTS and is deliberately generous: the
# spine, baseline and ticks are all DETECTED inside it, never assumed. `labels`
# is left-to-right as the panel's own category axis prints them.
# ---------------------------------------------------------------------------
PANELS = [
    dict(key="S1b", page=2, clip=(443, 150, 522, 228),
         decades=[8, 7, 6, 5, 4, 3],
         labels=["SM-102", "C3-K2-A14", "C3-K2-E14"],
         round_id="R1_hbond", cargo="unmodified Fluc mRNA",
         n_text="n = 4 hindlimbs", caption_page=2),
    dict(key="S2b", page=3, clip=(443, 137, 520, 214),
         decades=[6, 5, 4, 3],
         labels=["C3-K2-E14", "C3-K2-E10", "C3-K2-E15"],
         round_id="R2_tail_length", cargo="unmodified Fluc mRNA",
         n_text="n = 4 or 5 hindlimbs", caption_page=3),
    dict(key="S3b", page=3, clip=(321, 456, 381, 540),
         decades=[7, 6, 5, 4, 3],
         labels=["C3-K2-E10", "C3-K3-E10"],
         round_id="R3_linker", cargo="unmodified Fluc mRNA",
         n_text="n = 4 or 6 hindlimbs", caption_page=4),
    dict(key="S4b", page=4, clip=(132, 260, 237, 336),
         decades=[8, 7, 6, 5, 4],
         labels=["SM-102", "C3-K3-E10", "OC2-K3-E10", "OC3-K3-E10"],
         round_id="R4_core_unmod", cargo="unmodified Fluc mRNA",
         n_text="n = 8 or 10 hindlimbs", caption_page=4),
    dict(key="S4c", page=4, clip=(268, 260, 325, 336),
         decades=[8, 7, 6, 5, 4],
         labels=["SM-102", "OC2-K3-E10"],
         round_id="R5_core_modmRNA", cargo="modified Fluc mRNA",
         n_text="n = 4", caption_page=4),
    dict(key="S7a", page=7, clip=(188, 106, 290, 200),
         decades=[8, 7, 6, 5],
         labels=["DLin-MC3-DMA", "SM-102", "OC2-K3-E10", "OC2-K2-A10"],
         round_id="R6_NP5_benchmarks", cargo="unmodified Fluc mRNA",
         n_text="n = 8 or 12 hindlimbs", caption_page=7),
]


# ------------------------------------------------------------------ render --
def native_dpi(page, clip):
    """Resolution of the EMBEDDED image this clip falls inside.

    Needed because every panel is a flattened raster at ~190-221 dpi. Asking a
    600 dpi re-render for tick positions good to 1.5 px is asking for 0.55 of a
    REAL pixel, which no correct detection can deliver; axis.LOO_TOL_PX is
    therefore applied in the image's own pixels and converted here. This is a
    measured property of the file, not a tolerance anybody chose.
    """
    doc = pymupdf.open(PDF)
    pg = doc[page - 1]
    best = None
    for im in pg.get_images(full=True):
        for r in pg.get_image_rects(im[0]):
            if (r.x0 - 1 <= clip[0] and r.y0 - 1 <= clip[1]
                    and r.x1 + 1 >= clip[2] and r.y1 + 1 >= clip[3]):
                dx = im[2] / (r.x1 - r.x0) * 72.0
                dy = im[3] / (r.y1 - r.y0) * 72.0
                best = (dx + dy) / 2.0
    if best is None:
        raise DetectionError(
            f"native_dpi: no embedded image on page {page} contains clip "
            f"{clip}. Refusing to assume a resolution.")
    return best


def render(page, clip, dpi=DPI):
    doc = pymupdf.open(PDF)
    pg = doc[page - 1]
    pix = pg.get_pixmap(dpi=dpi, clip=pymupdf.Rect(*clip))
    a = np.frombuffer(pix.samples, dtype=np.uint8)
    a = a.reshape(pix.height, pix.width, pix.n)[:, :, :3].astype(int)
    return a


def bar_ink(rgb):
    """True where a pixel belongs to a BAR -- fill, stroke, dot or cap alike.

    Not white, and either coloured or bright. No per-bar colour, no tolerance:
    a chroma ceiling is what silently deleted a whole bar in a sibling folder,
    and a saturation floor is what split this paper's teal bar in two.
    """
    mx, mn = rgb.max(2), rgb.min(2)
    return (mn < WHITE_MIN) & (((mx - mn) > SPREAD) | (mn >= PALE_MIN))


def find_frame(rgb):
    """Spine column band and baseline row band, from the longest dark runs."""
    dark = rgb.max(2) < DARK
    H, W = dark.shape
    colsum, rowsum = dark.sum(0), dark.sum(1)

    cand = np.flatnonzero(colsum >= 0.55 * colsum.max())
    if cand.size == 0:
        raise DetectionError("find_frame: no vertical spine.")
    spine_cols = cand[cand <= cand.min() + 25]
    spine_l, spine_r = int(spine_cols.min()), int(spine_cols.max())

    rcand = np.flatnonzero(rowsum >= 0.55 * rowsum.max())
    rcand = rcand[rcand > H // 2]
    if rcand.size == 0:
        raise DetectionError("find_frame: no horizontal baseline below mid.")
    base_rows = rcand[rcand >= rcand.max() - 25]
    return dict(spine_l=spine_l, spine_r=spine_r,
                baseline_top=int(base_rows.min()),
                baseline_bot=int(base_rows.max()))


def ticks(rgb, fr):
    """Decade tick centres, from the dark strip LEFT of the spine."""
    gray = rgb.min(2)
    lo = max(0, fr["spine_l"] - 14)
    return sorted(axis.detect_ticks_1d(gray, "y", lo, fr["spine_l"],
                                       dark=140, min_len=2, max_len=30))


def bar_columns(ink, fr):
    """Contiguous column runs that carry a bar, right of the spine."""
    bot = fr["baseline_top"]
    count = ink[:bot, :].sum(0)
    on = count > 0.10 * bot
    on[:fr["spine_r"] + 2] = False
    runs, start = [], None
    for i, v in enumerate(on):
        if v and start is None:
            start = i
        elif not v and start is not None:
            runs.append((start, i)); start = None
    if start is not None:
        runs.append((start, len(on)))
    return [r for r in runs if r[1] - r[0] >= 12]


def modal(px):
    cols, cnt = np.unique(px, axis=0, return_counts=True)
    return tuple(int(v) for v in cols[int(np.argmax(cnt))])


def near(rgb, colour, tol):
    return np.abs(rgb - np.array(colour)[None, None, :]).max(2) <= tol


def walls(rgb, ink, dark, x0, x1, r0, r1):
    """Stroke width from the bar's own side walls, in a clean row band.

    Measured as the distance from the OUTER edge of the dark side wall to the
    first pale FILL pixel inside it. Anchoring on the dark wall rather than on
    the column run matters: the JPEG leaves a pale halo 1-3 px outside every
    bar, (255, 251, 232) against a (253, 228, 198) fill, and a first version
    that measured from the run boundary read that halo as a 3 px stroke on one
    side and 19 px on the other -- caught by the asymmetry guard below.
    """
    sub = rgb[r0:r1, x0:x1]
    mn = sub.min(2)
    m = ink[r0:r1, x0:x1]
    d = dark[r0:r1, x0:x1]
    pale = m & (mn >= PALE_MIN)
    if pale.sum() < 200:
        raise DetectionError(
            f"walls: only {int(pale.sum())} pale pixels in rows[{r0},{r1}) "
            f"cols[{x0},{x1}). Refusing.")
    fill_c = modal(sub[pale])
    fmask = near(sub, fill_c, 30)
    W = x1 - x0
    lefts, rights = [], []
    for r in range(fmask.shape[0]):
        dc = np.flatnonzero(d[r])
        if dc.size < 4:
            continue
        dl, dr = int(dc.min()), int(dc.max())
        fc = np.flatnonzero(fmask[r])
        fc = fc[(fc > dl) & (fc < dr)]
        if fc.size < 0.3 * W:
            continue
        lefts.append(int(fc.min()) - dl)
        rights.append(dr - int(fc.max()))
    if len(lefts) < 20:
        raise DetectionError("walls: too few clean rows to measure a stroke.")
    tl, tr = float(np.median(lefts)), float(np.median(rights))
    if abs(tl - tr) > WALL_ASYM_TOL_PX:
        raise DetectionError(
            f"walls: left stroke {tl:.1f} px, right {tr:.1f} px, asymmetry "
            f"{abs(tl-tr):.1f} px > {WALL_ASYM_TOL_PX}. A stroke centred on "
            "the drawn path is the same width on both walls; a difference "
            "this large means one wall is being read against a halo, a "
            "neighbour or a marker. Refusing (the F-031 shape).")
    return dict(fill_rgb=fill_c, t_left=round(tl, 2), t_right=round(tr, 2),
                t=round((tl + tr) / 2.0, 2), n_rows=len(lefts))


def trim_to_walls(dark, x0, x1, r0, r1):
    """Trim a column run to the bar's TRUE outer edges, from its dark walls.

    The JPEG leaves a pale halo 2-3 px outside every bar and bar_ink includes
    it, so the raw run is a few pixels wider than the bar. bars.bar_top then
    asks for 90 per cent coverage of an inflated width and the top rule, which
    covers the BAR, tops out at 0.92-0.95 -- enough for the full-width rule to
    pass while a flank fails. That is how a 3 px halo turns into a refusal.

    This is NOT the inner trimming that caused F-018. It moves the window
    OUTWARD to the bar's real edge, measured from the dark side walls in a
    clean band; it never narrows the window inside the bar.
    """
    firsts, lasts = [], []
    for r in range(r0, r1):
        c = np.flatnonzero(dark[r, x0:x1])
        if c.size < 4:
            continue
        firsts.append(int(c.min())); lasts.append(int(c.max()))
    if len(firsts) < 20:
        raise DetectionError("trim_to_walls: too few rows with a dark wall.")
    nx0 = x0 + int(round(float(np.median(firsts))))
    nx1 = x0 + int(round(float(np.median(lasts)))) + 1
    if nx1 - nx0 < 12 or nx0 < x0 or nx1 > x1:
        raise DetectionError(
            f"trim_to_walls: walls give [{nx0},{nx1}) from run [{x0},{x1}). "
            "Refusing.")
    return nx0, nx1


def rule_thickness(dark, x0, x1, top, bot, frac=0.90):
    """Contiguous rows from `top` whose dark coverage is at least `frac`."""
    n = 0
    r = top
    while r < bot and dark[r, x0:x1].mean() >= frac:
        n += 1; r += 1
    return n


def solid_top(ink, x0, x1, bot):
    """Topmost row from which the bar is SOLID all the way to the baseline.

    This is the paper-specific step, and it is what makes bars.bar_top usable
    here. Scanning for a top rule from row 0 walks into the replicate dots: in
    five of the six panels the dots sit on or above the bar's top corners, and
    two dots placed symmetrically satisfy bars.bar_top's flank test perfectly.
    What distinguishes them from a rule is underneath, not beside: a bar is
    continuous down to the baseline and a dot cluster is not. So the search
    band starts here, and the library then does the full-width and flank work
    on the real rule.
    """
    prof = ink[:bot - 2, x0:x1].mean(1)
    r = bot - 3
    if prof[r] < MIN_INSIDE_COVER:
        raise DetectionError(
            f"solid_top: the row just above the baseline is only {prof[r]:.3f} "
            "bar ink; this window is not sitting on a bar.")
    while r > 0 and prof[r - 1] >= MIN_INSIDE_COVER:
        r -= 1
    return int(r)


# -------------------------------------------------------------------- main --
def do_panel(p):
    rgb = render(p["page"], p["clip"])
    nat = native_dpi(p["page"], p["clip"])
    ink = bar_ink(rgb)
    fr = find_frame(rgb)
    tk = ticks(rgb, fr)

    want = len(p["decades"])
    if len(tk) != want:
        raise DetectionError(
            f"{p['key']}: detected {len(tk)} ticks at "
            f"{np.round(tk,1).tolist()}, expected {want} for decades "
            f"{p['decades']}. REFUSING -- do not drop a tick to fit.")

    # Ticks ascend in ROW; the panel's decades therefore DESCEND in value.
    tick_px = sorted(tk)
    values = [10.0 ** d for d in sorted(p["decades"], reverse=True)]

    # axis.LOO_TOL_PX is 1.5 px of the ORIGINAL raster. Expressed in the pixels
    # of this 600 dpi re-render of a `nat` dpi image that is 1.5*600/nat.
    loo_tol = axis.LOO_TOL_PX * DPI / nat
    cal = axis.calibrate(tick_px, values, kind="log10", axis="y",
                         loo_tol_px=loo_tol)
    loo = axis.loo_residuals(tick_px)

    runs = bar_columns(ink, fr)
    if len(runs) != len(p["labels"]):
        raise DetectionError(
            f"{p['key']}: found {len(runs)} bar column runs {runs}, expected "
            f"{len(p['labels'])} for {p['labels']}. REFUSING.")

    bot = fr["baseline_top"]
    r0 = int(bot * (1.0 - BAND_FRAC))
    r1 = bot - 6

    # bars.FLANK_TOL_PX is 2 px OF THE ORIGINAL RASTER. These panels are
    # flattened images at 188-221 dpi, so 2 px of them is 2*DPI/nat px of this
    # 600 dpi re-render. Asking for 2 render px would be asking the left and
    # right ends of a rule to agree to 0.7 of a real pixel, which JPEG blur on
    # a 188 dpi panel cannot deliver. Same reasoning, same source, as the LOO
    # tolerance above; both are reported in BOTH units in 05_digitize.json.
    flank_tol = bars.FLANK_TOL_PX * DPI / nat

    out = []
    dark = ink & (rgb.min(2) < PALE_MIN)
    for (rx0, rx1), lab in zip(runs, p["labels"]):
        x0, x1 = trim_to_walls(dark, rx0, rx1, r0, r1)
        band_cov = float(ink[r0:r1, x0:x1].mean())
        if band_cov < 0.95:
            raise DetectionError(
                f"{p['key']} / {lab}: the measuring band rows[{r0},{r1}) is "
                f"only {band_cov:.3f} bar ink. Refusing.")
        w = walls(rgb, ink, dark, x0, x1, r0, r1)
        t = w["t"]
        need = max(3, int(round(RULE_THICK_FRAC * t)))

        st = solid_top(ink, x0, x1, bot)

        # Find the TOP RULE, not the error-bar cap. bars.bar_top does the
        # full-width and flank work; this loop supplies the acceptance test it
        # cannot have: a top rule is as thick as the side wall that is drawn
        # with the same pen, and a whisker cap is not. Every rejection is
        # recorded, so a reviewer sees what was skipped and why.
        y = st
        rejected = []
        while True:
            res = bars.bar_top(dark, x0, x1, y, bot, baseline_row=bot,
                               flank_tol_px=flank_tol)
            thick = rule_thickness(dark, x0, x1, res["row"], bot)
            if thick >= need:
                break
            rejected.append({"row": int(res["row"]), "thick_px": int(thick),
                             "needed_px": need,
                             "why": "thinner than the side wall drawn with "
                                    "the same pen; this is the error-bar cap"})
            y = res["row"] + max(3, thick)
            if len(rejected) > 6:
                raise DetectionError(
                    f"{p['key']} / {lab}: no candidate rule reached "
                    f"{need} px thickness after {len(rejected)} tries "
                    f"{rejected}. REFUSING.")
        top = res["row"]

        prof = ink[top:bot - 2, x0:x1].mean(1)
        worst = float(prof.min()) if prof.size else 0.0
        if worst < MIN_INSIDE_COVER:
            bad = int(top + int(np.argmin(prof)))
            raise DetectionError(
                f"{p['key']} / {lab}: mask coverage falls to {worst:.3f} at "
                f"row {bad} between the detected top {top} and the baseline "
                f"{bot}. A real bar is solid all the way down. REFUSING "
                "(error class 32).")

        row = top + t / 2.0
        val = float(cal.to_value(row))
        out.append(dict(
            label=lab, run=[int(rx0), int(rx1)], x0=int(x0), x1=int(x1),
            width_px=int(x1 - x0), fill_rgb=w["fill_rgb"],
            solid_top_row=st, outer_row=int(top),
            rule_thickness_px=int(rule_thickness(dark, x0, x1, top, bot)),
            rule_thickness_needed_px=need,
            rejected_candidates=rejected,
            stroke_t_px=t, stroke_t_left_px=w["t_left"],
            stroke_t_right_px=w["t_right"], wall_rows_used=w["n_rows"],
            flank_delta_px=res["flank_delta_px"],
            flank_delta_native_px=round(res["flank_delta_px"] * nat / DPI, 3),
            coverage=res["coverage"],
            band_cover=round(band_cov, 4), min_inside_cover=round(worst, 4),
            value_row=round(row, 2), value=val,
            log10_value=round(math.log10(val), 4)))

    return dict(
        key=p["key"], page=p["page"], round_id=p["round_id"],
        cargo=p["cargo"], n_text=p["n_text"], clip_pt=p["clip"],
        render_dpi=DPI, native_image_dpi=round(nat, 2),
        frame=fr, tick_rows=[round(v, 2) for v in tick_px],
        tick_values=values, calibration=cal.qc(),
        loo_residual_px=[round(v, 3) for v in loo],
        loo_worst_px=round(max(loo), 3),
        loo_tol_px_used=round(loo_tol, 3),
        loo_worst_native_px=round(max(loo) * nat / DPI, 3),
        loo_tol_native_px=axis.LOO_TOL_PX,
        flank_tol_px_used=round(flank_tol, 3),
        flank_tol_native_px=bars.FLANK_TOL_PX,
        px_per_decade=round(cal.px_per_unit, 3),
        bars=out)


def main():
    res, fails = [], []
    for p in PANELS:
        try:
            r = do_panel(p)
            res.append(r)
            print(f"\n=== {r['key']}  ({r['round_id']})  "
                  f"{r['px_per_decade']:.2f} px/decade  "
                  f"native {r['native_image_dpi']} dpi  "
                  f"resid {r['calibration']['tick_fit_residual_px']} px  "
                  f"LOO worst {r['loo_worst_px']} px "
                  f"({r['loo_worst_native_px']} native, tol "
                  f"{r['loo_tol_native_px']}) ===")
            for b in r["bars"]:
                print(f"   {b['label']:14s} {b['value']:12.4g} p/s   "
                      f"log10 {b['log10_value']:7.4f}   "
                      f"stroke {b['stroke_t_px']:5.1f}px "
                      f"(L {b['stroke_t_left_px']:4.1f} R {b['stroke_t_right_px']:4.1f})  "
                      f"solid>={b['min_inside_cover']:.3f}")
        except Exception as e:
            fails.append({"panel": p["key"], "error": str(e)})
            print(f"\n=== {p['key']} REFUSED ===\n{e}")

    doc = {"_what": "Stage 05 geometry and per-bar reads. Tier 5 (pixels).",
           "dpi": DPI, "panels": res, "refusals": fails}
    with open(os.path.join(HERE, "05_digitize.json"), "w") as f:
        json.dump(doc, f, indent=1)
    print(f"\npanels read {len(res)} / {len(PANELS)}; refusals {len(fails)}")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
