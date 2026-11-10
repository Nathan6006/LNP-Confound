#!/usr/bin/env python3
"""Stage 02 -- which amine head groups this SI actually determines, MEASURED.

The library is 4 thiolactones x 4 pyridyl disulfides x 18 amines (Scheme S2,
SI p.7). The SI enumerates the thiolactones and the disulfides in full. IT NEVER
ENUMERATES A1-A18: that panel is in the main article, which this curation does
not hold (00_acquisition.json). What the SI gives is

  A1   NAMED in Scheme S1c ("N1,N1-diethyl ethylene diamine"), DRAWN there in
       vector, and pinned a third time by the MALDI formula of T18U-PY16-A1.
  A5   DRAWN only -- Figure S1 top scheme and the Figure S1b inset, both raster.
  A10  DRAWN only -- Figure S1 top scheme and the Figure S1d inset, both raster.
  A8   molecular FORMULA only, C7H17N3, solved from two MALDI formulae.

Every count below is MEASURED: vector bond coordinates for A1, raster path
geometry for A5, raster polygon geometry for A10. CLAUDE.md section 5 / error
class 48 -- a chain counted by eye is this project's documented +/-1 failure
mode, and "N,N-diethylETHYLENEdiamine" vs "N,N-diethylPROPYLENEdiamine" is
exactly a +/-1 on a linker, with A1 and A5 as the two answers.

Where a measurement does not settle the question this file SAYS SO and returns
no number. A10's ring size and A5's two N-ethyls are both in that category.

Run:  python3 -I stage02_amines.py
"""
import json
import os
import sys

import numpy as np
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))

import pymupdf                                              # noqa: E402
import chains                                               # noqa: E402

SI = os.path.join(HERE, "source", "bc8b00007_si_001.pdf")
FIGS1_XREF = 179          # the single raster image object of SI p.15
DR = 10                   # red-minus-green / red-minus-blue floor
MAXSUM = 730              # darker than this is ink, not paper


def red_mask(rgb):
    """Pixels of the Figure S1 structures' red ink.

    The insets are drawn in red on a JPEG-compressed raster, so saturation is
    low -- measured ink values run from (255,204,207) down to (173,157,158).
    R-G and R-B both >= DR separates the red drawing from the black axes, the
    black data curve and the black 'A5'/'A10' labels in the same window.
    """
    a = rgb.astype(int)
    return ((a[:, :, 0] - a[:, :, 1] >= DR)
            & (a[:, :, 0] - a[:, :, 2] >= DR)
            & (a.sum(2) < MAXSUM))


def _path_along_x(mask, x0, x1, y0, y1):
    """Ink centreline per column, by MEDIAN.

    The mean is wrong here: this raster carries isolated stray red pixels
    (measured, e.g. a single pixel 14 rows off the stroke at x=706) and one of
    them moves a column's mean by more than half a bond's amplitude. The median
    of a 2-5 pixel stroke is the stroke.
    """
    out = []
    for x in range(x0, x1):
        ys = np.nonzero(mask[y0:y1, x])[0]
        if len(ys):
            out.append((float(x), float(np.median(ys) + y0)))
    return out


def measure_A1(pg):
    """Vector. The three substituent lengths around A1's tertiary N."""
    segs = chains.vector_bonds(pg, (235, 215, 305, 260))
    # The amine N's text box on this page is x 269.4-273.3, y 234.9-239.7.
    # A bond of the amine has exactly one endpoint within 4 pt of it. The 47-pt
    # reaction arrow in the same window is removed by the length filter.
    NBOX = (269.4, 234.9, 273.3, 239.7)

    def near_n(p):
        return (NBOX[0] - 4 <= p[0] <= NBOX[2] + 4
                and NBOX[1] - 4 <= p[1] <= NBOX[3] + 4)

    short = [s for s in segs
             if ((s[2] - s[0]) ** 2 + (s[3] - s[1]) ** 2) ** 0.5 < 12]
    first = [s for s in short if near_n((s[0], s[1])) != near_n((s[2], s[3]))]
    if len(first) != 3:
        raise SystemExit(
            f"stage02: A1's amine N has {len(first)} bonds in the drawing, "
            "not 3. Refusing to proceed on a tertiary-amine assumption.")

    hetero = [tuple(round(v, 1) for v in w[:4]) for w in pg.get_text("words")
              if w[4] in ("NH2", "NH", "N", "O", "S")
              and 215 < w[1] < 260 and 235 < w[0] < 310]

    branches = []
    for s in first:
        far = (s[2], s[3]) if near_n((s[0], s[1])) else (s[0], s[1])
        chain, used = [far], {id(s)}
        while True:
            nxt = None
            for t in short:
                if id(t) in used:
                    continue
                for a, b in (((t[0], t[1]), (t[2], t[3])),
                             ((t[2], t[3]), (t[0], t[1]))):
                    if (abs(a[0] - chain[-1][0]) < 1.2
                            and abs(a[1] - chain[-1][1]) < 1.2):
                        nxt = (t, b)
                        break
                if nxt:
                    break
            if not nxt:
                break
            used.add(id(nxt[0]))
            chain.append(nxt[1])
        # A branch that ENDS ON A HETEROATOM GLYPH has one fewer carbon than it
        # has vertices. A1's terminal NH2 is a text box on this page, so the
        # test is positional and checkable rather than an assumption about
        # which of the three branches is the linker.
        tip = chain[-1]
        on_glyph = any(b[0] - 4 <= tip[0] <= b[2] + 4
                       and b[1] - 4 <= tip[1] <= b[3] + 4 for b in hetero)
        branches.append({"vertices": len(chain),
                         "tip": [round(tip[0], 1), round(tip[1], 1)],
                         "tip_is_heteroatom_glyph": on_glyph,
                         "carbons": len(chain) - (1 if on_glyph else 0)})

    return {"method": "vector bond coordinates (curation/lib/chains."
                      "vector_bonds), SI p.2 Scheme S1c",
            "n_bonds_at_N": len(first),
            "heteroatom_glyph_boxes": hetero,
            "branches": branches,
            "branch_carbon_counts": sorted(b["carbons"] for b in branches)}


def measure_A5(mask):
    """Raster. The LINKER is measured; the two N-ethyls are NOT. Both reported.

    The linker is the +/-1-critical quantity -- a two-carbon linker would make
    A5 the same molecule as A1 -- so it is measured by a route with no tuned
    parameter fitted to the expected answer. The turning points of the zig-zag
    are the sign changes of the centreline's slope, and the only free number --
    the lag over which that slope is measured -- is SWEPT rather than chosen,
    with the modal count taken and the whole sweep reported. 7 of the 8 lags
    agree; the one that does not (lag 5) splits a single turn into three on
    sub-pixel jitter, which is exactly what a sweep is for.

    The two N-substituents are NOT measurable at this resolution and this
    function returns no number for them. Their ink sits at the JPEG noise
    floor: the red mask finds nothing at y=318-319 where a grey-level dump
    shows a stroke, so the recovered path is non-monotone and Douglas-Peucker
    returns 4 vertices for a one-corner path at every eps from 1.2 to 5.0.
    Reporting "2 carbons" from that would be a guess dressed as a measurement.
    """
    link = _path_along_x(mask, 639, 706, 345, 375)   # strictly between the
    # 'N' glyph (last ink column 635) and the 'NH2' glyph (first ink column 709)
    xs = np.array([p[0] for p in link])
    ys = np.array([p[1] for p in link])
    sweep, places = {}, {}
    for lag in range(3, 11):
        sign, cur = [], 0
        for i in range(len(ys) - lag):
            d = ys[i + lag] - ys[i]
            if d:
                cur = 1 if d > 0 else -1
            sign.append(cur)
        turns = [round(float(xs[i]), 1) for i in range(1, len(sign))
                 if sign[i] != sign[i - 1] and sign[i - 1] != 0]
        sweep[str(lag)] = len(turns)
        places[str(lag)] = turns
    counts = list(sweep.values())
    modal = max(set(counts), key=counts.count)
    agree = counts.count(modal)
    if agree < 6:
        raise SystemExit(
            f"stage02: A5's linker turning-point count is unstable across the "
            f"lag sweep ({sweep}). Refusing to pick one. A chain length is an "
            "integer and an unstable count means the measurement, not the "
            "molecule, is the problem.")
    return {
        "method": "raster: turning points of the linker's y(x) centreline, "
                  "counted as sign changes of y(x+lag)-y(x) with the lag SWEPT "
                  "from 3 to 10 px. No threshold is fitted to the expected "
                  "answer (error class 27) and no single lag is trusted "
                  "(error class 28): the count is the modal one and the sweep "
                  "is reported so a reviewer sees the spread.",
        "linker_path_points": len(link),
        "linker_amplitude_px": round(float(ys.max() - ys.min()), 2),
        "lag_sweep_turning_point_counts": sweep,
        "lag_sweep_turning_point_positions": places,
        "lags_agreeing": f"{agree}/8",
        "linker_carbons": modal,
        "linker_rdp_vertices_by_eps": {
            str(e): chains.raster_chain_length(link, eps_px=e)
            for e in (1.2, 2.0, 3.0, 4.0, 5.0)},
        "N_substituents_measured": False,
        "N_substituents_note":
            "NOT MEASURED -- below the raster's noise floor (see docstring). "
            "Read from the drawing as two ethyls, and flagged on every "
            "emitted row as structure_assumed = 1 naming this feature.",
    }


def measure_A10_ring(mask):
    """Raster. Ring size from the drawn polygon's own aspect ratio.

    A ring is not a chain, so raster_chain_length does not apply. The
    discriminator is the polygon's width-to-vertical-side ratio, fixed by the
    polygon and by nothing else: a regular n-gon drawn with vertical left and
    right sides has width/side = cos(pi/n)/sin(pi/n) x ... -> 1.6180 (n=5),
    1.7321 (n=6), 1.8019 (n=7).

    The vertical side's extent depends on how near the global extreme a row
    must be, and AT THIS NOISE FLOOR THAT TOLERANCE CHANGES THE ANSWER.
    Sweeping it is the honest form of the measurement, and the sweep is what
    makes this inconclusive rather than confidently wrong.
    """
    ys, xs = np.nonzero(mask[578:662, 495:582])
    ys, xs = ys + 578, xs + 495
    width = int(xs.max() - xs.min())
    cand = {"pentagon_5": 1.6180, "hexagon_6": 1.7321, "heptagon_7": 1.8019}
    sweep = {}
    for tol in (1, 2, 3, 4, 5):
        rows = sorted(y for y in set(ys) if xs[ys == y].min() <= xs.min() + tol)
        side = (rows[-1] - rows[0]) if len(rows) > 1 else 0
        r = (width / side) if side else None
        sweep[str(tol)] = {
            "side_px": int(side),
            "width_over_side": round(r, 4) if r else None,
            "nearest": (min(cand, key=lambda k: abs(cand[k] - r))
                        if r else None)}
    answers = sorted({v["nearest"] for v in sweep.values() if v["nearest"]})
    verdict = (f"{answers[0]} at every tolerance" if len(answers) == 1 else
               "INCONCLUSIVE -- the tolerance sweep returns "
               + "/".join(answers)
               + ", so the drawn ring size is NOT decided by this measurement")
    return {"method": "raster polygon aspect ratio (width / vertical side), "
                      "tolerance swept",
            "width_px": width, "candidates": cand,
            "tolerance_sweep": sweep, "verdict": verdict}


def main():
    doc = pymupdf.open(SI)
    pix = pymupdf.Pixmap(doc, FIGS1_XREF)
    img = np.array(Image.frombytes("RGB", (pix.width, pix.height),
                                   pix.samples))
    mask = red_mask(img)

    out = {
        "_what": "What this SI determines about the 18 amine head groups.",
        "n_amines_in_library": 18,
        "figS1_image": {"xref": FIGS1_XREF, "w": pix.width, "h": pix.height,
                        "red_ink_px": int(mask.sum())},
        "A1": measure_A1(doc[1]),
        "A5": measure_A5(mask),
        "A10": measure_A10_ring(mask),
    }
    a1 = out["A1"]["branch_carbon_counts"]
    if a1 != [2, 2, 2]:
        raise SystemExit(
            f"stage02: A1's measured branches are {a1}; "
            "N,N-diethylethylenediamine requires [2, 2, 2]. Refusing to emit.")
    out["A1"].update(
        name_printed="N1,N1-diethyl ethylene diamine (Scheme S1c, SI p.2)",
        SMILES="CCN(CC)CCN", enumerated=True, confidence="high",
        third_anchor="MALDI of T18U-PY16-A1, C44H88N3O2S2 [M+H]+ 754.62. "
                     "Subtracting T18U (C22H39NO2S) and the hexadecyl arm "
                     "(C16H33S, replacing the thiol hydrogen) leaves C6H16N2 "
                     "for the amine -- this structure's formula.")
    if out["A5"]["linker_carbons"] != 3:
        raise SystemExit(
            f"stage02: A5's measured linker is "
            f"{out['A5']['linker_carbons']} carbons, not 3. Refusing to emit.")
    out["A5"].update(
        name_printed=None, SMILES="CCN(CC)CCCN", enumerated=True,
        confidence="medium", third_anchor=None,
        why="The linker is measured (3 carbons, so A5 is NOT A1). The two "
            "N-ethyls are read from the drawing and could not be measured at "
            "this raster's resolution, and the paper gives A5 no name, no "
            "formula and no mass, so there is no anchor orthogonal to the "
            "drawing at all. Mass cannot see connectivity -- and here there "
            "is no mass either.")
    out["A10"].update(
        name_printed=None, SMILES=None, enumerated=False, confidence=None,
        candidate_SMILES="CCN1CCNCC1",
        why="NOT ENUMERATED. Drawn only, and the one geometric measurement "
            "available at this raster resolution is INCONCLUSIVE between a "
            "5-, 6- and 7-membered ring once its tolerance is swept. "
            "1-ethylpiperazine is recorded as a CANDIDATE, not as a "
            "structure. Nothing is lost: A10 appears in NO row of Table S2. "
            "DECISIONS.md D-3.")
    out["A8"] = {
        "SMILES": None, "enumerated": False,
        "molecular_formula": "C7H17N3",
        "method": "solved from two independent MALDI formulae -- "
                  "T18U-PY12-A8 = C41H81N4O2S2 and T14-PY12-A8 = "
                  "C37H75N4O2S2. Their difference, C4H6, is exactly "
                  "T18U - T14, which is the cross-check that the solve is "
                  "self-consistent.",
        "confidence": None,
        "why_not_enumerated":
            "A FORMULA IS NOT A STRUCTURE. C7H17N3 admits "
            "1-(2-aminoethyl)-4-methylpiperazine, N-(3-aminopropyl)-N,N'-"
            "dimethylethylenediamine and many other triamines; nothing in "
            "the SI chooses between them, and A8 is the second-best and "
            "third-best lipidoid in the paper. DECISIONS.md D-2.",
    }
    out["n_enumerable_amines"] = sum(
        1 for k in ("A1", "A5", "A10") if out[k].get("enumerated"))
    out["summary"] = {
        "A1_branch_carbons": a1,
        "A5_linker_carbons": out["A5"]["linker_carbons"],
        "A5_lag_sweep": out["A5"]["lag_sweep_turning_point_counts"],
        "A5_lags_agreeing": out["A5"]["lags_agreeing"],
        "A5_N_substituents_measured": out["A5"]["N_substituents_measured"],
        "A10_verdict": out["A10"]["verdict"],
        "amines_enumerable": out["n_enumerable_amines"],
        "amines_in_library": 18,
    }
    with open(os.path.join(HERE, "01_amines.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps(out["summary"], indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
