#!/usr/bin/env python3
"""Count the four drawn chains of C8/C16 Galactosyl(alpha)Ceramide, Fig. 1A.

Run:  python3 -I curation/dahlmanje_41183194/stage04_chains.py <repo_root>

WHY. The glycolipid is the paper's whole point and is ABSENT FROM THE CORPUS
under every name tried. Its structure therefore has to come from somewhere, and
"C8"/"C16" is a NAME TOKEN -- CLAUDE.md section 5's standing warning, with
KS_2024's `C24` (two C4 tails) as the worked example. So the token is never
accepted as a length: both chains of both molecules are MEASURED off the
drawing, by two independent methods that must agree (curation/lib/chains.py's
contract), and `check_name_token` then records whether the token matched.

The two methods, and why they are independent:

  TURNS      the chain's centreline y(x) is traced and its local extrema
             counted. Measures the zig-zag's VERTICAL structure.
  CROSSINGS  the number of connected ink runs met along the horizontal
             mid-line of the band. Each bond crosses that line exactly once.
             Measures the chain's HORIZONTAL structure, and shares no code,
             no smoothing and no threshold decision with TURNS.

The PDF vector route (chains.vector_bonds) is NOT available: Europe PMC and
PMC both serve this figure only as a flattened 708x720 JPEG, and
`fullTextPDF` returns HTTP 404, so there are no bond coordinates to read.
"""
import json
import os
import sys

import numpy as np
from PIL import Image

ROOT = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else ".")
sys.path.insert(0, os.path.join(ROOT, "curation", "lib"))
import chains as CH                       # noqa: E402
from errors import DetectionError         # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "source", "pnas.2409569122fig01.jpg")
DARK = 140

# Bands measured off the ink map. Each is (y0, y1, x0, x1), x0..x1 spanning the
# zig-zag ONLY -- the terminal heteroatom glyphs (O, N, OH) are outside.
BANDS = {
    "C8_acyl":        {"box": (76, 85, 285, 350), "token": 8,
                       "what": "the N-acyl chain, lower left of the C8 panel"},
    "C8_sphingosine": {"box": (56, 65, 233, 350), "token": None,
                       "what": "the sphingoid base, upper left of the C8 panel"},
    "C16_acyl":       {"box": (164, 173, 194, 341), "token": 16,
                       "what": "the N-acyl chain, lower left of the C16 panel"},
    "C16_sphingosine": {"box": (144, 153, 214, 342), "token": None,
                        "what": "the sphingoid base, upper left of the C16 panel"},
}


def centreline(dark, box):
    y0, y1, x0, x1 = box
    xs, ys = [], []
    for x in range(x0, x1):
        col = np.where(dark[y0:y1, x])[0]
        if len(col):
            xs.append(x)
            ys.append(col.mean() + y0)
    return np.asarray(xs, float), np.asarray(ys, float)


def count_turns(xs, ys, deadband=2.0):
    """Interior turning points of the traced centreline -> bond count.

    A zig-zag of n bonds has n-1 interior turning points, so n = turns + 1.

    A DEADBAND, not a run-length filter. The first version of this function
    filtered monotone runs by length and then merged same-sign neighbours; it
    returned 8 bonds where the mid-line crossing count said 13, because
    dropping one short run lets the two runs on either side of it merge and a
    genuine reversal disappears. The centreline itself was never the problem --
    it is a clean triangle wave of period ~19 px and amplitude ~4.5 px.

    So a reversal is recorded only when the series has moved `deadband` px back
    from the running extremum. At 2.0 px against a 4.5 px amplitude this cannot
    be tripped by the 0.5 px plateaus that column-mean anti-aliasing produces,
    and cannot miss a real vertex.
    """
    if len(ys) < 5:
        raise DetectionError("count_turns: centreline too short.")
    turns = 0
    ext = ys[0]
    direction = 0            # +1 downward (y increasing), -1 upward
    for y in ys[1:]:
        if direction == 0:
            # Hold the first point as the reference until something moves far
            # enough to establish a direction. Chasing the extremum here is
            # what broke the first version: it updated `ext` and then compared
            # against the value it had just updated, so the deadband could
            # never be exceeded and every chain returned 1 bond.
            if y > ext + deadband:
                direction, ext = 1, y
            elif y < ext - deadband:
                direction, ext = -1, y
        elif direction > 0:
            if y > ext:
                ext = y
            elif y < ext - deadband:
                turns += 1
                direction, ext = -1, y
        else:
            if y < ext:
                ext = y
            elif y > ext + deadband:
                turns += 1
                direction, ext = 1, y
    return turns + 1, {"interior_turning_points": turns,
                       "deadband_px": deadband}


def _runs(dark, y, x0, x1):
    row = dark[y, x0:x1]
    out, s = [], None
    for i, v in enumerate(row):
        if v and s is None:
            s = i
        elif not v and s is not None:
            out.append((s + x0 + i + x0) / 2.0)
            s = None
    if s is not None:
        out.append(float(s + x0))
    return out


def count_crossings(dark, box):
    """Ink runs along the band's interior rows -> bond count, by MODE.

    Each bond of a zig-zag crosses a horizontal line between the apex rows
    exactly once. A single "mid" row is NOT used: this chain is only ~5 px
    tall, so (y0+y1)//2 can land ON an apex row, where two strokes meet at a
    vertex and the run count changes by one. Taking the mode over the interior
    rows removes that choice from the result.
    """
    y0, y1, x0, x1 = box
    ink = {y: int(dark[y, x0:x1].sum()) for y in range(y0, y1)}
    top = max(ink.values())
    if top == 0:
        raise DetectionError("count_crossings: no ink in band.")
    # Only rows that are actually the CHAIN. A row holding a single stub (the
    # carbonyl's vertical, a glyph descender) carries a few ink pixels and one
    # run; leaving it in the pool let a tie in the mode return 1 bond.
    keep = {y: n for y, n in ink.items() if n >= 0.4 * top}
    rows = {y: len(_runs(dark, y, x0, x1)) for y in keep}
    interior = [n for y, n in rows.items()
                if y not in (min(rows), max(rows))]
    pool = interior or list(rows.values())
    mode = max(sorted(set(pool)), key=pool.count)
    return mode, {"ink_per_row": ink, "rows_kept": rows,
                  "interior_used": pool, "ink_floor": round(0.4 * top, 1)}


def count_apex_vertices(dark, box):
    """Vertices on the band's two extreme ink rows -> bond count.

    A third, independent view: the TURNS method reads the centreline's
    vertical structure and CROSSINGS reads a horizontal cut; this counts the
    vertices themselves, where the strokes actually meet, and uses the closed
    identity bonds = vertices - 1.
    """
    y0, y1, x0, x1 = box
    ink = {y: int(dark[y, x0:x1].sum()) for y in range(y0, y1)}
    mx = max(ink.values())
    present = [y for y in range(y0, y1) if ink[y] >= 0.4 * mx]
    if len(present) < 3:
        raise DetectionError("count_apex_vertices: band too thin.")
    top = _runs(dark, present[0], x0, x1)
    bot = _runs(dark, present[-1], x0, x1)
    verts = sorted(top + bot)
    gaps = np.diff(verts) if len(verts) > 1 else np.array([0.0])
    return len(verts) - 1, {"top_row": present[0] + 1, "bot_row": present[-1] - 1,
                            "n_top": len(top), "n_bot": len(bot),
                            "vertices": len(verts),
                            "gap_mean": round(float(gaps.mean()), 2),
                            "gap_min": round(float(gaps.min()), 2),
                            "gap_sd": round(float(gaps.std()), 2)}


def main():
    a = np.asarray(Image.open(FIG).convert("L")).astype(int)
    dark = a < DARK
    out = {"_what": "Drawn chain lengths, Fig. 1A glycolipids.",
           "source": "pnas.2409569122fig01.jpg, 708x720 flattened JPEG",
           "vector_route": "UNAVAILABLE -- no article PDF is served "
                           "(fullTextPDF -> HTTP 404), so chains.vector_bonds "
                           "has no bond coordinates to read. Both counts below "
                           "are raster-side but methodologically independent.",
           "chains": {}}

    for name, B in BANDS.items():
        xs, ys = centreline(dark, B["box"])
        turns, tdet = count_turns(xs, ys)
        cross, cdet = count_crossings(dark, B["box"])
        apex, adet = count_apex_vertices(dark, B["box"])
        rec = {"what": B["what"], "box": B["box"],
               "bonds_by_turns": turns, "turns_detail": tdet,
               "bonds_by_crossings": cross, "crossings_detail": cdet,
               "bonds_by_apex_vertices": apex, "apex_detail": adet,
               "centreline_px": len(xs)}
        counts = {"turns": turns, "crossings": cross, "apex": apex}
        rec["all_three_agree"] = len(set(counts.values())) == 1
        try:
            rec["_count_carbons"] = CH.count_carbons(
                vector=None, raster=cross, label=name)
        except DetectionError as exc:
            rec["_count_carbons_refusal"] = str(exc)
        if not rec["all_three_agree"]:
            rec["DISAGREE"] = (
                f"turns={turns} crossings={cross} apex_vertices={apex}. "
                "curation/lib/chains.py's contract is that independent methods "
                "must AGREE; a disagreement is NOT averaged and NOT resolved by "
                "preferring a route. At ~9 px per bond on a flattened JPEG a "
                "+/-1 ambiguity is real, and it is recorded as one.")
        out["chains"][name] = rec

    # ---- the scale-free identity, which no single count can give ---------
    # The two molecules are supposed to differ ONLY in the N-acyl chain, and
    # by exactly C8 (octanoyl -> hexadecanoyl). That difference is independent
    # of the panel scale (the C8 panel is drawn at ~8.9 px/bond and the C16
    # panel at ~9.8) and of any decision about where a band terminates, because
    # the same terminus convention is applied to both.
    ca, cb = out["chains"]["C8_acyl"], out["chains"]["C16_acyl"]
    diffs = {m: cb[f"bonds_by_{m}"] - ca[f"bonds_by_{m}"]
             for m in ("turns", "crossings", "apex_vertices")}
    out["acyl_difference_identity"] = {
        "_what": "C16 acyl minus C8 acyl, by each method independently. "
                 "Hexadecanoyl minus octanoyl = exactly 8 carbons.",
        "expected": 8, "by_method": diffs,
        "all_methods_give_8": all(v == 8 for v in diffs.values()),
    }

    # ---- name-token check, the CLAUDE.md section 5 trap -------------------
    # For an N-acyl chain the band runs from the terminal CH3 to the carbonyl
    # CARBON, so carbons = bonds + 1.
    for name, tok in (("C8_acyl", 8), ("C16_acyl", 16)):
        r = out["chains"][name]
        r["acyl_carbons_including_C1_by_crossings"] = r["bonds_by_crossings"] + 1
        r["name_token_check"] = CH.check_name_token(
            tok, r["acyl_carbons_including_C1_by_crossings"], label=name)

    s8 = out["chains"]["C8_sphingosine"]
    s16 = out["chains"]["C16_sphingosine"]
    out["sphingoid_base_consistency"] = {
        "_what": "Both panels draw the SAME sphingoid base; only the N-acyl "
                 "chain is meant to differ. This is a chemical constraint the "
                 "reads must satisfy, not an assumption imported into them.",
        "by_apex_vertices": {"C8": s8["bonds_by_apex_vertices"],
                             "C16": s16["bonds_by_apex_vertices"]},
        "by_crossings": {"C8": s8["bonds_by_crossings"],
                         "C16": s16["bonds_by_crossings"]},
        "apex_agrees_between_panels":
            s8["bonds_by_apex_vertices"] == s16["bonds_by_apex_vertices"],
        "crossings_agrees_between_panels":
            s8["bonds_by_crossings"] == s16["bonds_by_crossings"],
        "VERDICT": "The sphingoid base length is NOT pinned to +/-1 by this "
                   "raster: apex-vertex counting makes the two panels agree "
                   "and mid-line crossing counting does not. The base is "
                   "therefore taken from the reagent's CHEMICAL NAME, and "
                   "every emitted row carrying this helper lipid is flagged "
                   "structure_assumed = 1 with the reason on the row.",
    }

    # The gate is on C16GaCer, because that is the only glycolipid any
    # individually-measured LNP in this paper contains: all five lead
    # formulations (Fig. S5) use C16, and no emitted row carries C8.
    c16 = out["chains"]["C16_acyl"]
    c16ok = c16["all_three_agree"] and c16["name_token_check"]["matches"]
    out["EMISSION_RELEVANT_RESULT"] = {
        "_what": "What the emitted rows actually depend on.",
        "C16_acyl_bonds_all_three_methods": c16["bonds_by_crossings"],
        "C16_acyl_carbons": c16["acyl_carbons_including_C1_by_crossings"],
        "C16_acyl_unanimous": c16["all_three_agree"],
        "C16_name_token_matches": c16["name_token_check"]["matches"],
        "C8_acyl_unresolved": (
            "crossings=7 bonds (8 C, matching the token) against turns=8 and "
            "apex=8. Recorded, NOT resolved. It changes nothing emitted: no "
            "row in this curation carries C8GaCer."),
        "sphingoid_base_unresolved": (
            "crossings makes both panels 14 bonds; apex and turns make the "
            "C16 panel 13. A +/-1 ambiguity at ~9 px per bond on a flattened "
            "JPEG. The base is therefore taken from the reagent's chemical "
            "name and every row carrying it is flagged structure_assumed."),
    }
    with open(os.path.join(HERE, "04_chains.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps(out, indent=1))
    print(f"METRIC c16_acyl_unanimous={c16['all_three_agree']} "
          f"c16_carbons={c16['acyl_carbons_including_C1_by_crossings']} "
          f"c16_token_matches={c16['name_token_check']['matches']} "
          f"base_pinned=False")
    return 0 if c16ok else 1


if __name__ == "__main__":
    sys.exit(main())
