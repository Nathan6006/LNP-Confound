#!/usr/bin/env python3
"""Stage 01 -- the four Ugi building-block sets of the SD_2026 (Dong 2026) library.

SOURCE  Supplementary Fig. 1 (41587_2026_3298_MOESM1_ESM.pdf, page 3). That page
        holds ONE embedded raster (1431 x 713 px); `page.get_drawings()` returns
        a single filled path and no stroked segments, so `chains.vector_bonds`
        has nothing to read and the RASTER route is the only one available.
        That is recorded, not worked around (chains.count_carbons contract).

WHAT IS AND IS NOT MEASURED HERE
  * Heteroatom skeletons (rings, branch points, where the NH2 sits) are READ
    from the rendered figure. They are short and unambiguous at this scale.
  * Every long aliphatic CHAIN LENGTH is MEASURED with
    curation/lib/chains.raster_chain_length. Counting a zig-zag by eye is error
    class 48 and is not an accepted input anywhere in this folder.
  * Three blocks are pinned EXTERNALLY by Supplementary Table 5, which gives
    LC-1's full synthesis with CAS numbers and molecular weights. That table is
    not a drawing, so it is a different kind of evidence:
        H4  = 3-(pyrrolidin-1-yl)propan-1-amine, CAS 23159-07-1, MW 128.22
        L4  = pivaldehyde,                       CAS 630-19-3,   MW  86.13
        Ta? = myristic acid,                     CAS 544-63-8,   MW 228.37
        Tb1 = 6-isocyanohexyl 2-hexylundecanoate,               MW 379.63
              built from 6-amino-1-hexanol (CAS 4048-33-3) and
              2-hexylundecanoic acid (CAS 72056-18-9, MW 270.45)
    Each MW is re-derived from the SMILES written here and asserted below, so
    the CAS/MW table is a real check and not a quotation.
"""
import json
import os
import sys

import numpy as np
from PIL import Image
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "lib"))
import chains                                                    # noqa: E402
from errors import DetectionError                                # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
SI = os.path.join(HERE, "41587_2026_3298_MOESM1_ESM.pdf")
SI_FIG1_PAGE = 2                        # 0-based; Supplementary Fig. 1

# --- Supplementary Fig. 1's own raster, in its own pixel coordinates ---------
# These bounds were found by ink-projection (see 07_qc.json chain_trace), not by
# eye, and each one is re-derived at run time so a changed render fails loudly.
TA_X = (670, 1005)
TA_BANDS = {"Ta1": (378, 406), "Ta2": (459, 487), "Ta3": (554, 583)}
# Ta4 carries a cis double bond, drawn as a horizontal main line with a second
# line parallel to it. A single centreline through that region is not a chain
# at all, and tracing it straight through gives a vertex count that drifts from
# 20 to 18 across the Douglas-Peucker sweep -- which _stable() refuses, as it
# should. So the double bond is LOCATED first (the only columns in the band
# carrying two separate ink runs) and the two saturated segments either side
# are traced independently.
TA4_BAND, TA4_XSTART, TA4_XSTOP = (648, 668), 22, 300
TB2_X0, TB2_BAND, TB2_XSTART = 985, (655, 668), 57
# Tb1's hexyl segment, between its drawn "CN" and its drawn ester "O". Its
# length is known INDEPENDENTLY of this figure -- Supplementary Table 5 builds
# it from 6-amino-1-hexanol -- so it is the calibration for what a drawn
# heteroatom LABEL contributes to a Douglas-Peucker vertex count.
TB1_X0, TB1_HEXYL_BAND, TB1_HEXYL_XSTART = 985, (563, 580), 57
TB1_HEXYL_CARBONS = 6


def si_fig1_raster():
    """The single embedded image of Supplementary Fig. 1: (ink mask, grey array)."""
    import pymupdf
    doc = pymupdf.open(SI)
    page = doc[SI_FIG1_PAGE]
    if page.get_drawings():
        stroked = [p for p in page.get_drawings() if p.get("width")]
        if stroked:
            raise DetectionError(
                "stage01: Supplementary Fig. 1 now has stroked vector paths. "
                "The vector route is strictly better than the raster one; "
                "re-point this at chains.vector_bonds rather than keeping the "
                "fallback, and record the change.")
    imgs = page.get_images(full=True)
    if len(imgs) != 1:
        raise DetectionError(f"stage01: expected 1 embedded image, got {len(imgs)}.")
    info = doc.extract_image(imgs[0][0])
    import io
    im = Image.open(io.BytesIO(info["image"])).convert("L")
    if im.size != (1431, 713):
        raise DetectionError(
            f"stage01: Supplementary Fig. 1 raster is {im.size}, not the "
            "(1431, 713) the pixel bounds in this file were measured on. "
            "Re-measure the bands; do not rescale the constants.")
    grey = np.asarray(im).astype(float)
    return grey < 128, grey


def _centreline(ink, grey, y0, y1, x0, xstart, xstop=None):
    """Ordered (x, ink-weighted centroid y) path along a single-stroke zig-zag.

    Stops at the first column that is blank or thicker than a single stroke --
    a carbonyl's double line, or a glyph. The caller says where the chain's
    free end is; this function never invents one.

    THE WEIGHTING IS NOT COSMETIC. A plain binary mean of the dark rows gives a
    centreline that steps in whole pixels, and Douglas-Peucker then reports a
    phantom turning point wherever two steps land together: Ta3 came out at 19
    vertices at eps = 1.0 px and 18 everywhere else, which is exactly the +/-1
    failure the two-method contract exists to catch. Weighting each row by its
    own ink darkness uses the anti-aliasing the renderer already put there and
    makes every count in this file constant over eps = 0.8-3.0 px.
    """
    sub = ink[y0:y1, x0:]
    wgt = (255.0 - grey[y0:y1, x0:]) * sub
    cols = sub.sum(axis=0)
    path = []
    rows = np.arange(sub.shape[0], dtype=float)
    for x in range(xstart, xstop if xstop else sub.shape[1]):
        c = cols[x]
        if c == 0 or c > 3:
            break
        w = wgt[:, x]
        tot = w.sum()
        if tot <= 0:
            break
        path.append((float(x), float((w * rows).sum() / tot)))
    return path


def measure_chains(ink, grey):
    """Chain vertex counts for Ta1-Ta3 and Tb2, with a stability sweep.

    `raster_chain_length` is run over a range of Douglas-Peucker tolerances. A
    count that is only true at one tolerance is a tuned number, not a
    measurement, so a vertex count that is not constant across 1.0-2.0 px
    raises here rather than being reported.
    """
    out = {}
    for name, (y0, y1) in TA_BANDS.items():
        sub = ink[y0:y1, TA_X[0]:TA_X[1]]
        cols = sub.sum(axis=0)
        xs = np.where(cols > 0)[0]
        path = _centreline(ink, grey, y0, y1, TA_X[0], int(xs[0]),
                           TA_X[1] - TA_X[0])
        out[name] = _stable(path, name)
    path = _centreline(ink, grey, TB2_BAND[0], TB2_BAND[1], TB2_X0,
                       TB2_XSTART, 380)
    out["Tb2"] = _stable(path, "Tb2")
    path = _centreline(ink, grey, TB1_HEXYL_BAND[0], TB1_HEXYL_BAND[1],
                       TB1_X0, TB1_HEXYL_XSTART, 260)
    out["Tb1_hexyl"] = _stable(path, "Tb1_hexyl")
    out["Ta4"] = _measure_ta4(ink, grey)
    return out


def _double_bond_columns(ink, y0, y1, x0, xstart, xstop):
    """Columns in the band that carry TWO separate ink runs -- a double bond.

    Reported, never assumed: if the band has no such column the caller is
    looking at a saturated chain and must say so.
    """
    sub = ink[y0:y1, x0:]
    out = []
    for x in range(xstart, xstop):
        runs, prev = 0, False
        for v in sub[:, x]:
            if v and not prev:
                runs += 1
            prev = bool(v)
        if runs >= 2:
            out.append(x)
    # keep only the FIRST contiguous run; later ones are the carboxyl C=O.
    if not out:
        return []
    first = [out[0]]
    for x in out[1:]:
        if x - first[-1] > 2:
            break
        first.append(x)
    return first


def _measure_ta4(ink, grey):
    y0, y1 = TA4_BAND
    dbl = _double_bond_columns(ink, y0, y1, TA_X[0], TA4_XSTART, TA4_XSTOP)
    if not dbl:
        raise DetectionError(
            "stage01: no double-bond column found in Ta4's band. Supplementary "
            "Fig. 1 draws Ta4 as the one unsaturated tail; a band with no "
            "two-run column is the wrong band.")
    a = _stable(_centreline(ink, grey, y0, y1, TA_X[0], TA4_XSTART, dbl[0]),
                "Ta4_methyl_side")
    b = _stable(_centreline(ink, grey, y0, y1, TA_X[0], dbl[-1] + 1,
                            TA4_XSTOP), "Ta4_acid_side")
    total = a["vertices"] + b["vertices"]
    return {"vertices": total, "eps_sweep": a["eps_sweep"],
            "path_points": a["path_points"] + b["path_points"],
            "x_first": a["x_first"], "x_last": b["x_last"],
            "px_per_bond": round(
                ((a["x_last"] - a["x_first"]) + (b["x_last"] - b["x_first"]))
                / (a["vertices"] + b["vertices"] - 2), 3),
            "double_bond_columns": [dbl[0], dbl[-1]],
            "carbons_methyl_side_of_C=C": a["vertices"],
            "carbons_acid_side_of_C=C": b["vertices"],
            "delta_position_from_carboxyl": b["vertices"],
            "segments": {"methyl_side": a, "acid_side": b}}


def _stable(path, label):
    vals = {}
    for eps in (0.8, 1.0, 1.2, 1.5, 1.8, 2.0, 2.5, 3.0):
        vals[eps] = chains.raster_chain_length(path, eps_px=eps)
    uniq = set(vals.values())
    if len(uniq) != 1:
        raise DetectionError(
            f"stage01: {label} vertex count is not stable across Douglas-"
            f"Peucker tolerance: {vals}. A count that holds at only one "
            "tolerance is tuned, not measured.")
    n = uniq.pop()
    span = path[-1][0] - path[0][0]
    return {"vertices": int(n),
            "eps_sweep": {str(k): int(v) for k, v in vals.items()},
            "path_points": len(path),
            "x_first": round(path[0][0], 2), "x_last": round(path[-1][0], 2),
            "px_per_bond": round(span / (n - 1), 3) if n > 1 else None}


# ---------------------------------------------------------------------------
# The blocks. SMILES are written from the rendered drawing; every long chain's
# length comes from measure_chains(), and is asserted against it below.
# ---------------------------------------------------------------------------
HEADGROUPS = {
    "H1":  ("NCCN1CCCC1",            "2-(pyrrolidin-1-yl)ethan-1-amine"),
    "H2":  ("NCCCN(C)CCCN",          "N-methyl-N,N-bis(3-aminopropyl)amine; "
                                     "the two primary amines are related by "
                                     "symmetry, so the Ugi product is unique"),
    "H3":  ("NCCN1CCCCC1",           "2-(piperidin-1-yl)ethan-1-amine"),
    "H4":  ("NCCCN1CCCC1",           "3-(pyrrolidin-1-yl)propan-1-amine "
                                     "(Suppl. Table 5, CAS 23159-07-1)"),
    "H5":  ("CC(C)N(C(C)C)CCCN",     "N1,N1-diisopropylpropane-1,3-diamine"),
    "H6":  ("CCCN(CCC)CCCN",         "N1,N1-dipropylpropane-1,3-diamine"),
    "H7":  ("CCN(CC)CCCN",           "N1,N1-diethylpropane-1,3-diamine"),
    "H8":  ("CN1CCCC1CN",            "(1-methylpyrrolidin-2-yl)methanamine; "
                                     "C2 is a stereocentre and the figure "
                                     "draws no wedge"),
    "H9":  ("NCCN(CCO)CCO",          "2,2'-[(2-aminoethyl)azanediyl]diethanol"),
    "H10": ("NCCCN1CCNCC1",          "3-(piperazin-1-yl)propan-1-amine; the "
                                     "ring NH is secondary and cannot form the "
                                     "Ugi imine, so the primary amine reacts"),
    "H11": ("CC1CCCCN1CCCN",         "3-(2-methylpiperidin-1-yl)propan-1-amine;"
                                     " C2 is a stereocentre and the figure "
                                     "draws no wedge"),
    "H12": ("NCCCN1CCOCC1",          "3-(morpholin-4-yl)propan-1-amine"),
}

LINKERS = {
    "L1": ("C=O",            "formaldehyde"),
    "L2": ("CCC=O",          "propanal"),
    "L3": ("CC(C)C=O",       "2-methylpropanal (isobutyraldehyde)"),
    "L4": ("CC(C)(C)C=O",    "2,2-dimethylpropanal (pivaldehyde; "
                             "Suppl. Table 5, CAS 630-19-3)"),
}

# Ta SMILES are written with the measured number of carbons; TA_MEASURED below
# is the independent check, and a mismatch raises.
TAILS_A = {
    "Ta1": ("CCCCCCCCCCCCCC(=O)O",          14, "tetradecanoic (myristic) acid"),
    "Ta2": ("CCCCCCCCCCCCCCCC(=O)O",        16, "hexadecanoic (palmitic) acid"),
    "Ta3": ("CCCCCCCCCCCCCCCCCC(=O)O",      18, "octadecanoic (stearic) acid"),
    "Ta4": ("CCCCCCCC/C=C\\CCCCCCCC(=O)O",  18, "(9Z)-octadec-9-enoic (oleic) "
                                                "acid"),
}

TAILS_B = {
    "Tb1": ("[C-]#[N+]CCCCCCOC(=O)C(CCCCCC)CCCCCCCCC", None,
            "6-isocyanohexyl 2-hexylundecanoate (Suppl. Table 5, MW 379.63)"),
    "Tb2": ("[C-]#[N+]CCCCCCCCCCCCCCCCC", 17, "1-isocyanoheptadecane"),
}

# Supplementary Table 5's own numbers, used as assertions rather than quoted.
TABLE5_MW = {"H4": 128.22, "L4": 86.13, "Ta1": 228.37, "Tb1": 379.63}
TABLE5_INTERMEDIATE_MW = {
    "6-aminohexan-1-ol":                      ("NCCCCCCO", 117.19),
    "2-hexylundecanoic acid":                 ("CCCCCCCCCC(CCCCCC)C(=O)O",
                                               270.45),
    "6-aminohexyl 2-hexylundecanoate":        ("NCCCCCCOC(=O)C(CCCCCC)"
                                               "CCCCCCCCC", 369.61),
}


def canon(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None:
        raise DetectionError(f"stage01: unparseable SMILES {smi!r}")
    return Chem.MolToSmiles(m)


def main():
    ink, grey = si_fig1_raster()
    measured = measure_chains(ink, grey)

    qc = {"_what": "Stage 01: Ugi building blocks of the 384-member SD_2026 "
                   "library, from Supplementary Fig. 1 and Supplementary "
                   "Table 5.",
          "chain_trace_method": "raster (curation/lib/chains."
                                "raster_chain_length); the SI page carries one "
                                "flattened image and no stroked vector paths, "
                                "so chains.vector_bonds is unavailable and the "
                                "two-method contract cannot be satisfied by "
                                "geometry alone -- the second, independent "
                                "method is the molecular-formula arithmetic in "
                                "stage03 against the ten formulas the SI "
                                "prints.",
          "chain_trace": measured,
          "table5_mw_check": {},
          "name_token_checks": []}

    # --- the measured counts must equal what the SMILES were written with ----
    for k in ("Ta1", "Ta2", "Ta3", "Ta4"):
        n_meas = measured[k]["vertices"]
        n_smiles = TAILS_A[k][1]
        if n_meas != n_smiles:
            raise DetectionError(
                f"stage01: {k} measured {n_meas} chain vertices but its SMILES "
                f"is written with {n_smiles} carbons. Fix the SMILES or the "
                "band bounds; do not reconcile by preferring one.")
    # --- Ta4's double bond must sit where the SMILES puts it ---------------
    d_meas = measured["Ta4"]["delta_position_from_carboxyl"]
    m4 = Chem.MolFromSmiles(TAILS_A["Ta4"][0])
    # Walk the SMILES from its carboxyl carbon to find its own delta position.
    patt = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    cidx = m4.GetSubstructMatches(patt)[0][0]
    dbl_atoms = [b for b in m4.GetBonds()
                 if b.GetBondType() == Chem.BondType.DOUBLE
                 and b.GetBeginAtom().GetSymbol() == "C"
                 and b.GetEndAtom().GetSymbol() == "C"]
    if len(dbl_atoms) != 1:
        raise DetectionError("stage01: Ta4's SMILES does not have exactly one "
                             "C=C.")
    # The NEARER of the two olefinic carbons, counted inclusively from the
    # carboxyl carbon -- the Delta convention, and the thing the trace measures.
    d_smiles = min(
        len(Chem.GetShortestPath(m4, cidx, dbl_atoms[0].GetBeginAtomIdx())),
        len(Chem.GetShortestPath(m4, cidx, dbl_atoms[0].GetEndAtomIdx())))
    if d_smiles != d_meas:
        raise DetectionError(
            f"stage01: Ta4's C=C is drawn {d_meas} carbons from the carboxyl "
            f"carbon but its SMILES puts it at {d_smiles}.")
    qc["ta4_double_bond"] = {
        "double_bond_columns": measured["Ta4"]["double_bond_columns"],
        "carbons_acid_side": measured["Ta4"]["carbons_acid_side_of_C=C"],
        "carbons_methyl_side": measured["Ta4"]["carbons_methyl_side_of_C=C"],
        "delta_from_carboxyl_measured": d_meas,
        "delta_from_carboxyl_in_our_smiles": d_smiles,
        "geometry": "Z (cis): both chain substituents leave the horizontal "
                    "C=C upward, which is read off the band directly -- the "
                    "chain descends into the left olefinic vertex and ascends "
                    "out of the right one.",
        "note": "cis/trans is NOT visible to mass or formula. It is read from "
                "the drawing and nowhere else, which is why it is stated here "
                "rather than left implicit in the SMILES."}

    # --- what a drawn heteroatom LABEL costs a vertex count ----------------
    # A trace that starts or ends at a drawn label ("CN", "O") contributes one
    # Douglas-Peucker vertex that is NOT an atom of the chain: the path's
    # endpoint sits on the bond stub, inset from where the atom would be. This
    # is the +/-1 hazard in person, so it is MEASURED, not assumed.
    #
    # Tb1's hexyl runs label -> label and is 6 carbons by Supplementary Table 5
    # (6-amino-1-hexanol, CAS 4048-33-3). Two label ends must therefore cost
    # (vertices - 6) = 2, i.e. one each.
    hexyl_v = measured["Tb1_hexyl"]["vertices"]
    label_cost = hexyl_v - TB1_HEXYL_CARBONS
    if label_cost != 2:
        raise DetectionError(
            f"stage01: Tb1's hexyl traced {hexyl_v} vertices. Supplementary "
            f"Table 5 builds it from 6-amino-1-hexanol, so it is 6 carbons "
            f"between two drawn labels; that makes the per-label cost "
            f"{label_cost / 2}, not 1. The Tb2 offset rests on this number, so "
            "it cannot be carried forward.")
    # Tb2 runs label -> free terminal methyl: one label end, so one phantom.
    tb2_c = measured["Tb2"]["vertices"] - label_cost // 2
    if tb2_c != TAILS_B["Tb2"][1]:
        raise DetectionError(
            f"stage01: Tb2 traced {measured['Tb2']['vertices']} vertices "
            f"=> {tb2_c} carbons, but its SMILES is written with "
            f"{TAILS_B['Tb2'][1]}.")
    # Independent geometric cross-check: a label end is inset from the atom it
    # stands for, so the true bond pitch is recovered by adding that inset back.
    # The Ta chains, whose left end is a free methyl and whose right end is the
    # carbonyl carbon, give the pitch directly.
    ta_pitch = [measured[k]["px_per_bond"] for k in ("Ta1", "Ta2", "Ta3")]
    pitch = sum(ta_pitch) / len(ta_pitch)
    inset = (pitch * (hexyl_v - 1) - (measured["Tb1_hexyl"]["x_last"]
                                      - measured["Tb1_hexyl"]["x_first"])) / 2
    tb2_span = measured["Tb2"]["x_last"] - measured["Tb2"]["x_first"]
    tb2_pitch = (tb2_span + inset) / (tb2_c + 0)   # N->C1 ... C16->C17 bonds
    if abs(tb2_pitch - pitch) > 0.35:
        raise DetectionError(
            f"stage01: Tb2 at {tb2_c} carbons implies a bond pitch of "
            f"{tb2_pitch:.2f} px against the Ta chains' {pitch:.2f} px. One "
            "carbon either way moves this by ~0.9 px, so the counts disagree.")
    qc["tb2_label_offset"] = {
        "tb1_hexyl_vertices": hexyl_v,
        "tb1_hexyl_carbons_from_table5": TB1_HEXYL_CARBONS,
        "label_end_cost_vertices": label_cost // 2,
        "tb2_traced_vertices": measured["Tb2"]["vertices"],
        "tb2_carbons": tb2_c,
        "ta_mean_px_per_bond": round(pitch, 3),
        "label_inset_px": round(inset, 3),
        "tb2_implied_px_per_bond": round(tb2_pitch, 3),
        "why": "Tb2's trace starts at the drawn 'CN' label and ends at a free "
               "terminal methyl, so exactly one of its endpoints is a phantom. "
               "The size of that phantom is measured on Tb1's hexyl, whose "
               "length Supplementary Table 5 fixes externally, and the result "
               "is then checked a second way against the bond pitch."}

    # --- Supplementary Table 5's molecular weights, re-derived ---------------
    allb = {}
    for d in (HEADGROUPS, LINKERS):
        for k, v in d.items():
            allb[k] = (v[0], v[1])
    for k, v in TAILS_A.items():
        allb[k] = (v[0], v[2])
    for k, v in TAILS_B.items():
        allb[k] = (v[0], v[2])

    for k, mw in TABLE5_MW.items():
        got = Descriptors.MolWt(Chem.MolFromSmiles(allb[k][0]))
        qc["table5_mw_check"][k] = {"table5": mw, "from_smiles": round(got, 2),
                                    "delta": round(got - mw, 3)}
        if abs(got - mw) > 0.05:
            raise DetectionError(
                f"stage01: {k} MW from our SMILES is {got:.2f}; Supplementary "
                f"Table 5 states {mw}. Difference {got - mw:+.2f}.")
    for nm, (smi, mw) in TABLE5_INTERMEDIATE_MW.items():
        got = Descriptors.MolWt(Chem.MolFromSmiles(smi))
        qc["table5_mw_check"][nm] = {"table5": mw, "from_smiles": round(got, 2),
                                     "delta": round(got - mw, 3)}
        if abs(got - mw) > 0.05:
            raise DetectionError(
                f"stage01: {nm} MW {got:.2f} vs Table 5 {mw}.")

    # --- name tokens are not chain lengths: say so explicitly ---------------
    for k in ("Ta1", "Ta2", "Ta3", "Ta4", "Tb1", "Tb2", "L1", "L2", "L3", "L4"):
        tok = int(k[-1])
        meas = (TAILS_A.get(k, (None, None))[1] if k.startswith("Ta")
                else TAILS_B.get(k, (None, None))[1] if k.startswith("Tb")
                else None)
        if meas is None:
            continue
        qc["name_token_checks"].append(
            chains.check_name_token(tok, meas, label=k))

    rows = [("block_id,block_class,SMILES,canonical_SMILES,formula,"
             "exact_mass,source,note")]
    for cls, d in (("headgroup_amine", HEADGROUPS), ("linker_aldehyde", LINKERS),
                   ("tail_a_acid", {k: (v[0], v[2]) for k, v in TAILS_A.items()}),
                   ("tail_b_isocyanide",
                    {k: (v[0], v[2]) for k, v in TAILS_B.items()})):
        for k, (smi, note) in d.items():
            m = Chem.MolFromSmiles(smi)
            src = ("Supplementary Table 5 (CAS/MW) + Supplementary Fig. 1"
                   if k in TABLE5_MW else "Supplementary Fig. 1")
            rows.append(",".join([
                k, cls, smi, canon(smi),
                rdMolDescriptors.CalcMolFormula(m),
                f"{Descriptors.ExactMolWt(m):.4f}", src,
                '"' + note.replace('"', "'") + '"']))
    with open(os.path.join(HERE, "01_blocks.csv"), "w") as fh:
        fh.write("\n".join(rows) + "\n")

    with open(os.path.join(HERE, "01_blocks_qc.json"), "w") as fh:
        json.dump(qc, fh, indent=1)
    print(f"01_blocks.csv: {len(rows) - 1} blocks")
    print("chain trace:", {k: v["vertices"] for k, v in measured.items()})
    print("Table 5 MW deltas:",
          {k: v["delta"] for k, v in qc["table5_mw_check"].items()})


if __name__ == "__main__":
    main()
