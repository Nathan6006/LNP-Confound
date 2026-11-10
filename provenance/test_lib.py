#!/usr/bin/env python3
"""Seeded-mutation tests for curation/lib. EXITS NON-ZERO on failure.

    python3 -I provenance/test_lib.py

CLAUDE.md section 6: "A NEW CHECK IS NOT DONE UNTIL IT HAS BEEN MADE TO FAIL ON
PURPOSE." Three incidents paid for that rule -- a re-copied C10 anchor that
could only pass, three ZT_2026 marginal checks that scored a literal 0 for
twelve days while the summary beside them said "ALL FOUR EXACT", and two
versions of C16 that reported "all derived" while examining nothing. Every one
printed a clean line.

So this file is the other half of curation/lib. Two kinds of test:

  MUTATION  seed the exact defect an invariant exists for; the module MUST
            refuse. A test that cannot fail proves nothing.
  CONTROL   the correct input; the module must NOT refuse. Without controls a
            module that refuses everything passes every mutation.

  REPLAY    the historical defect itself, reconstructed -- not a stand-in. The
            F-013 replay rebuilds a sterol ester through RWMol.RemoveAtom on an
            atom adjacent to C3 and confirms the product is epicholesterol:
            same formula, same InChIKey SKELETON block, different stereo block.
            A seeded mutation proves a check CAN fire; a replay proves it fires
            on the thing that actually happened.

Add a test here for every new invariant. Numbering is P-<n> for a replay and
M-<n> for a mutation; take max+1 and check it is unused (the rule from the
class 18/19 collision).
"""
import csv
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
LIB = os.path.join(HERE, "..", "curation", "lib")
sys.path.insert(0, LIB)

import axis            # noqa: E402
import bars            # noqa: E402
import heatmap         # noqa: E402
import chains          # noqa: E402
import assemble        # noqa: E402
import reagents        # noqa: E402
import emit            # noqa: E402
from errors import (AmbiguityError, CalibrationError,  # noqa: E402,F401
                    DetectionError, EmitError, StructureError,
                    ValidationTargetError)

from rdkit import Chem  # noqa: E402

CHOL = ("CC(C)CCC[C@@H](C)[C@H]1CC[C@H]2[C@@H]3CC=C4C[C@@H](O)CC[C@]4(C)"
        "[C@H]3CC[C@]12C")

RESULTS = []


def check(tid, what, fn, expect, note=""):
    """expect: an exception class to require, or None to require success."""
    try:
        val = fn()
    except Exception as e:                       # noqa: BLE001
        if expect is None:
            RESULTS.append((tid, what, "FAIL",
                            f"refused a CORRECT input: {type(e).__name__}: "
                            f"{str(e).splitlines()[0]}"))
            return None
        if isinstance(e, expect):
            RESULTS.append((tid, what, "PASS",
                            f"refused: {str(e).splitlines()[0][:88]}"))
            return None
        RESULTS.append((tid, what, "FAIL",
                        f"wrong exception {type(e).__name__} "
                        f"(wanted {expect.__name__}): "
                        f"{str(e).splitlines()[0][:70]}"))
        return None
    if expect is None:
        RESULTS.append((tid, what, "PASS", note or "accepted"))
        return val
    RESULTS.append((tid, what, "FAIL",
                    f"DID NOT REFUSE -- returned {str(val)[:80]!r}. A check "
                    "that cannot fail is the F-038 shape."))
    return val


# ===================================================================== axis
def t_axis():
    good_px = [2911.5, 2760.5, 2610.5, 2459.5, 2308.5]
    good_v = [1e3, 1e4, 1e5, 1e6, 1e7]

    cal = check("C-1", "axis.calibrate accepts real QL Fig 3c ticks",
                lambda: axis.calibrate(good_px, good_v), None)
    if cal:
        # The published floor of that axis is ~500 RLU at the panel box bottom
        # (y~2959). An independent landmark the calibration was not fitted to.
        y500 = cal.to_pixel(500.0)
        ok = abs(y500 - 2959) < 5
        RESULTS.append(("C-2", "axis: 500 RLU lands on the panel box bottom",
                        "PASS" if ok else "FAIL",
                        f"y(500 RLU) = {y500:.1f}, box bottom 2959, "
                        f"delta {y500 - 2959:+.1f} px"))

    # M-1  error class 20: glyph centroids instead of tick marks. Label
    # centroids drift with the text, so their spacing is not uniform.
    check("M-1", "axis: glyph-centroid spacing (class 20)",
          lambda: axis.calibrate([2911.5, 2759.0, 2604.0, 2455.0, 2308.5],
                                 good_v),
          CalibrationError)

    # M-2  error class 27 / F-006: the fit input IS the validation target.
    check("M-2", "axis: fit input is the validation target (class 27, F-006)",
          lambda: axis.calibrate(good_px, good_v, forbid_target=good_px),
          ValidationTargetError)

    # M-3  F-019: a colorbar anchored on the two outer edges of its own border.
    # Two points define a line through anything, so the border cannot be
    # distinguished from the ticks.
    check("M-3", "axis: 2-point border anchoring (F-019)",
          lambda: axis.calibrate([100.0, 900.0], [1e2, 1e6]),
          CalibrationError)

    # M-4  F-033: an origin maintained separately from the scale. There is no
    # API for it, so the uniform-step form that would allow it must refuse.
    check("M-4", "axis: scale without a fitted origin (F-033)",
          lambda: axis.calibrate(good_px, 1.0), CalibrationError)

    # M-5  tick VALUES unequally spaced -- a single scale cannot describe it.
    check("M-5", "axis: non-uniform tick values",
          lambda: axis.calibrate(good_px, [1e3, 1e4, 1e5, 1e6, 1e8]),
          CalibrationError)

    # M-27 / C-23: a SINGLE bad tick is invisible to the max-residual test.
    # Measured 2026-11-10 on a clean 5-tick 100px ladder: moving one END tick
    # by +6 px gives a max residual of 2.40 px and PASSES a 2.5 px tolerance,
    # because that tick has leverage 0.2 and four fifths of its error is
    # absorbed into the fit. Leave-one-out has no such masking.
    bent = [500.0 + 6.0, 400.0, 300.0, 200.0, 100.0]
    bent_vals = [1e1, 1e2, 1e3, 1e4, 1e5]
    RESULTS.append(("C-23", "axis: a +6px single-tick error DOES pass the "
                            "max-residual test (the weakness)", "PASS"
                    if axis.calibrate(bent, bent_vals).residual_px <= 2.5
                    else "FAIL",
                    f"max residual "
                    f"{axis.calibrate(bent, bent_vals).residual_px:.2f} px "
                    f"<= 2.5 -- documented, which is why LOO exists"))
    check("M-27", "axis: leave-one-out catches the single bad tick",
          lambda: axis.calibrate(bent, bent_vals,
                                 loo_tol_px=axis.LOO_TOL_PX),
          CalibrationError)
    RESULTS.append(("C-24", "axis: leave-one-out accepts a clean ladder",
                    "PASS" if max(axis.loo_residuals(
                        [500., 400., 300., 200., 100.])) < 1e-9 else "FAIL",
                    "all 5 LOO residuals 0.000 px"))

    check("C-3", "axis: linear axis accepted",
          lambda: axis.calibrate([100., 200., 300., 400.], [0., 25., 50., 75.],
                                 kind="linear"), None)


# ===================================================================== bars
def _panel():
    """A bar with a full-width top rule at row 60, plus a replicate dot
    cluster at rows 40-48 in the central strip. The F-018/F-020 geometry."""
    ink = np.zeros((220, 120), bool)
    ink[60:150, 10:40] = True
    ink[40:48, 22:28] = True
    return ink


def t_bars():
    ink = _panel()
    r = check("C-4", "bars.bar_top reads the rule, not the dots",
              lambda: bars.bar_top(ink, 10, 40, 0, 160), None)
    if r:
        RESULTS.append(("C-5", "bars: dot cluster at row 40 was skipped",
                        "PASS" if r["row"] == 60 else "FAIL",
                        f"row={r['row']} (true top 60, dots at 40); "
                        f"flank delta {r['flank_delta_px']} px"))

    # M-6  F-018/F-020: the dots alone must not become a bar top.
    dots = np.zeros((220, 120), bool)
    dots[40:48, 22:28] = True
    check("M-6", "bars: replicate dots only (class 32, F-018/F-020)",
          lambda: bars.bar_top(dots, 10, 40, 0, 160), DetectionError)

    # M-7  the inner-trimmed window that let F-018 through. Three overlapping
    # dots can exceed 85% of a trimmed window; the full-width rule rejects them.
    wide = np.zeros((220, 120), bool)
    wide[40:50, 20:31] = True          # ~37% of the full 30px bar
    check("M-7", "bars: wide dot cluster, trimmed-window shape (F-018)",
          lambda: bars.bar_top(wide, 10, 40, 0, 160), DetectionError)

    # M-8  class 32 via flank disagreement: a sloped edge, not a rule.
    slope = np.zeros((220, 120), bool)
    for i, x in enumerate(range(10, 40)):
        slope[60 + i:150, x] = True
    check("M-8", "bars: sloped edge fails the flank test",
          lambda: bars.bar_top(slope, 10, 40, 0, 160), DetectionError)

    # M-9  the first scanned row is already the axis rule.
    base = np.zeros((220, 120), bool)
    base[150:153, 0:120] = True
    check("M-9", "bars: baseline read as a zero-height bar",
          lambda: bars.bar_top(base, 10, 40, 140, 200, baseline_row=151),
          DetectionError)

    # M-10 too narrow for a flank test at all.
    check("M-10", "bars: bar too narrow to flank-test",
          lambda: bars.bar_top(ink, 10, 14, 0, 160), DetectionError)


# ================================================================== heatmap
def t_heatmap():
    sw = np.array([[255, 255, 255], [255, 200, 200], [200, 0, 0]], float)
    cells = np.array([[254, 255, 255], [255, 201, 199], [200, 1, 0]], float)
    m = check("C-6", "heatmap.match accepts clean swatches",
              lambda: heatmap.match(cells, sw), None)
    if m:
        RESULTS.append(("C-7", "heatmap: max distance recorded", "PASS",
                        f"max nearest-swatch distance "
                        f"{m['max_distance']:.3f} (tol {heatmap.MAX_DIST})"))

    # M-11 the geometry diagnostic CLAUDE.md section 5 says has caught an
    # off-by-one header row twice.
    check("M-11", "heatmap: cell colour in no swatch (geometry wrong)",
          lambda: heatmap.match(np.array([[128, 128, 128]], float), sw),
          DetectionError)

    # M-12 QL E-1: argmin resolving an exact tie by LUT row order.
    dup = np.array([[255, 255, 255], [200, 0, 0], [200, 0, 0]], float)
    check("M-12", "heatmap: exact tie resolved by row order (QL E-1)",
          lambda: heatmap.match(np.array([[200, 0, 0]], float), dup),
          DetectionError)

    # F-032: a band's stored value must be its midpoint, not an edge.
    inv = heatmap.invert([1], [(None, 1e4), (1e4, 1e5), (1e5, None)])
    mid = inv[0]["value"]
    ok = abs(mid - 31622.776) < 1.0
    RESULTS.append(("C-8", "heatmap.invert returns the band MIDPOINT (F-032)",
                    "PASS" if ok else "FAIL",
                    f"band [1e4,1e5] -> {mid:.1f} (geometric midpoint "
                    f"31622.8), bounds carried: lo={inv[0]['lo']} "
                    f"hi={inv[0]['hi']}"))

    # F-040: floor and ceiling flagged by COLOUR IDENTITY, symmetrically.
    both = heatmap.invert([0, 2], [(None, 1e4), (1e4, 1e5), (1e5, None)],
                          floor_index=0, ceiling_index=2)
    ok = both[0]["bound"] == "upper" and both[1]["bound"] == "lower"
    RESULTS.append(("C-9", "heatmap.invert flags floor and ceiling "
                           "symmetrically (F-040)",
                    "PASS" if ok else "FAIL",
                    f"floor -> {both[0]['bound']!r}, "
                    f"ceiling -> {both[1]['bound']!r}"))

    # C-10 (a control, not a mutation): a middle band must NOT be flagged as a
    # bound. F-040's defect was a value threshold straddling the floor band and
    # a band measurably not white, so "does not over-flag" needs a control.
    midband = heatmap.invert([1], [(None, 1e4), (1e4, 1e5), (1e5, None)],
                             floor_index=0)
    RESULTS.append(("C-10", "heatmap.invert does not flag a middle band "
                            "(F-040)",
                    "PASS" if midband[0]["bound"] == "" else "FAIL",
                    f"middle band bound = {midband[0]['bound']!r}"))

    check("M-14", "heatmap: LUT and band list disagree",
          lambda: heatmap.invert([5], [(None, 1e4)]), DetectionError)


# =================================================================== chains
def t_chains():
    zig = [(0, 0, 1, 1), (1, 1, 2, 0), (2, 0, 3, 1), (3, 1, 4, 0),
           (4, 0, 5, 1), (5, 1, 6, 0)]
    n = chains.vector_chain_length(zig, tol=0.1)
    RESULTS.append(("C-11", "chains.vector_chain_length counts a zig-zag",
                    "PASS" if n == 7 else "FAIL",
                    f"6 bonds -> {n} vertices (want 7)"))

    dbl = [(0, 0, 10, 0), (1, 0.15, 9, 0.15), (10, 0, 20, 0)]
    n2 = chains.vector_chain_length(dbl, tol=0.5, double_offset_pt=0.5)
    RESULTS.append(("C-12", "chains: double bond's inner line dropped",
                    "PASS" if n2 == 3 else "FAIL",
                    f"2 bonds + 1 inner line -> {n2} vertices (want 3)"))

    check("C-13", "chains.count_carbons accepts two agreeing methods",
          lambda: chains.count_carbons(vector=18, raster=18, label="OF-02"),
          None, "vector=18, raster=18 -> 18")

    # M-15 F-031 / class 48: OF-2 read as C17 where the reagent is C18.
    check("M-15", "chains: vector and raster disagree by one (F-031)",
          lambda: chains.count_carbons(vector=17, raster=18, label="OF-2"),
          DetectionError)

    # M-16 counting by eye, i.e. supplying nothing measured.
    check("M-16", "chains: no measured count supplied (class 48)",
          lambda: chains.count_carbons(label="by eye"), DetectionError)

    tok = chains.check_name_token(24, 4, "KS_2024 C24")
    RESULTS.append(("C-14", "chains.check_name_token reports, never sets",
                    "PASS" if tok["matches"] is False and tok["note"]
                    else "FAIL",
                    "C24 vs measured 4 flagged for adjudication, not "
                    "corrected"))


# ================================================================= assemble
def t_assemble():
    good = check("C-15", "assemble.esterify builds a 3-beta sterol ester",
                 lambda: assemble.esterify("CCCCC(=O)O", CHOL), None)
    if good:
        fp = assemble.stereo_fingerprint(good)
        want = assemble.stereo_fingerprint(CHOL)
        RESULTS.append(("C-16", "assemble: all 8 sterol centres preserved",
                        "PASS" if fp == want else "FAIL",
                        f"product {fp} vs cholesterol {want}"))

    # ---- P-1: the F-013 REPLAY, not a stand-in -----------------------
    A = Chem.MolFromSmiles("CCCCC(=O)O")
    S = Chem.MolFromSmiles(CHOL)
    cidx, _dbl, aoh = A.GetSubstructMatches(
        Chem.MolFromSmarts("[CX3](=O)[OX2H1]"))[0]
    oh = S.GetSubstructMatches(Chem.MolFromSmarts("[OX2H1]"))[0][0]
    c3 = [n.GetIdx() for n in S.GetAtomWithIdx(oh).GetNeighbors()][0]
    rw = Chem.RWMol(Chem.CombineMols(A, S))
    off = A.GetNumAtoms()
    rw.AddBond(aoh, c3 + off, Chem.BondType.SINGLE)
    rw.RemoveAtom(oh + off)          # <- deletes an atom ADJACENT to C3
    bad_mol = rw.GetMol()
    Chem.SanitizeMol(bad_mol)
    bad = Chem.MolToSmiles(bad_mol)

    from rdkit.Chem import rdMolDescriptors as rdd
    same_formula = (rdd.CalcMolFormula(Chem.MolFromSmiles(bad))
                    == rdd.CalcMolFormula(Chem.MolFromSmiles(good)))
    skel_bad = Chem.MolToInchiKey(Chem.MolFromSmiles(bad)).split("-")[0]
    skel_good = Chem.MolToInchiKey(Chem.MolFromSmiles(good)).split("-")[0]
    RESULTS.append(("P-1a", "F-013 replay is genuinely invisible downstream",
                    "PASS" if (same_formula and skel_bad == skel_good
                               and bad != good) else "FAIL",
                    f"formula identical: {same_formula}; InChIKey skeleton "
                    f"identical: {skel_bad == skel_good} ({skel_bad}); "
                    f"SMILES differ: {bad != good}"))
    # C-25: the blindness, as a NUMBER rather than a sentence. A sibling
    # curation measured the same thing on a real library on 2026-11-10:
    # seeding epicholesterol changed 15 of 18 FULL InChIKeys and 0 of 18
    # SKELETON InChIKeys. This is why C11 missed F-013 twice -- it compared
    # skeleton blocks -- and why any sterol check must use the full key.
    sk_same = (Chem.MolToInchiKey(Chem.MolFromSmiles(bad)).split("-")[0]
               == Chem.MolToInchiKey(Chem.MolFromSmiles(good)).split("-")[0])
    full_diff = (Chem.MolToInchiKey(Chem.MolFromSmiles(bad))
                 != Chem.MolToInchiKey(Chem.MolFromSmiles(good)))
    RESULTS.append(("C-25", "F-013: skeleton InChIKey is BLIND, full key is not",
                    "PASS" if (sk_same and full_diff) else "FAIL",
                    f"skeleton identical: {sk_same}; full key differs: "
                    f"{full_diff}; exact mass identical to 4 dp "
                    f"({rdd.CalcExactMolWt(Chem.MolFromSmiles(good)):.4f})"))

    check("P-1b", "F-013 replay caught by check_stereo_preserved",
          lambda: assemble.check_stereo_preserved(bad, ["CCCCC(=O)O", CHOL],
                                                  label="F-013 replay"),
          StructureError)

    # C-22 (control): a bond that CREATES a stereocentre must be ACCEPTED.
    # The original rule compared the whole fingerprint for equality and refused
    # these, firing on 40 legitimate assemblies in andersondg_41203968.
    frag, tail = "C1C[NH:1]C2CCCC2C1", "[CH2:2]CCCCC"
    try:
        prodc = assemble.attach(frag, tail)
        r = assemble.check_stereo_preserved(prodc, [frag.replace("[NH:1]", "N"),
                                                    "CCCCCC"], label="gains")
        RESULTS.append(("C-22", "assemble: creating a stereocentre is allowed",
                        "PASS", f"accepted; centres_created="
                                f"{r.get('centres_created')}"))
    except Exception as e:                                   # noqa: BLE001
        RESULTS.append(("C-22", "assemble: creating a stereocentre is allowed",
                        "FAIL", f"refused a legitimate build: "
                                f"{str(e).splitlines()[0][:70]}"))

    # M-26 losing an assigned centre must still be refused.
    check("M-26", "assemble: an assigned centre LOST",
          lambda: assemble.check_stereo_preserved(
              "CCCCC(=O)OC1CCC2(C)C(=CCC3C2CCC2(C)C(C(C)CCCC(C)C)CCC32)C1",
              ["CCCCC(=O)O", CHOL], label="stripped stereo"),
          StructureError)

    # C-26 / M-28: a CIP RELABEL is not an inversion, and counts cannot tell
    # them apart. Measured: (S)-glycidol's C2 is S, its decanoate ester is R,
    # with no atom moved -- the free OH leaves the C1/C3 comparison tied and
    # the ester breaks the tie. The counts check REFUSED that correct build
    # (reported from chamolstad_41691524, 2026-11-10). Parity is
    # numbering- and substituent-independent and gets both cases right.
    gly = "OC[C@H]1CO1"
    ester = "CCCCCCCCCC(=O)OC[C@H]1CO1"
    ok_relabel, det = assemble.parity_preserved(ester, gly)
    RESULTS.append(("C-26", "assemble.parity_preserved accepts a CIP RELABEL",
                    "PASS" if ok_relabel else "FAIL",
                    f"(S)-glycidol -> its ester: CIP S->R with no spatial "
                    f"change; parity ok={ok_relabel}, "
                    f"centres={det.get('centres_compared')}"))
    ok_epimer, det2 = assemble.parity_preserved(bad, CHOL)
    RESULTS.append(("M-28", "assemble.parity_preserved still REJECTS the "
                            "F-013 epimer", "PASS" if not ok_epimer else "FAIL",
                    f"parity ok={ok_epimer} (must be False); "
                    f"flipped={det2.get('flipped')}"))

    # C-27 / M-29: MAX_DIST is calibrated for matching cells against swatches
    # from THE SAME raster. Across panels it is dangerous. Measured on
    # peerd_31999380: S3 "Lipid 13" (127,4,127) vs S4 "Lipid 12" (128,0,129) are
    # 4.583 apart, BELOW MAX_DIST 5.0 -- so a cross-panel match assigns one
    # lipid's values to another and never trips the guard.
    pa = {"Lipid 12": (0, 112, 192), "Lipid 13": (127, 4, 127)}
    pb = {"Lipid 12": (128, 0, 129), "Lipid 13": (0, 110, 190)}
    coll = heatmap.cross_panel_collisions(pa, pb)
    RESULTS.append(("M-29", "heatmap.cross_panel_collisions FINDS the real "
                            "cross-panel swap", "PASS" if len(coll) == 2 else "FAIL",
                    f"{len(coll)} different-label pairs below MAX_DIST "
                    f"{heatmap.MAX_DIST}; closest {coll[0]['distance'] if coll else None}"))
    clean = heatmap.cross_panel_collisions({"L1": (255, 0, 0), "L2": (0, 255, 0)},
                                           {"L1": (254, 1, 1), "L2": (1, 254, 1)})
    RESULTS.append(("C-27", "heatmap.cross_panel_collisions stays quiet on "
                            "well-separated panels", "PASS" if not clean else "FAIL",
                    f"{len(clean)} collisions (must be 0)"))

    # C-28 / M-30: dpi-aware tolerance. A flattened SI raster has its own
    # resolution; rendering the page at 600 dpi upsamples rather than adding
    # information. Measured on khanof_36652629: panels at 188-221 dpi, where
    # this module's 1.5 px is 0.55 of a REAL pixel. native_tol expresses the
    # same physical tolerance in the panel's units -- and must REFUSE a guessed
    # dpi, because a tolerance rescaled by a guess looks derived.
    t204 = axis.native_tol(axis.LOO_TOL_PX, 204)
    ok_scale = abs(t204 - 1.5 * 600 / 204) < 1e-9 and \
        axis.native_tol(axis.LOO_TOL_PX, 600) == axis.LOO_TOL_PX
    RESULTS.append(("C-28", "axis.native_tol rescales by dpi and is identity at "
                            "600", "PASS" if ok_scale else "FAIL",
                    f"204 dpi -> {t204:.3f} px; 600 dpi -> "
                    f"{axis.native_tol(axis.LOO_TOL_PX, 600)} px"))
    refused = False
    for bad in (0, None, -1):
        try:
            axis.native_tol(1.5, bad)
        except Exception:
            refused = True
        else:
            refused = False
            break
    RESULTS.append(("M-30", "axis.native_tol REFUSES a missing or absurd dpi",
                    "PASS" if refused else "FAIL",
                    "a tolerance rescaled by a guessed dpi looks derived"))

    # M-17 the Siegwart radical carbon.
    check("M-17", "assemble: radical atom in the product",
          lambda: assemble.check_no_radicals("CCC[CH2]"), StructureError)

    # M-18 bonding to an unmapped fragment.
    check("M-18", "assemble: attachment atom not atom-mapped",
          lambda: assemble.attach("CCC[N:1]", "CCC"), StructureError)

    check("C-17", "assemble.attach makes a secondary amine, not a radical",
          lambda: assemble.attach("CCCC[NH:1]", "[CH2:2]CCO"), None,
          "CCCC[NH:1] + [CH2:2]CCO -> CCCCNCCCO")

    # M-19 two candidate acid groups: which reacts is not ours to choose.
    check("M-19", "assemble.esterify: two free acids",
          lambda: assemble.esterify("OC(=O)CCCC(=O)O", CHOL), StructureError)


# ================================================================= reagents
def t_reagents():
    RESULTS.append(("C-18", "reagents: MC3 and KC2 skeletons differ",
                    "PASS" if (reagents.skeleton(
                        reagents.registered("Dlin-MC3-DMA"))
                        != reagents.skeleton(reagents.registered("KC2")))
                    else "FAIL",
                    "same formula and exact mass 641.6111 Da, different "
                    "skeleton InChIKey -- the pair that motivated C11"))

    check("M-20", "reagents.smi: unregistered name with no corpus",
          lambda: reagents.smi("Nonesuch-PEG"), AmbiguityError)

    # M-21 the class-46 cross-check on a synthetic OF-2 / OF-02 collision,
    # placed in a GATED column so it must be gated rather than reported.
    of18 = ("CCCCCCCCCCCCCCCCCCOC(=O)CCN(CCC(=O)OCCCCCCCCCCCCCCCCCC)"
            "CCC(=O)OCCCCCCCCCCCCCCCCCC")
    of17 = of18.replace("CCCCCCCCCCCCCCCCCC", "CCCCCCCCCCCCCCCCC")
    rows = [{"Experiment_ID": "LL_2025", "HL_name": "OF-02",
             "HL_SMILES": of18},
            {"Experiment_ID": "ZT_2026", "HL_name": "OF-02",
             "HL_SMILES": of17}]
    res = reagents.cross_check(rows)
    gated = len(res["gated"])
    RESULTS.append(("M-21", "reagents.cross_check gates an OF-2 style "
                            "collision (class 46, F-031)",
                    "PASS" if gated == 1 else "FAIL",
                    f"{gated} gated entry(ies); "
                    + (reagents.explain(res['gated'][0]).splitlines()[0]
                       if gated else "nothing gated")))

    # C-19 a generic combinatorial IL label must be REPORTED, never gated.
    gen = [{"Experiment_ID": "JC_2023", "IL_name": "A1_B1_C1",
            "IL_SMILES": "CCCCN(C)C"},
           {"Experiment_ID": "YX_2024", "IL_name": "A1_B1_C1",
            "IL_SMILES": "CCCCCN(C)C"}]
    r2 = reagents.cross_check(gen)
    RESULTS.append(("C-19", "reagents.cross_check does not gate a generic "
                            "IL label",
                    "PASS" if (not r2["gated"] and len(r2["reported"]) == 1)
                    else "FAIL",
                    f"gated={len(r2['gated'])} reported={len(r2['reported'])} "
                    "-- 288 IL names legitimately collide (CLAUDE.md s3)"))


# ===================================================================== emit
def t_emit():
    recs = {"L1": {"site": "ASSUMED: no ring N and the paper draws no "
                           "structure"},
            "L2": {"site": "read from Fig 2a"}}
    check("C-20", "emit.emit_rows derives confidence per item",
          lambda: emit.emit_rows(
              ["L1", "L2"], recs,
              confidence=lambda i, r: ("low" if r["site"].startswith("ASSUMED")
                                       else "high"),
              note=lambda i, r: (r["site"]
                                 if r["site"].startswith("ASSUMED") else ""),
              label="control"),
          None, "L1 -> low (its record says ASSUMED), L2 -> high")

    # M-22 F-011 / F-025 / F-046: a hard-coded confidence literal.
    check("M-22", "emit: literal confidence (class 29; F-011/F-025/F-046)",
          lambda: emit.emit_rows(["L1"], recs, confidence="high",
                                 note=lambda i, r: "x"), EmitError)

    # M-23 F-024's confidence_without_a_note, made unwritable.
    check("M-23", "emit: caveat with an empty note (F-024)",
          lambda: emit.emit_rows(["L1"], recs,
                                 confidence=lambda i, r: "low",
                                 note=lambda i, r: ""), EmitError)

    # M-24 a value outside the schema vocabulary.
    check("M-24", "emit: value outside schema_quality.json's vocabulary",
          lambda: emit.emit_rows(["L1"], recs,
                                 confidence=lambda i, r: "quite good",
                                 note=lambda i, r: "n"), EmitError)

    # M-25 a flag with no enumeration record behind it.
    check("M-25", "emit: item with no enumeration record",
          lambda: emit.emit_rows(["L9"], recs,
                                 confidence=lambda i, r: "high",
                                 note=lambda i, r: ""), EmitError)

    # C-21 the live-corpus cross-validation of a number that had no
    # definition in code and moved four times (F-024).
    rows = []
    for f in ("data_vitro.csv", "data_vivo.csv"):
        p = os.path.join(HERE, "..", "data", f)
        if os.path.exists(p):
            with open(p) as fh:
                rows += list(csv.DictReader(fh))
    if rows:
        n = emit.caveat_count(rows)
        RESULTS.append(("C-21", "emit.caveat_count reproduces the gated "
                                "rows_with_a_caveat (F-024)",
                        "PASS" if n == 1424 else "FAIL",
                        f"independent re-implementation gives {n}; "
                        "CLAUDE.md s6 and audit.py C13 state 1,424"))


def main():
    for fn in (t_axis, t_bars, t_heatmap, t_chains, t_assemble, t_reagents,
               t_emit):
        try:
            fn()
        except Exception as e:                   # noqa: BLE001
            RESULTS.append((fn.__name__, "test group crashed", "FAIL",
                            f"{type(e).__name__}: {e}"))

    width = max(len(t) for t, *_ in RESULTS)
    npass = sum(1 for *_, s, _ in ((a, b, c, d) for a, b, c, d in RESULTS)
                if s == "PASS")
    print("=" * 78)
    print("curation/lib seeded-mutation suite")
    print("=" * 78)
    for tid, what, status, detail in RESULTS:
        mark = "ok  " if status == "PASS" else "FAIL"
        print(f"[{mark}] {tid:<{width}}  {what}")
        if detail:
            print(f"         {detail}")
    nfail = len(RESULTS) - npass
    print("-" * 78)
    print(f"{npass} passed, {nfail} failed, {len(RESULTS)} total")
    mut = sum(1 for t, *_ in RESULTS if t.startswith("M-"))
    rep = sum(1 for t, *_ in RESULTS if t.startswith("P-"))
    ctl = sum(1 for t, *_ in RESULTS if t.startswith("C-"))
    print(f"{mut} seeded mutations, {rep} historical replays, {ctl} controls")
    if nfail:
        print("\nFAILED -- an invariant in curation/lib is not enforced, or a "
              "module refused a correct input.")
        return 1
    print("\nAll invariants fire on their seeded defect and accept their "
          "controls.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
