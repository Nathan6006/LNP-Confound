"""Figure S1 - ESI-MS base-peak m/z by AXIS GEOMETRY, not by reading glyphs.

WHY THIS FILE EXISTS, AND WHY IT DID NOT BEFORE. The first pass measured the
base-peak labels at 5-7 px tall, found them undecidable between 784.73 / 794.75
/ 764.73, and refused - correctly, because with no hypothesis a 4-7 Da/px axis
cannot pin a formula.

That changed when head 9322 was recovered from the corpus (YR_2024, 22 rows,
IL_head_SMILES `Cc1nccn1CCCN`). With a known head the question stops being
"what is this number" and becomes "which of three candidate numbers is it",
and those candidates are 20 and 30 Da apart. A ~6.8 Da/px axis separates them.
A measurement that cannot pin a value can still decide between hypotheses.

The tick VALUES are not read either. They are fixed by a geometric signature:
ten evenly spaced ticks whose labels are three digits wide for the first four
and four digits wide for the last six. The only ladder that does that on an
instrument whose own header states a [150-2000] scan range is 200,400,...,2000.
"""
import json, os, sys
import numpy as np
from PIL import Image
from rdkit import Chem, RDLogger
from rdkit.Chem.Descriptors import ExactMolWt
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "lib")))
import axis                                               # noqa: E402

RASTER = os.path.join(HERE, "source", "figS1_page4_xref95.png")
PROTON = 1.0072765

# head 9322, from the corpus via the guarded lookup (see stage01_blocks.py).
HEAD_9322 = "Cc1nccn1CCCN"

# the six Figure S1 panels, in drawn order
PANELS = [(0, 0, "9322-O17C"), (0, 1, "9322-O17O"), (0, 2, "9322-O17S"),
          (1, 0, "9322-O16B"), (1, 1, "9322-N16B"), (1, 2, "9322-EC18")]

# candidate tail structures. O16B / N16B are HYPOTHESES (D-3); EC18 is from
# reference [3]; O17C has no candidate at all (D-4).
TAILS = {
    "O17O":  ("C=CC(=O)OCCOCCCCCCCCCCCCCC",    "SI 1.1, built"),
    "O17S":  ("C=CC(=O)OCCSCCCCCCCCCCCCCC",    "SI 1.2, built"),
    "O16B":  ("C=CC(=O)OCCSSCCCCCCCCCCCC",     "HYPOTHESIS: O16B = O16S-S (D-3)"),
    "N16B":  ("C=CC(=O)NCCSSCCCCCCCCCCCC",     "HYPOTHESIS: N16B = N16S-S (D-3)"),
    "EC18":  ("CCCCCCCCCCCCCCCCC1CO1",         "HYPOTHESIS: 1,2-epoxyoctadecane, ref [3]"),
    "O17C":  (None,                            "NO CANDIDATE - D-4 is open"),
}
# my three candidate glyph reads of the 9322-O17O base-peak label, first pass
FIRST_PASS_CANDIDATES = [784.73, 794.75, 764.73]


def predicted_MH(tail_smi):
    """head + 2 x tail, no atom loss (Michael addition / epoxide opening)."""
    h = Chem.MolFromSmiles(HEAD_9322)
    t = Chem.MolFromSmiles(tail_smi)
    m = ExactMolWt(h) + 2 * ExactMolWt(t)
    return m + PROTON


def panel_box(a, r, c):
    H, W, _ = a.shape
    return a[r * H // 2:(r + 1) * H // 2, c * W // 3:(c + 1) * W // 3]


def axis_row(dark):
    """The x-axis rule: the row of the LOWER HALF with the most dark pixels.

    Not the longest contiguous run. Measured 2026-11-10: on these Xcalibur
    panels the rule is broken by the tick notches, so a longest-run search
    picked a panel-title underline (20 px) in one panel and a fragment of the
    rule (74 px) in another, and the tick labels were then read off the wrong
    rows. A dark-pixel COUNT is insensitive to the notches.
    """
    h = dark.shape[0]
    counts = dark[h // 2:].sum(axis=1)
    i = int(np.argmax(counts)) + h // 2
    return i, int(counts.max())


def label_boxes(dark, y0, y1, gap=3, min_w=5):
    """Bounding boxes of the tick-label glyph groups in rows [y0, y1)."""
    prof = dark[y0:y1].sum(axis=0)
    on = prof > 0
    out, start, hole = [], None, 0
    for j, v in enumerate(on):
        if v:
            if start is None:
                start = j
            hole = 0
        elif start is not None:
            hole += 1
            if hole > gap:
                if j - hole - start >= min_w:
                    out.append((start, j - hole))
                start, hole = None, 0
    if start is not None:
        out.append((start, len(on)))
    return [b for b in out if b[1] - b[0] >= min_w]


def main():
    a = np.array(Image.open(RASTER).convert("RGB")).astype(int)
    out = {"_what": "base-peak m/z by axis geometry; tick VALUES fixed by a label-width "
                    "signature, not by reading them", "panels": {}}
    for r, c, name in PANELS:
        sub = panel_box(a, r, c)
        g = sub.mean(2)
        dark = g < 140
        arow, alen = axis_row(dark)
        # The glyph gap inside a four-digit label and the gap between two
        # labels differ per panel (the rasters are at slightly different
        # scales). Search the gap, but ACCEPT only on criteria that are
        # independent of the answer: exactly 10 boxes, spacing uniform to
        # <= 4 % of the pitch, and the 3-digit/4-digit width signature. If no
        # gap satisfies all three the panel is refused -- the search cannot
        # manufacture a pass.
        boxes, chosen_gap, tried = None, None, []
        for gp in range(2, 13):
            cand = label_boxes(dark, arow + 2, arow + 12, gap=gp)
            cen = [(aa + b) / 2.0 for aa, b in cand]
            ok10 = len(cand) == 10
            unif = None
            if ok10:
                d = np.diff(cen)
                unif = float(np.max(np.abs(d - d.mean())) / d.mean())
            wsig = None
            if ok10:
                ww = [b - aa for aa, b in cand]
                wsig = float(np.mean(ww[4:]) / np.mean(ww[:4]))
            tried.append(dict(gap=gp, n=len(cand),
                              spacing_nonuniformity=None if unif is None else round(unif, 4),
                              width_ratio=None if wsig is None else round(wsig, 3)))
            if ok10 and unif <= 0.04 and 1.15 <= wsig <= 1.55:
                boxes, chosen_gap = cand, gp
                break
        rec = dict(axis_row=int(arow), axis_len_px=int(alen),
                   gap_search=tried, gap_chosen=chosen_gap,
                   n_label_boxes=None if boxes is None else len(boxes))
        if boxes is None:
            rec["REFUSED"] = ("no glyph gap in 2..12 gives 10 tick labels that are "
                              "uniformly spaced to 4 % AND carry the 3-digit/4-digit "
                              "width signature")
            out["panels"][name] = rec
            continue
        rec["label_widths"] = [b - aa for aa, b in boxes]
        w = rec["label_widths"]
        narrow, wide = np.mean(w[:4]), np.mean(w[4:])
        rec["label_width_signature"] = dict(
            mean_first_four=round(float(narrow), 2), mean_last_six=round(float(wide), 2),
            ratio=round(float(wide / narrow), 3), expected_ratio_4_over_3=1.333,
            verdict=("three-digit labels 200-800 then four-digit 1000-2000"
                     if 1.15 <= wide / narrow <= 1.55 else "SIGNATURE FAILS"))
        centres = [(aa + b) / 2.0 for aa, b in boxes]
        vals = [200.0 * (k + 1) for k in range(10)]
        cal = axis.calibrate(centres, vals, kind="linear", axis="x",
                             residual_tol_px=3.0)
        rec["calibration"] = cal.qc()
        loo = axis.loo_residuals(centres)
        rec["loo_worst_px"] = round(max(loo), 3)
        rec["Da_per_px"] = round(1.0 / cal.px_per_unit, 3)

        # the base peak: the tallest ink column inside the plot area
        plot = dark[:arow - 2, :]
        col = plot.sum(axis=0)
        # ignore the y-axis furniture on the left and the panel title glyphs
        lo = int(round(cal.to_pixel(150.0)))
        col[:max(0, lo)] = 0
        jmax = int(np.argmax(col))
        win = slice(max(0, jmax - 2), jmax + 3)
        wts = col[win].astype(float)
        xs = np.arange(win.start, win.stop, dtype=float)
        centroid = float((xs * wts).sum() / wts.sum())
        rec["base_peak_col_px"] = round(centroid, 2)
        rec["base_peak_height_px"] = int(col[jmax])
        rec["measured_mz"] = round(float(cal.to_value(centroid)), 1)

        smi, prov = TAILS[name.split("-", 1)[1]]
        rec["tail_provenance"] = prov
        if smi:
            pred = predicted_MH(smi)
            rec["predicted_MplusH"] = round(pred, 3)
            rec["formula_neutral"] = CalcMolFormula(Chem.MolFromSmiles(
                Chem.MolToSmiles(Chem.MolFromSmiles(HEAD_9322))))
            rec["deviation_Da"] = round(rec["measured_mz"] - pred, 1)
            rec["deviation_px"] = round(rec["deviation_Da"] / rec["Da_per_px"], 2)
            # THE DISCRIMINATION: does this panel pick the 2-methyl head over the
            # plain-imidazole head 93? The two differ by exactly one CH2 (14.016).
            alt = pred - ExactMolWt(Chem.MolFromSmiles(HEAD_9322)) \
                       + ExactMolWt(Chem.MolFromSmiles("NCCCn1ccnc1"))
            rec["head_discrimination"] = dict(
                predicted_with_head_9322_2_methylimidazole=round(pred, 3),
                predicted_with_head_93_plain_imidazole=round(alt, 3),
                separation_Da=round(pred - alt, 3),
                separation_px=round((pred - alt) / rec["Da_per_px"], 2),
                deviation_from_9322_px=rec["deviation_px"],
                deviation_from_93_px=round((rec["measured_mz"] - alt) /
                                           rec["Da_per_px"], 2),
                picks=("9322 (2-methylimidazole)"
                       if abs(rec["measured_mz"] - pred) < abs(rec["measured_mz"] - alt)
                       else "93 (plain imidazole)"))
        else:
            rec["predicted_MplusH"] = None
            # no candidate tail: DERIVE the per-arm mass instead of asserting one
            hm = ExactMolWt(Chem.MolFromSmiles(HEAD_9322))
            per_arm = (rec["measured_mz"] - PROTON - hm) / 2.0
            rec["derived_tail_mass_Da"] = round(per_arm, 2)
            rec["derived_tail_uncertainty_Da"] = round(rec["Da_per_px"] / 2.0, 2)
            cands = {"heptadecyl acrylate C20H38O2 (plain C17 alkyl)":
                     ExactMolWt(Chem.MolFromSmiles("C=CC(=O)OCCCCCCCCCCCCCCCCC")),
                     "2-(hexadecyloxycarbonyloxy)ethyl acrylate (a carbonate, C22H40O5)":
                     ExactMolWt(Chem.MolFromSmiles(
                         "C=CC(=O)OCCOC(=O)OCCCCCCCCCCCCCCCC")),
                     "O17O (ether, for scale)":
                     ExactMolWt(Chem.MolFromSmiles("C=CC(=O)OCCOCCCCCCCCCCCCCC"))}
            rec["derived_tail_vs_candidates_Da"] = {
                k: round(per_arm - v, 2) for k, v in cands.items()}
        out["panels"][name] = rec

    # the decisive comparison: which first-pass candidate does the geometry pick?
    o17o = out["panels"].get("9322-O17O", {})
    if "measured_mz" in o17o:
        out["_first_pass_candidate_test"] = {
            "candidates_read_as_glyphs": FIRST_PASS_CANDIDATES,
            "predicted_from_the_recovered_head": o17o["predicted_MplusH"],
            "measured_by_geometry": o17o["measured_mz"],
            "Da_per_px": o17o["Da_per_px"],
            "distance_to_each_candidate_Da": {str(k): round(o17o["measured_mz"] - k, 1)
                                              for k in FIRST_PASS_CANDIDATES},
            "distance_to_each_candidate_px": {
                str(k): round((o17o["measured_mz"] - k) / o17o["Da_per_px"], 2)
                for k in FIRST_PASS_CANDIDATES},
        }
    with open(os.path.join(HERE, "09_figS1_mass_test.json"), "w") as f:
        json.dump(out, f, indent=1)
    for n, rec in out["panels"].items():
        if "REFUSED" in rec:
            print("%-12s REFUSED %s" % (n, rec["REFUSED"]))
        else:
            print("%-12s meas %7.1f  pred %s  dev %s Da (%s px)  [%.2f Da/px, loo %.2f]"
                  % (n, rec["measured_mz"], rec["predicted_MplusH"],
                     rec.get("deviation_Da"), rec.get("deviation_px"),
                     rec["Da_per_px"], rec["loo_worst_px"]))
    if "_first_pass_candidate_test" in out:
        print("candidate test:", json.dumps(out["_first_pass_candidate_test"]
                                            ["distance_to_each_candidate_Da"]))

if __name__ == "__main__":
    main()
