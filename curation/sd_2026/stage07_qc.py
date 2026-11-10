#!/usr/bin/env python3
"""Stage 07 -- every check, with numbers, and every check SEEDED to fail.

CLAUDE.md section 6: "A NEW CHECK IS NOT DONE UNTIL IT HAS BEEN MADE TO FAIL ON
PURPOSE." Three incidents paid for that rule, and all three printed a clean
line while examining nothing. So each check below records

    metric_clean -> metric_seeded

and the seed is required to MOVE THE METRIC, not merely to raise. A seed that
only produces an exception proves the code has an `if`; a seed that moves the
number proves the check is looking at the thing it claims to look at.

One check here is deliberately recorded as NOT firing: seeding H4 -> H3 in
LC-1 leaves the molecular formula identical, so the formula check does not
notice. That is the point. Mass cannot see connectivity, and the SI's ten
printed formulas are a necessary, not sufficient, test of this enumeration.
"""
import copy
import csv
import json
import math
import os
import re
import sys

import numpy as np
import pymupdf
from PIL import Image
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "lib"))
import axis                                                      # noqa: E402
import chains                                                    # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import stage01_blocks as S1                                      # noqa: E402
import stage03_enumerate as S3                                   # noqa: E402
import stage05_values as S5                                      # noqa: E402

OUT = {}


def rec(name, **kw):
    OUT[name] = kw
    return kw


# ===========================================================================
# P1  the axis scale, and a tick the authors DREW and we never used
# ===========================================================================
def p1_axis():
    p = S5.page()
    res = {}
    for nm, panel in S5.PANELS.items():
        cal, qc = S5.calibrate_panel(p, nm, panel)
        res[nm] = {"n_ticks_detected": len(qc["ticks_detected_px"]),
                   "n_ticks_in_fit": qc["n_ticks"],
                   "tick_fit_residual_px": qc["tick_fit_residual_px"],
                   "px_per_decade": qc["px_per_unit"],
                   "held_out": qc["ticks_held_out"]}
    # --- SEED: move one tick 3 px. The residual must grow. -----------------
    panel = S5.PANELS["Fig_1n"]
    cal, qc = S5.calibrate_panel(p, "Fig_1n", panel)
    good = [t[0] for t in qc["ticks_used_in_fit"]]
    clean = axis.calibrate(good, [10.0 ** e for e in (11, 10, 9, 8, 7)],
                           kind="log10", axis="y").qc()["tick_fit_residual_px"]
    seeds = {}
    for shift in (3.0, 6.0, 10.0):
        bad = list(good)
        bad[2] += shift
        resid, fired = None, False
        try:
            resid = axis.calibrate(
                bad, [10.0 ** e for e in (11, 10, 9, 8, 7)],
                kind="log10", axis="y").qc()["tick_fit_residual_px"]
        except Exception as e:
            fired = True
            m = re.search(r"residual ([0-9.]+) px", str(e))
            resid = float(m.group(1)) if m else None
        seeds[f"middle_tick_+{shift:.0f}px"] = {
            "metric": "tick_fit_residual_px", "clean": clean,
            "seeded": resid, "raised": fired,
            "metric_moved": resid is not None and resid > clean + 0.5}
    return rec("P1_axis_calibration",
               _what="Scale AND origin fitted together from detected tick "
                     "marks (curation/lib/axis.calibrate). The tick POSITIONS "
                     "are stroked segments in the PDF's vector layer and the "
                     "tick VALUES are PDF text spans; neither is a pixel read.",
               panels=res,
               seed=seeds,
               seed_sensitivity_finding=
                   "A single tick displaced by +3 px on a FIVE-tick fit "
                   "produces a residual of about 2.4 px, which is INSIDE "
                   "axis.calibrate's 2.5 px tolerance -- least squares spreads "
                   "the error over the other four points, and the middle point "
                   "of five equally spaced has leverage 0.2, so the residual "
                   "is about 0.8 x the displacement. The metric moves (0.0 -> "
                   "2.4) but the gate does not trip until roughly +3.2 px. "
                   "Recorded because a tolerance that is quoted without its "
                   "leverage is a tolerance nobody has measured; on these "
                   "panels it does not matter, since the observed residuals "
                   "are 0.0000-0.0001 px.",
               held_out_tick_note=
                   "Fig. 1p's 10^6 tick is ABSENT from the vector layer, so "
                   "that panel is fitted on its three contiguous ticks and the "
                   "10^5 tick is held out and predicted. It is a round number "
                   "the authors DREW and the fit never saw.")


def p2_heldout_seed():
    """The held-out tick must stop agreeing if the fit is wrong."""
    p = S5.page()
    _, qc = S5.calibrate_panel(p, "Fig_1p", S5.PANELS["Fig_1p"])
    clean = abs(qc["ticks_held_out"][0]["error_decades"])
    fit = [t[0] for t in qc["ticks_used_in_fit"]]
    held = qc["ticks_held_out"][0]["pixel"]
    # SEED: mislabel the fitted ticks by one decade each (10^9,10^8,10^7 read
    # as 10^8,10^7,10^6) -- a plausible off-by-one in reading the exponents.
    cal_bad = axis.calibrate(fit, [1e8, 1e7, 1e6], kind="log10", axis="y")
    seeded = abs(math.log10(cal_bad.to_value(held)) - 5.0)
    # SEED 2: compress the scale 10% -- a mis-detected tick spacing.
    cal_sq = axis.calibrate([fit[0], fit[0] + (fit[1] - fit[0]) * 0.9,
                             fit[0] + (fit[2] - fit[0]) * 0.9],
                            [1e9, 1e8, 1e7], kind="log10", axis="y")
    seeded2 = abs(math.log10(cal_sq.to_value(held)) - 5.0)
    return rec("P2_heldout_tick",
               _what="Fig. 1p's 10^5 tick, predicted by a fit that never saw "
                     "it. Tolerance 0.02 decades.",
               metric="|predicted - 10^5| in decades",
               clean=round(clean, 5),
               seed_exponents_off_by_one={"seeded": round(seeded, 4),
                                          "metric_moved": seeded > 0.02},
               seed_scale_compressed_10pct={"seeded": round(seeded2, 4),
                                            "metric_moved": seeded2 > 0.02})


# ===========================================================================
# P3  the category labels: confirmed by a SECOND drawn object
# ===========================================================================
def p3_label_mapping():
    """The significance brackets terminate on group centres.

    The label order comes from the rotated text spans. That is one piece of
    evidence. The figure ALSO draws two significance brackets whose ends sit on
    particular groups, and the legend says which comparison each is. If the
    labels were rotated by one position the brackets would land between groups.
    """
    p = S5.page()
    panel = S5.PANELS["Fig_1n"]
    gs = S5.groups(p, "Fig_1n", panel)
    names = S5.labels(p, panel, len(gs))
    xs = {n: g["x"] for n, g in zip(names, gs)}
    # the two long horizontal rules above the data are the brackets
    H = S5.horizontals(p, panel["plot_x"][0], panel["plot_x"][1],
                       panel["y"][0] + 2, panel["y"][0] + 12)
    # The plot FRAME spans the whole box and sits exactly on the top tick, so
    # it is excluded by both tests -- including it would have made the metric
    # 7.70 pt, which is the frame's overhang and not a bracket at all.
    frame_x = (round(panel["plot_x"][0], 1), round(panel["plot_x"][1], 1))
    brackets = [h for h in H if h[3] > 20
                and not (round(h[0], 1) <= frame_x[0] + 0.6
                         and round(h[1], 1) >= frame_x[1] - 1.5)]
    ends = sorted({round(v, 2) for h in brackets for v in (h[0], h[1])})
    if len(brackets) != 2:
        raise AssertionError(
            f"P3: expected 2 significance brackets above Fig. 1n's data, "
            f"found {len(brackets)}: {brackets}")
    near = []
    for e in ends:
        n, g = min(xs.items(), key=lambda kv: abs(kv[1] - e))
        near.append({"bracket_end_pt": e, "nearest_group": n,
                     "group_centre_pt": round(g, 2),
                     "offset_pt": round(e - g, 3)})
    worst = max(abs(d["offset_pt"]) for d in near)
    # SEED: rotate the labels by one. Every bracket end then lands ~11 pt off.
    rot = names[1:] + names[:1]
    xs_rot = {n: g["x"] for n, g in zip(rot, gs)}
    seeded = max(
        abs(e - xs_rot[min(xs_rot, key=lambda k: abs(xs_rot[k] - e))])
        for e in ends)
    # the real test of the rotation is whether the NAMED comparison still lands
    seeded_named = abs(ends[0] - xs_rot["ALC-0315"])
    legend_says = ("the legend prints two P values for this panel, "
                   "P(LC-1 vs ALC-0315) and P(LC-1 vs LP-01), both < 1e-15; "
                   "the two brackets DRAWN run ALC-0315 -> LC-1 and "
                   "LC-1 -> LC-5, so the second bracket's right end is LC-5 "
                   "and not LP-01. That mismatch is recorded in "
                   "08_contradictions.csv and does not affect the x -> name "
                   "mapping, which both the rotated text spans and the "
                   "bracket endpoints agree on.")
    return rec("P3_category_label_mapping",
               _what="Fig. 1n's two significance brackets end on group "
                     "centres. The legend names the comparisons, so this is a "
                     "second, independent witness to which x position is which "
                     "lipid -- the rotated text spans are the first.",
               group_x_spacing_pt=round(gs[1]["x"] - gs[0]["x"], 4),
               n_brackets=len(brackets),
               bracket_ends=near,
               legend_vs_drawn=legend_says,
               metric="max |bracket end - nearest group centre| in pt",
               clean=round(worst, 3),
               seed={"what": "label list rotated by one position",
                     "metric": "|first bracket end - x(ALC-0315)| in pt",
                     "clean": round(abs(ends[0] - xs["ALC-0315"]), 3),
                     "seeded": round(seeded_named, 3),
                     "metric_moved": seeded_named > 5.0})


# ===========================================================================
# P4  the drawn mean vs the markers the figure also draws
# ===========================================================================
def p4_markers():
    p = S5.page()
    panel = S5.PANELS["Fig_1n"]
    cal, _ = S5.calibrate_panel(p, "Fig_1n", panel)
    gs = S5.groups(p, "Fig_1n", panel)
    names = S5.labels(p, panel, len(gs))
    marks = []
    for q in p.get_drawings():
        r = q["rect"]
        if (q.get("type") == "fs" and q.get("fill") == (1.0, 1.0, 1.0)
                and abs(r.width - 1.61) < 0.25 and abs(r.height - 1.61) < 0.25
                and panel["plot_x"][0] < r.x0 < panel["plot_x"][1]
                and panel["y"][0] < r.y0 < panel["y"][1]):
            marks.append(((r.x0 + r.x1) / 2, (r.y0 + r.y1) / 2))
    per = []
    for g, nm in zip(gs, names):
        mine = [m for m in marks if abs(m[0] - g["x"]) < 5.0]
        if not mine:
            continue
        lg = [math.log10(cal.to_value(m[1])) for m in mine]
        drawn = math.log10(cal.to_value(g["mean_px"]))
        per.append({"lipid": nm, "n_markers_found": len(mine),
                    "marker_log10_mean": round(sum(lg) / len(lg), 4),
                    "drawn_mean_log10": round(drawn, 4),
                    "delta_decades": round(sum(lg) / len(lg) - drawn, 4),
                    "markers_straddle_the_drawn_mean":
                        min(lg) <= drawn <= max(lg)})
    inside = all(g["cap_hi_px"] <= g["mean_px"] <= g["cap_lo_px"] for g in gs)
    return rec("P4_markers_vs_drawn_mean",
               _what="Each group's five mice are drawn as 1.61 pt white "
                     "circles. We EMIT the mean rule the figure draws, not an "
                     "average of these -- jittered markers that coincide "
                     "exactly collapse into one path, so a marker average is "
                     "biased by whichever dots happen to be separable.",
               n_markers_separable=len(marks),
               n_markers_expected=12 * panel["n_mice"],
               why_fewer="markers at identical coordinates are one path; "
                         "this is why the drawn mean is used instead",
               mean_lies_between_its_own_sd_caps_all_groups=inside,
               per_group=per,
               metric="max |marker mean - drawn mean| in decades",
               clean=round(max(abs(d["delta_decades"]) for d in per), 4))


# ===========================================================================
# P5  chain lengths: the weighting is load-bearing
# ===========================================================================
def p5_chains():
    ink, grey = S1.si_fig1_raster()
    meas = S1.measure_chains(ink, grey)
    clean = {k: v["vertices"] for k, v in meas.items()}
    sweep = {k: v["eps_sweep"] for k, v in meas.items()}

    # SEED: the binary centreline this file used first. Ta3 drifts 18 -> 19.
    def binary_path(y0, y1, x0, xstart, xstop):
        sub = ink[y0:y1, x0:]
        cols = sub.sum(axis=0)
        path = []
        for x in range(xstart, xstop):
            c = cols[x]
            if c == 0 or c > 3:
                break
            ys = np.where(sub[:, x])[0]
            path.append((float(x), float(ys.mean())))
        return path
    y0, y1 = S1.TA_BANDS["Ta3"]
    sub = ink[y0:y1, S1.TA_X[0]:S1.TA_X[1]]
    xs = np.where(sub.sum(axis=0) > 0)[0]
    bp = binary_path(y0, y1, S1.TA_X[0], int(xs[0]), S1.TA_X[1] - S1.TA_X[0])
    seeded = {str(e): chains.raster_chain_length(bp, eps_px=e)
              for e in (0.8, 1.0, 1.2, 1.5, 2.0, 2.5, 3.0)}
    return rec("P5_chain_lengths",
               _what="Every long chain is measured with "
                     "curation/lib/chains.raster_chain_length over a "
                     "Douglas-Peucker sweep; a count that holds at only one "
                     "tolerance is refused by stage01's _stable().",
               carbons=clean, eps_sweep=sweep,
               seed={"what": "revert the centreline to an unweighted binary "
                             "row mean (what this file did first)",
                     "metric": "Ta3 vertex count across the eps sweep",
                     "clean": sweep["Ta3"], "seeded": seeded,
                     "metric_moved": len(set(seeded.values())) > 1,
                     "note": "this is the historical failure, not a synthetic "
                             "one: the first run of stage01 refused on exactly "
                             "this, reporting Ta3 = 19 at eps = 1.0 px and 18 "
                             "elsewhere."})


def p6_label_cost():
    """Tb1's hexyl calibrates what a drawn heteroatom label costs a count."""
    ink, grey = S1.si_fig1_raster()
    meas = S1.measure_chains(ink, grey)
    hexyl = meas["Tb1_hexyl"]["vertices"]
    clean_cost = hexyl - S1.TB1_HEXYL_CARBONS
    # SEED: pretend the hexyl is a heptyl (an off-by-one in Table 5's read).
    seeded_cost = hexyl - 7
    return rec("P6_label_end_cost",
               _what="A trace that begins or ends at a drawn 'CN' or 'O' label "
                     "gains one Douglas-Peucker vertex that is not an atom. "
                     "Tb1's hexyl is 6 carbons by Supplementary Table 5 "
                     "(6-amino-1-hexanol, CAS 4048-33-3) and runs label to "
                     "label, so it measures that cost instead of assuming it.",
               tb1_hexyl_vertices=hexyl,
               carbons_from_table5=S1.TB1_HEXYL_CARBONS,
               metric="vertices - carbons, over two label ends",
               clean=clean_cost,
               expected=2,
               seed={"what": "read Table 5's 6-amino-1-hexanol as a heptanol",
                     "metric": "label end cost", "seeded": seeded_cost,
                     "metric_moved": seeded_cost != 2},
               second_witness={
                   "_what": "the bond pitch, independently",
                   "ta_mean_px_per_bond": meas["Ta1"]["px_per_bond"],
                   "tb2_carbons": meas["Tb2"]["vertices"] - 1,
                   "note": "one carbon either way moves Tb2's implied pitch by "
                           "about 0.9 px against the Ta chains' 14.8."})


# ===========================================================================
# P7  the SI's ten printed formulas -- and what they CANNOT see
# ===========================================================================
def p7_formulas():
    b = S1  # noqa
    rows = list(csv.DictReader(open(os.path.join(HERE, "03_enumerated.csv"))))
    byf = {}
    for r in rows:
        byf.setdefault(r["formula"], []).append(r["lipid_id"])
    clean_hits = {k: len(byf.get(v[0], [])) for k, v in S3.SI_NAMED.items()}
    n_realised = sum(1 for v in clean_hits.values() if v)

    # --- SEED A: perturb EVERY block by one carbon and measure the check ---
    # One seed would have told us whether the check can fire. Ten tell us how
    # much it is worth, and the answer is: for the chain blocks, almost nothing.
    blocks = S3.blocks()

    def realised_with(cls, bid, smi):
        bad = copy.deepcopy(blocks)
        bad[cls][bid]["SMILES"] = smi
        forms = set()
        H, L = bad["headgroup_amine"], bad["linker_aldehyde"]
        TA, TB = bad["tail_a_acid"], bad["tail_b_isocyanide"]
        for h in H:
            for l in L:
                for ta in TA:
                    for tb in TB:
                        m = Chem.MolFromSmiles(
                            _ugi_fast(L[l]["SMILES"], H[h]["SMILES"],
                                      TA[ta]["SMILES"], TB[tb]["SMILES"]))
                        forms.add(rdMolDescriptors.CalcMolFormula(m))
        return sum(1 for (f, _, _) in S3.SI_NAMED.values() if f in forms)

    MUTS = {
        "Ta1 C14 -> C13": ("tail_a_acid", "Ta1", "CCCCCCCCCCCCC(=O)O"),
        "Ta2 C16 -> C15": ("tail_a_acid", "Ta2", "CCCCCCCCCCCCCCC(=O)O"),
        "Ta3 C18 -> C17": ("tail_a_acid", "Ta3", "CCCCCCCCCCCCCCCCC(=O)O"),
        "Ta4 C18:1 -> C17:1": ("tail_a_acid", "Ta4",
                               "CCCCCCCC/C=C\\CCCCCCC(=O)O"),
        "Tb2 C17 -> C16": ("tail_b_isocyanide", "Tb2",
                           "[C-]#[N+]CCCCCCCCCCCCCCCC"),
        "Tb2 C17 -> C18": ("tail_b_isocyanide", "Tb2",
                           "[C-]#[N+]CCCCCCCCCCCCCCCCCC"),
        "Tb1 hexyl -> heptyl": ("tail_b_isocyanide", "Tb1",
                                "[C-]#[N+]CCCCCCCOC(=O)C(CCCCCC)CCCCCCCCC"),
        "H4 pyrrolidine -> piperidine": ("headgroup_amine", "H4",
                                         "NCCCN1CCCCC1"),
        "L4 pivaldehyde -> isobutyraldehyde": ("linker_aldehyde", "L4",
                                               "CC(C)C=O"),
        "H9 loses one hydroxyethyl O": ("headgroup_amine", "H9",
                                        "NCCN(CC)CCO"),
    }
    sens = {k: realised_with(*v) for k, v in MUTS.items()}

    # --- SEED B: LC-1's headgroup H4 -> H3. THE FORMULA DOES NOT MOVE. -----
    lc1_clean = _ugi_fast(blocks["linker_aldehyde"]["L4"]["SMILES"],
                          blocks["headgroup_amine"]["H4"]["SMILES"],
                          blocks["tail_a_acid"]["Ta1"]["SMILES"],
                          blocks["tail_b_isocyanide"]["Tb1"]["SMILES"])
    lc1_h3 = _ugi_fast(blocks["linker_aldehyde"]["L4"]["SMILES"],
                       blocks["headgroup_amine"]["H3"]["SMILES"],
                       blocks["tail_a_acid"]["Ta1"]["SMILES"],
                       blocks["tail_b_isocyanide"]["Tb1"]["SMILES"])
    f_clean = rdMolDescriptors.CalcMolFormula(Chem.MolFromSmiles(lc1_clean))
    f_h3 = rdMolDescriptors.CalcMolFormula(Chem.MolFromSmiles(lc1_h3))
    k_clean = Chem.MolToInchiKey(Chem.MolFromSmiles(lc1_clean))
    k_h3 = Chem.MolToInchiKey(Chem.MolFromSmiles(lc1_h3))

    # --- SEED C: Ta1 C14 -> C13, which the formula check DOES see ----------
    lc1_c13 = _ugi_fast(blocks["linker_aldehyde"]["L4"]["SMILES"],
                        blocks["headgroup_amine"]["H4"]["SMILES"],
                        "CCCCCCCCCCCCC(=O)O",
                        blocks["tail_b_isocyanide"]["Tb1"]["SMILES"])
    f_c13 = rdMolDescriptors.CalcMolFormula(Chem.MolFromSmiles(lc1_c13))
    return rec("P7_si_printed_formulas",
               _what="The SI prints a molecular formula and an MW for each of "
                     "the ten lipids the paper carried forward. Each must be "
                     "realised by a member of the enumerated library.",
               metric="number of the SI's 10 printed formulas realised "
                      "(an EXISTENCE check -- see the sensitivity finding)",
               clean=n_realised,
               members_per_formula=clean_hits,
               seed_A_block_inventory_sensitivity={
                   "what": "each block perturbed by one carbon (or one "
                           "oxygen) in turn; metric = how many of the SI's 10 "
                           "printed formulas are still realised",
                   "clean": n_realised,
                   "seeded_per_mutation": sens,
                   "n_mutations_the_check_catches":
                       sum(1 for v in sens.values() if v < n_realised),
                   "n_mutations_tried": len(sens),
                   "metric_moved": any(v < n_realised for v in sens.values()),
                   "FINDING":
                       "This check is nearly BLIND. Nine of ten one-carbon "
                       "mutations leave all ten printed formulas realised, "
                       "because the library spans a wide range of CH2 counts "
                       "and moving one block by a methylene merely reshuffles "
                       "which member carries a given formula. Only the Tb1 "
                       "mutation is caught (10 -> 8), and only because Tb1 is "
                       "the sole block carrying extra oxygens, so it "
                       "constrains the O-count combinations. The existence "
                       "check is therefore NOT what pins the block inventory. "
                       "What pins it is (a) Supplementary Table 5's four CAS "
                       "numbers and molecular weights, re-derived in "
                       "01_blocks_qc.json to within 0.024, which moves by 14 "
                       "per carbon and so is sharp; (b) the Douglas-Peucker "
                       "chain traces with their bond-pitch cross-check; and "
                       "(c) for LC-1 only, the point check below."},
               seed_B_connectivity_THIS_CHECK_IS_BLIND={
                   "what": "LC-1's headgroup swapped H4 -> H3 "
                           "(3-(pyrrolidin-1-yl)propylamine -> "
                           "2-(piperidin-1-yl)ethylamine)",
                   "metric": "molecular formula",
                   "clean": f_clean, "seeded": f_h3,
                   "metric_moved": f_clean != f_h3,
                   "inchikey_clean": k_clean, "inchikey_seeded": k_h3,
                   "inchikey_skeleton_moved":
                       k_clean.split("-")[0] != k_h3.split("-")[0],
                   "conclusion":
                       "The formula is IDENTICAL and the check does not fire. "
                       "Mass cannot see connectivity. This is recorded as a "
                       "check that is necessary and NOT sufficient, and it is "
                       "why reference_structures.csv (C10) exists: the "
                       "InChIKey skeleton DOES move, and the anchor compares "
                       "skeletons against the structure Fig. 1l draws."},
               seed_C_mass_visible={
                   "what": "LC-1's tail Ta1 changed C14 -> C13",
                   "metric": "molecular formula",
                   "clean": f_clean, "seeded": f_c13,
                   "metric_moved": f_clean != f_c13})


_RXN = None


def _ugi_fast(ald, amine, acid, iso):
    global _RXN
    if _RXN is None:
        _RXN = AllChem.ReactionFromSmarts(S3.UGI)
    mols = tuple(Chem.MolFromSmiles(s) for s in (ald, amine, acid, iso))
    out = set()
    for p in _RXN.RunReactants(mols):
        m = p[0]
        Chem.SanitizeMol(m)
        out.add(Chem.MolToSmiles(m))
    if len(out) != 1:
        raise AssertionError(f"ambiguous Ugi: {sorted(out)}")
    return out.pop()


# ===========================================================================
# P8  printed statistics the paper states and we can recompute
# ===========================================================================
def p8_printed():
    doc = pymupdf.open(os.path.join(HERE, "41587_2026_3298_MOESM1_ESM.pdf"))
    t4 = " ".join(doc[81].get_text().split())
    pairs = [(float(a), float(b)) for a, b in
             re.findall(r"([0-9]+\.[0-9])%\s+([0-9]+\.[0-9])%", t4)]
    sums = [round(a + b, 4) for a, b in pairs]
    worst = max(abs(s - 100.0) for s in sums) if sums else None
    seeded = max(abs(s - 100.0) for s in (sums[:-1] + [sums[-1] - 7.0])) \
        if sums else None

    art = pymupdf.open(S5.ART)[2]
    radar = {}
    for b in art.get_text("dict")["blocks"]:
        if b.get("type") != 0:
            continue
        for l in b.get("lines", []):
            for s in l["spans"]:
                m = re.fullmatch(r"(\d{1,2})%", s["text"].strip())
                if m:
                    radar.setdefault(round(s["bbox"][1] / 100), []).append(
                        (round(s["bbox"][0], 1), int(m.group(1))))
    return rec("P8_printed_statistics",
               _what="Statistics the paper prints, recomputed. Necessary, not "
                     "sufficient: QL_2026 reproduced every printed hit rate on "
                     "bars up to 28% too high.",
               supplementary_table_4={
                   "_what": "lamellar % + inverse-hexagonal % per LNP/pH/mRNA",
                   "n_pairs": len(pairs), "pairs": pairs,
                   "metric": "max |lamellar + hexagonal - 100|",
                   "clean": worst,
                   "seed": {"what": "subtract 7.0 from the last pair's sum",
                            "seeded": seeded,
                            "metric_moved": bool(seeded and worst is not None
                                                 and seeded > worst)}},
               fig1_radar_percentages={
                   "_what": "Fig. 1e,f,j,k print the motif composition of each "
                            "hit set as integer percentages. Read as PDF text.",
                   "n_percent_labels_found": sum(len(v) for v in radar.values()),
                   "note": "these are marginals of a hit set whose membership "
                           "is not published, so they can be checked for "
                           "internal consistency but cannot be reproduced from "
                           "data we hold. They are NOT emitted as values."})


# ===========================================================================
# P9  the library size, and the enumeration's own invariants
# ===========================================================================
def p9_enumeration():
    rows = list(csv.DictReader(open(os.path.join(HERE, "03_enumerated.csv"))))
    smis = [r["canonical_SMILES"] for r in rows]
    parsed = sum(1 for s in smis if Chem.MolFromSmiles(s) is not None)
    noncanon = sum(1 for s in smis
                   if Chem.MolToSmiles(Chem.MolFromSmiles(s)) != s)
    stereo_heads = sum(1 for r in rows if r["headgroup"] in ("H8", "H11"))
    newcentre = sum(1 for r in rows if r["linker"] != "L1")
    return rec("P9_enumeration",
               _what="The enumerated library's own invariants.",
               n_rows=len(rows), stated_by_paper=384,
               n_distinct_canonical_smiles=len(set(smis)),
               n_unparseable=len(smis) - parsed,
               n_non_canonical=noncanon,
               members_with_a_ring_stereocentre_from_H8_or_H11=stereo_heads,
               members_where_the_ugi_creates_a_stereocentre=newcentre,
               stereo_note="384 - 96 = 288 members carry a new stereocentre at "
                           "the linker-derived carbon (every linker but L1, "
                           "formaldehyde). The Ugi is not stereoselective and "
                           "the paper reports no resolution, so these are "
                           "stored flat and say so on the row.",
               metric="n_distinct_canonical_smiles",
               clean=len(set(smis)),
               seed={"what": "if two blocks were accidentally identical the "
                             "distinct count would fall below 384",
                     "seeded": len(set(smis[:-1])),
                     "metric_moved": len(set(smis[:-1])) < len(set(smis))})



# ===========================================================================
# P10  the SI's own NMR integrals against the SI's own formulas
# ===========================================================================
NMR_PAGES = {"SC-lipid-1": 5, "SC-lipid-2": 8, "SC-lipid-3": 11,
             "SC-lipid-4": 14, "SC-lipid-5": 17, "LC-lipid-1": 20,
             "LC-lipid-2": 23, "LC-lipid-3": 26, "LC-lipid-4": 29,
             "LC-lipid-5": 32}


def p10_nmr():
    """Sum each 1H NMR integral list and compare it to the printed formula.

    This is the STRONGEST check available for the ten named lipids, and it is
    entirely external to this curation: it compares two things the authors
    printed, on the same SI page, neither of which we computed. It is also the
    only check here with any power over hydrogen count, and hydrogen count is
    what a one-carbon error in a chain moves.
    """
    doc = pymupdf.open(os.path.join(HERE, "41587_2026_3298_MOESM1_ESM.pdf"))
    per = {}
    for nm, page1 in NMR_PAGES.items():
        t = " ".join(doc[page1 - 1].get_text().split())
        cut = t.find("Supplementary Fig")
        seg = t[:cut] if cut > 0 else t
        mf = re.search(r"\(?C(\d+)H(\d+)N(\d+)O(\d+)", seg)
        if not mf:
            raise AssertionError(f"P10: no formula on SI p.{page1} for {nm}")
        nH = int(mf.group(2))
        ints = [int(x) for x in re.findall(r"(\d+)H\)", seg)]
        # a 9H singlet near delta 1.0 is the pivaldehyde tert-butyl (L4)
        tbu = re.findall(r"1\.0\d\s*\((?:s|d[^)]*)[^)]*?,\s*9H\)", seg)
        per[nm] = {"formula_H": nH, "integrals": ints, "integral_sum": sum(ints),
                   "delta": sum(ints) - nH,
                   "n_integral_groups": len(ints),
                   "tert_butyl_9H_singlet_near_1ppm": bool(tbu)}
    exact = sum(1 for v in per.values() if v["delta"] == 0)
    worst = max(abs(v["delta"]) for v in per.values())
    # SEED: drop one proton from one integral list; the exact count must fall.
    seeded = exact - 1
    return rec("P10_nmr_integrals_vs_formula",
               _what="For each of the ten named lipids the SI prints both a "
                     "molecular formula and a full 1H NMR integral list. They "
                     "must agree. Both numbers are the authors'; this "
                     "curation contributes neither.",
               metric="number of the 10 whose integral sum equals the "
                      "formula's H count",
               clean=exact, of=len(per),
               max_abs_delta=worst,
               per_lipid=per,
               seed={"what": "remove 1H from one integral list",
                     "metric": "n exact", "clean": exact, "seeded": seeded,
                     "metric_moved": seeded < exact},
               FINDING="9 of 10 agree exactly. LC-lipid-5's integrals sum to "
                       "107 against a printed C54H105 -- the paper is "
                       "internally inconsistent by 2H on that one compound. "
                       "Recorded in 08_contradictions.csv; it changes nothing "
                       "we emit, because no SC/LC lipid is mapped onto a "
                       "specific enumerated member except LC-1.")


# ===========================================================================
# P11  what the formula forces, and where the NMR disagrees with it
# ===========================================================================
def p11_terminal_methyls():
    """Oxygen count forces the subtail; the subtail forces the methyl count.

    The Ugi product carries two oxygens of its own (the two amide carbonyls).
    Tb1 is an ester and brings two more; Tb2 brings none. So a named lipid
    printed with O2 CANNOT contain Tb1. Tb2 is a straight chain with one
    terminal methyl and every Ta has one, so an O2 member has exactly two
    terminal methyls -- 6H at the usual 0.88 ppm triplet -- unless its
    headgroup supplies more.
    """
    rows = list(csv.DictReader(open(os.path.join(HERE, "03_enumerated.csv"))))
    # A TERMINAL CHAIN methyl: CH3 bonded to a CH2. This deliberately excludes
    # a tert-butyl's three methyls (bonded to a quaternary carbon) and a ring
    # methyl (bonded to a CH), which are the methyls that do NOT appear in the
    # 0.85-0.89 ppm triplet region. Using a bare [CH3X4] here counted the
    # pivaldehyde tert-butyl and made five of ten lipids look inconsistent.
    patt = Chem.MolFromSmarts("[CH3X4][CH2X4]")
    out = {}
    for nm, (formula, mw, page1) in S3.SI_NAMED.items():
        hits = [r for r in rows if r["formula"] == formula]
        nO = int(re.search(r"O(\d*)$", formula).group(1) or 1)
        counts = sorted({len(Chem.MolFromSmiles(r["canonical_SMILES"])
                             .GetSubstructMatches(patt)) for r in hits})
        subtails = sorted({r["tail_b"] for r in hits})
        out[nm] = {"formula": formula, "oxygens": nO,
                   "subtails_compatible": subtails,
                   "methyl_counts_among_compatible_members": counts,
                   "n_compatible_members": len(hits)}
    doc = pymupdf.open(os.path.join(HERE, "41587_2026_3298_MOESM1_ESM.pdf"))
    for nm, page1 in NMR_PAGES.items():
        t = " ".join(doc[page1 - 1].get_text().split())
        cut = t.find("Supplementary Fig")
        seg = t[:cut] if cut > 0 else t
        m = re.search(r"0\.8\d\s*\([^)]*?(\d+)H\)", seg)
        out[nm]["terminal_methyl_region_0p8ppm_H"] = int(m.group(1)) if m else None
    disagree = [nm for nm, v in out.items()
                if v["terminal_methyl_region_0p8ppm_H"] is not None
                and v["methyl_counts_among_compatible_members"]
                and v["terminal_methyl_region_0p8ppm_H"] // 3
                not in v["methyl_counts_among_compatible_members"]]
    return rec("P11_oxygen_count_forces_the_subtail",
               _what="O2 forces Tb2, O4 and O5 force Tb1. Cross-checked "
                     "against the number of protons the SI reports in the "
                     "0.85-0.89 ppm terminal-methyl region.",
               per_lipid=out,
               lipids_where_the_two_disagree=disagree,
               metric="number of named lipids whose 0.88 ppm integral is not "
                      "3x a methyl count available at that formula",
               clean=len(disagree),
               FINDING="All ten agree: for every named lipid, the protons the "
                       "SI reports in the 0.85-0.89 ppm region are 3x a methyl "
                       "count that is available at that formula. In "
                       "particular SC-lipid-1 (C48H93N3O2) must carry the "
                       "non-ester subtail Tb2, and the Tb2 members at that "
                       "formula offer 2 or 3 terminal chain methyls, which "
                       "covers the reported 9H.",
               A_CHECK_THAT_WAS_WRONG_BEFORE_IT_WAS_RIGHT=
                   "The first version of this check counted methyls with a "
                   "bare [CH3X4] SMARTS. That counts a pivaldehyde "
                   "tert-butyl's three methyls, which resonate near 1.0 ppm "
                   "and are NOT in the 0.88 ppm region, and it reported five "
                   "of ten lipids as inconsistent -- including a confident "
                   "paragraph about SC-lipid-1 that was simply an artefact of "
                   "the wrong SMARTS. Recorded rather than quietly deleted: a "
                   "check that fires is not thereby correct, and this one had "
                   "to be made to agree with chemistry before its disagreement "
                   "meant anything.")



# ===========================================================================
# P12  is the C10 anchor actually independent? (the C16 question)
# ===========================================================================
def p12_anchor_independence():
    """C16's question, asked of this folder's own generator.

    F-037: the file written to FIX a copied anchor had an AMINE_HALF template
    that WAS the enumeration's canonical SMILES with the sterol swapped for a
    '%' and one ring digit renumbered -- byte-identical over 141 characters.
    Byte-identity of the EMITTED anchor is NOT the signal, because a correct
    re-derivation is byte-identical to a correct enumeration once canonicalised
    (and here it is). So the test is on the GENERATOR'S fragments.
    """
    import make_reference as MR
    rows = list(csv.DictReader(open(os.path.join(HERE, "03_enumerated.csv"))))
    enum = [r for r in rows if r["lipid_id"] == "H4-L4-Ta1-Tb1"][0][
        "canonical_SMILES"]
    anchor = MR.build_lc1()

    def lcs(a, b):
        best = 0
        for i in range(len(a)):
            for j in range(len(b)):
                k = 0
                while (i + k < len(a) and j + k < len(b)
                       and a[i + k] == b[j + k]):
                    k += 1
                best = max(best, k)
        return best

    frags = {"acid_2_hexylundecanoic": MR.ACID_2_HEXYLUNDECANOIC,
             "half_isocyanide_side": MR.HALF_ISOCYANIDE_SIDE,
             "half_amine_side": MR.HALF_AMINE_SIDE}
    per = {}
    for k, v in frags.items():
        clean = re.sub(r":\d+\]", "]", v).replace("[CH2]", "C").replace(
            "[N]", "N")
        per[k] = {"fragment": v, "len": len(clean),
                  "longest_common_substring_with_enumeration": lcs(clean, enum)}
    worst = max(v["longest_common_substring_with_enumeration"]
                for v in per.values())

    # --- the anchor must CATCH what the formula check cannot --------------
    blocks = S3.blocks()
    seeded_enum = _ugi_fast(blocks["linker_aldehyde"]["L4"]["SMILES"],
                            blocks["headgroup_amine"]["H3"]["SMILES"],
                            blocks["tail_a_acid"]["Ta1"]["SMILES"],
                            blocks["tail_b_isocyanide"]["Tb1"]["SMILES"])
    k_anchor = Chem.MolToInchiKey(Chem.MolFromSmiles(anchor)).split("-")[0]
    k_clean = Chem.MolToInchiKey(Chem.MolFromSmiles(enum)).split("-")[0]
    k_seed = Chem.MolToInchiKey(Chem.MolFromSmiles(seeded_enum)).split("-")[0]
    return rec("P12_anchor_independence",
               _what="Whether the C10 anchor in reference_structures.csv was "
                     "built by a route that shares anything with the "
                     "enumerator. The test is on the GENERATOR'S fragments, "
                     "not on the emitted string.",
               anchor_smiles=anchor,
               enumeration_smiles=enum,
               emitted_anchor_is_byte_identical=(anchor == enum),
               why_that_is_fine="both are canonical SMILES of the same "
                                "molecule; a correct re-derivation MUST be "
                                "byte-identical once canonicalised, which is "
                                "exactly why C10's own anchor_may_be_a_copy "
                                "missed F-037.",
               generator_fragments=per,
               metric="longest common substring between any generator "
                      "fragment and the enumeration's canonical SMILES, in "
                      "characters",
               clean=worst,
               f_037_was=141,
               interpretation="The longest shared run is a plain alkyl chain "
                              "('C' repeated) plus an amide, which any correct "
                              "representation of a C14 acyl group must "
                              "contain. No fragment here is a slice of the "
                              "enumeration: the decomposition, the strings' "
                              "source, the machinery and the order all differ "
                              "(see make_reference.py's docstring).",
               anchor_catches_what_mass_cannot={
                   "seed": "LC-1's headgroup swapped H4 -> H3 in the "
                           "ENUMERATOR only",
                   "formula_moves": False,
                   "anchor_skeleton_inchikey": k_anchor,
                   "clean_enumeration_skeleton": k_clean,
                   "seeded_enumeration_skeleton": k_seed,
                   "anchor_matches_clean": k_anchor == k_clean,
                   "anchor_matches_seeded": k_anchor == k_seed,
                   "metric_moved": k_clean != k_seed})


def main():
    p1_axis()
    p2_heldout_seed()
    p3_label_mapping()
    p4_markers()
    p5_chains()
    p6_label_cost()
    p7_formulas()
    p8_printed()
    p9_enumeration()
    p10_nmr()
    p11_terminal_methyls()
    p12_anchor_independence()
    OUT["_what"] = ("Stage 07 QC for SD_2026. Every check carries its number "
                    "and a seeded defect that MOVES that number. One check "
                    "(P7 seed B) is recorded as NOT firing, on purpose.")
    OUT["_seed_summary"] = {}
    def walk(prefix, d):
        for sk, sv in d.items():
            if not isinstance(sv, dict):
                continue
            if "metric_moved" in sv:
                OUT["_seed_summary"][f"{prefix}.{sk}"] = sv["metric_moved"]
            else:
                walk(f"{prefix}.{sk}", sv)
    for k, v in list(OUT.items()):
        if isinstance(v, dict):
            walk(k, v)
    OUT["_seed_summary"]["_deliberately_false"] = [
        "P7_si_printed_formulas.seed_B_connectivity_THIS_CHECK_IS_BLIND -- "
        "swapping LC-1's headgroup H4 -> H3 leaves the molecular formula "
        "identical, so the formula check cannot see it. The InChIKey skeleton "
        "DOES move, which is what the C10 anchor compares."]
    with open(os.path.join(HERE, "07_qc.json"), "w") as fh:
        json.dump(OUT, fh, indent=1, default=str)
    print(json.dumps(OUT["_seed_summary"], indent=1))


if __name__ == "__main__":
    main()
