#!/usr/bin/env python3
"""XQ_2020 stage 07 -- every check, with numbers, and a gate.

    python3 -I stage07_qc.py              # run the checks, write 07_qc.json
    python3 -I stage07_qc.py --seed LIST  # seed defects, prove the checks fire
    python3 -I stage07_qc.py --self-test  # seed each one in turn, report

EXIT CODE 1 if any GATING check fails. Verify it directly with `echo $?`, not
through a pipe: `cmd | tail` reports the exit status of `tail`.

WHY A --seed MODE EXISTS. CLAUDE.md section 6: "A NEW CHECK IS NOT DONE UNTIL
IT HAS BEEN MADE TO FAIL ON PURPOSE." Three incidents paid for that rule, one
of which was three ZT_2026 marginal checks scoring a literal 0 for twelve days
while the summary beside them said "ALL FOUR EXACT". So every check here has a
seedable defect, and the self-test credits a check only when the AGREEMENT
COUNT MOVES -- not merely when something prints.
"""
import argparse
import csv
import json
import math
import os
import re
import sys
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))

import numpy as np                                                  # noqa: E402
from rdkit import Chem, RDLogger                                    # noqa: E402
from rdkit.Chem.Descriptors import ExactMolWt                       # noqa: E402
from rdkit.Chem.rdMolDescriptors import CalcMolFormula              # noqa: E402

RDLogger.DisableLog("rdApp.*")
PROTON = 1.00727646

# ---------------------------------------------------------------------------
# Printed statistics. Every one is a PDF/JATS TEXT SPAN, never a pixel.
# ---------------------------------------------------------------------------
PRINTED_POINTS = {
    # (amine, suffix): (printed %, where it is printed)
    ("113", ""):    (29, "RESULTS: '113-O16CBA and 306-O16CBA, induced 29% and 44% GFP+ cells respectively'"),
    ("306", ""):    (44, "same sentence"),
    ("400", ""):    (68, "RESULTS: 'the most active 400-O16CBA induced 68% GFP+ cells'"),
    ("81", "-F1"):  (93, "RESULTS: '81-O16CBA-F1 being the highest in producing ~93% GFP+ cells'"),
    ("93", "-F1"):  (2,  "RESULTS: 'all LNPs, except for 93-O16CBA-F1 (~2% GFP+ cells)'"),
    ("81", "-F2"):  (78, "RESULTS: '81-O16CBA-F2 was determined to be the most efficient ... ~78% GFP+ cells'"),
}
PRINTED_LPF2K = (91, "RESULTS: 'Lpf2k was highly efficient ... ~91% GFP+ cells were recorded'")

# ESI-MS of the acid-degradation products, read as SI text spans.
# R-O16CBA-1  = one benzylidene acetal hydrolysed (lose one HexDMBA, gain H2O)
# R-O16CBA-0  = both hydrolysed
DEGRADATION = {
    ("75", 1): (713.55, "RESULTS, Figure 3E"),
    ("76", 1): (711.64, "SI p. S6, Figure S5 label"),
    ("76", 0): (463.36, "SI p. S6, Figure S5 label"),
    ("77", 1): (725.64, "SI p. S7, Figure S6 label"),
    ("77", 0): (477.45, "SI p. S7, Figure S6 label"),
}
HEXDMBA = "O=Cc1c(OC)cc(OCCCCCC)cc1OC"
WATER = ExactMolWt(Chem.MolFromSmiles("O"))

# The ion-trap tolerance. The Experimental Section names a Finnigan LTQ, a
# linear ion trap: unit-resolution, NOT an accurate-mass instrument. Their own
# calculated-vs-observed spread (up to 0.15 Da) sets the scale, so the
# tolerance is taken from their data rather than chosen.
MZ_TOL = 0.30


def load():
    vals = list(csv.DictReader(open(os.path.join(HERE, "05_values.csv"))))
    enum = list(csv.DictReader(open(os.path.join(HERE, "03_enumerated.csv"))))
    b01 = json.load(open(os.path.join(HERE, "01_blocks_qc.json")))
    q03 = json.load(open(os.path.join(HERE, "03_enumerated_qc.json")))
    q05 = json.load(open(os.path.join(HERE, "05_values_qc.json")))
    return vals, enum, b01, q03, q05


def grid(vals):
    g = {}
    for r in vals:
        suf = "" if r["formulation_suffix"] == "(none)" else r["formulation_suffix"]
        g[(r["amine_head"], suf)] = float(r["value_pct_GFP_positive"])
    return g


def build(seed=()):
    vals, enum, b01, q03, q05 = load()
    g = grid(vals)
    enum_by = {e["head_amine"]: dict(e) for e in enum}

    # ------------------------------------------------------------------ seeds
    if "value_shift" in seed:          # a 5-point bias on every bar
        g = {k: v + 5.0 for k, v in g.items()}
    if "swap_groups" in seed:          # the x-axis label order is wrong
        for suf in ("", "-F1", "-F2"):
            g[("113", suf)], g[("306", suf)] = g[("306", suf)], g[("113", suf)]
    if "wrong_arm_count" in seed:      # 113 built with 2 arms instead of 3
        e = dict(enum_by["113"])
        m = Chem.MolFromSmiles(e["IL_SMILES"])
        e["IL_SMILES"] = q03["members"][8]["head_SMILES"]
        enum_by["113"] = e
    if "wrong_head" in seed:           # amine 80 given 78's structure
        e = dict(enum_by["80"]); e["IL_SMILES"] = enum_by["78"]["IL_SMILES"]
        enum_by["80"] = e
    if "bad_tick" in seed:             # one y tick displaced by 6 px
        q05 = json.loads(json.dumps(q05))
        q05["y_axis"]["tick_rows_px"][1] += 6.0
    if "reagent_mw" in seed:           # the hexyl chain called an octyl
        b01 = json.loads(json.dumps(b01))
        for r in b01["named_reagent_mw_crosscheck"]["rows"]:
            if r.get("reagent") == "1-bromohexane":
                r["MW_of_the_structure_we_assign"] = 193.12
                r["deviation_pct"] = round(
                    abs(r["MW_implied_by_their_own_g_per_mmol"] - 193.12)
                    / 193.12 * 100, 2)
    if "their_calc_typo" in seed:      # a transcription slip in Fig 2F's Cal.
        enum_by = dict(enum_by)
        e = dict(enum_by["78"]); e["MplusH_their_calc"] = "933.97"
        enum_by["78"] = e
        enum = [e if r["head_amine"] == "78" else r for r in enum]
    if "head75" in seed:               # amine 75 given 78's structure
        e = dict(enum_by["75"]); e["IL_SMILES"] = enum_by["78"]["IL_SMILES"]
        enum_by["75"] = e
    if "f1_drop" in seed:              # the F1 leader is no longer the leader
        g = dict(g); g[("81", "-F1")] = 48.0
    if "f2_swap" in seed:
        g = dict(g); g[("82", "-F2")] = 99.0
    if "plain_swap" in seed:
        g = dict(g); g[("306", "")] = 99.0
    if "sign_flip" in seed:            # break the printed pairwise signs
        g = dict(g); g[("77", "-F1")] = 1.0; g[("76", "-F2")] = 99.0
    if "minor_tick" in seed:
        q05 = json.loads(json.dumps(q05))
        q05["y_axis"]["tick_length_probe"]["probe_rows"][0] += 4.0
    if "baseline_off" in seed:
        q05 = json.loads(json.dumps(q05))
        q05["y_axis"]["zero_tick_lands_on_the_drawn_baseline"] = False
        q05["y_axis"]["baseline_vs_zero_tick_px"] = 7.4
    if "flank" in seed:                # an error-bar cap read as a bar top
        q05 = json.loads(json.dumps(q05))
        q05["series_geometry"]["red"]["worst_flank_delta_px"] = 9.0
    if "lost_bar" in seed:
        q05 = json.loads(json.dumps(q05))
        q05["series_geometry"]["blue"]["n_windows"] = 10
        q05["series_geometry"]["blue"]["pitch_px"][3] = 146
    if "ambiguous" in seed:
        q03 = json.loads(json.dumps(q03))
        q03["members"][8]["n_distinct_products"] = 2
    if "acetal_sym" in seed:
        q03 = json.loads(json.dumps(q03))
        q03["acetal_regiochemistry"]["all_three_equivalent"] = False
        q03["acetal_regiochemistry"]["symmetry_ranks_of_the_three_oxygens"] = [1, 2, 3]
    if "arm_routes" in seed:
        q03 = json.loads(json.dumps(q03))
        q03["arm_fragment"]["agree"] = False
    if "rebuild" in seed:
        q03 = json.loads(json.dumps(q03))
        q03["o16cba_rebuilt_from_HexDMBAH_plus_acrylic_acid"]["match"] = False
    if "isomer" in seed:               # 77 and 82 collapsed onto one structure
        e = dict(enum_by["82"]); e["IL_SMILES"] = enum_by["77"]["IL_SMILES"]
        enum_by["82"] = e
    if "corpus_head" in seed:          # our amine 75 disagrees with QL_2026's
        q03 = json.loads(json.dumps(q03))
        for m in q03["members"]:
            if m["id"] == "75":
                m["head_SMILES"] = "CCN(CC)CCCN"

    checks = []

    def add(cid, what, gating, measured, expected, ok, detail=None):
        checks.append({"id": cid, "what": what, "gating": gating,
                       "measured": measured, "expected": expected,
                       "pass": bool(ok), "detail": detail})

    # ---------------------------------------------------------------- P1 mass
    # THE THREE-WAY SPLIT. Collapsing these into one comparison is how a
    # published typo gets absorbed into a passing check.
    mass = []
    for e in enum:
        hid = e["head_amine"]
        m = Chem.MolFromSmiles(enum_by[hid]["IL_SMILES"])
        ours_formula = CalcMolFormula(m)
        ours_mz = ExactMolWt(m) + PROTON
        their_formula = e["MplusH_their_calc"], e["MplusH_their_obs"]
        theirs_cal = float(e["MplusH_their_calc"])
        theirs_obs = float(e["MplusH_their_obs"])
        # their FORMULA, independently re-weighed by us
        ionf = [r for r in q03["members"] if r["id"] == hid][0]["their_formula_MplusH"]
        c = Counter()
        for el, n in re.findall(r"([A-Z][a-z]?)(\d*)", ionf):
            if el:
                c[el] += int(n) if n else 1
        w = {"C": 12.0, "H": 1.00782503, "N": 14.0030740, "O": 15.9949146}
        their_formula_mz = sum(w[el] * n for el, n in c.items()) - 0.000548579
        mass.append({
            "lipidoid": e["lipidoid"],
            "A_ours_vs_their_formula_mz": round(ours_mz - their_formula_mz, 4),
            "B_their_formula_vs_their_calc_mz": round(
                their_formula_mz - theirs_cal, 4),
            "C_ours_vs_their_observed_mz": round(ours_mz - theirs_obs, 4),
            "our_formula_matches_theirs_atom_for_atom":
                CalcMolFormula(Chem.AddHs(m)).replace("H", "H") is not None
                and _formula_plus_proton(ours_formula) == ionf,
        })
    add("P1a", "our structure's [M+H]+ vs the m/z implied by the molecular "
               "FORMULA the authors print in Figure 2F", True,
        f"max |delta| = {max(abs(r['A_ours_vs_their_formula_mz']) for r in mass):.4f} Da",
        f"<= 0.01 Da (same formula => same mass)",
        max(abs(r["A_ours_vs_their_formula_mz"]) for r in mass) <= 0.01, mass)
    add("P1b", "the authors' own printed FORMULA vs their own printed "
               "CALCULATED m/z -- a check on THEIR arithmetic, not ours", True,
        f"max |delta| = {max(abs(r['B_their_formula_vs_their_calc_mz']) for r in mass):.4f} Da",
        "<= 0.02 Da",
        max(abs(r["B_their_formula_vs_their_calc_mz"]) for r in mass) <= 0.02)
    add("P1c", "our structure's [M+H]+ vs the authors' OBSERVED m/z", True,
        f"max |delta| = {max(abs(r['C_ours_vs_their_observed_mz']) for r in mass):.4f} Da",
        f"<= {MZ_TOL} Da (Finnigan LTQ linear ion trap, unit resolution)",
        max(abs(r["C_ours_vs_their_observed_mz"]) for r in mass) <= MZ_TOL)
    add("P1d", "our molecular formula equals the authors' printed formula, "
               "atom for atom, for all eleven lipidoids", True,
        f"{sum(1 for r in mass if r['our_formula_matches_theirs_atom_for_atom'])}/11",
        "11/11",
        sum(1 for r in mass if r["our_formula_matches_theirs_atom_for_atom"]) == 11)

    # --------------------------------------------- P2 named-reagent MW ladder
    cr = b01["named_reagent_mw_crosscheck"]
    worst = max(r["deviation_pct"] for r in cr["rows"] if "deviation_pct" in r)
    n_ok = sum(1 for r in cr["rows"]
               if "deviation_pct" in r and r["deviation_pct"] <= 2.5)
    add("P2", "every reagent named in the Experimental Section: the MW implied "
              "by its own printed grams/millimoles vs the MW of the structure "
              "we assign it. Pins the hexyl chain with no pixel measurement.",
        True, f"{n_ok}/9 within 2.5%, worst {worst:.2f}%", "9/9",
        n_ok == 9, cr["rows"])

    # ------------------------------------- P3 degradation fragments (arm count)
    deg = []
    hex_m = ExactMolWt(Chem.MolFromSmiles(HEXDMBA))
    for (hid, n_intact), (mz, where) in DEGRADATION.items():
        m = Chem.MolFromSmiles(enum_by[hid]["IL_SMILES"])
        n_lost = 2 - n_intact
        pred = ExactMolWt(m) - n_lost * hex_m + n_lost * WATER + PROTON
        deg.append({"fragment": f"{hid}-O16CBA-{n_intact}", "printed_mz": mz,
                    "our_mz": round(pred, 3), "delta": round(pred - mz, 3),
                    "where": where})
    add("P3", "the acid-degradation fragments. Each hydrolysed acetal releases "
              "one HexDMBA and takes up one water, so the fragment masses are "
              "an INDEPENDENT measure of how many tails a lipidoid has -- "
              "evidence the parent-ion mass alone cannot give.", True,
        f"{sum(1 for d in deg if abs(d['delta']) <= MZ_TOL)}/5 within "
        f"{MZ_TOL} Da, worst {max(abs(d['delta']) for d in deg):.3f}",
        "5/5", all(abs(d["delta"]) <= MZ_TOL for d in deg), deg)

    # ------------------------------------------ P4 printed point values (Fig 5)
    pts = []
    for (hid, suf), (printed, where) in PRINTED_POINTS.items():
        ours = g[(hid, suf)]
        pts.append({"item": f"{hid}-O16CBA{suf}", "printed": printed,
                    "ours": round(ours, 2), "delta": round(ours - printed, 2),
                    "where": where})
    lp = q05["controls_right_of_the_axis_break"]["measured"]
    lpf = [v for k, v in lp.items() if "value_pct" in v]
    pts.append({"item": "Lpf2k (positive control)", "printed": PRINTED_LPF2K[0],
                "ours": lpf[0]["value_pct"] if lpf else None,
                "delta": round(lpf[0]["value_pct"] - PRINTED_LPF2K[0], 2)
                if lpf else None, "where": PRINTED_LPF2K[1]})
    worst_pt = max(abs(p["delta"]) for p in pts if p["delta"] is not None)
    add("P4", "our digitised bar heights against the seven values the RESULTS "
              "text prints in words. The prose rounds to whole percent and "
              "writes '~' on three of them, so 1.0 is the tightest honest "
              "tolerance.", True,
        f"{sum(1 for p in pts if p['delta'] is not None and abs(p['delta']) <= 1.0)}/7 "
        f"within 1.0 pp, worst {worst_pt:.2f} pp", "7/7",
        worst_pt <= 1.0, pts)

    # ------------------------------------------------- P5 printed COUNTS
    below15 = [k for k in g if k[1] == "" and g[k] < 15]
    above50 = [k for k in g if k[1] == "-F1" and g[k] > 50]
    add("P5a", "'eight of the LNP formulations (75-O16CBA and 76-O16CBA etc.) "
               "had efficiencies below 15%' -- counted over the eleven "
               "helper-lipid-free bars", True,
        len(below15), 8, len(below15) == 8,
        sorted(k[0] for k in below15))
    add("P5b", "'Nine of the R-O16CBA-F1 formulations had >50% delivery "
               "efficiency'", True,
        len(above50), 9, len(above50) == 9, sorted(k[0] for k in above50))

    # ------------------------------------------------- P6 printed RANGES
    f2_band = [g[(h, "-F2")] for h in ("113", "306", "400")]
    add("P6a", "'113-, 306-, and 400-O16CBA-F2 were still among the top LNPs, "
               "as their delivery efficacies were determined to be 57-73%'. "
               "BOTH endpoints are tested, not just containment.", True,
        f"[{min(f2_band):.2f}, {max(f2_band):.2f}]", "[57, 73], both ends",
        56.5 <= min(f2_band) <= 57.5 and 72.5 <= max(f2_band) <= 73.5,
        {h: round(g[(h, '-F2')], 2) for h in ("113", "306", "400")})

    f1_band = {h: g[(h, "-F1")] for h in ("76", "77", "113", "306", "400")}
    in_band = [h for h, v in f1_band.items() if 77.5 <= v <= 83.5]
    add("P6b", "'Other top formulations such as 76-, 77-, 113-, 306-, and "
               "400-O16CBA-F1 induced delivery efficacies of 78-83%'",
        False,                      # see 08_contradictions.csv, row C-1
        f"{len(in_band)}/5 inside the printed band", "5/5",
        len(in_band) == 5,
        {h: round(v, 2) for h, v in f1_band.items()})

    # ------------------------------------------------- P7 printed MAXIMA
    add("P7a", "'81-O16CBA-F1 being the highest'", True,
        max(((h, round(v, 2)) for (h, s), v in g.items() if s == "-F1"),
            key=lambda t: t[1])[0], "81",
        max(((h, v) for (h, s), v in g.items() if s == "-F1"),
            key=lambda t: t[1])[0] == "81")
    add("P7b", "'81-O16CBA-F2 was determined to be the most efficient in this "
               "formulated library'", True,
        max(((h, round(v, 2)) for (h, s), v in g.items() if s == "-F2"),
            key=lambda t: t[1])[0], "81",
        max(((h, v) for (h, s), v in g.items() if s == "-F2"),
            key=lambda t: t[1])[0] == "81")
    add("P7c", "'the most active 400-O16CBA' among the helper-lipid-free LNPs",
        True,
        max(((h, round(v, 2)) for (h, s), v in g.items() if s == ""),
            key=lambda t: t[1])[0], "400",
        max(((h, v) for (h, s), v in g.items() if s == ""),
            key=lambda t: t[1])[0] == "400")

    # ------------------------------- P8 the three printed SIGN statements
    heads = [e["head_amine"] for e in enum]
    s1_exc = [h for h in heads if not g[(h, "-F1")] > g[(h, "")]]
    s2_exc = [h for h in heads if not g[(h, "-F2")] > g[(h, "")]]
    s3_exc = [h for h in heads if not g[(h, "-F2")] < g[(h, "-F1")]]
    add("P8a", "'all LNPs, except for 93-O16CBA-F1, showed higher efficacies "
               "than their corresponding R-O16CBA formulations' -- 11 pairwise "
               "signs", True, s1_exc, ["93"], s1_exc == ["93"])
    add("P8b", "'except for the 93-O16CBA-F2 and 400-O16CBA-F2, all other LNPs "
               "had higher efficacies than R-O16CBA'", True,
        sorted(s2_exc), sorted(["93", "400"]), sorted(s2_exc) == ["400", "93"])
    add("P8c", "'the delivery efficacies of R-O16CBA-F2 were slightly lower "
               "(except for 75- and 78-O16CBA-F2)'", True,
        sorted(s3_exc), sorted(["75", "78"]), sorted(s3_exc) == ["75", "78"])

    # ---------------------------------------------- P9 geometry and axis
    y = q05["y_axis"]
    ticks = y["tick_rows_px"]
    import axis as axmod
    loo = axmod.loo_residuals(ticks)
    A = np.column_stack([np.arange(len(ticks), dtype=float),
                         np.ones(len(ticks))])
    coef, *_ = np.linalg.lstsq(A, np.sort(np.asarray(ticks, float)),
                               rcond=None)
    resid = float(np.max(np.abs(np.sort(np.asarray(ticks, float))
                                - A @ coef)))
    add("P9a", "y-axis tick fit: maximum residual of the four labelled ticks",
        True, f"{resid:.3f} px", f"<= {axmod.TICK_RESIDUAL_PX} px",
        resid <= axmod.TICK_RESIDUAL_PX)
    add("P9b", "y-axis LEAVE-ONE-OUT residuals. calibrate's max-residual test "
               "is weak against a single bad tick -- an end tick has leverage "
               "0.2, so four fifths of its error is absorbed into the fit.",
        True, f"worst {max(loo):.3f} px", f"<= {axmod.LOO_TOL_PX} px",
        max(loo) <= axmod.LOO_TOL_PX, [round(v, 3) for v in loo])

    # THE ROUND NUMBER THE AUTHORS DREW AND NEVER PRINTED.
    cal_px_per = (ticks[-1] - ticks[0]) / 90.0
    zero_px = max(ticks)
    minors = [27.41, 119.5, 211.64, 303.64]        # from stage05's probe
    minors = q05["y_axis"]["tick_length_probe"]["probe_rows"]
    minor_rows = [p for p in minors
                  if min(abs(p - t) for t in ticks) > 2.0]
    minor_vals = [(zero_px - p) / cal_px_per for p in minor_rows]
    nearest = [round(v / 15.0) * 15.0 for v in minor_vals]
    dev = [abs(v - n) for v, n in zip(minor_vals, nearest)]
    add("P9c", "THE UNPRINTED ROUND NUMBER. The y axis carries four MINOR "
               "ticks that bear no label and were excluded from the fit: the "
               "scale is set only by the four long, labelled ticks. Predicting "
               "the minor ticks' values is therefore external to the "
               "calibration in the way a recomputed marginal is not. They "
               "should land on 15, 45, 75 and 105.",
        True,
        f"{[round(v, 2) for v in minor_vals]}, worst deviation "
        f"{max(dev):.3f} percentage points",
        "15.0, 45.0, 75.0, 105.0", max(dev) <= 0.25,
        {"minor_tick_rows_px": minor_rows,
         "predicted_values": [round(v, 3) for v in minor_vals]})

    add("P9d", "the 0 tick predicted by the fit lands on the DRAWN baseline "
               "rule, which took no part in the fit", True,
        f"{y['baseline_vs_zero_tick_px']} px", "< 2 px",
        y["zero_tick_lands_on_the_drawn_baseline"])

    geo = q05["series_geometry"]
    add("P9e", "bars.bar_top left/right FLANK AGREEMENT on all 33 library "
               "bars. The error bars on this figure are drawn in the BAR'S OWN "
               "COLOUR, so the colour mask alone does not separate them; the "
               "flank test is what makes an error-bar cap unreadable as a bar "
               "top (error class 32 / F-018).", True,
        f"worst {max(v['worst_flank_delta_px'] for v in geo.values())} px",
        f"<= {2.0} px",
        max(v["worst_flank_delta_px"] for v in geo.values()) <= 2.0,
        {k: v["flank_deltas_px"] for k, v in geo.items()})
    add("P9f", "every library bar detected, none refused", True,
        sum(v["n_windows"] for v in geo.values()), 33,
        sum(v["n_windows"] for v in geo.values()) == 33
        and all(v["n_refused"] == 0 for v in geo.values()))
    add("P9g", "bar pitch and width uniform within each series", True,
        {k: {"pitch": sorted(set(v["pitch_px"])),
             "width": sorted(set(v["widths_px"]))} for k, v in geo.items()},
        "pitch spread <= 2 px, width spread <= 2 px",
        all(max(v["pitch_px"]) - min(v["pitch_px"]) <= 2
            and max(v["widths_px"]) - min(v["widths_px"]) <= 2
            for v in geo.values()))

    # ------------------------------------------- P10 structure-level checks
    add("P10a", "the regiochemistry of every lipidoid is FORCED: exhaustive "
                "enumeration of all ways to place the formula-derived number "
                "of Michael arms on the head's N-H sites returns exactly one "
                "distinct constitution", True,
        f"{sum(1 for m in q03['members'] if m['n_distinct_products'] == 1)}/11",
        "11/11",
        all(m["n_distinct_products"] == 1 for m in q03["members"]),
        {m["id"]: m["products_per_step"] for m in q03["members"]})
    add("P10b", "the acetal ring closure is not a choice: the three CH2OH "
                "groups of 1,1,1-tris(hydroxymethyl)ethane have identical "
                "canonical rank, so all three pairings give one constitution",
        True, q03["acetal_regiochemistry"]["symmetry_ranks_of_the_three_oxygens"],
        "all equal", q03["acetal_regiochemistry"]["all_three_equivalent"])
    add("P10c", "the Michael arm built by two routes that share no code "
                "(reduce O16CBA's alkene / esterify propanoic acid onto "
                "HexDMBAH) gives the same fragment", True,
        q03["arm_fragment"]["agree"], True, q03["arm_fragment"]["agree"])
    add("P10d", "O16CBA rebuilt from its own named precursors equals the tail "
                "SMILES written by hand in stage 01", True,
        q03["o16cba_rebuilt_from_HexDMBAH_plus_acrylic_acid"]["match"], True,
        q03["o16cba_rebuilt_from_HexDMBAH_plus_acrylic_acid"]["match"])
    add("P10e", "77-O16CBA and 82-O16CBA are constitutional isomers sharing "
                "one formula and one printed m/z -- our two SMILES must "
                "DIFFER while their masses agree. Mass cannot see "
                "connectivity.", True,
        {"same_mass": abs(
            ExactMolWt(Chem.MolFromSmiles(enum_by["77"]["IL_SMILES"]))
            - ExactMolWt(Chem.MolFromSmiles(enum_by["82"]["IL_SMILES"]))) < 1e-6,
         "different_structure":
            enum_by["77"]["IL_SMILES"] != enum_by["82"]["IL_SMILES"]},
        "same mass, different structure",
        abs(ExactMolWt(Chem.MolFromSmiles(enum_by["77"]["IL_SMILES"]))
            - ExactMolWt(Chem.MolFromSmiles(enum_by["82"]["IL_SMILES"]))) < 1e-6
        and enum_by["77"]["IL_SMILES"] != enum_by["82"]["IL_SMILES"])

    # ------------------------------- P11 four heads corroborated by the corpus
    corpus_heads = {"75": "CCN(CC)CCN", "80": "CN(C)CCCN",
                    "113": "CN(CCN)CCN", "400": "CN(C)CCCNCCCN"}
    ours_heads = {m["id"]: m["head_SMILES"] for m in q03["members"]}
    agree = [h for h, s in corpus_heads.items()
             if Chem.CanonSmiles(s) == Chem.CanonSmiles(ours_heads[h])]
    add("P11", "four of the eleven amine heads are also used by QL_2026, a "
               "SEPARATELY curated paper that took them from a vendor "
               "catalogue by chemical name rather than from a drawing. Our "
               "structures, read off Figure 1C, must equal theirs.", True,
        f"{len(agree)}/4", "4/4", len(agree) == 4,
        {"source": "curation/ql_2026/liu2026_adhm202502225_blocks.csv rows "
                   "H2, H3, H4, H5", "agreeing": agree})

    return checks, {"vals": vals, "grid": g, "enum": enum, "q03": q03,
                    "q05": q05, "mass": mass, "deg": deg, "points": pts}


def _formula_plus_proton(neutral):
    c = Counter()
    for el, n in re.findall(r"([A-Z][a-z]?)(\d*)", neutral):
        if el:
            c[el] += int(n) if n else 1
    c["H"] += 1
    order = ["C", "H", "N", "O"]
    return "".join(f"{el}{c[el] if c[el] != 1 else ''}"
                   for el in order if c.get(el))


SEEDS = ["value_shift", "swap_groups", "wrong_arm_count", "wrong_head",
         "bad_tick", "reagent_mw", "their_calc_typo", "head75", "f1_drop",
         "f2_swap", "plain_swap", "sign_flip", "minor_tick", "baseline_off",
         "flank", "lost_bar", "ambiguous", "acetal_sym", "arm_routes",
         "rebuild", "isomer", "corpus_head"]


def failures(checks, gating_only=True):
    return [c["id"] for c in checks
            if not c["pass"] and (c["gating"] or not gating_only)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", default="")
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args()

    if a.self_test:
        base, _ = build()
        n_base = sum(1 for c in base if c["pass"])
        report = {"baseline_agreements": n_base,
                  "baseline_gating_failures": failures(base)}
        rows = []
        for s in SEEDS:
            ch, _ = build(seed=(s,))
            n = sum(1 for c in ch if c["pass"])
            moved = n_base - n
            rows.append({
                "seed": s,
                "agreements_without_the_seed": n_base,
                "agreements_with_the_seed": n,
                "agreement_count_moved_by": moved,
                "checks_that_fired": [c["id"] for c in ch if not c["pass"]
                                      and c["id"] not in failures(base, False)],
                "credited": moved > 0})
        report["seeds"] = rows
        report["all_seeds_move_the_agreement_count"] = all(
            r["credited"] for r in rows)
        print(json.dumps(report, indent=1))
        return 0 if report["all_seeds_move_the_agreement_count"] else 1

    seeds = tuple(x for x in a.seed.split(",") if x)
    checks, ctx = build(seed=seeds)
    bad = failures(checks)

    out = {
        "_what": "stage 07 QC for XQ_2020. Every check carries its NUMBERS. "
                 "`gating` checks set the exit code; a non-gating check is one "
                 "whose failure is a statement about the PAPER, and it names "
                 "the row of 08_contradictions.csv that records it.",
        "_how_to_verify_the_gate": "python3 -I stage07_qc.py; echo $?  -- "
                                   "directly, never through a pipe, because "
                                   "`cmd | tail` reports tail's exit status.",
        "seeds_applied": list(seeds),
        "n_checks": len(checks),
        "n_pass": sum(1 for c in checks if c["pass"]),
        "n_gating_failures": len(bad),
        "gating_failures": bad,
        "non_gating_failures": [c["id"] for c in checks
                                if not c["pass"] and not c["gating"]],
        "checks": checks,
        "three_way_mass_split": ctx["mass"],
        "degradation_fragments": ctx["deg"],
        "printed_point_values": ctx["points"],
        "threshold_sensitivity": {
            "_what": "stage05 moved the black/grey colour split from 120 to "
                     "85 so the Lpf2k control bar could be measured.",
            "effect_on_the_33_library_values": "none -- the 05_values.csv "
                     "columns formulation..log10_value were byte-identical "
                     "before and after, verified by diff on 2026-11-10.",
        },
        "what_the_automated_checks_CANNOT_see": [
            "CONNECTIVITY. Every mass check above passes on any isomer. The "
            "regiochemistry is settled by P10a's enumeration (there is only "
            "one constitution) and by the synthesis the paper describes, not "
            "by any mass.",
            "STEREOCHEMISTRY, twice over. (i) The cyclic benzylidene acetal "
            "has a C2/C5 cis-trans relationship that the paper never draws and "
            "never discusses, and the stored SMILES is flat. (ii) WORSE: "
            "assemble.stereocentres(), which is what assemble.attach's "
            "check_stereo_preserved relies on, returns {} for this tail, while "
            "Chem.FindPotentialStereo returns two unspecified centres. The "
            "library's stereo guard is blind to the one stereochemical feature "
            "this paper's chemistry turns on, so its passing is not evidence "
            "here. DECISIONS.md D-2.",
            "WHETHER THE X-AXIS LABELS WERE READ CORRECTLY. The group order is "
            "read off a 1050 px raster by eye. P4/P5/P7/P8 constrain it hard "
            "-- five groups are identified by name in the prose and the "
            "remaining order is pinned by eleven sign statements -- but no "
            "check reads the glyphs.",
            "WHETHER FIGURE 5 IS THE SAME DATA AS THE PROSE. If the authors "
            "typed a number that their own figure does not show, every check "
            "here would split rather than agree, which is exactly what P6b "
            "does.",
        ],
    }
    with open(os.path.join(HERE, "07_qc.json"), "w") as f:
        json.dump(out, f, indent=1)

    print(f"stage07: {out['n_pass']}/{out['n_checks']} checks pass; "
          f"{len(bad)} GATING failure(s) {bad}; "
          f"non-gating failures {out['non_gating_failures']}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
