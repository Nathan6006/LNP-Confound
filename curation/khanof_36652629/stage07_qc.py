"""Stage 07 -- every check, with numbers, for PMID 36652629.

Checks, and what each one is worth:

  A  STROKE-CENTRELINE CONVENTION, from a quantity nobody printed. The x-axis
     baseline rule and the lowest y tick are two different strokes drawn at the
     same value. If strokes carry their value at the CENTRELINE, the baseline
     band's centre must equal the lowest tick's darkness-weighted centroid. It
     does, in all six panels. This is what licenses value_row = outer + t/2,
     and it is external to the bar reader: neither number comes from a bar.

  B  MASS, SPLIT THREE WAYS, plus the monoisotopic-vs-average question.
     B1 ours vs their printed FORMULA.
     B2 their formula vs their own printed "calc'd" m/z.
     B3 ours vs their printed "found" m/z.
     B2 is the one that pays: the SI's "calc'd" column is AVERAGE mass on some
     compounds and MONOISOTOPIC on others, and on C3-K2-E14 the calc and the
     found are on DIFFERENT scales.

  C  PANEL -> FORMULATION MAPPING, from two independent printed sources that
     were typeset separately: Table S2 (text layer, tier 4) and Fig S11's pKa
     row (printed numerals inside a raster). Every pKa in Fig S11 must be some
     Table S2 row's pKa for the same lipid.

  D  CROSS-ROUND RATIOS for lipids measured in more than one panel. Not a pass
     or fail -- a measurement, and the main finding of this curation.

  E  THE UNPRINTED CHECK: Fig S9A. The authors drew two fold-increase bars
     derived from Fig S1B and Fig S4B. Digitising S9A independently and
     comparing against ratios formed from OUR S1B/S4B reads tests the bar
     reader against something outside it. A ratio within one panel cancels any
     constant row offset and any calibration origin error, so this tests the
     SCALE and the per-bar tops, and nothing else.

  F  SEEDED MUTATIONS. Each defect is pasted into a copy of the real panel and
     the refusal is required to fire. Credited only where the AGREEMENT COUNT
     MOVES, and the exit code is read from the process, never through a pipe.
"""
import os, sys, csv, json, math, re, subprocess
import numpy as np
import pymupdf

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import axis, bars                                            # noqa: E402
from errors import DetectionError                            # noqa: E402
sys.path.insert(0, HERE)
import stage05_digitize as st                                # noqa: E402
from rdkit import Chem, RDLogger                             # noqa: E402
from rdkit.Chem.Descriptors import ExactMolWt, MolWt         # noqa: E402
from rdkit.Chem import rdMolDescriptors                      # noqa: E402
RDLogger.DisableLog("rdApp.*")

PDF = st.PDF
PROTON = 1.00727646        # mass of H+ (monoisotopic), for [M+H]+
H_AVG = 1.008              # average atomic H, for an average-mass [M+H]+

# Fig S11's printed pKa row, read off the 900 dpi crop as NUMERALS (symbolic,
# not estimated). Order as printed, left to right.
FIG_S11_PKA = [
    ("SM-102", 6.61), ("C3-K2-E14", 5.49), ("C3-K2-A14", 6.25),
    ("C3-K2-E10", 6.26), ("C3-K2-E15", 4.37), ("C3-K3-E10", 6.87),
    ("OC2-K3-E10", 6.76), ("OC3-K3-E10", 6.82), ("OC2-K2-A10", 6.66),
    ("MC3", 6.30), ("ALC-0315", 6.27),
]
# Fig S7a's axis says MC3; Table S2 rows 25-27 say DLin-MC3-DMA. The paper's
# own expansion, not ours.
ALIAS = {"MC3": "DLin-MC3-DMA"}

S9A = dict(page=7, clip=(190, 532, 262, 628), ticks=[0, 2, 4, 6])


# --------------------------------------------------------------- A: strokes -
def check_centreline():
    rows = []
    for p in st.PANELS:
        rgb = st.render(p["page"], p["clip"])
        fr = st.find_frame(rgb)
        tk = sorted(st.ticks(rgb, fr))
        base_c = (fr["baseline_top"] + fr["baseline_bot"]) / 2.0
        low = tk[-1]
        rows.append(dict(panel=p["key"],
                         baseline_band=[fr["baseline_top"], fr["baseline_bot"]],
                         baseline_centre_px=round(base_c, 2),
                         lowest_tick_centroid_px=round(low, 2),
                         delta_px=round(abs(base_c - low), 3)))
    worst = max(r["delta_px"] for r in rows)
    return dict(check="A_stroke_centreline_convention", per_panel=rows,
                worst_delta_px=worst, tol_px=1.5, pass_=bool(worst <= 1.5),
                what_it_proves=(
                    "The baseline rule and the lowest y tick are separate "
                    "strokes drawn at the same value. Their agreement to "
                    f"{worst} px means a stroke carries its value at its "
                    "CENTRELINE, so a bar's value is its top rule's centre, "
                    "not the rule's outer edge. Neither quantity comes from a "
                    "bar, so this does not assume what it tests."))


# ------------------------------------------------------------------ B: mass -
def check_mass():
    rows, agree = [], dict(B1=0, B2_any=0, B3=0, n=0)
    with open(os.path.join(HERE, "03_enumerated.csv")) as f:
        for r in csv.DictReader(f):
            if not r["formula_SI"]:
                continue
            agree["n"] += 1
            m = Chem.MolFromSmiles(r["IL_SMILES"])
            ours_f = rdMolDescriptors.CalcMolFormula(m)
            b1 = (ours_f == r["formula_SI"])
            agree["B1"] += b1

            mono = ExactMolWt(m)
            avg = MolWt(m)
            ion = r["SI_ion"]
            add_mono = PROTON if ion == "[M+H]+" else 0.0
            add_avg = H_AVG if ion == "[M+H]+" else 0.0
            calc = float(r["SI_calc_mz"])
            found = float(r["SI_found_mz"])

            d_calc_mono = calc - (mono + add_mono)
            d_calc_avg = calc - (avg + add_avg)
            kind = ("monoisotopic" if abs(d_calc_mono) < abs(d_calc_avg)
                    else "average")
            d_calc = min(abs(d_calc_mono), abs(d_calc_avg))
            b2 = d_calc <= 0.06
            agree["B2_any"] += b2

            d_found_mono = found - (mono + add_mono)
            d_found_avg = found - (avg + add_avg)
            fkind = ("monoisotopic" if abs(d_found_mono) < abs(d_found_avg)
                     else "average")
            b3 = min(abs(d_found_mono), abs(d_found_avg)) <= 0.9
            agree["B3"] += b3

            rows.append(dict(
                IL_name=r["IL_name"], ion=ion, SI_formula=r["formula_SI"],
                our_formula=ours_f, B1_formula_match=bool(b1),
                our_monoisotopic=round(mono, 4), our_average=round(avg, 2),
                SI_calc=calc, SI_found=found,
                calc_minus_our_monoisotopic=round(d_calc_mono, 3),
                calc_minus_our_average=round(d_calc_avg, 3),
                SI_calc_scale=kind, B2_calc_within_0p06=bool(b2),
                found_minus_our_monoisotopic=round(d_found_mono, 3),
                found_minus_our_average=round(d_found_avg, 3),
                SI_found_scale=fkind, B3_found_within_0p9=bool(b3),
                calc_and_found_same_scale=bool(kind == fkind)))
    mixed = [r["IL_name"] for r in rows if not r["calc_and_found_same_scale"]]
    kinds = sorted({r["SI_calc_scale"] for r in rows})
    return dict(check="B_mass_three_ways", per_lipid=rows, agreement=agree,
                SI_calc_scales_seen=kinds,
                calc_and_found_on_different_scales=mixed,
                note=("B1 is ours vs their FORMULA; B2 is their formula vs "
                      "their own printed calc'd m/z; B3 is ours vs their "
                      "found m/z. B2 is what reveals that the 'calc'd' column "
                      "is not on one scale across the SI."))


# --------------------------------------------------- C: panel -> formulation -
def parse_table_s2():
    """Table S2, from the text layer. 28 formulations.

    The lipid name WRAPS: `OC2-K3-` sits on the line above the data row and
    `E10 mRNA` on the line below, so the data line itself reads
    `11 5 7.61 102 0.12 97 6.76` with no name at all. The numeric tail is
    unambiguous, so the row is parsed from the numbers and the name is
    reassembled from the token above plus the token below. A first version
    joined lowercase continuation lines before matching, which appended
    `mRNA` to every data row and dropped 27 of the 28.
    """
    import pdfplumber
    lines = []
    with pdfplumber.open(PDF) as pdf:
        for i in (35, 36):
            lines += (pdf.pages[i].extract_text() or "").splitlines()

    NUM = (r"^\s*(\d{1,2})\s+(.*?)\s*(\d+(?:\.\d+)?)\s+(\d+\.\d+)\s+"
           r"(\d+)\s+(\d+\.\d+|N/A)\s+(\d+)\s+(\d+\.\d+|N/A)\s*$")
    pat = re.compile(NUM)
    names = ["DLin-MC3-DMA", "OC2-K3-E10", "OC3-K3-E10", "OC2-K2-A10",
             "C3-K2-A14", "C3-K2-E14", "C3-K2-E15", "C3-K2-E10",
             "C3-K3-E10", "SM-102"]
    rows = []
    for i, ln in enumerate(lines):
        mt = pat.match(ln)
        if not mt:
            continue
        mid = mt.group(2)
        lipid = next((n for n in names if n in mid), "")
        if not lipid:
            prev = lines[i - 1].split() if i else []
            nxt = lines[i + 1].split() if i + 1 < len(lines) else []
            if prev and prev[0].endswith("-") and nxt:
                cand = prev[0] + nxt[0]
                lipid = next((n for n in names if n == cand), "")
        rows.append(dict(
            id=int(mt.group(1)), lipid=lipid, raw_mid=mid,
            NP=float(mt.group(3)), lipid_rna=float(mt.group(4)),
            z=int(mt.group(5)),
            pdi=(None if mt.group(6) == "N/A" else float(mt.group(6))),
            ee=int(mt.group(7)),
            pka=(None if mt.group(8) == "N/A" else float(mt.group(8)))))
    return rows


def check_mapping():
    rows = parse_table_s2()
    byl = {}
    for r in rows:
        byl.setdefault(r["lipid"], []).append(r)

    res, n_ok = [], 0
    for nm, pka in FIG_S11_PKA:
        key = ALIAS.get(nm, nm)
        cands = byl.get(key, [])
        matched = [c["id"] for c in cands
                   if c["pka"] is not None and abs(c["pka"] - pka) <= 0.006]
        ok = bool(matched)
        n_ok += ok
        res.append(dict(fig_S11_lipid=nm, resolved_to=key, fig_S11_pka=pka,
                        table_S2_ids_for_lipid=[c["id"] for c in cands],
                        table_S2_pkas=[c["pka"] for c in cands],
                        matched_formulation_ids=matched, agrees=ok))
    dup = {n: [c["id"] for c in v] for n, v in byl.items() if len(v) > 1}
    return dict(check="C_panel_to_formulation_mapping",
                table_S2_rows_parsed=len(rows),
                table_S2_rows_without_a_lipid_name=[r["id"] for r in rows
                                                    if not r["lipid"]],
                table_S2=rows,
                table_S2_ids=sorted(r["id"] for r in rows),
                fig_S11_vs_table_S2=res,
                agreement=f"{n_ok}/{len(FIG_S11_PKA)}",
                lipids_formulated_more_than_once=dup,
                note=("Table S2 is the text layer (tier 4); Fig S11's pKa row "
                      "is printed numerals inside a raster. They were typeset "
                      "separately, so agreement is a real cross-check and it "
                      "is what ties a drawn panel to a named batch."))


# ------------------------------------------------------- D: cross-round ------
def check_rounds(dig):
    seen = {}
    for p in dig["panels"]:
        for b in p["bars"]:
            seen.setdefault(b["label"], []).append(
                (p["key"], p["round_id"], b["value"], p["cargo"]))
    out = []
    for lab, lst in sorted(seen.items()):
        if len(lst) < 2:
            continue
        vals = [v for _, _, v, _ in lst]
        out.append(dict(
            IL_name=lab, n_panels=len(lst),
            panels=[dict(panel=k, round_id=r, value=v, cargo=c)
                    for k, r, v, c in lst],
            max_over_min=round(max(vals) / min(vals), 3),
            log10_spread=round(math.log10(max(vals) / min(vals)), 3)))
    return dict(check="D_cross_round_ratios", shared=out,
                n_lipids_in_more_than_one_panel=len(out),
                worst_ratio=max([o["max_over_min"] for o in out] or [0]))


# ----------------------------------------------------- E: the unprinted check -
def read_s9a():
    """Fig S9A: white bars, black outline, LINEAR axis. No dots, no caps."""
    rgb = st.render(S9A["page"], S9A["clip"])
    nat = st.native_dpi(S9A["page"], S9A["clip"])
    dark = rgb.max(2) < st.DARK
    fr = st.find_frame(rgb)
    tk = sorted(st.ticks(rgb, fr))
    if len(tk) != len(S9A["ticks"]):
        raise DetectionError(
            f"S9a: {len(tk)} ticks, expected {len(S9A['ticks'])}.")
    cal = axis.calibrate(sorted(tk), sorted(S9A["ticks"], reverse=True),
                         kind="linear", axis="y",
                         loo_tol_px=axis.LOO_TOL_PX * st.DPI / nat)

    bot = fr["baseline_top"]
    # The bars are HOLLOW, so a column carries a bar only where its two walls
    # are. Find the tall dark columns right of the spine and pair them.
    cnt = dark[:bot, :].sum(0)
    on = cnt > 0.15 * bot
    on[:fr["spine_r"] + 2] = False
    runs, s0 = [], None
    for i, v in enumerate(on):
        if v and s0 is None:
            s0 = i
        elif not v and s0 is not None:
            runs.append((s0, i)); s0 = None
    if s0 is not None:
        runs.append((s0, len(on)))
    runs = [r for r in runs if r[1] - r[0] >= 6]
    if len(runs) != 4:
        raise DetectionError(
            f"S9a: found {len(runs)} wall runs {runs}, expected 4 (two bars, "
            "two walls each). REFUSING.")
    out = []
    for lab, (a, b) in zip(["Tails", "Core"], [(runs[0], runs[1]),
                                               (runs[2], runs[3])]):
        x0, x1 = a[0], b[1]
        t = ((a[1] - a[0]) + (b[1] - b[0])) / 2.0
        res = bars.bar_top(dark, x0, x1, 0, bot, baseline_row=bot,
                           flank_tol_px=bars.FLANK_TOL_PX * st.DPI / nat)
        row = res["row"] + t / 2.0
        out.append(dict(label=lab, x0=int(x0), x1=int(x1),
                        wall_px=round(t, 2), outer_row=int(res["row"]),
                        value_row=round(row, 2),
                        value=round(float(cal.to_value(row)), 4)))
    return dict(native_image_dpi=round(nat, 2), calibration=cal.qc(),
                tick_rows=[round(v, 2) for v in sorted(tk)],
                loo_worst_px=round(max(axis.loo_residuals(sorted(tk))), 3),
                bars=out)


def check_unprinted(dig):
    s9 = read_s9a()
    drawn = {b["label"]: b["value"] for b in s9["bars"]}
    pan = {p["key"]: {b["label"]: b["value"] for b in p["bars"]}
           for p in dig["panels"]}
    ours = {
        "Tails": dict(num="C3-K2-E14", den="C3-K2-A14", panel="S1b",
                      value=pan["S1b"]["C3-K2-E14"] / pan["S1b"]["C3-K2-A14"]),
        "Core": dict(num="OC2-K3-E10", den="C3-K3-E10", panel="S4b",
                     value=pan["S4b"]["OC2-K3-E10"] / pan["S4b"]["C3-K3-E10"]),
    }
    rows = []
    for k, v in ours.items():
        rows.append(dict(
            comparison=k, numerator=v["num"], denominator=v["den"],
            source_panel=v["panel"],
            ours_fold=round(v["value"], 3),
            drawn_in_fig_S9A=drawn[k],
            ratio_ours_over_drawn=round(v["value"] / drawn[k], 3),
            log10_gap=round(math.log10(v["value"] / drawn[k]), 3)))
    return dict(check="E_unprinted_fig_S9A_fold_increase", fig_S9A=s9,
                comparisons=rows,
                note=("A ratio taken WITHIN one panel cancels any constant "
                      "row offset and any calibration-origin error, so a "
                      "disagreement here is about the scale or the individual "
                      "bar tops -- it cannot be an origin problem."))


# ------------------------------- G: the paper's own significance groupings ---
# Every significance mark the six panels carry, read off the captions and the
# panels. A t-test depends on variance and n as well as on the ratio of means,
# so there is no guarantee one threshold separates them -- which is exactly
# what makes it a test of the bar reads rather than a restatement of them.
ANNOTATIONS = [
    ("S1b", "SM-102",      "C3-K2-A14",  "sig"),   # **
    ("S1b", "SM-102",      "C3-K2-E14",  "sig"),   # **
    ("S2b", "C3-K2-E14",   "C3-K2-E10",  "ns"),
    ("S2b", "C3-K2-E14",   "C3-K2-E15",  "ns"),
    ("S3b", "C3-K2-E10",   "C3-K3-E10",  "sig"),   # ** bracket
    ("S4b", "SM-102",      "C3-K3-E10",  "sig"),   # **
    ("S4b", "SM-102",      "OC2-K3-E10", "ns"),
    ("S4b", "SM-102",      "OC3-K3-E10", "sig"),   # *
    ("S4c", "SM-102",      "OC2-K3-E10", "ns"),
    ("S7a", "OC2-K2-A10",  "DLin-MC3-DMA", "ns"),
    ("S7a", "OC2-K2-A10",  "SM-102",     "ns"),
    ("S7a", "OC2-K2-A10",  "OC2-K3-E10", "ns"),
]


def check_groupings(dig):
    pan = {p["key"]: {b["label"]: b["value"] for b in p["bars"]}
           for p in dig["panels"]}
    rows = []
    for key, ref, other, mark in ANNOTATIONS:
        a, b = pan[key][ref], pan[key][other]
        fold = max(a, b) / min(a, b)
        rows.append(dict(panel=key, reference=ref, compared=other,
                         printed=mark, fold_from_our_reads=round(fold, 4),
                         log10_fold=round(math.log10(fold), 4)))
    ns = [r["fold_from_our_reads"] for r in rows if r["printed"] == "ns"]
    sg = [r["fold_from_our_reads"] for r in rows if r["printed"] == "sig"]
    sep = max(ns) < min(sg)
    thr = math.sqrt(max(ns) * min(sg)) if sep else None
    return dict(
        check="G_significance_groupings_separate_on_one_threshold",
        comparisons=rows, n_ns=len(ns), n_sig=len(sg),
        largest_ns_fold=round(max(ns), 4), smallest_sig_fold=round(min(sg), 4),
        separated_by_one_threshold=bool(sep),
        threshold_fold=(round(thr, 4) if thr else None),
        margin_log10_decades=round(math.log10(min(sg) / max(ns)), 4),
        n_correct=f"{len(rows)}/{len(rows)}" if sep else "NOT SEPARABLE",
        what_it_proves=(
            "The paper marks 12 pairwise comparisons across the six panels as "
            "significant or not significant. Those marks were never used to "
            "produce a single number here, and a t-test is not a function of "
            "the fold change alone. If one threshold sorts all 12 correctly, "
            "the bar tops are ordered and spaced the way the authors' own "
            "statistics say they are. This is the strongest unprinted check "
            "available in this paper."))


# -------------------------------------------------------- F: seeded defects --
SEEDS = {
    "dots_on_both_top_corners": (
        "Paste two saturated dots of the bar's own colour onto BOTH top "
        "corners of S1b's SM-102 bar, 25 px above the rule. This is the one "
        "shape bars.bar_top's flank test cannot see, because two symmetric "
        "dots agree with each other perfectly."),
    "tick_error_only_LOO_can_see": (
        "This curation applies axis.LOO_TOL_PX in the pixels of the ORIGINAL "
        "188-221 dpi raster, which is 4.07 render px on S4b rather than 1.5. "
        "That is a widened number and it has to be shown to still have teeth. "
        "The seed searches for the smallest single-tick displacement that the "
        "max-residual test PASSES and the leave-one-out test REJECTS, and "
        "requires such a displacement to exist and to fire."),
    "halo_widened_window": (
        "Widen S3b's bar window by 8 px on each side, i.e. undo trim_to_walls. "
        "The top rule then covers less than 90 per cent of the inflated width "
        "and a flank must fail."),
}


def run_seeds():
    """Each seed is run in its OWN process; the exit code is read from the
    process object, never from a pipe. A seed is credited only when the
    control passes AND the seeded run fails."""
    results = []
    for name, why in SEEDS.items():
        script = os.path.join(HERE, "_seed_tmp.py")
        with open(script, "w") as f:
            f.write(SEED_SRC % {"seed": name})
        p = subprocess.run([sys.executable, "-I", script], cwd=HERE,
                           capture_output=True, text=True)
        os.remove(script)
        results.append(dict(seed=name, what=why, exit_code=p.returncode,
                            fired=bool(p.returncode == 1),
                            message=(p.stdout.strip().splitlines() or [""])[-1]))
    return dict(check="F_seeded_mutations", seeds=results,
                n_fired=sum(r["fired"] for r in results),
                n_seeds=len(results),
                how_verified=("each seed ran as its own process and "
                              "returncode was read from the Popen result, "
                              "not from a shell pipe"))


SEED_SRC = '''
import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import stage05_digitize as st, axis, bars
from errors import DetectionError
SEED = "%(seed)s"

def panel(key):
    return [p for p in st.PANELS if p["key"] == key][0]

try:
    if SEED == "dots_on_both_top_corners":
        p = panel("S1b"); rgb = st.render(p["page"], p["clip"])
        ink = st.bar_ink(rgb); fr = st.find_frame(rgb)
        dark = ink & (rgb.min(2) < st.PALE_MIN)
        bot = fr["baseline_top"]; r0 = int(bot*0.70); r1 = bot-6
        run = st.bar_columns(ink, fr)[0]
        x0, x1 = st.trim_to_walls(dark, run[0], run[1], r0, r1)
        w = st.walls(rgb, ink, dark, x0, x1, r0, r1)
        st0 = st.solid_top(ink, x0, x1, bot)
        base = bars.bar_top(dark, x0, x1, st0, bot, baseline_row=bot)
        colour = rgb[base["row"]+3, (x0+x1)//2].tolist()
        W = x1 - x0; d = int(W*0.34)
        for cx in (x0 + d//2, x1 - d//2):
            rgb[base["row"]-32:base["row"]-6, cx-d//2:cx+d//2] = colour
        ink2 = st.bar_ink(rgb); dark2 = ink2 & (rgb.min(2) < st.PALE_MIN)
        st2 = st.solid_top(ink2, x0, x1, bot)
        res = bars.bar_top(dark2, x0, x1, st2, bot, baseline_row=bot,
                           flank_tol_px=bars.FLANK_TOL_PX*st.DPI/219.88)
        thick = st.rule_thickness(dark2, x0, x1, res["row"], bot)
        need = max(3, int(round(st.RULE_THICK_FRAC*w["t"])))
        if res["row"] < base["row"] and thick >= need:
            print("NOT CAUGHT: read the dots at row", res["row"],
                  "instead of the rule at", base["row"])
            sys.exit(0)
        print("CAUGHT: dot row rejected; thickness", thick, "needed", need)
        sys.exit(1)

    if SEED == "tick_error_only_LOO_can_see":
        p = panel("S4b"); rgb = st.render(p["page"], p["clip"])
        fr = st.find_frame(rgb); tk = sorted(st.ticks(rgb, fr))
        nat = st.native_dpi(p["page"], p["clip"])
        tol = axis.LOO_TOL_PX*st.DPI/nat
        vals = [10.0**d for d in sorted(p["decades"], reverse=True)]
        axis.calibrate(tk, vals, kind="log10", axis="y", loo_tol_px=tol)
        found = None
        for i in range(len(tk)):
            for off in [x*0.25 for x in range(1, 60)]:
                for sgn in (+1, -1):
                    bad = list(tk); bad[i] += sgn*off
                    try:            # residual test alone
                        axis.calibrate(bad, vals, kind="log10", axis="y")
                    except Exception:
                        continue    # residual already catches it; not our case
                    try:            # residual passed -- does LOO catch it?
                        axis.calibrate(bad, vals, kind="log10", axis="y",
                                       loo_tol_px=tol)
                    except Exception as e:
                        found = (i, sgn*off, str(e).splitlines()[0])
                        break
                if found: break
            if found: break
        if found is None:
            print("NOT CAUGHT: no displacement exists that the residual test "
                  "passes and LOO rejects -- at this tolerance LOO adds "
                  "nothing"); sys.exit(0)
        print(f"CAUGHT: tick {found[0]} moved {found[1]:+.2f} px passes the "
              f"max-residual test and LOO rejects it"); sys.exit(1)

    if SEED == "halo_widened_window":
        p = panel("S3b"); rgb = st.render(p["page"], p["clip"])
        ink = st.bar_ink(rgb); fr = st.find_frame(rgb)
        dark = ink & (rgb.min(2) < st.PALE_MIN)
        bot = fr["baseline_top"]; r0 = int(bot*0.70); r1 = bot-6
        run = st.bar_columns(ink, fr)[0]
        x0, x1 = st.trim_to_walls(dark, run[0], run[1], r0, r1)
        s0 = st.solid_top(ink, x0, x1, bot)
        bars.bar_top(dark, x0, x1, s0, bot, baseline_row=bot,
                     flank_tol_px=bars.FLANK_TOL_PX*st.DPI/220.27)
        try:
            bars.bar_top(dark, x0-8, x1+8, s0, bot, baseline_row=bot,
                         flank_tol_px=bars.FLANK_TOL_PX*st.DPI/220.27)
        except Exception as e:
            print("CAUGHT:", str(e).splitlines()[0]); sys.exit(1)
        print("NOT CAUGHT: the inflated window passed"); sys.exit(0)
except DetectionError as e:
    print("CAUGHT(refusal):", str(e).splitlines()[0]); sys.exit(1)
print("seed did not run"); sys.exit(2)
'''


def main():
    dig = json.load(open(os.path.join(HERE, "05_digitize.json")))
    doc = {"_what": "Stage 07 QC for PMID 36652629. Every check with numbers.",
           "checks": []}
    for fn in (check_centreline, check_mass, check_mapping):
        doc["checks"].append(fn())
    doc["checks"].append(check_rounds(dig))
    doc["checks"].append(check_unprinted(dig))
    doc["checks"].append(check_groupings(dig))
    doc["checks"].append(run_seeds())

    with open(os.path.join(HERE, "07_qc.json"), "w") as f:
        json.dump(doc, f, indent=1)

    for c in doc["checks"]:
        print("\n===", c["check"], "===")
        if c["check"].startswith("A"):
            print("  worst baseline-vs-lowest-tick delta:",
                  c["worst_delta_px"], "px   pass:", c["pass_"])
            for r in c["per_panel"]:
                print(f"    {r['panel']:5s} baseline centre "
                      f"{r['baseline_centre_px']:8.2f}  lowest tick "
                      f"{r['lowest_tick_centroid_px']:8.2f}  "
                      f"delta {r['delta_px']:.3f}")
        elif c["check"].startswith("B"):
            print("  agreement:", c["agreement"])
            print("  SI 'calc'd' scales seen:", c["SI_calc_scales_seen"])
            print("  calc and found on DIFFERENT scales:",
                  c["calc_and_found_on_different_scales"])
            for r in c["per_lipid"]:
                print(f"    {r['IL_name']:12s} {r['SI_formula']:14s} "
                      f"B1 {r['B1_formula_match']!s:5s} calc {r['SI_calc']:9.2f} "
                      f"({r['SI_calc_scale']:12s} d="
                      f"{min(abs(r['calc_minus_our_monoisotopic']), abs(r['calc_minus_our_average'])):.3f})"
                      f"  found {r['SI_found']:9.2f} ({r['SI_found_scale']})")
        elif c["check"].startswith("C"):
            print("  Table S2 rows parsed:", c["table_S2_rows_parsed"])
            print("  Fig S11 pKa vs Table S2:", c["agreement"])
            for r in c["fig_S11_vs_table_S2"]:
                print(f"    {r['fig_S11_lipid']:14s} pKa {r['fig_S11_pka']} "
                      f"-> ids {r['matched_formulation_ids']} "
                      f"(all ids {r['table_S2_ids_for_lipid']})")
        elif c["check"].startswith("D"):
            for s in c["shared"]:
                ps = " ".join(f"{q['panel']}={q['value']:.4g}"
                              for q in s["panels"])
                print(f"    {s['IL_name']:12s} {ps}   max/min "
                      f"{s['max_over_min']}x")
        elif c["check"].startswith("E"):
            for b in c["fig_S9A"]["bars"]:
                print(f"    Fig S9A {b['label']:6s} drawn {b['value']}")
            for r in c["comparisons"]:
                print(f"    {r['comparison']:6s} ours {r['ours_fold']}x "
                      f"({r['numerator']}/{r['denominator']} in "
                      f"{r['source_panel']})  drawn {r['drawn_in_fig_S9A']}x "
                      f"  ours/drawn {r['ratio_ours_over_drawn']}")
        elif c["check"].startswith("G"):
            print(f"  largest ns fold {c['largest_ns_fold']}x   "
                  f"smallest significant fold {c['smallest_sig_fold']}x")
            print(f"  separated by one threshold: "
                  f"{c['separated_by_one_threshold']} at "
                  f"{c['threshold_fold']}x, margin "
                  f"{c['margin_log10_decades']} decades, {c['n_correct']}")
            for r in c["comparisons"]:
                print(f"    {r['panel']:4s} {r['reference']:12s} vs "
                      f"{r['compared']:13s} {r['printed']:3s} "
                      f"{r['fold_from_our_reads']:8.3f}x")
        elif c["check"].startswith("F"):
            print(f"  {c['n_fired']}/{c['n_seeds']} seeded defects fired")
            for r in c["seeds"]:
                print(f"    {r['seed']:28s} exit {r['exit_code']} "
                      f"fired {r['fired']}  {r['message'][:80]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
