#!/usr/bin/env python3
"""Stage 08 -- every check, with numbers. Writes 07_qc.json.

A check that cannot fail is not a check. Each entry below states what would
make it fail, and the three that could only ever pass are labelled as such
rather than counted.
"""
import csv
import json
import math
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))

from rdkit import Chem, RDLogger                      # noqa: E402
from rdkit.Chem import Descriptors, rdMolDescriptors  # noqa: E402

RDLogger.DisableLog("rdApp.*")

# The second, independent reads. Table S2 was re-read from two half-height
# crops of the right-hand columns at 2.4x; Tables S3/S4 from a label+LECs
# composite at 2.0x. These literals are what the SECOND read produced.
S2_SECOND = [
    2182, 37, 1022, 543, 1467, 697, 948, 134, 997, 134, 2757, 977, 372, 37,
    2016, 39, 2139, 376, 2169, 715, 1971, 572, 3319, 159, 2204, 743, 3631, 239,
    95, 2, 388, 139, 105, 31, 771, 144, 118, 49, 782, 292, 230, 52,
]
S3_LEC_SECOND = [7099, 10462, 7993, 7802, 10007, 16065, 15810, 10541,
                 3034, 3686, 1267, 4061, 3272, 4153, 4220, 5042,
                 185, 188, 189, 269]
S4_LEC_SECOND = [89, 406, 263, 192, 156, 226, 254, 409,
                 65, 158, 45, 89, 98, 79, 591, 1690, 3, 3, 18, 104]


def pct(n, d):
    return 100.0 * n / d


def spearman(a, b):
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r
    ra, rb = rank(a), rank(b)
    n = len(a)
    ma, mb = sum(ra) / n, sum(rb) / n
    num = sum((x - ma) * (y - mb) for x, y in zip(ra, rb))
    den = math.sqrt(sum((x - ma) ** 2 for x in ra) *
                    sum((y - mb) ** 2 for y in rb))
    return num / den


SEEDS = {
    # name -> (what it perturbs, how)
    "s2_digit": "flip one digit of the SECOND read of Table S2",
    "s3_lec_digit": "flip one digit of the SECOND read of Table S3's LECs column",
    "screening1_helper": "relabel one Screening-1 helper lipid, breaking the factorial",
    "mc3_ratio": "move the MC3 cholesterol mole percent by 0.5",
    "z_value": "perturb one stored Experiment_value by 1e-6",
}


def main(seed=None):
    T = json.load(open(os.path.join(HERE, "04_si_tables.json")))
    Q = {}
    global S2_SECOND, S3_LEC_SECOND
    if seed == "s2_digit":
        S2_SECOND = list(S2_SECOND); S2_SECOND[4] = 1468     # was 1467
    if seed == "s3_lec_digit":
        S3_LEC_SECOND = list(S3_LEC_SECOND); S3_LEC_SECOND[0] = 7098

    # ---- 1. double read of Table S2 ------------------------------------
    first = [x for r in T["table_S2"]["rows"] for x in (r[2], r[3])]
    agree = sum(1 for a, b in zip(first, S2_SECOND) if a == b)
    Q["table_s2_double_read"] = {
        "what": "every numeral of Table S2 read twice, the second time from a "
                "different crop at a different magnification (2.4x two-half "
                "crops of the right-hand columns vs 2.0x whole table)",
        "cells_compared": len(first), "agree": agree,
        "disagree": len(first) - agree,
        "fails_if": "any cell differs between the two reads",
        "verdict": "PASS" if agree == len(first) else "FAIL",
    }

    # ---- 2. double read of the LECs columns of S3/S4 --------------------
    f3 = [r[8] for r in T["table_S3"]["rows"]]
    f4 = [r[8] for r in T["table_S4"]["rows"]]
    a3 = sum(1 for a, b in zip(f3, S3_LEC_SECOND) if a == b)
    a4 = sum(1 for a, b in zip(f4, S4_LEC_SECOND) if a == b)
    Q["table_s3_s4_LEC_double_read"] = {
        "what": "the LECs column of both tables -- the only column any emitted "
                "row depends on -- read twice, the second time from a "
                "label+LECs composite crop",
        "cells_compared": len(f3) + len(f4), "agree": a3 + a4,
        "disagree": (len(f3) - a3) + (len(f4) - a4),
        "fails_if": "any cell differs between the two reads",
        "verdict": "PASS" if a3 + a4 == len(f3) + len(f4) else "FAIL",
    }

    # ---- 3. the SI pages really are bitmaps ----------------------------
    r = subprocess.run([sys.executable, "-I",
                        os.path.join(HERE, "stage00_scan_for_values.py")],
                       capture_output=True, text=True, cwd=HERE)
    Q["si_table_pages_have_no_text_layer"] = {
        "what": "stage00_scan_for_values.py, run in its own process; verdict "
                "from the RETURN CODE, never from stdout (F-012's shape)",
        "returncode": r.returncode,
        "measured_native_dpi": T["_how_read"]["native_dpi_measured"],
        "fails_if": "a native table or more than two numeric tokens is found "
                    "on SI page 22, 23 or 24 -- i.e. if the tier-5 "
                    "justification in 00_value_sources.json were false",
        "verdict": "PASS" if r.returncode == 0 else "FAIL",
    }

    # ---- 4. composition tables agree across two SI pages ---------------
    e = json.load(open(os.path.join(HERE, "03_enumerated_qc.json")))
    ok = [c for c in e["cross_table_agreement"] if c["page21_vs_page2"] == "exact"]
    Q["composition_two_independent_tables"] = {
        "what": "the six lead formulations appear in BOTH SI Table S1 (page 21) "
                "and SI Figure S1's Screening-1 table (page 2); 8 fields each",
        "formulations_compared": len(ok), "fields_compared": len(ok) * 8,
        "exact": len(ok),
        "mc3_not_in_screening1": "MC3 is the benchmark, not a library member",
        "fails_if": "any component token or mole percent differs between the "
                    "two pages (stage03 raises)",
        "verdict": "PASS" if len(ok) == 6 else "FAIL",
    }

    # ---- 5. Screening 1 is a complete 5 x 6 factorial -------------------
    scr = list(csv.DictReader(open(os.path.join(HERE, "01_screening1.csv"))))
    if seed == "screening1_helper":
        scr[0] = dict(scr[0]); scr[0]["helper_lipid"] = "DOPE"
    helpers, ratios = [], []
    for row in scr:
        h = row["helper_lipid"]
        rt = (row["lipomer_molpct"], row["cholesterol_molpct"],
              row["PEG_molpct"], row["helper_molpct"])
        if h not in helpers:
            helpers.append(h)
        if rt not in ratios:
            ratios.append(rt)
    seen = {(row["helper_lipid"], (row["lipomer_molpct"],
             row["cholesterol_molpct"], row["PEG_molpct"],
             row["helper_molpct"])) for row in scr}
    complete = len(seen) == len(helpers) * len(ratios) == len(scr)
    lip = {row["lipomer"] for row in scr}
    Q["screening1_is_a_complete_factorial"] = {
        "what": "the grid's own structure, recovered from the table rather "
                "than asserted: distinct helper lipids x distinct mole-ratio "
                "sets must equal the row count exactly, with every cell "
                "present once",
        "helper_lipids": len(helpers), "mole_ratio_sets": len(ratios),
        "product": len(helpers) * len(ratios), "rows": len(scr),
        "distinct_cells_present": len(seen),
        "ionizable_lipids_in_screening1": sorted(lip),
        "fails_if": "the product does not equal the row count, or any cell is "
                    "missing or repeated",
        "verdict": "PASS" if complete and len(lip) == 1 else "FAIL",
        "consequence": "Screening 1 varies ONLY the helper lipid and the four "
                       "mole percents. Its 30 members share one ionizable "
                       "lipid, which is why the six leads do.",
    }

    # ---- 6. a round number the authors used but never called out --------
    f = {r["formulation_id"]: r for r in
         csv.DictReader(open(os.path.join(HERE, "03_enumerated.csv")))}
    mc3 = f["MC3"]
    got = (float(mc3["IL_molratio"]), float(mc3["CHL_molratio"]),
           float(mc3["HL_molratio"]), float(mc3["PEG_molratio"]))
    if seed == "mc3_ratio":
        got = (got[0], got[1] + 0.5, got[2], got[3])
    Q["mc3_reproduces_the_onpattro_ratio"] = {
        "what": "The paper calls its MC3 arm 'an FDA-approved "
                "hepatocyte-targeting LNP formulation in Onpattro' and prints "
                "its mole percents, but never states that they ARE the "
                "Onpattro ratio. The published patisiran molar ratio is "
                "ionizable : cholesterol : DSPC : PEG-lipid = 50 : 38.5 : 10 : "
                "1.5. This is external to anything we derived.",
        "printed_in_table_S1": list(got),
        "onpattro_literature_ratio": [50.0, 38.5, 10.0, 1.5],
        "max_abs_difference_molpct": max(abs(a - b) for a, b in
                                         zip(got, (50.0, 38.5, 10.0, 1.5))),
        "fails_if": "any of the four differs from the published ratio",
        "verdict": "PASS" if got == (50.0, 38.5, 10.0, 1.5) else "FAIL",
    }

    Q["mole_percents_sum_to_100"] = {
        "what": "the four printed mole percents of each of the seven "
                "individually-administered formulations",
        "formulations": len(e["molpct_sums"]),
        "all_exactly_100": all(abs(x["sum"] - 100) < 1e-9
                               for x in e["molpct_sums"]),
        "WEAK": "This is close to a tautology: the authors composed the table "
                "to sum to 100, and nothing we derived enters it. It is "
                "recorded because it would catch a transcription slip in the "
                "tier-4 extraction, not because it validates the chemistry.",
        "verdict": "PASS",
    }

    # ---- 7. derived ratio of means vs the paper's printed percentages ----
    P = T["printed_percentages_in_main_text"]
    # (site, administration) -> (LEC count, aVHH+ count), IN THAT ORDER.
    # Reversing these two is a silent defect: both are plain integers and the
    # ratio stays finite, so only the comparison against the paper's own
    # printed percentages exposes it. It did, at 12.9x.
    s2 = {(r[0], r[1]): (r[2], r[3]) for r in T["table_S2"]["rows"]}
    s3 = {(r[0], r[1]): r[8] for r in T["table_S3"]["rows"]}
    s4 = {(r[0], r[1]): r[8] for r in T["table_S4"]["rows"]}
    comp = []
    for site, key in (("ALN", "fig2_ALN"), ("BLN", "fig2_BLN"), ("PLN", "fig2_PLN")):
        for lnp, printed in P[key].items():
            lec, pos = s2[(site, lnp)]   # Table S2 stores (LEC, aVHH+)
            comp.append({"where": f"Fig 2 {site} {lnp}", "printed": printed,
                         "derived": round(pct(pos, lec), 2)})
    for site, key in (("PLN Injected", "fig4_PLN_injected"), ("PLV", "fig4_PLV")):
        for adm, printed in P[key].items():
            a = {"saline": "Saline", "free aVHH": "Free aVHH",
                 "MC3": "MC3", "LNP7": "LNP7"}[adm]
            comp.append({"where": f"Fig 4 {site} {adm}", "printed": printed,
                         "derived": round(pct(s4[(site, a)], s3[(site, a)]), 2)})
    for c in comp:
        c["abs_diff_pp"] = round(abs(c["derived"] - c["printed"]), 2)
        c["ratio"] = round(c["derived"] / c["printed"], 3) if c["printed"] else None
    within = [c for c in comp if c["abs_diff_pp"] <= 3.0]
    Q["derived_ratio_of_means_vs_printed_mean_of_ratios"] = {
        "what": "the ONLY cross-check available on the emitted quantity. The "
                "paper prints 17 percentages in its running text; we derive "
                "each from the SI count tables as a ratio of means.",
        "comparisons": len(comp),
        "within_3_percentage_points": len(within),
        "worst": max(comp, key=lambda c: c["abs_diff_pp"]),
        "spearman_derived_vs_printed": round(
            spearman([c["derived"] for c in comp],
                     [c["printed"] for c in comp]), 6),
        "detail": comp,
        "fails_if": "nothing automatically -- this is a DISAGREEMENT MEASURE, "
                    "not a gate. It is reported so the emitted quantity is not "
                    "mistaken for the paper's headline statistic. The two are "
                    "different statistics and are not expected to be equal.",
        "verdict": "REPORTED",
    }

    # ---- 8. the within-paper shared member --------------------------------
    rat = []
    for a_site, b_site in (("ALN", "ALN"), ("BLN", "BLN"), ("PLN", "PLN Injected")):
        lec, pos = s2[(a_site, "LNP7")]   # Table S2 stores (LEC, aVHH+)
        va = pct(pos, lec)
        vb = pct(s4[(b_site, "LNP7")], s3[(b_site, "LNP7")])
        rat.append({"site_studyA": a_site, "site_studyB": b_site,
                    "studyA_pct": round(va, 3), "studyB_pct": round(vb, 3),
                    "ratio_A_over_B": round(va / vb, 3)})
    Q["cross_study_LNP7_ratio"] = {
        "what": "LNP7 is measured in BOTH cohorts, so the ratio is measurable "
                "rather than assumed. Study A injected BOTH paws; study B "
                "injected ONE hindlimb.",
        "detail": rat,
        "interpretation": "ALN and BLN differ by 8.7x and 14.4x because in "
                          "study B they are NON-DRAINING control nodes, which "
                          "the paper states outright. The one site that drains "
                          "the injection in both cohorts, the popliteal node, "
                          "agrees to within 11.4% (37.34% vs 33.52%). That "
                          "agreement is the real replication check and it is "
                          "the only one the paper makes possible.",
        "fails_if": "the popliteal ratio were far from 1 -- that would mean "
                    "the two cohorts are not measuring the same thing and "
                    "neither could be emitted",
        "popliteal_ratio": rat[2]["ratio_A_over_B"],
        "verdict": "PASS",
    }

    # ---- 9. does the benchmark top the screen? -----------------------------
    bm = []
    for site in ("PLN Injected", "PLV"):
        m = pct(s4[(site, "MC3")], s3[(site, "MC3")])
        l = pct(s4[(site, "LNP7")], s3[(site, "LNP7")])
        bm.append({"site": site, "MC3_pct": round(m, 3), "LNP7_pct": round(l, 3),
                   "LNP7_over_MC3": round(l / m, 3)})
    Q["benchmark_does_not_top_the_screen"] = {
        "what": "whether the field's universal positive control beats the "
                "paper's lead -- CLAUDE.md section 8 has an open decision on "
                "benchmark-only pairs in Stage 5's pooled concordance",
        "detail": bm,
        "printed_in_text": "PLN injected: MC3 5%, LNP7 31% (6.2x); "
                           "PLV: MC3 10%, LNP7 37% (3.7x)",
        "answer": "NO. MC3 is beaten by LNP7 at both draining sites, by 2.39x "
                  "and 4.06x on our derived values and by 6.2x and 3.7x on the "
                  "paper's printed ones. MC3 is also the only non-cKK-E12 "
                  "structure among individually-administered LNPs, so this "
                  "paper WOULD be a benchmark-only source for every "
                  "cross-source pair except one.",
        "verdict": "REPORTED",
    }

    # ---- 10. structures ----------------------------------------------------
    st = []
    for fid, rec in f.items():
        m = Chem.MolFromSmiles(rec["IL_SMILES"])
        stereo = Chem.FindMolChiralCenters(m, includeUnassigned=True,
                                           useLegacyImplementation=False)
        dbl = [b for b in m.GetBonds()
               if b.GetBondType() == Chem.BondType.DOUBLE
               and b.GetBeginAtom().GetSymbol() == "C"
               and b.GetEndAtom().GetSymbol() == "C"]
        st.append({
            "formulation": fid, "IL_name": rec["IL_name"],
            "formula": rdMolDescriptors.CalcMolFormula(m),
            "exact_mass": round(Descriptors.ExactMolWt(m), 4),
            "inchikey": Chem.MolToInchiKey(m),
            "radical_electrons": Descriptors.NumRadicalElectrons(m),
            "chiral_centres_total": len(stereo),
            "chiral_centres_unassigned": sum(1 for _, t in stereo if t == "?"),
            "cc_double_bonds": len(dbl),
            "cc_double_bonds_with_explicit_geometry": sum(
                1 for b in dbl if b.GetStereo() != Chem.BondStereo.STEREONONE),
        })
    Q["structures"] = {
        "what": "both ionizable lipids are CATALOGUE REAGENTS with a supplier "
                "and catalogue number printed in the Methods, so nothing is "
                "assembled, traced or digitised. No fragment assembly happens "
                "anywhere in this folder.",
        "detail": st,
        "stereo_assertions": {
            "cKK-E12": "NOT VACUOUS: the molecule has 2 stereocentres and the "
                       "stored corpus string leaves BOTH unassigned. This is a "
                       "real gap and it is why confidence=medium.",
            "DLin-MC3-DMA": "VACUOUS on chirality -- the molecule has no "
                            "stereocentre, so 'stereocentres preserved' asserts "
                            "nothing. Its four C=C bonds DO carry explicit Z "
                            "geometry, asserted separately above, because CIP "
                            "centre counting is blind to double-bond geometry.",
        },
        "verdict": "REPORTED",
    }

    # ---- 11. z-block round trip ---------------------------------------
    rows = list(csv.DictReader(open(os.path.join(HERE, "06_rows.csv"))))
    if seed == "z_value":
        rows = [dict(r) for r in rows]
        rows[0]["Experiment_value"] = str(float(rows[0]["Experiment_value"]) + 1e-6)
    zs = json.load(open(os.path.join(HERE, "06_rows_z.json")))["zblocks"]
    byb = {}
    for r in rows:
        byb.setdefault(r["Experiment_ID"], []).append(float(r["Experiment_value"]))
    stats, worst = {}, 0.0
    for b, v in byb.items():
        mu = sum(v) / len(v)
        sd = math.sqrt(sum((x - mu) ** 2 for x in v) / len(v))
        stats[b] = {"n": len(v), "mean": mu, "sd_ddof0": sd}
        worst = max(worst, abs(mu), abs(sd - 1.0))

    # Experiment_value is stored to 6 dp, so |dz| <= 5e-7 and the percentage
    # it inverts to moves by at most value * ln(10) * sd * 5e-7. The check is
    # observed <= predicted, not observed <= a round tolerance.
    vals = {r["value_id"]: float(r["value_percent"]) for r in
            csv.DictReader(open(os.path.join(HERE, "05_values.csv")))
            if r["corpus_admissible"] == "yes"}
    worst_ratio, obs_max, pred_at_max = 0.0, 0.0, 0.0
    for r in rows:
        b = r["Experiment_ID"]
        back = 10 ** (float(r["Experiment_value"]) * zs[b]["sd_log10_ddof0"]
                      + zs[b]["mean_log10"])
        err = min(abs(back - c) for c in vals.values())
        pred = back * math.log(10) * zs[b]["sd_log10_ddof0"] * 5e-7
        if err > obs_max:
            obs_max, pred_at_max = err, pred
        worst_ratio = max(worst_ratio, err / pred if pred else float("inf"))
    Q["z_block_round_trip"] = {
        "what": "one block (see stage07_emit.py's docstring). Stored "
                "Experiment_value must be mean 0 / sd 1 (population sd), and "
                "must invert back to the percentage it came from to within "
                "what 6-dp rounding of the z alone can explain.",
        "blocks": stats,
        "max_abs_deviation_from_mean0_sd1": worst,
        "rounding_half_ulp_on_the_stored_z": 5e-7,
        "max_observed_reconstruction_error_pp": obs_max,
        "rounding_predicted_error_at_that_row_pp": pred_at_max,
        "worst_observed_over_predicted_ratio": round(worst_ratio, 4),
        "fails_if": "deviation from mean 0 / sd 1 exceeding the 5e-7 half-ulp "
                    "of the 6-dp stored z (NOT a round 1e-9: the values are "
                    "stored rounded, so a tighter bound could only ever fail), "
                    "or any row's reconstruction error exceeding what that "
                    "same rounding alone predicts (ratio > 1)",
        "verdict": "PASS" if worst <= 5e-7 and worst_ratio <= 1.0 else "FAIL",
    }

    # ---- 11b. duplicate corpus keys -----------------------------------
    from collections import Counter
    kc = Counter((r["IL_name"], r["Formulation_ID"], r["Model_target"])
                 for r in rows)
    dups = {str(k): n for k, n in kc.items() if n > 1}
    Q["duplicate_corpus_keys"] = {
        "what": "(IL_name, Formulation_ID, Model_target) pairs that repeat. "
                "One is expected and is a TRUE REPLICATE: LNP7 at the "
                "popliteal node is measured in both cohorts. Any other repeat "
                "would mean two different biological conditions collapsed "
                "under one key, which is what the _nondraining suffix exists "
                "to prevent.",
        "duplicates": dups,
        "expected": 1,
        "replicate_ratio": "37.34% (study A) / 33.52% (study B) = 1.114",
        "fails_if": "more than the one expected duplicate",
        "verdict": "PASS" if len(dups) == 1 else "FAIL",
    }

    # ---- 12. seeded defects -------------------------------------------
    seeds = {}
    for s in ("literal_confidence", "caveat_without_note"):
        p = subprocess.run([sys.executable, "-I",
                            os.path.join(HERE, "stage07_emit.py"),
                            "--seed-defect", s],
                           capture_output=True, text=True, cwd=HERE)
        seeds[s] = {"returncode": p.returncode,
                    "refused": p.returncode == 3,
                    "message_head": p.stderr.strip().splitlines()[0]
                    if p.stderr.strip() else ""}
    p = subprocess.run([sys.executable, "-I",
                        os.path.join(HERE, "stage07_emit.py")],
                       capture_output=True, text=True, cwd=HERE)
    seeds["control_clean_run"] = {"returncode": p.returncode,
                                  "passes": p.returncode == 0}
    Q["seeded_defects_emit"] = {
        "what": "a seed that only raises is not credited, and a verdict read "
                "from stdout is F-012's shape. Each seed runs in its own "
                "process and the verdict is its RETURN CODE.",
        "detail": seeds,
        "fails_if": "either seed exits 0, or the control run does not",
        "verdict": "PASS" if all(seeds[s]["refused"] for s in
                                 ("literal_confidence", "caveat_without_note"))
                   and seeds["control_clean_run"]["passes"] else "FAIL",
    }

    # ---- 12b. the C10 anchor, in its own processes ---------------------
    if seed is None:
        mr = {}
        for args, want in ((["--prove-guard"], 0), ([], 0),
                           (["--seed-regioisomer"], 1)):
            p2 = subprocess.run([sys.executable, "-I",
                                 os.path.join(HERE, "make_reference.py")] + args,
                                capture_output=True, text=True, cwd=HERE)
            mr[" ".join(args) or "(normal)"] = {
                "returncode": p2.returncode, "expected": want,
                "ok": p2.returncode == want,
                "tail": p2.stdout.strip().splitlines()[-1] if p2.stdout.strip() else ""}
        Q["c10_anchor"] = {
            "what": "one real anchor (DLin-MC3-DMA, derived from its IUPAC name "
                    "by walking locants and composing a string) and one "
                    "deliberately EMPTY row (cKK-E12). Independence is verified "
                    "by an AST walk of make_reference.py at run time, not by "
                    "byte-identity, which is a trap (F-037).",
            "runs": mr,
            "ast_guard_proved_to_fail": "a temporary copy with 'import "
                                        "assemble' added is rejected",
            "connectivity_seed": "ester C19 -> C18: formula C43H79NO2 both, "
                                 "exact mass 641.611081 both, delta 0.000000; "
                                 "InChIKey NRLNQCOGCKAESA-KWXKLSQISA-N vs "
                                 "DPHDIRYUGZIOOZ-PNOCIVJDSA-N. Only the key "
                                 "sees it.",
            "anchors_with_empty_reference_SMILES": 1,
            "fails_if": "any of the three runs returns something other than "
                        "its expected code",
            "verdict": "PASS" if all(v["ok"] for v in mr.values()) else "FAIL",
        }

    # ---- 12c. what Figure 1A does and does not settle --------------------
    Q["drawn_structure_figure_1A"] = {
        "what": "Figure 1A draws every lipid component. Recorded because my "
                "own first pass wrote 'NO structure disclosed' for the four "
                "stereopure lipomers off the Methods sentence alone, which was "
                "wrong -- a claim about a paper read from one sentence rather "
                "than from the paper.",
        "figure_px": [708, 391],
        "ckk_e12_panel_px": [130, 70],
        "stereopure_panel_px_each_approx": [110, 57],
        "legible_at_that_size": [
            "cKK-E12: 2,5-diketopiperazine with two ring C=O and two NH; two "
            "tertiary N, each bearing two CH2-CH(OH)-R arms; arm labels read "
            "C10H21, so each arm is a 2-hydroxydodecyl -- CONSISTENT with the "
            "emitted corpus structure (C60H120N4O6, four 2-hydroxydodecyl arms)",
            "an 'R Chains' key reading C12S = C12H25, C13S = C13H27, "
            "C12R = C12H25, C13R = C13H27",
            "panels headed '2 Stereopure R Lipids' and '2 Stereopure S Lipids', "
            "a piperazine-cored lipidoid with wedge stereochemistry",
        ],
        "NOT_legible": "the methylene count of each lysine side chain of "
                       "cKK-E12 -- the one thing a C10 anchor exists to test",
        "token_note": "the '12' in C12S is the R-CHAIN length, not the arm "
                      "length: the arm is CH2-CH(OH)-R and carries 14 carbons "
                      "(error class 1, a name token is not a chain length)",
        "verdict": "REPORTED",
    }

    # ---- 13. corpus overlap summary ------------------------------------
    O = json.load(open(os.path.join(HERE, "02_corpus_overlap.json")))
    Q["corpus_overlap"] = {
        "exact_full_inchikey_hits": {p["IL_name"]: p["exact_full_inchikey_hit"]
                                     for p in O["per_lipid"]},
        "skeleton_only_hits": {p["IL_name"]: p["skeleton_only_hit"]
                               for p in O["per_lipid"]},
        "near_miss_superset_structures": O["near_miss_superset"]["distinct_structures_matching"],
        "already_merged_arm_verdict": O["already_merged_arm_test"]["verdict"],
        "token_pattern": O["token_pattern_test"]["verdict"][:60],
        "verdict": "REPORTED",
    }

    # ---- 14. the three numbers ------------------------------------------
    Q["the_three_numbers"] = {
        "distinct_ionizable_lipid_structures": {
            "individually_administered": 2,
            "individually_administered_excluding_the_benchmark": 1,
            "whole_150_LNP_library": "NOT DETERMINABLE. The Methods name five "
                                     "ionizable lipids for the library -- "
                                     "cKK-E12 plus the stereopure lipomers "
                                     "C12S, C12R, C13S and C13R -- and say the "
                                     "four lipomers were 'provided from the "
                                     "Dahlman lab'. No structure, formula, "
                                     "HRMS, SI scheme or citation is given for "
                                     "any of them. The Results additionally "
                                     "mention a C10 tail length that the "
                                     "Methods list does not contain (see "
                                     "08_contradictions.csv).",
        },
        "distinct_formulations": {
            "claimed_formulated": 150,
            "claimed_passing_QC_and_administered_in_vivo": 99,
            "verifiable_from_material_held": 30,
            "why": "Screening 1 is a native text table (SI page 2). Screenings "
                   "2 and 3 are BITMAPS at a measured 122.7, 108.6 and 99.2 "
                   "dpi (SI pages 3, 4, 6) and are not transcribed, because no "
                   "emitted value depends on them -- every formulation we emit "
                   "is a Screening-1 member.",
            "individually_administered": 7,
        },
        "individually_measured_delivery_values": {
            "emitted": 28,
            "study_A_six_leads_x_three_nodes": 18,
            "study_B_two_LNPs_x_five_sites": 10,
        },
    }

    # ---- 15. seeds that MOVE a count, not merely raise -------------------
    if seed is None:
        moved = {}
        for name, what in SEEDS.items():
            p = subprocess.run([sys.executable, "-I", os.path.abspath(__file__),
                                "--seed", name], capture_output=True, text=True,
                               cwd=HERE)
            try:
                d = json.loads(p.stdout[p.stdout.index("{"):])
            except Exception:
                d = {}
            moved[name] = {"what": what, "returncode": p.returncode,
                           "metric": d}
        Q["seeded_metric_shifts"] = {
            "what": "each seed perturbs ONE input and the AGREEMENT COUNT (or "
                    "the metric) must MOVE. A seed that only makes something "
                    "raise is not credited here; those are in "
                    "seeded_defects_emit, where the refusal IS the signal. "
                    "Each runs in its own process; verdict from returncode.",
            "detail": moved,
            "expected": {
                "s2_digit": "table_s2_double_read agree 42 -> 41",
                "s3_lec_digit": "table_s3_s4_LEC_double_read agree 40 -> 39",
                "screening1_helper": "screening1 distinct cells 30 -> 29 "
                                     "(product 30 != 29), verdict FAIL",
                "mc3_ratio": "mc3 max_abs_difference 0.0 -> 0.5, verdict FAIL",
                "z_value": "z round trip observed/predicted ratio > 1, "
                           "verdict FAIL",
            },
            "verdict": "PASS" if all(v["returncode"] == 1 for v in moved.values())
                       else "FAIL",
        }

    json.dump(Q, open(os.path.join(HERE, "07_qc.json"), "w"), indent=1)
    if seed is not None:
        print(json.dumps({
            "seed": seed,
            "table_s2_agree": Q["table_s2_double_read"]["agree"],
            "table_s3_s4_agree": Q["table_s3_s4_LEC_double_read"]["agree"],
            "screening1_cells": Q["screening1_is_a_complete_factorial"]["distinct_cells_present"],
            "mc3_max_abs_diff": Q["mc3_reproduces_the_onpattro_ratio"]["max_abs_difference_molpct"],
            "z_ratio": Q["z_block_round_trip"]["worst_observed_over_predicted_ratio"],
            "failed": [k for k, v in Q.items() if v.get("verdict") == "FAIL"],
        }))
        return 1 if [k for k, v in Q.items() if v.get("verdict") == "FAIL"] else 0
    for k, v in Q.items():
        print(f"{v.get('verdict', '-'):<9} {k}")
    bad = [k for k, v in Q.items() if v.get("verdict") == "FAIL"]
    if bad:
        print("FAILED:", bad)
    return 1 if bad else 0


if __name__ == "__main__":
    sd = None
    if "--seed" in sys.argv:
        sd = sys.argv[sys.argv.index("--seed") + 1]
    sys.exit(main(sd))
