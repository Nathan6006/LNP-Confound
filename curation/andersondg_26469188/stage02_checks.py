#!/usr/bin/env python3
"""KK_2015 -- paper-internal checks. Writes _checks.json for 07_qc.json.

Every check reports NUMBERS. Where a check is a tautology of our own
derivation it says so and is not credited (trap 7).
"""
import json
import math
import os
import statistics
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.environ.get("KK2015_OUTDIR", HERE)
P = json.load(open(os.path.join(OUTDIR, "_parsed.json")))

S1 = {r["code"]: r for r in P["table_s1"]}
DES = P["designs"]
S2 = P["table_s2"]
S3 = {r["code"]: r for r in P["table_s3"]}
LIB = {c: r for c, r in S1.items() if c not in ("PBS", "Original")}


def libof(code):
    return code[0]


# --------------------------------------------------------------- CHECK 1 ----
def check_molpct_sums():
    """Four printed molar percentages must sum to 100.

    NOT CREDITED as a design check: cholesterol is the balance, so the identity
    holds by the AUTHORS' construction. It is credited only as a check on OUR
    TRANSCRIPTION of four independently printed numbers -- a misread digit in
    any of the four breaks it.
    """
    res, worst = [], 0.0
    for c, r in S1.items():
        if c == "PBS":
            continue
        tot = r["il_molpct"] + r["ph_molpct"] + r["chol_molpct"] + r["peg_molpct"]
        d = abs(tot - 100.0)
        worst = max(worst, d)
        res.append((c, round(tot, 10), d < 1e-9))
    for t in (P["table_s3"], [x for x in P["table_s4"] if x["code"] != "PBS"]):
        for r in t:
            tot = r["il_molpct"] + r["ph_molpct"] + r["chol_molpct"] + r["peg_molpct"]
            worst = max(worst, abs(tot - 100.0))
            res.append((r["code"], round(tot, 10), abs(tot - 100.0) < 1e-9))
    return {"name": "mol% of the four components sums to 100",
            "n_formulations_tested": len(res),
            "n_pass": sum(1 for x in res if x[2]),
            "max_abs_deviation_pct": worst,
            "credited_as": "transcription check only",
            "NOT_credited_as": "a design check -- cholesterol is printed as the "
                               "balance to 100, so the identity is a tautology of "
                               "the authors' own derivation (trap 7)",
            "offenders": [x[0] for x in res if not x[2]]}


# --------------------------------------------------------------- CHECK 2 ----
def check_cross_parameterisation():
    """Table S1 vs the two design tables on pp10-11: the SAME formulations
    printed twice, in two different places, by two different table layouts.
    This is the strong check -- it is NOT a tautology of anything we derived.
    Cholesterol is absent from the design tables, so it is NOT tested here.
    """
    fields = ["wt_ratio", "il_molpct", "ph_molpct", "peg_molpct"]
    cells = 0
    mism = []
    for code, d in sorted(DES.items()):
        s = S1[code]
        for f in fields:
            cells += 1
            if abs(s[f] - d[f]) > 1e-9:
                mism.append((code, f, s[f], d[f]))
        cells += 1
        if s["phospholipid"] != d["phospholipid"]:
            mism.append((code, "phospholipid", s["phospholipid"], d["phospholipid"]))
    return {"name": "row-wise identity, Table S1 vs the DSD/FFD design tables",
            "n_formulations": len(DES),
            "n_cells_compared": cells,
            "n_exact": cells - len(mism),
            "mismatches": mism,
            "note": "5 fields x 32 formulations. Cholesterol mol%% is NOT in the "
                    "design tables, so it has no independent attestation anywhere "
                    "in the SI -- stated rather than quietly skipped.",
            "independent": True}


# --------------------------------------------------------------- CHECK 3 ----
def check_coded_levels():
    """The coded design strings (+/0/-) are a THIRD parameterisation.

    For every factor position the character -> level map must be a well-defined
    bijection across all rows of that library, and for a quantitative factor it
    must be order-preserving with '-' < '0' < '+'. A transposed column or a
    swapped row breaks the bijection.
    """
    out = []
    specs = [("A", 6, ["head", "tail", "wt_ratio", "il_molpct", "ph_molpct", "peg_molpct"]),
             ("B", 5, ["phospholipid", "wt_ratio", "il_molpct", "ph_molpct", "peg_molpct"])]
    total = ok = 0
    for lib, width, factors in specs:
        rows = {c: d for c, d in DES.items() if d["library"] == lib}
        for r in rows.values():
            assert len(r["design"]) == width, (r, width)
        for pos, fac in enumerate(factors):
            m = {}
            bad = []
            for code, d in sorted(rows.items()):
                ch, lvl = d["design"][pos], d[fac]
                total += 1
                if ch in m and m[ch] != lvl:
                    bad.append((code, ch, lvl, m[ch]))
                else:
                    m.setdefault(ch, lvl)
                    ok += 1
            quant = fac not in ("head", "tail", "phospholipid")
            mono = None
            if quant and set(m) >= {"-", "+"}:
                seq = [m[c] for c in "-0+" if c in m]
                mono = seq == sorted(seq)
            out.append({"library": lib, "position": pos, "factor": fac,
                        "char_to_level": {k: m[k] for k in sorted(m)},
                        "bijection": not bad, "violations": bad,
                        "quantitative": quant, "order_preserving": mono})
    return {"name": "coded design string (+/0/-) is a consistent third parameterisation",
            "n_character_assertions": total,
            "n_consistent": ok,
            "n_positions": len(out),
            "all_bijective": all(o["bijection"] for o in out),
            "all_quantitative_order_preserving":
                all(o["order_preserving"] for o in out if o["quantitative"]),
            "positions": out,
            "independent": True}


# --------------------------------------------------------------- CHECK 4 ----
def check_table_s2_statistics():
    """Printed regression statistics (Table S2) recomputed from the Table S1
    Library B values. External to our extraction in the strong sense: these are
    numbers the authors printed and never meant us to recompute."""
    b = {c: r for c, r in LIB.items() if libof(c) == "B"}
    vals = [r["epo_mean"] for r in b.values()]
    grand = statistics.fmean(vals)
    dope = [r["epo_mean"] for r in b.values() if r["phospholipid"] == "DOPE"]
    dspc = [r["epo_mean"] for r in b.values() if r["phospholipid"] == "DSPC"]
    half = (statistics.fmean(dope) - statistics.fmean(dspc)) / 2.0
    return {"name": "Table S2 printed statistics recomputed from Table S1",
            "n_library_B_formulations": len(b),
            "printed_n": S2["n"],
            "n_times_3": len(b) * 3,
            "n_matches_3_replicates": len(b) * 3 == S2["n"],
            "printed_Mean": S2["mean"],
            "recomputed_mean_of_18_group_means": round(grand, 4),
            "mean_abs_error": round(abs(grand - S2["mean"]), 4),
            "mean_rel_error_pct": round(100 * abs(grand - S2["mean"]) / S2["mean"], 4),
            "printed_Phospholipid[DOPE]_estimate": S2["dope_estimate"],
            "n_DOPE": len(dope), "n_DSPC": len(dspc),
            "design_balanced_on_phospholipid": len(dope) == len(dspc),
            "mean_DOPE": round(statistics.fmean(dope), 1),
            "mean_DSPC": round(statistics.fmean(dspc), 1),
            "half_difference": round(half, 1),
            "half_diff_vs_printed_estimate_abs": round(abs(half - S2["dope_estimate"]), 1),
            "half_diff_vs_printed_estimate_rel_pct":
                round(100 * abs(half - S2["dope_estimate"]) / S2["dope_estimate"], 2),
            "interpretation": "For a 2-level nominal factor coded +/-1 in a balanced "
                              "orthogonal design the least-squares estimate IS half the "
                              "group difference. The residual gap is what the other four "
                              "terms in the reduced model absorb, so agreement is expected "
                              "to be close but not exact.",
            "independent": True}


# --------------------------------------------------------------- CHECK 5 ----
def check_abstract_fold():
    """The abstract's headline '7-fold' claim, against Table S1."""
    orig = S1["Original"]["epo_mean"]
    best = max(LIB.values(), key=lambda r: r["epo_mean"] or -1)
    c35 = S1["C-35"]["epo_mean"]
    return {"name": "abstract's 7-fold potency claim",
            "Original_EPO_ng_per_mL": orig,
            "C-35_EPO_ng_per_mL": c35,
            "C-35_over_Original": round(c35 / orig, 3),
            "max_of_screen_code": best["code"] if "code" in best else None,
            "max_EPO_ng_per_mL": best["epo_mean"],
            "max_over_Original": round(best["epo_mean"] / orig, 3),
            "claim_reproduced_on_C-35": abs(c35 / orig - 7.0) < 0.5,
            "note": "7.34x on C-35, the formulation the paper calls optimized and "
                    "carries into Tables S3/S4. The screen MAXIMUM is C-36 at 7.85x, "
                    "so '7-fold' is the claim about C-35 specifically, not about the "
                    "best formulation measured."}


# --------------------------------------------------------------- CHECK 6 ----
def check_library_c_and_B26():
    """Library C holds molar composition FIXED and varies only the weight
    ratio. B-26 (Library B, a different design and a different batch) has
    IDENTICAL molar composition at weight ratio 12.5, which Library C does not
    sample -- so B-26 is an INDEPENDENT interpolation test of the C ladder.
    """
    cs = sorted((r["wt_ratio"], c, r["epo_mean"])
                for c, r in LIB.items() if libof(c) == "C")
    comp = [(r["il_molpct"], r["ph_molpct"], r["chol_molpct"], r["peg_molpct"],
             r["phospholipid"]) for c, r in LIB.items() if libof(c) == "C"]
    const = len(set(comp)) == 1
    b26 = LIB["B-26"]
    b26comp = (b26["il_molpct"], b26["ph_molpct"], b26["chol_molpct"],
               b26["peg_molpct"], b26["phospholipid"])
    same = b26comp == comp[0]
    lo = [x for x in cs if x[0] < b26["wt_ratio"]][-1]
    hi = [x for x in cs if x[0] > b26["wt_ratio"]][0]
    bracketed = lo[2] <= b26["epo_mean"] <= hi[2]
    return {"name": "Library C weight-ratio ladder, and the independent B-26 interpolation",
            "library_C_composition_constant": const,
            "library_C_composition": list(comp[0]),
            "n_cells_asserted_constant": 5 * len(comp),
            "library_C_ladder_wt_ratio_to_EPO": [[x[0], x[1], x[2]] for x in cs],
            "ladder_strictly_increasing": [x[2] for x in cs] == sorted(x[2] for x in cs),
            "ladder_peak_at_wt_ratio": max(cs, key=lambda x: x[2])[0],
            "B-26_composition_identical_to_library_C": same,
            "B-26_wt_ratio": b26["wt_ratio"],
            "B-26_EPO": b26["epo_mean"],
            "bracketing_neighbours": [[lo[0], lo[1], lo[2]], [hi[0], hi[1], hi[2]]],
            "B-26_falls_between_them": bracketed,
            "why_this_is_strong": "B-26 was run in a different library, from a "
                                  "different DOE, as a different batch, and nothing "
                                  "in the paper presents it as a replicate. Its value "
                                  "landing inside the bracket its Library C neighbours "
                                  "draw is an independent corroboration of both.",
            "independent": True}


# --------------------------------------------------------------- CHECK 7 ----
def check_s3_s4_formulation_identity():
    """Original and C-35 are re-specified in Tables S3 and S4. Row-wise identity."""
    fields = ["wt_ratio", "il_molpct", "ph_molpct", "chol_molpct", "peg_molpct"]
    cells = 0
    mism = []
    s4 = {}
    for r in P["table_s4"]:
        if r["code"] != "PBS":
            s4.setdefault(r["code"], r)
    for tname, tbl in (("S3", S3), ("S4", s4)):
        for code, r in tbl.items():
            base = S1[code]
            for f in fields:
                cells += 1
                if abs(r[f] - base[f]) > 1e-9:
                    mism.append((tname, code, f, r[f], base[f]))
            cells += 1
            if r["phospholipid"] != base["phospholipid"]:
                mism.append((tname, code, "phospholipid", r["phospholipid"],
                             base["phospholipid"]))
    return {"name": "Original and C-35 re-specified identically in Tables S3 and S4",
            "n_cells_compared": cells, "n_exact": cells - len(mism),
            "mismatches": mism, "independent": True}


# --------------------------------------------------------------- CHECK 8 ----
def check_batch_spread():
    """THE SAME NOMINAL FORMULATION, three independent batches (EPO-mRNA in
    Table S1, Luc-mRNA in S3, FVII-siRNA in S4). Trap 4: when two sources share
    a member, measure the ratio. Here the shared members are the formulations
    and the shared quantities are EE / size / PDI."""
    s4 = {}
    for r in P["table_s4"]:
        if r["code"] != "PBS":
            s4.setdefault(r["code"], r)
    out = {}
    for code in ("Original", "C-35"):
        b = {"S1_EPO_mRNA": S1[code], "S3_Luc_mRNA": S3[code], "S4_FVII_siRNA": s4[code]}
        rec = {}
        for q in ("ee_pct", "size_nm", "pdi"):
            vals = {k: v[q] for k, v in b.items()}
            rec[q] = {"values": vals,
                      "min": min(vals.values()), "max": max(vals.values()),
                      "max_over_min": round(max(vals.values()) / min(vals.values()), 3)}
        out[code] = rec
    return {"name": "batch-to-batch spread of the SAME nominal formulation across 3 tables",
            "per_formulation": out,
            "headline": "'Original' particle size is 152 / 113 / 65 nm across the three "
                        "batches, a 2.34x spread, and its encapsulation efficiency is "
                        "24 / 27 / 59 %, a 2.46x spread -- same lab, same nominal "
                        "recipe, three cargoes. Reported because it is a measured "
                        "within-lab reproducibility figure, not because the paper "
                        "presents it as one.",
            "independent": True}


# --------------------------------------------------------------- CHECK 9 ----
def check_ranks():
    """Does the benchmark top the screen? Does the declared optimum?"""
    scored = sorted(((r["epo_mean"], c) for c, r in S1.items()
                     if c != "PBS" and r["epo_mean"] is not None), reverse=True)
    rank = {c: i + 1 for i, (_, c) in enumerate(scored)}
    n = len(scored)
    return {"name": "rank of the benchmark and of the declared optimum",
            "n_ranked_LNPs": n,
            "benchmark_is": "'Original' -- the C12-200 formulation previously "
                            "optimized for siRNA (50/10/38.5/1.5 DSPC, 5:1 wt), i.e. "
                            "the field-standard starting point this paper sets out to beat",
            "Original_rank": rank["Original"], "Original_EPO": S1["Original"]["epo_mean"],
            "benchmark_tops_the_screen": rank["Original"] == 1,
            "declared_optimum": "C-35",
            "C-35_rank": rank["C-35"], "C-35_EPO": S1["C-35"]["epo_mean"],
            "declared_optimum_tops_the_screen": rank["C-35"] == 1,
            "top5": [[c, v, rank[c]] for v, c in scored[:5]],
            "beat_the_benchmark": sum(1 for v, c in scored
                                      if v > S1["Original"]["epo_mean"]),
            "note": "The declared 'optimized formulation' C-35 is NOT the maximum of "
                    "the screen: C-36 (7548) and B-26 (7485) both exceed it. C-35 was "
                    "selected as the Library B regression optimum and Library C was the "
                    "confirmation ladder run around it, so this is a model optimum "
                    "reported as the optimized formulation, not a selection-into-a-"
                    "second-screen tautology."}


# -------------------------------------------------------------- CHECK 10 ----
def check_phospholipid_effect():
    """The paper's central structural finding: DOPE beats DSPC. Measured."""
    out = {}
    for lib in ("A", "B"):
        rows = [r for c, r in LIB.items()
                if libof(c) == lib and r["epo_mean"] is not None]
        by = {}
        for r in rows:
            by.setdefault(r["phospholipid"], []).append(r["epo_mean"])
        out[lib] = {k: {"n": len(v), "mean": round(statistics.fmean(v), 1),
                        "geomean": round(math.exp(statistics.fmean(
                            [math.log(x) for x in v])), 1)}
                    for k, v in sorted(by.items())}
    return {"name": "phospholipid identity effect, measured per library",
            "per_library": out,
            "note": "The only chemical identity that varies in this paper, and it is "
                    "a HELPER lipid, not the ionizable lipid."}


def main():
    checks = {
        "C1_molpct_sums": check_molpct_sums(),
        "C2_cross_parameterisation": check_cross_parameterisation(),
        "C3_coded_levels": check_coded_levels(),
        "C4_table_s2_statistics": check_table_s2_statistics(),
        "C5_abstract_fold_claim": check_abstract_fold(),
        "C6_library_c_and_B26": check_library_c_and_B26(),
        "C7_s3_s4_identity": check_s3_s4_formulation_identity(),
        "C8_batch_spread": check_batch_spread(),
        "C9_ranks": check_ranks(),
        "C10_phospholipid_effect": check_phospholipid_effect(),
    }
    with open(os.path.join(OUTDIR, "_checks.json"), "w") as fh:
        json.dump(checks, fh, indent=1, sort_keys=True, default=str)
    for k, v in checks.items():
        print("=" * 72)
        print(k, "--", v["name"])
        for kk, vv in v.items():
            if kk == "name":
                continue
            s = json.dumps(vv, default=str)
            print("   %-46s %s" % (kk, s if len(s) < 300 else s[:300] + " ..."))
    return 0


if __name__ == "__main__":
    sys.exit(main())
