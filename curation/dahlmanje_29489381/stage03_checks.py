#!/usr/bin/env python3
"""Checks on the two tier-4 extractions, each reported as a NUMBER.

Every check here can be made to fail on purpose by stage04_seed.py, and the
seed asserts that the check's own METRIC MOVES -- not merely that something
raised.  A seed that only proves an exception can be thrown proves nothing
about what the check measures.
"""
import csv
import json
import os
import re
import sys
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))

STATED = {"1": (144, 112), "2": (120, 105), "3": (156, 64), "4": (144, 85)}
# (formulated, stable)  -- main text P9/P11/P15, SI Fig S3h/i, Fig S5a/b


def norm(s):
    """Compare compositions across two tables printed with different spacing.

    'PEI 600C15 Epoxy' (Fig S4e) and 'PEI 600' + 'C15 Epoxy' (Fig S2a) are the
    same composition; 'C12-Epoxy' and 'C12 Epoxy' are the same lipid.
    """
    return re.sub(r"[^a-z0-9]", "", s.lower())


def load(path):
    with open(os.path.join(HERE, path), newline="") as fh:
        return list(csv.DictReader(fh))


def g(r, *names):
    for n in names:
        v = r.get(n, "")
        if v not in ("", None):
            return v
    return ""


def design_key(r):
    """The composition fields a Fig S4 rank row also prints."""
    return (norm(g(r, "Lipomer")),
            norm(g(r, "Amine")) + norm(g(r, "Lipid")),
            norm(g(r, "Mass Ratio", "Mass")),
            norm(g(r, "PEG Chain", "PEG")),
            norm(g(r, "PEG MW")),
            norm(g(r, "Lipomer Mole %", "Lipomer Mole")),
            norm(g(r, "Cholesterol Mole %", "Cholesterol Mole")),
            norm(g(r, "PEG Mole %")),
            norm(g(r, "Helper Lipid Mole %")))


def rank_key(r):
    return (norm(g(r, "Lipomer")),
            norm(g(r, "Amine Lipid")) or (norm(g(r, "Amine"))
                                          + norm(g(r, "Lipid"))),
            norm(g(r, "Mass Ratio", "Mass RatioPEG Chain")),
            norm(g(r, "PEG Chain", "Mass RatioPEG Chain #2")),
            norm(g(r, "PEG MW")),
            norm(g(r, "Lipomer Mole %")),
            norm(g(r, "Cholesterol Mole %")),
            norm(g(r, "PEG Mole %")),
            norm(g(r, "Helper Lipid Mole %")))


RANK_COLS = ["iMAEC", "RAWs", "Heart Ecs", "Lung Ecs", "Marrow Ecs",
             "Heart Macs", "Lung Macs", "Marrow Macs"]
INVIVO = ["Heart Ecs", "Lung Ecs", "Marrow Ecs",
          "Heart Macs", "Lung Macs", "Marrow Macs"]


def run(design_csv="01_design_tables.csv", rank_csv="02_rank_tables.csv"):
    des = load(design_csv)
    rank = load(rank_csv)
    out = {}

    # ---- C-A  library sizes against the stated counts -------------------
    n_by_lib = Counter(r["library"] for r in des)
    out["A_library_row_counts"] = {
        "extracted": dict(sorted(n_by_lib.items())),
        "stated": {k: v[0] for k, v in STATED.items()},
        "n_matching": sum(1 for k, v in STATED.items()
                          if n_by_lib.get(k) == v[0]),
        "n_libraries": len(STATED),
    }

    # ---- C-B  mole percents sum to exactly 100 on every row -------------
    sums = Counter()
    for r in des:
        s = 0.0
        for c in ("Lipomer Mole %", "Lipomer Mole", "Cholesterol Mole %",
                  "Cholesterol Mole", "PEG Mole %", "Helper Lipid Mole %"):
            v = r.get(c, "")
            if v not in ("", None):
                s += float(v)
        sums[round(s, 6)] += 1
    out["B_mole_percent_sums"] = {
        "distribution": {str(k): v for k, v in sorted(sums.items())},
        "n_rows": len(des),
        "n_exactly_100": sums.get(100.0, 0),
    }

    # ---- C-B2  the same invariant on the Fig S4 RANK tables --------------
    # This is the check that decides a disagreement between two printings of
    # one particle: the four mole percents summing to 100 is a constraint the
    # authors imposed and never printed, so it is external to both sides.
    rsums, offenders = Counter(), []
    for r in rank:
        s2 = 0.0
        for c in ("Lipomer Mole %", "Cholesterol Mole %", "PEG Mole %",
                  "Helper Lipid Mole %"):
            v = r.get(c, "")
            if v not in ("", None):
                s2 += float(v)
        rsums[round(s2, 6)] += 1
        if round(s2, 6) != 100.0:
            offenders.append({"subtable": r.get("subtable"),
                              "si_page": r.get("si_page"),
                              "lipomer": r.get("Lipomer"),
                              "lipomer_mole_pct": r.get("Lipomer Mole %"),
                              "chol_mole_pct": r.get("Cholesterol Mole %"),
                              "peg_mole_pct": r.get("PEG Mole %"),
                              "helper_mole_pct": r.get("Helper Lipid Mole %"),
                              "diameter_nm": r.get("Diameter (nm)"),
                              "sum": round(s2, 6)})
    out["B2_rank_table_mole_percent_sums"] = {
        "distribution": {str(k): v for k, v in sorted(rsums.items())},
        "n_rows": len(rank),
        "n_exactly_100": rsums.get(100.0, 0),
        "offenders": offenders,
    }

    # ---- C-C  stable-pool reconciliation against the printed 66.9% ------
    formulated = sum(STATED[k][0] for k in "123")
    stable = sum(STATED[k][1] for k in "123")
    out["C_stable_pool"] = {
        "formulated_libs_1_3": formulated,
        "stable_libs_1_3": stable,
        "tracker_and_abstract_say": 281,
        "matches_281": stable == 281,
        "pct_included": round(100.0 * stable / formulated, 1),
        "fig_S3i_prints_included_pct": 66.9,
        "pct_removed": round(100.0 * (formulated - stable) / formulated, 1),
        "fig_S3i_prints_removed_pct": 33.1,
    }

    # ---- C-D  block inventory -------------------------------------------
    blocks = defaultdict(set)
    for r in des:
        blocks[r["library"]].add(g(r, "Lipomer"))
    out["D_distinct_lipomers_per_library"] = {
        k: len(v) for k, v in sorted(blocks.items())}
    out["D_stated_per_library"] = {"1": 2, "2": 15, "3": 2, "4": 3}
    ids = set()
    for v in blocks.values():
        ids |= v
    out["D_distinct_lipomer_ids_all_libraries"] = sorted(ids)
    out["D_n_distinct_lipomer_ids"] = len(ids)
    amines = {norm(g(r, "Amine")) for r in des if g(r, "Amine")}
    tails = {norm(g(r, "Lipid")) for r in des if g(r, "Lipid")}
    out["D_distinct_amines"] = len(amines)
    out["D_distinct_tails"] = len(tails)

    # ---- C-E  every Fig S4 rank row matches a design row ----------------
    dkeys = defaultdict(list)
    for r in des:
        dkeys[design_key(r)].append(r)
    hit = miss = 0
    misses = []
    for r in rank:
        k = rank_key(r)
        if k in dkeys:
            hit += 1
        else:
            miss += 1
            if len(misses) < 6:
                misses.append({c: r.get(c, "") for c in
                               ("subtable", "Lipomer", "Amine", "Lipid",
                                "Amine Lipid", "PEG Chain", "PEG MW",
                                "Cholesterol Mole %")})
    out["E_rank_rows_matched_to_a_design_row"] = {
        "n_rank_rows": len(rank), "matched": hit, "unmatched": miss,
        "example_misses": misses,
    }

    # ---- C-F  rank columns are permutations without repeats -------------
    per = defaultdict(list)
    for r in rank:
        for c in RANK_COLS:
            v = r.get(c, "")
            if v not in ("", None):
                per[(r["subtable"], r["si_page"], c)].append(int(v))
    dupes, invivo_over_10, n_cols = 0, 0, 0
    detail = {}
    for k, vs in sorted(per.items()):
        n_cols += 1
        d = len(vs) - len(set(vs))
        dupes += d
        over = sum(1 for v in vs if k[2] in INVIVO and v > 10)
        invivo_over_10 += over
        detail["%s p%s %s" % k] = {"n": len(vs), "min": min(vs),
                                   "max": max(vs), "repeats": d,
                                   "gt_10_in_vivo": over}
    out["F_rank_column_integrity"] = {
        "n_rank_columns": n_cols,
        "repeated_ranks_total": dupes,
        "in_vivo_ranks_above_10": invivo_over_10,
        "per_column": detail,
    }

    # ---- C-H  in vitro ranks must fit inside their own stable pool ------
    # An independent handle on the rank scale: the worst in vitro rank listed
    # for a library cannot exceed that library's stable pool, and for library
    # 3 it should land ON it, because the last-coloured particle there is the
    # worst in the pool.
    pool = {"Library 1": STATED["1"][1], "Library 2": STATED["2"][1],
            "Library 3": STATED["3"][1]}
    hcheck = {}
    for lab, n in pool.items():
        for col in ("iMAEC", "RAWs"):
            vs = [int(r[col]) for r in rank
                  if r.get("subtable") == lab and r.get(col)]
            if vs:
                hcheck["%s %s" % (lab, col)] = {
                    "max_rank_listed": max(vs), "stable_pool": n,
                    "fits": max(vs) <= n, "equals_pool": max(vs) == n}
    out["H_in_vitro_rank_vs_stable_pool"] = {
        "per_column": hcheck,
        "n_columns": len(hcheck),
        "n_fitting": sum(1 for v in hcheck.values() if v["fits"]),
        "n_hitting_the_pool_exactly": sum(
            1 for v in hcheck.values() if v["equals_pool"]),
    }

    # ---- C-I  each library is a COMPLETE factorial ----------------------
    fact = {}
    for lib in sorted(n_by_lib):
        sub = [r for r in des if r["library"] == lib]
        lev = {}
        for col in ("Lipomer", "PEG Chain", "PEG", "PEG MW"):
            vals = {g(r, col) for r in sub if g(r, col)}
            if vals:
                lev[col] = len(vals)
        chol = {(g(r, "Lipomer Mole %", "Lipomer Mole"),
                 g(r, "Cholesterol Mole %", "Cholesterol Mole"),
                 g(r, "PEG Mole %"), g(r, "Helper Lipid Mole %"))
                for r in sub}
        prod = 1
        for v in lev.values():
            prod *= v
        fact[lib] = {"factor_levels": lev, "distinct_mole_pct_sets": len(chol),
                     "product": prod * len(chol), "n_rows": len(sub),
                     "is_complete_factorial": prod * len(chol) == len(sub)}
    out["I_complete_factorial"] = {
        "per_library": fact,
        "n_libraries_complete": sum(
            1 for v in fact.values() if v["is_complete_factorial"]),
    }

    # ---- C-G  how much chemistry is actually varied ---------------------
    pairs = Counter()
    for r in des:
        pairs[g(r, "Lipomer")] += 1
    out["G_formulations_per_distinct_lipomer_chemistry"] = {
        "n_distinct_lipomers": len(pairs),
        "n_formulations": len(des),
        "median_formulations_each": sorted(pairs.values())[len(pairs) // 2],
        "max_formulations_on_one_chemistry": max(pairs.values()),
    }
    return out


def main():
    rep = run()
    with open(os.path.join(HERE, "03_checks.json"), "w") as fh:
        json.dump(rep, fh, indent=1)
    print(json.dumps(rep, indent=1)[:4000])
    return 0


if __name__ == "__main__":
    sys.exit(main())
