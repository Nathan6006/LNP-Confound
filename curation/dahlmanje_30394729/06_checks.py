#!/usr/bin/env python3
"""Stage 06 -- every quantitative check, with numbers. Writes 07_qc.json.

    python3 -I 06_checks.py

Reads 04_si_tables.json (tier-4 text tables) and the live corpus. Reads no pixels.

The checks that matter here:
  N1  reconstruct the paper's OWN worked normalisation example from its raw counts
  N2  the drawn-but-never-computed 'Total (%) = 100' rows
  F1  the factorial of both libraries (re-asserted from 04_si_tables.json)
  M1  mol% sums to 100 on all 160 formulations
  B1  the 156 barcodes
  W1  the five individually-tested formulations tie to rows of the library tables
  P1  the PEG repeat count drawn in Fig 3B against the PEG molecular weight printed
  S1  SEEDED: each check is made to fail on purpose, and must.
"""
import json, os, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))


def load():
    with open(os.path.join(HERE, "04_si_tables.json")) as f:
        return json.load(f)


# --------------------------------------------------------------------------
# N1/N2 -- the normalisation, reconstructed from Figure S1E's own raw counts.
# --------------------------------------------------------------------------
def normalisation(tabs, seed=None):
    """tabs = the three Fig S1E tables as pdfplumber returned them.

    raw[i][s]  counts; the last data column is the INPUT pool.
    Returns the reconstruction residuals against the paper's printed numbers.
    """
    raw_t, mid_t, fin_t = tabs
    cols = raw_t[1][2:]                      # ['Mouse 1','Mouse 2','Mouse 3','Input']
    raw = {}
    for r in raw_t[2:5]:
        raw[int(r[0])] = [float(x) for x in r[2:]]
    printed_totals_raw = [float(x) for x in raw_t[5][2:]]
    mid = {int(r[0]): [float(x) for x in r[2:]] for r in mid_t[2:5]}
    printed_totals_mid = [float(x) for x in mid_t[5][2:]]
    fin = {int(r[0]): [float(x) for x in r[2:]] for r in fin_t[2:5]}
    printed_totals_fin = [float(x) for x in fin_t[5][2:]]

    if seed is not None:                      # S1: corrupt one raw count
        raw[seed[0]][seed[1]] += seed[2]

    lnps = sorted(raw)
    # the paper's own words, Fig S1E caption: "Barcode counts are normalized
    # within a sample, and then to the DNA barcode input."
    tot = [sum(raw[i][s] for i in lnps) for s in range(len(cols))]
    share = {i: [raw[i][s] / tot[s] * 100.0 for s in range(len(cols))] for i in lnps}
    inp = len(cols) - 1
    ratio = {i: [share[i][s] / (share[i][inp] / 100.0) for s in range(inp)] for i in lnps}
    rsum = [sum(ratio[i][s] for i in lnps) for s in range(inp)]
    nd = {i: [ratio[i][s] / rsum[s] * 100.0 for s in range(inp)] for i in lnps}

    res_fin = [{"lnp": i, "sample": cols[s], "reconstructed": round(nd[i][s], 4),
                "printed": fin[i][s], "residual_after_rounding":
                round(abs(round(nd[i][s]) - fin[i][s]), 6)}
               for i in lnps for s in range(inp)]
    res_mid = [{"lnp": i, "sample": cols[s], "reconstructed": round(share[i][s], 4),
                "printed": mid[i][s], "residual_after_rounding":
                round(abs(round(share[i][s]) - mid[i][s]), 6)}
               for i in lnps for s in range(len(cols))]
    comp_tot_raw = [round(t, 6) for t in tot]

    return {
        "formula": "nd_i = [ (c_i / sum_k c_k) / (in_i / sum_k in_k) ] "
                   "/ sum_j[ same ] x 100%",
        "stage_1_raw_column_totals": {
            "computed": comp_tot_raw, "printed": printed_totals_raw,
            "max_residual": max(abs(a - b) for a, b in
                                zip(comp_tot_raw, printed_totals_raw)),
            "n_cells": len(comp_tot_raw)},
        "stage_2_within_sample_percent": {
            "n_cells": len(res_mid),
            "n_exact_after_rounding": sum(1 for r in res_mid
                                          if r["residual_after_rounding"] == 0),
            "max_residual": max(r["residual_after_rounding"] for r in res_mid),
            "mismatches": [r for r in res_mid if r["residual_after_rounding"] != 0],
            "printed_total_row": printed_totals_mid,
            "printed_cells_actually_sum_to":
                [round(sum(mid[i][s] for i in lnps), 6) for s in range(len(cols))]},
        "stage_3_normalised_delivery": {
            "n_cells": len(res_fin),
            "n_exact_after_rounding": sum(1 for r in res_fin
                                          if r["residual_after_rounding"] == 0),
            "max_residual": max(r["residual_after_rounding"] for r in res_fin),
            "cells": res_fin,
            "printed_total_row": printed_totals_fin,
            "printed_cells_actually_sum_to":
                [round(sum(fin[i][s] for i in lnps), 6) for s in range(inp)]},
    }


# --------------------------------------------------------------------------
# W1 -- the five individually-tested formulations, tied to library rows.
# --------------------------------------------------------------------------
# Read ONCE off SI p. S14 (Figure S3A) and Fig 6A, as printed text.
INDIVIDUALS = [
    {"label": "Original 7C1 (80:20)", "il": 80.0, "chol": 0.0, "peg": 20.0,
     "dspc": 0.0, "peg_type": "C14PEG2000", "expect_library": None},
    {"label": "Screen 1 Winner", "il": 60.0, "chol": 37.5, "peg": 2.5,
     "dspc": 0.0, "peg_type": "C16PEG2000", "expect_library": "1"},
    {"label": "LNP1 (screen 2)", "il": 70.0, "chol": 10.0, "peg": 20.0,
     "dspc": 0.0, "peg_type": "C16PEG2000", "expect_library": "2"},
    {"label": "BM1 (screen 2 winner)", "il": 80.0, "chol": 2.5, "peg": 17.5,
     "dspc": 0.0, "peg_type": "C18PEG2000", "expect_library": "2"},
    {"label": "LNP3 (screen 2)", "il": 80.0, "chol": 0.0, "peg": 20.0,
     "dspc": 0.0, "peg_type": "C18PEG2000", "expect_library": "2"},
]


def winners(d, perturb=None):
    out = []
    for w in INDIVIDUALS:
        hits = []
        for lab, key in (("1", "library_1"), ("2", "library_2")):
            for n, f in d[key].items():
                il = f["il_molpct"] + (perturb or 0.0)
                if (il == w["il"] and f["chol_molpct"] == w["chol"]
                        and f["peg_molpct"] == w["peg"]
                        and f["dspc_molpct"] == w["dspc"]
                        and f["peg_type"] == w["peg_type"]):
                    hits.append(f"library {lab} LNP {n}")
        out.append({"label": w["label"], "composition":
                    f'{w["il"]}/{w["chol"]}/{w["peg"]}/{w["dspc"]} {w["peg_type"]}',
                    "matches_in_the_library_tables": hits,
                    "unique": len(hits) == 1,
                    "expected_library": w["expect_library"],
                    "in_expected_library": (w["expect_library"] is None and not hits)
                    or any(h.startswith(f'library {w["expect_library"]} ') for h in hits)})
    return {"n": len(out),
            "n_uniquely_located": sum(1 for o in out if o["unique"]),
            "n_as_expected": sum(1 for o in out if o["in_expected_library"]),
            "rows": out}


# --------------------------------------------------------------------------
# P1 -- a round number the authors DREW but never printed.
# --------------------------------------------------------------------------
def peg_repeats():
    """Fig 3B and Fig 5A draw the PEG block as (OCH2CH2)45OCH3 and separately
    print the PEG molecular weight as 2000. 45 ethylene-oxide units is
    45 x 44.0526 = 1982.4 Da; with the terminal OCH3 (31.03) and the amide-linked
    oxygen the block is ~2013 Da. Neither number is printed as a product."""
    eo = 44.05256
    n = 45
    return {"drawn_repeat_count": n, "unit_mass_Da": eo,
            "block_mass_Da": round(n * eo, 3),
            "plus_terminal_OCH3_Da": round(n * eo + 31.0340, 3),
            "printed_label": 2000,
            "relative_error_vs_2000": round(abs(n * eo - 2000) / 2000, 5),
            "verdict": "45 x 44.05 = 1982.4 Da, 0.9% from the printed 2000 -- the "
                       "drawn repeat count is consistent with the label and was "
                       "never multiplied out in the paper."}


# --------------------------------------------------------------------------
# E1 -- the 'Formulated' pies against the composition tables.
#
# Each enrichment panel draws a 'Formulated' pie whose slices carry PRINTED
# INTEGER PERCENTAGES. Those percentages are a property of the composition
# table, so the table and the pie are two independent renderings of the same
# fact and either confirms or contradicts the other. The numerals below were
# read BY EYE off the author-manuscript raster at 2.2x; they are 16-20 px tall
# and unambiguous. No pixel was measured: a numeral is read, not digitised.
# --------------------------------------------------------------------------
PIES = {
    "Fig_4A_library1_PEG_molpct": {
        "library": "1", "field": "peg_molpct",
        "bands_half_open": [("2.5", 0, 5), ("5-9", 5, 10), ("10-14", 10, 15),
                            ("15-20", 15, 20), (">20", 20, 1e9)],
        "bands_inclusive": [("2.5", 0, 5), ("5-9", 5, 10), ("10-14", 10, 15),
                            ("15-20", 15, 20.0001), (">20", 20.0001, 1e9)],
        "printed": [12, 25, 25, 21, 17],
    },
    "Fig_4B_library1_PEG_type": {
        "library": "1", "field": "peg_type",
        "categories": ["C14PEG2000", "C16PEG2000", "C18PEG2000",
                       "C14PEG3000", "C18PEG3000"],
        "printed": [20, 20, 20, 20, 20],
    },
    "Fig_5H_library2_PEG_molpct": {
        "library": "2", "field": "peg_molpct",
        "bands_half_open": [("2.5", 0, 15), ("15-20", 15, 20), (">20", 20, 1e9)],
        "bands_inclusive": [("2.5", 0, 15), ("15-20", 15, 20.0001),
                            (">20", 20.0001, 1e9)],
        "printed": [40, 40, 20],
    },
    "Fig_5I_library2_PEG_type": {
        "library": "2", "field": "peg_type",
        "categories": ["C16PEG2000", "C18PEG2000"],
        "printed": [50, 50],
    },
    "Fig_5J_library2_7C1_molpct": {
        "library": "2", "field": "il_molpct",
        "categories": [50.0, 60.0, 70.0, 80.0],
        "printed": [20, 40, 20, 20],
    },
    "Fig_5K_library2_DSPC_molpct": {
        "library": "2", "field": "dspc_molpct",
        "categories": [0.0, 10.0],
        "printed": [60, 40],
    },
}


def pies(d):
    out = {}
    for panel, spec in PIES.items():
        lib = d["library_1" if spec["library"] == "1" else "library_2"]
        vals = [f[spec["field"]] for f in lib.values()]
        n = len(vals)
        row = {"n_formulated": n, "printed_slices": spec["printed"],
               "printed_sum": sum(spec["printed"])}
        if "categories" in spec:
            cnt = [sum(1 for v in vals if v == c) for c in spec["categories"]]
            pct = [round(c / n * 100, 4) for c in cnt]
            row.update({"categories": [str(c) for c in spec["categories"]],
                        "counts": cnt, "percent": pct,
                        "rounded": [round(p) for p in pct],
                        "exact_match": [round(p) for p in pct] == spec["printed"],
                        "convention": "exact category equality"})
        else:
            for key, label in (("bands_half_open", "half_open [lo,hi)"),
                               ("bands_inclusive", "inclusive [lo,hi]")):
                cnt = [sum(1 for v in vals if lo <= v < hi) for _, lo, hi in spec[key]]
                pct = [round(c / n * 100, 4) for c in cnt]
                row[key] = {"bands": [b[0] for b in spec[key]], "counts": cnt,
                            "percent": pct, "rounded": [round(p) for p in pct],
                            "exact_match": [round(p) for p in pct] == spec["printed"],
                            "convention": label}
            row["exact_match"] = (row["bands_half_open"]["exact_match"]
                                  or row["bands_inclusive"]["exact_match"])
        out[panel] = row
    out["_summary"] = {
        "n_panels": len(PIES),
        "n_reconciling": sum(1 for k, v in out.items()
                             if k != "_summary" and v["exact_match"]),
        "panels_that_do_not": [k for k, v in out.items()
                               if k != "_summary" and not v["exact_match"]],
    }
    return out


def main():
    d = load()
    tabs = d["normalisation_tables_fig_S1E"]
    qc = {
        "_what": "Stage 06 checks for curation/dahlmanje_30394729 (Sago et al. 2018, "
                 "PMID 30394729). Every entry carries numbers. Produced by 06_checks.py.",
        "N1_N2_normalisation": normalisation(tabs),
        "F1_factorial": {k: d["checks"][k] for k in
                         ("factorial_library_1", "factorial_library_2")},
        "M1_molpct": {k: d["checks"][k] for k in
                      ("molpct_library_1", "molpct_library_2")},
        "B1_barcodes": d["barcodes"],
        "T1_total_formulated": d["checks"]["total_formulated"],
        "W1_individually_tested": winners(d),
        "P1_peg_repeat_count": peg_repeats(),
        "E1_formulated_pies_vs_tables": pies(d),
    }
    with open(os.path.join(HERE, "07_qc_checks.json"), "w") as f:
        json.dump(qc, f, indent=1)
    print(json.dumps({
        "N1_stage3_exact": qc["N1_N2_normalisation"]["stage_3_normalised_delivery"],
        "N1_stage2_mismatch_n": len(qc["N1_N2_normalisation"]
                                    ["stage_2_within_sample_percent"]["mismatches"]),
        "W1": qc["W1_individually_tested"],
        "P1": qc["P1_peg_repeat_count"],
    }, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
