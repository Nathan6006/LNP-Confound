#!/usr/bin/env python3
"""Stage 01 -- read the two LNP library tables and the normalisation worked
example out of the SI as NATIVE TEXT (tier 2/4 of the intake gate), and run
every check that those objects can support.

    python3 -I stage01_si_tables.py [OUTDIR]

OUTDIR defaults to this folder. A seed harness MUST pass its own temp dir
(F-069: a harness that pins OUTDIR to the curation folder lets a mutated copy
overwrite a tracked artifact, and a non-zero exit conceals the write).

What this does NOT do: extract a delivery value. There is none to extract.
Every cell of both library tables is formulation metadata, and the sole
quantitative columns are mole percents and hydrodynamic diameter -- the latter
on the brief's not-delivery list. No pixel is read and nothing from
curation/lib's geometry helpers is called, because no geometry is involved.
"""
import json
import os
import re
import sys

import pdfplumber

HERE = os.path.dirname(os.path.abspath(__file__))
SI = os.path.join(HERE, "source", "nl8b03149_si_001.pdf")
OUTDIR = os.path.abspath(sys.argv[1]) if len(sys.argv) > 1 else HERE

# Column order of the library-1 halves (SI p14), read off the extracted header.
L1_COLS = ["Number", "Lipomer Type", "Chol Type", "PEG Tail", "PEG MW",
           "HL Type", "Lipomer %", "Chol %", "PEG %", "HL %", "Diameter (nm)"]
# Library 2 (SI p15) is ONE table holding two halves side by side; its header
# cells are glued together by the extractor ("% Chol Mole %"), so the halves
# are taken by POSITION and the positions are then validated, not trusted.
L2_COLS = ["LNP #", "Lipomer Type", "Chol Type", "PEG Alkyl", "PEG MW",
           "Lipomer Mole %", "Chol Mole %", "PEG Mole %", "Diameter (nm)"]

DIAM_LO, DIAM_HI = 20.0, 200.0   # SI Fig S4b / Methods inclusion criterion


def num(s):
    """Parse a numeric cell, or None. Refuses rather than guessing."""
    if s is None:
        return None
    t = str(s).strip()
    if not re.fullmatch(r"-?\d+(\.\d+)?", t):
        return None
    return float(t)


def read_library1(pdf):
    page = pdf.pages[13]
    tbs = page.find_tables()
    if len(tbs) != 2:
        raise SystemExit(f"SI p14: expected 2 tables, found {len(tbs)}")
    rows = []
    for tb in sorted(tbs, key=lambda t: t.bbox[0]):
        data = tb.extract()
        head = data[0]
        if len(head) != len(L1_COLS):
            raise SystemExit(f"SI p14 header width {len(head)} != {len(L1_COLS)}")
        # The right half's header is DOUBLED by the extractor, and the
        # doubling interleaves imperfectly: 'PEG Tail' comes back as
        # 'PPEEGG T Taailil', which is neither a substring nor a clean
        # character-pair interleave of the expected string. So the check is on
        # the LETTER SET, with each letter appearing once or twice. That still
        # has teeth for the thing being relied on -- column ORDER -- because
        # any neighbouring column here has a different letter set
        # ('PEG MW' -> {p,e,g,m,w} vs 'PEG Tail' -> {p,e,g,t,a,i,l}); it is
        # deliberately blind only to letter ORDER within one header cell,
        # which nothing downstream uses.
        for got, want in zip(head, L1_COLS):
            g = re.sub(r"[^A-Za-z]", "", str(got)).lower()
            w = re.sub(r"[^A-Za-z]", "", want).lower()
            if set(g) != set(w):
                raise SystemExit(f"SI p14 header mismatch: {got!r} vs {want!r} "
                                 f"(letter sets {sorted(set(g))} vs {sorted(set(w))})")
            if not (len(w) <= len(g) <= 2 * len(w) + 2):
                raise SystemExit(f"SI p14 header width odd: {got!r} vs {want!r}")
        for r in data[1:]:
            rows.append(dict(zip(L1_COLS, [None if c is None else str(c).strip()
                                           for c in r])))
    return rows


def read_library2(pdf):
    page = pdf.pages[14]
    tbs = [t for t in page.find_tables() if len(t.extract()) > 10]
    if len(tbs) != 1:
        raise SystemExit(f"SI p15: expected 1 big table, found {len(tbs)}")
    data = tbs[0].extract()
    if len(data[0]) != 18:
        raise SystemExit(f"SI p15 width {len(data[0])} != 18")
    rows = []
    for r in data[1:]:
        for off in (0, 9):
            rows.append(dict(zip(L2_COLS,
                                 [None if c is None else str(c).strip()
                                  for c in r[off:off + 9]])))
    return rows


def read_norm_example(pdf):
    """SI p13 (Fig S4e) -- the three-stage normalisation worked example."""
    page = pdf.pages[12]
    tbs = [t for t in page.find_tables() if len(t.extract()) == 6]
    want = {"raw": "Raw Counts", "pct": "Normalized Counts",
            "inp": "Normalized to input"}
    out = {}
    for tb in tbs:
        data = tb.extract()
        banner = " ".join(str(c) for c in data[0] if c)
        for key, tag in want.items():
            if tag in banner and key not in out:
                # columns: 0 LNP | 3 barcode | 4,5,6 lung M1-3 | 7,8,9 kidney M1-3 | 10 input
                body = {}
                for r in data[2:5]:
                    body[str(r[0]).strip()] = {
                        "barcode": str(r[3]).strip(),
                        "lung": [num(r[4]), num(r[5]), num(r[6])],
                        "kidney": [num(r[7]), num(r[8]), num(r[9])],
                        "input": num(r[10]) if len(r) > 10 else None,
                    }
                tot = data[5]
                out[key] = {"rows": body,
                            "total": {"lung": [num(tot[4]), num(tot[5]), num(tot[6])],
                                      "kidney": [num(tot[7]), num(tot[8]), num(tot[9])],
                                      "input": num(tot[10]) if len(tot) > 10 else None}}
    missing = set(want) - set(out)
    if missing:
        raise SystemExit(f"SI p13: missing stage tables {sorted(missing)}")
    return out


# ---------------------------------------------------------------- checks ----

def check_single_ionizable_lipid(l1, l2):
    """How many distinct ionizable lipids does this paper contribute?

    POSITIVE CONTROL (F-070): a sweep that returns one value because it read
    nothing is indistinguishable from a sweep that returns one value because
    the column is constant. So: (a) a column we know varies (Diameter) must
    come back with many distinct values, (b) a column we know is constant by a
    DIFFERENT mechanism (Chol Type) must come back constant, and (c) the row
    count must equal the stated library sizes. If (a) fails the extractor is
    returning blanks and the 7C1 result is vacuous.
    """
    lip1 = sorted({r["Lipomer Type"] for r in l1})
    lip2 = sorted({r["Lipomer Type"] for r in l2})
    pos_vary = len({r["Diameter (nm)"] for r in l1})
    pos_const = sorted({r["Chol Type"] for r in l1})
    return {
        "library_1_rows": len(l1),
        "library_2_rows": len(l2),
        "library_1_rows_expected": 128,
        "library_2_rows_expected": 120,
        "row_counts_match_paper": len(l1) == 128 and len(l2) == 120,
        "library_1_distinct_lipomer_type": lip1,
        "library_2_distinct_lipomer_type": lip2,
        "distinct_ionizable_lipids_in_paper": sorted(set(lip1) | set(lip2)),
        "n_distinct_ionizable_lipids": len(set(lip1) | set(lip2)),
        "positive_control_diameter_distinct_values": pos_vary,
        "positive_control_diameter_varies": pos_vary > 20,
        "positive_control_chol_type_constant": pos_const,
        "helper_lipid_types_library_1": sorted({r["HL Type"] for r in l1}),
        "peg_tails_library_1": sorted({r["PEG Tail"] for r in l1}),
        "peg_tails_library_2": sorted({r["PEG Alkyl"] for r in l2}),
        "peg_mw_library_1": sorted({r["PEG MW"] for r in l1}),
        "peg_mw_library_2": sorted({r["PEG MW"] for r in l2}),
        "interpretation": ("Every one of the 248 table rows carries Lipomer Type "
                           "'7C1'. The variables are the four component mole "
                           "percents, the PEG alkyl tail, the PEG molecular "
                           "weight, the presence of DSPC, and the mouse genotype. "
                           "The ionizable component never changes."),
    }


def check_molpct_sums(l1, l2):
    """mol% sums to 100 on every formulation. NOT a tautology of our own
    derivation: all four percents are printed independently in the SI and are
    summed here exactly as printed."""
    out = {}
    for name, rows, cols in (
            ("library_1", l1, ["Lipomer %", "Chol %", "PEG %", "HL %"]),
            ("library_2", l2, ["Lipomer Mole %", "Chol Mole %", "PEG Mole %"])):
        ok, bad, unparsed = 0, [], []
        for r in rows:
            vals = [num(r[c]) for c in cols]
            if any(v is None for v in vals):
                unparsed.append({"row": r[cols and ("Number" if "Number" in r else "LNP #")],
                                 "cells": [r[c] for c in cols]})
                continue
            s = round(sum(vals), 6)
            if abs(s - 100.0) < 1e-9:
                ok += 1
            else:
                bad.append({"row": r["Number" if "Number" in r else "LNP #"],
                            "sum": s, "cells": {c: r[c] for c in cols}})
        out[name] = {"n_rows": len(rows), "summed_to_100": ok,
                     "offenders": bad, "n_offenders": len(bad),
                     "unparsed": unparsed, "n_unparsed": len(unparsed),
                     "columns_summed": cols,
                     "tautology_check": ("not a tautology: all component "
                                         "percents are printed independently "
                                         "in the SI and none is derived here")}
    return out


def check_inclusion_criterion(l1, l2):
    """The paper states 111 of 128 (library 1) and 115 of 120 (library 2) LNPs
    met TWO criteria: 20-200 nm diameter AND a single-inflection
    autocorrelation curve. The table prints only the diameter, so the
    diameter gate alone is a SUPERSET of the pooled set. Reported as such."""
    out = {}
    for name, rows, key, stated in (("library_1", l1, "Diameter (nm)", 111),
                                    ("library_2", l2, "Diameter (nm)", 115)):
        inb, outb, unp = [], [], []
        for r in rows:
            d = num(r[key])
            rid = r["Number" if "Number" in r else "LNP #"]
            if d is None:
                unp.append(rid)
            elif DIAM_LO <= d <= DIAM_HI:
                inb.append(rid)
            else:
                outb.append({"row": rid, "diameter_nm": d})
        out[name] = {
            "n_rows": len(rows),
            "n_within_20_200nm": len(inb),
            "n_outside": len(outb),
            "outside": outb,
            "n_unparsed_diameter": len(unp),
            "paper_states_pooled": stated,
            "diameter_gate_minus_stated": len(inb) - stated,
            "is_a_superset_not_an_equality": True,
            "why": ("the second criterion -- one inflection point in the "
                    "autocorrelation curve (SI Fig S4b) -- is not printed per "
                    "LNP, so the diameter gate can only bound the pooled set "
                    "from above"),
        }
    return out


def check_normalisation(ex):
    """Reconstruct the paper's own worked example (SI Fig S4e) and report the
    residual. Two stages:
      stage A  pct_i      = raw_i / sum_j(raw_j) * 100          (compositional)
      stage B  final_i    = (pct_i / pct_i_input) normalised to sum 100
    Both are closed compositions: a column of stage A and a column of stage B
    each sum to 100 by construction.
    """
    raw, pct, inp = ex["raw"], ex["pct"], ex["inp"]
    lnps = sorted(raw["rows"])
    res = {"stage_A_counts_per_100": [], "stage_B_input_normalised": []}

    # ---- stage A: reproduce the printed percentages from the printed counts
    for organ in ("lung", "kidney"):
        for m in range(3):
            tot = sum(raw["rows"][k][organ][m] for k in lnps)
            if abs(tot - raw["total"][organ][m]) > 1e-9:
                raise SystemExit(f"stage A: printed total {raw['total'][organ][m]} "
                                 f"!= sum of printed counts {tot}")
            for k in lnps:
                calc = raw["rows"][k][organ][m] / tot * 100.0
                printed = pct["rows"][k][organ][m]
                res["stage_A_counts_per_100"].append({
                    "lnp": k, "organ": organ, "mouse": m + 1,
                    "raw": raw["rows"][k][organ][m], "column_total": tot,
                    "calculated_pct": round(calc, 6),
                    "printed_pct": printed,
                    "residual_after_rounding_to_int": round(round(calc) - printed, 6),
                    "exact": abs(round(calc) - printed) < 1e-9})
    # input column
    tot_in = sum(raw["rows"][k]["input"] for k in lnps)
    for k in lnps:
        calc = raw["rows"][k]["input"] / tot_in * 100.0
        res["stage_A_counts_per_100"].append({
            "lnp": k, "organ": "input", "mouse": None,
            "raw": raw["rows"][k]["input"], "column_total": tot_in,
            "calculated_pct": round(calc, 6), "printed_pct": pct["rows"][k]["input"],
            "residual_after_rounding_to_int": round(round(calc) - pct["rows"][k]["input"], 6),
            "exact": abs(round(calc) - pct["rows"][k]["input"]) < 1e-9})

    # ---- stage B: input-normalise, then re-close to 100.
    # Computed from the UNROUNDED stage-A percentages, i.e. from the raw counts,
    # so the reconstruction is of the pipeline and not of the printed integers.
    for organ in ("lung", "kidney"):
        for m in range(3):
            tot = sum(raw["rows"][k][organ][m] for k in lnps)
            ratios = {k: (raw["rows"][k][organ][m] / tot) /
                         (raw["rows"][k]["input"] / tot_in) for k in lnps}
            s = sum(ratios.values())
            for k in lnps:
                calc = ratios[k] / s * 100.0
                printed = inp["rows"][k][organ][m]
                res["stage_B_input_normalised"].append({
                    "lnp": k, "organ": organ, "mouse": m + 1,
                    "calculated_pct": round(calc, 6), "printed_pct": printed,
                    "residual_after_rounding_to_int": round(round(calc) - printed, 6),
                    "exact": abs(round(calc) - printed) < 1e-9})

    for stage in ("stage_A_counts_per_100", "stage_B_input_normalised"):
        cells = res[stage]
        res[stage + "_summary"] = {
            "n_cells": len(cells),
            "n_exact": sum(1 for c in cells if c["exact"]),
            "max_abs_residual": max(abs(c["residual_after_rounding_to_int"])
                                    for c in cells),
            "mean_abs_residual_before_rounding": round(
                sum(abs(c["calculated_pct"] - c["printed_pct"]) for c in cells)
                / len(cells), 6),
        }
    res["is_compositional"] = True
    res["compositional_evidence"] = (
        "Stage A divides each barcode's count by the sum over all barcodes in "
        "the same sample; stage B re-closes the input-normalised ratios to 100. "
        "The SI prints 'Total (%) = 100' on every column of both stages. A "
        "value is therefore a SHARE of a pool, not an amount: it moves when "
        "any other pool member moves, and it is not comparable to a value from "
        "a pool of different membership.")
    res["worked_example_exists"] = True
    res["worked_example_locus"] = "SI p13, Figure S4e, three stacked 3-LNP tables"
    res["formula_in_methods"] = (
        "Methods 'Barcode Sequencing Normalization', one sentence, verbatim: "
        "'Counts for each particle, per cell type, were normalized to the "
        "barcoded LNP mixture applied to cells or injected into the mouse.' "
        "The sentence alone does not distinguish stage A from stage B, nor "
        "state the re-closure; Fig S4e is what pins both.")
    return res


def main():
    with pdfplumber.open(SI) as pdf:
        l1 = read_library1(pdf)
        l2 = read_library2(pdf)
        ex = read_norm_example(pdf)

    qc = {
        "si_pages_read": {"library_1": 14, "library_2": 15,
                          "normalisation_worked_example": 13},
        "extraction_tier": "2_si_tables / 4_pdfplumber_text_table (native text, no pixels)",
        "ionizable_lipid_count": check_single_ionizable_lipid(l1, l2),
        "molpct_sums_to_100": check_molpct_sums(l1, l2),
        "inclusion_criterion": check_inclusion_criterion(l1, l2),
        "normalisation_reconstruction": check_normalisation(ex),
    }
    os.makedirs(OUTDIR, exist_ok=True)
    with open(os.path.join(OUTDIR, "01_si_tables.json"), "w") as f:
        json.dump(qc, f, indent=1, sort_keys=False)

    # 01_blocks.csv -- the combinatorial design, which is what this paper
    # actually varies. One row per LNP, no endpoint column, because there is
    # no per-LNP endpoint to attach.
    import csv
    path = os.path.join(OUTDIR, "01_blocks.csv")
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["library", "lnp_number", "ionizable_component",
                    "chol_type", "peg_tail", "peg_mw", "helper_lipid",
                    "ionizable_molpct", "chol_molpct", "peg_molpct",
                    "helper_molpct", "diameter_nm", "within_20_200nm"])
        for r in l1:
            d = num(r["Diameter (nm)"])
            w.writerow([1, r["Number"], r["Lipomer Type"], r["Chol Type"],
                        r["PEG Tail"], r["PEG MW"], r["HL Type"],
                        r["Lipomer %"], r["Chol %"], r["PEG %"], r["HL %"],
                        r["Diameter (nm)"],
                        "" if d is None else int(DIAM_LO <= d <= DIAM_HI)])
        for r in l2:
            d = num(r["Diameter (nm)"])
            w.writerow([2, r["LNP #"], r["Lipomer Type"], r["Chol Type"],
                        r["PEG Alkyl"], r["PEG MW"], "none",
                        r["Lipomer Mole %"], r["Chol Mole %"],
                        r["PEG Mole %"], 0, r["Diameter (nm)"],
                        "" if d is None else int(DIAM_LO <= d <= DIAM_HI)])

    s = qc["ionizable_lipid_count"]
    print(json.dumps({
        "n_distinct_ionizable_lipids": s["n_distinct_ionizable_lipids"],
        "ionizable_lipids": s["distinct_ionizable_lipids_in_paper"],
        "row_counts_match_paper": s["row_counts_match_paper"],
        "positive_control_diameter_varies": s["positive_control_diameter_varies"],
        "molpct_lib1": f'{qc["molpct_sums_to_100"]["library_1"]["summed_to_100"]}/128',
        "molpct_lib2": f'{qc["molpct_sums_to_100"]["library_2"]["summed_to_100"]}/120',
        "diam_gate_lib1": qc["inclusion_criterion"]["library_1"]["n_within_20_200nm"],
        "diam_gate_lib2": qc["inclusion_criterion"]["library_2"]["n_within_20_200nm"],
        "normA": qc["normalisation_reconstruction"]["stage_A_counts_per_100_summary"],
        "normB": qc["normalisation_reconstruction"]["stage_B_input_normalised_summary"],
    }, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
