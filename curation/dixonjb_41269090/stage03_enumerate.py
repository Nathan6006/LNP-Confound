#!/usr/bin/env python3
"""Stage 03 -- the formulations that were administered INDIVIDUALLY.

Only seven LNPs in this paper were ever administered on their own: the six
leads (LNP1, LNP2, LNP3, LNP4, LNP7, LNP11) and the MC3 benchmark. Every other
member of the 150-LNP library was only ever injected as part of a pool, and the
pooled barcode values are not published (see 00_value_sources.json).

The composition table is read from the SI page-21 NATIVE table (Table S1,
tier 4). It is then cross-checked row by row against the independent
SI page-2 native table (Fig S1, Screening 1) for the six cKK-E12 leads -- two
different tables on two different pages that must agree.

THE HEADLINE OF THIS STAGE: all six leads carry the SAME ionizable lipid. The
'library' of individually-measured LNPs has exactly ONE distinct ionizable-lipid
structure (cKK-E12); only the helper lipid and the four mole percentages vary.
MC3 is the seventh formulation and the second structure, and it is a benchmark.
"""
import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))

import pdfplumber                                   # noqa: E402
from rdkit import Chem, RDLogger                    # noqa: E402
from rdkit.Chem import Descriptors, rdMolDescriptors  # noqa: E402

RDLogger.DisableLog("rdApp.*")
SI = os.path.join(HERE, "source", "nn5c12080_si_001.pdf")

# Table S1's printed component tokens -> the component key in 01_blocks.csv.
TOKEN = {"cKK-E12": "cKK-E12", "MC3": "DLin-MC3-DMA",
         "Cholesterol": "Cholesterol", "C14PEG2K": "C14PEG2K",
         "DSPC": "DSPC", "DOPE": "DOPE"}


def blocks():
    with open(os.path.join(HERE, "01_blocks.csv")) as fh:
        return {r["component"]: r for r in csv.DictReader(fh)}


def native(page_idx, expect_rows):
    with pdfplumber.open(SI) as doc:
        t = doc.pages[page_idx].find_tables()
        if len(t) != 1:
            raise SystemExit(f"stage03: expected one native table on SI page "
                             f"{page_idx + 1}, found {len(t)}.")
        g = t[0].extract()
    body = [[c.strip() for c in r] for r in g[1:]]
    if len(body) != expect_rows:
        raise SystemExit(f"stage03: SI page {page_idx + 1} gave {len(body)} "
                         f"rows, expected {expect_rows}.")
    return body


def main():
    blk = blocks()
    s1 = native(20, 7)     # Table S1, page 21: the 6 leads + MC3
    scr1 = {r[0]: r for r in native(1, 30)}   # Fig S1 Screening 1, page 2

    qc = {"table_s1_rows": len(s1), "screening1_rows": len(scr1),
          "cross_table_agreement": [], "molpct_sums": [],
          "distinct_ionizable_structures": None}

    out = []
    for r in s1:
        name, lip, chol, peg, helper = r[0], r[1], r[2], r[3], r[4]
        pct = [float(x) for x in r[5:9]]
        tot = round(sum(pct), 6)
        qc["molpct_sums"].append({"LNP": name, "sum": tot})
        if abs(tot - 100.0) > 1e-9:
            raise SystemExit(f"stage03: {name} mole percents sum to {tot}, "
                             f"not 100. Refusing.")

        # cross-table check against the page-2 Screening-1 table
        key = name.replace("LNP", "LNP")
        if key in scr1:
            o = scr1[key]
            same = ([o[1], o[2], o[3], o[4]] == [lip, chol, peg, helper]
                    and [float(x) for x in o[5:9]] == pct)
            qc["cross_table_agreement"].append(
                {"LNP": name, "page21_vs_page2": "exact" if same else "MISMATCH",
                 "page2": o[1:], "page21": r[1:]})
            if not same:
                raise SystemExit(f"stage03: {name} differs between SI page 2 "
                                 f"and SI Table S1. Refusing to pick one.")
        else:
            qc["cross_table_agreement"].append(
                {"LNP": name, "page21_vs_page2": "not_in_screening1_table",
                 "why": "MC3 is the benchmark, not a library member"})

        il = TOKEN[lip]
        rec = {
            "formulation_id": name,
            "IL_name": il,
            "IL_SMILES": blk[il]["SMILES"],
            "IL_structure_source": blk[il]["structure_source"],
            "IL_molratio": pct[0],
            "CHL_name": TOKEN[chol], "CHL_SMILES": blk[TOKEN[chol]]["SMILES"],
            "CHL_molratio": pct[1],
            "PEG_name_as_printed": peg,
            "PEG_name_corpus": "C14-PEG2000",
            "PEG_SMILES": blk[TOKEN[peg]]["SMILES"],
            "PEG_molratio": pct[2],
            "HL_name": TOKEN[helper], "HL_SMILES": blk[TOKEN[helper]]["SMILES"],
            "HL_molratio": pct[3],
            "role": "benchmark" if name == "MC3" else "lead_candidate",
            "source": "SI Table S1 (page 21), cross-checked against SI Fig S1 "
                      "Screening 1 (page 2)",
        }
        m = Chem.MolFromSmiles(rec["IL_SMILES"])
        if m is None:
            raise SystemExit(f"stage03: {il} SMILES does not parse.")
        rec["IL_formula"] = rdMolDescriptors.CalcMolFormula(m)
        rec["IL_exact_mass"] = round(Descriptors.ExactMolWt(m), 4)
        rec["IL_inchikey"] = Chem.MolToInchiKey(m)
        rec["IL_skeleton"] = rec["IL_inchikey"].split("-")[0]
        out.append(rec)

    ils = sorted({r["IL_inchikey"] for r in out})
    leads = sorted({r["IL_inchikey"] for r in out if r["role"] == "lead_candidate"})
    qc["distinct_ionizable_structures"] = {
        "all_individually_measured": len(ils),
        "excluding_the_MC3_benchmark": len(leads),
        "inchikeys": ils,
        "statement": "The six lead LNPs are ONE ionizable lipid (cKK-E12) in "
                     "six different formulations. The only second structure "
                     "among individually-administered LNPs is the MC3 "
                     "benchmark.",
    }

    cols = list(out[0])
    with open(os.path.join(HERE, "03_enumerated.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(out)
    with open(os.path.join(HERE, "03_enumerated_qc.json"), "w") as fh:
        json.dump(qc, fh, indent=1)

    print(f"03_enumerated.csv   {len(out)} individually-administered "
          f"formulations")
    print(f"  distinct ionizable structures: {len(ils)} "
          f"({len(leads)} excluding the MC3 benchmark)")
    print(f"  mole percents sum to 100 on {len(qc['molpct_sums'])}/"
          f"{len(out)} formulations")
    ok = sum(1 for c in qc["cross_table_agreement"]
             if c["page21_vs_page2"] == "exact")
    print(f"  SI page 21 vs SI page 2: {ok}/6 leads exact "
          f"(MC3 is not a Screening-1 member)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
