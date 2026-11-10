#!/usr/bin/env python3
"""Emit the rows in the corpus's 79-column schema. Flags go ON THE ROW.

    python3 -I stage07_emit.py

THE ROWS ARE NOT MERGEABLE YET, BY DESIGN. DECISIONS.md carries OPEN entries
-- above all D-1, the identity of the assay behind Fig. 1C and Fig. 1E -- and
`build_dataset.py` refuses a folder with an OPEN entry. This file exists so
that everything that IS settled is written down in the corpus's own shape, and
so the only thing standing between here and a merge is a question a human has
to answer, not work that has not been done.

WHY EVERY ROW CARRIES `confidence = medium`. Not because any single feature is
doubtful, but because the whole library's connectivity rests on two drawings
and this paper publishes nothing that could test them: no HRMS, no NMR table,
no SMILES/InChI, no deposited structure file. The two bonds at risk --
where the thiolactone opens and which epoxide carbon takes the sulfur -- are
invisible to molecular formula, so a reader who sees only the formula would
have no way to know. `confidence` is derived per lipid from its own
enumeration record, never written as a literal: a hard-coded literal is error
class 29, found three times in this project (F-011, F-025, F-046), and in
F-011 it inverted the labelling so that the most uncertain rows carried the
strongest claim.

WHY NO ROW CARRIES `value_is_bound`. The values are the authors' own printed
numbers, not a digitised heatmap: there is no colormap floor to sit on and no
axis to be censored past. A bound flag here would be false.
"""
import csv
import json
import math
import os
import statistics
import sys

LIB = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib")
sys.path.insert(0, LIB)
import emit  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
CORPUS_HEADER = os.path.join(ROOT, "data", "data_vitro.csv")

EXPERIMENT_ID = "QG_2026"
DOI = "10.1021/acsnano.6c00706"
PMID = "42130331"
LINK = "https://doi.org/10.1021/acsnano.6c00706"

# D-2 is OPEN. The default below is the corpus convention -- ONE z-score per
# Experiment_ID -- and it is the choice that does NOT delete the ~50x
# difference between the two libraries, which is the paper's central
# structure-activity claim. Splitting Library 1 and Library 2 into two
# experiments would standardise that difference away, which is what happened
# in QL_2026 before the 1.1116x offset was found. The raw values stay in
# 05_values.csv so the other choice is one line.
POOL_LIBRARIES = True


def header():
    with open(CORPUS_HEADER) as fh:
        return next(csv.reader(fh))


def main():
    vals = list(csv.DictReader(open(os.path.join(HERE, "05_values.csv"))))
    enum = {r["IL_name"]: r for r in csv.DictReader(
        open(os.path.join(HERE, "03_enumerated.csv")))}

    groups = {}
    for v in vals:
        key = "" if POOL_LIBRARIES else v["library"]
        groups.setdefault(key, []).append(v)

    z = {}
    norm = {}
    for key, rows in groups.items():
        logs = [math.log10(float(r["value"])) for r in rows]
        mu = statistics.fmean(logs)
        sd = statistics.pstdev(logs)          # population sd, ddof=0
        if sd == 0:
            raise SystemExit("zero spread; refusing to divide by it.")
        for r, lg in zip(rows, logs):
            z[r["IL_name"]] = (lg - mu) / sd
        norm[key] = {"n": len(rows), "mean_log10": mu, "sd_log10_pop": sd}

    names = [v["IL_name"] for v in vals]
    records = {n: enum[n] for n in names}

    quality = emit.emit_rows(
        names, records,
        confidence=lambda n, rec: rec["confidence"],
        note=lambda n, rec: (
            "beta-hydroxy thioether lipid enumerated from the Fig. S1A scheme "
            "and the Fig. S1B / Fig. S2 block drawings; the paper publishes no "
            "HRMS, NMR table, SMILES or structure file, so the thiolactone "
            "ring-opening regiochemistry and the thiol-epoxide regiochemistry "
            "are read from the drawing and have no independent check. Value is "
            "the authors' own printed number from SI Table "
            + ("S1 (Fig. 1C)" if rec["library"] == "1" else "S2 (Fig. 1E)")
            + "; the assay behind that figure is not stated anywhere in the SI "
              "and the main text is paywalled -- see DECISIONS.md D-1."),
        structure_assumed=lambda n, rec: "0",
        positional_ambiguity=lambda n, rec: "0",
        value_is_bound=lambda n, rec: "",
        value_normalization_unit=lambda n, rec: "experiment",
        within_model=lambda n, rec: round(z[n], 6),
        id_key="IL_name", label="QG_2026")

    cols = header()
    out = []
    for v, q in zip(vals, quality):
        row = {c: "" for c in cols}
        row.update({
            "Experiment_ID": EXPERIMENT_ID,
            "IL_name": v["IL_name"],
            "IL_SMILES": enum[v["IL_name"]]["IL_SMILES"],
            "Model": "",            # D-1: in vitro vs in vivo not stated
            "Model_type": "",       # D-1
            "Model_target": "",     # D-1
            "Cargo": "",            # D-1
            "Cargo_type": "",
            "Experiment_method": "",
            "Experiment_value": round(z[v["IL_name"]], 6),
            "Publication_link": LINK,
            "Publication_DOI": DOI,
            "Publication_PMID": PMID,
        })
        row.update(q)
        out.append(row)

    path = os.path.join(HERE, f"07_{EXPERIMENT_ID}_rows.csv")
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(out)

    summary = {
        "experiment_id": EXPERIMENT_ID,
        "rows": len(out),
        "columns": len(cols),
        "pooled_z_across_both_libraries": POOL_LIBRARIES,
        "normalisation": norm,
        "rows_with_a_caveat": emit.caveat_count(out),
        "z_mean": round(statistics.fmean(
            float(r["Experiment_value"]) for r in out), 12),
        "z_sd_pop": round(statistics.pstdev(
            float(r["Experiment_value"]) for r in out), 12),
        "blank_Model_fields": "Model, Model_type, Model_target, Cargo and "
                              "Experiment_method are blank because the SI "
                              "never states the assay behind Fig. 1C/1E and "
                              "the main text is paywalled (DECISIONS.md D-1). "
                              "A blank is 'not stated', not 'clean'.",
    }
    with open(os.path.join(HERE, "07_emit_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1)
    for k, v in summary.items():
        print(f"{k}: {v}")
    print(f"\nwrote {path}")


if __name__ == "__main__":
    main()
