#!/usr/bin/env python3
"""Turn the decoded table cells into 05_values.csv, and check them.

    python3 -I stage05_values.py

Table S1 is the source data for Fig. 1C (Library 1, 10 amines x 5 epoxides).
Table S2 is the source data for Fig. 1E (Library 2, the same 50 lipids x 5
acyl chlorides). Both are the AUTHORS' OWN numbers -- the captions say so --
so nothing here is estimated from a colour or a bar height. What could still
go wrong is the mapping from a cell to a compound, and that is what the checks
below are for.

THE CHECK THIS CURATION LEANS ON. The paper names 113-AA-C8C14 as its lead
compound, in its title-level claim and throughout the SI figure captions. That
name decomposes to Table S2's column `113-AA-C8` and row `Cl14`. If our grid
mapping were off by a row or a column, that cell would be an ordinary value.
It is the MAXIMUM of all 250, and the maximum is unique. We did not fit
anything to make that true, and nothing in the decoding knew which cell was
supposed to win.
"""
import csv
import json
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))

# Table S2's row labels are `Cl<n>`; the acyl block in a Library 2 name is
# `C<n>`. The paper writes both, e.g. Fig. S3B's axis reads `70-AA-C8C6`
# while Table S2's row reads `Cl6`.
ROW_TO_ACYL = {"Cl6": "C6", "Cl8": "C8", "Cl10": "C10",
               "Cl12": "C12", "Cl14": "C14"}


def load():
    with open(os.path.join(HERE, "04_glyphs.json")) as fh:
        return json.load(fh)["cells"]


def num(s):
    s = s.strip()
    if not re.fullmatch(r"\d+(\.\d+)?([Ee]\+\d+)?", s):
        raise SystemExit(f"cell {s!r} is not a number; refusing to coerce it.")
    return float(s)


def main():
    cells = load()
    rows = []

    # --- Table S1 -------------------------------------------------------
    amines = [cells[f"S1|0|{c}"] for c in range(1, 11)]
    for r in range(1, 6):
        epox = cells[f"S1|{r}|0"]
        for c in range(1, 11):
            raw = cells[f"S1|{r}|{c}"]
            rows.append({
                "library": 1, "table": "S1", "figure_panel": "Fig 1C",
                "cell": f"row {r} ({epox}) x col {c} ({amines[c - 1]})",
                "IL_name": f"{amines[c - 1]}-AA-{epox}",
                "amine": amines[c - 1], "epoxide": epox, "acyl": "",
                "value_as_printed": raw, "value": num(raw),
            })

    # --- Table S2 -------------------------------------------------------
    s2_rows = sorted({int(k.split("|")[1]) for k in cells
                      if k.startswith("S2|")})
    headers = [r for r in s2_rows if f"S2|{r}|0" not in cells]
    for h in headers:
        cols = [cells[f"S2|{h}|{c}"] for c in range(1, 11)]
        for r in range(h + 1, h + 6):
            lab = cells[f"S2|{r}|0"]
            if lab not in ROW_TO_ACYL:
                raise SystemExit(f"unexpected Table S2 row label {lab!r}")
            for c in range(1, 11):
                raw = cells[f"S2|{r}|{c}"]
                rows.append({
                    "library": 2, "table": "S2", "figure_panel": "Fig 1E",
                    "cell": f"row {r} ({lab}) x col {c} ({cols[c - 1]})",
                    "IL_name": cols[c - 1] + ROW_TO_ACYL[lab],
                    "amine": cols[c - 1].split("-")[0],
                    "epoxide": cols[c - 1].split("-")[-1],
                    "acyl": ROW_TO_ACYL[lab],
                    "value_as_printed": raw, "value": num(raw),
                })

    # --- checks ---------------------------------------------------------
    enum = {r["IL_name"]: r for r in csv.DictReader(
        open(os.path.join(HERE, "03_enumerated.csv")))}
    names = [r["IL_name"] for r in rows]
    missing = sorted(set(names) - set(enum))
    extra = sorted(set(enum) - set(names))
    dup = len(names) - len(set(names))

    lib2 = [r for r in rows if r["library"] == 2]
    top = max(lib2, key=lambda r: r["value"])
    ties = [r for r in lib2 if r["value"] == top["value"]]

    checks = {
        "n_values": len(rows),
        "n_library_1": sum(1 for r in rows if r["library"] == 1),
        "n_library_2": len(lib2),
        "duplicate_names": dup,
        "values_with_no_enumerated_structure": missing,
        "enumerated_structures_with_no_value": extra,
        "lead_compound_named_by_the_paper": "113-AA-C8C14",
        "library_2_maximum": {"name": top["IL_name"],
                              "value": top["value"],
                              "cell": top["cell"],
                              "is_unique": len(ties) == 1},
        "lead_is_the_maximum": (top["IL_name"] == "113-AA-C8C14"
                                and len(ties) == 1),
        "library_1_range": [min(r["value"] for r in rows
                                if r["library"] == 1),
                            max(r["value"] for r in rows
                                if r["library"] == 1)],
        "library_2_range": [min(r["value"] for r in lib2),
                            max(r["value"] for r in lib2)],
        "library_1_median": sorted(r["value"] for r in rows
                                   if r["library"] == 1)[25],
        "library_2_median": sorted(r["value"] for r in lib2)[125],
    }
    for k in ("values_with_no_enumerated_structure",
              "enumerated_structures_with_no_value"):
        if checks[k]:
            raise SystemExit(f"{k}: {checks[k]}")
    if dup:
        raise SystemExit(f"{dup} duplicate IL_names among the values.")
    if not checks["lead_is_the_maximum"]:
        raise SystemExit(
            "The paper's own lead compound, 113-AA-C8C14, is NOT the unique "
            f"maximum of Table S2 -- the maximum is {top['IL_name']} at "
            f"{top['value']}. Either the grid mapping is off by a row or a "
            "column, or the lead was not chosen on this readout. Stop and "
            "find out which; do not emit rows on an unverified mapping.")

    with open(os.path.join(HERE, "05_values.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    with open(os.path.join(HERE, "05_value_checks.json"), "w") as fh:
        json.dump(checks, fh, indent=1)
    for k, v in checks.items():
        print(f"{k}: {v}")


if __name__ == "__main__":
    main()
