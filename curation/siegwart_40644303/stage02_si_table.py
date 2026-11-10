#!/usr/bin/env python3
"""SI Figure S3 (SI p3): the only TEXT table in this paper. Tier 4.

    python3 -I stage02_si_table.py

It carries Diameter (nm), PDI, pK_a and encapsulation efficiency for all 42
in-vivo formulations. None of those is a corpus delivery endpoint -- CLAUDE.md's
in-vivo definition excludes diameter and zeta -- so nothing here is emitted as
an Experiment_value. It is read because it is the ONLY tier-4 artefact in the
paper and it does two jobs no figure can:

  1. It names all 42 lipids and their order, independently of any panel.
  2. Its pK_a column is the x vector that identifies the unlabelled markers in
     Figure 5a, which is how the 42 bar readings get a second, independent
     reading (stage06_scatter_check.py).

PARSING. pymupdf's text layer emits the page column-major: a run of lipid
names, then the literal "Diameter (nm)" followed by one value per lipid, then
"PDI", "pKa", "E.E. (%)" likewise. Two such blocks, 18 lipids then 24. The
parser asserts the counts rather than zipping whatever it found: a short row
would silently shift every value onto the next lipid.
"""
import csv
import os
import re

import pymupdf

HERE = os.path.dirname(os.path.abspath(__file__))
SI = os.path.join(HERE, "source", "ab5c00463_si_001.pdf")
OUT = os.path.join(HERE, "02_si_table_S3.csv")

PAGE = 2                      # 0-based; SI p3
FIELDS = ["Diameter (nm)", "PDI", "pKa", "E.E. (%)"]
NAME_RE = re.compile(r"^\d[A-Z]\d(?:V\d?)?-SC\d{1,2}$")


def main():
    doc = pymupdf.open(SI)
    lines = [ln.strip() for ln in doc[PAGE].get_text().split("\n")]
    lines = [ln for ln in lines if ln]

    rows, i = [], 0
    while i < len(lines):
        if not NAME_RE.match(lines[i]):
            i += 1
            continue
        names = []
        while i < len(lines) and NAME_RE.match(lines[i]):
            names.append(lines[i])
            i += 1
        cols = {}
        for field in FIELDS:
            if i >= len(lines) or lines[i] != field:
                raise SystemExit(
                    f"stage02: after {len(names)} names the next line is "
                    f"{lines[i]!r}, expected {field!r}. REFUSING to zip "
                    "values onto names by position when the layout is not "
                    "what this parser was written against.")
            i += 1
            vals = []
            while i < len(lines) and re.fullmatch(r"-?\d+(\.\d+)?", lines[i]):
                vals.append(lines[i])
                i += 1
            if len(vals) != len(names):
                raise SystemExit(
                    f"stage02: {field} has {len(vals)} values for "
                    f"{len(names)} lipids. REFUSING: a short column shifts "
                    "every value onto the wrong lipid, which is the F-004 "
                    "shape and passed a range check there.")
            cols[field] = vals
        for k, nm in enumerate(names):
            rows.append({"lipid": nm,
                         "diameter_nm": cols["Diameter (nm)"][k],
                         "PDI": cols["PDI"][k],
                         "pKa": cols["pKa"][k],
                         "encapsulation_efficiency_pct": cols["E.E. (%)"][k]})

    if len(rows) != 42:
        raise SystemExit(f"stage02: parsed {len(rows)} lipids, expected 42.")
    seen = [r["lipid"] for r in rows]
    if len(set(seen)) != 42:
        raise SystemExit("stage02: duplicate lipid names in the SI table.")

    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    # the paper's own claims about this table, checked here
    ee = [float(r["encapsulation_efficiency_pct"]) for r in rows]
    dia = [float(r["diameter_nm"]) for r in rows]
    pka = [float(r["pKa"]) for r in rows]
    print(f"wrote {len(rows)} rows to {OUT}")
    print(f"  E.E.  min {min(ee):.2f}  (paper: 'greater than 80% across all "
          f"formulations') -> {'PASS' if min(ee) > 80 else 'FAIL'}")
    print(f"  size  {min(dia):.0f}-{max(dia):.0f} nm (paper: 'ranged from 100 "
          f"to 200 nm') -> {sum(1 for d in dia if 100 <= d <= 200)}/42 inside")
    print(f"  pKa   {min(pka):.2f}-{max(pka):.2f} (paper: acceptable range "
          f"'6.2-7.4')")


if __name__ == "__main__":
    main()
