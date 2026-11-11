#!/usr/bin/env python3
"""Build data/excluded_endpoints.csv -- the rows the corpus deliberately omits.

WHY THIS EXISTS (policy P6, 2026-11-10). The corpus holds one kind of number:
whether the cargo reached the cytosol and acted. An a-priori, metric-blind rule
of 2026-07-23 therefore excluded `uptake`, `diameter` and `zeta_potential` --
**671 rows, all of `LX_2024`** -- and those rows have had no home since. They
are the whole of the `-585` in-vivo delta between `Original_LNPDB.csv` (2,388
in-vivo rows) and `data/data_vivo.csv` (1,803), and until now the only record
of them was a sentence in a changelog.

The user's instruction was to keep rather than delete them: *"keep them in a
separate file in case they are ever useful in the future."* Two reasons that is
right beyond convenience:

  * **An exclusion nobody can inspect is an assertion, not a method.** Methods
    will state the delivery-only rule; a reader should be able to see exactly
    what it removed, and count it.
  * **A deleted row cannot be re-examined when the rule changes.** `uptake` is
    excluded because it measures where the particle went rather than whether it
    worked -- which is a judgement about the research question, not about data
    quality, and a later paper may want precisely that.

This file is NOT part of the corpus. Nothing reads it, `audit.py` does not gate
it, and it must never be concatenated with `data_vitro.csv` or
`data_vivo.csv` -- every row in it is here because it is not comparable to
them. The `excluded_because` column says which rule removed each row.

Run: python3 -I provenance/build_excluded.py
"""
import csv, json, os, sys
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

# The a-priori delivery-only rule of 2026-07-23, by Experiment_method.
# Each entry must say WHY, in the terms of the research question.
RULES = {
    "uptake": "measures where the particle went, not whether the cargo acted "
              "-- a cell can internalise an LNP whose cargo never escapes the "
              "endosome (a-priori rule, 2026-07-23)",
    "diameter": "a physicochemical property of the particle, not a delivery "
                "outcome (a-priori rule, 2026-07-23)",
    "zeta_potential": "a physicochemical property of the particle, not a "
                      "delivery outcome (a-priori rule, 2026-07-23)",
    "hemolysis_percent": "a toxicity readout, not delivery (a-priori rule, "
                         "2026-07-23; 90 rows, excluded at the original trim)",
}


def main():
    orig = os.path.join(ROOT, "data", "Original_LNPDB.csv")
    vitro = os.path.join(ROOT, "data", "data_vitro.csv")
    if not os.path.exists(orig):
        print(f"missing {orig}", file=sys.stderr)
        return 1

    with open(vitro, newline="") as f:
        cols = next(csv.reader(f))

    rows_in = list(csv.DictReader(open(orig)))
    out, why = [], Counter()
    for r in rows_in:
        m = (r.get("Experiment_method") or "").strip()
        if m not in RULES:
            continue
        row = {c: r.get(c, "") for c in cols}
        # Descriptors are NOT recomputed. These rows are not corpus rows and
        # recomputing would imply they had been through the corpus pipeline.
        row["excluded_because"] = RULES[m]
        row["excluded_rule"] = m
        out.append(row)
        why[m] += 1

    if not out:
        print("NO EXCLUDED ROWS FOUND. Refusing to write an empty file as if "
              "it were a result -- the delivery-only rule removed 761 rows "
              "(671 in-vivo + 90 hemolysis) and finding zero means the method "
              "names or the input have changed.", file=sys.stderr)
        return 1

    hdr = cols + ["excluded_rule", "excluded_because"]
    p = os.path.join(ROOT, "data", "excluded_endpoints.csv")
    with open(p, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=hdr, extrasaction="ignore")
        w.writeheader()
        w.writerows(out)
    print(f"wrote {os.path.relpath(p, ROOT)}: {len(out)} rows x {len(hdr)} cols")
    for k, n in why.most_common():
        print(f"  {k:20s} {n:5d}")
    eids = Counter(r["Experiment_ID"] for r in out)
    print(f"  experiments: {dict(eids)}")
    # The arithmetic this file is supposed to make inspectable.
    print(f"\n  Original_LNPDB in-vivo rows            2388")
    print(f"  minus the {sum(n for k, n in why.items() if k != 'hemolysis_percent')}"
          f" non-delivery in-vivo rows   -> 1717")
    print(f"  plus 86 added by our own curations      -> 1803 = data_vivo.csv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
