#!/usr/bin/env python3
"""Collapse the FALCON curation into the single block that merges into data/.

Decision of 2026-10-29, revised the same day. ALL 2,544 measurements merge, on
all seven cell lines, under ONE Experiment_ID `WT_2026`: the corpus treats a
source as a paper, and different cell lines from the same paper are one source.

An earlier version held back the five cell lines new to the corpus (Ramos,
THP-1, INS-1E, NIH-3T3, C2C12) to avoid adding Model_type categories. That was
dropped as indefensible -- a new cell line is new information, not a reason to
discard 1,851 real measurements, and Model_type is a free-text category that
already spans eleven values.

Z-SCORE. The SEVEN contributing blocks -- one per source file, corrected from
fourteen (per file x cell type) on 2026-10-30, review finding F-003 -- were
normalized by the authors against their own per-file internal controls, and
`dataset_info.txt` warns that values "cannot be cross compared across files".
Pooling the raw Ln(nLE) and z-scoring once would fold each file's normalization
offset into the values as if it were signal. Instead each block keeps the
z-score computed inside its own normalization unit, and those are concatenated.

That still satisfies the corpus convention exactly, and not by approximation:
for blocks each having mean 0 and population sd 1, the pooled mean is
sum(n_i * 0)/N = 0 and the pooled second moment is sum(n_i * 1)/N = 1. So the
combined WT_2026 block sits at mean 0 / sd 1 with no cross-file offset mixed in.
Asserted below rather than assumed.

BUT NOTE WHAT THAT ASSERTION IS AND IS NOT. It is the arithmetic above, so it
holds for ANY set of standardized blocks, including wrongly chosen ones: it held
just as well when there were fourteen blocks and every within-file
cross-cell-type difference had been deleted. It is a statement about
concatenation, not evidence that the normalization unit is right. The check that
actually tests the unit is in `build_rows.py`, where the stored z must reproduce
the authors' printed `nLnLE_Selectivity` column up to one constant per file --
that one fails if the unit is ever changed back.
"""
import csv
import json
import os
import statistics as st

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "01_WT_2026_rows.csv")
OUT = os.path.join(HERE, "05_WT_2026_merge.csv")

EID = "WT_2026"


def main():
    rows = list(csv.DictReader(open(SRC)))
    hdr = list(rows[0])

    keep, held = rows, []

    contributing = {}
    for r in keep:
        contributing.setdefault(r["Experiment_ID"], []).append(r)

    # each contributing block must already be standardized
    for eid, block in contributing.items():
        v = [float(r["Experiment_value"]) for r in block]
        assert abs(st.mean(v)) < 1e-9, (eid, st.mean(v))
        assert abs(st.pstdev(v) - 1) < 1e-9, (eid, st.pstdev(v))

    for r in keep:
        r["Experiment_ID"] = EID
        r["Formulation_ID"] = f"{EID}_{r['Formulation_ID']}"

    pooled = [float(r["Experiment_value"]) for r in keep]
    assert abs(st.mean(pooled)) < 1e-9, st.mean(pooled)
    assert abs(st.pstdev(pooled) - 1) < 1e-9, st.pstdev(pooled)

    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=hdr)
        w.writeheader()
        w.writerows(keep)

    qc = {
        "merged_rows": len(keep),
        "experiment_id": EID,
        "cells_merged": sorted({r["Model_type"] for r in keep}),
        "contributing_blocks": {k: len(v) for k, v in contributing.items()},
        "pooled_z": {"mean": st.mean(pooled), "pop_sd": st.pstdev(pooled)},
        "doses_ug": sorted({(r["Model_type"], r["Dose_ug_nucleicacid"]) for r in keep}),
        "rows_per_cell": {c: sum(1 for r in keep if r["Model_type"] == c)
                          for c in sorted({r["Model_type"] for r in keep})},
        "new_model_types": sorted({"Ramos", "THP-1", "INS-1E", "NIH-3T3", "C2C12"}),
    }
    json.dump(qc, open(os.path.join(HERE, "06_merge_subset_qc.json"), "w"), indent=1)

    print(f"{os.path.basename(OUT)}: {len(keep)} rows, one block {EID}")
    print(f"  cells {qc['cells_merged']}  doses {qc['doses_ug']}")
    print(f"  from {qc['contributing_blocks']}")
    print(f"  pooled z mean {st.mean(pooled):.2e} pop sd {st.pstdev(pooled):.12f}")
    print(f"  rows per cell {qc['rows_per_cell']}")


if __name__ == "__main__":
    main()
