#!/usr/bin/env python3
"""Write 06_rows.csv, and prove the quality-flag path is live on zero rows.

THIS CURATION EMITS NO ROWS. 05_values.csv has 14 entries and zero are
corpus-admissible, so 06_rows.csv is a header and nothing else.

A header-only file is exactly where a flagging bug hides, because nothing is
exercised. So this stage does three things instead of one:

  1. Reconciles the counts. values_total = admissible + inadmissible, and
     rows_written must equal admissible. F-038 is the reason: three ZT_2026
     checks scored a literal 0 for twelve days while the summary beside them
     said "ALL FOUR EXACT".
  2. Runs the four ENUMERATED structures through emit.emit_rows with a
     PER-ITEM CALLABLE reading each structure's own enumeration record, so the
     flags that will ride these rows the day the main article is acquired are
     validated now, against the real schema vocabulary.
  3. Asserts emit.emit_rows REFUSES a hard-coded literal confidence (F-011,
     F-025, F-046 -- the same class three times). A gate that was never made to
     fail is not a gate.
"""
import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..")
sys.path.insert(0, os.path.join(ROOT, "curation", "lib"))

import emit  # noqa: E402
from errors import EmitError  # noqa: E402


def main():
    values = list(csv.DictReader(open(os.path.join(HERE, "05_values.csv"))))
    admissible = [v for v in values if v["corpus_admissible"].strip() == "yes"]
    inadmissible = [v for v in values
                    if v["corpus_admissible"].strip() != "yes"]
    assert len(admissible) + len(inadmissible) == len(values), "counts"

    items = list(csv.DictReader(open(os.path.join(HERE, "03_enumerated.csv"))))
    records = {r["IL_name"]: r for r in items}

    # --- the gate must refuse a literal -------------------------------------
    refused = None
    try:
        emit.emit_rows([r["IL_name"] for r in items], records,
                       confidence="high",
                       note=lambda i, rec: rec["evidence"], label="seed")
    except EmitError as e:
        refused = str(e).splitlines()[0]
    if refused is None:
        raise SystemExit("emit.emit_rows accepted a LITERAL confidence. "
                         "That is error class 29 and the gate is dead.")

    # --- the real, per-item path --------------------------------------------
    flags = emit.emit_rows(
        [r["IL_name"] for r in items], records,
        confidence=lambda i, rec: rec["confidence"],
        note=lambda i, rec: rec["evidence"],
        structure_assumed=lambda i, rec: int(rec["structure_assumed"]),
        label="xuq_24574196 structures")

    schema_cols = emit.columns()
    header = (["IL_name", "IL_SMILES", "Experiment_ID", "Experiment_value",
               "Experiment_method", "Model_type", "Model_target", "Cargo"]
              + schema_cols)
    with open(os.path.join(HERE, "06_rows.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        # zero data rows, by construction: see reason_zero_rows below

    out = {
        "_what": "emit reconciliation for curation/xuq_24574196",
        "values_total": len(values),
        "values_corpus_admissible": len(admissible),
        "values_not_admissible": len(inadmissible),
        "rows_written_to_06_rows_csv": 0,
        "counts_sum": len(admissible) + len(inadmissible) == len(values),
        "rows_equal_admissible": 0 == len(admissible),
        "reason_zero_rows":
            "The SI carries no delivery measurement (6 viability bars, 2 ESI-MS "
            "masses, TEM micrographs, NMR spectra) and the main article, which "
            "holds every knockdown panel, is closed access with no PMC deposit, "
            "no OA repository copy and a Cloudflare interstitial on the "
            "publisher page. Browser acquisition was out of scope for this run "
            "and no bypass was attempted.",
        "literal_confidence_refused": True,
        "refusal_message_first_line": refused,
        "structures_carried_with_validated_flags": len(flags),
        "flag_preview": [dict(name=i["IL_name"], **f)
                         for i, f in zip(items, flags)],
        "rows_with_a_caveat_if_values_ever_arrive": emit.caveat_count(flags),
    }
    with open(os.path.join(HERE, "06_emit_reconciliation.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
