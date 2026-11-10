#!/usr/bin/env python3
"""Stage 06 -- emit the corpus-admissible rows. There are none.

06_rows.csv is written with its header and ZERO data rows. An empty-but-headed
file is the same posture as an empty reference_SMILES: it makes the zero
VISIBLE. A missing file would be indistinguishable from a stage that was never
run.

WHY ZERO, in one sentence each:
  * the delivery values (EGFP from plasmid DNA and from mRNA in HeLa, plus
    fibroblasts and cancer lines) are all in the MAIN ARTICLE, which is behind
    the ACS subscription wall with no PMC deposit -- see
    00_value_sources.json -> acquisition.unobtainable for the checkable
    negatives;
  * every number in the SI we do have is viability, nucleic-acid binding or a
    micrograph, and none of those is in data/'s delivery-only
    Experiment_method vocabulary;
  * and even if a value existed, no structure does: zero of the 16 lipidoids
    are enumerable, so a row would have a value and no molecule.

Any ONE of those three is sufficient. All three hold.

BUT THE EMIT GUARDS ARE STILL EXERCISED HERE, against synthetic items, because
"we would have flagged correctly" is not a claim this folder gets to make for
free. Two guards must RAISE:
  * a hard-coded literal `confidence`  (F-011 / F-025 / F-046, error class 29)
  * a caveat with an empty note        (F-024)

Run:  python3 -I stage06_emit.py
"""
import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))

import emit                                                # noqa: E402
from errors import EmitError                               # noqa: E402

# Identity columns a merge would need, followed by the quality columns, which
# are read FROM THE SCHEMA and never restated here (F-001).
IDENTITY = ["Experiment_ID", "IL_name", "IL_SMILES", "Model_type",
            "Model_target", "Experiment_method", "Experiment_value"]


def guard_literal_confidence():
    """emit_rows must REFUSE a literal confidence."""
    try:
        emit.emit_rows(["x"], {"x": {"conf": "low"}},
                       confidence="high",                 # the defect
                       note=lambda i, r: "a note",
                       label="guard_literal")
    except EmitError as e:
        return True, str(e).splitlines()[0]
    return False, "emit_rows ACCEPTED a literal confidence"


def guard_caveat_without_note():
    """emit_rows must REFUSE a caveat with no note."""
    try:
        emit.emit_rows(["x"], {"x": {"conf": "medium"}},
                       confidence=lambda i, r: r["conf"],
                       note=lambda i, r: "",               # the defect
                       label="guard_note")
    except EmitError as e:
        return True, str(e).splitlines()[0]
    return False, "emit_rows ACCEPTED a caveat with an empty note"


def guard_accepts_a_correct_call():
    """The control. A guard suite with no passing case cannot distinguish
    'refuses everything' from 'refuses the right things'."""
    rows = emit.emit_rows(
        ["x"], {"x": {"conf": "medium", "why": "the reason, on the row"}},
        confidence=lambda i, r: r["conf"],
        note=lambda i, r: r["why"],
        label="guard_control")
    ok = (len(rows) == 1 and rows[0]["confidence"] == "medium"
          and rows[0]["curation_note"] == "the reason, on the row")
    return ok, f"control returned {rows}"


def main():
    admissible = []          # measured in stage 05: nothing is admissible

    cols = IDENTITY + emit.columns()
    with open(os.path.join(HERE, "06_rows.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in admissible:            # loop runs zero times, by measurement
            w.writerow(r)

    guards = {}
    allok = True
    for name, fn in (("refuses_a_literal_confidence", guard_literal_confidence),
                     ("refuses_a_caveat_with_no_note", guard_caveat_without_note),
                     ("accepts_a_correct_per_item_call",
                      guard_accepts_a_correct_call)):
        ok, msg = fn()
        guards[name] = {"fired_as_required": bool(ok), "detail": msg}
        allok = allok and ok

    rec = {
        "rows_emitted": len(admissible),
        "rows_with_a_caveat": emit.caveat_count([]),
        "quality_columns_from_schema": emit.columns(),
        "quality_columns_restated_in_this_file": 0,
        "why_zero_rows": {
            "1_no_admissible_value_in_hand": "every number in the SI is "
                "viability (Fig S3), nucleic-acid binding (Fig S2) or a "
                "micrograph (Fig S1). 0 of 18 rows in 05_values.csv are "
                "corpus_admissible = yes.",
            "2_the_delivery_values_are_unobtainable": "they are in the main "
                "article. ACS subscription; elink gives only pubmed_pmc_refs "
                "with no pubmed_pmc self-link; esearch db=pmc on the DOI "
                "returns count=0; no author repository.",
            "3_no_structure_exists_either": "0 of 16 lipidoids are "
                "enumerable; the SI names, draws and formulates no amine.",
            "any_one_is_sufficient": True,
        },
        "emit_guards": guards,
        "guards_note": "exercised against SYNTHETIC items. This folder emits "
            "nothing, so these demonstrate the guards are live; they do not "
            "demonstrate that this folder's (non-existent) flags are right.",
    }
    with open(os.path.join(HERE, "06_emit_reconciliation.json"), "w") as fh:
        json.dump(rec, fh, indent=2)

    print(f"06_rows.csv: header + {len(admissible)} data rows "
          f"({len(cols)} columns, {len(emit.columns())} from the schema)")
    for k, v in guards.items():
        print(f"  guard {k}: {'OK' if v['fired_as_required'] else 'DID NOT FIRE'}")
    print(f"n_agree={sum(1 for v in guards.values() if v['fired_as_required'])} "
          f"n_total={len(guards)}")
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main())
