#!/usr/bin/env python3
"""Stage 07 -- emit corpus rows. Writes 06_rows.csv.

    python3 -I 07_emit.py

THIS PAPER EMITS ZERO ROWS, so 06_rows.csv is written HEADER-ONLY, which makes
the zero visible rather than absent. The header is the live 79-column corpus
schema, read off data/data_vivo.csv so it cannot drift from the corpus, with the
seven quality columns checked against provenance/schema_quality.json.

Why zero: the paper contributes exactly ONE ionizable lipid, 7C1, which is a
polydisperse branched-PEI/C15-epoxide polymer with no single correct SMILES, and
the corpus has 0 rows with a blank IL_SMILES in either pool. DECISIONS D-1 is
OPEN, which by itself blocks the merge.

The emit path is still EXERCISED here, on the five formulations that would have
been emitted had D-1 gone the other way, so that (a) the flags they would carry
are written down now rather than reconstructed later, and (b) the two guards
that matter are made to fire on purpose: a literal `confidence` must raise
(F-011/F-025/F-046, error class 29) and a caveat without a note must raise
(F-024). A guard that is never exercised is a guard nobody has tested.
"""
import csv, json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "curation", "lib"))

import emit                                                    # noqa: E402
from errors import EmitError                                   # noqa: E402

VIVO = os.path.join(ROOT, "data", "data_vivo.csv")

# The five formulations individually measured in BMECs (Fig S3A/S3B, Fig 6E).
# Each record is what the paper itself states, so every flag below is resolved
# FROM THE RECORD and never from a literal.
WOULD_HAVE_EMITTED = {
    "original_7C1": {"il": "7C1", "il_molpct": 80.0, "chol": 0.0, "peg": 20.0,
                     "peg_type": "C14PEG2000", "library": None,
                     "il_structure": None, "value_source": "Fig S3B bar 1 (=1.0 by "
                     "definition: the normaliser)"},
    "screen1_winner": {"il": "7C1", "il_molpct": 60.0, "chol": 37.5, "peg": 2.5,
                       "peg_type": "C16PEG2000", "library": 1,
                       "il_structure": None,
                       "value_source": "Fig S3B bar 2; text '2.2-fold more potent'"},
    "LNP1": {"il": "7C1", "il_molpct": 70.0, "chol": 10.0, "peg": 20.0,
             "peg_type": "C16PEG2000", "library": 2, "il_structure": None,
             "value_source": "Fig S3B bar 3 (figure only)"},
    "BM1": {"il": "7C1", "il_molpct": 80.0, "chol": 2.5, "peg": 17.5,
            "peg_type": "C18PEG2000", "library": 2, "il_structure": None,
            "value_source": "Fig S3B bar 4; Fig 6E; text '37% silencing', '4.8x'"},
    "LNP3": {"il": "7C1", "il_molpct": 80.0, "chol": 0.0, "peg": 20.0,
             "peg_type": "C18PEG2000", "library": 2, "il_structure": None,
             "value_source": "Fig S3B bar 5 (figure only)"},
}


def confidence_of(item_id, rec):
    """From the record: no structure was resolved for this lipid at all."""
    return "low" if rec["il_structure"] is None else "high"


def note_of(item_id, rec):
    return (f"{rec['il']} is a polydisperse branched-PEI(Mn~600)/C15-epoxide "
            f"polymer; the paper draws only one representative arm, so no single "
            f"SMILES is correct. Formulation {rec['il_molpct']}/{rec['chol']}/"
            f"{rec['peg']} mol% with {rec['peg_type']}; value from "
            f"{rec['value_source']}.")


def structure_assumed_of(item_id, rec):
    return "1" if rec["il_structure"] is None else "0"


def norm_unit_of(item_id, rec):
    # One z per experiment, pooling marrow / lung / heart; the organ difference
    # belongs in Model_target, not in a second standardisation (error class 22).
    return "experiment"


def header():
    with open(VIVO) as f:
        cols = next(csv.reader(f))
    q = [c["name"] for c in emit.schema()["columns"]]
    missing = [c for c in q if c not in cols]
    if missing:
        raise RuntimeError(f"quality columns missing from the corpus header: {missing}")
    return cols, q


def guard_selftests():
    """Make the two emit guards fail on purpose. Both must raise."""
    items = list(WOULD_HAVE_EMITTED)
    out = {}

    try:
        emit.emit_rows(items, WOULD_HAVE_EMITTED, confidence="low",
                       note=note_of, structure_assumed=structure_assumed_of,
                       value_normalization_unit=norm_unit_of, label="seed-literal")
        out["literal_confidence_raises"] = {"raised": False, "credited": False}
    except EmitError as e:
        out["literal_confidence_raises"] = {
            "raised": True, "credited": "LITERAL" in str(e),
            "message_head": str(e).splitlines()[0]}

    try:
        emit.emit_rows(items, WOULD_HAVE_EMITTED, confidence=confidence_of,
                       note=lambda i, r: "", structure_assumed=structure_assumed_of,
                       value_normalization_unit=norm_unit_of, label="seed-no-note")
        out["caveat_without_a_note_raises"] = {"raised": False, "credited": False}
    except EmitError as e:
        out["caveat_without_a_note_raises"] = {
            "raised": True, "credited": "curation_note is empty" in str(e),
            "message_head": str(e).splitlines()[0]}

    rows = emit.emit_rows(items, WOULD_HAVE_EMITTED, confidence=confidence_of,
                          note=note_of, structure_assumed=structure_assumed_of,
                          value_normalization_unit=norm_unit_of,
                          label="control")
    out["control_passes"] = {
        "n_rows": len(rows), "rows_with_a_caveat": emit.caveat_count(rows),
        "credited": len(rows) == 5 and emit.caveat_count(rows) == 5,
        "example": rows[0]}
    return out


def main():
    cols, q = header()
    with open(os.path.join(HERE, "06_rows.csv"), "w", newline="") as f:
        csv.writer(f).writerow(cols)

    rep = {
        "_what": "Stage 07. 06_rows.csv is HEADER-ONLY: this paper emits zero "
                 "corpus rows. Produced by 07_emit.py.",
        "rows_emitted": 0,
        "header_columns": len(cols),
        "header_source": "data/data_vivo.csv, read live",
        "quality_columns_present": q,
        "why_zero": [
            "1 distinct ionizable lipid (7C1), which is a polydisperse polymer "
            "with no single correct SMILES; data/ has 0 blank IL_SMILES rows.",
            "The 146 pooled-and-measured formulations carry no attributable "
            "value: Fig 3E/5D are unlabelled scatters and the stable subsets "
            "(115 of 120, 31 of 40) are never identified.",
            "The pooled readout is DNA barcode recovery, not expression; the "
            "corpus's one barcoded experiment (LR_2023) is a translated "
            "peptide-barcode readout.",
            "DECISIONS.md D-1 is OPEN, which blocks the merge on its own.",
        ],
        "would_have_emitted_if_D1_resolved": {
            "n_formulations": len(WOULD_HAVE_EMITTED),
            "pool": "in vivo",
            "Experiment_ID": "CS_2018",
            "Model_target": "bone_marrow / lung / heart, carried in Model_target "
                            "with ONE z over the experiment (error class 22)",
            "Experiment_method": "protein_adundance_normalized (Fig S3B, 6E, 6G) "
                                 "and editing_efficiency_normalized (Fig 6H)",
            "Experiment_batching": "individual",
        },
        "guard_selftests": guard_selftests(),
    }
    with open(os.path.join(HERE, "06_rows_report.json"), "w") as f:
        json.dump(rep, f, indent=1)
    print(json.dumps(rep, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
