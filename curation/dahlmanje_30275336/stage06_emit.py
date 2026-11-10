"""Emit the corpus rows for this paper -- of which there are ZERO -- and PROVE
that the emit guards fire.

A header-only 06_rows.csv is indistinguishable from an emitter that was never
run. So this stage does two things:

  1. writes 06_rows.csv with the real schema header and no data rows;
  2. SEEDS the two defects emit.emit_rows exists to make unwritable, and
     requires each to RAISE. A check that cannot fail is not a check, and a
     clean line from an emitter that emitted nothing is the cheapest possible
     false reassurance.
"""
import sys, os, csv, json
sys.path.insert(0, os.path.join(sys.argv[1], "curation", "lib"))
import emit                                          # noqa: E402
from errors import EmitError                         # noqa: E402

qc = {"n_rows_emitted": 0,
      "why_zero": "No compound in this paper has a determinable structure "
                  "(03_enumerated.csv: 8 compounds, 0 with a SMILES) and the "
                  "corpus has 0 blank IL_SMILES in 30,998 rows. See D-2.",
      "guard_tests": {}}

# A real record, so the guard tests are about the GUARD and not about a missing
# record. This item would be the 7C3 splenic-EC value if 7C1 had a structure.
REC = {"7C3": {"compound": "7C1", "amine": "PEI600",
               "why": "polydisperse branched polymer, degree of substitution undisclosed"}}

def conf_from_record(item_id, record):
    """Per-item, derived from the record -- never a literal (F-011/F-025/F-046)."""
    return "low" if "polydisperse" in record["why"] else "high"

def note_from_record(item_id, record):
    return ("Structure not emitted: %s is a %s conjugate and %s."
            % (record["compound"], record["amine"], record["why"]))

# --- test 1: the emitter WORKS when driven correctly -----------------------
ok = emit.emit_rows(["7C3"], REC, confidence=conf_from_record,
                    note=note_from_record, label="smoke")
qc["guard_tests"]["callable_confidence_accepted"] = {
    "raised": False, "confidence": ok[0]["confidence"],
    "note_nonempty": bool(ok[0]["curation_note"].strip()),
    "caveat_count": emit.caveat_count(ok)}

# --- test 2: a LITERAL confidence must raise (F-011) -----------------------
try:
    emit.emit_rows(["7C3"], REC, confidence="high", note=note_from_record,
                   label="seed_literal_confidence")
    qc["guard_tests"]["literal_confidence_raises"] = {"raised": False,
        "VERDICT": "FAIL - the guard did not fire"}
except EmitError as e:
    qc["guard_tests"]["literal_confidence_raises"] = {"raised": True,
        "exception": "EmitError", "message_head": str(e).split("\n")[0],
        "VERDICT": "PASS"}

# --- test 3: a caveat with an EMPTY note must raise (F-024) ----------------
try:
    emit.emit_rows(["7C3"], REC, confidence=conf_from_record,
                   note=lambda i, r: "", label="seed_caveat_without_note")
    qc["guard_tests"]["caveat_without_note_raises"] = {"raised": False,
        "VERDICT": "FAIL - the guard did not fire"}
except EmitError as e:
    qc["guard_tests"]["caveat_without_note_raises"] = {"raised": True,
        "exception": "EmitError", "message_head": str(e).split("\n")[0],
        "VERDICT": "PASS"}

# --- test 4: an out-of-vocabulary flag must raise (C13) --------------------
try:
    emit.emit_rows(["7C3"], REC, confidence=lambda i, r: "probably",
                   note=note_from_record, label="seed_bad_vocabulary")
    qc["guard_tests"]["out_of_vocabulary_raises"] = {"raised": False,
        "VERDICT": "FAIL - the guard did not fire"}
except EmitError as e:
    qc["guard_tests"]["out_of_vocabulary_raises"] = {"raised": True,
        "exception": "EmitError", "message_head": str(e).split("\n")[0],
        "VERDICT": "PASS"}

# --- the actual (empty) emission -------------------------------------------
CORE = ["Experiment_ID", "IL_name", "IL_SMILES", "HL_name", "HL_molratio",
        "CHL_name", "CHL_molratio", "PEG_name", "PEG_molratio",
        "Model", "Model_type", "Model_target", "Route_of_administration",
        "Cargo", "Cargo_type", "Dose_ug_nucleicacid", "Experiment_method",
        "Experiment_batching", "Experiment_value",
        "Publication_DOI", "Publication_PMID"]
header = CORE + emit.columns()
with open("06_rows.csv", "w", newline="") as fh:
    csv.writer(fh).writerow(header)

qc["06_rows_header"] = header
qc["n_quality_columns_from_schema"] = len(emit.columns())
print(json.dumps(qc, indent=1))
