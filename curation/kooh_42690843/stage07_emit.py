#!/usr/bin/env python3
"""Quality columns for the nine structures, via curation/lib emit.emit_rows.

WHY THERE IS NO CORPUS ROW HERE. A corpus row needs an `Experiment_value`, and
this paper's values are all in the main article, which we do not hold
(00_value_sources.json). So this stage does the half it CAN do honestly: it
produces the seven quality/provenance columns of
provenance/schema_quality.json for each of the nine structures, written by the
one audited writer, so that when a value source arrives the flags are already
attached to the structure rather than sitting in a QC file.

"A flag that lives only in a per-folder QC file is not a flag, it is a memo"
(F-001, 1,424 corpus rows). This file exists so that sentence is not true of
this folder.

emit.emit_rows REFUSES a literal confidence (error class 29, found three times:
F-011, F-025, F-046), so `confidence` and `curation_note` below are callables
reading each lipid's OWN record out of 03_enumerated.csv.

Run:  python3 -I stage07_emit.py                writes 07_rows_quality.csv
      python3 -I stage07_emit.py --seed-defect  passes a literal "high" and
                                                confirms emit_rows refuses it.
"""
import csv
import os
import sys

LIB = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib")
sys.path.insert(0, LIB)
import emit                                        # noqa: E402
from errors import EmitError                       # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "03_enumerated.csv")
OUT = os.path.join(HERE, "07_rows_quality.csv")


def records():
    with open(SRC) as fh:
        return {r["lipid"]: r for r in csv.DictReader(fh)}


def main():
    recs = records()
    items = list(recs)

    def confidence(item_id, rec):
        # the lipid's OWN assessment, carried from its enumeration record
        return rec["confidence"]

    def note(item_id, rec):
        n = int(rec["unassigned_stereocentres"])
        return ("Alpha carbon%s formed by the 427 nm photo-Petasis coupling "
                "%s %s unassigned stereocentre%s: the reaction is not "
                "stereoselective and the SI reports no resolution, so the "
                "stored SMILES is the constitution only. Regiochemistry "
                "(alkyl on the alpha carbon, not the nitrogen) is pinned by "
                "the authors' own IUPAC name and the 1H NMR N-CH methine, "
                "not by HRMS, which cannot see it."
                % ("s" if n != 1 else "", "are" if n != 1 else "is",
                   n, "s" if n != 1 else ""))

    def structure_assumed(item_id, rec):
        # Nothing was assumed: every block is pinned by its printed reagent
        # amount and every structure matches the drawn bond-line count.
        return "0"

    def positional_ambiguity(item_id, rec):
        # A18 is symmetric and both nitrogens are primary, so there is no
        # choice of which N bears which arm. A2 has one primary N and one
        # tertiary N; the tertiary N cannot condense with an aldehyde.
        return "0"

    def value_is_bound(item_id, rec):
        return ""          # no value held at all -- not a bound either

    def normalisation(item_id, rec):
        return ""          # no value held; blank means "not stated"

    rows = emit.emit_rows(
        items, recs,
        confidence=confidence, note=note,
        positional_ambiguity=positional_ambiguity,
        structure_assumed=structure_assumed,
        value_is_bound=value_is_bound,
        value_normalization_unit=normalisation,
        within_model=None, id_key="lipid", label="HK_2026")

    print("emit.emit_rows accepted %d rows." % len(rows))
    print("rows_with_a_caveat (the F-024 definition, the number to quote): %d"
          % emit.caveat_count(rows))

    if "--seed-defect" in sys.argv:
        print("\n-- SEEDED: pass a LITERAL confidence, as F-011/F-025/F-046 "
              "each did --")
        try:
            emit.emit_rows(items, recs, confidence="high", note=note,
                           label="SEEDED literal")
        except EmitError as e:
            print("   REFUSED: %s" % str(e).split("\n")[0])
        else:
            raise SystemExit("emit_rows accepted a literal confidence. Error "
                             "class 29 is writable again.")
        print("\n-- SEEDED: a caveat with an empty note --")
        try:
            emit.emit_rows(items, recs,
                           confidence=lambda i, r: "low",
                           note=lambda i, r: "", label="SEEDED empty note")
        except EmitError as e:
            print("   REFUSED: %s" % str(e).split("\n")[0])
        else:
            raise SystemExit("emit_rows accepted a caveat with no note.")
        return 0

    with open(OUT, "w", newline="") as fh:
        cols = ["lipid"] + emit.columns()
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for it, r in zip(items, rows):
            r = dict(r)
            r["lipid"] = it
            w.writerow(r)
    print("wrote %s" % OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
