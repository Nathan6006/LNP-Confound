#!/usr/bin/env python3
"""Attach the 79-column quality block to every candidate row, via lib/emit.

Why this runs even though nothing is being merged: F-001 is that four of five
curations COMPUTED their flags, wrote them into a folder QC file, and emitted
rows without them -- the flags died one stage before the merge. Building the
block here, keyed on each row's OWN enumeration/measurement record, is what
makes it survivable. "A flag that lives only in a per-folder QC file is not a
flag, it is a memo."

`confidence` is derived by a RULE over the row's own evidence, never written
as a literal; emit.emit_rows refuses a literal (error class 29, found three
times: F-011, F-025, F-046).
"""
import csv, json, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "lib"))
import emit                                                  # noqa: E402
HERE = os.path.dirname(os.path.abspath(__file__))

# The structures whose SI synthesis paragraph names an aldehyde that three
# independent lines of the authors' own evidence contradict. See
# 08_contradictions.csv C-03/C-04 and 04_structure_qc.json P4.
ADJUDICATED = {"A4B4-S3", "A4B4-S4"}


def confidence_of(item_id, rec):
    """From the row's own record. High only with >=2 independent external
    anchors AND no unresolved conflict in the source."""
    if rec["library_status"].startswith("benchmark"):
        return "high"          # provenance/reagent_registry.py, from the name
    if rec["IL_name"] in ADJUDICATED:
        return "medium"
    return "high"              # drawn in SI Figs S1-S9 + Table S1 formula + MW


def note_of(item_id, rec):
    bits = []
    if rec["IL_name"] in ADJUDICATED:
        bits.append(
            "The SI synthesis paragraph names the aldehyde '6-oxohexyl "
            "dodecanoate' (B3); the product's own IUPAC name, its MS formula "
            "and its printed reagent mass/mol all give 6-oxohexyl "
            "2-hexyldecanoate (B4), which is what is stored")
    if rec["value_is_bound"] == "lower":
        bits.append(
            "the bar top is above the axis break and the upper axis segment "
            "carries one labelled tick, so the panel gives only a LOWER "
            "BOUND of %s and no point value" % rec["lower_bound"])
    if rec["value_is_bound"] == "upper":
        bits.append(
            "the bar is drawn with no resolvable height, so %s is an UPPER "
            "BOUND at the panel's detection floor, not a measurement of zero"
            % rec["upper_bound"])
    if rec["IL_name"] == "H2A4-Ole":
        bits.append(
            "H2A4-Ole has the same full InChIKey as library member H2A4B2; "
            "the paper presents it as newly synthesised")
    if rec["refused"] and not bits:
        bits.append("no value: " + rec["refused"].split(".")[0])
    if not bits:
        bits.append(
            "value digitised from %s (tier 5); structure from the drawn SI "
            "library figure, formula and MW confirmed against Table S1"
            % rec["panel"])
    return "; ".join(bits) + "."


def main():
    rows = list(csv.DictReader(open(os.path.join(HERE, "05_values.csv"))))
    ids = ["%s|%s" % (r["panel"], r["IL_name"]) for r in rows]
    recs = dict(zip(ids, rows))

    block = emit.emit_rows(
        ids, recs,
        confidence=confidence_of,
        note=note_of,
        structure_assumed=lambda i, r: "0",
        positional_ambiguity=lambda i, r: "0",
        value_is_bound=lambda i, r: r["value_is_bound"],
        value_normalization_unit=lambda i, r: "",   # nothing is z-scored here
        label="li_39883839")

    out = os.path.join(HERE, "05_values_flagged.csv")
    cols = list(rows[0].keys()) + emit.columns()
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r, q in zip(rows, block):
            rr = dict(r); rr.update(q); w.writerow(rr)

    n_cav = emit.caveat_count(block)
    print("05_values_flagged.csv: %d rows, %d columns" % (len(rows), len(cols)))
    print("rows_with_a_caveat (the F-024 metric, the one to quote):", n_cav)
    print("by confidence:", {c: sum(1 for b in block if b["confidence"] == c)
                             for c in ("high", "medium", "low", "")})
    print("value_is_bound=upper:",
          sum(1 for b in block if b["value_is_bound"] == "upper"))
    json.dump({"n_rows": len(rows), "rows_with_a_caveat": n_cav,
               "by_confidence": {c: sum(1 for b in block
                                        if b["confidence"] == c)
                                 for c in ("high", "medium", "low", "")},
               "value_is_bound_upper": sum(1 for b in block
                                           if b["value_is_bound"] == "upper")},
              open(os.path.join(HERE, "05_flag_summary.json"), "w"), indent=1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
