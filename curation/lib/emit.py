#!/usr/bin/env python3
"""THE single writer of the 79-column schema's quality columns.

WHY THIS FILE EXISTS. Two defects, and they are opposite halves of one problem.

  F-001  Four of five curations COMPUTED their flags, wrote them into their own
         folder's QC JSON, and emitted rows without them. The flags died one
         stage BEFORE the merge, so no change to build_dataset.py could have
         rescued them and nobody reading the builder would have found it.
         1,424 corpus rows carry a convention. "A flag that lives only in a
         per-folder QC file is not a flag, it is a memo."

  F-011  A hard-coded `"high"` in build_new_curation.py relabelled the 24
         lipids whose enumeration record reads "ASSUMED: no ring N and the
         paper draws no structure" as the project's highest confidence. F-025
         then found the same class surviving the round meant to clear it, and
         F-046 found it a third time. The labelling was INVERTED: the most
         uncertain rows carried the strongest claim.

So `emit_rows` does two things no per-curation emitter reliably did:

  1. It REFUSES A LITERAL CONFIDENCE. `confidence` must arrive as a callable or
     a per-item mapping keyed on the item's own id. Passing the string "high"
     raises. This makes error class 29 unwritable rather than discouraged.
  2. It refuses a caveat with no note. `confidence_without_a_note` was added as
     an audit metric in F-024 and found 6 rows immediately; here the condition
     cannot be emitted in the first place.

The column names and their allowed values come from
provenance/schema_quality.json and are NOT restated here. Restating them is
why they drifted and died: that file's own `_what` says so.
"""
import json
import os

from errors import EmitError

HERE = os.path.dirname(os.path.abspath(__file__))
SCHEMA = os.path.join(HERE, "..", "..", "provenance", "schema_quality.json")

# A caveat is a flag that QUALIFIES the row. `confidence = high` is an
# assurance, not a caveat -- counting it as one is what made F-024's 1,652
# differ from the honest 1,424 by exactly the 495 rows that carry only that.
CAVEAT_CONFIDENCE = {"medium", "low"}


def schema():
    """The quality columns, read from the single source of truth."""
    with open(SCHEMA) as fh:
        return json.load(fh)


def columns():
    return [c["name"] for c in schema()["columns"]]


def allowed(colname):
    for c in schema()["columns"]:
        if c["name"] == colname:
            return c.get("allowed")
    raise EmitError(f"emit: {colname!r} is not a quality column.")


def _resolve(spec, item_id, record, field):
    """Per-item value from a callable or a mapping. A literal raises."""
    if callable(spec):
        return spec(item_id, record)
    if isinstance(spec, dict):
        if item_id not in spec:
            raise EmitError(
                f"emit: no {field} for item {item_id!r}. A mapping must cover "
                "every item; a missing key is not a blank, it is an untested "
                "assumption.")
        return spec[item_id]
    if spec is None:
        return None
    raise EmitError(
        f"emit: {field}={spec!r} is a LITERAL.\n"
        f"Every row would receive {spec!r} regardless of its own record. That "
        "is error class 29, and it has now been found three times (F-011, "
        "F-025, F-046). In F-011 it relabelled 24 lipids whose record reads "
        "'ASSUMED: no ring N and the paper draws no structure' as the "
        "project's highest confidence -- the labelling was inverted, so the "
        "most uncertain rows carried the strongest claim.\n"
        f"Pass a callable (item_id, record) -> value, or a dict keyed on "
        "item id, so the value is traceable to the item's own enumeration "
        "record.")


def emit_rows(items, records, confidence, note, positional_ambiguity=None,
              structure_assumed=None, value_is_bound=None,
              value_normalization_unit=None, within_model=None,
              id_key="id", label=""):
    """Build the quality-column dict for each item. Validates, never defaults.

    items      iterable of item ids, in emission order
    records    {item_id: the enumeration record that item came from}
    confidence callable or dict -> "high"|"medium"|"low"|""   (NOT a literal)
    note       callable or dict -> one sentence                (NOT a literal)
    the rest   callable, dict, or None (meaning "this curation did not assess
               it", which the schema defines as an absence of information)

    Returns [{column: value}], one dict per item, ready to splice into the row.
    """
    sch = {c["name"]: c for c in schema()["columns"]}
    out = []
    for it in items:
        rec = records.get(it)
        if rec is None:
            raise EmitError(
                f"emit_rows({label!r}): item {it!r} has no enumeration record. "
                "A flag with no record behind it cannot be traced, which is "
                "the whole defect F-011 was.")
        row = {}
        conf = _resolve(confidence, it, rec, "confidence")
        nt = _resolve(note, it, rec, "curation_note")
        for col, spec in (("positional_ambiguity", positional_ambiguity),
                          ("structure_assumed", structure_assumed),
                          ("value_is_bound", value_is_bound),
                          ("value_normalization_unit",
                           value_normalization_unit),
                          ("Experiment_value_within_model", within_model)):
            row[col] = ("" if spec is None
                        else _resolve(spec, it, rec, col))
        row["confidence"] = "" if conf is None else str(conf)
        row["curation_note"] = "" if nt is None else str(nt)

        # --- vocabulary, from the schema, not from a list in this file ---
        for col, val in row.items():
            allow = sch[col].get("allowed")
            if allow is not None and str(val) not in allow:
                raise EmitError(
                    f"emit_rows({label!r}) item {it!r}: {col}={val!r} is not "
                    f"in the schema's allowed set {allow}. The vocabulary "
                    "lives in provenance/schema_quality.json; audit.py C13 "
                    "enforces the same set, so an out-of-vocabulary value "
                    "would fail the build after the merge instead of here.")

        # --- a caveat must say what it is ---
        caveat = (row["confidence"] in CAVEAT_CONFIDENCE
                  or str(row["positional_ambiguity"]) == "1"
                  or str(row["structure_assumed"]) == "1"
                  or row["value_is_bound"] in ("upper", "lower"))
        if caveat and not row["curation_note"].strip():
            raise EmitError(
                f"emit_rows({label!r}) item {it!r}: this row carries a caveat "
                f"(confidence={row['confidence']!r}, "
                f"positional_ambiguity={row['positional_ambiguity']!r}, "
                f"structure_assumed={row['structure_assumed']!r}, "
                f"value_is_bound={row['value_is_bound']!r}) and its "
                "curation_note is empty.\n"
                "A flag that cannot say what it means gets ignored. "
                "audit.py C13's confidence_without_a_note was added for "
                "exactly this and found 6 rows on its first run; writing the "
                "note here is cheaper than rediscovering the reason later -- "
                "in that case it took re-reading 03_enumerated.csv to learn "
                "the arm was a branched isomer mixture while the QC file "
                "blamed the wrong component.")
        out.append(row)
    return out


def caveat_count(rows):
    """rows_with_a_caveat, by the F-024 definition. The number to quote.

    Deliberately NOT rows_carrying_a_flag: that counts `confidence = high`,
    which is an assurance. The two differ by the 495 rows carrying only that,
    and conflating them is how the headline count moved four times.
    """
    n = 0
    for r in rows:
        if (r.get("confidence") in CAVEAT_CONFIDENCE
                or str(r.get("positional_ambiguity") or "") == "1"
                or str(r.get("structure_assumed") or "") == "1"
                or (r.get("value_is_bound") or "") in ("upper", "lower")):
            n += 1
    return n
