#!/usr/bin/env python3
"""Emit rows in the 79-column schema -> 07_emit_preview.csv  (NOT a merge file)

Run from this folder:   python3 -I stage07_emit.py

THIS FOLDER EMITS NOTHING INTO data/. DECISIONS.md D-1 is OPEN -- whether an
antibacterial MIC against *Bacillus subtilis* belongs in a delivery-only corpus
at all is a corpus-scope question that this curation is not entitled to answer,
and `build_dataset.py` refuses a folder with an OPEN entry. So what this stage
produces is a PREVIEW: the quality columns, computed per member, written to
07_emit_preview.csv, with the row's identity and value beside them.

WHY IT RUNS AT ALL. F-001: four of five curations computed their flags, wrote
them into their own folder's QC JSON, and emitted rows without them. 1,424
corpus rows carry a convention and say so only because the repair went back and
rebuilt them. "A flag that lives only in a per-folder QC file is not a flag, it
is a memo." The flags for this paper are therefore produced HERE, by
curation/lib/emit.emit_rows, in the shape a merge would consume -- so that if
D-1 ever resolves to merge, the flags are already attached to the rows and not
to a summary of them.

AND THE CONFIDENCE IS DERIVED, NEVER WRITTEN. emit_rows refuses a literal
(error class 29: F-011, F-025, F-046, three findings for one mistake). Both
`confidence` and `curation_note` below are callables that read the member's own
record out of 03_enumerated.csv.
"""
import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))

import emit                                            # noqa: E402

OUT = os.path.join(HERE, "07_emit_preview.csv")
SIDE = os.path.join(HERE, "07_emit_preview_summary.json")

SCORE_MEANING_BOUND = {
    0: ("upper",
        "Antimicrobial score 0 means 'no inhibition at any studied "
        "concentration' (SI Table S2), i.e. the lowest concentration at which "
        "this lipidoid acts is above the highest tested 100 uM. The stored "
        "value is a ceiling on the activity, not a measurement of it."),
    6: ("lower",
        "Antimicrobial score 6 means 'delayed or no growth at 25 uM' (SI Table "
        "S2), which is the LOWEST concentration the 104-member screen tested, "
        "so the true potency is at or below 25 uM and the score is a floor. "
        "For the 27 members re-tested at 1-12.5 uM the MIC in main-text Table "
        "1 resolves it; 05_values.csv carries that MIC where it exists."),
}


def main():
    enum = {r["IL_name"]: r for r in
            csv.DictReader(open(os.path.join(HERE, "03_enumerated.csv")))}
    vals = list(csv.DictReader(open(os.path.join(HERE, "05_values.csv"))))

    # one record per member: the enumeration record plus its value row, so
    # every flag below is traceable to the member's own record.
    records = {}
    for v in vals:
        r = dict(enum[v["IL_name"]])
        r["antimicrobial_score"] = int(v["antimicrobial_score"])
        r["MIC_uM"] = v["MIC_uM"]
        records[v["IL_name"]] = r
    items = [v["IL_name"] for v in vals]

    def confidence(item_id, rec):
        """medium where the structure rests on a convention, else high.

        The only thing that lowers confidence here is the tail being a
        commercial isomer mixture (DECISIONS.md D-5). Nothing about the VALUE
        lowers it: every value in this curation is a printed numeral read out
        of a text table, not a pixel.
        """
        return "medium" if rec["structure_caveat_kind"] else "high"

    def note(item_id, rec):
        bits = []
        if rec["structure_caveat_kind"]:
            bits.append(rec["structure_caveat_why"])
        b = SCORE_MEANING_BOUND.get(rec["antimicrobial_score"])
        if b:
            bits.append(b[1])
        return " ".join(bits)

    def structure_assumed(item_id, rec):
        return "1" if rec["structure_caveat_kind"] == "structure_assumed" \
            else "0"

    def value_is_bound(item_id, rec):
        b = SCORE_MEANING_BOUND.get(rec["antimicrobial_score"])
        return b[0] if b else ""

    def positional_ambiguity(item_id, rec):
        """0 everywhere, and that is a measured 0, not a blank.

        Every nitrogen in every member is exhaustively substituted and then
        exhaustively methylated, so there is no 'which nitrogen bears which
        arm' choice to make anywhere in this library -- unlike LL_2026, where
        the regiochemistry was the bug. Recorded as an explicit 0 so the
        absence of ambiguity is a statement rather than a gap.
        """
        return "0"

    flags = emit.emit_rows(
        items, records,
        confidence=confidence, note=note,
        structure_assumed=structure_assumed,
        positional_ambiguity=positional_ambiguity,
        value_is_bound=value_is_bound,
        value_normalization_unit=None,   # no Experiment_value is produced:
        within_model=None,               # DECISIONS.md D-3 is OPEN.
        label="pabst_37602460")

    cols = (["IL_name", "headgroup", "tail", "antimicrobial_score", "MIC_uM",
             "IL_SMILES", "parent_SMILES"] + emit.columns())
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for name, fl in zip(items, flags):
            r = records[name]
            row = {"IL_name": name, "headgroup": r["headgroup"],
                   "tail": r["tail"],
                   "antimicrobial_score": r["antimicrobial_score"],
                   "MIC_uM": r["MIC_uM"], "IL_SMILES": r["IL_SMILES"],
                   "parent_SMILES": r["parent_SMILES"]}
            row.update(fl)
            w.writerow(row)

    summary = {
        "rows_previewed": len(flags),
        "rows_merged_into_data": 0,
        "why_zero_merged": "DECISIONS.md D-1 is OPEN: the endpoint is "
                           "antibacterial growth inhibition of Bacillus "
                           "subtilis, not nucleic-acid delivery, and "
                           "build_dataset.py refuses a folder with an OPEN "
                           "entry.",
        "rows_with_a_caveat": emit.caveat_count(flags),
        "rows_carrying_any_flag": sum(
            1 for f in flags if any(str(v).strip() for v in f.values())),
        "by_flag": {
            "confidence_medium": sum(1 for f in flags
                                     if f["confidence"] == "medium"),
            "confidence_high": sum(1 for f in flags
                                   if f["confidence"] == "high"),
            "structure_assumed_1": sum(1 for f in flags
                                       if f["structure_assumed"] == "1"),
            "value_is_bound_upper": sum(1 for f in flags
                                        if f["value_is_bound"] == "upper"),
            "value_is_bound_lower": sum(1 for f in flags
                                        if f["value_is_bound"] == "lower"),
            "positional_ambiguity_1": sum(1 for f in flags
                                          if f["positional_ambiguity"] == "1"),
        },
        "bound_rows_resolved_by_a_measured_MIC": sum(
            1 for name, f in zip(items, flags)
            if f["value_is_bound"] and records[name]["MIC_uM"]),
    }
    with open(SIDE, "w") as fh:
        json.dump(summary, fh, indent=1)

    print(f"wrote {len(flags)} PREVIEW rows -> {OUT}")
    print(f"  merged into data/        0  (DECISIONS.md D-1 is OPEN)")
    print(f"  rows_with_a_caveat     {summary['rows_with_a_caveat']}")
    for k, v in summary["by_flag"].items():
        print(f"  {k:28s} {v}")
    print(f"  of the bounded rows, resolved by a measured MIC: "
          f"{summary['bound_rows_resolved_by_a_measured_MIC']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
