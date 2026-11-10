#!/usr/bin/env python3
"""Stage 06 -- emit rows, with the quality flags ON THE ROW.

F-001: four of five curations computed their flags, wrote them into their own
folder's QC JSON and emitted rows without them. "A flag that lives only in a
per-folder QC file is not a flag, it is a memo." So every flag here goes
through curation/lib/emit.emit_rows, which refuses a literal confidence
(error class 29, found three times: F-011, F-025, F-046) and refuses a caveat
with no note (F-024).

WHAT IS EMITTABLE. 45 lipidoids carry a value. Five of them have a structure:
the 13 amine head groups used in Table S2 include only A1 and A5, and the other
eleven are defined nowhere in this SI (00_acquisition.json, DECISIONS.md D-1).
A row with no IL_SMILES is not a corpus row -- C1 would reject it -- so the
other 40 are written here with `emittable = False` and their reason, so the
gap is visible rather than absent.

Experiment_value is deliberately LEFT BLANK. The corpus convention is a
per-experiment z of log10, and Table S2 reports only the 45 lipidoids at or
above 40% of 288 measured. Any mean or sd computed here is a moment of a
truncated top tail, not of the experiment. DECISIONS.md D-4 queues that;
this stage does not pre-empt it.

Run:  python3 -I stage06_emit.py
"""
import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))

import emit                                                 # noqa: E402

EXPERIMENT_ID = "PL_2018"
RACEMIC = ("the alpha carbon is racemic -- the core is DL-homocysteine "
           "thiolactone (SI p.2), so the stored SMILES carries no stereo "
           "descriptor there by design")


def main():
    vals = list(csv.DictReader(open(os.path.join(HERE, "05_values.csv"))))
    amines = json.load(open(os.path.join(HERE, "01_amines.json")))

    items, records = [], {}
    for r in vals:
        items.append(r["lipidoid"])
        records[r["lipidoid"]] = r

    def confidence(i, rec):
        if rec["structure_determined"] != "True":
            return ""
        return amines[rec["amine"]]["confidence"]

    def note(i, rec):
        if rec["structure_determined"] != "True":
            return ""
        if rec["amine"] == "A5":
            return ("A5's three-carbon linker is measured from the Figure S1b "
                    "drawing but its two N-ethyl branches are read, not "
                    "measured (the raster is at its noise floor); the paper "
                    "gives A5 no name, formula or mass, so there is no anchor "
                    "orthogonal to the drawing; " + RACEMIC)
        return RACEMIC

    def structure_assumed(i, rec):
        if rec["structure_determined"] != "True":
            return ""
        return "1" if rec["amine"] == "A5" else "0"

    def positional(i, rec):
        if rec["structure_determined"] != "True":
            return ""
        # Not ambiguous: the thiolactone has one carbonyl and releases one
        # thiol, and Scheme S1c plus Figure S1 draw both regiochemistries
        # explicitly. The crude methanol-wash workup affects purity, not which
        # atom bears which substituent.
        return "0"

    quality = emit.emit_rows(
        items, records, confidence=confidence, note=note,
        positional_ambiguity=positional, structure_assumed=structure_assumed,
        value_is_bound=None,              # measurements, not bounds
        value_normalization_unit=None,    # not standardised here; see D-4
        within_model=None, label="PL_2018")

    rows = []
    for r, q in zip(vals, quality):
        det = r["structure_determined"] == "True"
        row = {
            "Experiment_ID": EXPERIMENT_ID,
            "IL_name": r["structure_name_canonical"],
            "given_name": r["lipidoid"],
            "IL_SMILES": r["IL_SMILES"],
            "Model_type": "HEK293T",
            "Model_target": "eGFP plasmid DNA",
            "Experiment_value": "",
            "value_raw_percent_TE": r["value"],
            "readout": r["readout"], "unit": r["unit"],
            "emittable": det,
            "not_emittable_because": "" if det else (
                f"head group {r['amine']} is not defined anywhere in this SI"
                if r["amine"] not in ("A8", "A10") else
                "head group A8 is known only as the molecular formula C7H17N3"
                if r["amine"] == "A8" else
                "head group A10's ring size measurement is inconclusive"),
        }
        row.update(q)
        rows.append(row)

    with open(os.path.join(HERE, "06_rows.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    emittable = [r for r in rows if r["emittable"]]
    summary = {
        "rows_written": len(rows),
        "rows_emittable_into_the_corpus": len(emittable),
        "rows_with_a_caveat (emit.caveat_count, F-024 definition)":
            emit.caveat_count(quality),
        "rows_carrying_any_flag":
            sum(1 for q in quality if any(str(v).strip() for v in q.values())),
        "emittable": [{k: r[k] for k in
                       ("given_name", "IL_name", "value_raw_percent_TE",
                        "confidence", "structure_assumed")}
                      for r in emittable],
        "under_100_usable_rows": len(emittable) < 100,
    }
    with open(os.path.join(HERE, "06_emit_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1)
    print(json.dumps(summary, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
