#!/usr/bin/env python3
"""Join structures to values and write the emittable rows, flags ON THE ROW.

    python3 -I stage05_emit.py      ->  06_emitted_rows.csv

Every quality column is written by curation/lib/emit.emit_rows, which refuses a
literal `confidence` and refuses a caveat with an empty note. The confidence of
each row is a function of THAT row's own record -- whether the SI draws the
structure, whether its HRMS m/z was reproduced, and whether re-assembling the
blocks gave back the drawn molecule -- never a constant.

NOT MERGED. This writes the curation's own file; nothing here touches data/.
The experiment key, the z-score convention and the grouping of the four in vivo
panels are left open in DECISIONS.md, because each of them would change a
published number and none of them is this script's to settle.
"""
import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "lib"))

import emit                                                # noqa: E402

PANEL_KEY = {
    "Extended Figure 1": "DG_2025",
    "Figure 1f": "DG_2025_vivo_screen3",
    "2a": "DG_2025_vivo_gen2_tails",
    "2b": "DG_2025_vivo_gen3_heads",
    "2c": "DG_2025_vivo_benchmarks",
    "Extended Figure 2b": "DG_2025_vivo_gen4_tails",
}
# Benchmarks on the comparison panel are not members of this library.
NOT_LIBRARY = {"SM-102", "ckk-e12", "cKK-E12", "LP01", "Lipid5", "PBS"}


def main():
    enum = {r["label"]: r for r in csv.DictReader(
        open(os.path.join(HERE, "03_enumerated.csv")))}
    struct = {}
    for r in csv.DictReader(open(os.path.join(HERE, "02_structures.csv"))):
        struct.setdefault(r["label"], r)
    values = list(csv.DictReader(open(os.path.join(HERE, "05_values.csv"))))

    items, records, rows, skipped = [], {}, [], []
    for v in values:
        lab = v["lipid_label"]
        if lab in NOT_LIBRARY:
            skipped.append((lab, v["panel"], "benchmark, not a library member"))
            continue
        e = enum.get(lab)
        if e is None or not e["SMILES"]:
            skipped.append((lab, v["panel"],
                            "no structure established" if e is None
                            else e["formula"]))
            continue
        s = struct.get(lab, {})
        item = f'{v["panel"]}|{lab}'
        items.append(item)
        records[item] = {
            "label": lab, "panel": v["panel"], "SMILES": e["SMILES"],
            "provenance": e["provenance"],
            "drawn_in_SI": e["drawn_in_SI"] == "True",
            "reassembly_matches_drawing": e["matches_drawing"] == "True",
            "hrms_status": s.get("status", "NO_STRUCTURE_ROW"),
            "hrms_ppm": s.get("ppm", ""),
            "mean": float(v["mean"]), "value": v,
        }

    def confidence(item, rec):
        if rec["drawn_in_SI"] and rec["hrms_status"] == "PASS" \
                and rec["reassembly_matches_drawing"]:
            return "high"
        return "medium"

    def note(item, rec):
        bits = []
        if rec["drawn_in_SI"]:
            bits.append("structure read from the vector coordinates of the "
                        "SI's own drawing")
            if rec["hrms_status"] == "PASS":
                bits.append(f'its printed HRMS m/z reproduced to '
                            f'{rec["hrms_ppm"]} ppm')
            else:
                bits.append("the SI prints no HRMS beside it, so the only "
                            "external mass check is absent")
            if rec["reassembly_matches_drawing"]:
                bits.append("and re-assembling the headgroup and tail blocks "
                            "reproduces it exactly")
        else:
            bits.append("the SI draws no structure for this member; it is "
                        "built from the headgroup and tail blocks, each of "
                        "which was derived independently from other lipids "
                        "in the same series")
        if float(rec["mean"]) == 0.0:
            bits.append("the published value is exactly 0 RFU, which is the "
                        "assay floor and not a measurement of zero")
        return "; ".join(bits) + "."

    def structure_assumed(item, rec):
        return "0" if rec["drawn_in_SI"] else "1"

    def value_is_bound(item, rec):
        return "upper" if float(rec["mean"]) == 0.0 else ""

    quality = emit.emit_rows(
        items, records, confidence=confidence, note=note,
        structure_assumed=structure_assumed, value_is_bound=value_is_bound,
        positional_ambiguity=lambda i, r: "0",
        value_normalization_unit=lambda i, r: "",
        label="andersondg_41203968")

    cols = ["Experiment_ID_proposed", "IL_name", "IL_SMILES", "panel", "model",
            "biological_system", "readout", "dose", "timepoint",
            "value_as_published", "sd_as_published", "n_values", "value_tier",
            "workbook", "sheet", "cells", "structure_provenance",
            "hrms_status", "hrms_ppm"] + emit.columns()
    out = []
    for item, q in zip(items, quality):
        rec = records[item]
        v = rec["value"]
        row = {
            "Experiment_ID_proposed": PANEL_KEY[v["panel"]],
            "IL_name": rec["label"], "IL_SMILES": rec["SMILES"],
            "panel": v["panel"], "model": v["model"],
            "biological_system": v["biological_system"],
            "readout": v["readout"], "dose": v["dose"],
            "timepoint": v["timepoint"], "value_as_published": v["mean"],
            "sd_as_published": v["sd_as_published"],
            "n_values": v["n_values"], "value_tier": v["value_tier"],
            "workbook": v["workbook"], "sheet": v["sheet"], "cells": v["cells"],
            "structure_provenance": rec["provenance"],
            "hrms_status": rec["hrms_status"], "hrms_ppm": rec["hrms_ppm"]}
        row.update(q)
        out.append(row)

    with open(os.path.join(HERE, "06_emitted_rows.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in out:
            w.writerow({k: r.get(k, "") for k in cols})

    from collections import Counter
    print("rows emitted:", len(out))
    for k, n in Counter(r["Experiment_ID_proposed"] for r in out).items():
        print(f"   {k:<30} {n}")
    print("rows with a caveat (emit.caveat_count):", emit.caveat_count(quality))
    print("confidence:", dict(Counter(r["confidence"] for r in out)))
    print("not emitted:", len(skipped))
    for s in sorted(set(skipped)):
        print("   ", s)
    json.dump({"skipped": sorted(set(skipped))},
              open(os.path.join(HERE, "06_not_emitted.json"), "w"), indent=1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
