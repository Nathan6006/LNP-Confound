#!/usr/bin/env python3
"""Stage 05 -- the values. Table S2, SI p.17, transcribed and checked.

READOUT, RECORDED AS PRINTED AND NOT CONVERTED (CLAUDE.md section 3 /
00_value_sources.json):

    column header            "%TE"
    caption                  "transfection efficiency >= 40% in HEK293T cells"
    assay                    eGFP plasmid DNA (100% eGFP) lipoplex, HEK293T,
                             24 h, PI counterstain, flow cytometry
                             (Guava easyCyte HT)
    unit                     percent of cells GFP-positive
    NOT                      luminescence; NOT knockdown; NOT a fold change

No log, no z-score, no renormalisation is applied here. The corpus convention
is a per-experiment z of log10, and applying it to THIS table would be wrong
for a reason that is on the face of the source: only the 45 lipidoids at or
above 40% are reported, so any moment computed here is a moment of a truncated
top tail, not of the experiment. DECISIONS.md D-4.

TRANSCRIPTION PROVENANCE. The table is a raster image: 600 dpi render, borders
located by ink projection (rules at y=3297/3403/5756 and x=814/1343/2017/2522
left, 2635/3163/3766/4306 right), rows cropped at those detected boundaries.
Read TWICE with different crop geometry -- once as a half-scale whole-table
image, once as eight native-resolution chunks -- and the two reads agree on
45/45 names and 45/45 values.

Run:  python3 -I stage05_values.py
"""
import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

# (given name, structure name AS PRINTED, %TE)
TABLE_S2 = [
    ("L1",  "T18U-PY16-A1",  95), ("L2",  "T18U-PY12-A8", 89),
    ("L3",  "T14-PY12-A8",   87), ("L4",  "T14-PY16-A14", 79),
    ("L5",  "T14-PY16-A13",  75), ("L6",  "T18u-PY18-A4", 75),
    ("L7",  "T14-PY14-A13",  75), ("L8",  "T18-PY12-A3",  72),
    ("L9",  "T14-PY16-A17",  70), ("L10", "T14-PY14-A17", 70),
    ("L11", "T16-PY12-A1",   68), ("L12", "T16-PY12-A17", 65),
    ("L13", "T14-PY14-A4",   60), ("L14", "T18U-PY12-A17", 60),
    ("L15", "T16-PY12-A4",   58), ("L16", "T14-PY12-A13", 57),
    ("L17", "T16-PY18-A13",  57), ("L18", "T16-PY18-A13", 57),
    ("L19", "T18-PY16-A3",   57), ("L20", "T18-PY18-A3",  57),
    ("L21", "T16-PY18-A9",   55), ("L22", "T14-PY14-A3",  52),
    ("L23", "T16-PY12-A13",  52), ("L24", "T18-PY12-A6",  52),
    ("L25", "T18U-PY14-A7",  52), ("L26", "T16-PY18-A13", 50),
    ("L27", "T16-PY12-A8",   50), ("L28", "T16-PY12-A14", 50),
    ("L29", "T16-PY14-A3",   50), ("L30", "T18-PY16-A11", 50),
    ("L31", "T18-PY18-A5",   50), ("L32", "T18U-PY18-A15", 50),
    ("L33", "T14-PY18-A9",   48), ("L34", "T18-PY18-A4",  48),
    ("L35", "T18-PY18-A1",   47), ("L36", "T16-PY14-A9",  46),
    ("L37", "T18-PY16-A4",   45), ("L38", "T16-PY18-A14", 44),
    ("L39", "T14-PY12-A17",  44), ("L40", "T18U-PY18-A13", 42),
    ("L41", "T18U-PY16-A6",  40), ("L42", "T18U-PY16-A8", 40),
    ("L43", "T18U-PY14-A13", 40), ("L44", "T18U-PY16-A5", 40),
    ("L45", "T16-PY18-A3",   40),
]
# SI p.2 writes the oleoyl thiolactone as T18U everywhere; Table S2 row L6
# writes "T18u". Same block, a case slip in the table, recorded not silently
# normalised -- 08_contradictions.csv C-1.
CASE_SLIPS = {"T18u": "T18U"}


def main():
    enum = {r["lipidoid"]: r for r in csv.DictReader(
        open(os.path.join(HERE, "03_enumerated.csv")))}

    if len(TABLE_S2) != 45:
        raise SystemExit(f"stage05: {len(TABLE_S2)} rows, Figure S3's caption "
                         "says 45.")
    if [t[0] for t in TABLE_S2] != [f"L{i}" for i in range(1, 46)]:
        raise SystemExit("stage05: the L-numbering is not 1..45 in order.")
    vals = [t[2] for t in TABLE_S2]
    if any(vals[i] < vals[i + 1] for i in range(44)):
        raise SystemExit(
            "stage05: %TE is not non-increasing down the table, but Figure "
            "S3's caption says the list starts from the top performer. A "
            "break here is a transcription error, not a finding.")
    if min(vals) < 40:
        raise SystemExit(
            f"stage05: a value of {min(vals)} is below the caption's stated "
            "inclusion threshold of 40%.")

    rows = []
    for name, printed, te in TABLE_S2:
        canon = printed
        for bad, good in CASE_SLIPS.items():
            canon = canon.replace(bad, good)
        t, p, a = canon.split("-")
        e = enum.get(canon)
        if e is None:
            raise SystemExit(
                f"stage05: {canon!r} is not one of the 288 enumerated labels. "
                "A value whose lipidoid is not in the library grid is a "
                "transcription error on one side or the other.")
        rows.append({
            "lipidoid": name,
            "structure_name_as_printed": printed,
            "structure_name_canonical": canon,
            "thiolactone": t, "disulfide": p, "amine": a,
            "readout": "%TE (percent transfection efficiency)",
            "unit": "percent of GFP-positive cells",
            "value": te,
            "value_transformed": "",
            "cell_line": "HEK293T",
            "payload": "eGFP plasmid DNA, 100% eGFP",
            "timepoint": "24 h",
            "instrument": "Guava easyCyte HT flow cytometer",
            "n_replicates": "not stated",
            "source_locus": "SI p.17, Table S2, column '%TE'; the same 45 "
                            "lipidoids are drawn in Figure S3 on the same page",
            "source_tier": "5_pixels (printed numerals in a raster table)",
            "structure_determined": e["structure_determined"],
            "IL_SMILES": e["IL_SMILES"],
        })

    with open(os.path.join(HERE, "05_values.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    dup_name = {}
    for r in rows:
        dup_name.setdefault(r["structure_name_canonical"], []).append(
            (r["lipidoid"], r["value"]))
    repeats = {k: v for k, v in dup_name.items() if len(v) > 1}
    det = [r for r in rows if r["structure_determined"] == "True"]
    summary = {
        "rows": len(rows),
        "distinct_structure_names": len(dup_name),
        "repeated_structure_names": repeats,
        "rows_with_a_determined_structure": len(det),
        "determined_lipidoids": [(r["lipidoid"],
                                  r["structure_name_canonical"], r["value"])
                                 for r in det],
        "value_range": [min(vals), max(vals)],
        "amines_used": sorted({r["amine"] for r in rows},
                              key=lambda s: int(s[1:])),
        "thiolactones_used": sorted({r["thiolactone"] for r in rows}),
        "disulfides_used": sorted({r["disulfide"] for r in rows}),
    }
    with open(os.path.join(HERE, "05_values_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1)
    print(json.dumps(summary, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
