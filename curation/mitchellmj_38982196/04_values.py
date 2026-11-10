#!/usr/bin/env python3
"""Stage 04 -- the per-lipid values, straight out of the publisher Source Data.

TIER 1. No pixel is read for any value in this curation (00_value_sources.json).
The two sheets of MOESM4 hold the SAME 200 numbers in two layouts:

  'Fig. 2b'  columns 1..100, an EXPERIMENT number, with two replicates
  'Fig. 2c'  a 25-row x 4-tail block, the row being the AMINE number

so 'Fig. 2c' is what makes an experiment number resolvable to a molecule. This
script does not assume the order; it DERIVES the map by matching replicate
pairs, then checks the derived map three ways:

  * every experiment number resolves to exactly one (amine, tail) cell
  * Supplementary Figs. S3 and S4, which are labelled by AMINE number and were
    not used to build the map, must reproduce rows 11-14
  * the derived map must be a bijection

Run: python3 -I 04_values.py
"""
import csv
import json
import os
import sys

import openpyxl

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import emit                                                # noqa: E402
from errors import CurationLibError as CurationError        # noqa: E402

MOESM4 = os.path.join(HERE, "source", "41557_2024_1557_MOESM4_ESM.xlsx")
MOESM3 = os.path.join(HERE, "source", "41557_2024_1557_MOESM3_ESM.xlsx")
ENUM = os.path.join(HERE, "03_enumerated.csv")
MAP_JSON = os.path.join(HERE, "02_index_map.json")
VALUES = os.path.join(HERE, "05_values.csv")

TAILS = ("O9", "O12", "O14", "O16")


def sheet(path, name):
    wb = openpyxl.load_workbook(path, data_only=True)
    if name not in wb.sheetnames:
        raise CurationError(f"04_values: {os.path.basename(path)} has no "
                            f"sheet {name!r}; it has {wb.sheetnames}.")
    return wb[name]


def read_fig2b():
    ws = sheet(MOESM4, "Fig. 2b")
    hdr = {}
    for c in ws[6]:
        if c.value is not None and c.column_letter != "A":
            hdr[c.column_letter] = int(c.value)
    out = {}
    for col, no in hdr.items():
        r1 = ws[f"{col}7"].value
        r2 = ws[f"{col}8"].value
        if r1 is None or r2 is None:
            raise CurationError(f"04_values: Fig. 2b column {col} is blank.")
        out[no] = (float(r1), float(r2), f"{col}6:{col}8")
    return out


def read_fig2c():
    ws = sheet(MOESM4, "Fig. 2c")
    cols = {}
    for c in ws[7]:
        if c.value in TAILS:
            cols[c.value] = c.column_letter
    if sorted(cols) != sorted(TAILS):
        raise CurationError(f"04_values: Fig. 2c tail headers are {cols}.")
    # Each tail header sits over a PAIR of columns (two replicates).
    letters = sorted(cols.values(), key=lambda L: openpyxl.utils
                     .column_index_from_string(L))
    out = {}
    for row in range(8, 33):
        amine = ws[f"C{row}"].value
        if amine is None:
            continue
        amine = int(amine)
        for tail in TAILS:
            L = cols[tail]
            i = openpyxl.utils.column_index_from_string(L)
            a = ws.cell(row=row, column=i).value
            b = ws.cell(row=row, column=i + 1).value
            if a is None or b is None:
                raise CurationError(
                    f"04_values: Fig. 2c cell for amine {amine}, {tail} is "
                    "blank.")
            out[(amine, tail)] = (
                float(a), float(b),
                f"{L}{row}:{openpyxl.utils.get_column_letter(i + 1)}{row}")
    if len(out) != 100:
        raise CurationError(f"04_values: Fig. 2c gave {len(out)} cells, not "
                            "100.")
    return out, letters


def derive_map(b, c):
    """experiment no. -> (amine, tail), by matching replicate PAIRS.

    Matching the pair rather than the mean is deliberate: a mean collides far
    more easily, and a pair that matched in the wrong order would be a
    different experiment. Any collision raises.
    """
    index = {}
    for key, (r1, r2, _cells) in c.items():
        index.setdefault((r1, r2), []).append(key)
    fwd, unmatched = {}, []
    for no, (r1, r2, _cells) in sorted(b.items()):
        hit = index.get((r1, r2))
        if not hit:
            unmatched.append(no)
            continue
        if len(hit) > 1:
            raise CurationError(
                f"04_values: experiment {no} ({r1}, {r2}) matches "
                f"{len(hit)} different Fig. 2c cells {hit}. The map would be "
                "ambiguous; refusing to pick one.")
        fwd[no] = hit[0]
    if unmatched:
        raise CurationError(
            f"04_values: {len(unmatched)} Fig. 2b experiment(s) have no "
            f"Fig. 2c cell with the same replicate pair: {unmatched[:10]}.")
    if len(set(fwd.values())) != len(fwd):
        raise CurationError("04_values: the derived map is not a bijection.")
    return fwd


def check_si(c):
    """Supplementary Figs. S3 and S4 are labelled by AMINE. Independent check."""
    res = []
    ws = sheet(MOESM3, "Fig. S3")
    hdr = {}
    for cell in ws[6]:
        if cell.value is not None:
            hdr[cell.column_letter] = int(cell.value)
    for row in range(7, 11):
        tail = ws[f"B{row}"].value
        for L, amine in hdr.items():
            i = openpyxl.utils.column_index_from_string(L)
            a = float(ws.cell(row=row, column=i).value)
            b = float(ws.cell(row=row, column=i + 1).value)
            want = c[(amine, tail)]
            res.append({"sheet": "Fig. S3", "amine": amine, "tail": tail,
                        "si": [a, b], "fig2c": [want[0], want[1]],
                        "match": (a, b) == (want[0], want[1])})
    ws = sheet(MOESM3, "Fig. S4")
    hdr = {}
    for cell in ws[5]:
        if cell.value is not None:
            hdr[cell.column_letter] = int(cell.value)
    for row in range(6, 10):
        tail = ws[f"B{row}"].value
        for L, amine in hdr.items():
            i = openpyxl.utils.column_index_from_string(L)
            a = float(ws.cell(row=row, column=i).value)
            b = float(ws.cell(row=row, column=i + 1).value)
            want = c[(amine, tail)]
            res.append({"sheet": "Fig. S4", "amine": amine, "tail": tail,
                        "si": [a, b], "fig2c": [want[0], want[1]],
                        "match": (a, b) == (want[0], want[1])})
    bad = [r for r in res if not r["match"]]
    if bad:
        raise CurationError(
            f"04_values: {len(bad)} of {len(res)} Supplementary Fig. S3/S4 "
            f"cells disagree with Fig. 2c. First: {bad[0]}")
    return res


def main():
    b = read_fig2b()
    c, _letters = read_fig2c()
    if len(b) != 100:
        raise CurationError(f"04_values: Fig. 2b has {len(b)} experiments.")
    fwd = derive_map(b, c)
    si = check_si(c)

    with open(MAP_JSON, "w") as fh:
        json.dump({
            "_what": "Fig. 2b's 'experiment no.' 1..100 resolved to (amine, "
                     "tail). DERIVED by matching replicate pairs between the "
                     "two Source Data sheets, not assumed from a pattern.",
            "n": len(fwd),
            "map": {str(k): list(v) for k, v in sorted(fwd.items())},
            "shape_that_emerged":
                "Experiments 1-60 run amines 1-15 in blocks of 15, one block "
                "per tail in the order O9, O12, O14, O16; experiments 61-100 "
                "run amines 16-25 in blocks of 10, same tail order. That "
                "two-block layout is an OUTPUT of the matching, not an input "
                "to it, which is why it is recorded here rather than coded.",
            "independent_confirmation":
                f"{len(si)} cells of Supplementary Figs. S3 and S4, which are "
                "labelled by AMINE number and play no part in building the "
                "map, reproduce Fig. 2c exactly.",
            "si_check": si}, fh, indent=1)
    print(f"wrote {MAP_JSON}: {len(fwd)} experiments mapped; "
          f"{len(si)} SI cells agree")

    # ---------------- 05_values.csv ---------------------------------------
    with open(ENUM) as fh:
        enum = {r["lipid_name"]: r for r in csv.DictReader(fh)}

    recs, order = {}, []
    for no in sorted(fwd):
        amine, tail = fwd[no]
        name = f"{amine}T-{tail}"
        if name not in enum:
            raise CurationError(f"04_values: {name} is not in 03_enumerated.")
        r1, r2, cells_b = b[no]
        _a, _b2, cells_c = c[(amine, tail)]
        recs[name] = {
            "lipid_name": name, "amine": amine, "tail": tail,
            "experiment_no": no,
            "replicate_1": r1, "replicate_2": r2,
            "mean": round((r1 + r2) / 2.0, 3),
            "n_replicates": 2,
            "Model_target": "luminescence_normalized",
            "Model_type": "HepG2",
            "Model": "in_vitro",
            "readout": "firefly luciferase luminescence (a.u.), mLuc",
            "dose": "15 ng mRNA per 5,000 cells",
            "timepoint_h": 24,
            "source_file": "41557_2024_1557_MOESM4_ESM.xlsx",
            "sheet_fig2b": "Fig. 2b", "cells_fig2b": cells_b,
            "sheet_fig2c": "Fig. 2c", "cells_fig2c": cells_c,
            "value_tier": "1_publisher_source_data",
            "IL_SMILES": enum[name]["IL_SMILES"],
            "n_arms": int(enum[name]["n_arms"]),
            "confidence_from_structure": enum[name]["confidence"],
            "note_from_structure": enum[name]["curation_note"],
            "structure_assumed_from_structure": enum[name]["structure_assumed"],
            "positional_ambiguity_from_structure":
                enum[name]["positional_ambiguity"]}
        order.append(name)

    if len(order) != 100:
        raise CurationError(f"04_values: {len(order)} rows, expected 100.")

    q = emit.emit_rows(
        order, recs,
        confidence=lambda i, r: r["confidence_from_structure"],
        note=lambda i, r: r["note_from_structure"],
        structure_assumed=lambda i, r: r["structure_assumed_from_structure"],
        positional_ambiguity=lambda i, r: r[
            "positional_ambiguity_from_structure"],
        # Nothing here is a bound: both replicates are printed numbers in the
        # publisher's own workbook, including the two the paper calls
        # "background". A background-level READING is a measurement.
        value_is_bound=lambda i, r: "",
        # Left blank on purpose. value_normalization_unit describes the unit
        # over which Experiment_value was z-scored, and no z-scoring has
        # happened here: this file carries raw luminescence. The corpus
        # convention ("experiment") is applied by build_dataset.py at merge.
        value_normalization_unit=lambda i, r: "",
        label="MW_2024:values")

    cols = [k for k in recs[order[0]] if not k.endswith("_from_structure")]
    with open(VALUES, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols + emit.columns())
        w.writeheader()
        for name, extra in zip(order, q):
            w.writerow({**{k: recs[name][k] for k in cols}, **extra})
    print(f"wrote {VALUES}: {len(order)} rows, "
          f"{emit.caveat_count(q)} carrying a caveat")


if __name__ == "__main__":
    main()
