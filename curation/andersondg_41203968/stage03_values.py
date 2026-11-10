#!/usr/bin/env python3
"""Per-lipid values, read from the PUBLISHER SOURCE DATA workbooks (tier 1).

    python3 -I stage03_values.py

Reads  source/41565_2025_2044_MOESM{3,4,8,9}_ESM.xlsx
Writes 05_values.csv and 05_printed_statistics.json

No pixel is read anywhere in this curation. Every value below is a number the
publisher released as Source Data, and each row records the workbook, the sheet
and the exact cell range it came from, so a reviewer can open the file and look
at the same cells.

WHAT IS AND IS NOT A PER-LIPID VALUE HERE. Figure 1d and 1e are BATCH screens:
each bar is a pool of ten (1d) or three-to-four (1e) lipids injected together,
so no bar belongs to a lipid. They are not emitted. Figure 1f (screen 3) is six
lipids injected individually and is emitted. Extended Figure 1 is the whole
100-member library, one well per lipid, and is the only panel that covers the
library.
"""
import csv
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "source")
sys.path.insert(0, HERE)

from openpyxl import load_workbook                         # noqa: E402

F = {n: os.path.join(SRC, f"41565_2025_2044_MOESM{n}_ESM.xlsx")
     for n in (3, 4, 5, 6, 7, 8, 9, 10, 11, 12)}


def sheet(n, name):
    wb = load_workbook(F[n], read_only=True, data_only=True)
    ws = wb[name]
    rows = [list(r) for r in ws.iter_rows(values_only=True)]
    wb.close()
    return rows


def col_letter(i):
    s = ""
    i += 1
    while i:
        i, r = divmod(i - 1, 26)
        s = chr(65 + r) + s
    return s


def num(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool)


# ---------------------------------------------------------------- panels ---
OUT = []


def add(panel, wb_n, sh, cells, label, model, system, readout, dose, t, vals,
        note="", sd=None, n_stated=None):
    vals = [float(v) for v in vals if num(v)]
    OUT.append({
        "panel": panel, "workbook": os.path.basename(F[wb_n]), "sheet": sh,
        "cells": cells, "lipid_label": label, "model": model,
        "biological_system": system, "readout": readout, "dose": dose,
        "timepoint": t, "n_values": len(vals),
        "n_stated_in_legend": "" if n_stated is None else n_stated,
        "replicates": ";".join(f"{v:.6g}" for v in vals),
        "mean": sum(vals)/len(vals) if vals else "",
        "sd_as_published": "" if sd is None else sd,
        "value_tier": "1_publisher_source_data", "note": note})


def extended_figure_1():
    """The 100-member library, C2C12 myoblasts. Rows G1-G10 x columns Am1-Am10.

    The sheet gives RFU, Stdev and n per cell, so the mean and its spread are
    the authors' own, not something recomputed from bar heights.
    """
    rows = sheet(8, "Extended Figure 1")
    hdr_am, hdr_sub = rows[11], rows[12]
    am_at = {}
    for j, c in enumerate(hdr_am):
        if isinstance(c, str) and re.fullmatch(r"Am\d+", c.strip()):
            am_at[j] = c.strip()
    assert len(am_at) == 10, f"expected 10 Am headers, found {len(am_at)}"
    for j in am_at:
        assert [str(hdr_sub[j+k]).strip() for k in range(3)] == \
               ["RFU", "Stdev", "n"], f"column {j} is not RFU/Stdev/n"
    seen = set()
    for i in range(13, 23):
        g = str(rows[i][0]).strip()
        assert re.fullmatch(r"G\d+", g), f"row {i+1} label {g!r}"
        for j, am in am_at.items():
            rfu, sd, n = rows[i][j], rows[i][j+1], rows[i][j+2]
            lab = f"AMG{am[2:]}{g[1:]}"
            seen.add(lab)
            add("Extended Figure 1", 8, "Extended Figure 1",
                f"{col_letter(j)}{i+1}:{col_letter(j+2)}{i+1}", lab,
                "in_vitro", "C2C12 myoblast",
                "firefly luciferase, RLU (relative luminescence units)",
                "40 ng mRNA/well", "24 h", [rfu], sd=sd, n_stated=n,
                note="library-wide in vitro screen; the sheet gives the mean "
                     "and Stdev, not the wells")
    assert len(seen) == 100, f"{len(seen)} distinct lipids, expected 100"
    return seen


def figure_1f():
    """Screen 3: six lipids injected individually (Am8-Am10 x G4, G5)."""
    rows = sheet(3, "Figure 1")
    hdr = rows[25]
    assert str(hdr[0]).strip() == "Figure 1f", hdr[0]
    tails = [(2, str(hdr[2]).strip()), (5, str(hdr[5]).strip())]
    assert [t for _, t in tails] == ["Tail 4", "Tail 5"], tails
    for i in range(26, 29):
        head = str(rows[i][1]).strip()
        assert re.fullmatch(r"Head \d+", head), head
        for j, tail in tails:
            vals = rows[i][j:j+3]
            lab = f"AMG{head.split()[1]}{tail.split()[1]}"
            add("Figure 1f", 3, "Figure 1",
                f"{col_letter(j)}{i+1}:{col_letter(j+2)}{i+1}", lab,
                "in_vivo", "BALB/c mouse, whole-body IVIS, IM injection",
                "firefly luciferase, total flux [p/s]", "1 ug mRNA", "24 h",
                vals, n_stated=3)


def _block_panel(wb_n, sh, tag, model, system, readout, dose, t, hdr_row,
                 first_row, last_row, skip=()):
    rows = sheet(wb_n, sh)
    hdr = rows[hdr_row]
    assert str(hdr[0]).strip() == tag, (tag, hdr[0])
    for j, c in enumerate(hdr):
        if j == 0 or not isinstance(c, str) or not c.strip():
            continue
        lab = c.strip().replace(" ", "")
        if lab in skip:
            continue
        vals = [rows[i][j] for i in range(first_row, last_row + 1)
                if j < len(rows[i])]
        add(tag, wb_n, sh,
            f"{col_letter(j)}{first_row+1}:{col_letter(j)}{last_row+1}",
            lab, model, system, readout, dose, t, vals)


def main():
    extended_figure_1()
    figure_1f()
    # Figure 2a -- 2nd-generation tails on Am10
    _block_panel(4, "Figure 2", "2a", "in_vivo",
                 "BALB/c mouse, whole-body IVIS, IM injection",
                 "firefly luciferase, total flux [p/s]", "5 ug mRNA", "24 h",
                 0, 1, 3)
    # Figure 2b -- 3rd-generation headgroups, all carrying G41
    _block_panel(4, "Figure 2", "2b", "in_vivo",
                 "BALB/c mouse, whole-body IVIS, IM injection",
                 "firefly luciferase, total flux [p/s]", "1 ug mRNA", "6 h",
                 6, 7, 11)
    # Figure 2c -- comparison with state-of-the-art controls
    _block_panel(4, "Figure 2", "2c", "in_vivo",
                 "BALB/c mouse, whole-body IVIS, IM injection",
                 "firefly luciferase, total flux [p/s]", "1 ug mRNA", "6 h",
                 13, 14, 20)
    # Extended Figure 2 -- 4th-generation tails on Am15
    rows = sheet(9, "Extended Figure 2")
    for j, c in enumerate(rows[0]):
        if not isinstance(c, str) or not c.strip():
            continue
        add("Extended Figure 2b", 9, "Extended Figure 2",
            f"{col_letter(j)}2:{col_letter(j)}6", c.strip().replace(" ", ""),
            "in_vivo", "BALB/c mouse, whole-body IVIS, IM injection",
            "firefly luciferase, total flux [p/s]", "0.1 ug mRNA", "6 h",
            [rows[i][j] for i in range(1, 6) if j < len(rows[i])])

    cols = ["panel", "lipid_label", "model", "biological_system", "readout",
            "dose", "timepoint", "mean", "sd_as_published", "n_values",
            "n_stated_in_legend", "replicates", "value_tier", "workbook",
            "sheet", "cells", "note"]
    with open(os.path.join(HERE, "05_values.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in OUT:
            w.writerow({k: r.get(k, "") for k in cols})

    from collections import Counter
    print("values written:", len(OUT))
    for p, n in Counter(r["panel"] for r in OUT).items():
        print(f"   {p:<22} {n}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
