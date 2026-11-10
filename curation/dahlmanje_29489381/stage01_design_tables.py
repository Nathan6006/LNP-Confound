#!/usr/bin/env python3
"""Extract the four LNP library DESIGN tables from the SI text layer (tier 4).

These are real text tables -- pdfplumber returns every cell as a positioned
word, so no pixel is read anywhere in this file.  Columns are recovered from
the header word positions on the same page, never from a hard-coded x.

Pages (0-based): 6 = Fig S2a (library 1), 11 = Fig S3d (library 2),
12 = Fig S3f (library 3), 22 = Fig S5a (library 4).

Chemistry glyphs from the reaction scheme sit on the same pages, so the
extractor is restricted to the x-range of the table's own header words and
rows are keyed on the LNP # column being a plain integer.
"""
import json
import os
import sys

import pdfplumber

HERE = os.path.dirname(os.path.abspath(__file__))
SI = os.path.join(HERE, "source", "nl8b00432_si_001.pdf")

# (page0, library, n_rows_expected, header tokens in order)
PAGES = [
    (6, 1, 144),
    (11, 2, 120),
    (12, 3, 156),
    (22, 4, 144),
]


def header_top(words):
    """Top coordinate of the row containing the literal header cell 'LNP #'.

    Found, not assumed: the panel letter ('a') also sits above the table, and
    hard-coding a y band is exactly the class of defect that a detected
    anchor avoids.
    """
    cands = [w for w in words if w["text"] == "LNP"]
    for w in cands:
        nbr = [v for v in words
               if abs(v["top"] - w["top"]) < 1.0 and 0 <= v["x0"] - w["x1"] < 4.0]
        if any(v["text"] == "#" for v in nbr):
            return w["top"]
    raise SystemExit("header cell 'LNP #' not found on page")


def header_cols(words, htop):
    """Column boundaries from the header words' own x positions."""
    hdr = [w for w in words if abs(w["top"] - htop) < 1.0]
    hdr.sort(key=lambda w: w["x0"])
    # group header words into column labels by x gap
    cols, cur = [], []
    for w in hdr:
        if cur and w["x0"] - cur[-1]["x1"] > 3.0:
            cols.append(cur)
            cur = []
        cur.append(w)
    if cur:
        cols.append(cur)
    out = []
    for c in cols:
        out.append({
            "name": " ".join(w["text"] for w in c),
            "x0": min(w["x0"] for w in c),
            "x1": max(w["x1"] for w in c),
        })
    return out


GLUE_PT = 0.6   # below this inter-word gap the glyphs abut: one token


def join_cell(ws):
    """Rebuild a cell from positioned words.

    Some rows of these tables are emitted glyph-by-glyph by the PDF producer,
    so '104' arrives as three 'words' with zero gap while 'PEI' and '600' sit
    ~1.0 pt apart on every row.  Joining on the measured gap recovers both
    without a per-row special case.
    """
    ws = sorted(ws, key=lambda w: w["x0"])
    out = ws[0]["text"]
    for a, b in zip(ws, ws[1:]):
        out += ("" if b["x0"] - a["x1"] < GLUE_PT else " ") + b["text"]
    return out.strip()


CELL_GAP_PT = 3.0   # inter-cell gaps measured >= 9 pt; intra-cell <= 1.1 pt


def row_cells(ws):
    """Group one row's words into CELLS on the measured inter-word gap.

    Not on the header midpoints: cell text is centred and can be wider than
    its own header label, so 'C15 Epoxy' spills its final glyph past the
    midpoint into the next column.  Grouping first, then matching each cell
    to the nearest header centre, never splits a cell.
    """
    ws = sorted(ws, key=lambda w: w["x0"])
    groups, cur = [], [ws[0]]
    for a, b in zip(ws, ws[1:]):
        if b["x0"] - a["x1"] >= CELL_GAP_PT:
            groups.append(cur)
            cur = []
        cur.append(b)
    groups.append(cur)
    return [(join_cell(g),
             (min(w["x0"] for w in g) + max(w["x1"] for w in g)) / 2.0)
            for g in groups]


def nearest(xc, cols):
    centres = [(c["x0"] + c["x1"]) / 2.0 for c in cols]
    return min(range(len(cols)), key=lambda i: abs(centres[i] - xc))


def extract(page0, library, n_expected):
    with pdfplumber.open(SI) as pdf:
        pg = pdf.pages[page0]
        words = pg.extract_words()
    htop = header_top(words)
    cols = header_cols(words, htop)
    if not cols or not cols[0]["name"].startswith("LNP"):
        raise SystemExit(
            "p%d: header not found; first column is %r"
            % (page0 + 1, cols[0]["name"] if cols else None))
    xlo, xhi = cols[0]["x0"] - 6, cols[-1]["x1"] + 6

    # group body words into rows by their top coordinate
    body = [w for w in words
            if w["top"] > htop + 1.0 and xlo <= w["x0"] and w["x1"] <= xhi]
    # Cluster words into rows on y with a tolerance, NOT on a rounded top.
    # Measured row pitch is 4.7-5.1 pt and individual cells carry up to 0.4 pt
    # of baseline jitter, so rounding to 0.1 pt splits one row in two and the
    # halves then each lose cells to the other.
    ROW_TOL = 1.5
    body.sort(key=lambda w: w["top"])
    rows, cur = [], []
    for w in body:
        if cur and w["top"] - cur[0]["top"] > ROW_TOL:
            rows.append(cur)
            cur = []
        cur.append(w)
    if cur:
        rows.append(cur)

    recs, rejected = [], 0
    for grp in rows:
        ws = sorted(grp, key=lambda w: w["x0"])
        cells = {}
        collisions = 0
        for txt, xc in row_cells(ws):
            i = nearest(xc, cols)
            if i in cells:
                collisions += 1
                cells[i] = cells[i] + " " + txt
            else:
                cells[i] = txt
        first = cells.get(0, "")
        if not first.isdigit():
            rejected += 1
            continue
        rec = {"library": library, "si_page": page0 + 1}
        for i, c in enumerate(cols):
            rec[c["name"]] = cells.get(i, "")
        rec["_cell_collisions"] = collisions
        recs.append(rec)
    return cols, recs, rejected


def main():
    out = {}
    allrecs = []
    for page0, lib, n_exp in PAGES:
        cols, recs, rej = extract(page0, lib, n_exp)
        nums = [int(r["LNP #"]) for r in recs]
        out["library_%d" % lib] = {
            "si_page": page0 + 1,
            "columns": [c["name"] for c in cols],
            "n_rows_extracted": len(recs),
            "n_rows_stated": n_exp,
            "lnp_numbers_are_1_to_n": nums == list(range(1, len(nums) + 1)),
            "rows_rejected_non_integer_first_cell": rej,
            "rows_with_two_cells_on_one_column": sum(
                1 for r in recs if r["_cell_collisions"]),
        }
        allrecs.extend(recs)
    with open(os.path.join(HERE, "01_design_tables.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    import csv
    keys = []
    for r in allrecs:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(os.path.join(HERE, "01_design_tables.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        for r in allrecs:
            w.writerow(r)
    print(json.dumps(out, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
