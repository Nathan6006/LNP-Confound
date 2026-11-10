#!/usr/bin/env python3
"""Stage 05 -- read Tables S3, S4, S5, S6 out of the SI's TEXT LAYER.

TIER 4, NOT TIER 5. Every value this writes is a PDF text object. Nothing here
digitises anything; there is no axis, no colour ramp and no bar, so none of the
geometric error classes apply and curation/lib's pixel modules are deliberately
not imported.

WHY THIS FILE EXISTS AT ALL, GIVEN find_tables() RETURNS 0.
pdfplumber's default table strategy looks for `line` objects. This SI draws its
rules as FILLED RECTS 0.36-0.48 pt tall, so the default finds no edges and
reports zero tables on pages 4, 5, 18 and 19. The text is all there. The right
response is to hand pdfplumber the rects as explicit edges, which is what this
file does -- NOT to reach for pixels.

THE DEFECT THIS FILE IS BUILT TO AVOID.
Table S4 (p.18) carries FOUR (efficiency, SD) column pairs but only THREE
labelled lipidoid:DOPE ratios. A row-wise read that simply takes the numbers in
order would assign the L2K control's "14.0 1.8" to the 2:1 ratio. MEASURED: the
L2K row's two numbers sit at x0 = 525.8 and 557.2, inside the FOURTH, unlabelled
band (506.16-545.64 and 546.00-572.28), while every lipidoid row's first number
sits near x0 = 211.7, inside the FIRST band (133.08-226.68). The fourth pair is
an L2K-only column. So every cell here is placed by its MEASURED x against
edges taken from the drawn rules, and a word that straddles an edge raises.

Run:  python3 -I stage05_values.py
"""
import csv
import json
import os
import sys

import pdfplumber

HERE = os.path.dirname(os.path.abspath(__file__))
PDF = os.path.join(HERE, "source", "bc400158w_si_001.pdf")

# ---------------------------------------------------------------------------
# Table geometry. Page index (0-based) and the y of the table's TOP RULE, which
# is the rect run whose segmentation defines the columns. Both are read off the
# PDF, never guessed: see 07_qc.json check Q-1.
# ---------------------------------------------------------------------------
TABLES = [
    dict(id="S3", page=17, top_rule=92.64, body=(107.64, 290.88),
         caption="corresponds to Figure 2A",
         cols=["name", "eff", "sd"], n_rows=13),
    dict(id="S4", page=17, top_rule=330.24, body=(376.80, 561.24),
         caption="corresponds to Figure 3",
         cols=["name", "eff_2:1", "sd_2:1", "eff_1:1", "sd_1:1",
               "eff_1:2", "sd_1:2", "eff_ctrl", "sd_ctrl"], n_rows=13),
    # Table S5 spans a page break: the 75 ng block starts on p.18 and the rest
    # is on p.19. Its two page-parts share one column geometry, asserted equal.
    dict(id="S5a", page=17, top_rule=597.96, body=(612.60, 792.0),
         caption="corresponds to Figure 6",
         cols=["dna", "formulation", "eff", "sd"], n_rows=8),
    dict(id="S5b", page=18, top_rule=376.92, body=(0.0, 376.0),
         caption="corresponds to Figure 6 (continuation)",
         cols=["dna", "formulation", "eff", "sd"], n_rows=22),
    dict(id="S6", page=18, top_rule=416.28, body=(430.80, 529.92),
         caption="corresponds to Figure 8",
         cols=["name", "eff", "sd"], n_rows=7),
]

RULE_MAX_H = 1.0      # a rule is a rect thinner than this
RULE_MIN_W = 5.0      # a column band is wider than this; narrower = a gap
Y_TOL = 2.0           # rects within this of the stated rule y belong to it
ROW_TOL = 3.0         # words within this in `top` are one row
WATERMARK_X = 580.0   # the ACS download stamp runs down the right margin


def column_bands(page, rule_y):
    """Column x-ranges, from the drawn rule's own segmentation.

    The rule is a run of filled rects: wide ones are the column bands, the
    0.36-0.48 pt ones between them are the vertical-rule gaps. Refuses rather
    than returning a guess.
    """
    segs = []
    for r in page.rects:
        if (r["y1"] - r["y0"]) > RULE_MAX_H:
            continue
        if abs(r["top"] - rule_y) > Y_TOL:
            continue
        if (r["x1"] - r["x0"]) < RULE_MIN_W:
            continue                      # a vertical-rule gap, not a column
        segs.append((round(r["x0"], 2), round(r["x1"], 2)))
    segs.sort()
    if not segs:
        raise RuntimeError(
            f"column_bands: no rule rects within {Y_TOL} pt of y={rule_y}. "
            "This file will not fall back to a hard-coded column layout: "
            "error class 28 is exactly a detector miss falling through to a "
            "literal.")
    return segs


def assign(page, bands, y0, y1, ncols, tid):
    """Words in (y0, y1] -> rows of `ncols` cells, placed by measured x.

    A word whose box straddles a band edge raises: that is the only way a cell
    can be silently misassigned, and it is the defect in the docstring.
    """
    if len(bands) != ncols:
        raise RuntimeError(
            f"{tid}: the drawn rule segments into {len(bands)} column bands "
            f"but the declared layout has {ncols}. The rule is the authority; "
            "fix the declaration, do not pad the bands.")
    ws = [w for w in page.extract_words()
          if y0 < w["top"] <= y1 and w["x0"] < WATERMARK_X]
    rows = {}
    for w in ws:
        hit = [i for i, (bx0, bx1) in enumerate(bands)
               if w["x0"] >= bx0 - 1.0 and w["x1"] <= bx1 + 1.0]
        if len(hit) != 1:
            # Not a straddle if it simply falls outside every band (a stray
            # caption word); a straddle is a word overlapping two bands.
            over = [i for i, (bx0, bx1) in enumerate(bands)
                    if w["x1"] > bx0 and w["x0"] < bx1]
            if len(over) > 1:
                raise RuntimeError(
                    f"{tid}: word {w['text']!r} at x=[{w['x0']:.1f},"
                    f"{w['x1']:.1f}] straddles column edges {over}. Refusing "
                    "to place it by nearest centre -- that is how the L2K "
                    "control would land under a DOPE ratio it was never "
                    "measured at.")
            continue
        key = None
        for k in rows:
            if abs(k - w["top"]) <= ROW_TOL:
                key = k
                break
        if key is None:
            key = w["top"]
            rows[key] = [[] for _ in range(ncols)]
        rows[key][hit[0]].append((w["x0"], w["text"]))
    out = []
    for k in sorted(rows):
        cells = []
        for cell in rows[k]:
            cell.sort()
            cells.append(" ".join(t for _, t in cell))
        out.append((round(k, 2), cells))
    return out


def main():
    doc = pdfplumber.open(PDF)
    raw = {}
    geom = {}
    for t in TABLES:
        pg = doc.pages[t["page"]]
        bands = column_bands(pg, t["top_rule"])
        rows = assign(pg, bands, t["body"][0], t["body"][1], len(t["cols"]),
                      t["id"])
        if len(rows) != t["n_rows"]:
            raise RuntimeError(
                f"{t['id']}: read {len(rows)} body rows, declared "
                f"{t['n_rows']}. A row-count mismatch is not cosmetic -- a "
                "dropped row on a page break is how Table S5's 1C12-D 75 ng "
                "pair would vanish silently.")
        raw[t["id"]] = rows
        geom[t["id"]] = dict(page=t["page"] + 1, top_rule=t["top_rule"],
                             bands=bands, n_rows=len(rows), cols=t["cols"])
    doc.close()

    # --- S5's two page-parts must share one column geometry ----------------
    ba, bb = geom["S5a"]["bands"], geom["S5b"]["bands"]
    drift = max(abs(a[0] - b[0]) for a, b in zip(ba, bb))
    drift = max(drift, max(abs(a[1] - b[1]) for a, b in zip(ba, bb)))
    if drift > 1.5:
        raise RuntimeError(
            f"S5: the two page-parts' column bands differ by {drift:.2f} pt. "
            "They are one table and must share one geometry.")
    geom["S5_page_break_band_drift_pt"] = round(drift, 3)

    # --- flatten to 05_values.csv ------------------------------------------
    recs = []

    def num(s):
        s = s.strip()
        if s == "":
            return ""
        return s                    # kept AS PRINTED; never reformatted

    for ypos, c in raw["S3"]:
        recs.append(dict(
            value_id=f"S3-{c[0]}", source_table="Table S3",
            source_page=18, source_row_y=ypos,
            lipidoid=c[0], formulation_lipidoid_DOPE="not stated",
            dna_ng="not printed (see Q-7)",
            readout="average transfection efficiency (percentage of transfected cells)",
            unit="percent", value_as_printed=num(c[1]), sd_as_printed=num(c[2]),
            corresponds_to="Figure 2A"))

    ratio_cols = [("2:1", 1, 2), ("1:1", 3, 4), ("1:2", 5, 6)]
    for ypos, c in raw["S4"]:
        name = c[0]
        if name == "L2K":
            # Measured into the fourth, UNLABELLED band -- not a DOPE ratio.
            recs.append(dict(
                value_id="S4-L2K", source_table="Table S4", source_page=18,
                source_row_y=ypos, lipidoid="L2K",
                formulation_lipidoid_DOPE="n/a (control; its numbers sit in "
                "the fourth unlabelled column pair, x0=525.8/557.2)",
                dna_ng="not printed (see Q-7)",
                readout="average transfection efficiency (percentage of transfected cells)",
                unit="percent", value_as_printed=num(c[7]),
                sd_as_printed=num(c[8]), corresponds_to="Figure 3"))
            for i in (1, 2, 3, 4, 5, 6):
                if c[i].strip():
                    raise RuntimeError(
                        "S4 L2K: a ratio column is populated; the control was "
                        "supposed to be measured at no DOPE ratio.")
            continue
        for lab, ie, isd in ratio_cols:
            recs.append(dict(
                value_id=f"S4-{name}-{lab}", source_table="Table S4",
                source_page=18, source_row_y=ypos, lipidoid=name,
                formulation_lipidoid_DOPE=lab,
                dna_ng="not printed (see Q-7)",
                readout="average transfection efficiency (percentage of transfected cells)",
                unit="percent", value_as_printed=num(c[ie]),
                sd_as_printed=num(c[isd]), corresponds_to="Figure 3"))
        if c[7].strip() or c[8].strip():
            raise RuntimeError(f"S4 {name}: the control column is populated.")

    # -----------------------------------------------------------------------
    # Table S5's DNA amount is a MERGED CELL spanning its whole 10-row block,
    # so the label glyph sits on whichever body row is nearest the block's
    # VERTICAL CENTRE -- measured: row 6 of 10, every time.
    #
    # THE FIRST VERSION OF THIS BLOCK WAS WRONG AND THE ERROR WAS SILENT.
    # It grouped rows by "flush when a DNA label appears", which treats the
    # label as the block's LAST row. That shifted every block four rows: the
    # 1C12-D 75 ng values (32.6 / 40.3 / 3.4) were stamped 150 ng, the real
    # 150 ng rows were stamped 200 ng, and the output still looked like a
    # clean 30-row table. It was caught by reading the emitted CSV against the
    # page, not by anything the parser checked -- so the grouping is now
    # STRUCTURAL and asserted.
    #
    # A block runs from '1C16-D:DOPE=2:1' to 'L2K' inclusive. Each must be
    # exactly 10 rows and carry exactly one DNA label, on its 6th row. The
    # dose is the single most load-bearing metadata field in this table, and a
    # four-row shift silently moves 9 of 30 values onto the wrong dose.
    # -----------------------------------------------------------------------
    s5 = ([(18, y, c) for y, c in raw["S5a"]]
          + [(19, y, c) for y, c in raw["S5b"]])
    blocks, cur = [], []
    for item in s5:
        c = item[2]
        if c[1].strip() == "1C16-D:DOPE=2:1":
            if cur:
                raise RuntimeError(
                    "S5: a new block opened before the previous one closed "
                    "on an L2K row.")
            cur = [item]
        else:
            if not cur:
                raise RuntimeError(
                    f"S5: row {c[1]!r} appears before any block opened.")
            cur.append(item)
            if c[1].strip() == "L2K":
                blocks.append(cur)
                cur = []
    if cur:
        raise RuntimeError("S5: the last block never closed on an L2K row.")
    if len(blocks) != 3:
        raise RuntimeError(f"S5: {len(blocks)} dose blocks, expected 3.")
    s5_blocks = []
    for bi, blk in enumerate(blocks):
        if len(blk) != 10:
            raise RuntimeError(
                f"S5 block {bi}: {len(blk)} rows, expected 10 "
                "(3 lipidoids x 3 DOPE ratios + L2K).")
        labs = [(i, it[2][0].strip()) for i, it in enumerate(blk)
                if it[2][0].strip()]
        if len(labs) != 1:
            raise RuntimeError(
                f"S5 block {bi}: {len(labs)} DNA labels {labs}, expected 1. "
                "The dose is a merged cell; two labels or none means the "
                "block boundaries are wrong.")
        idx, dna = labs[0]
        if idx != 5:
            raise RuntimeError(
                f"S5 block {bi}: the DNA label {dna!r} sits on row {idx + 1} "
                "of 10, not row 6. The merged cell's glyph is vertically "
                "centred in its block, so row 6 is the signature of correct "
                "block boundaries; anything else means this parser has "
                "regrouped the table.")
        s5_blocks.append(dict(block=bi, dna=dna, label_row_of_10=idx + 1,
                              rows=[it[2][1] for it in blk]))
        for pn, yy, cc in blk:
            recs.append(dict(
                value_id=f"S5-{dna.replace(' ', '')}-{cc[1]}",
                source_table="Table S5", source_page=pn, source_row_y=yy,
                lipidoid=cc[1].split(":DOPE=")[0],
                formulation_lipidoid_DOPE=(
                    cc[1].split(":DOPE=")[1] if ":DOPE=" in cc[1]
                    else "n/a (control)"),
                dna_ng=dna,
                readout="average transfection efficiency (percentage of transfected cells)",
                unit="percent", value_as_printed=num(cc[2]),
                sd_as_printed=num(cc[3]), corresponds_to="Figure 6"))
    if [b["dna"] for b in s5_blocks] != ["75 ng", "150 ng", "200 ng"]:
        raise RuntimeError(
            f"S5: dose order {[b['dna'] for b in s5_blocks]} is not "
            "75/150/200 ng as printed.")

    for ypos, c in raw["S6"]:
        recs.append(dict(
            value_id=f"S6-{c[0]}", source_table="Table S6", source_page=19,
            source_row_y=ypos, lipidoid=c[0],
            formulation_lipidoid_DOPE="not stated for this panel",
            dna_ng="not printed (see Q-7)",
            readout="average transfection efficiency (percentage of transfected cells)",
            unit="percent", value_as_printed=num(c[1]), sd_as_printed=num(c[2]),
            corresponds_to="Figure 8"))

    cols = ["value_id", "source_table", "source_page", "source_row_y",
            "lipidoid", "formulation_lipidoid_DOPE", "dna_ng", "readout",
            "unit", "value_as_printed", "sd_as_printed", "corresponds_to"]
    with open(os.path.join(HERE, "05_values.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in recs:
            w.writerow(r)

    # every value must be a parseable number as printed
    bad = [r["value_id"] for r in recs
           if not r["value_as_printed"].replace(".", "", 1).isdigit()]
    if bad:
        raise RuntimeError(f"non-numeric values: {bad}")

    summary = dict(
        _what="Stage 05. Where every value came from, read from the PDF text "
              "layer with column edges taken from the drawn rules.",
        tier="4_pdf_text_table -- zero pixels, zero digitisation",
        n_values=len(recs),
        by_table={t: sum(1 for r in recs if r["source_table"] ==
                         ("Table S5" if t.startswith("S5") else f"Table {t}"))
                  for t in ["S3", "S4", "S5", "S6"]},
        geometry=geom,
        S5_dose_blocks=s5_blocks,
        distinct_lipidoid_labels=sorted({r["lipidoid"] for r in recs}))
    with open(os.path.join(HERE, "05_values_geometry.json"), "w") as fh:
        json.dump(summary, fh, indent=1)

    print(f"OK  {len(recs)} values written to 05_values.csv")
    print(f"    S5 page-break band drift = "
          f"{geom['S5_page_break_band_drift_pt']} pt")
    for t in TABLES:
        print(f"    {t['id']}: {geom[t['id']]['n_rows']} rows, "
              f"{len(geom[t['id']]['bands'])} column bands")
    return 0


if __name__ == "__main__":
    sys.exit(main())
