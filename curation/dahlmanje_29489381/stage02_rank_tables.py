#!/usr/bin/env python3
"""Extract the six RANK tables of Fig S4e/f/g from the SI text layer (tier 4).

These are the ONLY per-LNP-identified readouts anywhere in the Supporting
Information.  They are RANKS, not delivery values:

  * the in vitro column (iMAEC or RAWs) is that LNP's rank within its own
    library's stable pool;
  * the three in vivo columns are ranks WITHIN THE TOP 10% of that bed, so a
    blank is 'not in the top 10%', not 'no measurement'.

Only LNPs that reached the top 10% in vivo somewhere are listed, so every
table is a selected subset of its library.  No pixel is read anywhere.

COLUMN GEOMETRY.  Columns are NOT taken from the header: two header labels
abut ('Mass Ratio' / 'PEG Chain', measured gap 1.3 pt, below the 3.0 pt
inter-cell gap) so a header-only split silently merges them and shifts every
later column by one.  Instead the column spans are the connected components
of the union of the BODY cells' own x-intervals, and the header text is then
mapped onto those spans.  The header is a label; the body is the evidence.
"""
import csv
import json
import os
import sys

import pdfplumber

HERE = os.path.dirname(os.path.abspath(__file__))
SI = os.path.join(HERE, "source", "nl8b00432_si_001.pdf")
PAGES = [19, 20]            # 0-based: Fig S4e+f on p20, Fig S4g on p21
GLUE_PT = 0.6               # below this, glyphs abut: one token
CELL_GAP_PT = 3.0           # intra-cell word gaps measured <= 1.6 pt
ROW_TOL = 1.5               # row pitch measured 4.7-7.4 pt


def join_cell(ws):
    ws = sorted(ws, key=lambda w: w["x0"])
    out = ws[0]["text"]
    for a, b in zip(ws, ws[1:]):
        out += ("" if b["x0"] - a["x1"] < GLUE_PT else " ") + b["text"]
    return out.strip()


def _groups(ws, gap=CELL_GAP_PT):
    ws = sorted(ws, key=lambda w: w["x0"])
    groups, cur = [], [ws[0]]
    for a, b in zip(ws, ws[1:]):
        if b["x0"] - a["x1"] >= gap:
            groups.append(cur)
            cur = []
        cur.append(b)
    groups.append(cur)
    return groups


def cells_of(ws):
    return [(join_cell(g), min(w["x0"] for w in g), max(w["x1"] for w in g))
            for g in _groups(ws)]


def cluster_rows(ws, tol=ROW_TOL):
    ws = sorted(ws, key=lambda w: w["top"])
    rows, cur = [], []
    for w in ws:
        if cur and w["top"] - cur[0]["top"] > tol:
            rows.append(cur)
            cur = []
        cur.append(w)
    if cur:
        rows.append(cur)
    return rows


def column_spans(body_rows):
    """Connected components of the union of every body cell's x-interval."""
    ivs = []
    for grp in body_rows:
        ivs.extend((a, b) for _, a, b in cells_of(grp))
    ivs.sort()
    out = []
    for a, b in ivs:
        if out and a <= out[-1][1]:
            out[-1][1] = max(out[-1][1], b)
        else:
            out.append([a, b])
    return out


def truncate_at_gap(body):
    """Cut the sub-table where its own row pitch breaks.

    The last sub-table on a page has no following header to bound it, so
    without this it swallows whatever is drawn beneath -- on p21 that is Fig
    S4h's category labels, which added 13 phantom rows.  The cut is made on
    the MEASURED pitch of the rows already accepted, not on a literal y.
    """
    if len(body) < 3:
        return body
    tops = [min(w["top"] for w in g) for g in body]
    gaps = sorted(b - a for a, b in zip(tops, tops[1:]))
    med = gaps[len(gaps) // 2]
    for i, (a, b) in enumerate(zip(tops, tops[1:])):
        if b - a > 2.0 * med:
            return body[:i + 1]
    return body


def label(spans, header_ws):
    names = []
    hc = cells_of(header_ws)
    for a, b in spans:
        hit = [t for t, x0, x1 in hc if not (x1 < a or x0 > b)]
        names.append(" ".join(hit))
    return names


def unique(names):
    """Column keys must be distinct: one header label can straddle two body
    columns ('Mass RatioPEG Chain'), and a dict would then drop one of them.
    """
    out, seen = [], {}
    for i, n in enumerate(names):
        n = n or ("col%d" % i)
        seen[n] = seen.get(n, 0) + 1
        out.append(n if seen[n] == 1 else "%s #%d" % (n, seen[n]))
    return out


def which(xc, spans):
    for i, (a, b) in enumerate(spans):
        if a - 0.5 <= xc <= b + 0.5:
            return i
    return min(range(len(spans)),
               key=lambda i: min(abs(spans[i][0] - xc), abs(spans[i][1] - xc)))


def main():
    tables, rows_out = [], []
    for page0 in PAGES:
        with pdfplumber.open(SI) as pdf:
            words = pdf.pages[page0].extract_words()
        # The ACS download stamp is set VERTICALLY in the right margin.
        # Dropped on its orientation, not on an x threshold: the library-2
        # tables run to x=584 and a margin guess at 560 silently ate their
        # last column ('Marrow Ecs').
        words = [w for w in words if w.get("upright", 1)]
        rows = cluster_rows(words)
        hidx = [i for i, g in enumerate(rows)
                if any(w["text"] == "Diameter" for w in g)]
        for k, hi in enumerate(hidx):
            end = hidx[k + 1] - 1 if k + 1 < len(hidx) else len(rows)
            sub = rows[hi + 1:end]
            body = [g for g in sub
                    if join_cell(g) not in ("Library 1", "Library 2",
                                            "Library 3", "f", "g")]
            body = truncate_at_gap(body)
            if not body:
                continue
            spans = column_spans(body)
            names = label(spans, rows[hi])
            uniq = unique(names)
            lab = join_cell(rows[hi - 1]) if hi else ""
            recs = []
            for grp in body:
                cells, dup = {}, 0
                for txt, a, b in cells_of(grp):
                    i = which((a + b) / 2.0, spans)
                    if i in cells:
                        dup += 1
                        cells[i] += " " + txt
                    else:
                        cells[i] = txt
                rec = {"si_page": page0 + 1, "subtable": lab}
                for i, n in enumerate(uniq):
                    rec[n] = cells.get(i, "")
                rec["_dup_cells"] = dup
                recs.append(rec)
            tables.append({"si_page": page0 + 1, "label": lab,
                           "columns": uniq, "n_columns": len(names),
                           "n_rows": len(recs),
                           "rows_with_two_cells_on_one_column":
                               sum(1 for r in recs if r["_dup_cells"])})
            rows_out.extend(recs)

    keys = []
    for r in rows_out:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(os.path.join(HERE, "02_rank_tables.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, restval="")
        w.writeheader()
        for r in rows_out:
            w.writerow(r)
    rep = {"tables": tables, "n_rows_total": len(rows_out)}
    with open(os.path.join(HERE, "02_rank_tables.json"), "w") as fh:
        json.dump(rep, fh, indent=1)
    print(json.dumps(rep, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
