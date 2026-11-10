#!/usr/bin/env python3
"""Stage-00 triage: does this source contain an extractable LIBRARY SCREEN?

    python3 -I triage.py <folder> <file> [file...]

Reports, per file: page count, how much real text, how many compound-like
labels, whether any page looks like a data TABLE of library members, and how
many numeric values sit near those labels. The point is to decide, cheaply and
before any extraction effort, whether a paper can yield >=100 usable rows --
the run's stop threshold.
"""
import os, re, sys, json, collections

LABEL = re.compile(r'\b(?:[A-Z]{1,3}\d{1,3}[-_][A-Za-z0-9]{1,8}|'
                   r'[A-Z]{1,4}-?\d{1,4}[A-Za-z]?|'
                   r'\d{1,2}[A-Z]\d{1,2}[A-Za-z0-9]*)\b')
NUM = re.compile(r'\d+\.\d+|\d{3,}')
HDRISH = re.compile(r'(compound|lipid|entry|yield|structure|MW|m/z|HRMS|'
                    r'formulation|RLU|luminescen|transfect|efficien|IC50|'
                    r'%\s*encaps|size|PDI|zeta|pKa)', re.I)


def pdf_report(path):
    import pymupdf
    d = pymupdf.open(path)
    pages = d.page_count
    txt_pages, labels, nums, tabley = 0, collections.Counter(), 0, []
    for i, p in enumerate(d):
        t = p.get_text()
        if len(t.strip()) > 200:
            txt_pages += 1
        ls = LABEL.findall(t)
        labels.update(ls)
        nums += len(NUM.findall(t))
        if len(set(ls)) >= 8 and HDRISH.search(t):
            tabley.append(i + 1)
    d.close()
    out = {"pages": pages, "text_pages": txt_pages,
           "distinct_labels": len(labels), "numbers": nums,
           "table_pages": tabley[:25], "n_table_pages": len(tabley),
           "top_labels": [l for l, _ in labels.most_common(8)]}
    out["numbers_regex"] = NUM.pattern   # see native_and_ruled() docstring
    nat, ruled = native_and_ruled(path)
    out["native_tables"] = nat
    out["ruled_tables"] = ruled
    return out


def native_and_ruled(path, max_pages=120):
    """Both table answers, because MY TEXT HEURISTIC HAS BEEN WRONG BOTH WAYS.

    Measured 2026-11-10 on `adma201906128-sup-0001-SuppMat.pdf`: this file
    reported **"0 table-like pages"** while `pdfplumber.find_tables()` succeeds
    on two of them -- p47 a 12x14 and p48 a 12x17, matching a curation agent's
    hand count of Tables S1/S2 exactly. The heuristic above needs >=8 distinct
    compound-like labels AND a header-ish word on the same page, so a table of
    toxicology rows passes neither test and vanishes. I put that false zero into
    a reviewer's brief.

    The opposite error is in `ruled_tables` below: `find_tables()` returning 0 on
    a page that is nothing but table. So the honest report is BOTH numbers plus
    the text heuristic, each labelled, and none of them trusted alone.

    The token count has the same disease. This file reported "1,745 numeric
    tokens"; the same text gives **2,533 / 3,123 / 6,958** under three other
    perfectly reasonable regexes. A count whose definition is not printed beside
    it is not quotable -- the project has already been through this with the
    flagged-row count that moved four times (F-024). `numbers_regex` now ships
    with the number.
    """
    try:
        import pdfplumber
    except Exception:
        return {"available": False}, {"available": False}
    native, found = [], []
    try:
        with pdfplumber.open(path) as pdf:
            for i, pg in enumerate(pdf.pages[:max_pages]):
                try:
                    t = pg.find_tables()
                    if t:
                        native.append({"page": i + 1, "n": len(t),
                                       "shapes": [
                                           f"{len(x.rows)}x"
                                           f"{len(x.rows[0].cells) if x.rows else 0}"
                                           for x in t[:4]]})
                        continue
                    found += _ruled_page(pg, i)
                except Exception:
                    continue
    except Exception as e:
        err = {"available": True, "error": f"{type(e).__name__}: {e}"}
        return err, err
    return ({"available": True, "pages": native, "n": len(native),
             "how": "pdfplumber.find_tables() succeeded on these pages"},
            {"available": True, "pages": found, "n": len(found),
             "how": ("find_tables() returned 0 here, but the rules are FILLED "
                     "RECTS, not line objects. Re-run with "
                     "horizontal_strategy='explicit' and the rect runs as "
                     "explicit_horizontal_lines.") if found else
                    "no page has a ruled table that find_tables() missed"})


def _ruled_page(pg, i):
    thin = [r for r in pg.rects if 0 < (r["y1"] - r["y0"]) < 1.0]
    if len(thin) < 6:
        return []
    tops = sorted({round(r["y0"], 1) for r in thin})
    if len(tops) < 3:
        return []
    try:
        t = pg.find_tables({"horizontal_strategy": "explicit",
                            "explicit_horizontal_lines": tops,
                            "vertical_strategy": "text"})
    except Exception:
        return []
    if not t:
        return []
    rows = len(t[0].rows)
    cols = len(t[0].rows[0].cells) if t[0].rows else 0
    try:
        cells = t[0].extract()
    except Exception:
        return []
    filled = sum(1 for r in cells for c in r if c and str(c).strip())
    numeric = sum(1 for r in cells for c in r
                  if c and any(ch.isdigit() for ch in str(c)))
    if filled < max(8, rows) or numeric < 4:
        return []
    return [{"page": i + 1, "find_tables": 0, "thin_rects": len(thin),
             "recovered": f"{rows}x{cols}", "filled_cells": filled,
             "numeric_cells": numeric,
             "rule_heights": sorted({round(r["y1"] - r["y0"], 2)
                                     for r in thin})[:4]}]


def ruled_tables(path, max_pages=120):
    """Pages holding a RULED table that `pdfplumber.find_tables()` reports as 0.

    WHY THIS EXISTS. Measured on `bc400158w_si_001.pdf` (LL_2013), 2026-11-10:
    pages 18-19 carry every delivery value in the paper, as pure text objects,
    and `find_tables()` returns **0** on both. The rules are not `line` objects
    at all -- they are **filled rectangles 0.36-0.48 pt tall**, so the default
    "lines" strategy finds no edges. Passing the rect runs as explicit
    horizontal lines recovers a 7x13 and a 3x7 table.

        p18: find_tables()=0  words=245  images=0  thin rects=80
             with rect runs as explicit h-lines: 1 table, 7 rows x 13 cols

    Two false zeros compound here and both are silent:
      * `find_tables()` says "no table" on a page that is nothing but a table;
      * that SI's own List of Contents stops at Figure S2 and names none of the
        four value tables, so a triage that trusts the contents page agrees.

    A scanner that cannot see a structure reports a clean run -- the same shape
    as `find_screens` hiding a 12-lipid panel in column W, and `batch_gate`
    checking zero structures in four folders. So this does not replace the text
    heuristic above, it runs beside it and reports what that one cannot see.
    """
    try:
        import pdfplumber
    except Exception:
        return {"available": False}
    found = []
    try:
        with pdfplumber.open(path) as pdf:
            for i, pg in enumerate(pdf.pages[:max_pages]):
                try:
                    native = len(pg.find_tables())
                    thin = [r for r in pg.rects if 0 < (r["y1"] - r["y0"]) < 1.0]
                    if native or len(thin) < 6:
                        continue
                    tops = sorted({round(r["y0"], 1) for r in thin})
                    if len(tops) < 3:
                        continue
                    t = pg.find_tables({"horizontal_strategy": "explicit",
                                        "explicit_horizontal_lines": tops,
                                        "vertical_strategy": "text"})
                    if not t:
                        continue
                    rows = len(t[0].rows)
                    cols = len(t[0].rows[0].cells) if t[0].rows else 0
                    # A RECOVERED GEOMETRY IS NOT A TABLE. The first version of
                    # this check reported `sd_2026` p19 as a 31x11 table from 197
                    # thin rects; extracting it gave 31 rows of entirely EMPTY
                    # cells. Those rects are figure graphics, not rules. Scoring
                    # the grid and not its contents is the same defect this
                    # function exists to catch -- a scanner printing a plausible
                    # number -- so require real text in real cells.
                    try:
                        cells = t[0].extract()
                    except Exception:
                        continue
                    filled = sum(1 for r in cells for c in r
                                 if c and str(c).strip())
                    numeric = sum(1 for r in cells for c in r
                                  if c and any(ch.isdigit() for ch in str(c)))
                    if filled < max(8, rows) or numeric < 4:
                        continue
                    found.append({"page": i + 1, "find_tables": native,
                                  "thin_rects": len(thin),
                                  "recovered": f"{rows}x{cols}",
                                  "filled_cells": filled,
                                  "numeric_cells": numeric,
                                  "rule_heights": sorted(
                                      {round(r["y1"] - r["y0"], 2)
                                       for r in thin})[:4]})
                except Exception:
                    continue
    except Exception as e:
        return {"available": True, "error": f"{type(e).__name__}: {e}"}
    return {"available": True, "pages": found, "n": len(found),
            "how": ("pdfplumber find_tables() returned 0 on these pages, but "
                    "their rules are FILLED RECTS, not line objects. Re-run "
                    "with horizontal_strategy='explicit' and the rect runs as "
                    "explicit_horizontal_lines.") if found else
                   "no page has a ruled table that find_tables() missed"}


def docx_report(path):
    import zipfile
    with zipfile.ZipFile(path) as z:
        xml = z.read("word/document.xml").decode("utf-8", "replace")
    txt = re.sub(r'<[^>]+>', ' ', xml)
    ls = LABEL.findall(txt)
    return {"kind": "docx", "tables": xml.count("<w:tbl>"),
            "distinct_labels": len(set(ls)), "numbers": len(NUM.findall(txt)),
            "top_labels": [l for l, _ in collections.Counter(ls).most_common(8)]}


def discover(folder):
    """Every PDF/DOCX under `folder`, including a `source/` subdirectory.

    `new_curation.py` scaffolds sources into `<folder>/source/`, so a caller who
    passes only the folder used to get `{"files": {}}` -- a clean, empty,
    entirely silent report. That is the same false-zero shape this file's
    `ruled_tables` exists to catch, committed by the triage tool itself.
    """
    out = []
    for root, _dirs, files in os.walk(folder):
        for f in sorted(files):
            if f.lower().endswith((".pdf", ".docx")):
                out.append(os.path.join(root, f))
    return out


def main():
    folder = sys.argv[1]
    paths = sys.argv[2:]
    if not paths:
        if os.path.isdir(folder):
            paths = discover(folder)
        elif folder.lower().endswith((".pdf", ".docx")):
            paths, folder = [folder], os.path.dirname(folder)  # a file, not a dir
    if not paths:
        print(json.dumps({"folder": folder, "files": {},
                          "error": "NO PDF OR DOCX FOUND. Refusing to print an "
                                   "empty report as if it were a result -- pass "
                                   "the files explicitly, or check that the "
                                   "sources are really on disk."}))
        return 1
    out = {"folder": folder, "files": {}}
    for p in paths:
        name = os.path.basename(p)
        try:
            if p.lower().endswith(".pdf"):
                out["files"][name] = pdf_report(p)
            elif p.lower().endswith(".docx"):
                out["files"][name] = docx_report(p)
        except Exception as e:
            out["files"][name] = {"error": str(e)[:100]}
    print(json.dumps(out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
