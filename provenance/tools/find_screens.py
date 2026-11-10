#!/usr/bin/env python3
"""Flag spreadsheet sheets that look like a LIPID LIBRARY SCREEN.

Why: cell count is not library size. A 1,975-row Source Data sheet in
41551_2026_1640_MOESM11 turned out to be a molecular-dynamics time series
(column header "Time (ns)"), not a screen. A sheet is interesting here only if
it has many DISTINCT COMPOUND-LIKE LABELS, not merely many numbers.
"""
import re, sys, os
from openpyxl import load_workbook

TIMEISH = re.compile(r'\b(time|ns\)|min\)|hour|sec|wavelength|nm\)|cycle|frame|'
                     r'temperature|pH\b|size \(|zeta|PDI|diameter)', re.I)
LIPIDISH = re.compile(r'^(?:[A-Z]{1,4}[-_ ]?\d{1,4}[A-Za-z0-9\-_]*|'
                      r'[A-Za-z]{1,3}\d{1,3}[-_][A-Za-z0-9]{1,6}.*|'
                      r'Lipid\s*\d+|L\d+|\d+[A-Z]\d+.*)$')

def classify(ws):
    rows = list(ws.iter_rows(values_only=True))
    if not rows: return None
    strs = {}
    for r in rows:
        for j, c in enumerate(r):
            if isinstance(c, str) and 1 < len(c.strip()) < 30:
                strs.setdefault(j, []).append(c.strip())
    # Report EVERY qualifying column, not just the best one, and at a lower
    # threshold. The first version kept only the single richest column per
    # sheet at >=15 labels, and so missed `Fig. 3i` -- the largest in-vivo
    # panel in 41551_2026_1640_MOESM9 -- because its 12 lipid labels sit in
    # column W of an 89-column sheet while a different column held 16. One
    # spreadsheet sheet routinely holds SEVERAL panels side by side; a
    # best-column scanner silently hides all but one of them.
    cols = []
    for j, vals in strs.items():
        uniq = {v for v in vals if LIPIDISH.match(v)}
        if len(uniq) >= 8:
            cols.append((j, len(uniq), sorted(uniq)[:6]))
    cols.sort(key=lambda c: -c[1])
    best = cols[0] if cols else None
    hdr = " | ".join(str(c)[:18] for c in (rows[0] or [])[:8] if c)
    timey = bool(TIMEISH.search(hdr)) or any(
        isinstance(c, str) and TIMEISH.search(c) for r in rows[:3] for c in (r or []) if c)
    return {"rows": len(rows), "header": hdr[:90], "timeseries_like": timey,
            "label_col": best, "label_cols": cols}

for p in sys.argv[1:]:
    try:
        wb = load_workbook(p, read_only=True, data_only=True)
    except Exception as e:
        print(f"  !! {os.path.basename(p)}: {e}"); continue
    for ws in wb.worksheets:
        c = classify(ws)
        if not c: continue
        if c["label_cols"]:
            for j, n, sample in c["label_cols"]:
                print(f"  SCREEN? {os.path.basename(p)[:38]:<38} "
                      f"{ws.title[:22]:<22} {c['rows']:>5}r  {n:>4} labels "
                      f"col{j}  {sample}")
        elif c["rows"] > 40 and not c["timeseries_like"]:
            print(f"  maybe   {os.path.basename(p)[:38]:<38} {ws.title[:22]:<22} "
                  f"{c['rows']:>5}r  hdr: {c['header'][:48]}")
    wb.close()
