#!/usr/bin/env python3
"""Parse SI Table S1 (pp. 19-22) into 02_tableS1.csv.

Tier 4: the table is live text in the PDF. No pixels are read.
REFUSES rather than guessing: every row must match the strict pattern
  <Name> <Formula> <MW> <Z-Ave> <PDI> <EE%>
and the row count must be exactly 144.
"""
import csv, os, re, sys
import pymupdf

HERE = os.path.dirname(os.path.abspath(__file__))
SI = os.path.join(HERE, "source", "pnas.2409572122.sapp.pdf")
PAGES = (18, 19, 20, 21)          # 0-based: SI pp. 19-22

NAME = re.compile(r"^H[1-9]A[1-4]B[1-4]$")
FORM = re.compile(r"^C(\d+)H(\d+)N(\d+)O(\d+)$")


def tokens():
    doc = pymupdf.open(SI)
    out = []
    for p in PAGES:
        for ln in doc[p].get_text().splitlines():
            t = ln.strip()
            if t:
                out.append(t)
    return out


def main():
    tk = tokens()
    rows, i = [], 0
    while i < len(tk):
        if NAME.match(tk[i]):
            if i + 5 >= len(tk):
                raise SystemExit("Table S1: truncated row at %r" % tk[i])
            name, form, mw, z, pdi, ee = tk[i:i + 6]
            if not FORM.match(form):
                raise SystemExit("Table S1: %s formula %r unparseable" % (name, form))
            if not ee.endswith("%"):
                raise SystemExit("Table S1: %s EE %r not a percent" % (name, ee))
            rows.append(dict(Name=name, Formula=form, MW=float(mw),
                             Z_Ave_nm=float(z), PDI=float(pdi),
                             EE_pct=float(ee.rstrip("%"))))
            i += 6
        else:
            i += 1

    if len(rows) != 144:
        raise SystemExit("Table S1: parsed %d rows, expected 144" % len(rows))
    names = [r["Name"] for r in rows]
    if len(set(names)) != 144:
        raise SystemExit("Table S1: %d distinct names" % len(set(names)))
    expect = ["H%dA%dB%d" % (h, a, b)
              for h in range(1, 10) for a in range(1, 5) for b in range(1, 5)]
    if names != expect:
        bad = [(x, y) for x, y in zip(names, expect) if x != y]
        raise SystemExit("Table S1: name order departs from H1..H9 x A1..A4 x "
                         "B1..B4 at %r" % bad[:5])

    with open(os.path.join(HERE, "02_tableS1.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print("02_tableS1.csv: %d rows, names in exact H1..H9 x A1..A4 x B1..B4 order"
          % len(rows))
    print("MW range %.2f - %.2f ; Z-Ave %.1f - %.1f ; PDI %.3f - %.3f ; EE %.2f - %.2f"
          % (min(r["MW"] for r in rows), max(r["MW"] for r in rows),
             min(r["Z_Ave_nm"] for r in rows), max(r["Z_Ave_nm"] for r in rows),
             min(r["PDI"] for r in rows), max(r["PDI"] for r in rows),
             min(r["EE_pct"] for r in rows), max(r["EE_pct"] for r in rows)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
