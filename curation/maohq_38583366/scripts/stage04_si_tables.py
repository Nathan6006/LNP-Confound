#!/usr/bin/env python3
"""Parse Supplementary Tables S1-S6 out of source/mmc1.docx -> 01_si_tables.csv.

READS THE DOCX IN THIS FOLDER, not an unpacked copy in a session scratchpad.
Scratchpads are session-scoped and get cleared, so a tracked output that depends
on one has no recoverable lineage -- that is F-023, and it is the second
instance of the class in this project. The archive is read with zipfile in
memory; nothing from it is written to disk and nothing in it is executed.

The table is read as a real w:tbl -> w:tr -> w:tc tree and every header cell is
asserted before a single data row is taken. Counting <row> tags gave 51 and 49
where reading the cells gave 25 on another paper this run; three tag-counting
methods, three wrong answers.
"""
import csv, os, re, sys, zipfile
from xml.etree import ElementTree as ET

W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"
HERE = os.path.dirname(os.path.abspath(__file__))
DOCX = os.environ.get("STAGE04_DOCX",
                      os.path.join(HERE, "..", "source", "mmc1.docx"))
OUT = os.environ.get("STAGE04_OUT",
                     os.path.join(HERE, "..", "01_si_tables.csv"))


def ptext(p):
    return "".join(t.text or "" for t in p.iter(W + "t"))


def main():
    with zipfile.ZipFile(DOCX) as zf:
        if zf.read("[Content_Types].xml")[:5] != b"<?xml":
            raise SystemExit("stage04: unexpected archive content")
        root = ET.fromstring(zf.read("word/document.xml"))
    body = root.find(W + "body")

    tables, pending = [], None
    for el in body:
        if el.tag == W + "p":
            s = ptext(el).strip()
            if re.match(r"^Supplementary Table \d+\.", s):
                pending = s
        elif el.tag == W + "tbl":
            rows = [[" ".join(ptext(p).strip() for p in tc.findall(W + "p")).strip()
                     for tc in tr.findall(W + "tc")]
                    for tr in el.findall(W + "tr")]
            tables.append((pending, rows))
            pending = None
    if len(tables) != 6:
        raise SystemExit("stage04: expected 6 tables, got %d" % len(tables))

    recs = []
    for cap, rows in tables:
        m = re.search(r"top 32 (\S+)-based", cap)
        if not m:
            raise SystemExit("stage04: caption does not name a helper: %r" % cap)
        helper = m.group(1)
        h0, h1 = rows[0], rows[1]
        if h0[0] != "LNP ID" or h0[2] != "RLU in HepG2":
            raise SystemExit("stage04: unexpected header %r" % (h0,))
        if h1[1:6] != [helper, "DLin-MC3", "Chol.", "DMG-PEG", "pDNA"]:
            raise SystemExit("stage04: unexpected subheader %r" % (h1,))
        data = rows[2:]
        if len(data) != 32:
            raise SystemExit("stage04: %s has %d data rows, not 32"
                             % (cap, len(data)))
        for r in data:
            if len(r) != 9:
                raise SystemExit("stage04: row of width %d: %r" % (len(r), r))
            recs.append(dict(si_table=cap.split(".")[0], helper=helper,
                             lnp_id=r[0], wt_helper=r[1], wt_mc3=r[2],
                             wt_chol=r[3], wt_peg=r[4], wt_pdna=r[5],
                             rlu_hepg2=r[6], z_average=r[7], pdi=r[8]))
    if len(recs) != 192 or len({r["lnp_id"] for r in recs}) != 192:
        raise SystemExit("stage04: %d rows, %d distinct ids"
                         % (len(recs), len({r["lnp_id"] for r in recs})))
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(recs[0].keys()))
        w.writeheader(); w.writerows(recs)
    print("stage04: %d rows, %d distinct LNP ids, 6 tables x 32"
          % (len(recs), len({r["lnp_id"] for r in recs})))


if __name__ == "__main__":
    main()
