#!/usr/bin/env python3
"""Stage 01 -- parse the SI's own TEXT LAYER for every FIL record.

TIER 4 (`4_pdf_text_table`), not tier 5. `ja5c13982_si_001.pdf` carries a full
text layer over sections 3.1-3.2: every FIL has a numbered synthesis heading, a
one-sentence statement of WHICH bromide/amine it was made from, an
`MS (ESI, m/z) of [M + H]+ calculated for <formula>, <m/z>; Found <m/z>` line,
and a 1H NMR integral list. None of that is read from pixels.

What this file does NOT do is decide anything. It records the SI's own strings
verbatim, including the ones that contradict each other, so that stage 03 can
reconcile them against the formula and the series pattern and stage 08 can
record the disagreement. A parser that silently normalises its input has given
up the ability to detect non-canonical input (error class 37).

    python3 -I stage01_parse_si.py            # writes 01_si_records.json
"""
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SI = os.path.join(HERE, "source", "ja5c13982_si_001.pdf")

# The ACS watermark is injected into every page's text layer as one word per
# line. Cutting at it is safe because it is always the page's last block.
WATERMARK = "Downloaded"

HEAD = re.compile(r"^(3\.2(?:\.\d+)+)\s+Synthesis of\s+(\S+)\s*$", re.M)
# Formula + the calculated m/z values that follow it. The SI sometimes prints
# two (the [M+H]+ and the [M+Na]+), sometimes one, and once "846.52.62".
# Whitespace-tolerant on purpose: the SI's hard wraps split "calculated for"
# and split formulas mid-token ("C44H84F3\nNO7"). An intolerant pattern silently
# dropped 3 of 42 records on the first run, which is error class 28 wearing a
# parser's clothes -- a miss that looks like an absence.
MS = re.compile(
    r"calculated\s*for\s+(C\s*\d+\s*H\s*\d+(?:\s*F\s*\d*)?\s*N\s*\d*"
    r"\s*O\s*\d*)\s*,\s*([0-9.,\s]+?)\s*(?:;|\.\s*1H|Found)", re.S)
FOUND = re.compile(r"Found\s*,?\s*([0-9.,\s]+?)\s*(?:\.\s*1H|1H|$)", re.S)
MATERIAL = re.compile(
    r"(?:The|Specifically,)?\s*([^.]*?)\s+(?:was|were) used(?:\s+as\s+"
    r"(?:the\s+)?(?:starting\s+)?materials?)?", re.S)


def page_text(doc, i):
    t = doc.pages[i].extract_text() or ""
    return t.split(WATERMARK)[0]


def main():
    import pdfplumber
    doc = pdfplumber.open(SI)
    # Sections 3.1-3.2 live on SI PDF pages 4-38 (1-based). Take the whole
    # range rather than a hand-picked set: a heading that moves must still be
    # found, and an extra page costs nothing.
    pages = {i + 1: page_text(doc, i) for i in range(3, 38)}
    body = "".join(pages[k] for k in sorted(pages))
    # Collapse the SI's hard line wraps, which split every formula and name.
    flat = re.sub(r"\n(?=[a-z0-9(,])", "", body)
    flat = re.sub(r"[ \t]+", " ", flat)

    heads = list(HEAD.finditer(body))
    recs = []
    for n, h in enumerate(heads):
        sec, name = h.group(1), h.group(2)
        start = h.end()
        end = heads[n + 1].start() if n + 1 < len(heads) else len(body)
        blk = body[start:end]
        blk_flat = re.sub(r"\s+", " ", blk)
        ms = MS.search(blk_flat)
        fo = FOUND.search(blk_flat)
        mat = MATERIAL.search(blk_flat)
        # which SI page the heading sits on
        off, pg = 0, None
        for k in sorted(pages):
            off2 = off + len(pages[k])
            if off <= h.start() < off2:
                pg = k
                break
            off = off2
        recs.append({
            "section": sec,
            "name_in_heading": name,
            "si_pdf_page": pg,
            "formula_stated": (re.sub(r"\s+", "", ms.group(1))
                               if ms else None),
            "mz_calc_stated": (ms.group(2).strip() if ms else None),
            "mz_found_stated": (fo.group(1).strip() if fo else None),
            "material_sentence": (mat.group(1).strip() if mat else None),
            "block": blk_flat.strip(),
        })

    out = os.path.join(HERE, "01_si_records.json")
    with open(out, "w") as fh:
        json.dump({
            "_what": "Verbatim SI section-3.2 records, one per numbered "
                     "synthesis heading. Nothing here is reconciled; the "
                     "contradictions are preserved on purpose.",
            "source_pdf": os.path.basename(SI),
            "n_records": len(recs),
            "records": recs,
        }, fh, indent=1)
    print(f"{len(recs)} records -> {out}")
    missing = [r["section"] for r in recs if not r["formula_stated"]]
    if missing:
        print("NO FORMULA PARSED for:", missing)
    return 0


if __name__ == "__main__":
    sys.exit(main())
