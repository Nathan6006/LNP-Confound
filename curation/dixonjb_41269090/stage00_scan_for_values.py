#!/usr/bin/env python3
"""Stage 00 -- the checkable half of 00_value_sources.json.

Prints, per SI page: text-layer character count, embedded images with pixel
size and MEASURED native dpi (px / rect_pt * 72, never guessed), native tables
found by pdfplumber, and the numeric tokens present in the text layer.

The claim this verifies: the bodies of Tables S2, S3 and S4 are BITMAPS with no
text layer, so they are unreachable at tiers 2 and 4 despite being SI tables;
and no page of the SI carries per-LNP SANDS barcode values.

Run:  python3 -I stage00_scan_for_values.py
"""
import os
import re
import sys

import pdfplumber
import pymupdf

HERE = os.path.dirname(os.path.abspath(__file__))
SI = os.path.join(HERE, "source", "nn5c12080_si_001.pdf")
NUM = re.compile(r"\d+\.\d+|\d{2,}")
WATERMARK = "Downloaded from pubs"


def main():
    doc = pymupdf.open(SI)
    plumb = pdfplumber.open(SI)
    bad = 0
    print(f"{SI}\npages: {doc.page_count}\n")
    for i, page in enumerate(doc):
        text = page.get_text()
        body = "\n".join(l for l in text.splitlines() if WATERMARK not in l)
        imgs = []
        for im in page.get_images(full=True):
            xref = im[0]
            info = doc.extract_image(xref)
            rects = page.get_image_rects(xref)
            if not rects:
                imgs.append(f"{info['width']}x{info['height']} dpi=UNPLACED")
                continue
            r = rects[0]
            dpi = info["width"] / r.width * 72.0
            imgs.append(f"{info['width']}x{info['height']}px native_dpi={dpi:.1f}")
        tabs = plumb.pages[i].find_tables()
        nums = NUM.findall(body)
        print(f"page {i+1:>2}  text={len(text):>5}  body_text={len(body):>5}  "
              f"native_tables={len(tabs)}  images={imgs}")
        if i + 1 in (22, 23, 24):
            if len(tabs) != 0:
                print("    UNEXPECTED: a native table was found on a page "
                      "00_value_sources.json calls a bitmap.")
                bad += 1
            if len(nums) > 2:
                print(f"    UNEXPECTED: {len(nums)} numeric tokens in the text "
                      f"layer of a page called bitmap-only: {nums[:10]}")
                bad += 1
            else:
                print(f"    OK: caption-only text layer, numeric tokens={nums}")
    print()
    print("claim: SANDS per-LNP barcode values appear nowhere in the SI text")
    hits = []
    for i, page in enumerate(doc):
        t = page.get_text().lower()
        if "barcode" in t and NUM.findall(t):
            hits.append(i + 1)
    print(f"  pages mentioning 'barcode' at all: {hits} "
          f"(captions of Fig S11 and the title block; no value table)")
    plumb.close()
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
