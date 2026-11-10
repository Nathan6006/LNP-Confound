#!/usr/bin/env python3
"""Back the intake gate's claim: the SI PDF prints no per-peptoid delivery value.

Run from this folder:   python3 -I stage00_scan_for_values.py

00_value_sources.json asserts that every per-member expression value in this
paper lives in a main-article figure we do not hold, and that the SI's tables
carry particle properties only. That is a negative claim, and a negative claim
nobody can re-run is an impression. This makes it one command.

It does TWO passes, because either alone would be easy to fool:

  A. vocabulary -- pages whose text layer mentions luciferase / flux / photon /
     radiance / p/s / RLU / expression at all.
  B. magnitude  -- numeric tokens anywhere in the SI that fall in the
     1e6 .. 1e11 band a murine whole-organ photon flux occupies.

A page that trips BOTH is a candidate value source and is printed in full for a
human to adjudicate. The expected answer is that the only trips are figure
CAPTIONS (which describe the main-article panels) and axis TICK LABELS drawn as
text inside SI figures -- never a table cell keyed to a peptoid number.
"""
import os
import re
import sys

import pdfplumber

HERE = os.path.dirname(os.path.abspath(__file__))
SI = os.path.join(HERE, "source", "nn4c05513_si_001.pdf")

VOCAB = re.compile(
    r"luciferase|fluc|flux|photon|radiance|\bp/s\b|\bRLU\b|luminescen|expression",
    re.I)

# A photon-flux magnitude, in the two forms this SI writes numbers:
#   plain        1234567890
#   sci / drawn  4x10 9  /  4*10^9  /  1e9
PLAIN = re.compile(r"\b\d{7,11}\b")
SCI = re.compile(r"\b\d(?:\.\d+)?\s*[x*×]\s*10\s*\^?\s*(\d{1,2})\b")


def main():
    if not os.path.exists(SI):
        raise SystemExit(f"missing {SI}")
    both, vocab_only, mag_only = [], [], []
    with pdfplumber.open(SI) as pdf:
        pages = []
        for i, pg in enumerate(pdf.pages, start=1):
            t = pg.extract_text() or ""
            v = bool(VOCAB.search(t))
            mags = [m for m in PLAIN.findall(t)]
            mags += [f"1e{e}" for e in SCI.findall(t) if 6 <= int(e) <= 11]
            m = bool(mags)
            pages.append((i, t, v, mags))
            if v and m:
                both.append(i)
            elif v:
                vocab_only.append(i)
            elif m:
                mag_only.append(i)

    print(f"SI pages scanned: {len(pages)}")
    print(f"  vocabulary AND flux-magnitude : {both}")
    print(f"  vocabulary only               : {vocab_only}")
    print(f"  flux-magnitude only           : {mag_only}")
    print()
    for i, t, v, mags in pages:
        if i in both:
            print(f"--- page {i}: magnitudes {sorted(set(mags))}")
            for line in t.split("\n"):
                if VOCAB.search(line) or PLAIN.search(line) or SCI.search(line):
                    print("    " + line.strip()[:150])
    print()
    print("ADJUDICATION: a page is a VALUE SOURCE only if a flux magnitude sits")
    print("in a table cell keyed to a peptoid number. A magnitude that is an")
    print("axis tick label, or a caption describing a main-article panel, is")
    print("not. See 00_value_sources.json si_tables_inventory for the verdict")
    print("on every SI table.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
