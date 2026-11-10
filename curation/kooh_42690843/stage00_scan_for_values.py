#!/usr/bin/env python3
"""Make the intake gate's central claim CHECKABLE: the SI holds no value.

00_value_sources.json says every transfection measurement in this paper lives
in the main article, which we do not hold. "Not available" is not a sentence a
reviewer can check, so this script is the sentence: it scans all 49 SI pages
for the vocabulary a transfection readout is written in, and for the figure
captions that would carry one, and prints what it finds.

Run:  python3 -I stage00_scan_for_values.py

A NEW CHECK IS NOT DONE UNTIL IT HAS BEEN MADE TO FAIL ON PURPOSE
(CLAUDE.md section 6). Run with --seed to inject a synthetic page whose text
reads like a screening table; the scan must then report a hit. 07_qc.json
records that it did.
"""
import re
import sys
import os

import pymupdf

HERE = os.path.dirname(os.path.abspath(__file__))
SI = os.path.join(HERE, "source", "ja6c13674_si_001.pdf")

# Vocabulary a per-lipid transfection readout is written in.
VALUE_WORDS = [
    r"\bRLU\b", r"\bp/s\b", r"photons?/s", r"radiance",
    r"luminescen\w+", r"luciferase activity", r"relative expression",
    r"normali[sz]ed (?:to|expression|luminescence)", r"fold[- ]change",
    r"transfection efficien\w+ of", r"%\s*GFP", r"MFI",
    r"total flux", r"average radiance",
]
# A table of 275 members would have to name members. These are the library's
# own labels.
MEMBER_LABEL = re.compile(r"\bA(?:[1-9]|1\d|2[0-6])B[1-8]C(?:[1-9]|1[0-4])\b")


def scan(pages):
    hits = {"value_words": [], "member_labels": {}, "pages_with_both": []}
    for pno, text in pages:
        vw = []
        for pat in VALUE_WORDS:
            for m in re.finditer(pat, text, re.I):
                vw.append((pat, text[max(0, m.start() - 40):m.end() + 40]
                           .replace("\n", " ")))
        labels = sorted(set(MEMBER_LABEL.findall(text)))
        if vw:
            hits["value_words"].append((pno, vw))
        if labels:
            hits["member_labels"][pno] = labels
        if vw and labels:
            hits["pages_with_both"].append(pno)
    return hits


def main():
    doc = pymupdf.open(SI)
    pages = [(i + 1, p.get_text()) for i, p in enumerate(doc)]
    seeded = "--seed" in sys.argv
    if seeded:
        pages.append((999, "Table S99. Screening results.\n"
                           "A18B8C14 luminescence 1.24e7 RLU normalized to "
                           "untreated; A2B1C1 3.1e5 RLU; fold-change 40.1"))

    h = scan(pages)
    print("SI pages scanned: %d%s" % (len(pages),
                                      "  (+1 SEEDED page)" if seeded else ""))
    print("\n-- pages whose text uses transfection-readout vocabulary --")
    for pno, vw in h["value_words"]:
        print("  S%-3d %d hit(s)" % (pno, len(vw)))
        for pat, ctx in vw[:3]:
            print("        %-28s ...%s..." % (pat, ctx.strip()[:90]))
    if not h["value_words"]:
        print("  (none)")

    print("\n-- pages naming library members (AxByCz) --")
    for pno, labels in sorted(h["member_labels"].items()):
        print("  S%-3d %2d label(s): %s" % (pno, len(labels),
                                            ", ".join(labels[:8])))

    both = h["pages_with_both"]
    print("\n-- pages carrying BOTH a member label and a readout word --")
    print("  %s" % (both if both else "NONE"))
    n_members = len(set(x for v in h["member_labels"].values() for x in v))
    print("\nDistinct member labels anywhere in the SI: %d "
          "(the paper states 275 were made)." % n_members)

    # Every candidate page must be adjudicated HERE, by page number, with the
    # reason it is not a transfection value. An UNADJUDICATED candidate fails.
    # A bare "no hits" verdict would have been a lie: S29 is a real hit.
    ADJUDICATED = {
        29: ("Figure S7 caption. The readout is Cy5-mRNA CELLULAR UPTAKE by "
             "flow cytometry, normalised to untreated controls -- not a "
             "transfection readout -- and the caption prints no number; the "
             "values are bars in a raster panel of the SI figure, for ONE "
             "lipid (A18B8C14), with no per-library-member content."),
    }
    unadjudicated = [p for p in both if p not in ADJUDICATED]
    for p in both:
        if p in ADJUDICATED:
            print("  S%-3d ADJUDICATED: %s" % (p, ADJUDICATED[p]))
    print("\nVERDICT: %s" % (
        "NO page pairs a library member with a transfection readout. The SI "
        "holds no per-lipid transfection value."
        if not unadjudicated else
        "UNADJUDICATED candidate page(s) %s -- look at them." % unadjudicated))
    if seeded:
        if not unadjudicated:
            raise SystemExit(
                "SEEDED RUN DID NOT FIRE. The scan cannot see the thing it "
                "exists to see, so a clean run from it means nothing. "
                "CLAUDE.md section 6.")
        print("SEEDED RUN FIRED CORRECTLY on page(s) %s -- the scan has teeth."
              % unadjudicated)
        return 0
    if unadjudicated:
        raise SystemExit(1)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
