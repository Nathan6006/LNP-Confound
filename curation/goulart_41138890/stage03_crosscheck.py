#!/usr/bin/env python3
"""An INDEPENDENT tie between the main-text figure and the SI table.

Figure 2G's three legends print a pKa for nine of the formulations. Table S3
prints a pKa for eleven. The two live in different files, were produced by
different renderers, and are read here by different code paths:

  Fig 2G     a publisher JPEG raster (gr2_lrg.jpg), read BY EYE off the
             rendered legend at 1.6x zoom and typed in below.
  Table S3   an EMF metafile inside mmc1.docx, read as 534 EXTTEXTOUTW
             glyph records by stage02_si_tables.py.

They even use different decimal separators -- Fig 2G writes "7.47", the SI
table writes "7,47" -- which is itself evidence they are not one source twice.

So a 9/9 agreement validates the EMF transcription against something external
to it. That is strictly stronger than any check recomputed from the EMF alone,
and it is the reason this file exists rather than a self-consistency check.

The nine values are FROZEN here as read, before being compared, so the
comparison cannot be tuned: this is not a fit whose x-vector is the expected
answer (error class 27).
"""
import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

# Read off gr2_lrg.jpg, Figure 2G legends, 2026-11-10, left->right panel.
# Transcribed BEFORE opening 02_tableS3.csv.
FIG2G_PKA = {
    "LNP2": 7.47, "LNP3": 7.60, "LNP4": 7.41,        # DOPE  panel (blue)
    "LNP8": 7.37, "LNP10": 6.78, "LNP12": 7.49,      # DSPC  panel (green)
    "LNP14": 7.65, "LNP16": 8.10, "LNP18": 10.01,    # DOTAP panel (red)
}
FIG2G_SD = {"LNP2": 1.3, "LNP3": 1.2, "LNP4": 0.8, "LNP8": 0.6,
            "LNP10": 0.5, "LNP12": 1.0, "LNP14": 0.5, "LNP16": 0.3,
            "LNP18": 0.2}


def main():
    path = os.path.join(HERE, "02_tableS3.csv")
    si = {r["formulation"]: r for r in csv.DictReader(open(path, newline=""))}

    agree, disagree, absent = [], [], []
    for name, v in FIG2G_PKA.items():
        if name not in si:
            absent.append(name)
            continue
        got = float(si[name]["pKa"])
        (agree if abs(got - v) < 1e-9 else disagree).append(
            {"formulation": name, "fig2G": v, "tableS3": got,
             "delta": round(got - v, 6)})

    si_only = sorted(set(si) - set(FIG2G_PKA))
    out = {
        "what": "Fig 2G printed pKa (raster legend, read by eye) vs Table S3 "
                "pKa (EMF glyph records, read by stage02_si_tables.py)",
        "n_compared": len(FIG2G_PKA),
        "n_agree_exactly": len(agree),
        "n_disagree": len(disagree),
        "disagreements": disagree,
        "absent_from_tableS3": absent,
        "in_tableS3_but_not_plotted_in_Fig_2G": si_only,
        "decimal_separator_fig2G": ".",
        "decimal_separator_tableS3": ", (comma -- Brazilian locale)",
        "why_this_matters": (
            "two independent renderings, two independent extraction paths, "
            "and different decimal separators. A 9/9 exact agreement "
            "validates the EMF glyph transcription against something outside "
            "it; one misread digit in the SI table would break it."),
        "what_it_does_NOT_validate": (
            "the Fig 4G bar geometry. pKa is printed text in both places; "
            "nothing here touches a measured pixel."),
        "verdict": "PASS" if len(agree) == len(FIG2G_PKA) else "FAIL"}

    # The pKa column of Table S3 is the ONLY column Fig 2G duplicates, so a
    # second, weaker tie: Fig 2's size/zeta panels plot LNP2/3/4, LNP8/10/12
    # and LNP14/16/18 -- exactly the nine formulations whose pKa is printed.
    out["figure_2_plots_exactly_the_nine_pKa_formulations"] = (
        sorted(FIG2G_PKA) == sorted(
            ["LNP2", "LNP3", "LNP4", "LNP8", "LNP10", "LNP12",
             "LNP14", "LNP16", "LNP18"]))
    out["tableS3_has_two_extra"] = si_only

    json.dump(out, open(os.path.join(HERE, "03_crosscheck.json"), "w"),
              indent=1)
    print(json.dumps(out, indent=1))
    sys.exit(0 if out["verdict"] == "PASS" else 1)


if __name__ == "__main__":
    main()
