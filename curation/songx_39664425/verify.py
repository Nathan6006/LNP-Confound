#!/usr/bin/env python3
"""Gate for curation/songx_39664425.

This curation emits ZERO corpus rows. The claim that needs a gate is therefore
not "the values are right" but "there are no values, and the article is a
Highlight commentary whose only quantitative figure is a licensed reprint of a
paper already in the corpus".

Every check below is re-derived from the archived source files in ./source/ and
reports an AGREEMENT COUNT. Seed a defect into a copy of a source file and the
agreement count must MOVE and the process must exit 1 (see 07_qc.json, which
records the seeding runs and their observed exit codes).

Usage:  python3 -I verify.py [source_dir]
"""
import hashlib
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

# --- facts asserted by this curation, transcribed BY HAND (C15 posture:
# --- the expected value is not parsed out of the prose it is checking).
EXPECT = {
    "pmid": "39664425",
    "pmcid": "PMC11628822",
    "article_type": "letter",
    "subject_heading": "Highlight",
    "n_table_wrap": 0,
    "n_html_table": 0,
    "n_supplementary_material": 0,
    "has_suppl": "N",
    "n_fig": 1,
    "n_graphic": 1,
    # Markup-level, not prose-level: the caption's "Ref." is followed by an
    # <xref> that POINTS AT bib1. That mechanically ties the reprinted figure to
    # reference 1 (PMID 38740955), which a prose string match would not.
    # (\xa0 is a real non-breaking space in the source; \s covers it.)
    "reprint_xref_pattern": r'Reprinted with the permission from Ref\.\s*<xref rid="bib1"',
    "primary_paper_pmid": "38740955",
    "figure_md5": "8cca629044a5e7ffee0f9d91f7c320c6",
    # A primary research article has at least one of these as a section title.
    "forbidden_section_titles": [
        "Methods", "Materials and methods", "Results", "Results and discussion",
        "Experimental", "Statistical analysis", "Data availability",
    ],
}


def main(argv):
    src = argv[1] if len(argv) > 1 else os.path.join(HERE, "source")
    xml = open(os.path.join(src, "europepmc_PMC11628822_fullTextXML.xml"),
               encoding="utf-8").read()
    html = open(os.path.join(src, "PMC11628822.html"),
                encoding="utf-8", errors="replace").read()
    core = json.load(open(os.path.join(src, "europepmc_core_39664425.json"),
                          encoding="utf-8"))["resultList"]["result"][0]
    figb = open(os.path.join(src, "gr1.jpg"), "rb").read()

    checks = []   # (name, ok, observed, expected)

    def chk(name, observed, expected):
        checks.append((name, observed == expected, observed, expected))

    # 1-2. identity, from the independent core record
    chk("pmid", str(core.get("pmid")), EXPECT["pmid"])
    chk("pmcid", str(core.get("pmcid")), EXPECT["pmcid"])

    # 3-4. article type and PMC subject heading -> this is a commentary
    m = re.search(r'<article\b[^>]*article-type="([^"]+)"', xml)
    chk("article_type", m.group(1) if m else None, EXPECT["article_type"])
    m = re.search(r'<subj-group subj-group-type="heading">\s*<subject>([^<]+)</subject>', xml)
    chk("subject_heading", m.group(1).strip() if m else None, EXPECT["subject_heading"])

    # 5-7. no tabular data anywhere, in EITHER representation
    chk("n_table_wrap", len(re.findall(r"<table-wrap\b", xml)), EXPECT["n_table_wrap"])
    chk("n_html_table", len(re.findall(r"<table\b", html, re.I)), EXPECT["n_html_table"])
    chk("n_supplementary_material",
        len(re.findall(r"<supplementary-material\b", xml)), EXPECT["n_supplementary_material"])

    # 8. publisher-side confirmation that no SI exists
    chk("has_suppl", str(core.get("hasSuppl")), EXPECT["has_suppl"])

    # 9-10. exactly one figure, one raster
    chk("n_fig", len(re.findall(r"<fig\b", xml)), EXPECT["n_fig"])
    chk("n_graphic", len(re.findall(r"<graphic\b", xml)), EXPECT["n_graphic"])

    # 11. the figure is a licensed reprint, not original data -- and the thing
    #     it is reprinted from is bib1, asserted through the xref target.
    figm = re.search(r"<fig\b.*?</fig>", xml, re.S)
    figxml = figm.group(0) if figm else ""
    chk("caption_reprint_xref_to_bib1",
        bool(re.search(EXPECT["reprint_xref_pattern"], figxml)), True)

    # 12. the paper it reprints from, by PMID, read out of reference bib1
    m = re.search(r'<ref id="bib1">.*?<pub-id pub-id-type="pmid">(\d+)</pub-id>', xml, re.S)
    chk("primary_paper_pmid", m.group(1) if m else None, EXPECT["primary_paper_pmid"])

    # 13. no primary-research section exists
    titles = {re.sub(r"<[^>]+>", "", t).strip()
              for t in re.findall(r"<title>(.*?)</title>", xml, re.S)}
    found = sorted(t for t in titles if t in EXPECT["forbidden_section_titles"])
    chk("primary_research_sections_present", found, [])

    # 14. the raster we reasoned over is the one archived
    chk("figure_md5", hashlib.md5(figb).hexdigest(), EXPECT["figure_md5"])

    # 15. nothing was emitted
    for f in ("03_enumerated.csv", "05_values.csv", "06_rows.csv"):
        p = os.path.join(HERE, f)
        n = 0
        if os.path.exists(p):
            n = max(0, sum(1 for _ in open(p, encoding="utf-8")) - 1)
        chk("emitted_rows:" + f, n, 0)

    n_ok = sum(1 for _, ok, _, _ in checks if ok)
    n = len(checks)
    for name, ok, obs, exp in checks:
        print("%-4s %-38s observed=%r expected=%r" % ("PASS" if ok else "FAIL",
                                                      name, obs, exp))
    print("\nAGREEMENT COUNT: %d/%d" % (n_ok, n))
    if n_ok != n:
        print("GATE: FAIL")
        return 1
    print("GATE: PASS  (article is a Highlight commentary; 0 rows emittable)")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
