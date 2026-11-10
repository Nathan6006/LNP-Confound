#!/usr/bin/env python3
"""Flatten the SI .docx to text, with every embedded image named in place.

Run from this folder:  python3 -I stage01_si_text.py  -> 00_pages/si_text.txt
                                                      -> 00_pages/si_media.json

WHY THIS STAGE EXISTS AT ALL, and why its output is tracked in the folder
rather than left in a scratchpad: CLAUDE.md section 7 -- no script in this repo
may READ an input from a session scratchpad, because scratchpads are
session-scoped and get cleared, and a tracked artifact that depends on one has
no recoverable lineage. That has already happened twice (lit_review/build_data.py,
and build_manifest.py as F-023). stage03 parses the 39 printed ESI-MS values out
of 00_pages/si_text.txt, so this script is what makes that file regenerable from
source/ alone.

WHAT THE SI TURNED OUT TO BE. The briefing for this curation warned that a
<w:tbl> count of 0 can mean "the tables are EMF images", because a sibling
curation in this run found exactly that and recovered three tables from
EMR_EXTTEXTOUTW glyph records. That is NOT what is going on here, and this
script is how that was settled rather than assumed:

  * <w:tbl> element count = 0, AND
  * the string "Table S" occurs ZERO times in 959,060 characters of
    word/document.xml -- the SI has no numbered tables to have missed, and
  * all 83 files in word/media were typed: 48 PNG + 28 SVG (Schemes S1-S12 and
    Figures S1-S29: synthesis schemes, 1H NMR and ESI-MS spectra), 4 TIFF
    (Figures S30-S33), 2 EMF (Figures S34-S35) and the Figure S37 pair.
  * The two EMF metafiles DO carry a text layer and were opened. They are bar
    charts of cell viability and of endosomal escape. Neither is a table and
    neither is a delivery endpoint.

So every value in this paper is in a figure, and the structures are in prose
plus drawn schemes.
"""
import json
import os
import zipfile
from xml.etree import ElementTree as ET

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "source", "anie74215-sup-0001-suppmat.docx")
OUTDIR = os.path.join(HERE, "00_pages")

W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"
R = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"


def para_text(p, rels):
    out = []
    for n in p.iter():
        if n.tag == W + "t":
            out.append(n.text or "")
        elif n.tag == W + "tab":
            out.append("\t")
        elif n.tag == W + "br":
            out.append("\n")
        elif n.tag.endswith("}blip"):
            e = n.get(R + "embed")
            out.append("[[IMG %s]]" % rels.get(e, e))
        elif n.tag.endswith("}imagedata"):
            e = n.get(R + "id")
            out.append("[[IMGV %s]]" % rels.get(e, e))
        elif n.tag.endswith("}svgBlip"):
            e = n.get(R + "embed")
            out.append("[[SVG %s]]" % rels.get(e, e))
    return "".join(out)


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    z = zipfile.ZipFile(SRC)
    rels = {c.get("Id"): c.get("Target")
            for c in ET.fromstring(z.read("word/_rels/document.xml.rels"))}
    doc = z.read("word/document.xml").decode("utf-8", "replace")
    root = ET.fromstring(doc)
    body = root.find(W + "body")

    ntbl = doc.count("<w:tbl>")
    ntable_caption = doc.count("Table S")

    lines, i = [], 0
    for el in body:
        if el.tag == W + "p":
            i += 1
            t = para_text(el, rels)
            if t.strip():
                lines.append("%04d| %s" % (i, t))
        elif el.tag == W + "tbl":
            i += 1
            lines.append("%04d| <<<TABLE>>>" % i)
            for row in el.findall(W + "tr"):
                cells = [" ".join(para_text(p, rels) for p in c.findall(W + "p"))
                         for c in row.findall(W + "tc")]
                lines.append("     | " + " || ".join(cells))

    out = os.path.join(OUTDIR, "si_text.txt")
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")

    media = {}
    for n in z.namelist():
        if n.startswith("word/media/"):
            media[os.path.basename(n)] = {
                "bytes": z.getinfo(n).file_size,
                "ext": os.path.splitext(n)[1].lower().lstrip("."),
            }
    ext = {}
    for v in media.values():
        ext[v["ext"]] = ext.get(v["ext"], 0) + 1
    with open(os.path.join(OUTDIR, "si_media.json"), "w") as fh:
        json.dump({
            "_what": "Every embedded object in the SI, typed. This is the "
                     "check that the SI's tables are not images: there are no "
                     "tables at all.",
            "w_tbl_element_count": ntbl,
            "occurrences_of_the_string_Table_S": ntable_caption,
            "document_xml_chars": len(doc),
            "media_file_count": len(media),
            "media_by_extension": ext,
            "emf_files_opened_and_identified": {
                "image79.emf": "Figure S34, cell viability bar chart - not a "
                               "table, not a delivery endpoint",
                "image80.emf": "Figure S35, endosomal-escape quantification - "
                               "not a table, not a delivery endpoint",
            },
            "media": media,
        }, fh, indent=1)

    print(f"wrote {out} ({len(lines)} lines)")
    print(f"<w:tbl> elements: {ntbl}")
    print(f"occurrences of 'Table S' in document.xml: {ntable_caption}")
    print(f"document.xml characters: {len(doc)}")
    print(f"media files: {len(media)} -> {ext}")
    if ntbl == 0 and ntable_caption == 0:
        print("=> the SI has NO tables, images or otherwise. Every value in "
              "this paper is in a figure.")
    else:
        raise SystemExit(
            "stage01: this SI DOES contain table content. Re-read it before "
            "digitising anything -- a transcribed table beats a digitised "
            "panel on every axis.")


if __name__ == "__main__":
    main()
