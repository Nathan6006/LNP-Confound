#!/usr/bin/env python3
"""The measured values, from the tier 00_value_sources.json names -> 05_values.csv

Run from this folder:   python3 -I stage04_values.py

NO PIXEL IS READ. Three value sources, all text:

  antimicrobial_score   SI Table S3, page S10 of am3c05481_si_001.pdf, read by
                        pdfplumber word coordinates. 104 cells. Confirmed cell
                        by cell against a 600 dpi render of the same page.
  MIC_uM                main-text Table 1, from the PMC JATS <table> markup.
                        27 cells.
  TI and HC50_uM        main-text Table 2 (27 cells) and HC50 = TI x MIC, which
                        lands on the paper's own discrete hemolysis ladder for
                        27 of 27 compounds. See 07_qc.json P7.

WHAT THE ENDPOINT IS. Not transfection. Growth inhibition of Bacillus subtilis
in Luria Broth, read as optical density at 420-580 nm over 24 h, with "active"
defined as >99.9% inhibition of growth. There is no nucleic-acid cargo, no
helper lipid, no PEG-lipid, no nanoparticle formulation and no cell line: the
methylated lipidoid is dispersed in growth medium by itself. DECISIONS.md D-1
is the question of whether anything here belongs in a delivery corpus at all,
and it is OPEN, which is why this stage writes values and emits no corpus row.
"""
import csv
import json
import os
import re
import sys
from xml.etree import ElementTree as ET

import pdfplumber

HERE = os.path.dirname(os.path.abspath(__file__))
SI = os.path.join(HERE, "source", "am3c05481_si_001.pdf")
JATS = os.path.join(HERE, "source", "PMC10472336.jats.xml")
OUT = os.path.join(HERE, "05_values.csv")

# Table S3's grid, page S10. The eight headgroup columns sit at these x0
# positions in PDF points; the row label sits left of x=190. These are the
# page's OWN word coordinates, not a guess: every one of the 104 digits lands
# within 1 pt of one of them, which is asserted below.
S3_PAGE = 9
S3_COL_X = [194, 219, 245, 272, 298, 329, 363, 395]
S3_HEADS = ["2N", "3N", "3N'", "4N", "4N'", "3NB", "3NB'", "4NC"]
S3_TAILS = ["6B", "8", "8CB", "9", "9C", "10", "10B", "10B'", "11", "12",
            "13B", "14", "16"]
S3_TOTAL_X = 436          # the printed "Tail total" column
S3_Y_RANGE = (320, 520)   # the data band: below the header row, above "Head total"

# What the paper says about the assay, read from the Results and the
# Experimental Section. Recorded as READINGS of the text, with the sentence
# each came from, never as inferred facts.
ASSAY = {
    "organism": "Bacillus subtilis",
    "organism_kind": "Gram-positive bacterium, model strain (NOT a cell line)",
    "cargo": "none -- the lipidoid is tested as the compound itself; there is "
             "no nucleic acid, no helper lipid, no PEG-lipid and no "
             "nanoparticle in this assay",
    "formulation": "methylated (quaternised) lipidoid dispersed directly in "
                   "Luria Broth growth medium; DLS shows 100-1000 nm "
                   "aggregates for most members (Figure S2)",
    "route": "in vitro, liquid broth growth-inhibition assay (bioscreen OD)",
    "readout_as_printed": "optical density at 420-580 nm (OD420-580) of "
                          "planktonic bacteria",
    "activity_definition": ">99.9% inhibition of growth, evidenced either by "
                           "complete lack of growth over 24 h or by a "
                           "significantly delayed onset of the log phase; the "
                           "150 min delay threshold is derived in Figure S3",
    "inoculum_screen": "1e6 CFU/mL",
    "inoculum_PI_assay": "1e7 CFU/mL",
    "timepoint": "24 h",
    "screen_doses_uM": "25, 50, 100",
    "mic_doses_uM": "1 - 12.5 (the 27-member follow-up)",
    "hemolysis_doses_uM": "up to 200 (human erythrocytes)",
}


def norm(s):
    """Prime, curly quote and NBSP all mean the same apostrophe here."""
    return (s.replace("′", "'").replace("‘", "'")
             .replace("’", "'").replace(" ", " ").strip())


def read_table_s3():
    """The 104 printed scores, by word coordinate. Returns (grid, tail_totals).

    Refuses rather than guessing: a row that does not yield exactly 8 digits on
    the 8 known column positions raises, because a silently short row is how a
    grid read turns into a plausible wrong number.
    """
    with pdfplumber.open(SI) as pdf:
        words = pdf.pages[S3_PAGE].extract_words()
    rows = {}
    for w in words:
        if not (S3_Y_RANGE[0] <= w["top"] <= S3_Y_RANGE[1]):
            continue
        if not re.fullmatch(r"\d{1,2}", w["text"]):
            continue
        # Keep only the grid band. Left of 190 pt are the tail ROW LABELS
        # ("6", "8", "10", "13"...), which are digits too and would otherwise
        # be read as scores; right of 460 pt is the ACS download watermark,
        # which runs vertically down the margin and contains "10" and "2026".
        if not (190 <= w["x0"] <= 460):
            continue
        rows.setdefault(round(w["top"]), []).append(w)

    grid, totals, used = {}, {}, 0
    for i, key in enumerate(sorted(rows)):
        line = sorted(rows[key], key=lambda w: w["x0"])
        cells, total = {}, None
        for w in line:
            for j, x in enumerate(S3_COL_X):
                if abs(w["x0"] - x) <= 1.5:
                    cells[S3_HEADS[j]] = int(w["text"])
            if abs(w["x0"] - S3_TOTAL_X) <= 1.5:
                total = int(w["text"])
        if len(cells) != 8:
            raise SystemExit(
                f"Table S3 row at y={key}: found {len(cells)} of 8 score "
                f"cells on the known column positions {S3_COL_X}. Refusing to "
                "emit a partial row -- a short row read as complete is how a "
                "grid becomes a plausible wrong number.")
        tail = S3_TAILS[used]
        grid[tail] = cells
        totals[tail] = total
        used += 1
    if used != 13:
        raise SystemExit(f"Table S3: read {used} tail rows, expected 13.")
    return grid, totals


def read_jats_table(label):
    root = ET.parse(JATS).getroot()
    for tw in root.iter("table-wrap"):
        lb = tw.find("label")
        if lb is None or "".join(lb.itertext()).strip() != label:
            continue
        rows = [[norm("".join(c.itertext())) for c in tr]
                for tr in tw.iter("tr")]
        heads = rows[1][1:]
        out = {}
        for r in rows[2:]:
            out[norm(r[0])] = dict(zip(heads, r[1:]))
        return heads, out
    raise SystemExit(f"{label} not found in the JATS record.")


def main():
    grid, tail_totals = read_table_s3()
    mic_heads, mic = read_jats_table("Table 1")
    ti_heads, ti = read_jats_table("Table 2")
    if mic_heads != ti_heads:
        raise SystemExit(f"Table 1 heads {mic_heads} != Table 2 {ti_heads}")

    enum = {r["IL_name"]: r for r in
            csv.DictReader(open(os.path.join(HERE, "03_enumerated.csv")))}

    rows = []
    for tail in S3_TAILS:
        for head in S3_HEADS:
            name = f"{head}{tail}"
            if name not in enum:
                raise SystemExit(f"{name} is in Table S3 but not enumerated.")
            score = grid[tail][head]
            m = mic.get(tail, {}).get(head)
            t = ti.get(tail, {}).get(head)
            hc50 = ""
            if m and t:
                bound = t.startswith("≥")
                hc50 = round(float(t.lstrip("≥")) * float(m), 4)
                hc50 = f"{'>200' if bound else hc50}"
            rows.append({
                "IL_name": name, "headgroup": head, "tail": tail,
                "antimicrobial_score": score,
                "score_source": "SI Table S3 (page S10), printed integer",
                "score_meaning_as_printed": {
                    0: "No inhibition at any studied concentration",
                    1: "Delayed growth at 100 uM only",
                    2: "Delayed growth > 150 mins at 100 uM only",
                    3: "No growth at 100 uM only",
                    4: "Delayed growth at 50 uM",
                    5: "No growth at 50 uM",
                    6: "Delayed or no growth at 25 uM",
                }[score],
                "MIC_uM": m or "",
                "MIC_source": "main-text Table 1 (JATS <table>)" if m else
                              "not measured -- this member is outside the "
                              "27-compound follow-up set",
                "therapeutic_index_HC50_over_MIC": t or "",
                "TI_source": "main-text Table 2 (JATS <table>)" if t else "",
                "HC50_uM_derived": hc50,
                "HC50_derivation": ("TI x MIC; a bounded TI means HC50 was "
                                    "above the highest tested 200 uM")
                                   if hc50 else "",
                "readout_as_printed": ASSAY["readout_as_printed"],
                "units": "score: dimensionless 0-6 ordinal. MIC, HC50: "
                         "micromolar (uM). TI: dimensionless ratio.",
                "cargo": ASSAY["cargo"],
                "cell_or_organ": ASSAY["organism"]
                + " (" + ASSAY["organism_kind"] + ")",
                "route": ASSAY["route"],
                "dose": "screen 25/50/100 uM; MIC follow-up 1-12.5 uM",
                "timepoint": ASSAY["timepoint"],
                "in_27_member_subset": int(bool(m)),
                "parent_SMILES": enum[name]["parent_SMILES"],
                "methylated_cation_SMILES":
                    enum[name]["methylated_cation_SMILES"],
            })

    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    side = os.path.join(HERE, "05_values_context.json")
    with open(side, "w") as fh:
        json.dump({"assay": ASSAY,
                   "printed_tail_totals_tableS3": tail_totals,
                   "note": "printed_tail_totals_tableS3 is recorded because it "
                           "CONTRADICTS the cells it totals; see "
                           "08_contradictions.csv C-1."}, fh, indent=1)

    print(f"wrote {len(rows)} value rows -> {OUT}")
    print(f"  with a score  {sum(1 for r in rows if r['antimicrobial_score'] != '')}")
    print(f"  with an MIC   {sum(r['in_27_member_subset'] for r in rows)}")
    print(f"  with an HC50  {sum(1 for r in rows if r['HC50_uM_derived'])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
