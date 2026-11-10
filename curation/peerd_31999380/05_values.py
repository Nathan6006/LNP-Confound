#!/usr/bin/env python3
"""Every NUMBER this SI yields, with where it came from. Writes 05_values.csv.

RUN:  python3 -I 05_values.py

NONE OF THESE IS CORPUS-ADMISSIBLE. That is the paper's answer, not a shortfall
in the extraction, and the `corpus_admissible` / `reason_not_admissible`
columns carry it per row so nobody has to re-derive it from prose. The reason
is in 00_value_sources.json: data/'s Experiment_method vocabulary is delivery
only (luminescence_normalized, luminescence_discretized_normalized,
LRP6_knockdown_normalized, protein_adundance_normalized,
editing_efficiency_normalized, luminescence_relative_to_Spikevax) and nothing
below is any of those. The combinatorial-library delivery SCREEN is in main-text
Figures 1-4, which are not in source/.

TIERS. Tables S1 and S2 come out of the PDF text layer through
pdfplumber.find_tables() -- tier 4, exact, no pixels touched. Figures S2 and S4
are tier 5, and 00_value_sources.json states why tiers 1-4 fail for them.
"""
import csv
import json
import os
import re
import sys

import pdfplumber

HERE = os.path.dirname(os.path.abspath(__file__))
PDF = os.path.join(HERE, "source", "adma201906128-sup-0001-SuppMat.pdf")

# Shared metadata for the two in-vivo toxicology tables (captions, p.47/48).
VIVO = {"cargo": "siRNA", "cargo_detail": "unspecified siRNA",
        "model_system": "in_vivo", "cell_line_or_species": "mouse",
        "route": "intravenous", "dose": "1 mg/kg siRNA",
        "n_as_printed": "not stated in the caption"}

# Shared metadata for the two U266 viability figures (captions, p.38/40).
VITRO = {"cargo": "siRNA", "model_system": "in_vitro",
         "cell_line_or_species": "U266 (human multiple myeloma)",
         "organ": "", "route": "cell culture", "timepoint": "72 h",
         "assay": "XTT", "n_as_printed": "triplicate experiments"}


def split_mean_sd(cell):
    """'5.580 ±1.82' -> (5.580, 1.82). Returns the raw text on failure --
    never a guessed number."""
    if cell is None:
        return "", "", ""
    txt = " ".join(cell.split())
    m = re.match(r"^(-?[\d.]+)\s*±\s*([\d.]+)$", txt)
    if m:
        # [\d.]+ happily matches "0.2.3", which is Table S1's Lipid 10 / 2 h
        # Hematocrit SD. A character class is not a number test; float() is.
        # Without this the malformed cell was silently stored as an SD.
        try:
            float(m.group(1))
            float(m.group(2))
        except ValueError:
            return "", "", (f"UNPARSED: not a number -- {txt!r} "
                            "(printed this way in the SI)")
        return m.group(1), m.group(2), ""
    return "", "", f"UNPARSED as mean+-SD: {txt!r}"


def table_rows(page_index, table_label, caption_lipids, readout_kind, vid):
    rows = []
    with pdfplumber.open(PDF) as pdf:
        tbl = pdf.pages[page_index].find_tables()[0].extract()
    header = [" ".join((h or "").split()) for h in tbl[0]]
    for r in tbl[1:]:
        group = " ".join((r[0] or "").split())
        tp = "2 h" if group.endswith("2h") else (
            "24 h" if group.endswith("24h") else "n/a (untreated)")
        lipid = group.replace("2h", "").replace("24h", "").strip() or "UT"
        for col, cell in zip(header[1:], r[1:]):
            mean, sd, note = split_mean_sd(cell)
            vid += 1
            rows.append({
                "value_id": f"V{vid:04d}",
                "source": table_label, "page": page_index + 1, "panel": "",
                "lipid_or_group": lipid,
                "readout_as_printed": col.split(" ")[0],
                "units_as_printed": " ".join(col.split(" ")[1:]),
                "value": mean, "sd": sd, "value_raw_text":
                    " ".join((cell or "").split()),
                "timepoint": tp, "assay": readout_kind,
                "organ": "blood" if "S1" in table_label else "serum",
                "tier": "4_pdfplumber_text",
                "corpus_admissible": "no",
                "reason_not_admissible":
                    "haematology / clinical-chemistry toxicology; data/ has "
                    "no such Experiment_method",
                "extraction_note": note,
                **{k: v for k, v in VIVO.items()},
            })
    return rows, vid, caption_lipids


def figure_rows(vid):
    rows = []
    fp = os.path.join(HERE, "06_figure_values.json")
    if not os.path.exists(fp):
        return rows, vid
    with open(fp) as f:
        fig = json.load(f)
    for key, src, page, dose, cargo_detail in (
            ("S2", "Figure S2", 38, "0.06 uM",
             "NC-siRNA (negative control)"),
            ("S4", "Figure S4", 40, "0.06 uM",
             "siCon / siNC (negative control)")):
        for label, d in fig.get(key, {}).get("bars", {}).items():
            if "value_pct" not in d:
                continue
            vid += 1
            rows.append({
                "value_id": f"V{vid:04d}", "source": src, "page": page,
                "panel": ("whole panel" if key == "S2"
                          else "0.06 uM dose group"),
                "lipid_or_group": label,
                "readout_as_printed": "U266 cell viability",
                "units_as_printed": "%",
                "value": d["value_pct"], "sd": "",
                "value_raw_text": "",
                "timepoint": VITRO["timepoint"], "assay": VITRO["assay"],
                "dose": dose, "cargo_detail": cargo_detail,
                "tier": "5_pixels",
                "corpus_admissible": "no",
                "reason_not_admissible":
                    "cell viability under a NEGATIVE-CONTROL siRNA is a "
                    "safety readout, not delivery",
                "extraction_note":
                    f"bar top row {d['row']}, width {d['width_px']} px, "
                    f"flank delta {d['flank_delta_px']} px; "
                    "axis calibrated from 7 detected ticks",
                **{k: v for k, v in VITRO.items() if k != "timepoint"
                   and k != "assay"},
            })
    return rows, vid


def main():
    vid = 0
    all_rows = []
    r1, vid, _ = table_rows(46, "Table S1", "2, 6, 8 and 10",
                            "complete blood count", vid)
    all_rows += r1
    r2, vid, _ = table_rows(47, "Table S2", "2, 6, 8 and 10",
                            "serum chemistry", vid)
    all_rows += r2
    r3, vid = figure_rows(vid)
    all_rows += r3

    cols = ["value_id", "source", "page", "panel", "lipid_or_group",
            "readout_as_printed", "units_as_printed", "value", "sd",
            "value_raw_text", "cargo", "cargo_detail", "model_system",
            "cell_line_or_species", "organ", "route", "dose", "timepoint",
            "assay", "n_as_printed", "tier", "corpus_admissible",
            "reason_not_admissible", "extraction_note"]
    with open(os.path.join(HERE, "05_values.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in all_rows:
            w.writerow({c: r.get(c, "") for c in cols})

    unparsed = [r for r in all_rows if r.get("extraction_note", "").startswith(
        "UNPARSED")]
    print(f"{len(all_rows)} values -> 05_values.csv")
    for src in ("Table S1", "Table S2", "Figure S2", "Figure S4"):
        print(f"  {src}: {sum(1 for r in all_rows if r['source'] == src)}")
    print(f"  corpus-admissible: "
          f"{sum(1 for r in all_rows if r['corpus_admissible'] == 'yes')}")
    print(f"  cells that would not parse as mean+-SD: {len(unparsed)}")
    for r in unparsed:
        print("   ", r["value_id"], r["source"], r["lipid_or_group"],
              r["readout_as_printed"], "|", r["extraction_note"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
