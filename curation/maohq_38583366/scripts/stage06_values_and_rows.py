#!/usr/bin/env python3
"""05_values.csv (every value, with corpus_admissible and a reason),
06_rows.csv (header-only -- the zero has to be visible) and 07_qc.json.

The 79-column header for 06_rows.csv is read LIVE from data/data_vitro.csv, not
restated here.
"""
import csv, json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..", "..")
F = lambda n: os.path.join(HERE, "..", n)

COLS = ["value_id", "source_location", "quantity_as_printed", "units_as_printed",
        "tier", "digitised", "value", "ionizable_lipid", "formulation_or_pool",
        "pool_size", "cargo", "cell_or_organ", "route", "dose", "timepoint",
        "assay", "corpus_admissible", "reason"]


def main():
    si = list(csv.DictReader(open(F("01_si_tables.csv"))))
    f4b = json.load(open(F("01_fig4b.json")))
    f4d = json.load(open(F("02_fig4d.json")))
    rows = []

    for r in si:
        rows.append(dict(
            value_id="SI_" + r["lnp_id"], source_location=r["si_table"],
            quantity_as_printed="RLU in HepG2",
            units_as_printed="header says RLU; the values are ln(RLU)",
            tier="2_si_numeric_table", digitised="no", value=r["rlu_hepg2"],
            ionizable_lipid="DLin-MC3-DMA", formulation_or_pool=r["lnp_id"],
            pool_size=1, cargo="fLuc pDNA", cell_or_organ="HepG2",
            route="in_vitro", dose="1 ug/mL pDNA equivalent", timepoint="24 h",
            assay="luciferase, luminometer vs a standard curve",
            corpus_admissible="NO",
            reason="ALREADY MERGED: this value is already in data_vitro.csv as "
                   "YZ_2022 (PMID 35879315), which this paper cites as ref [16] "
                   "for this screen. 192/192 verified -- see "
                   "02_corpus_overlap.json."))

    for b in f4b["bars"]:
        rows.append(dict(
            value_id="F4B_" + b["label"], source_location="Fig 4B",
            quantity_as_printed="Total Bioluminescence Flux from Liver",
            units_as_printed="photon/s", tier="5_pixels", digitised="yes",
            value="%.4g" % b["value_photons_per_s"],
            ionizable_lipid="DLin-MC3-DMA", formulation_or_pool=b["label"],
            pool_size=1, cargo="fLuc pDNA", cell_or_organ="liver",
            route="intraduodenal", dose="100 ug pDNA per mouse",
            timepoint="not stated for this panel",
            assay="ex vivo IVIS bioluminescence, 5 min exposure",
            corpus_admissible="ADMISSIBLE IN KIND, NOT EMITTED",
            reason="A translated luciferase reporter, individually measured, "
                   "so the readout qualifies. Not emitted: the paper "
                   "contributes ONE ionizable lipid, already present in the "
                   "in-vivo pool from six sources and in liver from three of "
                   "them, so these 20 values are a FORMULATION rank of a "
                   "single benchmark molecule, not a chemistry rank. See "
                   "DECISIONS.md D-2."))

    for b in f4d["bars"]:
        rows.append(dict(
            value_id="F4D_%s_%s" % (b["formulation"], b["organ"]),
            source_location="Fig 4D",
            quantity_as_printed="Total Bioluminescence Flux",
            units_as_printed="photons/s", tier="5_pixels", digitised="yes",
            value="%.4g" % b["value_photons_per_s"],
            ionizable_lipid="DLin-MC3-DMA",
            formulation_or_pool=b["formulation"], pool_size=1,
            cargo="fLuc pDNA", cell_or_organ=b["organ"].lower(),
            route="intraduodenal", dose="not stated for this panel",
            timepoint="one week (per Sec 3.3 for CII-b3)",
            assay="ex vivo IVIS bioluminescence",
            corpus_admissible="ADMISSIBLE IN KIND, NOT EMITTED",
            reason="As Fig 4B. Read here to test the method against the six "
                   "liver-share percentages printed on the pie charts above "
                   "the panel."))

    # panels inventoried but not digitised -- one row per measurement slot, so
    # the counts sum and the absence is visible.
    slots = [
        ("Fig 3A", 24, 8, "intrahepatic", "liver",
         "POOLED, NOT ATTRIBUTABLE: the injected dose is a physical mixture of "
         "8 formulations with no barcode, so the flux belongs to no single "
         "LNP. This is NOT the LR_2023 barcoded precedent, which resolves "
         "members."),
        ("Fig 3B", 9, 8, "intravenous", "liver", "POOLED, NOT ATTRIBUTABLE."),
        ("Fig 3C", 10, 4, "intraduodenal", "liver", "POOLED, NOT ATTRIBUTABLE."),
        ("Fig 3D", 25, 4, "intraduodenal", "5 organs",
         "POOLED, NOT ATTRIBUTABLE."),
        ("Fig 4E", 12, 1, "intraduodenal", "liver",
         "ADMISSIBLE IN KIND, NOT EMITTED: a dose series (100/500/1000 ug) of "
         "four formulations already inventoried from Fig 4B; same single "
         "lipid."),
        ("Fig 4G", 12, 1, "intraduodenal and intravenous", "6 organs",
         "ADMISSIBLE IN KIND, NOT EMITTED: organ distribution of flux for "
         "CII-b3 only."),
        ("Fig 5B/5C", 2, 1, "intraduodenal", "liver",
         "ADMISSIBLE IN KIND, NOT EMITTED: Cas9 excision -> tdTomato is "
         "downstream of translation (tdTomato+ cells per unit area; 16.9% by "
         "FACS). One formulation."),
        ("Fig S7", 4, 1, "in_vitro", "HEK293T-Ai14",
         "ADMISSIBLE IN KIND, NOT EMITTED: % tdTomato+ at day 5, a "
         "recombination readout. The panel's formulation count is not stated "
         "in the caption."),
        ("Fig S9", 6, 1, "intraduodenal", "liver",
         "ADMISSIBLE IN KIND, NOT EMITTED: TIDE indel % on 6 liver samples, "
         "one formulation."),
        ("Fig S3", 6, 1, "intraduodenal", "6 organs",
         "NOT DELIVERY: Cy5-labelled pDNA is cargo/particle biodistribution, "
         "measured without any translation step."),
        ("Fig 6B/6C, S10, S11", 4, 1, "intraduodenal", "serum",
         "NOT DELIVERY: serum LDL-C, total cholesterol and ALT are downstream "
         "pharmacodynamic and toxicity endpoints."),
        ("Fig S1", 4, 1, "in_vitro", "HepG2",
         "NOT DELIVERY for the size/PDI part; the transfection part is the "
         "lyophilisation/SIF/SGF comparison of ONE formulation."),
        ("Fig S2, S4, S5, S6, S8", 5, 1, "various", "various",
         "NOT DELIVERY: survival, body weight, cytokine panel, confocal "
         "images."),
    ]
    for loc, n, pool, route, organ, reason in slots:
        for i in range(n):
            adm = ("NO" if reason.startswith(("POOLED", "NOT DELIVERY"))
                   else "ADMISSIBLE IN KIND, NOT EMITTED")
            rows.append(dict(
                value_id="%s_%02d" % (loc.replace(" ", "").replace(",", "")
                                      .replace("/", ""), i + 1),
                source_location=loc, quantity_as_printed="see reason",
                units_as_printed="see reason", tier="5_pixels",
                digitised="no", value="", ionizable_lipid="DLin-MC3-DMA",
                formulation_or_pool="pool" if pool > 1 else "single",
                pool_size=pool, cargo="pDNA", cell_or_organ=organ,
                route=route, dose="", timepoint="", assay="",
                corpus_admissible=adm, reason=reason))

    with open(F("05_values.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS)
        w.writeheader(); w.writerows(rows)

    # 06_rows.csv -- header-only, header read LIVE from the corpus
    with open(os.path.join(ROOT, "data", "data_vitro.csv")) as fh:
        header = fh.readline().rstrip("\n")
    with open(F("06_rows.csv"), "w") as fh:
        fh.write(header + "\n")

    counts = {}
    for r in rows:
        counts[r["corpus_admissible"]] = counts.get(r["corpus_admissible"], 0) + 1
    print(json.dumps({"values_inventoried": len(rows), "by_admissibility":
                      counts, "rows_emitted": 0,
                      "06_rows_columns": len(header.split(","))}, indent=1))
    return rows


if __name__ == "__main__":
    main()
