#!/usr/bin/env python3
"""Stage 05 -- the catalogue of every value-bearing panel. Writes 05_values.csv.

    python3 -I 05_values.py

Hand-assembled from the figure captions, the Results text and the SI, each row
citing where it lives. The arithmetic (the counts, and whether they sum) is done
here rather than asserted, so the totals in 07_qc.json cannot drift from the
table.

Two separate judgements are kept in separate columns on purpose:
  readout_is_delivery   does the QUANTITY map to data/'s Experiment_method
                        vocabulary (luminescence / protein-abundance /
                        editing-efficiency / knockdown)?
  corpus_admissible     could this value become a corpus ROW?
A value can pass the first and fail the second, and most of this paper's do.
"""
import csv, json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))

COLS = ["panel", "source_location", "quantity_as_printed", "units_as_printed",
        "readout_class", "cargo", "biological_system", "organ_or_cell", "route",
        "dose", "timepoint", "assay", "n_values", "n_distinct_formulations",
        "attributable", "pooled", "pool_identity", "pool_size", "compositional",
        "readout_is_delivery", "corpus_admissible", "reason"]

MOUSE = "C57BL/6J mouse (SpCas9-expressing B6J, Jax 026179, for editing)"
IL_BLOCK = ("BLOCKED UPSTREAM OF THE VALUE: the only ionizable lipid is 7C1, a "
            "polydisperse branched-PEI/C15-epoxide polymer with no single correct "
            "SMILES, and data/ has 0 rows with a blank IL_SMILES. See DECISIONS D-1.")

R = [
    dict(panel="Figure 1B + SI Fig S1A (p. S3)",
         source_location="SI p. S3, three native text tables; Fig 1B table",
         quantity_as_printed="dose required to silence a target gene in a vascular bed",
         units_as_printed="'65% @ 1.25 mg / kg' etc.; mRNA or protein as marked",
         readout_class="literature review of EIGHT OTHER PAPERS",
         cargo="siRNA", biological_system="various, as published elsewhere",
         organ_or_cell="lung/heart/kidney/liver/spleen/pancreas/lymph node/brain",
         route="systemic", dose="as published elsewhere", timepoint="various",
         assay="as published elsewhere", n_values=24, n_distinct_formulations=0,
         attributable="yes, to the cited paper", pooled="no", pool_identity="",
         pool_size="", compositional="no", readout_is_delivery="n/a",
         corpus_admissible="no",
         reason="NOT THIS PAPER'S DATA. Eight other laboratories' published "
                "numbers, transcribed into a review table. Emitting them under "
                "this Experiment_ID would fabricate a laboratory."),
    dict(panel="Figure 1C",
         source_location="Fig 1C, 2 bars + SEM, 'ns'",
         quantity_as_printed="ICAM-2 Protein Expression (MFI, Normalized)",
         units_as_printed="% of siLuc control, axis 0-125",
         readout_class="protein_abundance_knockdown", cargo="siRNA (siICAM-2/siLuc)",
         biological_system=MOUSE, organ_or_cell="bone marrow endothelial cells (CD31+CD45-)",
         route="intravenous", dose="1.0 mg/kg siRNA", timepoint="3 days",
         assay="flow cytometry MFI", n_values=2, n_distinct_formulations=1,
         attributable="yes - the 'original' 80:20 7C1", pooled="no",
         pool_identity="", pool_size="", compositional="no",
         readout_is_delivery="yes", corpus_admissible="no",
         reason=IL_BLOCK + " Also note the Results text claims lung and heart "
                "silencing 'as expected' and cites Fig 1C for it; Fig 1C shows "
                "ONE bar pair and the caption assigns it to BMECs, so the "
                "lung/heart comparison is shown nowhere."),
    dict(panel="Figure 2B, 2C",
         source_location="Fig 2C, 4 bars + SEM, ** and ****",
         quantity_as_printed="ICAM-2 Protein Silencing (MFI, Normalized)",
         units_as_printed="axis 0-60",
         readout_class="protein_abundance_knockdown", cargo="siRNA + DNA barcode",
         biological_system=MOUSE, organ_or_cell="lung endothelial cells",
         route="intravenous",
         dose="0.16 / 0.5 / 1.5 mg/kg total nucleic acid (+ PBS)",
         timepoint="72 h", assay="flow cytometry MFI", n_values=4,
         n_distinct_formulations=1, attributable="yes - the 'original' 80:20 7C1",
         pooled="no", pool_identity="", pool_size="", compositional="no",
         readout_is_delivery="yes", corpus_admissible="no",
         reason=IL_BLOCK + " One formulation at three doses; the corpus key is "
                "the lipid, so this is one lipid however many doses."),
    dict(panel="Figure 2D",
         source_location="Fig 2D, 5 bars",
         quantity_as_printed="ICAM-2 Low Cells (%)", units_as_printed="%, axis 0-40",
         readout_class="FACS gate frequency", cargo="siRNA + DNA barcode",
         biological_system=MOUSE, organ_or_cell="lung endothelial cells",
         route="intravenous", dose="0.16 / 0.5 / 1.5 mg/kg", timepoint="72 h",
         assay="flow cytometry", n_values=5, n_distinct_formulations=1,
         attributable="yes", pooled="no", pool_identity="", pool_size="",
         compositional="no", readout_is_delivery="no", corpus_admissible="no",
         reason="A gate frequency, not a potency: the fraction of cells below a "
                "threshold, which depends on where the gate sits. Not in data/'s "
                "Experiment_method vocabulary."),
    dict(panel="Figure 3E, 3F, 3G (and Fig S2I)",
         source_location="Fig 3E one-column scatter, 115 grey circles + 2 blue; "
                         "3F/3G the same values against diameter",
         quantity_as_printed="Normalized DNA Delivery (%)",
         units_as_printed="%, log axis 0.001 - 100",
         readout_class="pooled DNA barcode recovery from a functionally sorted "
                       "population",
         cargo="siICAM-2 + DNA barcode (siRNA:barcode mass ratio 10:1)",
         biological_system=MOUSE,
         organ_or_cell="ICAM-2-LOW bone marrow endothelial cells, FACS-sorted",
         route="intravenous", dose="1.5 mg/kg total nucleic acid, pooled",
         timepoint="3 days", assay="Illumina MiniSeq barcode counts",
         n_values=117, n_distinct_formulations=115,
         attributable="NO - 115 identical unlabelled circles; and which 115 of the "
                      "120 tabulated LNPs were pooled is never stated",
         pooled="yes", pool_identity="library 1, the stable subset", pool_size=115,
         compositional="YES - every value is a share of a total that sums to 100% "
                       "across the pool",
         readout_is_delivery="no", corpus_admissible="no",
         reason="Two independent blockers. (1) The quantity is DNA barcode "
                "abundance, not expression: the functional element is the FACS "
                "gate, not the measurement. The corpus's only barcoded "
                "experiment, LR_2023 (407 of 1,803 in-vivo rows), carries "
                "Cargo_type = peptide_barcode with a TRANSLATED readout "
                "(protein_adundance_normalized); pooling is not the objection. "
                "(2) Not one point is attributable. " + IL_BLOCK),
    dict(panel="Figure 4C",
         source_location="Fig 4C, paired before/after plot, ~15 line pairs, *",
         quantity_as_printed="Normalized DNA Delivery (%)",
         units_as_printed="%, axis 0.0 - 7.5",
         readout_class="pooled DNA barcode recovery",
         cargo="siICAM-2 + DNA barcode", biological_system=MOUSE,
         organ_or_cell="ICAM-2-LOW BMECs", route="intravenous",
         dose="1.5 mg/kg total nucleic acid, pooled", timepoint="3 days",
         assay="Illumina MiniSeq", n_values=30, n_distinct_formulations=30,
         attributable="NO - the pairing is by identical 7C1/chol/PEG ratio, but no "
                      "line carries a label",
         pooled="yes", pool_identity="library 1, C18PEG2000 vs C18PEG3000 pairs",
         pool_size=115, compositional="yes", readout_is_delivery="no",
         corpus_admissible="no",
         reason="Same two blockers as Fig 3E. The paired design is the paper's "
                "strongest internal argument for PEG structure, and it still "
                "yields no attributable value."),
    dict(panel="Figure 5C",
         source_location="Fig 5C, 2 bars + SEM",
         quantity_as_printed="ICAM-2 Protein Silencing (MFI, Normalized)",
         units_as_printed="axis 0-20",
         readout_class="protein_abundance_knockdown", cargo="siICAM-2 + barcode",
         biological_system=MOUSE, organ_or_cell="bone marrow endothelial cells",
         route="intravenous", dose="1.5 mg/kg total nucleic acid", timepoint="3 days",
         assay="flow cytometry MFI", n_values=2, n_distinct_formulations=0,
         attributable="NO - each value is the whole POOL (library 1, library 2), "
                      "not any LNP",
         pooled="yes", pool_identity="library 1 (115) and library 2 (31)",
         pool_size="115 / 31", compositional="no", readout_is_delivery="yes",
         corpus_admissible="no",
         reason="A pool-level effect. There is no formulation to attach it to: "
                "the mice received 115 or 31 LNPs at once."),
    dict(panel="Figure 5D, 5E, 5F (and Fig S2J)",
         source_location="Fig 5D one-column scatter, 31 grey circles + 2 blue",
         quantity_as_printed="Normalized Barcode Count (%)",
         units_as_printed="%, log axis 0.1 - 100",
         readout_class="pooled DNA barcode recovery",
         cargo="siICAM-2 + DNA barcode", biological_system=MOUSE,
         organ_or_cell="ICAM-2-LOW BMECs", route="intravenous",
         dose="1.5 mg/kg total nucleic acid, pooled", timepoint="3 days",
         assay="Illumina MiniSeq", n_values=33, n_distinct_formulations=31,
         attributable="NO - unlabelled; which 31 of the 40 tabulated LNPs were "
                      "pooled is never stated",
         pooled="yes", pool_identity="library 2, the stable subset", pool_size=31,
         compositional="yes", readout_is_delivery="no", corpus_admissible="no",
         reason="Same two blockers as Fig 3E. Note the axis says 'Normalized "
                "Barcode Count (%)' while the caption says 'Normalized DNA "
                "delivery' -- see 08_contradictions.csv."),
    dict(panel="Figure 4A, 4B, 5H, 5I, 5J, 5K, 5L, Fig S2B",
         source_location="pie pairs + enrichment bars",
         quantity_as_printed="Enrichment in Top 10%",
         units_as_printed="% (pies), unitless ratio x100 (bars)",
         readout_class="derived summary over the pool", cargo="siICAM-2 + barcode",
         biological_system=MOUSE, organ_or_cell="ICAM-2-LOW BMECs",
         route="intravenous", dose="1.5 mg/kg pooled", timepoint="3 days",
         assay="computed from the barcode counts", n_values=46,
         n_distinct_formulations=0, attributable="no - a property of a BAND, not an LNP",
         pooled="yes", pool_identity="library 1 / library 2", pool_size="115 / 31",
         compositional="yes", readout_is_delivery="no", corpus_admissible="no",
         reason="A derived statistic over a band of formulations. USED HERE AS A "
                "CHECK INSTEAD: the 'Formulated' pies are a second rendering of "
                "the composition tables, and 5 of 6 reproduce them exactly "
                "(07_qc.json, check E1)."),
    dict(panel="Figure 6E, 6F",
         source_location="Fig 6E, 2 bars + SEM, **; Results text '37% protein "
                         "silencing', '4.8x increase in potency'",
         quantity_as_printed="ICAM-2 Protein Silencing (MFI, Normalized)",
         units_as_printed="%, axis 0-50",
         readout_class="protein_abundance_knockdown", cargo="siICAM-2 (vs siLuc)",
         biological_system=MOUSE, organ_or_cell="bone marrow endothelial cells",
         route="intravenous", dose="1.0 mg/kg siRNA", timepoint="3 days",
         assay="flow cytometry MFI", n_values=2, n_distinct_formulations=2,
         attributable="YES - 'original' 80:20 7C1 and BM1", pooled="no",
         pool_identity="", pool_size="", compositional="no",
         readout_is_delivery="yes", corpus_admissible="no", reason=IL_BLOCK),
    dict(panel="Figure 6G",
         source_location="Fig 6G, 4 bars + SEM",
         quantity_as_printed="ICAM-2 Protein Silencing (MFI, Normalized)",
         units_as_printed="%, axis 0-80",
         readout_class="protein_abundance_knockdown", cargo="siICAM-2",
         biological_system=MOUSE, organ_or_cell="lung and heart endothelial cells",
         route="intravenous", dose="1.0 mg/kg siRNA", timepoint="3 days",
         assay="flow cytometry MFI", n_values=4, n_distinct_formulations=2,
         attributable="YES - 'original' 7C1 and BM1, two organs each", pooled="no",
         pool_identity="", pool_size="", compositional="no",
         readout_is_delivery="yes", corpus_admissible="no",
         reason=IL_BLOCK + " Had rows been emittable, these two organs and the "
                "marrow of Fig 6E would belong in ONE experiment with Model_target "
                "carrying the difference and a single z (error class 22)."),
    dict(panel="Figure 6H",
         source_location="Fig 6H, 2 bars + SEM, *; Results text '15% indel'",
         quantity_as_printed="ICAM-2 Indel (%)", units_as_printed="%, axis 0-20",
         readout_class="editing_efficiency", cargo="sgICAM-2",
         biological_system="SpCas9-expressing B6J mouse (Jax 026179)",
         organ_or_cell="bone marrow endothelial cells", route="intravenous",
         dose="1.0 mg/kg sgRNA", timepoint="5 days", assay="TIDE decomposition",
         n_values=2, n_distinct_formulations=1, attributable="YES - BM1 (and PBS)",
         pooled="no", pool_identity="", pool_size="", compositional="no",
         readout_is_delivery="yes", corpus_admissible="no", reason=IL_BLOCK),
    dict(panel="Figure S2D, S2E",
         source_location="SI p. S10; Results text '16% reduction', '2.2-fold more potent'",
         quantity_as_printed="ICAM-2 Protein Expression (MFI, Normalized); Fold "
                             "ICAM-2 protein silencing (Normalized to 7C1)",
         units_as_printed="axis 0-125; fold, axis 0-3",
         readout_class="protein_abundance_knockdown", cargo="siICAM-2 (vs siLuc)",
         biological_system=MOUSE, organ_or_cell="bone marrow endothelial cells",
         route="intravenous", dose="1.0 mg/kg siRNA", timepoint="3 days",
         assay="flow cytometry MFI", n_values=4, n_distinct_formulations=2,
         attributable="YES - the screen-1 winner (60/37.5/2.5 C16PEG2000 = library 1 "
                      "LNP 25) and the original 7C1",
         pooled="no", pool_identity="", pool_size="", compositional="no",
         readout_is_delivery="yes", corpus_admissible="no", reason=IL_BLOCK),
    dict(panel="Figure S3A, S3B",
         source_location="SI p. S14; 5 bars + SEM, **",
         quantity_as_printed="Fold ICAM-2 protein silencing (Normalized to 7C1)",
         units_as_printed="fold, axis 0-6",
         readout_class="protein_abundance_knockdown", cargo="siICAM-2",
         biological_system=MOUSE, organ_or_cell="bone marrow endothelial cells",
         route="intravenous", dose="1.0 mg/kg siICAM-2", timepoint="3 days",
         assay="flow cytometry MFI", n_values=5, n_distinct_formulations=5,
         attributable="YES - original 7C1, screen-1 winner, LNP1, BM1, LNP3; S3A "
                      "prints each one's molar composition",
         pooled="no", pool_identity="", pool_size="", compositional="no",
         readout_is_delivery="yes", corpus_admissible="no",
         reason=IL_BLOCK + " THIS IS THE PANEL THAT WOULD HAVE BEEN EMITTED: five "
                "distinct formulations of one lipid, individually dosed, one "
                "organ, one timepoint."),
    dict(panel="Figure S3H, S3I, S3J",
         source_location="SI p. S15, QUANT ddPCR",
         quantity_as_printed="Barcode Counts (AU); Normalized Nucleic Acid Delivery "
                             "(fold increase compared to 7C1)",
         units_as_printed="AU, axis 0-1500; fold, axis 0-7",
         readout_class="biodistribution", cargo="DNA barcode only",
         biological_system=MOUSE,
         organ_or_cell="heart / lung / marrow ECs and CD34+ HSPCs",
         route="intravenous", dose="0.5 mg/kg", timepoint="4 h",
         assay="droplet digital PCR", n_values=12, n_distinct_formulations=2,
         attributable="yes - 7C1 and BM1", pooled="no", pool_identity="",
         pool_size="", compositional="no", readout_is_delivery="no",
         corpus_admissible="no",
         reason="BIODISTRIBUTION -- where the drug went, not whether it worked. "
                "The paper's own thesis is that this is a different quantity "
                "from functional delivery."),
    dict(panel="Figure S3L",
         source_location="SI p. S15, 2 bars",
         quantity_as_printed="ICAM-2 Indel (%)", units_as_printed="%, axis 0-20",
         readout_class="editing_efficiency", cargo="sgICAM-2",
         biological_system="SpCas9-expressing B6J mouse",
         organ_or_cell="CD34+ haematopoietic stem and progenitor cells",
         route="intravenous", dose="1.0 mg/kg sgRNA", timepoint="5 days",
         assay="TIDE", n_values=2, n_distinct_formulations=1,
         attributable="yes - BM1 (and PBS)", pooled="no", pool_identity="",
         pool_size="", compositional="no", readout_is_delivery="yes",
         corpus_admissible="no",
         reason=IL_BLOCK + " The text also calls the result 'undetectable', so "
                "the value is at or below the assay floor (value_is_bound = upper) "
                "even if it were emittable."),
    dict(panel="Figure 6C, 6D, Fig S3F, S3G",
         source_location="Fig 6C/6D bars; S3F/S3G histograms",
         quantity_as_printed="Alexa647 MFI; Normalized Alexa647 MFI",
         units_as_printed="MFI, axis 0-200000; % of vehicle, axis 0-125",
         readout_class="cellular uptake", cargo="AlexaFluor647-siGFP",
         biological_system="immortalised mouse aortic endothelial cells (iMAECs)",
         organ_or_cell="iMAECs in vitro", route="in vitro, 20 nM siRNA",
         dose="20 nM siRNA/well", timepoint="1 h", assay="flow cytometry",
         n_values=8, n_distinct_formulations=2, attributable="yes - 7C1 and BM1",
         pooled="no", pool_identity="", pool_size="", compositional="no",
         readout_is_delivery="no", corpus_admissible="no",
         reason="CELLULAR UPTAKE, explicitly excluded by the endpoint vocabulary: "
                "fluorescence associated with cells at 1 h includes surface-bound "
                "and endosomal particle."),
    dict(panel="Figure 3D, 3H, 5G, 6B, Fig S1H, S2G, S3C, S3D, S3E",
         source_location="scatter/bar panels and the Fig 6B table",
         quantity_as_printed="Diameter (nm); PDI; pKa; fluorescence intensity vs pH",
         units_as_printed="nm; unitless; pH units",
         readout_class="physicochemical characterisation", cargo="n/a",
         biological_system="n/a", organ_or_cell="n/a", route="n/a", dose="n/a",
         timepoint="n/a", assay="DLS; TNS assay", n_values=320,
         n_distinct_formulations=160, attributable="partly", pooled="no",
         pool_identity="", pool_size="", compositional="no",
         readout_is_delivery="no", corpus_admissible="no",
         reason="Size, PDI and pKa are not delivery endpoints. Note the diameters "
                "are also the gate that decided which LNPs were pooled, and the "
                "per-LNP values are unlabelled, so they cannot even be used to "
                "identify the 115 and 31."),
    dict(panel="Figure S2H, S2K, S2L, S3M, S3N",
         source_location="SI pp. S12, S15",
         quantity_as_printed="mouse age/sex/weight; % initial body weight",
         units_as_printed="g; %", readout_class="tolerability",
         cargo="n/a", biological_system=MOUSE, organ_or_cell="whole animal",
         route="intravenous", dose="1.0-1.5 mg/kg", timepoint="0-72 h",
         assay="weighing", n_values=40, n_distinct_formulations=0,
         attributable="no", pooled="mixed", pool_identity="", pool_size="",
         compositional="no", readout_is_delivery="no", corpus_admissible="no",
         reason="Body weight is a tolerability readout, excluded by the endpoint "
                "vocabulary."),
]


def main():
    with open(os.path.join(HERE, "05_values.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLS)
        w.writeheader()
        for r in R:
            w.writerow({c: r.get(c, "") for c in COLS})

    summary = {
        "n_panel_groups": len(R),
        "n_values_catalogued": sum(r["n_values"] for r in R),
        "by_corpus_admissible": {
            "yes": sum(r["n_values"] for r in R if r["corpus_admissible"] == "yes"),
            "no": sum(r["n_values"] for r in R if r["corpus_admissible"] == "no")},
        "by_readout_is_delivery": {
            "yes": sum(r["n_values"] for r in R if r["readout_is_delivery"] == "yes"),
            "no": sum(r["n_values"] for r in R if r["readout_is_delivery"] == "no"),
            "n/a": sum(r["n_values"] for r in R if r["readout_is_delivery"] == "n/a")},
        "delivery_values_that_are_ATTRIBUTABLE": sum(
            r["n_values"] for r in R
            if r["readout_is_delivery"] == "yes" and r["attributable"].startswith("YES")),
        "delivery_values_from_POOLED_unattributable_panels": sum(
            r["n_values"] for r in R
            if r["readout_is_delivery"] == "yes" and r["attributable"].startswith("NO")),
        "distinct_ionizable_lipids": 1,
        "distinct_formulations_formulated": 160,
        "distinct_formulations_pooled_and_measured": 146,
        "distinct_formulations_individually_measured": 5,
        "design_supported_per_LNP_values": 146,
        "attributed_per_LNP_values": 5,
        "counts_sum": True,
    }
    summary["counts_sum"] = (summary["by_corpus_admissible"]["yes"]
                             + summary["by_corpus_admissible"]["no"]
                             == summary["n_values_catalogued"])
    assert summary["counts_sum"]
    assert (sum(summary["by_readout_is_delivery"].values())
            == summary["n_values_catalogued"])
    with open(os.path.join(HERE, "05_values_summary.json"), "w") as f:
        json.dump(summary, f, indent=1)
    print(json.dumps(summary, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
