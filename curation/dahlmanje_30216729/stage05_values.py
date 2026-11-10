#!/usr/bin/env python3
"""Stage 05 -- the complete panel inventory, with corpus_admissible + reason on
EVERY row, and the counts made to sum.

    python3 -I stage05_values.py [OUTDIR]

This is the whole-paper sweep: every panel, table and Methods-only readout that
could conceivably carry a value, judged against the corpus endpoint vocabulary.
It is written as an inventory rather than as extracted values because the
verdict on every row is the same and it is reached before extraction: nothing
in this paper is measured downstream of translation.

Writes 05_values.csv and 06_rows.csv. 06_rows.csv is HEADER-ONLY, with the
79-column header read LIVE off data/data_vivo.csv so the zero is visible in the
corpus's own schema rather than in a remembered one.
"""
import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
OUTDIR = os.path.abspath(sys.argv[1]) if len(sys.argv) > 1 else HERE

# Shared experimental context for the in vivo screens, as printed.
IV = dict(route="intravenous (lateral tail vein)", dose="0.5 mg/kg total DNA",
          cargo="single-stranded DNA barcode, 91 nt, phosphorothioate-capped",
          cargo_type="DNA barcode", species="C57BL/6J (WT) and Cav1-/- (JAX 007083)")

NOT_TRANSLATED = ("readout is the abundance of the administered DNA barcode "
                  "itself (ddPCR and/or Illumina); nothing downstream of "
                  "translation is measured")

FIELDS = ["panel", "locus", "readout_as_printed", "units_as_printed",
          "assay", "cargo", "cargo_type", "cell_type", "organ", "route",
          "dose", "timepoint", "genotype", "pooled", "pool_identity",
          "pool_size", "n_values_design_supported", "n_values_attributable",
          "corpus_admissible", "admissible_reason"]


def R(**kw):
    row = {f: "" for f in FIELDS}
    row.update(kw)
    return row


ROWS = [
    # ---------------- Figure 1: assay characterisation, no biology ----------
    R(panel="Fig 1a-b", locus="main text Fig 1",
      readout_as_printed="barcode schematic", units_as_printed="none",
      assay="schematic", corpus_admissible="no",
      admissible_reason="schematic, carries no measurement"),
    R(panel="Fig 1c", locus="main text Fig 1",
      readout_as_printed="Standard curve of QUANT barcodes diluted in TE buffer",
      units_as_printed="aM (750 aM - 12 fM), R2 = 1.00",
      assay="ddPCR", cargo=IV["cargo"], cargo_type=IV["cargo_type"],
      cell_type="none (buffer)", corpus_admissible="no",
      admissible_reason="instrument calibration in buffer; no biological "
                        "system, so not a delivery measurement of any kind"),
    R(panel="Fig 1d", locus="main text Fig 1",
      readout_as_printed="barcodes identified above background at 300 aM",
      units_as_printed="aM", assay="ddPCR", cell_type="none (buffer)",
      corpus_admissible="no",
      admissible_reason="limit of detection; assay sensitivity, not delivery"),
    R(panel="Fig 1e", locus="main text Fig 1",
      readout_as_printed="An in vitro standard curve",
      units_as_printed="zeptogram/well (60-16,000 zg), R2 = 0.91",
      assay="ddPCR after Lipofectamine 2000 transfection",
      cargo=IV["cargo"], cargo_type=IV["cargo_type"],
      cell_type="iMAEC (immortalised mouse aortic endothelial)",
      timepoint="24 h", corpus_admissible="no",
      admissible_reason=NOT_TRANSLATED + "; also a transfection-reagent dose "
                        "response, not an LNP screen"),
    R(panel="Fig 1f", locus="main text Fig 1",
      readout_as_printed="QUANT barcode readouts immediately after DNA was "
                         "isolated, or after storage at -20C for 20 or 31 days",
      units_as_printed="ddPCR counts", assay="ddPCR",
      cell_type="lung endothelial cells", organ="lung", route=IV["route"],
      dose=IV["dose"], timepoint="24 h", corpus_admissible="no",
      admissible_reason="sample-stability / operator repeatability control; "
                        + NOT_TRANSLATED),
    R(panel="Fig 1g", locus="main text Fig 1",
      readout_as_printed="comparison of fluorescence and QUANT methods",
      units_as_printed="none", assay="schematic", corpus_admissible="no",
      admissible_reason="schematic, carries no measurement"),

    # ---------------- Figure 2: biodistribution, 13 cell types --------------
    R(panel="Fig 2a", locus="main text Fig 2", readout_as_printed="schematic",
      units_as_printed="none", assay="schematic", corpus_admissible="no",
      admissible_reason="schematic, carries no measurement"),
    R(panel="Fig 2b", locus="main text Fig 2",
      readout_as_printed="Cumulative biodistribution measured by QUANT or "
                         "fluorescence in liver and non-liver cell types",
      units_as_printed="% of total signal (87% fluorescence vs 56% QUANT in liver)",
      assay="ddPCR and Alexa-647 MFI by flow cytometry",
      cargo=IV["cargo"], cargo_type=IV["cargo_type"],
      cell_type="13 sorted cell types", organ="liver, heart, lung, spleen, kidney",
      route=IV["route"], dose=IV["dose"], timepoint="1 h",
      genotype="WT", n_values_design_supported="26",
      n_values_attributable="2", corpus_admissible="no",
      admissible_reason="BIODISTRIBUTION, named as such by the paper; on the "
                        "brief's explicit not-delivery list. " + NOT_TRANSLATED),
    R(panel="Fig 2c-d", locus="main text Fig 2",
      readout_as_printed="Cumulative biodistribution within the 5 tissues "
                         "examined by QUANT and fluorescence",
      units_as_printed="% of total signal",
      assay="ddPCR and Alexa-647 MFI", cell_type="13 sorted cell types",
      organ="liver, heart, lung, spleen, kidney", route=IV["route"],
      dose=IV["dose"], timepoint="1 h", genotype="WT",
      n_values_design_supported="10", n_values_attributable="0",
      corpus_admissible="no",
      admissible_reason="biodistribution as a share of total signal; "
                        + NOT_TRANSLATED),
    R(panel="Fig 2e", locus="main text Fig 2",
      readout_as_printed="Comparison of biodistribution in the 13 cell types",
      units_as_printed="normalised to Kupffer cells",
      assay="ddPCR and Alexa-647 MFI", cell_type="13 sorted cell types",
      organ="liver, heart, lung, spleen, kidney", route=IV["route"],
      dose=IV["dose"], timepoint="1 h", genotype="WT",
      n_values_design_supported="26", n_values_attributable="0",
      corpus_admissible="no",
      admissible_reason="biodistribution, re-expressed relative to one cell "
                        "type; " + NOT_TRANSLATED),

    # ---------------- Figure 3: pharmacokinetics ----------------------------
    R(panel="Fig 3a", locus="main text Fig 3", readout_as_printed="schematic",
      units_as_printed="none", assay="schematic", corpus_admissible="no",
      admissible_reason="schematic, carries no measurement"),
    R(panel="Fig 3b", locus="main text Fig 3",
      readout_as_printed="Relative nanoparticle biodistribution (normalized to "
                         "maximal signal in any cell type)",
      units_as_printed="relative, max-normalised",
      assay="ddPCR and Alexa-647 MFI",
      cell_type="liver ECs, Kupffer, hepatocytes, lung ECs, lung macrophages",
      organ="liver, lung", route=IV["route"], dose=IV["dose"],
      timepoint="0.4, 0.75, 1.25, 12, 24, 36 h", genotype="WT",
      n_values_design_supported="60", n_values_attributable="0",
      corpus_admissible="no",
      admissible_reason="biodistribution time course; " + NOT_TRANSLATED),
    R(panel="Fig 3c", locus="main text Fig 3",
      readout_as_printed="Comparisons of area under the curve",
      units_as_printed="AUC", assay="ddPCR and Alexa-647 MFI",
      cell_type="5 cell types", organ="liver, lung", route=IV["route"],
      dose=IV["dose"], timepoint="0.4-36 h integrated", genotype="WT",
      n_values_design_supported="10", n_values_attributable="0",
      corpus_admissible="no",
      admissible_reason="integral of a biodistribution curve; " + NOT_TRANSLATED),
    R(panel="Fig 3d", locus="main text Fig 3",
      readout_as_printed="Peak DNA delivery (normalized to liver ECs)",
      units_as_printed="relative to liver ECs", assay="ddPCR and Alexa-647 MFI",
      cell_type="5 cell types", organ="liver, lung", route=IV["route"],
      dose=IV["dose"], timepoint="peak over 0.4-36 h", genotype="WT",
      n_values_design_supported="10", n_values_attributable="0",
      corpus_admissible="no",
      admissible_reason="peak barcode abundance; the word 'delivery' in the "
                        "axis label does not change what was measured. "
                        + NOT_TRANSLATED),
    R(panel="Fig 3e", locus="main text Fig 3",
      readout_as_printed="R2 analysis of QUANT absolute counts across two "
                         "experiments",
      units_as_printed="R2 = 0.98", assay="ddPCR", cell_type="5 cell types",
      corpus_admissible="no",
      admissible_reason="between-experiment reproducibility statistic, not a "
                        "per-LNP measurement"),

    # ---------------- Figure 4: the screens set up --------------------------
    R(panel="Fig 4a-c", locus="main text Fig 4",
      readout_as_printed="LNP library 1 was synthesized with the amine 7C1, "
                         "Cholesterol, DSPC, and PEG compounds at variable "
                         "molar ratios; 128 different LNPs",
      units_as_printed="none (design)", assay="schematic / design",
      corpus_admissible="no",
      admissible_reason="design description, no measurement. Quoted here "
                        "because it is the paper's own statement that ONE "
                        "amine is used and only the molar ratios vary"),
    R(panel="Fig 4d", locus="main text Fig 4",
      readout_as_printed="The diameter of each LNP in screen 1 was measured "
                         "individually",
      units_as_printed="nm", assay="dynamic light scattering",
      n_values_design_supported="128", n_values_attributable="128",
      corpus_admissible="no",
      admissible_reason="hydrodynamic diameter is on the brief's explicit "
                        "not-delivery list (size/PDI). Fully attributable and "
                        "fully extracted anyway, into 01_blocks.csv, because "
                        "it supports the inclusion-criterion check"),
    R(panel="Fig 4e", locus="main text Fig 4",
      readout_as_printed="The average normalized delivery from all LNPs and "
                         "the naked barcode (negative control) from screen 1",
      units_as_printed="Normalized DNA Delivery (%)",
      assay="Illumina deep sequencing of recovered barcodes",
      cargo=IV["cargo"], cargo_type=IV["cargo_type"],
      cell_type="10 sorted cell types", organ="liver, lung, kidney, heart",
      route=IV["route"], dose=IV["dose"], timepoint="24 h",
      genotype="WT and Cav1-/-", pooled="yes",
      pool_identity="library 1 (7C1 / Chol / DSPC / PEG), co-injected",
      pool_size="111 LNPs + 1 naked-barcode control",
      n_values_design_supported="111", n_values_attributable="0",
      corpus_admissible="no",
      admissible_reason=NOT_TRANSLATED + "; AND compositional (a share of the "
                        "pool, sums to 100%); AND pool membership is not "
                        "published, so no value is attributable to a named LNP"),
    R(panel="Fig 4f", locus="main text Fig 4",
      readout_as_printed="the naked barcode was delivered less efficiently "
                         "than barcodes carried by LNPs in every cell type",
      units_as_printed="Normalized DNA Delivery (%)",
      assay="Illumina deep sequencing", cell_type="10 sorted cell types",
      route=IV["route"], dose=IV["dose"], timepoint="24 h",
      genotype="WT and Cav1-/-", pooled="yes",
      pool_identity="library 1", pool_size="111 LNPs + naked control",
      n_values_design_supported="20", n_values_attributable="20",
      corpus_admissible="no",
      admissible_reason="negative-control contrast (naked DNA vs LNP mean); "
                        + NOT_TRANSLATED),

    # ---------------- Figure 5: the Cav1 result -----------------------------
    R(panel="Fig 5a", locus="main text Fig 5",
      readout_as_printed="The total ddPCR counts in all tested cell types - "
                         "which are equal to the area of the circle",
      units_as_printed="ddPCR counts (encoded as pie area)", assay="ddPCR",
      cell_type="10 sorted cell types", organ="whole body, liver, heart, lung, kidney",
      route=IV["route"], dose=IV["dose"], timepoint="24 h",
      genotype="WT and Cav1-/-", pooled="yes", pool_identity="library 1",
      pool_size="111", n_values_design_supported="10",
      n_values_attributable="0", corpus_admissible="no",
      admissible_reason="total barcode counts per organ, i.e. biodistribution; "
                        + NOT_TRANSLATED),
    R(panel="Fig 5b-d", locus="main text Fig 5",
      readout_as_printed="The total ddPCR counts were determined in different "
                         "cell-types from the liver, lung and kidney",
      units_as_printed="ddPCR counts (liver -93%, lung -43%, kidney -27% in Cav1-/-)",
      assay="ddPCR", cell_type="ECs, macrophages/Kupffer, immune, hepatocytes",
      organ="liver, lung, kidney", route=IV["route"], dose=IV["dose"],
      timepoint="24 h", genotype="WT and Cav1-/-", pooled="yes",
      pool_identity="library 1", pool_size="111",
      n_values_design_supported="18", n_values_attributable="0",
      corpus_admissible="no",
      admissible_reason="organ- and cell-level barcode counts; the paper's "
                        "central biological result, and still not delivery in "
                        "the corpus sense. " + NOT_TRANSLATED),
    R(panel="Fig 5e", locus="main text Fig 5",
      readout_as_printed="normalized library 1 nanoparticle biodistribution "
                         "within the liver cell-types",
      units_as_printed="% Normalized Nanoparticle Uptake (to WT Kupffer Cells)",
      assay="Illumina deep sequencing", cell_type="Kupffer, liver ECs, hepatocytes",
      organ="liver", route=IV["route"], dose=IV["dose"], timepoint="24 h",
      genotype="WT and Cav1-/-", pooled="yes", pool_identity="library 1",
      pool_size="111", n_values_design_supported="6",
      n_values_attributable="0", corpus_admissible="no",
      admissible_reason="the axis label says UPTAKE; uptake is on the brief's "
                        "explicit not-delivery list. " + NOT_TRANSLATED),
    R(panel="Fig 5f", locus="main text Fig 5",
      readout_as_printed="Combined sequencing data and ddPCR results shows the "
                         "absolute delivery of 111 nanoparticles for each LNP "
                         "in the liver",
      units_as_printed="Absolute Delivery (AU), log10 axis 10 to 1,000,000",
      assay="Illumina normalized delivery x ddPCR counts",
      cargo=IV["cargo"], cargo_type=IV["cargo_type"],
      cell_type="Kupffer, liver ECs, hepatocytes", organ="liver",
      route=IV["route"], dose=IV["dose"], timepoint="24 h",
      genotype="WT and Cav1-/-", pooled="yes", pool_identity="library 1",
      pool_size="111", n_values_design_supported="666",
      n_values_attributable="0", corpus_admissible="no",
      admissible_reason="THE per-LNP panel, and the one a corpus would want. "
                        "Inadmissible on endpoint (" + NOT_TRANSLATED + "), "
                        "and separately unattributable: 666 points are drawn "
                        "with no per-point label and the identity of the 111 "
                        "pooled members is never published"),

    # ---------------- Figure 6: macrophages and host phenotype --------------
    R(panel="Fig 6a", locus="main text Fig 6",
      readout_as_printed="Nanoparticle biodistribution in macrophages isolated "
                         "from multiple tissues",
      units_as_printed="ddPCR counts", assay="ddPCR",
      cell_type="Kupffer, lung macrophages, kidney macrophages",
      organ="liver, lung, kidney", route=IV["route"], dose=IV["dose"],
      timepoint="24 h", genotype="WT and Cav1-/-", pooled="yes",
      pool_identity="library 1", pool_size="111",
      n_values_design_supported="6", n_values_attributable="0",
      corpus_admissible="no", admissible_reason="biodistribution; " + NOT_TRANSLATED),
    R(panel="Fig 6b", locus="main text Fig 6",
      readout_as_printed="Absolute nanoparticle delivery to wild-type and "
                         "Cav1-/- macrophages in the liver, lung, and kidney",
      units_as_printed="Absolute Delivery (AU)",
      assay="Illumina normalized delivery x ddPCR counts",
      cell_type="Kupffer, lung macrophages, kidney macrophages",
      organ="liver, lung, kidney", route=IV["route"], dose=IV["dose"],
      timepoint="24 h", genotype="WT and Cav1-/-", pooled="yes",
      pool_identity="library 1", pool_size="111",
      n_values_design_supported="666", n_values_attributable="0",
      corpus_admissible="no",
      admissible_reason=NOT_TRANSLATED + "; same attributability failure as Fig 5f"),
    R(panel="Fig 6c", locus="main text Fig 6",
      readout_as_printed="The percentage of Kupffer cells (CD68+ CD45+) within "
                         "the immune cell population (CD45+)",
      units_as_printed="% of CD45+", assay="flow cytometry",
      cell_type="Kupffer cells", organ="liver",
      dose="1 mg/kg siLuciferase", genotype="WT and Cav1-/-",
      corpus_admissible="no",
      admissible_reason="host cell frequency, a control for the Cav1 "
                        "interpretation; not a delivery readout"),
    R(panel="Fig 6d-e", locus="main text Fig 6",
      readout_as_printed="Phenotype variations determined by MFI of CD86 and CD206",
      units_as_printed="MFI", assay="flow cytometry",
      cell_type="Kupffer cells", organ="liver",
      dose="1 mg/kg siLuciferase", genotype="WT and Cav1-/-",
      corpus_admissible="no",
      admissible_reason="macrophage activation markers. MFI of a HOST protein "
                        "that the cargo does not target is not protein "
                        "abundance of a delivered construct; this is the "
                        "cytokine/host-response class on the not-delivery list"),

    # ---------------- SI ----------------------------------------------------
    R(panel="SI Fig S1b", locus="SI p1-2",
      readout_as_printed="Comparison of MFI of barcodes conjugated with "
                         "Alexa-488 and Alexa-647 at 200 ng and 50 ng per well",
      units_as_printed="MFI (Above Background)", assay="flow cytometry",
      cell_type="iMAEC", corpus_admissible="no",
      admissible_reason="fluorophore brightness comparison; a reagent choice, "
                        "not a delivery measurement"),
    R(panel="SI Fig S1c-e", locus="SI p1-3",
      readout_as_printed="FACS marker definitions and representative gating",
      units_as_printed="none", assay="flow cytometry gating",
      corpus_admissible="no",
      admissible_reason="gating metadata, carries no measurement"),
    R(panel="SI Fig S1f", locus="SI p3",
      readout_as_printed="MFI of 12 cell types in WT and CavKO",
      units_as_printed="647 MFI (Above Background)", assay="flow cytometry",
      cell_type="12 sorted cell types", organ="liver, heart, lung, kidney",
      route=IV["route"], dose=IV["dose"], timepoint="1 h",
      genotype="WT and Cav1-/-", n_values_design_supported="24",
      n_values_attributable="0", corpus_admissible="no",
      admissible_reason="fluorophore biodistribution; the paper's own argument "
                        "is that this readout does NOT track the nucleic acid"),
    R(panel="SI Fig S2a-h", locus="SI p4-7",
      readout_as_printed="primer screens, gels, annealing-temperature and probe "
                         "optimisation, scrambled-probe control",
      units_as_printed="Signal / Noise (Au); gel images", assay="PCR / ddPCR",
      corpus_admissible="no",
      admissible_reason="assay development; no biological delivery measured"),
    R(panel="SI Fig S3a", locus="SI p8-9",
      readout_as_printed="Alexa-647 fluorescence 24 h after labeled barcodes "
                         "administered in vitro with Lipofectamine 2000",
      units_as_printed="MFI (above background) vs g/well", assay="flow cytometry",
      cell_type="iMAEC", timepoint="24 h", corpus_admissible="no",
      admissible_reason="in vitro uptake of a transfection reagent; uptake is "
                        "on the not-delivery list and " + NOT_TRANSLATED),
    R(panel="SI Fig S4a", locus="SI p10-11",
      readout_as_printed="156 8-nucleotide barcode sequences",
      units_as_printed="none", assay="sequence table", corpus_admissible="no",
      admissible_reason="barcode sequence registry; metadata"),
    R(panel="SI Fig S4b", locus="SI p12",
      readout_as_printed="DLS inclusion criteria: autocorrelation curves with 1 "
                         "inflection point and diameters between 20 and 200 nm",
      units_as_printed="nm; autocorrelation intensity", assay="DLS",
      corpus_admissible="no",
      admissible_reason="particle characterisation; on the not-delivery list"),
    R(panel="SI Fig S4c-d", locus="SI p12",
      readout_as_printed="Average normalized delivery of each LNP from library "
                         "2; naked barcode delivered less efficiently",
      units_as_printed="Normalized DNA Delivery (%)",
      assay="Illumina deep sequencing", cell_type="10 sorted cell types",
      organ="liver, lung, kidney, heart", route=IV["route"], dose=IV["dose"],
      timepoint="24 h", genotype="WT and Cav1-/-", pooled="yes",
      pool_identity="library 2 (7C1 / Chol / PEG, no DSPC), co-injected",
      pool_size="115 LNPs + 1 naked-barcode control",
      n_values_design_supported="115", n_values_attributable="0",
      corpus_admissible="no",
      admissible_reason=NOT_TRANSLATED + "; compositional; pool membership "
                        "unpublished"),
    R(panel="SI Fig S4e", locus="SI p13",
      readout_as_printed="The following example illustrates how our deep "
                         "sequencing data was normalized",
      units_as_printed="raw counts; Normalized Counts (%); Normalized to input "
                       "Counts (%)",
      assay="worked example", corpus_admissible="no",
      admissible_reason="ILLUSTRATIVE numbers for 3 fictitious LNPs, not "
                        "measurements. Extracted and reconstructed anyway "
                        "because it is the only place the normalisation is "
                        "pinned; see 01_si_tables.json"),
    R(panel="SI Fig S4f-k", locus="SI p10",
      readout_as_printed="FACS gating and per-LNP DLS for both libraries",
      units_as_printed="nm; gating", assay="DLS / flow cytometry",
      corpus_admissible="no",
      admissible_reason="characterisation and gating metadata"),
    R(panel="SI Fig S5a-e", locus="SI p16-17",
      readout_as_printed="Total ddPCR barcode counts for library 2 across "
                         "organs and liver cell types",
      units_as_printed="ddPCR counts; % Normalized Nanoparticle Uptake",
      assay="ddPCR / Illumina", cell_type="10 sorted cell types",
      organ="whole body, liver, lung, kidney", route=IV["route"],
      dose=IV["dose"], timepoint="24 h", genotype="WT and Cav1-/-",
      pooled="yes", pool_identity="library 2", pool_size="115",
      n_values_design_supported="28", n_values_attributable="0",
      corpus_admissible="no",
      admissible_reason="biodistribution / uptake; " + NOT_TRANSLATED),
    R(panel="SI Fig S5f", locus="SI p17",
      readout_as_printed="absolute delivery of 115 nanoparticles for each LNP "
                         "in the liver",
      units_as_printed="Absolute Delivery (AU)", assay="Illumina x ddPCR",
      cell_type="Kupffer, liver ECs, hepatocytes", organ="liver",
      route=IV["route"], dose=IV["dose"], timepoint="24 h",
      genotype="WT and Cav1-/-", pooled="yes", pool_identity="library 2",
      pool_size="115", n_values_design_supported="690",
      n_values_attributable="0", corpus_admissible="no",
      admissible_reason=NOT_TRANSLATED + "; unattributable scatter, 115 "
                        "unlabelled points per panel"),
    R(panel="SI Fig S6a-b", locus="SI p18-19",
      readout_as_printed="Normalized nanoparticle biodistribution across two "
                         "screens (226 LNPs) in liver, lung, heart, kidney ECs",
      units_as_printed="Normalized Nanoparticle Biodistribution (%)",
      assay="Illumina deep sequencing", cell_type="endothelial cells",
      organ="liver, lung, heart, kidney", route=IV["route"], dose=IV["dose"],
      timepoint="24 h", genotype="WT and Cav1-/-", pooled="yes",
      pool_identity="libraries 1 and 2 (separate injections, pooled analysis)",
      pool_size="226 (111 + 115)", n_values_design_supported="8",
      n_values_attributable="0", corpus_admissible="no",
      admissible_reason="biodistribution, explicitly named; " + NOT_TRANSLATED),
    R(panel="SI Fig S6c-e", locus="SI p18-19",
      readout_as_printed="Combined sequencing data and ddPCR results for each "
                         "LNP in lung, heart and kidney endothelial cells",
      units_as_printed="Nanoparticle Accumulation (Au)",
      assay="Illumina x ddPCR", cell_type="endothelial cells",
      organ="lung, heart, kidney", route=IV["route"], dose=IV["dose"],
      timepoint="24 h", genotype="WT and Cav1-/-", pooled="yes",
      pool_identity="libraries 1 and 2", pool_size="226",
      n_values_design_supported="1356", n_values_attributable="0",
      corpus_admissible="no",
      admissible_reason="the axis says ACCUMULATION. " + NOT_TRANSLATED
                        + "; 226 unlabelled points per panel"),
    R(panel="SI p14 table", locus="SI p14",
      readout_as_printed="Library 1 formulation table: Lipomer Type, Chol Type, "
                         "PEG Tail, PEG MW, HL Type, four mole percents, Diameter",
      units_as_printed="mole %; nm", assay="formulation record",
      n_values_design_supported="128", n_values_attributable="128",
      corpus_admissible="no",
      admissible_reason="formulation metadata and diameter; no endpoint. "
                        "Extracted in full to 01_blocks.csv"),
    R(panel="SI p15 table", locus="SI p15",
      readout_as_printed="Library 2 formulation table: LNP #, Lipomer Type, "
                         "Chol Type, PEG Alkyl, PEG MW, three mole percents, Diameter",
      units_as_printed="mole %; nm", assay="formulation record",
      n_values_design_supported="120", n_values_attributable="120",
      corpus_admissible="no",
      admissible_reason="formulation metadata and diameter; no endpoint. "
                        "Extracted in full to 01_blocks.csv"),

    # ---------------- the one readout class that WOULD be admissible --------
    R(panel="Methods 'Endothelial RNA interference' (NO FIGURE)",
      locus="main text Methods",
      readout_as_printed="7C1 carrying 1 mg/kg siCTRL (siLuc); 72 hours after "
                         "injection, tissues were isolated and protein "
                         "expression was determined via flow cytometry. ICAM2 "
                         "MFI in siLuc-treated mice (for each background) was "
                         "normalized to 100 percent",
      units_as_printed="ICAM2 MFI, % of siLuc",
      assay="flow cytometry (protein abundance)",
      cargo="siLuciferase (non-targeting control)", cargo_type="siRNA",
      cell_type="endothelial cells", route=IV["route"],
      dose="1 mg/kg siRNA", timepoint="72 h", genotype="WT and Cav1-/-",
      n_values_design_supported="0", n_values_attributable="0",
      corpus_admissible="no",
      admissible_reason="THE ONE ADMISSIBLE READOUT CLASS IN THE PAPER -- "
                        "protein-silencing MFI -- and it has no data. Only the "
                        "non-targeting control siLuc is ever injected, and its "
                        "ICAM2 MFI is DEFINED as 100%, so no knockdown is "
                        "measured. No figure, table or sentence anywhere "
                        "reports an ICAM2 silencing value; 'ICAM' occurs twice "
                        "in the whole article XML, both inside this one "
                        "Methods paragraph, and zero times in the SI. Recorded "
                        "as a contradiction (see 08_contradictions.csv)"),
]


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    with open(os.path.join(OUTDIR, "05_values.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        w.writeheader()
        for r in ROWS:
            w.writerows([r])

    # 06_rows.csv -- header-only, header read LIVE off the corpus
    vivo = os.path.join(REPO, "data", "data_vivo.csv")
    with open(vivo) as f:
        header = next(csv.reader(f))
    vitro = os.path.join(REPO, "data", "data_vitro.csv")
    with open(vitro) as f:
        header_vitro = next(csv.reader(f))
    if header != header_vitro:
        raise SystemExit("vitro and vivo headers differ; CLAUDE.md says they "
                         "are identical -- refusing to guess which to use")
    with open(os.path.join(OUTDIR, "06_rows.csv"), "w", newline="") as f:
        csv.writer(f).writerow(header)

    n_yes = sum(1 for r in ROWS if r["corpus_admissible"] == "yes")
    n_no = sum(1 for r in ROWS if r["corpus_admissible"] == "no")
    if n_yes + n_no != len(ROWS):
        raise SystemExit("counts do not sum: every row must carry yes or no")
    design = sum(int(r["n_values_design_supported"] or 0) for r in ROWS)
    attrib = sum(int(r["n_values_attributable"] or 0) for r in ROWS)
    summary = {
        "n_panels_inventoried": len(ROWS),
        "n_corpus_admissible_yes": n_yes,
        "n_corpus_admissible_no": n_no,
        "counts_sum": n_yes + n_no == len(ROWS),
        "n_values_design_supported_total": design,
        "n_values_attributable_total": attrib,
        "n_values_attributable_that_are_a_delivery_endpoint": 0,
        "rows_emitted": 0,
        "corpus_header_columns": len(header),
        "corpus_header_identical_vitro_vivo": True,
    }
    with open(os.path.join(OUTDIR, "05_values_summary.json"), "w") as f:
        json.dump(summary, f, indent=1)
    print(json.dumps(summary, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
