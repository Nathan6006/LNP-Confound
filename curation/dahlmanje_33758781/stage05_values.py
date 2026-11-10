#!/usr/bin/env python3
"""ZG_2020 -- every value locus in the paper, with an admissibility verdict.

Run:  python3 -I stage05_values.py   ->  05_values.csv, 06_rows.csv

NO VALUE IS DIGITISED HERE, ON PURPOSE. Every locus below ends with
corpus_admissible = no, for a reason recorded on the row, so digitising would
produce numbers nothing consumes. The only figures obtainable are 709 px web
JPEGs with no page rect (dpi not derivable), and CLAUDE.md section 7 forbids
inventing precision we cannot state. `value_as_printed` therefore records the
READOUT AND UNITS as the paper prints them, and `value` records why no number
was taken.

06_rows.csv is written HEADER-ONLY, from the live 79-column header of
data/data_vivo.csv, so that the zero is visible rather than implied.
"""
import csv, json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))

# Outputs go to OUTDIR, which is HERE for a normal run. The seed harnesses
# redirect it to a temp directory: a MUTATED copy must never overwrite the
# real curation outputs. It did, and the determinism check caught it --
# reference_structures.csv was left holding seed G3's REGIOISOMER anchors.
OUTDIR = HERE
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))

POOL = ("pooled: one intravenous injection of all 109 barcoded LNPs per "
        "mouse, pool size 109, total nucleic acid 1.0 mg/kg "
        "(0.009 mg/kg per particle on average)")
NOT_POOLED = "not pooled: LNP1 injected on its own at 0.5 mg/kg"

CELLS_FIG3B = [
    "Liver Kupffer Cells", "Liver ECs", "Liver Hepatocytes", "Spleen B Cells",
    "Spleen DCs", "Spleen Monocytes", "Spleen T Cells", "Lung ECs",
    "Lung Dentritic Cells", "Lung EPs", "Bone Marrow ECs",
    "Bone Marrow Monocytes", "Bone Marrow HSCs", "Thymus Cytotoxic T cells",
    "Thymus Helper T Cells", "Thymus Other Thymocytes",
]
PHOSPHOLIPIDS = ["A-9", "A-10", "A-11", "A-12", "A-13", "A-14", "A-15",
                 "A-16", "A-17", "A-17-2Z"]
STEROLS = ["Cholesterol", "20alpha-OH Cholesterol", "beta-Sitosterol"]

# reason codes, written once so the counts below cannot drift from the prose
R_ONE_IL = ("NO - every row would carry the SAME ionizable lipid (cKK-E15); "
            "the varying factor is the helper phospholipid / sterol / molar "
            "ratio, so a corpus keyed on ionizable-lipid chemistry gains "
            "nothing from additional rows")
R_POOL_LEVEL = ("NO - this is a POOL-level measurement of all 109 LNPs "
                "together, not a per-LNP value; it cannot be attributed to "
                "any formulation")
R_UNATTRIBUTED = ("NO - the 109 individual values exist in the authors' hands "
                  "but are plotted UNLABELLED, so no point can be attributed "
                  "to an LNP; and the 11 of 120 formulations dropped at QC "
                  "are never identified, so the 109 members cannot even be "
                  "enumerated")
R_COLLAPSED = ("NO - a mean over 3 sterols x 4 molar ratios x 16 cell types "
               "(~11 LNPs per bar), so CHL_name, the four *_molratio columns "
               "and Model_target are all unattributable")
R_NOT_DELIVERY = "NO - not a delivery endpoint under the project's vocabulary"
R_DERIVED = ("NO - a derived enrichment statistic (top-10% minus bottom-10% "
             "frequency), not a delivery measurement")
R_UNAVAILABLE = ("NO - panel not served with the article; Appendix S1 contains "
                 "Methods, NMR and HRMS only, and Supplementary Figures S1-S3 "
                 "are absent from every copy we can obtain")

ROWS = []


def add(**kw):
    ROWS.append(kw)


def main():
    common = dict(
        cargo="Cre mRNA (co-encapsulated with a single-stranded DNA barcode "
              "at a 10:1 mRNA:barcode mass ratio)",
        cargo_type="mRNA", organism="Mouse", strain="Ai14 (LSL-tdTomato) "
        "on C57BL/6", route="intravenous, lateral tail vein",
        timepoint="72 h post-injection", n_mice="N = 3-5 per group",
        ionizable_lipid="cKK-E15 (the SAME on every LNP in this paper)",
    )

    # --- Fig 3b : pool-level functional delivery, 16 cell types ------------
    for c in CELLS_FIG3B:
        add(locus="Figure 3b", item=c,
            readout_as_printed="tdTom+Cells above background",
            units_as_printed="% of cells (axis 0-70)",
            assay="flow cytometry of tdTomato fluorescence in FACS-gated "
                  "populations; Cre mRNA -> Cre protein -> Lox-Stop-Lox "
                  "excision -> tdTomato",
            requires_expression="yes", compositional="no",
            organ=c.split()[0], cell_type=c, dose="1.0 mg/kg total nucleic acid",
            pool=POOL, value="not digitised - no row emitted",
            value_as_printed="bar + SEM, 709 px raster, not transcribed",
            corpus_admissible="no",
            reason=f"{R_POOL_LEVEL}. {R_ONE_IL}.", **common)

    # --- Fig 3c : the 109 individual normalized-delivery points ------------
    add(locus="Figure 3c", item="109 individual LNPs (unlabelled points)",
        readout_as_printed="Normalized Delivery",
        units_as_printed="dimensionless, log axis 0.1 - 10, centred on 1",
        assay="next-generation sequencing of DNA barcodes recovered from "
              "tdTomato+ FACS-sorted cells",
        requires_expression="yes",
        compositional="yes - a share of the pooled barcode total, then "
                      "rescaled against the injected input mixture",
        organ="all", cell_type="all 16 populations",
        dose="1.0 mg/kg total nucleic acid", pool=POOL,
        value="not digitised - no row emitted",
        value_as_printed="~109 overlapping grey points, no labels",
        corpus_admissible="no",
        reason=f"{R_UNATTRIBUTED}. {R_ONE_IL}.", **common)
    add(locus="Figure 3c", item="CTRL (free, unencapsulated DNA barcode)",
        readout_as_printed="Normalized Delivery",
        units_as_printed="dimensionless, log axis 0.1 - 10",
        assay="same sequencing readout", requires_expression="yes",
        compositional="yes", organ="all", cell_type="all 16 populations",
        dose="1.0 mg/kg total nucleic acid", pool=POOL,
        value="not digitised - no row emitted",
        value_as_printed="green point, below the LNP cloud",
        corpus_admissible="no",
        reason="NO - a negative control (naked DNA), not an LNP formulation",
        **common)

    # --- Fig 3d / 3e : spread and biodistribution --------------------------
    for organ_cell in ["Kupffer Cells", "ECs", "Hepatocytes"]:
        add(locus="Figure 3d", item=organ_cell,
            readout_as_printed="Standard Deviation",
            units_as_printed="dimensionless (axis 0.0-0.6)",
            assay="spread of normalized delivery across mouse replicates",
            requires_expression="yes", compositional="n/a",
            organ="Liver", cell_type=organ_cell,
            dose="1.0 mg/kg total nucleic acid", pool=POOL,
            value="not digitised - no row emitted",
            value_as_printed="bar + SEM", corpus_admissible="no",
            reason=f"{R_NOT_DELIVERY} - it is a reproducibility statistic "
                   "about the measurement, not a measurement of delivery",
            **common)
        add(locus="Figure 3e", item=organ_cell,
            readout_as_printed="Barcode Counts (AU)",
            units_as_printed="arbitrary units (axis 0-10000)",
            assay="droplet digital PCR (QUANT) of absolute barcode in sorted "
                  "cells",
            requires_expression="no", compositional="no",
            organ="Liver", cell_type=organ_cell,
            dose="1.0 mg/kg total nucleic acid", pool=POOL,
            value="not digitised - no row emitted",
            value_as_printed="bar + SEM", corpus_admissible="no",
            reason=f"{R_NOT_DELIVERY} - this is PARTICLE BIODISTRIBUTION "
                   "(cargo abundance measured by ddPCR, before and "
                   "independent of translation), which the project's endpoint "
                   "vocabulary excludes. It is also pool-level.",
            **common)

    # --- Fig 3f / 3g : the group means -------------------------------------
    for p in PHOSPHOLIPIDS:
        add(locus="Figure 3f", item=f"phospholipid {p}",
            readout_as_printed="Normalized Delivery",
            units_as_printed="dimensionless (axis 0.00-1.25)",
            assay="barcode sequencing from tdTomato+ sorted cells, averaged",
            requires_expression="yes", compositional="yes",
            organ="all", cell_type="mean over all cell types",
            dose="1.0 mg/kg total nucleic acid", pool=POOL,
            value="not digitised - no row emitted",
            value_as_printed="bar + SEM", corpus_admissible="no",
            reason=f"{R_COLLAPSED}. {R_ONE_IL}. The paper's own reading of "
                   "this panel is that it shows no substantial difference.",
            **common)
    for s in STEROLS:
        add(locus="Figure 3g", item=f"sterol {s}",
            readout_as_printed="Normalized Delivery",
            units_as_printed="dimensionless (axis 0.00-1.25)",
            assay="barcode sequencing from tdTomato+ sorted cells, averaged",
            requires_expression="yes", compositional="yes",
            organ="all", cell_type="mean over all cell types",
            dose="1.0 mg/kg total nucleic acid", pool=POOL,
            value="not digitised - no row emitted",
            value_as_printed="bar + SEM", corpus_admissible="no",
            reason=f"{R_COLLAPSED} (here over 10 phospholipids x 4 ratios). "
                   f"{R_ONE_IL}.", **common)

    # --- Fig 3h : enrichment -----------------------------------------------
    for p in PHOSPHOLIPIDS:
        add(locus="Figure 3h", item=f"phospholipid {p}",
            readout_as_printed="Enrichment",
            units_as_printed="fold enrichment change (axis -3 to +3)",
            assay="frequency in the top 10% minus frequency in the bottom 10%",
            requires_expression="yes", compositional="n/a",
            organ="Liver", cell_type="Kupffer / EC / hepatocyte",
            dose="1.0 mg/kg total nucleic acid", pool=POOL,
            value="SIGN ONLY, by colour identity: see 04_figure_checks.json",
            value_as_printed="teal bar (enriched) or grey bar (depleted)",
            corpus_admissible="no", reason=f"{R_DERIVED}. {R_ONE_IL}.",
            **common)

    # --- Fig 3i : the one individually-injected LNP ------------------------
    for c in ["Kupffer Cells", "ECs", "Hepatocytes"]:
        add(locus="Figure 3i", item=f"LNP1 (phospholipid A-11) -> {c}",
            readout_as_printed="tdTom+Cells above background",
            units_as_printed="% of cells (axis 0-100)",
            assay="flow cytometry, same Cre->tdTomato readout",
            requires_expression="yes", compositional="no",
            organ="Liver", cell_type=c, dose="0.5 mg/kg", pool=NOT_POOLED,
            value="not digitised - no row emitted",
            value_as_printed="bar + SEM, **** p<.0001 two-way ANOVA",
            corpus_admissible="no",
            reason="NO - this IS an individually measured, individually "
                   "injected LNP and it is the closest this paper comes to an "
                   "emittable row, but its FORMULATION IS INCOMPLETE: the "
                   "paper names only the phospholipid (A-11) and never says "
                   "which of the 3 sterols or which of the 4 molar ratios "
                   "LNP1 uses, so CHL_name/CHL_SMILES and all four "
                   "*_molratio columns would have to be guessed. "
                   f"{R_ONE_IL}.", **common)

    # --- Fig 2d/2e/2f, S1, S2 ---------------------------------------------
    for loc, it, rd in [("Figure 2d", "diameter by phospholipid (10 bars)",
                         "Diameter (nm)"),
                        ("Figure 2e", "diameter by sterol (3 bars)",
                         "Diameter (nm)"),
                        ("Figure 2f", "diameter of the 109-LNP pool",
                         "Diameter (nm)")]:
        add(locus=loc, item=it, readout_as_printed=rd,
            units_as_printed="nm", assay="dynamic light scattering",
            requires_expression="no", compositional="no", organ="n/a",
            cell_type="n/a", dose="n/a", pool="n/a (in vitro QC)",
            value="not digitised - no row emitted",
            value_as_printed="bar + SD", corpus_admissible="no",
            reason=f"{R_NOT_DELIVERY} - particle size", **common)
    add(locus="Figure S1c / S2b", item="mouse body weight",
        readout_as_printed="body weight vs PBS", units_as_printed="n/a",
        assay="weighing", requires_expression="no", compositional="no",
        organ="n/a", cell_type="n/a", dose="1.0 / 0.5 mg/kg",
        pool="both", value="panel unavailable",
        value_as_printed="n/a", corpus_admissible="no",
        reason=f"{R_NOT_DELIVERY} - tolerability. {R_UNAVAILABLE}", **common)
    add(locus="Figure S2a", item="LNP1 -> spleen T cells / B cells / monocytes",
        readout_as_printed="tdTom+ cells", units_as_printed="% of cells",
        assay="flow cytometry", requires_expression="yes", compositional="no",
        organ="Spleen", cell_type="T cells, B cells, monocytes",
        dose="0.5 mg/kg", pool=NOT_POOLED,
        value="panel unavailable; the text states no tdTomato+ cells",
        value_as_printed="n/a", corpus_admissible="no",
        reason=f"{R_UNAVAILABLE}. Even with the panel, {R_ONE_IL}.", **common)

    # ---------------------------------------------------------------------
    fields = list(ROWS[0].keys())
    with open(os.path.join(OUTDIR, "05_values.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(ROWS)

    yes = sum(1 for r in ROWS if r["corpus_admissible"] == "yes")
    no = sum(1 for r in ROWS if r["corpus_admissible"] == "no")
    assert yes + no == len(ROWS), "admissibility counts must sum"

    # 06_rows.csv: HEADER ONLY, read live from the corpus so the schema cannot
    # drift, and so the zero is visible instead of implied.
    with open(os.path.join(ROOT, "data", "data_vivo.csv")) as fh:
        header = fh.readline().rstrip("\n")
    with open(os.path.join(OUTDIR, "06_rows.csv"), "w") as fh:
        fh.write(header + "\n")
    n_cols = len(header.split(","))

    summary = {
        "value_loci_recorded": len(ROWS),
        "corpus_admissible_yes": yes,
        "corpus_admissible_no": no,
        "counts_sum": yes + no == len(ROWS),
        "rows_emitted": 0,
        "06_rows_csv": f"header only, {n_cols} columns read live from "
                       f"data/data_vivo.csv",
        "distinct_ionizable_lipids_contributed": 1,
    }
    print("METRIC values_summary=" + json.dumps(summary))
    with open(os.path.join(OUTDIR, "05_values_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1)
    if n_cols != 79:
        raise SystemExit(f"GATE FAILED: expected 79 corpus columns, got {n_cols}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
