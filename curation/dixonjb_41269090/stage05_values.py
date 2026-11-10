#!/usr/bin/env python3
"""Stage 05 -- every value in the paper, classified, with corpus_admissible on
every row and the counts summing.

The quantity. For a given site and administration the paper reports two average
counts: LEC (or other cell type) count, and aVHH+ cell count. The delivery
readout is their ratio -- the fraction of that cell population expressing the
protein translated from the delivered mRNA.

WHAT THIS RATIO IS NOT. It is a RATIO OF MEANS. The figures plot one point per
independent experiment and the mean of those per-animal PERCENTAGES, which is a
MEAN OF RATIOS. They are not the same statistic and they disagree by up to 2.5x
on one cell of this paper. The ratio of means is the only per-formulation number
derivable from material we hold; it is emitted with that stated on every row.

Two cohorts, and they must not be pooled:

  study_A  Figure 2 / Table S2.  Six leads + saline, each administered on its
           own, intradermally into BOTH paws. Sites ALN, BLN, PLN.
  study_B  Figure 4 / Tables S3-S4.  Saline, free aVHH mRNA, MC3 and LNP7, ONE
           hindlimb injected, so ALN, BLN and PLN Contralateral are off-target
           control sites while PLN Injected and PLV drain the injection.

LNP7 is in both, which gives a within-paper shared member and therefore a
measurable ratio (07_qc.json, check cross_study_LNP7_ratio).
"""
import csv
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

CELLS = ["Macrophages", "Monocytes", "cDC2", "cDC1", "DNs", "FRCs",
         "LECs", "BECs", "VECs"]
LNPS = {"LNP1", "LNP2", "LNP3", "LNP4", "LNP7", "LNP11", "MC3"}

SITE_TARGET = {
    "ALN": "lymph_node_axillary",
    "BLN": "lymph_node_brachial",
    "PLN": "lymph_node_popliteal",
    "PLN Contralateral": "lymph_node_popliteal_contralateral",
    "PLN Injected": "lymph_node_popliteal",
    "PLV": "lymphatic_vessel_popliteal",
}

NOT_LNP = ("not an LNP: {} carries no ionizable lipid and no formulation, so "
           "there is no structure for a corpus row to key on. Recorded here so "
           "the assay background is visible rather than dropped.")

BYSTANDER = ("cargo-expression measurement, but in a BYSTANDER cell population "
             "({}), not the targeted lymphatic endothelium. The corpus "
             "Model_target vocabulary is organ-level (liver, muscle, spleen, "
             "ear, brain, lung_epithelium, bone_marrow, whole_body) and has no "
             "term for a FACS-sorted subset inside an organ; emitting nine "
             "subsets per site would invent that convention on the strength of "
             "two ionizable lipids, and would give one formulation nine "
             "correlated values inside one z-block. See DECISIONS.md D-3.")

ADMISSIBLE = ("delivery endpoint: fraction of lymphatic endothelial cells "
              "(CD45-/CD31+/PDPN+) expressing the protein translated from the "
              "delivered aVHH mRNA, 12-16 h after intradermal administration "
              "of a single named formulation. Cargo expression, not particle "
              "accumulation.")

NON_DELIVERY = {
    "fig1D_hydrodynamic_diameter": "particle size, not a delivery endpoint",
    "fig1E_PDI": "polydispersity, not a delivery endpoint",
    "fig2A_diameter_PDI": "particle characterisation, not a delivery endpoint",
    "figS4a_zeta_potential": "zeta potential, not a delivery endpoint",
    "figS4a_encapsulation_efficiency": "encapsulation efficiency, not a delivery endpoint",
    "figS9_characterisation": "diameter/PDI/EE/mRNA concentration, not a delivery endpoint",
    "fig3A_model_MSE_R2": "model performance, not a measurement",
    "fig3B_SHAP": "model attribution, not a measurement",
    "figS5_coefficients_importances": "model attribution, not a measurement",
    "fig5B_Cy3_radiant_efficiency": ("PARTICLE/CARGO-LABEL BIODISTRIBUTION: IVIS "
                                     "fluorescence of Cy3-conjugated miRNA with no "
                                     "expression step. This is the 'lymphatic targeting "
                                     "measured as accumulation' case; it is present here "
                                     "alongside, not instead of, the expression readout"),
    "fig6C_PDPN_DAPI_area": "lymphatic vessel area, a therapeutic phenotype",
    "fig6D_EdU_PROX1_counts": "LEC proliferation, a therapeutic phenotype",
    "fig7_VEGFC_ELISA": ("protein product of the delivered mRNA, so cargo expression in "
                         "kind, but excluded: the values are plotted only and never "
                         "printed, the panels are not held, the quantity is an optical "
                         "density NORMALISED TO A DIFFERENT-CARGO CONTROL GROUP rather "
                         "than a per-formulation delivery value, and it adds no "
                         "formulation beyond LNP7 and MC3 already emitted"),
    "fig8_fig9_NIR_function_oedema": "lymphatic pump function and tail oedema, therapeutic phenotype",
    "figS2_per_screen_LEC_counts": ("per-SCREEN aggregate over the whole pooled library; "
                                    "no value attaches to a formulation"),
    "figS3a_SANDS_histogram": ("histogram of the pooled-barcode delivery distribution: bin "
                               "heights are counts of LNPs, so no per-member value exists "
                               "even in principle"),
    "figS3b_UMAP": "embedding coordinates of composition, unitless",
    "figS7_viability": "Alamar Blue / Live-Dead viability, not a delivery endpoint",
    "figS10_IV_vs_ID": ("route comparison DURING THE POOLED SCREEN: the administered "
                        "article is the pool, so no value attaches to a formulation"),
    "SANDS_99_LNP_barcode_values": ("POOLED BARCODE SCREEN. Values are normalised barcode "
                                    "counts, relative and competitive -- each depends on the "
                                    "other members of the pool. Moot in any case: they are "
                                    "published nowhere. See DECISIONS.md D-1"),
}


def pct(num, den):
    return 100.0 * num / den


def main():
    T = json.load(open(os.path.join(HERE, "04_si_tables.json")))
    out = []

    # ---- study A, Table S2 (LECs only; the table has no other cell type) ----
    for site, adm, lec, pos in T["table_S2"]["rows"]:
        v = pct(pos, lec)
        adm_is_lnp = adm in LNPS
        out.append({
            "value_id": f"A|{site}|{adm}|LECs",
            "study": "study_A",
            "source_location": "SI Table S2, page 22 (bitmap; read as printed numerals)",
            "figure_equivalent": {"ALN": "Figure 2B", "BLN": "Figure 2C",
                                  "PLN": "Figure 2D"}[site],
            "administration": adm,
            "site": site,
            "cell_type": "LECs",
            "readout_as_printed": "Average LEC count; Average aVHH+ Cell Count",
            "units_as_printed": "cells (counts)",
            "derived_quantity": "percent of LECs that are aVHH+",
            "numerator_aVHH_pos": pos,
            "denominator_cells": lec,
            "value_percent": round(v, 4),
            "cargo": "aVHH mRNA (anchored VHH anti-RSV single-domain antibody, DAF-GPI anchored)",
            "cargo_type": "aVHH",
            "model": "C57BL/6 mouse, female",
            "route": "intradermal, both paws",
            "dose_mg_per_kg": None,
            "dose_note": "not stated for the single-LNP validation studies; "
                         "1.5 mg/kg is stated only for the pooled SANDS screen",
            "timepoint": "12-16 h",
            "assay": "flow cytometry, Live/Dead then CD45-/CD31+/PDPN+ LEC gate, aVHH+ fraction",
            "pool_identity": "single LNP, administered individually",
            "corpus_admissible": "yes" if adm_is_lnp else "no",
            "admissibility_reason": ADMISSIBLE if adm_is_lnp else NOT_LNP.format(adm),
        })

    # ---- study B, Tables S3 and S4, every cell type ------------------------
    s3 = {(r[0], r[1]): r[2:] for r in T["table_S3"]["rows"]}
    s4 = {(r[0], r[1]): r[2:] for r in T["table_S4"]["rows"]}
    for site, adm in s3:
        tot, pos_row = s3[(site, adm)], s4[(site, adm)]
        for i, cell in enumerate(CELLS):
            if tot[i] is None or pos_row[i] is None:
                continue
            v = pct(pos_row[i], tot[i])
            adm_is_lnp = adm in LNPS
            is_lec = cell == "LECs"
            if not adm_is_lnp:
                ok, why = "no", NOT_LNP.format(adm)
            elif not is_lec:
                ok, why = "no", BYSTANDER.format(cell)
            else:
                ok, why = "yes", ADMISSIBLE
            out.append({
                "value_id": f"B|{site}|{adm}|{cell}",
                "study": "study_B",
                "source_location": "SI Tables S3 and S4, pages 23-24 "
                                   "(bitmaps; read as printed numerals)",
                "figure_equivalent": ("Figure 4A" if (is_lec and site in
                                      ("PLN Injected", "PLV")) else
                                      "Figure 4B / Figure S6"),
                "administration": adm,
                "site": site,
                "cell_type": cell,
                "readout_as_printed": "Average cell count (S3); Average aVHH+ cell count (S4)",
                "units_as_printed": "cells (counts)",
                "derived_quantity": f"percent of {cell} that are aVHH+",
                "numerator_aVHH_pos": pos_row[i],
                "denominator_cells": tot[i],
                "value_percent": round(v, 4),
                "cargo": "aVHH mRNA",
                "cargo_type": "aVHH",
                "model": "C57BL/6 mouse, female",
                "route": "intradermal, ONE hindlimb",
                "dose_mg_per_kg": None,
                "dose_note": "not stated for the single-LNP validation studies",
                "timepoint": "12-16 h",
                "assay": "12-channel flow cytometry, gating in SI Figs S12/S13",
                "pool_identity": "single LNP, administered individually",
                "corpus_admissible": ok,
                "admissibility_reason": why,
            })

    # ---- every non-delivery or unextractable panel, recorded explicitly -----
    for panel, why in NON_DELIVERY.items():
        out.append({
            "value_id": f"X|{panel}",
            "study": "whole_paper",
            "source_location": panel,
            "figure_equivalent": panel,
            "administration": "", "site": "", "cell_type": "",
            "readout_as_printed": "", "units_as_printed": "",
            "derived_quantity": "", "numerator_aVHH_pos": "",
            "denominator_cells": "", "value_percent": "",
            "cargo": "", "cargo_type": "", "model": "", "route": "",
            "dose_mg_per_kg": "", "dose_note": "", "timepoint": "",
            "assay": "", "pool_identity": "",
            "corpus_admissible": "no",
            "admissibility_reason": why,
        })

    cols = list(out[0])
    with open(os.path.join(HERE, "05_values.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(out)

    yes = [r for r in out if r["corpus_admissible"] == "yes"]
    no = [r for r in out if r["corpus_admissible"] == "no"]
    assert len(yes) + len(no) == len(out)
    print(f"05_values.csv   {len(out)} rows = {len(yes)} admissible + "
          f"{len(no)} inadmissible")
    from collections import Counter
    print("  admissible by study:", dict(Counter(r["study"] for r in yes)))
    c = Counter()
    for r in no:
        c[r["admissibility_reason"].split(":")[0].split(".")[0][:60]] += 1
    for k, n in c.most_common():
        print(f"  inadmissible {n:>3}  {k}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
