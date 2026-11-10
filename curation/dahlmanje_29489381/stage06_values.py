#!/usr/bin/env python3
"""05_values.csv -- every per-LNP readout this paper makes available to us.

There are 87 of them and NOT ONE IS A DELIVERY VALUE.  They are ranks, and the
file says so on every row, because a rank column named 'Heart Ecs' sitting
next to a composition is exactly the thing a later reader would mistake for a
measurement.

Each row carries its own source locus (SI figure panel and page), the readout
and units exactly as the panel prints them, and the experimental context
(cargo, organ/cell, route, dose, timepoint) from Methods.  Where Methods does
not state something for THIS panel, the field says 'not stated' and names what
the surrounding panels do state -- it is never filled by inference.
"""
import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

PANEL = {("Library 1", "20"): "Fig S4e",
         ("Library 2", "20"): "Fig S4f",
         ("Library 3", "21"): "Fig S4g"}

INVITRO = {"iMAEC": ("immortalised mouse aortic endothelial cell line "
                     "(iMAEC)", "in_vitro"),
           "RAWs": ("RAW 264.7 mouse macrophage cell line", "in_vitro")}

INVIVO = {
    "Heart Ecs": ("heart endothelial cells (CD31+CD45-)", "in_vivo"),
    "Lung Ecs": ("lung endothelial cells (CD31+CD45-)", "in_vivo"),
    "Marrow Ecs": ("bone marrow endothelial cells (CD31+)", "in_vivo"),
    "Heart Macs": ("heart macrophages (CD31-CD45+CD11b+)", "in_vivo"),
    "Lung Macs": ("lung macrophages (CD31-CD45+CD11b+)", "in_vivo"),
    "Marrow Macs": ("bone marrow macrophages (CD31-CD45+CD11b+)", "in_vivo"),
}

CTX_VITRO = {
    "cargo": "56-nt single-stranded DNA barcode, phosphorothioate-capped",
    "route": "added to culture medium, 24-well plate, 40,000 cells/well",
    "dose": ("20 ng total DNA per well -- Methods: 'cells were treated with "
             "20 ng total DNA in all other experiments'. NOT printed on the "
             "Fig S4 panel itself."),
    "timepoint": ("72 h -- the library's main in vitro panels (Fig S3j-l, "
                  "S3m-o) are labelled '72 hours'. NOT printed on the Fig S4 "
                  "panel itself; library 2 also has 4 h and 48 h datasets."),
    "animal": "",
}
CTX_VIVO = {
    "cargo": "56-nt single-stranded DNA barcode, phosphorothioate-capped",
    "route": "intravenous, lateral tail vein",
    "dose": ("0.5 mg/kg TOTAL DNA for the whole pooled library; the per-LNP "
             "dose is that divided by the number of pooled LNPs (Methods "
             "gives 0.005 mg/kg for a 100-particle pool)"),
    "timepoint": ("72 h -- Methods: 'In all other cases, tissues and cells "
                  "were isolated 72 hours after injection'. NOT printed on "
                  "the Fig S4 panel."),
    "animal": "female C57BL/6J, 4-12 weeks, N = 3-5 per group",
}

READOUT = {
    "in_vitro": ("rank within this library's stable pool, best = 1 "
                 "(unitless; derived by the authors from normalized "
                 "delivery, % of pooled barcode reads)"),
    "in_vivo": ("rank WITHIN THE TOP 10% of this cell type, best = 1 "
                "(unitless, CENSORED: only top-10% particles are listed, so "
                "a blank means 'outside the top 10%', not 'not measured')"),
}


def main():
    with open(os.path.join(HERE, "02_rank_tables.csv"), newline="") as fh:
        rank = list(csv.DictReader(fh))

    out = []
    for r in rank:
        panel = PANEL[(r["subtable"], r["si_page"])]
        comp = {
            "lipomer_id": r.get("Lipomer", ""),
            "amine": r.get("Amine", "") or r.get("Amine Lipid", ""),
            "tail_lipid": r.get("Lipid", ""),
            "helper_lipid": r.get("Phospholipid", ""),
            "mass_ratio": r.get("Mass Ratio", "")
                          or r.get("Mass RatioPEG Chain", ""),
            "peg_chain": r.get("PEG Chain", "")
                         or r.get("Mass RatioPEG Chain #2", ""),
            "peg_mw": r.get("PEG MW", ""),
            "lipomer_mol_pct": r.get("Lipomer Mole %", ""),
            "cholesterol_mol_pct": r.get("Cholesterol Mole %", ""),
            "peg_mol_pct": r.get("PEG Mole %", ""),
            "helper_lipid_mol_pct": r.get("Helper Lipid Mole %", ""),
            "diameter_nm": r.get("Diameter (nm)", ""),
        }
        for col, (tissue, model) in list(INVITRO.items()) + list(INVIVO.items()):
            v = r.get(col, "")
            if v in ("", None):
                continue
            ctx = CTX_VITRO if model == "in_vitro" else CTX_VIVO
            row = {"library": r["subtable"].replace("Library ", ""),
                   "source_figure": panel,
                   "source_page_si_pdf": r["si_page"],
                   "source_column_as_printed": col,
                   "model": model,
                   "organ_or_cell": tissue,
                   "readout_as_printed": col,
                   "readout_quantity": READOUT[model],
                   "value": int(v),
                   "value_is_a_rank": "yes",
                   "value_is_a_delivery_measurement": "no"}
            row.update(comp)
            row.update(ctx)
            row["emittable_to_corpus"] = "no"
            row["not_emittable_because"] = (
                "the lipomer has no determinable SMILES (see "
                "03_enumerated.csv) AND a censored rank over an "
                "outcome-selected subset is not an Experiment_value")
            out.append(row)

    with open(os.path.join(HERE, "05_values.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]))
        w.writeheader()
        for r in out:
            w.writerow(r)

    rep = {
        "n_rank_table_rows": len(rank),
        "n_value_cells": len(out),
        "n_in_vitro_cells": sum(1 for r in out if r["model"] == "in_vitro"),
        "n_in_vivo_cells": sum(1 for r in out if r["model"] == "in_vivo"),
        "n_delivery_values": 0,
        "n_emittable_to_corpus": 0,
        "by_panel": {},
    }
    for r in out:
        k = "%s %s" % (r["source_figure"], r["source_column_as_printed"])
        rep["by_panel"][k] = rep["by_panel"].get(k, 0) + 1
    with open(os.path.join(HERE, "05_values_qc.json"), "w") as fh:
        json.dump(rep, fh, indent=1)
    print(json.dumps(rep, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
