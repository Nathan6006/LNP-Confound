#!/usr/bin/env python3
"""peer_39912611 stage 06 -- 05_values.csv and the corpus-schema rows.

Nothing is merged. This writes what a merge WOULD take, so a reviewer can see
every value with its source location beside it, and so the quality flags are on
the rows rather than in a QC file (F-001: a flag that lives only in a per-folder
QC file is not a flag, it is a memo).
"""
import csv
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "lib"))
import emit                                                   # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
EXP_VITRO = "GN_2025"        # first author Gonna Somu Naidu, 2025. PROPOSED.
EXP_VIVO = "GN_2025"

PKA = {"Lipid 14": 6.44, "Lipid 32": 5.92, "Lipid 33": 6.21, "Lipid 34": 6.24,
       "Lipid 35": 7.67, "Lipid 36": 7.14, "Lipid 37": 6.08, "Lipid 38": 6.17,
       "Lipid 39": 5.93, "Lipid 40": 6.86, "Lipid 41": 7.20, "Lipid 42": 6.90}

COMMON = {
    "IL_molratio": 40.0, "HL_name": "DSPC", "HL_molratio": 10.0,
    "CHL_name": "Cholesterol", "CHL_molratio": 48.5,
    "PEG_name": "PEG-DMG2000", "PEG_molratio": 1.5,
    "IL_to_nucleicacid_chargeratio": 9.0,
    "Aqueous_buffer": "25 mM citrate pH 4.5",
    "Dialysis_buffer": "1x PBS",
    "Mixing_method": "microfluidic_NanoAssemblr_FRR3to1_TFR12",
    "Cargo": "FLuc", "Cargo_type": "mRNA",
    "Publication_DOI": "10.1021/acsnano.4c18636",
    "Publication_PMID": "39912611",
    "Publication_link": "https://doi.org/10.1021/acsnano.4c18636",
}


def read(name):
    with open(os.path.join(HERE, name)) as fh:
        return list(csv.DictReader(fh))


def main():
    struct = {r["lipid_name"]: r for r in read("03_enumerated.csv")}
    f2 = read("04_fig2_values.csv")
    f3 = read("05_fig3_values.csv")

    vals = []
    for r in f2:
        vals.append({
            "lipid_name": r["lipid_name"],
            "pool": "in_vitro",
            "source_location": (f"{r['panel']} ({r['cell_line']}), bar "
                                f"{r['bar_index']} of 36, x=[{r['x0']},"
                                f"{r['x1']}), top row {r['top_row_px']} px"),
            "model_type": r["cell_line"],
            "model_target": "in_vitro",
            "route": "in_vitro",
            "cargo": "firefly luciferase mRNA (mLuc), chemically modified",
            "dose": r["dose_ug_per_mL"],
            "dose_unit": "ug/mL (AXIS LABEL -- contradicts Methods, see DECISIONS D-1)",
            "timepoint": "24 h after treatment",
            "assay": ("luciferase luminescence, Promega Luciferase Assay "
                      "System on a GloMax reader, 10 s exposure, 30 uL lysate"),
            "replicates": "n=3 per caption / quadruplicate per Methods (08_contradictions row 4)",
            "value": r["luminescence"],
            "value_unit": "relative luminescence units",
            "log10": r["log10"],
            "status": r["status"],
            "value_is_bound": "",
        })
    for r in f3:
        if r["lipid_name"] == "Untreated":
            continue
        vals.append({
            "lipid_name": r["lipid_name"],
            "pool": "in_vivo",
            "source_location": (f"{r['panel']} ({r['organ']}), category "
                                f"{r['cat_index']} of 13, top row "
                                f"{r['top_row_px']} px"),
            "model_type": "Mouse_C57BL6J_female_8to10wk",
            "model_target": r["organ"],
            "route": "intravenous_retro_orbital",
            "cargo": "firefly luciferase mRNA (mLuc), chemically modified",
            "dose": "10",
            "dose_unit": "ug mRNA per animal",
            "timepoint": "6 h post injection (D-luciferin IP 5 min before)",
            "assay": ("ex vivo IVIS Lumina, average radiance in an ad hoc ROI,"
                      " Living Image 4.1"),
            "replicates": "n=4 mice per group",
            "value": r["avg_radiance_p_s_cm2_sr"],
            "value_unit": "p/s/cm2/sr",
            "log10": r["log10"],
            "status": r["status"],
            "value_is_bound": r["value_is_bound"],
        })

    with open(os.path.join(HERE, "05_values.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(vals[0].keys()))
        w.writeheader()
        w.writerows(vals)

    # ---- corpus-schema rows -------------------------------------------
    # One Experiment_ID per pool; Experiment_value is the per-experiment
    # z-score of log10, population sd, which is the section-3 convention.
    rows = []
    for pool, expid in (("in_vitro", EXP_VITRO), ("in_vivo", EXP_VIVO)):
        sub = [v for v in vals if v["pool"] == pool and v["value"]]
        lg = np.array([float(v["log10"]) for v in sub])
        z = (lg - lg.mean()) / lg.std(ddof=0)
        for v, zz in zip(sub, z):
            v["_z"] = float(zz)
            rows.append((pool, expid, v))

    # ---- the quality columns, through the ONE writer -------------------
    items, records = [], {}
    for n, (pool, expid, v) in enumerate(rows):
        key = f"{pool}|{v['lipid_name']}|{v['model_target']}|{v['dose']}"
        items.append(key)
        records[key] = {
            "lipid": v["lipid_name"], "pool": pool,
            "status": v["status"], "bound": v["value_is_bound"],
            "structure_evidence": struct[v["lipid_name"]]["route"],
            "iupac": struct[v["lipid_name"]]["iupac_from_SI"],
        }

    def confidence(_k, rec):
        # Derived from the item's OWN record, never a literal (error class 29).
        if rec["bound"] == "upper":
            return "low"
        if not rec["iupac"]:
            # lipid 14: no IUPAC name in this SI, structure read off the
            # drawn scheme and corroborated against ZT_2026 in the corpus.
            return "medium"
        return "medium"

    def note(_k, rec):
        bits = []
        if rec["bound"] == "upper":
            bits.append("no bar is drawn for this lipid in its panel, so the "
                        "value is at or below the axis floor; stored at the "
                        "floor as an UPPER BOUND, not a measurement")
        if not rec["iupac"]:
            bits.append("structure has no IUPAC name in this SI (lipid 14 is "
                        "from the authors' ref 41); it is read off the drawn "
                        "scheme Supp Fig 15 and corroborated against ZT_2026's "
                        "Lipid 14-CB1/2/3 in the corpus")
        bits.append("value digitised from a 700 px publisher raster (the only "
                    "one reachable without a browser); 1 px is "
                    "0.02 decades in Fig 2A and 0.04 in Fig 2B/2C")
        return "; ".join(bits)

    q = emit.emit_rows(items, records, confidence=confidence, note=note,
                       value_is_bound=lambda _k, rec: rec["bound"],
                       value_normalization_unit=lambda _k, _r: "experiment",
                       label="peer_39912611")

    out = []
    for (pool, expid, v), qq in zip(rows, q):
        s = struct[v["lipid_name"]]
        d = dict(COMMON)
        d.update({
            "Experiment_ID": expid,
            "Formulation_ID": f"{expid}_F1",
            "IL_name": v["lipid_name"],
            "IL_SMILES": s["SMILES"],
            "IL_linker_name": s["linker"],
            "IL_head_name": s["head"],
            "Model": pool,
            "Model_type": v["model_type"],
            "Model_target": v["model_target"],
            "Route_of_administration": v["route"],
            "Dose_ug_nucleicacid": v["dose"] if pool == "in_vivo" else "",
            "Experiment_method": ("luminescence_normalized" if pool == "in_vitro"
                                  else "radiance_normalized"),
            "Experiment_batching": "individual",
            "Experiment_value": round(v["_z"], 9),
            "Experiment_value_within_model": "",
            "apparent_pKa_of_LNP": PKA[v["lipid_name"]],
            "raw_value": v["value"],
            "raw_value_unit": v["value_unit"],
        })
        d.update(qq)
        out.append(d)

    cols = sorted({k for r in out for k in r})
    with open(os.path.join(HERE, "06_rows.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(out)

    print(json.dumps({
        "values_written": len(vals),
        "rows_written": len(out),
        "in_vitro_rows": sum(1 for r in out if r["Model"] == "in_vitro"),
        "in_vivo_rows": sum(1 for r in out if r["Model"] == "in_vivo"),
        "rows_with_a_caveat": emit.caveat_count(out),
        "upper_bounds": sum(1 for r in out if r["value_is_bound"] == "upper"),
        "z_check_vitro_mean_sd": [
            round(float(np.mean([r["Experiment_value"] for r in out
                                 if r["Model"] == "in_vitro"])), 12),
            round(float(np.std([r["Experiment_value"] for r in out
                                if r["Model"] == "in_vitro"], ddof=0)), 9)],
        "z_check_vivo_mean_sd": [
            round(float(np.mean([r["Experiment_value"] for r in out
                                 if r["Model"] == "in_vivo"])), 12),
            round(float(np.std([r["Experiment_value"] for r in out
                                if r["Model"] == "in_vivo"], ddof=0)), 9)],
    }, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
