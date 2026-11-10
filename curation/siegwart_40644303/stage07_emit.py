#!/usr/bin/env python3
"""Emit the corpus rows in the 79-column schema. NOT merged -- see DECISIONS.md.

    python3 -I stage07_emit.py

Writes 09_rows_vitro.csv and 09_rows_vivo.csv. `build_dataset.py` is NOT
pointed at this folder and must not be while DECISIONS.md carries an OPEN
entry: the merge gate refuses it, and it should.

THE QUALITY COLUMNS ARE WRITTEN BY curation/lib/emit.emit_rows AND NOTHING
ELSE. It refuses a literal confidence (error class 29, found three times:
F-011, F-025, F-046) and refuses a caveat with an empty note. Every value here
comes from the item's own enumeration or digitisation record.

NORMALISATION. The corpus convention is a per-experiment z-score of log10,
population sd. For a DISCRETIZED source the quantity z-scored is the BIN
ORDINAL, not a midpoint (CLAUDE.md section 3, verified there against SL_2021),
so the in-vitro rows z-score the 1/2/3 ordinal. The in-vivo rows z-score
log10(radiance).

THE IN-VIVO ROWS ARE SPLIT INTO TWO EXPERIMENTS, PROVISIONALLY. Figure 2b tops
out at 1e6 and Figure 3b at 6e7; they share no lipid and no benchmark control,
so the ~60x offset between them is not recoverable from anything in the paper.
Pooling them into one z would delete an offset we cannot measure, which is
error class 22 and is what F-016 cost on 5,040 rows. Splitting them is also a
choice, so it is marked provisional on every row and queued as D-1. Whichever
way D-1 goes, every in-vivo Experiment_value here changes.
"""
import csv
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))

import emit   # noqa: E402

ENUM = os.path.join(HERE, "03_enumerated.csv")
VALUES = os.path.join(HERE, "05_values.csv")
SI = os.path.join(HERE, "02_si_table_S3.csv")
OUT_VITRO = os.path.join(HERE, "09_rows_vitro.csv")
OUT_VIVO = os.path.join(HERE, "09_rows_vivo.csv")

HEADER = [
    "Index", "LNP_ID", "Experiment_ID", "Formulation_ID", "IL_name",
    "IL_SMILES", "IL_protonated_SMILES", "IL_head_name", "IL_head_SMILES",
    "IL_linker_name", "IL_linker_SMILES", "IL_tail1_name", "IL_tail1_SMILES",
    "IL_tail2_name", "IL_tail2_SMILES", "IL_tail3_name", "IL_tail3_SMILES",
    "IL_tail4_name", "IL_tail4_SMILES", "IL_molratio",
    "IL_to_nucleicacid_massratio", "IL_to_nucleicacid_chargeratio",
    "HL_name", "HL_SMILES", "HL_molratio", "CHL_name", "CHL_SMILES",
    "CHL_molratio", "PEG_name", "PEG_SMILES", "PEG_molratio",
    "fifthcomponent_name", "fifthcomponent_SMILES", "fifthcomponent_molratio",
    "Aqueous_buffer", "Dialysis_buffer", "Mixing_method", "Model",
    "Model_type", "Model_target", "Route_of_administration", "Cargo",
    "Cargo_type", "Dose_ug_nucleicacid", "Experiment_method",
    "Experiment_batching", "Experiment_value", "Publication_link",
    "Publication_DOI", "Publication_PMID",
] + ["Heavy.Atoms", "Rings", "Aromatic.Rings", "Rotatable.Bonds",
     "ExactMolWt", "Topological.Polar.Surface.Area", "Hydrogen.Bond.Donors",
     "Hydrogen.Bond.Acceptors", "LogP", "Molar.Refractivity",
     "Fraction.sp3.Carbons", "sp3.Carbons", "Nitrogen.Count",
     "Molecular.Weight", "has_ester", "has_carbonate", "has_disulfide",
     "Num_carbon_in_tail", "num_protonatable_nitrogens",
     "num_permanent_cationic_N", "formal_net_charge",
     "num_unsaturated_cc_bonds"] + emit.columns()

# Methods, "In Vitro/In Vivo Nanoparticle Formulations": ionizable lipid, DOPE,
# cholesterol and DMG-PEG2k at 15:15:30:3 mol ratio = 23.81 / 23.81 / 47.62 /
# 4.76 mol%. Luc mRNA in 10 mM citrate pH 4.4, combined 3:1 aqueous:organic.
FORM = {
    "IL_molratio": "23.81", "HL_name": "DOPE", "HL_molratio": "23.81",
    "CHL_name": "Cholesterol", "CHL_molratio": "47.62",
    "PEG_name": "DMG-PEG2000", "PEG_molratio": "4.76",
    "Aqueous_buffer": "citrate", "Mixing_method": "handmixed",
}
PUB = {"Publication_link": "https://doi.org/10.1021/acsbiomaterials.5c00463",
       "Publication_DOI": "10.1021/acsbiomaterials.5c00463",
       "Publication_PMID": "40644303"}


def zscore(vals):
    """Population sd (ddof=0), the corpus convention."""
    m = sum(vals) / len(vals)
    sd = math.sqrt(sum((v - m) ** 2 for v in vals) / len(vals))
    if sd == 0:
        raise SystemExit("stage07: zero variance; cannot z-score.")
    return [(v - m) / sd for v in vals]


def main():
    enum = {r["lipid_name"]: r for r in csv.DictReader(open(ENUM))}
    si = {r["lipid"]: r for r in csv.DictReader(open(SI))}
    vals = [r for r in csv.DictReader(open(VALUES)) if r["emitted"] == "yes"]

    vitro = [r for r in vals if r["record"] == "fig1a_mRNA"]
    vivo = {"SL_2025": [r for r in vals if r["record"] == "fig2b"],
            "SL_2025_2": [r for r in vals if r["record"] == "fig3b"]}

    def quality(recs, items, bound_of, extra_note):
        """One call, one writer. emit refuses a literal and a silent caveat."""
        def conf(i, rec):
            return rec["confidence"]

        def note(i, rec):
            parts = [p for p in (rec["note"], extra_note(i, rec)) if p]
            return "; ".join(parts)

        return emit.emit_rows(
            items, recs,
            confidence=conf, note=note,
            positional_ambiguity=lambda i, rec: rec["positional_ambiguity"],
            structure_assumed=lambda i, rec: "0",
            value_is_bound=bound_of,
            value_normalization_unit=lambda i, rec: "experiment",
            label="SL_2025")

    out_vitro, out_vivo = [], []

    # ---------------- in vitro: 465 mRNA ordinals ----------------------
    zs = zscore([float(r["raw_value"]) for r in vitro])
    items = [r["lipid"] for r in vitro]
    recs = {r["lipid"]: enum[r["lipid"]] for r in vitro}
    bound = {r["lipid"]: r["value_is_bound"] for r in vitro}

    def vitro_note(i, rec):
        v = bound[i]
        s = ("Figure 1a is a three-level heat map; this cell's band is the "
             "open top band '> log2', so the value is a LOWER bound."
             if v == "lower" else
             "value is the ordinal of a three-level band read from Figure 1a, "
             "not a continuous measurement")
        return s

    q = quality(recs, items, lambda i, rec: bound[i], vitro_note)
    for k, r in enumerate(vitro):
        e = enum[r["lipid"]]
        row = dict.fromkeys(HEADER, "")
        row.update(FORM)
        row.update(PUB)
        row.update({
            "Experiment_ID": "SL_2025",
            "Formulation_ID": "SL_2025_Fvitro",
            "IL_name": r["lipid"], "IL_SMILES": e["smiles"],
            "IL_tail1_name": e["tail"], "IL_head_name": e["core"],
            "Model": "in_vitro", "Model_type": "IGROV1",
            "Model_target": "in_vitro",
            "Route_of_administration": "in_vitro",
            "Cargo": "mRNA", "Cargo_type": "FLuc",
            "Dose_ug_nucleicacid": "0.025",
            "Experiment_method": "luminescence_discretized_normalized",
            "Experiment_batching": "individual",
            "Experiment_value": f"{zs[k]:.9f}",
        })
        row.update(q[k])
        out_vitro.append(row)

    # ---------------- in vivo: two panels, two experiments --------------
    for exp, recs_in in vivo.items():
        items = [r["lipid"] for r in recs_in]
        recs = {r["lipid"]: enum[r["lipid"]] for r in recs_in}
        zs = zscore([math.log10(float(r["raw_value"])) for r in recs_in])
        fig = "Figure 2b" if exp == "SL_2025" else "Figure 3b"

        def vivo_note(i, rec, fig=fig):
            s = si.get(i, {})
            return (f"liver radiance digitised from {fig}; PROVISIONAL "
                    "normalisation -- Figures 2b and 3b are drawn on scales "
                    "~60x apart, share no lipid and no benchmark control, so "
                    "they are z-scored separately pending DECISIONS.md D-1 "
                    f"(SI Fig S3: {s.get('diameter_nm','?')} nm, PDI "
                    f"{s.get('PDI','?')}, pKa {s.get('pKa','?')}, E.E. "
                    f"{s.get('encapsulation_efficiency_pct','?')}%)")

        q = quality(recs, items, lambda i, rec: "", vivo_note)
        for k, r in enumerate(recs_in):
            e = enum[r["lipid"]]
            row = dict.fromkeys(HEADER, "")
            row.update(FORM)
            row.update(PUB)
            row.update({
                "Experiment_ID": exp,
                "Formulation_ID": exp + "_Fvivo",
                "IL_name": r["lipid"], "IL_SMILES": e["smiles"],
                "IL_tail1_name": e["tail"], "IL_head_name": e["core"],
                "Dialysis_buffer": "PBS",
                "Model": "in_vivo", "Model_type": "Mouse_B6",
                "Model_target": "liver",
                "Route_of_administration": "intravenous",
                "Cargo": "mRNA", "Cargo_type": "FLuc",
                "Dose_ug_nucleicacid": "5",
                "Experiment_method": "luminescence_normalized",
                "Experiment_batching": "individual",
                "Experiment_value": f"{zs[k]:.9f}",
            })
            row.update(q[k])
            out_vivo.append(row)

    for path, rows in ((OUT_VITRO, out_vitro), (OUT_VIVO, out_vivo)):
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=HEADER)
            w.writeheader()
            w.writerows(rows)

    summary = {
        "vitro_rows": len(out_vitro), "vivo_rows": len(out_vivo),
        "total": len(out_vitro) + len(out_vivo),
        "rows_with_a_caveat": emit.caveat_count(out_vitro + out_vivo),
        "rows_with_value_is_bound": sum(
            1 for r in out_vitro + out_vivo if r["value_is_bound"]),
        "rows_with_positional_ambiguity": sum(
            1 for r in out_vitro + out_vivo
            if str(r["positional_ambiguity"]) == "1"),
        "experiments": sorted({r["Experiment_ID"]
                               for r in out_vitro + out_vivo}),
        "z_check": {
            # grouped by FILE x Experiment_ID, because SL_2025 names both the
            # in-vitro screen and the Figure 2b in-vivo cohort -- the same
            # convention ZT_2026 uses -- and a mean taken across both blocks
            # would be the mean of two different z-scores and would read clean
            # for the wrong reason.
            f"{tag}/{e}": {
                "n": sum(1 for r in rs if r["Experiment_ID"] == e),
                "mean": round(sum(float(r["Experiment_value"]) for r in rs
                                  if r["Experiment_ID"] == e)
                              / sum(1 for r in rs
                                    if r["Experiment_ID"] == e), 12),
                "sd": round((sum((float(r["Experiment_value"])) ** 2
                                 for r in rs if r["Experiment_ID"] == e)
                             / sum(1 for r in rs
                                   if r["Experiment_ID"] == e)) ** 0.5, 9),
            }
            for tag, rs in (("vitro", out_vitro), ("vivo", out_vivo))
            for e in sorted({r["Experiment_ID"] for r in rs})},
        "descriptor_columns": "left blank on purpose: build_dataset.py "
                              "recomputes all 22 from IL_SMILES, and a "
                              "curation that writes its own is how 676 rows "
                              "came to carry descriptors inconsistent with "
                              "their own structure.",
    }
    with open(os.path.join(HERE, "09_emit_qc.json"), "w") as fh:
        json.dump(summary, fh, indent=1)
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
