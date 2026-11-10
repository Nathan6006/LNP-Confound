#!/usr/bin/env python3
"""05_values.csv -- the 49 corpus rows, in the merge schema.

ONLY 49 OF THE 580 LIBRARY MEMBERS HAVE A DELIVERY MEASUREMENT. The other 531
have a TNS fluorescence value (06_tns.csv) and nothing else, and TNS is a
surface-ionization measurement of the formulation, not a delivery outcome -- the
same class as diameter, zeta and uptake, which the corpus's own definition
excludes. So the honest row count for this paper is 49, not 580, and
00_value_sources.json says so in the intake gate rather than at the end.

FORMULATION, from Methods "LNP fabrication" (the 16PE formulation used for the
whole screen): lipidoid : cholesterol : DOPE : PEG2000 = 35 : 46.5 : 16 : 2.5
mol%, lipidoid-to-mRNA 10:1 w/w, 5moU-modified mRNA in 10 mM citrate pH 4,
ethanol/citrate 9:1, pipette-mixed then vortexed, dialysed 1 h against PBS.

ENDPOINT, from Methods "LNP efficacy in murine organs": 6-8 week old C57BL/6
female mice, tail vein, 0.75 mg/kg total mRNA, firefly luciferase; at 3 h
D-luciferin IP, 15 min later sacrifice; liver, spleen, lungs, heart, pancreas
and kidneys excised and imaged on an IVIS. Figure 2B plots the TOTAL.

Run: python3 -I stage06_rows.py   (after stage03 and stage04)
"""
import csv
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))

import emit                                                     # noqa: E402
import reagents                                                 # noqa: E402
from errors import AmbiguityError                               # noqa: E402

EXPERIMENT_ID = "KW_2024"
DOI = "10.1002/adhm.202400225"
PMID = "38888972"
LINK = "https://onlinelibrary.wiley.com/doi/10.1002/adhm.202400225"

# The 67-column LNPDB header plus the 7 quality columns, i.e. the schema
# curation/new_curation.csv uses. The 5 ML-change-log descriptors and all 17
# RDKit descriptor columns are recomputed by build_dataset.py's
# merge_curations(), so they are emitted blank rather than restated here.
COLS = [
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
    "Heavy.Atoms", "Rings", "Aromatic.Rings", "Rotatable.Bonds", "ExactMolWt",
    "Topological.Polar.Surface.Area", "Hydrogen.Bond.Donors",
    "Hydrogen.Bond.Acceptors", "LogP", "Molar.Refractivity",
    "Fraction.sp3.Carbons", "sp3.Carbons", "Nitrogen.Count",
    "Molecular.Weight", "has_ester", "has_carbonate", "has_disulfide",
    "Experiment_value_within_model", "value_normalization_unit", "confidence",
    "positional_ambiguity", "structure_assumed", "value_is_bound",
    "curation_note",
]


def peg():
    """DMPE-PEG2000, and the reason it cannot be resolved by the shared guard.

    Methods: "PEG2000 (14:0 PEG2000 PE) from Avanti Polar Lipids", i.e.
    1,2-dimyristoyl-sn-glycero-3-phosphoethanolamine-N-[methoxy(PEG)-2000] --
    DMPE-PEG2000, which HAS a phosphoethanolamine head.

    `reagents.smi` refuses this name, and the refusal is a real finding rather
    than an obstacle: the corpus carries TWO structures under it. 3,611 rows
    hold a phosphoethanolamine (P and N present); 1,873 `JW_2024` rows hold the
    DMG-PEG2000 diglyceride string, which has no P and no N. That is open
    finding F-049. We take the structure the CHEMICAL NAME requires, from the
    3,611-row majority, and record the collision -- we do not repair F-049 here,
    because a 41-to-1 count is not evidence about which half is wrong
    (the AC_2025 lesson).
    """
    rows = list(csv.DictReader(open(os.path.join(
        HERE, "..", "..", "data", "data_vivo.csv"))))
    rows += list(csv.DictReader(open(os.path.join(
        HERE, "..", "..", "data", "data_vitro.csv"))))
    seen = {}
    for r in rows:
        if r["PEG_name"] == "DMPE-PEG2000" and r["PEG_SMILES"]:
            seen[r["PEG_SMILES"]] = seen.get(r["PEG_SMILES"], 0) + 1
    guard = None
    try:
        reagents.smi("DMPE-PEG2000", rows)
    except (AmbiguityError, SystemExit) as e:
        # corpus_components.smi raises SystemExit, reagents.smi AmbiguityError;
        # both are the same refusal and both are recorded rather than bypassed.
        guard = str(e).split("\n")[0]
    if guard is None:
        raise SystemExit(
            "peg(): the shared guard did NOT refuse 'DMPE-PEG2000'. Either "
            "F-049 has been repaired since this was written -- in which case "
            "take the registry's answer -- or the guard has stopped working. "
            "Refusing either way.")
    cand = [s for s in seen if "P" in s and "N" in s]
    if len(cand) != 1:
        raise SystemExit(f"peg(): {len(cand)} phosphoethanolamine candidates")
    return cand[0], {k: v for k, v in seen.items()}, guard


def main():
    enum = {r["IL_name"]: r for r in csv.DictReader(
        open(os.path.join(HERE, "03_enumerated.csv")))}
    bars = [r for r in csv.DictReader(
        open(os.path.join(HERE, "04_fig2b_bars.csv")))
        if r["is_control"] == "False"]
    if len(bars) != 49:
        raise SystemExit(f"stage06: {len(bars)} library bars, expected 49")
    missing = [r["IL_name"] for r in bars if r["IL_name"] not in enum]
    if missing:
        raise SystemExit(
            f"stage06: {len(missing)} Figure 2B labels are not members of the "
            f"enumerated 580: {missing}. Refusing -- a name on the x axis that "
            "the library does not contain means the label transcription or "
            "the enumeration is wrong, and guessing which would put measured "
            "values on molecules that were never made.")

    dope = reagents.registered("DOPE")
    chol = reagents.registered("Cholesterol")
    peg_smi, peg_seen, peg_guard = peg()
    if not dope or not chol:
        raise SystemExit("stage06: DOPE/Cholesterol are not in the registry")

    # --- the corpus z-score: log10, per experiment, POPULATION sd ------
    logv = np.array([float(r["log10_flux"]) for r in bars])
    z = (logv - logv.mean()) / logv.std(ddof=0)

    names = [r["IL_name"] for r in bars]
    recs = {n: enum[n] for n in names}
    VALUE_NOTE = ("in-vivo total luminescent flux digitised from Figure 2B "
                  "(log axis, 6 decade ticks, tick-fit residual 0.12 px; the "
                  "drawn 10^8 reference line lands at 0.997x nominal)")

    def note(i, rec):
        base = (rec["curation_note"] or "").strip()
        return (base + " " if base else "") + f"Value: {VALUE_NOTE}."

    flags = emit.emit_rows(
        names, recs,
        confidence=lambda i, r: r["confidence"],
        note=note,
        structure_assumed=lambda i, r: r["structure_assumed"],
        positional_ambiguity=lambda i, r: r["positional_ambiguity"],
        value_is_bound=lambda i, r: "",
        value_normalization_unit=lambda i, r: "experiment",
        label=f"{EXPERIMENT_ID} rows")

    out = []
    for bar, zz, fl in zip(bars, z, flags):
        e = enum[bar["IL_name"]]
        row = dict.fromkeys(COLS, "")
        row.update({
            "Experiment_ID": EXPERIMENT_ID,
            "Formulation_ID": f"{EXPERIMENT_ID}_F1",
            "IL_name": e["IL_name"], "IL_SMILES": e["IL_SMILES"],
            "IL_head_name": e["head"], "IL_head_SMILES": e["head_SMILES"],
            "IL_tail1_name": e["tail"],
            "IL_tail1_SMILES": e["tail_acrylate_SMILES"],
            "IL_molratio": "35", "IL_to_nucleicacid_massratio": "10",
            "HL_name": "DOPE", "HL_SMILES": dope, "HL_molratio": "16",
            "CHL_name": "Cholesterol", "CHL_SMILES": chol,
            "CHL_molratio": "46.5",
            "PEG_name": "DMPE-PEG2000", "PEG_SMILES": peg_smi,
            "PEG_molratio": "2.5",
            "Aqueous_buffer": "citrate", "Dialysis_buffer": "PBS",
            "Mixing_method": "handmixed",
            "Model": "in_vivo", "Model_type": "Mouse_B6",
            "Model_target": "whole_body",
            "Route_of_administration": "intravenous",
            "Cargo": "mRNA", "Cargo_type": "FLuc",
            "Dose_ug_nucleicacid": "",
            "Experiment_method": "luminescence_normalized",
            "Experiment_batching": "individual",
            "Experiment_value": repr(float(zz)),
            "Publication_link": LINK, "Publication_DOI": DOI,
            "Publication_PMID": PMID,
        })
        row.update(fl)
        out.append(row)

    with open(os.path.join(HERE, "05_values.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS)
        w.writeheader()
        for r in out:
            w.writerow(r)

    qc = {
        "Experiment_ID": EXPERIMENT_ID,
        "n_rows": len(out),
        "n_library_members_enumerated": len(enum),
        "n_members_with_a_delivery_value": len(out),
        "all_49_are_members_of_the_580": True,
        "z": {"mean": round(float(z.mean()), 12),
              "sd_ddof0": round(float(z.std(ddof=0)), 12),
              "min": round(float(z.min()), 5),
              "max": round(float(z.max()), 5)},
        "raw_flux_p_per_s": {
            "min": float(f"{10 ** logv.min():.6g}"),
            "max": float(f"{10 ** logv.max():.6g}"),
            "geometric_mean": float(f"{10 ** logv.mean():.6g}")},
        "flags": {
            "rows_with_a_caveat": emit.caveat_count(flags),
            "confidence_high": sum(1 for f in flags
                                   if f["confidence"] == "high"),
            "confidence_medium": sum(1 for f in flags
                                     if f["confidence"] == "medium"),
            "confidence_low": sum(1 for f in flags
                                  if f["confidence"] == "low"),
            "structure_assumed": sum(1 for f in flags
                                     if f["structure_assumed"] == "1"),
            "value_is_bound": sum(1 for f in flags if f["value_is_bound"]),
        },
        "PEG": {
            "name_in_paper": "PEG2000 (14:0 PEG2000 PE), Avanti; the same "
                             "Methods also calls it C14-PEG2000 in the TNS "
                             "section",
            "resolved_to": "DMPE-PEG2000",
            "shared_guard_refused": peg_guard,
            "corpus_structures_under_this_name": {
                s[:48] + "...": n for s, n in sorted(
                    peg_seen.items(), key=lambda kv: -kv[1])},
            "_why": "F-049. Two structures share the name in the corpus; the "
                    "phosphoethanolamine one is what '14:0 PEG2000 PE' means, "
                    "so that is used. F-049 is NOT repaired here.",
        },
        "unrecorded_metadata": {
            "Dose_ug_nucleicacid": "BLANK. The dose is stated as 0.75 mg/kg "
                                   "and no body weight is given, so micrograms "
                                   "are not derivable. DECISIONS.md D-5.",
            "IL_to_nucleicacid_chargeratio": "not stated",
            "Model_target": "'whole_body' is the nearest corpus term. The "
                            "measurement is the SUM over six EXCISED organs "
                            "imaged ex vivo, not a whole-animal IVIS image.",
        },
    }
    with open(os.path.join(HERE, "05_values_qc.json"), "w") as fh:
        json.dump(qc, fh, indent=1, sort_keys=True)
    print(json.dumps(qc, indent=1, sort_keys=True))


if __name__ == "__main__":
    main()
