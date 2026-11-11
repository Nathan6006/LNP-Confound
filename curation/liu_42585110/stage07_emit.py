#!/usr/bin/env python3
"""Emit KJ_2026 in the live 79-column schema. -> 06_rows.csv

Run from this folder:  python3 -I stage07_emit.py

NOTE: this folder is NOT merged. DECISIONS.md carries OPEN entries and
build_dataset.py refuses a folder with one, which is the intended behaviour --
D-1 (an All-Rights-Reserved source against CLAUDE.md section 8) and D-2 (what
Figure 1a's unitless -1..+1 scale actually is) are questions this curation is
not entitled to answer. The rows are emitted anyway so the work is reviewable
and so the decisions are made against real numbers.

THE HEADER IS READ LIVE from data/data_vitro.csv. It is not restated here: the
schema has moved from 66 to 72 to 79 columns, and a restated header is how a
curation's flags died one stage before the merge (F-001).

ONE Z PER EXPERIMENT, POOLING ORGANS. The 82 in-vivo values are standardised
together across liver AND spleen, not per organ. A per-organ z is error class
22 and would delete exactly the spleen-versus-liver difference this paper is
about -- the same mistake that forced WT_2026's published selectivity column to
identically zero (F-016/F-003). The organ lives in Model_target, as the
convention requires.

THE TWO ARMS SHARE ONE Experiment_ID AND ARE STANDARDISED SEPARATELY, because
they live in different corpus files on different scales: data_vitro.csv carries
the DC2.4 screen and data_vivo.csv the organ radiances. ZT_2026 and XH_2024 are
both already in the corpus this way.
"""
import csv
import json
import os
import sys

import numpy as np

LIB = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib")
sys.path.insert(0, LIB)
import emit                                                      # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "provenance"))
sys.path.insert(0, os.path.join(ROOT, "curation"))
import reagent_registry as RR                                    # noqa: E402
from corpus_components import load as cc_load, smi as cc_smi     # noqa: E402

sys.path.insert(0, HERE)
from stage05_formulation import CONTROL_MOLPCT, rail_molpct, META  # noqa: E402

OUT = os.path.join(HERE, "06_rows.csv")
EXP = "KJ_2026"
PMID = "42585110"
DOI = "10.1002/anie.3823572"
LINK = "https://doi.org/10.1002/anie.3823572"
BENCH = "DLin-MC3-DMA"      # the corpus/registry spelling; the paper writes
                            # "Dlin-MC3", a THIRD spelling that matches 0
                            # corpus rows by name while the molecule is
                            # present 3,098 times (02_corpus_overlap.json)

NOTE_VITRO = (
    "Digitised from Fig 1a, a continuous red-white-blue heat map on an "
    "UNITLESS, unlabelled -1.0 to +1.0 colour scale; the paper never states "
    "what the scale is a transform of. Experiment_value is the per-experiment "
    "z of that scaled quantity, which equals the corpus convention (z of "
    "log10 signal) if and only if the plotted quantity is log-scaled. See "
    "curation/liu_42585110/DECISIONS.md D-2."
)
NOTE_VIVO = (
    "Digitised from Fig 1c (Avg Radiance, x1e8 p/s/cm2/sr); one z per "
    "experiment over liver and spleen together, organ in Model_target."
)
NOTE_BOUND = (
    "UPPER BOUND, not a measurement: the Fig 1c bar's entire body lies within "
    "4 px of the drawn baseline and is not separable from the axis rule's own "
    "thickness, so the true Avg Radiance is at or below 0.036e8 p/s/cm2/sr. "
) + NOTE_VIVO


def header():
    with open(os.path.join(ROOT, "data", "data_vitro.csv")) as fh:
        return next(csv.reader(fh))


def components():
    corpus = list(csv.DictReader(open(os.path.join(ROOT, "data",
                                                   "data_vitro.csv"))))
    comp, amb = cc_load(corpus)
    return {
        # DSPC and cholesterol come from provenance/reagent_registry.py, which
        # states the structure FROM THE CHEMICAL NAME and so does not inherit
        # whatever the corpus happens to hold.
        "HL_name": "DSPC", "HL_SMILES": RR.REGISTRY["DSPC"][0],
        "CHL_name": "Cholesterol", "CHL_SMILES": RR.REGISTRY["Cholesterol"][0],
        # The registry deliberately holds NO PEG-lipid: a PEG block length is
        # a distribution, so no single SMILES is right. This one comes from
        # the corpus through the guarded loader, which RAISES on an ambiguous
        # name. The name match is not a coincidence of spelling: the SI says
        # "methoxy-PEG2000-carbamate-dimyristoyl glycerol", and the corpus's
        # DMG-C-PEG2000 is C119H235NO49 with two myristoyl esters and one
        # carbamate -- which is that compound. The corpus's DMG-PEG2000, by
        # contrast, is the ether-linked C122H242O50 and would be the wrong
        # molecule; this is the C14-PEG2000 trap named in EXTRACTION_BRIEF 5.
        "PEG_name": "DMG-C-PEG2000",
        "PEG_SMILES": cc_smi(comp, amb, "DMG-C-PEG2000"),
    }


def main():
    cols = header()
    if len(cols) != 79:
        raise SystemExit(f"stage07: the live header has {len(cols)} columns.")
    comp = components()
    rail = rail_molpct()
    enum = {r["IL_name"]: r for r in csv.DictReader(
        open(os.path.join(HERE, "03_enumerated.csv")))}
    a = list(csv.DictReader(open(os.path.join(HERE, "04_fig1a_cells.csv"))))
    c = list(csv.DictReader(open(os.path.join(HERE, "04_fig1c_bars.csv"))))

    def structure(name):
        if name == BENCH:
            return RR.REGISTRY["DLin-MC3-DMA"][0]
        return enum[name]["IL_SMILES"]

    def molpct(name):
        return CONTROL_MOLPCT if name == BENCH else rail

    # ---------------- in vitro: 41 rows ----------------------------------
    means = {}
    for r in a:
        means.setdefault(r["IL_name"], []).append(float(r["scaled_value"]))
    order = [r["IL_name"] for r in a]
    order = list(dict.fromkeys(order))
    names = [BENCH if n == "Dlin-MC3" else n for n in order]
    vals = np.array([np.mean(means[n]) for n in order])
    z_vitro = (vals - vals.mean()) / vals.std(ddof=0)

    # ---------------- in vivo: 82 rows -----------------------------------
    vv = [r for r in c]
    raw = np.array([float(r["value_avg_radiance_e8"]) for r in vv])
    lg = np.log10(raw)
    z_vivo = (lg - lg.mean()) / lg.std(ddof=0)

    # ---------------- flags, through emit.emit_rows ----------------------
    recs = {}
    for n in names:
        recs[("vitro", n)] = enum.get(n, {"IL_name": n, "confidence": "high",
                                          "confidence_reason":
                                          "benchmark; structure from "
                                          "provenance/reagent_registry.py"})
    for i, r in enumerate(vv):
        n = BENCH if r["IL_name"] == "Dlin-MC3" else r["IL_name"]
        recs[("vivo", i)] = dict(enum.get(n, {"IL_name": n,
                                              "confidence": "high"}),
                                 **{"organ": r["organ"],
                                    "bound": r["value_is_bound"]})

    vitro_items = [("vitro", n) for n in names]
    vivo_items = [("vivo", i) for i in range(len(vv))]

    def conf(item, rec):
        return rec.get("confidence", "")

    def note_vitro(item, rec):
        base = NOTE_VITRO
        if rec.get("confidence") == "medium":
            base = ("Structure has NO ESI-MS in the SI (the only one of the "
                    "40 without one), so its molecular formula is "
                    "unconfirmed. ") + base
        return base

    def note_vivo(item, rec):
        base = NOTE_BOUND if rec.get("bound") == "upper" else NOTE_VIVO
        if rec.get("confidence") == "medium":
            base = ("Structure has NO ESI-MS in the SI (the only one of the "
                    "40 without one), so its molecular formula is "
                    "unconfirmed. ") + base
        return base

    q_vitro = emit.emit_rows(
        vitro_items, recs, confidence=conf, note=note_vitro,
        positional_ambiguity=lambda i, r: "0",
        structure_assumed=lambda i, r: "0",
        value_is_bound=lambda i, r: "",
        value_normalization_unit=lambda i, r: "experiment",
        within_model=lambda i, r: "", label="KJ_2026 vitro")
    q_vivo = emit.emit_rows(
        vivo_items, recs, confidence=conf, note=note_vivo,
        positional_ambiguity=lambda i, r: "0",
        structure_assumed=lambda i, r: "0",
        value_is_bound=lambda i, r: r.get("bound", ""),
        value_normalization_unit=lambda i, r: "experiment",
        within_model=lambda i, r: "", label="KJ_2026 vivo")

    rows = []
    for k, n in enumerate(names):
        mp = molpct(n)
        row = {col: "" for col in cols}
        row.update({
            "Experiment_ID": EXP,
            "Formulation_ID": f"{EXP}_{'MC3control' if n == BENCH else 'RAIL'}",
            "IL_name": n, "IL_SMILES": structure(n),
            "IL_molratio": f"{mp['IL']:.4g}",
            "HL_name": comp["HL_name"], "HL_SMILES": comp["HL_SMILES"],
            "HL_molratio": f"{mp['HL']:.4g}",
            "CHL_name": comp["CHL_name"], "CHL_SMILES": comp["CHL_SMILES"],
            "CHL_molratio": f"{mp['CHL']:.5g}",
            "PEG_name": comp["PEG_name"], "PEG_SMILES": comp["PEG_SMILES"],
            "PEG_molratio": f"{mp['PEG']:.4g}",
            "Aqueous_buffer": META["Aqueous_buffer"],
            "Dialysis_buffer": META["Dialysis_buffer"],
            "Mixing_method": META["Mixing_method"],
            "Model": "in_vitro", "Model_type": "DC2.4",
            "Model_target": "in_vitro", "Route_of_administration": "in_vitro",
            "Cargo": "mRNA", "Cargo_type": "FLuc",
            "Dose_ug_nucleicacid": "",
            "Experiment_method": "luminescence_normalized",
            "Experiment_batching": "individual",
            "Experiment_value": repr(float(z_vitro[k])),
            "Publication_link": LINK, "Publication_DOI": DOI,
            "Publication_PMID": PMID,
        })
        q = dict(q_vitro[k])
        q["Experiment_value_within_model"] = repr(float(z_vitro[k]))
        q["curation_note"] = (
            q["curation_note"] + " Dose: 0.5 ug/mL mRNA for 6 h; the well "
            "volume is never stated, so the absolute microgram dose this "
            "column wants is not recoverable and is left blank.")
        row.update(q)
        rows.append(row)

    for i, r in enumerate(vv):
        n = BENCH if r["IL_name"] == "Dlin-MC3" else r["IL_name"]
        mp = molpct(n)
        row = {col: "" for col in cols}
        row.update({
            "Experiment_ID": EXP,
            "Formulation_ID": f"{EXP}_{'MC3control' if n == BENCH else 'RAIL'}",
            "IL_name": n, "IL_SMILES": structure(n),
            "IL_molratio": f"{mp['IL']:.4g}",
            "HL_name": comp["HL_name"], "HL_SMILES": comp["HL_SMILES"],
            "HL_molratio": f"{mp['HL']:.4g}",
            "CHL_name": comp["CHL_name"], "CHL_SMILES": comp["CHL_SMILES"],
            "CHL_molratio": f"{mp['CHL']:.5g}",
            "PEG_name": comp["PEG_name"], "PEG_SMILES": comp["PEG_SMILES"],
            "PEG_molratio": f"{mp['PEG']:.4g}",
            "Aqueous_buffer": META["Aqueous_buffer"],
            "Dialysis_buffer": META["Dialysis_buffer"],
            "Mixing_method": META["Mixing_method"],
            "Model": "in_vivo", "Model_type": "BALB/c",
            "Model_target": r["organ"],
            "Route_of_administration": "intravenous",
            "Cargo": "mRNA", "Cargo_type": "FLuc",
            "Dose_ug_nucleicacid": "",
            "Experiment_method": "luminescence_normalized",
            "Experiment_batching": "individual",
            "Experiment_value": repr(float(z_vivo[i])),
            "Publication_link": LINK, "Publication_DOI": DOI,
            "Publication_PMID": PMID,
        })
        q = dict(q_vivo[i])
        q["Experiment_value_within_model"] = repr(float(z_vivo[i]))
        q["curation_note"] = (
            q["curation_note"] + " Dose: 0.1 mg/kg mRNA IV; this column is in "
            "micrograms and the animal mass is never stated, so the absolute "
            "dose is not recoverable and is left blank (F-004's posture).")
        row.update(q)
        rows.append(row)

    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)

    nv = len(names)
    print(f"wrote {len(rows)} rows to {OUT}")
    print(f"  in vitro {nv} (DC2.4)   z mean {z_vitro.mean():+.3e} "
          f"sd {z_vitro.std(ddof=0):.12f}")
    print(f"  in vivo  {len(vv)} (liver+spleen, ONE z)  z mean "
          f"{z_vivo.mean():+.3e} sd {z_vivo.std(ddof=0):.12f}")
    print(f"  distinct IL structures emitted: "
          f"{len({r['IL_SMILES'] for r in rows})} "
          f"(40 RAILs + the Dlin-MC3 benchmark)")
    print(f"  rows_with_a_caveat (emit.caveat_count): "
          f"{emit.caveat_count(q_vitro + q_vivo)}")
    for q in (q_vitro, q_vivo):
        pass
    print(f"  value_is_bound = upper: "
          f"{sum(1 for r in rows if r['value_is_bound'] == 'upper')}")
    print(f"  confidence = medium: "
          f"{sum(1 for r in rows if r['confidence'] == 'medium')} "
          "(AT3-Bex, in vitro and both organs)")


if __name__ == "__main__":
    main()
