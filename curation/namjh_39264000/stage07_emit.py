#!/usr/bin/env python3
"""Stage 07 -- emit rows in the corpus's 74-column curation schema.

    python3 -I stage07_emit.py      ->  06_rows.csv

TWO EXPERIMENTS, NOT ONE. The brief for a vaccine paper is right to insist on
this. The per-lipid SCREEN readout and the VACCINE readout are different
experiments here, and a third set is not emitted at all:

  SB_2025    72 rows. Renilla luciferase in the harvested EAR, 6 h after an
             intradermal injection of 5 ug RLuc mRNA, ICR mouse. Fig S2c (LNP
             type 1) and Fig S2e (LNP type 4).
  SB_2025_2  36 rows. Human EPO in SERUM by ELISA, 6 h after an intramuscular
             injection of 10 ug hEPO mRNA, ICR mouse, LNP type 4 only. Fig 4b.
  NOT EMITTED  total IgG, IgG1/IgG2a, IFN-gamma ELISpot, TNF-alpha, IL-2
             (Fig 4f, 4g, Fig S4). These are the actual VACCINE readouts and
             they are a third experiment: different cargo (HPV E6/E7 mRNA, not
             a reporter), different schedule (prime/boost a week apart, read a
             week after the boost, against a single injection read at 6 h),
             different mouse strain, and units that are antibody titre and
             spot-forming cells. Nothing is converted and nothing is pooled.
             See DECISIONS.md D-3.

UNITS ARE AS PRINTED. Both emitted readouts are FOLD CHANGE versus the Nil
group of the experiment in question (Fig S2 caption, Fig 4 caption). The
corpus convention is then applied -- `Experiment_value` is a per-experiment
z-score of log10, population sd -- and nothing else is converted.

THE TWO LNP TYPES ARE Z-SCORED TOGETHER, ON PURPOSE. F-016 is 5,040 rows whose
between-assay offset was deleted by standardising each sub-block separately.
Here the type-1 and type-4 panels are already on a common scale by
construction, both being fold change over their own Nil control, and they
differ by a factor of two at the top of the range. Pooling preserves that;
standardising each panel to mean 0 would delete it. `value_normalization_unit`
says `experiment` on every row.

FLAGS GO ON THE ROW, through `emit.emit_rows`, which refuses a literal
`confidence` (error class 29, found three times) and refuses a caveat with an
empty note. The confidence of a member is the WEAKEST of its own blocks, which
is computed in stage03 from 01_blocks.csv and read back here -- never written
as a constant.
"""
import csv
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import emit                                                 # noqa: E402
import reagents                                             # noqa: E402

DOI = "10.1002/smll.202405618"
PMID = "39264000"
LINK = "https://onlinelibrary.wiley.com/doi/10.1002/smll.202405618"

# Table S2 of the SI, as molar PERCENT, which is how the corpus stores
# IL/HL/CHL/PEG molratio (SB_2024 carries 50.0 / 10.0 / 38.5 / 1.5).
LNP_TYPE = {
    1: {"IL_molratio": 50.0, "HL_name": "DSPC", "HL_molratio": 10.0,
        "CHL_molratio": 38.5, "PEG_molratio": 1.5, "NP": 12.0},
    4: {"IL_molratio": 25.0, "HL_name": "DOPE", "HL_molratio": 25.0,
        "CHL_molratio": 48.5, "PEG_molratio": 1.5, "NP": 6.0},
}
# The excipient structures are taken from provenance/reagent_registry.json
# THROUGH curation/lib `reagents`, which REFUSES an ambiguous name. They are
# not written out here: error class 30 is an unguarded name -> structure
# lookup, it was "repaired" in round 1 and was still live in three emitters
# when F-044 found it, and 292 corpus names carry more than one structure.
# DMG-PEG2000 is deliberately absent from the registry -- a PEG block length
# is a distribution, so no single SMILES is right -- and PEG_SMILES stays
# blank for that reason, not by oversight.
CHL_NAME = "Cholesterol"
PEG_NAME = "DMG-PEG2000"

EXPERIMENTS = {
    "Fig S2c": {"Experiment_ID": "SB_2025", "lnp_type": 1,
                "Model_target": "ear", "Route": "intradermal",
                "Cargo_type": "RLuc", "Dose": 5.0,
                "Experiment_method": "luminescence_normalized"},
    "Fig S2e": {"Experiment_ID": "SB_2025", "lnp_type": 4,
                "Model_target": "ear", "Route": "intradermal",
                "Cargo_type": "RLuc", "Dose": 5.0,
                "Experiment_method": "luminescence_normalized"},
    "Fig 4b": {"Experiment_ID": "SB_2025_2", "lnp_type": 4,
               "Model_target": "serum", "Route": "intramuscular",
               "Cargo_type": "hEPO", "Dose": 10.0,
               "Experiment_method": "protein_adundance_normalized"},
}

COLUMNS = [
    "Index", "LNP_ID", "Experiment_ID", "Formulation_ID", "IL_name",
    "IL_SMILES", "IL_protonated_SMILES", "IL_head_name", "IL_head_SMILES",
    "IL_linker_name", "IL_linker_SMILES", "IL_tail1_name", "IL_tail1_SMILES",
    "IL_tail2_name", "IL_tail2_SMILES", "IL_tail3_name", "IL_tail3_SMILES",
    "IL_tail4_name", "IL_tail4_SMILES", "IL_molratio",
    "IL_to_nucleicacid_massratio", "IL_to_nucleicacid_chargeratio", "HL_name",
    "HL_SMILES", "HL_molratio", "CHL_name", "CHL_SMILES", "CHL_molratio",
    "PEG_name", "PEG_SMILES", "PEG_molratio", "fifthcomponent_name",
    "fifthcomponent_SMILES", "fifthcomponent_molratio", "Aqueous_buffer",
    "Dialysis_buffer", "Mixing_method", "Model", "Model_type", "Model_target",
    "Route_of_administration", "Cargo", "Cargo_type", "Dose_ug_nucleicacid",
    "Experiment_method", "Experiment_batching", "Experiment_value",
    "Publication_link", "Publication_DOI", "Publication_PMID", "Heavy.Atoms",
    "Rings", "Aromatic.Rings", "Rotatable.Bonds", "ExactMolWt",
    "Topological.Polar.Surface.Area", "Hydrogen.Bond.Donors",
    "Hydrogen.Bond.Acceptors", "LogP", "Molar.Refractivity",
    "Fraction.sp3.Carbons", "sp3.Carbons", "Nitrogen.Count",
    "Molecular.Weight", "has_ester", "has_carbonate", "has_disulfide",
] + emit.columns()


def load(name):
    with open(os.path.join(HERE, name)) as fh:
        return list(csv.DictReader(fh))


def descriptors(smiles):
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors
    RDLogger.DisableLog("rdApp.*")
    m = Chem.MolFromSmiles(smiles)
    return {
        "Heavy.Atoms": m.GetNumHeavyAtoms(),
        "Rings": rdMolDescriptors.CalcNumRings(m),
        "Aromatic.Rings": rdMolDescriptors.CalcNumAromaticRings(m),
        "Rotatable.Bonds": rdMolDescriptors.CalcNumRotatableBonds(m),
        "ExactMolWt": round(Descriptors.ExactMolWt(m), 6),
        "Topological.Polar.Surface.Area": round(rdMolDescriptors.CalcTPSA(m), 4),
        "Hydrogen.Bond.Donors": rdMolDescriptors.CalcNumHBD(m),
        "Hydrogen.Bond.Acceptors": rdMolDescriptors.CalcNumHBA(m),
        "LogP": round(Crippen.MolLogP(m), 4),
        "Molar.Refractivity": round(Crippen.MolMR(m), 4),
        "Fraction.sp3.Carbons": round(rdMolDescriptors.CalcFractionCSP3(m), 6),
        "sp3.Carbons": sum(1 for a in m.GetAtoms()
                           if a.GetSymbol() == "C"
                           and a.GetHybridization().name == "SP3"),
        "Nitrogen.Count": sum(1 for a in m.GetAtoms()
                              if a.GetSymbol() == "N"),
        "Molecular.Weight": round(Descriptors.MolWt(m), 4),
        "has_ester": bool(m.HasSubstructMatch(
            Chem.MolFromSmarts("[CX3](=O)[OX2H0][#6]"))),
        "has_carbonate": bool(m.HasSubstructMatch(
            Chem.MolFromSmarts("[OX2][CX3](=O)[OX2]"))),
        "has_disulfide": bool(m.HasSubstructMatch(
            Chem.MolFromSmarts("[#16X2][#16X2]"))),
    }


def main():
    blocks = {r["block_id"]: r for r in load("01_blocks.csv")}
    enum = {r["IL_name"]: r for r in load("03_enumerated.csv")}
    values = load("05_values.csv")

    # ---- the enumeration record every flag must be traceable to ---------
    records = {}
    for v in values:
        key = (v["panel"], v["IL_name"])
        e = enum[v["IL_name"]]
        records[key] = {"enum": e, "value": v,
                        "head": blocks[e["FH"]], "bridge": blocks[e["EB"]],
                        "tail": blocks[e["T"]]}

    def confidence(key, rec):
        # The member's own confidence, computed in stage03 as the weakest of
        # its three blocks. Read back, never written here: emit.emit_rows
        # REFUSES a literal, because a hard-coded "high" once relabelled 24
        # lipids whose record reads "ASSUMED..." as the project's highest
        # confidence (F-011), and the same class was found twice more.
        return rec["enum"]["confidence"]

    def note(key, rec):
        e, v = rec["enum"], rec["value"]
        bits = []
        weak = [b["block_id"] for b in (rec["head"], rec["bridge"],
                                        rec["tail"])
                if b["confidence"] != "high"]
        if weak:
            bits.append(
                "Block(s) " + "/".join(weak) + " are drawn only in Fig 3a and "
                "are named in no synthetic procedure, no HRMS and no NMR "
                "assignment; they were measured from that drawing by a tracer "
                "calibrated on the blocks the SI does name "
                "(2-hexyl-1-decanol, 2-octyl-1-dodecanol, glutaric anhydride, "
                "3,5-dihydroxybenzoic acid)")
        if int(e["n_unassigned_stereocentres"]):
            bits.append(
                f"{e['n_unassigned_stereocentres']} sp3 stereocentres are "
                "stored unassigned because the paper assigns none: the "
                "glycerol centre comes from racemic 3-iodopropane-1,2-diol "
                "and the tail branch carbons from racemic branched alcohols")
        if v["value_is_bound"] == "upper":
            bits.append(
                f"the value is an UPPER BOUND: this cell matches the "
                f"colourbar's lightest swatch, so any value in "
                f"[{v['band_lo']}, {v['band_hi']}] fold renders identically "
                "(Fig 4c draws this curation's lowest cell at 12.1 fold "
                "against a stored 53.1, inside the band)")
        elif v["value_is_bound"] == "lower":
            bits.append("the value is a LOWER BOUND: this cell matches the "
                        "colourbar's darkest swatch")
        bits.append(
            f"value digitised from {v['panel']}'s colour ramp; the raw "
            f"reading is {v['value']} fold change vs the Nil group, band "
            f"[{v['band_lo']}, {v['band_hi']}]")
        return ". ".join(bits) + "."

    items = list(records)
    q = emit.emit_rows(
        items, records,
        confidence=confidence, note=note,
        positional_ambiguity=lambda k, r: "0",
        structure_assumed=lambda k, r: "0",
        value_is_bound=lambda k, r: r["value"]["value_is_bound"],
        value_normalization_unit=lambda k, r: "experiment",
        within_model=None,
        label="SB_2025")
    quality = dict(zip(items, q))

    # ---- per-experiment z-score of log10, population sd -----------------
    by_exp = {}
    for v in values:
        by_exp.setdefault(EXPERIMENTS[v["panel"]]["Experiment_ID"],
                          []).append(v)
    zmap, zqc = {}, {}
    for exp, vs in by_exp.items():
        logs = [math.log10(float(v["value"])) for v in vs]
        mu = sum(logs) / len(logs)
        sd = math.sqrt(sum((x - mu) ** 2 for x in logs) / len(logs))
        for v, x in zip(vs, logs):
            zmap[(v["panel"], v["IL_name"])] = (x - mu) / sd
        zqc[exp] = {"n": len(vs), "mean_log10": round(mu, 6),
                    "population_sd_log10": round(sd, 6),
                    "panels": sorted(set(v["panel"] for v in vs))}

    desc_cache = {}
    rows = []
    for i, v in enumerate(values, start=1):
        name = v["IL_name"]
        e = enum[name]
        spec = EXPERIMENTS[v["panel"]]
        t = LNP_TYPE[spec["lnp_type"]]
        if e["IL_SMILES"] not in desc_cache:
            desc_cache[e["IL_SMILES"]] = descriptors(e["IL_SMILES"])
        d = desc_cache[e["IL_SMILES"]]
        head, bridge, tail = blocks[e["FH"]], blocks[e["EB"]], blocks[e["T"]]
        r = {c: "" for c in COLUMNS}
        r.update({
            "Index": i, "LNP_ID": f"{spec['Experiment_ID']}_{i:03d}",
            "Experiment_ID": spec["Experiment_ID"],
            "Formulation_ID": f"LNP_type_{spec['lnp_type']}",
            "IL_name": name, "IL_SMILES": e["IL_SMILES"],
            "IL_head_name": "" if e["FH"] == "FH-3" else e["FH"],
            "IL_head_SMILES": head["SMILES"],
            "IL_linker_name": e["EB"], "IL_linker_SMILES": bridge["SMILES"],
            "IL_tail1_name": e["T"], "IL_tail1_SMILES": tail["SMILES"],
            "IL_tail2_name": e["T"], "IL_tail2_SMILES": tail["SMILES"],
            "IL_molratio": t["IL_molratio"],
            "IL_to_nucleicacid_chargeratio": t["NP"],
            "HL_name": t["HL_name"], "HL_SMILES": reagents.registered(t["HL_name"]),
            "HL_molratio": t["HL_molratio"],
            "CHL_name": CHL_NAME,
            "CHL_SMILES": reagents.registered(CHL_NAME),
            "CHL_molratio": t["CHL_molratio"],
            "PEG_name": PEG_NAME, "PEG_molratio": t["PEG_molratio"],
            "Aqueous_buffer": "citrate",
            "Mixing_method": "microfluidics",
            "Model": "in_vivo", "Model_type": "Mouse_ICR",
            "Model_target": spec["Model_target"],
            "Route_of_administration": spec["Route"],
            "Cargo": "mRNA", "Cargo_type": spec["Cargo_type"],
            "Dose_ug_nucleicacid": spec["Dose"],
            "Experiment_method": spec["Experiment_method"],
            "Experiment_batching": "individual",
            "Experiment_value": round(zmap[(v["panel"], name)], 9),
            "Publication_link": LINK, "Publication_DOI": DOI,
            "Publication_PMID": PMID,
        })
        r.update(d)
        r.update(quality[(v["panel"], name)])
        rows.append(r)

    with open(os.path.join(HERE, "06_rows.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS)
        w.writeheader()
        w.writerows(rows)

    qc = {"_what": "What stage07 emitted, and the z-score it applied.",
          "n_rows": len(rows),
          "experiments": zqc,
          "rows_with_a_caveat": emit.caveat_count(
              [quality[k] for k in items]),
          "rows_carrying_any_flag": sum(
              1 for k in items if any(str(quality[k].get(c) or "")
                                      for c in ("confidence",
                                                "positional_ambiguity",
                                                "structure_assumed",
                                                "value_is_bound"))),
          "confidence_breakdown": {
              c: sum(1 for k in items if quality[k]["confidence"] == c)
              for c in ("high", "medium", "low")},
          "value_is_bound_breakdown": {
              b: sum(1 for k in items if quality[k]["value_is_bound"] == b)
              for b in ("upper", "lower", "")},
          "NOT_EMITTED": {
              "the other 177 of the paper's 213 LNPs":
                  "no structure is recoverable -- the 14 functional heads and "
                  "14 EB-T chains of the full screen are drawn nowhere in the "
                  "paper or its SI",
              "total IgG, IgG1/IgG2a, IFN-gamma, TNF-alpha, IL-2":
                  "the vaccine readouts; a third experiment with different "
                  "cargo, schedule, strain and units (DECISIONS.md D-3)",
              "Table S3 zeta potential, Fig 3c encapsulation efficiency, "
              "Fig S3 size and PDI":
                  "physicochemical properties, not delivery measurements; "
                  "the corpus's in-vivo definition excludes diameter, zeta "
                  "and uptake",
              "Fig 6 expression-score ranks and multi-compartmental ratios":
                  "derived scores over a subset, not measurements"},
          }
    with open(os.path.join(HERE, "06_emit_qc.json"), "w") as fh:
        json.dump(qc, fh, indent=1)
    print(f"06_rows.csv: {len(rows)} rows, {len(COLUMNS)} columns")
    for exp, z in zqc.items():
        print(f"  {exp}: n={z['n']} panels={z['panels']}")
    print("  rows_with_a_caveat:", qc["rows_with_a_caveat"],
          " confidence:", qc["confidence_breakdown"],
          " bounds:", qc["value_is_bound_breakdown"])


if __name__ == "__main__":
    main()
