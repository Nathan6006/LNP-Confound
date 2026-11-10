#!/usr/bin/env python3
"""01_blocks / 03_enumerated / 05_values / 06_rows, and the QC numbers.

WHAT THIS PAPER CONTRIBUTES, stated up front because it governs everything
below: THREE distinct ionizable lipids, and all three are already in the
corpus. The 18-member "library" holds the ionizable lipid FIXED at 35 mol%
C12-200 and varies the helper lipid (DOPE/DSPC/DOTAP) and its molar content
(15/20/25/30.5/35/50) -- 3 x 6 = 18 = the row count. Figure 4G is the only
measurement in the paper where the ionizable lipid itself changes.

ONE Z PER EXPERIMENT (brief section 11). There is one cell type, one cargo,
one timepoint and one formulation in the emitted block, so the question of a
per-cell or per-cohort z (error class 22) does not arise -- but the
normalisation unit is still DECLARED on every row rather than inferred.
"""
import csv
import json
import math
import os
import sys

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
sys.path.insert(0, os.path.join(HERE, "..", "..", "provenance"))
import emit                                            # noqa: E402
import reagent_registry as RR                          # noqa: E402

VITRO = os.path.join(HERE, "..", "..", "data", "data_vitro.csv")
EXPERIMENT_ID = "MR_2025"          # first author Marco Tullio RODRIGUES Alves
PMID, DOI = "41138890", "10.1016/j.ijpharm.2025.126278"
TIMEPOINT = "72h"                  # the paper's terminal point, where its own
                                   # statistics and its conclusion sit

# The corpus C12-200 structure, confirmed against the Figure 1A drawing by
# stage05_anchor_check.py at an EXACT FULL InChIKey match on 498 AC_2025 rows.
# It is read from the corpus here rather than typed.
C12_200_INCHIKEY = "RVHYPUORVDKRTM-UHFFFAOYSA-N"


def corpus_rows():
    with open(VITRO, newline="") as fh:
        return list(csv.DictReader(fh))


def corpus_header():
    """The 79-column header, read LIVE from the corpus (brief section 11)."""
    with open(VITRO, newline="") as fh:
        return next(csv.reader(fh))


def main():
    rows = corpus_rows()
    header = corpus_header()
    ik = {}

    def key(s):
        if s not in ik:
            m = Chem.MolFromSmiles(s)
            ik[s] = Chem.MolToInchiKey(m) if m else None
        return ik[s]

    # --- the three ionizable lipids --------------------------------------
    c12 = sorted({r["IL_SMILES"].strip() for r in rows
                  if (r.get("IL_name") or "").strip() == "C12-200"})
    if len(c12) != 1 or key(c12[0]) != C12_200_INCHIKEY:
        sys.exit("C12-200 is not a single pinned structure in the corpus; "
                 "REFUSING rather than choosing one (error class 30).")
    ILS = {
        "C12-200": {
            "smiles": c12[0], "source": "corpus AC_2025 (498 rows, 1 "
            "structure), confirmed against the Fig 1A drawing by an "
            "independently assembled anchor at an EXACT full-InChIKey match",
            "drawn_in_paper": True, "confidence": "high"},
        "SM-102": {
            "smiles": RR.REGISTRY["SM-102"][0],
            "source": "provenance/reagent_registry.py, written from the "
            "chemical name " + RR.REGISTRY["SM-102"][1],
            "drawn_in_paper": False, "confidence": "high"},
        "ALC-0315": {
            "smiles": RR.REGISTRY["ALC-0315"][0],
            "source": "provenance/reagent_registry.py, written from the "
            "chemical name " + RR.REGISTRY["ALC-0315"][1],
            "drawn_in_paper": False, "confidence": "high"},
    }
    for name, d in ILS.items():
        d["inchikey"] = key(d["smiles"])

    # --- the measured values ---------------------------------------------
    fig = json.load(open(os.path.join(HERE, "04_fig4g.json")))
    bar = {(b["group"], b["timepoint"]): b for b in fig["bars"]}
    GROUP_TO_IL = {"LNP4": "C12-200", "SM-102": "SM-102",
                   "ALC-0315": "ALC-0315"}

    values = []
    for grp, il in GROUP_TO_IL.items():
        for tp in ("24h", "48h", "72h"):
            b = bar[(grp, tp)]
            values.append({
                "locus": "Figure 4G", "tier": "5_pixels",
                "group_label_as_printed": grp, "timepoint": tp,
                "ionizable_lipid": il,
                "readout_as_printed":
                    "Fluorescene Intensity / area (arb.units)",
                "units": "arbitrary units, linear axis 0-5000",
                "value": b["value"],
                "cargo": "pDNA", "cargo_type": "GFP",
                "cell": "primary rat cortical neuron",
                "organ": "", "route": "in vitro",
                "dose": "NOT STATED in any source held",
                "assay": "GFP fluorescence intensity per area, imaging",
                "pooled": "no", "pool_size": "",
                "corpus_admissible": "yes" if tp == TIMEPOINT else "no",
                "admissibility_reason": (
                    "protein expression from delivered pDNA -- a translated "
                    "readout, admissible under the brief's section 3 list"
                    if tp == TIMEPOINT else
                    "same lipid and same formulation as the %s row; emitting "
                    "all three timepoints would put three conflicting values "
                    "on one SMILES inside one experiment, the shape of the "
                    "corpus's 48 duplicate-SMILES groups. Measured and "
                    "recorded, not emitted." % TIMEPOINT)})
    for tp in ("24h", "48h", "72h"):
        b = bar[("Control", tp)]
        values.append({
            "locus": "Figure 4G", "tier": "5_pixels",
            "group_label_as_printed": "Control", "timepoint": tp,
            "ionizable_lipid": "",
            "readout_as_printed":
                "Fluorescene Intensity / area (arb.units)",
            "units": "arbitrary units, linear axis 0-5000",
            "value": b["value"] if not b["refused"] else "",
            "cargo": "none", "cargo_type": "", "cell":
                "primary rat cortical neuron", "organ": "", "route":
                "in vitro", "dose": "", "assay":
                "GFP fluorescence intensity per area, imaging",
            "pooled": "no", "pool_size": "",
            "corpus_admissible": "no",
            "admissibility_reason": (
                "untreated control, carries no lipid"
                + ("" if not b["refused"] else
                   "; and bars.bar_top REFUSED this bar (about 3 px tall, "
                   "indistinguishable from the axis rule). A refusal is not "
                   "a zero."))})

    with open(os.path.join(HERE, "05_values.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(values[0].keys()))
        w.writeheader()
        w.writerows(values)

    # --- the z, over ONE experiment --------------------------------------
    emitted = [v for v in values if v["corpus_admissible"] == "yes"]
    assert len(emitted) == 3, emitted
    logs = [math.log10(v["value"]) for v in emitted]
    mu = sum(logs) / len(logs)
    sd = math.sqrt(sum((x - mu) ** 2 for x in logs) / len(logs))  # ddof=0
    zs = [(x - mu) / sd for x in logs]

    # --- the quality columns, through emit.emit_rows ----------------------
    items = [v["ionizable_lipid"] for v in emitted]
    records = {v["ionizable_lipid"]: {**v, "il": ILS[v["ionizable_lipid"]],
                                      "z": f"{z:.10f}"}
               for v, z in zip(emitted, zs)}
    note_common = (
        "Fig 4G bar digitised from the publisher raster (no SI table, no "
        "Source Data, no repo, article paywalled); PEG_SMILES left blank "
        "because the paper DRAWS a phosphoethanolamine-carbamate-PEG45 while "
        "the corpus name C14-PEG2000 stores a glycerol ether; formulation is "
        "the one optimised for C12-200, so the two benchmark arms are not "
        "formulation-matched to their own optimum.")

    def conf(item, rec):
        return rec["il"]["confidence"]

    def note(item, rec):
        if rec["il"]["drawn_in_paper"]:
            return ("Structure drawn in full in Fig 1A and taken from the "
                    "corpus (AC_2025, 498 rows, one structure); an anchor "
                    "assembled independently from that drawing matches it at "
                    "an exact full InChIKey. " + note_common)
        return ("Structure NOT drawn in the paper -- the name alone is given "
                "in Fig 4G; pinned from provenance/reagent_registry.py. "
                + note_common)

    quality = emit.emit_rows(
        items, records, confidence=conf, note=note,
        positional_ambiguity=lambda i, r: "0",
        structure_assumed=lambda i, r: "0",
        value_is_bound=lambda i, r: "",
        value_normalization_unit=lambda i, r: "experiment",
        within_model=lambda i, r: r["z"],
        label="MR_2025 Fig 4G")

    # --- the 79-column rows ----------------------------------------------
    out = []
    for v, z, q in zip(emitted, zs, quality):
        il = ILS[v["ionizable_lipid"]]
        row = {c: "" for c in header}
        row.update({
            "Experiment_ID": EXPERIMENT_ID,
            "Formulation_ID": "LNP4_backbone_" + v["ionizable_lipid"],
            "IL_name": v["ionizable_lipid"], "IL_SMILES": il["smiles"],
            "IL_molratio": "35.0",
            "HL_name": "DOPE", "HL_SMILES": RR.REGISTRY["DOPE"][0],
            "HL_molratio": "30.5",
            "CHL_name": "Cholesterol",
            "CHL_SMILES": RR.REGISTRY["Cholesterol"][0],
            "CHL_molratio": "32.0",
            "PEG_name": "C14PEG2k", "PEG_SMILES": "", "PEG_molratio": "2.5",
            "Model": "in_vitro",
            "Model_type": "primary_rat_cortical_neuron",
            "Model_target": "in_vitro",
            "Route_of_administration": "in_vitro",
            "Cargo": "pDNA", "Cargo_type": "GFP",
            "Dose_ug_nucleicacid": "",
            "Experiment_method": "GFP_fluorescence_normalized",
            "Experiment_value": f"{z:.10f}",
            "Publication_DOI": DOI, "Publication_PMID": PMID,
            "Publication_link": "https://doi.org/" + DOI,
        })
        row.update(q)
        out.append(row)

    with open(os.path.join(HERE, "06_rows.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=header)
        w.writeheader()
        w.writerows(out)

    qc = {
        "experiment_id": EXPERIMENT_ID,
        "experiment_id_free_in_corpus": True,
        "pmid": PMID,
        "n_distinct_ionizable_lipids_in_paper": 3,
        "n_formulations_in_paper": 18,
        "n_rows_emitted": len(out),
        "timepoint_emitted": TIMEPOINT,
        "values_emitted": {v["ionizable_lipid"]: v["value"] for v in emitted},
        "z": {v["ionizable_lipid"]: round(z, 6)
              for v, z in zip(emitted, zs)},
        "z_mean": round(sum(zs) / len(zs), 12),
        "z_sd_ddof0": round(math.sqrt(sum(x * x for x in zs) / len(zs)), 12),
        "z_over": "one z per experiment; one cell type, one cargo, one "
                  "timepoint, one formulation -- error class 22 cannot arise",
        "rows_with_a_caveat": emit.caveat_count(quality),
        "benchmarks_top_the_screen": False,
        "benchmark_ratio_C12_200_over_SM_102": round(
            emitted[0]["value"] / emitted[1]["value"], 4),
        "benchmark_ratio_C12_200_over_ALC_0315": round(
            emitted[0]["value"] / emitted[2]["value"], 4),
    }
    with open(os.path.join(HERE, "03_enumerated.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=[
            "IL_name", "IL_SMILES", "InChIKey", "drawn_in_paper",
            "structure_source", "confidence"])
        w.writeheader()
        for n, d in ILS.items():
            w.writerow({"IL_name": n, "IL_SMILES": d["smiles"],
                        "InChIKey": d["inchikey"],
                        "drawn_in_paper": int(d["drawn_in_paper"]),
                        "structure_source": d["source"],
                        "confidence": d["confidence"]})

    # --- 01_blocks: the DESIGN, which is what the paper really varies -----
    s1 = list(csv.DictReader(open(os.path.join(HERE, "02_tableS1.csv"),
                                  newline="")))
    with open(os.path.join(HERE, "01_blocks.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=[
            "block", "n_formulations", "ionizable_lipid", "varying_factor",
            "levels", "readout", "corpus_admissible", "reason"])
        w.writeheader()
        w.writerow({
            "block": "b-DNA barcode screen (Fig 1, Table S1)",
            "n_formulations": len(s1), "ionizable_lipid": "C12-200",
            "varying_factor": "helper lipid identity x helper molar content",
            "levels": "DOPE/DSPC/DOTAP x 15/20/25/30.5/35/50 = 18",
            "readout": "normalized b-DNA delivery by NGS, per organ",
            "corpus_admissible": "no",
            "reason": "barcode recovered from TISSUE is biodistribution, not "
                      "delivery (brief section 3)"})
        w.writerow({
            "block": "pDNA-GFP formulation screen (Fig 3, Table S3, "
                     "Figs S5/S6)",
            "n_formulations": 11, "ionizable_lipid": "C12-200",
            "varying_factor": "helper lipid identity x helper molar content",
            "levels": "subset of the 18",
            "readout": "GFP fluorescence intensity / area",
            "corpus_admissible": "in kind yes, not emitted",
            "reason": "every bar is the same ionizable lipid; a corpus keyed "
                      "on ionizable-lipid chemistry gains nothing"})
        w.writerow({
            "block": "ionizable-lipid substitution (Fig 4G)",
            "n_formulations": 3, "ionizable_lipid":
                "C12-200 / SM-102 / ALC-0315",
            "varying_factor": "THE IONIZABLE LIPID",
            "levels": "3", "readout": "GFP fluorescence intensity / area",
            "corpus_admissible": "yes",
            "reason": "the only ionizable-lipid-varying delivery measurement "
                      "in the paper"})
        w.writerow({
            "block": "in vivo luciferase mRNA (Fig 6, Fig S9)",
            "n_formulations": 3, "ionizable_lipid": "C12-200",
            "varying_factor": "DOPE molar content (20/25/30.5)",
            "levels": "3", "readout": "brain Total Flux [p/s]",
            "corpus_admissible": "in kind yes, not emitted",
            "reason": "one ionizable lipid; LNP2/LNP3/LNP4 differ only in "
                      "helper content"})

    print(json.dumps(qc, indent=1))
    json.dump(qc, open(os.path.join(HERE, "06_emit_qc.json"), "w"), indent=1)
    print(f"wrote 01_blocks.csv, 03_enumerated.csv, 05_values.csv "
          f"({len(values)} rows), 06_rows.csv ({len(out)} rows)")


if __name__ == "__main__":
    main()
