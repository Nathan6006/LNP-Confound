#!/usr/bin/env python3
"""Stage 04 -- 05_values.csv (every value, admissible or not) and 06_rows.csv
(the admissible ones, in the 79-column corpus schema).

RUN: python3 -I scripts/stage04_values.py <1080csv> <validation.xlsx> <vitro.csv> <outdir>

ADMISSIBILITY. data/'s Experiment_method is delivery-only: luminescence,
protein abundance, editing efficiency, knockdown. Size, zeta potential and PDI
are physicochemical and are NOT delivery endpoints, so they are recorded in
05_values.csv with corpus_admissible = no and never reach 06_rows.csv. That is
three quarters of this paper's 'targets': an ML dataset has several target
columns and only some are delivery.

The second, larger exclusion is specific to this paper and is the finding that
matters: RLU_HepG2 is ALREADY IN THE CORPUS as YZ_2022. Not similar -- the same
numbers. 1,080 of 1,080 formulations join one-to-one and the z-scored values
agree at Spearman 1.000000 with a maximum residual of 3.2e-4, which is the
3-decimal rounding of the repository file. Emitting it would double-count a
published screen.

NORMALISATION. Experiment_value is a per-experiment z-score of the log signal
(CLAUDE.md section 3), population sd. The stored values are ln(mean raw RLU),
and a per-block z of ln IS the per-block z of log10 -- the same affine argument
schema_quality.json makes for Experiment_value_within_model -- so no rescaling
happens anywhere here. ONE z over the whole experiment, pooling all cell types,
so the between-cell-type offset survives: z-scoring per (experiment x cell
type) is error class 22 and it is what deleted WT_2026's published selectivity
(F-016). value_normalization_unit = 'experiment' accordingly.

FLAGS go on the row through curation/lib/emit.py with PER-ITEM CALLABLES. A
literal confidence raises there by construction (F-011/F-025/F-046), and so
does a caveat with an empty note (F-024).
"""
import csv
import json
import math
import os
import sys

from openpyxl import load_workbook

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "lib"))
import emit                                                 # noqa: E402

OHE = ["18PG", "14PA", "DOPE", "DSPC", "DOTAP", "DDAB"]

# column label -> (Model_type as the corpus writes it, how we know)
CELLS = {
    "RLU_HepG2": ("HepG2", "ATCC HB-8065, Methods"),
    "RLU_N2a": ("Neuro2a", "ATCC CCL-131, Methods; the column label is N2a"),
    "RLU_ARPE19": ("ARPE-19", "ATCC CRL-2302, Methods; column label ARPE19"),
    "RLU_HEK293": ("HEK293T", "ATCC CRL-3216 IS HEK293T, per Methods, although "
                              "the column is labelled RLU_HEK293; the corpus "
                              "already uses HEK293T"),
    "RLU_B16": ("B16F10", "ATCC CRL-6475, Methods; column label B16"),
    "RLU_PC3": ("PC3", "ATCC CRL-1435, Methods"),
}
PHYS = {"Size": "nm, z-average by DLS",
        "Zeta": "mV, zeta potential by DLS in PBS",
        "PDI": "dimensionless polydispersity index by DLS"}

DOSE = "0.1"          # ug pDNA per well: 4 uL of 25 ug/mL into 100 uL medium
TIMEPOINT = "48 h"

PMID = "39375194"
DOI = "10.1021/acsnano.4c07615"
LINK = "https://pubs.acs.org/doi/10.1021/acsnano.4c07615"

NOTE_COMMON = (
    "Values are the authors' own ln(mean RLU) from "
    "github.com/MaoResearchGroup/ml_lnp_design_rules@82cee58 "
    "Raw_Data/, not digitised. Dose 0.1 ug pDNA/well and 48 h come from this "
    "paper's Methods, which describe the same screen that Nat Commun 2022 "
    "(PMID 35879315) describes at 0.2 ug/well and 24 h -- see "
    "08_contradictions.csv C-1.")
NOTE_18PG = (
    "Helper lipid '18PG' is emitted as DIOLEOYL PG (18:1/18:1), which is what "
    "SI Table 2 draws and what SI Table 3's Double_bonds = 2 states; the "
    "corpus registry pins the name '18PG' to SOPG (18:0/18:1) on the "
    "authority of PMID 35879315, which describes THIS SAME library. "
    "Unresolved: DECISIONS.md D-2.")
NOTE_DOTAP = (
    "DOTAP is emitted without a stereocentre, matching the reagent registry "
    "and 2,564 existing corpus rows; SI Table 2 draws a defined wedge at the "
    "propane C2, so the emitted stereochemistry is a convention, not a "
    "reading. See 08_contradictions.csv C-3.")


def f(x):
    return None if x is None or x == "" else float(x)


def zscore(v):
    n = len(v)
    m = sum(v) / n
    sd = math.sqrt(sum((x - m) ** 2 for x in v) / n)
    return [(x - m) / sd for x in v], m, sd


def load(repo_csv, xlsx):
    with open(repo_csv, newline="", encoding="utf-8-sig") as fh:
        tr = list(csv.DictReader(fh))
    wb = load_workbook(xlsx, data_only=True)
    ws = wb.worksheets[0]
    it = ws.iter_rows(values_only=True)
    hdr = [str(h).strip() for h in next(it)]
    va = [dict(zip(hdr, r)) for r in it if not all(v is None for v in r)]
    return tr, va


def molpct(r):
    a, h, p = (float(r["(IL+HL)"]), float(r["HL_(IL+HL)"]),
               float(r["PEG_(Chol+PEG)"]))
    return (a * (100 - h) / 100.0, a * h / 100.0,
            (100 - a) * (100 - p) / 100.0, (100 - a) * p / 100.0)


def main(repo_csv, xlsx, vitro, outdir):
    tr, va = load(repo_csv, xlsx)
    blocks = {r["published_name"]: r
              for r in csv.DictReader(open(os.path.join(outdir,
                                                        "01_blocks.csv")))}
    enum = list(csv.DictReader(open(os.path.join(outdir,
                                                 "03_enumerated.csv"))))[0]
    with open(vitro, newline="") as fh:
        corpus = list(csv.DictReader(fh))
    yz = [r for r in corpus if r["Experiment_ID"] == "YZ_2022"]
    peg_smiles = sorted(set(r["PEG_SMILES"] for r in yz))
    assert len(peg_smiles) == 1, peg_smiles

    # ------------------------------------------------- 05_values.csv
    vals = []

    def add(src_file, locus, lib, fid, hl, readout, units, value,
            admissible, reason, model_type="", assay=""):
        vals.append({
            "source_file": src_file,
            "source_locus": locus,
            "library": lib,
            "formula_label": fid,
            "helper_lipid": hl,
            "readout_as_printed": readout,
            "units_as_printed": units,
            "value": "" if value is None else value,
            "cargo": "pDNA (firefly luciferase pDNA, Aldevron)",
            "cell_line": model_type,
            "organ": "",
            "route": "in_vitro",
            "dose": DOSE + " ug pDNA per well (1 ug pDNA/mL final)"
                    if admissible == "yes" or readout.startswith("RLU") else "",
            "timepoint": TIMEPOINT if readout.startswith("RLU") else "at formulation",
            "assay": assay,
            "corpus_admissible": admissible,
            "admissibility_reason": reason,
        })

    REASON_DUP = ("NO -- delivery endpoint, but it is ALREADY IN THE CORPUS. "
                  "All 1,080 formulations join YZ_2022 one-to-one and the "
                  "z-scored values agree at Spearman 1.000000, max residual "
                  "3.2e-4 (the file's 3-decimal rounding). Emitting would "
                  "double-count PMID 35879315.")
    REASON_OK = ("YES -- luminescence from a reporter transgene, the corpus's "
                 "luminescence_normalized endpoint; this cell line is not in "
                 "the corpus for this library.")
    REASON_BLANK = ("NO -- never measured. 14PA and DDAB were not run in "
                    "Neuro2a or ARPE-19; the authors' own "
                    "utilities.py::select_input_params drops those two "
                    "one-hot columns for exactly these two cell types.")
    REASON_PHYS = ("NO -- physicochemical nanoparticle property, not a "
                   "delivery endpoint. data/'s Experiment_method vocabulary "
                   "is luminescence / protein abundance / editing efficiency "
                   "/ knockdown only.")
    REASON_VAL = ("YES -- luminescence from a reporter transgene, measured on "
                  "a 72-formulation de-novo library whose composition levels "
                  "are disjoint from the training library and absent from the "
                  "corpus.")

    for r in tr:
        hl = [c for c in OHE if r[c].strip() == "1"][0]
        for col, (mt, _) in CELLS.items():
            v = r[col].strip()
            if v == "":
                add("1080_LNP_OHE.csv", col, "training_1080",
                    r["Formula label"], hl, col, "ln(mean RLU)", None,
                    "no", REASON_BLANK, mt, "firefly luciferase, luminometer")
            elif col == "RLU_HepG2":
                add("1080_LNP_OHE.csv", col, "training_1080",
                    r["Formula label"], hl, col, "ln(mean RLU)", float(v),
                    "no", REASON_DUP, mt, "firefly luciferase, luminometer")
            else:
                add("1080_LNP_OHE.csv", col, "training_1080",
                    r["Formula label"], hl, col, "ln(mean RLU)", float(v),
                    "yes", REASON_OK, mt, "firefly luciferase, luminometer")
        for col, units in PHYS.items():
            add("1080_LNP_OHE.csv", col, "training_1080",
                r["Formula label"], hl, col, units, float(r[col]),
                "no", REASON_PHYS, "", "dynamic light scattering, Nano ZS90")

    for r in va:
        hl = str(r["Helper_lipid"]).strip()
        for col, (mt, _) in CELLS.items():
            v = f(r.get(col))
            add("Compiled_validation.xlsx", "Sheet1!" + col,
                "validation_72", str(r["Formula label"]), hl, col,
                "ln(mean RLU)", v,
                "yes" if v is not None else "no",
                REASON_VAL if v is not None else REASON_BLANK,
                mt, "firefly luciferase, luminometer")
        for col, units in PHYS.items():
            add("Compiled_validation.xlsx", "Sheet1!" + col, "validation_72",
                str(r["Formula label"]), hl, col, units, f(r.get(col)),
                "no", REASON_PHYS, "", "dynamic light scattering, Nano ZS90")

    vpath = os.path.join(outdir, "05_values.csv")
    with open(vpath, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(vals[0].keys()))
        w.writeheader()
        w.writerows(vals)

    # ------------------------------------------------- 06_rows.csv
    SCHEMA = [c.strip() for c in
              open(os.path.join(HERE, "..", "..", "..", "data",
                                "data_vitro.csv")).readline().strip().split(",")]
    DESC = [c for c in SCHEMA if c in enum]

    def build(lib_rows, lib, exp_id, source_file):
        items, records = [], {}
        for r in lib_rows:
            if lib == "training_1080":
                hl = [c for c in OHE if r[c].strip() == "1"][0]
                fid = r["Formula label"]
                np_ratio = float(r["NP_ratio"])
                cols = CELLS
            else:
                hl = str(r["Helper_lipid"]).strip()
                fid = str(r["Formula label"])
                np_ratio = float(r["NP_ratio"])
                cols = CELLS
            il, hlp, chl, peg = molpct(r)
            for col, (mt, _) in cols.items():
                v = (f(r[col]) if lib == "training_1080"
                     and r[col].strip() != "" else
                     (f(r.get(col)) if lib != "training_1080" else None))
                if v is None:
                    continue
                if lib == "training_1080" and col == "RLU_HepG2":
                    continue            # already in the corpus as YZ_2022
                iid = "%s|%s|%s" % (exp_id, fid, col)
                items.append(iid)
                records[iid] = {"hl": hl, "fid": fid, "cell": mt, "ln": v,
                                "np": np_ratio, "il": il, "hlp": hlp,
                                "chl": chl, "peg": peg, "col": col,
                                "lib": lib, "source_file": source_file}
        return items, records

    out_rows, summary = [], {}
    for lib, src, exp_id, data in (
            ("training_1080", "1080_LNP_OHE.csv", "LC_2024", tr),
            ("validation_72", "Compiled_validation.xlsx", "LC_2024_2", va)):
        items, records = build(data, lib, exp_id, src)
        z, mean, sd = zscore([records[i]["ln"] for i in items])
        for i, zz in zip(items, z):
            records[i]["z"] = zz

        def conf(item_id, rec):
            # PER-ITEM, from the row's own helper lipid. A literal raises.
            return "medium" if rec["hl"] == "18PG" else "high"

        def note(item_id, rec):
            parts = [NOTE_COMMON]
            if rec["hl"] == "18PG":
                parts.insert(0, NOTE_18PG)
            if rec["hl"] == "DOTAP":
                parts.insert(0, NOTE_DOTAP)
            return " ".join(parts)

        def assumed(item_id, rec):
            return "1" if rec["hl"] == "DOTAP" else "0"

        flags = emit.emit_rows(
            items, records,
            confidence=conf, note=note,
            positional_ambiguity=lambda i, r: "0",
            structure_assumed=assumed,
            value_is_bound=lambda i, r: "",
            value_normalization_unit=lambda i, r: "experiment",
            within_model=None, label=exp_id)

        for iid, fl in zip(items, flags):
            rec = records[iid]
            row = {c: "" for c in SCHEMA}
            row.update({
                "Experiment_ID": exp_id,
                "Formulation_ID": "%s_F%s" % (exp_id, rec["fid"]),
                "IL_name": enum["IL_name"],
                "IL_SMILES": enum["IL_SMILES"],
                "IL_molratio": round(rec["il"], 4),
                "IL_to_nucleicacid_chargeratio": rec["np"],
                "HL_name": rec["hl"],
                "HL_SMILES": blocks[rec["hl"]]["SMILES"],
                "HL_molratio": round(rec["hlp"], 4),
                "CHL_name": "Cholesterol",
                "CHL_SMILES": blocks["Cholesterol"]["SMILES"],
                "CHL_molratio": round(rec["chl"], 4),
                "PEG_name": "DMG-PEG2000",
                "PEG_SMILES": peg_smiles[0],
                "PEG_molratio": round(rec["peg"], 4),
                "Aqueous_buffer": "acetate",
                "Dialysis_buffer": "",
                "Mixing_method": "handmixed",
                "Model": "in_vitro",
                "Model_type": rec["cell"],
                "Model_target": "in_vitro",
                "Route_of_administration": "in_vitro",
                "Cargo": "pDNA",
                "Cargo_type": "FLuc",
                "Dose_ug_nucleicacid": DOSE,
                "Experiment_method": "luminescence_normalized",
                "Experiment_batching": "individual",
                "Experiment_value": rec["z"],
                "Publication_link": LINK,
                "Publication_DOI": DOI,
                "Publication_PMID": PMID,
            })
            for c in DESC:
                row[c] = enum[c]
            row.update(fl)
            out_rows.append(row)
        summary[exp_id] = {
            "library": lib, "source_file": src, "rows": len(items),
            "cells": sorted(set(records[i]["cell"] for i in items)),
            "ln_mean": round(mean, 6), "ln_pop_sd": round(sd, 6),
            "z_mean": round(sum(z) / len(z), 12),
            "z_pop_sd": round(math.sqrt(sum(x * x for x in z) / len(z)), 12),
            "caveat_rows": emit.caveat_count(
                [r for r in out_rows if r["Experiment_ID"] == exp_id]),
        }

    rpath = os.path.join(outdir, "06_rows.csv")
    with open(rpath, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=SCHEMA)
        w.writeheader()
        w.writerows(out_rows)

    # ------------------------------------------------- counts must sum
    n_yes = sum(1 for v in vals if v["corpus_admissible"] == "yes")
    n_no = sum(1 for v in vals if v["corpus_admissible"] == "no")
    counts = {
        "values_rows_total": len(vals),
        "admissible_yes": n_yes,
        "admissible_no": n_no,
        "yes_plus_no_equals_total": n_yes + n_no == len(vals),
        "emitted_rows_total": len(out_rows),
        "emitted_equals_admissible": len(out_rows) == n_yes,
        # NOTE: this dict previously used the key "  of which emitted" TWICE,
        # under training and again under validation. Python keeps the last, so
        # the training figure was silently overwritten by the validation one
        # and the breakdown reported 432 where 4,680 was true -- while every
        # total still summed correctly, because the totals are computed
        # separately. Keys are unique now and the sub-counts are asserted
        # against the totals below.
        "breakdown": {
            "training_RLU_cells": 1080 * 6,
            "training_RLU_blank_never_run": sum(
                1 for v in vals if v["library"] == "training_1080"
                and v["readout_as_printed"].startswith("RLU")
                and v["value"] == ""),
            "training_RLU_HepG2_duplicate_of_YZ_2022": sum(
                1 for v in vals if v["library"] == "training_1080"
                and v["readout_as_printed"] == "RLU_HepG2"),
            "training_RLU_emitted": sum(
                1 for v in vals if v["library"] == "training_1080"
                and v["corpus_admissible"] == "yes"),
            "training_physicochemical": 1080 * 3,
            "validation_RLU_cells": 72 * 6,
            "validation_RLU_blank_never_run": sum(
                1 for v in vals if v["library"] == "validation_72"
                and v["readout_as_printed"].startswith("RLU")
                and v["value"] == ""),
            "validation_RLU_emitted": sum(
                1 for v in vals if v["library"] == "validation_72"
                and v["corpus_admissible"] == "yes"),
            "validation_physicochemical": 72 * 3,
        },
        "per_experiment": summary,
    }
    json.dump(counts, open(os.path.join(outdir, "05_value_counts.json"), "w"),
              indent=1)
    print(json.dumps(counts, indent=1))
    b = counts["breakdown"]
    counts["breakdown_sums_to_total"] = (
        b["training_RLU_cells"] + b["training_physicochemical"]
        + b["validation_RLU_cells"] + b["validation_physicochemical"]
        == counts["values_rows_total"])
    counts["training_RLU_parts_sum"] = (
        b["training_RLU_blank_never_run"]
        + b["training_RLU_HepG2_duplicate_of_YZ_2022"]
        + b["training_RLU_emitted"] == b["training_RLU_cells"])
    counts["emitted_parts_sum"] = (
        b["training_RLU_emitted"] + b["validation_RLU_emitted"]
        == counts["emitted_rows_total"])
    ok = (counts["yes_plus_no_equals_total"]
          and counts["emitted_equals_admissible"]
          and counts["breakdown_sums_to_total"]
          and counts["training_RLU_parts_sum"]
          and counts["emitted_parts_sum"])
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:5]))
