#!/usr/bin/env python3
"""Stage 05 -- validate 06_rows.csv, sweep the corpus for near misses, and
assemble 07_qc.json from every stage's measurements.

RUN: python3 -I scripts/stage05_qc.py <vitro.csv> <vivo.csv> <validation.xlsx> <outdir>
"""
import csv
import json
import math
import os
import sys

from openpyxl import load_workbook
from rdkit import Chem
from rdkit import DataStructs
from rdkit import RDLogger
from rdkit.Chem import rdFingerprintGenerator

RDLogger.DisableLog("rdApp.*")

OHE = ["18PG", "14PA", "DOPE", "DSPC", "DOTAP", "DDAB"]
GEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


def fp(smi, cache={}):
    if smi not in cache:
        m = Chem.MolFromSmiles(smi)
        cache[smi] = None if m is None else GEN.GetFingerprint(m)
    return cache[smi]


def main(vitro, vivo, xlsx, outdir):
    Q = {}
    rows = list(csv.DictReader(open(os.path.join(outdir, "06_rows.csv"))))
    schema = open(vitro).readline().strip().split(",")
    got = list(rows[0].keys())

    # ------------------------------------------------ 1. schema conformance
    Q["schema"] = {
        "corpus_columns": len(schema),
        "emitted_columns": len(got),
        "identical_header_in_order": schema == got,
        "missing": [c for c in schema if c not in got],
        "extra": [c for c in got if c not in schema],
    }

    # ------------------------------------------------ 2. row-level integrity
    REQUIRED = ["Experiment_ID", "Formulation_ID", "IL_name", "IL_SMILES",
                "HL_name", "HL_SMILES", "CHL_name", "CHL_SMILES", "PEG_name",
                "PEG_SMILES", "Model", "Model_type", "Cargo", "Cargo_type",
                "Dose_ug_nucleicacid", "Experiment_method", "Experiment_value",
                "Publication_PMID", "value_normalization_unit"]
    blank = {c: sum(1 for r in rows if not r[c].strip()) for c in REQUIRED}
    dup = {}
    for r in rows:
        k = (r["Experiment_ID"], r["Formulation_ID"], r["Model_type"])
        dup[k] = dup.get(k, 0) + 1
    unparse = sum(1 for s in set(r["IL_SMILES"] for r in rows)
                  if Chem.MolFromSmiles(s) is None)
    noncanon = sum(1 for s in set(r["IL_SMILES"] for r in rows)
                   if Chem.MolToSmiles(Chem.MolFromSmiles(s)) != s)
    Q["row_integrity"] = {
        "rows": len(rows),
        "required_columns_blank": {k: v for k, v in blank.items() if v},
        "all_required_populated": all(v == 0 for v in blank.values()),
        "duplicate_(exp,formulation,cell)_keys":
            sum(1 for v in dup.values() if v > 1),
        "unparseable_IL_SMILES": unparse,
        "non_canonical_IL_SMILES": noncanon,
        "padded_cells": sum(1 for r in rows for c in r
                            if r[c] != r[c].strip()),
    }

    # ------------------------------------------------ 3. molar percents
    worst, bad = 0.0, 0
    for r in rows:
        s = sum(float(r[c]) for c in ("IL_molratio", "HL_molratio",
                                      "CHL_molratio", "PEG_molratio"))
        worst = max(worst, abs(s - 100.0))
        if abs(s - 100.0) > 1e-6:
            bad += 1
    Q["molar_percents"] = {
        "rows_checked": len(rows),
        "rows_off_100_by_more_than_1e-6": bad,
        "max_abs_deviation_from_100": worst,
        "CREDITED_AS_A_CHECK": False,
        "why_not": ("IL+HL = a and CHL+PEG = 100-a by construction, so the "
                    "sum is 100 identically for any input, including garbage. "
                    "It is reported because its ABSENCE would look like an "
                    "oversight, and because a non-zero value here would mean "
                    "a rounding bug in the emitter -- but it cannot detect a "
                    "wrong composition. The checks that CAN are in "
                    "01_checks.json: the row-wise pct == 100/(1+ratio) "
                    "identity against the file's second parameterisation, and "
                    "the levels against SI Table 1."),
    }

    # ------------------------------------------------ 4. z-score convention
    zq = {}
    for exp in sorted(set(r["Experiment_ID"] for r in rows)):
        v = [float(r["Experiment_value"]) for r in rows
             if r["Experiment_ID"] == exp]
        m = sum(v) / len(v)
        sd = math.sqrt(sum((x - m) ** 2 for x in v) / len(v))
        zq[exp] = {"n": len(v), "mean": m, "pop_sd": sd,
                   "mean_within_1e-9": abs(m) < 1e-9,
                   "sd_within_1e-9": abs(sd - 1) < 1e-9,
                   "cells_pooled": sorted(set(
                       r["Model_type"] for r in rows
                       if r["Experiment_ID"] == exp))}
    Q["zscore"] = {
        "per_experiment": zq,
        "unit": "experiment",
        "why_not_per_cell_type": (
            "Pooling all cell types into one z keeps the between-cell-type "
            "offset. Standardising per (experiment x Model_type) is error "
            "class 22 and is what deleted WT_2026's own published selectivity "
            "(F-016). Experiment_value_within_model is left blank for "
            "build_dataset.py to compute."),
        "ln_vs_log10": (
            "The stored source values are ln(mean RLU) and the corpus "
            "convention is a z of log10. A per-block z of an affine function "
            "of a log signal IS the per-block z of the signal, so no "
            "conversion is applied or needed; this is the same argument "
            "schema_quality.json makes for Experiment_value_within_model."),
    }

    # ------------------------------------------------ 5. validation vs corpus
    corpus = list(csv.DictReader(open(vitro)))
    corpus += list(csv.DictReader(open(vivo)))
    wb = load_workbook(xlsx, data_only=True)
    ws = wb.worksheets[0]
    it = ws.iter_rows(values_only=True)
    hdr = [str(h).strip() for h in next(it)]
    va = [dict(zip(hdr, r)) for r in it if not all(v is None for v in r)]

    def key(hl, il, h, c, p, npr):
        return (hl, round(il, 1), round(h, 1), round(c, 1), round(p, 1),
                round(npr, 2))
    ck = set()
    for r in corpus:
        try:
            ck.add(key(r["HL_name"].strip(), float(r["IL_molratio"]),
                       float(r["HL_molratio"]), float(r["CHL_molratio"]),
                       float(r["PEG_molratio"]),
                       float(r["IL_to_nucleicacid_chargeratio"] or -1)))
        except ValueError:
            continue
    hits = 0
    for r in va:
        a, h, p = (float(r["(IL+HL)"]), float(r["HL_(IL+HL)"]),
                   float(r["PEG_(Chol+PEG)"]))
        k = key(str(r["Helper_lipid"]).strip(), a * (100 - h) / 100.0,
                a * h / 100.0, (100 - a) * (100 - p) / 100.0,
                (100 - a) * p / 100.0, float(r["NP_ratio"]))
        if k in ck:
            hits += 1
    Q["validation_library_vs_corpus"] = {
        "validation_formulations": len(va),
        "matching_a_corpus_formulation": hits,
        "note": ("The 72-formulation de-novo library uses composition levels "
                 "the training library never used, so it cannot collide with "
                 "YZ_2022/YZ_2024, and no other corpus experiment carries "
                 "this formulation grid. Checked against ALL 30,998 corpus "
                 "rows, not just the two obvious experiments."),
    }

    # ------------------------------------------------ 6. near-miss sweep
    # Honestly labelled a SUPERSET: Tanimoto neighbours are candidates for
    # 'the same compound under another name', never a determination.
    target = rows[0]["IL_SMILES"]
    tfp = fp(target)
    seen, near = {}, []
    for r in corpus:
        s = (r.get("IL_SMILES") or "").strip()
        if not s or s in seen:
            continue
        f = fp(s)
        seen[s] = 0.0 if f is None else DataStructs.TanimotoSimilarity(tfp, f)
    for s, t in seen.items():
        if t >= 0.70:
            ex = sorted(set(r["Experiment_ID"] for r in corpus
                            if (r.get("IL_SMILES") or "").strip() == s))
            n = sum(1 for r in corpus
                    if (r.get("IL_SMILES") or "").strip() == s)
            m = Chem.MolFromSmiles(s)
            near.append({"tanimoto": round(t, 4), "rows": n,
                         "experiments": ex,
                         "exact_same_molecule": t == 1.0,
                         "inchikey": Chem.MolToInchiKey(m),
                         "names": sorted(set(
                             r["IL_name"] for r in corpus
                             if (r.get("IL_SMILES") or "").strip() == s))[:4]})
    near.sort(key=lambda x: -x["tanimoto"])
    Q["near_miss_sweep"] = {
        "query": "DLin-MC3-DMA, the paper's only ionizable lipid",
        "fingerprint": "Morgan radius 2, 2048 bits (CLAUDE.md section 7)",
        "threshold_tanimoto": 0.70,
        "distinct_corpus_IL_structures_scanned": len(seen),
        "n_neighbours": len(near),
        "neighbours": near[:12],
        "THIS_IS_A_SUPERSET": True,
        "how_to_read_it": (
            "A Tanimoto neighbour is a LEAD, not a match. A sibling's first "
            "SMARTS returned 332 under a label claiming exactness when the "
            "true figure was 128. Only the rows with exact_same_molecule = "
            "true are the same compound; everything else is reported for a "
            "human to read and is acted on nowhere."),
    }

    # ------------------------------------------------ 7. fold in the stages
    for name, path in (("stage01_file_checks", "01_checks.json"),
                       ("stage01_seeded_mutations", "01_seed_results.json"),
                       ("stage02_corpus_overlap", "02_corpus_overlap.json"),
                       ("stage03_structures", "03_structures_qc.json"),
                       ("stage03b_c10_anchors", "04_c10_anchor_check.json"),
                       ("stage04_emit_seeds", "04_emit_seed_results.json"),
                       ("stage04_value_counts", "05_value_counts.json")):
        p = os.path.join(outdir, path)
        if os.path.exists(p):
            Q[name] = json.load(open(p))

    # drop the bulky per-item bodies so 07_qc.json stays readable
    if "stage02_corpus_overlap" in Q:
        so = Q["stage02_corpus_overlap"].get("structure_overlap", {})
        for v in so.values():
            v.pop("SKELETON_only_hits_stereo_differs", None) if not v.get(
                "n_rows_skeleton_only") else None

    fails = []
    if not Q["schema"]["identical_header_in_order"]:
        fails.append("emitted header is not data_vitro.csv's header")
    if not Q["row_integrity"]["all_required_populated"]:
        fails.append("a required column is blank")
    if Q["row_integrity"]["duplicate_(exp,formulation,cell)_keys"]:
        fails.append("duplicate (experiment, formulation, cell) keys")
    if Q["row_integrity"]["unparseable_IL_SMILES"]:
        fails.append("unparseable IL_SMILES")
    if Q["row_integrity"]["non_canonical_IL_SMILES"]:
        fails.append("non-canonical IL_SMILES")
    if Q["row_integrity"]["padded_cells"]:
        fails.append("padded cells (F-022 / error class 37)")
    for e, v in zq.items():
        if not (v["mean_within_1e-9"] and v["sd_within_1e-9"]):
            fails.append("%s is not at mean 0 / sd 1" % e)
    if Q["validation_library_vs_corpus"]["matching_a_corpus_formulation"]:
        fails.append("a validation formulation already exists in the corpus")
    Q["FAILURES"] = fails
    json.dump(Q, open(os.path.join(outdir, "07_qc.json"), "w"), indent=1)
    print(json.dumps({"schema": Q["schema"],
                      "row_integrity": Q["row_integrity"],
                      "molar_percents": {k: v for k, v in
                                         Q["molar_percents"].items()
                                         if k != "why_not"},
                      "zscore": {k: {kk: vv for kk, vv in v.items()}
                                 for k, v in zq.items()},
                      "validation_vs_corpus":
                          Q["validation_library_vs_corpus"]
                          ["matching_a_corpus_formulation"],
                      "near_miss": {"n": len(near),
                                    "top": near[:5]},
                      "FAILURES": fails}, indent=1))
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:5]))
