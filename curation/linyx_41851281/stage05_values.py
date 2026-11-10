#!/usr/bin/env python3
"""Measured values -> 05_values.csv, with the quality columns via lib/emit.

EVERY VALUE HERE IS MEASURED. The paper's model outputs live in
MOESM10 'Fig. 4' column C ('Predicted') and are not read by this file at all;
00_value_sources.json entry V6 says why. No predicted number is emitted.

VALUE SOURCES, one per row, named to the sheet and cell range:
  A  MOESM4 'Supplementary Fig. 10'  A1:D113   lipids 1-110 + 3 benchmarks
  B  MOESM10 'Fig. 4'                E2:F23    P1-P6 + ALC-0315, in vitro
  C  MOESM10 'Fig. 4'                E24:F45   P1-P6 + ALC-0315, in vivo
  D  MOESM16 'Extended Data Fig. 2'  A2:B46    Lipid T1, Lipid T2 (2c)
  E  MOESM9  'Fig. 3'                 W2:X43    Fig. 3i, in vivo: PBS,
                                                ALC-0315 and 12 library lipids

Source A and source B give ALC-0315 the SAME three numbers
(957212.2 / 948389.9 / 816056.5), so Fig. 4c's in-vitro panel is the same
measurement as the screen, not a repeat of it -- which is what licenses
carrying the screen's assay annotation onto the P rows (DECISIONS.md D-2).
"""
import csv
import json
import math
import os
import re
import sys

import openpyxl

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import emit                                                # noqa: E402
import reagents                                            # noqa: E402
from rdkit import Chem, RDLogger                           # noqa: E402

RDLogger.DisableLog("rdApp.*")
SRC = os.path.join(HERE, "source")

# Supplementary Table 3, MOESM1 p.39. Every one of the 8 formulations in that
# table sums to exactly 100 mol% -- checked in 07_qc.json, and it is an
# invariant the authors never printed.
FORM_LIBRARY = "ionisable lipid:DOPE:Chol:DMG-PEG2000 = 45:15:38.5:1.5"
FORM_BENCH = {
    "DLin-MC3-DMA": "MC3:DSPC:Chol:DMG-PEG2000 = 50:10:38.5:1.5",
    "ALC-0315": "ALC-0315:DSPC:Chol:ALC-0159 = 46.3:9.4:42.7:1.6",
    "SM-102": "SM-102:DSPC:Chol:DMG-PEG2000 = 50:10:38.5:1.5",
}
BENCH_NAME = {"MC3": "DLin-MC3-DMA", "ALC-0315": "ALC-0315",
              "SM-102": "SM-102"}
# Structures for the three commercial benchmarks come from the corpus reagent
# registry (C11), never typed into this file -- `18PG` became a symmetric
# 18:0/18:0 that way when it is SOPG. reagents.smi refuses an ambiguous name.
BENCH_SMILES = {n: reagents.smi(n) for n in BENCH_NAME.values()}
BENCH_KEY = {n: Chem.MolToInchiKey(Chem.MolFromSmiles(s))
             for n, s in BENCH_SMILES.items()}

VITRO = dict(
    assay="luciferase luminescence (RLU)", cargo="mFluc mRNA",
    system="HEK293T", dose="1 ug/mL mRNA", timepoint="24 h",
    route="", n=3,
    evidence="Supplementary Fig. 10 legend, MOESM1 p.14")
VIVO = dict(
    assay="whole-animal luciferase total flux (p/s)", cargo="Fluc mRNA",
    system="mouse", dose="", timepoint="", route="", n=3,
    evidence="main-text Fig. 4c legend NOT IN HAND; strain, route, dose and "
             "timepoint are not stated anywhere in the SI or the peer-review "
             "file. Left blank. DECISIONS.md D-3.")


def sheet(book, name):
    wb = openpyxl.load_workbook(os.path.join(SRC, book), data_only=True)
    return list(wb[name].iter_rows(values_only=True))


def main():
    enum = {r["ml_id"]: r for r in
            csv.DictReader(open(os.path.join(HERE, "03_enumerated.csv")))}
    by_name = {r["lipid_name"]: r for r in enum.values()}
    names = json.load(open(os.path.join(HERE, "00_names.json")))

    items, recs, vals = [], {}, {}

    # ---- A: the screen -------------------------------------------------
    for i, r in enumerate(sheet("41551_2026_1640_MOESM4_ESM.xlsx",
                                "Supplementary Fig. 10"), 1):
        if r[0] is None:
            continue
        reps = [r[1], r[2], r[3]]
        key = str(r[0])
        if key in BENCH_NAME:
            iid = "A:%s" % BENCH_NAME[key]
            b = BENCH_NAME[key]
            rec = dict(kind="benchmark", name=b, ml_id="", lipid_name=b,
                       smiles=BENCH_SMILES[b], inchikey=BENCH_KEY[b],
                       formulation=FORM_BENCH[b])
        else:
            e = enum[str(int(float(key)))]
            iid = "A:%s" % e["ml_id"]
            rec = dict(kind="library", name="Lipid %s" % e["ml_id"],
                       ml_id=e["ml_id"], lipid_name=e["lipid_name"],
                       smiles=e["smiles"], inchikey=e["inchikey"],
                       formulation=FORM_LIBRARY, enum=e)
        rec.update(pool="in_vitro", source="A", book="MOESM4",
                   sheet_name="Supplementary Fig. 10",
                   cells="A%d:D%d" % (i, i), **VITRO)
        items.append(iid)
        recs[iid] = rec
        vals[iid] = reps

    # ---- B/C: Fig. 4c --------------------------------------------------
    rows4 = sheet("41551_2026_1640_MOESM10_ESM.xlsx", "Fig. 4")
    block = None
    cur, start = None, None
    for i, r in enumerate(rows4, 1):
        lab, v = r[4], r[5]
        if lab == "In virto":
            block = "in_vitro"
            continue
        if lab == "In vivo":
            block = "in_vivo"
            continue
        if lab:
            cur, start = str(lab).strip(), i
        if v is None or cur is None or block is None:
            continue
        key = (block, cur)
        vals.setdefault(key, []).append(v)
        if key not in recs:
            nm = cur.replace("Lipid ", "")
            if nm in ("ALC-0315",):
                rec = dict(kind="benchmark", name="ALC-0315", ml_id="",
                           lipid_name="ALC-0315",
                           smiles=BENCH_SMILES["ALC-0315"],
                           inchikey=BENCH_KEY["ALC-0315"],
                           formulation=FORM_BENCH["ALC-0315"])
            else:
                hit = [e for e in enum.values() if e["p_label"] == nm][0]
                rec = dict(kind="library_AI_designed", name="Lipid %s" % nm,
                           ml_id=hit["ml_id"], lipid_name=hit["lipid_name"],
                           smiles=hit["smiles"], inchikey=hit["inchikey"],
                           formulation=FORM_LIBRARY, enum=hit)
            rec.update(pool=block, source="B" if block == "in_vitro" else "C",
                       book="MOESM10", sheet_name="Fig. 4",
                       cells="E%d:F%d" % (start, start + 2),
                       **(VITRO if block == "in_vitro" else VIVO))
            recs[key] = rec
            items.append(key)

    # ALC-0315's in-vitro triplicate appears in BOTH source A and source B.
    # They are not two measurements: assert they are the same three numbers
    # and keep only source A, so no lipid is emitted twice.
    a = vals["A:ALC-0315"]
    b = vals[("in_vitro", "ALC-0315")]
    if [float(x) for x in a] != [float(x) for x in b]:
        raise SystemExit("ALC-0315 differs between Supp Fig 10 %s and "
                         "Fig 4c %s -- they are then two measurements and "
                         "both must be emitted." % (a, b))
    items.remove(("in_vitro", "ALC-0315"))
    del recs[("in_vitro", "ALC-0315")], vals[("in_vitro", "ALC-0315")]

    # ---- D: Lipid T1 / T2 ----------------------------------------------
    rowsE = sheet("41551_2026_1640_MOESM16_ESM.xlsx", "Extended Data Fig. 2")
    cur, start = None, None
    seen = {}
    for i, r in enumerate(rowsE, 1):
        if r[0]:
            cur, start = str(r[0]).strip(), i
        if r[1] is None or cur is None:
            continue
        seen.setdefault(cur, (start, []))[1].append(r[1])
    for nm in ("Lipid T1", "Lipid T2"):
        start, reps = seen[nm]
        iid = "D:%s" % nm
        recs[iid] = dict(kind="external_ML_test", name=nm, ml_id="",
                         lipid_name=nm, smiles="", inchikey="",
                         formulation=FORM_LIBRARY, pool="in_vitro", source="D",
                         book="MOESM16", sheet_name="Extended Data Fig. 2",
                         cells="A%d:B%d" % (start, start + 2), **VITRO)
        items.append(iid)
        vals[iid] = reps

    # ---- E: Fig. 3i, in vivo ------------------------------------------
    rows3 = sheet("41551_2026_1640_MOESM9_ESM.xlsx", "Fig. 3")
    cur, start = None, None
    seen3 = {}
    for i, r in enumerate(rows3, 1):
        if r[21]:
            cur, start = str(r[21]).strip(), i
        if r[22] is None or cur is None:
            continue
        seen3.setdefault(cur, [start, []])[1].append(r[22])
    # Fig. 3i and Fig. 4c share ALC-0315's three in-vivo numbers, so they are
    # one experiment, not two. Assert it, then emit ALC-0315 once.
    a3 = [float(x) for x in seen3["ALC-0315"][1]]
    a4 = [float(x) for x in vals[("in_vivo", "ALC-0315")]]
    if a3 != a4:
        raise SystemExit("ALC-0315 in vivo differs: Fig 3i %s vs Fig 4c %s"
                         % (a3, a4))
    for nm, (start, reps) in seen3.items():
        if nm in ("PBS", "ALC-0315"):
            continue                      # PBS is a vehicle, not a lipid
        mlid = nm.replace("Lipid ", "")
        e = enum[mlid]
        iid = "E:%s" % mlid
        recs[iid] = dict(kind="library", name="Lipid %s" % mlid, ml_id=mlid,
                         lipid_name=e["lipid_name"], smiles=e["smiles"],
                         inchikey=e["inchikey"], formulation=FORM_LIBRARY,
                         enum=e, pool="in_vivo", source="E", book="MOESM9",
                         sheet_name="Fig. 3",
                         cells="W%d:X%d" % (start, start + 2), **VIVO)
        items.append(iid)
        vals[iid] = reps

    # ---- quality columns, through the one audited writer ---------------
    def confidence(iid, rec):
        if rec["kind"] == "external_ML_test":
            return "low"
        if rec["kind"] == "benchmark":
            return "high"
        if rec["enum"]["positional_ambiguity"] == "1":
            return "medium"
        if rec["enum"]["hrms_dev_Da"] in ("", None):
            return "medium"
        return "high"

    def note(iid, rec):
        if rec["kind"] == "external_ML_test":
            return ("Structure NOT derived: this lipid is outside the H-T "
                    "combinatorial library and is drawn only in Supplementary "
                    "Fig. 12. Its arm formula is pinned by HRMS (see "
                    "07_qc.json) but the connectivity was not re-derived here. "
                    "DECISIONS.md D-5.")
        if rec["kind"] == "benchmark":
            return ("Commercial benchmark; structure taken from the corpus "
                    "reagent registry, not from this paper. Formulated "
                    "DIFFERENTLY from the library lipids: %s, against the "
                    "library's %s." % (rec["formulation"], FORM_LIBRARY))
        e = rec["enum"]
        if e["positional_ambiguity"] == "1":
            return ("%d arms on a head with %s N-H sites: %s distinct "
                    "regiochemistries fit the formula and the paper neither "
                    "draws nor states which. The maximal-dispersion "
                    "convention was used. DECISIONS.md D-1."
                    % (int(e["n_tails"]), e["n_NH_sites"],
                       e["n_distinct_isomers"]))
        if e["hrms_dev_Da"] in ("", None):
            return ("No printed exact mass for this lipid, so the enumerated "
                    "structure is unconfirmed by HRMS; head, tail and arm "
                    "count all come from the name plus Supplementary Fig. 1.")
        return ("Enumerated structure reproduces the paper's printed exact "
                "mass to %s Da." % e["hrms_dev_Da"])

    qual = emit.emit_rows(
        items, recs, confidence=confidence, note=note,
        positional_ambiguity=lambda i, r: (
            r["enum"]["positional_ambiguity"] if "enum" in r else "0"),
        structure_assumed=lambda i, r: (
            "1" if r["kind"] == "external_ML_test" else "0"),
        value_is_bound=None, value_normalization_unit=None,
        within_model=None, label="linyx_41851281")

    out = []
    for iid, q in zip(items, qual):
        r = recs[iid]
        reps = [float(v) for v in vals[iid]]
        mean = sum(reps) / len(reps)
        row = {k: r[k] for k in
               ("kind", "name", "ml_id", "lipid_name", "smiles", "inchikey",
                "pool", "assay", "cargo", "system", "route", "dose",
                "timepoint", "formulation", "book", "sheet_name", "cells",
                "evidence")}
        row["measured_or_predicted"] = "measured"
        row["n_replicates"] = len(reps)
        for j, v in enumerate(reps, 1):
            row["rep%d" % j] = "%.6g" % v
        row["mean"] = "%.6g" % mean
        row["log10_mean"] = "%.5f" % math.log10(mean)
        row["source_tier"] = "1_publisher_source_data"
        row.update(q)
        out.append(row)

    cols = list(out[0].keys())
    with open(os.path.join(HERE, "05_values.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in out:
            w.writerow(r)
    print("05_values.csv: %d rows (%d in vitro, %d in vivo), caveats %d"
          % (len(out), sum(r["pool"] == "in_vitro" for r in out),
             sum(r["pool"] == "in_vivo" for r in out),
             emit.caveat_count(out)))
    return out


if __name__ == "__main__":
    main()
