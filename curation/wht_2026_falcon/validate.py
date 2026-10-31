#!/usr/bin/env python3
"""Independent check of 01_WT_2026_rows.csv. Reads the OUTPUT and the SOURCE and
compares them; it does not import build_rows, so a bug there cannot hide here.

The external anchor is Supplementary Table 9, which prints the initial
grid-search library size for each experiment. Iteration 0 of each dataset must
contain exactly that many formulations -- a number the authors published and the
repo never states, so agreement is real evidence the files were read correctly.
"""
import csv
import json
import os
import random
import statistics as st

from rdkit import Chem, RDLogger
from rdkit.Chem import Crippen, Descriptors as D, rdMolDescriptors as rd

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "repo", "datasets", "FALCON_Manuscript_Datasets")
ROWS = os.path.join(HERE, "01_WT_2026_rows.csv")

# Supplementary Table 9, "Initial Library Size" column
TABLE9_LIBRARY = {"0": 81, "1": 108, "2": 36, "3": 108, "4": 108}
FILES = {"0": "0_Norm_RAMOS_THP1_Single_Dual_Objective_4ITER_MC3_DSPC.csv",
         "1": "1_Norm_RAMOS_HEPG2_Single_Dual_Objective_4ITER_SM102_18PG.csv",
         "2": "2_Norm_THP1_HEPG2_Dual_Objective_5ITER_SM102_DSPC.csv",
         "3": "3_Norm_INS_Single_Objective_3ITER_SM102_DOTAP.csv",
         "4": "4_Norm_RAMOS_C2C12_DC24_3T3_Quad_Objective_3ITER_SM102_18PG.csv"}

fails, notes = [], []


def check(ok, msg):
    (notes if ok else fails).append(("PASS" if ok else "FAIL") + "  " + msg)


# Derived, not a literal: 67 curation columns + 5 recomputed descriptors + the
# quality columns in provenance/schema_quality.json. A frozen count here is the
# error class F-002 was logged for -- it fired on a correct file the moment the
# schema legitimately grew, which trains whoever hits it to edit the number.
NCOL = 72 + len(json.load(open(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "..",
    "provenance", "schema_quality.json")))["columns"])


def main():
    rows = list(csv.DictReader(open(ROWS)))
    check(len(rows) == 2544, f"row count {len(rows)} (expect 2544)")
    check(len(rows[0]) == NCOL, f"column count {len(rows[0])} (expect {NCOL})")

    # --- external anchor: published initial-library sizes -------------------
    for tag, want in TABLE9_LIBRARY.items():
        src = list(csv.DictReader(open(os.path.join(SRC, FILES[tag]), encoding="utf-8-sig")))
        got = sum(1 for r in src if r["Iter"] == "0")
        check(got == want,
              f"dataset {tag} iteration-0 library = {got}, Suppl Table 9 prints {want}")

    # --- z-score blocks, recomputed from the output ------------------------
    blocks = {}
    for r in rows:
        blocks.setdefault(r["Experiment_ID"], []).append(float(r["Experiment_value"]))
    # SEVEN blocks, one per source file, corrected from fourteen on 2026-10-30
    # (review finding F-003). Per-(file x cell type) standardisation forced the
    # authors' own nLnLE_Selectivity to identically zero.
    check(len(blocks) == 7, f"{len(blocks)} experiment blocks, one per source "
                            f"file (expect 7; it was 14 before F-003)")
    worst_mu = max(abs(st.mean(v)) for v in blocks.values())
    worst_sd = max(abs(st.pstdev(v) - 1) for v in blocks.values())
    check(worst_mu < 1e-9, f"max |block mean| = {worst_mu:.2e}")
    check(worst_sd < 1e-9, f"max |block sd - 1| = {worst_sd:.2e}")
    check(all(len(v) >= 2 for v in blocks.values()), "every block has >1 row")

    # --- composition -------------------------------------------------------
    # Dataset 1 prints its percentages to one decimal (e.g. 44.4/33.3/21.1/1.1),
    # so its rows sum to 99.9 or 100.1. That is the authors' rounding carried
    # through unchanged, not a curation error; anything beyond it would be.
    sums = [sum(float(r[f"{p}_molratio"]) for p in ("IL", "HL", "CHL", "PEG")) for r in rows]
    worst = max(abs(s - 100) for s in sums)
    check(worst <= 0.1 + 1e-9,
          f"molar percents sum to 100, max deviation {worst:.4f} (source rounds to 0.1)")
    np_r = [float(r["IL_to_nucleicacid_chargeratio"]) for r in rows]
    check(min(np_r) > 0, f"N/P range {min(np_r):.2f}-{max(np_r):.2f}")

    # --- structures --------------------------------------------------------
    uniq = {(p, r[f"{p}_name"], r[f"{p}_SMILES"])
            for r in rows for p in ("IL", "HL", "CHL", "PEG")}
    bad = []
    for p, name, smi in uniq:
        m = Chem.MolFromSmiles(smi)
        if m is None:
            bad.append((name, "unparseable")); continue
        if Chem.MolToSmiles(m) != smi:
            bad.append((name, "not canonical"))
        if len(Chem.GetMolFrags(m)) != 1:
            bad.append((name, "disconnected"))
        if sum(a.GetNumRadicalElectrons() for a in m.GetAtoms()):
            bad.append((name, "radical"))
    check(not bad, f"{len(uniq)} distinct components parse, canonical, single fragment, "
                   f"no radicals {bad if bad else ''}")
    # one name must not carry two structures
    collide = {}
    for p, name, smi in uniq:
        collide.setdefault((p, name), set()).add(smi)
    multi = {k: v for k, v in collide.items() if len(v) > 1}
    check(not multi, f"no component name maps to two SMILES {list(multi)}")

    # --- descriptors recomputed independently ------------------------------
    cache, mism = {}, 0
    for r in random.Random(42).sample(rows, 300):
        s = r["IL_SMILES"]
        if s not in cache:
            m = Chem.MolFromSmiles(s)
            cache[s] = {"Heavy.Atoms": m.GetNumHeavyAtoms(),
                        "Rotatable.Bonds": rd.CalcNumRotatableBonds(m),
                        "Nitrogen.Count": sum(1 for a in m.GetAtoms() if a.GetSymbol() == "N"),
                        "LogP": round(Crippen.MolLogP(m), 6),
                        "ExactMolWt": round(D.ExactMolWt(m), 6)}
        for k, v in cache[s].items():
            got = float(r[k]) if k in ("LogP", "ExactMolWt") else int(float(r[k]))
            if abs(got - v) > 1e-5:
                mism += 1
    check(mism == 0, f"descriptors agree with their own SMILES on 300 sampled rows "
                     f"({mism} mismatches)")

    # --- round-trip a sample of values back to the source ------------------
    src_cache = {}
    bad_rt = 0
    tagmap = {"URS": "0_Norm_RAMOS_Single_Objective_UniformRandomSampling_Ctrl_MC3_DSPC.csv",
              "E5": "5_Norm_RAMOS_HEPG2_Dual_Objective_4ITER_SM102_DSPC.csv", **FILES}
    colmap = {"RAMOS": "nLnLE_RAMOS", "THP1": "nLnLE_THP1", "HEPG2": "nLnLE_HEPG2",
              "INS1E": "nLnLE_INS", "DC24": "nLnLE_DC", "3T3": "nLnLE_3T3",
              "C2C12": "nLnLE_C2C12"}
    # Which cell-type columns each file contributes. Needed now that a block is
    # a FILE rather than a (file x cell type) pair: the z is pooled over all of
    # a file's columns, so the round-trip has to pool them the same way.
    OUTCOLS = {"URS": ["nLnLE_RAMOS"],
               "E0": ["nLnLE_RAMOS", "nLnLE_THP1"],
               "E1": ["nLnLE_RAMOS", "nLnLE_HEPG2"],
               "E2": ["nLnLE_THP1", "nLnLE_HEPG2"],
               "E3": ["nLnLE_INS"],
               "E4": ["nLnLE_RAMOS", "nLnLE_DC", "nLnLE_3T3", "nLnLE_C2C12"],
               "E5": ["nLnLE_HEPG2", "nLnLE_RAMOS"]}
    for eid, vals in blocks.items():
        _, _, tag = eid.split("_", 2)
        f = tagmap[tag[1:] if tag.startswith("E") and tag != "E5" else tag]
        if f not in src_cache:
            src_cache[f] = list(csv.DictReader(open(os.path.join(SRC, f),
                                                    encoding="utf-8-sig")))
        src = src_cache[f]
        raw = [float(r[c]) for c in OUTCOLS[tag] for r in src]
        mu, sd = st.mean(raw), st.pstdev(raw)
        want = sorted((v - mu) / sd for v in raw)
        if len(want) != len(vals) or max(abs(a - b) for a, b
                                         in zip(want, sorted(vals))) > 1e-12:
            bad_rt += 1
    check(bad_rt == 0, f"all 7 blocks round-trip to the source CSV ({bad_rt} differ)")

    # THE CHECK F-003 EXISTS FOR, re-done here against the emitted rows rather
    # than inside the builder: the stored z difference between two cell types in
    # one file must reproduce the authors' printed Selectivity column up to a
    # single per-file constant. Under per-column z-scoring no such constant
    # exists, so this fails if the normalisation unit is ever changed back.
    sel_bad, sel_done = 0, 0
    byfid = {}
    for r in rows:
        byfid.setdefault((r["Experiment_ID"], r["Formulation_ID"]), {})[
            r["Model_type"]] = float(r["Experiment_value"])
    CELL = {"nLnLE_RAMOS": "Ramos", "nLnLE_THP1": "THP-1", "nLnLE_HEPG2": "HepG2",
            "nLnLE_INS": "INS-1E", "nLnLE_DC": "DC2.4", "nLnLE_3T3": "NIH-3T3",
            "nLnLE_C2C12": "C2C12"}
    for tag, cols in OUTCOLS.items():
        if len(cols) != 2:
            continue
        f = tagmap[tag[1:] if tag.startswith("E") and tag != "E5" else tag]
        if f not in src_cache:
            src_cache[f] = list(csv.DictReader(open(os.path.join(SRC, f),
                                                    encoding="utf-8-sig")))
        src = src_cache[f]
        selcol = next((c for c in src[0] if c.endswith("Selectivity")), None)
        if not selcol:
            continue
        raw = [float(r[c]) for c in cols for r in src]
        sd = st.pstdev(raw)
        a, b = CELL[cols[0]], CELL[cols[1]]
        pairs = [(v[a] - v[b]) for (e, _), v in byfid.items()
                 if e == f"WT_2026_{tag}" and a in v and b in v]
        printed = [float(r[selcol]) for r in src]
        if len(pairs) != len(printed):
            sel_bad += 1
            continue
        # The printed column's SIGN convention is the authors', not ours:
        # DATASETS lists file 5's columns HepG2-first while its Selectivity is
        # Ramos - HepG2. So accept whichever orientation matches, as the
        # builder's own solver does -- what is being tested is that a SINGLE
        # per-file constant recovers the difference at all, which is false
        # under per-column z-scoring in either orientation.
        best = min(max(abs(x - y) for x, y in zip(sorted(sgn * p * sd for p in pairs),
                                                  sorted(printed)))
                   for sgn in (1, -1))
        if best > 1e-6:
            sel_bad += 1
        sel_done += 1
    check(sel_bad == 0, f"stored z differences reproduce the authors' printed "
                        f"nLnLE_Selectivity up to one constant per file "
                        f"({sel_done} files checked, {sel_bad} fail)")

    # --- things to be aware of, not failures -------------------------------
    comp_key = lambda r: (r["IL_name"], r["HL_name"],
                          r["IL_molratio"], r["HL_molratio"], r["CHL_molratio"],
                          r["PEG_molratio"], r["IL_to_nucleicacid_chargeratio"])
    per_exp = {}
    for r in rows:
        per_exp.setdefault(r["Experiment_ID"], []).append(comp_key(r))
    dups = {k: len(v) - len(set(v)) for k, v in per_exp.items() if len(v) != len(set(v))}

    out = {"failures": fails, "passes": notes,
           "repeated_compositions_within_an_experiment": dups,
           "known_duplicate_in_the_source_table": {
               "rows": ["F1109 (iter 1)", "F1120 (iter 2)"],
               "file": "1_Norm_RAMOS_HEPG2_..._SM102_18PG.csv",
               "what": "same composition; nLnLE_RAMOS byte-identical "
                       "(2.483377939) while nLnLE_HEPG2 differs (3.237 vs "
                       "4.414). Published that way by the authors - the Ramos "
                       "value looks carried forward between iterations. Kept "
                       "rather than dropped; it is their data, not a read "
                       "error. Shows up as audit C4's one exact duplicate row.",
           },
           "note": "repeats are real: FALCON retests reference/control formulations "
                   "across iterations, so a composition can appear more than once "
                   "with different measured values."}
    json.dump(out, open(os.path.join(HERE, "04_validation.json"), "w"), indent=1)

    for line in notes + fails:
        print(" ", line)
    print(f"\nrepeated compositions within an experiment: {dups or 'none'}")
    print(f"\n{len(fails)} failures, {len(notes)} passes")
    return 1 if fails else 0


if __name__ == "__main__":
    raise SystemExit(main())
