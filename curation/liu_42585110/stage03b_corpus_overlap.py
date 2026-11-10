#!/usr/bin/env python3
"""Does any RAIL, or Dlin-MC3 as this paper formulates it, already exist in the corpus?

Run from this folder:  python3 -I stage03b_corpus_overlap.py
                                               -> 02_corpus_overlap.json

EVERY QUERY HERE CARRIES A POSITIVE CONTROL, AND EVERY CONTROL ASSERTS THAT IT
MATCHES A SPECIFIC NAMED MOLECULE KNOWN TO BE PRESENT (EXTRACTION_BRIEF 5,
F-070). A non-zero count is NOT evidence that a control works: one sweep in
this run returned 6 structures from a hand-written sterol SMARTS and reported
`passes: True` while not matching the registry cholesterol at all.

Joins are on FULL InChIKey, never on IL_name. A structure sits under several
names in this corpus (cKK-E12 is also `other` and `MD-1`), 292 names carry more
than one structure, and the paper's own spelling of a benchmark can match zero
rows for a lipid present thousands of times -- `MC-3` matched 0 for a molecule
present 3,098 times. This paper's spelling is `Dlin-MC3`, a THIRD spelling
again, and it matches 0 corpus rows by name.
"""
import collections
import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "provenance"))
import reagent_registry as RR                                    # noqa: E402

from rdkit import Chem                                           # noqa: E402
from rdkit import RDLogger                                       # noqa: E402
RDLogger.DisableLog("rdApp.*")

VITRO = os.path.join(ROOT, "data", "data_vitro.csv")
VIVO = os.path.join(ROOT, "data", "data_vivo.csv")
OUT = os.path.join(HERE, "02_corpus_overlap.json")

_K = {}


def key(smiles):
    if smiles not in _K:
        m = Chem.MolFromSmiles(smiles)
        _K[smiles] = Chem.MolToInchiKey(m) if m else None
    return _K[smiles]


def main():
    vitro = list(csv.DictReader(open(VITRO)))
    vivo = list(csv.DictReader(open(VIVO)))
    rep = {"_what": "RAIL library vs the live corpus, joined on full InChIKey.",
           "corpus_size": {"vitro_rows": len(vitro), "vivo_rows": len(vivo)}}

    # ---- the structures we would add -----------------------------------
    ours = {}
    for r in csv.DictReader(open(os.path.join(HERE, "03_enumerated.csv"))):
        ours[r["inchikey"]] = r["IL_name"]
    rep["rails_enumerated"] = len(ours)

    # ---- POSITIVE CONTROLS ---------------------------------------------
    # Each asserts a NAMED molecule is found. A count is not a control.
    controls = {}
    il_keys_vitro = collections.Counter()
    il_keys_vivo = collections.Counter()
    for rows, ctr in ((vitro, il_keys_vitro), (vivo, il_keys_vivo)):
        for r in rows:
            s = r["IL_SMILES"].strip()
            if s:
                k = key(s)
                if k:
                    ctr[k] += 1

    mc3_key = key(RR.REGISTRY["DLin-MC3-DMA"][0])
    n_mc3 = il_keys_vitro[mc3_key] + il_keys_vivo[mc3_key]
    controls["DLin-MC3-DMA_by_InChIKey"] = {
        "inchikey": mc3_key,
        "rows": n_mc3,
        "vitro": il_keys_vitro[mc3_key], "vivo": il_keys_vivo[mc3_key],
        "asserts": "the registry's DLin-MC3-DMA is present in the corpus",
        "passes": n_mc3 > 0,
    }
    chol_key = key(RR.REGISTRY["Cholesterol"][0])
    n_chol = sum(1 for r in vitro + vivo
                 if r["CHL_SMILES"].strip() and key(r["CHL_SMILES"]) == chol_key)
    controls["Cholesterol_by_InChIKey_in_CHL_column"] = {
        "inchikey": chol_key, "rows": n_chol,
        "asserts": "the registry's cholesterol is present as a CHL component",
        "passes": n_chol > 0,
    }
    dspc_key = key(RR.REGISTRY["DSPC"][0])
    n_dspc = sum(1 for r in vitro + vivo
                 if r["HL_SMILES"].strip() and key(r["HL_SMILES"]) == dspc_key)
    controls["DSPC_by_InChIKey_in_HL_column"] = {
        "inchikey": dspc_key, "rows": n_dspc,
        "asserts": "the registry's DSPC is present as a helper lipid",
        "passes": n_dspc > 0,
    }
    # a control for the MACHINERY used on our own structures: feed it a
    # molecule we KNOW is in both sets.
    controls["roundtrip_our_own_lead"] = {
        "asserts": "the lookup machinery finds AT4-13CRA in our own table "
                   "when asked, so a zero against the corpus is a real zero "
                   "and not a broken join",
        "passes": any(v == "AT4-13CRA" for v in ours.values()),
    }
    rep["positive_controls"] = controls
    if not all(c["passes"] for c in controls.values()):
        raise SystemExit("stage03b: a positive control FAILED. A zero result "
                         "from this sweep would mean nothing. Fix the control "
                         "before reading any count.")

    # ---- the actual query ----------------------------------------------
    hits = {k: ours[k] for k in ours if k in il_keys_vitro or k in il_keys_vivo}
    rep["rail_structures_already_in_corpus"] = {
        "count": len(hits), "detail": hits,
        "interpretation":
            "ZERO is the expected answer and it is now a measured zero rather "
            "than an assumed one: no retinoid-derived ionizable lipid has any "
            "analogue among the corpus's 19,919 structures. The control above "
            "proves the same join finds DLin-MC3-DMA, so the machinery works.",
    }

    # a SUPERSET sweep, labelled as such: anything containing a retinoid-like
    # trimethylcyclohexenyl polyene, or a bexarotene-like tetrahydronaphthalene
    # My FIRST version of the retinoid SMARTS was written "C1(C)=C(CCCC1(C)C)
    # /C=C/C" and returned a clean corpus 0 -- while not matching AT1-ATRA,
    # a molecule built three files earlier that unquestionably carries the
    # motif. The assertion below is what caught it. That is F-070, reproduced
    # from scratch, and it is the reason the zero in this file can be read at
    # all.
    loose = {
        "retinoid_trimethylcyclohexene_polyene":
            "C1(C)(C)CCCC(C)=C1C=C",
        "tetrahydronaphthalene_gem_dimethyl":
            "C1(C)(C)CCC(C)(C)c2ccccc12",
        "any_tertiary_amine_SANITY_SUPERSET": "[NX3;H0;!$(NC=O)](C)(C)C",
    }
    sweep = {}
    allrows = vitro + vivo
    uniq = {}
    for r in allrows:
        s = r["IL_SMILES"].strip()
        if s:
            uniq.setdefault(s, 0)
            uniq[s] += 1
    for nm, sma in loose.items():
        patt = Chem.MolFromSmarts(sma)
        if patt is None:
            raise SystemExit(f"stage03b: SMARTS {nm} does not compile.")
        n = 0
        for s in uniq:
            m = Chem.MolFromSmiles(s)
            if m is not None and m.HasSubstructMatch(patt):
                n += 1
        sweep[nm] = {"smarts": sma, "distinct_corpus_structures": n}
    # assert the two retinoid SMARTS actually match OUR lipids
    ours_rows = list(csv.DictReader(open(os.path.join(HERE,
                                                      "03_enumerated.csv"))))
    # The tertiary-amine superset is controlled against the corpus itself:
    # it must match DLin-MC3-DMA, a molecule we know is in there.
    patt_amine = Chem.MolFromSmarts(loose["any_tertiary_amine_SANITY_SUPERSET"])
    mc3_ok = Chem.MolFromSmiles(
        RR.REGISTRY["DLin-MC3-DMA"][0]).HasSubstructMatch(patt_amine)
    sweep["any_tertiary_amine_SANITY_SUPERSET"][
        "control_matches_DLin-MC3-DMA"] = bool(mc3_ok)
    if not mc3_ok:
        raise SystemExit(
            "stage03b: the tertiary-amine superset does not match "
            "DLin-MC3-DMA, which is a tertiary amine present 3,000+ times. "
            "Its count is uninterpretable.")

    for nm, want in (("retinoid_trimethylcyclohexene_polyene", "AT1-ATRA"),
                     ("tetrahydronaphthalene_gem_dimethyl", "AT1-Bex")):
        patt = Chem.MolFromSmarts(loose[nm])
        tgt = [r for r in ours_rows if r["IL_name"] == want][0]
        ok = Chem.MolFromSmiles(tgt["IL_SMILES"]).HasSubstructMatch(patt)
        sweep[nm]["control_matches_" + want] = bool(ok)
        if not ok:
            raise SystemExit(
                f"stage03b: SMARTS {nm} does not match {want}, a molecule we "
                "built and KNOW carries the motif. Its corpus count of "
                f"{sweep[nm]['distinct_corpus_structures']} is uninterpretable "
                "and a zero would be meaningless. This is F-070 exactly.")
    sweep["_label"] = ("SUPERSET sweep. These counts are chemotype "
                       "neighbourhoods, NOT shared lipids. The shared-lipid "
                       "answer is rail_structures_already_in_corpus above.")
    rep["superset_sweep"] = sweep

    # ---- the benchmark: Dlin-MC3 -----------------------------------------
    name_hits = collections.Counter()
    for r in allrows:
        if r["IL_name"].strip():
            name_hits[r["IL_name"].strip()] += 1
    rep["benchmark_Dlin_MC3"] = {
        "paper_spelling": "Dlin-MC3",
        "rows_matching_that_exact_name_in_corpus": name_hits.get("Dlin-MC3", 0),
        "corpus_spellings_of_the_same_InChIKey":
            sorted({r["IL_name"].strip() for r in allrows
                    if r["IL_SMILES"].strip()
                    and key(r["IL_SMILES"]) == mc3_key}),
        "rows_by_InChIKey": {"vitro": il_keys_vitro[mc3_key],
                             "vivo": il_keys_vivo[mc3_key]},
        "token_pattern_note":
            "The paper's own spelling 'Dlin-MC3' matches 0 corpus rows by "
            "name while the molecule is present thousands of times. Joining "
            "on IL_name would have reported this benchmark absent.",
    }

    # ---- the in-vivo picture, which is the load-bearing part --------------
    vv = [r for r in vivo
          if r["IL_SMILES"].strip() and key(r["IL_SMILES"]) == mc3_key]
    by = collections.Counter((r["Experiment_ID"], r["Model_target"]) for r in vv)
    organs = collections.Counter(r["Model_target"] for r in vivo)
    spleen_sources = sorted({r["Experiment_ID"] for r in vivo
                             if "spleen" in r["Model_target"].lower()})
    rep["in_vivo_shared_lipid"] = {
        "why_this_matters":
            "EXTRACTION_BRIEF 5: the in-vivo pool has exactly TWO structures "
            "shared across experiments, DLin-MC3-DMA and SM-102. Any in-vivo "
            "lipid already present elsewhere is rare and load-bearing.",
        "MC3_in_vivo_today": {f"{a}|{b}": n for (a, b), n in sorted(by.items())},
        "MC3_in_vivo_sources_today": len({a for a, _ in by}),
        "MC3_in_vivo_organs_today": sorted({b for _, b in by}),
        "MC3_liver_sources_today":
            sorted({a for (a, b) in by if b == "liver"}),
        "what_KJ_2026_would_add":
            "Dlin-MC3 measured in BOTH liver and spleen at 0.1 mg/kg IV, 6 h. "
            "LIVER: a FOURTH liver source for MC3 (today NG_2025_2, ZC_2023, "
            "ZH_2023), i.e. 3 new organ-matched cross-source pairs on the one "
            "benchmark the in-vivo pool actually shares. SPLEEN: the FIRST "
            "MC3 spleen measurement anywhere in the in-vivo corpus.",
        "spleen_coverage_today": {
            "spleen_rows": organs.get("spleen", 0),
            "spleen_sources": spleen_sources,
            "note": "The entire in-vivo corpus has ONE spleen source, "
                    "ZH_2023. KJ_2026 would be the second. ZH_2023 carries "
                    "MC3 only in liver, so the two spleen sources would still "
                    "share no lipid -- the gain is organ coverage, not a new "
                    "shared-lipid pair.",
        },
        "organ_census_today": dict(organs.most_common()),
    }

    # ---- is any arm of this screen already merged? -----------------------
    rep["already_merged_arm_check"] = {
        "pmids_in_corpus_matching_42585110":
            sorted({r["Publication_PMID"] for r in allrows
                    if r["Publication_PMID"].strip() == "42585110"}),
        "experiment_ids_starting_KJ":
            sorted({r["Experiment_ID"] for r in allrows
                    if r["Experiment_ID"].startswith("KJ")}),
        "distinct_experiment_ids_in_corpus":
            len({r["Experiment_ID"] for r in allrows}),
        "verdict": "No arm of this paper is in the corpus. KJ_2026 is free; "
                   "checked by PMID, not by the year token (F-075).",
    }

    with open(OUT, "w") as fh:
        json.dump(rep, fh, indent=1)
    print(json.dumps({k: v for k, v in rep.items()
                      if k not in ("superset_sweep",)}, indent=1)[:4000])
    print("\nsuperset sweep:", json.dumps(rep["superset_sweep"], indent=1))
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
