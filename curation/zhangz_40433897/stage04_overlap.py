#!/usr/bin/env python3
"""Stage 04 - overlap with data/ : full InChIKey, skeleton key, IL_name, and a
substructure / near-miss search.

Nothing in this paper has a PINNED structure (see DECISIONS.md D-1), so the
exact-key arms of this check are run on the PROPOSED structures and are
reported as such. The arms that do not depend on our proposal - the IL_name
search, the tocopherol substructure search, and the benchmark check - are the
ones that carry weight.

The benchmark arm is the free check the brief describes: SM-102 appears in
Figure S8b, so if OUR SM-102 string reproduces the corpus on all three keys at
once, our chemistry handling is right independently of this paper.
"""
import json
import os
import sys

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import inchi
from rdkit.Chem import DataStructs
from rdkit.Chem import rdFingerprintGenerator

RDLogger.DisableLog("rdApp.*")
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "..", "data")

# Tocopherol (vitamin E) core, as the substructure arm's query.
TOCOPHEROL_CORE = "CC1=C(C)C2=C(CCC(C)(CCCC(C)CCCC(C)CCCC(C)C)O2)C(C)=C1O"
CHROMANOL = "CC1(CCc2c(C)c(O)c(C)c(C)c2O1)C"

# Benchmarks that should already be in the corpus.
#
# These are NOT written from memory. My first attempt was, and all three were
# wrong: SM-102 came out C43H85NO5 against the true C44H87NO5, MC3 was one CH2
# short on a tail, and ALC-0315 had the ester connectivity reversed. All three
# parsed, canonicalised and carried no radicals - only the InChIKey saw it.
# They are taken instead from provenance/reagent_registry.py, the project's own
# pinned name -> structure authority (C11).
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "provenance"))
from reagent_registry import REGISTRY           # noqa: E402

BENCH = {k: REGISTRY[k][0] for k in ("SM-102", "DLin-MC3-DMA", "ALC-0315")}

PROPOSED = json.load(open(os.path.join(HERE, "03_proposed_structures.json")))


def keys(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None, None
    k = inchi.MolToInchiKey(m)
    return k, k.split("-")[0]


def main():
    vit = pd.read_csv(os.path.join(DATA, "data_vitro.csv"), low_memory=False)
    viv = pd.read_csv(os.path.join(DATA, "data_vivo.csv"), low_memory=False)
    corp = pd.concat([vit.assign(_pool="vitro"), viv.assign(_pool="vivo")],
                     ignore_index=True)

    uniq = corp[["IL_SMILES", "IL_name", "Experiment_ID", "_pool"]].dropna(
        subset=["IL_SMILES"]).drop_duplicates("IL_SMILES")
    out = {"corpus": {"vitro_rows": len(vit), "vivo_rows": len(viv),
                      "unique_IL_SMILES": int(uniq.shape[0]),
                      "experiments_vitro": int(vit.Experiment_ID.nunique()),
                      "experiments_vivo": int(viv.Experiment_ID.nunique())}}

    full, skel = {}, {}
    bad = 0
    for smi, name, exp, pool in uniq.itertuples(index=False):
        f, s = keys(smi)
        if f is None:
            bad += 1
            continue
        full.setdefault(f, []).append((name, exp, pool))
        skel.setdefault(s, []).append((name, exp, pool))
    out["corpus"]["unparseable_IL_SMILES"] = bad

    # ---- benchmark arm -------------------------------------------------
    bench = {}
    for nm, smi in BENCH.items():
        f, s = keys(smi)
        hits_f = full.get(f, [])
        hits_s = skel.get(s, [])
        nm_hits = sorted(set(
            corp.loc[corp.IL_name.astype(str).str.replace("-", "", regex=False)
                     .str.upper() == nm.replace("-", "").upper(),
                     "Experiment_ID"]))
        bench[nm] = {
            "full_inchikey": f,
            "full_key_hit": bool(hits_f),
            "n_experiments_full_key": len(set(e for _, e, _ in hits_f)),
            "skeleton_key_hit": bool(hits_s),
            "n_experiments_skeleton": len(set(e for _, e, _ in hits_s)),
            "IL_name_hit_experiments": nm_hits,
            "three_of_three": bool(hits_f and hits_s and nm_hits)}
    out["benchmarks"] = bench
    out["benchmarks_summary"] = "%d/3 hit on all three keys at once" % sum(
        1 for b in bench.values() if b["three_of_three"])

    # ---- this paper's proposed structures ------------------------------
    paper = {}
    for rec in PROPOSED["lipids"]:
        e = {"formula_from_HRMS": rec["formula_ion"],
             "core_formula": rec["core_formula"],
             "structure_pinned": rec["structure_pinned"],
             "proposed_IL_SMILES": rec.get("proposed_IL_SMILES")}
        smi = rec.get("proposed_IL_SMILES")
        if smi:
            f, s = keys(smi)
            e["full_inchikey"] = f
            e["full_key_hit"] = bool(full.get(f))
            e["skeleton_key_hit"] = bool(skel.get(s))
            e["skeleton_hits"] = sorted(set(x[1] for x in skel.get(s, [])))
        e["IL_name_hit"] = sorted(set(
            corp.loc[corp.IL_name.astype(str).str.upper() ==
                     rec["name"].upper(), "Experiment_ID"]))
        paper[rec["name"]] = e
    out["this_paper"] = paper

    # ---- substructure / near-miss arm ----------------------------------
    sub = {}
    for label, q in (("tocopherol_full", TOCOPHEROL_CORE),
                     ("chromanol_ring_system", CHROMANOL)):
        patt = Chem.MolFromSmiles(q)
        n = 0
        hits = []
        for smi, name, exp, pool in uniq.itertuples(index=False):
            m = Chem.MolFromSmiles(smi)
            if m is not None and m.HasSubstructMatch(patt):
                n += 1
                hits.append({"IL_name": name, "experiment": exp, "pool": pool})
        sub[label] = {"n_corpus_lipids_containing_it": n, "hits": hits[:20]}
    out["substructure_search"] = sub

    # nearest Tanimoto neighbour of the proposed VN2C2, as a near-miss probe
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    probe = next((r for r in PROPOSED["lipids"]
                  if r.get("proposed_IL_SMILES")), None)
    if probe:
        pm = Chem.MolFromSmiles(probe["proposed_IL_SMILES"])
        pf = gen.GetFingerprint(pm)
        best = []
        for smi, name, exp, pool in uniq.itertuples(index=False):
            m = Chem.MolFromSmiles(smi)
            if m is None:
                continue
            t = DataStructs.TanimotoSimilarity(pf, gen.GetFingerprint(m))
            best.append((t, name, exp, pool))
        best.sort(reverse=True)
        out["near_miss"] = {
            "probe": probe["name"],
            "probe_smiles": probe["proposed_IL_SMILES"],
            "fingerprint": "Morgan r=2, 2048 bits",
            "top10": [{"tanimoto": round(t, 4), "IL_name": n,
                       "experiment": e, "pool": p} for t, n, e, p in best[:10]]}

    json.dump(out, open(os.path.join(HERE, "02_corpus_overlap.json"), "w"),
              indent=1)
    print(json.dumps({k: v for k, v in out.items()
                      if k != "substructure_search"}, indent=1)[:4000])
    print("\nsubstructure:", json.dumps(
        {k: v["n_corpus_lipids_containing_it"]
         for k, v in out["substructure_search"].items()}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
