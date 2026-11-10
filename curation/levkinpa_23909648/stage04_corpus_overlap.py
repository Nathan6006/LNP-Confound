#!/usr/bin/env python3
"""Stage 04 -- is any of this already in the corpus?

WHY. Siblings in this run found 360/489, 246/580, 15/116 and 11/104 of their
members already present by full InChIKey, under KZ_2016, NC_2024, KW_2014 and
AA_2008. Twice this became the strongest external check in a folder, because
the corpus's own lipidoid names encode chain length, so an overlap is an
independent witness on the head/tail assignment.

REPORT, DO NOT ACT. This stage writes a report and changes nothing.

Three comparisons, because they fail differently:
  (a) FULL InChIKey -- the same compound, stereo included.
  (b) SKELETON InChIKey block -- the same connectivity, stereo ignored. This is
      the comparison that was blind to epicholesterol twice (F-013), so a
      skeleton hit with no full hit is reported as a STEREO DIFFERENCE, loudly.
  (c) IL_name -- a label collision, which is NOT a structural claim: 288 IL
      names already carry different structures in different experiments
      (CLAUDE.md section 3), so a name hit is a flag to check, never a merge.

Plus a substructure sweep for the 2-(pyrrolidin-1-yl)ethylamine head, which
would catch the same chemistry entering the corpus under a different label.

Run:  python3 -I stage04_corpus_overlap.py
"""
import csv
import json
import os
import sys

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import inchi

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.abspath(os.path.join(HERE, "..", "..", "data"))

HEAD_SMARTS = Chem.MolFromSmarts("[NX3;!R]CC[NX3;R]1CCCC1")


def main():
    ours = list(csv.DictReader(open(os.path.join(HERE, "03_enumerated.csv"))))
    crude = list(csv.DictReader(
        open(os.path.join(HERE, "03_not_determinable.csv"))))

    rep = dict(
        _what="Stage 04. Corpus-overlap report. READ-ONLY; nothing is acted "
              "on here.",
        corpus_files={}, our_structures=[], name_collisions=[],
        head_substructure_sweep={})

    frames = {}
    for pool, fn in [("vitro", "data_vitro.csv"), ("vivo", "data_vivo.csv")]:
        df = pd.read_csv(os.path.join(DATA, fn), low_memory=False)
        frames[pool] = df
        rep["corpus_files"][pool] = dict(
            file=fn, rows=int(df.shape[0]), cols=int(df.shape[1]),
            experiments=int(df["Experiment_ID"].nunique()))

    # --- build the corpus's InChIKey index once ---------------------------
    key_index = {}          # full key -> {(pool, experiment_id, IL_name)}
    skel_index = {}
    smiles_seen = {}
    heads = []
    for pool, df in frames.items():
        sub = df[["Experiment_ID", "IL_name", "IL_SMILES"]].dropna(
            subset=["IL_SMILES"])
        for smi, grp in sub.groupby("IL_SMILES"):
            if smi in smiles_seen:
                k, sk, has_head = smiles_seen[smi]
            else:
                m = Chem.MolFromSmiles(smi)
                if m is None:
                    smiles_seen[smi] = (None, None, False)
                    continue
                k = inchi.MolToInchiKey(m)
                sk = k.split("-")[0] if k else None
                has_head = bool(HEAD_SMARTS) and m.HasSubstructMatch(
                    HEAD_SMARTS)
                smiles_seen[smi] = (k, sk, has_head)
            if not k:
                continue
            for _, r in grp.iterrows():
                tag = (pool, str(r["Experiment_ID"]), str(r["IL_name"]))
                key_index.setdefault(k, set()).add(tag)
                skel_index.setdefault(sk, set()).add(tag)
                if has_head:
                    heads.append(dict(pool=pool,
                                      experiment=str(r["Experiment_ID"]),
                                      IL_name=str(r["IL_name"]),
                                      IL_SMILES=smi))

    rep["corpus_files"]["unique_IL_SMILES_indexed"] = len(
        [s for s in smiles_seen if smiles_seen[s][0]])

    # --- (a) and (b) ------------------------------------------------------
    for r in ours:
        full = r["inchikey_full"]
        skel = r["inchikey_skeleton"]
        fh = sorted(key_index.get(full, []))
        sh = sorted(skel_index.get(skel, []))
        rep["our_structures"].append(dict(
            IL_name=r["IL_name"], IL_SMILES=r["IL_SMILES"],
            inchikey_full=full, inchikey_skeleton=skel,
            full_key_hits=[list(t) for t in fh],
            n_full_key_hits=len(fh),
            skeleton_only_hits=[list(t) for t in sh
                                if t not in set(fh)],
            verdict=("ALREADY IN CORPUS" if fh else
                     "SKELETON MATCH ONLY -- a stereo difference, which is "
                     "the F-013 blind spot; investigate before merging"
                     if sh else "NOT IN CORPUS")))

    # --- (c) name collisions, over every label this paper names ----------
    our_labels = ([r["IL_name"] for r in ours]
                  + [r["IL_name"] for r in crude] + ["L2K"])
    for pool, df in frames.items():
        names = set(df["IL_name"].dropna().astype(str))
        for lab in our_labels:
            if lab in names:
                hit = df[df["IL_name"].astype(str) == lab]
                rep["name_collisions"].append(dict(
                    pool=pool, label=lab,
                    experiments=sorted(set(hit["Experiment_ID"].astype(str))),
                    n_rows=int(hit.shape[0]),
                    distinct_SMILES=sorted(
                        set(hit["IL_SMILES"].dropna().astype(str)))[:3],
                    meaning="A NAME COLLISION IS NOT A STRUCTURAL CLAIM. 288 "
                            "IL names already carry different structures in "
                            "different experiments (CLAUDE.md section 3). "
                            "Never key on IL_name across experiments."))

    # --- head substructure sweep -----------------------------------------
    rep["head_substructure_sweep"] = dict(
        smarts="[NX3;!R]CC[NX3;R]1CCCC1",
        what="an N-alkylated 2-(pyrrolidin-1-yl)ethylamine, i.e. this paper's "
             "amine-1 head after alkylation. A hit would mean the same "
             "chemistry is in the corpus under a different label.",
        n_hits=len(heads),
        hits=heads[:40])

    with open(os.path.join(HERE, "04_corpus_overlap.json"), "w") as fh:
        json.dump(rep, fh, indent=1)

    print("OK  corpus-overlap report written")
    for pool, d in rep["corpus_files"].items():
        if isinstance(d, dict):
            print(f"    {pool}: {d['rows']} rows, {d['experiments']} experiments")
    print(f"    unique IL_SMILES indexed: "
          f"{rep['corpus_files']['unique_IL_SMILES_indexed']}")
    for s in rep["our_structures"]:
        print(f"    {s['IL_name']:<7} full-key hits {s['n_full_key_hits']}, "
              f"skeleton-only {len(s['skeleton_only_hits'])}  -> {s['verdict']}")
    print(f"    name collisions: {len(rep['name_collisions'])}")
    for c in rep["name_collisions"]:
        print(f"      {c['pool']}/{c['label']}: {c['n_rows']} rows in "
              f"{c['experiments']}")
    print(f"    amine-1 head substructure in corpus: "
          f"{rep['head_substructure_sweep']['n_hits']} rows")
    return 0


if __name__ == "__main__":
    sys.exit(main())
