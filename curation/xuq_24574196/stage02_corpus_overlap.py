#!/usr/bin/env python3
"""Corpus overlap for the four enumerated lipidoids.

Four independent keys, because each is blind to something the others see:
  full InChIKey      identity including stereo (there is none here)
  skeleton InChIKey  connectivity only -- the block that F-013's epicholesterol
                     shared with cholesterol, so it is a NEAR-MISS detector,
                     not an identity test
  IL_name            the corpus's own label. A LEADING NUMBER TOKEN IS NOT AN
                     IDENTITY: `1-O16B` vs a corpus `1_O16B` would be a
                     coincidence of combinatorial labelling until the SMILES
                     agree. Reported, never trusted.
  substructure       the amine-1 core and the arm, searched separately, so a
                     corpus lipid built on the SAME head or the SAME tail
                     surfaces even when the whole molecule differs. This is
                     what gave a sibling ground truth its SI had withheld.
"""
import csv
import json
import os
import sys

from rdkit import Chem, RDLogger
from rdkit.Chem import inchi

RDLogger.DisableLog("rdApp.*")
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..")
CORPUS = [os.path.join(ROOT, "data", "data_vitro.csv"),
          os.path.join(ROOT, "data", "data_vivo.csv")]

# amine 1's core as it appears IN a lipidoid: a dimethylaminopropyl on a
# tertiary N bearing two 3-oxopropyl arms.
Q_HEAD = Chem.MolFromSmarts("CN(C)CCCN(CCC(=O)O[#6])CCC(=O)O[#6]")
Q_DISULFIDE_ARM = Chem.MolFromSmarts("C(=O)OCCSS[CX4]")
Q_ANY_SS = Chem.MolFromSmarts("[#6]SS[#6]")


def keys(smiles):
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None, None
    k = inchi.MolToInchiKey(m)
    return k, k.split("-")[0]


def main():
    ours = list(csv.DictReader(open(os.path.join(HERE, "03_enumerated.csv"))))
    ours_keys = {}
    for r in ours:
        full, skel = keys(r["IL_SMILES"])
        ours_keys[r["IL_name"]] = (full, skel)

    seen_full, seen_skel = {}, {}
    name_tokens = {}
    sub_head, sub_ssarm, sub_ss = [], [], []
    n_rows = 0
    uniq = {}
    for path in CORPUS:
        pool = os.path.basename(path)
        for row in csv.DictReader(open(path)):
            n_rows += 1
            smi = (row.get("IL_SMILES") or "").strip()
            nm = (row.get("IL_name") or "").strip()
            exp = row.get("Experiment_ID", "")
            if nm.split("_")[0] in {"1"} or nm.startswith("1-"):
                name_tokens.setdefault(nm, set()).add(exp)
            if not smi or smi in uniq:
                if smi in uniq:
                    uniq[smi][1].add(exp)
                continue
            uniq[smi] = (nm, {exp}, pool)

    # structural passes over UNIQUE structures only
    parsed = 0
    for smi, (nm, exps, pool) in uniq.items():
        m = Chem.MolFromSmiles(smi)
        if m is None:
            continue
        parsed += 1
        k = inchi.MolToInchiKey(m)
        seen_full.setdefault(k, []).append((nm, sorted(exps), pool))
        seen_skel.setdefault(k.split("-")[0], []).append((nm, sorted(exps), pool))
        if m.HasSubstructMatch(Q_HEAD):
            sub_head.append((nm, sorted(exps), pool, smi))
        if m.HasSubstructMatch(Q_DISULFIDE_ARM):
            sub_ssarm.append((nm, sorted(exps), pool, smi))
        if m.HasSubstructMatch(Q_ANY_SS):
            sub_ss.append((nm, sorted(exps), pool))

    out = {
        "_what": "corpus overlap for curation/xuq_24574196, four keys",
        "corpus_rows_scanned": n_rows,
        "unique_IL_SMILES": len(uniq),
        "unique_parsed": parsed,
        "exact_full_inchikey_hits": {},
        "skeleton_inchikey_hits": {},
        "name_token_hits_NOT_AN_IDENTITY": {k: sorted(v) for k, v in
                                           sorted(name_tokens.items())},
        "substructure_amine1_head": sub_head,
        "substructure_disulfide_ethyl_ester_arm": sub_ssarm,
        "substructure_any_disulfide_n": len(sub_ss),
        "substructure_any_disulfide_examples": sub_ss[:15],
    }
    for name, (full, skel) in ours_keys.items():
        out["exact_full_inchikey_hits"][name] = seen_full.get(full, [])
        out["skeleton_inchikey_hits"][name] = seen_skel.get(skel, [])
    with open(os.path.join(HERE, "02_corpus_overlap.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps({k: v for k, v in out.items()
                      if k != "substructure_any_disulfide_examples"},
                     indent=1)[:4000])


if __name__ == "__main__":
    main()
