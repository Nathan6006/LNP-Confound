#!/usr/bin/env python3
"""Is any of these 14 lipids already in data/? Three independent keys.

RUN:  python3 -I 04_corpus_overlap.py

WHY THREE KEYS. CLAUDE.md section 3: 288 IL names carry DIFFERENT structures in
different experiments, so a name match is not a structure match; and the F-013
epicholesterol class means a SKELETON match is not a stereochemistry match.
So compare on all three and report where they disagree:

  full InChIKey      identity including stereochemistry
  skeleton block     connectivity only -- a hit here with a miss on the full
                     key is exactly the F-013 signature (same skeleton, same
                     formula, same exact mass, different molecule)
  IL_name            "Lipid 1" .. "Lipid 14" are generic combinatorial labels,
                     so a hit here is expected to be a COLLISION, not a match

READ-ONLY. This script opens data/ and writes only into this folder.
"""
import csv
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "lib"))
import reagents                                            # noqa: E402
from rdkit import Chem, RDLogger                           # noqa: E402

RDLogger.DisableLog("rdApp.*")
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))


def load_corpus():
    rows = []
    for pool, fn in (("vitro", "data_vitro.csv"), ("vivo", "data_vivo.csv")):
        with open(os.path.join(ROOT, "data", fn), newline="") as f:
            for r in csv.DictReader(f):
                r["_pool"] = pool
                rows.append(r)
    return rows


def main():
    with open(os.path.join(HERE, "03_enumerated.csv"), newline="") as f:
        ours = list(csv.DictReader(f))
    corpus = load_corpus()

    # Index the corpus three ways. reagents.full_key / .skeleton canonicalise
    # through RDKit, so a corpus SMILES written differently still matches.
    by_full, by_skel, by_name = {}, {}, {}
    bad = 0
    for r in corpus:
        smi = (r.get("IL_SMILES") or "").strip()
        nm = (r.get("IL_name") or "").strip()
        if nm:
            by_name.setdefault(nm, set()).add(r["Experiment_ID"])
        if not smi:
            continue
        try:
            fk = reagents.full_key(smi)
            sk = reagents.skeleton(smi)
        except Exception:
            bad += 1
            continue
        by_full.setdefault(fk, set()).add(r["Experiment_ID"])
        by_skel.setdefault(sk, set()).add(r["Experiment_ID"])

    out, n_full, n_skel, n_name = [], 0, 0, 0
    for o in ours:
        mol = Chem.MolFromSmiles(o["IL_SMILES"])
        fk = Chem.MolToInchiKey(mol)
        sk = fk.split("-")[0]
        f_hit = sorted(by_full.get(fk, ()))
        s_hit = sorted(by_skel.get(sk, ()))
        n_hit = sorted(by_name.get(o["IL_name"], ()))
        n_full += bool(f_hit)
        n_skel += bool(s_hit)
        n_name += bool(n_hit)
        out.append({
            "lipid_id": o["lipid_id"], "IL_name": o["IL_name"],
            "inchikey": fk, "skeleton": sk,
            "full_key_hit_experiments": ";".join(f_hit),
            "skeleton_hit_experiments": ";".join(s_hit),
            "IL_name_hit_experiments": ";".join(n_hit),
            "skeleton_hit_without_full_hit": int(bool(s_hit) and not f_hit),
        })

    with open(os.path.join(HERE, "04_corpus_overlap.csv"), "w",
              newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0].keys()))
        w.writeheader()
        w.writerows(out)

    summary = {
        "corpus_rows_scanned": len(corpus),
        "corpus_unique_full_inchikeys": len(by_full),
        "corpus_unique_skeletons": len(by_skel),
        "corpus_unparseable_IL_SMILES_skipped": bad,
        "our_lipids": len(ours),
        "matches_on_full_inchikey": n_full,
        "matches_on_skeleton_block": n_skel,
        "matches_on_IL_name": n_name,
        "free_external_check_available": n_full > 0,
        "interpretation": (
            "0 / 0 / 0 means this paper contributes 14 structures none of "
            "which the corpus has ever seen, and -- the part that costs "
            "something -- that NO free external check is available from "
            "overlap. A sibling got 21/21 independent confirmation of its "
            "head/tail assignment because corpus lipidoid names encode chain "
            "length; there is no such gift here. The ESI-MS anchors and the "
            "1H NMR integrals are the only external checks this paper has, "
            "and 07_qc.json leans on them accordingly."),
    }
    with open(os.path.join(HERE, "04_corpus_overlap.json"), "w") as f:
        json.dump(summary, f, indent=1)
    print(json.dumps(summary, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
