#!/usr/bin/env python3
"""Trap 2 of the brief: are these members already in the corpus?

Compared three ways, because each misses something the others catch:
  FULL InChIKey      stereo included -- the strictest identity.
  SKELETON block     connectivity only; catches a corpus row that lost a Z.
  IL_name            catches a re-screened library whose names were reused
                     (the KZ_2016 case, 465 of 465 names matched).
REPORTS. Does not act. CLAUDE.md: "Report it; do not act on it."
"""
import csv, json, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "lib"))
import reagents                                              # noqa: E402
HERE = os.path.dirname(os.path.abspath(__file__))
CORPUS = os.path.join(HERE, "..", "..", "data")


def main():
    enum = list(csv.DictReader(open(os.path.join(HERE, "03_enumerated.csv"))))
    ours_full = {r["InChIKey"]: r["IL_name"] for r in enum}
    ours_skel = {}
    for r in enum:
        ours_skel.setdefault(r["InChIKey"].split("-")[0], []).append(r["IL_name"])
    ours_names = {r["IL_name"] for r in enum}

    rep = {"n_enumerated": len(enum),
           "n_distinct_full_InChIKey": len(ours_full),
           "pools": {}}
    for pool, fn in (("vitro", "data_vitro.csv"), ("vivo", "data_vivo.csv")):
        hits_full, hits_skel, hits_name = {}, {}, {}
        n_rows = 0
        seen = {}
        for row in csv.DictReader(open(os.path.join(CORPUS, fn))):
            n_rows += 1
            smi = row.get("IL_SMILES") or ""
            key = seen.get(smi)
            if key is None:
                key = reagents.full_key(smi) or ""
                seen[smi] = key
            if key and key in ours_full:
                hits_full.setdefault(ours_full[key], set()).add(row["Experiment_ID"])
            if key and key.split("-")[0] in ours_skel:
                for nm in ours_skel[key.split("-")[0]]:
                    hits_skel.setdefault(nm, set()).add(row["Experiment_ID"])
            nm = (row.get("IL_name") or "").strip()
            if nm in ours_names:
                hits_name.setdefault(nm, set()).add(row["Experiment_ID"])
        rep["pools"][pool] = {
            "corpus_rows_scanned": n_rows,
            "n_our_structures_present_by_FULL_InChIKey": len(hits_full),
            "by_full_InChIKey": {k: sorted(v) for k, v in hits_full.items()},
            "n_our_structures_present_by_SKELETON": len(hits_skel),
            "by_skeleton": {k: sorted(v) for k, v in hits_skel.items()},
            "n_our_names_present": len(hits_name),
            "by_name": {k: sorted(v) for k, v in hits_name.items()}}
    rep["same_lab_experiments_already_in_corpus"] = {
        "BL_2023": "PMID 36997680 (Li lab, A*-T* amino-lipid library)",
        "BL_2023_3CR": "PMID 37679571 (Li lab, Ugi 3CR, A*_I*_T* names)",
        "BL_2024": "PMID 38740955 (Li lab, numeric names)",
        "note": "Three Bowen Li experiments are already in the corpus. None is "
                "PMID 39883839. The Ugi 3CR paper is this paper's own ref 1 "
                "and is the closest relative: a different multicomponent "
                "reaction giving an alpha-AMINO amide, not an alpha-ACYLOXY "
                "amide."}
    json.dump(rep, open(os.path.join(HERE, "04_corpus_overlap.json"), "w"),
              indent=1)
    for p, v in rep["pools"].items():
        print(p, "rows", v["corpus_rows_scanned"],
              "| full-key hits", v["n_our_structures_present_by_FULL_InChIKey"],
              "| skeleton hits", v["n_our_structures_present_by_SKELETON"],
              "| name hits", v["n_our_names_present"])
        if v["by_skeleton"]:
            print("   ", v["by_skeleton"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
