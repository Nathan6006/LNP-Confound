"""Stage 04 -- is any of this paper's chemistry already in the corpus?

Compared on THREE keys, because they fail differently:
  * full InChIKey   -- identical molecule including stereochemistry
  * skeleton block  -- the first InChIKey block: same connectivity, stereo
                       ignored. This is the one that matters here, because
                       these lipids carry four racemic carbinols and nobody's
                       SMILES records them.
  * IL_name         -- a label match, which is evidence about the LABEL and
                       not about the molecule (288 corpus names already carry
                       more than one structure).

REPORT ONLY. Nothing here changes an emitted row.
"""
import os, sys, csv, json, collections

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import reagents                                              # noqa: E402
from rdkit import Chem, RDLogger                             # noqa: E402
RDLogger.DisableLog("rdApp.*")

CORPUS = {
    "vitro": os.path.join(HERE, "..", "..", "data", "data_vitro.csv"),
    "vivo": os.path.join(HERE, "..", "..", "data", "data_vivo.csv"),
}


def keys(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None, None
    k = Chem.MolToInchiKey(m)
    return k, k.split("-")[0]


def main():
    ours = []
    with open(os.path.join(HERE, "03_enumerated.csv")) as f:
        for r in csv.DictReader(f):
            full, skel = keys(r["IL_SMILES"])
            ours.append(dict(name=r["IL_name"], smiles=r["IL_SMILES"],
                             full=full, skel=skel, origin=r["origin"]))

    report = {"_what": "Stage 04 corpus overlap, three keys. Report only.",
              "n_ours": len(ours), "pools": {}}

    for pool, path in CORPUS.items():
        by_full = collections.defaultdict(set)
        by_skel = collections.defaultdict(set)
        by_name = collections.defaultdict(set)
        nrow = 0
        with open(path) as f:
            for r in csv.DictReader(f):
                nrow += 1
                nm = (r.get("IL_name") or "").strip()
                exp = r.get("Experiment_ID", "")
                if nm:
                    by_name[nm].add(exp)
                smi = r.get("IL_SMILES") or ""
                if not smi:
                    continue
                full, skel = keys(smi)
                if full:
                    by_full[full].add(exp)
                    by_skel[skel].add(exp)

        hits = []
        for o in ours:
            h = {"IL_name": o["name"],
                 "full_inchikey_hit": sorted(by_full.get(o["full"], [])),
                 "skeleton_hit": sorted(by_skel.get(o["skel"], [])),
                 "name_hit": sorted(by_name.get(o["name"], []))}
            if any(h[k] for k in ("full_inchikey_hit", "skeleton_hit",
                                  "name_hit")):
                hits.append(h)
        report["pools"][pool] = {
            "file": os.path.relpath(path, HERE), "rows": nrow,
            "unique_full_inchikeys": len(by_full),
            "n_ours_with_any_hit": len(hits), "hits": hits}

    with open(os.path.join(HERE, "04_corpus_overlap.json"), "w") as f:
        json.dump(report, f, indent=1)

    for pool, d in report["pools"].items():
        print(f"\n{pool}: {d['rows']} rows, {d['unique_full_inchikeys']} "
              f"unique lipids; {d['n_ours_with_any_hit']} of {len(ours)} of "
              f"ours hit something")
        for h in d["hits"]:
            print(f"   {h['IL_name']:14s} full={h['full_inchikey_hit']} "
                  f"skel={h['skeleton_hit']} name={h['name_hit']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
