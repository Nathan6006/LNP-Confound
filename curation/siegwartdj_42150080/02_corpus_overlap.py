#!/usr/bin/env python3
"""Stage 02 - corpus overlap. Run: python3 -I 02_corpus_overlap.py

Brief's trap 2: "Members may already be in the corpus ... After enumerating,
compare full InChIKeys and IL_names against data/data_vitro.csv and
data/data_vivo.csv; report overlap and the Experiment_ID. Report, do not act."

This paper's SI exposes exactly one ionizable library member by name, 6A1-SC8,
so the comparison is small but decisive: it is a published Siegwart dendrimer
lipid and the corpus already carries it.

Comparison is by FULL InChIKey, not by name and not by skeleton block. The
project has been bitten twice by a skeleton-only comparison (F-013: the
cholesterol/epicholesterol epimers share the skeleton block RWTQCZGAMKTBRV) and
once by keying on IL_name across experiments (288 corpus names carry more than
one structure).

READ-ONLY. Writes 02_corpus_overlap.json inside this folder and nothing else.
"""
import json
import os

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "..", "data")

# The only ionizable library member this paper's SI names and characterises.
# Fig S1 gives its 1H-NMR and LC/MS; the SI draws no structure for it.
PAPER_IL_NAMES = ["6A1-SC8"]

# The three helper lipids of the titration series, by their Avanti shorthand as
# printed on the Fig S2/S3 axes. Looked up in the corpus's excipient columns.
HELPERS = {"18:1 TAP": "DOTAP", "18:1 DAP": "DODAP", "18:1 PA": "DOPA"}


def key(smiles):
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None, None, None
    return Chem.MolToInchiKey(m), CalcMolFormula(m), Chem.MolToSmiles(m)


def main():
    out = {
        "_what": "Corpus overlap for Sun 2026 PNAS (PMID 42150080). "
                 "Report only - nothing in data/ was read for writing and "
                 "nothing was merged.",
        "corpus": {},
        "il_name_overlap": [],
        "helper_lipid_overlap": {},
        "verdict": None,
    }

    frames = {}
    for pool in ("vitro", "vivo"):
        p = os.path.join(DATA, f"data_{pool}.csv")
        d = pd.read_csv(p, low_memory=False)
        frames[pool] = d
        out["corpus"][pool] = {"path": os.path.relpath(p, HERE),
                               "rows": int(len(d)),
                               "experiments": int(d["Experiment_ID"].nunique())}

    for name in PAPER_IL_NAMES:
        rec = {"paper_IL_name": name, "hits": []}
        for pool, d in frames.items():
            m = d["IL_name"].astype(str).str.strip() == name
            for _, r in d[m].iterrows():
                ik, formula, canon = key(r["IL_SMILES"])
                rec["hits"].append({
                    "pool": pool,
                    "Experiment_ID": r["Experiment_ID"],
                    "Publication_PMID": (None if pd.isna(r["Publication_PMID"])
                                         else str(r["Publication_PMID"])),
                    "Model_target": r["Model_target"],
                    "Model_type": r["Model_type"],
                    "Route_of_administration": (
                        None if pd.isna(r["Route_of_administration"])
                        else r["Route_of_administration"]),
                    "Cargo": None if pd.isna(r["Cargo"]) else r["Cargo"],
                    "Experiment_value": float(r["Experiment_value"]),
                    "inchikey_full": ik,
                    "inchikey_skeleton": ik.split("-")[0] if ik else None,
                    "formula": formula,
                    "canonical_SMILES": canon,
                })
        # every structure stored under this name, across experiments
        iks = {h["inchikey_full"] for h in rec["hits"]}
        rec["n_hits"] = len(rec["hits"])
        rec["n_distinct_full_inchikeys_under_this_name"] = len(iks)
        rec["name_is_unambiguous_in_corpus"] = (len(iks) <= 1)
        out["il_name_overlap"].append(rec)

    # siblings of the same dendrimer series, for context
    sib = {}
    for pool, d in frames.items():
        n = d["IL_name"].astype(str).str.strip()
        m = n.str.match(r"^6A1-S")
        if m.any():
            sib[pool] = sorted(set(n[m]))
    out["same_series_names_in_corpus"] = sib

    for short, common in HELPERS.items():
        hit = {}
        for pool, d in frames.items():
            for col in ("HL_name", "CHL_name", "fifthcomponent_name", "IL_name"):
                s = d[col].astype(str).str.strip()
                n = int((s.str.upper() == common).sum()
                        + (s.str.upper() == short.upper()).sum())
                if n:
                    hit.setdefault(pool, {})[col] = n
        out["helper_lipid_overlap"][f"{short} ({common})"] = hit or "absent from the corpus"

    overlapped = [r for r in out["il_name_overlap"] if r["n_hits"]]
    out["verdict"] = (
        f"{len(overlapped)}/{len(PAPER_IL_NAMES)} named ionizable library "
        "member(s) of this paper are ALREADY in the corpus. "
        + "; ".join(f"{r['paper_IL_name']} -> "
                    + ",".join(sorted({h['Experiment_ID'] for h in r['hits']}))
                    for r in overlapped)
        + ". Reported, not acted on: no row was added, changed or removed.")

    path = os.path.join(HERE, "02_corpus_overlap.json")
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps(out, indent=1))
    print("wrote", path)


if __name__ == "__main__":
    main()
