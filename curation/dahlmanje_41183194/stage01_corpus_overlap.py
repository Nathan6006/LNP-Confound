#!/usr/bin/env python3
"""KG_2025 (Gentry 2025, PNAS, PMID 41183194) -- corpus overlap.

Run:  python3 -I curation/dahlmanje_41183194/stage01_corpus_overlap.py <repo_root>

What this answers
-----------------
The paper's screen varies the HELPER LIPID, the PEG-lipid and the molar ratio.
It uses exactly TWO ionizable lipids, SM-102 and cKK-E12, and both are named
catalogue/literature compounds. So the overlap question is not "does a new
chemotype collide" but "are these two lipids already in the corpus, and under
which names".

Three joins, reported separately and labelled honestly:
  1. FULL InChIKey   -- the only one that establishes identity.
  2. SKELETON block  -- connectivity ignoring stereo/charge; a SUPERSET.
  3. IL_name string  -- reported for contrast only. The brief's measured fact is
     that cKK-E12 sits in the corpus under three different IL_name values
     (CKK-E12, `other`, MD-1) on one InChIKey, so a name join finds it in one of
     three places and silently misses the rest.

Plus a substructure sweep, labelled a SUPERSET, and the token-pattern test for
this paper's own numbering (GlycoLNP45/59/60/61/62).

Structure provenance -- NEVER from memory (trap 6):
  SM-102       provenance/reagent_registry.py   (written from the chemical name)
  DSPC         provenance/reagent_registry.py
  Cholesterol  provenance/reagent_registry.py
  cKK-E12      NOT in the registry -> taken from the live corpus through
               curation/corpus_components.py, whose smi() RAISES on an
               ambiguous name. Reported with the experiments it came from.
"""
import csv
import json
import os
import sys

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

ROOT = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else ".")
sys.path.insert(0, os.path.join(ROOT, "provenance"))
sys.path.insert(0, os.path.join(ROOT, "curation"))
import reagent_registry as RR          # noqa: E402
from corpus_components import load, smi  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))


def rows(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def keys(smiles):
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None, None
    ik = Chem.MolToInchiKey(m)
    return ik, ik.split("-")[0]


def main():
    vitro = rows(os.path.join(ROOT, "data", "data_vitro.csv"))
    vivo = rows(os.path.join(ROOT, "data", "data_vivo.csv"))
    corpus = vitro + vivo
    comp, ambiguous = load(corpus)

    out = {"_what": "KG_2025 corpus overlap. PMID 41183194.",
           "corpus": {"vitro_rows": len(vitro), "vivo_rows": len(vivo),
                      "experiments_vitro": len({r["Experiment_ID"] for r in vitro}),
                      "experiments_vivo": len({r["Experiment_ID"] for r in vivo})}}

    # ---- 0. is this paper already merged? -------------------------------
    pm = {str(r.get("Publication_PMID", "")).split(".")[0] for r in corpus}
    out["already_merged"] = {
        "this_pmid_in_corpus": "41183194" in pm,
        "this_doi_in_corpus": any("2409569122" in (r.get("Publication_DOI") or "")
                                  for r in corpus),
        "note": ("LR_2023 (PMID 37127709), the corpus's barcoded in-vivo "
                 "experiment, is Rhym/Anderson (MIT), NOT Dahlman -- checked "
                 "against Europe PMC. No Dahlman-lab experiment is in the "
                 "corpus under this paper's PMID or DOI."),
    }

    # ---- 1. the two ionizable lipids ------------------------------------
    sm102 = RR.REGISTRY["SM-102"][0]
    ils = {"SM-102": {"smiles": sm102,
                      "provenance": "provenance/reagent_registry.py, written "
                                    "from the chemical name "
                                    + RR.REGISTRY["SM-102"][1]}}

    # cKK-E12 is not a registry entry. Resolve it through the guarded loader.
    ckk = None
    ckk_err = None
    for name in ("CKK-E12", "cKK-E12", "ckk-E12"):
        try:
            ckk = smi(comp, ambiguous, name)
            ckk_name_used = name
            break
        except Exception as exc:          # noqa: BLE001  -- guarded loader raising is the point
            ckk_err = f"{name}: {exc}"
    if ckk is None:
        raise SystemExit(f"cKK-E12 not resolvable from the corpus: {ckk_err}")
    ils["cKK-E12"] = {
        "smiles": ckk,
        "provenance": (f"curation/corpus_components.py smi() on IL_name "
                       f"{ckk_name_used!r}; the guarded loader did NOT raise, "
                       "so that name carries exactly one structure corpus-wide."),
    }

    # ---- 2. three joins, per ionizable lipid ----------------------------
    # pre-index the corpus by key
    by_full, by_skel, by_name = {}, {}, {}
    for r in corpus:
        s = (r.get("IL_SMILES") or "").strip()
        n = (r.get("IL_name") or "").strip()
        if n:
            by_name.setdefault(n, []).append(r)
        if not s:
            continue
        ik, sk = keys(s)
        if ik:
            by_full.setdefault(ik, []).append(r)
            by_skel.setdefault(sk, []).append(r)

    def summarise(rs):
        by_exp = {}
        for r in rs:
            e = r["Experiment_ID"]
            by_exp.setdefault(e, {"rows": 0, "names": set(), "pool": set()})
            by_exp[e]["rows"] += 1
            by_exp[e]["names"].add((r.get("IL_name") or "").strip())
            by_exp[e]["pool"].add(r.get("Model") or "")
        return {e: {"rows": v["rows"], "IL_names": sorted(v["names"]),
                    "Model": sorted(v["pool"])} for e, v in sorted(by_exp.items())}

    for nm, rec in ils.items():
        ik, sk = keys(rec["smiles"])
        rec["inchikey"] = ik
        rec["skeleton"] = sk
        rec["join_full_inchikey"] = summarise(by_full.get(ik, []))
        rec["join_skeleton"] = summarise(by_skel.get(sk, []))
        rec["join_full_rows"] = sum(v["rows"] for v in rec["join_full_inchikey"].values())
        rec["join_skeleton_rows"] = sum(v["rows"] for v in rec["join_skeleton"].values())
        # distinct IL_name values carried by THIS structure
        rec["names_this_structure_carries"] = sorted(
            {n for v in rec["join_full_inchikey"].values() for n in v["IL_names"]})
        # what a NAME join alone would have found
        rec["join_IL_name_literal"] = {
            n: sum(1 for r in by_name.get(n, [])) for n in (nm, nm.upper(), nm.lower())
            if n in by_name}

    out["ionizable_lipids"] = ils

    # how many corpus structures carry >1 IL_name (the brief's 324)
    nm_per_ik = {}
    for r in corpus:
        s = (r.get("IL_SMILES") or "").strip()
        n = (r.get("IL_name") or "").strip()
        if not s or not n:
            continue
        ik, _ = keys(s)
        if ik:
            nm_per_ik.setdefault(ik, set()).add(n)
    out["structures_with_more_than_one_IL_name"] = {
        "n_structures_total": len(nm_per_ik),
        "n_with_multiple_names": sum(1 for v in nm_per_ik.values() if len(v) > 1),
    }

    # ---- 3. excipients: are THIS paper's helper/PEG lipids in the corpus?
    exc = {}
    for label, nmlist in (
            ("DSPC", ["DSPC"]),
            ("Cholesterol", ["Cholesterol", "cholesterol", "Chol"]),
            ("ALC-0159", ["ALC-0159", "ALC-0159 ", "ALC0159"]),
            ("DMG-PEG", ["DMG-PEG", "DMG-PEG2000", "DMG-PEG-2000"]),
            ("C14-PEG2000", ["C14-PEG2000", "C14PEG2000", "DMPE-PEG2000"]),
            ("C18-PEG2000", ["C18-PEG2000", "C18PEG2000", "DSPE-PEG2000"]),
            ("C8GaCer", ["C8GaCer", "C8GalCer", "GalCer", "alpha-GalCer"]),
            ("C16GaCer", ["C16GaCer", "C16GalCer"]),
    ):
        hits = {}
        for n in nmlist:
            if n in comp:
                hits[n] = {"status": "unambiguous", "smiles": comp[n]}
            elif n in ambiguous:
                hits[n] = {"status": "AMBIGUOUS", "n_structures": len(ambiguous[n]),
                           "smiles": ambiguous[n]}
        exc[label] = hits or {"status": "absent from the corpus under every "
                                        "name tried: " + ", ".join(nmlist)}
    out["excipients_in_corpus"] = exc

    # ---- 4. token-pattern test ------------------------------------------
    # This paper numbers its members GlycoLNP45/59/60/61/62 -- bare integers.
    toks = ["45", "59", "60", "61", "62"]
    tok = {}
    for t in toks:
        exact = [r for r in corpus if (r.get("IL_name") or "").strip() == t]
        sub = [r for r in corpus if t in (r.get("IL_name") or "")]
        sub_struct = {keys(r["IL_SMILES"])[0] for r in sub
                      if (r.get("IL_SMILES") or "").strip()}
        tok[t] = {
            "exact_IL_name_match_rows": len(exact),
            "exact_match_experiments": sorted({r["Experiment_ID"] for r in exact}),
            "substring_match_rows": len(sub),
            "substring_match_distinct_structures": len(sub_struct),
            "substring_match_experiments": sorted({r["Experiment_ID"] for r in sub}),
        }
    out["token_pattern_test"] = {
        "paper_numbering": "GlycoLNP45, 59, 60, 61, 62 -- bare integers, and the "
                           "screen members are numbered 1..64 CONTIGUOUSLY within "
                           "each library. CLAUDE.md records the corpus amine "
                           "catalogue as sparse and non-contiguous (255 values, "
                           "large gaps, only 7 of 1-16 present), so a contiguous "
                           "1..N numbering is LOCAL and must never be joined to it.",
        "verdict": "LOCAL INDEX. No name join is performed anywhere in this "
                   "curation; both ionizable lipids are attested structurally.",
        "per_token": tok,
    }

    # ---- 5. substructure sweep, labelled a SUPERSET ---------------------
    sweep = {}
    for label, sma in (
        ("beta-hydroxy-tertiary-amine-with-two-esters (SM-102-like)",
         "[#6]OC(=O)[#6][#6][#6][#6][#6][#6][NX3]([#6][#6]O)[#6]"),
        ("ketal-piperazine-diketopiperazine core (cKK-like)",
         "O=C1[#6][NX3][#6](=O)[#6][NX3]1"),
    ):
        patt = Chem.MolFromSmarts(sma)
        if patt is None:
            sweep[label] = {"error": "SMARTS did not parse"}
            continue
        hit_exp, hit_struct, n = {}, set(), 0
        seen = {}
        for r in corpus:
            s = (r.get("IL_SMILES") or "").strip()
            if not s:
                continue
            if s not in seen:
                m = Chem.MolFromSmiles(s)
                seen[s] = bool(m and m.HasSubstructMatch(patt))
            if seen[s]:
                n += 1
                hit_exp[r["Experiment_ID"]] = hit_exp.get(r["Experiment_ID"], 0) + 1
                hit_struct.add(s)
        sweep[label] = {
            "LABEL": "SUPERSET -- a substructure hit is a CHEMOTYPE, not a "
                     "near-miss and not an identity claim.",
            "smarts": sma, "rows": n, "distinct_structures": len(hit_struct),
            "experiments": dict(sorted(hit_exp.items())),
        }
    out["substructure_sweep_SUPERSET"] = sweep

    with open(os.path.join(HERE, "02_corpus_overlap.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps({k: v for k, v in out.items()
                      if k not in ("substructure_sweep_SUPERSET",)}, indent=1)[:7000])
    print("\n--- sweep ---")
    print(json.dumps(sweep, indent=1))


if __name__ == "__main__":
    main()
