#!/usr/bin/env python3
"""Stage 02 -- overlap against the LIVE corpus, including already-merged arms.

Four separate questions, answered separately because they have different
failure modes:

 1  EXACT      full InChIKey of each enumerated ionizable lipid vs every
                distinct IL_SMILES in data_vitro.csv + data_vivo.csv.
 2  SKELETON   InChIKey first block (connectivity, stereo-blind). A skeleton
                hit with no full hit means a stereoisomer is already present.
 3  NAME       IL_name match, case-folded. Reported SEPARATELY from structure,
                because a name match is not a structural attestation.
 4  NEAR-MISS  a deliberately LOOSE SMARTS sweep, reported as a SUPERSET. It is
                not evidence of anything on its own; its purpose is to make the
                zero in question 1 checkable rather than asserted.

Plus two things the brief insists on:

 *  the ALREADY-MERGED ARM test. AC_2025 (the COMET dataset) contains 216
    cKK-E12 rows that vary the SAME four formulation parameters this paper
    varies. If any of our seven formulations were also measured there, and the
    readouts joined, re-emitting them under a new Experiment_ID would
    manufacture exactly the source-separability signal this project measures.
    So: match formulation-by-formulation on (IL, helper, 4 mole percents).

 *  the TOKEN-PATTERN test. 'LNP1'...'LNP150' are contiguous local indices, not
    members of the shared lipidoid amine catalogue. A contiguous 1..N numbering
    is the AC_2025 failure mode, so no corpus name lookup is performed on them.
"""
import csv
import json
import os
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))

from rdkit import Chem, RDLogger   # noqa: E402

RDLogger.DisableLog("rdApp.*")

# Loose on purpose: a 2-hydroxyalkyl on a tertiary amine, the cKK-E12 motif.
NEAR_MISS_SMARTS = "[NX3;!$(N=*)](CC(O)[CH2])"


def corpus():
    rows = []
    for f in ("data_vitro.csv", "data_vivo.csv"):
        p = os.path.join(ROOT, "data", f)
        with open(p) as fh:
            for r in csv.DictReader(fh):
                r["_file"] = f
                rows.append(r)
    return rows


def main():
    rows = corpus()
    with open(os.path.join(HERE, "03_enumerated.csv")) as fh:
        mine = list(csv.DictReader(fh))

    # index the corpus by structure
    by_full, by_skel, by_name = defaultdict(set), defaultdict(set), defaultdict(set)
    names_of_structure = defaultdict(lambda: defaultdict(int))
    smiles_seen = {}
    for r in rows:
        s = (r.get("IL_SMILES") or "").strip()
        n = (r.get("IL_name") or "").strip()
        e = r["Experiment_ID"]
        if n:
            by_name[n.casefold()].add(e)
        if not s:
            continue
        if s not in smiles_seen:
            m = Chem.MolFromSmiles(s)
            smiles_seen[s] = Chem.MolToInchiKey(m) if m else None
        k = smiles_seen[s]
        if k:
            by_full[k].add(e)
            by_skel[k.split("-")[0]].add(e)
            names_of_structure[k][(e, n)] += 1

    res = {
        "_what": "overlap of this paper's individually-administered ionizable "
                 "lipids against the live corpus (29,195 vitro + 1,803 vivo "
                 "rows as held on 2026-11-10)",
        "corpus_rows": len(rows),
        "corpus_distinct_IL_structures": len(by_full),
        "per_lipid": [],
    }
    for rec in mine:
        k, sk, nm = rec["IL_inchikey"], rec["IL_skeleton"], rec["IL_name"]
        res["per_lipid"].append({
            "formulation": rec["formulation_id"],
            "IL_name": nm,
            "inchikey": k,
            "exact_full_inchikey_hit": sorted(by_full.get(k, [])),
            "skeleton_only_hit": sorted(by_skel.get(sk, []) - by_full.get(k, set())),
            "IL_name_hit": sorted(by_name.get(nm.casefold(), [])),
            "names_the_corpus_gives_this_structure": sorted(
                [{"experiment": e, "IL_name": n, "rows": c}
                 for (e, n), c in names_of_structure.get(k, {}).items()],
                key=lambda d: (d["experiment"], d["IL_name"])),
        })

    # --- near-miss superset ------------------------------------------------
    patt = Chem.MolFromSmarts(NEAR_MISS_SMARTS)
    hits = defaultdict(set)
    for s, k in smiles_seen.items():
        m = Chem.MolFromSmiles(s)
        if m is not None and m.HasSubstructMatch(patt):
            for r in rows:
                if (r.get("IL_SMILES") or "").strip() == s:
                    hits[k].add(r["Experiment_ID"])
    exact = {p["inchikey"] for p in res["per_lipid"]}
    res["near_miss_superset"] = {
        "smarts": NEAR_MISS_SMARTS,
        "what_it_is": "A SUPERSET, deliberately loose. Any corpus lipid "
                      "bearing a 2-hydroxyalkyl arm on a non-amide nitrogen "
                      "matches, which is a whole chemotype, not a near-miss "
                      "of cKK-E12. It is reported so the exact-match answer "
                      "can be seen not to rest on a silent zero.",
        "distinct_structures_matching": len(hits),
        "experiments_touched": sorted({e for v in hits.values() for e in v}),
        "of_these_are_exact_matches_to_ours": len(exact & set(hits)),
    }

    # --- already-merged arm: AC_2025's cKK-E12 formulations ----------------
    ck = [p for p in res["per_lipid"] if p["IL_name"] == "cKK-E12"][0]
    arm_exps = ck["exact_full_inchikey_hit"]
    arm_rows = [r for r in rows if r["Experiment_ID"] in arm_exps
                and (r.get("IL_name") or "").strip().casefold() == "ckk-e12"]
    corpus_forms = set()
    for r in arm_rows:
        corpus_forms.add((
            (r.get("HL_name") or "").strip(),
            round(float(r["IL_molratio"]), 3) if r.get("IL_molratio") else None,
            round(float(r["CHL_molratio"]), 3) if r.get("CHL_molratio") else None,
            round(float(r["PEG_molratio"]), 3) if r.get("PEG_molratio") else None,
            round(float(r["HL_molratio"]), 3) if r.get("HL_molratio") else None))
    mine_forms = {}
    for rec in mine:
        if rec["IL_name"] != "cKK-E12":
            continue
        mine_forms[rec["formulation_id"]] = (
            rec["HL_name"], round(float(rec["IL_molratio"]), 3),
            round(float(rec["CHL_molratio"]), 3),
            round(float(rec["PEG_molratio"]), 3),
            round(float(rec["HL_molratio"]), 3))
    collide = {k: v for k, v in mine_forms.items() if v in corpus_forms}
    res["already_merged_arm_test"] = {
        "why": "AC_2025 (COMET) holds cKK-E12 rows that vary the SAME four "
               "formulation parameters. Re-emitting a formulation the corpus "
               "already measures, under a new Experiment_ID, would manufacture "
               "source separability.",
        "corpus_experiments_holding_cKK-E12": arm_exps,
        "corpus_cKK-E12_rows": len(arm_rows),
        "corpus_distinct_cKK-E12_formulations": len(corpus_forms),
        "our_cKK-E12_formulations": len(mine_forms),
        "formulation_collisions": sorted(collide),
        "verdict": ("NO formulation collision" if not collide
                    else "COLLISION -- do not emit without resolving"),
        "endpoint_is_different_anyway": {
            "corpus_arm": sorted({(r["Experiment_ID"], r["Model"],
                                   r["Model_type"], r["Experiment_method"])
                                  for r in arm_rows}),
            "ours": "in_vivo / Mouse_B6 / fraction of lymph-node LECs "
                    "expressing the delivered aVHH protein",
            "note": "Different organism-level model and a different readout, "
                    "so no value could be a re-publication of a corpus value. "
                    "The formulation test is the stronger one and it is clean.",
        },
    }

    # --- token-pattern test -------------------------------------------------
    labels = [rec["formulation_id"] for rec in mine if rec["formulation_id"] != "MC3"]
    allnums = sorted(int(x[3:]) for x in
                     [r[0] for r in csv.reader(open(os.path.join(HERE, "01_screening1.csv")))][1:])
    res["token_pattern_test"] = {
        "tokens": labels,
        "library_numbering": f"LNP1..LNP{max(allnums)} contiguous over "
                             f"Screening 1 ({len(allnums)} members, "
                             f"min {min(allnums)}, max {max(allnums)}, "
                             f"gaps {sorted(set(range(min(allnums), max(allnums) + 1)) - set(allnums))})",
        "verdict": "CONTIGUOUS LOCAL INDICES. A contiguous 1..N numbering is "
                   "the AC_2025 failure mode: the shared lipidoid amine "
                   "catalogue is sparse and non-contiguous, so a corpus name "
                   "lookup on 'LNP7' would return a plausible wrong molecule. "
                   "NO corpus name lookup is performed on these tokens. They "
                   "are carried in Formulation_ID only, never in IL_name.",
        "structural_attestation_instead": "Both ionizable lipids are named "
                                          "catalogue reagents with a supplier "
                                          "and catalogue number printed in the "
                                          "Methods (BroadPharm BP-29590; MedKoo "
                                          "555308), and their structures come "
                                          "from the registry / the guarded "
                                          "corpus lookup, not from a number.",
    }

    res["one_structure_many_names"] = {
        "_what": "C17 pins one NAME to one structure. The reverse is not "
                 "checked, and it matters here: our cKK-E12 is already in the "
                 "corpus under names that do not look like each other.",
        "detail": [{"IL_name": p["IL_name"],
                    "names": p["names_the_corpus_gives_this_structure"]}
                   for p in res["per_lipid"] if p["IL_name"] == "cKK-E12"],
        "consequence": "Keying on IL_name across experiments would miss this "
                       "lipid in two of the three experiments that already "
                       "hold it. Our rows must carry the SMILES the corpus "
                       "already uses, not a fresh one.",
    }

    with open(os.path.join(HERE, "02_corpus_overlap.json"), "w") as fh:
        json.dump(res, fh, indent=1)

    for p in res["per_lipid"]:
        print(f"{p['formulation']:<6} {p['IL_name']:<14} full={p['exact_full_inchikey_hit']} "
              f"skel_only={p['skeleton_only_hit']} name={p['IL_name_hit']}")
    print()
    print("near-miss superset:", res["near_miss_superset"]["distinct_structures_matching"],
          "structures in", res["near_miss_superset"]["experiments_touched"])
    print("already-merged arm:", res["already_merged_arm_test"]["verdict"],
          f"({res['already_merged_arm_test']['corpus_distinct_cKK-E12_formulations']} "
          f"corpus formulations vs our {res['already_merged_arm_test']['our_cKK-E12_formulations']})")
    print("token pattern:", res["token_pattern_test"]["library_numbering"])
    print("one structure, many names:",
          res["one_structure_many_names"]["detail"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
