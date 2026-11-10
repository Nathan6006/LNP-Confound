#!/usr/bin/env python3
"""Corpus overlap for the 306 family, three ways plus a near-miss sweep.

WHY THIS IS THE PAPER'S MAIN STRUCTURE CHECK. The SI prints no HRMS, no
molecular formula, no molecular weight and no synthetic scheme beyond one
sentence; the main text is behind a subscription. So there is no in-paper
anchor for any structure. What there IS: two prior, independently built
curations already in data/ that store these exact lipidoid names --
NC_2024 (A306_O10..O14, A306_Oi10, mRNA, BeWo_b30) and KW_2014
(A306_C10..C14, siRNA, HeLa). Comparing against them tests connectivity
against an external ground truth built by a different code path, which is
what C10 exists to do and what mass can never do.

Reported three ways, per CLAUDE.md: FULL InChIKey, SKELETON InChIKey block and
IL_name. The name half matters on its own: CLAUDE.md records 288 IL names that
carry different structures in different experiments, so a name match with a
structure mismatch is the thing to look for, not a reassurance.

The near-miss sweep is here because a sibling's near miss supplied the ground
truth its SI had withheld. Reported, never acted on.
"""
import csv
import json
import os
import sys
from collections import defaultdict

from rdkit import Chem, RDLogger
from rdkit.Chem import inchi

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "lib"))
import reagents                                                  # noqa: E402

RDLogger.DisableLog("rdApp.*")
HERE = os.path.dirname(os.path.abspath(__file__))
CORPUS = os.path.join(HERE, "..", "..", "data")


def load_corpus():
    rows = []
    for pool, fn in (("vitro", "data_vitro.csv"), ("vivo", "data_vivo.csv")):
        with open(os.path.join(CORPUS, fn), newline="") as fh:
            for r in csv.DictReader(fh):
                r["_pool"] = pool
                rows.append(r)
    return rows


def main():
    with open(os.path.join(HERE, "03_enumerated.csv"), newline="") as fh:
        mine = list(csv.DictReader(fh))
    corpus = load_corpus()

    # --- index the corpus by full key, skeleton block and name -----------
    by_full, by_skel, by_name = defaultdict(set), defaultdict(set), defaultdict(set)
    smi_of = {}
    for r in corpus:
        smi = (r.get("IL_SMILES") or "").strip()
        nm = (r.get("IL_name") or "").strip()
        exp = r["Experiment_ID"]
        if nm:
            by_name[nm].add(exp)
        if not smi:
            continue
        try:
            fk = reagents.full_key(smi)
            sk = reagents.skeleton(smi)
        except Exception:
            continue
        by_full[fk].add(exp)
        by_skel[sk].add(exp)
        smi_of.setdefault(fk, smi)

    out = {
        "_what": ("Overlap of this paper's 11 enumerated 306-family lipidoids "
                  "with data/data_vitro.csv + data/data_vivo.csv, by FULL "
                  "InChIKey, SKELETON InChIKey block and IL_name. Reported, "
                  "not acted on."),
        "corpus_rows_scanned": len(corpus),
        "corpus_distinct_full_keys": len(by_full),
        "per_lipid": [],
    }
    n_full = n_skel = n_name = 0
    for m in mine:
        smi = m["IL_SMILES"]
        fk, sk = reagents.full_key(smi), reagents.skeleton(smi)
        name_hits = sorted(by_name.get(m["IL_name"], ()))
        full_hits = sorted(by_full.get(fk, ()))
        skel_hits = sorted(by_skel.get(sk, ()))
        byte_identical = (fk in smi_of
                          and Chem.CanonSmiles(smi) == Chem.CanonSmiles(smi_of[fk]))
        rec = {
            "IL_name": m["IL_name"], "paper_name": m["paper_name"],
            "full_inchikey": fk, "skeleton_block": sk,
            "corpus_full_key_match": full_hits,
            "corpus_skeleton_match": skel_hits,
            "corpus_name_match": name_hits,
            "corpus_smiles_canonically_identical": bool(byte_identical),
            "name_matched_but_structure_did_not":
                bool(name_hits) and not full_hits,
        }
        out["per_lipid"].append(rec)
        n_full += bool(full_hits)
        n_skel += bool(skel_hits)
        n_name += bool(name_hits)

    out["summary"] = {
        "enumerated": len(mine),
        "already_in_corpus_by_full_inchikey": n_full,
        "already_in_corpus_by_skeleton_block": n_skel,
        "already_in_corpus_by_IL_name": n_name,
        "novel_to_the_corpus": len(mine) - n_full,
        "name_match_with_structure_mismatch":
            sum(1 for r in out["per_lipid"]
                if r["name_matched_but_structure_did_not"]),
    }

    # --- near-miss sweep -------------------------------------------------
    # TIER A, exact: the tetra-substituted 306 core itself -- N-methyl, both
    # 3-aminopropyl branches, all four N-CH2CH2-C(=O)O-alkyl arms spelled out.
    # TIER B, loose: any lipid with four acrylate-ester arms on a bridged
    # diamine, i.e. the SAME CHEMISTRY under a different amine core. Tier B is
    # reported as a COUNT PER EXPERIMENT and not as a row list, because my
    # first pattern for it was loose enough to return 332 rows while its label
    # said "same 306 core" -- a scanner whose stated scope is wider than what
    # it measured reports a reassuring number about the wrong set.
    core306 = Chem.MolFromSmarts(
        "[CH3][NX3]([CH2][CH2][CH2][NX3]([CH2][CH2]C(=O)O[CX4])"
        "[CH2][CH2]C(=O)O[CX4])[CH2][CH2][CH2][NX3]"
        "([CH2][CH2]C(=O)O[CX4])[CH2][CH2]C(=O)O[CX4]")
    arm = Chem.MolFromSmarts("[NX3][CH2][CH2]C(=O)O[CX4]")
    ours_keys = {reagents.full_key(x["IL_SMILES"]) for x in mine}
    near, loose = [], defaultdict(int)
    seen = set()
    for r in corpus:
        smi = (r.get("IL_SMILES") or "").strip()
        nm = (r.get("IL_name") or "").strip()
        if not smi or (nm, smi, r["_pool"]) in seen:
            continue
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        seen.add((nm, smi, r["_pool"]))
        if mol.HasSubstructMatch(core306):
            near.append({
                "Experiment_ID": r["Experiment_ID"], "pool": r["_pool"],
                "IL_name": nm, "IL_SMILES": smi,
                "full_inchikey": reagents.full_key(smi),
                "is_one_of_ours": reagents.full_key(smi) in ours_keys,
                "Model_target": r.get("Model_target", ""),
                "Cargo": r.get("Cargo", ""),
                "Experiment_value": r.get("Experiment_value", ""),
            })
        elif len(mol.GetSubstructMatches(arm)) >= 4:
            loose[r["Experiment_ID"]] += 1

    out["near_miss_tierA_exact_306_core"] = {
        "_what": ("Every corpus lipid CONTAINING the tetra-substituted 306 "
                  "core as a SUBSTRUCTURE (N-methyl-3,3'-dipropyl, four "
                  "N-CH2CH2-C(=O)O-alkyl arms). A substructure match is a "
                  "SUPERSET: a larger polyamine lipidoid that embeds this "
                  "motif matches too, which is why the decisive figure here "
                  "is n_identical_to_one_of_ours and not n. This is the "
                  "paper's only "
                  "external structure anchor: KW_2014 stores these same "
                  "structures under 'C<n>' names and NC_2024 under 'O<n>' "
                  "names, so the token's meaning is attested twice, "
                  "independently, instead of assumed from the label."),
        "n": len(near),
        "n_identical_to_one_of_ours": sum(1 for x in near if x["is_one_of_ours"]),
        "rows": sorted(near, key=lambda x: (x["Experiment_ID"], x["IL_name"])),
    }
    out["near_miss_tierB_same_chemistry_other_amine_core"] = {
        "_what": ("Count per experiment of corpus lipids with >=4 acrylate-"
                  "ester Michael arms but NOT the 306 core -- the same "
                  "lipidoid chemistry built on a different amine. Context "
                  "only: it bounds how much of the corpus is the same "
                  "chemical neighbourhood. Reported, not acted on."),
        "per_experiment": dict(sorted(loose.items())),
        "total": sum(loose.values()),
    }

    with open(os.path.join(HERE, "02_corpus_overlap.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps(out["summary"], indent=1))
    print("\nper-lipid:")
    for r in out["per_lipid"]:
        print(f"  {r['IL_name']:12s} full={str(r['corpus_full_key_match']):34s} "
              f"name={str(r['corpus_name_match']):24s} "
              f"canon_identical={r['corpus_smiles_canonically_identical']}")
    print(f"\nTIER A exact 306-core rows in corpus: {len(near)} "
          f"({out['near_miss_tierA_exact_306_core']['n_identical_to_one_of_ours']}"
          " identical to one of ours)")
    for r in out["near_miss_tierA_exact_306_core"]["rows"]:
        print(f"  {r['Experiment_ID']:10s} {r['pool']:6s} {r['IL_name']:12s} "
              f"{r['Model_target']:12s} {r['Cargo']:6s} ours={r['is_one_of_ours']}")
    print("\nTIER B same chemistry, other amine core, per experiment:")
    print("  ", out["near_miss_tierB_same_chemistry_other_amine_core"]["per_experiment"])


if __name__ == "__main__":
    main()
