#!/usr/bin/env python3
"""MY_2026 (PMID 42569578) -- corpus overlap for the four ionizable lipids.

Structures come from provenance/reagent_registry.py (written from the chemical
NAME, so independent of whatever the corpus contains) for ALC-0315, SM-102 and
DLin-MC3-DMA. cKK-E12 is NOT in the registry, so it is taken from the corpus
through corpus_components.smi(), which RAISES on an ambiguous name.

A corpus join taken FROM the corpus cannot also VERIFY the corpus, so the
cKK-E12 match is reported as NOT INDEPENDENT.

Every query carries a positive control that names the molecule it must match.
Run: python3 -I curation/zhangn_42569578/stage03_corpus_overlap.py
"""
import json
import os
import sys

from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors

RDLogger.DisableLog("rdApp.*")

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "provenance"))
sys.path.insert(0, os.path.join(ROOT, "curation"))
import reagent_registry as REG                      # noqa: E402
from corpus_components import load, smi             # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))


def ikey(smiles):
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        raise SystemExit(f"unparseable SMILES: {smiles}")
    return Chem.MolToInchiKey(m)


def read_corpus():
    import csv
    out = {}
    for tag, fn in (("vitro", "data_vitro.csv"), ("vivo", "data_vivo.csv")):
        with open(os.path.join(ROOT, "data", fn), newline="", encoding="utf8") as f:
            out[tag] = list(csv.DictReader(f))
    return out


def main():
    corpus = read_corpus()
    allrows = corpus["vitro"] + corpus["vivo"]
    comp, ambig = load(allrows)
    rep = {"paper": "PMID 42569578", "experiment_id_proposed": "MY_2026",
           "corpus": {k: len(v) for k, v in corpus.items()},
           "ambiguous_component_names_in_corpus": len(ambig)}

    # ---- IL_SMILES -> InChIKey index over the whole corpus -------------------
    il = {}
    for tag in ("vitro", "vivo"):
        for r in corpus[tag]:
            s = (r.get("IL_SMILES") or "").strip()
            if not s:
                continue
            il.setdefault(s, []).append((tag, r["Experiment_ID"],
                                         (r.get("IL_name") or "").strip()))
    keyidx = {}
    for s in il:
        m = Chem.MolFromSmiles(s)
        if m is None:
            continue
        keyidx.setdefault(Chem.MolToInchiKey(m), []).append(s)
    rep["corpus_distinct_IL_SMILES"] = len(il)
    rep["corpus_distinct_IL_inchikeys"] = len(keyidx)

    # ---- POSITIVE CONTROLS: each asserts a SPECIFIC named molecule ----------
    # A non-zero count is not evidence the control works (F-070 amendment).
    ctrl = []

    def control(label, named, structure_smiles, expect_min):
        """Assert the query matches the molecule `named`, by full InChIKey."""
        k = ikey(structure_smiles)
        hits = keyidx.get(k, [])
        rows = sum(len(il[s]) for s in hits)
        names = sorted({n for s in hits for _, _, n in il[s] if n})
        ok = rows >= expect_min and named.lower() in " | ".join(names).lower()
        ctrl.append({"control": label, "asserts_molecule": named,
                     "inchikey": k, "corpus_rows": rows,
                     "distinct_SMILES_strings": len(hits),
                     "IL_names_carrying_it": names,
                     "expect_min_rows": expect_min,
                     "matched_the_named_molecule": named.lower() in
                     " | ".join(names).lower(), "passes": ok})
        return ok

    # control 1: DLin-MC3-DMA from the registry must be found, under its name.
    control("inchikey_join_registry_MC3", "DLin-MC3-DMA",
            REG.REGISTRY["DLin-MC3-DMA"][0], 1)
    # control 2: a name scan for an excipient known present (DSPC).
    dspc_rows = sum(1 for r in allrows
                    if (r.get("HL_name") or "").strip() == "DSPC")
    ctrl.append({"control": "name_scan_DSPC", "asserts_molecule": "DSPC",
                 "corpus_rows": dspc_rows, "expect_min_rows": 1,
                 "matched_the_named_molecule": dspc_rows > 0,
                 "passes": dspc_rows > 0})
    # control 3: the InChIKey index itself must contain cholesterol under CHL.
    chol_k = ikey(REG.REGISTRY["Cholesterol"][0])
    chl_rows = sum(1 for r in allrows
                   if (r.get("CHL_SMILES") or "").strip() and
                   ikey(r["CHL_SMILES"]) == chol_k) if False else None
    # (full-corpus per-row InChIKey on CHL is slow; do it on distinct strings)
    chl = {}
    for r in allrows:
        s = (r.get("CHL_SMILES") or "").strip()
        if s:
            chl.setdefault(s, 0)
            chl[s] += 1
    chl_hit = sum(n for s, n in chl.items() if ikey(s) == chol_k)
    chl_names = sorted({(r.get("CHL_name") or "").strip() for r in allrows
                        if (r.get("CHL_SMILES") or "").strip() and
                        ikey(r["CHL_SMILES"]) == chol_k}) if chl_hit else []
    ctrl.append({"control": "inchikey_join_registry_cholesterol",
                 "asserts_molecule": "Cholesterol", "inchikey": chol_k,
                 "corpus_rows": chl_hit, "CHL_names_carrying_it": chl_names,
                 "expect_min_rows": 1,
                 "matched_the_named_molecule": any("chol" in n.lower()
                                                   for n in chl_names),
                 "passes": chl_hit > 0 and any("chol" in n.lower()
                                               for n in chl_names)})
    rep["positive_controls"] = ctrl
    if not all(c["passes"] for c in ctrl):
        rep["CONTROLS_FAILED"] = True

    # ---- the paper's four ionizable lipids ----------------------------------
    lipids = []
    for paper_name, src, s in (
            ("ALC-0315", "reagent_registry", REG.REGISTRY["ALC-0315"][0]),
            ("SM-102", "reagent_registry", REG.REGISTRY["SM-102"][0]),
            ("MC-3 (DLin-MC3-DMA)", "reagent_registry",
             REG.REGISTRY["DLin-MC3-DMA"][0]),
            ("CKK-E12 (cKK-E12)", "corpus_components.smi(CKK-E12)",
             smi(comp, ambig, "CKK-E12")),
    ):
        m = Chem.MolFromSmiles(s)
        k = Chem.MolToInchiKey(m)
        hits = keyidx.get(k, [])
        occ = [(tag, eid, nm) for ss in hits for tag, eid, nm in il[ss]]
        by_exp = {}
        for tag, eid, nm in occ:
            by_exp.setdefault(f"{eid}:{tag}", {"rows": 0, "names": set()})
            by_exp[f"{eid}:{tag}"]["rows"] += 1
            if nm:
                by_exp[f"{eid}:{tag}"]["names"].add(nm)
        names = sorted({nm for _, _, nm in occ if nm})
        ent = {
            "paper_name": paper_name,
            "structure_source": src,
            "structure_is_independent_of_corpus":
                src == "reagent_registry",
            "canonical_SMILES": Chem.MolToSmiles(m),
            "inchikey": k, "skeleton_block": k.split("-")[0],
            "formula": Chem.rdMolDescriptors.CalcMolFormula(m),
            "exact_mass": round(Descriptors.ExactMolWt(m), 4),
            "corpus_rows_by_full_inchikey": len(occ),
            "corpus_distinct_SMILES_strings": len(hits),
            "corpus_IL_names": names,
            "corpus_experiments": {k2: {"rows": v["rows"],
                                        "names": sorted(v["names"])}
                                   for k2, v in sorted(by_exp.items())},
            "n_corpus_experiments": len(by_exp),
            "vivo_experiments": sorted({e.split(":")[0] for e in by_exp
                                        if e.endswith(":vivo")}),
        }
        # name join vs inchikey join, for the three-names hazard
        nj = sum(1 for _, _, nm in occ if nm)
        ent["name_join_would_find"] = {
            n: sum(1 for _, _, nm in occ if nm == n) for n in names}
        ent["inchikey_join_finds_rows"] = len(occ)
        ent["rows_a_single_name_join_would_MISS"] = (
            len(occ) - max(ent["name_join_would_find"].values())
            if names else len(occ))
        lipids.append(ent)
    rep["paper_lipids"] = lipids

    # ---- token-pattern test -------------------------------------------------
    toks = {}
    for t in ("E12", "R93", "S2", "S1", "D04", "102", "0315", "MC-3", "MC3"):
        rows = 0
        strucs = set()
        names = set()
        for ss, occ in il.items():
            for tag, eid, nm in occ:
                if nm and t.lower() in nm.lower():
                    rows += 1
                    strucs.add(ss)
                    names.add(nm)
        toks[t] = {"substring_rows": rows, "distinct_structures": len(strucs),
                   "distinct_IL_names": sorted(names)[:25],
                   "n_distinct_IL_names": len(names)}
    rep["token_pattern_test"] = toks

    # ---- names carrying >1 structure (the three-names hazard) --------------
    nm2s = {}
    for ss, occ in il.items():
        for tag, eid, nm in occ:
            if nm:
                nm2s.setdefault(nm, set()).add(ss)
    rep["IL_names_with_more_than_one_structure"] = sum(
        1 for v in nm2s.values() if len(v) > 1)
    s2nm = {}
    for ss, occ in il.items():
        for tag, eid, nm in occ:
            if nm:
                s2nm.setdefault(ss, set()).add(nm)
    rep["IL_structures_with_more_than_one_name"] = sum(
        1 for v in s2nm.values() if len(v) > 1)

    out = os.path.join(HERE, "02_corpus_overlap.json")
    with open(out, "w", encoding="utf8") as f:
        json.dump(rep, f, indent=2, sort_keys=False)
    print(json.dumps(rep, indent=2)[:400])
    print("...\nwrote", out)
    return 0 if all(c["passes"] for c in ctrl) else 1


if __name__ == "__main__":
    sys.exit(main())
