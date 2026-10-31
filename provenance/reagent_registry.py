#!/usr/bin/env python3
"""Tier-1 check C11: named commercial reagents must mean one fixed molecule.

WHY. Two 2026-10-29 bugs shared a shape that no existing check could see:

  * `AC_2025` held the DLin-KC2-DMA dioxolane under the name `DLin-MC3-DMA` and
    the MC3 ester under `KC2`. Both are C43H79NO2 at 641.6111 Da, so mass,
    formula and descriptor checks were blind to it.
  * `18PG` was stored as a symmetric 18:0/18:0 lipid across 360 rows. Avanti's
    18PG is 18:0/18:1 (SOPG), as YZ_2022 spells out in its own Materials.

C10 cannot cover these: it checks lipids a curation folder ENUMERATED against a
structure a paper DREW. Catalogue reagents arrive already named, from no figure.

WHAT THIS DOES. Each entry below pins a reagent name to a structure written from
its chemical name, with the authority for that name recorded next to it. The
check then asserts, over the whole corpus, that:

  1. every occurrence of a registered name carries the registered molecule
     (compared by InChIKey SKELETON, so stereo descriptors do not cause noise), and
  2. no registered name carries two different structures anywhere.

An entry is only worth having if its structure is genuinely fixed and public.
Polymer-bearing reagents (DMG-PEG2000 and friends) are deliberately absent: the
PEG block length is a distribution, so there is no single right SMILES and a
registry entry would assert something false.

Disagreement is a FINDING, not proof the corpus is wrong -- when this was first
run, one disagreement was the registry's own typo and the rest were real. Look
at every one.
"""
import csv
import json
import os
import sys
from collections import defaultdict

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
OUT = os.path.join(HERE, "reagent_registry.json")

# name -> (SMILES written from the chemical name, authority for the name)
REGISTRY = {
    "DSPC": (
        "CCCCCCCCCCCCCCCCCC(=O)OC[C@H](COP(=O)([O-])OCC[N+](C)(C)C)OC(=O)CCCCCCCCCCCCCCCCC",
        "1,2-distearoyl-sn-glycero-3-phosphocholine; 18:0/18:0"),
    "DOPC": (
        "CCCCCCCC/C=C\\CCCCCCCC(=O)OC[C@H](COP(=O)([O-])OCC[N+](C)(C)C)OC(=O)CCCCCCC/C=C\\CCCCCCCC",
        "1,2-dioleoyl-sn-glycero-3-phosphocholine; 18:1/18:1 cis-9"),
    "DOPE": (
        "CCCCCCCC/C=C\\CCCCCCCC(=O)OC[C@H](COP(=O)([O-])OCC[NH3+])OC(=O)CCCCCCC/C=C\\CCCCCCCC",
        "1,2-dioleoyl-sn-glycero-3-phosphoethanolamine; 18:1/18:1 cis-9"),
    "DOTAP": (
        "CCCCCCCC/C=C\\CCCCCCCC(=O)OCC(C[N+](C)(C)C)OC(=O)CCCCCCC/C=C\\CCCCCCCC",
        "1,2-dioleoyl-3-trimethylammonium-propane"),
    "DDAB": (
        "CCCCCCCCCCCCCCCCCC[N+](C)(C)CCCCCCCCCCCCCCCCCC",
        "dimethyldioctadecylammonium"),
    "Cholesterol": (
        "CC(C)CCC[C@@H](C)[C@H]1CC[C@H]2[C@@H]3CC=C4C[C@@H](O)CC[C@]4(C)[C@H]3CC[C@]12C",
        "cholest-5-en-3beta-ol"),
    "18PG": (
        "CCCCCCCCCCCCCCCCCC(=O)OC[C@H](COP(=O)([O-])OCC(O)CO)OC(=O)CCCCCCC/C=C\\CCCCCCCC",
        "1-stearoyl-2-oleoyl-sn-glycero-3-phospho-(1'-rac-glycerol), SOPG; named "
        "verbatim in YZ_2022 Materials (Nat Commun 2022, PMC9310361). NOT 18:0/18:0."),
    "14PA": (
        "CCCCCCCCCCCCCC(=O)OC[C@H](COP(=O)([O-])O)OC(=O)CCCCCCCCCCCCC",
        "1,2-dimyristoyl-sn-glycero-3-phosphate, DMPA; same YZ_2022 sentence"),
    "DOPG": (
        "CCCCCCCC/C=C\\CCCCCCCC(=O)OC[C@H](COP(=O)([O-])OC[C@@H](O)CO)OC(=O)CCCCCCC/C=C\\CCCCCCCC",
        "1,2-dioleoyl-sn-glycero-3-phospho-(1'-sn-glycerol)"),
    "DLin-MC3-DMA": (
        "CCCCC/C=C\\C/C=C\\CCCCCCCCC(OC(=O)CCCN(C)C)CCCCCCCC/C=C\\C/C=C\\CCCCC",
        "(6Z,9Z,28Z,31Z)-heptatriaconta-6,9,28,31-tetraen-19-yl "
        "4-(dimethylamino)butanoate; ACYCLIC ESTER, no ring"),
    "KC2": (
        "CCCCC/C=C\\C/C=C\\CCCCCCCCC1(CCCCCCCC/C=C\\C/C=C\\CCCCC)OCC(CCN(C)C)O1",
        "DLin-KC2-DMA, 2,2-dilinoleyl-4-(2-dimethylaminoethyl)-1,3-dioxolane; "
        "RING, no ester. Same formula and exact mass as MC3 -- the pair that "
        "motivated this check."),
    "SM-102": (
        "CCCCCCCCC(CCCCCCCC)OC(=O)CCCCCCCN(CCO)CCCCCC(=O)OCCCCCCCCCCC",
        "heptadecan-9-yl 8-[(2-hydroxyethyl)(6-oxo-6-undecoxyhexyl)amino]octanoate"),
    "ALC-0315": (
        "CCCCCCCCC(CCCCCC)C(=O)OCCCCCCN(CCCCO)CCCCCCOC(=O)C(CCCCCC)CCCCCCCC",
        "((4-hydroxybutyl)azanediyl)bis(hexane-6,1-diyl) bis(2-hexyldecanoate)"),
}

ALIASES = {"Dlin-MC3-DMA": "DLin-MC3-DMA", "DLin-KC2-DMA": "KC2"}
PREFIXES = ("IL", "HL", "CHL", "PEG", "fifthcomponent")


def skeleton(smiles):
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None
    return Chem.MolToInchiKey(m).split("-")[0]


def main():
    want = {}
    for name, (smi, why) in REGISTRY.items():
        sk = skeleton(smi)
        if sk is None:
            print(f"  REGISTRY ENTRY DOES NOT PARSE: {name}")
            return 1
        want[name] = (sk, Chem.MolToSmiles(Chem.MolFromSmiles(smi)), why)

    seen = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    for fn in ("data_vitro.csv", "data_vivo.csv"):
        path = os.path.join(DATA, fn)
        if not os.path.exists(path):
            continue
        for r in csv.DictReader(open(path)):
            for pre in PREFIXES:
                raw = (r.get(f"{pre}_name") or "").strip()
                smi = r.get(f"{pre}_SMILES") or ""
                name = ALIASES.get(raw, raw)
                if name in want and smi:
                    seen[name][smi][r["Experiment_ID"]] += 1

    results, n_bad, n_ok, n_split = [], 0, 0, 0
    for name, (sk, canon_smi, why) in sorted(want.items()):
        if name not in seen:
            results.append({"reagent": name, "status": "absent_from_corpus"})
            continue
        variants = []
        for smi, exps in seen[name].items():
            ok = skeleton(smi) == sk
            variants.append({"smiles": smi, "matches_registry": ok,
                             "rows": sum(exps.values()),
                             "experiments": sorted(exps)})
        bad = [v for v in variants if not v["matches_registry"]]
        entry = {"reagent": name, "authority": why,
                 "registry_smiles": canon_smi, "variants": variants}
        if bad:
            entry["status"] = "STRUCTURE_MISMATCH"
            n_bad += sum(v["rows"] for v in bad)
        elif len(variants) > 1:
            entry["status"] = "same_skeleton_multiple_strings"
            n_split += 1
        else:
            entry["status"] = "ok"
            n_ok += 1
        results.append(entry)

    json.dump({"what": "named commercial reagents vs their fixed structures",
               "registered": len(REGISTRY), "rows_mismatched": n_bad,
               "results": results}, open(OUT, "w"), indent=1)

    print(f"named-reagent registry: {len(REGISTRY)} reagents")
    for e in results:
        tag = {"ok": "OK  ", "STRUCTURE_MISMATCH": "FAIL",
               "same_skeleton_multiple_strings": "SPLIT",
               "absent_from_corpus": "-   "}[e["status"]]
        rows = sum(v["rows"] for v in e.get("variants", []))
        print(f"  {tag:6} {e['reagent']:16} {rows:6} rows   {e['status']}")
        for v in e.get("variants", []):
            if not v["matches_registry"]:
                print(f"         wrong structure in {v['experiments']} "
                      f"({v['rows']} rows): {v['smiles'][:72]}")
    print(f"\nreagents clean: {n_ok}   mismatched rows: {n_bad}")
    return 1 if n_bad else 0


if __name__ == "__main__":
    sys.exit(main())
