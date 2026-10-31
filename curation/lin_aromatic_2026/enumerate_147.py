#!/usr/bin/env python3
"""Enumerate the 147 initial-screen lipids and emit them in the corpus schema.

Reuses the reaction machinery in enumerate_library.py unchanged, so the 147 and
the 840 are built by identical chemistry; only the block sets differ (all twelve
amines 1A-12A here, EP10 as the sole tail).

Experiment_ID LL_2026_2, following the corpus convention for a second screen
from one paper (LX_2024 / LX_2024_2 / LX_2024_3, XH_2024 / XH_2024_2). Kept
separate from LL_2026 because the formulation differs: 15/20/25/2
(24.2/32.3/40.3/3.2 mol%) at IL:mRNA 20:1, the initial formulation, versus
15/15/25/2 at 10:1 for the 840.
"""
import csv
import json
import os
import statistics as st
import sys

# run with python -I, which drops the script directory from sys.path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from rdkit import Chem, RDLogger

import enumerate_library as EL
import build_rows as BR

RDLogger.DisableLog("rdApp.*")

EID = "LL_2026_2"
MOLPCT = {"IL": 24.2, "HL": 32.3, "CHL": 40.3, "PEG": 3.2}
OUT = "12_LL_2026_2_rows.csv"


def main():
    cells = list(csv.DictReader(open("04b_fig2a_147.csv")))
    rows, fails = [], []
    for c in cells:
        smi, n_tails, n_sites = EL.build(c["amine"], c["aromatic"], "EP10",
                                         int(c["equiv"]))
        if smi is None:
            fails.append(c["lipid_name"])
            continue
        rows.append({**c, "IL_SMILES": smi, "n_NH_sites": n_sites,
                     "n_tail_arms": n_tails,
                     "positional_ambiguity": int(int(c["equiv"]) < n_sites)})
    assert not fails, f"failed to build: {fails}"

    with open("13_enumerated_147.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    corpus = list(csv.DictReader(open(BR.CORPUS)))
    hdr = list(corpus[0].keys())
    for c in (x["name"] for x in json.load(open(os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            "..", "..", "provenance", "schema_quality.json")))["columns"]):
        if c not in hdr:
            hdr.append(c)
    comp = {}
    for r in corpus:
        for pre in ("HL", "PEG", "CHL"):
            if r[f"{pre}_name"] and r[f"{pre}_SMILES"]:
                comp.setdefault(r[f"{pre}_name"], r[f"{pre}_SMILES"])
    start = 1 + max(int(r["Index"]) for r in corpus if r["Index"].isdigit())

    ordinals = [int(r["bin_index"]) - 1 for r in rows]
    mu, sd = st.mean(ordinals), st.pstdev(ordinals)
    zs = [(o - mu) / sd for o in ordinals]

    out = []
    for i, (r, z) in enumerate(zip(rows, zs)):
        row = dict.fromkeys(hdr, "")
        row.update({
            "Index": start + i, "LNP_ID": f"LNP_{start + i:07d}",
            "Experiment_ID": EID, "Formulation_ID": f"{EID}_F1",
            "IL_name": r["lipid_name"], "IL_SMILES": r["IL_SMILES"],
            "IL_molratio": MOLPCT["IL"], "IL_to_nucleicacid_massratio": 20.0,
            "HL_name": "DOPE", "HL_SMILES": comp.get("DOPE", ""),
            "HL_molratio": MOLPCT["HL"],
            "CHL_name": "Cholesterol", "CHL_SMILES": comp.get("Cholesterol", ""),
            "CHL_molratio": MOLPCT["CHL"],
            "PEG_name": "DMG-PEG2000", "PEG_SMILES": comp.get("DMG-PEG2000", ""),
            "PEG_molratio": MOLPCT["PEG"],
            "Aqueous_buffer": "citrate", "Mixing_method": "handmixed",
            "Model": "in_vitro", "Model_type": "IGROV1", "Model_target": "in_vitro",
            "Route_of_administration": "in_vitro",
            "Cargo": "mRNA", "Cargo_type": "FLuc",
            "Dose_ug_nucleicacid": 0.025,
            "Experiment_method": "luminescence_discretized_normalized",
            "Experiment_batching": "individual",
            "Experiment_value": z,
            "Publication_link": BR.LINK, "Publication_DOI": BR.DOI,
            "Publication_PMID": "41689530",
        })
        row.update({
            "IL_head_name": r.get("amine", ""),
            "IL_tail1_name": r.get("aromatic", ""),
            "IL_tail2_name": "EP10",
        })
        row.update(BR.descriptors(r["IL_SMILES"]))
        # Same flags as the 840-row slice (F-001). `confidence` is set here
        # rather than carried, because this file has no confidence column of its
        # own: the 840-row rule assigns `medium` only to the Ac8-1/Ac10-1 tails
        # (commercial branched isomer mixtures) and `high` otherwise, and every
        # lipid here is built on the single defined EP10 tail. So it is derived
        # by that same rule, not invented -- stated because the sibling file
        # states it for identical chemistry and a silent blank would read as
        # "not assessed".
        row.update(BR.quality({**r, "confidence": "high"}))
        out.append(row)

    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=hdr)
        w.writeheader()
        w.writerows(out)

    big = {r["IL_SMILES"] for r in csv.DictReader(open("06_merged.csv"))}
    novel_amines = sorted({r["amine"] for r in rows} - {"8A", "9A", "10A", "11A", "12A"})
    qc = {
        "experiment_id": EID, "rows": len(out), "columns": len(hdr),
        "source": "Lin et al., JACS 2026, Fig 2a (147-lipid initial screen)",
        "molratio_sum": round(sum(MOLPCT.values()), 2),
        "z": {"mean": round(st.mean(zs), 12), "pop_sd": round(st.pstdev(zs), 12),
              "distinct": len({round(x, 9) for x in zs})},
        "unique_SMILES": len({r["IL_SMILES"] for r in rows}),
        "structures_also_in_the_840": len({r["IL_SMILES"] for r in rows} & big),
        "amine_head_groups_new_to_the_corpus": novel_amines,
        "why_separate_experiment_id": "initial formulation 15/20/25/2 at IL:mRNA "
                                      "20:1, vs 15/15/25/2 at 10:1 for the 840",
        "legend_bins": "3 (no lipid fell below 1e4); indexed 2-4 on the same RLU "
                       "scale as Fig 2e so the two screens stay comparable",
    }
    json.dump(qc, open("14_LL_2026_2_qc.json", "w"), indent=1)
    print(f"{OUT}: {len(out)} rows x {len(hdr)} cols")
    print(f"  unique SMILES {qc['unique_SMILES']}, "
          f"{qc['structures_also_in_the_840']} also in the 840")
    print(f"  new amine head groups: {novel_amines}")
    print(f"  molratios sum {qc['molratio_sum']}, z mean {qc['z']['mean']:.2e} "
          f"sd {qc['z']['pop_sd']}")


if __name__ == "__main__":
    main()
