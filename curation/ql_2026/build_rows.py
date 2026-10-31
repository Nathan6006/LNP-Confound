#!/usr/bin/env python3
"""Emit the 344 Liu et al. lipids in the 72-column corpus schema.
Merged into data/ on 2026-10-30; see provenance/CHANGELOG.md section 11.

Experiment_ID QL_2026 (Qiang Liu, first author). Adv Healthc Mater 2026, 15,
e02225, DOI 10.1002/adhm.202502225.

VALUE SOURCE, per lipid:
  library 1 (320)  Figure 2c where it plots the lipid - a log axis, so precise
                   and uncensored - otherwise the Figure 2b heatmap.
  library 2 (24)   Figure 3c, the only place these appear.

HOW THE TWO LIBRARIES ARE PUT ON ONE SCALE. Four benchmark LNPs (DLin-MC3-DMA,
SM102, ALC-0315, Lipid 5) are plotted in BOTH Fig 2c and Fig 3c, so the ratio
between their readings measures the offset between the two screens. Library 2 is
divided by that ratio and all 344 are then z-scored ONCE, together.

THE OFFSET IS NOT DISTINGUISHABLE FROM 1.0. Corrected 2026-10-30, review
finding F-018. The measured offset is **1.0071, sd 0.0247** (DLin 1.011,
SM102 1.035, ALC-0315 0.967, Lipid 5 1.015) -- the two screens agree to within
+/-2.5%, which is the spread of the four benchmarks themselves.

The earlier "SCALAR ~1.1x" (1.1116, later 1.1069) was an artifact of the bar-top
detector reading replicate dots instead of the bar top: ALC-0315 was
dot-contaminated in BOTH panels (+24 px in Fig 2c, +16 px in Fig 3c) and Lipid 5
in Fig 3c (+8 px), and those are three of the four benchmarks the offset is
measured from. So the correction is not a refinement of the offset, it removes
it: the quantity was mostly measuring our own detector error.

The division is kept rather than dropped, because 1.0071 is the honest estimate
and hard-coding 1.0 would assert an agreement the figures do not quite show. The
sd is reported beside it so a reader can see it covers 1.0.

An earlier version standardised the two libraries separately. That was wrong
twice over. The offset it cited, 1.54x, was an artifact of a miscalibrated Fig
3c axis; and separate standardisation does not merely remove an offset, it
forces library 2's spread onto library 1's and deletes the ~40x gap between
them - which is the paper's central finding, that six-tail lipids beat
three-tail ones. It inverted the raw ordering of 36% of cross-library pairs,
while every automated check still passed, because per-block mean 0 / sd 1 is a
tautology for concatenated standardised blocks.

FLAGS carried on the rows:
  value_censored   the Fig 2b cell is saturated at the colorbar maximum and
                   Fig 2c does not plot it, so only a lower bound is known (3).
  value_conflict   Fig 2b and Fig 2c disagree about the lipid (1: T3A-L8-C14,
                   cream in the heatmap at 7.9e3 but plotted at 5.8e4 in 2c).
                   The HEATMAP value is carried: the cell is visibly cream and
                   the heatmap is the complete screen. An earlier version said
                   so in prose while the code silently took the Fig 2c value.
  value_floored    the Fig 2b cell sits at the colorbar floor, where the
                   colormap is degenerate (runs of identical RGB span up to
                   ~980 RLU), so the value is an upper bound within its band
                   (113 of 320).
"""
import csv
import json
import math
import os
import statistics as st

from rdkit import Chem, RDLogger
from rdkit.Chem import Crippen, Descriptors as D, rdMolDescriptors as rd

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
CORPUS = os.path.join(HERE, "..", "..", "data", "data_vitro.csv")
OUT = os.path.join(HERE, "11_QL_2026_rows.csv")

EID = "QL_2026"
DOI = "10.1002/adhm.202502225"
LINK = "https://doi.org/10.1002/adhm.202502225"
# Verified 2026-10-29 against NCBI esummary: PMID 41145968 carries DOI
# 10.1002/adhm.202502225 and first author Liu Q, so the PMID, the DOI and the
# initials-and-year convention all agree. An earlier version left this blank on
# the belief the article was not yet indexed; it is.
PMID = "41145968"

MOLPCT = {"IL": 35.0, "HL": 16.0, "CHL": 46.5, "PEG": 2.5}

ESTER = Chem.MolFromSmarts("[CX3](=O)[OX2H0][#6]")
CARBONATE = Chem.MolFromSmarts("[OX2][CX3](=O)[OX2]")
DISULFIDE = Chem.MolFromSmarts("[#16X2][#16X2]")
AMIDE = Chem.MolFromSmarts("[NX3][CX3]=[OX1]")


def n_protonatable(m):
    am = {t[0] for t in m.GetSubstructMatches(AMIDE)}
    n = 0
    for a in m.GetAtoms():
        if a.GetSymbol() != "N":
            continue
        hv = sum(1 for nb in a.GetNeighbors() if nb.GetAtomicNum() > 1)
        if a.GetIsAromatic():
            if a.GetTotalNumHs() == 0 and hv == 2:
                n += 1
        elif a.GetFormalCharge() == 0 and a.GetIdx() not in am and hv < 4:
            n += 1
    return n


def longest_chain(m):
    idx = [a.GetIdx() for a in m.GetAtoms() if a.GetSymbol() == "C" and not a.IsInRing()]
    S = set(idx)
    adj = {i: [nb.GetIdx() for nb in m.GetAtomWithIdx(i).GetNeighbors() if nb.GetIdx() in S]
           for i in idx}
    best = 0
    for start in idx:
        stack = [(start, frozenset([start]))]
        while stack:
            u, vis = stack.pop()
            best = max(best, len(vis))
            for v in adj[u]:
                if v not in vis:
                    stack.append((v, vis | {v}))
    return best


def descriptors(smiles):
    m = Chem.MolFromSmiles(smiles)
    return {
        "Heavy.Atoms": m.GetNumHeavyAtoms(), "Rings": rd.CalcNumRings(m),
        "Aromatic.Rings": rd.CalcNumAromaticRings(m),
        "Rotatable.Bonds": rd.CalcNumRotatableBonds(m),
        "ExactMolWt": D.ExactMolWt(m),
        "Topological.Polar.Surface.Area": round(rd.CalcTPSA(m), 2),
        "Hydrogen.Bond.Donors": D.NumHDonors(m),
        "Hydrogen.Bond.Acceptors": D.NumHAcceptors(m),
        "LogP": Crippen.MolLogP(m), "Molar.Refractivity": Crippen.MolMR(m),
        "Fraction.sp3.Carbons": rd.CalcFractionCSP3(m),
        "sp3.Carbons": sum(1 for a in m.GetAtoms()
                           if a.GetSymbol() == "C" and str(a.GetHybridization()) == "SP3"),
        "Nitrogen.Count": sum(1 for a in m.GetAtoms() if a.GetSymbol() == "N"),
        "Molecular.Weight": round(D.MolWt(m), 2),
        "has_ester": int(m.HasSubstructMatch(ESTER)),
        "has_carbonate": int(m.HasSubstructMatch(CARBONATE)),
        "has_disulfide": int(m.HasSubstructMatch(DISULFIDE)),
        "Num_carbon_in_tail": longest_chain(m),
        "num_protonatable_nitrogens": n_protonatable(m),
        "num_permanent_cationic_N": sum(
            1 for a in m.GetAtoms() if a.GetSymbol() == "N" and a.GetFormalCharge() == 1
            and sum(1 for nb in a.GetNeighbors() if nb.GetAtomicNum() > 1) == 4),
        "formal_net_charge": sum(a.GetFormalCharge() for a in m.GetAtoms()),
        "num_unsaturated_cc_bonds": sum(
            1 for b in m.GetBonds()
            if b.GetBondType() == Chem.BondType.DOUBLE and not b.GetIsAromatic()
            and b.GetBeginAtom().GetSymbol() == "C" and b.GetEndAtom().GetSymbol() == "C"),
    }


def row_flags(r):
    """The per-row quality flags for one library member.

    `value_is_bound` is `upper` on a floored cell because the whole reason the
    cell is flagged is that the colormap bottoms out: identical RGB spans about
    980 RLU there, so every true value from the ramp's white point up to the
    assigned number renders the same colour and the assigned number is the TOP
    of the band the truth could occupy. It is `lower` on a censored cell for the
    mirror reason -- the cell is saturated at the colorbar maximum, so the truth
    is at or above what we stored.

    113 of these 344 rows are floored, i.e. a third of library 1. That is not a
    small caveat and it is the dominant uncertainty in the slice; see F-019.
    """
    out = {}
    if r["floored"]:
        out["value_is_bound"] = "upper"
        out["curation_note"] = (
            "Fig 2b heatmap cell sits on the colormap floor. The ramp reaches "
            "pure white at about 948 RLU with no gradient below it, so every "
            "true value from 0 to there renders the same colour: the stored "
            "value is the UPPER BOUND of a [0, ~982] RLU band, not a "
            "measurement. Adjacent LUT steps near the floor differ by 1.7-3.5 "
            "RGB units while each spans 349-1045 RLU, so the limit is the "
            "colormap's own quantisation, not image noise. 113 of 344 rows are "
            "in this condition and 56% of the slice carries one of only six "
            "values, all from the bottom six LUT steps. See F-019.")
    elif r["censored"]:
        out["value_is_bound"] = "lower"
        out["curation_note"] = (
            "Fig 2b heatmap cell is saturated at the colorbar maximum and the "
            "lipid is absent from Fig 2c, so the stored value is a lower bound.")
    if r["conflict"]:
        out["curation_note"] = (
            "Fig 2b and Fig 2c disagree about this lipid by 7.3x. The heatmap "
            "value is stored, because Fig 2c's bar is dot-contaminated and the "
            "heatmap is the complete 320-member panel; see PROVENANCE.md.")
    return out


def main():
    enum = {r["lipid_name"]: r for r in csv.DictReader(open("09_enumerated_344.csv"))}
    blocks = {r["name_in_paper"]: r["SMILES"]
              for r in csv.DictReader(open("liu2026_adhm202502225_blocks.csv"))}
    heat = {r["lipid_name"]: r for r in csv.DictReader(open("01_fig2b_320.csv"))}
    bars2c = {r["lipid_name"]: r for r in csv.DictReader(open("07_fig2c_bars.csv"))
              if r["is_benchmark"] == "0"}
    bars3c = {r["lipid_name"]: r for r in csv.DictReader(open("05_fig3c_24.csv"))
              if r["is_benchmark"] == "0"}

    recs = []
    for name, e in enum.items():
        if e["library"] == "1":
            h = heat[name]
            conflict = int(name in bars2c and h["saturated"] == "0"
                           and float(bars2c[name]["RLU"]) / max(float(h["RLU"]), 1) > 3)
            if name in bars2c and not conflict:
                value, src = float(bars2c[name]["RLU"]), "Fig 2c"
            else:
                # conflict rows fall back to the heatmap, as documented
                value, src = float(h["RLU"]), "Fig 2b"
            censored = int(h["saturated"] == "1" and name not in bars2c)
            floored = int(h["floored"] == "1")
        else:
            value, src = float(bars3c[name]["RLU"]), "Fig 3c"
            censored = conflict = floored = 0
        recs.append({"name": name, "lib": e["library"], "smiles": e["IL_SMILES"],
                     "value": value, "src": src,
                     "censored": censored, "conflict": conflict,
                     "floored": floored})

    # cross-figure offset, measured on the four benchmarks plotted in both
    bm2 = {r["lipid_name"]: float(r["RLU"]) for r in
           csv.DictReader(open("07_fig2c_bars.csv")) if r["is_benchmark"] == "1"}
    bm3 = {r["lipid_name"]: float(r["RLU"]) for r in
           csv.DictReader(open("05_fig3c_24.csv")) if r["is_benchmark"] == "1"}
    ratios = [bm3[k] / bm2[k] for k in bm2]
    offset = st.mean(ratios)

    # put library 2 on library 1's scale, then z-score all 344 together
    for r in recs:
        r["value_scaled"] = r["value"] / offset if r["lib"] == "2" else r["value"]
    logs = [math.log10(r["value_scaled"]) for r in recs]
    mu, sd = st.mean(logs), st.pstdev(logs)
    for r, L in zip(recs, logs):
        r["z"] = (L - mu) / sd

    corpus = list(csv.DictReader(open(CORPUS)))
    hdr = list(corpus[0].keys())
    for c in (x["name"] for x in json.load(open(os.path.join(
            HERE, "..", "..", "provenance", "schema_quality.json")))["columns"]):
        if c not in hdr:
            hdr.append(c)
    comp = {}
    for r in corpus:
        for pre in ("HL", "CHL", "PEG"):
            n = (r[f"{pre}_name"] or "").strip()
            if n and r[f"{pre}_SMILES"]:
                comp.setdefault(n, r[f"{pre}_SMILES"])
    start = 1 + max(int(r["Index"]) for r in corpus if r["Index"].isdigit())

    rows = []
    for i, r in enumerate(recs):
        row = dict.fromkeys(hdr, "")
        row.update({
            "Index": start + i, "LNP_ID": f"LNP_{start + i:07d}",
            "Experiment_ID": EID, "Formulation_ID": f"{EID}_lib{r['lib']}",
            "IL_name": r["name"], "IL_SMILES": r["smiles"],
            "IL_molratio": MOLPCT["IL"],
            # The paper says "20:1" but calls it "lipid/mRNA" in the Results and
            # "LNP/mRNA" in the Methods, and never defines it as ionizable
            # lipid to mRNA. If it is TOTAL lipid, the true IL:mRNA runs 4.4-12.5
            # across this library (the IL mass fraction varies 0.22-0.63 as the
            # lipid MW varies), so writing 20.0 here would assert something the
            # paper does not say. Left blank; see PROVENANCE.md.
            "IL_to_nucleicacid_massratio": "",
            "HL_name": "DOPE", "HL_SMILES": comp.get("DOPE", ""),
            "HL_molratio": MOLPCT["HL"],
            "CHL_name": "Cholesterol", "CHL_SMILES": comp.get("Cholesterol", ""),
            "CHL_molratio": MOLPCT["CHL"],
            # The SI Materials list names a phosphoethanolamine-PEG2000 from
            # Sigma, while the main text, Methods and every figure say
            # DMG-PEG 2000. The corpus distinguishes the two, so the figures'
            # name is used and the contradiction recorded in PROVENANCE.md.
            "PEG_name": "DMG-PEG2000", "PEG_SMILES": comp.get("DMG-PEG2000", ""),
            "PEG_molratio": MOLPCT["PEG"],
            "Aqueous_buffer": "citrate", "Mixing_method": "handmixed",
            "Model": "in_vitro", "Model_type": "HeLa", "Model_target": "in_vitro",
            "Route_of_administration": "in_vitro",
            "Cargo": "mRNA", "Cargo_type": "FLuc",
            "Dose_ug_nucleicacid": 0.020,
            "Experiment_method": "luminescence_normalized",
            "Experiment_batching": "individual",
            "Experiment_value": r["z"],
            "Publication_link": LINK, "Publication_DOI": DOI,
            "Publication_PMID": PMID,
        })
        # Quality flags onto the ROW, 2026-10-30 (review finding F-001). These
        # three counts were already computed and already written into
        # 12_QL_2026_qc.json -- and then dropped, because the 72-column schema
        # had nowhere to put them. A consumer of data_vitro.csv could not tell
        # a colormap floor from a measurement. See provenance/schema_quality.json.
        row.update(row_flags(r))
        # decomposition: this is a clean three-component library, so record it
        e = enum[r["name"]]
        row.update({
            "IL_head_name": e["amine"], "IL_head_SMILES": blocks[e["amine"]],
            "IL_linker_name": e["bl"], "IL_linker_SMILES": blocks[e["bl"]],
            "IL_tail1_name": e["tail"], "IL_tail1_SMILES": blocks[e["tail"]],
        })
        row.update(descriptors(r["smiles"]))
        rows.append(row)

    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=hdr)
        w.writeheader()
        w.writerows(rows)

    allz = [r["z"] for r in recs]
    qc = {
        "experiment_id": EID, "rows": len(rows), "columns": len(hdr),
        "paper": "Liu et al., Adv Healthc Mater 2026, 15, e02225",
        "doi": DOI,
        "library1": sum(1 for r in recs if r["lib"] == "1"),
        "library2": sum(1 for r in recs if r["lib"] == "2"),
        "value_source": {s: sum(1 for r in recs if r["src"] == s)
                         for s in ("Fig 2b", "Fig 2c", "Fig 3c")},
        "cross_figure_offset": {
            "benchmark_ratios_fig3c_over_fig2c": {k: round(bm3[k] / bm2[k], 3)
                                                  for k in bm2},
            "applied": round(offset, 4), "sd": round(st.pstdev(ratios), 4)},
        "flags": {"value_censored": sum(r["censored"] for r in recs),
                  "value_conflict": sum(r["conflict"] for r in recs),
                  "value_floored": sum(r["floored"] for r in recs),
                  "censored_lipids": [r["name"] for r in recs if r["censored"]],
                  "conflict_lipids": [r["name"] for r in recs if r["conflict"]]},
        "geometric_mean_RLU": {
            lib: round(10 ** st.mean([math.log10(r["value"]) for r in recs
                                      if r["lib"] == lib]), 1)
            for lib in ("1", "2")},
        "z": {"mean": round(st.mean(allz), 12), "pop_sd": round(st.pstdev(allz), 12),
              "pooled_once_over_all_344": True,
              "per_library_mean": {lib: round(st.mean([r["z"] for r in recs
                                                       if r["lib"] == lib]), 4)
                                   for lib in ("1", "2")}},
        "molratio_sum": round(sum(MOLPCT.values()), 2),
        "conditions": "HeLa, 1e4 cells/well, 96-well; mFluc mRNA 20 ng/well; 24 h; "
                      "n=3; IL/DOPE/Chol/DMG-PEG2000 = 35/16/46.5/2.5 mol%; "
                      "lipid:mRNA 20:1 w/w; citrate pH 4.0; pipette mixing",
        "merged": "2026-10-30 into data/data_vitro.csv (28,541 -> 28,885 "
                  "rows, 41 -> 42 experiments); see provenance/CHANGELOG.md "
                  "section 11. Round-1 review findings open against this "
                  "slice: provenance/review/round_2026-10-30/FINDINGS.md",
    }
    json.dump(qc, open("12_QL_2026_qc.json", "w"), indent=1)

    print(f"{os.path.basename(OUT)}: {len(rows)} rows x {len(hdr)} cols")
    print(f"  library 1 {qc['library1']}, library 2 {qc['library2']}")
    print(f"  value source {qc['value_source']}")
    print(f"  cross-figure offset applied to library 2: "
          f"{qc['cross_figure_offset']['applied']}x "
          f"(sd {qc['cross_figure_offset']['sd']})")
    print(f"  geometric mean RLU: library 1 {qc['geometric_mean_RLU']['1']}, "
          f"library 2 {qc['geometric_mean_RLU']['2']}")
    print(f"  flags censored {qc['flags']['value_censored']}, floored "
          f"{qc['flags']['value_floored']}, conflict "
          f"{qc['flags']['value_conflict']} {qc['flags']['conflict_lipids']}")
    print(f"  per-library mean z: {qc['z']['per_library_mean']}")
    print(f"  z overall mean {qc['z']['mean']:.2e} sd {qc['z']['pop_sd']:.12f}")
    print(f"  molratios sum {qc['molratio_sum']}")


if __name__ == "__main__":
    main()
