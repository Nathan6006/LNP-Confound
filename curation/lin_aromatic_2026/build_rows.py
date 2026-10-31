#!/usr/bin/env python3
"""Emit the Lin et al. 840-lipid screen in the 67-column corpus schema.

Experiment_ID LL_2026 (Lixin Lin, first author). No collision: the corpus has
LL_2012 in vitro and LL_2025 in vivo.

Conditions, all from the SI Experimental Procedures rather than inferred:
  cells        IGROV1, 96-well, 1.5e4 per well, 24 h, n = 3
  cargo        firefly luciferase mRNA, 25 ng per well  (= 0.025 ug)
  formulation  lipid / DOPE / cholesterol / DMG-PEG2000 = 15/15/25/2
               (26.3 / 26.3 / 43.9 / 3.5 mol%)
  IL : mRNA    10 : 1 by mass
  buffer       citrate 10 mM pH 4.4, rapid 1:3 volumetric mixing (handmixed)

Experiment_value follows the corpus convention, reverse-engineered and verified
against SL_2021 in provenance/build_dataset.py: a per-experiment z-score, and
for a DISCRETIZED source the quantity z-scored is the BIN ORDINAL, not a
midpoint. Fig 2e is 4-binned, so ordinals 0-3 are z-scored with the population
sd.
"""
import csv
import json
import os
import statistics as st

from rdkit import Chem, RDLogger
from rdkit.Chem import Crippen, Descriptors as D, rdMolDescriptors as rd

RDLogger.DisableLog("rdApp.*")

CORPUS = os.path.join("..", "..", "data", "data_vitro.csv")
OUT = "10_LL_2026_rows.csv"

EID = "LL_2026"
LINK = "https://pubs.acs.org/doi/10.1021/jacs.5c16971"
DOI = "10.1021/jacs.5c16971"

MOLPCT = {"IL": 26.3, "HL": 26.3, "CHL": 43.9, "PEG": 3.5}

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


def quality(r):
    """Carry this lipid's own flags onto the row (F-001), plus the arm
    decomposition that 03_enumerated.csv has always held (F-010).

    `positional_ambiguity` is 1 on 987 of these rows and it is not a technicality:
    the library was made in one pot and purified without chromatography, so WHICH
    nitrogen bears which arm is a convention we chose, not something the paper
    measured. `confidence = medium` on the 140 Ac8-1/Ac10-1 rows because the SI
    names commercial iso-octyl and iso-decyl acrylate -- branched isomer
    MIXTURES sold under one name -- so the single stored SMILES is one member of
    a mixture. Both facts were computed by enumerate_library.py, written into
    03_enumerated.csv, and then dropped at the emit step.
    """
    out = {}
    if str(r.get("positional_ambiguity")) == "1":
        out["positional_ambiguity"] = "1"
        out["curation_note"] = (
            "One-pot library, no chromatographic separation of regioisomers: "
            "which ring nitrogen bears which arm is a curation convention, not "
            "a measurement.")
    if r.get("confidence"):
        out["confidence"] = r["confidence"]
    if r.get("confidence") == "medium":
        out["structure_assumed"] = "1"
        out["curation_note"] = (
            (out.get("curation_note", "") + " ") if out.get("curation_note") else ""
        ) + ("Tail reagent is a commercial branched isomer mixture (iso-octyl / "
             "iso-decyl acrylate); the stored SMILES is one isomer of that mixture.")
    return out


def main():
    corpus = list(csv.DictReader(open(CORPUS)))
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

    lipids = list(csv.DictReader(open("06_merged.csv")))
    ordinals = [int(r["bin_index"]) - 1 for r in lipids]
    mu, sd = st.mean(ordinals), st.pstdev(ordinals)
    zs = [(o - mu) / sd for o in ordinals]

    rows = []
    for i, (r, z) in enumerate(zip(lipids, zs)):
        row = dict.fromkeys(hdr, "")
        row.update({
            "Index": start + i, "LNP_ID": f"LNP_{start + i:07d}",
            "Experiment_ID": EID, "Formulation_ID": f"{EID}_F1",
            "IL_name": r["lipid_name"], "IL_SMILES": r["IL_SMILES"],
            "IL_molratio": MOLPCT["IL"], "IL_to_nucleicacid_massratio": 10.0,
            "HL_name": "DOPE", "HL_SMILES": comp.get("DOPE", ""),
            "HL_molratio": MOLPCT["HL"],
            "CHL_name": "Cholesterol", "CHL_SMILES": comp.get("Cholesterol", ""),
            "CHL_molratio": MOLPCT["CHL"],
            "PEG_name": "DMG-PEG2000", "PEG_SMILES": comp.get("DMG-PEG2000", ""),
            "PEG_molratio": MOLPCT["PEG"],
            "Aqueous_buffer": "citrate", "Dialysis_buffer": "",
            "Mixing_method": "handmixed",
            "Model": "in_vitro", "Model_type": "IGROV1", "Model_target": "in_vitro",
            "Route_of_administration": "in_vitro",
            "Cargo": "mRNA", "Cargo_type": "FLuc",
            "Dose_ug_nucleicacid": 0.025,
            "Experiment_method": "luminescence_discretized_normalized",
            "Experiment_batching": "individual",
            "Experiment_value": z,
            "Publication_link": LINK, "Publication_DOI": DOI,
            "Publication_PMID": "41689530",
        })
        # arm decomposition: 03_enumerated.csv has carried it all along and the
        # corpus has the columns, but they were left blank on all 987 rows while
        # 16,471 of 28,885 other vitro rows do carry them (F-010).
        row.update({
            "IL_head_name": r.get("amine_name") or r.get("amine", ""),
            "IL_tail1_name": r.get("aromatic_name") or r.get("aromatic", ""),
            "IL_tail2_name": r.get("tail_name") or r.get("tail", ""),
        })
        row.update(descriptors(r["IL_SMILES"]))
        row.update(quality(r))
        rows.append(row)

    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=hdr)
        w.writeheader()
        w.writerows(rows)

    tot = sum(MOLPCT.values())
    qc = {
        "experiment_id": EID, "rows": len(rows), "columns": len(hdr),
        "paper": "Lin et al., JACS 2026, 148, 6960-6969",
        "molratio_sum": round(tot, 2),
        "z": {"mean": round(st.mean(zs), 12), "pop_sd": round(st.pstdev(zs), 12),
              "distinct": len({round(x, 9) for x in zs})},
        "conditions": "IGROV1 96-well 1.5e4/well 24 h n=3; Fluc mRNA 25 ng/well; "
                      "IL:mRNA 10:1 w/w; citrate 10 mM pH 4.4; 1:3 rapid mixing",
        "validation": {
            "HRMS_anchors": "4/4 exact to +/-0.0001 Da "
                            "(10A2-3B-EP10, 9A2-3B-EP10, 9A2-4B-EP10, 10A2-3B-EP12)",
            "fig2fgh_hit_rates": "48/48 printed hit rates reproduced from a TICK-ANCHORED scale, max |delta| 1.0 pp, mean -0.52 pp. Corrected 2026-10-30 (F-006): until then the pixels-per-percent scale was least-squares fitted onto the Fig 2e rates this panel validates, so the agreement could not detect a scale error AND was absorbing a real one -- a systematic ~1% top-edge bias from the colour mask. Sub-pixel top edges plus a tick-only scale; the residual -0.52 pp one-sided bias is recorded, not tuned away",
            "colour_match": "max distance 2.24, zero ambiguous cells",
        },
    }
    json.dump(qc, open("11_LL_2026_qc.json", "w"), indent=1)
    print(f"{OUT}: {len(rows)} rows x {len(hdr)} cols")
    print(f"  molratios sum to {tot}")
    print(f"  z mean {qc['z']['mean']:.2e} pop sd {qc['z']['pop_sd']}  "
          f"{qc['z']['distinct']} distinct values")


if __name__ == "__main__":
    main()
