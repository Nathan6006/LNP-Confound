#!/usr/bin/env python3
"""Emit the FALCON datasets (Toh et al., bioRxiv 2026) in the 72-column corpus schema.

SOURCE. github.com/MaoResearchGroup/falcon-lnp-optimization, cloned to ./repo.
`datasets/FALCON_Manuscript_Datasets/*.csv` are the authors' own tables, so there
is no figure digitization here and no reading error to bound. 1,218 formulations
carry 2,544 cell-type measurements between them.

WHAT THE VALUE IS. The repo's `dataset_info.txt` and the paper's Methods agree:
`nLnLE_<CELL>` is Ln(nLE), where

    LE   = treatment luciferase signal / no-treatment signal   (a fold-change)
    nLE  = LE after a per-batch additive correction in log space, computed from
           six internal control formulations (Methods Eq. 2)

The batch correction is an ADDITIVE SHIFT in log space, so Ln(nLE) is an affine
function of the underlying log signal. The corpus convention is a per-experiment
z-score of log10; z-scoring Ln(nLE) gives exactly that, because ln and log10
differ by a constant factor that cancels. Min-max scaling and the output floor
mentioned in the paper are applied inside the ML pipeline, NOT to these files:
`run_FALCON.py` sets `nLnLE_floor = -100` with the comment "nLnLE values are
already blank normalized, so not needed". The published numbers are untouched.

ONE EXPERIMENT PER SOURCE FILE -- corrected 2026-10-30, review finding F-003.

The earlier version used one z-score per (dataset x cell type): fourteen blocks,
on the argument that each cell type was dosed differently and normalized against
its own controls, so a cell type is the largest unit on one measurement scale.
That was wrong, and the authors' own file says so.

`dataset_info.txt` warns that values "cannot be cross compared ACROSS FILES".
It says nothing against comparing cell types WITHIN a file -- and the authors do
exactly that themselves: every multi-output file carries an `nLnLE_<X>_Selectivity`
column, and the Methods define it as on-target minus off-target Ln(nLE). The
solver below verifies that identity to 1e-15 on every file. So a within-file
difference between two cell types is a quantity the authors computed and
published.

Standardizing each cell-type column separately forces every column to mean 0 and
population sd 1, which makes that difference IDENTICALLY ZERO -- not shrunk,
deleted. On file 5 the raw HepG2 - Ramos mean is +4.0342 ln (a factor of 56.5)
and 173 of 180 formulations are HepG2 > Ramos; under per-column z-scoring the
stored mean difference was exactly 0.0000, and Spearman(raw difference, stored
difference) was -0.517, i.e. ANTI-correlated. That is error class 22, and the old
version's own assertion (each block at mean 0 / sd 1) could never catch it
because it is a tautology of standardizing.

One z-score per FILE, pooled over that file's cell-type columns, honours the
across-file warning and keeps the within-file structure. It is also checkable,
which per-column was not: with a single per-file (mu, sd),

    z_a - z_b = (nLnLE_a - nLnLE_b) / sd = Selectivity / sd

exactly -- one constant per file. Asserted below against the authors' printed
column, max residual 1.8e-15. That assertion is the real check; it fails if the
normalization unit is ever changed back.

NOT CURATED as rows: the `*_Selectivity` columns themselves. Selectivity is a
function of Ln(nLE) values already present as rows, so it is derived, not an
independent measurement. It is used here as validation instead.
"""
import csv
import json
import os
import statistics as st

from rdkit import Chem, RDLogger
from rdkit.Chem import Crippen, Descriptors as D, rdMolDescriptors as rd

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "repo", "datasets", "FALCON_Manuscript_Datasets")
CORPUS = os.path.join(HERE, "..", "..", "data", "data_vitro.csv")
OUT = os.path.join(HERE, "01_WT_2026_rows.csv")

LINK = "https://www.biorxiv.org/content/10.64898/2026.06.13.731774v1"
DOI = "10.64898/2026.06.13.731774"
PMID = ""                                  # preprint, not in PubMed as of 2026-10-29

# Per-cell doses, Methods: "150, 25, 5, and 10 ng/well for Ramos B, THP-1 cells,
# HepG2, and INS-1E/DC2.4/C2C12/NIH-3T3, respectively", with one stated
# exception: "For the experiments described in Fig. 6, THP-1 and HepG2 cells were
# treated with mRNA-LNPs at 10 and 2.5 ng/well". Supplementary Table 9 identifies
# dataset 2 as Fig. 6, so dataset 2 alone takes the lower pair.
DOSE_UG = {"Ramos": 0.150, "THP-1": 0.025, "HepG2": 0.005,
           "INS-1E": 0.010, "DC2.4": 0.010, "C2C12": 0.010, "NIH-3T3": 0.010}
DOSE_FIG6 = {"THP-1": 0.010, "HepG2": 0.0025}

# Component structures. Ionizable lipids and cholesterol/PEG/DSPC/DOTAP are taken
# from the corpus so an identical reagent keeps one SMILES across experiments --
# the premise of manuscript Result 3. 18PG is the exception, see NOTE below.
SOPG = ("CCCCCCCCCCCCCCCCCC(=O)OC[C@H](COP(=O)([O-])OCC(O)CO)"
        "OC(=O)CCCCCCC/C=C\\CCCCCCCC")
# NOTE on 18PG. Avanti's "18PG" is NOT a symmetric 18:0 lipid. The paper that
# introduced this reagent shorthand into the corpus (Nat Commun 2022,
# PMC9310361, curated as YZ_2022/YZ_2024) spells it out: "1-stearoyl-2-oleoyl-
# sn-glycero-3-phospho-(1'-rac-glycerol) (18PG)" -- sn-1 stearoyl 18:0, sn-2
# oleoyl 18:1 cis-9. FALCON cites the same vendor and shorthand. The corpus
# currently stores 18PG as 18:0/18:0 (DSPG) on 360 rows, which is wrong; these
# new rows use the correct mixed-chain structure and the discrepancy is reported
# rather than silently harmonized in either direction.

IL_ALIAS = {"Dlin-MC3-DMA": "DLin-MC3-DMA", "SM-102": "SM-102"}

# (file, dataset tag, {column: cell type})
DATASETS = [
    ("0_Norm_RAMOS_Single_Objective_UniformRandomSampling_Ctrl_MC3_DSPC.csv",
     "URS", {"nLnLE_RAMOS": "Ramos"}),
    ("0_Norm_RAMOS_THP1_Single_Dual_Objective_4ITER_MC3_DSPC.csv",
     "E0", {"nLnLE_RAMOS": "Ramos", "nLnLE_THP1": "THP-1"}),
    ("1_Norm_RAMOS_HEPG2_Single_Dual_Objective_4ITER_SM102_18PG.csv",
     "E1", {"nLnLE_RAMOS": "Ramos", "nLnLE_HEPG2": "HepG2"}),
    ("2_Norm_THP1_HEPG2_Dual_Objective_5ITER_SM102_DSPC.csv",
     "E2", {"nLnLE_THP1": "THP-1", "nLnLE_HEPG2": "HepG2"}),
    ("3_Norm_INS_Single_Objective_3ITER_SM102_DOTAP.csv",
     "E3", {"nLnLE_INS": "INS-1E"}),
    ("4_Norm_RAMOS_C2C12_DC24_3T3_Quad_Objective_3ITER_SM102_18PG.csv",
     "E4", {"nLnLE_RAMOS": "Ramos", "nLnLE_DC": "DC2.4",
            "nLnLE_3T3": "NIH-3T3", "nLnLE_C2C12": "C2C12"}),
    ("5_Norm_RAMOS_HEPG2_Dual_Objective_4ITER_SM102_DSPC.csv",
     "E5", {"nLnLE_HEPG2": "HepG2", "nLnLE_RAMOS": "Ramos"}),
]
EID_CELL = {"Ramos": "RAMOS", "THP-1": "THP1", "HepG2": "HEPG2",
            "INS-1E": "INS1E", "DC2.4": "DC24", "NIH-3T3": "3T3", "C2C12": "C2C12"}

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


def pct_columns(hdr):
    """(IL, HL, Chol, PEG) molar-percent column names, which differ per file."""
    pct = [c for c in hdr if c.endswith("PcntMolar")]
    chol = next(c for c in pct if c.startswith("Cholesterol"))
    peg = next(c for c in pct if "PEG" in c)
    rest = [c for c in pct if c not in (chol, peg)]
    il = next(c for c in rest if c.split("_")[0] in ("MC3", "SM102", "SM-102"))
    hl = next(c for c in rest if c != il)
    return il, hl, chol, peg


def main():
    corpus = list(csv.DictReader(open(CORPUS)))
    hdr = list(corpus[0].keys())
    # The quality/provenance columns, read from their single source rather than
    # restated here. Needed explicitly only on a build that runs before the
    # corpus has them; once it does, the line above already picks them up.
    for c in (x["name"] for x in json.load(open(os.path.join(
            HERE, "..", "..", "provenance", "schema_quality.json")))["columns"]):
        if c not in hdr:
            hdr.append(c)
    smiles = {}
    for r in corpus:
        for pre in ("IL", "HL", "CHL", "PEG"):
            n = (r[f"{pre}_name"] or "").strip()
            if n and r[f"{pre}_SMILES"]:
                smiles.setdefault((pre, n), r[f"{pre}_SMILES"])
    # MC3 carries two different SMILES in the corpus (AC_2025 stores a dioxolane
    # that is not MC3); take the one on the 2,194-row majority explicitly.
    smiles[("IL", "DLin-MC3-DMA")] = (
        "CCCCC/C=C\\C/C=C\\CCCCCCCCC(CCCCCCCC/C=C\\C/C=C\\CCCCC)OC(=O)CCCN(C)C")
    # SOPG is hand-written from the chemical name, so canonicalise it before use;
    # the corpus stores canonical SMILES and a non-canonical string would make the
    # same molecule look like two.
    smiles[("HL", "18PG")] = Chem.MolToSmiles(Chem.MolFromSmiles(SOPG))

    index = 1 + max(int(r["Index"]) for r in corpus if r["Index"].isdigit())
    rows, blocks, sel_check = [], [], []

    for fname, tag, outputs in DATASETS:
        src = list(csv.DictReader(open(os.path.join(SRC, fname), encoding="utf-8-sig")))
        head = list(src[0])
        ilc, hlc, cholc, pegc = pct_columns(head)
        selcol = next((c for c in head if c.endswith("Selectivity")), None)

        # ONE z-score per FILE, pooled over this file's cell-type columns, so a
        # within-file cross-cell-type difference survives (F-003). mu and sd are
        # computed once here and used for every column in the file.
        eid = f"WT_2026_{tag}"
        cells = []
        pooled_raw = [float(r[col]) for col in outputs for r in src]
        mu, sd = st.mean(pooled_raw), st.pstdev(pooled_raw)
        zcol = {col: [(float(r[col]) - mu) / sd for r in src] for col in outputs}

        for col, cell in outputs.items():
            vals = [float(r[col]) for r in src]
            zs = zcol[col]
            dose = DOSE_FIG6.get(cell, DOSE_UG[cell]) if tag == "E2" else DOSE_UG[cell]

            for r, v, z in zip(src, vals, zs):
                il_name = IL_ALIAS[r["Ionizable_lipid"]]
                hl_name = r["Helper_lipid"]
                fid = r.get("FALCON_ID") or f"{tag}_{r['Formula_label']}"
                row = dict.fromkeys(hdr, "")
                row.update({
                    "Index": index, "LNP_ID": f"LNP_{index:07d}",
                    "Experiment_ID": eid, "Formulation_ID": fid,
                    "IL_name": il_name, "IL_SMILES": smiles[("IL", il_name)],
                    "IL_molratio": float(r[ilc]),
                    "IL_to_nucleicacid_chargeratio": float(r["NP_ratio"]),
                    "HL_name": hl_name, "HL_SMILES": smiles[("HL", hl_name)],
                    "HL_molratio": float(r[hlc]),
                    "CHL_name": "Cholesterol", "CHL_SMILES": smiles[("CHL", "Cholesterol")],
                    "CHL_molratio": float(r[cholc]),
                    "PEG_name": "DMG-PEG2000", "PEG_SMILES": smiles[("PEG", "DMG-PEG2000")],
                    "PEG_molratio": float(r[pegc]),
                    "Aqueous_buffer": "Unknown", "Mixing_method": "handmixed",
                    "Model": "in_vitro", "Model_type": cell, "Model_target": "in_vitro",
                    "Route_of_administration": "in_vitro",
                    "Cargo": "mRNA", "Cargo_type": "FLuc",
                    "Dose_ug_nucleicacid": dose,
                    "Experiment_method": "luminescence_normalized",
                    "Experiment_batching": "individual",
                    "Experiment_value": z,
                    "Publication_link": LINK, "Publication_DOI": DOI,
                    "Publication_PMID": PMID,
                    # the unit this z was computed over, stated on the row so a
                    # consumer never has to infer it (F-001/F-016)
                    "value_normalization_unit": "source_file",
                })
                row.update(descriptors(row["IL_SMILES"]))
                rows.append(row)
                index += 1

            cells.append({
                "cell": cell, "column": col, "n": len(vals), "dose_ug": dose,
                "raw_LnnLE": {"min": round(min(vals), 4), "max": round(max(vals), 4),
                              "mean": round(st.mean(vals), 4),
                              "pop_sd": round(st.pstdev(vals), 4)},
                # NOT 0 and NOT 1, by design: these are the per-cell-type
                # statistics INSIDE the file's single z-score, and the fact that
                # they differ from each other is the between-assay offset that
                # per-column standardizing used to delete.
                "z_within_file": {"mean": round(st.mean(zs), 6),
                                  "pop_sd": round(st.pstdev(zs), 6)},
            })

        blocks.append({
            "experiment_id": eid, "source_file": fname,
            "n": sum(c["n"] for c in cells), "cells": cells,
            "pooled_raw_LnnLE": {"mean": round(mu, 4), "pop_sd": round(sd, 4)},
            "z": {"mean": round(st.mean([z for c in outputs for z in zcol[c]]), 12),
                  "pop_sd": round(st.pstdev([z for c in outputs for z in zcol[c]]), 12)},
        })
        # the file's own z IS standardized -- that one is not a tautology,
        # it is the definition of the normalization unit
        assert abs(blocks[-1]["z"]["mean"]) < 1e-9, (eid, blocks[-1]["z"])
        assert abs(blocks[-1]["z"]["pop_sd"] - 1) < 1e-9, (eid, blocks[-1]["z"])

        # Selectivity is a function of columns already curated as rows, not an
        # independent measurement. Rather than assume which form, solve for it:
        # every ordered pair difference, plus "on-target minus the mean of the
        # rest" for the quad-objective file. The winning form must reproduce the
        # printed column to ~1e-9, or the assumption that it is derived fails.
        if selcol and len(outputs) >= 2:
            cols = list(outputs)
            cands = {f"{a} - {b}": (lambda r, a=a, b=b: float(r[a]) - float(r[b]))
                     for a in cols for b in cols if a != b}
            if len(cols) > 2:
                for a in cols:
                    rest = [c for c in cols if c != a]
                    cands[f"{a} - mean(rest)"] = (
                        lambda r, a=a, rest=rest:
                        float(r[a]) - sum(float(r[c]) for c in rest) / len(rest))
            scored = {k: max(abs(f(r) - float(r[selcol])) for r in src)
                      for k, f in cands.items()}
            best = min(scored, key=scored.get)
            entry = {"file": fname, "form": best,
                     "max_abs_residual": round(scored[best], 12),
                     "derived": scored[best] < 1e-6}

            # THE CHECK THAT F-003 EXISTS FOR. Under one z-score per file,
            # z_a - z_b = (raw_a - raw_b)/sd = Selectivity/sd, a single constant
            # per file. So the stored values must reproduce the authors' printed
            # column exactly, up to that constant. This is a real test of the
            # normalization unit: it FAILS if the z is ever taken per cell type
            # again, because then each column carries its own mu and sd and no
            # single constant can recover the difference.
            #
            # It is deliberately checked against the PRINTED column rather than
            # against our own difference of raw columns -- the authors' number is
            # independent of how we recover it. Only the pairwise forms are
            # testable this way; "a - mean(rest)" (file 4) is linear in the z as
            # well, so it is checked by the same identity.
            if " - mean(rest)" in best:
                a = best.split(" - ")[0]
                rest = [c for c in cols if c != a]
                zsel = [zcol[a][i] - sum(zcol[c][i] for c in rest) / len(rest)
                        for i in range(len(src))]
            else:
                a, b = best.split(" - ")
                zsel = [zcol[a][i] - zcol[b][i] for i in range(len(src))]
            printed = [float(r[selcol]) for r in src]
            resid = max(abs(z * sd - p) for z, p in zip(zsel, printed))
            entry["stored_z_recovers_printed_selectivity"] = {
                "scale": round(sd, 6), "max_abs_residual": resid}
            assert resid < 1e-6, (fname, "stored z does not reproduce the "
                                  "authors' Selectivity column", resid)
            sel_check.append(entry)

    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=hdr)
        w.writeheader()
        w.writerows(rows)

    # structures for the C10 drawn-structure check
    with open(os.path.join(HERE, "03_enumerated.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["lipid_name", "IL_SMILES"])
        for n in ("DLin-MC3-DMA", "SM-102"):
            w.writerow([n, smiles[("IL", n)]])

    molsum = [round(r["IL_molratio"] + r["HL_molratio"] + r["CHL_molratio"]
                    + r["PEG_molratio"], 4) for r in rows]
    qc = {
        "paper": "Toh, Cheng et al., bioRxiv 2026 (FALCON)", "doi": DOI,
        "source": "authors' own CSVs; github.com/MaoResearchGroup/falcon-lnp-optimization",
        "formulations": sum(len(list(csv.DictReader(
            open(os.path.join(SRC, f), encoding="utf-8-sig")))) for f, _, _ in DATASETS),
        "rows": len(rows), "columns": len(hdr),
        "experiments": len(blocks),
        "unique_IL": sorted({r["IL_name"] for r in rows}),
        "unique_HL": sorted({r["HL_name"] for r in rows}),
        "cell_types": sorted({r["Model_type"] for r in rows}),
        "molratio_sum": {"min": min(molsum), "max": max(molsum)},
        "normalization_unit": "source_file (one z per file, pooled over that "
                              "file's cell-type columns) -- see F-003",
        "z_block_max_abs_mean": max(abs(b["z"]["mean"]) for b in blocks),
        "z_block_sd_range": [min(b["z"]["pop_sd"] for b in blocks),
                             max(b["z"]["pop_sd"] for b in blocks)],
        # the between-assay offset per-column z-scoring used to delete, now
        # preserved and reported: the spread of per-cell-type means INSIDE each
        # file's single z. All zero would mean the bug is back.
        "within_file_cell_offsets": {
            b["experiment_id"]: {c["cell"]: c["z_within_file"]["mean"]
                                 for c in b["cells"]} for b in blocks},
        "selectivity_is_derived": sel_check,
        "blocks": blocks,
    }
    json.dump(qc, open(os.path.join(HERE, "02_WT_2026_qc.json"), "w"), indent=1)

    print(f"{os.path.basename(OUT)}: {len(rows)} rows x {len(hdr)} cols, "
          f"{len(blocks)} experiments")
    print(f"  ILs {qc['unique_IL']}  HLs {qc['unique_HL']}")
    print(f"  cells {qc['cell_types']}")
    print(f"  molratio sums {qc['molratio_sum']}")
    print(f"  z blocks: |mean| <= {qc['z_block_max_abs_mean']:.2e}, "
          f"sd in {qc['z_block_sd_range']}")
    for s in sel_check:
        print(f"  selectivity {s['file'][:26]:28} = {s['form']:30} "
              f"residual {s['max_abs_residual']:.1e}  derived={s['derived']}")


if __name__ == "__main__":
    main()
