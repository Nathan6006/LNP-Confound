"""
Stage 1 — Curation, canonicalization, and the frozen analysis set.

All exclusion/inclusion decisions here are frozen a priori and are metric-blind:
none of them were chosen by looking at how they affect any downstream
correlation, accuracy, or separability number. They rest only on (a) whether a
value is present/parseable and (b) whether a row is a byte-for-byte duplicate
of another.

Outputs:
  results/curated.parquet        — frozen LNPDB (vitro+vivo) analysis table
  results/curated_molea.parquet  — frozen MOLEA analysis table (kept separate)
  results/exclusions.json        — every exclusion rule, reason, and row count
  results/environments.csv       — endpoint-comparable environment groups
  reports/01_curation.md         — human-readable summary of every decision
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import rdkit
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
REPORTS = ROOT / "reports"
RESULTS = ROOT / "results"
REPORTS.mkdir(exist_ok=True)
RESULTS.mkdir(exist_ok=True)

MIN_ROWS_PER_ENVIRONMENT = 20  # a priori threshold for a "comparable environment" to be individually analyzable

with open(RESULTS / "schema.json") as f:
    SCHEMA = json.load(f)


def canonicalize_series(series, label, exclusions_log):
    """Canonicalize a pandas Series of SMILES strings with RDKit.
    Returns (canonical_series, n_failed). Non-null inputs that fail to parse
    become NaN in the output; failures are counted but the caller decides
    whether to exclude the row (only done for the primary IL SMILES column)."""
    def _canon(smi):
        if pd.isna(smi):
            return np.nan
        mol = Chem.MolFromSmiles(str(smi))
        if mol is None:
            return None  # parse failure, distinct from missing input
        return Chem.MolToSmiles(mol, canonical=True)

    canon = series.map(_canon)
    n_input_present = series.notna().sum()
    n_failed = int((canon.isna() & series.notna()).sum())
    exclusions_log.append(
        {
            "column": label,
            "n_input_present": int(n_input_present),
            "n_parse_failed": n_failed,
            "pct_parse_failed": round(100 * n_failed / n_input_present, 3) if n_input_present else 0.0,
        }
    )
    return canon, n_failed


def main():
    report = []
    exclusions = []  # list of {rule, dataset, reason, rows_affected}
    smiles_parse_log = []

    report.append("# Stage 1 — Curation, Canonicalization, and the Frozen Analysis Set\n")
    report.append(f"RDKit version: **{rdkit.__version__}**\n")
    report.append(
        "All exclusion criteria below were fixed before any downstream metric (separability, "
        "concordance, protocol comparison) was computed. Nothing here was chosen to make a "
        "result look better or worse.\n"
    )

    # ================= Load =================
    vitro = pd.read_csv(DATA / "LNPDB_vitro.csv", low_memory=False)
    vivo = pd.read_csv(DATA / "LNPDB_vivo.csv", low_memory=False)
    molea = pd.read_csv(DATA / "molea.csv", low_memory=False)

    n_vitro_raw, n_vivo_raw, n_molea_raw = len(vitro), len(vivo), len(molea)
    report.append(f"## 1. Raw row counts\n")
    report.append(f"- `LNPDB_vitro.csv`: {n_vitro_raw:,} rows")
    report.append(f"- `LNPDB_vivo.csv`: {n_vivo_raw:,} rows")
    report.append(f"- `molea.csv`: {n_molea_raw:,} rows\n")

    vitro["_source_dataset"] = "LNPDB_vitro"
    vivo["_source_dataset"] = "LNPDB_vivo"
    lnpdb = pd.concat([vitro, vivo], ignore_index=True)
    lnpdb["_orig_row_index"] = lnpdb.groupby("_source_dataset").cumcount()

    # ================= SMILES canonicalization =================
    report.append("## 2. SMILES canonicalization\n")
    report.append("Canonicalized with RDKit `Chem.MolFromSmiles` → `Chem.MolToSmiles(canonical=True)`. "
                   "Parse failures on the primary ionizable-lipid column are excluded (Rule E1 below); "
                   "failures on formulation SMILES columns are logged but do not drop rows (formulation "
                   "identity is not the primary key).\n")

    lnpdb["canonical_IL_SMILES"], n_il_failed_lnpdb = canonicalize_series(lnpdb["IL_SMILES"], "LNPDB.IL_SMILES", smiles_parse_log)
    molea["canonical_smiles"], n_il_failed_molea = canonicalize_series(molea["smiles"], "molea.smiles", smiles_parse_log)

    for col, label in [
        ("IL_protonated_SMILES", "LNPDB.IL_protonated_SMILES"),
        ("HL_SMILES", "LNPDB.HL_SMILES"),
        ("CHL_SMILES", "LNPDB.CHL_SMILES"),
        ("PEG_SMILES", "LNPDB.PEG_SMILES"),
        ("fifthcomponent_SMILES", "LNPDB.fifthcomponent_SMILES"),
    ]:
        canon, _ = canonicalize_series(lnpdb[col], label, smiles_parse_log)
        lnpdb["canonical_" + col] = canon

    report.append("| column | n present | n parse-failed | % failed |")
    report.append("|---|---|---|---|")
    for row in smiles_parse_log:
        report.append(f"| {row['column']} | {row['n_input_present']:,} | {row['n_parse_failed']:,} | {row['pct_parse_failed']}% |")
    report.append("\n")

    # ================= Exclusion Rule E1: IL SMILES parse failure =================
    mask_il_fail_lnpdb = lnpdb["IL_SMILES"].notna() & lnpdb["canonical_IL_SMILES"].isna()
    n_e1_lnpdb = int(mask_il_fail_lnpdb.sum())
    exclusions.append({
        "rule": "E1_il_smiles_parse_failure", "dataset": "LNPDB",
        "reason": "IL_SMILES did not parse as a valid molecule under RDKit",
        "rows_affected": n_e1_lnpdb,
    })
    lnpdb = lnpdb[~mask_il_fail_lnpdb].copy()

    mask_il_fail_molea = molea["smiles"].notna() & molea["canonical_smiles"].isna()
    n_e1_molea = int(mask_il_fail_molea.sum())
    exclusions.append({
        "rule": "E1_il_smiles_parse_failure", "dataset": "molea",
        "reason": "smiles did not parse as a valid molecule under RDKit",
        "rows_affected": n_e1_molea,
    })
    molea = molea[~mask_il_fail_molea].copy()

    # ================= Exclusion Rule E2: missing delivery readout =================
    mask_missing_value = lnpdb["Experiment_value"].isna()
    n_e2 = int(mask_missing_value.sum())
    exclusions.append({
        "rule": "E2_missing_delivery_readout", "dataset": "LNPDB",
        "reason": "Experiment_value (delivery readout) is missing",
        "rows_affected": n_e2,
    })
    lnpdb = lnpdb[~mask_missing_value].copy()

    # molea: both endpoint columns must be checked (a row missing BOTH is unusable;
    # a row missing one endpoint is still usable for the other)
    mask_molea_both_missing = molea["mTP_(Log2)_Hepatocyte"].isna() & molea["mTP_(Log2)_Chondrocyte"].isna()
    n_e2_molea = int(mask_molea_both_missing.sum())
    exclusions.append({
        "rule": "E2_missing_delivery_readout", "dataset": "molea",
        "reason": "Both mTP_(Log2)_Hepatocyte and mTP_(Log2)_Chondrocyte are missing",
        "rows_affected": n_e2_molea,
    })
    molea = molea[~mask_molea_both_missing].copy()

    # ================= Exclusion Rule E3: exact duplicate rows =================
    id_like_cols = {"Index", "LNP_ID", "Formulation_ID", "_source_dataset", "_orig_row_index"}
    content_cols = [c for c in lnpdb.columns if c not in id_like_cols and not c.startswith("canonical_")]
    dup_mask = lnpdb.duplicated(subset=content_cols, keep="first")
    n_e3 = int(dup_mask.sum())
    exclusions.append({
        "rule": "E3_exact_duplicate_row", "dataset": "LNPDB",
        "reason": "Row is a byte-for-byte content duplicate of an earlier row (all fields incl. readout value identical)",
        "rows_affected": n_e3,
    })
    lnpdb = lnpdb[~dup_mask].copy()

    dup_mask_molea = molea.duplicated(subset=[c for c in molea.columns if c not in ("canonical_smiles",)], keep="first")
    n_e3_molea = int(dup_mask_molea.sum())
    exclusions.append({
        "rule": "E3_exact_duplicate_row", "dataset": "molea",
        "reason": "Row is a byte-for-byte content duplicate of an earlier row",
        "rows_affected": n_e3_molea,
    })
    molea = molea[~dup_mask_molea].copy()

    report.append("## 3. Frozen exclusion rules applied (in order)\n")
    report.append("| rule | dataset | reason | rows affected |")
    report.append("|---|---|---|---|")
    for e in exclusions:
        report.append(f"| {e['rule']} | {e['dataset']} | {e['reason']} | {e['rows_affected']:,} |")
    report.append("\n")

    n_vitro_final = int((lnpdb["_source_dataset"] == "LNPDB_vitro").sum())
    n_vivo_final = int((lnpdb["_source_dataset"] == "LNPDB_vivo").sum())
    n_molea_final = len(molea)

    report.append("## 4. Post-exclusion row counts\n")
    report.append(f"- LNPDB in-vitro: {n_vitro_raw:,} → **{n_vitro_final:,}** ({n_vitro_raw - n_vitro_final:,} excluded)")
    report.append(f"- LNPDB in-vivo: {n_vivo_raw:,} → **{n_vivo_final:,}** ({n_vivo_raw - n_vivo_final:,} excluded)")
    report.append(f"- molea: {n_molea_raw:,} → **{n_molea_final:,}** ({n_molea_raw - n_molea_final:,} excluded)\n")

    # ================= Source label & endpoint definition =================
    report.append("## 5. Source-label and endpoint definitions\n")
    report.append(
        "- **LNPDB**: source label = `Experiment_ID` (native, 1:1 with `Publication_PMID`). Endpoint = "
        "the row itself; no further split needed since `Experiment_ID` + `Experiment_method` + "
        "`Model_type` already uniquely identify an assay context.\n"
        "- **MOLEA**: no source/publication column exists in the raw file. **DECISION (frozen, matches "
        "CLAUDE.md §6, 2026-07-23/24):** treat MOLEA as **one source, two endpoints** — "
        "`source_label='MOLEA'`, `endpoint in {hepatocyte, chondrocyte}`. This is a single combinatorial "
        "screen (`Markush code`) from one lab/library measuring two cell types, not two independent "
        "publications, so splitting it into two 'sources' would misrepresent its provenance. "
        "`mTP_Difference` is a derived selectivity contrast between the two endpoints, not an "
        "independently poolable readout, and is carried through unchanged for reference but not treated "
        "as a third endpoint.\n"
    )

    # ================= Ionizable lipid vs full formulation =================
    report.append("## 6. Ionizable lipid vs. full formulation\n")
    report.append(
        "Downstream chemical analyses (source separability, shared-lipid concordance) use "
        "`canonical_IL_SMILES` only. Formulation identity (`HL_*`, `CHL_*`, `PEG_*`, "
        "`fifthcomponent_*`, each canonicalized where present) is retained as separate columns for the "
        "Stage 4 excipient-invariance analysis and is never mixed into the lipid identity key.\n"
    )

    # ================= Endpoint-comparable environments =================
    report.append("## 7. Endpoint-comparable environments\n")
    report.append(
        f"An **environment** = a unique combination of (in-vitro/in-vivo, `Experiment_method` "
        f"[readout modality], `Model_type` [cell line or animal strain — the biological system], "
        f"`Route_of_administration`). An environment is **comparable/analyzable** if it has "
        f"**≥ {MIN_ROWS_PER_ENVIRONMENT} rows** (frozen threshold, chosen for basic statistical "
        f"stability, not to hit a target count). This is stricter than counting raw sources: two "
        f"sources using the identical readout/cell-line/route combination are the same environment, "
        f"and a source split across multiple readouts contributes multiple environments.\n"
    )

    lnpdb["environment_key"] = (
        lnpdb["_source_dataset"].astype(str) + " | " +
        lnpdb["Experiment_method"].astype(str) + " | " +
        lnpdb["Model_type"].astype(str) + " | " +
        lnpdb["Route_of_administration"].astype(str)
    )
    env_counts = lnpdb.groupby(["_source_dataset", "environment_key"]).size().reset_index(name="n_rows")
    env_counts["comparable"] = env_counts["n_rows"] >= MIN_ROWS_PER_ENVIRONMENT
    env_counts.to_csv(RESULTS / "environments.csv", index=False)

    for ds in ["LNPDB_vitro", "LNPDB_vivo"]:
        sub = env_counts[env_counts["_source_dataset"] == ds]
        n_raw_sources = lnpdb.loc[lnpdb["_source_dataset"] == ds, "Experiment_ID"].nunique()
        n_envs = len(sub)
        n_comparable = int(sub["comparable"].sum())
        report.append(f"**{ds}**: {n_raw_sources} raw sources → {n_envs} distinct environments, "
                       f"of which **{n_comparable} are comparable** (≥{MIN_ROWS_PER_ENVIRONMENT} rows).\n")
        report.append("| environment | n_rows | comparable |")
        report.append("|---|---|---|")
        for _, r in sub.sort_values("n_rows", ascending=False).iterrows():
            report.append(f"| {r['environment_key']} | {r['n_rows']:,} | {'yes' if r['comparable'] else 'no'} |")
        report.append("\n")

    lnpdb = lnpdb.merge(env_counts[["environment_key", "comparable"]].drop_duplicates("environment_key"),
                         on="environment_key", how="left")
    lnpdb = lnpdb.rename(columns={"comparable": "environment_comparable"})

    # ================= molea: canonical component structure sanity =================
    report.append("## 8. MOLEA structure notes\n")
    n_markush_dup = int(molea["Markush code"].duplicated().sum())
    n_smiles_dup = int(molea["canonical_smiles"].duplicated().sum())
    report.append(f"- Distinct `Markush code` values: {molea['Markush code'].nunique():,} (of {len(molea):,} rows; "
                   f"{n_markush_dup:,} duplicate codes)")
    report.append(f"- Distinct canonical SMILES: {molea['canonical_smiles'].nunique():,} "
                   f"({n_smiles_dup:,} rows share a canonical SMILES with another row)\n")

    # ================= Provenance & write outputs =================
    lnpdb["provenance"] = lnpdb["_source_dataset"] + ":row" + lnpdb["_orig_row_index"].astype(str) + \
        ":" + lnpdb["LNP_ID"].astype(str)
    molea["provenance"] = "molea:row" + molea.index.astype(str)

    keep_cols_lnpdb = [
        "provenance", "_source_dataset", "LNP_ID", "Formulation_ID",
        "Experiment_ID", "Publication_PMID", "Publication_link",
        "IL_name", "IL_SMILES", "canonical_IL_SMILES", "IL_protonated_SMILES", "canonical_IL_protonated_SMILES",
        "Experiment_value", "Experiment_method", "Experiment_batching",
        "Model", "Model_type", "Model_target", "Route_of_administration",
        "Cargo", "Cargo_type", "Dose_ug_nucleicacid",
        "HL_name", "HL_SMILES", "canonical_HL_SMILES", "HL_molratio",
        "CHL_name", "CHL_SMILES", "canonical_CHL_SMILES", "CHL_molratio",
        "PEG_name", "PEG_SMILES", "canonical_PEG_SMILES", "PEG_molratio",
        "fifthcomponent_name", "fifthcomponent_SMILES", "canonical_fifthcomponent_SMILES", "fifthcomponent_molratio",
        "IL_molratio", "environment_key", "environment_comparable",
        "IL_head_name", "IL_head_SMILES", "IL_linker_name", "IL_linker_SMILES",
        "IL_tail1_name", "IL_tail1_SMILES", "IL_tail2_name", "IL_tail2_SMILES",
    ]
    curated_lnpdb = lnpdb[keep_cols_lnpdb].reset_index(drop=True)
    curated_lnpdb.to_parquet(RESULTS / "curated.parquet", index=False)

    keep_cols_molea = [
        "provenance", "canonical_smiles", "smiles", "Markush code",
        "mTP_(Log2)_Hepatocyte", "mTP_(Log2)_Chondrocyte", "mTP_Difference",
    ]
    curated_molea = molea[keep_cols_molea].reset_index(drop=True)
    curated_molea.to_parquet(RESULTS / "curated_molea.parquet", index=False)

    with open(RESULTS / "exclusions.json", "w") as f:
        json.dump({
            "min_rows_per_environment": MIN_ROWS_PER_ENVIRONMENT,
            "rdkit_version": rdkit.__version__,
            "smiles_parse_log": smiles_parse_log,
            "exclusion_rules": exclusions,
            "raw_counts": {"LNPDB_vitro": n_vitro_raw, "LNPDB_vivo": n_vivo_raw, "molea": n_molea_raw},
            "final_counts": {"LNPDB_vitro": n_vitro_final, "LNPDB_vivo": n_vivo_final, "molea": n_molea_final},
        }, f, indent=2)

    report.append("## 9. Frozen outputs\n")
    report.append(f"- `results/curated.parquet`: {len(curated_lnpdb):,} rows (LNPDB vitro+vivo, post-exclusion)")
    report.append(f"- `results/curated_molea.parquet`: {len(curated_molea):,} rows")
    report.append(f"- `results/exclusions.json`, `results/environments.csv`\n")

    (REPORTS / "01_curation.md").write_text("\n".join(str(x) for x in report))

    print("\n--- STAGE 1 SUMMARY ---")
    print(f"LNPDB in-vitro: {n_vitro_raw:,} -> {n_vitro_final:,}")
    print(f"LNPDB in-vivo:  {n_vivo_raw:,} -> {n_vivo_final:,}")
    print(f"molea:          {n_molea_raw:,} -> {n_molea_final:,}")
    print(f"Comparable environments (vitro): {int(env_counts[(env_counts._source_dataset=='LNPDB_vitro') & (env_counts.comparable)].shape[0])}")
    print(f"Comparable environments (vivo):  {int(env_counts[(env_counts._source_dataset=='LNPDB_vivo') & (env_counts.comparable)].shape[0])}")
    print("Wrote reports/01_curation.md, results/curated.parquet, results/curated_molea.parquet, results/exclusions.json, results/environments.csv")


if __name__ == "__main__":
    main()
