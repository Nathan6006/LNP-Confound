"""
Stage 2 — Corpus overview (Table 1).

Descriptive backbone: per-source and aggregate counts for the frozen curated
corpus from Stage 1. Reads only results/curated.parquet and
results/environments.csv; computes nothing that touches a downstream metric.

Outputs:
  results/table1_corpus.csv
  reports/02_corpus.md
  registry keys: corpus.n_sources_vitro, corpus.n_sources_vivo,
                 corpus.n_measurements_vitro, corpus.n_measurements_vivo,
                 corpus.n_comparable_environments_vitro/vivo, corpus.n_unique_lipids_*
"""
from common import load_curated, load_curated_molea, registry_set, RESULTS, REPORTS

import pandas as pd


def main():
    df = load_curated()
    molea = load_curated_molea()
    env = pd.read_csv(RESULTS / "environments.csv")

    report = ["# Stage 2 — Corpus Overview (Table 1)\n"]

    # ---- Per-source table ----
    rows = []
    for (ds, exp_id), g in df.groupby(["_source_dataset", "Experiment_ID"]):
        env_key_mode = g["environment_key"].mode().iloc[0]
        rows.append({
            "source_dataset": ds,
            "Experiment_ID": exp_id,
            "n_measurements": len(g),
            "n_unique_lipids": g["canonical_IL_SMILES"].nunique(),
            "readout_methods": ";".join(sorted(g["Experiment_method"].unique())),
            "cell_or_animal_models": ";".join(sorted(g["Model_type"].astype(str).unique())),
            "routes": ";".join(sorted(g["Route_of_administration"].astype(str).unique())),
            "cargo": ";".join(sorted(g["Cargo"].astype(str).unique())),
            "n_environments": g["environment_key"].nunique(),
            "n_comparable_environments": g.loc[g["environment_comparable"], "environment_key"].nunique(),
        })
    table1 = pd.DataFrame(rows).sort_values(["source_dataset", "n_measurements"], ascending=[True, False])
    table1.to_csv(RESULTS / "table1_corpus.csv", index=False)

    report.append("## Per-source summary\n")
    report.append(table1.to_markdown(index=False))
    report.append("\n")

    # ---- Aggregate counts ----
    report.append("## Aggregate counts\n")
    for ds in ["LNPDB_vitro", "LNPDB_vivo"]:
        sub = df[df["_source_dataset"] == ds]
        n_sources = sub["Experiment_ID"].nunique()
        n_meas = len(sub)
        n_lipids = sub["canonical_IL_SMILES"].nunique()
        env_sub = env[env["_source_dataset"] == ds]
        n_envs = len(env_sub)
        n_comparable = int(env_sub["comparable"].sum())

        report.append(f"**{ds}**")
        report.append(f"- Sources (raw): {n_sources}")
        report.append(f"- Measurements: {n_meas:,}")
        report.append(f"- Unique ionizable lipids: {n_lipids:,}")
        report.append(f"- Environments: {n_envs} total, **{n_comparable} comparable** (>= threshold in Stage 1)\n")

        registry_set(f"corpus.n_sources_{ds}", n_sources, stage="stage2")
        registry_set(f"corpus.n_measurements_{ds}", n_meas, stage="stage2")
        registry_set(f"corpus.n_unique_lipids_{ds}", n_lipids, stage="stage2")
        registry_set(f"corpus.n_environments_{ds}", n_envs, stage="stage2")
        registry_set(f"corpus.n_comparable_environments_{ds}", n_comparable, stage="stage2")

    report.append(f"**molea** (kept separate; single source, two endpoints)")
    report.append(f"- Rows: {len(molea):,}")
    report.append(f"- Unique canonical SMILES: {molea['canonical_smiles'].nunique():,}")
    report.append(f"- Hepatocyte values present: {molea['mTP_(Log2)_Hepatocyte'].notna().sum():,}")
    report.append(f"- Chondrocyte values present: {molea['mTP_(Log2)_Chondrocyte'].notna().sum():,}\n")
    registry_set("corpus.n_rows_molea", len(molea), stage="stage2")

    total_meas = len(df)
    total_vitro = (df["_source_dataset"] == "LNPDB_vitro").sum()
    total_vivo = (df["_source_dataset"] == "LNPDB_vivo").sum()
    report.append("## Totals\n")
    report.append(f"- Total LNPDB measurements: {total_meas:,} ({total_vitro:,} in-vitro / {total_vivo:,} in-vivo)")
    report.append(f"- Total LNPDB sources: {df['Experiment_ID'].nunique()}")
    registry_set("corpus.total_measurements", int(total_meas), stage="stage2")

    (REPORTS / "02_corpus.md").write_text("\n".join(str(x) for x in report))

    print("--- STAGE 2 SUMMARY ---")
    print(f"vitro: {total_vitro:,} measurements, {df[df._source_dataset=='LNPDB_vitro'].Experiment_ID.nunique()} sources, "
          f"{df[df._source_dataset=='LNPDB_vitro'].canonical_IL_SMILES.nunique():,} unique lipids")
    print(f"vivo:  {total_vivo:,} measurements, {df[df._source_dataset=='LNPDB_vivo'].Experiment_ID.nunique()} sources, "
          f"{df[df._source_dataset=='LNPDB_vivo'].canonical_IL_SMILES.nunique():,} unique lipids")
    print("Wrote results/table1_corpus.csv, reports/02_corpus.md")


if __name__ == "__main__":
    main()
