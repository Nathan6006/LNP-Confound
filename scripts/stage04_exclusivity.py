"""
Stage 4 — Source-exclusivity and excipient invariance.

Explains Stage 3's separability structurally: how much of it is "a lipid only
ever appears in one source" (trivial lookup) vs. genuine chemical clustering,
and whether the surrounding formulation (helper lipid / cholesterol / PEG) is
too invariant across the corpus to explain cross-lab disagreement itself.

Outputs:
  results/stage4_exclusivity.csv
  figures/stage4_*.png (+ .svg)
  reports/04_exclusivity.md
  registry keys: exclusivity.{vitro,vivo}.*, formulation.{vitro,vivo}.*
"""
import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency
from sklearn.metrics import normalized_mutual_info_score

from common import load_curated, registry_set, RESULTS, REPORTS, FIGURES, smiles_to_mol
from rdkit import Chem


def canonical_or_none(smi):
    if pd.isna(smi):
        return None
    mol = smiles_to_mol(smi)
    if mol is None:
        return None
    return Chem.MolToSmiles(mol, canonical=True)


def exclusivity_table(df, id_col, source_col="Experiment_ID"):
    """For each distinct value of id_col, how many distinct sources does it appear in?"""
    g = df.dropna(subset=[id_col]).groupby(id_col)[source_col].nunique()
    return g


def cramers_v(confusion_matrix):
    chi2, p, dof, exp = chi2_contingency(confusion_matrix)
    n = confusion_matrix.sum().sum()
    r, k = confusion_matrix.shape
    denom = n * (min(r, k) - 1)
    if denom <= 0:
        return np.nan, p
    return np.sqrt(chi2 / denom), p


def concentration_stats(series):
    vc = series.value_counts(dropna=True)
    n = vc.sum()
    top1 = vc.iloc[0] / n if len(vc) else np.nan
    top3 = vc.iloc[:3].sum() / n if len(vc) else np.nan
    herfindahl = ((vc / n) ** 2).sum() if len(vc) else np.nan
    return {"n_distinct": len(vc), "n_total": int(n), "top1_frac": top1, "top3_frac": top3,
            "herfindahl_index": herfindahl, "top_values": vc.head(5).to_dict()}


def run_for_dataset(df_all, dataset_name, report, csv_rows):
    df = df_all[df_all["_source_dataset"] == dataset_name].copy()
    report.append(f"## `{dataset_name}`\n")

    # ---- 1. Lipid exclusivity ----
    lipid_src_counts = exclusivity_table(df, "canonical_IL_SMILES")
    frac_single = (lipid_src_counts == 1).mean()
    report.append(f"### Ionizable-lipid exclusivity\n")
    report.append(f"- Distinct lipids: {len(lipid_src_counts):,}")
    report.append(f"- Appear in exactly 1 source: {(lipid_src_counts==1).sum():,} ({frac_single:.1%})")
    dist = lipid_src_counts.value_counts().sort_index()
    report.append(f"- Distribution of #sources-per-lipid: {dist.to_dict()}\n")
    registry_set(f"exclusivity.{dataset_name}.lipid.frac_single_source", frac_single,
                 n=len(lipid_src_counts), stage="stage4")
    csv_rows.append({"dataset": dataset_name, "level": "ionizable_lipid", "n_distinct": len(lipid_src_counts),
                      "frac_single_source": frac_single})

    # ---- component-level exclusivity (DB-provided head/linker/tail decomposition) ----
    report.append("### Component-level exclusivity (DB-provided head/linker/tail1/tail2 decomposition)\n")
    report.append("Using the pre-existing `IL_head/linker/tail*_SMILES` decomposition columns already in "
                   "LNPDB, rather than a fresh BRICS/fragmentation pass — these are the DB curators' "
                   "chemically-informed split, which is more reliable than an automatic fragmentation "
                   "run cold on already-decomposed lipids. Coverage (non-missing) varies by component and "
                   "is reported for each.\n")
    for comp, col in [("head", "IL_head_SMILES"), ("linker", "IL_linker_SMILES"),
                       ("tail1", "IL_tail1_SMILES"), ("tail2", "IL_tail2_SMILES")]:
        if col not in df.columns:
            report.append(f"- **{comp}**: column `{col}` not available in curated table — skipped.\n")
            continue
        present = df[col].notna()
        coverage = present.mean()
        canon_col = df.loc[present, col].map(canonical_or_none)
        sub = df.loc[present].copy()
        sub["_canon_comp"] = canon_col
        sub = sub.dropna(subset=["_canon_comp"])
        comp_counts = exclusivity_table(sub, "_canon_comp")
        frac_single_comp = (comp_counts == 1).mean() if len(comp_counts) else np.nan
        report.append(f"- **{comp}** (`{col}`): coverage={coverage:.1%} ({present.sum():,}/{len(df):,} rows), "
                      f"{len(comp_counts):,} distinct components, {frac_single_comp:.1%} single-source")
        registry_set(f"exclusivity.{dataset_name}.component_{comp}.frac_single_source", frac_single_comp,
                     n=len(comp_counts), stage="stage4", confidence="OK" if coverage > 0.3 else "LOW-CONFIDENCE",
                     note=f"coverage={coverage:.1%}")
        csv_rows.append({"dataset": dataset_name, "level": f"component_{comp}", "n_distinct": len(comp_counts),
                          "frac_single_source": frac_single_comp, "coverage": coverage})
    report.append("\n")

    # ---- 2. Excipient invariance ----
    report.append("### Excipient invariance\n")
    for comp, col in [("helper_lipid", "HL_name"), ("peg_lipid", "PEG_name"), ("cholesterol", "CHL_name")]:
        stats = concentration_stats(df[col])
        report.append(f"- **{comp}** (`{col}`): {stats['n_distinct']} distinct identities across {stats['n_total']:,} rows, "
                      f"top-1 = {stats['top1_frac']:.1%}, top-3 = {stats['top3_frac']:.1%}, "
                      f"Herfindahl = {stats['herfindahl_index']:.3f} — top values: {stats['top_values']}")
        registry_set(f"formulation.{dataset_name}.{comp}.top1_frac", stats["top1_frac"],
                     n=stats["n_total"], stage="stage4")
        registry_set(f"formulation.{dataset_name}.{comp}.herfindahl", stats["herfindahl_index"],
                     n=stats["n_total"], stage="stage4")
    report.append(
        "\nInterpretation: if the excipient identities are highly concentrated (few distinct HL/PEG/CHL "
        "choices shared across nearly all sources — high top-1/top-3 fraction, high Herfindahl), the "
        "formulation backbone is nearly invariant across the corpus and **cannot itself explain** "
        "cross-lab disagreement in Stage 5 (there is little formulation variance to correlate with the "
        "disagreement). If instead formulation is a strong source signature, it becomes an alternative "
        "explanation for Stage 3 separability that competes with lipid-chemistry exclusivity.\n"
    )

    # ---- 3. Formulation-source association ----
    report.append("### Formulation-source association (Cramer's V, NMI)\n")
    fmt = lambda x: "NA" if pd.isna(x) else round(float(x), 3)
    df["_formulation_key"] = (
        df["HL_name"].astype(str) + "|" + df["CHL_name"].astype(str) + "|" + df["PEG_name"].astype(str) + "|" +
        df["HL_molratio"].round(1).astype(str) + "|" + df["CHL_molratio"].round(1).astype(str) + "|" +
        df["PEG_molratio"].round(1).astype(str)
    )
    ct = pd.crosstab(df["Experiment_ID"], df["_formulation_key"])
    v, p = cramers_v(ct)
    nmi = normalized_mutual_info_score(df["Experiment_ID"], df["_formulation_key"])
    report.append(f"- Formulation key = HL/CHL/PEG identity + molar ratios (rounded to 0.1)")
    report.append(f"- Distinct formulation keys: {df['_formulation_key'].nunique():,}")
    report.append(f"- Cramer's V(source, formulation) = **{fmt(v)}** (chi2 p={p:.2e})")
    report.append(f"- Normalized mutual information(source, formulation) = **{fmt(nmi)}**\n")
    registry_set(f"formulation.{dataset_name}.cramers_v_source", v, n=len(df), stage="stage4",
                 extra={"chi2_p": float(p)})
    registry_set(f"formulation.{dataset_name}.nmi_source", nmi, n=len(df), stage="stage4")

    # ---- figure: lipid exclusivity histogram ----
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    dist.plot(kind="bar", ax=axes[0])
    axes[0].set_xlabel("# sources a lipid appears in"); axes[0].set_ylabel("# lipids")
    axes[0].set_title(f"Lipid exclusivity — {dataset_name}")
    for comp, col in [("HL", "HL_name"), ("PEG", "PEG_name"), ("CHL", "CHL_name")]:
        vc = df[col].value_counts(normalize=True).head(8)
        axes[1].plot(range(len(vc)), vc.values, marker="o", label=comp)
    axes[1].set_xlabel("rank"); axes[1].set_ylabel("fraction of rows"); axes[1].legend()
    axes[1].set_title(f"Excipient concentration — {dataset_name}")
    fig.tight_layout()
    fig.savefig(FIGURES / f"stage4_exclusivity_{dataset_name}.png", dpi=300)
    fig.savefig(FIGURES / f"stage4_exclusivity_{dataset_name}.svg")
    plt.close(fig)
    report.append(f"Figure: `figures/stage4_exclusivity_{dataset_name}.png` (+ .svg)\n")


def main():
    df = load_curated()
    report = ["# Stage 4 — Source-Exclusivity and Excipient Invariance\n"]
    csv_rows = []

    for ds in ["LNPDB_vitro", "LNPDB_vivo"]:
        run_for_dataset(df, ds, report, csv_rows)

    pd.DataFrame(csv_rows).to_csv(RESULTS / "stage4_exclusivity.csv", index=False)
    (REPORTS / "04_exclusivity.md").write_text("\n".join(str(x) for x in report))

    print("--- STAGE 4 SUMMARY ---")
    for r in csv_rows:
        if r["level"] == "ionizable_lipid":
            print(f"{r['dataset']}: {r['frac_single_source']:.1%} of {r['n_distinct']:,} lipids are single-source")
    print("Wrote results/stage4_exclusivity.csv, reports/04_exclusivity.md, figures/stage4_exclusivity_*.png")


if __name__ == "__main__":
    main()
