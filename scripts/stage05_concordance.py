"""
Stage 5 — Do identical lipids rank the same across labs? (CENTERPIECE)

Holds chemistry fixed (identical canonical_IL_SMILES) and asks whether two
independent sources agree on its relative delivery rank. This is the paper's
load-bearing result, so every step is deliberately conservative: normalized
ranks only (never raw readouts across sources), explicit per-pair n next to
every pooled number, a computed (not asserted) chemical-bias caveat, and a
real assay-noise ceiling from replicate measurements rather than an assumed
one.

Outputs:
  results/stage5_concordance.csv (per-pair + pooled + ceiling + anchoring)
  figures/stage5_concordance_{vitro,vivo}.png (+ .svg)
  reports/05_concordance.md
  registry keys: concordance.{vitro,vivo}.pooled_spearman[...], .ceiling.*, .anchoring.*, .molea.*
"""
from itertools import combinations

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, kendalltau, pearsonr

from common import (
    GLOBAL_SEED, load_curated, load_curated_molea, normalized_rank, bootstrap_ci,
    bemis_murcko_scaffold, registry_set, RESULTS, REPORTS, FIGURES,
)

REFERENCE_LIPID_KEYWORDS = [
    "mc3", "sm-102", "sm102", "alc-0315", "alc0315", "c12-200", "c12200",
    "ckk-e12", "306oi10", "306o-i10", "98n12", "dlin-dma", "dodap",
]
MIN_REPLICATE_GROUPS_FOR_CEILING = 20


def add_within_source_endpoint_rank(df):
    df = df.copy()
    df["_rank"] = df.groupby(["Experiment_ID", "Experiment_method"])["Experiment_value"].transform(
        lambda s: normalized_rank(s).values
    )
    return df


def lipid_source_table(df):
    """One row per (Experiment_ID, canonical_IL_SMILES): mean within-source-endpoint
    normalized rank, and how many raw rows were averaged."""
    g = df.groupby(["Experiment_ID", "canonical_IL_SMILES"]).agg(
        rank=("_rank", "mean"), n_rows=("_rank", "size"),
    ).reset_index()
    return g


def build_pairwise_points(lipid_src, min_sources=2):
    """For every canonical_IL_SMILES shared by >=2 sources, emit one
    (source_A, source_B, rank_A, rank_B, smiles) row per pair of sources sharing it."""
    rows = []
    for smi, g in lipid_src.groupby("canonical_IL_SMILES"):
        if g["Experiment_ID"].nunique() < min_sources:
            continue
        for (idx_a, ra), (idx_b, rb) in combinations(g.set_index("Experiment_ID")["rank"].items(), 2):
            rows.append({"smiles": smi, "source_a": idx_a, "source_b": idx_b, "rank_a": ra, "rank_b": rb})
    return pd.DataFrame(rows)


def pooled_concordance(points_df, seed=GLOBAL_SEED):
    if len(points_df) < 2:
        return {"spearman": np.nan, "kendall": np.nan, "n": len(points_df), "ci_lo": np.nan, "ci_hi": np.nan}
    a, b = points_df["rank_a"].values, points_df["rank_b"].values
    rho = spearmanr(a, b).statistic
    tau = kendalltau(a, b).statistic
    ci = bootstrap_ci((a, b), lambda x, y: spearmanr(x, y).statistic, seed=seed)
    return {"spearman": rho, "kendall": tau, "n": len(points_df), "ci_lo": ci[0], "ci_hi": ci[1]}


def weighted_pair_average(pair_table):
    """Sample-size-weighted average of per-pair Spearman rho (robustness cross-check
    against the naive pooled-points correlation)."""
    valid = pair_table.dropna(subset=["spearman"])
    if valid["n_shared"].sum() == 0:
        return np.nan
    return float((valid["spearman"] * valid["n_shared"]).sum() / valid["n_shared"].sum())


def replicate_noise_ceiling(df, seed=GLOBAL_SEED):
    """Assay-noise ceiling from condition-complete replicate rows (same source, lipid,
    formulation identity+ratios, dose, model, readout method, cargo — differing only in
    which physical well/LNP_ID was measured). For groups with >=2 such rows and genuine
    value variation, randomly pick 2 as (repA, repB), look up their already-computed
    within-source-endpoint ranks, and pool the (rankA, rankB) pairs."""
    key_cols = [
        "Experiment_ID", "canonical_IL_SMILES", "canonical_HL_SMILES", "HL_molratio",
        "canonical_CHL_SMILES", "CHL_molratio", "canonical_PEG_SMILES", "PEG_molratio",
        "Dose_ug_nucleicacid", "Model_type", "Experiment_method", "Cargo_type",
    ]
    d = df.copy()
    d["_key"] = d[key_cols].astype(str).agg("|".join, axis=1)
    rng = np.random.default_rng(seed)
    pairs = []
    per_source_counts = {}
    for key, g in d.groupby("_key"):
        if len(g) < 2 or g["Experiment_value"].nunique() < 2:
            continue
        idx = rng.choice(len(g), size=2, replace=False)
        r = g["_rank"].values[idx]
        pairs.append((r[0], r[1]))
        src = g["Experiment_ID"].iloc[0]
        per_source_counts[src] = per_source_counts.get(src, 0) + 1
    if not pairs:
        return None
    arr = np.array(pairs)
    result = pooled_concordance(pd.DataFrame({"rank_a": arr[:, 0], "rank_b": arr[:, 1]}), seed=seed)
    result["per_source_group_counts"] = per_source_counts
    return result


def reference_anchoring_test(df, lipid_src, pair_table, seed=GLOBAL_SEED):
    """Pick the canonical SMILES shared across the most sources; re-express each
    eligible source's readout relative to its own measurement of that anchor lipid;
    recompute pooled concordance on the SAME restricted source-pair set (both members
    must contain the anchor) and the same shared-lipid set, excluding the anchor
    itself, for a fair before/after comparison."""
    counts = lipid_src.groupby("canonical_IL_SMILES")["Experiment_ID"].nunique().sort_values(ascending=False)
    if len(counts) == 0 or counts.iloc[0] < 2:
        return None
    anchor_smiles = counts.index[0]
    eligible_sources = set(lipid_src.loc[lipid_src["canonical_IL_SMILES"] == anchor_smiles, "Experiment_ID"])
    if len(eligible_sources) < 2:
        return None

    d = df[df["Experiment_ID"].isin(eligible_sources)].copy()
    anchor_value = d[d["canonical_IL_SMILES"] == anchor_smiles].groupby("Experiment_ID")["Experiment_value"].mean()
    d["_anchored_value"] = d["Experiment_value"] - d["Experiment_ID"].map(anchor_value)
    d["_anchored_rank"] = d.groupby(["Experiment_ID", "Experiment_method"])["_anchored_value"].transform(
        lambda s: normalized_rank(s).values
    )
    anchored_lipid_src = d.groupby(["Experiment_ID", "canonical_IL_SMILES"]).agg(
        rank=("_anchored_rank", "mean")
    ).reset_index()
    anchored_lipid_src = anchored_lipid_src[anchored_lipid_src["canonical_IL_SMILES"] != anchor_smiles]

    anchored_points = build_pairwise_points(anchored_lipid_src)
    anchored_result = pooled_concordance(anchored_points, seed=seed)

    # baseline restricted to the SAME source-pair set and SAME excluded-anchor lipid set
    restricted_pairs = pair_table[
        pair_table["source_a"].isin(eligible_sources) & pair_table["source_b"].isin(eligible_sources) &
        (pair_table["smiles"] != anchor_smiles)
    ]
    baseline_result = pooled_concordance(restricted_pairs, seed=seed)

    return {
        "anchor_smiles": anchor_smiles, "eligible_sources": sorted(eligible_sources),
        "baseline": baseline_result, "anchored": anchored_result,
    }


def molea_cross_endpoint(molea):
    hep = molea["mTP_(Log2)_Hepatocyte"]
    chon = molea["mTP_(Log2)_Chondrocyte"]
    both = hep.notna() & chon.notna()
    n = int(both.sum())
    if n < 2:
        return None
    hep_rank = normalized_rank(hep[both])
    chon_rank = normalized_rank(chon[both])
    rho = spearmanr(hep_rank, chon_rank).statistic
    pear = pearsonr(hep[both], chon[both]).statistic
    ci = bootstrap_ci((hep_rank.values, chon_rank.values), lambda x, y: spearmanr(x, y).statistic, seed=GLOBAL_SEED)
    return {"n": n, "spearman_rank": rho, "pearson_raw": pear, "ci_lo": ci[0], "ci_hi": ci[1]}


def run_for_dataset(df_all, dataset_name, report, csv_rows):
    report.append(f"## `{dataset_name}`\n")
    df = df_all[df_all["_source_dataset"] == dataset_name].copy()
    df = add_within_source_endpoint_rank(df)

    lipid_src = lipid_source_table(df)
    shared_counts = lipid_src.groupby("canonical_IL_SMILES")["Experiment_ID"].nunique()
    n_shared_lipids = int((shared_counts >= 2).sum())
    report.append(f"- Distinct lipids: {lipid_src['canonical_IL_SMILES'].nunique():,}; "
                  f"shared across >=2 sources: **{n_shared_lipids}**")

    pair_table_points = build_pairwise_points(lipid_src)
    if len(pair_table_points) == 0:
        report.append("- **No shared lipids across sources — cross-lab concordance cannot be computed "
                      "for this pool.** (LOW-CONFIDENCE / not computable, not fabricated.)\n")
        registry_set(f"concordance.{dataset_name}.pooled_spearman", None, n=0, stage="stage5",
                     confidence="LOW-CONFIDENCE", note="no shared lipids across sources")
        return

    n_source_pairs = pair_table_points[["source_a", "source_b"]].drop_duplicates().shape[0]
    report.append(f"- Source-pairs sharing >=1 lipid: **{n_source_pairs}**")
    report.append(f"- Total pooled (source-pair, lipid) points: **{len(pair_table_points)}**\n")

    # ---- per-pair table ----
    per_pair = []
    for (sa, sb), g in pair_table_points.groupby(["source_a", "source_b"]):
        n = len(g)
        rho = spearmanr(g["rank_a"], g["rank_b"]).statistic if n >= 2 else np.nan
        tau = kendalltau(g["rank_a"], g["rank_b"]).statistic if n >= 2 else np.nan
        per_pair.append({"source_a": sa, "source_b": sb, "n_shared": n, "spearman": rho, "kendall": tau})
    pair_table = pd.DataFrame(per_pair).sort_values("n_shared", ascending=False)
    # re-attach per-point smiles/source for the anchoring baseline restriction
    pair_table_points_named = pair_table_points

    report.append("**Per-pair concordance** (top by n; full table in results/stage5_concordance.csv):\n")
    report.append("| source_a | source_b | n_shared | spearman | kendall |")
    report.append("|---|---|---|---|---|")
    for _, r in pair_table.head(12).iterrows():
        report.append(f"| {r.source_a} | {r.source_b} | {r.n_shared} | "
                      f"{r.spearman:.3f} | {r.kendall:.3f} |" if not pd.isna(r.spearman)
                      else f"| {r.source_a} | {r.source_b} | {r.n_shared} | n<2 | n<2 |")
    report.append("\n")

    # ---- pooled concordance (primary: pooled points; robustness: n-weighted pair average) ----
    pooled = pooled_concordance(pair_table_points)
    weighted_avg = weighted_pair_average(pair_table)
    confidence = "OK" if pooled["n"] >= 30 else "LOW-CONFIDENCE"
    report.append(f"**Pooled concordance (primary — all cross-source points concatenated)**: "
                  f"Spearman = **{pooled['spearman']:.3f}** [{pooled['ci_lo']:.3f}, {pooled['ci_hi']:.3f}] "
                  f"(95% bootstrap CI), Kendall = {pooled['kendall']:.3f}, n = {pooled['n']} — "
                  f"**{confidence}**")
    report.append(f"**Robustness cross-check (n-weighted average of per-pair Spearman)**: {weighted_avg:.3f}\n")

    registry_set(f"concordance.{dataset_name}.pooled_spearman", pooled["spearman"], n=pooled["n"],
                 ci=[pooled["ci_lo"], pooled["ci_hi"]], confidence=confidence, stage="stage5",
                 extra={"kendall": pooled["kendall"], "n_weighted_pair_avg_spearman": weighted_avg,
                        "n_source_pairs": int(n_source_pairs), "n_shared_lipids": n_shared_lipids})

    # ---- chemical-bias caveat ----
    shared_smiles = shared_counts[shared_counts >= 2].index.tolist()
    shared_names = df[df["canonical_IL_SMILES"].isin(shared_smiles)].groupby("canonical_IL_SMILES")["IL_name"].first()
    is_reference = shared_names.astype(str).str.lower().apply(lambda n: any(k in n for k in REFERENCE_LIPID_KEYWORDS))
    n_reference = int(is_reference.sum())
    scaffolds = {smi: bemis_murcko_scaffold(smi) for smi in shared_smiles}
    n_acyclic = sum(1 for v in scaffolds.values() if v == "")
    # an empty Murcko scaffold means "no ring system" (fully acyclic lipid, common for
    # lipidoid libraries) — NOT a shared scaffold. Treating "" as a real scaffold value
    # would falsely lump every acyclic lipid into one giant fake "family", so acyclic
    # lipids are excluded from the congeneric count and reported separately instead.
    ring_scaffolds = {smi: s for smi, s in scaffolds.items() if s}
    scaffold_counts = pd.Series(list(ring_scaffolds.values())).value_counts()
    n_congeneric = int(sum(1 for smi, s in ring_scaffolds.items() if scaffold_counts.get(s, 0) > 1))
    n_ring = len(ring_scaffolds)
    ring_frac_str = f"{n_congeneric/n_ring:.1%}" if n_ring else "NA"
    report.append(f"**Chemical-bias caveat**: of {len(shared_smiles)} shared lipids, "
                  f"{n_reference} ({n_reference/len(shared_smiles):.1%}) match a known reference-lipid name "
                  f"(MC3/SM-102/ALC-0315/etc.); {n_acyclic} ({n_acyclic/len(shared_smiles):.1%}) are fully "
                  f"acyclic (no Bemis-Murcko scaffold to compare — reported separately, not counted as "
                  f"'congeneric' since an empty scaffold is not a real shared feature); of the "
                  f"{n_ring} ring-containing shared lipids, {n_congeneric} ({ring_frac_str} of "
                  f"those, or {n_congeneric/len(shared_smiles):.1%} of all shared lipids) share a real "
                  f"Bemis-Murcko scaffold with another shared lipid (congeneric family member). "
                  f"Top single pair contributes {pair_table['n_shared'].max()} / {pooled['n']} "
                  f"({pair_table['n_shared'].max()/pooled['n']:.1%}) of pooled points.\n")
    registry_set(f"concordance.{dataset_name}.shared_set_frac_reference_lipid",
                 n_reference / len(shared_smiles), n=len(shared_smiles), stage="stage5")

    # ---- assay-noise ceiling ----
    ceiling = replicate_noise_ceiling(df)
    if ceiling is None or ceiling["n"] < MIN_REPLICATE_GROUPS_FOR_CEILING:
        n_ceiling = ceiling["n"] if ceiling else 0
        report.append(f"**Assay-noise ceiling**: only {n_ceiling} condition-complete replicate groups found "
                      f"(threshold {MIN_REPLICATE_GROUPS_FOR_CEILING}) — **LOW-CONFIDENCE / not reliably "
                      f"estimable** for this pool.\n")
        registry_set(f"concordance.{dataset_name}.ceiling_spearman", ceiling["spearman"] if ceiling else None,
                     n=n_ceiling, stage="stage5", confidence="LOW-CONFIDENCE")
    else:
        dominant_source = max(ceiling["per_source_group_counts"], key=ceiling["per_source_group_counts"].get)
        dominant_frac = ceiling["per_source_group_counts"][dominant_source] / sum(ceiling["per_source_group_counts"].values())
        report.append(f"**Assay-noise ceiling** (within-source replicate split-half): Spearman = "
                      f"**{ceiling['spearman']:.3f}** [{ceiling['ci_lo']:.3f}, {ceiling['ci_hi']:.3f}], "
                      f"n = {ceiling['n']} replicate-group pairs across {len(ceiling['per_source_group_counts'])} "
                      f"sources (dominant source `{dominant_source}` contributes {dominant_frac:.1%} of groups — "
                      f"noted, not hidden). **This is the ceiling the cross-lab concordance above should be "
                      f"compared against.**\n")
        registry_set(f"concordance.{dataset_name}.ceiling_spearman", ceiling["spearman"], n=ceiling["n"],
                     ci=[ceiling["ci_lo"], ceiling["ci_hi"]], stage="stage5",
                     extra={"dominant_source": dominant_source, "dominant_source_frac": dominant_frac})

    # ---- reference-lipid anchoring test ----
    anchoring = reference_anchoring_test(df, lipid_src, pair_table_points_named)
    if anchoring is None:
        report.append("**Reference-lipid anchoring test**: no lipid is shared across enough sources to run "
                      "this test for this pool.\n")
    else:
        b, a = anchoring["baseline"], anchoring["anchored"]
        report.append(f"**Reference-lipid anchoring test**: anchor = `{anchoring['anchor_smiles'][:40]}...` "
                      f"shared by {len(anchoring['eligible_sources'])} sources "
                      f"({', '.join(anchoring['eligible_sources'])}). Baseline pooled Spearman (same restricted "
                      f"source-pairs, anchor excluded) = {b['spearman']:.3f} (n={b['n']}) -> after anchoring = "
                      f"**{a['spearman']:.3f}** (n={a['n']}). "
                      f"{'Anchoring improved concordance.' if (not pd.isna(a['spearman']) and not pd.isna(b['spearman']) and a['spearman'] > b['spearman']) else ('Anchoring did not improve concordance.' if not pd.isna(a['spearman']) and not pd.isna(b['spearman']) else 'Comparison inconclusive (insufficient n).')}\n")
        registry_set(f"concordance.{dataset_name}.anchoring.baseline_spearman", b["spearman"], n=b["n"], stage="stage5")
        registry_set(f"concordance.{dataset_name}.anchoring.anchored_spearman", a["spearman"], n=a["n"], stage="stage5",
                     extra={"anchor_smiles": anchoring["anchor_smiles"], "eligible_sources": anchoring["eligible_sources"]})

    # ---- figure ----
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(5.5, 5.5))
    ax.scatter(pair_table_points["rank_a"], pair_table_points["rank_b"], s=14, alpha=0.5)
    ax.plot([0, 1], [0, 1], "--", color="gray", lw=1)
    ax.set_xlabel("within-source normalized rank (source A)")
    ax.set_ylabel("within-source normalized rank (source B)")
    ax.set_title(f"Cross-lab shared-lipid concordance — {dataset_name}\n"
                 f"pooled Spearman={pooled['spearman']:.2f} [{pooled['ci_lo']:.2f},{pooled['ci_hi']:.2f}], n={pooled['n']}")
    fig.tight_layout()
    fig.savefig(FIGURES / f"stage5_concordance_{dataset_name}.png", dpi=300)
    fig.savefig(FIGURES / f"stage5_concordance_{dataset_name}.svg")
    plt.close(fig)
    report.append(f"Figure: `figures/stage5_concordance_{dataset_name}.png` (+ .svg)\n")

    pair_table["dataset"] = dataset_name
    csv_rows.append(pair_table)


def main():
    df = load_curated()
    molea = load_curated_molea()
    report = ["# Stage 5 — Cross-Lab Shared-Lipid Concordance (Centerpiece)\n",
              "Within-source normalized rank is computed within (Experiment_ID, Experiment_method) groups; "
              "raw readouts are never compared across sources directly. See module docstring for the pooling "
              "method (all cross-source-pair points concatenated) and the assay-noise-ceiling construction.\n"]
    csv_rows = []

    for ds in ["LNPDB_vitro", "LNPDB_vivo"]:
        run_for_dataset(df, ds, report, csv_rows)

    report.append("## MOLEA cross-endpoint check (independent corroborating data point)\n")
    molea_result = molea_cross_endpoint(molea)
    if molea_result is None:
        report.append("Not computable (insufficient paired rows).\n")
    else:
        report.append(f"Within MOLEA (single lab, identical chemistry, only the assayed cell type changes): "
                      f"Hepatocyte-vs-Chondrocyte rank Spearman = **{molea_result['spearman_rank']:.3f}** "
                      f"[{molea_result['ci_lo']:.3f}, {molea_result['ci_hi']:.3f}], raw-value Pearson = "
                      f"{molea_result['pearson_raw']:.3f}, n = {molea_result['n']}. This is a same-lab, "
                      f"single-variable-changed data point — a middle rung between within-experiment "
                      f"learnable agreement and cross-lab shared-lipid agreement above.\n")
        registry_set("concordance.molea.hepatocyte_vs_chondrocyte_spearman", molea_result["spearman_rank"],
                     n=molea_result["n"], ci=[molea_result["ci_lo"], molea_result["ci_hi"]], stage="stage5",
                     extra={"pearson_raw": molea_result["pearson_raw"]})

    if csv_rows:
        pd.concat(csv_rows, ignore_index=True).to_csv(RESULTS / "stage5_concordance.csv", index=False)
    (REPORTS / "05_concordance.md").write_text("\n".join(str(x) for x in report))

    print("--- STAGE 5 SUMMARY ---")
    from common import registry_get_all
    reg = registry_get_all()
    for ds in ["LNPDB_vitro", "LNPDB_vivo"]:
        k = f"concordance.{ds}.pooled_spearman"
        if k in reg:
            e = reg[k]
            print(f"{ds}: pooled Spearman={e['value']}, n={e['n']}, CI={e['ci']}, confidence={e['confidence']}")
        kc = f"concordance.{ds}.ceiling_spearman"
        if kc in reg:
            e = reg[kc]
            print(f"{ds}: ceiling Spearman={e['value']}, n={e['n']}, confidence={e['confidence']}")
    if molea_result:
        print(f"molea: hep-vs-chon Spearman={molea_result['spearman_rank']:.3f}, n={molea_result['n']}")
    print("Wrote results/stage5_concordance.csv, reports/05_concordance.md, figures/stage5_concordance_*.png")


if __name__ == "__main__":
    main()
