"""
Stage 7 — What stays learnable, and how fast it decays.

Two distinct questions from Stage 6's grid, both needing dedicated analysis:

1. Within-experiment learnability: train and test a model *entirely inside one
   source* (scaffold-disjoint split), i.e. does structure-activity hold at all
   for that source's own local chemistry, independent of cross-lab pooling.
   This differs from Stage 6's "within_experiment_spearman" cells, which score
   a pooled-trained model's performance measured inside each source — here the
   model itself never sees any other source's data.

2. Decay curve: reuse Stage 6's leave-one-source-out folds and ask whether the
   held-out source's Spearman score falls off smoothly or as a cliff as a
   function of its chemical distance (median Tanimoto NN distance) to the
   pooled training support.

Outputs:
  results/stage7_learnable.csv (per-source learnability + decay points)
  figures/stage7_decay_{vitro,vivo}.png (+ .svg)
  reports/07_learnable.md
  registry keys: learnable.{vitro,vivo}.size_weighted_spearman, decay.{vitro,vivo}.trend_slope
"""
import warnings

import numpy as np
import pandas as pd
from rdkit import DataStructs
from scipy.stats import spearmanr, linregress
from sklearn.model_selection import GroupKFold
from xgboost import XGBRegressor

warnings.filterwarnings("ignore")

from common import (
    GLOBAL_SEED, MORGAN_RADIUS, MORGAN_NBITS, load_curated, registry_set,
    RESULTS, REPORTS, FIGURES, batch_morgan_fingerprints, morgan_fingerprint_bitvect,
    butina_clusters,
)
from stage03_separability import rdkit_descriptor_matrix
from stage05_concordance import add_within_source_endpoint_rank, lipid_source_table

MIN_LIPIDS_FOR_WITHIN_SOURCE_CV = 30
BUTINA_CUTOFF = 0.4
XGB_PARAMS = dict(n_estimators=200, max_depth=4, learning_rate=0.1, random_state=GLOBAL_SEED,
                   n_jobs=4, verbosity=0)


def within_source_learnability(X, y, cluster_ids, seed=GLOBAL_SEED):
    """Scaffold-disjoint (Butina cluster) k-fold entirely within one source's own
    data. Returns (spearman, n_test, reason) — reason is None on success, else why
    it could not be computed (too few lipids vs. too chemically homogeneous to
    scaffold-split, which are different findings and should not look the same)."""
    if len(X) < MIN_LIPIDS_FOR_WITHIN_SOURCE_CV:
        return np.nan, 0, f"fewer than {MIN_LIPIDS_FOR_WITHIN_SOURCE_CV} lipids"
    n_unique_clusters = len(set(cluster_ids))
    n_splits = min(5, n_unique_clusters)
    if n_splits < 2:
        return np.nan, 0, f"only {n_unique_clusters} Butina cluster(s) — too chemically homogeneous to scaffold-split"
    gkf = GroupKFold(n_splits=n_splits)
    all_true, all_pred = [], []
    for train_idx, test_idx in gkf.split(X, y, groups=cluster_ids):
        if len(train_idx) == 0 or len(test_idx) == 0:
            continue
        model = XGBRegressor(**XGB_PARAMS)
        model.fit(X[train_idx], y[train_idx])
        pred = model.predict(X[test_idx])
        all_true.extend(y[test_idx].tolist())
        all_pred.extend(pred.tolist())
    if len(set(all_true)) < 2:
        return np.nan, len(all_true), "degenerate test set (constant true values)"
    return spearmanr(all_true, all_pred).statistic, len(all_true), None


def median_nn_tanimoto_distance(held_out_bitvects, train_bitvects):
    dists = []
    for bv in held_out_bitvects:
        sims = DataStructs.BulkTanimotoSimilarity(bv, train_bitvects)
        dists.append(1 - max(sims))
    return float(np.median(dists))


def run_for_dataset(df_all, dataset_name, report, csv_rows):
    report.append(f"## `{dataset_name}`\n")
    df = df_all[df_all["_source_dataset"] == dataset_name].copy()
    df = add_within_source_endpoint_rank(df)
    lipid_src = lipid_source_table(df).dropna(subset=["rank"]).reset_index(drop=True)

    smiles = lipid_src["canonical_IL_SMILES"].tolist()
    sources = lipid_src["Experiment_ID"].values
    y_all = lipid_src["rank"].values

    fps, valid_fp = batch_morgan_fingerprints(smiles, MORGAN_RADIUS, MORGAN_NBITS, cache_name=f"stage6_{dataset_name}")
    X_desc, valid_desc, _ = rdkit_descriptor_matrix(smiles)
    valid = valid_fp & valid_desc
    X_all = np.concatenate([fps[valid], X_desc[valid]], axis=1)
    y_all = y_all[valid]
    sources = sources[valid]
    smiles_valid = np.array([s for s, v in zip(smiles, valid) if v])
    bitvects_all = [morgan_fingerprint_bitvect(s) for s in smiles_valid]

    # ---- 1. Within-experiment learnability, per source ----
    report.append(f"### Within-experiment learnability (train+test inside one source)\n")
    report.append(f"Minimum {MIN_LIPIDS_FOR_WITHIN_SOURCE_CV} lipids/source to attempt a local "
                  f"scaffold-disjoint split; smaller sources reported as not-computable, not skipped silently.\n")
    per_source_rows = []
    for src in sorted(set(sources)):
        mask = sources == src
        n_src = int(mask.sum())
        if n_src < MIN_LIPIDS_FOR_WITHIN_SOURCE_CV:
            per_source_rows.append({"source": src, "n_lipids": n_src, "spearman": np.nan, "n_test": 0,
                                     "computable": False,
                                     "reason": f"fewer than {MIN_LIPIDS_FOR_WITHIN_SOURCE_CV} lipids"})
            continue
        X_src, y_src = X_all[mask], y_all[mask]
        bv_src = [bitvects_all[i] for i in np.where(mask)[0]]
        cluster_ids_src, _ = butina_clusters(bv_src, cutoff=BUTINA_CUTOFF)
        rho, n_test, reason = within_source_learnability(X_src, y_src, cluster_ids_src)
        per_source_rows.append({"source": src, "n_lipids": n_src, "spearman": rho, "n_test": n_test,
                                 "computable": not np.isnan(rho), "reason": reason})

    per_source_df = pd.DataFrame(per_source_rows)
    computable = per_source_df[per_source_df["computable"]]
    if len(computable):
        size_weighted = float(np.average(computable["spearman"], weights=computable["n_lipids"]))
    else:
        size_weighted = np.nan
    report.append(f"- Sources with enough lipids to test: {len(computable)} / {len(per_source_df)}")
    report.append(f"- Size-weighted within-experiment scaffold-disjoint Spearman: **{size_weighted:.3f}**\n")
    report.append("| source | n_lipids | spearman | computable | reason if not |")
    report.append("|---|---|---|---|---|")
    for _, r in per_source_df.sort_values("n_lipids", ascending=False).iterrows():
        sp = f"{r.spearman:.3f}" if r.computable else "NA"
        reason = "" if r.computable else r.get("reason", "")
        report.append(f"| {r.source} | {r.n_lipids} | {sp} | {r.computable} | {reason} |")
    report.append("\n")
    registry_set(f"learnable.{dataset_name}.size_weighted_spearman", size_weighted,
                 n=int(computable["n_lipids"].sum()) if len(computable) else 0, stage="stage7",
                 extra={"n_sources_computable": int(len(computable))})
    per_source_df["dataset"] = dataset_name
    per_source_df["metric"] = "within_experiment_learnability"

    # ---- 2. Decay curve (reuse Stage 6 LOSO folds' chemical distance vs. performance) ----
    report.append("### Decay curve: LOSO fold performance vs. chemical distance to training support\n")
    decay_rows = []
    for src in sorted(set(sources)):
        test_mask = sources == src
        train_mask = ~test_mask
        if test_mask.sum() == 0 or train_mask.sum() == 0:
            continue
        bv_test = [bitvects_all[i] for i in np.where(test_mask)[0]]
        bv_train = [bitvects_all[i] for i in np.where(train_mask)[0]]
        dist = median_nn_tanimoto_distance(bv_test, bv_train)

        model = XGBRegressor(**XGB_PARAMS)
        model.fit(X_all[train_mask], y_all[train_mask])
        pred = model.predict(X_all[test_mask])
        y_test = y_all[test_mask]
        rho = spearmanr(y_test, pred).statistic if len(set(y_test)) > 1 and len(set(pred)) > 1 else np.nan
        decay_rows.append({"dataset": dataset_name, "metric": "decay", "source": src,
                            "n_lipids": int(test_mask.sum()), "median_nn_tanimoto_distance": dist,
                            "spearman": rho})

    decay_df = pd.DataFrame(decay_rows)
    valid_decay = decay_df.dropna(subset=["spearman"])
    if len(valid_decay) >= 3:
        lr = linregress(valid_decay["median_nn_tanimoto_distance"], valid_decay["spearman"])
        report.append(f"- Linear trend: Spearman_fold = {lr.intercept:.3f} + {lr.slope:.3f} * distance "
                      f"(r={lr.rvalue:.3f}, p={lr.pvalue:.3f}, n={len(valid_decay)} sources)")
        report.append(f"- {'Negative slope: performance degrades as held-out chemistry gets further from training support (smooth decay), consistent with the slope direction expected if the model is extrapolating, not just failing categorically.' if lr.slope < 0 else 'Slope is non-negative — performance does not show the expected smooth decay with distance in this pool; a cliff/threshold pattern or no relationship is at least as consistent with the data.'}\n")
        registry_set(f"decay.{dataset_name}.trend_slope", lr.slope, n=len(valid_decay), stage="stage7",
                     extra={"r": lr.rvalue, "p": lr.pvalue, "intercept": lr.intercept})
    else:
        report.append("- Too few LOSO folds with a computable Spearman to fit a trend.\n")

    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    ax.scatter(decay_df["median_nn_tanimoto_distance"], decay_df["spearman"], s=40)
    for _, r in decay_df.iterrows():
        ax.annotate(r["source"], (r["median_nn_tanimoto_distance"], r["spearman"]), fontsize=6, alpha=0.7)
    if len(valid_decay) >= 3:
        xs = np.linspace(valid_decay["median_nn_tanimoto_distance"].min(), valid_decay["median_nn_tanimoto_distance"].max(), 50)
        ax.plot(xs, lr.intercept + lr.slope * xs, "--", color="gray")
    ax.axhline(0, color="black", lw=0.5)
    ax.set_xlabel("median NN Tanimoto distance (held-out source -> training support)")
    ax.set_ylabel("LOSO fold Spearman")
    ax.set_title(f"Decay curve — {dataset_name}")
    fig.tight_layout()
    fig.savefig(FIGURES / f"stage7_decay_{dataset_name}.png", dpi=300)
    fig.savefig(FIGURES / f"stage7_decay_{dataset_name}.svg")
    plt.close(fig)
    report.append(f"Figure: `figures/stage7_decay_{dataset_name}.png` (+ .svg)\n")

    csv_rows.append(per_source_df)
    csv_rows.append(decay_df)


def main():
    df = load_curated()
    report = ["# Stage 7 — Within-Experiment Learnability and Transfer Decay\n"]
    csv_rows = []

    for ds in ["LNPDB_vitro", "LNPDB_vivo"]:
        run_for_dataset(df, ds, report, csv_rows)

    pd.concat(csv_rows, ignore_index=True).to_csv(RESULTS / "stage7_learnable.csv", index=False)
    (REPORTS / "07_learnable.md").write_text("\n".join(str(x) for x in report))

    print("--- STAGE 7 SUMMARY ---")
    from common import registry_get_all
    reg = registry_get_all()
    for ds in ["LNPDB_vitro", "LNPDB_vivo"]:
        k = f"learnable.{ds}.size_weighted_spearman"
        if k in reg:
            print(f"{ds}: within-experiment learnability Spearman = {reg[k]['value']:.3f} (n_sources={reg[k]['note'] or reg[k].get('n_sources_computable')})")
        kd = f"decay.{ds}.trend_slope"
        if kd in reg:
            print(f"{ds}: decay trend slope = {reg[kd]['value']:.3f}")
    print("Wrote results/stage7_learnable.csv, reports/07_learnable.md, figures/stage7_decay_*.png")


if __name__ == "__main__":
    main()
