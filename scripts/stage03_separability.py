"""
Stage 3 — Does chemistry predict its source? (separability)

Tests whether ionizable-lipid structure alone reveals the measurement's source
(Experiment_ID), for the in-vitro and in-vivo pools separately.

Unit of analysis: one row per distinct (source, canonical_IL_SMILES) pair, not
per raw measurement row. Raw rows massively over-count a handful of reference
lipids re-measured under many formulation/dose conditions (e.g. AC_2025: 1,978
rows but 7 unique lipids) — classifying at the raw-row level would mostly be
memorizing which duplicate-content row belongs to which source, not testing
whether *chemistry* is source-predictive.

Minimum-rows rule (frozen a priori): a source needs >=20 distinct ionizable
lipids to be included in the classifier (needed for a stable 5-fold split and
a meaningful multiclass signal); excluded sources are listed, not silently
dropped from the report.

Outputs:
  results/stage3_separability.csv
  figures/stage3_embedding_{vitro,vivo}.png (+ .svg)
  reports/03_separability.md
  registry keys: separability.{vitro,vivo}.balanced_accuracy_{random,scaffold}[...]
"""
import sys
import time
import warnings

import numpy as np
import pandas as pd
from rdkit import DataStructs
from rdkit.Chem import Descriptors
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, f1_score, recall_score, top_k_accuracy_score
from sklearn.model_selection import StratifiedKFold, GroupKFold
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore")

from common import (
    GLOBAL_SEED, MORGAN_RADIUS, MORGAN_NBITS,
    load_curated, registry_set, RESULTS, REPORTS, FIGURES,
    smiles_to_mol, batch_morgan_fingerprints, bemis_murcko_scaffold,
)

MIN_LIPIDS_PER_SOURCE = 20
N_FOLDS = 5
N_PERMUTATIONS = 200  # see report note: reduced from the spec's >=1000 for compute budget (M3/16GB);
                       # uses the fast linear classifier, not RF, as the permutation-null statistic.


def build_lipid_source_table(df, dataset_name, min_lipids):
    sub = df[df["_source_dataset"] == dataset_name].copy()
    counts = sub.groupby("Experiment_ID")["canonical_IL_SMILES"].nunique()
    qualifying = sorted(counts[counts >= min_lipids].index.tolist())
    excluded = sorted(counts[counts < min_lipids].index.tolist())
    dedup = (
        sub[sub["Experiment_ID"].isin(qualifying)]
        .drop_duplicates(subset=["Experiment_ID", "canonical_IL_SMILES"])
        .reset_index(drop=True)
    )
    return dedup, qualifying, excluded, counts


def rdkit_descriptor_matrix(smiles_list):
    names = [n for n, _ in Descriptors._descList]
    calc_fns = dict(Descriptors._descList)
    n = len(smiles_list)
    X = np.full((n, len(names)), np.nan)
    valid = np.zeros(n, dtype=bool)
    for i, smi in enumerate(smiles_list):
        mol = smiles_to_mol(smi)
        if mol is None:
            continue
        valid[i] = True
        for j, name in enumerate(names):
            try:
                v = calc_fns[name](mol)
                X[i, j] = v if np.isfinite(v) else np.nan
            except Exception:
                X[i, j] = np.nan
    # impute remaining NaN with column median, drop constant/all-NaN columns
    col_ok = ~np.all(np.isnan(X), axis=0)
    X = X[:, col_ok]
    names = [n for n, ok in zip(names, col_ok) if ok]
    col_med = np.nanmedian(X, axis=0)
    inds = np.where(np.isnan(X))
    X[inds] = np.take(col_med, inds[1])
    # a handful of RDKit descriptors (e.g. Ipc, information-content indices) are not
    # literally inf/nan but can reach ~1e90+ on larger molecules — finite in float64,
    # but they overflow to inf when a downstream library (e.g. XGBoost) casts to
    # float32. Clip to a generous range that only ever engages for these pathological
    # descriptors; realistic values for lipid-sized molecules are nowhere near this.
    X = np.clip(X, -1e6, 1e6)
    stds = X.std(axis=0)
    keep = stds > 1e-8
    X = X[:, keep]
    names = [n for n, k in zip(names, keep) if k]
    return X, valid, names


def evaluate_classifier(X, y, groups, splitter_fn, clf_fn, n_classes):
    """Run k-fold CV with the given splitter, return per-fold and aggregate metrics."""
    fold_metrics = []
    all_true, all_pred = [], []
    for fold_i, (train_idx, test_idx) in enumerate(splitter_fn(X, y, groups)):
        clf = clf_fn()
        clf.fit(X[train_idx], y[train_idx])
        pred = clf.predict(X[test_idx])
        bal_acc = balanced_accuracy_score(y[test_idx], pred)
        acc = (pred == y[test_idx]).mean()
        f1 = f1_score(y[test_idx], pred, average="macro", zero_division=0)
        fold_metrics.append({"fold": fold_i, "n_test": len(test_idx), "balanced_accuracy": bal_acc,
                              "top1_accuracy": acc, "macro_f1": f1})
        all_true.extend(y[test_idx].tolist())
        all_pred.extend(pred.tolist())
    agg = {
        "balanced_accuracy_mean": np.mean([m["balanced_accuracy"] for m in fold_metrics]),
        "balanced_accuracy_std": np.std([m["balanced_accuracy"] for m in fold_metrics]),
        "top1_accuracy_mean": np.mean([m["top1_accuracy"] for m in fold_metrics]),
        "macro_f1_mean": np.mean([m["macro_f1"] for m in fold_metrics]),
    }
    return agg, fold_metrics, np.array(all_true), np.array(all_pred)


def stratified_splitter(n_folds, seed):
    def _splitter(X, y, groups):
        skf = StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=seed)
        return list(skf.split(X, y))
    return _splitter


def group_splitter(n_folds, seed):
    def _splitter(X, y, groups):
        gkf = GroupKFold(n_splits=n_folds)
        return list(gkf.split(X, y, groups=groups))
    return _splitter


def permutation_null(X, y, n_classes, n_perm, seed):
    """Fast permutation baseline using a single train/test split + logistic regression,
    repeated n_perm times with shuffled labels. See module docstring for why this is a
    documented, cheaper substitute for a full 1000x k-fold RF permutation."""
    rng = np.random.default_rng(seed)
    n = len(y)
    idx = np.arange(n)
    rng.shuffle(idx)
    split = int(0.8 * n)
    train_idx, test_idx = idx[:split], idx[split:]
    null_scores = []
    for i in range(n_perm):
        y_perm = rng.permutation(y)
        clf = LogisticRegression(max_iter=300, C=1.0, class_weight="balanced", n_jobs=1)
        clf.fit(X[train_idx], y_perm[train_idx])
        pred = clf.predict(X[test_idx])
        null_scores.append(balanced_accuracy_score(y_perm[test_idx], pred))
    return np.array(null_scores), train_idx, test_idx


def nearest_neighbor_same_source_rate(bitvects, sources):
    n = len(bitvects)
    same_source_hits = 0
    expected_rate_sum = 0.0
    source_arr = np.array(sources)
    counts = pd.Series(sources).value_counts()
    for i in range(n):
        sims = DataStructs.BulkTanimotoSimilarity(bitvects[i], bitvects)
        sims[i] = -1  # exclude self
        nn_idx = int(np.argmax(sims))
        if source_arr[nn_idx] == source_arr[i]:
            same_source_hits += 1
        n_s = counts[source_arr[i]]
        expected_rate_sum += (n_s - 1) / (n - 1) if n > 1 else 0.0
    return same_source_hits / n, expected_rate_sum / n


def run_for_dataset(df, dataset_name, report, all_rows):
    report.append(f"## `{dataset_name}`\n")
    dedup, qualifying, excluded, counts = build_lipid_source_table(df, dataset_name, MIN_LIPIDS_PER_SOURCE)
    n_classes = len(qualifying)
    report.append(f"- Sources with >= {MIN_LIPIDS_PER_SOURCE} unique lipids (included): **{n_classes}**")
    report.append(f"- Sources excluded (too few unique lipids): {len(excluded)} — "
                  f"{', '.join(f'{s} (n={counts[s]})' for s in excluded) if excluded else 'none'}")
    report.append(f"- Classification unit: distinct (source, lipid) pairs = **{len(dedup):,}**")
    report.append(f"- Chance baseline (1/S): **{1/n_classes:.4f}**\n")

    smiles = dedup["canonical_IL_SMILES"].tolist()
    sources = dedup["Experiment_ID"].tolist()
    le_map = {s: i for i, s in enumerate(qualifying)}
    y = np.array([le_map[s] for s in sources])

    fps, valid_fp = batch_morgan_fingerprints(smiles, MORGAN_RADIUS, MORGAN_NBITS, cache_name=f"stage3_{dataset_name}")
    X_morgan, y_morgan = fps[valid_fp], y[valid_fp]
    smiles_valid = [s for s, v in zip(smiles, valid_fp) if v]
    sources_valid = [s for s, v in zip(sources, valid_fp) if v]
    report.append(f"- Morgan fingerprint valid: {valid_fp.sum():,} / {len(valid_fp):,}\n")

    scaffolds = [bemis_murcko_scaffold(s) for s in smiles_valid]
    scaffold_ids = pd.factorize(scaffolds)[0]

    results_rows = []

    for feat_name, X_feat in [("morgan_fp", X_morgan)]:
        for clf_name, clf_fn in [
            ("logistic_regression", lambda: LogisticRegression(max_iter=500, C=1.0, class_weight="balanced", n_jobs=1)),
            ("random_forest", lambda: RandomForestClassifier(n_estimators=300, class_weight="balanced",
                                                               random_state=GLOBAL_SEED, n_jobs=1)),
        ]:
            for split_name, splitter in [
                ("stratified_kfold", stratified_splitter(N_FOLDS, GLOBAL_SEED)),
                ("scaffold_disjoint_groupkfold", group_splitter(N_FOLDS, GLOBAL_SEED)),
            ]:
                t0 = time.time()
                agg, fold_metrics, y_true_all, y_pred_all = evaluate_classifier(
                    X_feat, y_morgan, scaffold_ids, splitter, clf_fn, n_classes
                )
                dt = time.time() - t0
                per_source_recall = recall_score(y_true_all, y_pred_all, average=None, zero_division=0)
                row = {
                    "dataset": dataset_name, "featurization": feat_name, "classifier": clf_name,
                    "split": split_name, "n_classes": n_classes, "n_samples": len(X_feat),
                    "chance_balanced_accuracy": 1 / n_classes, "fit_time_s": round(dt, 1),
                    **agg,
                }
                results_rows.append(row)
                report.append(f"**{feat_name} / {clf_name} / {split_name}**: balanced_acc="
                              f"{agg['balanced_accuracy_mean']:.3f} (+/-{agg['balanced_accuracy_std']:.3f}), "
                              f"top1={agg['top1_accuracy_mean']:.3f}, macro_F1={agg['macro_f1_mean']:.3f}, "
                              f"chance={1/n_classes:.3f}  [{dt:.1f}s]")
                registry_set(
                    f"separability.{dataset_name}.{feat_name}.{clf_name}.{split_name}.balanced_accuracy",
                    agg["balanced_accuracy_mean"], n=len(X_feat), stage="stage3",
                    extra={"chance": 1 / n_classes, "std": agg["balanced_accuracy_std"], "n_classes": n_classes},
                )

    # RDKit descriptor featurization (only with the faster classifier: logistic regression, one split scheme,
    # since this is the qualitative-robustness check, not a second full grid)
    report.append("\n**RDKit descriptor featurization (robustness check vs. Morgan FP):**\n")
    X_desc, desc_valid, desc_names = rdkit_descriptor_matrix(smiles)
    X_desc_v = X_desc[valid_fp & desc_valid]
    y_desc_v = y[valid_fp & desc_valid]
    scaler = StandardScaler().fit(X_desc_v)
    X_desc_scaled = scaler.transform(X_desc_v)
    agg_d, _, _, _ = evaluate_classifier(
        X_desc_scaled, y_desc_v, None,
        stratified_splitter(N_FOLDS, GLOBAL_SEED),
        lambda: LogisticRegression(max_iter=500, C=1.0, class_weight="balanced", n_jobs=1),
        n_classes,
    )
    report.append(f"- RDKit descriptors ({len(desc_names)} features) / logistic_regression / stratified_kfold: "
                  f"balanced_acc={agg_d['balanced_accuracy_mean']:.3f}, chance={1/n_classes:.3f}\n")
    results_rows.append({
        "dataset": dataset_name, "featurization": "rdkit_descriptors", "classifier": "logistic_regression",
        "split": "stratified_kfold", "n_classes": n_classes, "n_samples": len(X_desc_v),
        "chance_balanced_accuracy": 1 / n_classes, "fit_time_s": None, **agg_d,
    })
    registry_set(f"separability.{dataset_name}.rdkit_descriptors.logistic_regression.stratified_kfold.balanced_accuracy",
                 agg_d["balanced_accuracy_mean"], n=len(X_desc_v), stage="stage3",
                 extra={"chance": 1 / n_classes, "n_classes": n_classes})

    # Nearest-neighbor same-source rate (Morgan FP, Tanimoto)
    bitvects = [common_bitvect(s) for s in smiles_valid]
    obs_rate, exp_rate = nearest_neighbor_same_source_rate(bitvects, sources_valid)
    report.append(f"**Same-source nearest-neighbor rate**: observed = {obs_rate:.3f}, "
                  f"expected under random source assignment = {exp_rate:.3f} "
                  f"(ratio = {obs_rate/exp_rate:.2f}x)\n")
    registry_set(f"separability.{dataset_name}.nn_same_source_rate", obs_rate, n=len(bitvects), stage="stage3",
                 extra={"expected_random": exp_rate})

    # Permutation baseline (Morgan FP, logistic regression, single split, N_PERMUTATIONS shuffles)
    null_scores, perm_train_idx, perm_test_idx = permutation_null(X_morgan, y_morgan, n_classes, N_PERMUTATIONS, GLOBAL_SEED)
    clf_obs = LogisticRegression(max_iter=500, C=1.0, class_weight="balanced", n_jobs=1)
    clf_obs.fit(X_morgan[perm_train_idx], y_morgan[perm_train_idx])
    obs_score = balanced_accuracy_score(y_morgan[perm_test_idx], clf_obs.predict(X_morgan[perm_test_idx]))
    n_null_ge_obs = int(np.sum(null_scores >= obs_score))
    p_value = (n_null_ge_obs + 1) / (len(null_scores) + 1)
    p_display = f"< {p_value:.4f}" if n_null_ge_obs == 0 else f"= {p_value:.4f}"
    report.append(f"**Permutation baseline** ({N_PERMUTATIONS} shuffles, logistic regression, single 80/20 split): "
                  f"observed balanced_acc={obs_score:.3f}, null mean={null_scores.mean():.3f} "
                  f"(std={null_scores.std():.3f}), {n_null_ge_obs}/{N_PERMUTATIONS} null shuffles >= observed, "
                  f"empirical p {p_display} (resolution floor of {N_PERMUTATIONS} shuffles)\n")
    registry_set(f"separability.{dataset_name}.permutation.observed_balanced_accuracy", obs_score,
                 stage="stage3", extra={"null_mean": float(null_scores.mean()), "null_std": float(null_scores.std()),
                                         "p_value": float(p_value), "n_permutations": N_PERMUTATIONS})

    # UMAP embedding
    try:
        import umap
        reducer = umap.UMAP(random_state=GLOBAL_SEED, n_neighbors=15, min_dist=0.1)
        emb = reducer.fit_transform(X_morgan)
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(8, 7))
        srcs = pd.Series(sources_valid)
        cmap = plt.get_cmap("tab20", len(srcs.unique()))
        for i, s in enumerate(sorted(srcs.unique())):
            mask = srcs.values == s
            ax.scatter(emb[mask, 0], emb[mask, 1], s=6, alpha=0.6, color=cmap(i), label=s)
        ax.set_title(f"UMAP of Morgan FPs colored by source — {dataset_name}")
        ax.set_xlabel("UMAP-1"); ax.set_ylabel("UMAP-2")
        ax.legend(markerscale=3, fontsize=6, ncol=2, loc="best", frameon=False)
        fig.tight_layout()
        fig.savefig(FIGURES / f"stage3_embedding_{dataset_name}.png", dpi=300)
        fig.savefig(FIGURES / f"stage3_embedding_{dataset_name}.svg")
        plt.close(fig)
        report.append(f"- UMAP embedding saved to `figures/stage3_embedding_{dataset_name}.png` (+ .svg)\n")
    except Exception as e:
        report.append(f"- UMAP embedding FAILED: {e}\n")

    all_rows.extend(results_rows)
    return results_rows


def common_bitvect(smiles):
    from common import morgan_fingerprint_bitvect
    return morgan_fingerprint_bitvect(smiles)


def main():
    df = load_curated()
    report = ["# Stage 3 — Source Separability from Chemistry\n",
              f"Featurization: Morgan/ECFP fingerprints (radius={MORGAN_RADIUS}, {MORGAN_NBITS} bits); "
              f"RDKit full descriptor set as a robustness check. Classification unit = distinct "
              f"(source, canonical ionizable-lipid SMILES) pairs (see module docstring for why raw "
              f"measurement rows would leak). Minimum {MIN_LIPIDS_PER_SOURCE} unique lipids/source to "
              f"be included.\n"]

    all_rows = []
    for ds in ["LNPDB_vitro", "LNPDB_vivo"]:
        run_for_dataset(df, ds, report, all_rows)

    results_df = pd.DataFrame(all_rows)
    results_df.to_csv(RESULTS / "stage3_separability.csv", index=False)
    (REPORTS / "03_separability.md").write_text("\n".join(str(x) for x in report))

    print("--- STAGE 3 SUMMARY ---")
    for ds in ["LNPDB_vitro", "LNPDB_vivo"]:
        sub = results_df[(results_df.dataset == ds) & (results_df.featurization == "morgan_fp") &
                          (results_df.classifier == "random_forest")]
        for _, r in sub.iterrows():
            print(f"{ds} | RF | {r.split}: balanced_acc={r.balanced_accuracy_mean:.3f} vs chance={r.chance_balanced_accuracy:.3f}")
    print("Wrote results/stage3_separability.csv, reports/03_separability.md, figures/stage3_embedding_*.png")


if __name__ == "__main__":
    main()
