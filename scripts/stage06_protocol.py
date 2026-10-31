"""
Stage 6 — Does the evaluation protocol inflate performance? (3x3 grid)

A self-contained model (Morgan FP + RDKit descriptors -> XGBoost) is trained
under three objectives:
  - pooled_regression: plain regression on pooled within-source-normalized rank.
  - within_source_pairwise: XGBRanker rank:pairwise, pairs formed only within a source.
  - within_source_lambdarank: graded-relevance LambdaRank (see scripts/lambdarank.py),
    ported from a sibling project's within_exp_lambdarank2.py. Same within-source
    pairing discipline as the pairwise objective, but pairs are weighted by a
    graded relevance gap (rel in {0,1,2,3}, rel=3 = top-decile-within-source "hit")
    via NDCG-style lambda weights, with budgeted stratified pair sampling that
    oversamples hit-involving pairs. See lambdarank.py's docstring for exactly how
    `rel` is derived here (there is no curated-hit label in this corpus, unlike the
    source project) and why 10% was chosen over that project's ~3%.

...and evaluated under three split schemes (random k-fold, Butina-cluster
scaffold-disjoint k-fold, leave-one-source-out). The point is the *relative*
change across the grid, not state-of-the-art accuracy — this is a from-scratch
ablation, not a reproduction of any published model.

Every cell reports two metrics: pooled Spearman (predictions vs. true pooled
within-source-normalized rank, concatenated across all test folds) and
within-experiment Spearman (Spearman computed separately inside each source's
test predictions, then size-weighted across sources — gauge-invariant to
per-source relabeling). LOSO additionally gets a per-fold n and a per-fold
permutation baseline (fast Ridge regressor on shuffled labels, not full
XGBoost refits — a documented compute-budget substitute, same rationale as
Stage 3's permutation test).

Each cell ALSO gets top-k ranking-quality metrics vs. a random baseline (graded
NDCG@k_e / NDCG_full / NDCG@{3,5,10} and hit_rate@{3,5,10}, each paired with a
1000-permutation random-ordering floor; EF@k reported as a ratio of means) —
ported from the sibling project's analyze.py evaluation pattern (see
scripts/lambdarank.py). Unit of analysis is an (experiment, fold) fragment:
LOSO folds give each source intact (primary home for hit metrics); random/
scaffold-disjoint k-fold splits fragment a source's rows across multiple
folds' test sets (secondary, partial admissibility — reported honestly, not
hidden).

Outputs:
  results/stage6_protocol_grid.csv
  results/stage6_topk_metrics.csv (per (objective, split, experiment, fold) fragment)
  figures/stage6_protocol_grid_{vitro,vivo}.png (+ .svg)
  reports/06_protocol.md
  registry keys: protocol.{vitro,vivo}.{objective}.{split}.pooled_spearman / within_experiment_spearman
                 protocol.{vitro,vivo}.{objective}.{split}.topk.{metric}
"""
import warnings

import numpy as np
import pandas as pd
import xgboost as xgb
from scipy.stats import spearmanr
from sklearn.linear_model import Ridge
from sklearn.model_selection import KFold, GroupKFold
from xgboost import XGBRegressor, XGBRanker

warnings.filterwarnings("ignore")

from common import (
    GLOBAL_SEED, MORGAN_RADIUS, MORGAN_NBITS, load_curated, registry_set,
    RESULTS, REPORTS, FIGURES, batch_morgan_fingerprints, morgan_fingerprint_bitvect,
    butina_clusters,
)
from stage03_separability import rdkit_descriptor_matrix
from stage05_concordance import add_within_source_endpoint_rank, lipid_source_table
from lambdarank import (
    compute_relevance_grades, make_lambdarank_objective,
    per_experiment_fold_metrics, aggregate_topk_metrics,
)

N_FOLDS = 5
BUTINA_CUTOFF = 0.4
N_PERMUTATIONS_PER_LOSO_FOLD = 20  # fast Ridge-based null, see module docstring
XGB_PARAMS = dict(n_estimators=200, max_depth=4, learning_rate=0.1, random_state=GLOBAL_SEED,
                   n_jobs=4, verbosity=0)
LAMBDARANK_NUM_BOOST_ROUND = 200  # matches XGB_PARAMS n_estimators for comparability
LAMBDARANK_BOOSTER_PARAMS = dict(max_depth=4, eta=0.1, tree_method="hist", seed=GLOBAL_SEED, nthread=4)
LAMBDARANK_HIT_FRAC = 0.10
TOPK_N_PERM = 1000  # random-baseline permutations, fixed a priori (never a knob) — matches source project
TOPK_K_FRAC, TOPK_K_MIN, TOPK_K_MAX = 0.10, 5, 50
TOPK_MIN_N = 3
OBJECTIVES = ["pooled_regression", "within_source_pairwise", "within_source_lambdarank"]


def within_experiment_spearman(sources, y_true, y_pred, min_n=3):
    """Size-weighted average of per-source Spearman(y_true, y_pred), sources with
    < min_n test points excluded (too few to correlate)."""
    df = pd.DataFrame({"source": sources, "y_true": y_true, "y_pred": y_pred})
    rhos, weights = [], []
    for src, g in df.groupby("source"):
        if len(g) < min_n or g["y_true"].nunique() < 2 or g["y_pred"].nunique() < 2:
            continue
        rho = spearmanr(g["y_true"], g["y_pred"]).statistic
        if not np.isnan(rho):
            rhos.append(rho)
            weights.append(len(g))
    if not rhos:
        return np.nan, 0
    return float(np.average(rhos, weights=weights)), len(rhos)


def fit_predict_pooled(X_train, y_train, X_test):
    model = XGBRegressor(**XGB_PARAMS)
    model.fit(X_train, y_train)
    return model.predict(X_test)


def fit_predict_ranker(X_train, y_train, groups_train, X_test):
    order = np.argsort(groups_train, kind="stable")
    X_sorted, y_sorted = X_train[order], y_train[order]
    _, group_sizes = np.unique(groups_train[order], return_counts=True)
    model = XGBRanker(objective="rank:pairwise", **XGB_PARAMS)
    model.fit(X_sorted, y_sorted, group=group_sizes)
    return model.predict(X_test)


def fit_predict_lambdarank(X_train, y_train, sources_train, rel_train, X_test, seed=GLOBAL_SEED):
    """rel_train must already be sliced from the GLOBALLY precomputed rel array (see
    module docstring / lambdarank.py: relevance grades are frozen once on the full
    per-source rank distribution, not recomputed on each training subset)."""
    obj_fn, structs, dropped = make_lambdarank_objective(
        sources_train, rel_train, beta=1.0, budget_B=1500, top_frac=0.25,
        top_rel_threshold=3, base_seed=seed,
    )
    dtrain = xgb.DMatrix(X_train, label=y_train)
    booster = xgb.train(LAMBDARANK_BOOSTER_PARAMS, dtrain, num_boost_round=LAMBDARANK_NUM_BOOST_ROUND, obj=obj_fn)
    pred = booster.predict(xgb.DMatrix(X_test))
    return pred, dropped


def random_kfold_splits(n, seed):
    kf = KFold(n_splits=N_FOLDS, shuffle=True, random_state=seed)
    return list(kf.split(np.arange(n)))


def butina_kfold_splits(cluster_ids, seed):
    gkf = GroupKFold(n_splits=N_FOLDS)
    return list(gkf.split(np.arange(len(cluster_ids)), groups=cluster_ids))


def loso_splits(sources):
    src_arr = np.array(sources)
    folds = []
    for s in sorted(set(src_arr)):
        test_idx = np.where(src_arr == s)[0]
        train_idx = np.where(src_arr != s)[0]
        folds.append((train_idx, test_idx, s))
    return folds


def loso_permutation_baseline(X_train, y_train, X_test, y_test, n_perm, seed):
    rng = np.random.default_rng(seed)
    null_scores = []
    for _ in range(n_perm):
        y_perm = rng.permutation(y_train)
        model = Ridge(alpha=1.0)
        model.fit(X_train, y_perm)
        pred = model.predict(X_test)
        if len(set(y_test)) < 2 or len(set(pred)) < 2:
            continue
        null_scores.append(spearmanr(y_test, pred).statistic)
    return np.array(null_scores)


def run_for_dataset(df_all, dataset_name, report, csv_rows, topk_csv_rows_out):
    report.append(f"## `{dataset_name}`\n")
    df = df_all[df_all["_source_dataset"] == dataset_name].copy()
    df = add_within_source_endpoint_rank(df)
    lipid_src = lipid_source_table(df)
    lipid_src = lipid_src.dropna(subset=["rank"]).reset_index(drop=True)

    smiles = lipid_src["canonical_IL_SMILES"].tolist()
    sources = lipid_src["Experiment_ID"].values
    y = lipid_src["rank"].values

    fps, valid_fp = batch_morgan_fingerprints(smiles, MORGAN_RADIUS, MORGAN_NBITS, cache_name=f"stage6_{dataset_name}")
    X_desc, valid_desc, _ = rdkit_descriptor_matrix(smiles)
    valid = valid_fp & valid_desc
    X = np.concatenate([fps[valid], X_desc[valid]], axis=1)
    y = y[valid]
    sources = sources[valid]
    smiles_valid = [s for s, v in zip(smiles, valid) if v]
    n = len(X)
    report.append(f"- Unit: distinct (source, lipid) pairs, n = {n:,}, features = Morgan({MORGAN_NBITS}) + "
                  f"RDKit descriptors({X_desc.shape[1]}) = {X.shape[1]} dims\n")

    # Relevance grades for the LambdaRank objective, frozen ONCE on the full per-source
    # rank distribution (see lambdarank.py docstring) — never recomputed per CV fold.
    rel_all = compute_relevance_grades(y, sources, hit_frac=LAMBDARANK_HIT_FRAC)
    n_hits = int((rel_all == 3).sum())
    report.append(f"- LambdaRank relevance grades: hit_frac={LAMBDARANK_HIT_FRAC} -> {n_hits:,} 'hits' "
                  f"(rel=3) of {n:,} lipids; level counts {np.bincount(rel_all).tolist()}\n")

    bitvects = [morgan_fingerprint_bitvect(s) for s in smiles_valid]
    cluster_ids, _ = butina_clusters(bitvects, cutoff=BUTINA_CUTOFF)
    n_clusters = len(set(cluster_ids))
    report.append(f"- Butina clusters (Tanimoto cutoff={BUTINA_CUTOFF}): {n_clusters}\n")

    splitters = {
        "random_kfold": lambda: random_kfold_splits(n, GLOBAL_SEED),
        "butina_scaffold_disjoint_kfold": lambda: butina_kfold_splits(cluster_ids, GLOBAL_SEED),
        "leave_one_source_out": lambda: [(tr, te) for tr, te, _ in loso_splits(sources)],
    }
    loso_fold_names = [s for _, _, s in loso_splits(sources)]

    grid_rows = []
    for objective in OBJECTIVES:
        dropped_log = []
        for split_name, split_fn in splitters.items():
            folds = split_fn()
            all_true, all_pred, all_src = [], [], []
            per_fold_records = []
            topk_fragment_rows = []
            for fold_i, fold in enumerate(folds):
                train_idx, test_idx = fold[0], fold[1]
                if len(train_idx) == 0 or len(test_idx) == 0:
                    continue
                if objective == "pooled_regression":
                    pred = fit_predict_pooled(X[train_idx], y[train_idx], X[test_idx])
                elif objective == "within_source_pairwise":
                    pred = fit_predict_ranker(X[train_idx], y[train_idx], sources[train_idx], X[test_idx])
                else:
                    pred, dropped = fit_predict_lambdarank(
                        X[train_idx], y[train_idx], sources[train_idx], rel_all[train_idx], X[test_idx],
                        seed=GLOBAL_SEED + fold_i,
                    )
                    dropped_log.extend(dropped)
                all_true.extend(y[test_idx].tolist())
                all_pred.extend(pred.tolist())
                all_src.extend(sources[test_idx].tolist())

                # Top-k metrics vs. random baseline, one (experiment, fold) fragment at a
                # time: LOSO gives each source its full test_idx (intact); random/butina
                # k-fold test_idx spans many sources, each only partially represented here.
                test_sources = sources[test_idx]
                test_rel = rel_all[test_idx]
                test_y = y[test_idx]
                for src in np.unique(test_sources):
                    m = test_sources == src
                    frag = per_experiment_fold_metrics(
                        test_y[m], pred[m], test_rel[m], src, fold_i,
                        min_n=TOPK_MIN_N, k_frac=TOPK_K_FRAC, k_min=TOPK_K_MIN, k_max=TOPK_K_MAX,
                        n_perm=TOPK_N_PERM, seed=GLOBAL_SEED,
                    )
                    topk_fragment_rows.append(frag)

                fold_rho = spearmanr(y[test_idx], pred).statistic if len(set(y[test_idx])) > 1 and len(set(pred)) > 1 else np.nan
                fold_record = {"objective": objective, "split": split_name, "fold": fold_i,
                                "n_test": len(test_idx), "fold_spearman": fold_rho}

                if split_name == "leave_one_source_out":
                    held_out_source = loso_fold_names[fold_i]
                    null_scores = loso_permutation_baseline(X[train_idx], y[train_idx], X[test_idx], y[test_idx],
                                                              N_PERMUTATIONS_PER_LOSO_FOLD, GLOBAL_SEED + fold_i)
                    null_mean = float(np.nanmean(null_scores)) if len(null_scores) else np.nan
                    null_std = float(np.nanstd(null_scores)) if len(null_scores) else np.nan
                    at_chance = (not np.isnan(fold_rho)) and (not np.isnan(null_std)) and \
                                (fold_rho <= null_mean + 1.645 * null_std)
                    fold_record.update({"source": held_out_source, "null_mean": null_mean, "null_std": null_std,
                                        "at_or_below_chance": bool(at_chance)})
                per_fold_records.append(fold_record)
                grid_rows.append(fold_record | {"dataset": dataset_name})

            pooled_rho = spearmanr(all_true, all_pred).statistic if len(set(all_true)) > 1 and len(set(all_pred)) > 1 else np.nan
            we_rho, we_n_sources = within_experiment_spearman(all_src, all_true, all_pred)

            report.append(f"**{objective} / {split_name}**: pooled Spearman = **{pooled_rho:.3f}**, "
                          f"within-experiment Spearman = **{we_rho:.3f}** (n_sources={we_n_sources}), "
                          f"n_test_total={len(all_true)}")
            registry_set(f"protocol.{dataset_name}.{objective}.{split_name}.pooled_spearman", pooled_rho,
                         n=len(all_true), stage="stage6")
            registry_set(f"protocol.{dataset_name}.{objective}.{split_name}.within_experiment_spearman", we_rho,
                         n=we_n_sources, stage="stage6")

            if split_name == "leave_one_source_out":
                n_at_chance = sum(1 for r in per_fold_records if r.get("at_or_below_chance"))
                n_folds_total = len(per_fold_records)
                report.append(f"  - LOSO folds at-or-below permutation-null chance level: "
                              f"**{n_at_chance}/{n_folds_total}**")
                registry_set(f"protocol.{dataset_name}.{objective}.leave_one_source_out.n_folds_at_chance",
                             n_at_chance, n=n_folds_total, stage="stage6")

            # ---- Top-k metrics vs. random baseline (aggregated over this cell's fragments) ----
            topk_df = pd.DataFrame(topk_fragment_rows)
            if topk_df.empty:
                topk_df = pd.DataFrame(columns=["n_vals", "n_rel_levels", "has_hits"])
            topk_df["objective"], topk_df["split"], topk_df["dataset"] = objective, split_name, dataset_name
            topk_csv_rows_out.append(topk_df)
            agg = aggregate_topk_metrics(topk_df).set_index("metric") if len(topk_df) else \
                pd.DataFrame(columns=["mean", "sqrt_n_mean", "n_weighted_mean", "n_frag"]).rename_axis("metric")

            n_considered = int(((topk_df["n_vals"] >= TOPK_MIN_N) & (topk_df["n_rel_levels"] >= 2)).sum()) if len(topk_df) else 0
            n_admissible_hits = int(topk_df["has_hits"].sum()) if "has_hits" in topk_df.columns and len(topk_df) else 0
            role = "intact experiments" if split_name == "leave_one_source_out" else "row-fragmented, partial admissibility"

            def _m(name):
                return agg.loc[name, "mean"] if name in agg.index else float("nan")

            report.append(
                f"  - top-k vs random ({role}): NDCG@k_e={_m('ndcg@k_e'):.3f} "
                f"(random={_m('ndcg@k_e_rand'):.3f}); hit_rate@5={_m('hit_rate@5'):.3f} "
                f"(random={_m('hit_rate@5_rand'):.3f}, EF@5={_m('ef@5'):.2f}x); "
                f"hit-admissible fragments={n_admissible_hits}/{n_considered}"
            )
            registry_set(f"protocol.{dataset_name}.{objective}.{split_name}.topk.ndcg_at_k_e",
                         _m("ndcg@k_e"), n=int(agg.loc["ndcg@k_e", "n_frag"]) if "ndcg@k_e" in agg.index else 0,
                         stage="stage6", extra={"random_baseline": _m("ndcg@k_e_rand")})
            registry_set(f"protocol.{dataset_name}.{objective}.{split_name}.topk.hit_rate_at_5",
                         _m("hit_rate@5"), n=n_admissible_hits, stage="stage6",
                         extra={"random_baseline": _m("hit_rate@5_rand"), "ef_at_5": _m("ef@5"),
                                "n_hit_admissible": n_admissible_hits, "n_considered": n_considered})

        if objective == "within_source_lambdarank" and dropped_log:
            unique_dropped = sorted(set(dropped_log))
            report.append(f"  - LambdaRank dropped {len(unique_dropped)} distinct source(s) across all folds "
                          f"(IDCG=0 or single relevance level present in that training subset): "
                          f"{unique_dropped[:15]}{' ...' if len(unique_dropped) > 15 else ''}")
    report.append("\n")

    grid_df = pd.DataFrame(grid_rows)
    csv_rows.append(grid_df)

    # ---- figure: 3x3 grid of pooled vs within-experiment spearman ----
    import matplotlib.pyplot as plt
    cells = []
    for objective in OBJECTIVES:
        for split_name in splitters:
            key_pooled = f"protocol.{dataset_name}.{objective}.{split_name}.pooled_spearman"
            key_we = f"protocol.{dataset_name}.{objective}.{split_name}.within_experiment_spearman"
            from common import registry_get_all
            reg = registry_get_all()
            cells.append({"objective": objective, "split": split_name,
                          "pooled": reg[key_pooled]["value"], "within_experiment": reg[key_we]["value"]})
    cell_df = pd.DataFrame(cells)
    split_order = list(splitters.keys())  # random -> scaffold-disjoint -> LOSO, so the
                                           # inflation trend reads left-to-right in the figure
    fig, axes = plt.subplots(1, 2, figsize=(11, 5.5))
    for ax, metric in zip(axes, ["pooled", "within_experiment"]):
        pivot = cell_df.pivot(index="objective", columns="split", values=metric).reindex(
            index=OBJECTIVES, columns=split_order
        )
        im = ax.imshow(pivot.values, vmin=-0.3, vmax=1.0, cmap="RdYlGn", aspect="auto")
        ax.set_xticks(range(len(pivot.columns))); ax.set_xticklabels(pivot.columns, rotation=30, ha="right", fontsize=8)
        ax.set_yticks(range(len(pivot.index))); ax.set_yticklabels(pivot.index, fontsize=8)
        for i in range(pivot.shape[0]):
            for j in range(pivot.shape[1]):
                v = pivot.values[i, j]
                ax.text(j, i, f"{v:.2f}" if not pd.isna(v) else "NA", ha="center", va="center", fontsize=9)
        ax.set_title(f"{metric} Spearman")
        fig.colorbar(im, ax=ax, fraction=0.046)
    fig.suptitle(f"Protocol comparison grid — {dataset_name}")
    fig.tight_layout()
    fig.savefig(FIGURES / f"stage6_protocol_grid_{dataset_name}.png", dpi=300)
    fig.savefig(FIGURES / f"stage6_protocol_grid_{dataset_name}.svg")
    plt.close(fig)
    report.append(f"Figure: `figures/stage6_protocol_grid_{dataset_name}.png` (+ .svg)\n")


def main():
    df = load_curated()
    report = ["# Stage 6 — Protocol Comparison (3x3 Grid): Does the Evaluation Split Inflate Performance?\n",
              f"Model: XGBoost (Morgan FP {MORGAN_NBITS}-bit + RDKit descriptors), {N_FOLDS}-fold CV. "
              f"Three objectives (pooled regression / within-source pairwise ranking / within-source "
              f"graded-relevance LambdaRank, ported from a sibling project — see scripts/lambdarank.py) "
              f"x three splits (random, Butina-cluster scaffold-disjoint, leave-one-source-out). This is a "
              f"from-scratch ablation built for this comparison, not a reproduction of any published "
              f"architecture — the point is the *relative* change across cells.\n"]
    csv_rows = []
    topk_csv_rows = []

    for ds in ["LNPDB_vitro", "LNPDB_vivo"]:
        run_for_dataset(df, ds, report, csv_rows, topk_csv_rows)

    pd.concat(csv_rows, ignore_index=True).to_csv(RESULTS / "stage6_protocol_grid.csv", index=False)
    pd.concat(topk_csv_rows, ignore_index=True).to_csv(RESULTS / "stage6_topk_metrics.csv", index=False)
    (REPORTS / "06_protocol.md").write_text("\n".join(str(x) for x in report))

    print("--- STAGE 6 SUMMARY ---")
    from common import registry_get_all
    reg = registry_get_all()
    for ds in ["LNPDB_vitro", "LNPDB_vivo"]:
        for obj in OBJECTIVES:
            for split in ["random_kfold", "butina_scaffold_disjoint_kfold", "leave_one_source_out"]:
                k = f"protocol.{ds}.{obj}.{split}.pooled_spearman"
                if k in reg:
                    print(f"{ds} | {obj} | {split}: pooled={reg[k]['value']:.3f}")
    print("Wrote results/stage6_protocol_grid.csv, reports/06_protocol.md, figures/stage6_protocol_grid_*.png")


if __name__ == "__main__":
    main()
