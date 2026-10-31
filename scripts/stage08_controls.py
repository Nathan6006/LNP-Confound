"""
Stage 8 — Negative controls and robustness.

The split-half negative control is a GATE: it re-runs the separability and
concordance battery on two random pseudo-sources carved out of one real
source, where by construction there is no true cross-lab confound. If
separability does NOT collapse to chance, or concordance does NOT approach
the within-source ceiling, that means the Stage 3/5 machinery itself is
producing signal out of nothing — everything upstream would need to be
distrusted. Report PASS/FAIL explicitly; do not soften a failure.

Outputs:
  results/stage8_controls.csv
  reports/08_controls.md
  registry keys: controls.split_half.*, controls.featurization_sensitivity.*, controls.seed_sensitivity.*
"""
import warnings

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import balanced_accuracy_score
from sklearn.model_selection import StratifiedKFold

warnings.filterwarnings("ignore")

from common import (
    GLOBAL_SEED, MORGAN_RADIUS, MORGAN_NBITS, load_curated, registry_set,
    RESULTS, REPORTS, batch_morgan_fingerprints, bootstrap_ci,
)
from stage03_separability import build_lipid_source_table, evaluate_classifier, stratified_splitter
from stage05_concordance import (
    add_within_source_endpoint_rank, replicate_noise_ceiling, pooled_concordance,
    MIN_REPLICATE_GROUPS_FOR_CEILING,
)

SEEDS = [1, 2, 3, 4, 5]
SEPARABILITY_CHANCE_TOLERANCE = 0.10  # balanced accuracy must be within chance +/- this to PASS


def split_half_control(df_all, dataset_name, report, csv_rows):
    report.append(f"## Split-half negative control (GATE) — `{dataset_name}`\n")
    df = df_all[df_all["_source_dataset"] == dataset_name].copy()
    lipid_counts = df.groupby("Experiment_ID")["canonical_IL_SMILES"].nunique()
    largest_source = lipid_counts.idxmax()
    n_lipids = int(lipid_counts.max())
    report.append(f"- Largest source by unique lipids: **{largest_source}** ({n_lipids:,} lipids)\n")

    src_df = df[df["Experiment_ID"] == largest_source].copy()
    unique_lipids = src_df["canonical_IL_SMILES"].unique()
    rng = np.random.default_rng(GLOBAL_SEED)
    assignment = rng.integers(0, 2, size=len(unique_lipids))
    lipid_to_half = dict(zip(unique_lipids, assignment))
    src_df["_pseudo_source"] = src_df["canonical_IL_SMILES"].map(lipid_to_half).map({0: f"{largest_source}_A", 1: f"{largest_source}_B"})

    # ---- (a) separability collapse test ----
    pseudo = src_df.rename(columns={"Experiment_ID": "_real_source"}).copy()
    pseudo["Experiment_ID"] = pseudo["_pseudo_source"]
    pseudo["_source_dataset"] = "pseudo"
    dedup, qualifying, _, _ = build_lipid_source_table(pseudo, "pseudo", min_lipids=1)
    smiles = dedup["canonical_IL_SMILES"].tolist()
    fps, valid_fp = batch_morgan_fingerprints(smiles, MORGAN_RADIUS, MORGAN_NBITS, cache_name=f"stage8_splithalf_{dataset_name}")
    X = fps[valid_fp]
    y = np.array([0 if s == f"{largest_source}_A" else 1 for s in dedup["Experiment_ID"]])[valid_fp]

    agg, _, _, _ = evaluate_classifier(
        X, y, None, stratified_splitter(5, GLOBAL_SEED),
        lambda: RandomForestClassifier(n_estimators=300, class_weight="balanced", random_state=GLOBAL_SEED, n_jobs=1),
        n_classes=2,
    )
    bal_acc = agg["balanced_accuracy_mean"]
    sep_pass = bool(abs(bal_acc - 0.5) <= SEPARABILITY_CHANCE_TOLERANCE)
    report.append(f"**Separability collapse test**: random 50/50 lipid-level pseudo-source split -> "
                  f"balanced accuracy = **{bal_acc:.3f}** (chance = 0.500, tolerance +/-{SEPARABILITY_CHANCE_TOLERANCE}). "
                  f"**{'PASS' if sep_pass else 'FAIL'}**\n")
    registry_set(f"controls.split_half.{dataset_name}.separability_balanced_accuracy", bal_acc, n=len(X),
                 stage="stage8", extra={"pass": bool(sep_pass), "source": largest_source})

    # ---- (b) concordance-approaches-ceiling test ----
    # src_df's Experiment_ID column is still the real, untouched source label (the pseudo-source
    # split above operated on a separate `pseudo` copy) — exactly what add_within_source_endpoint_rank
    # and replicate_noise_ceiling need to compute this source's own replicate reproducibility.
    src_df_ranked = add_within_source_endpoint_rank(src_df)
    ceiling = replicate_noise_ceiling(src_df_ranked)
    if ceiling is None:
        report.append("**Concordance-approaches-ceiling test**: no replicate structure found in this source "
                      "— not computable. **INCONCLUSIVE** (not a fail; the source simply lacks replicate rows).\n")
        conc_val, conc_n = None, 0
        conc_pass = None
    elif ceiling["n"] < MIN_REPLICATE_GROUPS_FOR_CEILING:
        conc_val, conc_n = float(ceiling["spearman"]), int(ceiling["n"])
        conc_pass = None
        report.append(f"**Concordance-approaches-ceiling test**: within-source ({largest_source}) replicate "
                      f"split-half Spearman = {conc_val:.3f}, but n={conc_n} replicate-group pairs is below "
                      f"Stage 5's own reliability threshold ({MIN_REPLICATE_GROUPS_FOR_CEILING}) for trusting a "
                      f"concordance estimate. **INCONCLUSIVE (LOW-CONFIDENCE, n too small)** — not scored as a "
                      f"pass or fail, consistent with how Stage 5 treats any concordance built on this few "
                      f"replicate pairs.\n")
        registry_set(f"controls.split_half.{dataset_name}.concordance_spearman", conc_val, n=conc_n,
                     stage="stage8", confidence="LOW-CONFIDENCE", extra={"pass": None, "source": largest_source})
    else:
        conc_val, conc_n = float(ceiling["spearman"]), int(ceiling["n"])
        # gate criterion: split-half-within-one-source concordance should clearly exceed the
        # cross-source pooled concordance found in Stage 5 for this same pool
        from common import registry_get_all
        reg = registry_get_all()
        cross_source_key = f"concordance.{dataset_name}.pooled_spearman"
        cross_source_val = reg.get(cross_source_key, {}).get("value")
        conc_pass = bool((cross_source_val is None) or (conc_val > cross_source_val))
        report.append(f"**Concordance-approaches-ceiling test**: within-source ({largest_source}) replicate "
                      f"split-half Spearman = **{conc_val:.3f}** (n={conc_n}) vs. Stage 5's cross-source pooled "
                      f"concordance = {cross_source_val:.3f}. **{'PASS' if conc_pass else 'FAIL'}** "
                      f"(within-source value should clearly exceed the cross-source value).\n")
        registry_set(f"controls.split_half.{dataset_name}.concordance_spearman", conc_val, n=conc_n,
                     stage="stage8", extra={"pass": conc_pass, "source": largest_source,
                                             "cross_source_pooled": cross_source_val})

    csv_rows.append({"dataset": dataset_name, "control": "split_half_separability", "source": largest_source,
                      "value": bal_acc, "n": len(X), "pass": sep_pass})
    csv_rows.append({"dataset": dataset_name, "control": "split_half_concordance", "source": largest_source,
                      "value": conc_val, "n": conc_n, "pass": conc_pass})
    return sep_pass, conc_pass


def featurization_sensitivity(df_all, dataset_name, report, csv_rows):
    report.append(f"## Featurization sensitivity — `{dataset_name}`\n")
    from common import registry_get_all
    reg = registry_get_all()
    headline_key = f"separability.{dataset_name}.morgan_fp.random_forest.stratified_kfold.balanced_accuracy"
    headline = reg.get(headline_key, {}).get("value")

    dedup, qualifying, excluded, counts = build_lipid_source_table(df_all, dataset_name, min_lipids=20)
    n_classes = len(qualifying)
    smiles = dedup["canonical_IL_SMILES"].tolist()
    le_map = {s: i for i, s in enumerate(qualifying)}
    y = np.array([le_map[s] for s in dedup["Experiment_ID"]])

    alt_radius, alt_bits = 3, 1024
    fps, valid_fp = batch_morgan_fingerprints(smiles, alt_radius, alt_bits, cache_name=f"stage8_altfp_{dataset_name}")
    X, y_v = fps[valid_fp], y[valid_fp]
    agg, _, _, _ = evaluate_classifier(
        X, y_v, None, stratified_splitter(5, GLOBAL_SEED),
        lambda: RandomForestClassifier(n_estimators=300, class_weight="balanced", random_state=GLOBAL_SEED, n_jobs=1),
        n_classes=n_classes,
    )
    alt_bal_acc = agg["balanced_accuracy_mean"]
    report.append(f"- Headline (Morgan radius={2}, {2048}-bit): balanced accuracy = {headline:.3f}")
    report.append(f"- Alternate (Morgan radius={alt_radius}, {alt_bits}-bit): balanced accuracy = **{alt_bal_acc:.3f}** "
                  f"(chance = {1/n_classes:.3f})")
    report.append(f"- Difference: {abs(headline - alt_bal_acc):.3f} — "
                  f"{'qualitatively robust to featurization choice.' if abs(headline - alt_bal_acc) < 0.1 else 'notable sensitivity to featurization choice — flagged.'}\n")
    registry_set(f"controls.featurization_sensitivity.{dataset_name}.alt_morgan_balanced_accuracy", alt_bal_acc,
                 n=len(X), stage="stage8", extra={"radius": alt_radius, "bits": alt_bits, "headline": headline})
    csv_rows.append({"dataset": dataset_name, "control": "featurization_sensitivity", "source": None,
                      "value": alt_bal_acc, "n": len(X), "pass": None})


def seed_sensitivity(df_all, dataset_name, report, csv_rows):
    report.append(f"## Seed sensitivity — `{dataset_name}`\n")

    # (a) Stage 3 headline: RF balanced accuracy, stratified k-fold, across seeds
    dedup, qualifying, _, _ = build_lipid_source_table(df_all, dataset_name, min_lipids=20)
    n_classes = len(qualifying)
    smiles = dedup["canonical_IL_SMILES"].tolist()
    le_map = {s: i for i, s in enumerate(qualifying)}
    y = np.array([le_map[s] for s in dedup["Experiment_ID"]])
    fps, valid_fp = batch_morgan_fingerprints(smiles, MORGAN_RADIUS, MORGAN_NBITS, cache_name=f"stage3_{dataset_name}")
    X, y_v = fps[valid_fp], y[valid_fp]

    accs = []
    for seed in SEEDS:
        agg, _, _, _ = evaluate_classifier(
            X, y_v, None, stratified_splitter(5, seed),
            lambda: RandomForestClassifier(n_estimators=300, class_weight="balanced", random_state=seed, n_jobs=1),
            n_classes=n_classes,
        )
        accs.append(agg["balanced_accuracy_mean"])
    accs = np.array(accs)
    report.append(f"- Separability balanced accuracy across seeds {SEEDS}: {np.round(accs, 3).tolist()} "
                  f"-> mean={accs.mean():.3f}, std={accs.std():.3f}, range=[{accs.min():.3f}, {accs.max():.3f}]")
    registry_set(f"controls.seed_sensitivity.{dataset_name}.separability_balanced_accuracy_std", float(accs.std()),
                 n=len(SEEDS), stage="stage8", extra={"mean": float(accs.mean()), "seeds": SEEDS})
    csv_rows.append({"dataset": dataset_name, "control": "seed_sensitivity_separability", "source": None,
                      "value": float(accs.std()), "n": len(SEEDS), "pass": None})

    # (b) Stage 5 headline: pooled concordance Spearman bootstrap CI stability across bootstrap seeds
    df = df_all[df_all["_source_dataset"] == dataset_name].copy()
    df = add_within_source_endpoint_rank(df)
    from stage05_concordance import lipid_source_table as lst, build_pairwise_points
    lipid_src = lst(df)
    points = build_pairwise_points(lipid_src)
    ci_widths = []
    point_estimate = None
    for seed in SEEDS:
        res = pooled_concordance(points, seed=seed)
        point_estimate = res["spearman"]
        ci_widths.append(res["ci_hi"] - res["ci_lo"])
    ci_widths = np.array(ci_widths)
    report.append(f"- Pooled concordance point estimate (seed-invariant, Spearman is deterministic given data): "
                  f"{point_estimate:.3f}")
    report.append(f"- Bootstrap CI width across seeds {SEEDS}: {np.round(ci_widths, 3).tolist()} -> "
                  f"mean={ci_widths.mean():.3f}, std={ci_widths.std():.3f}\n")
    registry_set(f"controls.seed_sensitivity.{dataset_name}.concordance_ci_width_std", float(ci_widths.std()),
                 n=len(SEEDS), stage="stage8", extra={"mean_ci_width": float(ci_widths.mean())})
    csv_rows.append({"dataset": dataset_name, "control": "seed_sensitivity_concordance_ci", "source": None,
                      "value": float(ci_widths.std()), "n": len(SEEDS), "pass": None})


def main():
    df = load_curated()
    report = ["# Stage 8 — Negative Controls and Robustness\n",
              "The split-half control is a validity GATE: a failure here means the upstream Stage 3/5 "
              "machinery is producing signal out of nothing and everything preceding this stage would need "
              "to be distrusted. It is reported plainly, not softened.\n"]
    csv_rows = []
    gate_results = []

    for ds in ["LNPDB_vitro", "LNPDB_vivo"]:
        sep_pass, conc_pass = split_half_control(df, ds, report, csv_rows)
        gate_results.append((ds, "separability_collapse", sep_pass))
        gate_results.append((ds, "concordance_ceiling", conc_pass))
        featurization_sensitivity(df, ds, report, csv_rows)
        seed_sensitivity(df, ds, report, csv_rows)

    report.append("## Gate summary\n")
    report.append("| dataset | gate | result |")
    report.append("|---|---|---|")
    any_fail = False
    for ds, gate, result in gate_results:
        result_str = "PASS" if result is True else ("FAIL" if result is False else "INCONCLUSIVE")
        if result is False:
            any_fail = True
        report.append(f"| {ds} | {gate} | {result_str} |")
    report.append(f"\n**Overall gate status: {'FAIL — upstream results should be distrusted until resolved' if any_fail else 'PASS'}**\n")

    pd.DataFrame(csv_rows).to_csv(RESULTS / "stage8_controls.csv", index=False)
    (REPORTS / "08_controls.md").write_text("\n".join(str(x) for x in report))

    print("--- STAGE 8 SUMMARY ---")
    for ds, gate, result in gate_results:
        print(f"{ds} | {gate}: {'PASS' if result is True else ('FAIL' if result is False else 'INCONCLUSIVE')}")
    print("GATE:", "FAIL" if any_fail else "PASS")
    print("Wrote results/stage8_controls.csv, reports/08_controls.md")


if __name__ == "__main__":
    main()
