"""
lambdarank.py - Within-source hit-status LambdaRank objective + metrics.

Ported from a sibling project's `within_exp_lambdarank2.py`
(~/Downloads/LNP_ML/scripts/within_exp_lambdarank2.py), which implements a
graded-relevance LambdaRank objective for within-experiment ranking. Loss
math, pair sampler, and gradient/hessian derivation are UNCHANGED from that
source; only the relevance-label construction differs (see
`compute_relevance_grades` below).

ADAPTATION FOR THIS PROJECT: the source implementation expects `rel` in
{0,1,2,3} to be precomputed externally by a `label_rel.py` script tied to
that project's wet-lab-curated "hit" flags (rel==3 = a curated hit, ~3% of
rows there). That script does not exist in this project — our LNPDB corpus
is a literature aggregation with no curated-hit column. `rel` is instead
derived here from the within-source-and-endpoint normalized rank (the same
quantity Stages 5-7 already compute): the **top 10% within each source**
becomes rel=3 ("hit"), and the remaining 90% is split into three
equal-sized tiers (rel=0/1/2). 10% (not the source project's ~3%) is a
documented choice, not a default we inherited — our per-source lipid counts
are much smaller, and 3% would leave most sources with zero qualifying
"hits" and be dropped outright by `compute_experiment_relevance` below.

Loss:
    rel_i in {0,1,2,3}, gain_i = 2^rel_i - 1        (exponential gain)
    d_i   = 1 / log2(1 + pos_i)                     (position discount)
    IDCG_e = sum_m gain_(m) / log2(1+m)             (per-experiment, per-split)
    For every eligible pair (k,j) with rel_k > rel_j in experiment e:
        dZ_kj  = (gain_k - gain_j) * |d_k - d_j| / IDCG_e
        L_kj   = w_k w_j dZ_kj * [-log sigma(beta (s_k - s_j))]
    grad/hess follow the standard LambdaRank pairwise-logistic derivation.

Sampling, per experiment per boosting round:
    - Exhaustive if n*(n-1)/2 <= 2*budget_B.
    - Otherwise a budgeted, stratified sample: `top_frac` of the budget is an
      "anchor" stream (one endpoint drawn from the hit set rel>=top_rel_threshold,
      the other uniformly), the rest is a uniform stream (both endpoints
      uniform). Each stream over-draws 3x then filters to eligible
      (rel_a != rel_b) pairs, since tied-relevance pairs have dZ=0.
    - Resampled every boosting round, seeded off `base_seed + round` so a run
      is reproducible but each round sees a different pair subset (variance
      reduction over training).
"""
import numpy as np
import pandas as pd
from scipy.special import expit


# =============================================================================
# Relevance-grade construction (the one piece NOT ported — see module docstring)
# =============================================================================

def compute_relevance_grades(rank, source_ids, hit_frac=0.10):
    """Derive rel in {0,1,2,3} from a continuous within-source-normalized rank
    (already in [0,1], comparable only within a source — never across sources).
    Top `hit_frac` fraction within each source -> rel=3 ("hit"); the rest split
    into three ~equal tiers by rank -> rel=0/1/2. Monotone in rank by construction.
    """
    rank = np.asarray(rank, dtype=np.float64)
    source_ids = np.asarray(source_ids)
    rel = np.zeros(len(rank), dtype=np.int64)
    for src in np.unique(source_ids):
        idx = np.where(source_ids == src)[0]
        r = rank[idx]
        n = len(idx)
        order = np.argsort(r)  # ascending: lowest rank first
        ranks_asc = np.empty(n, dtype=np.int64)
        ranks_asc[order] = np.arange(n)
        n_hit = max(1, int(round(hit_frac * n))) if n >= 4 else 0
        n_rest = n - n_hit
        grade = np.zeros(n, dtype=np.int64)
        if n_hit > 0:
            hit_cut = n_rest
            is_hit = ranks_asc >= hit_cut
            tier_size = max(1, n_rest // 3)
            grade[~is_hit & (ranks_asc < tier_size)] = 0
            grade[~is_hit & (ranks_asc >= tier_size) & (ranks_asc < 2 * tier_size)] = 1
            grade[~is_hit & (ranks_asc >= 2 * tier_size)] = 2
            grade[is_hit] = 3
        rel[idx] = grade
    return rel


# =============================================================================
# Internal utilities (verbatim math from the source implementation)
# =============================================================================

def _ideal_dcg(gain, k=None):
    g = np.sort(np.asarray(gain, dtype=np.float64))[::-1]
    if k is not None:
        g = g[:k]
    disc = 1.0 / np.log2(np.arange(2, g.size + 2))
    return float(np.sum(g * disc))


def _positions_discount(scores):
    n = scores.size
    order = np.argsort(-scores, kind="stable")
    pos = np.empty(n, dtype=np.float64)
    pos[order] = np.arange(1, n + 1, dtype=np.float64)
    disc = 1.0 / np.log2(1.0 + pos)
    return pos, disc


def _sample_pairs(rng, n, top_local, rel, budget_B, top_frac):
    """Return (a, b) local index arrays of ELIGIBLE (rel_a != rel_b) pairs, |pairs| <= B."""
    n_pairs_all = n * (n - 1) // 2
    if n_pairs_all <= 2 * budget_B:
        ii, jj = np.triu_indices(n, k=1)
        elig = rel[ii] != rel[jj]
        return ii[elig], jj[elig]

    OVER = 3
    n_top = int(round(top_frac * budget_B))
    n_uni = budget_B - n_top
    a_parts, b_parts = [], []
    if n_top > 0 and top_local.size > 0:
        a = rng.choice(top_local, size=n_top * OVER)
        b = rng.integers(0, n, size=n_top * OVER)
        keep = rel[a] != rel[b]
        a_parts.append(a[keep][:n_top]); b_parts.append(b[keep][:n_top])
    if n_uni > 0:
        a = rng.integers(0, n, size=n_uni * OVER)
        b = rng.integers(0, n, size=n_uni * OVER)
        keep = rel[a] != rel[b]
        a_parts.append(a[keep][:n_uni]); b_parts.append(b[keep][:n_uni])
    if not a_parts:
        return np.empty(0, np.int64), np.empty(0, np.int64)
    return np.concatenate(a_parts), np.concatenate(b_parts)


def _accumulate_experiment(grad, hess, s_full, struct, gweight, beta, a, b):
    """LambdaRank grad/hess for one experiment's eligible pairs."""
    idx = struct["idx"]
    if a.size == 0:
        return 0.0
    s = s_full[idx]
    _, disc = _positions_discount(s)
    rel, gain, idcg = struct["rel"], struct["gain"], struct["idcg"]

    a_is_k = rel[a] > rel[b]
    k = np.where(a_is_k, a, b)
    j = np.where(a_is_k, b, a)

    dgain = gain[k] - gain[j]
    dZ = dgain * np.abs(disc[k] - disc[j]) / idcg
    w_pair = gweight[idx][k] * gweight[idx][j] * dZ

    ds = s[k] - s[j]
    P = expit(beta * ds)
    g = w_pair * beta * (P - 1.0)
    h = w_pair * (beta * beta) * P * (1.0 - P)

    gk, gj = idx[k], idx[j]
    np.add.at(grad, gk, g)
    np.add.at(grad, gj, -g)
    np.add.at(hess, gk, h)
    np.add.at(hess, gj, h)

    return float(np.sum(w_pair * -np.log(np.clip(P, 1e-12, 1.0))))


def group_indices(source_ids, min_size=2):
    groups = {}
    for i, e in enumerate(source_ids):
        groups.setdefault(e, []).append(i)
    return {e: np.asarray(ix, dtype=np.int64) for e, ix in groups.items() if len(ix) >= min_size}


def compute_experiment_relevance(precomputed_rel, source_ids, top_rel_threshold=3, min_size=2):
    """Freeze per-source relevance structures. rel/gain are precomputed once on the
    full dataset (see module docstring); IDCG is computed from rows present in the
    CURRENT split only (so a LOSO/k-fold training subset gets its own IDCG)."""
    precomputed_rel = np.asarray(precomputed_rel, dtype=np.int64)
    groups = group_indices(source_ids, min_size=min_size)
    out, dropped = {}, []
    for e, idx in groups.items():
        rel = precomputed_rel[idx]
        gain = 2.0 ** rel - 1.0
        idcg = _ideal_dcg(gain)
        n_levels = int(np.unique(rel).size)
        if idcg <= 0.0 or n_levels < 2:
            dropped.append((e, "idcg=0" if idcg <= 0.0 else "single_level"))
            continue
        out[e] = {
            "idx": np.asarray(idx, dtype=np.int64),
            "rel": rel, "gain": gain, "idcg": idcg,
            "top_local": np.where(rel >= top_rel_threshold)[0].astype(np.int64),
            "n_levels": n_levels,
        }
    return out, dropped


# =============================================================================
# Objective factory (for xgboost's low-level xgb.train(obj=...) API)
# =============================================================================

def make_lambdarank_objective(
    source_ids, precomputed_rel,
    beta=1.0, budget_B=1500, top_frac=0.25, top_rel_threshold=3,
    base_seed=0, hess_floor=1e-6,
):
    """Build an xgboost custom objective closed over one training set's source_ids
    and precomputed rel. Pass to `xgb.train(params, dtrain, obj=this)`."""
    structs, dropped = compute_experiment_relevance(precomputed_rel, source_ids, top_rel_threshold, min_size=2)
    n_valid = len(structs)
    state = {"round": 0}

    def obj(predt, dtrain):
        predt = np.asarray(predt, dtype=np.float64)
        n_rows = predt.size
        w = dtrain.get_weight()
        w = np.ones(n_rows) if (w is None or len(w) == 0) else np.asarray(w, dtype=np.float64)

        grad = np.zeros(n_rows, dtype=np.float64)
        hess = np.zeros(n_rows, dtype=np.float64)
        if n_valid == 0:
            return grad.astype(np.float32), np.full(n_rows, hess_floor, np.float32)

        rng = np.random.default_rng(base_seed + state["round"])
        state["round"] += 1

        for struct in structs.values():
            n = struct["idx"].size
            a, b = _sample_pairs(rng, n, struct["top_local"], struct["rel"], budget_B, top_frac)
            _accumulate_experiment(grad, hess, predt, struct, w, beta, a, b)

        valid = hess > 0
        mean_h = hess[valid].mean() if valid.any() else 0.0
        if mean_h > 0:
            grad /= mean_h
            hess /= mean_h
        hess = np.maximum(hess, hess_floor)
        return grad.astype(np.float32), hess.astype(np.float32)

    return obj, structs, dropped


# =============================================================================
# Top-k ranking-quality metrics vs. a random baseline (verbatim from the source
# implementation's within_exp_lambdarank2.py + the evaluation pattern in its
# sibling analyze.py — see module docstring for what's ported vs. adapted).
# =============================================================================

def k_for_n(n, k_frac=0.10, k_min=5, k_max=50):
    return int(min(n, max(k_min, min(k_max, int(round(k_frac * n))))))


def ndcg_at_k_graded(rel, gain, scores, k):
    """Graded exponential-gain NDCG@k for one experiment (fragment)."""
    n = scores.size
    k = int(min(k, n))
    if k <= 0:
        return float("nan")
    order = np.argsort(-np.asarray(scores, dtype=np.float64), kind="stable")[:k]
    disc = 1.0 / np.log2(np.arange(2, k + 2))
    dcg = float(np.sum(gain[order] * disc))
    idcg = _ideal_dcg(gain, k)
    return 0.0 if idcg <= 0 else dcg / idcg


def hit_rate_at_k(rel, scores, k):
    """Fraction of the model's top-k rows that are true hits (rel == 3)."""
    rel = np.asarray(rel, dtype=np.int64)
    scores = np.asarray(scores, dtype=np.float64)
    k = int(min(k, len(scores)))
    if k <= 0:
        return float("nan")
    order = np.argsort(-scores, kind="stable")[:k]
    return float(np.mean(rel[order] == 3))


def enrichment_factor_at_k(rel, scores, k):
    """hit_rate@k / base_hit_rate for ONE experiment. NaN if no hits present.
    NOTE: for aggregate reporting, analyze.py deliberately uses a ratio of MEANS
    (mean hit_rate@k across experiments / mean random@k across experiments) rather
    than the mean of this per-experiment ratio, since the latter is unstable when
    base hit-rates vary a lot in size across experiments. This function is kept for
    single-experiment diagnostics; see `aggregate_topk_metrics` for the headline EF."""
    rel = np.asarray(rel, dtype=np.int64)
    base = float(np.mean(rel == 3))
    return float("nan") if base == 0 else hit_rate_at_k(rel, scores, k) / base


def pooled_hit_recovery_at_k(rel_score_pairs, k):
    """POOLED: sum(hits recovered in top_k) / sum(total hits available), across a
    list of (rel_array, scores_array) per admissible experiment(-fold)."""
    total_hits = total_recovered = 0
    for rel, ss in rel_score_pairs:
        rel = np.asarray(rel, dtype=np.int64)
        ss = np.asarray(ss, dtype=np.float64)
        n_hits = int(np.sum(rel == 3))
        if n_hits == 0:
            continue
        order = np.argsort(-ss, kind="stable")[:k]
        total_hits += n_hits
        total_recovered += int(np.sum(rel[order] == 3))
    return total_recovered / total_hits if total_hits > 0 else float("nan")


def random_baseline_ndcg(gain, rel, ks_by_name, hit_ks=(), n_perm=1000, seed=0):
    """Expected NDCG and hit_rate of a RANDOM ordering (labels/gains fixed) — the
    random-baseline floor every top-k metric here is compared against.
    ks_by_name: {name: k} for NDCG metrics. hit_ks: k values for hit_rate@k.
    Returns {name: (mean, std)}."""
    gain = np.asarray(gain, dtype=np.float64)
    rel = np.asarray(rel, dtype=np.int64)
    n = gain.size
    rng = np.random.default_rng(seed)
    disc = 1.0 / np.log2(np.arange(2, n + 2))
    idcg = {name: _ideal_dcg(gain, min(int(k), n)) for name, k in ks_by_name.items()}

    ndcg_acc = {name: np.empty(n_perm) for name in ks_by_name}
    hit_acc = {k: np.empty(n_perm) for k in hit_ks}

    for p in range(n_perm):
        perm = rng.permutation(n)
        cdcg = np.cumsum(gain[perm] * disc)
        for name, k in ks_by_name.items():
            kk = min(int(k), n)
            ndcg_acc[name][p] = 0.0 if idcg[name] <= 0 else cdcg[kk - 1] / idcg[name]
        for k in hit_ks:
            kk = min(int(k), n)
            hit_acc[k][p] = float(np.mean(rel[perm[:kk]] == 3))

    result = {name: (float(v.mean()), float(v.std())) for name, v in ndcg_acc.items()}
    for k in hit_ks:
        result[f"hit_rate@{k}"] = (float(hit_acc[k].mean()), float(hit_acc[k].std()))
    return result


def gain_weighted_pair_accuracy(precomputed_rel, scores, source_ids, min_n=8):
    """Pooled within-source pairwise accuracy weighted by |gain_i - gain_j|. A
    top-focused (hit-vs-nonhit pairs weigh 2^3-1=7, rel1-vs-rel0 weighs 1,
    rel0-vs-rel0 weighs 0) but smooth diagnostic — useful as a training-progress /
    selection signal since it pools over many hit-involving pairs rather than
    hinging on the exact position of 1-2 hit rows the way NDCG@k does."""
    rel = np.asarray(precomputed_rel, dtype=np.int64)
    scores = np.asarray(scores, dtype=np.float64)
    gain = 2.0 ** rel - 1.0
    num = den = 0.0
    for idx in group_indices(source_ids, min_size=min_n).values():
        gg, ss = gain[idx], scores[idx]
        n = idx.size
        ii, jj = np.triu_indices(n, k=1)
        w = np.abs(gg[ii] - gg[jj])
        m = w > 0
        if not np.any(m):
            continue
        sd = ss[ii] - ss[jj]
        gd = gg[ii] - gg[jj]
        correct = (np.sign(sd[m]) == np.sign(gd[m])).astype(np.float64)
        num += float(np.sum(w[m] * correct))
        den += float(np.sum(w[m]))
    return num / den if den > 0 else float("nan")


# =============================================================================
# Per-(experiment, fold) evaluation + aggregation (adapted from the sibling
# project's analyze.py: per_experiment_metrics / aggregate). Unit of analysis is
# an (experiment, fold) FRAGMENT, not a whole experiment — LOSO folds give each
# experiment intact (primary home for hit metrics); random/scaffold-disjoint
# k-fold splits fragment an experiment's rows across multiple folds' test sets
# (secondary, partial admissibility) — both are reported honestly, never conflated.
# =============================================================================

NDCG_METRIC_NAMES = ["ndcg@k_e", "ndcg_full", "ndcg@3", "ndcg@5", "ndcg@10"]
HIT_KS = (3, 5, 10)


def per_experiment_fold_metrics(y, scores, rel, source_id, fold_id, min_n=3,
                                 k_frac=0.10, k_min=5, k_max=50, n_perm=1000, seed=0):
    """One (experiment, fold) fragment's NDCG family + hit-rate family, each paired
    with its random-ordering baseline. Returns a dict; metrics are NaN (not
    silently omitted) when n < min_n or fewer than 2 relevance levels are present."""
    y = np.asarray(y, dtype=np.float64)
    scores = np.asarray(scores, dtype=np.float64)
    rel = np.asarray(rel, dtype=np.int64)
    n = rel.size
    gain = 2.0 ** rel - 1.0
    n_levels = int(np.unique(rel).size)
    has_hits = bool(np.any(rel == 3))
    row = {"source": source_id, "fold": fold_id, "n_vals": n, "n_rel_levels": n_levels, "has_hits": has_hits}

    if n >= min_n and n_levels >= 2:
        ks = {"ndcg@k_e": k_for_n(n, k_frac, k_min, k_max), "ndcg_full": n, "ndcg@3": 3, "ndcg@5": 5, "ndcg@10": 10}
        for name, k in ks.items():
            row[name] = ndcg_at_k_graded(rel, gain, scores, k)
        hit_ks_for_rand = HIT_KS if has_hits else ()
        base = random_baseline_ndcg(gain, rel, ks, hit_ks=hit_ks_for_rand, n_perm=n_perm, seed=seed)
        for name in ks:
            row[name + "_rand"] = base[name][0]
    else:
        for name in NDCG_METRIC_NAMES:
            row[name] = float("nan")
            row[name + "_rand"] = float("nan")
        base = None

    if n >= min_n and n_levels >= 2 and has_hits:
        for k in HIT_KS:
            row[f"hit_rate@{k}"] = hit_rate_at_k(rel, scores, k)
            row[f"hit_rate@{k}_rand"] = base[f"hit_rate@{k}"][0]
    else:
        for k in HIT_KS:
            row[f"hit_rate@{k}"] = float("nan")
            row[f"hit_rate@{k}_rand"] = float("nan")

    return row


def aggregate_topk_metrics(per_exp_fold_df):
    """Equal-weight / sqrt(n)-weighted / n-weighted means over (experiment, fold)
    rows with a defined metric, plus EF@k as a RATIO OF MEANS (not mean of
    per-fragment ratios — see enrichment_factor_at_k docstring for why)."""
    out = []
    n = per_exp_fold_df["n_vals"].to_numpy(dtype=np.float64)
    metric_cols = NDCG_METRIC_NAMES + [m + "_rand" for m in NDCG_METRIC_NAMES] + \
        [f"hit_rate@{k}" for k in HIT_KS] + [f"hit_rate@{k}_rand" for k in HIT_KS]
    for col in metric_cols:
        if col not in per_exp_fold_df.columns:
            continue
        vals = per_exp_fold_df[col].to_numpy(dtype=np.float64)
        finite = np.isfinite(vals)
        if finite.sum() == 0:
            out.append({"metric": col, "mean": np.nan, "sqrt_n_mean": np.nan, "n_weighted_mean": np.nan, "n_frag": 0})
            continue
        w_sqrt = np.sqrt(n[finite])
        w_n = n[finite]
        out.append({
            "metric": col,
            "mean": float(np.mean(vals[finite])),
            "sqrt_n_mean": float(np.sum(vals[finite] * w_sqrt) / np.sum(w_sqrt)),
            "n_weighted_mean": float(np.sum(vals[finite] * w_n) / np.sum(w_n)),
            "n_frag": int(finite.sum()),
        })
    agg_df = pd.DataFrame(out).set_index("metric")

    ef_rows = []
    for k in HIT_KS:
        hr_col, rand_col = f"hit_rate@{k}", f"hit_rate@{k}_rand"
        if hr_col not in agg_df.index or rand_col not in agg_df.index:
            continue
        mean_hr, mean_rand = agg_df.loc[hr_col, "mean"], agg_df.loc[rand_col, "mean"]
        ef = mean_hr / mean_rand if (mean_rand and mean_rand > 0) else float("nan")
        ef_rows.append({"metric": f"ef@{k}", "mean": ef, "sqrt_n_mean": np.nan, "n_weighted_mean": np.nan,
                        "n_frag": int(agg_df.loc[hr_col, "n_frag"])})
    if ef_rows:
        agg_df = pd.concat([agg_df, pd.DataFrame(ef_rows).set_index("metric")])
    return agg_df.reset_index()


# =============================================================================
# Self-test (5 correctness checks, ported from the source implementation)
# =============================================================================

def _self_test(seed=0):
    print("Running lambdarank.py self-test...")
    rng = np.random.default_rng(seed)

    y = np.array([-2.0, -1.0, 0.1, 0.6, 1.5, 2.5])
    rel = np.array([0, 0, 1, 2, 3, 3], dtype=np.int64)
    src = np.array(["E"] * 6)
    w = np.array([1.0, 1.3, 0.7, 1.1, 0.9, 1.2])
    beta = 1.0

    structs, dropped = compute_experiment_relevance(rel, src, top_rel_threshold=3, min_size=2)
    assert "E" in structs
    struct = structs["E"]
    print(f"  rel={struct['rel'].tolist()} gain={struct['gain'].tolist()} idcg={struct['idcg']:.4f} "
          f"top_local={struct['top_local'].tolist()}")

    s0 = rng.normal(size=6)
    ii, jj = np.triu_indices(6, 1)
    keep = struct["rel"][ii] != struct["rel"][jj]
    a_fix, b_fix = ii[keep], jj[keep]

    def grad_hess_L(s):
        g, h = np.zeros(6), np.zeros(6)
        L = _accumulate_experiment(g, h, s, struct, w, beta, a_fix.copy(), b_fix.copy())
        return g, h, L

    g0, h0, _ = grad_hess_L(s0)

    _, disc0 = _positions_discount(s0)
    rel_s, gain_s, idcg_s = struct["rel"], struct["gain"], struct["idcg"]
    a_is_k = rel_s[a_fix] > rel_s[b_fix]
    kk = np.where(a_is_k, a_fix, b_fix)
    jjk = np.where(a_is_k, b_fix, a_fix)
    dZ0 = (gain_s[kk] - gain_s[jjk]) * np.abs(disc0[kk] - disc0[jjk]) / idcg_s
    wp0 = w[kk] * w[jjk] * dZ0

    def L_frozen(s):
        ds = s[kk] - s[jjk]
        P = expit(beta * ds)
        return float(np.sum(wp0 * -np.log(np.clip(P, 1e-12, 1.0))))

    def grad_frozen(s):
        g = np.zeros(6)
        ds = s[kk] - s[jjk]
        P = expit(beta * ds)
        gk = wp0 * beta * (P - 1.0)
        np.add.at(g, kk, gk); np.add.at(g, jjk, -gk)
        return g

    eps = 1e-5
    fd = np.zeros(6)
    for r in range(6):
        sp, sm = s0.copy(), s0.copy()
        sp[r] += eps; sm[r] -= eps
        fd[r] = (L_frozen(sp) - L_frozen(sm)) / (2 * eps)
    gf = grad_frozen(s0)
    assert np.allclose(fd, gf, atol=1e-5), f"FD mismatch: fd={fd} grad={gf}"
    print(f"  (1) finite-difference == analytic grad: max|diff|={np.max(np.abs(fd-gf)):.2e} PASS")

    obj, _, _ = make_lambdarank_objective(src, rel, beta=beta, budget_B=50)

    class _DM:
        def __init__(self, w): self._w = w
        def get_weight(self): return self._w

    g_s, h_s = obj(s0.astype(np.float32), _DM(w))
    assert abs(float(g_s.sum())) < 1e-5, f"grad sum not ~0: {g_s.sum()}"
    print(f"  (2) grad sum ~0 (sampling on): {float(g_s.sum()):.2e} PASS")

    c = 2.7
    gA, hA, LA = grad_hess_L(s0)
    gB, hB, LB = grad_hess_L(s0 + c)
    assert np.allclose(gA, gB, atol=1e-8) and np.allclose(hA, hB, atol=1e-8) and abs(LA - LB) < 1e-8
    print(f"  (3) gauge invariant +{c}: max|dgrad|={np.max(np.abs(gA-gB)):.2e} PASS")

    assert (hA >= 0).all() and (h_s >= 0).all()
    print(f"  (4) hessian >= 0: min={float(min(hA.min(), h_s.min())):.2e} PASS")

    _, disc = _positions_discount(s0)

    def dz(i, j):
        return abs(gain_s[i] - gain_s[j]) * abs(disc[i] - disc[j]) / idcg_s

    floor = np.where(rel_s == 0)[0]
    top = np.where(rel_s == 3)[0]
    ff = dz(floor[0], floor[1]) if floor.size >= 2 else 0.0
    ft = dz(floor[0], top[0]) if floor.size and top.size else float("nan")
    assert ff == 0.0
    assert ft > 0
    print(f"  (5) dZ floor-floor={ff:.4f}(==0) floor-top={ft:.4f} PASS")

    # Bonus: relevance-grade construction sanity
    rank_test = np.linspace(0, 1, 40)
    src_test = np.array(["S"] * 40)
    rel_test = compute_relevance_grades(rank_test, src_test, hit_frac=0.10)
    assert (rel_test == 3).sum() == 4, "expected 10% of 40 = 4 hits"
    assert np.all(np.diff(rel_test[np.argsort(rank_test)]) >= 0), "rel must be monotone in rank"
    print(f"  (6) relevance-grade construction: {np.bincount(rel_test)} (monotone in rank) PASS")

    # (7) Top-k metrics: perfect ranker should crush its random baseline
    # struct has exactly 2 rel==3 rows (indices 4,5), so hit_rate@2 (k == n_hits) is the
    # clean "should be 1.0" check; hit_rate@3 necessarily dilutes to 2/3 since k > n_hits.
    perfect_scores = gain_s.astype(np.float64)  # gain-order == perfect ranking
    hr2 = hit_rate_at_k(rel_s, perfect_scores, 2)
    hr3 = hit_rate_at_k(rel_s, perfect_scores, 3)
    ef3 = enrichment_factor_at_k(rel_s, perfect_scores, 3)
    assert hr2 == 1.0, f"perfect ranker should have hit_rate@2=1.0 (2 hits, k=2), got {hr2}"
    assert abs(hr3 - 2 / 3) < 1e-9, f"perfect ranker hit_rate@3 should be 2/3 (2 hits, k=3), got {hr3}"
    assert ef3 > 1.0, f"perfect ranker EF@3 should exceed 1x chance, got {ef3}"
    ndcg_perfect = ndcg_at_k_graded(rel_s, gain_s, perfect_scores, 6)
    assert abs(ndcg_perfect - 1.0) < 1e-9, f"perfect ranker full NDCG should be 1.0, got {ndcg_perfect}"
    worst_scores = -gain_s.astype(np.float64)  # exactly inverted ranking
    ndcg_worst = ndcg_at_k_graded(rel_s, gain_s, worst_scores, 6)
    assert ndcg_worst < ndcg_perfect, "inverted ranking should score below perfect ranking"
    base = random_baseline_ndcg(gain_s, rel_s, {"ndcg_full": 6}, hit_ks=(3,), n_perm=2000, seed=1)
    assert base["ndcg_full"][0] < ndcg_perfect, "random-order NDCG should be below perfect NDCG"
    assert 0.0 <= base["hit_rate@3"][0] <= 1.0
    print(f"  (7) top-k vs random baseline: hit_rate@2(perfect)={hr2:.3f} hit_rate@3(perfect)={hr3:.3f} EF@3={ef3:.2f}x  "
          f"NDCG(perfect)={ndcg_perfect:.3f} NDCG(worst)={ndcg_worst:.3f} NDCG(random floor)={base['ndcg_full'][0]:.3f} "
          f"hit_rate@3(random floor)={base['hit_rate@3'][0]:.3f}  PASS")

    # (8) pooled hit recovery + gain-weighted pairwise accuracy sanity
    pooled = pooled_hit_recovery_at_k([(rel_s, perfect_scores)], k=3)
    assert pooled == 1.0, f"perfect ranker pooled hit recovery should be 1.0, got {pooled}"
    gw = gain_weighted_pair_accuracy(rel_s, perfect_scores, np.array(["E"] * 6), min_n=2)
    assert gw == 1.0, f"perfect ranker gain-weighted pairwise accuracy should be 1.0, got {gw}"
    print(f"  (8) pooled_hit_recovery@3(perfect)={pooled:.3f}  gain_weighted_pair_acc(perfect)={gw:.3f}  PASS")

    print("lambdarank.py self-test PASSED.\n")


if __name__ == "__main__":
    _self_test()
