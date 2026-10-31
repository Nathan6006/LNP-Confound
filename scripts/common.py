"""
Shared utilities for the Confound_analysis pipeline: paths, seeds, the results
registry, featurization (Morgan fingerprints, RDKit descriptors, scaffolds,
Butina clustering), rank normalization, and bootstrap CI.

Every stage script imports from here rather than re-implementing these
primitives, so the rank-normalization and Tanimoto logic — which the whole
paper's centerpiece result rests on — exists in exactly one place.
"""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import rdkit
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.ML.Cluster import Butina

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
RESULTS = ROOT / "results"
REPORTS = ROOT / "reports"
FIGURES = ROOT / "figures"
CACHE = ROOT / "cache"
ENVIRONMENT = ROOT / "environment"
for d in (RESULTS, REPORTS, FIGURES, CACHE, ENVIRONMENT):
    d.mkdir(exist_ok=True)

GLOBAL_SEED = 42
MORGAN_RADIUS = 2
MORGAN_NBITS = 2048


# ============================== Data loading ==============================

def load_curated():
    return pd.read_parquet(RESULTS / "curated.parquet")


def load_curated_molea():
    return pd.read_parquet(RESULTS / "curated_molea.parquet")


def load_schema():
    with open(RESULTS / "schema.json") as f:
        return json.load(f)


# ============================== Results registry ==============================

_REGISTRY_PATH = RESULTS / "registry.json"


def registry_set(key, value, n=None, ci=None, confidence="OK", stage=None, note=None, extra=None):
    """Write one entry to the single machine-readable results registry.
    confidence: "OK" or "LOW-CONFIDENCE" (never silently fabricate a number —
    if n is too small to trust, mark it here instead of omitting the flag)."""
    reg = {}
    if _REGISTRY_PATH.exists():
        with open(_REGISTRY_PATH) as f:
            reg = json.load(f)

    def _jsonable(v):
        if isinstance(v, (np.floating,)):
            return float(v)
        if isinstance(v, (np.integer,)):
            return int(v)
        if isinstance(v, np.ndarray):
            return v.tolist()
        return v

    entry = {
        "value": _jsonable(value),
        "n": _jsonable(n),
        "ci": [_jsonable(c) for c in ci] if ci is not None else None,
        "confidence": confidence,
        "stage": stage,
        "note": note,
    }
    if extra:
        entry.update({k: _jsonable(v) for k, v in extra.items()})
    reg[key] = entry
    with open(_REGISTRY_PATH, "w") as f:
        json.dump(reg, f, indent=2)
    return entry


def registry_get_all():
    if not _REGISTRY_PATH.exists():
        return {}
    with open(_REGISTRY_PATH) as f:
        return json.load(f)


# ============================== Caching ==============================

def cache_key(*parts):
    s = "||".join(str(p) for p in parts)
    return hashlib.sha256(s.encode()).hexdigest()[:16]


def cache_path(name, *key_parts):
    return CACHE / f"{name}_{cache_key(*key_parts)}.npz"


# ============================== Featurization ==============================

_morgan_gen_cache = {}


def _get_morgan_gen(radius, n_bits):
    key = (radius, n_bits)
    if key not in _morgan_gen_cache:
        _morgan_gen_cache[key] = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=n_bits)
    return _morgan_gen_cache[key]


def smiles_to_mol(smiles):
    if pd.isna(smiles):
        return None
    return Chem.MolFromSmiles(str(smiles))


def morgan_fingerprint(smiles, radius=MORGAN_RADIUS, n_bits=MORGAN_NBITS):
    """Return an (n_bits,) uint8 numpy array, or None if the SMILES doesn't parse."""
    mol = smiles_to_mol(smiles)
    if mol is None:
        return None
    gen = _get_morgan_gen(radius, n_bits)
    fp = gen.GetFingerprint(mol)
    arr = np.zeros((n_bits,), dtype=np.uint8)
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr


def morgan_fingerprint_bitvect(smiles, radius=MORGAN_RADIUS, n_bits=MORGAN_NBITS):
    """Return the raw RDKit ExplicitBitVect (fast for BulkTanimotoSimilarity), or None."""
    mol = smiles_to_mol(smiles)
    if mol is None:
        return None
    gen = _get_morgan_gen(radius, n_bits)
    return gen.GetFingerprint(mol)


def batch_morgan_fingerprints(smiles_list, radius=MORGAN_RADIUS, n_bits=MORGAN_NBITS, cache_name="morgan"):
    """Vectorized + cached Morgan fingerprinting for a list of canonical SMILES.
    Returns (fp_matrix [n, n_bits] uint8, valid_mask [n] bool). Invalid/unparseable
    entries get an all-zero row and valid_mask=False at that position."""
    smiles_list = list(smiles_list)
    key = cache_path(cache_name, tuple(smiles_list), radius, n_bits)
    if key.exists():
        z = np.load(key)
        return z["fps"], z["valid"]

    n = len(smiles_list)
    fps = np.zeros((n, n_bits), dtype=np.uint8)
    valid = np.zeros((n,), dtype=bool)
    for i, smi in enumerate(smiles_list):
        arr = morgan_fingerprint(smi, radius, n_bits)
        if arr is not None:
            fps[i] = arr
            valid[i] = True
    np.savez_compressed(key, fps=fps, valid=valid)
    return fps, valid


def bemis_murcko_scaffold(smiles):
    mol = smiles_to_mol(smiles)
    if mol is None:
        return None
    try:
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
        return Chem.MolToSmiles(scaffold, canonical=True)
    except Exception:
        return None


def butina_clusters(bitvects, cutoff=0.4):
    """Butina cluster a list of RDKit bit vectors. Returns a list of cluster-id
    ints, one per input (same order), and the list of cluster tuples RDKit returns.
    cutoff = Tanimoto DISTANCE threshold (1 - similarity) for two points to be
    in the same cluster; RDKit's Butina.ClusterData default convention."""
    n = len(bitvects)
    dists = []
    for i in range(1, n):
        sims = DataStructs.BulkTanimotoSimilarity(bitvects[i], bitvects[:i])
        dists.extend([1 - s for s in sims])
    clusters = Butina.ClusterData(dists, n, cutoff, isDistData=True)
    cluster_id = np.full(n, -1, dtype=int)
    for cid, members in enumerate(clusters):
        for m in members:
            cluster_id[m] = cid
    return cluster_id, clusters


# ============================== Rank normalization ==============================

def normalized_rank(values):
    """Convert a 1-D array-like of readout values to normalized ranks in [0, 1],
    average rank for ties. A constant-valued input (n==1, or all-equal values)
    returns all 0.5. This is the ONLY transform ever applied before comparing a
    delivery readout across two different sources — never compare raw
    Experiment_value across sources."""
    s = pd.Series(values, dtype=float)
    n = s.notna().sum()
    if n <= 1:
        return pd.Series(np.where(s.notna(), 0.5, np.nan), index=s.index)
    ranks = s.rank(method="average", na_option="keep")
    return (ranks - 1) / (n - 1)


# ============================== Bootstrap CI ==============================

def bootstrap_ci(data, stat_fn, n_boot=2000, seed=GLOBAL_SEED, alpha=0.05):
    """Percentile bootstrap CI for stat_fn applied to `data` (a sequence, or a
    tuple of parallel sequences resampled jointly by index — pass a tuple of
    same-length arrays and stat_fn(*arrays))."""
    rng = np.random.default_rng(seed)
    if isinstance(data, tuple):
        n = len(data[0])
        arrays = data
    else:
        n = len(data)
        arrays = (np.asarray(data),)
    if n == 0:
        return (np.nan, np.nan)
    stats = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        resampled = tuple(a[idx] if not isinstance(a, list) else [a[i] for i in idx] for a in arrays)
        try:
            stats.append(stat_fn(*resampled))
        except Exception:
            continue
    if not stats:
        return (np.nan, np.nan)
    lo = np.nanpercentile(stats, 100 * alpha / 2)
    hi = np.nanpercentile(stats, 100 * (1 - alpha / 2))
    return (float(lo), float(hi))


# ============================== Unit sanity checks ==============================

def _sanity_checks():
    # normalized_rank: monotone, spans [0,1], ties average
    r = normalized_rank([10, 20, 30, 40])
    assert np.allclose(r.values, [0.0, 1 / 3, 2 / 3, 1.0]), r.values
    r_tied = normalized_rank([5, 5, 10])
    assert np.isclose(r_tied.iloc[0], r_tied.iloc[1]), "ties should get equal rank"
    assert np.isclose(r_tied.iloc[2], 1.0)
    r_const = normalized_rank([7, 7, 7])
    assert np.allclose(r_const.values, 0.5), "constant input should map to 0.5, not raise"
    r_one = normalized_rank([42])
    assert np.allclose(r_one.values, 0.5)

    # Tanimoto sanity: identical molecule -> similarity 1; very different -> < 1
    fp1 = morgan_fingerprint_bitvect("CCCCCCCCCC")  # decane
    fp2 = morgan_fingerprint_bitvect("CCCCCCCCCC")
    fp3 = morgan_fingerprint_bitvect("c1ccccc1O")  # phenol
    sim_same = DataStructs.TanimotoSimilarity(fp1, fp2)
    sim_diff = DataStructs.TanimotoSimilarity(fp1, fp3)
    assert np.isclose(sim_same, 1.0), sim_same
    assert sim_diff < 1.0, sim_diff

    # canonical SMILES round-trip: two equivalent SMILES strings for the same
    # molecule must canonicalize to the same string
    a = Chem.MolToSmiles(Chem.MolFromSmiles("CCO"), canonical=True)
    b = Chem.MolToSmiles(Chem.MolFromSmiles("OCC"), canonical=True)
    assert a == b, (a, b)

    print("All common.py sanity checks passed.")


if __name__ == "__main__":
    _sanity_checks()
    print(f"rdkit version: {rdkit.__version__}")
