"""
Single entry point: regenerates the entire Confound_analysis pipeline from
raw data (Stages 0-9) and writes environment/manifest.json.

Usage: python3 run_all.py [--from-stage N]
"""
import argparse
import json
import platform
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "scripts"))

STAGES = [
    ("stage00_inspect", "Data inventory"),
    ("stage01_curation", "Curation and frozen analysis set"),
    ("stage02_corpus", "Corpus overview"),
    ("stage03_separability", "Source separability"),
    ("stage04_exclusivity", "Source exclusivity and excipient invariance"),
    ("stage05_concordance", "Cross-lab concordance (centerpiece)"),
    ("stage06_protocol", "Protocol comparison grid"),
    ("stage07_learnable", "Within-experiment learnability and decay"),
    ("stage08_controls", "Negative controls and robustness"),
    ("stage09_summary", "Aggregate summary and manifest"),
]


def git_commit_hash():
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, stderr=subprocess.DEVNULL).decode().strip()
    except Exception:
        return None


def package_versions():
    versions = {}
    for pkg in ["pandas", "numpy", "rdkit", "sklearn", "scipy", "xgboost", "umap", "matplotlib", "pyarrow"]:
        try:
            mod = __import__(pkg)
            versions[pkg] = getattr(mod, "__version__", "unknown")
        except ImportError:
            versions[pkg] = "not installed"
    return versions


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--from-stage", type=int, default=0, help="Resume from this stage index (0-9)")
    args = parser.parse_args()

    from common import GLOBAL_SEED, ENVIRONMENT

    run_start = time.time()
    stage_timings = {}

    for i, (module_name, description) in enumerate(STAGES):
        if i < args.from_stage:
            continue
        print(f"\n{'='*70}\nStage {i}: {description} ({module_name})\n{'='*70}")
        t0 = time.time()
        mod = __import__(module_name)
        mod.main()
        stage_timings[module_name] = round(time.time() - t0, 1)

    manifest = {
        "timestamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "global_seed": GLOBAL_SEED,
        "git_commit": git_commit_hash(),
        "os": platform.platform(),
        "python_version": platform.python_version(),
        "package_versions": package_versions(),
        "per_stage_runtime_seconds": stage_timings,
        "total_runtime_seconds": round(time.time() - run_start, 1),
    }
    with open(ENVIRONMENT / "manifest.json", "w") as f:
        json.dump(manifest, f, indent=2)
    print(f"\nWrote environment/manifest.json")
    print(f"Total runtime: {manifest['total_runtime_seconds']}s")


if __name__ == "__main__":
    main()
