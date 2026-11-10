#!/usr/bin/env python3
"""Stage 09 -- make every check fail on purpose. Writes 09_seedgate.json.

    python3 -I 09_seeds.py            # parent: runs one CHILD PROCESS per seed
    python3 -I 09_seeds.py <name>     # one child

A check that cannot fail is not a check. A seed that only raises is not credited:
each child must PRINT its metric on a line of the form

    METRIC <check> <value>

and the parent refuses to credit any seed whose metric line is absent, however
the child exited -- that is the failure shape where a child dies in a TypeError
before reaching the thing under test. The child exits NON-ZERO when its check
fires; the parent requires, for each seed: control exit 0 with a metric, seeded
exit non-zero with a metric, and the two metrics DIFFERENT.
"""
import copy, json, os, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def load_checks():
    """06_checks.py is not an importable module name; load it by path."""
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "stage06", os.path.join(HERE, "06_checks.py"))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def data():
    with open(os.path.join(HERE, "04_si_tables.json")) as f:
        return json.load(f)


# ---------------------------------------------------------------- the seeds
def seed_none(d):
    return d


def seed_norm_raw(d):
    """N1: move ONE raw barcode count by +10. Subtle: the count stays plausible,
    the column total is recomputed from the counts, and only the reconstruction
    residual against the paper's printed numbers can see it."""
    t = d["normalisation_tables_fig_S1E"]
    t[0][2][2] = str(int(t[0][2][2]) + 10)        # LNP 1, Mouse 1: 200 -> 210
    return d


def seed_molpct(d):
    """M1: one library-1 cholesterol value off by 2.5, so its row sums to 97.5."""
    d["library_1"]["7"]["chol_molpct"] = 20.0      # was 22.5
    return d


def seed_factorial(d):
    """F1: retype ONE LNP's PEG, which breaks 'every variant carries the same
    ratio set' without changing the row count."""
    d["library_1"]["7"]["peg_type"] = "C16PEG2000"
    return d


def seed_barcode(d):
    """B1: drop one barcode from the count."""
    d["barcodes"] = dict(d["barcodes"])
    d["barcodes"]["n_barcodes"] = 155
    d["barcodes"]["n_distinct_sequences"] = 155
    d["barcodes"]["matches_claim"] = False
    return d


def seed_winner_gross(d):
    """W1 gross: shift every ionizable mole % by +0.5 so BM1's composition
    matches nothing. Credited by the located-count collapsing, not by a raise."""
    return d          # applied via the perturb argument, see run_child


def seed_pie(d):
    """E1: move the 60/22.5/17.5/0 ratio to 60/17.5/22.5/0 in ALL FIVE PEG blocks.

    Deliberately specific. Each row still sums to 100 (cholesterol absorbs it), so
    M1 cannot see it; the change is identical across the five blocks, so the
    factorial F1 still holds; the row count and every PEG type are untouched. The
    ONLY thing that moves is the Fig 4A band census -- five LNPs cross from
    [15,20) into >=20 -- so this seed can only be caught by E1.

    An earlier version of this seed moved ONE row from 20 to 22.5 and was a NO-OP,
    because the half-open convention already puts 20 in the '>20' band. It is
    recorded here because a seed that moves nothing is the failure this gate is
    for: it exited non-zero on the factorial while E1's metric never moved, which
    is why 'credited' requires the metric to move and not merely a non-zero exit.
    """
    n = 0
    for k, f in d["library_1"].items():
        if (f["il_molpct"], f["chol_molpct"], f["peg_molpct"],
                f["dspc_molpct"]) == (60.0, 22.5, 17.5, 0.0):
            f["chol_molpct"], f["peg_molpct"] = 17.5, 22.5
            n += 1
    assert n == 5, f"expected the ratio in all 5 PEG blocks, moved {n}"
    return d


SEEDS = {
    "control": (seed_none, None),
    "E1_formulated_pie": (seed_pie, None),
    "N1_normalisation_raw_count": (seed_norm_raw, None),
    "M1_molpct_sum": (seed_molpct, None),
    "F1_factorial_peg_type": (seed_factorial, None),
    "B1_barcode_count": (seed_barcode, None),
    "W1_winner_location": (seed_winner_gross, 0.5),
}


# ---------------------------------------------------------------- the child
def run_child(name):
    m = load_checks()
    fn, perturb = SEEDS[name]
    d = fn(copy.deepcopy(data()))

    norm = m.normalisation(d["normalisation_tables_fig_S1E"])
    n_exact = norm["stage_3_normalised_delivery"]["n_exact_after_rounding"]
    # recompute mol% and factorial from the (possibly seeded) tables
    lib1 = {int(k): v for k, v in d["library_1"].items()}
    viol = sum(1 for f in lib1.values()
               if abs(f["il_molpct"] + f["chol_molpct"] + f["peg_molpct"]
                      + f["dspc_molpct"] - 100.0) > 1e-9)
    by_peg = {}
    for f in lib1.values():
        by_peg.setdefault(f["peg_type"], []).append(
            (f["il_molpct"], f["chol_molpct"], f["peg_molpct"], f["dspc_molpct"]))
    first = sorted(by_peg[sorted(by_peg)[0]])
    fact_ok = (len(by_peg) == 5 and len(first) == 24
               and all(sorted(v) == first for v in by_peg.values()))
    nbc = d["barcodes"]["n_barcodes"]
    w = m.winners(d, perturb=perturb)
    pie = m.pies(d)
    n_pies_ok = pie["_summary"]["n_reconciling"]

    print(f"METRIC N1_stage3_exact_cells {n_exact}")
    print(f"METRIC M1_molpct_violations_lib1 {viol}")
    print(f"METRIC F1_factorial_holds_lib1 {int(fact_ok)}")
    print(f"METRIC B1_n_barcodes {nbc}")
    print(f"METRIC W1_uniquely_located {w['n_uniquely_located']}")
    print(f"METRIC E1_pies_reconciling {n_pies_ok}")

    fired = (n_exact != 9) or (viol != 0) or (not fact_ok) or (nbc != 156) \
        or (w["n_uniquely_located"] != 3) or (n_pies_ok != 5)
    return 1 if fired else 0


# --------------------------------------------------------------- the parent
def parse(out):
    d = {}
    for line in out.splitlines():
        if line.startswith("METRIC "):
            _, k, v = line.split()
            d[k] = v
    return d


METRIC_FOR = {
    "N1_normalisation_raw_count": "N1_stage3_exact_cells",
    "M1_molpct_sum": "M1_molpct_violations_lib1",
    "F1_factorial_peg_type": "F1_factorial_holds_lib1",
    "B1_barcode_count": "B1_n_barcodes",
    "W1_winner_location": "W1_uniquely_located",
    "E1_formulated_pie": "E1_pies_reconciling",
}


def main():
    if len(sys.argv) > 1:
        return run_child(sys.argv[1])

    runs = {}
    for name in SEEDS:
        r = subprocess.run([sys.executable, "-I", os.path.abspath(__file__), name],
                           capture_output=True, text=True, cwd=HERE)
        runs[name] = {"returncode": r.returncode, "metrics": parse(r.stdout),
                      "stderr_tail": r.stderr.strip().splitlines()[-1:]}

    ctrl = runs["control"]
    report = {
        "_what": "Seeded-mutation gate. One CHILD PROCESS per seed; the verdict is "
                 "the child's returncode, and a seed whose metric line is absent is "
                 "NOT credited however it exited.",
        "control": {"returncode": ctrl["returncode"], "metrics": ctrl["metrics"],
                    "passes_clean": ctrl["returncode"] == 0 and len(ctrl["metrics"]) == 6},
        "seeds": {},
    }
    for name, metric in METRIC_FOR.items():
        r = runs[name]
        has = metric in r["metrics"]
        moved = has and r["metrics"][metric] != ctrl["metrics"].get(metric)
        report["seeds"][name] = {
            "metric": metric,
            "control_value": ctrl["metrics"].get(metric),
            "seeded_value": r["metrics"].get(metric) if has else None,
            "metric_line_present": has,
            "metric_moved": moved,
            "returncode": r["returncode"],
            "credited": bool(has and moved and r["returncode"] != 0),
            "stderr_tail": r["stderr_tail"],
        }
    report["n_credited"] = sum(1 for v in report["seeds"].values() if v["credited"])
    report["n_seeds"] = len(report["seeds"])
    report["verdict"] = ("ALL SEEDS CREDITED" if report["n_credited"] == report["n_seeds"]
                         and report["control"]["passes_clean"] else "NOT CLEAN")
    with open(os.path.join(HERE, "09_seedgate.json"), "w") as f:
        json.dump(report, f, indent=1)
    print(json.dumps(report, indent=1))
    return 0 if report["verdict"] == "ALL SEEDS CREDITED" else 1


if __name__ == "__main__":
    sys.exit(main())
