#!/usr/bin/env python3
"""Seed harness -- make every check in this folder fail ON PURPOSE.

    python3 -I seed_checks.py

F-069: THE HARNESS WRITES ONLY TO A TEMP DIRECTORY. Each seed runs in its own
child process with its own OUTDIR under tempfile.mkdtemp(); nothing in the
curation folder is opened for writing by a seeded run. The folder's tracked
artifacts are hashed before and after the whole suite and the hashes must be
identical -- because a seeded run that exits non-zero can still have written,
and the non-zero exit conceals it.

A SEED THAT ONLY RAISES IS NOT CREDITED. Each seed must print a line
`METRIC <name> <value>` for the metric it is supposed to move, both clean and
seeded. If the metric line is absent from either run the seed is recorded as
NOT CREDITED, whatever the exit code -- that is how a seed that dies in a
TypeError before reaching its own metric gets caught.

Verdicts come from each child's returncode and stdout, one child per seed.
Never `python3 ... | tail`: that returns tail's exit code.
"""
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))

TRACKED = ["00_value_sources.json", "01_blocks.csv", "01_si_tables.json",
           "02_corpus_overlap.json", "03_drawn_checks.json",
           "03_enumerated.csv", "05_values.csv", "06_rows.csv",
           "07_qc.json", "08_contradictions.csv",
           "reference_structures.csv", "PROVENANCE.md", "DECISIONS.md"]

# Each seed: (name, metric_name, mutation source, expectation)
# The child receives MUT via env and applies it to the in-memory objects AFTER
# extraction, so what is being tested is the CHECK, not the extractor.
CHILD = r'''
import json, os, sys, csv
sys.path.insert(0, os.environ["FOLDER"])
import pdfplumber
import stage01_si_tables as s1
import stage03_drawn_checks as s3

MUT = os.environ.get("MUT", "none")
METRIC = os.environ["METRIC"]

with pdfplumber.open(s1.SI) as pdf:
    l1 = s1.read_library1(pdf)
    l2 = s1.read_library2(pdf)
    ex = s1.read_norm_example(pdf)
band = s3.read_structure_band()

# ---- mutations, applied to extracted objects -------------------------------
if MUT == "molpct_perturbed":
    l1[0]["Chol %"] = str(float(l1[0]["Chol %"]) + 1.0)
elif MUT == "molpct_unparseable":
    l1[5]["PEG %"] = "n/a"
elif MUT == "lipomer_type_overwritten":
    for r in l2[:40]:
        r["Lipomer Type"] = "MC3"
elif MUT == "diameter_column_blanked":
    for r in l1:
        r["Diameter (nm)"] = ""
elif MUT == "printed_pct_perturbed":
    ex["pct"]["rows"]["1"]["lung"][0] += 5
elif MUT == "printed_total_inconsistent":
    ex["raw"]["total"]["lung"][0] += 25
elif MUT == "stage_B_skips_input_normalisation":
    for k in ex["inp"]["rows"]:
        ex["inp"]["rows"][k]["lung"] = list(ex["pct"]["rows"][k]["lung"])
        ex["inp"]["rows"][k]["kidney"] = list(ex["pct"]["rows"][k]["kidney"])
elif MUT == "peg_repeat_shifted":
    band["label_to_repeat_pairs"][0]["repeat"] = 60
    band["polyether_repeat_counts_in_drawing_order"][0] = 60
elif MUT == "peg_repeat_mw_classes_merged":
    for p in band["label_to_repeat_pairs"]:
        p["repeat"] = 45
elif MUT == "drawn_peg_label_dropped":
    band["drawn_peg_lipid_labels"] = band["drawn_peg_lipid_labels"][:3]
    band["label_to_repeat_pairs"] = band["label_to_repeat_pairs"][:3]
elif MUT != "none":
    print("MUTATION_UNKNOWN", MUT); sys.exit(3)

# ---- recompute the metric under test --------------------------------------
val = None
if METRIC == "molpct_lib1_ok":
    val = s1.check_molpct_sums(l1, l2)["library_1"]["summed_to_100"]
elif METRIC == "molpct_lib1_unparsed":
    val = s1.check_molpct_sums(l1, l2)["library_1"]["n_unparsed"]
elif METRIC == "n_distinct_ionizable":
    val = s1.check_single_ionizable_lipid(l1, l2)["n_distinct_ionizable_lipids"]
elif METRIC == "diam_positive_control":
    c = s1.check_single_ionizable_lipid(l1, l2)
    val = int(bool(c["positive_control_diameter_varies"]))
elif METRIC == "normA_exact":
    val = s1.check_normalisation(ex)["stage_A_counts_per_100_summary"]["n_exact"]
elif METRIC == "normA_total_guard_raises":
    # The printed-total guard is a REFUSAL, not a number, so it is credited by
    # observing the refusal -- reported as a metric line so that a seed dying
    # for some unrelated reason still cannot be credited.
    try:
        s1.check_normalisation(ex)
        val = 0
    except SystemExit:
        val = 1
elif METRIC == "normB_exact":
    val = s1.check_normalisation(ex)["stage_B_input_normalised_summary"]["n_exact"]
elif METRIC == "peg_all_within_one":
    val = int(bool(s3.check_peg_repeats(band)["all_within_one_repeat"]))
elif METRIC == "peg_ratio_residual_pct":
    r = s3.check_peg_repeats(band).get("ratio_test")
    val = None if r is None else round(r["relative_residual_pct"], 4)
elif METRIC == "peg_n_drawn":
    val = s3.check_peg_repeats(band)["n_lipids_drawn"]
else:
    print("METRIC_UNKNOWN", METRIC); sys.exit(4)

print("METRIC", METRIC, val)
'''

SEEDS = [
    ("molpct_perturbed", "molpct_lib1_ok",
     "one Chol % cell +1.0 -> that formulation no longer sums to 100",
     "decrease"),
    ("molpct_unparseable", "molpct_lib1_unparsed",
     "one PEG % cell set to 'n/a' -> must be counted as unparsed, not silently skipped",
     "increase"),
    ("lipomer_type_overwritten", "n_distinct_ionizable",
     "40 library-2 rows relabelled 'MC3' -> the one-ionizable-lipid finding must break",
     "increase"),
    ("diameter_column_blanked", "diam_positive_control",
     "Diameter column emptied -> the sweep's own positive control must stop passing",
     "decrease"),
    ("printed_pct_perturbed", "normA_exact",
     "one PRINTED percentage +5 -> stage A can no longer reproduce it. The "
     "first version of this seed inflated a RAW count instead, which trips "
     "check_normalisation's printed-total guard and raises before the metric "
     "is reached; the guard is correct, so the seed was retargeted rather "
     "than the guard widened, and the guard is credited by its own seed below",
     "decrease"),
    ("printed_total_inconsistent", "normA_total_guard_raises",
     "a printed column total +25 -> the guard that cross-checks the printed "
     "total against the sum of the printed counts must refuse",
     "increase"),
    ("stage_B_skips_input_normalisation", "normB_exact",
     "stage-B printed values replaced by stage-A values -> input normalisation "
     "is what distinguishes them, so the reconstruction must stop matching",
     "decrease"),
    ("peg_repeat_shifted", "peg_all_within_one",
     "one drawn repeat 45 -> 60 -> must fall outside one repeat of 2000/44.05",
     "decrease"),
    ("peg_repeat_mw_classes_merged", "peg_ratio_residual_pct",
     "PEG3000's repeat forced to 45 -> the ratio test must move off 0.74%",
     "increase"),
    ("drawn_peg_label_dropped", "peg_n_drawn",
     "one drawn PEG-lipid removed -> the 4-drawn count must drop",
     "decrease"),
]


def sha(path):
    if not os.path.exists(path):
        return None
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


def run(child_path, metric, mut, outdir):
    env = dict(os.environ, FOLDER=HERE, METRIC=metric, MUT=mut,
               OUTDIR=outdir, PYTHONDONTWRITEBYTECODE="1")
    p = subprocess.run([sys.executable, "-I", child_path],
                       capture_output=True, text=True, env=env, cwd=outdir)
    val, seen = None, False
    for line in p.stdout.splitlines():
        if line.startswith("METRIC "):
            parts = line.split(None, 2)
            if len(parts) == 3 and parts[1] == metric:
                seen = True
                try:
                    val = float(parts[2])
                except ValueError:
                    val = parts[2]
    return {"returncode": p.returncode, "metric_line_present": seen,
            "value": val, "stderr_tail": p.stderr.strip()[-300:]}


def main():
    before = {f: sha(os.path.join(HERE, f)) for f in TRACKED}
    tmp = tempfile.mkdtemp(prefix="seed_30216729_")
    child_path = os.path.join(tmp, "seed_child.py")
    with open(child_path, "w") as f:
        f.write(CHILD)

    results = []
    for name, metric, what, direction in SEEDS:
        od_clean = os.path.join(tmp, name + "_clean")
        od_seed = os.path.join(tmp, name + "_seed")
        os.makedirs(od_clean, exist_ok=True)
        os.makedirs(od_seed, exist_ok=True)
        clean = run(child_path, metric, "none", od_clean)
        seeded = run(child_path, metric, name, od_seed)

        credited, why = False, ""
        if not clean["metric_line_present"]:
            why = ("NOT CREDITED: the clean run printed no metric line, so "
                   "there is no baseline to compare against")
        elif not seeded["metric_line_present"]:
            why = ("NOT CREDITED: the seeded run printed no metric line -- it "
                   "died before reaching the metric, so the seed demonstrates "
                   "nothing about the check")
        elif clean["value"] == seeded["value"]:
            why = (f"NOT CREDITED: metric did not move "
                   f"({clean['value']} -> {seeded['value']}); the check is "
                   f"blind to this defect")
        else:
            moved = ("increase" if seeded["value"] > clean["value"]
                     else "decrease")
            if moved != direction:
                why = (f"NOT CREDITED: metric moved the wrong way "
                       f"({clean['value']} -> {seeded['value']}, expected "
                       f"{direction})")
            else:
                credited = True
                why = (f"CREDITED: {clean['value']} -> {seeded['value']} "
                       f"({direction})")
        results.append({"seed": name, "metric": metric, "what": what,
                        "expected_direction": direction,
                        "clean": clean, "seeded": seeded,
                        "credited": credited, "verdict": why})

    after = {f: sha(os.path.join(HERE, f)) for f in TRACKED}
    untouched = all(before[f] == after[f] for f in TRACKED)

    out = {
        "_what": ("Seeded-mutation results for this folder's checks. A seed is "
                  "CREDITED only if the metric line appears in both runs and "
                  "the metric moves in the expected direction."),
        "n_seeds": len(results),
        "n_credited": sum(1 for r in results if r["credited"]),
        "all_credited": all(r["credited"] for r in results),
        "curation_folder_artifacts_untouched_by_the_suite": untouched,
        "f069_posture": ("every seeded run received its own OUTDIR under "
                         f"{tmp}; the tracked-artifact hashes above are "
                         "compared across the whole suite, which is the only "
                         "check that sees a write from a run that exited "
                         "non-zero"),
        "hashes_before": before,
        "hashes_after": after,
        "seeds": results,
    }
    with open(os.path.join(HERE, "07_seed_results.json"), "w") as f:
        json.dump(out, f, indent=1)
    shutil.rmtree(tmp, ignore_errors=True)

    for r in results:
        print(f'{"PASS" if r["credited"] else "FAIL"}  {r["seed"]:38s} '
              f'{r["metric"]:26s} {r["verdict"]}')
    print(f'\n{out["n_credited"]}/{out["n_seeds"]} credited; '
          f'folder artifacts untouched: {untouched}')
    return 0 if out["all_credited"] and untouched else 1


if __name__ == "__main__":
    sys.exit(main())
