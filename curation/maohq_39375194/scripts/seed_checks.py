#!/usr/bin/env python3
"""Seeded-mutation harness for stage01_invert.py.

RUN:  python3 -I scripts/seed_checks.py <1080_LNP_OHE.csv> <Compiled_validation.xlsx> <workdir> <out.json>

CLAUDE.md section 6: "A NEW CHECK IS NOT DONE UNTIL IT HAS BEEN MADE TO FAIL ON
PURPOSE." And the stronger form from this task's brief, which is the one that
matters: a seed is credited ONLY WHEN THE AGREEMENT COUNT MOVES. A sibling's
category check scored 11/11 clean AND 11/11 seeded because it measured spread
about the mean of its own set, which is shift-invariant -- exit code flipped on
nothing. So every seed below names a NUMERIC METRIC by path, and the seed is
credited only if that number differs from baseline. A seed that merely flips a
boolean, or merely raises the exit code, is reported as NOT CREDITED.

Each seed runs in its OWN PROCESS and its returncode is read from
CompletedProcess.returncode, never inferred from output.
"""
import csv
import json
import os
import shutil
import subprocess
import sys

from openpyxl import load_workbook

HERE = os.path.dirname(os.path.abspath(__file__))
STAGE01 = os.path.join(HERE, "stage01_invert.py")


def dig(d, path):
    cur = d
    for k in path.split("."):
        if isinstance(cur, list):
            cur = cur[int(k)]
        else:
            cur = cur[k]
    return cur


def read_csv(p):
    with open(p, newline="", encoding="utf-8-sig") as fh:
        rd = csv.DictReader(fh)
        return rd.fieldnames, list(rd)


def write_csv(p, fields, rows):
    with open(p, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


# ----------------------------------------------------------- the mutations
def m_onehot_two_ones(fields, rows, _):
    rows[4]["DOPE"] = "1"          # row already has 18PG = 1


def m_onehot_disagrees(fields, rows, _):
    rows[6]["Helper_lipid"] = "DSPC"   # one-hot still says 18PG


def m_factorial_duplicate(fields, rows, _):
    rows[2]["NP_ratio"] = "8"      # collides with rows[1]'s combination


def m_level_count(fields, rows, _):
    rows[10]["(IL+HL)"] = "25"     # a fifth (IL+HL) level the SI never prints


def m_class40_rowwise(fields, rows, _):
    rows[0]["HL_(IL+HL)"] = "0.5"  # ratio 200 implies 0.497512438


def m_class40_integrality(fields, rows, _):
    # keep the row-wise relation EXACT but break the integrality of the ratio
    rows[0]["Chol_DMG-PEG_ratio"] = "499.5"
    rows[0]["PEG_(Chol+PEG)"] = repr(100.0 / 500.5)


def m_descriptor_nonconstant(fields, rows, _):
    rows[0]["cLogP"] = "7.7"       # one 18PG row disagrees with its siblings


def m_descriptor_vs_si(fields, rows, _):
    for r in rows:
        if r["Helper_lipid"] == "18PG":
            r["Double_bonds"] = "1"   # SI Table 3 prints 2


def m_blank_a_measurement(fields, rows, _):
    rows[0]["RLU_HepG2"] = ""


def m_rescale_to_unit(fields, rows, _):
    v = [float(r["RLU_HEK293"]) for r in rows]
    lo, hi = min(v), max(v)
    for r in rows:
        r["RLU_HEK293"] = repr((float(r["RLU_HEK293"]) - lo) / (hi - lo))


def v_overlap_level(ws):
    ws.cell(row=2, column=21).value = 12     # NP_ratio 10 -> 12, a training level


def v_descriptor_mismatch(ws):
    ws.cell(row=2, column=6).value = 9.9     # cLogP of DOTAP


def v_shared_formulation(ws):
    # make validation row 1 equal training row 1's formulation exactly
    ws.cell(row=2, column=3).value = "18PG"
    ws.cell(row=2, column=21).value = 4            # NP_ratio
    ws.cell(row=2, column=22).value = 0.497512438  # HL_(IL+HL)
    ws.cell(row=2, column=23).value = 20           # (IL+HL)
    ws.cell(row=2, column=24).value = 0.199600798  # PEG_(Chol+PEG)


SEEDS = [
    ("onehot_two_ones", m_onehot_two_ones, None,
     "onehot_inversion.rows_whose_group_sums_to_exactly_1",
     "a row with two 1s in the helper-lipid group"),
    ("onehot_disagrees_with_categorical", m_onehot_disagrees, None,
     "onehot_inversion.rows_where_inversion_equals_the_file's_Helper_lipid_column",
     "one-hot and the categorical column name different lipids"),
    ("factorial_duplicate_cell", m_factorial_duplicate, None,
     "factorial.combinations_appearing_exactly_once",
     "two rows occupy one factorial cell"),
    ("level_not_in_SI_Table_1", m_level_count, None,
     "si_table1_levels.(IL+HL).n_levels_in_file",
     "a composition level the SI never prints"),
    ("class40_rowwise_disagreement", m_class40_rowwise, None,
     "error_class_40_two_parameterisations.HL_(IL+HL).rows_violating_at_1e-6",
     "percent column disagrees with its own ratio column on one row"),
    ("class40_ratio_not_integer", m_class40_integrality, None,
     "error_class_40_two_parameterisations.PEG_(Chol+PEG)."
     "ratio_levels_are_exact_integers",
     "the ratio is no longer a round integer, with the row-wise relation "
     "left EXACT -- isolates integrality from the derivation"),
    ("descriptor_varies_within_lipid", m_descriptor_nonconstant, None,
     "helper_descriptor_block.non_constant",
     "one 18PG row carries a different cLogP from its siblings"),
    ("descriptor_contradicts_SI_Table_3", m_descriptor_vs_si, None,
     "helper_descriptor_block.vs_si_table3.18PG.agrees",
     "every 18PG row says 1 double bond where SI Table 3 prints 2"),
    ("measurement_blanked", m_blank_a_measurement, None,
     "measurement_coverage.blank_pattern_is_all_or_nothing",
     "one luciferase reading removed"),
    ("values_already_minmaxed", m_rescale_to_unit, None,
     "scale_is_ln_not_minmax.cells_whose_values_all_fit_[0,1]",
     "a target column arrives already rescaled to [0,1]"),
    ("validation_level_overlaps_training", None, v_overlap_level,
     "validation_workbook.levels_vs_training.NP_ratio.intersection",
     "a validation level that the training library also used"),
    ("validation_descriptor_mismatch", None, v_descriptor_mismatch,
     "validation_workbook.helper_descriptor_mismatches_vs_training",
     "the two files disagree on one helper lipid's cLogP"),
    ("shared_formulation_appears", None, v_shared_formulation,
     "shared_formulations_between_files.n",
     "a formulation present in both files, which would make a ratio "
     "measurable and must not go unnoticed"),
]


def run(csv_p, xlsx_p, out_p):
    r = subprocess.run([sys.executable, "-I", STAGE01, csv_p, xlsx_p, out_p],
                       capture_output=True, text=True)
    data = json.load(open(out_p)) if os.path.exists(out_p) else None
    return r.returncode, data


def main(csv_src, xlsx_src, work, out_path):
    os.makedirs(work, exist_ok=True)
    base_csv = os.path.join(work, "base.csv")
    base_xlsx = os.path.join(work, "base.xlsx")
    shutil.copy(csv_src, base_csv)
    shutil.copy(xlsx_src, base_xlsx)

    base_rc, base = run(base_csv, base_xlsx, os.path.join(work, "base.json"))
    results = {"control": {"returncode": base_rc, "passes": base_rc == 0}}
    if base_rc != 0:
        print(json.dumps({"FATAL": "control run does not pass"}, indent=1))
        return 1

    seeds = []
    for name, cmut, vmut, metric, what in SEEDS:
        d = os.path.join(work, name)
        os.makedirs(d, exist_ok=True)
        c = os.path.join(d, "t.csv")
        x = os.path.join(d, "v.xlsx")
        if cmut:
            fields, rows = read_csv(csv_src)
            cmut(fields, rows, None)
            write_csv(c, fields, rows)
            shutil.copy(xlsx_src, x)
        else:
            shutil.copy(csv_src, c)
            wb = load_workbook(xlsx_src)
            vmut(wb.worksheets[0])
            wb.save(x)

        rc, got = run(c, x, os.path.join(d, "o.json"))
        try:
            before = dig(base, metric)
        except Exception as e:
            before = "METRIC_PATH_ERROR:%s" % e
        try:
            after = dig(got, metric) if got else "NO_OUTPUT"
        except Exception as e:
            after = "METRIC_PATH_ERROR:%s" % e
        moved = (before != after)
        seeds.append({
            "seed": name,
            "defect_seeded": what,
            "metric": metric,
            "baseline_value": before,
            "seeded_value": after,
            "AGREEMENT_COUNT_MOVED": moved,
            "returncode": rc,
            "gate_fired": rc != 0,
            "CREDITED": bool(moved and rc != 0),
        })

    results["seeds"] = seeds
    results["n_seeds"] = len(seeds)
    results["n_credited"] = sum(1 for s in seeds if s["CREDITED"])
    results["n_moved_but_gate_silent"] = sum(
        1 for s in seeds if s["AGREEMENT_COUNT_MOVED"] and not s["gate_fired"])
    results["n_gate_fired_but_no_movement"] = sum(
        1 for s in seeds if s["gate_fired"] and not s["AGREEMENT_COUNT_MOVED"])
    results["_rule"] = ("A seed is credited only when the named metric's value "
                        "DIFFERS from baseline AND the gate exits non-zero. A "
                        "gate that fires without the metric moving is a no-op "
                        "seed, which is the failure mode this harness exists "
                        "to catch.")
    with open(out_path, "w") as fh:
        json.dump(results, fh, indent=1)
    print(json.dumps({k: v for k, v in results.items() if k != "seeds"},
                     indent=1))
    for s in seeds:
        print(("  %-38s moved=%-5s rc=%-2d CREDITED=%s   %s -> %s"
               % (s["seed"], s["AGREEMENT_COUNT_MOVED"], s["returncode"],
                  s["CREDITED"], s["baseline_value"], s["seeded_value"]))[:200])
    return 0 if results["n_credited"] == len(seeds) else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]))
