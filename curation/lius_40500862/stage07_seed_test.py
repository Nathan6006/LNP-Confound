#!/usr/bin/env python3
"""Make every check in stage07_checks.py FAIL ON PURPOSE.

CLAUDE.md section 6: "A NEW CHECK IS NOT DONE UNTIL IT HAS BEEN MADE TO FAIL
ON PURPOSE." Three incidents paid for that rule, and all three printed a clean
line while examining nothing.

A check is credited here only if, with its defect seeded:
  1. its OWN agreement count DROPS, and
  2. the gate's EXIT CODE is non-zero, read from the process directly and not
     through a pipe.

A seed that merely makes the gate exit 1 proves nothing about the check it was
aimed at -- some other check may have caught it. Hence condition 1.

Everything runs on a COPY in a temporary directory; this script never writes
inside the curation folder.
"""
import json
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
LIB = os.path.abspath(os.path.join(HERE, "..", "lib"))


def seed_csv(d, col, name, value):
    p = os.path.join(d, "03_enumerated.csv")
    lines = open(p).read().splitlines()
    head = lines[0].split(",")
    i = head.index(col)
    for n, ln in enumerate(lines[1:], 1):
        cells = ln.split(",")
        if cells[0] == name:
            cells[i] = str(value)
            lines[n] = ",".join(cells)
    open(p, "w").write("\n".join(lines) + "\n")


def seed_json(d, fn, path, value):
    p = os.path.join(d, fn)
    obj = json.load(open(p))
    cur = obj
    for k in path[:-1]:
        cur = cur[k]
    cur[path[-1]] = value
    json.dump(obj, open(p, "w"), indent=1)


SEEDS = {
    "P1_hrms_ours_vs_si_calculated_ppm":
        lambda d: seed_csv(d, "ppm_ours_vs_si_calc", "A5-CE-C8", 12.5),
    "P2_hrms_si_found_vs_si_calculated_ppm":
        lambda d: seed_csv(d, "ppm_si_found_vs_si_calc", "A5-CE-C9", 950.0),
    "P3_hrms_ours_vs_si_found_ppm":
        lambda d: seed_csv(d, "ppm_ours_vs_si_found", "A5-CE-C10", -44.0),
    "P4_carbonate_minus_carboxylate_is_n_oxygens":
        lambda d: seed_csv(d, "si_calc_MH_plus", "A5-C10", 741.6451),
    "P5_axis_leave_one_out_tick_residual":
        lambda d: seed_json(d, "01_digitised.json",
                            ["qc", "S5", "calibration", "loo_worst_px"], 6.1),
    "P6_axis_rule_sits_on_the_extreme_tick":
        lambda d: seed_json(d, "01_digitised.json",
                            ["qc", "S23A", "calibration", "ticks_px"],
                            [75.6, 127.8, 180.1, 232.3, 285.0, 325.0]),
    "P7_bar_top_left_right_flank_agreement":
        lambda d: seed_json(d, "01_digitised.json",
                            ["qc", "S22C", "bars", 0, "flank_delta_px"], 9.0),
    "P8_S5_marker_centroid_equals_bbox_centre":
        lambda d: seed_json(d, "01_digitised.json",
                            ["qc", "S5", "points", 3,
                             "centroid_minus_bboxmid_px"], 7.4),
    "P9_corpus_overlap_full_inchikey":
        lambda d: seed_json(d, "04_corpus_overlap.json",
                            ["full_inchikey_matches"], 3),
    "P10_linkage_is_a_carbonate_not_an_ester_or_carbamate":
        lambda d: _seed_c13(d),
}


def _seed_c13(d):
    """Move the carbonyl from the carbonate band to the ester band.

    This is the defect the check exists for: a structure whose linkage is an
    ester rather than a carbonate. Every mass check (P1-P4) still passes on it.
    """
    p = os.path.join(d, "stage07_checks.py")
    s = open(p).read()
    assert "c13 = [155.16," in s
    open(p, "w").write(s.replace("c13 = [155.16,", "c13 = [173.40,"))


def run_gate(d):
    """Run the gate as a process and return (exit_code, parsed 07_qc.json)."""
    r = subprocess.run([sys.executable, "-I", "stage07_checks.py"],
                       cwd=d, capture_output=True, text=True)
    qc = None
    p = os.path.join(d, "07_qc.json")
    if os.path.exists(p):
        qc = json.load(open(p))
    return r.returncode, qc, r.stdout


def snapshot(d):
    shutil.copytree(HERE, d, dirs_exist_ok=True)
    # curation/lib is imported by a relative path, so mirror the layout.
    lib = os.path.join(os.path.dirname(d), "lib")
    if not os.path.exists(lib):
        shutil.copytree(LIB, lib)


def main():
    root = tempfile.mkdtemp(prefix="lius_seed_")
    work = os.path.join(root, "paper")

    # --- the control: an unmutated copy must still PASS and exit 0 -------
    snapshot(work)
    code, qc, _ = run_gate(work)
    if code != 0:
        print(f"CONTROL FAILED: clean copy exited {code}. Teeth with no "
              "baseline fail every clean run; stop here.")
        return 2
    base = {k: v["agree"] for k, v in qc["checks"].items()}
    print(f"control: exit 0, {len(base)} checks, counts "
          f"{sum(base.values())}/{sum(v['of'] for v in qc['checks'].values())}")

    results, bad = [], []
    for name, seed in SEEDS.items():
        d = os.path.join(root, "m_" + name[:3])
        snapshot(d)
        seed(d)
        code, qc, out = run_gate(d)
        if qc is None:
            results.append((name, code, None, None, "gate wrote no QC file"))
            bad.append(name)
            continue
        now = qc["checks"][name]["agree"]
        moved = now < base[name]
        ok = moved and code != 0
        results.append((name, code, base[name], now, "" if ok else
                        ("count did not move" if not moved
                         else "exit code was 0")))
        if not ok:
            bad.append(name)

    print()
    for name, code, b, n, why in results:
        tag = "CREDITED" if not why else "NOT CREDITED"
        print(f"{tag:<13} {name:<52} count {b} -> {n}  exit {code}  {why}")

    print(f"\n{len(SEEDS) - len(bad)}/{len(SEEDS)} checks credited "
          "(own count moved AND gate exited non-zero)")
    shutil.rmtree(root, ignore_errors=True)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
