#!/usr/bin/env python3
"""Seed one defect per check and require the CHECK'S OWN METRIC to move.

CLAUDE.md section 6: "A NEW CHECK IS NOT DONE UNTIL IT HAS BEEN MADE TO FAIL
ON PURPOSE." The rule here is stricter than 'an exception was raised': three
ZT_2026 marginal checks scored a literal 0 for twelve days while the summary
beside them said ALL FOUR EXACT, and two versions of C16 reported clean while
examining nothing. Every one printed a clean line. So each case below names
the metric it expects to move, and the control run must leave it alone.
"""
import json, subprocess, sys, os
HERE = os.path.dirname(os.path.abspath(__file__))

CASES = [
 ("P1", "P1_tableS1_formula", "n_mismatch",
  "swap the head groups of H1A1B1 and H2A1B1 in the enumeration"),
 ("P2", "P2_tableS1_MW", "max_abs_delta_Da",
  "set one enumerated MW to 999.99"),
 ("P4", "P4_reagent_stoichiometry", "n_lipids_matching_as_a_multiset",
  "believe the printed reagent NAME for A4B4-S3 (B3 instead of B4)"),
 ("P5", "P5_block_MS", "n_matching_neutral_monoisotopic_within_0.02",
  "corrupt B4's printed calc'd mass 354.31 -> 340.0"),
]


def run(mut=None):
    cmd = [sys.executable, "-I", os.path.join(HERE, "stage03_checks.py")]
    if mut:
        cmd += ["--mutate", mut]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if mut:
        return json.loads(r.stdout)
    return json.load(open(os.path.join(HERE, "04_structure_qc.json")))


def main():
    base = run()
    out, ok = [], True
    for mut, check, metric, what in CASES:
        got = run(mut)
        b, g = base[check][metric], got[check][metric]
        moved = (b != g)
        ok &= moved
        out.append({"seeded": what, "check": check, "metric": metric,
                    "clean_run": b, "with_defect": g, "metric_moved": moved})
        print("%-4s %-28s %-42s %s -> %s  %s"
              % (mut, check, metric, b, g, "MOVED" if moved else "*** DID NOT MOVE ***"))
    # the control: an unmutated run must reproduce the clean values exactly
    ctrl = run()
    same = all(ctrl[c][m] == base[c][m] for _, c, m, _ in CASES)
    print("control run reproduces every clean metric:", same)
    ok &= same
    json.dump({"cases": out, "control_reproduces": same, "all_moved": ok},
              open(os.path.join(HERE, "04_seeded_defects.json"), "w"), indent=1)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
