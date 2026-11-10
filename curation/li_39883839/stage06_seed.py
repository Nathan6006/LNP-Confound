#!/usr/bin/env python3
"""Seed a defect into the VALUE checks and require each metric to move."""
import json, os, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__))
CASES = [
 ("V1", "V1_fig3H_printed_ROI_containment", "n_inside",
  "class 20: Fig 3I origin anchored one tick out (+5e7 on every value)"),
 ("V2", "V2_H2A4B4_two_panels", "ratio",
  "read panel K's LOG-styled minor ticks as the scale (x22 on a mid bar)"),
]


def run(mut=None):
    cmd = [sys.executable, "-I", os.path.join(HERE, "stage06_crosschecks.py")]
    if mut:
        return json.loads(subprocess.run(cmd + ["--mutate", mut],
                                         capture_output=True, text=True).stdout)
    subprocess.run(cmd, capture_output=True, text=True)
    return json.load(open(os.path.join(HERE, "06_value_qc.json")))


def main():
    base, out, ok = run(), [], True
    for mut, chk, met, what in CASES:
        g = run(mut)[chk][met]
        b = base[chk][met]
        moved = b != g
        ok &= moved
        out.append({"seeded": what, "check": chk, "metric": met,
                    "clean_run": b, "with_defect": g, "metric_moved": moved})
        print("%-3s %-38s %-10s %s -> %s  %s" % (mut, chk, met, b, g,
              "MOVED" if moved else "*** DID NOT MOVE ***"))
    ctrl = run()
    same = all(ctrl[c][m] == base[c][m] for _, c, m, _ in CASES)
    print("control run reproduces every clean metric:", same)
    json.dump({"cases": out, "control_reproduces": same, "all_moved": ok and same},
              open(os.path.join(HERE, "06_seeded_defects.json"), "w"), indent=1)
    return 0 if (ok and same) else 1


if __name__ == "__main__":
    sys.exit(main())
