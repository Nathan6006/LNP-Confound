#!/usr/bin/env python3
"""Seeded-mutation gate for the Fig. 2B digitisation.

Run:  python3 -I curation/dahlmanje_41183194/stage03_seedgate.py <repo_root>

A CHECK THAT CANNOT FAIL IS NOT A CHECK, AND A SEED THAT ONLY RAISES IS NOT
CREDITED. Rules this harness enforces on itself:

  * every seed runs in its OWN PROCESS and the verdict is read from
    `returncode`. (`python3 ... | tail` returns tail's exit code.)
  * a seed is credited only if its declared gate actually MOVED. A seed whose
    process died before reaching its metric line is recorded `CRASH` and is
    NOT credited -- that is how a seed that exits on a TypeError gets mistaken
    for a catch.
  * a seed credited BY A REFUSAL must name the refusal it expected, and the
    refusal text must contain it. "Something raised" is not a result.
  * the control must PASS. Teeth with no clean control fail everything;
    a gate with no teeth passes everything. Both halves are required.
"""
import json
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else ".")
TARGET = os.path.join(HERE, "stage02_fig2b.py")

SEEDS = [
    {
        "name": "tick_jitter",
        "defect": "One y tick displaced +3 px in the tick strip; the bars and "
                  "the other ticks are untouched.",
        "gate": "axis.calibrate's UNIFORM-SPACING assertion, which fires "
                "before the leave-one-out test is ever reached. Recorded "
                "honestly: on a panel whose ticks are evenly spaced, any "
                "non-uniform displacement trips uniform-spacing first, so LOO "
                "is NOT independently exercised here and its measured 0.5 / "
                "0.0 px residual is evidence of fit quality, not a gate this "
                "curation proved has teeth. A uniform displacement, which LOO "
                "also cannot see, is seeded separately as `origin_shift`.",
        "credited_by": "refusal",
        "expect_in_stderr": "equally spaced",
    },
    {
        "name": "stretch_scale",
        "defect": "Every bar body moved 6 px up. The ticks are untouched, so "
                  "the calibration is still perfect and LOO stays clean.",
        "gate": "the printed Fig. 3C anchor (54 / 78), which took no part in "
                "any geometry. 6 px = 3.4 percentage points > the 2.0 tol.",
        "credited_by": "metric",
        "expect_metric": "anchor_verdict=FAIL",
    },
    {
        "name": "origin_shift",
        "defect": "All y ticks moved 8 px uniformly -- a pure origin error. "
                  "Shifted UP, not down, so that no tick is pushed out of the "
                  "detection window: a seed that changes how many ticks are "
                  "found crashes the calibration and proves nothing about the "
                  "gate it was aimed at.",
        "gate": "the printed Fig. 3C anchor. Spacing stays uniform and LOO "
                "stays clean, so no calibration gate can see this.",
        "credited_by": "metric",
        "expect_metric": "anchor_verdict=FAIL",
    },
    {
        "name": "fill_colour_mask",
        "defect": "HISTORICAL REPLAY of the defect this curation actually "
                  "hit: mask keyed on the bar's fill colour instead of its "
                  "ink, used without the solidity walk. The bar's darker top "
                  "stroke is missed and the replicate dots punch holes in the "
                  "fill, so the two flanks of liver bar 61 disagree.",
        "gate": "bars.bar_top's full-width-coverage and left/right flank "
                "rules. These are raised PER BAR and collected, so the run "
                "reaches its metric line with every bar refused rather than "
                "dying -- credited on the refusal TEXT, not on a crash.",
        "credited_by": "refusal_text",
        "expect_in_stderr": "flank",
    },
    {
        "name": "no_solidity",
        "defect": "The solidity walk dropped, so free-floating replicate "
                  "clusters re-enter the mask.",
        "gate": "the spleen-62 value itself, reported on the METRIC line.",
        "credited_by": "metric",
        "expect_metric_change": "spleen62",
    },
]


def run(seed=None):
    cmd = [sys.executable, "-I", TARGET, ROOT]
    if seed:
        cmd += ["--seed", seed]
    p = subprocess.run(cmd, capture_output=True, text=True)
    m = re.search(r"^METRIC (.*)$", p.stdout, re.M)
    metric = m.group(1) if m else None
    return {"returncode": p.returncode, "metric": metric,
            "stderr_tail": p.stderr.strip().splitlines()[-1] if p.stderr.strip()
            else "", "reached_metric": metric is not None}


def main():
    report = {"_what": "Seeded-mutation gate, KG_2025 Fig. 2B.",
              "rule": "A seed is credited only if its own declared gate moved. "
                      "A crash before the metric line is NOT a catch."}

    control = run(None)
    report["control"] = control
    control_ok = control["returncode"] == 0 and control["reached_metric"] \
        and "anchor_verdict=PASS" in (control["metric"] or "")
    report["control_verdict"] = "PASS" if control_ok else "FAIL"

    base_sp62 = None
    if control["metric"]:
        mm = re.search(r"spleen62=([\d.]+)", control["metric"])
        base_sp62 = float(mm.group(1)) if mm else None

    results = []
    for s in SEEDS:
        r = run(s["name"])
        rec = dict(s)
        rec.update(r)
        credited, why = False, ""

        if s["credited_by"] == "refusal_text":
            txt = (r["metric"] or "")
            nref = re.search(r"n_refusals=(\d+)", txt)
            nref = int(nref.group(1)) if nref else 0
            if nref == 0:
                why = ("NOT CREDITED: no bar was refused, so the rule this "
                       "seed targets did not fire.")
            elif s["expect_in_stderr"].lower() not in txt.lower():
                why = (f"NOT CREDITED: {nref} refusal(s), but none names "
                       f"{s['expect_in_stderr']!r}: {txt[:200]}")
            elif r["returncode"] == 0:
                why = "NOT CREDITED: exited 0."
            else:
                credited = True
                why = (f"CREDITED by refusal text: {nref} bar(s) refused and "
                       f"the refusal names {s['expect_in_stderr']!r}.")
        elif s["credited_by"] == "refusal":
            if r["reached_metric"]:
                why = ("NOT CREDITED: the run reached its metric line, so the "
                       "refusal this seed targets did not fire.")
            elif r["returncode"] == 0:
                why = "NOT CREDITED: exited 0."
            elif s["expect_in_stderr"].lower() not in r["stderr_tail"].lower():
                why = (f"NOT CREDITED: it raised, but not with the expected "
                       f"refusal {s['expect_in_stderr']!r}. Got: "
                       f"{r['stderr_tail'][:160]}")
            else:
                credited = True
                why = (f"CREDITED by refusal: non-zero exit and the refusal "
                       f"names {s['expect_in_stderr']!r}.")
        else:
            if not r["reached_metric"]:
                why = ("NOT CREDITED (CRASH): the process died before its "
                       "metric line, so nothing proves the gate moved.")
            elif "expect_metric" in s:
                if s["expect_metric"] in r["metric"] and r["returncode"] != 0:
                    credited = True
                    why = f"CREDITED by metric: {s['expect_metric']}"
                else:
                    why = (f"NOT CREDITED: expected {s['expect_metric']!r} and "
                           f"a non-zero exit; got rc={r['returncode']} "
                           f"metric={r['metric']!r}")
            else:
                mm = re.search(r"spleen62=([\d.]+)", r["metric"] or "")
                got = float(mm.group(1)) if mm else None
                if got is None or base_sp62 is None:
                    why = "NOT CREDITED: metric absent."
                elif abs(got - base_sp62) < 1.0:
                    why = (f"NOT CREDITED (NO-OP SEED): spleen62 moved "
                           f"{base_sp62} -> {got}, under 1 percentage point. "
                           "The mutation did not reach the quantity.")
                else:
                    credited = True
                    why = (f"CREDITED by metric: spleen62 moved {base_sp62} "
                           f"-> {got}, a {abs(got-base_sp62):.2f} point shift.")
        rec["credited"] = credited
        rec["why"] = why
        results.append(rec)

    report["seeds"] = results
    report["n_seeds"] = len(results)
    report["n_credited"] = sum(1 for r in results if r["credited"])
    report["verdict"] = ("PASS" if control_ok and
                         report["n_credited"] == len(results) else "FAIL")

    with open(os.path.join(HERE, "03_seedgate.json"), "w") as fh:
        json.dump(report, fh, indent=1)
    print(json.dumps(report, indent=1))
    return 0 if report["verdict"] == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
