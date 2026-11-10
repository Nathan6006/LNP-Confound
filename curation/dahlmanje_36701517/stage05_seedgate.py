"""Drive seed_runner.py, one `python3 -I` PROCESS PER SEED.

The verdict for every seed is `proc.returncode`. Nothing is read from stdout
except the `METRIC` line, and a seed whose METRIC line is ABSENT is NOT
CREDITED whatever it exited with -- that is the "exited non-zero from a
TypeError and never reached its metric" shape.

A seed is credited only when the AGREEMENT COUNT MOVES against the control:
either the number of bars read changes, or a value changes by more than
0.25 percentage points (0.6 px on this 2.378 px/pct axis). A bare non-zero
exit earns nothing, and is recorded as its own route when the refusal happens
upstream of any count.

Usage: python3 -I stage05_seedgate.py   (writes 09_seedgate.json, merges into
                                         07_qc.json, exits 1 if the control
                                         fails or no seed is credited)
"""
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
RUNNER = os.path.join(HERE, "seed_runner.py")
MOVE_PCT = 0.25

SEEDS = [
    ("control", None, "the unmutated panel; must read 6 of 6"),
    ("tick_shift", "6", "top y tick stub moved +6 px: a gross scale error"),
    ("tick_shift", "3", "top y tick stub moved +3 px"),
    ("tick_shift", "2", "top y tick stub moved +2 px"),
    ("tick_shift", "1", "top y tick stub moved +1 px: below the 1.5 px gate"),
    ("category_stub_deleted", None,
     "one x category stub erased: six labels, five slots"),
    ("marker", "10", "a 10 px wide replicate marker on an 18 px bar (56%)"),
    ("marker", "16", "a 16 px wide replicate marker (89%, just under the "
                     "0.90 coverage rule)"),
    ("marker", "17", "a 17 px wide replicate marker (94%)"),
    ("marker", "18", "a full-bar-width marker (100%)"),
    ("spine_erased", None,
     "the y-axis spine erased: removes the furniture the 100%-level check "
     "leans on"),
]


def run(name, arg):
    cmd = [sys.executable, "-I", RUNNER, name] + ([arg] if arg else [])
    p = subprocess.run(cmd, capture_output=True, text=True, cwd=HERE)
    metric = None
    for line in p.stdout.splitlines():
        if line.startswith("METRIC"):
            metric = line.split()
    return {"returncode": p.returncode, "metric_line_present": metric
            is not None, "metric": metric,
            "stdout_tail": p.stdout.strip().splitlines()[-2:],
            "stderr_head": p.stderr.strip().splitlines()[:2]}


def parse(metric):
    if not metric:
        return None, []
    n = int(metric[1])
    vals = [None if t == "None" else float(t) for t in metric[2:]]
    return n, vals


def main():
    out = {"_what": "seeded mutations, ONE PROCESS PER SEED, verdict from "
                    "returncode; a seed with no METRIC line is NOT CREDITED",
           "move_threshold_pct": MOVE_PCT,
           "exit_codes": {"0": "panel read", "2": "detector refused",
                          "3": "crashed before the metric"},
           "seeds": []}

    ctrl = run("control", None)
    if ctrl["returncode"] != 0 or not ctrl["metric_line_present"]:
        print("CONTROL FAILED -- a seed proves nothing against a broken "
              "control.", ctrl)
        sys.exit(1)
    cn, cvals = parse(ctrl["metric"])
    out["control"] = {"returncode": 0, "bars_read": cn, "values_pct": cvals}

    credited = 0
    for name, arg, why in SEEDS:
        if name == "control":
            continue
        r = run(name, arg)
        rec = {"seed": name + (f" {arg}" if arg else ""), "why": why,
               "returncode": r["returncode"],
               "metric_line_present": r["metric_line_present"]}
        if not r["metric_line_present"]:
            rec["credited"] = False
            rec["verdict"] = (
                "NOT CREDITED -- no METRIC line. The process exited "
                f"{r['returncode']} without ever producing its measurement, "
                "so nothing about the check was tested.")
            rec["stdout_tail"] = r["stdout_tail"]
            rec["stderr_head"] = r["stderr_head"]
            out["seeds"].append(rec)
            continue
        n, vals = parse(r["metric"])
        rec["bars_read"] = n
        rec["values_pct"] = vals
        moved_vals = sum(
            1 for a, b in zip(cvals, vals)
            if a is not None and b is not None and abs(a - b) > MOVE_PCT)
        count_moved = (n != cn)
        rec["bars_read_moved"] = f"{cn} -> {n}" if count_moved else "no"
        rec["values_that_moved"] = moved_vals
        if r["returncode"] == 2:
            rec["credited"] = count_moved
            rec["route"] = "refusal"
            rec["verdict"] = (
                f"CREDITED by REFUSAL: the detector refused and the bars-read "
                f"count moved {cn} -> {n}, so real output was removed."
                if count_moved else
                "NOT CREDITED: the detector refused but the count did not "
                "move, so the refusal happened upstream of anything being "
                "measured.")
            rec["refusal_text"] = r["stdout_tail"][0] if r["stdout_tail"] \
                else ""
        else:
            rec["credited"] = count_moved or moved_vals > 0
            rec["route"] = "value change"
            if count_moved and not moved_vals:
                rec["route"] = "silent bar loss"
                rec["verdict"] = (
                    f"CREDITED by BAR LOSS: no exception was raised, but the "
                    f"bars-read count moved {cn} -> {n} because the slot "
                    "picker found no bar of measurable width in that slot. "
                    "Output was removed, which is the count moving.")
            elif moved_vals:
                rec["verdict"] = (
                    f"CREDITED as a demonstrated LEAK, not as a check with "
                    f"teeth: nothing refused, the panel still reported {n} "
                    f"bars, and {moved_vals} value(s) changed by more than "
                    f"{MOVE_PCT} pct. This is a measured blind spot.")
            if not (count_moved or moved_vals):
                rec["harmless"] = True
                rec["verdict"] = (
                    "HARMLESS, NOT CREDITED: nothing refused and no value "
                    f"moved by more than {MOVE_PCT} pct. Correct behaviour, "
                    "but it is not evidence that any check fired.")
        credited += bool(rec["credited"])
        out["seeds"].append(rec)

    out["n_seeds"] = len(out["seeds"])
    out["n_credited"] = credited
    out["n_harmless"] = sum(1 for s in out["seeds"] if s.get("harmless"))
    out["n_leaks"] = sum(1 for s in out["seeds"]
                         if s.get("route") == "value change"
                         and s.get("values_that_moved"))
    with open(os.path.join(HERE, "09_seedgate.json"), "w") as fh:
        json.dump(out, fh, indent=1, sort_keys=True, default=str)
    qc_path = os.path.join(HERE, "07_qc.json")
    qc = json.load(open(qc_path))
    qc["seed_gate_one_process_per_seed"] = out
    json.dump(qc, open(qc_path, "w"), indent=1, sort_keys=True, default=str)

    for s in out["seeds"]:
        print("%-28s rc=%s credited=%-5s %s"
              % (s["seed"], s["returncode"], s["credited"],
                 s["verdict"][:88]))
    print("\ncontrol %d bars; credited %d/%d; harmless %d; leaks %d"
          % (cn, credited, out["n_seeds"], out["n_harmless"], out["n_leaks"]))
    sys.exit(0 if credited else 1)


if __name__ == "__main__":
    main()
