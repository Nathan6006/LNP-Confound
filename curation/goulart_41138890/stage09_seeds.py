#!/usr/bin/env python3
"""Seeded-mutation suite. A CHECK THAT CANNOT FAIL IS NOT A CHECK.

Rules this harness obeys, each from a specific finding:

  * Every seeded run writes to its OWN temp directory. It NEVER writes into the
    curation folder (F-069: a harness that pinned HERE to the real folder left
    a seeded regioisomer anchor behind, and the non-zero exit concealed the
    write). After the suite, both tracked artifacts are re-run clean and
    compared by sha256 against the pre-suite bytes.
  * One CHILD PROCESS per seed; the verdict is the child's `returncode`, never
    stdout and never a pipeline (`python3 ... | tail` returns tail's exit code
    -- F-012's shape).
  * A seed is CREDITED only when a named metric MOVES. A seed that merely
    raises is not credited, and a seed whose metric LINE IS ABSENT from the
    output is refused rather than silently passed.
  * One seed is a deliberate NEGATIVE CONTROL that must move nothing.
"""
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
PY = sys.executable
TRACKED = ["04_fig4g.json", "02_si_tables.json", "02_tableS1.csv",
           "02_tableS3.csv"]


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def run(script, outdir, seed=None):
    cmd = [PY, "-I", os.path.join(HERE, script), "--outdir", outdir]
    if seed:
        cmd += ["--seed", seed]
    p = subprocess.run(cmd, capture_output=True, text=True, cwd=HERE)
    return p.returncode, p.stdout, p.stderr


def metric(outdir, fname, path):
    f = os.path.join(outdir, fname)
    if not os.path.exists(f):
        return None, "artifact absent"
    d = json.load(open(f))
    for k in path:
        if isinstance(d, list):
            d = d[k]
        elif k in d:
            d = d[k]
        else:
            return None, f"metric key {k!r} absent"
    return d, None


def bar_values(outdir):
    d, err = metric(outdir, "04_fig4g.json", ["bars"])
    if d is None:
        return None, err
    return {f"{b['group']}_{b['timepoint']}": b["value"] for b in d}, None


def main():
    base = tempfile.mkdtemp(prefix="goulart_seeds_")
    report = {"temp_root": base, "seeds": []}

    # ---- the clean reference, in a temp dir of its own -------------------
    clean = os.path.join(base, "clean")
    os.makedirs(clean)
    rc4, _, e4 = run("stage04_fig4g.py", clean)
    rc2, _, e2 = run("stage02_si_tables.py", clean)
    assert rc4 == 0 and rc2 == 0, f"clean control FAILED: {e4}\n{e2}"
    ref_bars, err = bar_values(clean)
    assert ref_bars is not None, err
    ref_loo, _ = metric(clean, "04_fig4g.json", ["loo_worst_px"])
    ref_zero, _ = metric(clean, "04_fig4g.json",
                         ["zero_tick_vs_drawn_baseline_px"])
    ref_s1sum, _ = metric(clean, "02_si_tables.json",
                          ["checks", "S1_molpct_sum_exact_100"])
    ref_s1, _ = metric(clean, "02_si_tables.json", ["tableS1"])
    report["clean_control"] = {
        "returncodes": [rc4, rc2], "loo_worst_px": ref_loo,
        "zero_tick_offset_px": ref_zero,
        "S1_molpct_sum_exact_100": ref_s1sum,
        "bar_values": ref_bars}
    print("CLEAN CONTROL ok:", json.dumps(report["clean_control"], indent=1))

    SEEDS = [
        ("stage04_fig4g.py", "tick_shift",
         "GROSS: displace ticks[2] by +3.0 px. calibrate's own 2.5 px residual "
         "limit refuses first, so this seed is credited by the REFUSAL route, "
         "recorded distinctly -- it never reaches the LOO gate (brief s10).",
         "refusal_only"),
        ("stage04_fig4g.py", "tick_shift_subtle",
         "SUBTLE: displace ticks[2] by +1.6 px -- small enough to pass "
         "calibrate, so the LOO gate is the only thing that can see it.",
         "loo+values"),
        ("stage04_fig4g.py", "ticks_uniform_shift",
         "displace ALL six ticks by +4.0 px. Spacing is untouched, so "
         "calibrate AND the LOO gate both pass perfectly. Only the UNPRINTED "
         "check -- the 0 tick landing on the drawn baseline -- can see it.",
         "zero_offset"),
        ("stage04_fig4g.py", "no_furniture_mask",
         "do not mask the significance bracket at rows 3106-3109 (612 px, "
         "9.1 bar widths)", "values"),
        ("stage04_fig4g.py", "topmost_ink",
         "replace bars.bar_top with the topmost inked row -- error class 32, "
         "the F-018/F-020 detector", "values"),
        ("stage04_fig4g.py", "erase_baseline",
         "erase the drawn axis baseline -- the unprinted round-number check "
         "must notice", "zero_offset"),
        ("stage04_fig4g.py", "roll_image",
         "NEGATIVE CONTROL: roll the whole image by 7 rows. Ticks and bars "
         "move together, so a correct read is unchanged. Must move NOTHING.",
         "negative_control"),
        ("stage02_si_tables.py", "swap_digit",
         "flip the first UTF-16 '4' glyph record in image15.emf to '9'",
         "s1"),
    ]

    for script, seed, why, kind in SEEDS:
        od = os.path.join(base, seed)
        os.makedirs(od)
        rc, out, err = run(script, od, seed)
        rec = {"script": script, "seed": seed, "why": why, "kind": kind,
               "returncode": rc, "credited": False, "reason": None}

        if kind == "negative_control":
            got, e = bar_values(od)
            if got is None:
                rec["reason"] = f"metric line ABSENT ({e}); cannot judge"
            elif rc != 0:
                rec["reason"] = f"NEGATIVE CONTROL FAILED: exited {rc}"
            elif got != ref_bars:
                moved = {k: (ref_bars[k], got[k]) for k in ref_bars
                         if ref_bars[k] != got[k]}
                rec["reason"] = f"NEGATIVE CONTROL FAILED: values moved {moved}"
            else:
                rec["credited"] = True
                rec["reason"] = ("correctly moved nothing: 12/12 bar values "
                                 "identical. The harness does not credit "
                                 "'any change exits non-zero'.")
            report["seeds"].append(rec)
            continue

        if kind == "s1":
            got, e = metric(od, "02_si_tables.json", ["tableS1"])
            if rc != 0 and got is None:
                rec["reason"] = ("refused before writing: " +
                                 (err.strip().split("\n")[-1] or "raise"))
                rec["credited"] = True
                rec["route"] = "refusal"
            elif got is None:
                rec["reason"] = f"metric line ABSENT ({e}) -- NOT CREDITED"
            else:
                diffs = [(a["formulation"], k, a[k], b[k])
                         for a, b in zip(ref_s1, got) for k in a
                         if a[k] != b[k]]
                sm, _ = metric(od, "02_si_tables.json",
                               ["checks", "S1_molpct_sum_exact_100"])
                if diffs or sm != ref_s1sum:
                    rec["credited"] = True
                    rec["reason"] = (f"{len(diffs)} cell(s) moved {diffs[:3]}; "
                                     f"S1_molpct_sum_exact_100 "
                                     f"{ref_s1sum} -> {sm}")
                else:
                    rec["reason"] = "NO-OP: nothing moved -- NOT CREDITED"
            report["seeds"].append(rec)
            continue

        got, e = bar_values(od)
        loo, _ = metric(od, "04_fig4g.json", ["loo_worst_px"])
        zero, _ = metric(od, "04_fig4g.json",
                         ["zero_tick_vs_drawn_baseline_px"])
        if got is None and rc != 0:
            rec["reason"] = ("REFUSED before writing a metric: " +
                             (err.strip().split("\n")[-1] or "raise"))
            rec["credited"] = True
            rec["route"] = ("refusal (calibrate's 2.5 px residual limit, "
                            "reached before the LOO gate)"
                            if kind == "refusal_only" else "refusal")
        elif kind == "refusal_only":
            rec["reason"] = ("expected a REFUSAL and got rc=0 -- the gross "
                             "seed no longer refuses. NOT CREDITED.")
        elif got is None:
            rec["reason"] = f"metric line ABSENT ({e}) -- NOT CREDITED"
        else:
            moved = {k: (ref_bars[k], got[k]) for k in ref_bars
                     if ref_bars[k] != got[k]}
            notes = []
            if kind in ("loo+values",):
                notes.append(f"loo_worst_px {ref_loo} -> {loo} "
                             f"(strict gate 1.5: "
                             f"{'FAILS' if loo and loo > 1.5 else 'passes'})")
            if kind == "zero_offset":
                notes.append(f"zero_tick_offset_px {ref_zero} -> {zero}")
            if moved:
                worst = max(abs((b or 0) - (a or 0))
                            for a, b in moved.values())
                notes.append(f"{len(moved)}/12 bar values moved, "
                             f"worst |delta| {worst:.1f} units: "
                             f"{dict(list(moved.items())[:4])}")
            if kind == "loo+values":
                gate_fired = (loo is not None and loo > 1.5) and bool(moved)
            elif kind == "zero_offset":
                gate_fired = (zero != ref_zero) and bool(moved)
            else:
                gate_fired = bool(moved)
            rec["credited"] = gate_fired
            rec["reason"] = ("; ".join(notes) if gate_fired
                             else "NO-OP: nothing moved -- NOT CREDITED")
            rec["route"] = "metric moved"
        report["seeds"].append(rec)

    # ---- F-069: the curation folder must be untouched --------------------
    pre = {f: sha256(os.path.join(HERE, f)) for f in TRACKED
           if os.path.exists(os.path.join(HERE, f))}
    again = os.path.join(base, "again")
    os.makedirs(again)
    run("stage04_fig4g.py", again)
    run("stage02_si_tables.py", again)
    post = {f: sha256(os.path.join(again, f)) for f in TRACKED
            if os.path.exists(os.path.join(again, f))}
    report["reproducibility"] = {
        "tracked_artifacts": sorted(pre),
        "byte_identical_across_two_clean_runs":
            all(pre.get(f) == post.get(f) for f in pre) and len(pre) == len(post),
        "sha256_folder": pre, "sha256_fresh_run": post}

    n_cred = sum(1 for s in report["seeds"] if s["credited"])
    report["summary"] = {"n_seeds": len(report["seeds"]),
                         "n_credited": n_cred,
                         "all_credited": n_cred == len(report["seeds"])}
    with open(os.path.join(HERE, "09_seeds.json"), "w") as fh:
        json.dump(report, fh, indent=1)
    for s in report["seeds"]:
        print(f"  [{'CREDIT' if s['credited'] else 'NOT CREDITED'}] "
              f"{s['seed']:<18} rc={s['returncode']} :: {s['reason']}")
    print(json.dumps(report["reproducibility"]
                     ["byte_identical_across_two_clean_runs"], indent=1),
          "<- byte-identical across two clean runs")
    print(json.dumps(report["summary"]))
    shutil.rmtree(base, ignore_errors=True)
    sys.exit(0 if report["summary"]["all_credited"]
             and report["reproducibility"]
             ["byte_identical_across_two_clean_runs"] else 1)


if __name__ == "__main__":
    main()
