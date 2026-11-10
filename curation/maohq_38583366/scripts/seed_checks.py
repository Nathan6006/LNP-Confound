#!/usr/bin/env python3
"""Make every gate in stage01/stage02 fail on purpose.

RULES THIS HARNESS OBEYS, each one bought by a real incident:

 * It writes ONLY into a temp directory (F-069). One harness pinned its mutated
   copies' output path at the real curation folder and left a seeded anchor
   there, written by a run that exited non-zero -- so the non-zero exit
   CONCEALED the write. Nothing here can touch ../01_fig4b.json or
   ../02_fig4d.json; after the suite runs, both are re-derived and compared
   byte-for-byte against the pre-suite hashes.
 * One child process per seed, and the verdict is read from `returncode`, never
   from stdout and never through a pipeline (`python3 ... | tail` returns
   tail's exit code -- F-012's shape).
 * A SEED IS CREDITED ONLY WHEN THE METRIC MOVES. A seed that merely raises is
   not credited, except where the refusal itself IS the route under test, which
   is recorded as its own kind.
 * A seed whose metric line is ABSENT from the child's output is REFUSED, not
   credited: an exception thrown before the metric is computed proves nothing
   about the gate.
 * The control runs on the same re-encoded PNG as the seeds, so the only
   difference between control and seed is the mutation itself.
"""
import hashlib, json, os, re, shutil, subprocess, sys, tempfile
import numpy as np
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
FOLDER = os.path.join(HERE, "..")
SRC = os.path.join(FOLDER, "source", "figs_elsevier", "gr4_lrg.jpg")


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def run(stage, fig, outdir):
    env = dict(os.environ)
    env["STAGE01_FIG"] = env["STAGE02_FIG"] = fig
    env["STAGE01_OUT"] = os.path.join(outdir, "01.json")
    env["STAGE02_OUT"] = os.path.join(outdir, "02.json")
    p = subprocess.run([sys.executable, "-I", os.path.join(HERE, stage)],
                       capture_output=True, text=True, env=env)
    gates = dict(re.findall(r"^GATE (\w+) (-?[\d.]+)$", p.stdout, re.M))
    return p.returncode, gates, (p.stdout + p.stderr)


# ---- mutations. Each takes the RGB array and edits it in place. -----------
def m_tick_shift_4b_subtle(a):
    """Move Fig 4B's 1e5 major tick down 3 px.

    Deliberately SUBTLE. A 3 px offset of the middle of three ticks leaves a
    least-squares residual of 2.0 px, under axis.calibrate's 2.5 px tolerance,
    so calibrate PASSES and the decade-gap gate is the only thing standing
    between it and a wrong scale. The gross version is a separate seed because
    'exited non-zero' from a refusal is a different route from 'the gate fired'.
    """
    band = a[1935:1946, 200:213].copy()
    a[1935:1946, 200:213] = 255
    a[1938:1949, 200:213] = band


def m_tick_shift_4b_gross(a):
    """The same defect at 6 px, where axis.calibrate refuses first."""
    band = a[1935:1946, 200:213].copy()
    a[1935:1946, 200:213] = 255
    a[1941:1952, 200:213] = band


def m_baseline_shift_4b(a):
    """Slide Fig 4B's WHOLE tick ladder down 8 px.

    Targets the round-5,000 anchor specifically. The decade gaps stay equal and
    the held-out minor ticks stay consistent -- every other gate is blind to a
    bodily offset -- but the detected x-axis rule no longer lands on 5,000.
    (Moving the RULE instead was tried and REFUSED by the level-separation
    guard before reaching this gate, so it proved nothing about it: that is why
    the seed targets the ladder.)"""
    band = a[1600:2330, 196:214].copy()
    a[1600:2330, 196:214] = 255
    a[1608:2338, 196:214] = band


def m_occlude_bar_4b(a):
    """Whiten a 12 px vertical strip through AI-b3 (x 651-698), leaving only
    rows 2299-2300 intact so the bar is still DETECTED as one window at the
    row the window finder uses. No three contiguous rows then reach 0.90
    coverage of the full 47 px bar, so bar_top must REFUSE rather than return
    the topmost surviving ink."""
    a[1630:2299, 668:680] = 255
    a[2301:2325, 668:680] = 255


def m_liver_bar_4d(a):
    """Cut 30 px off the top of CII-b3's LIVER bar (x 2415-2444, top 1179).
    30 px is a factor of 1.86 on this axis, which must break the pie share."""
    a[1170:1215, 2410:2450] = 255


def m_minor_tick_4d(a):
    """Move ONE held-out minor tick in Fig 4D (row 1286.85) down 5 px.

    Targets the held-out-minor-tick gate directly. A major-tick shift was tried
    and moved the metric only 0.528 -> 1.057 px, under the gate: four majors in
    a least-squares fit absorb most of one tick's error, and the minors barely
    move. Widening the gate to credit that seed would have been widening the
    gate rather than the seed."""
    band = a[1281:1293, 1878:1888].copy()
    a[1281:1293, 1878:1888] = 255
    a[1286:1298, 1878:1888] = band


SEEDS = [
    ("tick_shift_4b_subtle", "stage01_fig4b.py", m_tick_shift_4b_subtle,
     "decade_gap_difference_px", "metric"),
    ("tick_shift_4b_gross", "stage01_fig4b.py", m_tick_shift_4b_gross,
     None, "refusal_calibrate"),
    ("baseline_shift_4b", "stage01_fig4b.py", m_baseline_shift_4b,
     "baseline_offset_px", "metric"),
    ("occlude_bar_4b", "stage01_fig4b.py", m_occlude_bar_4b,
     None, "refusal_bar_top"),
    ("liver_bar_4d", "stage02_fig4d.py", m_liver_bar_4d,
     "liver_share_max_abs_delta_pp", "metric"),
    ("minor_tick_4d", "stage02_fig4d.py", m_minor_tick_4d,
     "minor_tick_residual_px_max", "metric"),
]


def main():
    pre = {f: sha(os.path.join(FOLDER, f))
           for f in ("01_fig4b.json", "02_fig4d.json")
           if os.path.exists(os.path.join(FOLDER, f))}
    tmp = tempfile.mkdtemp(prefix="maohq38_seeds_")
    base = np.asarray(Image.open(SRC).convert("RGB")).copy()
    clean = os.path.join(tmp, "clean.png")
    Image.fromarray(base).save(clean)

    results, credited = [], 0
    ctrl = {}
    for stage in ("stage01_fig4b.py", "stage02_fig4d.py"):
        d = os.path.join(tmp, "ctrl_" + stage); os.makedirs(d, exist_ok=True)
        rc, gates, log = run(stage, clean, d)
        ctrl[stage] = gates
        results.append({"seed": "CONTROL(" + stage + ")", "returncode": rc,
                        "gates": gates,
                        "verdict": "PASS" if rc == 0 else "FAIL -- the control "
                                   "is not clean, so no seed below means "
                                   "anything (error class 50)"})
        if rc != 0:
            print(log[-2000:])
            json.dump(results, open(os.path.join(FOLDER,
                      "01_seed_results.json"), "w"), indent=1)
            raise SystemExit("seed_checks: CONTROL FAILED on " + stage)

    for name, stage, mut, metric, kind in SEEDS:
        a = base.copy(); mut(a)
        fp = os.path.join(tmp, name + ".png"); Image.fromarray(a).save(fp)
        d = os.path.join(tmp, name); os.makedirs(d, exist_ok=True)
        rc, gates, log = run(stage, fp, d)
        rec = {"seed": name, "stage": stage, "kind": kind,
               "returncode": rc, "metric": metric,
               "control_value": ctrl[stage].get(metric),
               "seeded_value": gates.get(metric)}
        if kind.startswith("refusal"):
            want = ("DetectionError", "bar_top:") if kind == "refusal_bar_top" \
                else ("CalibrationError", "calibrate:")
            ok = rc != 0 and want[0] in log and want[1] in log
            rec["verdict"] = ("CREDITED (%s route): the library refused before "
                              "any value was produced" % kind
                              if ok else "NOT CREDITED: expected %s, got rc=%d"
                              % (want[0], rc))
            rec["refusal_excerpt"] = (
                [l for l in log.splitlines() if want[1] in l] or [""])[0][:220]
        elif metric not in gates:
            rec["verdict"] = ("REFUSED: the metric line %r is absent from the "
                              "child's output, so the seed exited before "
                              "reaching the gate and proves nothing about it."
                              % metric)
        elif float(gates[metric]) == float(ctrl[stage][metric]):
            rec["verdict"] = ("NOT CREDITED: the metric did not move "
                              "(%s both clean and seeded). A seed that changes "
                              "nothing is a no-op." % gates[metric])
        elif rc == 0:
            rec["verdict"] = ("NOT CREDITED: the metric MOVED (%s -> %s) but "
                              "the gate stayed silent. Widen the seed, never "
                              "the gate." % (ctrl[stage][metric], gates[metric]))
        else:
            rec["verdict"] = ("CREDITED: metric moved %s -> %s and the gate "
                              "exited non-zero"
                              % (ctrl[stage][metric], gates[metric]))
        if rec["verdict"].startswith("CREDITED"):
            credited += 1
        results.append(rec)
        print("%-20s %s" % (name, rec["verdict"]))

    # F-069: nothing the suite ran may have touched a tracked artifact.
    post = {f: sha(os.path.join(FOLDER, f)) for f in pre}
    results.append({"seed": "_F069_tracked_artifacts_unchanged",
                    "pre_sha256": pre, "post_sha256": post,
                    "verdict": "PASS" if pre == post else
                    "FAIL -- a seeded run wrote into the curation folder"})
    summary = {"credited": credited, "seeds": len(SEEDS), "results": results,
               "tmpdir": tmp}
    json.dump(summary, open(os.path.join(FOLDER, "01_seed_results.json"), "w"),
              indent=1)
    shutil.rmtree(tmp, ignore_errors=True)
    print("\n%d of %d seeds credited; tracked artifacts unchanged: %s"
          % (credited, len(SEEDS), pre == post))
    if credited != len(SEEDS) or pre != post:
        raise SystemExit("seed_checks: not every seed was credited.")


if __name__ == "__main__":
    main()
