#!/usr/bin/env python3
"""KK_2015 -- seeded-mutation suite. DOES THIS CURATION'S CHECKS HAVE TEETH?

F-069: EVERY SEED RUNS IN A CHILD PROCESS WITH ITS OWN TEMP OUTDIR. Nothing in
this file, and nothing any child writes, may land in the curation folder. A
sibling's harness pinned each mutated copy's HERE to the real folder, the copies
wrote their outputs there, and a seeded regioisomer anchor ended up in the
tracked reference_structures.csv -- written by a run that exited non-zero, so
the non-zero exit CONCEALED the write.

A SEED THAT ONLY RAISES IS NOT CREDITED. Each seed must name the METRIC it
moved; a seed whose metric line is absent from the child's output is reported as
NOT CREDITED even if the child exited non-zero, because an exception from a
TypeError three lines before the metric proves nothing.

Verdicts come from `returncode`, never from a pipeline's exit status.
"""
import json
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
LIB = os.path.join(HERE, "..", "lib")
PROV = os.path.join(HERE, "..", "..", "provenance")
CUR = os.path.join(HERE, "..")

TRACKED = ["00_value_sources.json", "01_blocks.csv", "02_corpus_overlap.json",
           "03_enumerated.csv", "05_values.csv", "06_rows.csv",
           "reference_structures.csv"]


def run(code, outdir, extra_env=None):
    env = dict(os.environ)
    env["KK2015_OUTDIR"] = outdir
    env["PYTHONPATH"] = os.pathsep.join([HERE, LIB, CUR, PROV])
    if extra_env:
        env.update(extra_env)
    p = subprocess.run([sys.executable, "-I", "-c", code], cwd=HERE, env=env,
                       capture_output=True, text=True)
    return p


PRELUDE = (
    "import os, sys\n"
    "sys.path[:0] = os.environ['PYTHONPATH'].split(os.pathsep)\n"
)

SEEDS = []


def seed(name, metric, why):
    def deco(fn):
        SEEDS.append((name, metric, why, fn))
        return fn
    return deco


# ---------------------------------------------------------------- emit guards --
@seed("emit_literal_confidence", "EmitError",
      "F-011/F-025/F-046: a hard-coded confidence relabelled 24 'ASSUMED' rows "
      "as the project's highest confidence. emit.emit_rows must REFUSE a literal.")
def s1(outdir):
    return run(PRELUDE + """
import emit
try:
    emit.emit_rows(["x"], {"x": {}}, confidence="high",
                   note=lambda i, r: "n", label="seed")
    print("METRIC exception=None")
except Exception as e:
    print("METRIC exception=%s" % type(e).__name__)
    print("MSG", str(e)[:120].replace(chr(10), " "))
""", outdir)


@seed("emit_caveat_without_note", "EmitError",
      "F-024: confidence_without_a_note found 6 rows on its first run. Here the "
      "condition must be unwritable, not merely audited later.")
def s2(outdir):
    return run(PRELUDE + """
import emit
try:
    emit.emit_rows(["x"], {"x": {}}, confidence=lambda i, r: "medium",
                   note=lambda i, r: "", label="seed")
    print("METRIC exception=None")
except Exception as e:
    print("METRIC exception=%s" % type(e).__name__)
""", outdir)


# ------------------------------------------------------- stage02 check teeth --
@seed("molpct_sum_offender", "n_pass",
      "Perturb one printed molar percentage and require C1 to lose a pass. A "
      "shift-invariant check scores the same clean and seeded.")
def s3(outdir):
    return run(PRELUDE + """
import json, os, shutil
od = os.environ['KK2015_OUTDIR']
import stage01_extract as s1
P = {"table_s1": s1.parse_table_s1(s1.si_pages()[1])}
pages = s1.si_pages()
des, fx = s1.parse_designs(pages[9], pages[10])
P["designs"] = des; P["table_s2"] = s1.parse_table_s2(pages[2])
P["table_s3"] = s1.parse_table_s3(pages[3])
import pdfplumber
with pdfplumber.open(s1.SI) as pdf:
    P["table_s4"], _ = s1.parse_table_s4(pdf.pages[4])
P["emdash_fixed"] = fx
base = None
for mode in ("clean", "seeded"):
    Q = json.loads(json.dumps(P))
    if mode == "seeded":
        for r in Q["table_s1"]:
            if r["code"] == "B-20":
                r["chol_molpct"] = r["chol_molpct"] + 1.0
    json.dump(Q, open(os.path.join(od, "_parsed.json"), "w"))
    for m in list(sys.modules):
        if m.startswith("stage02"): del sys.modules[m]
    import stage02_checks as c
    r = c.check_molpct_sums()
    print("METRIC %s n_pass=%d max_dev=%.1f offenders=%s"
          % (mode, r["n_pass"], r["max_abs_deviation_pct"], r["offenders"]))
    if mode == "clean": base = r["n_pass"]
    else: sys.exit(0 if r["n_pass"] == base - 1 else 1)
""", outdir)


@seed("cross_parameterisation_transpose", "n_exact",
      "Swap two formulations' weight ratios in the DESIGN table only. The two "
      "printed parameterisations must stop agreeing.")
def s4(outdir):
    return run(PRELUDE + """
import json, os
od = os.environ['KK2015_OUTDIR']
import stage01_extract as s1, pdfplumber
pages = s1.si_pages()
des, fx = s1.parse_designs(pages[9], pages[10])
P = {"table_s1": s1.parse_table_s1(pages[1]), "designs": des,
     "table_s2": s1.parse_table_s2(pages[2]), "table_s3": s1.parse_table_s3(pages[3]),
     "emdash_fixed": fx}
with pdfplumber.open(s1.SI) as pdf:
    P["table_s4"], _ = s1.parse_table_s4(pdf.pages[4])
base = None
for mode in ("clean", "seeded"):
    Q = json.loads(json.dumps(P))
    if mode == "seeded":
        a, b = Q["designs"]["B-15"], Q["designs"]["B-16"]
        a["wt_ratio"], b["wt_ratio"] = b["wt_ratio"], a["wt_ratio"]
    json.dump(Q, open(os.path.join(od, "_parsed.json"), "w"))
    for m in list(sys.modules):
        if m.startswith("stage02"): del sys.modules[m]
    import stage02_checks as c
    r = c.check_cross_parameterisation()
    print("METRIC %s n_exact=%d/%d mismatches=%d"
          % (mode, r["n_exact"], r["n_cells_compared"], len(r["mismatches"])))
    if mode == "clean": base = r["n_exact"]
    else: sys.exit(0 if r["n_exact"] == base - 2 else 1)
""", outdir)


@seed("coded_level_bijection", "n_consistent",
      "Corrupt one character of one coded design string. The third "
      "parameterisation must lose its bijection.")
def s5(outdir):
    return run(PRELUDE + """
import json, os
od = os.environ['KK2015_OUTDIR']
import stage01_extract as s1, pdfplumber
pages = s1.si_pages()
des, fx = s1.parse_designs(pages[9], pages[10])
P = {"table_s1": s1.parse_table_s1(pages[1]), "designs": des,
     "table_s2": s1.parse_table_s2(pages[2]), "table_s3": s1.parse_table_s3(pages[3]),
     "emdash_fixed": fx}
with pdfplumber.open(s1.SI) as pdf:
    P["table_s4"], _ = s1.parse_table_s4(pdf.pages[4])
base = None
for mode in ("clean", "seeded"):
    Q = json.loads(json.dumps(P))
    if mode == "seeded":
        d = Q["designs"]["A-01"]["design"]
        assert d[3] == "+", d
        Q["designs"]["A-01"]["design"] = d[:3] + "-" + d[4:]
    json.dump(Q, open(os.path.join(od, "_parsed.json"), "w"))
    for m in list(sys.modules):
        if m.startswith("stage02"): del sys.modules[m]
    import stage02_checks as c
    r = c.check_coded_levels()
    print("METRIC %s n_consistent=%d/%d all_bijective=%s"
          % (mode, r["n_consistent"], r["n_character_assertions"], r["all_bijective"]))
    if mode == "clean": base = r["n_consistent"]
    else: sys.exit(0 if r["n_consistent"] < base else 1)
""", outdir)


@seed("table_s2_statistic", "mean_abs_error",
      "Perturb one Library B EPO value. The printed Mean=2961.6 and the printed "
      "DOPE estimate=2307.1 must stop reproducing.")
def s6(outdir):
    return run(PRELUDE + """
import json, os
od = os.environ['KK2015_OUTDIR']
import stage01_extract as s1, pdfplumber
pages = s1.si_pages()
des, fx = s1.parse_designs(pages[9], pages[10])
P = {"table_s1": s1.parse_table_s1(pages[1]), "designs": des,
     "table_s2": s1.parse_table_s2(pages[2]), "table_s3": s1.parse_table_s3(pages[3]),
     "emdash_fixed": fx}
with pdfplumber.open(s1.SI) as pdf:
    P["table_s4"], _ = s1.parse_table_s4(pdf.pages[4])
base = None
for mode in ("clean", "seeded"):
    Q = json.loads(json.dumps(P))
    if mode == "seeded":
        for r in Q["table_s1"]:
            if r["code"] == "B-26": r["epo_mean"] = 1000.0
    json.dump(Q, open(os.path.join(od, "_parsed.json"), "w"))
    for m in list(sys.modules):
        if m.startswith("stage02"): del sys.modules[m]
    import stage02_checks as c
    r = c.check_table_s2_statistics()
    print("METRIC %s mean_abs_error=%.4f dope_abs=%.1f"
          % (mode, r["mean_abs_error"], r["half_diff_vs_printed_estimate_abs"]))
    if mode == "clean": base = r["mean_abs_error"]
    else: sys.exit(0 if r["mean_abs_error"] > base + 100 else 1)
""", outdir)


@seed("table_s4_geometry", "refusal",
      "Delete one anchor token from Table S4's word list. The coordinate-based "
      "dose assignment must refuse, not fall back to reading order.")
def s7(outdir):
    return run(PRELUDE + """
import os
import stage01_extract as s1, pdfplumber
class FakePage:
    def __init__(self, pg, drop): self._pg, self._drop = pg, drop
    def extract_words(self):
        return [w for w in self._pg.extract_words() if w["text"] != self._drop]
    def extract_text(self): return self._pg.extract_text()
with pdfplumber.open(s1.SI) as pdf:
    pg = pdf.pages[4]
    rows, geom = s1.parse_table_s4(pg)
    print("METRIC clean rows=%d bands=%d" % (len(rows), geom["n_dose_bands"]))
    try:
        s1.parse_table_s4(FakePage(pg, "DOPE"))
        print("METRIC seeded refusal=None")
        raise SystemExit(1)
    except ValueError as e:
        print("METRIC seeded refusal=ValueError %s" % str(e)[:70])
        raise SystemExit(0)
""", outdir)


@seed("emdash_width_guard", "refusal",
      "Hand the em-dash normaliser a string it cannot uniquely widen. It must "
      "refuse rather than pad to the expected width.")
def s8(outdir):
    return run(PRELUDE + """
import stage01_extract as s1
ok, was = s1._fix_emdash("+0\\u20140", 5)
print("METRIC clean expanded=%r changed=%s" % (ok, was))
bad = 0
for arg, want in (("+00", 5), ("+0\\u2014\\u20140", 5), ("+0\\u20140", 4)):
    try:
        s1._fix_emdash(arg, want)
    except ValueError:
        bad += 1
print("METRIC seeded refusals=%d/3" % bad)
raise SystemExit(0 if bad == 3 else 1)
""", outdir)


@seed("make_reference_ast_guard", "findings",
      "F-037. The anchor's independence guard must fire when the module reaches "
      "a banned module, and must NOT fire on its own ban list or docstring.")
def s9(outdir):
    p = subprocess.run([sys.executable, "-I",
                        os.path.join(HERE, "make_reference.py"),
                        "--seed-copied-template"],
                       cwd=HERE, capture_output=True, text=True,
                       env=dict(os.environ, KK2015_OUTDIR=outdir))
    if "BOTH GUARDS FIRED" in p.stdout:
        p.stdout += "\nMETRIC guards=BOTH_FIRED\n"
    return p


@seed("make_reference_stereo", "assertion",
      "The derived glycerol-C2=(R) assertion must fire when the backbone tag is "
      "flipped -- otherwise the stereo claim is vacuous.")
def s10(outdir):
    p = subprocess.run([sys.executable, "-I",
                        os.path.join(HERE, "make_reference.py"),
                        "--seed-wrong-stereo"],
                       cwd=HERE, capture_output=True, text=True,
                       env=dict(os.environ, KK2015_OUTDIR=outdir))
    if "ASSERTION FIRED" in p.stdout:
        p.stdout += "\nMETRIC stereo=ASSERTION_FIRED\n"
    return p


@seed("three_way_delta6", "all_three_agree",
      "A delta-6 oleoyl passes every internal assertion at IDENTICAL formula and "
      "IDENTICAL exact mass. Only the three-way InChIKey check can see it.")
def s11(outdir):
    p = subprocess.run([sys.executable, "-I",
                        os.path.join(HERE, "stage04_structures.py"),
                        "--seed-delta6"],
                       cwd=HERE, capture_output=True, text=True,
                       env=dict(os.environ, KK2015_OUTDIR=outdir))
    if "THREE-WAY CHECK FIRED" in p.stdout:
        p.stdout += "\nMETRIC three_way=FIRED\n"
    return p


@seed("make_reference_regioisomer", "refusal",
      "A 1,3-diacyl-2-phospho connectivity defect at delta-mass 0.000000 Da must "
      "be refused by the assembly route, proving the gate is not a mass check.")
def s12(outdir):
    p = subprocess.run([sys.executable, "-I",
                        os.path.join(HERE, "make_reference.py"),
                        "--seed-regioisomer"],
                       cwd=HERE, capture_output=True, text=True,
                       env=dict(os.environ, KK2015_OUTDIR=outdir))
    if "CONNECTIVITY DEFECT REFUSED" in p.stdout:
        p.stdout += "\nMETRIC regioisomer=REFUSED\n"
    return p


def digest(folder):
    import hashlib
    out = {}
    for f in TRACKED:
        p = os.path.join(folder, f)
        out[f] = (hashlib.sha256(open(p, "rb").read()).hexdigest()
                  if os.path.exists(p) else None)
    return out


def main():
    before = digest(HERE)
    results = []
    tmproot = tempfile.mkdtemp(prefix="kk2015_seeds_")
    try:
        for name, metric, why, fn in SEEDS:
            od = os.path.join(tmproot, name)
            os.makedirs(od, exist_ok=True)
            p = fn(od)
            lines = [l for l in (p.stdout or "").split("\n")
                     if l.startswith("METRIC")]
            credited = bool(lines) and p.returncode == 0
            results.append({"seed": name, "expected_metric": metric,
                            "returncode": p.returncode,
                            "metric_lines": lines,
                            "credited": credited,
                            "why": why,
                            "stderr_tail": (p.stderr or "").strip()[-200:]})
            print("%-32s rc=%-3d credited=%-5s %s"
                  % (name, p.returncode, credited,
                     " | ".join(l[7:] for l in lines)[:150]))
            if not credited:
                print("    NOT CREDITED -- %s"
                      % ("no METRIC line: the child may have died before "
                         "reaching its metric" if not lines
                         else "non-zero returncode"))
                if p.stderr:
                    print("    stderr: %s" % p.stderr.strip()[-300:])
    finally:
        shutil.rmtree(tmproot, ignore_errors=True)

    after = digest(HERE)
    untouched = before == after
    print()
    print("tracked artifacts byte-identical across the seed suite: %s" % untouched)
    if not untouched:
        for k in TRACKED:
            if before.get(k) != after.get(k):
                print("   CHANGED BY A SEED: %s" % k)

    n_ok = sum(1 for r in results if r["credited"])
    summary = {"n_seeds": len(results), "n_credited": n_ok,
               "all_credited": n_ok == len(results),
               "curation_folder_untouched_by_seeds": untouched,
               "seeds": results}
    with open(os.path.join(HERE, "_seeds.json"), "w") as fh:
        json.dump(summary, fh, indent=1, sort_keys=True, default=str)
    print("%d/%d seeds credited" % (n_ok, len(results)))
    return 0 if (n_ok == len(results) and untouched) else 1


if __name__ == "__main__":
    sys.exit(main())
