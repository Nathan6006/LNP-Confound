#!/usr/bin/env python3
"""Seeded-mutation harness for zhangn_42569578 (MY_2026).

RULES THIS HARNESS ENFORCES, each of which cost someone a finding:

  F-069  Every run writes to a TEMP OUTDIR. Never the curation folder. A
         sibling's seed copies wrote into the real folder and left a seeded
         regioisomer anchor behind; the run exited non-zero and the non-zero
         exit CONCEALED the write.
  A seed that only RAISES is not credited. A child that dies before its metric
         line is reached proves nothing about the metric, so a seed whose
         expected METRIC line is absent is reported UNCREDITED, not passed.
  Verdicts come from returncode of ONE CHILD PROCESS PER SEED. `python3 ... |
         tail` returns tail's exit code, so no pipes.
  Every seed is paired with a CLEAN CONTROL that must come back clean.

Run: python3 -I test_seeds.py
"""
import os
import re
import shutil
import subprocess
import sys
import tempfile
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
PY = [sys.executable, "-I"]


def run(argv, cwd=HERE):
    p = subprocess.run(argv, cwd=cwd, capture_output=True, text=True)
    return p.returncode, p.stdout + p.stderr


def metric(out, name):
    m = re.search(r"(?:METRIC )?%s=([^\s]+)" % re.escape(name), out)
    return m.group(1) if m else None


class Suite:
    def __init__(self):
        self.rows = []

    def record(self, name, kind, expect, rc, out, mname, mval_ok):
        got = metric(out, mname) if mname else "n/a"
        if mname and got is None:
            verdict = "UNCREDITED (metric line absent - the child died before "
            verdict += "reaching the check, so the check is not what fired)"
            ok = False
        else:
            ok = expect(rc, got)
            verdict = "PASS" if ok else "FAIL"
        self.rows.append((name, kind, rc, mname, got, verdict))
        return ok

    def report(self):
        print("\n%-44s %-8s %-4s %-34s %-10s %s"
              % ("seed/control", "kind", "rc", "metric", "value", "verdict"))
        print("-" * 135)
        bad = 0
        for r in self.rows:
            print("%-44s %-8s %-4s %-34s %-10s %s" % r)
            if not r[5].startswith("PASS"):
                bad += 1
        print("-" * 135)
        print("%d checks, %d not passing" % (len(self.rows), bad))
        return bad


def mutate_xlsx(src, dst, mutator):
    """Rewrite sheet1.xml of a COPY. Never touches `src`."""
    with zipfile.ZipFile(src) as z:
        names = z.namelist()
        data = {n: z.read(n) for n in names}
    data["xl/worksheets/sheet1.xml"] = mutator(
        data["xl/worksheets/sheet1.xml"].decode("utf8")).encode("utf8")
    with zipfile.ZipFile(dst, "w", zipfile.ZIP_DEFLATED) as z:
        for n in names:
            z.writestr(n, data[n])


def main():
    s = Suite()
    tmp = tempfile.mkdtemp(prefix="my2026_seeds_")
    print("TEMP OUTDIR:", tmp)
    src_x = os.path.join(HERE, "source", "mmc2.xlsx")
    before = sorted(os.listdir(HERE))

    # ---------------- CONTROL: clean stage01 ------------------------------
    o = os.path.join(tmp, "ctrl")
    rc, out = run(PY + ["stage01_stage.py", "--outdir", o])
    s.record("CONTROL stage01 clean", "control",
             lambda rc, g: rc == 0 and g == "0", rc, out,
             "design_level_violations_xlsx", None)
    s.record("CONTROL stage01 clean (bijection)", "control",
             lambda rc, g: rc == 0 and g == "True", rc, out,
             "S2_to_xlsx_bijection", None)

    # ---------------- SEED 1: an off-design factor level ------------------
    bad_x = os.path.join(tmp, "seed_level.xlsx")
    mutate_xlsx(src_x, bad_x,
                lambda t: t.replace('<c r="E2"><v>60</v></c>',
                                    '<c r="E2"><v>55</v></c>'))
    o = os.path.join(tmp, "seed1")
    rc, out = run(PY + ["stage01_stage.py", "--outdir", o, "--xlsx", bad_x])
    s.record("SEED ILR 60 -> 55 (off Table 1 levels)", "seed",
             lambda rc, g: g not in (None, "0"), rc, out,
             "design_level_violations_xlsx", None)

    # ---------------- SEED 2: break the S2<->xlsx bijection ---------------
    bad_x2 = os.path.join(tmp, "seed_bij.xlsx")
    mutate_xlsx(src_x, bad_x2,
                lambda t: t.replace('<c r="D2"><v>2.5</v></c>',
                                    '<c r="D2"><v>10</v></c>'))
    o = os.path.join(tmp, "seed2")
    rc, out = run(PY + ["stage01_stage.py", "--outdir", o, "--xlsx", bad_x2])
    s.record("SEED composition changed -> bijection breaks", "seed",
             lambda rc, g: g == "False", rc, out, "S2_to_xlsx_bijection", None)

    # ---------------- SEED 3: drop a formulation row ----------------------
    bad_x3 = os.path.join(tmp, "seed_drop.xlsx")
    mutate_xlsx(src_x, bad_x3,
                lambda t: re.sub(r'<row r="26".*?</row>', "", t, flags=re.S))
    o = os.path.join(tmp, "seed3")
    rc, out = run(PY + ["stage01_stage.py", "--outdir", o, "--xlsx", bad_x3])
    s.record("SEED one formulation row deleted", "seed",
             lambda rc, g: rc != 0, rc, out, None, None)

    # ---------------- CONTROL + SEEDS on make_reference -------------------
    o = os.path.join(tmp, "ref_ctrl")
    rc, out = run(PY + ["make_reference.py", "--outdir", o])
    s.record("CONTROL make_reference clean", "control",
             lambda rc, g: rc == 0 and g == "0", rc, out,
             "GUARD_A_ast_independence_violations", None)

    # SEED 4: guard A -- an import of the enumerator's own helper
    cp = os.path.join(tmp, "mr_astseed.py")
    txt = open(os.path.join(HERE, "make_reference.py"), encoding="utf8").read()
    open(cp, "w", encoding="utf8").write(
        txt.replace("import ast\n", "import ast\nimport corpus_components\n"))
    o = os.path.join(tmp, "seed4")
    rc, out = run(PY + ["make_reference.py", "--outdir", o, "--self", cp])
    s.record("SEED import corpus_components (guard A)", "seed",
             lambda rc, g: rc != 0 and g not in (None, "0"), rc, out,
             "GUARD_A_ast_independence_violations", None)

    # SEED 5: guard B -- a long literal lifted out of an enumerated SMILES
    import csv as _csv
    with open(os.path.join(HERE, "03_enumerated.csv"), newline="",
              encoding="utf8") as f:
        smis = [r["SMILES"] for r in _csv.DictReader(f)]
    frag = smis[0][:40]
    cp2 = os.path.join(tmp, "mr_tmplseed.py")
    open(cp2, "w", encoding="utf8").write(
        txt.replace("HEADER = [", "_STOLEN = %r\nHEADER = [" % frag))
    o = os.path.join(tmp, "seed5")
    rc, out = run(PY + ["make_reference.py", "--outdir", o, "--self", cp2])
    s.record("SEED 40-char literal from an emitted SMILES (guard B)", "seed",
             lambda rc, g: rc != 0 and g not in (None, "0"), rc, out,
             "GUARD_B_template_provenance_violations", None)

    # ---------------- emit.emit_rows: the two refusals --------------------
    probe = os.path.join(tmp, "emit_probe.py")
    open(probe, "w", encoding="utf8").write('''
import sys, os
ROOT = sys.argv[1]; mode = sys.argv[2]
sys.path.insert(0, os.path.join(ROOT, "curation", "lib"))
import emit
items = ["F1"]; recs = {"F1": {"why": "a record"}}
try:
    if mode == "literal_confidence":
        emit.emit_rows(items, recs, confidence="high",
                       note=lambda i, r: "n", label="probe")
    elif mode == "caveat_no_note":
        emit.emit_rows(items, recs, confidence=lambda i, r: "medium",
                       note=lambda i, r: "", label="probe")
    elif mode == "clean":
        rows = emit.emit_rows(items, recs,
                              confidence=lambda i, r: "medium",
                              note=lambda i, r: "a real sentence",
                              label="probe")
        print("EMIT_OK=1 caveats=%d" % emit.caveat_count(rows))
        sys.exit(0)
    print("EMIT_REFUSED=0")
    sys.exit(0)
except Exception as e:
    print("EMIT_REFUSED=1 %s" % type(e).__name__)
    sys.exit(3)
''')
    root = os.path.dirname(os.path.dirname(HERE))
    rc, out = run(PY + [probe, root, "clean"], cwd=tmp)
    s.record("CONTROL emit_rows accepts a proper caveat+note", "control",
             lambda rc, g: rc == 0 and g == "1", rc, out, "EMIT_OK", None)
    rc, out = run(PY + [probe, root, "literal_confidence"], cwd=tmp)
    s.record("SEED emit_rows literal confidence='high' (F-011)", "seed",
             lambda rc, g: rc == 3 and g == "1", rc, out, "EMIT_REFUSED", None)
    rc, out = run(PY + [probe, root, "caveat_no_note"], cwd=tmp)
    s.record("SEED emit_rows caveat with empty note (F-024)", "seed",
             lambda rc, g: rc == 3 and g == "1", rc, out, "EMIT_REFUSED", None)

    # ---------------- F-069: the folder must be untouched ------------------
    after = sorted(os.listdir(HERE))
    s.record("F-069 curation folder unchanged by the seed suite", "control",
             lambda rc, g: before == after, 0,
             "UNCHANGED=%s" % (before == after), "UNCHANGED", None)

    # ---------------- determinism AFTER the seeds --------------------------
    import hashlib
    r1 = os.path.join(tmp, "det1")
    r2 = os.path.join(tmp, "det2")
    run(PY + ["stage01_stage.py", "--outdir", r1])
    run(PY + ["stage01_stage.py", "--outdir", r2])
    run(PY + ["make_reference.py", "--outdir", r1])
    run(PY + ["make_reference.py", "--outdir", r2])
    tracked = ["01_blocks.csv", "03_enumerated.csv", "05_values.csv",
               "06_rows.csv", "07_qc_values.json", "reference_structures.csv"]
    ident = []
    for t in tracked:
        h = []
        for r in (r1, r2):
            p = os.path.join(r, t)
            h.append(hashlib.sha256(open(p, "rb").read()).hexdigest()
                     if os.path.exists(p) else "MISSING")
        ident.append(h[0] == h[1] and h[0] != "MISSING")
    s.record("two clean runs byte-identical on %d artifacts" % len(tracked),
             "control", lambda rc, g: g == "True", 0,
             "IDENTICAL=%s" % all(ident), "IDENTICAL", None)
    # and identical to what is committed in the folder
    same_as_folder = []
    for t in tracked:
        a = os.path.join(HERE, t)
        b = os.path.join(r1, t)
        if os.path.exists(a) and os.path.exists(b):
            same_as_folder.append(
                hashlib.sha256(open(a, "rb").read()).hexdigest() ==
                hashlib.sha256(open(b, "rb").read()).hexdigest())
    s.record("folder artifacts reproduce from a clean run", "control",
             lambda rc, g: g == "True", 0,
             "REPRO=%s" % all(same_as_folder), "REPRO", None)

    bad = s.report()
    print("\nTEMP OUTDIR (left for inspection):", tmp)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
