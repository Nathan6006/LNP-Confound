#!/usr/bin/env python3
"""Seed every check's own defect and confirm it FIRES. Writes 10_seeded.json.

CLAUDE.md section 6: "A NEW CHECK IS NOT DONE UNTIL IT HAS BEEN MADE TO FAIL
ON PURPOSE." Three incidents paid for that rule -- a re-copied C10 anchor that
could only pass, three ZT_2026 marginal checks that scored a literal 0 for
twelve days while the summary beside them read "ALL FOUR EXACT", and two
versions of C16 that reported "all derived" while examining nothing. Every one
of them printed a clean line.

Each mutation below is applied to a COPY of this folder (never to the folder
itself), the stage is re-run, and the run must fail AND the failure message
must contain the expected fragment. A mutation that fails for the wrong reason
is as bad as one that passes, so both are checked. A CONTROL run of each
mutated stage, unmutated, must still succeed.

Run:  python3 -I 10_seed_checks.py
"""
import json
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
LIB = os.path.join(HERE, "..", "lib")

# (id, stage file, old text, new text, expected fragment of the refusal)
MUTATIONS = [
    ("M1-chain-length", "stage01_chains.py",
     '"PY12": 17,', '"PY12": 19,',
     "geometry and HRMS disagree"),
    # The unsaturation POSITION. stage01 measures it and warns; the refusal
    # lives in make_reference.py, which will not write an anchor for a 9Z
    # oleoyl if the measurement says otherwise. So the mutation goes into
    # stage01 and the stage that must refuse is make_reference.
    ("M2-oleoyl-position", "make_reference.py",
     "    dbl = [n_left + 1, n_left + 2]", "    dbl = [n_left, n_left + 1]",
     "this builder writes 9Z", "stage01_chains.py"),
    ("M3-A1-branches", "stage02_amines.py",
     "    NBOX = (269.4, 234.9, 273.3, 239.7)",
     "    NBOX = (284.8, 234.9, 288.3, 239.7)",
     "bonds in the drawing"),
    ("M4-block-composition", "stage03_enumerate.py",
     '"T18": (18, "[CH3:6]" + "C" * 16),', '"T18": (18, "[CH3:6]" + "C" * 15),',
     "mDa from the SI's own printed calcd"),
    ("M5-lost-cis", "stage03_enumerate.py",
     '"[CH3:6]" + "C" * 6 + "/C=C\\\\" + "C" * 8',
     '"[CH3:6]" + "C" * 6 + "C=C" + "C" * 8',
     "stereo-defined double bonds"),
    ("M6-tick-set", "stage04_figS3.py",
     "AXIS_X - 33, AXIS_X - 25", "AXIS_X - 38, AXIS_X - 4",
     "major y ticks, expected"),
    ("M7-bar-count", "stage04_figS3.py",
     "N_BARS = 47", "N_BARS = 46",
     "bar columns, expected"),
    ("M8-class27", "stage04_figS3.py",
     "cal = axislib.calibrate(ticks, TICK_VALUES, kind=\"linear\",\n"
     "                            forbid_target=te, axis=\"y\")",
     "cal = axislib.calibrate(ticks, ticks, kind=\"linear\",\n"
     "                            forbid_target=ticks, axis=\"y\")",
     "the fit input IS the validation target"),
    ("M9-value-order", "stage05_values.py",
     '("L21", "T16-PY18-A9",   55)', '("L21", "T16-PY18-A9",   35)',
     "not non-increasing"),
    ("M10-unknown-label", "stage05_values.py",
     '("L24", "T18-PY12-A6",  52)', '("L24", "T18-PY12-A99", 52)',
     "not one of the 288 enumerated labels"),
    ("M11-literal-confidence", "stage06_emit.py",
     "quality = emit.emit_rows(\n        items, records, confidence=confidence,",
     "quality = emit.emit_rows(\n        items, records, confidence=\"high\",",
     "is a LITERAL"),
    ("M12-caveat-no-note", "stage06_emit.py",
     '        if rec["amine"] == "A5":\n            return ("A5\'s',
     '        if rec["amine"] == "A5":\n            return ("" or ("A5\'s',
     None),      # resolved below; see _m12
    ("M13-anchor-is-a-copy", "make_reference.py",
     "    built = (f[\"diethylamino_head\"] + f[\"ethylene_linker\"]",
     "    built = enum_copy() or (f[\"diethylamino_head\"] + f[\"ethylene_linker\"]",
     "byte-identical to the enumeration"),
    ("M14-wrong-anchor", "make_reference.py",
     '"hexadecyl": "C" * n_arm,', '"hexadecyl": "C" * (n_arm - 1),',
     "DIFFERENT molecules"),
]

# M12 needs a real empty note on a caveat row, which a one-line swap cannot
# express cleanly; it is written out in full here.
_M12_OLD = '''    def note(i, rec):
        if rec["structure_determined"] != "True":
            return ""'''
_M12_NEW = '''    def note(i, rec):
        return ""
        if rec["structure_determined"] != "True":
            return ""'''

_M13_HELPER = '''
def enum_copy():
    import csv as _csv, os as _os
    r = {x["lipidoid"]: x for x in _csv.DictReader(
        open(_os.path.join(HERE, "03_enumerated.csv")))}["T18U-PY16-A1"]
    return r["IL_SMILES"]
'''

ORDER = ["stage01_chains.py", "stage02_amines.py", "stage03_enumerate.py",
         "stage05_values.py", "stage04_figS3.py", "stage06_emit.py",
         "make_reference.py"]


def prepare(tmp):
    root = os.path.join(tmp, "curation")
    os.makedirs(root)
    shutil.copytree(LIB, os.path.join(root, "lib"))
    # curation/lib/emit.py resolves provenance/schema_quality.json relative to
    # itself, so the copied tree needs that sibling or every emit test fails
    # for a path reason instead of the reason it is testing.
    os.symlink(os.path.join(HERE, "..", "..", "provenance"),
               os.path.join(tmp, "provenance"))
    dst = os.path.join(root, "paper")
    os.makedirs(dst)
    for f in os.listdir(HERE):
        src = os.path.join(HERE, f)
        if os.path.isdir(src):
            os.symlink(src, os.path.join(dst, f))
        else:
            shutil.copy2(src, dst)
    return dst


def run(dst, stage):
    p = subprocess.run([sys.executable, "-I", stage], cwd=dst,
                       capture_output=True, text=True)
    return p.returncode, (p.stdout + p.stderr)


def build_all(dst, upto):
    """Run the pipeline in order so each stage has its inputs."""
    for st in ORDER:
        rc, out = run(dst, st)
        if st == upto:
            return rc, out
        if rc != 0:
            return rc, f"[prerequisite {st} failed] " + out
    raise AssertionError(upto)


def main():
    results = []
    for mut in MUTATIONS:
        mid, stage, old, new, want = mut[:5]
        mutate_in = mut[5] if len(mut) > 5 else stage
        with tempfile.TemporaryDirectory() as tmp:
            dst = prepare(tmp)
            path = os.path.join(dst, mutate_in)
            src = open(path).read()
            if mid == "M12-caveat-no-note":
                old, new, want = _M12_OLD, _M12_NEW, "curation_note is empty"
            if mid == "M13-anchor-is-a-copy":
                src = src.replace("def main():", _M13_HELPER + "\n\ndef main():", 1)
            if old not in src:
                results.append({"id": mid, "stage": stage, "FIRED": False,
                                "why": "the mutation's anchor text is no "
                                       "longer in the file -- the test is "
                                       "stale, which is itself a finding",
                                "mutate_in": mutate_in})
                continue
            open(path, "w").write(src.replace(old, new, 1))
            rc, out = build_all(dst, stage)
            fired = rc != 0 and (want in out)
            results.append({
                "id": mid, "stage": stage, "mutate_in": mutate_in,
                "FIRED": fired,
                "exit_code": rc,
                "expected_fragment": want,
                "message_head": out.strip().splitlines()[-1][:220]
                if out.strip() else ""})

    # CONTROL: the unmutated pipeline must run clean end to end.
    with tempfile.TemporaryDirectory() as tmp:
        dst = prepare(tmp)
        control = []
        for st in ORDER:
            rc, out = run(dst, st)
            control.append({"stage": st, "exit_code": rc,
                            "tail": out.strip().splitlines()[-1][:160]
                            if out.strip() else ""})

    report = {
        "_what": "Every check in this folder, seeded with the exact defect it "
                 "exists to catch. A check that has not been made to fail is "
                 "not a check (CLAUDE.md section 6).",
        "mutations": results,
        "n_mutations": len(results),
        "n_fired": sum(1 for r in results if r["FIRED"]),
        "all_fired": all(r["FIRED"] for r in results),
        "control_run_unmutated": control,
        "control_clean": all(c["exit_code"] == 0 for c in control),
    }
    with open(os.path.join(HERE, "10_seeded.json"), "w") as fh:
        json.dump(report, fh, indent=1)
    for r in results:
        print(("FIRED " if r["FIRED"] else "MISSED"), r["id"],
              "" if r["FIRED"] else r.get("message_head", r.get("why", "")))
    print("control clean:", report["control_clean"])
    return 0 if (report["all_fired"] and report["control_clean"]) else 1


if __name__ == "__main__":
    sys.exit(main())
