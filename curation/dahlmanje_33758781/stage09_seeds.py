#!/usr/bin/env python3
"""Seeded-mutation harness for stage01_enumerate.py's gates.

Run:  python3 -I stage09_seeds.py

A CHECK THAT CANNOT FAIL IS NOT A CHECK, AND A SEED THAT ONLY RAISES IS NOT
CREDITED. Every seed below is applied to a COPY of stage01_enumerate.py in a
temp directory and run in its OWN child process, so the verdict comes from
that child's returncode and never from a pipeline's. A seed is credited only
if the child actually REACHED the metric lines -- a seed that dies of a
TypeError on the way there proves nothing, so the harness refuses to credit a
run whose METRIC lines are absent (the shape a sibling got caught by today).

Three of the six seeds are deliberately expected NOT to fire. Those are the
point: they measure where this HRMS check is blind, and a blind spot you have
measured is worth more than a gate you have assumed. The regiochemistry seed
in particular moves the connectivity of the lead compound A-11 at delta mass
EXACTLY 0.000000 Da and delta formula EXACTLY zero.
"""
import os, re, shutil, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "stage01_enumerate.py")
SRC3 = os.path.join(HERE, "stage03_figure_checks.py")
SRC4 = os.path.join(HERE, "stage04_fig1c_and_yields.py")

# The ten products, correct build: adamantylacetyl on the PRIMARY (C3) OH of
# the glycerol, fatty acyl on the secondary (C2) OH -- what the route says.
# The regioisomer swaps them, which is how A-10 and A-11 are NAMED.
REGIO_SWAPPED_C4 = ('"[C@@H](COP(=O)([O-])OCC[N+](C)(C)C)'
                    '(OC(=O)CC12CC3CC(CC(C3)C1)C2)CO"')

SEEDS = [
    # (id, what it changes, find, replace, expectation)
    ("S1_chain_plus_one",
     "A-11's acyl donor becomes tridecanoic instead of lauric: the formula "
     "and the exact mass both move, so the HRMS gate MUST fire.",
     '("A-11",     "lauric acid",        "OC(=O)CCCCCCCCCCC",',
     '("A-11",     "lauric acid",        "OC(=O)CCCCCCCCCCCC",',
     "FIRE"),
    ("S2_enantiomer",
     "Invert the single glycerol stereocentre of intermediate 4 (R -> S). "
     "Formula and exact mass are IDENTICAL; only the CIP assertion can see "
     "it. This is the F-013 shape, reproduced from scratch on this paper.",
     'COMPOUND_4 = ("[C@@H](COP',
     'COMPOUND_4 = ("[C@H](COP',
     "FIRE"),
    ("S3_adduct_swapped",
     "Score A-9 as [M+H]+ while its printed calcd is the [M+Na]+ value: a "
     "21.98 Da error, the mass-convention-drift class (trap 2). MUST fire.",
     '"C30H54NNaO8P", "[M+Na]+", 610.3479, 610.3468)',
     '"C30H54NNaO8P", "[M+H]+", 610.3479, 610.3468)',
     "FIRE"),
    ("S4_formula_typo_unnamed",
     "Delete the Na from A-17's printed formula, the same defect A-16 really "
     "has but on a compound NOT named in the allow-list. MUST fire -- this "
     "proves the A-16 exemption is a named exception and not a blanket skip.",
     '"C38H70NNaO8P", "[M+Na]+", 722.4731, 722.4723)',
     '"C38H70NO8P", "[M+Na]+", 722.4731, 722.4723)',
     "FIRE"),
    ("S5_regiochemistry_IDENTICAL_MASS",
     "Swap the adamantylacetyl and the fatty acyl between the glycerol's "
     "primary and secondary oxygens, for ALL TEN compounds. Formula, exact "
     "mass and the CIP count are all unchanged. EXPECTED BLIND: this is the "
     "measurement that justifies resolving A-10/A-11's names from the "
     "synthesis rather than from HRMS.",
     'COMPOUND_4 = ("[C@@H](COP(=O)([O-])OCC[N+](C)(C)C)(O)"\n'
     '              "COC(=O)CC12CC3CC(CC(C3)C1)C2")',
     'COMPOUND_4 = ' + REGIO_SWAPPED_C4,
     "BLIND"),
    ("S6_adamantyl_to_noradamantyl_control",
     "Replace adamantyl with a CH2-shorter cage (noradamantyl) in compound 4: "
     "the mass moves by one CH2. A gross control, credited only because it "
     "fires for the RIGHT reason -- it confirms S5's silence is about "
     "connectivity and not about the gate being switched off.",
     'COC(=O)CC12CC3CC(CC(C3)C1)C2")',
     'COC(=O)CC12CC3CC(C1)CC3C2")',
     "FIRE"),
]

SEEDS3 = [
    ("S7_ratio_sum_broken",
     "Change one Fig 2b numeral (phospholipid 55 -> 54) so formulation 1 sums "
     "to 99. The sum-to-100 check MUST fire.",
     '"Phospholipid":    [55, 45, 35, 25],',
     '"Phospholipid":    [54, 45, 35, 25],',
     "FIRE"),
    ("S8_factorial_broken",
     "Claim 2 sterols instead of 3, so 10x2x4=80 != the paper's 120. A "
     "claimed factorial that does not multiply out MUST fire.",
     "N_PHOSPHOLIPID, N_STEROL, N_RATIO = 10, 3, 4",
     "N_PHOSPHOLIPID, N_STEROL, N_RATIO = 10, 2, 4",
     "FIRE"),
    ("S9_teal_mask_swallows_grey",
     "Loosen the teal mask (g > r + 18 becomes g > r - 40) so grey bars are "
     "read as enriched. The caption's five MUST stop matching.",
     "teal = (g > r + 18) & (g > b + 8) & (g < 235) & (g > 120)",
     "teal = (g > r - 40) & (g > b - 40) & (g < 235) & (g > 120)",
     "FIRE"),
    ("S10_slots_shifted_half_pitch",
     "Shift every panel-(h) bar slot by +11 px, half the 22.8 px pitch, so "
     "each window straddles two bars. This measures the SENSITIVITY LIMIT of "
     "the slot positions; whatever it does is recorded, not assumed.",
     "PANEL_H_CENTRES = [276, 300, 321, 344, 367, 390, 412, 435, 458, 481]",
     "PANEL_H_CENTRES = [287, 311, 332, 355, 378, 401, 423, 446, 469, 492]",
     "FIRE"),
]

TMP = [None]          # set by main() before any prepare()

METRIC = re.compile(r"^METRIC ([\w.]+)=(.*)$", re.M)


LIBDIR = os.path.abspath(os.path.join(HERE, "..", "lib"))
PATH_LINE = ('sys.path.insert(0, os.path.join(os.path.dirname('
             'os.path.abspath(__file__)),\n                                "..", "lib"))')


def prepare(src):
    """A copy that can run OUTSIDE the curation tree.

    `python3 -I` ignores PYTHONPATH on purpose, so the copy cannot inherit the
    library path from the environment; it has to be written in. This rewrite
    is applied to the control and to every seed IDENTICALLY, so it cannot be
    what makes a seed differ from the control.
    """
    if PATH_LINE in src:
        return src.replace(PATH_LINE, f"sys.path.insert(0, {LIBDIR!r})", 1)
    if "curation/lib" in src or "\"..\", \"lib\"" in src:
        raise SystemExit("prepare(): the sys.path line moved; fix the harness")
    # stage03 resolves its figure through HERE; a copy in a temp directory
    # must still read the REAL figure, so pin it rather than copy the JPEG.
    src = src.replace('HERE = os.path.dirname(os.path.abspath(__file__))',
                      f'HERE = {HERE!r}', 1)
    # and writes go to the temp directory, never to the real folder
    src = src.replace('OUTDIR = HERE', f'OUTDIR = {{TMP[0]!r}}'.format(TMP=TMP), 1)
    return src          # stage03 imports no project library


def run(path):
    r = subprocess.run([sys.executable, "-I", path],
                       capture_output=True, text=True, cwd=os.path.dirname(path),
                       timeout=600)
    return r.returncode, r.stdout + r.stderr


SEEDS4 = [
    ("S11_fig1c_formula_broken",
     "Corrupt one Fig 1c formula (A-11 C32 -> C33). The 10/10 formula match "
     "MUST fire.",
     '"A-11":    ("C32H58NO8P", 615.5, 2.30),',
     '"A-11":    ("C33H58NO8P", 615.5, 2.30),',
     "FIRE"),
    ("S12_mw_outlier_unnamed",
     "Move A-15's printed MW by +1.0 Da -- the same defect A-13 really has, "
     "but on a compound NOT in the known-typo list. MUST fire, which proves "
     "the A-13 exemption is a NAMED exception and not a widened threshold.",
     '"A-15":    ("C36H66NO8P", 671.5, 3.46),',
     '"A-15":    ("C36H66NO8P", 672.5, 3.46),',
     "FIRE"),
    ("S13_logp_not_monotone",
     "Swap A-12's and A-13's printed LogP so the series is no longer "
     "monotone in chain length. MUST fire.",
     '"A-12":    ("C33H60NO8P", 629.4, 2.63),',
     '"A-12":    ("C33H60NO8P", 629.4, 2.99),',
     "FIRE"),
    ("S14_yield_range_broken",
     "Change A-9's isolated yield from 92% to 90%, so the printed 77-92% "
     "range no longer reproduces. MUST fire -- this is the check that proves "
     "A-9 belongs to the series.",
     '"A-9": 92,',
     '"A-9": 90,',
     "FIRE"),
]

# One table, one loop. Writing the loop a third time is error class 49, and
# this harness is not exempt from it.
TARGETS = [
    ("stage01_enumerate",    SRC,  SEEDS,
     ("formula_mismatches_unexplained", "mass_failures_over_0.6mDa",
      "compound4_cip")),
    ("stage03_figure_checks", SRC3, SEEDS3, ("figure_check_failures",)),
    ("stage04_fig1c",         SRC4, SEEDS4, ("fig1c_failures",)),
]


def control_ok(keys, rc, m):
    if rc != 0 or any(k not in m for k in keys):
        return False, [k for k in keys if k not in m]
    for k in keys:
        v = m[k]
        if k == "compound4_cip":
            if v != "R":
                return False, []
        elif not (v.startswith("0") or v == "[]"):
            return False, []
    return True, []


def main():
    tmp = tempfile.mkdtemp(prefix="zg2020_seeds_")
    TMP[0] = tmp
    results = []
    for tname, src_path, seeds, keys in TARGETS:
        ctrl = os.path.join(tmp, f"control_{tname}.py")
        open(ctrl, "w").write(prepare(open(src_path).read()))
        rc, out = run(ctrl)
        cm = dict(METRIC.findall(out))
        ok, missing = control_ok(keys, rc, cm)
        print(f"\nCONTROL {tname:22s} rc={rc} "
              f"{ {k: cm.get(k) for k in keys} }")
        if missing:
            print(f"   missing metric keys: {missing}")
        print(f"CONTROL {tname:22s} VERDICT="
              f"{'PASS' if ok else 'FAIL -- the harness itself is broken'}")
        results.append((f"CONTROL_{tname}", "PASS" if ok else "FAIL",
                        rc, cm, ""))

        for sid, why, find, repl, expect in seeds:
            src = open(src_path).read()
            if find not in src:
                print(f"{sid:36s} SEED DID NOT APPLY (anchor text absent). "
                      "NOT CREDITED.")
                results.append((sid, "NOT_APPLIED", None, {}, why))
                continue
            mutated = prepare(src.replace(find, repl, 1))
            if mutated == prepare(src):
                print(f"{sid:36s} SEED IS A NO-OP. NOT CREDITED.")
                results.append((sid, "NO_OP", None, {}, why))
                continue
            path = os.path.join(tmp, sid + ".py")
            open(path, "w").write(mutated)
            rc, out = run(path)
            m = dict(METRIC.findall(out))
            # A seed that never reached a metric line proves NOTHING, however
            # loudly it exited.
            reached = any(k in m for k in keys) or "GATE FAILED" in out
            if not reached:
                verdict = "NOT_CREDITED_metric_absent"
            elif expect == "FIRE":
                verdict = "FIRED" if rc != 0 else "DID_NOT_FIRE"
            else:
                verdict = ("BLIND_as_predicted" if rc == 0
                           else "FIRED_unexpectedly")
            print(f"{sid:36s} rc={rc:<3} expect={expect:<6} "
                  f"VERDICT={verdict}")
            for l in out.splitlines():
                if ("GATE FAILED" in l or "expected exactly one R" in l
                        or any(f"METRIC {k}=" in l for k in keys)):
                    print(f"      {l.strip()[:150]}")
            results.append((sid, verdict, rc, m, why))

    shutil.rmtree(tmp, ignore_errors=True)
    bad = [r for r in results
           if r[1] in ("DID_NOT_FIRE", "FIRED_unexpectedly", "NO_OP",
                       "NOT_APPLIED", "NOT_CREDITED_metric_absent", "FAIL")]
    n_seeds = sum(len(t[2]) for t in TARGETS)
    print(f"\n{n_seeds} seeds + {len(TARGETS)} controls; "
          f"{len(bad)} uncredited/unexpected: {[b[0] for b in bad]}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
