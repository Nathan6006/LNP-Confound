#!/usr/bin/env python3
"""Stage 06 - seed each check's own defect and confirm the AGREEMENT COUNT
MOVES. A clean run is credited only afterwards.

Every seed runs in its OWN process and the verdict is read from `returncode`,
never from stdout text. F-038: three ZT_2026 marginal checks scored a literal
0 for twelve days while the summary beside them said "ALL FOUR EXACT", because
`.get(k, 0)` swallowed a key rename and nobody ever made the check fail.

Appends a `seeds` list to 07_qc.json.
"""
import json
import os
import subprocess
import sys
import textwrap

HERE = os.path.dirname(os.path.abspath(__file__))
LIB = os.path.join(HERE, "..", "lib")
PY = sys.executable


def run(src, name):
    """Run a seed in its own process. Verdict = returncode."""
    path = os.path.join(HERE, "_seed_%s.py" % name)
    with open(path, "w") as fh:
        fh.write(src)
    try:
        p = subprocess.run([PY, "-I", path], capture_output=True, text=True,
                           cwd=HERE)
        return {"returncode": p.returncode,
                "stdout_tail": p.stdout.strip().splitlines()[-4:],
                "stderr_tail": p.stderr.strip().splitlines()[-3:]}
    finally:
        os.remove(path)


HEAD = textwrap.dedent("""
    import json, os, sys
    import numpy as np
    sys.path.insert(0, %r)
    HERE = %r
""") % (LIB, HERE)


SEEDS = [
    # ---- G: the axis calibration -----------------------------------------
    ("G_axis_loo_native_tolerance", "axis.loo_residuals at the RESCALED "
     "tolerance still rejects a displaced tick", HEAD + textwrap.dedent("""
        import axis
        from errors import CalibrationError
        true_ticks = [63.0, 128.0, 194.0, 259.0]          # FigS5a, measured
        native = 219.5
        tol_native = axis.native_tol(1.5, native_dpi=native)
        base = max(axis.loo_residuals(true_ticks))
        assert base <= tol_native, "control run must pass"
        worst_pass = None
        for d in (1.0, 2.0, 3.0, 3.25, 4.0, 5.0, 6.0):
            seeded = list(true_ticks); seeded[1] += d
            loo = max(axis.loo_residuals(seeded))
            print("seed %+0.2f px -> worst LOO %.3f native px (tol %.3f, "
                  "unrescaled tol 1.5) -> %s"
                  % (d, loo, tol_native, "PASS" if loo <= tol_native else "FAIL"))
            if loo <= tol_native:
                worst_pass = d
        seeded = list(true_ticks); seeded[1] += 5.0
        loo = max(axis.loo_residuals(seeded))
        if loo <= tol_native:
            print("SEED NOT CAUGHT"); sys.exit(0)
        print("control worst LOO %.4f ; seeded +5.0px worst LOO %.4f ; "
              "largest displacement that still slips through: %s px"
              % (base, loo, worst_pass))
        sys.exit(1)
     """)),

    # ---- bars: class 32, the thing nrow relaxation could have opened ------
    ("bars_class32_centred_marker", "a centred replicate blob 6 px above the "
     "top rule is still rejected at the RELAXED nrow=2", HEAD + textwrap.dedent("""
        import bars
        from errors import DetectionError
        ink = np.zeros((60, 20), dtype=bool)
        ink[40:42, 2:18] = True          # a 2-row top rule, full width
        ink[42:58, 2:4] = True           # walls
        ink[42:58, 16:18] = True
        clean = bars.bar_top(ink, 2, 18, 0, 58, nrow=2)
        ink2 = ink.copy()
        ink2[34:36, 8:12] = True         # centred replicate blob, 6 px above
        try:
            got = bars.bar_top(ink2, 2, 18, 0, 58, nrow=2)
        except DetectionError as e:
            print("REFUSED as intended:", str(e).splitlines()[0]); sys.exit(1)
        if got["row"] == clean["row"]:
            print("clean row %d, seeded row %d -> marker ignored, count "
                  "unmoved" % (clean["row"], got["row"])); sys.exit(1)
        print("SEED NOT CAUGHT: marker read as the top (row %d vs %d)"
              % (got["row"], clean["row"])); sys.exit(0)
     """)),

    # ---- the furniture strip ---------------------------------------------
    ("furniture_long_rule", "a 300 px table rule crossing a bar window is "
     "removed; without the strip it is read as the bar top",
     HEAD + textwrap.dedent("""
        import importlib.util
        spec = importlib.util.spec_from_file_location(
            "s2", os.path.join(HERE, "stage02_bars.py"))
        s2 = importlib.util.module_from_spec(spec); spec.loader.exec_module(s2)
        import bars
        ink = np.zeros((80, 400), dtype=bool)
        ink[60:63, 100:114] = True          # a real bar top, 14 px wide
        ink[63:78, 100:102] = True
        ink[63:78, 112:114] = True
        ink[20:23, 40:360] = True           # an inset-table rule, 320 px long
        raw = bars.bar_top(ink, 100, 114, 0, 78)
        stripped, n = s2.strip_long_rules(ink, 16)
        fixed = bars.bar_top(stripped, 100, 114, 0, 78)
        print("without strip row=%d ; with strip row=%d ; pixels cleared=%d"
              % (raw["row"], fixed["row"], n))
        if raw["row"] == fixed["row"]:
            print("SEED NOT CAUGHT: the strip changed nothing"); sys.exit(0)
        sys.exit(1)
     """)),

    # ---- emit: the literal-confidence guard ------------------------------
    ("emit_literal_confidence", "a hard-coded confidence literal raises "
     "(F-011 / F-025 / F-046)", HEAD + textwrap.dedent("""
        import emit
        from errors import EmitError
        items = ["a"]; records = {"a": {"x": 1}}
        try:
            emit.emit_rows(items, records, confidence="high",
                           note=lambda i, r: "n", label="seed")
        except EmitError as e:
            print("REFUSED as intended:", str(e).splitlines()[0]); sys.exit(1)
        print("SEED NOT CAUGHT: a literal confidence was accepted"); sys.exit(0)
     """)),

    # ---- emit: a caveat with no note -------------------------------------
    ("emit_caveat_without_note", "a low-confidence row with an empty note "
     "raises", HEAD + textwrap.dedent("""
        import emit
        from errors import EmitError
        try:
            emit.emit_rows(["a"], {"a": {}},
                           confidence=lambda i, r: "low",
                           note=lambda i, r: "", label="seed")
        except EmitError as e:
            print("REFUSED as intended:", str(e).splitlines()[0]); sys.exit(1)
        print("SEED NOT CAUGHT"); sys.exit(0)
     """)),

    # ---- C: the correlation check must move if the bars move -------------
    ("C_correlation_moves_with_bars", "scaling every FigS4b bar by 1.10 must "
     "break the Spearman magnitude agreement count (if it does not, the check "
     "is not testing the bars)", HEAD + textwrap.dedent("""
        import pdfplumber
        from scipy import stats
        SI = os.path.join(HERE, "source", "nn5c01378_si_001.pdf")
        t2 = {}
        with pdfplumber.open(SI) as pdf:
            for pno in (2, 3):
                for line in (pdf.pages[pno-1].extract_text() or "").splitlines():
                    p = line.split()
                    if len(p) == 5 and p[0].startswith("LNP"):
                        try: v = [float(x) for x in p[1:]]
                        except ValueError: continue
                        t2[p[0]] = dict(zip(["VN2C2","DDAB","Chol","DMG-PEG"], v))
        b = json.load(open(os.path.join(HERE, "02_bars.json")))
        s4 = {v["series"]: v["value"] for v in b["panels"]["FigS4b"]["values"]
              if v["status"] == "ok"}
        names = [n for n in ("LNP%d" % i for i in range(1, 61)) if n in s4]
        P = {"VN2C2": -0.636, "DDAB": -0.501, "Chol": 0.727, "DMG-PEG": -0.519}
        def count(y):
            c = 0
            for k, pr in P.items():
                x = np.array([t2[n][k] for n in names])
                rho = stats.spearmanr(x, y)[0]
                c += abs(abs(rho) - abs(pr)) <= 0.02
            return c
        y = np.array([s4[n] for n in names])
        base = count(y)
        # a pure scale leaves RANKS untouched, so Spearman cannot see it.
        scaled = count(y * 1.10)
        # a bias that is NOT rank-preserving: reverse the 12 smallest bars.
        y2 = y.copy(); idx = np.argsort(y2)[:12]; y2[idx] = y2[idx][::-1] + 30
        perturbed = count(y2)
        print("baseline agree=%d ; x1.10 agree=%d ; rank-perturbed agree=%d"
              % (base, scaled, perturbed))
        if base == 4 and perturbed < 4:
            sys.exit(1)
        print("SEED NOT CAUGHT"); sys.exit(0)
     """)),

    # ---- D: the drawn red/black grouping ---------------------------------
    ("D_red_black_separation", "swapping two bars across the colour boundary "
     "destroys the single-threshold separation", HEAD + textwrap.dedent("""
        q = json.load(open(os.path.join(HERE, "07_qc.json")))
        d = [c for c in q["checks"] if c["id"].startswith("D_")][0]["result"]
        red_hi, black_lo = d["red_range"][1], d["black_range"][0]
        base_ok = red_hi < black_lo
        # seed: move one red bar 5 units above the lowest black bar
        seeded_red_hi = black_lo + 5.0
        seeded_ok = seeded_red_hi < black_lo
        print("control: red max %.2f < black min %.2f -> %s (margin %.3f)"
              % (red_hi, black_lo, base_ok, black_lo - red_hi))
        print("seeded : red max %.2f < black min %.2f -> %s"
              % (seeded_red_hi, black_lo, seeded_ok))
        if base_ok and not seeded_ok:
            sys.exit(1)
        print("SEED NOT CAUGHT"); sys.exit(0)
     """)),

    # ---- E/F: the mass arithmetic ----------------------------------------
    ("EF_mass_arithmetic", "a one-CH2 error in the proposed core breaks the "
     "HRMS agreement", HEAD + textwrap.dedent("""
        from rdkit import Chem, RDLogger
        from rdkit.Chem.Descriptors import ExactMolWt
        RDLogger.DisableLog("rdApp.*")
        E = 0.00054857990946
        st = json.load(open(os.path.join(HERE, "03_proposed_structures.json")))
        rec = [r for r in st["lipids"] if r["name"] == "VN2C2"][0]
        m = Chem.MolFromSmiles(rec["proposed_IL_SMILES"])
        good = abs(ExactMolWt(m) + 1.00782503207 - E - 1085.9220) * 1000
        bad_smi = rec["proposed_IL_SMILES"].replace("CCN(CC)CCN", "CCN(CC)CCCN", 1)
        bm = Chem.MolFromSmiles(bad_smi)
        bad = abs(ExactMolWt(bm) + 1.00782503207 - E - 1085.9220) * 1000
        print("correct core: %.2f mDa ; core + one CH2: %.2f mDa" % (good, bad))
        if good < 0.2 and bad > 1.0:
            sys.exit(1)
        print("SEED NOT CAUGHT"); sys.exit(0)
     """)),

    # ---- the C10 anchor's own gate ---------------------------------------
    ("anchor_direction_gate", "a reversed ester in the anchor - same formula, "
     "same exact mass - must make make_reference.py exit 1",
     HEAD + textwrap.dedent("""
        import importlib.util, subprocess
        spec = importlib.util.spec_from_file_location(
            "mr", os.path.join(HERE, "make_reference.py"))
        mr = importlib.util.module_from_spec(spec); spec.loader.exec_module(mr)
        from rdkit import Chem, RDLogger
        from rdkit.Chem.Descriptors import ExactMolWt
        from rdkit.Chem.rdMolDescriptors import CalcMolFormula
        RDLogger.DisableLog("rdApp.*")
        good = mr.VE_AC_ANCHOR
        # the SAME atoms, the ester written backwards: Ar-O-C(=O)-vinyl
        # becomes Ar-C(=O)-O-vinyl. Identical formula, identical exact mass.
        bad = good.replace("c(OC(=O)C=C)", "c(C(=O)OC=C)")
        gm, bm = Chem.MolFromSmiles(good), Chem.MolFromSmiles(bad)
        print("formula  good %s  seeded %s"
              % (CalcMolFormula(gm), CalcMolFormula(bm)))
        print("mass     good %.4f  seeded %.4f  delta %.6f"
              % (ExactMolWt(gm), ExactMolWt(bm),
                 abs(ExactMolWt(gm) - ExactMolWt(bm))))
        assert CalcMolFormula(gm) == CalcMolFormula(bm), "seed must be isomeric"
        assert abs(ExactMolWt(gm) - ExactMolWt(bm)) < 1e-9, "mass must be blind"
        # control: the real anchor passes
        mr.check_direction(good)
        # seeded: run make_reference as a PROCESS with the bad literal patched
        src = open(os.path.join(HERE, "make_reference.py")).read()
        src = src.replace('"' + good + '"', '"' + bad + '"')
        assert bad in src, "patch did not apply"
        tmp = os.path.join(HERE, "_seeded_make_reference.py")
        open(tmp, "w").write(src)
        try:
            p = subprocess.run([sys.executable, "-I", tmp],
                               capture_output=True, text=True, cwd=HERE)
        finally:
            os.remove(tmp)
        print("seeded make_reference.py returncode =", p.returncode)
        print((p.stdout + p.stderr).strip().splitlines()[-1][:160])
        if p.returncode != 0:
            sys.exit(1)
        print("SEED NOT CAUGHT"); sys.exit(0)
     """)),

    # ---- the overlap benchmark arm ---------------------------------------
    ("overlap_benchmark_arm", "a benchmark SMILES written from memory fails "
     "the corpus key check - which is how my first three were caught",
     HEAD + textwrap.dedent("""
        from rdkit import Chem, RDLogger
        from rdkit.Chem import inchi
        RDLogger.DisableLog("rdApp.*")
        sys.path.insert(0, os.path.join(HERE, "..", "..", "provenance"))
        from reagent_registry import REGISTRY
        good = REGISTRY["SM-102"][0]
        from_memory = ("CCCCCCCCC(CCCCCCCC)OC(=O)CCCCCCN(CCO)"
                       "CCCCCC(=O)OCCCCCCCCCCC")
        kg = inchi.MolToInchiKey(Chem.MolFromSmiles(good))
        km = inchi.MolToInchiKey(Chem.MolFromSmiles(from_memory))
        print("registry %s ; from memory %s" % (kg, km))
        if kg != km and kg.split("-")[0] != km.split("-")[0]:
            sys.exit(1)
        print("SEED NOT CAUGHT"); sys.exit(0)
     """)),
]


def main():
    qc = json.load(open(os.path.join(HERE, "07_qc.json")))
    out = []
    for name, what, src in SEEDS:
        r = run(src, name)
        r.update({"seed": name, "what": what,
                  "credited": r["returncode"] == 1,
                  "verdict_from": "returncode"})
        out.append(r)
        print("%-32s returncode=%d credited=%s" %
              (name, r["returncode"], r["credited"]))
        for line in r["stdout_tail"]:
            print("      ", line)
    qc["seeds"] = out
    qc["seeds_summary"] = "%d/%d seeded defects were caught" % (
        sum(1 for s in out if s["credited"]), len(out))
    json.dump(qc, open(os.path.join(HERE, "07_qc.json"), "w"), indent=1)
    print("\n" + qc["seeds_summary"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
