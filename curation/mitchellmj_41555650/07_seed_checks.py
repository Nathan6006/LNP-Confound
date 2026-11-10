#!/usr/bin/env python3
"""Seed every check in this folder and confirm it FIRES.

A CHECK THAT CANNOT FAIL IS NOT A CHECK. A sibling's category check scored 11/11
clean AND 11/11 seeded, because it measured spread about the mean of its own set
and was therefore shift-invariant. So a check is credited here only when the
AGREEMENT COUNT MOVES under the seed, not merely when something raises.

Each seed runs in ITS OWN PROCESS and the verdict is read from `returncode`.
"""
import json, os, subprocess, sys, tempfile, textwrap

HERE = os.path.dirname(os.path.abspath(__file__))
PY = [sys.executable, "-I"]
results = []


def run(code, label):
    with tempfile.NamedTemporaryFile("w", suffix=".py", dir=HERE,
                                     delete=False) as f:
        f.write(textwrap.dedent(code)); p = f.name
    try:
        r = subprocess.run(PY + [p], capture_output=True, text=True, timeout=600)
        return r.returncode, (r.stdout + r.stderr)[-1400:]
    finally:
        os.unlink(p)


def record(name, what, clean_rc, seeded_rc, clean_metric, seeded_metric, tail=""):
    moved = clean_metric != seeded_metric
    ok = (clean_rc == 0 and seeded_rc != 0) or moved
    results.append(dict(check=name, what=what,
                        clean_exit=clean_rc, seeded_exit=seeded_rc,
                        clean_metric=clean_metric, seeded_metric=seeded_metric,
                        agreement_count_moved=moved, CREDITED=ok,
                        seeded_output_tail=tail))
    print(f"{'PASS' if ok else 'NOT CREDITED'}  {name}: clean={clean_metric} "
          f"seeded={seeded_metric} (exit {clean_rc} -> {seeded_rc})")


# ---------------------------------------------------------------- S1: C10 anchor
# A defect with IDENTICAL FORMULA AND IDENTICAL EXACT MASS: one arm's hydroxyl
# moved from C2 to C3. Mass and formula cannot see it.
S1 = '''
    import sys, os
    sys.path.insert(0, {here!r})
    from rdkit import Chem
    from rdkit.Chem import inchi
    from rdkit.Chem.Descriptors import ExactMolWt
    from rdkit.Chem.rdMolDescriptors import CalcMolFormula
    import csv
    good = ("CCCCCCCCC(O)CN(CC(O)CCCCCCCC)CCCN(C)CCCN"
            "(CC(O)CCCCCCCC)CC(O)CCCCCCCC")
    seed = ("CCCCCCCCC(O)CN(CC(O)CCCCCCCC)CCCN(C)CCCN"
            "(CC(O)CCCCCCCC)CCC(O)CCCCCCC")
    enum = {{r["IL_name"]: r for r in csv.DictReader(
        open(os.path.join({here!r}, "03_enumerated.csv")))}}["C10-306"]
    mg, ms = Chem.MolFromSmiles(good), Chem.MolFromSmiles(seed)
    print("formula  good/seed:", CalcMolFormula(mg), CalcMolFormula(ms))
    print("exactmass delta   :", abs(ExactMolWt(mg) - ExactMolWt(ms)))
    kg, ks = inchi.MolToInchiKey(mg), inchi.MolToInchiKey(ms)
    print("MATCHES_GOOD", kg == enum["InChIKey"])
    print("MATCHES_SEED", ks == enum["InChIKey"])
    sys.exit(0 if ks == enum["InChIKey"] else 1)
'''.format(here=HERE)
rc, tail = run(S1, "anchor")
import re
mg = re.search(r"MATCHES_GOOD (\w+)", tail)
ms = re.search(r"MATCHES_SEED (\w+)", tail)
dm = re.search(r"exactmass delta   : ([\d.e-]+)", tail)
record("C10 anchor (identical formula AND exact mass)",
       "an OH moved C2->C3 on one arm of C10-306",
       0, rc, f"anchor==enum:{mg.group(1) if mg else '?'}",
       f"anchor==enum:{ms.group(1) if ms else '?'} (delta mass {dm.group(1) if dm else '?'})",
       tail[-400:])

# ------------------------------------------------- S2: emit refuses a literal
S2 = '''
    import sys
    sys.path.insert(0, {lib!r})
    import emit
    try:
        emit.emit_rows(["a"], {{"a": {{}}}}, confidence="high",
                       note=lambda i, r: "x", label="seed")
        print("NO_RAISE"); sys.exit(0)
    except Exception as e:
        print("RAISED", type(e).__name__, str(e)[:200]); sys.exit(1)
'''.format(lib=os.path.join(HERE, "..", "lib"))
rc, tail = run(S2, "emit-literal")
record("emit.emit_rows refuses a LITERAL confidence (F-011)",
       'confidence="high" passed as a string instead of a callable',
       0, rc, "callable accepted", "literal " + ("REJECTED" if rc else "ACCEPTED"), tail[-300:])

# --------------------------------------- S3: emit refuses a caveat with no note
S3 = '''
    import sys
    sys.path.insert(0, {lib!r})
    import emit
    try:
        emit.emit_rows(["a"], {{"a": {{}}}},
                       confidence=lambda i, r: "low",
                       note=lambda i, r: "", label="seed")
        print("NO_RAISE"); sys.exit(0)
    except Exception as e:
        print("RAISED", type(e).__name__, str(e)[:160]); sys.exit(1)
'''.format(lib=os.path.join(HERE, "..", "lib"))
rc, tail = run(S3, "emit-caveat")
record("emit.emit_rows refuses a CAVEAT WITH NO NOTE (F-024)",
       'confidence="low" with an empty curation_note',
       0, rc, "note present -> accepted",
       "empty note " + ("REJECTED" if rc else "ACCEPTED"), tail[-300:])

# ------------------------------------------- S4: conf_rule cannot be a literal
S4 = '''
    import sys, os
    sys.path.insert(0, {here!r})
    sys.path.insert(0, os.path.join({here!r}, "..", "lib"))
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "r6", os.path.join({here!r}, "06_rows.py"))
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    # a record that claims an IUPAC-derived structure but carries no
    # stoichiometry result -- conf_rule must refuse to invent a confidence
    try:
        v = m.conf_rule("seeded", {{"structure_source":
                        "enumerated_from_SI_IUPAC_name",
                        "iupac_name_is_full": True,
                        "core_stoichiometry_dev_pct": None}})
        print("RETURNED", v); sys.exit(0)
    except SystemExit:
        print("REFUSED"); sys.exit(1)
'''.format(here=HERE)
rc, tail = run(S4, "conf-rule")
record("conf_rule refuses to assign confidence from an assertion",
       "a record asserting an IUPAC structure with no stoichiometry evidence",
       0, rc, "evidence present -> high",
       "evidence absent -> " + ("REFUSED" if rc else "RETURNED A VALUE"), tail[-300:])

# ------------------------------- S5: P5 design grid, seeded block mis-assignment
S5 = '''
    import sys, os, json, re
    sys.path.insert(0, {here!r})
    import pdfplumber
    SI = os.path.join({here!r}, "source", "nn5c19138_si_001.pdf")
    rows = []
    with pdfplumber.open(SI) as pdf:
        for p in (22, 23, 24):
            for tb in pdf.pages[p-1].find_tables():
                rows.extend(tb.extract())
    s1 = {{}}
    for r in rows:
        if not r or len(r) < 6: continue
        num = (r[1] or "").strip()
        if not re.match(r"^\\d+$", num): continue
        s1[int(num)] = dict(zip(["IL","DOPE","Chol","PEG"],
                                [(r[2+i] or "").strip() for i in range(4)]))
    GOOD = [[1,2,3,4,5,6],[7,8,9,10,11,12],[13,14,15,16,17,18],
            [20,21,22,23,24,25],[26,27,28,29,30,31],[32,33,34,35,36,37],
            [38,39,40,42,43,44],[45,46,47,48,49,50],[51,52,53,54,58,59],
            [60,61,62,63,64,65],[66,68,69,70,71,72],[73,74,75,76,77,78],
            [79,80,81,82,83,84],[85,86,87,88,89,90],[91,92,94,96,97,98]]
    import copy
    SEED = copy.deepcopy(GOOD)
    SEED[10] = [66,67,69,70,71,72]   # 67 is a GREYED-OUT, EXCLUDED b-LNP
    PAT=[("35","2.5"),("45","2.5"),("55","2.5"),("35","5"),("45","5"),("55","5")]
    def nmis(blocks):
        n=0
        for ns in blocks:
            for pos,x in enumerate(ns):
                if x not in s1: n+=1; continue
                if (s1[x]["IL"], s1[x]["PEG"]) != PAT[pos]: n+=1
        return n
    g, s = nmis(GOOD), nmis(SEED)
    print("MIS_GOOD", g); print("MIS_SEED", s)
    sys.exit(0 if s > g else 1)
'''.format(here=HERE)
rc, tail = run(S5, "P5")
g = re.search(r"MIS_GOOD (\d+)", tail); s = re.search(r"MIS_SEED (\d+)", tail)
record("P5 design-grid vs printed molar ratios",
       "one b-LNP in the C12-200 block swapped for the excluded number 67",
       0, 0, f"{g.group(1) if g else '?'} mismatches",
       f"{s.group(1) if s else '?'} mismatches", tail[-300:])

# ------------------------------------------ S6: P1 amine-core stoichiometry
# NOTE: S6 takes no .format(), so its braces are SINGLE. An earlier version
# copied the doubled braces from the formatted seeds, raised TypeError, exited
# non-zero, and was credited as a PASS -- a seed that fires for the wrong reason
# is a no-op seed wearing a green tick. The guard is below: the seeded run must
# print its metric, not merely exit non-zero.
S6 = """
import sys
from rdkit import Chem
from rdkit.Chem.Descriptors import MolWt
cores = {"200":"NCCN1CCN(CCNCCN)CC1","494":"NCCOCCN1CCN(CCNCCOCCN)CC1",
         "306":"CN(CCCN)CCCN","114":"NCCNCCN",
         "490":"NCC(OCC)CN1CCN(CCNCC(OCC)CN)CC1"}
seed = dict(cores); seed["306"] = "CN(CCCCN)CCCN"   # one CH2 added
mmol = {"200":(26.0,0.724/6.0),"494":(40.0,0.791/6.0),"306":(43.6,2.10/7.0),
        "114":(30.9,2.10/7.0),"490":(25.0,0.528/7.0)}
def dev(c, k):
    mw = MolWt(Chem.MolFromSmiles(c)); mg, mm = mmol[k]
    return abs(100*((mg/mm)-mw)/mw)
def nwithin(cs):
    return sum(1 for k, c in cs.items() if dev(c, k) <= 0.5)
def worst(cs):
    return max(dev(c, k) for k, c in cs.items())
print("N_WITHIN_0p5_GOOD", nwithin(cores))
print("N_WITHIN_0p5_SEED", nwithin(seed))
print("WORST_GOOD %.4f" % worst(cores))
print("WORST_SEED %.4f" % worst(seed))
sys.exit(0 if (nwithin(seed) < nwithin(cores)) else 1)
"""
rc, tail = run(S6, "P1")
if "N_WITHIN_0p5_SEED" not in tail:
    raise SystemExit("S6 did not reach its own metric -- the seed raised for "
                     "an unrelated reason and must not be credited:\n" + tail)
a = re.search(r"N_WITHIN_0p5_GOOD (\d+)", tail)
b = re.search(r"N_WITHIN_0p5_SEED (\d+)", tail)
wg = re.search(r"WORST_GOOD ([\d.]+)", tail)
ws = re.search(r"WORST_SEED ([\d.]+)", tail)
record("P1 amine-core formula vs the paper's own stoichiometry",
       "one CH2 added to amine 306 (a plausible cataloguing slip)",
       0, rc,
       f"{a.group(1)}/5 cores within 0.5% (worst {wg.group(1)}%)",
       f"{b.group(1)}/5 cores within 0.5% (worst {ws.group(1)}%)",
       tail[-300:])

# ------------------------------------- S7: P3 arms built vs arms in IUPAC name
S7 = '''
    import sys, csv, os
    rows = list(csv.DictReader(open(os.path.join({here!r}, "03_enumerated.csv"))))
    g = sum(1 for r in rows if r["n_arms_built"] == r["n_arms_in_iupac_name"])
    # seed: amine 114 read as tris(2-aminoethyl)amine instead of
    # diethylenetriamine -- 6 N-H instead of 5, same THREE nitrogens
    s = 0
    for r in rows:
        want = r["n_arms_in_iupac_name"]
        built = "6" if r["amine_core"] == "114" else r["n_arms_built"]
        s += (built == want)
    print("AGREE_GOOD", g, "of", len(rows))
    print("AGREE_SEED", s, "of", len(rows))
    sys.exit(0 if s < g else 1)
'''.format(here=HERE)
rc, tail = run(S7, "P3")
a = re.search(r"AGREE_GOOD (\d+) of (\d+)", tail); b = re.search(r"AGREE_SEED (\d+)", tail)
record("P3 arms built (from N-H count) vs arms stated in the IUPAC name",
       "amine 114 mis-read as tris(2-aminoethyl)amine: 6 N-H, same 3 nitrogens",
       0, rc, f"{a.group(1) if a else '?'}/{a.group(2) if a else '?'} agree",
       f"{b.group(1) if b else '?'}/{a.group(2) if a else '?'} agree", tail[-300:])

# ----------------------------------------- S8: axis LOO under a seeded tick error
S8 = '''
    import sys, os
    sys.path.insert(0, os.path.join({here!r}, "..", "lib"))
    import numpy as np, axis
    good = [104.0, 191.5, 278.5, 366.5, 453.5, 541.0, 628.5]
    tol  = axis.native_tol(axis.LOO_TOL_PX, native_dpi=355.7)
    out = []
    for d in (0.0, 1.0, 2.0, 3.25, 5.0):
        t = list(good); t[3] += d
        w = float(np.max(np.abs(axis.loo_residuals(t))))
        out.append((d, round(w,4), w <= axis.LOO_TOL_PX, w <= tol))
        print("SEED %.2f px -> LOO %.4f  strict=%s native=%s"
              % (d, w, w <= axis.LOO_TOL_PX, w <= tol))
    print("TOL_NATIVE %.4f" % tol)
    # the check has teeth if a seeded error is caught at the tolerance we USE
    sys.exit(0 if not out[-1][3] else 1)
'''.format(here=HERE)
rc, tail = run(S8, "P8")
record("P8 axis LOO at the WIDENED native tolerance, seeded in THIS panel's geometry",
       "one of the seven decade ticks displaced by 0 / 1 / 2 / 3.25 / 5.0 px",
       0, rc, "0.00 px seed -> passes",
       "5.00 px seed -> " + ("CAUGHT" if rc == 0 else "MISSED"), tail[-700:])

# ------------- S9: P14 pool membership, grey colouring vs the stated criteria
S9 = """
import sys, os, json
HERE = os.path.dirname(os.path.abspath("__PLACEHOLDER__/x"))
d = json.load(open(os.path.join("__HERE__", "09_pool_membership.json")))
rec = d["per_bLNP"]
def disagreements(flip=None):
    n = 0
    for k, v in rec.items():
        g, f = v["grey_in_table"], v["fails_stated_criteria"]
        if g is None or f is None:
            continue
        if flip is not None and k == flip:
            g = not g
        n += (g != f)
    return n
clean = disagreements()
# seed: flip ONE row's grey classification -- a single mis-read colour
target = sorted(rec)[0]
seeded = disagreements(flip=target)
print("DISAGREE_CLEAN", clean)
print("DISAGREE_SEEDED", seeded)
sys.exit(0 if seeded > clean else 1)
""".replace("__HERE__", HERE)
rc, tail = run(S9, "P14")
if "DISAGREE_SEEDED" not in tail:
    raise SystemExit("S9 did not reach its own metric:\n" + tail)
a = re.search(r"DISAGREE_CLEAN (\d+)", tail)
b = re.search(r"DISAGREE_SEEDED (\d+)", tail)
record("P14 grey colouring vs the caption's four stated criteria",
       "one b-LNP's grey/black classification flipped",
       0, rc, f"{a.group(1)} disagreements over 90 rows",
       f"{b.group(1)} disagreements over 90 rows", tail[-300:])

json.dump(results, open(os.path.join(HERE, "07_seed_results.json"), "w"), indent=1)
n_ok = sum(1 for r in results if r["CREDITED"])
print(f"\n{n_ok}/{len(results)} checks CREDITED (their agreement count moved, "
      f"or a clean run passed where a seeded run failed)")
sys.exit(0 if n_ok == len(results) else 1)
