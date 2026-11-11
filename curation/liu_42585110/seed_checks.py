#!/usr/bin/env python3
"""Make every check in this folder fail on purpose. -> 09_seeds.json

    python3 -I seed_checks.py

CLAUDE.md section 6: A NEW CHECK IS NOT DONE UNTIL IT HAS BEEN MADE TO FAIL ON
PURPOSE. Three incidents paid for that rule, and every one printed a clean
line: a re-copied C10 anchor that could only pass, three ZT_2026 marginal
checks that scored a literal 0 for twelve days while the summary beside them
said "ALL FOUR EXACT", and two versions of C16 that reported "all derived"
while examining nothing.

TWO RULES THIS HARNESS OBEYS, both learned the expensive way:

  F-069. IT WRITES TO A TEMP DIRECTORY, NEVER THE CURATION FOLDER. One earlier
  harness pinned each mutated copy's HERE to the real folder so the copy could
  read its inputs, and the copies then wrote their OUTPUTS there -- leaving
  reference_structures.csv holding a seeded regioisomer anchor written by a
  run that exited non-zero, so the non-zero exit CONCEALED the write. A
  file-count or exit-code check passes that; only a cross-run hash comparison
  sees it. This harness copies the folder, mutates the copy, runs in the copy,
  and then verifies every tracked artifact in the REAL folder is byte-identical
  to what it was before the suite ran.

  F-012 / section 10. A SEED IS CREDITED ONLY WHEN ITS METRIC MOVES, and the
  verdict is read from the child's RETURNCODE, never from stdout or a
  pipeline. A seed that exits non-zero while moving nothing is not credited;
  nor is one whose metric line is ABSENT, because an exception before the
  metric is reached looks identical to a caught defect. Where a check's whole
  job is to REFUSE, the seed is credited by the refusal and recorded as a
  distinct route, `refusal`, with the message asserted.
"""
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "09_seeds.json")

TRACKED = ["01_blocks.csv", "02_corpus_overlap.json", "03_enumerated.csv",
           "04_fig1a_cells.csv", "04_fig1c_bars.csv", "04_geometry.json",
           "05_values.csv", "06_rows.csv", "07_qc.json",
           "reference_structures.csv"]

SEEDS = [
    # ---- the C10 anchor's independence guards -------------------------
    dict(id="S1_ast_guard", file="make_reference.py", script="make_reference.py",
         old="def _ast_guard(path):",
         # The seeded import sits inside a function that is NEVER CALLED. The
         # guard reads the module's AST, so it must see the import without
         # Python executing it -- my first version put a bare `import
         # assemble` at module level, which raised ModuleNotFoundError before
         # the guard ran, and the harness correctly refused to credit a seed
         # whose expected message never appeared.
         new="def _never_called():\n    import assemble  # SEEDED\n"
             "    return assemble.canonical(\"C\")\n\n\n"
             "def _ast_guard(path):",
         route="refusal",
         expect_in_stderr="NOT INDEPENDENT OF THE ENUMERATOR",
         why="If the anchor generator may import the enumerator's assembly "
             "library, the anchor can be built by the code it checks -- "
             "F-013 and F-037 in one line."),
    dict(id="S2_template_guard", file="make_reference.py",
         script="make_reference.py",
         old='"note":\n            "Built by writing one SMILES string, '
             'ring-end first.',
         new='"seeded_literal": "@@LONG@@",\n        "note":\n            '
             '"Built by writing one SMILES string, ring-end first.',
         route="refusal", literal_from_enum="AT1-c-ATRA",
         expect_in_stderr="is a SUBSTRING of an enumerated SMILES",
         why="F-037's AMINE_HALF was the enumeration's canonical SMILES with "
             "one token swapped, byte-identical over 141 characters. A long "
             "literal lifted out of the enumerator's output is the signal; "
             "byte-identity of the EMITTED anchor is not."),
    # ---- C10 itself, on a defect with IDENTICAL formula and mass -------
    dict(id="S3_anchor_regioisomer", file="make_reference.py",
         script="make_reference.py",
         old='"OCCCCCC"                       # the linker, one ester only\n'
             '            "N(CCCCO)"',
         new='"OCCCCCC"                       # the linker, one ester only\n'
             '            "N(CCC(C)O)"  # SEEDED regioisomer',
         route="metric",
         metric_re=r"exact=(True|False)", metric_kind="count_true",
         clean=2, seeded=1,
         isomass_of="AT4-Bex",
         why="The seeded head arm is 1-methyl-3-hydroxypropyl instead of "
             "4-hydroxybutyl: SAME molecular formula, SAME exact mass to "
             "0.000000 Da, different connectivity. Formula and mass checks "
             "are blind to it; only the full InChIKey sees it. That is what "
             "C10 exists for."),
    # ---- the double-bond geometry assertion ----------------------------
    dict(id="S4_retinoid_geometry", file="stage02_blocks.py",
         script="stage02_blocks.py",
         old='"13CRA": ("OC(=O)/C=C(C)\\\\C=C\\\\C=C(C)\\\\C=C\\\\'
             'C1=C(C)CCCC1(C)C",',
         new='"13CRA": ("OC(=O)/C=C(C)/C=C/C=C(C)/C=C/C1=C(C)CCCC1(C)C",',
         route="refusal",
         expect_in_stderr="double bonds, expected 3E/1Z",
         why="13-cis and all-trans retinoic acid are isobaric to 0.000000 Da "
             "and share an InChIKey skeleton block. Nothing else in this "
             "pipeline can tell them apart, so the E/Z assertion is the only "
             "thing standing between the corpus and three identical lipids "
             "under three different names."),
    dict(id="S5_z_bond_position", file="stage02_blocks.py",
         script="stage02_blocks.py",
         old='Z_BOND_DISTANCE_FROM_ACID = {"13CRA": 1, "9CRA": 5}',
         new='Z_BOND_DISTANCE_FROM_ACID = {"13CRA": 5, "9CRA": 1}',
         route="refusal",
         expect_in_stderr="bonds from the carboxyl carbon, expected",
         why="13CRA and 9CRA BOTH have 3E/1Z, so a count alone passes with "
             "the two isomers swapped. This seed swaps exactly that."),
    # ---- the HRMS gate --------------------------------------------------
    dict(id="S6_hrms_extra_CH2", file="stage02_blocks.py",
         script="stage03_enumerate.py",
         old='"AT1": ("[CH2:2]CCCCC(=O)OC(CCCCCCC)CCCCCCC",',
         new='"AT1": ("[CH2:2]CCCCCC(=O)OC(CCCCCCC)CCCCCCC",',
         route="metric",
         metric_re=r"HRMS agreement: (\d+)/\d+ within",
         metric_kind="int", clean=39, seeded=29,
         why="One extra CH2 in the AT1 tail. The 10 AT1 lipids must stop "
             "matching their printed exact masses; the other 29 must not "
             "move, which is what makes this a test of the gate and not of "
             "the parser."),
    # ---- the heat-map grid ---------------------------------------------
    dict(id="S7_heatmap_grid_shift", file="stage04_values.py",
         script="stage04_values.py",
         old='    xs_c = heatmap.grid(g["map_x0"], g["map_x1"] + 1, ncol)',
         new='    xs_c = heatmap.grid(g["map_x0"] + 24, g["map_x1"] + 25, '
             'ncol)  # SEEDED half-cell shift',
         route="refusal",
         expect_in_stderr="max nearest-swatch distance",
         why="A half-cell horizontal shift puts every sample on a cell "
             "boundary. heatmap.match's distance diagnostic is what catches "
             "an off-by-one grid, and it has caught a wrong header row twice "
             "elsewhere in this project. Credit is by refusal because the "
             "metric line is only printed on success -- which is itself the "
             "right design."),
    # ---- the bar-top flank test ----------------------------------------
    dict(id="S8_bar_flank_tolerance", file="stage04_values.py",
         script="stage04_values.py",
         old="                x0, x1 = det[i][0], det[i][1] + 1",
         new="                x0, x1 = det[i][0] - 6, det[i][1] + 1  # SEEDED",
         route="metric",
         metric_re=r"bars read (\d+)/82", metric_kind="int",
         clean=77, seeded_max=70,
         why="Widening every bar window 6 px to the left makes the left "
             "flank straddle white page. bars.bar_top's flank test is the "
             "part no single threshold can replace (error class 32, F-018 "
             "and F-020), and this seed is what shows it is live."),
]


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for ch in iter(lambda: fh.read(1 << 20), b""):
            h.update(ch)
    return h.hexdigest()


def snapshot():
    return {f: sha(os.path.join(HERE, f)) for f in TRACKED
            if os.path.exists(os.path.join(HERE, f))}


def make_copy(tmp):
    """A folder copy that keeps ../lib resolvable.

    The stage files do `sys.path.insert(0, <folder>/../lib)`, so a bare copy
    into a temp directory cannot import curation/lib at all and the clean
    control dies with ModuleNotFoundError -- which would have read as "the
    baseline is red" rather than "the harness is wrong". The copy therefore
    reproduces the two-level layout: <root>/lib beside <root>/folder.
    """
    root = tempfile.mkdtemp(prefix="run_", dir=tmp)
    lib_src = os.path.abspath(os.path.join(HERE, "..", "lib"))
    shutil.copytree(lib_src, os.path.join(root, "lib"),
                    ignore=shutil.ignore_patterns("__pycache__"))
    dst = os.path.join(root, "folder")
    os.makedirs(dst)
    for f in os.listdir(HERE):
        if f in ("source", "__pycache__", "09_seeds.json"):
            continue
        src = os.path.join(HERE, f)
        if os.path.isdir(src):
            shutil.copytree(src, os.path.join(dst, f))
        else:
            shutil.copy2(src, os.path.join(dst, f))
    return dst


def run(folder, script):
    """One child process per seed. The verdict is the RETURNCODE."""
    return subprocess.run([sys.executable, "-I", script], cwd=folder,
                          capture_output=True, text=True, timeout=900)


def metric(text, pattern, kind):
    ms = re.findall(pattern, text)
    if not ms:
        return None
    if kind == "int":
        return int(ms[0])
    if kind == "count_true":
        return sum(1 for m in ms if m == "True")
    raise SystemExit(f"seed_checks: unknown metric kind {kind!r}")


def main():
    before = snapshot()
    results = []

    with tempfile.TemporaryDirectory(prefix="kj2026_seeds_") as tmp:
        # ---------- the CLEAN CONTROL, first ---------------------------
        clean_folder = make_copy(tmp)
        clean_runs = {}
        for script in ("make_reference.py", "stage02_blocks.py",
                       "stage03_enumerate.py", "stage04_values.py"):
            r = run(clean_folder, script)
            clean_runs[script] = r
            if r.returncode != 0:
                raise SystemExit(
                    f"seed_checks: THE CLEAN CONTROL FAILED on {script} "
                    f"(returncode {r.returncode}). Every verdict below would "
                    "be meaningless -- a seed cannot be credited against a "
                    f"baseline that is already red.\n{r.stderr[-2000:]}")
        results.append({"id": "S0_clean_control", "route": "control",
                        "credited": True,
                        "detail": {k: v.returncode
                                   for k, v in clean_runs.items()},
                        "why": "A seed suite with no clean control credits "
                               "every seed, including the ones that fire for "
                               "the wrong reason."})

        enum_smiles = {}
        with open(os.path.join(HERE, "03_enumerated.csv")) as fh:
            import csv as _csv
            for row in _csv.DictReader(fh):
                enum_smiles[row["IL_name"]] = row["IL_SMILES"]

        for seed in SEEDS:
            folder = make_copy(tmp)
            path = os.path.join(folder, seed["file"])
            text = open(path, encoding="utf-8").read()
            new = seed["new"]
            if "literal_from_enum" in seed:
                new = new.replace("@@LONG@@",
                                  enum_smiles[seed["literal_from_enum"]][:60])
            if seed["old"] not in text:
                raise SystemExit(
                    f"seed_checks: {seed['id']} could not find its anchor "
                    f"text in {seed['file']}. A seed that does not apply is a "
                    "NO-OP that reads as a clean pass -- the exact shape "
                    "section 10 warns about. Fix the seed, never the guard.")
            open(path, "w", encoding="utf-8").write(
                text.replace(seed["old"], new, 1))

            r = run(folder, seed["script"])
            rec = {"id": seed["id"], "route": seed["route"],
                   "file": seed["file"], "script": seed["script"],
                   "returncode": r.returncode, "why": seed["why"]}

            if seed["route"] == "refusal":
                msg = (r.stderr or "") + (r.stdout or "")
                rec["credited"] = (r.returncode != 0
                                   and seed["expect_in_stderr"] in msg)
                rec["expected_message"] = seed["expect_in_stderr"]
                rec["message_seen"] = seed["expect_in_stderr"] in msg
                rec["note"] = ("Credited by REFUSAL, recorded as a distinct "
                               "route: this check's whole job is to refuse, "
                               "and its metric line is printed only on "
                               "success. The message is asserted so the seed "
                               "cannot be credited by an unrelated crash.")
            else:
                clean_txt = clean_runs[seed["script"]].stdout
                seeded_txt = (r.stdout or "")
                mc = metric(clean_txt, seed["metric_re"], seed["metric_kind"])
                ms = metric(seeded_txt, seed["metric_re"],
                            seed["metric_kind"])
                rec["clean_metric"] = mc
                rec["seeded_metric"] = ms
                if ms is None:
                    rec["credited"] = False
                    rec["note"] = (
                        "NOT CREDITED: the metric line is ABSENT from the "
                        "seeded run's stdout, so the run died before "
                        "reaching it. An exception on the way to a check "
                        "looks identical to the check firing, and crediting "
                        "it would be crediting a crash.")
                else:
                    moved = ms != mc
                    ok = moved
                    if "seeded" in seed:
                        ok = ok and ms == seed["seeded"]
                    if "seeded_max" in seed:
                        ok = ok and ms <= seed["seeded_max"]
                    rec["credited"] = bool(ok)
                    rec["metric_moved"] = bool(moved)
                    rec["note"] = ("Credited only because the AGREEMENT COUNT "
                                   f"MOVED, {mc} -> {ms}. A non-zero exit "
                                   "that moves nothing is not credited.")
            if "isomass_of" in seed:
                from rdkit import Chem
                from rdkit.Chem import Descriptors, rdMolDescriptors
                from rdkit import RDLogger
                RDLogger.DisableLog("rdApp.*")
                sys.path.insert(0, folder)
                for mod in ("make_reference",):
                    sys.modules.pop(mod, None)
                import importlib.util
                spec = importlib.util.spec_from_file_location(
                    "seeded_ref", os.path.join(folder, "make_reference.py"))
                mm = importlib.util.module_from_spec(spec)
                spec.loader.exec_module(mm)
                a = [x for x in mm.ANCHORS
                     if x["lipid_name"] == seed["isomass_of"]][0]
                seeded = Chem.MolFromSmiles(a["smiles"])
                real = Chem.MolFromSmiles(enum_smiles[seed["isomass_of"]])
                rec["isomass_check"] = {
                    "formula_seeded": rdMolDescriptors.CalcMolFormula(seeded),
                    "formula_real": rdMolDescriptors.CalcMolFormula(real),
                    "exact_mass_delta_Da": round(
                        Descriptors.ExactMolWt(seeded)
                        - Descriptors.ExactMolWt(real), 9),
                    "inchikey_seeded": Chem.MolToInchiKey(seeded),
                    "inchikey_real": Chem.MolToInchiKey(real),
                    "note": "The seeded defect is INVISIBLE to formula and to "
                            "exact mass and visible only to the InChIKey. "
                            "That is what proves the gate tests connectivity.",
                }
            results.append(rec)

    after = snapshot()
    changed = [f for f in after if before.get(f) != after[f]]
    report = {
        "_what": "Every check in curation/liu_42585110, made to fail on "
                 "purpose. Section 6: a new check is not done until it has.",
        "n_seeds": len(SEEDS),
        "n_credited": sum(1 for r in results
                          if r.get("credited") and r["route"] != "control"),
        "all_credited": all(r.get("credited") for r in results),
        "folder_untouched_by_the_suite": {
            "what": "F-069. Every mutated copy ran in a temp directory. These "
                    "sha256s are of the REAL folder's tracked artifacts, "
                    "before and after the suite.",
            "files_checked": len(before),
            "files_changed": changed,
            "clean": changed == [],
        },
        "seeds": results,
    }
    with open(OUT, "w") as fh:
        json.dump(report, fh, indent=1)

    for r in results:
        mark = "CREDITED" if r.get("credited") else "NOT CREDITED"
        extra = ""
        if r["route"] == "metric":
            extra = f"  metric {r.get('clean_metric')} -> {r.get('seeded_metric')}"
        print(f"  {mark:<12} {r['id']:<26} ({r['route']}){extra}")
    print(f"\n{report['n_credited']}/{len(SEEDS)} seeds credited; "
          f"folder untouched: {report['folder_untouched_by_the_suite']['clean']}")
    if not report["all_credited"]:
        raise SystemExit("seed_checks: a seed was NOT credited. Either the "
                         "check cannot fire, or the seed does not reach it. "
                         "Fix the seed; never weaken the guard.")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
