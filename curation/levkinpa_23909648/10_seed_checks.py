#!/usr/bin/env python3
"""Seeded-mutation suite. A GATE IS NOT DONE UNTIL IT HAS BEEN MADE TO FAIL.

CLAUDE.md section 6: "A NEW CHECK IS NOT DONE UNTIL IT HAS BEEN MADE TO FAIL
ON PURPOSE." Three incidents paid for that rule -- a re-copied C10 anchor that
could only pass, three ZT_2026 marginal checks that scored a literal 0 for
twelve days while the summary beside them said "ALL FOUR EXACT", and two
versions of C16 that reported "all derived" while examining nothing. Every one
printed a clean line.

TWO REQUIREMENTS PER SEED, not one:
  (1) the gate command EXITS NON-ZERO, checked on subprocess.returncode
      DIRECTLY -- never through a pipe, where the exit status is the pipe's;
  (2) the named gate's AGREEMENT COUNT MOVES. A non-zero exit alone is not
      enough: a sibling's +1-unit seed was a no-op under a loose tolerance and
      had to be replaced with x1.01, and the run still "failed" for an
      unrelated reason. So every numeric seed here is MULTIPLICATIVE and the
      count is compared before and after.

Each seed runs in a throwaway copy of the folder, so nothing in the real
folder is touched. Run:  python3 -I 10_seed_checks.py
"""
import json
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PY = [sys.executable, "-I"]


def baseline():
    with open(os.path.join(HERE, "07_qc.json")) as fh:
        return json.load(fh)["gates"]


BASE = baseline()


def run(cwd, script):
    """Run a stage. Returns (returncode, stdout+stderr). Exit code is read
    from returncode, not inferred from output."""
    p = subprocess.run(PY + [script], cwd=cwd, capture_output=True, text=True)
    return p.returncode, p.stdout + p.stderr


def sandbox():
    d = tempfile.mkdtemp(prefix="ll2013_seed_")
    work = os.path.join(d, "levkinpa_23909648")
    shutil.copytree(HERE, work, symlinks=True,
                    ignore=shutil.ignore_patterns("__pycache__"))
    # curation/lib must be importable at ../lib relative to the copy
    os.makedirs(os.path.join(d, "lib"), exist_ok=True)
    for f in os.listdir(os.path.join(ROOT, "curation", "lib")):
        if f.endswith(".py"):
            shutil.copy(os.path.join(ROOT, "curation", "lib", f),
                        os.path.join(d, "lib", f))
    return d, work


def patch_text(path, old, new):
    s = open(path).read()
    if old not in s:
        raise RuntimeError(f"seed target not found in {path}: {old[:70]!r}")
    open(path, "w").write(s.replace(old, new, 1))


def patch_json(path, fn):
    d = json.load(open(path))
    fn(d)
    json.dump(d, open(path, "w"), indent=1)


RESULTS = []


def seed(name, gate_name, mutate, stages=("stage07_qc.py",), expect_raise=False):
    """Apply `mutate` in a sandbox, run `stages`, require failure + movement."""
    d, work = sandbox()
    try:
        mutate(work)
        rc, out = 0, ""
        for s in stages:
            rc, out = run(work, s)
            if rc != 0:
                break
        moved, before, after = None, None, None
        if not expect_raise and os.path.exists(
                os.path.join(work, "07_qc.json")):
            try:
                g = json.load(open(os.path.join(work, "07_qc.json")))["gates"]
                if gate_name in g and gate_name in BASE:
                    before = BASE[gate_name]["agreement_count"]
                    after = g[gate_name]["agreement_count"]
                    moved = before != after
            except Exception:
                moved = None
        ok_exit = rc != 0
        ok_move = True if expect_raise else bool(moved)
        ok = ok_exit and ok_move
        RESULTS.append(dict(seed=name, gate=gate_name, exit_code=rc,
                            exit_nonzero=ok_exit, count_before=before,
                            count_after=after, count_moved=moved,
                            credited=ok,
                            mode=("refusal (stage raises; no count to move)"
                                  if expect_raise else "count movement"),
                            evidence=out.strip().splitlines()[-1][:200]
                            if out.strip() else ""))
        flag = "CREDITED" if ok else "NOT CREDITED"
        extra = ("" if expect_raise
                 else f", count {before} -> {after}")
        print(f"  [{flag}] {name:<26} gate {gate_name:<34} "
              f"exit={rc}{extra}")
        if not ok:
            print(f"             last line: "
                  f"{out.strip().splitlines()[-1][:160] if out.strip() else '(no output)'}")
    finally:
        shutil.rmtree(d, ignore_errors=True)


def main():
    print("Seeded-mutation suite for the LL_2013 gate.\n"
          "Each seed must (1) exit non-zero and (2) move its gate's "
          "agreement count.\n")

    # --- control: an UNMUTATED copy must still pass -----------------------
    d, work = sandbox()
    try:
        rc, out = run(work, "stage07_qc.py")
        print(f"  [CONTROL ] unmutated copy exit={rc} "
              f"({'OK' if rc == 0 else 'BROKEN -- the suite proves nothing'})")
        control_ok = rc == 0
    finally:
        shutil.rmtree(d, ignore_errors=True)

    # =====================================================================
    # Script-level seeds: the parser itself must refuse.
    # =====================================================================
    seed("S1_S4_declared_columns", "Q-1a_table_geometry",
         lambda w: patch_text(
             os.path.join(w, "stage05_values.py"),
             '"eff_1:2", "sd_1:2", "eff_ctrl", "sd_ctrl"], n_rows=13)',
             '"eff_1:2", "sd_1:2"], n_rows=13)'),
         stages=("stage05_values.py",), expect_raise=True)

    seed("S2_rule_y_detector_miss", "Q-1a_table_geometry",
         lambda w: patch_text(os.path.join(w, "stage05_values.py"),
                              "top_rule=330.24", "top_rule=350.24"),
         stages=("stage05_values.py",), expect_raise=True)

    seed("S3_S5_block_size", "Q-1b_S5_dose_blocks",
         lambda w: patch_text(os.path.join(w, "stage05_values.py"),
                              'if len(blk) != 10:', 'if len(blk) != 9:'),
         stages=("stage05_values.py",), expect_raise=True)

    seed("S4_S5_merged_cell_row", "Q-1b_S5_dose_blocks",
         lambda w: patch_text(os.path.join(w, "stage05_values.py"),
                              "if idx != 5:", "if idx != 3:"),
         stages=("stage05_values.py",), expect_raise=True)

    seed("S5_amine1_head_x_CH2", "Q-5e_tableS1_amine1_neutral_masses",
         lambda w: patch_text(os.path.join(w, "stage01_blocks.py"),
                              'smiles="NCCN1CCCC1"',
                              'smiles="NCCCN1CCCC1"'),
         stages=("stage01_blocks.py", "stage07_qc.py"))

    seed("S6_alkyl_chain_plus_one", "Q-11_two_build_routes_agree",
         lambda w: patch_text(os.path.join(w, "stage03_enumerate.py"),
                              'return "[CH2:2]" + "C" * (n - 1)',
                              'return "[CH2:2]" + "C" * n'),
         stages=("stage03_enumerate.py", "stage07_qc.py"))

    seed("S7_emit_literal_confidence", "Q-10_emission_reconciliation",
         lambda w: patch_text(os.path.join(w, "stage06_emit.py"),
                              "confidence=_confidence,",
                              'confidence="high",'),
         stages=("stage06_emit.py",), expect_raise=True)

    # =====================================================================
    # Artifact seeds. All MULTIPLICATIVE, so no seed can be absorbed by a
    # tolerance -- the x1.01 lesson.
    # =====================================================================
    def scale_hrms(w):
        patch_text(os.path.join(w, "stage07_qc.py"),
                   '("1C12-D", 12, 2, 451.4986, 451.4987, 63.20, 59,',
                   '("1C12-D", 12, 2, 455.99858, 451.4987, 63.20, 59,')
    seed("S8_their_calc_mz_x1.01", "Q-5a_ours_vs_their_calculated_mz",
         scale_hrms)

    seed("S9_their_obs_mz_x1.01", "Q-5c_ours_vs_their_observed_mz",
         lambda w: patch_text(
             os.path.join(w, "stage07_qc.py"),
             '("1C16-S", 16, 1, 339.3734, 339.3735, 44.45, 45,',
             '("1C16-S", 16, 1, 339.3734, 342.7672, 44.45, 45,'))

    seed("S10_unprinted_relation_x1.01", "Q-5d_unprinted_HRMS_relations",
         lambda w: patch_text(os.path.join(w, "stage07_qc.py"),
                              "d1 = 563.62378 - 339.3734",
                              "d1 = 563.62378 - 342.76713"))

    seed("S11_figure_integral_x1.10", "Q-6a_1H_figure_integral_sum",
         lambda w: patch_text(
             os.path.join(w, "stage07_qc.py"),
             '("1C16-D", 16, 2, 563.62378, 563.62344, 75.76, 78,',
             '("1C16-D", 16, 2, 563.62378, 563.62344, 83.34, 78,'))

    seed("S12_L2K_control_x1.01", "Q-7_unprinted_dose_inference",
         lambda w: patch_text(
             os.path.join(w, "05_values.csv"),
             "S3-L2K,Table S3,18,281.13,L2K,not stated,"
             "not printed (see Q-7),average transfection efficiency "
             "(percentage of transfected cells),percent,14.0,1.8,Figure 2A",
             "S3-L2K,Table S3,18,281.13,L2K,not stated,"
             "not printed (see Q-7),average transfection efficiency "
             "(percentage of transfected cells),percent,14.14,1.8,Figure 2A"))

    seed("S13_library_shape_x", "Q-2_library_shape",
         lambda w: patch_json(
             os.path.join(w, "01_blocks.json"),
             lambda d: d["library_shape"].update(
                 n_alkyl_halides=11, product=220)))

    seed("S14_mixture_evidence", "Q-3_members_are_mixtures",
         lambda w: patch_json(
             os.path.join(w, "01_blocks.json"),
             lambda d: d["multiple_alkylation_evidence"].update(
                 rows_with_more_than_one_tail_count=1)))

    seed("S15_tail_increment_count", "Q-4_unprinted_tail_increment",
         lambda w: patch_json(
             os.path.join(w, "01_blocks.json"),
             lambda d: d["unprinted_increment_check"].update(
                 n_within_tolerance=30)))

    seed("S16_corpus_overlap_hit", "Q-9_corpus_overlap",
         lambda w: patch_json(
             os.path.join(w, "04_corpus_overlap.json"),
             lambda d: d["our_structures"][0].update(
                 n_full_key_hits=1,
                 full_key_hits=[["vitro", "KW_2014", "FAKE"]])))

    seed("S17_emit_count", "Q-10_emission_reconciliation",
         lambda w: patch_json(
             os.path.join(w, "06_emit_summary.json"),
             lambda d: d.update(n_emitted=28)))

    seed("S18_structure_swap", "Q-11_two_build_routes_agree",
         lambda w: patch_text(
             os.path.join(w, "03_enumerated.csv"),
             "CCCCCCCCCCCCCCCCNCCN1CCCC1",
             "CCCCCCCCCCCCCCCCCNCCN1CCCC1"))

    # --- the three gates the first run of this suite did not exercise -----
    # Found by comparing gates_covered against the gate file, not by
    # inspection: 14 of 17 were seeded and the suite reported a clean
    # 14/14. Coverage of gates is the thing to measure, not seeds passed.
    seed("S19_their_observed_x1.01", "Q-5b_their_calc_vs_their_observed",
         lambda w: patch_text(
             os.path.join(w, "stage07_qc.py"),
             '("1C16-D", 16, 2, 563.62378, 563.62344, 75.76, 78,',
             '("1C16-D", 16, 2, 563.62378, 569.25967, 75.76, 78,'))

    seed("S20_text_listing_x1.09", "Q-6b_1H_text_listing_sum",
         lambda w: patch_text(
             os.path.join(w, "stage07_qc.py"),
             '("1C16-S", 16, 1, 339.3734, 339.3735, 44.45, 45,',
             '("1C16-S", 16, 1, 339.3734, 339.3735, 44.45, 49,'))

    seed("S21_cross_panel_row_removed", "Q-8_cross_panel_ratio",
         lambda w: patch_text(
             os.path.join(w, "05_values.csv"),
             "S6-1C16-S,Table S6,19,435.45,1C16-S,",
             "S6-1C16-S,Table S6,19,435.45,1C16-XX,"))

    # --- the drawn-structure legs (formula, neutral exact mass, MW) ------
    seed("S22_drawn_formula", "Q-5f_drawn_molecular_formula",
         lambda w: patch_text(os.path.join(w, "stage07_qc.py"),
                              '"C22H46N2", 338.36610, 338.61404',
                              '"C23H48N2", 338.36610, 338.61404'))

    seed("S23_drawn_exact_mass_x1.01", "Q-5g_drawn_neutral_exact_mass",
         lambda w: patch_text(os.path.join(w, "stage07_qc.py"),
                              '"C30H62N2", 450.49130, 450.82668',
                              '"C30H62N2", 454.99621, 450.82668'))

    seed("S24_drawn_mol_weight_x1.01", "Q-5h_drawn_molecular_weight",
         lambda w: patch_text(os.path.join(w, "stage07_qc.py"),
                              '"C38H78N2", 562.61650, 563.03932',
                              '"C38H78N2", 562.61650, 568.66971'))

    seed("S25_drawn_bond_count", "Q-12_drawn_bond_geometry",
         lambda w: patch_json(
             os.path.join(w, "02_drawn_geometry.json"),
             lambda d: d["structures"][1].update(
                 drawn_filled_bond_objects=26, bond_count_agrees=False,
                 closed_identity_holds=False)))

    # The F-037 replay: make the anchor generator's template a BYTE COPY of
    # the enumeration's canonical SMILES. The emitted anchor is unchanged and
    # still canonicalises correctly, so a check that compared only the
    # emitted anchor would see nothing -- which is exactly what happened to
    # C10's own `anchor_may_be_a_copy`. Q-13 examines the TEMPLATE, so it
    # fires.
    seed("S26_anchor_template_is_a_copy", "Q-13_C10_anchor_match_and_independence",
         lambda w: patch_text(
             os.path.join(w, "make_reference.py"),
             '        return HEAD + "N(" + tail + ")" + tail',
             '        return tail + "N(" + tail + ")CCN1CCCC1"'))

    n_cred = sum(1 for r in RESULTS if r["credited"])
    out = dict(
        _what="Seeded-mutation suite. Each gate in stage07_qc.py was made to "
              "fail on purpose. A seed is CREDITED only if the gate command "
              "exits non-zero (read from subprocess.returncode, never through "
              "a pipe) AND, where a count exists, that count moves.",
        control_unmutated_copy_passes=control_ok,
        n_seeds=len(RESULTS), n_credited=n_cred,
        gates_covered=sorted({r["gate"] for r in RESULTS}),
        n_gates_in_stage07=len(BASE),
        seeds=RESULTS)
    with open(os.path.join(HERE, "10_seeded.json"), "w") as fh:
        json.dump(out, fh, indent=1)

    print(f"\n{n_cred}/{len(RESULTS)} seeds credited; "
          f"{len(out['gates_covered'])} distinct gates exercised of "
          f"{len(BASE)} in the gate file")
    if not control_ok:
        print("CONTROL FAILED -- the suite proves nothing.")
        return 1
    return 0 if n_cred == len(RESULTS) else 1


if __name__ == "__main__":
    sys.exit(main())
