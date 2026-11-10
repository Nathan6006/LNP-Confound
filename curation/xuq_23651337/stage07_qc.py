#!/usr/bin/env python3
"""Stage 07 -- every check, with numbers, and every check SEEDED.

THE RULE THIS FILE ENFORCES. A check is credited only when a seeded defect
makes its AGREEMENT COUNT MOVE. A seed that merely makes the run exit non-zero
is not enough: a no-op seed can trip an unrelated check and look like proof.
So each seed is compared against the clean run count-by-count, and a seed that
moves nothing is reported as NO-OP and fails this stage.

Every seed runs in ITS OWN PROCESS and the verdict is read from
`subprocess.run(...).returncode`. Nothing is inferred from stdout text, and
nothing is inferred from a pipeline's exit status -- `python3 ... | tail` gives
tail's exit code, which is how a "FAILED" report can sit next to EXIT=0.
That is F-012's shape and it was re-encountered while building this folder.

Run:  python3 -I stage07_qc.py
"""
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
PY = [sys.executable, "-I"]

# script -> (clean args, {seed_name: seed args})
SEEDED = {
    "stage05_values.py": {
        "clean": ["--json"],
        "seeds": {s: ["--json", "--seed", s] for s in [
            "seed_major_tick", "seed_minor_tick", "seed_bartop_fill",
            "seed_errbar_origin", "seed_baseline", "seed_xtick",
            "seed_label", "seed_frame_top", "seed_grid",
            "seed_major_tick_gross"]},
        "counts": "json",
    },
    "stage03_enumerate.py": {
        "clean": [],
        "seeds": {s: ["--seed", s] for s in [
            "seed_chain_direction", "seed_alkene_stereo", "seed_chain_length"]},
        "counts": "nagree",
    },
    "stage02_corpus_overlap.py": {
        "clean": [],
        "seeds": {s: ["--seed", s] for s in [
            "seed_tail_structure", "seed_tokens_are_catalogue"]},
        "counts": "nagree",
    },
    "make_reference.py": {
        "clean": [],
        "seeds": {"reverse": ["--seed", "reverse"]},
        "counts": "nagree",
    },
}

NAGREE = re.compile(r"n_agree=(\d+)\s+n_total=(\d+)")


def run(script, args, cwd=HERE):
    p = subprocess.run(PY + [os.path.join(cwd, script)] + args,
                       cwd=cwd, capture_output=True, text=True)
    return p


def counts_of(p, mode):
    if mode == "json":
        for line in p.stdout.splitlines():
            line = line.strip()
            if line.startswith("{"):
                try:
                    d = json.loads(line)
                except json.JSONDecodeError:
                    continue
                return {k: v.get("n_agree") for k, v in d.items()}
        return None
    m = None
    for line in p.stdout.splitlines():
        mm = NAGREE.search(line)
        if mm:
            m = mm
    return None if m is None else {"n_agree": int(m.group(1))}


def ast_guard_selftest():
    """Prove assert_no_assembler() can FAIL. My first instinct was to trust a
    clean line from it; C16's first two versions reported 'all derived' while
    examining nothing, so the alarm is made to fire on purpose."""
    tmp = tempfile.mkdtemp(prefix="anchor_ast_")
    try:
        src = os.path.join(HERE, "make_reference.py")
        dst = os.path.join(tmp, "seeded_make_reference.py")
        text = open(src).read()
        # seed BOTH halves of the provenance claim
        text = text.replace("import ast\n", "import ast\nimport assemble\n", 1)
        text = text.replace("    prov = assert_no_assembler()",
                            "    _ = Chem.RWMol()\n"
                            "    prov = assert_no_assembler(__file__)", 1)
        open(dst, "w").write(text)
        p = subprocess.run(PY + [dst], cwd=tmp, capture_output=True, text=True)
        fired = p.returncode != 0
        why = (p.stderr or p.stdout).strip().splitlines()
        return {"seeded_file": "copy of make_reference.py with `import "
                               "assemble` and a `Chem.RWMol()` call added",
                "exit_code": p.returncode,
                "fired_as_required": bool(fired),
                "last_line": why[-1] if why else "",
                "note": "run in a temporary directory on a COPY. "
                        "make_reference.py itself is untouched."}
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def main():
    report = {"paper": {"pmid": "23651337", "doi": "10.1021/sb300023h"},
              "scripts": {}, "failures": []}

    for script, spec in SEEDED.items():
        clean = run(script, spec["clean"])
        base = counts_of(clean, spec["counts"])
        entry = {"clean_exit_code": clean.returncode,
                 "clean_counts": base, "seeds": {}}
        if clean.returncode != 0:
            report["failures"].append(f"{script}: CLEAN RUN FAILED "
                                      f"(exit {clean.returncode})")
        if base is None:
            report["failures"].append(f"{script}: no agreement count parsed "
                                      "from the clean run")
        for name, args in spec["seeds"].items():
            p = run(script, args)
            c = counts_of(p, spec["counts"])
            moved = []
            if base is not None and c is not None:
                for k in base:
                    if base.get(k) != c.get(k):
                        moved.append({"check": k, "clean": base.get(k),
                                      "seeded": c.get(k)})
            # A seed is credited by ONE of two routes, never by both being
            # blurred together:
            #   (1) an agreement count MOVED -- the check under test scored
            #       the defect;
            #   (2) curation/lib REFUSED, raising a named error before any
            #       count could be computed. That is the library behaving as
            #       designed (it refuses rather than returning a plausible
            #       number), and it is strictly stronger evidence than a count
            #       moving -- but it must be a NAMED refusal, not any crash.
            refusal = None
            for ln in (p.stderr or "").splitlines():
                m2 = re.match(r"^errors\.(\w+Error): (.*)$", ln.strip())
                if m2:
                    refusal = {"exception": m2.group(1),
                               "message": m2.group(2)[:200]}
            ok = p.returncode != 0 and (bool(moved) or refusal is not None)
            entry["seeds"][name] = {
                "exit_code": p.returncode,
                "exits_non_zero": p.returncode != 0,
                "agreement_counts_moved": moved,
                "credited_by": ("agreement_count_moved" if moved
                                else ("library_refusal" if refusal else None)),
                "library_refusal": refusal,
                "credited": bool(ok),
            }
            if p.returncode == 0:
                report["failures"].append(
                    f"{script} --seed {name}: exited 0; the defect was not "
                    "detected")
            elif not moved and refusal is None:
                report["failures"].append(
                    f"{script} --seed {name}: NO-OP SEED -- it exits non-zero "
                    "but no agreement count moved, so it credits nothing. "
                    "Replace it with a seed that lands inside the check it is "
                    "meant to test.")
        report["scripts"][script] = entry

    # the two stages with no seeds of their own, run for their exit code
    for script in ("stage01_blocks.py", "stage06_emit.py"):
        p = run(script, [])
        report["scripts"][script] = {"clean_exit_code": p.returncode,
                                     "seeds": {}}
        if p.returncode != 0:
            report["failures"].append(f"{script}: exit {p.returncode}")

    report["ast_independence_guard_selftest"] = ast_guard_selftest()
    if not report["ast_independence_guard_selftest"]["fired_as_required"]:
        report["failures"].append(
            "assert_no_assembler() did NOT fire on a seeded file that imports "
            "assemble and calls RWMol -- the independence claim is unverified")

    # ---- the headline numbers, measured not restated --------------------
    with open(os.path.join(HERE, "05_values_qc.json")) as fh:
        v = json.load(fh)
    with open(os.path.join(HERE, "02_corpus_overlap.json")) as fh:
        o = json.load(fh)
    with open(os.path.join(HERE, "01_blocks_qc.json")) as fh:
        b = json.load(fh)
    with open(os.path.join(HERE, "06_emit_reconciliation.json")) as fh:
        e = json.load(fh)

    n_seeds = sum(len(s["seeds"]) for s in report["scripts"].values())
    n_credited = sum(1 for s in report["scripts"].values()
                     for d in s["seeds"].values() if d.get("credited"))
    report["summary"] = {
        "endpoint_of_the_paper": "EGFP expression in HeLa after transfection "
            "with EGFP plasmid DNA and EGFP mRNA (plus fibroblasts and "
            "cancer lines) -- a delivery endpoint, ADMISSIBLE IN KIND",
        "where_those_values_are": "the MAIN ARTICLE, which is not in hand",
        "tier_of_what_IS_in_hand": "4b_pdf_vector_geometry",
        "native_dpi": "n/a -- vector; axis.native_tol() deliberately not "
                      "called, no tolerance widened",
        "distinct_structures_determinable": 0,
        "distinct_structures_existing":
            b["library_size_arithmetic"]["distinct_structures_that_exist"],
        "formulations": b["library_size_arithmetic"]["formulations"],
        "individually_measured_values_located":
            v["accounting"]["per_item_numeric_values_measured"],
        "corpus_admissible_values": v["accounting"]["corpus_admissible_yes"],
        "rows_emitted": e["rows_emitted"],
        "free_external_check":
            o["d_free_external_check"]["d1_tail_structure_attestation"][
                "liberate_our_octadecyl_acrylamide_exactly"],
        "catalogue_tokens_rejected":
            o["d_free_external_check"]["d2_catalogue_number_DISPROOF"][
                "tokens_that_are_NOT_catalogue_numbers"],
        "panel_checks": {k: c.get("agreement") for k, c in
                         v["checks"].items()},
        "seeds_run": n_seeds, "seeds_credited": n_credited,
    }

    with open(os.path.join(HERE, "07_qc.json"), "w") as fh:
        json.dump(report, fh, indent=2)

    for script, s in report["scripts"].items():
        print(f"{script}: clean exit {s['clean_exit_code']}, "
              f"{len(s['seeds'])} seeds")
        for name, d in s["seeds"].items():
            mv = ", ".join(f"{m['check']} {m['clean']}->{m['seeded']}"
                           for m in d["agreement_counts_moved"]) or "NOTHING"
            print(f"    {name}: exit {d['exit_code']}, moved: {mv}"
                  f"  {'[credited]' if d['credited'] else '[NOT CREDITED]'}")
    g = report["ast_independence_guard_selftest"]
    print(f"AST guard self-test: exit {g['exit_code']} "
          f"{'[fires]' if g['fired_as_required'] else '[DOES NOT FIRE]'}")
    print(f"seeds credited: {n_credited}/{n_seeds}")
    if report["failures"]:
        print("FAILURES:", *report["failures"], sep="\n  - ")
        return 1
    print("07_qc.json written; all checks pass and every seed is credited")
    return 0


if __name__ == "__main__":
    sys.exit(main())
