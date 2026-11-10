#!/usr/bin/env python3
"""Seed each check's own defect and require that THE CHECK'S METRIC MOVES.

Not 'an exception was raised'.  A harness that only proves something can throw
proves nothing about what the check measures -- a sibling curation in this run
found its seed harness was doing exactly that.  Every case below records the
metric's value on the clean input, the value on the mutated input, and asserts
the specific change.  A case that does not move is reported as NO TEETH.

Mutated copies are written to a temporary directory; nothing in the curation
folder is modified.
"""
import csv
import json
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import stage03_checks as C                                  # noqa: E402


def copy_csv(src, dst, mutate):
    with open(src, newline="") as fh:
        rdr = csv.DictReader(fh)
        rows = list(rdr)
        cols = rdr.fieldnames
    rows = mutate(rows)
    with open(dst, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, restval="")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def get(d, path):
    for p in path.split("."):
        d = d[p]
    return d


def main():
    des = os.path.join(HERE, "01_design_tables.csv")
    rnk = os.path.join(HERE, "02_rank_tables.csv")
    clean = C.run(des, rnk)
    tmp = tempfile.mkdtemp(prefix="dahlman_seed_")
    cases, failures = [], []

    def case(name, metric, mutate_des=None, mutate_rnk=None,
             expect=None, note=""):
        d2 = os.path.join(tmp, "d.csv")
        r2 = os.path.join(tmp, "r.csv")
        copy_csv(des, d2, mutate_des or (lambda rs: rs))
        copy_csv(rnk, r2, mutate_rnk or (lambda rs: rs))
        seeded = C.run(d2, r2)
        before, after = get(clean, metric), get(seeded, metric)
        ok = (after == expect) and (after != before)
        cases.append({"case": name, "metric": metric, "clean": before,
                      "seeded": after, "expected_when_seeded": expect,
                      "metric_moved": after != before, "PASS": ok,
                      "note": note})
        if not ok:
            failures.append(name)

    def drop_one_lib1(rs):
        out, dropped = [], False
        for r in rs:
            if not dropped and r["library"] == "1":
                dropped = True
                continue
            out.append(r)
        return out

    def break_a_mole_sum(rs):
        for r in rs:
            if r["library"] == "1":
                r["Cholesterol Mole %"] = str(float(r["Cholesterol Mole %"]) + 1)
                break
        return rs

    def break_a_rank_composition(rs):
        rs[0]["PEG MW"] = "9999"
        return rs

    def duplicate_a_rank(rs):
        seen = None
        for r in rs:
            if r.get("iMAEC"):
                if seen is None:
                    seen = r["iMAEC"]
                else:
                    r["iMAEC"] = seen
                    break
        return rs

    def push_invivo_rank_over_10(rs):
        for r in rs:
            if r.get("Heart Ecs"):
                r["Heart Ecs"] = "11"
                break
        return rs

    def repair_the_real_typo(rs):
        """The opposite direction: the one live offender is Fig S4f's
        'Lipomer Mole % = 36' on a particle whose other three percents are
        60/5/0.  Setting it to the 35 that Fig S3d and Fig S4f's own
        macrophage sub-table both print must take the offender count to 0."""
        for r in rs:
            if r.get("Lipomer Mole %") == "36" and r.get("Lipomer") == "58":
                r["Lipomer Mole %"] = "35"
        return rs

    case("A: drop one library-1 design row",
         "A_library_row_counts.n_matching",
         mutate_des=drop_one_lib1, expect=3,
         note="4 libraries match their stated size; dropping one row must "
              "take that to 3.")
    case("B: add 1 to one cholesterol mole %",
         "B_mole_percent_sums.n_exactly_100",
         mutate_des=break_a_mole_sum, expect=563,
         note="564 design rows sum to exactly 100; breaking one must give 563.")
    case("E: corrupt one rank row's PEG MW",
         "E_rank_rows_matched_to_a_design_row.matched",
         mutate_rnk=break_a_rank_composition, expect=85,
         note="86 of 87 rank rows match a design row; corrupting a matched "
              "row's composition must give 85.")
    case("F1: repeat a rank inside one in vitro column",
         "F_rank_column_integrity.repeated_ranks_total",
         mutate_rnk=duplicate_a_rank, expect=1,
         note="0 repeated ranks across 24 columns; one deliberate repeat "
              "must give 1.")
    case("F2: push an in vivo rank past the top-10 bound",
         "F_rank_column_integrity.in_vivo_ranks_above_10",
         mutate_rnk=push_invivo_rank_over_10, expect=1,
         note="0 in vivo ranks exceed 10; forcing an 11 must give 1.")
    case("B2 (reverse direction): repair the live 101% row",
         "B2_rank_table_mole_percent_sums.n_exactly_100",
         mutate_rnk=repair_the_real_typo, expect=87,
         note="86 of 87 rank rows sum to 100 BECAUSE of a real defect in the "
              "paper; repairing it must give 87. This proves the metric "
              "tracks the condition and is not stuck.")

    # ---- the structure checks, seeded in stage05's own tables ------------
    import stage05_blocks as B                               # noqa: E402
    saved = dict(B.AMINES)
    B.AMINES = dict(B.AMINES)
    B.AMINES["Spermidine"] = ("NCCCNCCCCNC", "C7H19N3", "amine", "SEEDED")
    out_before = os.path.join(HERE, "05_blocks_qc.json")
    with open(out_before) as fh:
        blocks_clean = json.load(fh)
    cwd = os.getcwd()
    try:
        os.chdir(tmp)
        B.HERE = tmp
        shutil.copy(os.path.join(HERE, "01_design_tables.csv"), tmp)
        rc = B.main()
        with open(os.path.join(tmp, "05_blocks_qc.json")) as fh:
            blocks_seeded = json.load(fh)
    finally:
        os.chdir(cwd)
        B.AMINES = saved
        B.HERE = HERE
    nf_clean = len(blocks_clean["failures"])
    nf_seeded = len(blocks_seeded["failures"])
    ok = nf_clean == 0 and nf_seeded == 1 and rc == 1
    cases.append({"case": "H: lengthen spermidine by one carbon",
                  "metric": "05_blocks_qc.failures (len) / exit code",
                  "clean": nf_clean, "seeded": nf_seeded,
                  "expected_when_seeded": 1, "metric_moved": nf_seeded != nf_clean,
                  "PASS": ok,
                  "note": "the name-derived formula check must catch an extra "
                          "CH2 that leaves the SMILES perfectly parseable; "
                          "exit code must also go 0 -> 1 (got %d)." % rc})
    if not ok:
        failures.append("H")

    # macrocycle ring walk
    pos_clean, size_clean = B.ring_nitrogen_positions("N1CCNCCCNCCCNCCC1")
    pos_bad, size_bad = B.ring_nitrogen_positions("N1CCNCCCNCCNCCCC1")
    # the claim that formula and mass are blind to this is MEASURED, not said
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    from rdkit.Chem import rdMolDescriptors as rdmd
    fm = [(rdmd.CalcMolFormula(Chem.MolFromSmiles(s)),
           round(Descriptors.ExactMolWt(Chem.MolFromSmiles(s)), 4))
          for s in ("N1CCNCCCNCCCNCCC1", "N1CCNCCCNCCNCCCC1")]
    ok = (pos_clean == [1, 4, 8, 12] and pos_bad != [1, 4, 8, 12]
          and size_clean == size_bad == 15 and fm[0] == fm[1])
    cases.append({"case": "I: move one macrocycle nitrogen (1,4,8,12 -> "
                          "1,4,8,11), same ring size, same formula",
                  "metric": "ring_nitrogen_positions",
                  "clean": pos_clean, "seeded": pos_bad,
                  "expected_when_seeded": "!= [1,4,8,12]",
                  "metric_moved": pos_clean != pos_bad, "PASS": ok,
                  "formula_and_mass_both_sides": fm,
                  "note": "both sides measured at %s / %s Da, so formula and "
                          "mass are blind to this; only the ring walk sees "
                          "it." % (fm[0][0], fm[0][1])})
    if not ok:
        failures.append("I")

    shutil.rmtree(tmp, ignore_errors=True)
    rep = {"n_cases": len(cases), "n_passed": len(cases) - len(failures),
           "n_without_teeth": len(failures), "failed_cases": failures,
           "cases": cases}
    with open(os.path.join(HERE, "04_seed_results.json"), "w") as fh:
        json.dump(rep, fh, indent=1)
    for c in cases:
        print("%-5s %-58s clean=%-6s seeded=%-6s" %
              ("PASS" if c["PASS"] else "NOTEETH", c["case"][:58],
               c["clean"], c["seeded"]))
    print("\n%d/%d cases have teeth" % (rep["n_passed"], rep["n_cases"]))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
