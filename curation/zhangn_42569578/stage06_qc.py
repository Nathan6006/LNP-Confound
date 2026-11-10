#!/usr/bin/env python3
"""Merge the stage outputs into 07_qc.json. Every check carries NUMBERS.

07_qc.json must be REGENERABLE, not hand-assembled: a QC file typed by hand is
a memo about checks rather than a record of them, and it drifts from the
artifacts the moment either changes.

Inputs : 07_qc_values.json (stage01), 02_corpus_overlap.json (stage03),
         05_values.csv, and a live run of test_seeds.py.
Run    : python3 -I stage06_qc.py [--outdir DIR] [--skip-seeds]
"""
import argparse
import csv
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=HERE)
    ap.add_argument("--skip-seeds", action="store_true")
    a = ap.parse_args()
    out = os.path.abspath(a.outdir)
    os.makedirs(out, exist_ok=True)

    v = json.load(open(os.path.join(HERE, "07_qc_values.json")))
    ov = json.load(open(os.path.join(HERE, "02_corpus_overlap.json")))
    rows = list(csv.DictReader(open(os.path.join(HERE, "05_values.csv"))))

    if a.skip_seeds:
        rc, nseed, ntot = None, 6, 14
    else:
        p = subprocess.run([sys.executable, "-I",
                            os.path.join(HERE, "test_seeds.py")],
                           capture_output=True, text=True, cwd=HERE)
        rc = p.returncode
        nseed = sum(1 for l in p.stdout.splitlines() if l.startswith("SEED "))
        ntot = sum(1 for l in p.stdout.splitlines()
                   if " PASS" in l or " FAIL" in l or "UNCREDITED" in l)

    qc = json.load(open(os.path.join(HERE, "_qc_template.json")))
    qc["CHECKS"].update({
        "xlsx_rows_with_data": v["xlsx_rows_with_data"],
        "xlsx_header_columns_verified": v["xlsx_header_columns_verified"],
        "xlsx_header_cells_resolved_to_strings":
            v["xlsx_header_cells_resolved_to_strings"],
        "si_tables_found": v["si_tables_found"],
        "si_table_formulation_counts": {"S1": v["si_S1_formulations"],
                                        "S2": v["si_S2_formulations"],
                                        "S3": v["si_S3_formulations"]},
        "design_level_violations_xlsx": v["design_level_violations_xlsx"],
        "design_level_violations_tableS1": v["design_level_violations_tableS1"],
        "design_space_combinations_recomputed": v["design_space_combinations"],
        "design_space_agrees": v["design_space_combinations"] == 864,
        "cholesterol_pct_range": [v["cholesterol_pct_min"],
                                  v["cholesterol_pct_max"]],
        "cholesterol_pct_nonpositive": v["cholesterol_pct_nonpositive"],
        "S2_to_xlsx_composition_bijection": v["S2_to_xlsx_bijection"],
        "S2_xlsx_same_index": v["S2_xlsx_rows_at_same_index"],
        "S2_xlsx_permuted": v["S2_xlsx_rows_permuted"],
        "values_total": len(rows),
        "values_admissible_yes": v["values_admissible_yes"],
        "values_admissible_no": v["values_admissible_no"],
        "values_counts_sum": (v["values_admissible_yes"] +
                              v["values_admissible_no"] == len(rows)),
        "rows_emitted": v["rows_emitted"],
        "rows_header_columns": v["rows_header_columns"],
    })
    qc["CHEMISTRY"]["distinct_ionizable_lipids"] = v["distinct_ionizable_lipids"]
    qc["CHEMISTRY"]["novel_ionizable_lipids"] = v["novel_ionizable_lipids"]
    qc["CORPUS_OVERLAP"]["positive_controls"] = [
        {"control": c["control"], "asserts": c["asserts_molecule"],
         "rows": c["corpus_rows"],
         "matched_the_named_molecule": c["matched_the_named_molecule"],
         "passes": c["passes"]} for c in ov["positive_controls"]]
    qc["CORPUS_OVERLAP"]["lipids"] = {
        l["paper_name"]: {
            "inchikey": l["inchikey"],
            "corpus_rows": l["corpus_rows_by_full_inchikey"],
            "corpus_names": l["corpus_IL_names"],
            "n_experiments": l["n_corpus_experiments"],
            "vivo_experiments": l["vivo_experiments"],
            "rows_a_single_name_join_would_miss":
                l["rows_a_single_name_join_would_MISS"],
            "structure_independent_of_corpus":
                l["structure_is_independent_of_corpus"]}
        for l in ov["paper_lipids"]}
    qc["CORPUS_OVERLAP"]["IL_names_with_more_than_one_structure"] = \
        ov["IL_names_with_more_than_one_structure"]
    qc["CORPUS_OVERLAP"]["IL_structures_with_more_than_one_name"] = \
        ov["IL_structures_with_more_than_one_name"]
    qc["SEED_HARNESS"].update({"ran": not a.skip_seeds, "exit_code": rc,
                               "all_pass": (rc == 0) if rc is not None else None,
                               "seeds": nseed,
                               "controls_and_invariants": ntot - nseed})
    qc["WOULD_BE_Z"] = v["would_be_z"]
    qc["D1_ILC_discriminator"] = v["ILC_discriminator"]
    qc["D1_S2_to_xlsx_map"] = v["S2_to_xlsx_map"]

    with open(os.path.join(out, "07_qc.json"), "w", encoding="utf8") as f:
        json.dump(qc, f, indent=2)
    print("wrote", os.path.join(out, "07_qc.json"))
    print("seed harness exit:", rc)
    for k in ("design_level_violations_xlsx", "S2_to_xlsx_composition_bijection",
              "values_admissible_yes", "values_admissible_no", "rows_emitted"):
        print("  %s = %s" % (k, qc["CHECKS"][k]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
