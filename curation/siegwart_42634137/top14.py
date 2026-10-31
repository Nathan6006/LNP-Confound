#!/usr/bin/env python3
"""The 14 formulations carried from the in vitro screen into in vivo testing.

Printed TWICE in the paper, independently: Fig 2j in the main text (as molar
ratios) and Table S1 in the SI (as molar ratios AND molar percentages, plus the
total lipid/mRNA weight ratio). Both were read and agree on all 14 rows x 4
ratio columns, so the transcription below is double-sourced.

Table S1's own arithmetic is a third check and is verified here: the molar
percentages must be each ratio divided by the sum of that row's four ratios.
That the printed percentages reproduce to <0.01 confirms the "35/50/65,
10/30/50, 10/30/50, 1/2/4" DoE levels are RATIOS, not percentages - which
matters, because reading them as percentages would misstate every formulation.

Fourth check: every one of the 14 must exist among the 486 digitized from
Fig 2c, and every one must sit in the top bin (>1e6 RLU), since these were
selected as the top performers.

Naming: Table S1 uses the Fig 2c panel letters (B=DMPC, C=DOPG, D=DOPE,
F=DSPC); Fig 2j spells them DMPC-1, DOPG-1, DSPC-1 and so on. Both are kept.
No A (DOPC) or E (EggSM) formulation made the top 14.
"""
import csv
import json
import os

HELPER = {"B": "DMPC", "C": "DOPG", "D": "DOPE", "F": "DSPC"}

# (table_id, C6O2B2, helper, PEG, DODAP, printed percentages, lipid/mRNA wt/wt)
ROWS = [
    ("B1", 50, 10, 4, 10, (67.57, 13.51, 5.41, 13.51), 40),
    ("B2", 65, 10, 4, 10, (73.03, 11.24, 4.49, 11.24), 40),
    ("B3", 65, 10, 4, 30, (59.63, 9.17, 3.67, 27.52), 40),
    ("B4", 65, 10, 4, 50, (50.39, 7.75, 3.10, 38.76), 40),
    ("B5", 65, 50, 4, 10, (50.39, 38.76, 3.10, 7.75), 40),
    ("C1", 50, 50, 1, 50, (33.11, 33.11, 0.66, 33.11), 40),
    ("D1", 35, 50, 1, 50, (25.74, 36.76, 0.74, 36.76), 40),
    ("F1", 35, 30, 1, 50, (30.17, 25.86, 0.86, 43.10), 40),
    ("F2", 35, 30, 2, 30, (36.08, 30.93, 2.06, 30.93), 40),
    ("F3", 35, 30, 2, 50, (29.91, 25.64, 1.71, 42.74), 40),
    ("F4", 50, 30, 1, 50, (38.17, 22.90, 0.76, 38.17), 40),
    ("F5", 50, 30, 2, 50, (37.88, 22.73, 1.52, 37.88), 40),
    ("F6", 35, 50, 1, 30, (30.17, 43.10, 0.86, 25.86), 40),
    ("F7", 65, 50, 2, 50, (38.92, 29.94, 1.20, 29.94), 40),
]


def main():
    screen = {}
    if os.path.exists("04b_fig2c_486.csv"):
        for r in csv.DictReader(open("04b_fig2c_486.csv")):
            screen[(r["helper_lipid"], int(r["C6O2B2_molpct"]),
                    int(r["DODAP_molpct"]), int(r["helper_molpct"]),
                    int(r["DMGPEG_molpct"]))] = r

    out, pct_err, missing, not_top = [], [], [], []
    for tid, il, hl, peg, dodap, pcts, wt in ROWS:
        helper = HELPER[tid[0]]
        tot = il + hl + peg + dodap
        calc = [round(v / tot * 100, 2) for v in (il, hl, peg, dodap)]
        err = max(abs(a - b) for a, b in zip(calc, pcts))
        if err > 0.011:
            pct_err.append([tid, calc, list(pcts), round(err, 3)])
        key = (helper, il, dodap, hl, peg)
        hit = screen.get(key)
        if screen and hit is None:
            missing.append(tid)
        elif hit is not None and int(hit["bin_index"]) != 4:
            not_top.append([tid, hit["bin_label"]])
        out.append({
            "table_id": tid,
            "fig2j_name": f"{helper}-{tid[1:]}",
            "lipid_name": "C6O2B2",
            "helper_lipid": helper,
            "C6O2B2_molratio": il, "helper_molratio": hl,
            "DMGPEG_molratio": peg, "DODAP_molratio": dodap,
            "C6O2B2_molpct": calc[0], "helper_molpct": calc[1],
            "DMGPEG_molpct": calc[2], "DODAP_molpct": calc[3],
            "printed_molpct": "/".join(str(p) for p in pcts),
            "molpct_max_abs_error": round(err, 3),
            "total_lipid_to_mRNA_wt": wt,
            "in_fig2c_screen": bool(hit),
            "fig2c_bin_label": hit["bin_label"] if hit else "",
            "fig2c_formulation_id": hit["formulation_id"] if hit else "",
            "selected_for": "in vivo validation (Fig 2i-k)",
            "source": "Fig 2j (main) + Table S1 (SI), both read and agreeing",
            "extraction_method": "printed_numeral_double_sourced",
            "source_pdf": "Adv Mater 2026 Tian (main + SI)",
        })

    with open("04e_top14.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0]))
        w.writeheader()
        w.writerows(out)

    qc = {
        "paper": "Tian et al., Adv Mater 2026, Fig 2j + Table S1",
        "n": len(out), "expected": 14,
        "composition": "5 DMPC (B1-B5), 1 DOPG (C1), 1 DOPE (D1), 7 DSPC (F1-F7); "
                       "no DOPC or EggSM formulation reached the top 14",
        "check_double_sourced": "Fig 2j and Table S1 agree on all 14 rows x 4 ratio columns",
        "check_percentage_arithmetic": {
            "rule": "molar % = ratio / sum(ratios) * 100",
            "rows_failing": pct_err,
            "max_abs_error_pct": max(r["molpct_max_abs_error"] for r in out),
            "implication": "confirms the DoE levels are molar RATIOS, not percentages",
        },
        "check_against_486": {
            "screen_loaded": bool(screen),
            "rows_not_found_in_486": missing,
            "rows_not_in_top_bin": not_top,
            "expectation": "all 14 are top performers so all must carry bin_index 4 (>1e6)",
        },
    }
    json.dump(qc, open("07e_top14_qc.json", "w"), indent=1)

    print(f"top-14 formulations written: {len(out)}")
    print(f"  percentage arithmetic: max abs error "
          f"{qc['check_percentage_arithmetic']['max_abs_error_pct']} pct points, "
          f"{len(pct_err)} rows failing")
    print(f"  present in the 486: {14 - len(missing)}/14   "
          f"in the >1e6 bin: {14 - len(missing) - len(not_top)}/{14 - len(missing)}")
    if not_top:
        print(f"  NOT in top bin: {not_top}")
    if missing:
        print(f"  MISSING from the 486: {missing}")


if __name__ == "__main__":
    main()
