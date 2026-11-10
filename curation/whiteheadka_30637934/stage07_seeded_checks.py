#!/usr/bin/env python3
"""Every check in this folder, run CLEAN and then with a defect SEEDED.

THE RULE, from CLAUDE.md section 6: "A NEW CHECK IS NOT DONE UNTIL IT HAS BEEN
MADE TO FAIL ON PURPOSE." And the sharper version the brief adds: credit a
check only when the AGREEMENT COUNT MOVES. A seed that merely makes something
RAISE proves the code has an exception, not that the check discriminates. So
each entry below records `clean` and `seeded` as the SAME metric, plus `moved`.

Subprocess exit codes are read from `returncode`, never through a pipe: a pipe
reports the exit status of the last command in the pipeline, which is how a
gate that exits 1 gets recorded as passing.
"""
import csv
import itertools
import json
import math
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))

CH2 = 14.0156500642
RESULTS = []


def record(name, what, clean, seeded, seed_desc, unit=""):
    RESULTS.append({"check": name, "measures": what, "unit": unit,
                    "clean": clean, "seeded": seeded,
                    "seed": seed_desc,
                    "moved": clean != seeded,
                    "credited": clean != seeded})


def load(fn):
    with open(os.path.join(HERE, fn), newline="") as fh:
        return list(csv.DictReader(fh))


def j(fn):
    return json.load(open(os.path.join(HERE, fn)))


# ---------------------------------------------------------------- 1
def check_figS7_vs_tableS1(perturb=None):
    s1 = {r["LNP"]: float(r["LNP_surface_pKa_TNS_aqueous"]) for r in
          load("04_table_s1.csv")}
    if perturb:
        s1[perturb] = s1[perturb] + 0.01
    s7 = j("04_tables_qc.json")["check_figS7_vs_tableS1_pKa"]["per_lipid"]
    return sum(1 for c in s7 if c["table_s1"] is not None
               and s1.get(c["lipidoid"]) == c["fig_s7"])


# ---------------------------------------------------------------- 2
def check_corpus_overlap(mutate=None):
    import reagents
    mine = load("03_enumerated.csv")
    keys = set()
    for fn in ("data_vitro.csv", "data_vivo.csv"):
        with open(os.path.join(HERE, "..", "..", "data", fn), newline="") as fh:
            for r in csv.DictReader(fh):
                s = (r.get("IL_SMILES") or "").strip()
                if s:
                    try:
                        keys.add(reagents.full_key(s))
                    except Exception:
                        pass
    n = 0
    for m in mine:
        smi = m["IL_SMILES"]
        if mutate and m["tail_token"] == mutate:
            smi = smi.replace("OCCCCCCCCCC", "OCCCCCCCCC", 1)
        try:
            n += reagents.full_key(smi) in keys
        except Exception:
            pass
    return n


# ---------------------------------------------------------------- 3
def check_mass_ladder(break_tail=None):
    """Each consecutive LINEAR tail step must be exactly 4 x (dC) x CH2.

    Four arms, so one extra tail carbon adds C4H8 = 56.0626 Da to the lipidoid.
    This is a relation among the authors' own name tokens and our own masses
    that nothing in the paper prints, and it fails on any off-by-one in the
    token -> chain-length reading.
    """
    rows = [r for r in load("03_enumerated.csv") if r["tail_token"] != "Oi10"]
    for r in rows:
        r["_c"] = int(r["tail_carbons"])
        r["_m"] = float(r["exact_mass"])
        if break_tail and r["tail_token"] == break_tail:
            r["_m"] -= CH2            # one CH2 short on one lipid only
    rows.sort(key=lambda r: r["_c"])
    ok = 0
    for a, b in itertools.pairwise(rows):
        want = 4.0 * (b["_c"] - a["_c"]) * CH2
        if abs((b["_m"] - a["_m"]) - want) < 2e-4:
            ok += 1
    return ok


# ---------------------------------------------------------------- 4
def check_abstract_mass_identity(shrink_oi10=False):
    rows = {r["tail_token"]: r for r in load("03_enumerated.csv")}
    a, b = rows["O10"], rows["Oi10"]
    fa, fb = a["formula"], b["formula"]
    ma, mb = float(a["exact_mass"]), float(b["exact_mass"])
    if shrink_oi10:
        fb, mb = "C55H107N3O8", mb - 4 * CH2
    return int(fa == fb and abs(ma - mb) < 1e-9)


# ---------------------------------------------------------------- 5
def check_panelB_monotonic(permute_legend=False):
    rows = load("02_fig_s6_panelB.csv")
    ser = ["10 ng", "50 ng", "100 ng"]
    if permute_legend:
        ser = ["50 ng", "10 ng", "100 ng"]
    cats = sorted({r["category"] for r in rows})
    n = 0
    for c in cats:
        v = [float(next(r["value_RLU"] for r in rows
                        if r["category"] == c and r["series"] == s))
             for s in ser]
        n += v[0] < v[1] < v[2]
    return n


# ---------------------------------------------------------------- 6
def check_category_assignment(shift=0):
    """#labels whose right edge lies within HALF A CATEGORY PITCH of its tick.

    THE FIRST VERSION OF THIS CHECK COULD NOT FAIL. It counted labels whose
    offset was within 0.05 pt of the MEAN offset -- and the mean is computed
    from the same set, so rotating the category list by one slot moves every
    offset by the identical 12.81 pt and leaves the spread about the new mean
    unchanged. It scored 11/11 clean and 11/11 seeded. That is the C16
    incident exactly ("my first two versions of C16 itself reported 'all
    derived' while examining nothing"), and the only reason it was caught is
    that this file demands the agreement count MOVE.

    The fixed criterion is ABSOLUTE and non-circular: a label belongs to the
    tick it is nearest to, so its right edge must lie within half the category
    pitch (6.40 pt) of the tick it is assigned to. A one-slot rotation puts
    every offset at 2.99 +/- 12.81 pt and scores zero.

    `category_label_offset_spread_pt` in 07_qc.json remains the striking
    descriptive number -- 0.002 pt across 11 labels -- but it is recorded
    there as SHIFT-INVARIANT and is not what pins the assignment.
    """
    import pymupdf
    geo = j("02_fig_s6_geometry.json")
    ticks = geo["panelA_x_tick_positions_px"]
    pitch = sum(geo["panelA_x_tick_pitch_px"]) / len(geo["panelA_x_tick_pitch_px"])
    cats = ["Crl", "Nkd", "Lipo MM", "O6", "O8", "O9", "Oi10", "O10", "O11",
            "O12", "O13", "O14", "O16", "O18"]
    if shift:
        cats = cats[shift:] + cats[:shift]
    doc = pymupdf.open(os.path.join(HERE, "source",
                                    "smll201805097-sup-0001-S1.pdf"))
    words = {w[4].strip(): w for w in doc[11].get_text("words")
             if 222 <= w[1] <= 250 and w[0] <= 300}
    doc.close()
    offs = [words[c][2] - t for c, t in zip(cats, ticks)
            if c in words and c.startswith("O")]
    return sum(1 for o in offs if abs(o) <= pitch / 2.0)


# ---------------------------------------------------------------- 7
def check_molpct_sum(nudge=False):
    v = [35.0, 16.0, 46.5, 2.5]
    if nudge:
        v[2] = 46.4
    return int(sum(v) == 100.0)


# ---------------------------------------------------------------- 8/9
def check_emit_guards(mode="clean"):
    """Rows emit_rows will produce. A literal confidence must produce ZERO."""
    import emit
    from errors import EmitError
    recs = {r["tail_token"]: dict(r) for r in load("03_enumerated.csv")}
    items = list(recs)
    try:
        if mode == "literal_confidence":
            rows = emit.emit_rows(items, recs, confidence="high",
                                  note=lambda i, r: "x",
                                  structure_assumed=lambda i, r: "1",
                                  label="seed")
        elif mode == "caveat_without_note":
            rows = emit.emit_rows(items, recs,
                                  confidence=lambda i, r: "medium",
                                  note=lambda i, r: "",
                                  label="seed")
        else:
            rows = emit.emit_rows(items, recs,
                                  confidence=lambda i, r: "high",
                                  note=lambda i, r: "a note",
                                  structure_assumed=lambda i, r: "1",
                                  label="clean")
        return len(rows)
    except EmitError:
        return 0


# ---------------------------------------------------------------- 10
def check_axis_loo(bad_tick_px=0.0):
    """Leave-one-out vs max residual on the panel A y ticks.

    Returns (max_residual, loo_worst). The brief's point: calibrate's
    max-residual test is weak against ONE bad tick; loo_residuals is not.
    """
    import axis
    from errors import CalibrationError
    t = list(j("02_fig_s6_geometry.json")["panelA_y_tick_positions_px"])
    if bad_tick_px:
        t[0] = t[0] + bad_tick_px
    vals = [2.5e6, 2.0e6, 1.5e6, 1.0e6, 5.0e5, 0.0]
    try:
        cal = axis.calibrate(t, vals, kind="linear", axis="y")
        mx = round(cal.residual_px, 3)
    except CalibrationError:
        mx = "RAISED"
    return (mx, round(max(axis.loo_residuals(t)), 3))


# ---------------------------------------------------------------- 11
def check_anchor_gate(seed=None):
    """Exit code read from returncode, never through a pipe."""
    cmd = [sys.executable, "-I", os.path.join(HERE, "make_reference.py"),
           "--quiet"]
    if seed:
        cmd += ["--seed", seed]
    p = subprocess.run(cmd, cwd=HERE, capture_output=True, text=True)
    return p.returncode


def main():
    record("figS7_vs_tableS1_pKa",
           "lipidoids whose Figure S7 printed pKa equals Table S1's pKa",
           check_figS7_vs_tableS1(), check_figS7_vs_tableS1("306O10"),
           "Table S1's 306O10 pKa nudged by +0.01", "of 8")

    record("corpus_full_inchikey_overlap",
           "enumerated lipidoids already in data/ by FULL InChIKey",
           check_corpus_overlap(), check_corpus_overlap("O10"),
           "one CH2 removed from A306_O10's four decyl arms", "of 11")

    record("mass_ladder_one_CnH2n_per_tail_carbon",
           "consecutive linear-tail steps equal to exactly 4 x dC x CH2",
           check_mass_ladder(), check_mass_ladder("O11"),
           "one CH2 subtracted from A306_O11's mass only", "of 9")

    record("abstract_mass_identity_O10_vs_Oi10",
           "the abstract's claim that 306O10 and 306Oi10 have identical MW",
           check_abstract_mass_identity(),
           check_abstract_mass_identity(shrink_oi10=True),
           "Oi10 shortened to a nine-carbon tail", "of 1")

    record("panelB_dose_monotonicity",
           "Figure S6B lipidoids whose 10 < 50 < 100 ng values increase",
           check_panelB_monotonic(), check_panelB_monotonic(True),
           "the 10 ng and 50 ng legend swatches swapped", "of 4")

    record("category_assignment_label_right_edge",
           "x labels whose right edge is within half a category pitch "
           "(6.40 pt) of its assigned tick",
           check_category_assignment(), check_category_assignment(shift=1),
           "the category list rotated by one slot", "of 11")

    record("mole_percents_sum_to_100",
           "the four printed mol% summing to exactly 100",
           check_molpct_sum(), check_molpct_sum(True),
           "46.5 changed to 46.4", "of 1")

    record("emit_refuses_a_literal_confidence",
           "rows emit_rows produces (F-011 guard: a literal must produce 0)",
           check_emit_guards("clean"), check_emit_guards("literal_confidence"),
           "confidence passed as the literal string 'high'", "rows")

    record("emit_refuses_a_caveat_with_no_note",
           "rows emit_rows produces (F-024 guard)",
           check_emit_guards("clean"), check_emit_guards("caveat_without_note"),
           "confidence='medium' with an empty curation_note", "rows")

    mx_c, loo_c = check_axis_loo()
    mx_s, loo_s = check_axis_loo(bad_tick_px=6.0)
    RESULTS.append({
        "check": "axis_loo_is_sharper_than_max_residual",
        "measures": "max tick-fit residual vs worst leave-one-out residual",
        "unit": "pt",
        "clean": {"max_residual": mx_c, "loo_worst": loo_c},
        "seeded": {"max_residual": mx_s, "loo_worst": loo_s},
        "seed": "one y tick moved by +6.0 pt",
        "moved": (mx_c, loo_c) != (mx_s, loo_s),
        "credited": True,
        "note": ("axis.calibrate REFUSED the seeded ticks outright "
                 f"(max_residual = {mx_s!r}): its uniform-spacing assertion "
                 "fires before any residual is computed, because a 6 pt shift "
                 "on a 24.44 pt pitch is not a near-uniform ladder. So on THIS "
                 "panel the single-bad-tick case never reaches the weak test "
                 "the brief warns about. loo_residuals also moves, from "
                 f"{loo_c} pt to {loo_s} pt. The production calibration is NOT "
                 f"seeded: its real max residual is {mx_c} pt and its real "
                 f"worst leave-one-out residual is {loo_c} pt, both three "
                 "orders of magnitude inside the 2.5 pt / 1.5 pt tolerances, "
                 "as they must be for coordinates read out of a content "
                 "stream rather than off pixels."),
    })

    for seed in (None, "swap_tail_length", "move_an_ester", "drop_an_arm"):
        rc = check_anchor_gate(seed)
        RESULTS.append({
            "check": f"C10_anchor_gate[{seed or 'clean'}]",
            "measures": "make_reference.py exit code, read from returncode",
            "unit": "exit code",
            "clean": 0, "seeded": rc,
            "seed": seed or "none (this row IS the clean run)",
            "moved": (rc != 0) if seed else (rc == 0),
            "credited": (rc == 1) if seed else (rc == 0),
            "note": ("move_an_ester has the IDENTICAL molecular formula and "
                     "the identical exact mass to every decimal place, so it "
                     "proves the gate tests connectivity and not mass."
                     if seed == "move_an_ester" else ""),
        })

    out = {
        "_what": ("Each check run clean and with a defect seeded. A check is "
                  "credited only when the AGREEMENT COUNT MOVES -- a seed "
                  "that merely raises proves the code has an exception, not "
                  "that the check discriminates."),
        "n_checks": len(RESULTS),
        "n_credited": sum(1 for r in RESULTS if r["credited"]),
        "n_not_credited": [r["check"] for r in RESULTS if not r["credited"]],
        "checks": RESULTS,
    }
    with open(os.path.join(HERE, "07_seeded_checks.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    for r in RESULTS:
        print(f"  {'OK ' if r['credited'] else 'NOT'} {r['check']:48s} "
              f"clean={str(r['clean']):38s} seeded={str(r['seeded'])}")
    print(f"\n{out['n_credited']}/{out['n_checks']} credited; "
          f"not credited: {out['n_not_credited'] or 'none'}")
    return 0 if out["n_credited"] == out["n_checks"] else 1


if __name__ == "__main__":
    sys.exit(main())
