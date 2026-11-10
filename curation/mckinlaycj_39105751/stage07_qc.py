#!/usr/bin/env python3
"""Every check, with numbers, AND a seeded defect for each. Writes 07_qc.json.

Run from this folder:   python3 -I stage07_qc.py

A NEW CHECK IS NOT DONE UNTIL IT HAS BEEN MADE TO FAIL ON PURPOSE (CLAUDE.md
section 6). Three incidents paid for that rule, and every one printed a clean
line: a re-copied C10 anchor that could only pass, three ZT_2026 marginal checks
that scored a literal 0 for twelve days while the summary beside them said "ALL
FOUR EXACT", and two versions of C16 that reported "all derived" while examining
nothing.

AND THE SEED MUST MOVE AN OUTCOME. It is not enough that a seeded defect raises
an exception -- an exception can come from the harness, from a typo in the seed,
or from a guard unrelated to the thing being tested. Every seed below reports a
METRIC BEFORE and a METRIC AFTER, and the check is only credited if the
agreement count actually drops. Where a seed is expected to raise, the metric is
still computed and reported, so "it raised" is never the whole evidence.

PASSING ALL OF THESE IS NECESSARY AND NOT SUFFICIENT. QL_2026 reproduced every
printed hit rate on bars up to 28% too high; XD_2026 landed 12/12 inside its
stated band on tops 0-21 px wrong. Nothing here digitises a value, which removes
that whole class -- but it removes it by not using the panels, not by being
careful with them.
"""
import copy
import csv
import json
import math
import os
import statistics
import sys
import xml.etree.ElementTree as ET

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
sys.path.insert(0, HERE)

import emit                                            # noqa: E402
import pymupdf                                         # noqa: E402
from rdkit import Chem, RDLogger                       # noqa: E402
from rdkit.Chem.Descriptors import ExactMolWt          # noqa: E402

import blocks                                          # noqa: E402
import make_reference                                  # noqa: E402
import seqtables                                       # noqa: E402
import stage01_blocks                                  # noqa: E402
import stage03_enumerate as enum3                      # noqa: E402
import stage05_values as val5                          # noqa: E402

RDLogger.DisableLog("rdApp.*")

PROTON = 1.00728


# ----------------------------------------------------------------- metrics

def m_exact_mass():
    """How many enumerated structures match the paper's calculated exact mass."""
    s9, order = seqtables.table_s9()
    ok = tot = 0
    for pid in order:
        if pid in enum3.SKIP:
            continue
        tot += 1
        smi = enum3.build(s9[pid]["monomers"], label=str(pid))
        mw = ExactMolWt(Chem.MolFromSmiles(smi))
        if abs(mw - s9[pid]["calc_exact_mass"]) <= enum3.MASS_TOL:
            ok += 1
    return ok, tot


def m_carbon_arithmetic():
    rows, agree = stage01_blocks.carbon_arithmetic()
    return agree, len(rows)


def m_anchor_order():
    """Anchors that reproduce the enumeration. Tests ORDER, which mass cannot."""
    order = make_reference.boxes()
    ok = tot = 0
    s9, _ = seqtables.table_s9()
    for pid, mons in order.items():
        tot += 1
        _, canon = make_reference.compose(mons)
        if canon == enum3.build(s9[pid]["monomers"], label=str(pid)):
            ok += 1
    return ok, tot


def m_s4_rejoin():
    """Table S4 rows whose rejoined sequence matches its own printed counts."""
    s4, _ = seqtables.table_s4()
    ok = sum(1 for p, r in s4.items()
             if len(r["monomers"]) == r["lipid_monomers"] + 1)
    return ok, len(s4)


def m_sup_exponent(flatten=False):
    """Table 1 values within a factor of 10 of the correct magnitude."""
    truth = {226: 2.46e7, 236: 1.08e9, 245: 2.93e9, 263: 1.06e9, 264: 1.90e9}
    root = ET.parse(val5.XML).getroot()
    tw = next(t for t in root.iter("table-wrap")
              if (t.find("label") is not None
                  and "".join(t.find("label").itertext()).strip() == "Table 1"))
    ok = tot = 0
    import re
    for tr in tw.iter("tr"):
        cells = [c for c in tr if c.tag in ("td", "th")]
        if len(cells) < 6 or cells[0].tag == "th":
            continue
        m = re.match(r"^(\d+)", "".join(cells[0].itertext()).strip())
        if not m or int(m.group(1)) not in truth:
            continue
        tot += 1
        v = cells[5]
        if flatten:
            # THE DEFECT: itertext() flattens 10<sup>7</sup> to "107".
            txt = "".join(v.itertext())
            nums = re.findall(r"[0-9.]+", txt)
            got = float(nums[0]) * float(nums[-1]) if len(nums) >= 2 else 0.0
        else:
            sup = v.find("sup")
            mm = val5.VALUE_RE.search((v.text or "").strip())
            got = float(mm.group(1)) * 10 ** int(sup.text) if (sup is not None
                                                               and mm) else 0.0
        t = truth[int(m.group(1))]
        if got > 0 and 0.1 <= got / t <= 10:
            ok += 1
    return ok, tot


def m_confidence_discriminates():
    """Distinct confidence values the rule produces over a probe set.

    A rule that returns the same string for every input is error class 29 with
    extra steps -- emit.emit_rows cannot see the difference between a derived
    constant and a literal. So the metric is DISCRIMINATION, not agreement.
    """
    probes = {
        "all_drawn": {"sequence": "Apd-Ehx-Dod-Ehx"},
        "has_name_only_monomer": {"sequence": "Apd-Ole-Dec-Ole"},
        "has_unresolved_monomer": {"sequence": "Dme-Dod-Dod-Dod-Dod"},
    }
    got = {k: val5._confidence(0, r) for k, r in probes.items()}
    return len(set(got.values())), len(probes), got


def m_found_minus_calc():
    """Rows whose printed 'found' mass is the protonated calculated mass.

    A REGULARITY THE AUTHORS DREW BUT NEVER PRINTED. SI Table S9 has columns
    'Calculated Exact Mass' and 'Found Exact Mass' and says nothing about what
    species was found. Across the table the difference is +1.01 Da with a
    spread of 0.02 -- the proton. That is external to our pipeline in a way a
    recomputed marginal is not, and the rows that MISS it turn out to be the
    paper's own data-entry errors, not ours.
    """
    s9, order = seqtables.table_s9()
    deltas = [(p, round(s9[p]["found_exact_mass"] - s9[p]["calc_exact_mass"], 3))
              for p in order]
    inband = [d for _, d in deltas if 0.95 <= d <= 1.10]
    out = [(p, d) for p, d in deltas if not 0.95 <= d <= 1.10]
    return {"n_rows": len(deltas), "n_at_plus_one_proton": len(inband),
            "mean_delta": round(statistics.mean(inband), 4),
            "sd_delta": round(statistics.pstdev(inband), 4),
            "proton_mass": PROTON,
            "outliers": out}


# ----------------------------------------------------------------- seeds

def seed(name, what, metric, mutate, restore):
    before = metric()
    raised = None
    try:
        mutate()
        after = metric()
    except Exception as e:
        after = None
        raised = f"{e.__class__.__name__}: {str(e)[:160]}"
    finally:
        restore()
    restored = metric()
    moved = (after is None) or (after[0] < before[0])
    return {"seed": name, "defect": what,
            "metric_before": f"{before[0]}/{before[1]}",
            "metric_after": (f"{after[0]}/{after[1]}" if after
                             else "check refused to run"),
            "raised": raised,
            "metric_after_restore": f"{restored[0]}/{restored[1]}",
            "outcome_moved": moved,
            "verdict": "the check has teeth" if moved else
                       "THE CHECK DID NOT MOVE -- it is not testing what it says"}


def run_seeds():
    out = []

    # 1. one CH2 removed from Dod. Composition error of 14.016 Da.
    orig = blocks.MONOMERS["Dod"]["sidechain"]
    out.append(seed(
        "dod_loses_one_carbon",
        "blocks.py Dod side chain shortened from C12 to C11 (-14.016 Da)",
        m_exact_mass,
        lambda: blocks.MONOMERS["Dod"].__setitem__("sidechain", "CCCCCCCCCCC"),
        lambda: blocks.MONOMERS["Dod"].__setitem__("sidechain", orig)))

    # 2. Hex carbon count wrong. Breaks Table S4's printed arithmetic.
    origc = blocks.MONOMERS["Hex"]["n_carbons"]
    out.append(seed(
        "hex_carbon_count_off_by_one",
        "blocks.py Hex n_carbons 6 -> 7",
        m_carbon_arithmetic,
        lambda: blocks.MONOMERS["Hex"].__setitem__("n_carbons", 7),
        lambda: blocks.MONOMERS["Hex"].__setitem__("n_carbons", origc)))

    # 3. THE IMPORTANT ONE. Transpose two monomers in the drawing-derived
    #    order. Composition is UNCHANGED, so the exact-mass check cannot see
    #    it -- this is the only check in the folder that can.
    real_boxes = make_reference.boxes

    def transposed():
        b = copy.deepcopy(real_boxes())
        for pid in b:
            # Swap the first ADJACENT PAIR THAT DIFFERS. Swapping two identical
            # monomers is a no-op, and the first version of this seed did
            # exactly that on peptoid 266 (whose positions 1 and 2 are both
            # Ehx), so the metric moved 2/2 -> 1/2 and looked like a partial
            # pass when half the seed had never been applied.
            seq = b[pid]
            i = next(k for k in range(len(seq) - 1) if seq[k] != seq[k + 1])
            seq[i], seq[i + 1] = seq[i + 1], seq[i]
        return b

    res = seed(
        "anchor_monomer_order_transposed",
        "two side chains swapped in the order read off SI p20's drawing; "
        "molecular formula and exact mass are IDENTICAL either way",
        m_anchor_order,
        lambda: setattr(make_reference, "boxes", transposed),
        lambda: setattr(make_reference, "boxes", real_boxes))
    # and show the mass check is blind to the very same defect
    s9, _ = seqtables.table_s9()
    a = enum3.build(s9[236]["monomers"], label="236")
    swapped = list(s9[236]["monomers"])
    j = next(k for k in range(1, len(swapped) - 1)
             if swapped[k] != swapped[k + 1])
    swapped[j], swapped[j + 1] = swapped[j + 1], swapped[j]
    b = enum3.build(swapped, label="236-transposed")
    res["mass_is_blind_to_this"] = {
        "same_structure": a == b,
        "exact_mass_correct_order": round(ExactMolWt(Chem.MolFromSmiles(a)), 4),
        "exact_mass_transposed": round(ExactMolWt(Chem.MolFromSmiles(b)), 4),
        "comment": "identical mass, different molecule. This is why the C10 "
                   "anchor exists and why 'HRMS passes' is not an answer."}
    out.append(res)

    # 4. the watermark column that defeated the row classifier (real bug).
    origx = seqtables.X_TABLE_MAX
    out.append(seed(
        "si_watermark_column_readmitted",
        "seqtables.X_TABLE_MAX raised so the right-margin download watermark "
        "rejoins the rows -- the ACTUAL defect found on this curation's first "
        "run, which made peptoid 246 parse as the single monomer 'Hex'",
        m_s4_rejoin,
        lambda: setattr(seqtables, "X_TABLE_MAX", 10000.0),
        lambda: setattr(seqtables, "X_TABLE_MAX", origx)))

    # 5. the superscript trap.
    before = m_sup_exponent(False)
    after = m_sup_exponent(True)
    out.append({
        "seed": "superscript_flattened",
        "defect": "Table 1 value cells read with itertext(), which turns "
                  "'2.46 +/- 0.55 x 10<sup>7</sup>' into '... x 107'",
        "metric_before": f"{before[0]}/{before[1]}",
        "metric_after": f"{after[0]}/{after[1]}",
        "raised": None,
        "metric_after_restore": f"{m_sup_exponent(False)[0]}/{before[1]}",
        "outcome_moved": after[0] < before[0],
        "verdict": "the check has teeth" if after[0] < before[0] else
                   "THE CHECK DID NOT MOVE",
        "note": "the metric is 'values within a factor of 10 of the truth', "
                "not 'it raised' -- the flattened read produces a perfectly "
                "well-formed number that is wrong by 14 orders of magnitude."})

    # 6. the confidence rule must DISCRIMINATE, not merely return a string.
    d_ok, d_tot, d_got = m_confidence_discriminates()
    real = val5._confidence
    try:
        val5._confidence = lambda p, r: "high"
        s_ok, s_tot, s_got = m_confidence_discriminates()
    finally:
        val5._confidence = real
    out.append({
        "seed": "confidence_rule_replaced_by_a_constant",
        "defect": "_confidence replaced by a function that returns 'high' for "
                  "every record -- a derived constant, which emit.emit_rows "
                  "CANNOT distinguish from a correct rule because it is a "
                  "callable",
        "metric_before": f"{d_ok}/{d_tot} distinct values over the probe set",
        "metric_after": f"{s_ok}/{s_tot} distinct values over the probe set",
        "raised": None,
        "metric_after_restore":
            f"{m_confidence_discriminates()[0]}/{d_tot} distinct values",
        "outcome_moved": s_ok < d_ok,
        "verdict": "the check has teeth" if s_ok < d_ok else
                   "THE CHECK DID NOT MOVE",
        "probe_results": d_got,
        "note": "emit.emit_rows refuses a LITERAL confidence (error class 29, "
                "found three times). It cannot refuse a callable that ignores "
                "its argument. This probe is what closes that gap."})

    # 7. emit.emit_rows still refuses an actual literal.
    recs = {1: {"sequence": "Apd-Dod-Dod-Dod"}}
    try:
        emit.emit_rows([1], recs, confidence="high", note=lambda p, r: "x",
                       label="seed")
        refused = False
    except Exception:
        refused = True
    nrows = len(emit.emit_rows([1], recs, confidence=val5._confidence,
                               note=lambda p, r: "x", label="seed"))
    out.append({
        "seed": "emit_literal_confidence",
        "defect": "confidence passed to emit.emit_rows as the string 'high'",
        "metric_before": f"{nrows}/1 rows emitted with a derived confidence",
        "metric_after": "0/1 rows emitted" if refused else "1/1 rows emitted",
        "raised": "EmitError" if refused else None,
        "metric_after_restore": f"{nrows}/1 rows emitted",
        "outcome_moved": refused,
        "verdict": "the check has teeth" if refused else
                   "THE CHECK DID NOT MOVE"})
    return out


def main():
    checks = {}

    ok, tot = m_exact_mass()
    checks["exact_mass_vs_paper_calculated"] = {
        "what": "enumerated structure's monoisotopic mass vs SI Table S9's "
                "printed Calculated Exact Mass",
        "result": f"{ok}/{tot}", "tolerance_Da": enum3.MASS_TOL,
        "blind_to": "monomer ORDER. 236/237/259/264/266 all read 1275.09."}

    ok, tot = m_carbon_arithmetic()
    checks["table_s4_carbon_arithmetic"] = {
        "what": "lipid monomer count, total lipid carbons and "
                "branched/unsaturated count recomputed from blocks.py and "
                "compared to SI Table S4's three printed columns",
        "result": f"{ok}/{tot} peptoids agree on all three columns",
        "why_it_matters": "34 equations over 6 unknown chain lengths. This, "
                          "not the three-letter code, is what pins "
                          "Hex/Oct/Dec/Dod at 6/8/10/12."}

    vec = stage01_blocks.vector_lengths()
    checks["vector_chain_lengths"] = {
        "what": "Dod and Ehx measured from the PDF's own stroked bond "
                "coordinates on SI page 20 (curation/lib chains)",
        "result": {c: f"{v['stroked_segments']} segments -> "
                      f"{v['merged_vertices']} vertices -> {v['carbons']} C"
                   for c, v in vec.items()},
        "agrees_with_blocks_py": all(
            v["carbons"] == blocks.MONOMERS[c]["n_carbons"]
            for c, v in vec.items()),
        "convention": "carbons = vertices - 1; the first vertex is the "
                      "attachment point. DERIVED: it is the only convention "
                      "that fits both chains, and they differ in length."}

    s9, _ = seqtables.table_s9()
    s4, _ = seqtables.table_s4()
    agree = sum(1 for p in s4 if s4[p]["monomers"] == s9[p]["monomers"])
    checks["table_s4_vs_s9_sequences"] = {
        "what": "the same 34 sequences printed in two different SI tables",
        "result": f"{agree}/{len(s4)} identical",
        "why_it_matters": "the only check on monomer ORDER that does not "
                          "depend on the drawing."}

    ok, tot = m_anchor_order()
    checks["c10_anchor"] = {
        "what": "peptoids 236 and 266 rebuilt from SI Figure S7A's DRAWING "
                "(fill colour + x centre of the side-chain boxes), C-terminus "
                "to N-terminus, by string composition, and compared to the "
                "Table-S9 enumeration",
        "result": f"{ok}/{tot} canonical SMILES identical",
        "independence": "different input (drawing geometry vs text table), "
                        "different traversal (C->N vs N->C), different "
                        "construction (string vs assemble.attach atom maps), "
                        "blocks.py not imported by make_reference.py",
        "box_spacing_pt": "the six boxes of each structure are evenly spaced "
                          "in x to within 2 pt, asserted before composing"}

    checks["found_minus_calculated_is_a_proton"] = m_found_minus_calc()
    checks["found_minus_calculated_is_a_proton"]["what"] = (
        "A ROUND NUMBER THE AUTHORS DREW BUT NEVER PRINTED. SI Table S9 never "
        "says what species 'Found Exact Mass' is. The difference is the "
        "proton, on 35 of 45 rows, to within 0.02 Da.")
    checks["found_minus_calculated_is_a_proton"]["what_the_misses_are"] = (
        "the paper's own data-entry errors, not ours -- 253's found mass is "
        "byte-identical to 245's and 255's, and 262's is printed as '2.00'. "
        "See 08_contradictions.csv. None of them touches a structure we "
        "emit, because we use the CALCULATED column.")

    b, t = m_sup_exponent(False)
    checks["table1_superscript"] = {
        "what": "all five Table 1 values parsed with an explicit <sup> "
                "exponent and within a factor of 10 of their magnitude",
        "result": f"{b}/{t}"}

    with open(os.path.join(HERE, "05_values.csv")) as fh:
        vrows = list(csv.DictReader(fh))
    logs = [math.log10(float(r["value_mean_p_per_s"])) for r in vrows]
    mu = sum(logs) / len(logs)
    sd = math.sqrt(sum((x - mu) ** 2 for x in logs) / len(logs))
    zs = [float(r["Experiment_value_z"]) for r in vrows]
    checks["z_score_round_trip"] = {
        "what": "the stored z recomputed from the stored raw value",
        "result": f"max |residual| = "
                  f"{max(abs(z - (l - mu) / sd) for z, l in zip(zs, logs)):.2e}",
        "mean_z": round(sum(zs) / len(zs), 12),
        "sd_z": round(math.sqrt(sum(z * z for z in zs) / len(zs)), 12),
        "n": len(zs),
        "caveat": "n = 5. A five-point z-score is a legitimate application of "
                  "the corpus convention and a thin one; see DECISIONS.md D-2."}

    d_ok, d_tot, d_got = m_confidence_discriminates()
    checks["flags_on_the_row"] = {
        "what": "quality columns written by curation/lib emit.emit_rows, "
                "never by this folder",
        "rows": len(vrows),
        "rows_with_a_caveat_F024": sum(
            1 for r in vrows
            if r["confidence"] in ("medium", "low")
            or r["positional_ambiguity"] == "1"
            or r["structure_assumed"] == "1"
            or r["value_is_bound"] in ("upper", "lower")),
        "confidence_rule_distinct_outputs": f"{d_ok}/{d_tot}",
        "confidence_rule_probe": d_got,
        "note": "0 caveats is the correct answer for these five: every "
                "monomer they use is drawn in Figure 1B, the sequences are "
                "confirmed by two tables, the masses by a third, and 236 by "
                "the drawing. The rule still discriminates -- see the probe "
                "and seed 6."}

    report = {
        "_what": "Every check this curation makes, with numbers, and a seeded "
                 "defect for each that is credited only if an agreement count "
                 "MOVES.",
        "checks": checks,
        "seeded_defects": run_seeds(),
    }
    report["all_seeds_moved"] = all(s["outcome_moved"]
                                    for s in report["seeded_defects"])

    out = os.path.join(HERE, "07_qc.json")
    with open(out, "w") as fh:
        json.dump(report, fh, indent=1)

    for k, v in checks.items():
        print(f"{k}: {v.get('result', '')}")
    print()
    for s in report["seeded_defects"]:
        print(f"SEED {s['seed']}: {s['metric_before']} -> "
              f"{s['metric_after']}  [{s['verdict']}]")
    print(f"\nall seeds moved an outcome: {report['all_seeds_moved']}")
    print(f"wrote {out}")
    return 0 if report["all_seeds_moved"] else 1


if __name__ == "__main__":
    sys.exit(main())
