#!/usr/bin/env python3
"""Stage 07 -- every check, with numbers, and every check SEEDED to confirm it
fires.

"A NEW CHECK IS NOT DONE UNTIL IT HAS BEEN MADE TO FAIL ON PURPOSE."
(CLAUDE.md section 6.) Three incidents paid for that rule, and one of them --
three ZT_2026 marginal checks scoring a literal 0 for twelve days while the
summary beside them said "ALL FOUR EXACT" -- is exactly what a QC file full of
clean lines looks like when nothing is actually running. So each check below
runs twice: once on the real data, and once on a copy with the defect it exists
to catch deliberately injected. A check whose seeded run does not fail is
reported as `teeth: false` and is worth nothing.

PASSING EVERY PRINTED STATISTIC IS NECESSARY, NOT SUFFICIENT. For this paper
the usual trap is inverted: there is no digitised quantity to be
multiplicatively wrong, because this curation emits no values at all. The
corresponding trap here is that a molecular formula is blind to connectivity --
every one of P1/P2/P3 would pass on a structure with the right atoms wired the
wrong way. P9 (the C10 anchors) is the only check below that tests topology
against something external, and it covers 2 of 74 lipids.

    python3 -I stage07_qc.py        # writes 07_qc.json; exit 1 on a failure
"""
import copy
import csv
import importlib
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "lib"))

import numpy as np                                            # noqa: E402
from PIL import Image                                         # noqa: E402
from rdkit import Chem                                        # noqa: E402

import stage03_enumerate as S3                                # noqa: E402
from errors import CurationLibError                           # noqa: E402

SI = os.path.join(HERE, "source", "ja5c13982_si_001.pdf")
FIG_S51_PAGE = 90          # 0-based page index; SI label "S91"


# ----------------------------------------------------------------- P1-P3 ----
def hrms_checks():
    rows = list(csv.DictReader(open(os.path.join(HERE, "03_enumerated.csv"))))
    si = [r for r in rows if r["formula_stated"]]
    f_ok = [r["name"] for r in si if r["formula_agrees"] == "1"]
    m_ok = [r["name"] for r in si if r["mz_agrees"] == "1"]
    both = [r["name"] for r in si
            if r["formula_agrees"] == "1" and r["mz_agrees"] == "1"]
    self_ok = [r["name"] for r in si if r["si_self_consistent"] == "1"]
    neither = [r["name"] for r in si
               if r["formula_agrees"] != "1" and r["mz_agrees"] != "1"]
    return {
        "n_lipids_enumerated": len(rows),
        "n_with_an_SI_HRMS_record": len(si),
        "n_without_any_SI_record": len(rows) - len(si),
        "P1_formula_agrees": f"{len(f_ok)}/{len(si)}",
        "P2_mh_plus_within_0.05_Da": f"{len(m_ok)}/{len(si)}",
        "P1_and_P2_both": f"{len(both)}/{len(si)}",
        "P3_SI_formula_agrees_with_SI_own_mz": f"{len(self_ok)}/{len(si)}",
        "P3_note": ("Entirely internal to the SI: the printed formula's "
                    "monoisotopic [M+H]+ against the printed [M+H]+. Where "
                    "these disagree the HRMS line is not evidence about "
                    "anything, which is how 5 of the 8 structural "
                    "disagreements below are adjudicated."),
        "lipids_failing_P1_or_P2": [r["name"] for r in si
                                    if r["name"] not in both],
        "lipids_where_NEITHER_SI_locus_confirms_us": neither,
        "max_abs_mh_delta_over_passing": round(max(
            float(r["mh_plus_delta_vs_SI"]) for r in si
            if r["name"] in both), 4),
    }


# -------------------------------------------------------------------- P4 ----
def _tick_run_widths(img, y_axis_band=(8, 22), thresh=100, cover=0.8):
    ink = img < thresh
    rowdark = ink.sum(1)
    ax = int(np.argmax(rowdark))
    band = ink[ax + y_axis_band[0]: ax + y_axis_band[1]]
    col = band.mean(0)
    on = col > cover
    runs, i = [], 0
    while i < len(on):
        if on[i]:
            j = i
            while j < len(on) and on[j]:
                j += 1
            runs.append((i, j))
            i = j
        else:
            i += 1
    return runs


def roster_from_figure(deface=False):
    """Count Figure S51's x-axis ticks: the library size the SI never prints.

    The SI states no library size anywhere -- the string "315" does not occur
    in it, nor does the word "library". So the count of bars the authors DREW
    is the only figure the SI supplies, and it is external to our text parse in
    the way CLAUDE.md section 5 asks for: a number drawn and never printed.
    """
    import pymupdf
    d = pymupdf.open(SI)
    pm = d[FIG_S51_PAGE].get_pixmap(dpi=600,
                                    clip=pymupdf.Rect(90, 230, 530, 250))
    img = np.array(Image.frombytes("RGB", (pm.width, pm.height),
                                   pm.samples).convert("L"))
    if deface:
        # SEED: erase one tick. A counter that is really counting must drop.
        runs = _tick_run_widths(img)
        full = [r for r in runs if r[1] - r[0] >= 6]
        a, b = full[len(full) // 2]
        ink_rows = (img < 100).sum(1)
        ax = int(np.argmax(ink_rows))
        img[ax + 6: ax + 30, a - 1: b + 1] = 255
    runs = _tick_run_widths(img)
    widths = [b - a for a, b in runs]
    full = [(a, b) for a, b in runs if b - a >= 6]
    cen = np.array([(a + b) / 2 for a, b in full])
    dif = np.diff(cen)
    return {"n_ticks_full_width": len(full),
            "n_runs_total": len(runs),
            "partial_runs_rejected": sorted(w for w in widths if w < 6),
            "tick_pitch_px_mean": round(float(dif.mean()), 3),
            "tick_pitch_px_sd": round(float(dif.std()), 4),
            "tick_pitch_cv": round(float(dif.std() / dif.mean()), 5)}


# -------------------------------------------------------------------- P7 ----
def fluorine_check():
    """The F subscript of the SI's printed formula vs the F count we built.

    Fluorine is the point of this paper, and a wrong fluoroalkyl is a defect
    that a chain-length count cannot see at all: 4-fluorobutyl, 4,4,4-
    trifluorobutyl, 3,3,4,4,4-pentafluorobutyl and 2,2,3,3,4,4,4-
    heptafluorobutyl are all FOUR carbons. They are distinguishable only by
    the F subscript, and swapping two of them changes the formula by
    H-2/F+2 -- a 35.98 Da shift that would otherwise look like an ordinary
    homologue step. So this is reported on its own rather than folded into
    P1, where a carbon error and a fluorine error would be one number.

    The POSITION of the fluorines is a separate matter and is NOT checked
    here: F5 as 3,3,4,4,4- and a hypothetical 2,2,3,3,4- are isomers with
    identical formula and identical exact mass. Position comes from the
    purchased reagent's own IUPAC name in the SI Materials list (01_blocks.csv),
    which is the only locus in the paper that states it.
    """
    import re as _re
    rows = list(csv.DictReader(open(os.path.join(HERE, "03_enumerated.csv"))))
    si = [r for r in rows if r["formula_stated"]]
    agree, disagree = 0, []
    for r in si:
        m = _re.search(r"F(\d*)", r["formula_stated"])
        stated = (int(m.group(1)) if m.group(1) else 1) if m else 0
        if stated == int(r["n_F"]):
            agree += 1
        else:
            disagree.append({"name": r["name"], "stated_F": stated,
                             "built_F": int(r["n_F"])})
    tokens = {t: sorted({int(r["n_F"]) for r in rows if r["fluoro"] == t})
              for t in ("F1", "F3", "F5", "F7")}
    return {"P7_F_subscript_agrees_with_built_F_count": f"{agree}/{len(si)}",
            "P7_disagreements": disagree,
            "P7_Fn_token_to_built_F_count": tokens,
            "P7_all_four_tails_are_four_carbons": True,
            "P7_what_it_cannot_see": ("fluorine POSITION -- 3,3,4,4,4- and "
                                      "2,2,3,3,4- pentafluorobutyl have the "
                                      "same formula and the same exact mass. "
                                      "Position rests on the Materials list's "
                                      "IUPAC reagent names alone.")}


# -------------------------------------------------------------- P5/P6/P8 ----
def name_token_report():
    nt = json.load(open(os.path.join(HERE, "03_name_token_checks.json")))
    cn = {k: v for k, v in nt.items() if v["name_token_length"] is not None}
    bad = [k for k, v in cn.items() if not v["matches"]]
    other = {}
    for k, v in nt.items():
        if v["name_token_length"] is None:
            other.setdefault(v["note"], []).append(k)
    return {"P5_Cn_token_equals_measured_chain":
            f"{len(cn) - len(bad)}/{len(cn)}",
            "P5_mismatches": bad,
            "P5_tokens_that_are_NOT_chain_lengths":
            {k: len(v) for k, v in other.items()}}


# -------------------------------------------------------------------- P9 ----
def anchors():
    enum = {r["name"]: r["smiles"]
            for r in csv.DictReader(open(os.path.join(HERE,
                                                      "03_enumerated.csv")))}
    out = []
    for r in csv.DictReader(open(os.path.join(HERE,
                                              "reference_structures.csv"))):
        a = Chem.MolFromSmiles(r["reference_SMILES"])
        b = Chem.MolFromSmiles(enum[r["lipid_name"]])
        out.append({"lipid": r["lipid_name"],
                    "inchikey_anchor": Chem.MolToInchiKey(a),
                    "inchikey_enumerated": Chem.MolToInchiKey(b),
                    "exact_match": Chem.MolToInchiKey(a) ==
                    Chem.MolToInchiKey(b),
                    "source": r["source_figure"]})
    return out


# ------------------------------------------------------------------- P10 ----
def duplicates_and_isomers():
    rows = list(csv.DictReader(open(os.path.join(HERE, "03_enumerated.csv"))))
    by_smi, by_formula = {}, {}
    for r in rows:
        by_smi.setdefault(r["smiles"], []).append(r["name"])
        by_formula.setdefault((r["formula"], r["exact_mass"]),
                              []).append(r["name"])
    dups = {k: v for k, v in by_smi.items() if len(v) > 1}
    iso = {f"{k[0]} / {k[1]}": v for k, v in by_formula.items() if len(v) > 1}
    return {"P10_unique_canonical_SMILES": f"{len(by_smi)}/{len(rows)}",
            "P10_duplicate_structures": dups,
            "P11_distinct_structures_sharing_formula_AND_exact_mass": iso,
            "P11_note": ("Mass cannot see connectivity, demonstrated inside "
                         "this one library: these pairs are different "
                         "molecules with identical formula and identical "
                         "exact mass to every decimal, so no HRMS line in the "
                         "SI can tell them apart. AL vs AS is the clean case "
                         "-- 2-hexyldecanoate against 2-heptylnonanoate -- "
                         "and the SI prints the SAME calculated m/z for both "
                         "members of each pair.")}


# ------------------------------------------------------- seeded mutations ---
def _restore():
    """Rebuild 03_enumerated.csv AND re-apply stage04's quality flags.

    stage03.main() rewrites 03_enumerated.csv with the structural columns only,
    so restoring after a seeded mutation by calling it alone SILENTLY DROPPED
    every quality flag -- the F-001 failure mode, reintroduced by the very
    machinery that exists to prove the checks work. Caught by the guard at the
    end of main(), which is why that guard is there and not merely implied.
    """
    S3.main()
    import stage04_flags
    importlib.reload(stage04_flags)
    stage04_flags.main()


def _agreement():
    """(formula_agrees, mz_agrees) counts, read off the artefact on disk."""
    rows = [r for r in csv.DictReader(open(os.path.join(HERE,
                                                        "03_enumerated.csv")))
            if r["formula_stated"]]
    return (sum(1 for r in rows if r["formula_agrees"] == "1"),
            sum(1 for r in rows if r["mz_agrees"] == "1"))


BASELINE_AGREEMENT = None          # set in main(), before anything is seeded


def seed(label, mutate, expect):
    """Inject the defect this check exists to catch, and prove it is caught.

    TEETH ARE NOT THE SAME AS AN EXCEPTION, and the first version of this
    function assumed they were. It required `mutate` to make stage03 RAISE,
    which meant:

      * the P1/P2 seed scored `teeth: true` for the wrong reason -- lengthening
        the C4 arm trips the LINKER_MAIN_CHAIN assertion, which raises, so the
        formula comparison was never shown to catch anything;
      * the P7 seed -- building every F5 member with the F3 tail, four carbons
        either way, two fluorines fewer -- scored `teeth: FALSE`, because
        stage03 reports a formula mismatch and exits 0. It does not gate, and
        deliberately so: eight real disagreements with the SI are open
        decisions (DECISIONS.md D-3) and must not block the build.

    That is the exact shape of the three ZT_2026 marginal checks that scored a
    literal 0 for twelve days beside a summary reading "ALL FOUR EXACT". So a
    seeded run now counts as caught if it raises OR if the agreement counts
    measured off the rebuilt artefact fall. Both outcomes are recorded.
    """
    importlib.reload(S3)
    saved = copy.deepcopy((S3.ARMS, S3.LINKER_MAIN_CHAIN, S3.TAILS))
    raised, msg, after = False, "", None
    try:
        mutate(S3)
        try:
            S3.build("SSC4F1")
            S3.measured_main_chain("SSC5F1")
            S3.main()
            after = _agreement()
        except Exception as e:
            raised, msg = True, f"{type(e).__name__}: {str(e)[:180]}"
    finally:
        S3.ARMS, S3.LINKER_MAIN_CHAIN, S3.TAILS = saved
        importlib.reload(S3)
        _restore()
    dropped = (after is not None and BASELINE_AGREEMENT is not None
               and (after[0] < BASELINE_AGREEMENT[0]
                    or after[1] < BASELINE_AGREEMENT[1]))
    return {"seeded": label, "expect": expect,
            "stage03_raised": raised, "message": msg,
            "agreement_before": list(BASELINE_AGREEMENT or ()),
            "agreement_after": list(after) if after else None,
            "agreement_dropped": dropped,
            "teeth": raised or dropped,
            "caught_by": ("an assertion in stage03" if raised else
                          "the formula/mass agreement counts falling"
                          if dropped else "NOTHING")}


def _m_arm(mod):
    mod.ARMS["C4"] = mod.ARMS["C4"].replace("OCCCCC(=O)", "OCCCCCC(=O)")


def _m_token(mod):
    mod.LINKER_MAIN_CHAIN["C5"] = 9


def _m_fluoro(mod):
    mod.TAILS["F5"]["ester"] = mod.TAILS["F3"]["ester"]    # F5 built as F3


def _m_ring(mod):
    mod.TAILS["F1"]["ester"] = "CCCC1CC1F"          # a ring digit in the tail


def main():
    qc = {
        "_what": ("Quality control for curation/liu_41162324 "
                  "(Liu et al., JACS 2025, PMID 41162324, "
                  "DOI 10.1021/jacs.5c13982). Structures only: this curation "
                  "emits NO Experiment_value rows. See 00_value_sources.json."),
        "paper": {"pmid": 41162324, "doi": "10.1021/jacs.5c13982",
                  "journal": "J Am Chem Soc 147(45), 2025",
                  "source_held": "ja5c13982_si_001.pdf only (100 pages); the "
                                 "article PDF was not supplied and is not "
                                 "open access"},
        "P1_P2_P3_hrms": hrms_checks(),
        "P4_roster_drawn_in_figure_S51": roster_from_figure(),
        "P5_name_tokens": name_token_report(),
        "P7_fluorine": fluorine_check(),
        "P9_C10_anchors": anchors(),
        "P10_P11_structures": duplicates_and_isomers(),
    }

    # --- the claim P4 exists to test ---
    n = qc["P4_roster_drawn_in_figure_S51"]["n_ticks_full_width"]
    qc["P4_roster_drawn_in_figure_S51"].update({
        "interpretation": f"{n} ticks = 1 benchmark (SM102) + {n - 1} FILs",
        "n_FILs_drawn": n - 1,
        "n_FILs_enumerated": qc["P1_P2_P3_hrms"]["n_lipids_enumerated"],
        "agrees_with_enumeration":
            n - 1 == qc["P1_P2_P3_hrms"]["n_lipids_enumerated"],
        "abstract_claims": 315,
        "reconciles_with_abstract": False,
        "why": ("The string '315' occurs nowhere in the 100-page SI, nor do "
                "'library', 'combinatorial' or 'total of'. The SI's own "
                "roster -- the x-axis of Figures S48-S51, the synthesis "
                "headings of section 3.2, and the NMR/MS figure captions "
                "S12-S47 -- is 74 FILs. The abstract's 315 cannot be "
                "reconciled from the material held; it needs the article "
                "text. Open: DECISIONS.md D-1."),
    })

    global BASELINE_AGREEMENT
    BASELINE_AGREEMENT = _agreement()

    qc["seeded_mutations"] = [
        seed("P1/P2: lengthen the C4 arm template by one CH2", _m_arm,
             "the C4 members' built formula and [M+H]+ stop matching the SI, "
             "and the main-chain assertion contradicts LINKER_MAIN_CHAIN"),
        seed("P5: set LINKER_MAIN_CHAIN['C5'] to 9", _m_token,
             "stage03 refuses: a hand-counted expectation and the measured "
             "structure disagree, and neither is allowed to win silently"),
        seed("P7: build every F5 member with the F3 tail (same 4 carbons, "
             "two fluorines fewer)", _m_fluoro,
             "the F5 members' formula and [M+H]+ stop matching the SI -- the "
             "defect a chain-length count is blind to"),
        seed("P6: put a ring-closure digit in the F1 tail", _m_ring,
             "_check_no_ring_digit_collision refuses before anything parses"),
    ]

    # P8: the SM-102 substructure assertion in make_reference.py.
    mr = os.path.join(HERE, "make_reference.py")
    src = open(mr).read()
    bad = src.replace('SM102_SHARED = "OCCN(',
                      'SM102_SHARED = "OCCCN(')          # one CH2 too many
    tmp = os.path.join(HERE, "_seeded_make_reference.py")
    try:
        open(tmp, "w").write(bad)
        r = subprocess.run([sys.executable, "-I", tmp], capture_output=True,
                           text=True, cwd=HERE)
        qc["seeded_mutations"].append({
            "seeded": "P8: lengthen the SM-102 shared fragment by one CH2",
            "expect": "make_reference.py refuses: the registry's SM-102 no "
                      "longer contains the fragment this library is built on",
            "fired": r.returncode != 0,
            "message": (r.stdout + r.stderr).strip()[-200:],
            "teeth": r.returncode != 0})
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)

    qc["seeded_mutations"].append({
        "seeded": "P4: erase one x-axis tick from the Figure S51 raster",
        "expect": "the tick count drops by exactly 1",
        "n_ticks_after": roster_from_figure(deface=True)["n_ticks_full_width"],
        "fired": roster_from_figure(deface=True)["n_ticks_full_width"] ==
                 n - 1,
        "teeth": roster_from_figure(deface=True)["n_ticks_full_width"] ==
                 n - 1})

    qc["seeded_mutations"].append({
        "seeded": "C16 (provenance/check_anchor_independence.py) on the first "
                  "version of make_reference.py",
        "expect": "C16 names the transcribed literal",
        "fired": True,
        "teeth": True,
        "message": ("Not a hypothetical: the first version of this folder's "
                    "make_reference.py typed SM102_SHARED nitrogen-tail-first "
                    "and C16 reported '1 TEMPLATE(S) APPEAR TO BE TRANSCRIBED "
                    "FROM THE ENUMERATION: curation/liu_41162324/"
                    "make_reference.py:SM102_SHARED'. Re-typed hydroxyl-first; "
                    "C16 now reports 'all derived'. The alarm was observed to "
                    "fire on this folder before it was trusted on it.")})

    qc["what_these_checks_CANNOT_see"] = [
        "Connectivity. P1/P2/P3 are atom counts; a molecule with the right "
        "formula and the wrong wiring passes all three. Only P9 tests "
        "topology against the paper, and it covers 2 of 74 lipids.",
        "Which nitrogen is alkylated -- there is only one N-H in each amine "
        "intermediate, so this is pinned by the SI's name for that "
        "intermediate rather than by any measurement here.",
        "The SP series entirely: 32 of 74 lipids have no formula, no m/z, no "
        "individual synthesis paragraph and no drawn scheme anywhere in the "
        "SI, and the SI names their core alcohol two different ways.",
        "Any transfection value. None is held; see 00_value_sources.json.",
        "Whether Figure S51's bar ORDER matches its label order -- the labels "
        "were read, the bars were not paired to them, because no value from "
        "that panel is emitted.",
    ]

    # GUARD: the seeded runs rebuild 03_enumerated.csv. Every quality column
    # must still be on it when this stage finishes, or the flags died one
    # stage before the merge -- which is precisely F-001.
    import emit as _emit
    hdr = next(csv.reader(open(os.path.join(HERE, "03_enumerated.csv"))))
    missing = [c for c in _emit.columns() if c not in hdr]
    if missing:
        print("FATAL: quality columns missing from 03_enumerated.csv:",
              missing)
        return 1
    rows = list(csv.DictReader(open(os.path.join(HERE,
                                                 "03_enumerated.csv"))))
    qc["flags_on_the_row"] = {
        "columns": _emit.columns(),
        "rows_with_a_caveat": _emit.caveat_count(rows),
        "confidence_counts": {c: sum(1 for r in rows if r["confidence"] == c)
                              for c in ("high", "medium", "low")},
        "structure_assumed_1": sum(1 for r in rows
                                   if r["structure_assumed"] == "1"),
        "note": ("Written by curation/lib/emit.emit_rows, which refuses a "
                 "literal confidence (error class 29) and refuses a caveat "
                 "with no note. No value row is emitted by this curation, so "
                 "value_is_bound and value_normalization_unit are left blank "
                 "-- 'not stated', which the schema defines as an absence of "
                 "information rather than an assurance."),
    }

    out = os.path.join(HERE, "07_qc.json")
    json.dump(qc, open(out, "w"), indent=1)
    toothless = [m["seeded"] for m in qc["seeded_mutations"]
                 if not m.get("teeth")]
    print(json.dumps({k: v for k, v in qc["P1_P2_P3_hrms"].items()
                      if not isinstance(v, list)}, indent=1))
    print("roster ticks:", qc["P4_roster_drawn_in_figure_S51"][
        "n_ticks_full_width"], "pitch cv",
        qc["P4_roster_drawn_in_figure_S51"]["tick_pitch_cv"])
    print("anchors:", [(a["lipid"], a["exact_match"])
                       for a in qc["P9_C10_anchors"]])
    for m in qc["seeded_mutations"]:
        print(f"  seeded[{m['teeth']}] {m['seeded']}")
    print(f"-> {out}")
    if toothless:
        print("TOOTHLESS CHECKS:", toothless)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
