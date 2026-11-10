#!/usr/bin/env python3
"""Stage 07 - QC with numbers. Run: python3 -I 07_qc.py   (writes 07_qc.json)

Exit code is the gate: 0 only if every check passes. VERIFY THE EXIT CODE
DIRECTLY, not through a pipe -- `python3 -I 07_qc.py | tail` reports tail's
status, which is how a sibling's gate was masked.

THE MASS CHECK IS SPLIT THREE WAYS, plus a fourth that the caption cannot
reach. Collapsing them into one boolean is how a published formula typo gets
adopted:

  Q1  our formula (from the structure the CORPUS stores for 6A1-SC8)
        vs the SI caption's printed formula string
  Q2  the SI caption's printed formula string
        vs the SI caption's own printed m/z
  Q3  our formula  vs the SI caption's printed m/z
  Q4  our formula  vs the ISOTOPE LABELS DRAWN INSIDE Fig S1B

Q4 is the valuable one and it is the "round number the authors DREW but never
printed" class. The caption says "found: 2170.3 Da"; the spectrum itself is
labelled 1086.1679, 1086.6692, 1087.1695, 1087.6700, 1088.1705. The 0.5 Da
spacing says z = 2 -- a fact the caption never states -- and 2170.3 is a
deconvoluted number, not a measurement. Q4 therefore tests the formula against
the instrument's own output, through a route the caption plays no part in.

WHAT THIS DOES NOT TEST. Mass cannot see connectivity. Q1-Q4 passing confirms
C110H203N5O24S6 and, with the NMR ratios in Q6, the per-arm composition. It
says nothing about which nitrogen of the polyamine core bears which arm. The SI
DRAWS NO STRUCTURE anywhere, so this paper supplies no connectivity ground
truth at all; see reference_structures.csv and D-4.
"""
import json
import os
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "lib"))

import numpy as np
import pymupdf
from rdkit import Chem, RDLogger
from rdkit.Chem.rdMolDescriptors import CalcExactMolWt, CalcMolFormula

import axis  # curation/lib/axis.py -- never a local re-implementation

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
PDF = os.path.join(HERE, "source", "pnas.2528477123.sapp.pdf")
PT = Chem.GetPeriodicTable()
PROTON = 1.00727646688  # mass of H+ (1H minus one electron)

# ---------------------------------------------------------------------------
# INPUTS, each tagged with where it was read and how.
# ---------------------------------------------------------------------------

# The structure the CORPUS already stores for 6A1-SC8. Not enumerated here:
# this paper enumerates nothing. Pulled from data/data_vivo.csv (KZ_2016) by
# 02_corpus_overlap.py and pasted verbatim; 02 re-derives it from the live file
# every run, so a drift would show up there.
CORPUS_SMILES = (
    "CCCCCCCCSCC(C)C(=O)OCCOC(=O)CCN(CCNCCN(CCC(=O)OCCOC(=O)C(C)CSCCCCCCCC)"
    "CCN(CCC(=O)OCCOC(=O)C(C)CSCCCCCCCC)CCN(CCC(=O)OCCOC(=O)C(C)CSCCCCCCCC)"
    "CCC(=O)OCCOC(=O)C(C)CSCCCCCCCC)CCC(=O)OCCOC(=O)C(C)CSCCCCCCCC")

# SI p.4 caption, PDF TEXT LAYER (word tokens at y=564.4): "C110H203N5O24S6:"
# "2170.31" "Da," "found:" "2170.3" "Da". Read symbolically, not off pixels.
SI_FORMULA_STRING = "C110H203N5O24S6"
SI_PRINTED_CALC_MZ = 2170.31
SI_PRINTED_FOUND_MZ = 2170.3

# Fig S1B isotope-cluster labels. These live inside an EMBEDDED RASTER, not in
# the text layer, so they were read by eye off a 900 dpi crop of the unscaled
# page (pymupdf get_pixmap dpi=900, clip in PDF points). A sixth label,
# 1088.6710, was legible only at 200 dpi and is deliberately EXCLUDED -- a
# number read at a resolution where the glyphs were not separable is not
# evidence.
MS_LABELS = [1086.1679, 1086.6692, 1087.1695, 1087.6700, 1088.1705]
MS_LABELS_READ_AT_DPI = 900

# SI p.4, Fig S1A. Integral values and the ppm-axis tick labels are BOTH in the
# PDF text layer, with coordinates, so the integral -> ppm assignment is
# measured by axis.calibrate rather than asserted. See q5/q6.
NMR_PAGE = 3  # 0-based

# Predicted proton counts for the corpus structure, per region. Derived from
# the SMILES by rdkit in q6 -- NOT typed here -- except for the ppm windows,
# which are the chemistry.
REGIONS = {
    "ester_glycol_OCH2CH2O": dict(ppm_lo=4.0, ppm_hi=4.6),
    "octyl_beta_to_S_CH2":   dict(ppm_lo=1.45, ppm_hi=1.70),
    "octyl_omega_CH3":       dict(ppm_lo=0.75, ppm_hi=1.00),
}


# ---------------------------------------------------------------------------
def parse_formula(s):
    out = {}
    el = ""
    num = ""
    for ch in s + "$":
        if ch.isupper() or ch == "$":
            if el:
                out[el] = out.get(el, 0) + (int(num) if num else 1)
            el, num = ch, ""
        elif ch.islower():
            el += ch
        elif ch.isdigit():
            num += ch
        else:
            raise ValueError(f"parse_formula: unexpected {ch!r} in {s!r}")
    out.pop("$", None)
    return out


def exact_mass_of_formula(d):
    return sum(PT.GetMostCommonIsotopeMass(el) * n for el, n in d.items())


def nmr_text(page):
    """Integral labels and ppm tick labels, with x centroids, from the text
    layer. Returns (integrals, ticks) as lists of (x_centre, value)."""
    words = page.get_text("words")
    integrals, ticks = [], []
    for x0, y0, x1, y1, t, *_ in words:
        xc = (x0 + x1) / 2.0
        if 309.0 < y0 < 314.0:                       # the integral row(s)
            integrals.append((xc, float(t)))
        elif 322.5 < y0 < 324.5 and t.count(".") == 1:  # the f1 (ppm) axis
            # The two outermost labels are clipped by the panel ("0." for 0.0
            # and a truncated "10.0" rendered "0.0" at the far left), so their
            # glyph centroids do not sit on their ticks. Excluded by value:
            # only the 0.5..9.5 ladder is used.
            try:
                v = float(t.rstrip("."))
            except ValueError:
                continue
            if 0.5 <= v <= 9.5:
                ticks.append((xc, v))
    return integrals, ticks


# ---------------------------------------------------------------------------
def run(smiles=CORPUS_SMILES, si_formula=SI_FORMULA_STRING,
        ms_labels=None, region_defs=None):
    """All checks. Arguments exist so the seeded mutations in `selftest` can
    perturb one input at a time and watch an AGREEMENT COUNT move."""
    ms_labels = list(MS_LABELS if ms_labels is None else ms_labels)
    regions = dict(REGIONS if region_defs is None else region_defs)
    R = {}

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise SystemExit("07_qc: the corpus SMILES does not parse.")
    ours_formula_str = CalcMolFormula(mol)
    ours = parse_formula(ours_formula_str)
    theirs = parse_formula(si_formula)
    ours_exact = CalcExactMolWt(mol)

    # -- Q1 --------------------------------------------------------------
    els = sorted(set(ours) | set(theirs))
    per_el = {e: {"ours": ours.get(e, 0), "SI_string": theirs.get(e, 0),
                  "agree": ours.get(e, 0) == theirs.get(e, 0)} for e in els}
    R["q1_our_formula_vs_SI_printed_string"] = {
        "_compares": "the structure the corpus stores for 6A1-SC8, against the "
                     "formula STRING printed in the Fig S1 caption",
        "ours": ours_formula_str,
        "SI_string": si_formula,
        "per_element": per_el,
        "n_elements": len(els),
        "agreement_count": sum(1 for e in els if per_el[e]["agree"]),
        "pass": all(per_el[e]["agree"] for e in els),
    }

    # -- Q2 --------------------------------------------------------------
    si_exact = exact_mass_of_formula(theirs)
    d2 = abs(si_exact - SI_PRINTED_CALC_MZ)
    R["q2_SI_string_vs_SI_printed_calc"] = {
        "_compares": "the SI's OWN formula string against the SI's OWN printed "
                     "'Calc.' m/z -- an internal consistency check on the "
                     "caption, in which we take no part",
        "monoisotopic_of_SI_string": round(si_exact, 4),
        "SI_printed_calc_mz": SI_PRINTED_CALC_MZ,
        "abs_delta_Da": round(d2, 4),
        "tol_Da": 0.005,
        "tol_why": "the caption prints 2 decimals, so |delta| <= 0.005 is pure "
                   "rounding and anything larger is a real disagreement",
        "pass": d2 <= 0.005,
    }

    # -- Q3 --------------------------------------------------------------
    d3 = abs(ours_exact - SI_PRINTED_FOUND_MZ)
    R["q3_our_formula_vs_SI_printed_found"] = {
        "_compares": "our monoisotopic mass against the SI's printed 'found' "
                     "value, which is a DECONVOLUTED neutral mass, not an "
                     "observed m/z (see q4)",
        "ours_monoisotopic": round(ours_exact, 4),
        "SI_printed_found_mz": SI_PRINTED_FOUND_MZ,
        "abs_delta_Da": round(d3, 4),
        "tol_Da": 0.05,
        "tol_why": "the caption prints 1 decimal for 'found', so the printed "
                   "value itself carries +/-0.05 Da of quantisation",
        "pass": d3 <= 0.05,
    }

    # -- Q4 ---------------------------------------------------------------
    sp = [round(b - a, 4) for a, b in zip(ms_labels, ms_labels[1:])]
    mean_sp = float(np.mean(sp)) if sp else float("nan")
    z = int(round(1.0 / mean_sp)) if sp and mean_sp else 0
    base = (ours_exact + z * PROTON) / z if z else float("nan")
    # 13C isotopologues step by (13C-12C)/z
    c13 = PT.GetMassForIsotope("C", 13) - PT.GetMostCommonIsotopeMass("C")
    pred = [base + i * c13 / z for i in range(len(ms_labels))] if z else []
    ppm_err = [abs(o - p) / p * 1e6 for o, p in zip(ms_labels, pred)]
    TOL_PPM = 10.0
    R["q4_our_formula_vs_the_labels_DRAWN_in_Fig_S1B"] = {
        "_compares": "our formula's predicted [M+2H]2+ isotope ladder against "
                     "the m/z labels printed inside the spectrum itself. The "
                     "caption plays no part in this check.",
        "labels_read_from": "embedded raster in SI p.4, read by eye at "
                            f"{MS_LABELS_READ_AT_DPI} dpi",
        "printed_labels": ms_labels,
        "label_spacings_Da": sp,
        "mean_spacing_Da": round(mean_sp, 4),
        "charge_state_implied": z,
        "charge_state_stated_in_caption": None,
        "charge_state_note": "z=2 is derived from the 0.5 Da spacing. The "
                             "caption's '2170.3 Da found' is therefore a "
                             "deconvoluted neutral mass; no 2170 ion was "
                             "observed.",
        "predicted_mz": [round(p, 4) for p in pred],
        "ppm_error_per_label": [round(e, 2) for e in ppm_err],
        "worst_ppm": round(max(ppm_err), 2) if ppm_err else None,
        "tol_ppm": TOL_PPM,
        "n_labels": len(ms_labels),
        "agreement_count": sum(1 for e in ppm_err if e <= TOL_PPM),
        "neutral_mass_back_derived_from_label_1": (
            round(ms_labels[0] * z - z * PROTON, 4) if z else None),
        "pass": bool(ppm_err) and max(ppm_err) <= TOL_PPM,
    }

    # -- Q5: the ppm axis, via curation/lib/axis -------------------------
    page = pymupdf.open(PDF)[NMR_PAGE]
    integrals, ticks = nmr_text(page)
    ticks.sort()
    tx = [t[0] for t in ticks]
    tv = [t[1] for t in ticks]
    cal = axis.calibrate(tx, tv, kind="linear", axis="x",
                         loo_tol_px=axis.LOO_TOL_PX)
    loo = axis.loo_residuals(tx)
    R["q5_nmr_ppm_axis"] = {
        "_compares": "nothing yet -- this is the calibration the integral "
                     "assignment in q6 rests on",
        "n_tick_labels_used": len(tx),
        "tick_values": tv,
        "excluded": "the two outermost labels, clipped by the panel frame, so "
                    "their glyph centroids are not on their ticks",
        "caveat": "these are text-layer LABEL centroids, not detected tick "
                  "marks. A label typeset off its tick is error class 20; the "
                  "LOO residuals below are what would expose it.",
        "calibration": cal.qc(),
        "loo_residual_px_worst": round(max(loo), 4),
        "loo_residual_px_all": [round(v, 4) for v in loo],
        "loo_tol_px": axis.LOO_TOL_PX,
        "pass": max(loo) <= axis.LOO_TOL_PX,
    }

    # -- Q6: integral ratios ---------------------------------------------
    # Predicted H counts come from the molecule, counted by rdkit.
    n_arms = sum(1 for a in mol.GetAtoms() if a.GetSymbol() == "S")
    patt = {
        "ester_glycol_OCH2CH2O": ("C(=O)OCCOC(=O)", 4),   # 4 H per arm
        "octyl_beta_to_S_CH2":   ("SCCCCCCCC", 2),        # 2 H per arm
        "octyl_omega_CH3":       ("CCCCCCCCS", 3),        # 3 H per arm
    }
    rows, ref = {}, None
    for name, win in regions.items():
        got = [(xc, v) for xc, v in integrals
               if win["ppm_lo"] <= cal.to_value(xc) <= win["ppm_hi"]]
        _, per_arm = patt[name]
        rows[name] = {
            "ppm_window": [win["ppm_lo"], win["ppm_hi"]],
            "integrals_found": [{"ppm": round(cal.to_value(xc), 3),
                                 "integral": v} for xc, v in got],
            "integral_sum": round(sum(v for _, v in got), 2),
            "predicted_H": per_arm * n_arms,
        }
        if name == "octyl_omega_CH3":
            ref = rows[name]

    TOL_FRAC = 0.05
    checks = []
    for name, r in rows.items():
        if name == "octyl_omega_CH3":
            continue
        obs_ratio = (r["integral_sum"] / ref["integral_sum"]
                     if ref["integral_sum"] else float("nan"))
        pred_ratio = r["predicted_H"] / ref["predicted_H"]
        rel = abs(obs_ratio - pred_ratio) / pred_ratio
        checks.append({"region": name,
                       "observed_ratio_to_omega_CH3": round(obs_ratio, 4),
                       "predicted_ratio": round(pred_ratio, 4),
                       "rel_error": round(rel, 4),
                       "within_tol": rel <= TOL_FRAC})
    R["q6_nmr_integral_ratios"] = {
        "_compares": "integral RATIOS against the corpus structure's own H "
                     "counts. Ratios, not absolute integrals: the 0.87 ppm "
                     "methyl is the spectrum's reference and is set to 18.00 "
                     "by the operator, so comparing it to 18 is circular and "
                     "is NOT counted as a check.",
        "n_arms_from_S_count": n_arms,
        "regions": rows,
        "checks": checks,
        "tol_rel": TOL_FRAC,
        "n_checks": len(checks),
        "agreement_count": sum(1 for c in checks if c["within_tol"]),
        "pass": all(c["within_tol"] for c in checks),
        "what_it_cannot_see": "arm REGIOCHEMISTRY. The ratios fix one glycol "
                             "diester and one octyl thioether per arm; they "
                             "are blind to which core nitrogen each arm is on.",
    }
    return R


# ---------------------------------------------------------------------------
def selftest():
    """Seed each check's defect and require the AGREEMENT COUNT to move.
    A check that cannot be made to fail has not been credited with anything."""
    base = run()
    out = []

    def rec(name, key, field, before, after, msg):
        moved = before != after
        out.append({"check": key, "mutation": name, field + "_clean": before,
                    field + "_seeded": after, "moved": moved, "note": msg})
        return moved

    ok = True

    # Q1: delete one carbon from our structure.
    mut = CORPUS_SMILES.replace("CCCCCCCCSCC(C)C(=O)OCCOC(=O)CC",
                                "CCCCCCCSCC(C)C(=O)OCCOC(=O)CC", 1)
    r = run(smiles=mut)
    ok &= rec("drop one CH2 from one tail", "q1", "agreement_count",
              base["q1_our_formula_vs_SI_printed_string"]["agreement_count"],
              r["q1_our_formula_vs_SI_printed_string"]["agreement_count"],
              "C and H counts must both stop matching")

    # Q2: a one-element typo in the SI's printed formula string.
    r = run(si_formula="C110H203N5O23S6")
    ok &= rec("typo O24 -> O23 in the printed string", "q2", "abs_delta_Da",
              base["q2_SI_string_vs_SI_printed_calc"]["abs_delta_Da"],
              r["q2_SI_string_vs_SI_printed_calc"]["abs_delta_Da"],
              "this is the published-typo case the three-way split exists for")

    # Q4a: shift every drawn label by one CH2/2 -- the ladder must break.
    r = run(ms_labels=[v + 7.0078 for v in MS_LABELS])
    ok &= rec("shift the drawn labels by +7.0078 (one CH2 at z=2)", "q4",
              "agreement_count",
              base["q4_our_formula_vs_the_labels_DRAWN_in_Fig_S1B"]["agreement_count"],
              r["q4_our_formula_vs_the_labels_DRAWN_in_Fig_S1B"]["agreement_count"],
              "5 -> 0")

    # Q4b: break the spacing so the charge state is misread.
    r = run(ms_labels=[MS_LABELS[0] + i * 1.0034 for i in range(5)])
    ok &= rec("1.0 Da spacing (z=1) instead of 0.5", "q4", "agreement_count",
              base["q4_our_formula_vs_the_labels_DRAWN_in_Fig_S1B"]["agreement_count"],
              r["q4_our_formula_vs_the_labels_DRAWN_in_Fig_S1B"]["agreement_count"],
              "charge state is derived, not assumed, so this must move")

    # Q6: point the glycol window at empty spectrum.
    reg = {k: dict(v) for k, v in REGIONS.items()}
    reg["ester_glycol_OCH2CH2O"] = dict(ppm_lo=5.8, ppm_hi=6.4)
    r = run(region_defs=reg)
    ok &= rec("move the glycol window to 5.8-6.4 ppm (empty)", "q6",
              "agreement_count",
              base["q6_nmr_integral_ratios"]["agreement_count"],
              r["q6_nmr_integral_ratios"]["agreement_count"], "2 -> 1")

    # 01_blocks.py: two seeded mutations, run as subprocesses, exit code read
    # directly (never through a pipe).
    blocks = []
    src = open(os.path.join(HERE, "01_blocks.py")).read()
    for label, old, new in [
        ("point the extractor at a page with no category axis",
         "PAGE = 4 ", "PAGE = 2 "),
        ("stop accepting the '50 mol%' level, splitting every axis",
         "(10|20|30|40|50)", "(10|20|30|40)"),
    ]:
        assert old in src, label
        p = os.path.join(HERE, ".seed_blocks.py")
        open(p, "w").write(src.replace(old, new))
        rc = subprocess.run([sys.executable, "-I", p], cwd=HERE,
                            capture_output=True, text=True)
        os.remove(p)
        blocks.append({"mutation": label, "exit_code": rc.returncode,
                       "stderr_first_line": rc.stderr.strip().splitlines()[0]
                       if rc.stderr.strip() else "",
                       "fired": rc.returncode != 0})
        ok &= rc.returncode != 0
    ctrl = subprocess.run([sys.executable, "-I",
                           os.path.join(HERE, "01_blocks.py")], cwd=HERE,
                          capture_output=True, text=True)
    blocks.append({"mutation": "CONTROL (unmutated)",
                   "exit_code": ctrl.returncode,
                   "fired": ctrl.returncode != 0})
    ok &= ctrl.returncode == 0

    # make_reference.py: its whole claim is "the paper draws nothing" and "no
    # enumerator exists to collide with". Both must be falsifiable.
    anchor = []
    asrc = open(os.path.join(HERE, "make_reference.py")).read()

    # (a) a caption that DOES announce drawn chemistry
    old = 'caps = caption_override if caption_override is not None else captions(doc)'
    assert old in asrc
    p = os.path.join(HERE, ".seed_anchor.py")
    open(p, "w").write(asrc.replace(
        old, old + '\n    caps = caps + ["Figure S12. Chemical structure of '
                   '6A1-SC8."] * 1'))
    rc = subprocess.run([sys.executable, "-I", p], cwd=HERE,
                        capture_output=True, text=True)
    os.remove(p)
    anchor.append({"mutation": "inject a caption announcing a drawn structure",
                   "exit_code": rc.returncode, "fired": rc.returncode != 0,
                   "stderr_first_line": rc.stderr.strip().splitlines()[0]
                   if rc.stderr.strip() else ""})
    ok &= rc.returncode != 0

    # (b) an enumerator appears in the folder
    fake = os.path.join(HERE, "03_enumerated.csv")
    created = not os.path.exists(fake)
    if created:
        open(fake, "w").write("name,SMILES\n")
    rc = subprocess.run([sys.executable, "-I",
                         os.path.join(HERE, "make_reference.py")], cwd=HERE,
                        capture_output=True, text=True)
    if created:
        os.remove(fake)
    anchor.append({"mutation": "03_enumerated.csv appears in the folder",
                   "exit_code": rc.returncode, "fired": rc.returncode != 0,
                   "stderr_first_line": rc.stderr.strip().splitlines()[0]
                   if rc.stderr.strip() else ""})
    ok &= rc.returncode != 0

    rc = subprocess.run([sys.executable, "-I",
                         os.path.join(HERE, "make_reference.py")], cwd=HERE,
                        capture_output=True, text=True)
    anchor.append({"mutation": "CONTROL (unmutated)",
                   "exit_code": rc.returncode, "fired": rc.returncode != 0})
    ok &= rc.returncode == 0

    return {"all_seeded_defects_moved_their_agreement_count": bool(ok),
            "value_checks": out, "01_blocks_subprocess": blocks,
            "make_reference_subprocess": anchor}


def main():
    R = run()
    R["_selftest"] = selftest()
    R["_meta"] = {
        "paper": "Sun et al., PNAS 2026, 10.1073/pnas.2528477123, PMID 42150080",
        "rows_emitted_to_the_corpus": 0,
        "why": "see 00_value_sources.json -- there is no per-member value "
               "source in any held file. These checks validate the ONE thing "
               "the SI does supply: an external mass/NMR anchor on 6A1-SC8, "
               "which the corpus already carries from KZ_2016.",
        "emit_rows_called": False,
        "emit_rows_why_not": "curation/lib/emit.emit_rows writes the quality "
                             "columns of corpus rows. There are no corpus rows "
                             "to write, so calling it would mean inventing "
                             "one. Skipped, and said so.",
    }
    fails = {k: v for k, v in R.items()
             if isinstance(v, dict) and v.get("pass") is False}
    R["_gate"] = {
        "checks_run": sum(1 for v in R.values()
                          if isinstance(v, dict) and "pass" in v),
        "failed": sorted(fails),
        "selftest_ok": R["_selftest"]["all_seeded_defects_moved_their_agreement_count"],
    }
    with open(os.path.join(HERE, "07_qc.json"), "w") as fh:
        json.dump(R, fh, indent=1)
    for k, v in R.items():
        if isinstance(v, dict) and "pass" in v:
            print(f"{'PASS' if v['pass'] else 'FAIL'}  {k}")
    print("selftest:", R["_gate"]["selftest_ok"])
    return 0 if (not fails and R["_gate"]["selftest_ok"]) else 1


if __name__ == "__main__":
    sys.exit(main())
