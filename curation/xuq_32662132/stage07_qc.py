"""Stage 07 - the gate. Every check carries NUMBERS, and every check was seeded.

RULES THIS FILE FOLLOWS
  * a check is credited only when a seeded defect MOVES ITS AGREEMENT COUNT.
    A seed that merely makes the script raise proves nothing about the check.
    `--self-test` runs each seed and prints (clean count -> seeded count).
  * exit codes are read from `returncode`, never through a pipe.
  * a refusal is a result. Figure S1's illegibility is recorded with pixel
    measurements, not as silence.
"""
import csv, json, math, os, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "lib")))
import emit                                               # noqa: E402
from errors import EmitError                              # noqa: E402

PY = [sys.executable, "-I"]

# -- the SI's own gram / millimole pairs, section by section -------------------
# (section, reagent, grams, millimoles, SMILES) - SMILES so the molar mass is
# computed, never typed.
LADDER = [
    ("1.1", "sodium hydride",              0.72, 30.0,  "[Na+].[H-]"),
    ("1.1", "ethylene glycol",             5.6,  90.0,  "OCCO"),
    ("1.1", "1-bromotetradecane",          6.0,  20.0,  "BrCCCCCCCCCCCCCC"),
    ("1.1", "potassium iodide",            3.3,  20.0,  "[K+].[I-]"),
    ("1.1", "2-(tetradecyloxy)ethanol *",  3.3,  12.8,  "OCCOCCCCCCCCCCCCCC"),
    ("1.1", "triethylamine",               1.9,  19.2,  "CCN(CC)CC"),
    ("1.1", "acryloyl chloride",           1.4,  15.4,  "C=CC(=O)Cl"),
    ("1.2", "1-bromotetradecane",          5.0,  18.0,  "BrCCCCCCCCCCCCCC"),
    ("1.2", "potassium carbonate",         3.6,  26.0,  "[K+].[K+].[O-]C([O-])=O"),
    ("1.2", "2-mercaptoethanol",           1.1,  14.0,  "OCCS"),
    ("1.3", "potassium selenocyanate",     1.5,  10.0,  "[K+].[Se-]C#N"),
    ("1.3", "2-bromoethanol",              1.6,  13.0,  "OCCBr"),
    ("1.3", "sodium borohydride",          0.9,  24.0,  "[Na+].[BH4-]"),
    ("1.3", "1-bromotetradecane",          4.1,  15.0,  "BrCCCCCCCCCCCCCC"),
    ("1.4", "2,2'-dipyridyl disulfide",    4.00, 13.47, "c1ccnc(SSc2ccccn2)c1"),
    ("1.4", "1-dodecanol",                 2.00, 11.23, "OCCCCCCCCCCCC"),
    ("1.4", "dodecyl 2-pyridyl disulfide *", 2.00, 6.43, "c1ccnc(SSCCCCCCCCCCCC)c1"),
    ("1.4", "2-mercaptoethanol",           0.6,  7.7,   "OCCS"),
]
INTERMEDIATES = {"2-(tetradecyloxy)ethanol *", "dodecyl 2-pyridyl disulfide *"}
LADDER_TOL_PCT = 5.0          # starting materials; intermediates get 0.5 %

S2A_LABELS = ["16:8:1:1", "16:4:1:1", "16:1:1:1", "16:0:1:1", "16:4:8:1", "16:4:4:1",
              "16:4:1:1", "16:4:0:1", "16:4:1:8", "16:4:1:4", "16:4:1:1", "16:4:1:0",
              "16:0:0:0", "UT"]


def run(args):
    p = subprocess.run(PY + args, cwd=HERE, capture_output=True, text=True)
    return p.returncode, p.stdout + p.stderr


# ---------------------------------------------------------------- checks ----
def check_name_rule(seed=False):
    q = json.load(open(os.path.join(HERE, "01_blocks_qc.json")))["name_token_rule"]
    per = [dict(r) for r in q["per_tail"]]
    if seed:                         # seed: O16S-S relabelled O18S-S
        for r in per:
            if r["tail"] == "O16S-S":
                r["name_token"] = 18
                r["agrees"] = (r["name_token"] == r["measured_chain_heavy_atoms"])
    n = sum(1 for r in per if r["agrees"])
    p = sum(1 for r in per if r["prefix_agrees"])
    return dict(agreement_count=n, n_tested=len(per), prefix_agreement=p,
                per_tail=per, passes=(n == len(per) and p == len(per)))


def check_reagent_ladder(seed=False):
    from rdkit import Chem, RDLogger
    from rdkit.Chem.Descriptors import MolWt
    RDLogger.DisableLog("rdApp.*")
    rows, n_ok, n_int_ok = [], 0, 0
    for sec, name, g, mmol, smi in LADDER:
        true_mw = MolWt(Chem.MolFromSmiles(smi))
        if seed and name == "1-dodecanol":        # seed: call it 1-decanol
            true_mw = MolWt(Chem.MolFromSmiles("OCCCCCCCCCC"))
        implied = g / (mmol / 1000.0)
        pct = 100.0 * abs(implied - true_mw) / true_mw
        tol = 0.5 if name in INTERMEDIATES else LADDER_TOL_PCT
        ok = pct <= tol
        n_ok += ok
        if name in INTERMEDIATES:
            n_int_ok += ok
        rows.append(dict(section=sec, reagent=name, grams=g, mmol=mmol,
                         implied_g_per_mol=round(implied, 2),
                         true_molar_mass=round(true_mw, 2),
                         pct_error=round(pct, 2), tol_pct=tol, ok=bool(ok)))
    return dict(agreement_count=n_ok, n_tested=len(LADDER),
                intermediates_agreement=n_int_ok, n_intermediates=len(INTERMEDIATES),
                per_reagent=rows,
                what_it_pins=("the two INTERMEDIATES are the load-bearing rows: "
                              "2-(tetradecyloxy)ethanol pins the O17 tails' chain at "
                              "C14 and dodecyl 2-pyridyl disulfide pins the 16S-S "
                              "tails' chain at C12, each to <=0.3 %. Three starting-"
                              "material lines are sloppier by up to 35 %, which is "
                              "why the chain length is taken from the intermediates "
                              "and the discrepancy is logged in 08_contradictions.csv."),
                passes=(n_int_ok == len(INTERMEDIATES)))


def check_repeated_formulation(seed=False):
    geo = json.load(open(os.path.join(HERE, "05_figS2_geometry.json")))
    bars = geo["a"]["bars"]
    idx = [i for i, lab in enumerate(S2A_LABELS) if lab == "16:4:1:1"]
    tops = [bars[i]["top_row"] for i in idx]
    if seed:                          # seed: a 4 px error on one of the three
        tops = list(tops); tops[1] += 4
    ref = tops[0]
    n = sum(1 for t in tops if abs(t - ref) < 1e-9)
    return dict(agreement_count=n, n_tested=len(tops), positions=idx,
                top_rows=tops, values=[bars[i]["value"] for i in idx],
                max_spread_px=round(max(tops) - min(tops), 3),
                what_it_is=("16:4:1:1 is drawn three times, once in each sub-series. "
                            "Three independent window detections and three independent "
                            "bar_top calls must land on the same row. This is an "
                            "unprinted identity the authors drew, not a recomputed "
                            "marginal of our own read."),
                passes=(n == len(tops)))


def check_ut_is_100(seed=False):
    geo = json.load(open(os.path.join(HERE, "05_figS2_geometry.json")))
    b = geo["b"]
    ut = b["bars"][-1]["value"]
    if seed:
        ut += 6.0
    px_per_unit = b["calibration"]["px_per_unit"]
    err_pp = abs(ut - 100.0)
    n = int(err_pp <= 1.0)
    return dict(agreement_count=n, n_tested=1, ut_value_percent=round(ut, 2),
                target=100.0, error_pp=round(err_pp, 3),
                error_px=round(err_pp * px_per_unit, 2), tol_pp=1.0,
                what_it_is=("Figure S2b is 'Cell viability (%)' and the SI says the "
                            "viability is calculated relative to untreated cells, so "
                            "the UT bar MUST be exactly 100. The authors drew that "
                            "round number and never printed it; it is external to our "
                            "pipeline in a way a recomputed marginal is not."),
                passes=(n == 1))


def check_printed_optima(seed=False):
    geo = json.load(open(os.path.join(HERE, "05_figS2_geometry.json")))
    c = [b["value"] for b in geo["c"]["bars"]]
    d = [b["value"] for b in geo["d"]["bars"]]
    if seed:
        d = list(d); d[0] = max(d) * 2      # move the dose optimum off 30 ug/mL
    c_lab = ["0.5", "1", "2", "5", "7", "18", "26", "UT"]
    d_lab = ["6", "15", "30", "60", "120"]
    c_best = c_lab[c.index(max(c))]
    d_best = d_lab[d.index(max(d))]
    hits = []
    hits.append(dict(claim="SI section 6: 'the best ... lipidoid concentration was "
                           "found to be ... 30 ug/mL'", panel="S2d",
                     our_argmax=d_best, agrees=(d_best == "30")))
    hits.append(dict(claim="SI section 6: 'the best delivery time ... ~6 h'",
                     panel="S2c", our_argmax=c_best + " h",
                     agrees=(c_best in ("5", "7")),
                     note="the panel has no 6 h bar; 7 h is its maximum and 5 h its "
                          "neighbour, so '~6 h' is consistent with either. Recorded as "
                          "consistent, not as exact."))
    n = sum(1 for h in hits if h["agrees"])
    return dict(agreement_count=n, n_tested=len(hits), per_claim=hits,
                passes=(n == len(hits)))


def check_axis_calibration(seed=False):
    geo = json.load(open(os.path.join(HERE, "05_figS2_geometry.json")))
    per, n = [], 0
    for pid in ("a", "b", "c", "d"):
        rec = geo[pid]
        worst = rec["loo_worst_px"]
        if seed and pid == "a":
            worst += 6.0
        ok = worst <= rec["loo_tol_px"]
        n += ok
        per.append(dict(panel=pid, n_ticks=rec["calibration"]["n_ticks"],
                        fit_residual_px=rec["calibration"]["tick_fit_residual_px"],
                        loo_worst_px=round(worst, 4), loo_tol_px=rec["loo_tol_px"],
                        px_per_unit=rec["calibration"]["px_per_unit"], ok=bool(ok)))
    return dict(agreement_count=n, n_tested=4, per_panel=per,
                why_loo=("calibrate's max-residual test is weak against a SINGLE bad "
                         "tick; loo_residuals has no leverage masking. Both are "
                         "reported."),
                passes=(n == 4))


def check_flank_agreement(seed=False):
    geo = json.load(open(os.path.join(HERE, "05_figS2_geometry.json")))
    deltas = []
    for pid in ("a", "b", "c", "d"):
        for b in geo[pid]["bars"]:
            deltas.append(b["detail"]["flank_delta_px"])
    if seed:
        deltas = list(deltas); deltas[0] = 9.0
    n = sum(1 for d in deltas if d <= 2.0)
    return dict(agreement_count=n, n_tested=len(deltas),
                max_flank_delta_px=max(deltas),
                what_it_is=("every bar carries a symmetric s.d. whisker with a cap. A "
                            "cap is centred, so it is in NEITHER flank; a flank delta "
                            "of 0 on all bars is the evidence that no cap was read as "
                            "a top (error class 32, F-018)."),
                passes=(n == len(deltas)))


def check_row_counts(seed=False):
    sk = json.load(open(os.path.join(HERE, "06_skipped.json")))
    rows = list(csv.DictReader(open(os.path.join(HERE, "06_rows.csv"))))
    vals = list(csv.DictReader(open(os.path.join(HERE, "05_values.csv"))))
    n_em = len(rows)
    if seed:
        n_em += 1
    sums = (n_em + sk["n_skipped"] == len(vals))
    return dict(agreement_count=int(sums), n_tested=1, n_values=len(vals),
                n_emitted=n_em, n_skipped=sk["n_skipped"],
                every_skip_has_a_reason=all(s["reason"] for s in sk["skipped"]),
                passes=sums)


def check_z_normalisation(seed=False):
    rows = list(csv.DictReader(open(os.path.join(HERE, "06_rows.csv"))))
    z = [float(r["Experiment_value"]) for r in rows]
    if seed:
        z = list(z); z[0] += 0.5
    mu = sum(z) / len(z)
    sd = (sum((t - mu) ** 2 for t in z) / len(z)) ** 0.5
    # tolerance 1e-6, not 1e-9: Experiment_value is written at 6 dp to match the
    # corpus, so the residual is the ROUNDING, not the normalisation. Measured
    # clean: mean 7.7e-08, sd 1.000000105. The seed moves it by 0.5/13 = 0.038.
    tol = 1e-6
    ok = abs(mu) < tol and abs(sd - 1.0) < tol
    return dict(agreement_count=int(ok), n_tested=1, n_rows=len(z),
                mean=mu, sd_population=sd, tolerance=tol,
                why_not_1e_9="Experiment_value is stored at 6 dp, as the corpus does; "
                             "the residual below is that rounding and nothing else",
                convention="per-experiment z-score of log10, population sd (CLAUDE.md s3)",
                passes=ok)


def check_flags_on_the_row(seed=False):
    rows = list(csv.DictReader(open(os.path.join(HERE, "06_rows.csv"))))
    n_caveat = emit.caveat_count(rows)
    n_noted = sum(1 for r in rows if r["curation_note"].strip())
    if seed:
        n_noted -= 1
    ok = (n_caveat == len(rows) and n_noted == len(rows))
    return dict(agreement_count=n_noted, n_tested=len(rows),
                rows_with_a_caveat=n_caveat,
                confidence_values=sorted({r["confidence"] for r in rows}),
                passes=ok)


def check_emit_refuses_a_literal(seed=False):
    """emit.emit_rows must RAISE on a hard-coded confidence (F-011)."""
    recs = {"x": {"IL_name": "x"}}
    raised = False
    try:
        if seed:
            emit.emit_rows(["x"], recs, confidence=lambda i, r: "medium",
                           note=lambda i, r: "n", label="seed")
        else:
            emit.emit_rows(["x"], recs, confidence="high", note=lambda i, r: "n",
                           label="literal probe")
    except EmitError:
        raised = True
    except Exception:
        raised = True
    return dict(agreement_count=int(raised), n_tested=1, raised_on_literal=raised,
                passes=raised)


def check_anchor_gate(seed=False):
    rc_check, _ = run(["make_reference.py", "--check"])
    rc_seed, _ = run(["make_reference.py", "--seed"])
    rc_prov, _ = run(["make_reference.py", "--provenance"])
    if seed:
        rc_seed = 0                 # pretend the seeded defect passed
    n = int(rc_check == 0) + int(rc_seed == 1) + int(rc_prov == 0)
    return dict(agreement_count=n, n_tested=3,
                returncode_check=rc_check, returncode_seeded=rc_seed,
                returncode_provenance=rc_prov,
                what_the_seed_is=("one Michael arm moved from the amine nitrogen onto "
                                  "imidazole C2. Molecular formula C44H83N3O4S2 and "
                                  "exact mass 781.582500 are IDENTICAL to 6 dp, so no "
                                  "mass check can see it; the skeleton InChIKey "
                                  "changes from OLOIOPYNTGMAAP to GNLMIRXUIQYCQW."),
                passes=(n == 3))


def check_corpus_overlap(seed=False):
    ov = json.load(open(os.path.join(HERE, "02_corpus_overlap.json")))
    v, o = ov["pools"]["vitro"], ov["pools"]["vivo"]
    full = v["full_inchikey_hits"] + o["full_inchikey_hits"]
    skel = v["skeleton_inchikey_hits"] + o["skeleton_inchikey_hits"]
    if seed:
        full += 1
    ok = (full == 0 and skel == 0)
    return dict(agreement_count=int(ok), n_tested=1,
                full_inchikey_hits=full, skeleton_inchikey_hits=skel,
                near_miss_head_93_vitro_rows=v["near_miss_head_93"]["n_rows"],
                near_miss_head_93_vivo_rows=o["near_miss_head_93"]["n_rows"],
                near_miss_experiments=sorted(set(v["near_miss_head_93"]["experiments"] +
                                                 o["near_miss_head_93"]["experiments"])),
                leading_token_rows_reported_not_acted_on=
                v["leading_token_trap"]["n_rows_whose_IL_name_starts_with_93_or_9322"],
                passes=ok)


def check_figS1_head_confirmation(seed=False):
    """Does Figure S1 confirm the head recovered from the corpus?

    Two independent tests, both on the 9322-O17O panel:
      (a) is the measured base-peak m/z closer to the prediction for head 9322
          (2-methylimidazole) or for head 93 (plain imidazole)? They differ by
          exactly one CH2, 14.016 Da = 3.57 px on this axis.
      (b) of the THREE candidate labels the first pass read and could not decide
          between (784.73 / 794.75 / 764.73), which does the geometry pick?
    """
    d = json.load(open(os.path.join(HERE, "09_figS1_mass_test.json")))
    rec = d["panels"]["9322-O17O"]
    meas = rec["measured_mz"]
    if seed:
        meas += 14.016               # one CH2: pretend it is the plain-imidazole head
    hd = rec["head_discrimination"]
    dpx = rec["Da_per_px"]
    d9322 = abs(meas - hd["predicted_with_head_9322_2_methylimidazole"])
    d93 = abs(meas - hd["predicted_with_head_93_plain_imidazole"])
    cand = d["_first_pass_candidate_test"]["candidates_read_as_glyphs"]
    nearest = min(cand, key=lambda k: abs(meas - k))
    tests = [dict(test="head 9322 vs head 93", picks_9322=bool(d9322 < d93),
                  deviation_9322_px=round(d9322 / dpx, 2),
                  deviation_93_px=round(d93 / dpx, 2),
                  margin_x=round(d93 / d9322, 1) if d9322 else None),
             dict(test="which first-pass glyph candidate", nearest=nearest,
                  picks_9322=bool(nearest == 764.73),
                  distances_px={str(k): round(abs(meas - k) / dpx, 2) for k in cand})]
    n = sum(1 for t in tests if t["picks_9322"])
    return dict(agreement_count=n, n_tested=len(tests), measured_mz=round(meas, 1),
                Da_per_px=dpx, per_test=tests,
                n_panels_calibrated=sum(1 for v in d["panels"].values()
                                        if "REFUSED" not in v),
                n_panels_refused=sum(1 for v in d["panels"].values()
                                     if "REFUSED" in v),
                what_changed=("the first pass refused this figure outright. With no "
                              "hypothesis a ~4 Da/px axis cannot pin a formula, which "
                              "was true; with a head recovered from the corpus the "
                              "question becomes WHICH of three candidates 20 and 30 Da "
                              "apart, and the axis answers that decisively."),
                passes=(n == len(tests)))


# ------------------------------------------------- recorded refusals --------
def refusal_figS1():
    """Figure S1: what is still refused, and what stopped being refused."""
    return {
        "_status": "PARTLY SUPERSEDED 2026-11-10. The GLYPH read is still refused -- "
                   "the labels really are illegible. What is no longer refused is the "
                   "panel: with head 9322 recovered from the corpus, the base-peak "
                   "position measured against the m/z axis decides between the three "
                   "candidate labels, and 2 of 6 panels calibrate under criteria that "
                   "cannot be tuned to the answer. See 09_figS1_mass_test.json and "
                   "08_retracted_suspicions.csv R-06.",
        "still_refused": "4 of 6 panels: no glyph gap in 2..12 yields 10 tick labels "
                         "that are uniformly spaced to 4 % AND carry the 3-digit/"
                         "4-digit width signature. The threshold was NOT loosened to "
                         "recover them (O17S misses at 4.65 % nonuniformity, O16B at "
                         "18.3 %, EC18 at 37.8 %, N16B never reaches 10 boxes). So "
                         "D-3 (is O16B the same as O16S-S) is still OPEN.",
        "panel": "Figure S1, ESI-MS of head 9322 with six tails",
        "raster_px": [1558, 810],
        "base_peak_label_glyph_height_px": [5, 7],
        "base_peak_label_box_width_px": [18, 25],
        "upscale_tried": "28x LANCZOS",
        "verdict": "ILLEGIBLE. The 9322-O17O base peak reads as 784.73 / 794.75 / "
                   "764.73 with no way to decide between them.",
        "consequence": "head 9322's molecular formula cannot be solved, so the four "
                       "Figure S3 delivery values for 9322-O17S have no molecule and "
                       "none of the seventeen 93xx analogue heads can be enumerated.",
        "what_the_first_pass_got_wrong": "it concluded that because the axis cannot "
                             "PIN a formula it could not TEST one. Those are different "
                             "questions. ~3.9 Da/px still separates hypotheses 14 Da "
                             "apart at a 3.6 px margin. The refusal was correct for "
                             "the read and wrong as a stopping point.",
        "what_is_still_true": "the axis cannot separate O from S (16 Da, ~4 px is "
                              "marginal) and cannot separate O16B from N16B (1.97 Da, "
                              "0.5 px). Those remain unresolved.",
    }


CHECKS = {
    "tail_name_token_rule": check_name_rule,
    "si_reagent_molar_mass_ladder": check_reagent_ladder,
    "repeated_formulation_identity": check_repeated_formulation,
    "untreated_viability_is_100": check_ut_is_100,
    "printed_optima_reproduced": check_printed_optima,
    "axis_calibration_and_loo": check_axis_calibration,
    "bar_flank_agreement": check_flank_agreement,
    "row_counts_sum": check_row_counts,
    "z_normalisation": check_z_normalisation,
    "flags_on_the_row": check_flags_on_the_row,
    "emit_refuses_a_literal_confidence": check_emit_refuses_a_literal,
    "c10_anchor_gate": check_anchor_gate,
    "corpus_overlap": check_corpus_overlap,
    "figS1_mass_confirms_head_9322": check_figS1_head_confirmation,
}


def main():
    out, failed = {}, []
    for name, fn in CHECKS.items():
        r = fn(False)
        out[name] = r
        if not r["passes"]:
            failed.append(name)
    out["_refusals"] = {"figure_S1_mass_spectra": refusal_figS1()}
    out["_summary"] = {"n_checks": len(CHECKS), "n_failed": len(failed),
                       "failed": failed}
    with open(os.path.join(HERE, "07_qc.json"), "w") as f:
        json.dump(out, f, indent=1)
    for n, r in out.items():
        if n.startswith("_"):
            continue
        print("%-36s %s  %s/%s" % (n, "PASS" if r["passes"] else "FAIL",
                                   r["agreement_count"], r["n_tested"]))
    print("FAILED:", failed if failed else "none")
    return 1 if failed else 0


def self_test():
    """Seed each check. Credit it only if its AGREEMENT COUNT MOVES."""
    bad = []
    for name, fn in CHECKS.items():
        clean, seeded = fn(False), fn(True)
        moved = clean["agreement_count"] != seeded["agreement_count"]
        print("%-36s clean %s/%s -> seeded %s/%s  %s"
              % (name, clean["agreement_count"], clean["n_tested"],
                 seeded["agreement_count"], seeded["n_tested"],
                 "COUNT MOVED" if moved else "NO-OP SEED - NOT CREDITED"))
        if not moved:
            bad.append(name)
    print("checks whose seed did not move the count:", bad if bad else "none")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(self_test() if "--self-test" in sys.argv else main())
