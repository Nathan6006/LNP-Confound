#!/usr/bin/env python3
"""Independent audit of data_vitro.csv / data_vivo.csv.

Deliberately NOT a rerun of build_dataset.py's own checks - this looks for
classes of error the builder does not know about, so that a clean build and a
clean audit are genuinely independent statements. Prints a finding table; writes
audit_findings.json.

Checks, and why each one is here:

 C1  SMILES wellformedness: parses, canonical, no disconnected fragments. A
     counter-ion or a stray '.' silently changes every descriptor.
 C2  Head-series mass monotonicity. Within one head group, mass must rise with
     tail length. This is how the 4A3-SC6 and KZ_2016 SC3 collisions were found,
     generalized and run over every series in the corpus - it needs no naming
     convention, only that a series exists.
 C3  Duplicate SMILES sharing one experiment, split by whether the colliding
     rows also disagree on the measured value. Same structure + same formulation
     + same cell line but different value is a contradiction; same value is
     merely redundancy.
 C4  Contradictory replicates: identical structure AND formulation AND model
     within an experiment, but different Experiment_value.
 C5  Formulation ratios: present components sum to 100, none negative.
 C6  Value sanity: finite, and |z| within a plausible bound.
 C7  Field consistency: Model/Route agreement, Cargo vs Cargo_type.
 C8  Descriptor agreement with the stored SMILES, recomputed from scratch.
 C9  Stereochemistry preserved where the source had it.
 C12 Z-score normalization unit: one z-score per experiment, not per
     (experiment x cell type). Error class 22, found four times by hand and
     never by a check before this one.
 C10 Enumerated structures vs the structure the PAPER DRAWS (structure_check.py).
     Runs last and is reported here so one command covers the whole Tier-1 suite.
     This is the only check that tests TOPOLOGY against an external ground truth;
     everything above validates composition, which is why three connectivity
     bugs once passed all of C1-C9.
"""
import csv
import json
import math
import os
import re
import statistics as st
import sys
from collections import Counter, defaultdict

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Crippen, Descriptors as D, rdMolDescriptors as rd

RDLogger.DisableLog("rdApp.*")
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")

findings = {}


def load(name):
    return list(csv.DictReader(open(os.path.join(DATA, name))))


def c1_smiles(rows):
    """Wellformedness, INCLUDING radical electrons.

    The radical check exists because a real bug shipped without it: deleting a
    hydroxyl from a SMILES string left a 2-coordinate carbon carrying an
    unpaired electron in 537 rows. It parsed, it canonicalised, and its
    descriptors were self-consistent, so every other check passed. A neutral
    closed-shell lipid should never carry a radical.
    """
    bad_parse, noncanon, frags, charged, radicals = [], [], [], [], []
    seen = {}
    for r in rows:
        s = r["IL_SMILES"]
        if s not in seen:
            m = Chem.MolFromSmiles(s)
            seen[s] = m
        m = seen[s]
        if m is None:
            bad_parse.append(r["LNP_ID"])
            continue
        if Chem.MolToSmiles(m) != s:
            noncanon.append(r["LNP_ID"])
        if "." in s:
            frags.append((r["Experiment_ID"], r["IL_name"]))
        if sum(a.GetNumRadicalElectrons() for a in m.GetAtoms()):
            radicals.append((r["Experiment_ID"], r["IL_name"]))
    return {"unparseable": len(bad_parse), "non_canonical": len(noncanon),
            "disconnected_fragments": len(frags),
            "fragment_examples": sorted(set(frags))[:6],
            "radical_electrons": len(radicals),
            "radical_by_experiment": dict(Counter(e for e, _ in radicals))}


def c2_monotonic(rows):
    """Within a HOMOLOGOUS series, heavier tail must mean heavier molecule.

    Restricted to series whose consecutive mass steps are whole multiples of
    CH2 (14.0157). Without that gate the check is meaningless: in BL_2023 the
    T1/T2/T3 token is a tail IDENTIFIER, not a length, and mass legitimately
    falls as unsaturation rises; in KS_2024 'C24' means two C4 tails, not 24
    carbons. Both produced large violation counts that were pure parser error.
    """
    series = defaultdict(list)
    for r in rows:
        n = r["IL_name"].strip()
        m = re.match(r"^(.*?)[-_]([A-Za-z]{0,2}\d{1,2})$", n)
        if not m:
            continue
        head, tail = m.group(1), m.group(2)
        tm = re.search(r"(\d+)$", tail)
        if not tm:
            continue
        try:
            n_tail = int(float(r["Num_carbon_in_tail"]))
        except (ValueError, TypeError):
            n_tail = None
        series[(r["Experiment_ID"], head, re.sub(r"\d+$", "", tail))].append(
            (int(tm.group(1)), float(r["Molecular.Weight"]), n, n_tail))
    # DOES THE TOKEN MEAN A CHAIN LENGTH? Decided PER SERIES, from the data,
    # rather than from a hand-maintained allowlist of experiments (F-010).
    #
    # This used to be `NAMING_IS_LENGTH = {"KW_2014"}`: one experiment, so 737 of
    # 972 series were skipped and the corpus line read "235 series checked, 0
    # violations" while the single largest curation contributed zero checked
    # series. A check that passes because it did not look is the worst failure
    # mode available, because it reads as a clean report.
    #
    # The test generalises KW_2014's own justification -- 1,271 of 1,273 of its
    # names agreed with the computed longest chain -- and applies it to each
    # series: the token is a length for THIS series only if it equals the
    # measured `Num_carbon_in_tail` for every member. That is self-verifying and
    # needs no list to maintain.
    #
    # It still correctly refuses the cases the allowlist existed for, and refuses
    # them for a stated reason rather than by name: KS_2024 'C24' is two C4
    # tails (token 24, measured 6), KZ_2016 'SC1/SC2/SC3' are tail identifiers
    # (isobutyl/isopentyl/neopentyl), YX_2024 'C1/C2/C3' are component indices.
    # All three pass the CH2 step gate, so that gate alone is NOT sufficient --
    # they contributed 380 spurious violations when it was tried on its own.
    #
    # It is also finer-grained than an experiment ever could be: in LL_2026 the
    # EP series pass (EP6->EP18 steps are exactly 4 CH2 per index unit = 4 arms
    # x 2 carbons) while the Ac8-1/Ac10-1 series correctly do not, because there
    # the trailing 1 is an isomer marker, not a length.
    #
    # Coverage: 235 -> 634 series in vitro across 11 experiments, 12 -> 22 in
    # vivo, with 0 violations.
    CH2 = 14.0157
    viol, homologous = [], 0
    skipped, skipped_steps = Counter(), Counter()
    for key, items in series.items():
        u = {}
        for idx, mw, nm, n_tail in items:
            u[idx] = (mw, nm, n_tail)
        if len(u) < 3:
            continue
        ks = sorted(u)
        steps = [(u[b][0] - u[a][0]) / (b - a) for a, b in zip(ks, ks[1:])]
        good = sum(1 for st_ in steps
                   if st_ > 0 and abs(st_ / CH2 - round(st_ / CH2)) < 0.02
                   and round(st_ / CH2) >= 1)
        if good < 0.7 * len(steps):
            skipped_steps[key[0]] += 1
            continue
        if not all(u[k][2] == k for k in ks):
            skipped[key[0]] += 1       # the token is not this series' chain length
            continue
        homologous += 1
        for (a, b), st_ in zip(zip(ks, ks[1:]), steps):
            if st_ <= 0:
                viol.append({"experiment": key[0], "head": key[1], "token": key[2],
                             "lower": u[a][1], "lower_mw": u[a][0],
                             "higher": u[b][1], "higher_mw": u[b][0]})
    return {"series_examined": len(series),
            "series_checked": homologous,
            "series_skipped_token_not_a_length": dict(skipped),
            "series_skipped_steps_not_whole_CH2": dict(skipped_steps),
            "violations": len(viol), "examples": viol[:8]}


def c3_dup_smiles(rows):
    by = defaultdict(list)
    for r in rows:
        by[(r["Experiment_ID"], r["IL_SMILES"])].append(r)
    contradictory, redundant = [], []
    for (eid, _), rs in by.items():
        names = {r["IL_name"].strip() for r in rs}
        if len(names) < 2:
            continue
        vals = {round(float(r["Experiment_value"]), 6) for r in rs}
        entry = {"experiment": eid, "names": sorted(names), "n_rows": len(rs),
                 "distinct_values": len(vals)}
        (contradictory if len(vals) > 1 else redundant).append(entry)
    return {"colliding_name_groups": len(contradictory) + len(redundant),
            "with_conflicting_values": len(contradictory),
            "with_identical_values": len(redundant),
            "conflicting_examples": contradictory[:6],
            "by_experiment": dict(Counter(e["experiment"]
                                          for e in contradictory + redundant))}


def c4_replicates(rows):
    """True replicates: identical in EVERY descriptive column.

    A partial key is useless here - an earlier version keyed on a handful of
    formulation fields and reported 1001 'contradictions', of which 746 simply
    differed in IL_to_nucleicacid_massratio, a real condition the key omitted.
    """
    cols = [c for c in rows[0] if c not in
            ("Index", "LNP_ID", "Formulation_ID", "Experiment_value")]
    by = defaultdict(list)
    for r in rows:
        by[tuple(r[c] for c in cols)].append(round(float(r["Experiment_value"]), 9))
    rep = {k: v for k, v in by.items() if len(v) > 1}
    exact_dup = {k: v for k, v in rep.items() if len(set(v)) == 1}
    return {"columns_in_key": len(cols),
            "replicate_groups": len(rep),
            "groups_that_are_true_replicates": len(rep) - len(exact_dup),
            "exact_duplicate_rows": sum(len(v) - 1 for v in exact_dup.values()),
            "by_experiment": dict(Counter(k[cols.index("Experiment_ID")]
                                          for k in rep).most_common(8))}


def c5_ratios(rows):
    bad_sum, negative = [], []
    for r in rows:
        tot, neg = 0.0, False
        for c in ("IL_molratio", "HL_molratio", "CHL_molratio", "PEG_molratio",
                  "fifthcomponent_molratio"):
            v = r.get(c, "")
            if v == "":
                continue
            try:
                f = float(v)
            except ValueError:
                continue
            tot += f
            if f < 0:
                neg = True
        if neg:
            negative.append(r["LNP_ID"])
        if tot and abs(tot - 100) > 0.5:
            bad_sum.append((r["Experiment_ID"], round(tot, 2)))
    return {"rows_with_negative_ratio": len(negative),
            "rows_not_summing_to_100": len(bad_sum),
            "offending_experiments": dict(Counter(e for e, _ in bad_sum))}


def c6_values(rows):
    nonfinite, extreme = [], []
    for r in rows:
        try:
            v = float(r["Experiment_value"])
        except ValueError:
            nonfinite.append(r["LNP_ID"])
            continue
        if not math.isfinite(v):
            nonfinite.append(r["LNP_ID"])
        elif abs(v) > 12:
            extreme.append((r["Experiment_ID"], r["IL_name"], round(v, 2)))
    return {"non_finite": len(nonfinite), "abs_z_over_12": len(extreme),
            "examples": extreme[:6]}


def c7_fields(rows, expect_model):
    bad_model = [r["LNP_ID"] for r in rows if r["Model"] != expect_model]
    bad_route = []
    for r in rows:
        if expect_model == "in_vitro" and r["Route_of_administration"] != "in_vitro":
            bad_route.append(r["LNP_ID"])
        if expect_model == "in_vivo" and r["Route_of_administration"] == "in_vitro":
            bad_route.append(r["LNP_ID"])
    cargo = Counter((r["Cargo"], r["Cargo_type"]) for r in rows)
    blank = Counter()
    for c in ("IL_name", "IL_SMILES", "Experiment_method", "Model_type",
              "Cargo", "Experiment_value", "Publication_PMID"):
        blank[c] = sum(1 for r in rows if not str(r.get(c, "")).strip())
    return {"wrong_Model": len(bad_model), "wrong_Route": len(bad_route),
            "cargo_pairs": {f"{a}|{b}": n for (a, b), n in cargo.most_common()},
            "blank_required_fields": {k: v for k, v in blank.items() if v}}


def c8_descriptors(rows):
    """Every published descriptor column, recomputed from the stored SMILES.

    This used to check 7 of the 14 computable columns. The other seven --
    LogP, TPSA, Molar.Refractivity, Fraction.sp3.Carbons, sp3.Carbons,
    Hydrogen.Bond.Acceptors, Aromatic.Rings -- were published and never
    verified, so a stale value in any of them passed. Found by the seeded
    mutation test in test_checks.py: a LogP of 42 on a real row exited 0.

    TOL is per column because the corpus stores them at different precisions.
    """
    TOL = {"ExactMolWt": 1e-6, "LogP": 0.02, "Molar.Refractivity": 0.02,
           "Topological.Polar.Surface.Area": 0.02,
           "Fraction.sp3.Carbons": 1e-4}

    def compute(m):
        return {
            "Molecular.Weight": round(D.MolWt(m), 2),
            "ExactMolWt": D.ExactMolWt(m),
            "Nitrogen.Count": sum(1 for a in m.GetAtoms() if a.GetSymbol() == "N"),
            "Heavy.Atoms": m.GetNumHeavyAtoms(),
            "Rings": rd.CalcNumRings(m),
            "Aromatic.Rings": rd.CalcNumAromaticRings(m),
            "Rotatable.Bonds": rd.CalcNumRotatableBonds(m),
            "Hydrogen.Bond.Donors": D.NumHDonors(m),
            "Hydrogen.Bond.Acceptors": D.NumHAcceptors(m),
            "LogP": Crippen.MolLogP(m),
            "Molar.Refractivity": Crippen.MolMR(m),
            "Topological.Polar.Surface.Area": rd.CalcTPSA(m),
            "Fraction.sp3.Carbons": rd.CalcFractionCSP3(m),
            "sp3.Carbons": sum(1 for a in m.GetAtoms() if a.GetSymbol() == "C"
                               and a.GetHybridization() == Chem.HybridizationType.SP3),
        }

    cache, bad, by_col = {}, 0, Counter()
    for r in rows:
        s = r["IL_SMILES"]
        if s not in cache:
            m = Chem.MolFromSmiles(s)
            cache[s] = None if m is None else compute(m)
        c = cache[s]
        if c is None:
            bad += 1
            continue
        hit = False
        for k, a in c.items():
            v = str(r.get(k, "")).strip()
            if v == "":
                continue                      # not populated is not disagreeing
            if abs(a - float(v)) > TOL.get(k, 0.02):
                by_col[k] += 1
                hit = True
        bad += hit
    return {"rows_checked": len(rows), "unique_lipids": len(cache),
            "rows_disagreeing": bad, "columns_checked": 14,
            "by_column": dict(by_col.most_common())}


def c9_stereo(rows):
    """Double bonds that COULD carry E/Z but do not.

    Uses FindPotentialStereoBonds, so a trisubstituted alkene bearing two
    identical methyls - the terpenoid tails in SL_2020 - is correctly not
    counted: there is no stereochemistry there to lose. C=N is reported
    separately because imine geometry is rarely specified or meaningful.
    """
    cache = {}

    def counts(s):
        if s not in cache:
            m = Chem.MolFromSmiles(s)
            Chem.FindPotentialStereoBonds(m)
            cc = cn = 0
            for b in m.GetBonds():
                if b.GetStereo() != Chem.BondStereo.STEREOANY:
                    continue
                pair = {b.GetBeginAtom().GetSymbol(), b.GetEndAtom().GetSymbol()}
                if pair == {"C"}:
                    cc += 1
                else:
                    cn += 1
            cache[s] = (cc, cn)
        return cache[s]

    cc_rows, cn_rows, mixed = Counter(), Counter(), Counter()
    for r in rows:
        cc, cn = counts(r["IL_SMILES"])
        if cc:
            cc_rows[r["Experiment_ID"]] += 1
            if "/" in r["IL_SMILES"] or "\\" in r["IL_SMILES"]:
                mixed[r["Experiment_ID"]] += 1
        if cn:
            cn_rows[r["Experiment_ID"]] += 1
    return {"rows_with_unassigned_stereogenic_CC": sum(cc_rows.values()),
            "by_experiment": dict(cc_rows.most_common()),
            "rows_mixing_assigned_and_unassigned_CC": dict(mixed),
            "rows_with_unassigned_CN_imine": sum(cn_rows.values()),
            "CN_by_experiment": dict(cn_rows)}


def c12_zscore_unit(rows):
    """ERROR CLASS 22: the z-score was taken per (experiment x cell type).

    The corpus convention is ONE z-score per experiment. When a curation
    standardizes each cell-type block separately instead, every block lands at
    mean 0 and population sd 1, which forces the between-assay offset to
    identically zero -- the quantity a reader would use to ask whether the same
    lipid is more potent in HepG2 than in Ramos. It is not shrunk, it is
    deleted, and the deletion is invisible to every other check because the
    pooled block still sits at mean 0 / sd 1 (that is arithmetic, not evidence:
    pooling blocks each at mean 0 / sd 1 gives mean 0 / sd 1).

    Found four times by hand before this check existed -- ql_2026, WT_2026
    (F-003), YX_2024 and LZ_2025 (F-016) -- and never once by a check, which is
    why it is here. The signature is unmistakable and has no false positives on
    the corpus: a correctly pooled experiment's cell-type blocks have no reason
    to sit at mean 0 / sd 1 individually, and will not.

    Reported against `value_normalization_unit`, so an experiment that DECLARES
    itself per-cell-type is accounted for rather than re-reported forever: the
    violation is a mismatch between what a row claims its normalization unit is
    and what its numbers say it is.
    """
    by_exp = defaultdict(list)
    for r in rows:
        by_exp[r["Experiment_ID"]].append(r)

    flagged, declared, checked = {}, {}, 0
    for eid, block in sorted(by_exp.items()):
        cells = defaultdict(list)
        for r in block:
            try:
                cells[r["Model_type"]].append(float(r["Experiment_value"]))
            except ValueError:
                pass
        usable = {c: v for c, v in cells.items() if len(v) >= 3}
        if len(usable) < 2:
            continue                      # one cell type: nothing to tell apart
        pooled = [v for vs in cells.values() for v in vs]
        if abs(st.mean(pooled)) > 5e-3 or abs(st.pstdev(pooled) - 1) > 0.05:
            continue                      # not on the z-score scale at all
        checked += 1
        # 1e-3 / 5e-3: the corpus stores values rounded, so an exactly
        # standardized block reads as mean ~1e-6 and sd ~1-1e-7, while a pooled
        # experiment's per-cell means are O(0.1-1). Three orders of margin.
        if all(abs(st.mean(v)) < 1e-3 and abs(st.pstdev(v) - 1) < 5e-3
               for v in usable.values()):
            unit = {r.get("value_normalization_unit", "") for r in block}
            flagged[eid] = {"rows": len(block),
                            "blocks": {c: len(v) for c, v in sorted(usable.items())},
                            "declared_unit": sorted(unit)}
            if unit == {"experiment_x_model_type"}:
                declared[eid] = len(block)
    undeclared = {k: v for k, v in flagged.items() if k not in declared}
    return {"experiments_checked": checked,
            "experiments_standardized_per_cell_type": len(flagged),
            "rows_affected": sum(v["rows"] for v in flagged.values()),
            "undeclared_violations": len(undeclared),
            "undeclared_rows": sum(flagged[k]["rows"] for k in undeclared),
            "detail": flagged}


STEROID = Chem.MolFromSmarts(
    "[#6]1~[#6]~[#6]~[#6]2~[#6]1~[#6]~[#6]~[#6]1~[#6]2~[#6]~[#6]~[#6]2~"
    "[#6]~[#6]~[#6]~[#6]~[#6]12")
ESTER_RX = AllChem.ReactionFromSmarts(
    "[C:1](=[O:2])[O:3][C:4]>>[C:1](=[O:2])[OH].[OH][C:4]")


def c14_embedded_sterol(rows):
    """A sterol ESTERIFIED INSIDE the ionizable lipid must still be cholesterol.

    This is the check that would have caught F-013 on day one, and the reason it
    is here rather than in reagent_registry.py is that the registry is blind to
    this by construction, twice over:

      * it compares the InChIKey SKELETON block only, deliberately, "so stereo
        descriptors do not cause noise" -- and 3alpha and 3beta cholesterol
        share their skeleton block exactly;
      * it inspects the CHL_* columns, which are EMPTY on all 537 ZT_2026 rows,
        because that paper's whole claim is cholesterol-free LNPs. The sterol
        lives inside IL_SMILES, where the registry never looks.

    Formula C117H218N5O12+, exact mass 1885.6597 and all 17 descriptor columns
    are bit-identical between the epimers, so nothing else in C1-C13 can see it
    either. The only discriminating measurement is the stereo block of the full
    InChIKey of the liberated alcohol.

    Corpus-wide on purpose, not scoped to ZT_2026: scoping a check to the
    experiment where a bug was found is how the same bug survives in the next
    curation. Pre-filtered on a steroid ring substructure so the reaction
    enumeration only runs on the few hundred molecules that could have one.

    ZERO CLEAVAGE SITES ON A STEROID-BEARING LIPID IS A FAILURE, not a pass.
    A check that silently examines nothing is the failure mode this whole round
    was called to find.
    """
    want = None
    for r in json.load(open(os.path.join(HERE, "reagent_registry.json")))["results"]:
        if r["reagent"] == "Cholesterol":
            want = Chem.MolToInchiKey(Chem.MolFromSmiles(r["registry_smiles"]))
    if want is None:
        return {"error": "no Cholesterol entry in reagent_registry.json",
                "rows_with_embedded_sterol": 0, "rows_wrong_epimer": -1}

    cache, by_exp, wrong, nosite = {}, Counter(), Counter(), Counter()
    for r in rows:
        s = r["IL_SMILES"]
        if s not in cache:
            m = Chem.MolFromSmiles(s)
            if m is None or not m.HasSubstructMatch(STEROID):
                cache[s] = None                 # no steroid: nothing to check
            else:
                keys = set()
                for prod in ESTER_RX.RunReactants((m,)):
                    for x in prod:
                        try:
                            Chem.SanitizeMol(x)
                        except Exception:
                            continue
                        if rd.CalcMolFormula(x) == "C27H46O":
                            keys.add(Chem.MolToInchiKey(x))
                cache[s] = keys
        keys = cache[s]
        if keys is None:
            continue
        by_exp[r["Experiment_ID"]] += 1
        if not keys:
            nosite[r["Experiment_ID"]] += 1     # steroid present, none cleavable
        elif keys != {want}:
            wrong[r["Experiment_ID"]] += 1
    return {"authority": want,
            "rows_with_embedded_sterol": sum(by_exp.values()),
            "by_experiment": dict(by_exp.most_common()),
            "rows_wrong_epimer": sum(wrong.values()),
            "wrong_by_experiment": dict(wrong),
            "rows_steroid_present_but_no_cleavable_ester": sum(nosite.values()),
            "nosite_by_experiment": dict(nosite)}


def c13_quality_vocabulary(rows):
    """Every quality column holds a value schema_quality.json allows.

    A provenance column whose vocabulary is not enforced drifts, and a drifted
    flag is worse than no flag: `med`, `medium` and `Medium` become three
    conditions and a consumer filtering on one of them silently gets a subset.
    `curation_note` and the derived numeric column are free-form by design --
    `allowed: null` in the schema says so -- and are skipped.

    Also counts rows flagged WITHOUT a note. A `structure_assumed` with nothing
    saying which feature was assumed is the pre-2026-10-30 situation in
    miniature: the reader learns that something is uncertain, not what.
    """
    schema = json.load(open(os.path.join(HERE, "schema_quality.json")))
    allowed = {c["name"]: c["allowed"] for c in schema["columns"]}
    bad, unexplained = [], 0
    for r in rows:
        for col, ok in allowed.items():
            if ok is None:
                continue
            v = str(r.get(col, ""))
            if v not in ok:
                bad.append((r["LNP_ID"], col, v))
        if (r.get("structure_assumed") == "1" or r.get("value_is_bound")) \
                and not r.get("curation_note"):
            unexplained += 1
    flagged = sum(1 for r in rows
                  if any(r.get(c) for c in ("confidence", "positional_ambiguity",
                                            "structure_assumed", "value_is_bound")))
    return {"bad_values": len(bad), "examples": bad[:8],
            "rows_carrying_a_flag": flagged,
            "flagged_without_a_note": unexplained}


def check_baseline(findings):
    """Compare every gated metric with audit_baseline.json. Returns failures.

    This is the whole of F-012. Before it, C1-C9 printed and main()'s return
    value was discarded, so the command exited 0 on a corpus with eight seeded
    violations in it. The reason it was print-only is a real one -- the corpus
    has known, deliberately accepted conditions and failing on anything printed
    would fail forever -- so the gate is a baseline rather than a zero.

    A DECREASE fails too, and that is the point: a known condition vanishing
    without anyone editing the baseline is how a check that stopped looking
    announces itself, and that failure mode reads as a clean report.
    """
    base = json.load(open(os.path.join(HERE, "audit_baseline.json")))
    fails = []

    def dig(f, path):
        # "C7_fields.blank_Publication_PMID" indexes blank_required_fields,
        # which OMITS a field entirely when its count is zero -- so a missing
        # key there means 0, not a missing metric.
        chk, _, leaf = path.partition(".")
        node = f.get(chk)
        if node is None:
            return None
        if leaf.startswith("blank_"):
            return node.get("blank_required_fields", {}).get(leaf[6:], 0)
        return node.get(leaf)

    for corpus in ("vitro", "vivo"):
        for path, spec in base[corpus].items():
            if path.startswith("_"):
                continue
            got = dig(findings[corpus], path)
            if got is None:
                fails.append(
                    f"{corpus} {path}: METRIC MISSING from the audit output. "
                    f"The check was renamed or removed, so this baseline is no "
                    f"longer enforcing anything.")
                continue
            want, mode = spec["value"], spec["mode"]
            why = spec["why"][:88]
            if mode == "exact" and got != want:
                what = "NEW" if got > want else "a known condition VANISHED"
                fails.append(f"{corpus} {path}: {got} != accepted {want} "
                             f"({what}) -- {why}")
            elif mode == "max" and got > want:
                fails.append(f"{corpus} {path}: {got} exceeds limit {want} "
                             f"-- {why}")
            elif mode == "min" and got < want:
                fails.append(f"{corpus} {path}: {got} is below the coverage "
                             f"floor {want}; the check is looking at LESS than "
                             f"it was -- {why}")
    return fails


def main():
    vit, viv = load("data_vitro.csv"), load("data_vivo.csv")
    for tag, rows, model in (("vitro", vit, "in_vitro"), ("vivo", viv, "in_vivo")):
        findings[tag] = {
            "rows": len(rows),
            "C1_smiles": c1_smiles(rows),
            "C2_mass_monotonicity": c2_monotonic(rows),
            "C3_duplicate_smiles": c3_dup_smiles(rows),
            "C4_contradictory_replicates": c4_replicates(rows),
            "C5_formulation_ratios": c5_ratios(rows),
            "C6_values": c6_values(rows),
            "C7_fields": c7_fields(rows, model),
            "C8_descriptors": c8_descriptors(rows),
            "C9_stereochemistry": c9_stereo(rows),
            "C12_zscore_unit": c12_zscore_unit(rows),
            "C13_quality_vocabulary": c13_quality_vocabulary(rows),
            "C14_embedded_sterol": c14_embedded_sterol(rows),
        }
    json.dump(findings, open(os.path.join(HERE, "audit_findings.json"), "w"), indent=1)

    for tag in ("vitro", "vivo"):
        f = findings[tag]
        print(f"\n===== {tag}  ({f['rows']} rows) =====")
        print(f"C1 SMILES      unparseable {f['C1_smiles']['unparseable']}  "
              f"non-canonical {f['C1_smiles']['non_canonical']}  "
              f"disconnected {f['C1_smiles']['disconnected_fragments']}  "
              f"RADICALS {f['C1_smiles']['radical_electrons']} "
              f"{f['C1_smiles']['radical_by_experiment'] or ''}")
        c2 = f["C2_mass_monotonicity"]
        print(f"C2 monotonic   {c2['series_checked']} verified-length series checked, "
              f"{c2['violations']} violations "
              f"({sum(c2['series_skipped_token_not_a_length'].values())} series skipped: "
              f"token is not a chain length)")
        c3 = f["C3_duplicate_smiles"]
        print(f"C3 dup SMILES  {c3['colliding_name_groups']} name-groups "
              f"({c3['with_conflicting_values']} with CONFLICTING values, "
              f"{c3['with_identical_values']} identical)")
        c4 = f["C4_contradictory_replicates"]
        print(f"C4 replicates  {c4['groups_that_are_true_replicates']} true replicate groups, "
              f"{c4['exact_duplicate_rows']} exact duplicate rows")
        c5 = f["C5_formulation_ratios"]
        print(f"C5 ratios      not-100 {c5['rows_not_summing_to_100']}  "
              f"negative {c5['rows_with_negative_ratio']}  {c5['offending_experiments']}")
        c6 = f["C6_values"]
        print(f"C6 values      non-finite {c6['non_finite']}  |z|>12 {c6['abs_z_over_12']}")
        c7 = f["C7_fields"]
        print(f"C7 fields      wrong Model {c7['wrong_Model']}  wrong Route {c7['wrong_Route']}  "
              f"blank required {c7['blank_required_fields'] or 'none'}")
        c8 = f["C8_descriptors"]
        print(f"C8 descriptors {c8['rows_disagreeing']} of {c8['rows_checked']} disagree "
              f"({c8['unique_lipids']} unique lipids)")
        c9 = f["C9_stereochemistry"]
        print(f"C9 stereo      {c9['rows_with_unassigned_stereogenic_CC']} rows with an "
              f"unassigned stereogenic C=C {c9['by_experiment']}")
        print(f"               mixing assigned+unassigned in one molecule: "
              f"{c9['rows_mixing_assigned_and_unassigned_CC']}")
        c12 = f["C12_zscore_unit"]
        print(f"C12 z-unit     {c12['undeclared_violations']} of "
              f"{c12['experiments_checked']} multi-cell experiments standardized "
              f"per cell type WITHOUT declaring it "
              f"({c12['undeclared_rows']} rows)")
        for eid, d in c12["detail"].items():
            print(f"               {eid}: {d['rows']} rows, blocks {d['blocks']}, "
                  f"declared {d['declared_unit']}")
        c13 = f["C13_quality_vocabulary"]
        print(f"C13 quality    {c13['rows_carrying_a_flag']} rows carry a flag, "
              f"{c13['bad_values']} values outside the schema vocabulary, "
              f"{c13['flagged_without_a_note']} flagged without a note")
        c14 = f["C14_embedded_sterol"]
        print(f"C14 sterol     {c14['rows_with_embedded_sterol']} rows have a sterol "
              f"esterified inside IL_SMILES; WRONG EPIMER "
              f"{c14['rows_wrong_epimer']} {c14['wrong_by_experiment'] or ''}  "
              f"no cleavable ester {c14['rows_steroid_present_but_no_cleavable_ester']}"
              f" {c14['nosite_by_experiment'] or ''}")

    fails = check_baseline(findings)
    print(f"\n===== baseline gate ({os.path.basename('audit_baseline.json')}) =====")
    if not fails:
        n = sum(1 for c in ("vitro", "vivo")
                for k in json.load(open(os.path.join(HERE, "audit_baseline.json")))[c]
                if not k.startswith("_"))
        print(f"{n} gated metrics, all at their accepted values.")
        return 0
    print(f"{len(fails)} DEPARTURE(S) from the accepted baseline:")
    for x in fails:
        print(f"  * {x}")
    print("\nEach is either a new problem or a check that stopped looking. "
          "Establish which, then update audit_baseline.json in the same commit "
          "with the reason -- do not edit the number to make this green.")
    return 1


def delegate(script, banner):
    """Run a sibling Tier-1 check so the whole suite stays one command."""
    import subprocess
    print(f"\n===== {banner} =====")
    r = subprocess.run([sys.executable, "-I", os.path.join(HERE, script)],
                       capture_output=True, text=True)
    print(r.stdout.rstrip() or r.stderr.rstrip())
    return r.returncode


if __name__ == "__main__":
    # main()'s return value used to be DISCARDED here, so C1-C9 could print
    # anything at all and the command still exited 0 (F-012). It is now the
    # baseline comparison, and it gates.
    rc = main()
    rc |= delegate("structure_check.py", "C10  drawn-structure anchors")
    rc |= delegate("reagent_registry.py", "C11  named commercial reagents")
    print("\nAUDIT PASS" if rc == 0 else "\nAUDIT FAIL")
    sys.exit(rc)
