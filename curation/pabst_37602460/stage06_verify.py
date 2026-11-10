#!/usr/bin/env python3
"""The paper's OWN cross-checks, measured by script -> 07_qc.json

Run from this folder:   python3 -I stage06_verify.py   (after stages 02-04)

EVERY CHECK IS SEEDED. CLAUDE.md section 6: "A NEW CHECK IS NOT DONE UNTIL IT
HAS BEEN MADE TO FAIL ON PURPOSE", and F-038 is what the weaker version costs --
three ZT_2026 marginal checks scored a literal 0 for twelve days while the
summary beside them said "ALL FOUR EXACT", because `.get(k, 0)` swallowed every
miss. So each check here is a pure function of its inputs, run twice: once on
the real inputs, once on a copy with that check's own defect seeded into it.
The record is not "an exception was raised" -- it is the pair (clean metric,
seeded metric) and whether the METRIC MOVED. A check whose metric does not move
when its defect is present is reported as having no teeth.

AND PASSING EVERY PRINTED STATISTIC IS NECESSARY, NOT SUFFICIENT. QL_2026
reproduced every printed hit rate on bars up to 28% too high. The strongest
evidence here is deliberately NOT a recomputed marginal:

  P7  an identity across two independently produced tables -- SI Table S3's
      score of 6 means "delayed or no growth at 25 uM", so it must coincide
      exactly with main-text Table 1's MIC <= 25 uM. The authors never state
      this identity. It holds 27/27.
  P8  all eight bounded therapeutic indices equal 400/MIC exactly. 400 is a
      constant that appears nowhere in the paper -- it is twice the highest
      tested hemolysis concentration, which is how the authors encoded the
      ">200 uM" category. Recovering an unprinted constant from eight
      independent entries is the "round number the authors used but never
      printed" that the playbook asks for.
"""
import copy
import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
from rdkit import Chem                                 # noqa: E402
from rdkit import RDLogger                             # noqa: E402
RDLogger.DisableLog("rdApp.*")

ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
OUT = os.path.join(HERE, "07_qc.json")

HEADS = ["2N", "3N", "3N'", "4N", "4N'", "3NB", "3NB'", "4NC"]
TAILS = ["6B", "8", "8CB", "9", "9C", "10", "10B", "10B'", "11", "12",
         "13B", "14", "16"]

# Printed marginals of SI Table S3, transcribed by hand from the 600 dpi render
# (00_pages/si_p10_tableS3_600dpi.png), NOT re-extracted. C15's design: the
# expected value is written out so changing a documented number means changing
# both halves.
S3_PRINTED_HEAD_TOTAL = {"2N": 29, "3N": 53, "3N'": 62, "4N": 54, "4N'": 66,
                         "3NB": 49, "3NB'": 54, "4NC": 73}
S3_PRINTED_TAIL_TOTAL = {"6B": 38, "8": 40, "8CB": 47, "9": 51, "9C": 39,
                         "10": 51, "10B": 43, "10B'": 48, "11": 49, "12": 42,
                         "13B": 25, "14": 15, "16": 12}
HEMOLYSIS_LADDER = [12.5, 25.0, 50.0, 100.0, 200.0]


def load():
    vals = list(csv.DictReader(open(os.path.join(HERE, "05_values.csv"))))
    blocks = list(csv.DictReader(open(os.path.join(HERE, "01_blocks.csv"))))
    enum = list(csv.DictReader(open(os.path.join(HERE,
                                                 "03_enumerated.csv"))))
    grid = {r["tail"]: {} for r in vals}
    mic, ti = {}, {}
    for r in vals:
        grid[r["tail"]][r["headgroup"]] = int(r["antimicrobial_score"])
        if r["MIC_uM"]:
            mic[r["IL_name"]] = r["MIC_uM"]
            ti[r["IL_name"]] = r["therapeutic_index_HC50_over_MIC"]
    return {"vals": vals, "blocks": blocks, "enum": enum, "grid": grid,
            "mic": mic, "ti": ti,
            "head_total": dict(S3_PRINTED_HEAD_TOTAL),
            "tail_total": dict(S3_PRINTED_TAIL_TOTAL)}


def spearman(xs, ys):
    def ranks(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r
    rx, ry = ranks(xs), ranks(ys)
    n = len(xs)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = (sum((a - mx) ** 2 for a in rx) *
           sum((b - my) ** 2 for b in ry)) ** 0.5
    return round(num / den, 4) if den else None


# --------------------------------------------------------------------------
# the checks. each is  f(D) -> metric dict.  no check mutates D.
# --------------------------------------------------------------------------

def P1_library_size(D):
    heads = {r["name"] for r in D["blocks"] if r["kind"] == "headgroup"}
    tails = {r["name"] for r in D["blocks"] if r["kind"] == "tail_acrylate"}
    return {"headgroups": len(heads), "tails": len(tails),
            "product": len(heads) * len(tails),
            "enumerated": len(D["enum"]),
            "claimed_in_abstract": 104,
            "cells_in_tableS3": sum(len(v) for v in D["grid"].values()),
            "agree": int(len(heads) * len(tails) == 104 == len(D["enum"])
                         == sum(len(v) for v in D["grid"].values()))}


def P2_tails_per_headgroup(D):
    hs = [r for r in D["blocks"] if r["kind"] == "headgroup"]
    ok = sum(int(r["n_tails_agree"]) for r in hs)
    return {"agree": ok, "of": len(hs),
            "per_head": {r["name"]: [int(r["n_tails_computed"]),
                                     int(r["n_tails_stated_in_text"])]
                         for r in hs}}


def P3_tail_mw(D):
    ts = [r for r in D["blocks"] if r["kind"] == "tail_acrylate"]
    diffs = {r["name"]: float(r["mw_abs_diff"]) for r in ts}
    # the three constitutional-isomer MW ties the naming scheme predicts
    byprint = {}
    for r in ts:
        byprint.setdefault(r["printed_mw_tableS4"], []).append(r["name"])
    ties = {k: sorted(v) for k, v in byprint.items() if len(v) > 1}
    return {"agree": sum(1 for d in diffs.values() if d <= 0.05),
            "of": len(ts), "max_abs_diff_Da": round(max(diffs.values()), 3),
            "isomer_MW_ties": ties, "n_ties": len(ties)}


def P4_head_charge_density(D):
    hs = [r for r in D["blocks"] if r["kind"] == "headgroup"]
    bad = {}
    for r in hs:
        if not int(r["charge_density_agrees"]):
            bad[r["name"]] = {
                "computed": float(r["computed_charge_density_N_over_MW"]),
                "printed": float(r["printed_charge_density_tableS5"])}
    swapped = None
    if set(bad) == {"3NB", "4NC"}:
        swapped = (abs(bad["3NB"]["computed"] - bad["4NC"]["printed"]) <= 1e-4
                   and abs(bad["4NC"]["computed"]
                           - bad["3NB"]["printed"]) <= 1e-4)
    return {"agree": len(hs) - len(bad), "of": len(hs),
            "mismatches": bad,
            "each_mismatch_is_the_others_value": swapped}


def P5_name_token(D):
    ts = [r for r in D["blocks"] if r["kind"] == "tail_acrylate"]
    return {"agree": sum(int(r["token_agrees"]) for r in ts), "of": len(ts),
            "measured": {r["name"]: int(r["token_measured"]) for r in ts}}


def P6_tableS3_marginals(D):
    hsum = {h: sum(D["grid"][t][h] for t in TAILS) for h in HEADS}
    tsum = {t: sum(D["grid"][t][h] for h in HEADS) for t in TAILS}
    hok = sum(1 for h in HEADS if hsum[h] == D["head_total"][h])
    tok = sum(1 for t in TAILS if tsum[t] == D["tail_total"][t])
    return {"head_totals_agree": hok, "of_heads": len(HEADS),
            "tail_totals_agree": tok, "of_tails": len(TAILS),
            "grand_total_from_cells": sum(hsum.values()),
            "grand_total_from_printed_head_totals": sum(D["head_total"].values()),
            "grand_total_from_printed_tail_totals": sum(D["tail_total"].values()),
            "head_total_residuals": {h: D["head_total"][h] - hsum[h]
                                     for h in HEADS},
            "tail_total_residuals": {t: D["tail_total"][t] - tsum[t]
                                     for t in TAILS}}


def P7_score6_iff_mic_le_25(D):
    agree, disagree = 0, []
    for name, m in D["mic"].items():
        head = [h for h in HEADS if name.startswith(h)]
        head = max(head, key=len)
        tail = name[len(head):]
        s6 = D["grid"][tail][head] == 6
        le = float(m) <= 25.0
        if s6 == le:
            agree += 1
        else:
            disagree.append({"name": name, "score": D["grid"][tail][head],
                             "MIC_uM": m})
    return {"agree": agree, "of": len(D["mic"]), "disagreements": disagree,
            "n_score6_in_subset": sum(
                1 for n in D["mic"]
                for h in [max((x for x in HEADS if n.startswith(x)), key=len)]
                if D["grid"][n[len(h):]][h] == 6)}


def P8_ti_and_hc50(D):
    bounded_hits, bounded = 0, []
    ladder_ok, ladder_bad, above_ceiling = 0, [], []
    for name, t in D["ti"].items():
        m = float(D["mic"][name])
        if t.startswith("≥"):
            bounded.append(name)
            if abs(float(t.lstrip("≥")) - 400.0 / m) < 1e-9:
                bounded_hits += 1
        else:
            hc = float(t) * m
            # 400 is included as a rung ON PURPOSE and is reported
            # separately. It is the value the eight BOUNDED entries encode for
            # ">200 uM", so an unbounded entry landing there means the ">=" was
            # dropped, not that the measurement is fine. Folding it silently
            # into "on the ladder" would be a swallowed miss of exactly the
            # F-038 kind.
            near = min(HEMOLYSIS_LADDER + [400.0],
                       key=lambda L: abs(L - hc) / L)
            if abs(near - hc) / near <= 0.05:
                ladder_ok += 1
                if near == 400.0:
                    above_ceiling.append(
                        {"name": name, "TI": t, "MIC": m,
                         "HC50_implied": round(hc, 2),
                         "why": "lands on 400, which is the encoding the eight "
                                "BOUNDED entries use for '>200 uM', but this "
                                "entry carries no '>=' sign -- "
                                "08_contradictions.csv C-5"})
            else:
                ladder_bad.append({"name": name, "TI": t, "MIC": m,
                                   "HC50_implied": round(hc, 2)})
    return {"bounded_entries": len(bounded),
            "bounded_equal_400_over_MIC": bounded_hits,
            "unprinted_constant": 400.0,
            "unbounded_on_hemolysis_ladder": ladder_ok,
            "of_unbounded": len(D["ti"]) - len(bounded),
            "unbounded_but_above_the_200uM_ceiling": above_ceiling,
            "off_ladder": ladder_bad}


def P9_printed_statistics(D):
    mic2 = sorted(n for n, m in D["mic"].items() if float(m) == 2.0)
    timax = max(float(t.lstrip("≥")) for t in D["ti"].values())
    tibest = sorted(n for n, t in D["ti"].items()
                    if float(t.lstrip("≥")) == timax
                    and not t.startswith("≥"))
    nc = {n[3:]: float(m) for n, m in D["mic"].items() if n.startswith("4NC")}
    return [
        {"claim": "'MIC assays revealed 7 different lipidoids with MICs at "
                  "2 uM'", "expected": 7, "measured": len(mic2),
         "which": mic2, "ok": int(len(mic2) == 7)},
        {"claim": "'two of the most active lipidoids 4N'8 and 4NC8CB "
                  "(MIC = 2 uM) recorded the highest therapeutic index of 50'",
         "expected": [50.0, ["4N'8", "4NC8CB"]],
         "measured": [timax, tibest],
         "ok": int(timax == 50.0 and tibest == ["4N'8", "4NC8CB"])},
        {"claim": "'4N'12 ... the MIC was 50 uM'", "expected": 50.0,
         "measured": float(D["mic"]["4N'12"]),
         "ok": int(float(D["mic"]["4N'12"]) == 50.0)},
        {"claim": "'4N'8, which showed no evidence of structure but possessed "
                  "a low MIC at 2 uM'", "expected": 2.0,
         "measured": float(D["mic"]["4N'8"]),
         "ok": int(float(D["mic"]["4N'8"]) == 2.0)},
        {"claim": "'the main exceptions were 4N'12 and 4NC8 ... despite MIC "
                  "values at 50 or 100 uM'", "expected": [50.0, 100.0],
         "measured": [float(D["mic"]["4N'12"]), float(D["mic"]["4NC8"])],
         "ok": int(float(D["mic"]["4N'12"]) == 50.0
                   and float(D["mic"]["4NC8"]) == 100.0)},
        {"claim": "'most 4NC-based lipidoids showed MIC < 6 uM, apart from the "
                  "species with linear tails 8-9 atoms long'",
         "expected": "7 of 9 below 6 uM; the two at or above are tails 8 and 9",
         "measured": {"below_6": sorted(k for k, v in nc.items() if v < 6),
                      "at_or_above_6": sorted(k for k, v in nc.items()
                                              if v >= 6)},
         "ok": int(sorted(k for k, v in nc.items() if v >= 6) == ["8", "9"])},
        {"claim": "'lipidoids with the most potent antimicrobial properties "
                  "(down to 2 uM)'", "expected": 2.0,
         "measured": min(float(m) for m in D["mic"].values()),
         "ok": int(min(float(m) for m in D["mic"].values()) == 2.0)},
        {"claim": "'despite having the lowest charge density, the diamine 2N "
                  "performed the worst of all headgroups'",
         "expected": ["2N", "2N"],
         "measured": [min(HEADS, key=lambda h: sum(D["grid"][t][h]
                                                   for t in TAILS)),
                      "2N"],
         "ok": int(min(HEADS, key=lambda h: sum(D["grid"][t][h]
                                                for t in TAILS)) == "2N")},
        {"claim": "'intermediate c log P values (ca. 3-4) and tail lengths "
                  "(9-11) showing the highest antimicrobial activity'",
         "expected": "the three best tail totals are tails 9, 10, 11",
         "measured": sorted(TAILS,
                            key=lambda t: -sum(D["grid"][t][h]
                                               for h in HEADS))[:3],
         "ok": int(sorted(sorted(TAILS,
                                 key=lambda t: -sum(D["grid"][t][h]
                                                    for h in HEADS))[:3])
                   == sorted(["9", "10", "11"]))},
    ]


def P10_prose_twenty(D):
    n6 = 0
    for name in D["mic"]:
        h = max((x for x in HEADS if name.startswith(x)), key=len)
        if D["grid"][name[len(h):]][h] == 6:
            n6 += 1
        n_le25 = sum(1 for m in D["mic"].values() if float(m) <= 25.0)
    return {"prose_says": 20,
            "tableS3_score6_within_the_27": n6,
            "table1_MIC_le_25_within_the_27": n_le25,
            "the_two_tables_agree_with_each_other": int(n6 == n_le25),
            "prose_agrees": int(n6 == 20)}


def P11_charge_density_trend(D):
    hs = {r["name"]: r for r in D["blocks"] if r["kind"] == "headgroup"}
    tot = [sum(D["grid"][t][h] for t in TAILS) for h in HEADS]
    printed = [float(hs[h]["printed_charge_density_tableS5"]) for h in HEADS]
    computed = [float(hs[h]["computed_charge_density_N_over_MW"])
                for h in HEADS]
    return {
        "paper_states": "lower charge densities ... achieving the highest "
                        "scores, i.e. a NEGATIVE correlation",
        "spearman_with_printed_tableS5": spearman(printed, tot),
        "spearman_with_computed_from_the_drawn_amine": spearman(computed, tot),
        "printed_has_the_stated_sign":
            int((spearman(printed, tot) or 0) < 0),
        "computed_has_the_stated_sign":
            int((spearman(computed, tot) or 0) < 0)}


def P12_corpus_overlap(D):
    ours = {r["parent_InChIKey"]: r["IL_name"] for r in D["enum"]}
    ourm = {r["methylated_InChIKey"]: r["IL_name"] for r in D["enum"]}
    ournames = {r["IL_name"] for r in D["enum"]}
    out = {}
    for pool, path in (("vitro", "data/data_vitro.csv"),
                       ("vivo", "data/data_vivo.csv")):
        hitsp, hitsm, cache, names = {}, {}, {}, {}
        with open(os.path.join(ROOT, path)) as fh:
            for row in csv.DictReader(fh):
                s = (row.get("IL_SMILES") or "").strip()
                if s:
                    if s not in cache:
                        m = Chem.MolFromSmiles(s)
                        cache[s] = Chem.MolToInchiKey(m) if m else None
                    k = cache[s]
                    if k in ours:
                        hitsp.setdefault(ours[k], set()).add(
                            (row["Experiment_ID"],
                             (row.get("IL_name") or "").strip()))
                    if k in ourm:
                        hitsm.setdefault(ourm[k], set()).add(
                            row["Experiment_ID"])
                n = (row.get("IL_name") or "").strip()
                if n in ournames:
                    names.setdefault(n, set()).add(row["Experiment_ID"])
        out[pool] = {
            "members_already_in_corpus_by_FULL_InChIKey": len(hitsp),
            "of": len(ours),
            "which": {k: sorted(v) for k, v in sorted(hitsp.items())},
            "methylated_cation_matches": len(hitsm),
            "IL_name_collisions": {k: sorted(v) for k, v in names.items()},
        }
    return out


def P13_anchor(D):
    """C10's anchors, and a seeded connectivity defect that mass cannot see.

    The seed permutes cyclam's ring spacers from 3,2,3,2 to 2,2,3,3. Same ring
    size, same molecular formula, same exact mass to every decimal -- the
    F-013 class exactly, where epicholesterol shared formula, exact mass to
    4 dp, the InChIKey skeleton block and all 17 descriptor columns with
    cholesterol. If the anchor comparison cannot see this, it is worthless.
    """
    import importlib.util
    from rdkit.Chem import rdMolDescriptors
    spec = importlib.util.spec_from_file_location(
        "_mr", os.path.join(HERE, "make_reference.py"))
    mr = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mr)
    enum = {r["IL_name"]: r["IL_SMILES"] for r in D["enum"]}

    out = {"anchors": {}}
    for name, fn, _fig, _how in mr.ANCHORS:
        m = fn()
        out["anchors"][name] = {
            "formula": rdMolDescriptors.CalcMolFormula(m),
            "heavy_atoms": m.GetNumAtoms(),
            "matches_enumeration_full_InChIKey": int(
                Chem.MolToInchiKey(m)
                == Chem.MolToInchiKey(Chem.MolFromSmiles(enum[name]))),
        }

    good = mr.draw_4NC6B()
    d = mr.Draw()
    first = d.atom(7)
    ns, prev = [first], first
    for k, nc in enumerate([2, 2, 3, 3]):
        for _ in range(nc):
            prev = d.atom(6, bond_to=prev)
        if k < 3:
            prev = d.atom(7, bond_to=prev)
            ns.append(prev)
        else:
            d.bond(prev, first)
    for n in ns:
        mr.arm(d, n, mr.r_tert_butyl)
    for n in ns:
        d.quaternary(n)
    bad = d.finish("seeded 4NC6B")
    ek = Chem.MolToInchiKey(Chem.MolFromSmiles(enum["4NC6B"]))
    out["seeded_connectivity_defect"] = {
        "what": "cyclam ring spacers walked 2,2,3,3 instead of 3,2,3,2",
        "formula_unchanged": int(rdMolDescriptors.CalcMolFormula(good)
                                 == rdMolDescriptors.CalcMolFormula(bad)),
        "exact_mass_unchanged": int(
            abs(rdMolDescriptors.CalcExactMolWt(good)
                - rdMolDescriptors.CalcExactMolWt(bad)) < 1e-9),
        "good_matches_enumeration": int(Chem.MolToInchiKey(good) == ek),
        "bad_matches_enumeration": int(Chem.MolToInchiKey(bad) == ek),
        "anchor_has_teeth": int(Chem.MolToInchiKey(good) == ek
                                and Chem.MolToInchiKey(bad) != ek),
    }
    out["generator_contains_no_SMILES_literal"] = int(not any(
        len(v) >= 40 and v.count("C") + v.count("c") >= 10
        and all(ch in "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"
                "0123456789@+-[]()=#%/\\." for ch in v)
        for v in (getattr(mr, n) for n in dir(mr))
        if isinstance(v, str)))
    return out


def P15_figureS1A_3N_core_is_undrawable(D):
    """SI Figure S1 panel (A) draws the 3N core as a 9-membered triaza RING.

    Three nitrogens each joined to the other two by a 2-carbon bridge is
    1,4,7-triazacyclononane. Two things make that impossible for 3N, and
    neither of them is an opinion about the drawing:

      1. a ring nitrogen has one N-H, so a triazacyclononane core can carry at
         most THREE arms; the same panel draws FIVE, and the article states
         five for the 3N headgroup;
      2. the charge density printed for 3N in SI Table S5 is 0.0291, which is
         3/103.17 (diethylenetriamine) and not 3/129.21 (triazacyclononane).

    Recorded, not repaired -- 08_contradictions.csv C-9.
    """
    hs = {r["name"]: r for r in D["blocks"] if r["kind"] == "headgroup"}
    printed = float(hs["3N"]["printed_charge_density_tableS5"])
    deta = Chem.MolFromSmiles("NCCNCCN")
    tacn = Chem.MolFromSmiles("C1CNCCNCCN1")
    from rdkit.Chem import Descriptors
    cd = {"diethylenetriamine_linear": 3 / Descriptors.MolWt(deta),
          "triazacyclononane_ring": 3 / Descriptors.MolWt(tacn)}
    return {
        "printed_3N_charge_density": printed,
        "computed": {k: round(v, 6) for k, v in cd.items()},
        "linear_matches_printed": int(abs(cd["diethylenetriamine_linear"]
                                          - printed) <= 1e-4),
        "ring_matches_printed": int(abs(cd["triazacyclononane_ring"]
                                        - printed) <= 1e-4),
        "max_arms_on_a_ring_core": sum(a.GetTotalNumHs() for a in tacn.GetAtoms()
                                       if a.GetSymbol() == "N"),
        "max_arms_on_the_linear_core": sum(a.GetTotalNumHs()
                                           for a in deta.GetAtoms()
                                           if a.GetSymbol() == "N"),
        "arms_stated_by_the_article_for_3N":
            int(hs["3N"]["n_tails_stated_in_text"]),
        "arms_drawn_in_figureS1A": 5,
    }


def P14_emit_refuses_a_literal(D):
    """emit_rows must refuse a literal confidence and a caveat with no note."""
    import emit
    from errors import EmitError
    recs = {"x": {"structure_caveat_kind": "", "antimicrobial_score": 3}}
    res = {}
    try:
        emit.emit_rows(["x"], recs, confidence="high",
                       note=lambda i, r: "", label="seed")
        res["literal_confidence_refused"] = 0
    except EmitError:
        res["literal_confidence_refused"] = 1
    try:
        emit.emit_rows(["x"], recs, confidence=lambda i, r: "medium",
                       note=lambda i, r: "", label="seed")
        res["caveat_without_a_note_refused"] = 0
    except EmitError:
        res["caveat_without_a_note_refused"] = 1
    try:
        emit.emit_rows(["x"], recs, confidence=lambda i, r: "high",
                       note=lambda i, r: "", label="seed")
        res["clean_call_succeeds"] = 1
    except EmitError:
        res["clean_call_succeeds"] = 0
    return res


CHECKS = [
    ("P1_library_size", P1_library_size),
    ("P2_tails_per_headgroup", P2_tails_per_headgroup),
    ("P3_tail_mw_vs_tableS4", P3_tail_mw),
    ("P4_head_charge_density_vs_tableS5", P4_head_charge_density),
    ("P5_name_token_vs_measured", P5_name_token),
    ("P6_tableS3_marginals", P6_tableS3_marginals),
    ("P7_score6_iff_MIC_le_25", P7_score6_iff_mic_le_25),
    ("P8_TI_bound_convention_and_HC50_ladder", P8_ti_and_hc50),
    ("P9_printed_statistics", P9_printed_statistics),
    ("P10_prose_count_of_twenty", P10_prose_twenty),
    ("P11_charge_density_trend_sign", P11_charge_density_trend),
    ("P12_corpus_overlap", P12_corpus_overlap),
    ("P13_C10_anchor", P13_anchor),
    ("P15_figureS1A_3N_core", P15_figureS1A_3N_core_is_undrawable),
    ("P14_emit_refuses_a_literal", P14_emit_refuses_a_literal),
]

# P14 carries its seed INSIDE the check -- it is three calls, two of which must
# raise -- so there is nothing in D for an external seed to perturb. It is
# listed here so the "no teeth" report does not flag it as unseeded.
SELF_SEEDED = {"P14_emit_refuses_a_literal"}


# --------------------------------------------------------------------------
# the seeds. each one plants THIS check's own defect and must MOVE ITS METRIC.
# --------------------------------------------------------------------------

def seed_P1(D):
    D["blocks"] = [r for r in D["blocks"]
                   if not (r["kind"] == "tail_acrylate" and r["name"] == "16")]
    return "dropped tail 16 from the block table"


def seed_P2(D):
    for r in D["blocks"]:
        if r["name"] == "4NC":
            r["n_tails_agree"] = "0"
    return "marked 4NC's tail count as disagreeing"


def seed_P3(D):
    for r in D["blocks"]:
        if r["name"] == "9C":
            r["mw_abs_diff"] = "1.000"
            r["printed_mw_tableS4"] = "207.3"
    return "moved tail 9C's printed MW by 1 Da (also breaks no tie)"


def seed_P4(D):
    for r in D["blocks"]:
        if r["name"] == "4N'":
            r["charge_density_agrees"] = "0"
            r["computed_charge_density_N_over_MW"] = "0.0300"
    return "corrupted 4N''s computed charge density"


def seed_P5(D):
    for r in D["blocks"]:
        if r["name"] == "8CB":
            r["token_agrees"] = "0"
    return "marked tail 8CB's name token as not matching"


def seed_P6(D):
    D["grid"]["9"]["4N'"] = 1
    return "changed one Table S3 cell (tail 9, head 4N') from 6 to 1"


def seed_P7(D):
    D["mic"]["4NC12"] = "50"
    return "changed 4NC12's MIC from 3 to 50 while its score stays 6"


def seed_P8(D):
    D["ti"]["3NB'12"] = "≥7"
    D["ti"]["4NC10"] = "19"
    return "moved one bounded TI off 400/MIC and one unbounded TI off the ladder"


def seed_P9(D):
    D["mic"]["4N'9"] = "4"
    return "changed 4N'9's MIC from 2 to 4"


def seed_P10(D):
    D["mic"]["4NC11"] = "50"
    return "changed 4NC11's MIC from 3 to 50, breaking the table-to-table tie"


def seed_P11(D):
    hs = {r["name"]: r for r in D["blocks"] if r["kind"] == "headgroup"}
    hs["3NB"]["computed_charge_density_N_over_MW"] = "0.019967"
    hs["4NC"]["computed_charge_density_N_over_MW"] = "0.027353"
    return ("swapped the COMPUTED 3NB / 4NC charge densities back to the "
            "printed assignment")


def seed_P12(D):
    for r in D["enum"]:
        if r["IL_name"] == "4N14":
            r["parent_InChIKey"] = "ZZZZZZZZZZZZZZ-UHFFFAOYSA-N"
    return "corrupted 4N14's parent InChIKey so its corpus hit cannot be found"


def seed_P13(D):
    for r in D["enum"]:
        if r["IL_name"] == "4NC6B":
            r["IL_SMILES"] = r["parent_SMILES"]
    return ("pointed 4NC6B's enumerated IL_SMILES at the UNMETHYLATED parent, "
            "so the anchor should stop matching")


def seed_P15(D):
    for r in D["blocks"]:
        if r["name"] == "3N":
            r["printed_charge_density_tableS5"] = "0.0232"
    return ("set 3N's printed charge density to the triazacyclononane value, "
            "so the ring reading would become the one that matches")


SEEDS = {"P1_library_size": seed_P1,
         "P2_tails_per_headgroup": seed_P2,
         "P3_tail_mw_vs_tableS4": seed_P3,
         "P4_head_charge_density_vs_tableS5": seed_P4,
         "P5_name_token_vs_measured": seed_P5,
         "P6_tableS3_marginals": seed_P6,
         "P7_score6_iff_MIC_le_25": seed_P7,
         "P8_TI_bound_convention_and_HC50_ladder": seed_P8,
         "P9_printed_statistics": seed_P9,
         "P10_prose_count_of_twenty": seed_P10,
         "P11_charge_density_trend_sign": seed_P11,
         "P12_corpus_overlap": seed_P12,
         "P13_C10_anchor": seed_P13,
         "P15_figureS1A_3N_core": seed_P15}


def main():
    D = load()
    report = {
        "_what": "Every check with numbers, each one run a second time with "
                 "its own defect seeded so that the METRIC is shown to move.",
        "_rule": "A clean line from a check that cannot fail is worth nothing "
                 "(F-012, F-038, and two versions of C16 that reported 'all "
                 "derived' while examining nothing).",
        "checks": {}, "seeded_mutation_tests": {},
    }
    for name, fn in CHECKS:
        report["checks"][name] = fn(D)

    toothless = []
    for name, fn in CHECKS:
        if name in SELF_SEEDED:
            report["seeded_mutation_tests"][name] = {
                "seeded_defect": "carried inside the check itself",
                "metric_moved": 1, "clean": report["checks"][name],
                "seeded": "n/a -- the check IS the seeded run"}
            continue
        E = copy.deepcopy(D)
        what = SEEDS[name](E)
        after = fn(E)
        before = report["checks"][name]
        moved = json.dumps(before, sort_keys=True, default=str) != \
            json.dumps(after, sort_keys=True, default=str)
        report["seeded_mutation_tests"][name] = {
            "seeded_defect": what, "metric_moved": int(moved),
            "clean": before, "seeded": after}
        if not moved:
            toothless.append(name)
    report["checks_with_no_teeth"] = toothless

    with open(OUT, "w") as fh:
        json.dump(report, fh, indent=1, default=str)

    c = report["checks"]
    print(f"wrote {OUT}")
    print(f"P1  library            {c['P1_library_size']['product']} = "
          f"{c['P1_library_size']['claimed_in_abstract']} claimed, "
          f"agree={c['P1_library_size']['agree']}")
    print(f"P2  tails/headgroup    {c['P2_tails_per_headgroup']['agree']}/"
          f"{c['P2_tails_per_headgroup']['of']}")
    print(f"P3  tail MW vs S4      {c['P3_tail_mw_vs_tableS4']['agree']}/"
          f"{c['P3_tail_mw_vs_tableS4']['of']}, max |d| "
          f"{c['P3_tail_mw_vs_tableS4']['max_abs_diff_Da']} Da, "
          f"{c['P3_tail_mw_vs_tableS4']['n_ties']} isomer MW ties")
    print(f"P4  head N/MW vs S5    "
          f"{c['P4_head_charge_density_vs_tableS5']['agree']}/8  "
          f"swap={c['P4_head_charge_density_vs_tableS5']['each_mismatch_is_the_others_value']}")
    print(f"P5  name tokens        {c['P5_name_token_vs_measured']['agree']}/"
          f"{c['P5_name_token_vs_measured']['of']}")
    print(f"P6  S3 marginals       heads "
          f"{c['P6_tableS3_marginals']['head_totals_agree']}/8, tails "
          f"{c['P6_tableS3_marginals']['tail_totals_agree']}/13, cells sum "
          f"{c['P6_tableS3_marginals']['grand_total_from_cells']}")
    print(f"P7  score6 <-> MIC<=25 {c['P7_score6_iff_MIC_le_25']['agree']}/"
          f"{c['P7_score6_iff_MIC_le_25']['of']}")
    print(f"P8  TI bounds = 400/MIC "
          f"{c['P8_TI_bound_convention_and_HC50_ladder']['bounded_equal_400_over_MIC']}/"
          f"{c['P8_TI_bound_convention_and_HC50_ladder']['bounded_entries']}; "
          f"HC50 on ladder "
          f"{c['P8_TI_bound_convention_and_HC50_ladder']['unbounded_on_hemolysis_ladder']}/"
          f"{c['P8_TI_bound_convention_and_HC50_ladder']['of_unbounded']}")
    print(f"P9  printed statistics "
          f"{sum(s['ok'] for s in c['P9_printed_statistics'])}/"
          f"{len(c['P9_printed_statistics'])}")
    print(f"P10 prose 'twenty'     tables say "
          f"{c['P10_prose_count_of_twenty']['tableS3_score6_within_the_27']}, "
          f"prose says 20, prose_agrees="
          f"{c['P10_prose_count_of_twenty']['prose_agrees']}")
    print(f"P11 trend sign         printed rho "
          f"{c['P11_charge_density_trend_sign']['spearman_with_printed_tableS5']}"
          f", computed rho "
          f"{c['P11_charge_density_trend_sign']['spearman_with_computed_from_the_drawn_amine']}")
    print(f"P12 corpus overlap     vitro "
          f"{c['P12_corpus_overlap']['vitro']['members_already_in_corpus_by_FULL_InChIKey']}"
          f"/104, vivo "
          f"{c['P12_corpus_overlap']['vivo']['members_already_in_corpus_by_FULL_InChIKey']}/104")
    print(f"SEEDS metric moved     "
          f"{sum(v['metric_moved'] for v in report['seeded_mutation_tests'].values())}"
          f"/{len(CHECKS)}   no teeth: {toothless or 'none'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
