#!/usr/bin/env python3
"""Every check, with numbers, and every check SEEDED -> 07_qc.json.

A check is credited only when the seeded defect MOVES ITS AGREEMENT COUNT.
A sibling curation found two of its seeds were no-ops -- one shadowed by a
default argument, one swapping two identical monomers -- so each seed here
records `clean`, `seeded` and the delta, and `credited` is false unless they
differ. A check whose seed does not move it is reported as NOT CREDITED rather
than quietly listed as passing.

Passing every printed statistic is necessary, not sufficient (QL_2026 did it on
bars 28% too high), so the checks that carry the weight are the UNPRINTED
invariants: 22 monomer masses with zero spread across 108 lipids, eight
formulations summing to exactly 100 mol%, and the same triplicate appearing in
two independent panels.
"""
import csv
import json
import math
import os
import re
import sys

import openpyxl

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
sys.path.insert(0, HERE)
SRC = os.path.join(HERE, "source")

from stage01_blocks import (HEADS, TAILS, exact, hrms_masses,   # noqa: E402
                            monomer_formula, names, MASS_EXCLUDE)


def sheet(book, name):
    wb = openpyxl.load_workbook(os.path.join(SRC, book), data_only=True)
    return list(wb[name].iter_rows(values_only=True))


def rel(a, b):
    return abs(a - b) / abs(b) if b else float("inf")


CHECKS = []


def record(cid, what, clean, seeded, total, detail, seed_desc, unprinted=False):
    CHECKS.append(dict(
        id=cid, what=what, agreements_clean=clean, total=total,
        agreements_with_seed=seeded, delta=clean - seeded,
        credited=bool(clean != seeded), seed=seed_desc,
        unprinted_invariant=unprinted, detail=detail))


# ------------------------------------------------------------------ P1
def chk_hrms(perturb=None):
    """Enumerated exact mass vs the paper's printed `Mass calculated m/z`."""
    rows = list(csv.DictReader(open(os.path.join(HERE, "03_enumerated.csv"))))
    ok, tot, worst = 0, 0, []
    for r in rows:
        if not r["hrms_printed"]:
            continue
        tot += 1
        calc = float(r["exact_mass"])
        if perturb and r["ml_id"] == perturb[0]:
            calc += perturb[1]
        d = calc - float(r["hrms_printed"])
        if abs(d) <= 5e-3:
            ok += 1
        else:
            worst.append((r["ml_id"], r["lipid_name"], round(d, 4)))
    return ok, tot, worst


# ------------------------------------------------------------------ P2
def chk_monomer(bump_head=None):
    """Per-(linker,tail) monomer mass solved from HRMS: spread and deviation."""
    masses, nm = hrms_masses(), names()
    heads = {k: exact(v[1]) for k, v in HEADS.items()}
    if bump_head:
        heads[bump_head] += 14.01565
    from collections import defaultdict
    g = defaultdict(list)
    for mlid, name in nm.items():
        m = re.match(r"H(\d+)-([ON])T(\d+)\((\d)\)$", name)
        if not m or mlid not in masses or mlid in MASS_EXCLUDE:
            continue
        h, lk, tl, n = ("H" + m.group(1), m.group(2), "T" + m.group(3),
                        int(m.group(4)))
        g[(lk, tl)].append((masses[mlid] - heads[h]) / n)
    ok, det = 0, []
    for k, v in sorted(g.items()):
        spread = max(v) - min(v)
        dev = sum(v) / len(v) - exact(monomer_formula(*k))
        good = spread < 1e-3 and abs(dev) < 5e-4
        ok += good
        det.append(dict(group="%s%s" % k, n_lipids=len(v),
                        spread_Da=round(spread, 5), dev_Da=round(dev, 5),
                        ok=good))
    return ok, len(g), det


# ------------------------------------------------------------------ P3
def chk_trace(bump_tail=None):
    """The raster trace's carbon count vs the HRMS formula's, per tail."""
    ok, det = 0, []
    for t, (smi, tf, geom, note) in TAILS.items():
        for lk in ("O", "N"):
            f = monomer_formula(lk, t)
            hrms_c = int(re.search(r"C(\d+)", f).group(1)) - 3
            g = geom + (1 if bump_tail == (lk, t) else 0)
            good = (g == hrms_c)
            ok += good
            det.append(dict(tail="%s%s" % (lk, t), trace_C=g, hrms_C=hrms_c,
                            ok=good))
    return ok, len(det), det


# ------------------------------------------------------------------ P4
def chk_fig3c(corrupt=None):
    """Fig. 3c log10 vs Supplementary Fig. 10 raw RLU, cell by cell."""
    raw = {}
    for r in sheet("41551_2026_1640_MOESM4_ESM.xlsx", "Supplementary Fig. 10"):
        if r[0] is None:
            continue
        raw[str(r[0])] = [r[1], r[2], r[3]]
    # Compare the triplet as a SORTED MULTISET. Fig. 3c lists lipids 107-110
    # with the three replicates ROTATED relative to Supplementary Fig. 10 --
    # same three numbers, different order. A positional comparison calls that
    # 12 mismatches and buries the one real disagreement (lipid 8) among them,
    # so the ordering difference is reported separately instead of scored.
    ok, tot, bad, reordered = 0, 0, [], []
    for r in sheet("41551_2026_1640_MOESM9_ESM.xlsx", "Fig. 3")[2:]:
        if r[0] is None:
            continue
        key = str(int(r[0])) if isinstance(r[0], (int, float)) else str(r[0])
        if key not in raw or any(x is None for x in r[1:4]):
            continue
        a = sorted(10 ** x for x in r[1:4])
        if corrupt and key == corrupt:
            a = sorted([a[0] * 1.3] + a[1:])
        b = sorted(raw[key])
        for x, y in zip(a, b):
            tot += 1
            if rel(x, y) < 2e-4:
                ok += 1
            else:
                bad.append((key, round(x), y))
        if ([round(10 ** x) for x in r[1:4]] != [round(y) for y in raw[key]]
                and all(rel(x, y) < 2e-4 for x, y in zip(a, b))):
            reordered.append(key)
    return ok, tot, dict(value_mismatches=bad, same_values_different_order=reordered)


# ------------------------------------------------------------------ P5
def chk_edfig2c(corrupt=False):
    """ED Fig. 2c's Lipid T3-T12 against the screen, cell by cell."""
    raw = {}
    for r in sheet("41551_2026_1640_MOESM4_ESM.xlsx", "Supplementary Fig. 10"):
        if r[0] is not None:
            raw[str(r[0])] = sorted([r[1], r[2], r[3]])
    MAP = {"Lipid T3": "107", "Lipid T4": "108", "Lipid T5": "109",
           "Lipid T6": "110", "Lipid T7": "1", "Lipid T8": "25",
           "Lipid T9": "33", "Lipid T10": "44", "Lipid T11": "53",
           "Lipid T12": "87"}
    cur, got = None, {}
    for r in sheet("41551_2026_1640_MOESM16_ESM.xlsx", "Extended Data Fig. 2"):
        if r[0]:
            cur = str(r[0]).strip()
        if r[1] is not None and cur:
            got.setdefault(cur, []).append(r[1])
    ok, tot, bad = 0, 0, []
    for nm, mlid in MAP.items():
        a, b = sorted(got[nm]), raw[mlid]
        if corrupt and nm == "Lipid T7":
            a = [a[0] + 1.0] + a[1:]
        for x, y in zip(a, b):
            tot += 1
            # 2e-4: the two sheets print the same number to different
            # significant figures (154595.1 vs 154595.1284).
            if rel(x, y) < 2e-4:
                ok += 1
            else:
                bad.append((nm, mlid, x, y))
    return ok, tot, bad


# ------------------------------------------------------------------ P6
def chk_label(corrupt=False):
    """label.txt (the ML target) against the mean of the three replicates."""
    import zipfile
    z = zipfile.ZipFile(os.path.join(SRC, "41551_2026_1640_MOESM7_ESM.zip"))
    nmf = [n for n in z.namelist() if n.endswith("data/label.txt")][0]
    lab = {}
    for line in z.read(nmf).decode().splitlines():
        if line.strip():
            a, b = line.split()
            lab[a] = float(b)
    raw = {}
    for r in sheet("41551_2026_1640_MOESM4_ESM.xlsx", "Supplementary Fig. 10"):
        if r[0] is not None:
            raw[str(r[0])] = [r[1], r[2], r[3]]
    ok, tot, isrep1, bad = 0, 0, [], []
    for mlid, v in sorted(lab.items(), key=lambda kv: int(kv[0])):
        if mlid not in raw:
            continue
        tot += 1
        reps = raw[mlid]
        if corrupt and mlid == "1":
            reps = [reps[0] * 2, reps[1], reps[2]]
        mean = sum(reps) / 3
        if rel(v, mean) < 2e-4:
            ok += 1
        elif rel(v, reps[0]) < 2e-4:
            isrep1.append(mlid)
        else:
            bad.append((mlid, v, round(mean, 2), reps[0]))
    return ok, tot, dict(equals_replicate_1=isrep1, neither=bad)


# ------------------------------------------------------------------ P7
def chk_shared_triplicates(corrupt=False):
    """The same triplicate appearing in two independent panels."""
    f10 = {}
    for r in sheet("41551_2026_1640_MOESM4_ESM.xlsx", "Supplementary Fig. 10"):
        if r[0] is not None:
            f10[str(r[0])] = [r[1], r[2], r[3]]
    f4 = sheet("41551_2026_1640_MOESM10_ESM.xlsx", "Fig. 4")
    alc_vitro = [f4[i][5] for i in range(2, 5)]
    f3 = sheet("41551_2026_1640_MOESM9_ESM.xlsx", "Fig. 3")
    alc_vivo_3i = [f3[i][22] for i in range(4, 7)]
    alc_vivo_4c = [f4[i][5] for i in range(24, 27)]
    ed = {}
    cur = None
    for r in sheet("41551_2026_1640_MOESM16_ESM.xlsx", "Extended Data Fig. 2"):
        if r[0]:
            cur = str(r[0]).strip()
        if r[1] is not None and cur:
            ed.setdefault(cur, []).append(r[1])
    pairs = [("ALC-0315 in vitro: SuppFig10 vs Fig4c", f10["ALC-0315"],
              alc_vitro),
             ("ALC-0315 in vivo: Fig3i vs Fig4c", alc_vivo_3i, alc_vivo_4c),
             ("MC3 in vitro: SuppFig10 vs EDFig2c", f10["MC3"], ed["MC3"]),
             ("ALC-0315 in vitro: SuppFig10 vs EDFig2c", f10["ALC-0315"],
              ed["ALC-0315"]),
             ("SM-102 in vitro: SuppFig10 vs EDFig2c", f10["SM-102"],
              ed["SM-102"])]
    ok, tot, bad = 0, 0, []
    for lbl, a, b in pairs:
        a = list(a)
        if corrupt and lbl.startswith("MC3"):
            a = [a[0] * 1.01] + a[1:]   # +1 RLU was below tolerance: a no-op
        for x, y in zip(sorted(a), sorted(b)):
            tot += 1
            # 1e-6, not 2e-4: these panels REUSE the same numbers, differing
            # only in printed significant figures, so the tolerance can be
            # tight here and the seed must be a real change.
            if rel(x, y) < 1e-6:
                ok += 1
            else:
                bad.append((lbl, x, y))
    return ok, tot, bad


# ------------------------------------------------------------------ P8
FORMULATIONS = {
    "In vitro screening": [45, 15, 38.5, 1.5],
    "DLin-MC3-DMA": [50, 10, 38.5, 1.5],
    "ALC-0315": [46.3, 9.4, 42.7, 1.6],
    "SM-102": [50, 10, 38.5, 1.5],
    "Lipid P1 coated IgM": [45, 15, 38.5, 1.5],
    "Lipid 11 coated IgM": [45, 15, 38, 1.5, 0.5],
    "SM-102 coated IgM": [25, 20, 52.5, 2, 0.5],
    "Lipid P1/mOVA": [40, 10, 48.5, 1.5],
}


def chk_formulations(corrupt=False):
    ok, bad = 0, []
    for k, v in FORMULATIONS.items():
        v = list(v)
        if corrupt and k == "Lipid P1/mOVA":
            v[0] = 41
        s = round(sum(v), 6)
        if s == 100.0:
            ok += 1
        else:
            bad.append((k, s))
    return ok, len(FORMULATIONS), bad


# ------------------------------------------------------------------ P9
def chk_unique(corrupt=False):
    rows = list(csv.DictReader(open(os.path.join(HERE, "03_enumerated.csv"))))
    keys = [r["inchikey"] for r in rows]
    if corrupt:
        keys[1] = keys[0]
    return len(set(keys)), len(keys), []


# ------------------------------------------------------------------ P10
def corpus_overlap():
    """Which of these lipids are ALREADY in the corpus, and under which source.

    Reported, never acted on. Amine + acrylate Michael chemistry is what the
    Akinc 2008 and Whitehead 2014 lipidoid libraries are made of, so overlap is
    expected; what matters for the corpus is that it is visible before a merge.
    """
    import pandas as pd
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    mine = {r["inchikey"]: r for r in
            csv.DictReader(open(os.path.join(HERE, "03_enumerated.csv")))}
    out = {}
    for pool, f in (("vitro", "data/data_vitro.csv"),
                    ("vivo", "data/data_vivo.csv")):
        path = os.path.join(HERE, "..", "..", f)
        d = pd.read_csv(path, low_memory=False)
        hit = {}
        for smi, exp, nm in zip(d["IL_SMILES"], d["Experiment_ID"],
                                d["IL_name"]):
            if not isinstance(smi, str):
                continue
            m = Chem.MolFromSmiles(smi)
            if m is None:
                continue
            k = Chem.MolToInchiKey(m)
            if k in mine:
                hit.setdefault(k, {}).setdefault(str(exp), set()).add(str(nm))
        out[pool] = dict(
            corpus_rows=len(d), matches=len(hit),
            detail=[dict(ml_id=mine[k]["ml_id"],
                         lipid_name=mine[k]["lipid_name"], inchikey=k,
                         corpus=[dict(experiment=e, il_names=sorted(v))
                                 for e, v in sorted(hit[k].items())])
                    for k in sorted(hit, key=lambda x: int(mine[x]["ml_id"]))])
        out[pool]["il_name_collisions"] = sorted(
            set(d["IL_name"].astype(str)) &
            {r["lipid_name"] for r in mine.values()})
    return out


# ------------------------------------------------------------------ P10
def chk_corpus_names(ov, swap=False):
    """The corpus's OWN names for the 15 shared lipids, against our assignment.

    This is the only EXTERNAL check on the head/tail assignment that does not
    come from this paper at all. The 15 lipids that already sit in the corpus
    arrived there from the Akinc 2008 and Whitehead 2014 lipidoid libraries,
    whose naming is `<amine code>_<linker><chain length>` (e.g. `80_O10`,
    `A113_C14`). Two things then have to hold, and neither needs any outside
    knowledge of what those codes mean:

      (a) the chain length printed in the corpus name equals the carbon count
          we derived for that tail from the trace and the HRMS formula;
      (b) our head and the corpus amine code are in BIJECTION across the 15 --
          one head never maps to two amine codes and vice versa.

    (a) is a third, independent witness on every tail length. (b) is a witness
    on the heads: five of our fourteen heads appear, each consistently.
    """
    tails = {t: TAILS[t][2] for t in TAILS}
    enum = {r["ml_id"]: r for r in
            csv.DictReader(open(os.path.join(HERE, "03_enumerated.csv")))}
    ok, tot, bad, pairs = 0, 0, [], {}
    for h in ov["vitro"]["detail"]:
        e = enum[h["ml_id"]]
        for c in h["corpus"]:
            for nm in c["il_names"]:
                m = re.match(r"A?(\d+)_([OCN])(\d+)$", nm)
                if not m:
                    bad.append((nm, "unparsed"))
                    tot += 1
                    continue
                amine, _, chain = m.group(1), m.group(2), int(m.group(3))
                tot += 1
                tl = e["tail"]
                want = tails[tl] + (1 if swap else 0)
                if chain == want:
                    ok += 1
                else:
                    bad.append((nm, e["lipid_name"], chain, want))
                pairs.setdefault(e["head"], set()).add(amine)
                pairs.setdefault("amine:" + amine, set()).add(e["head"])
    bij = all(len(v) == 1 for v in pairs.values())
    return ok, tot, dict(mismatches=bad, bijection_head_to_amine_code=bij,
                         pairs={k: sorted(v) for k, v in sorted(pairs.items())})


def main():
    c, t, d = chk_hrms()
    s, _, _ = chk_hrms(perturb=("1", 14.01565))
    record("P1", "enumerated exact mass vs the paper's printed HRMS", c, s, t,
           dict(failures=d), "add one CH2 (14.0157 Da) to lipid 1")

    c, t, d = chk_monomer()
    s, _, _ = chk_monomer(bump_head="H4")
    record("P2", "HRMS-solved monomer mass per (linker, tail): spread and "
                 "deviation from the formula", c, s, t, dict(groups=d),
           "add one CH2 to head H4's mass", unprinted=True)

    c, t, d = chk_trace()
    s, _, _ = chk_trace(bump_tail=("O", "T7"))
    record("P3", "raster trace carbon count vs the HRMS formula's, per tail",
           c, s, t, dict(tails=d), "add one carbon to the OT7 trace")

    c, t, d = chk_fig3c()
    s, _, _ = chk_fig3c(corrupt="1")
    record("P4", "Fig. 3c log10 values vs Supplementary Fig. 10 raw RLU",
           c, s, t, d, "scale lipid 1's smallest replicate by 1.3")

    c, t, d = chk_edfig2c()
    s, _, _ = chk_edfig2c(corrupt=True)
    record("P5", "Extended Data Fig. 2c Lipid T3-T12 vs the screen", c, s, t,
           dict(mismatches=d), "add 1.0 RLU to Lipid T7's first replicate",
           unprinted=True)

    c, t, d = chk_label()
    s, _, _ = chk_label(corrupt=True)
    record("P6", "label.txt (ML target) vs the mean of the three replicates",
           c, s, t, d, "double lipid 1's first replicate")

    c, t, d = chk_shared_triplicates()
    s, _, _ = chk_shared_triplicates(corrupt=True)
    record("P7", "the same triplicate re-used across independent panels",
           c, s, t, dict(mismatches=d), "scale MC3's first replicate by 1.01 in Supp Fig 10",
           unprinted=True)

    c, t, d = chk_formulations()
    s, _, _ = chk_formulations(corrupt=True)
    record("P8", "every Supplementary Table 3 formulation sums to 100 mol%",
           c, s, t, dict(failures=d), "change Lipid P1/mOVA's 40 to 41",
           unprinted=True)

    c, t, d = chk_unique()
    s, _, _ = chk_unique(corrupt=True)
    record("P9", "distinct InChIKeys among the enumerated lipids", c, s, t,
           {}, "duplicate one InChIKey")

    ov = corpus_overlap()
    c, t, d = chk_corpus_names(ov)
    sd, _, _ = chk_corpus_names(ov, swap=True)
    record("P10", "the corpus's own AA_2008/KW_2014 names for the 15 shared "
                  "lipids vs our head/tail assignment", c, sd, t, d,
           "add one carbon to every derived tail length", unprinted=True)

    out = dict(
        paper="Lin et al., Nat Biomed Eng 2026, PMID 41851281",
        written_utc="2026-11-10",
        library=dict(ml_ids=121, H_T_library_lipids=116,
                     commercial_benchmarks=3, external_ML_test_lipids=2,
                     distinct_structures_enumerated=116,
                     formulations_in_the_screen=116,
                     note="One formulation per lipid in the screen "
                          "(Supplementary Table 3), so formulations do not "
                          "exceed structures. The separate 18-formulation "
                          "sweep in Supplementary Fig. 26 is ONE lipid (P1) "
                          "and is not part of the screen."),
        checks=CHECKS,
        corpus_overlap=ov,
        credited=sum(1 for x in CHECKS if x["credited"]),
        not_credited=[x["id"] for x in CHECKS if not x["credited"]],
    )
    json.dump(out, open(os.path.join(HERE, "07_qc.json"), "w"), indent=1)
    for x in CHECKS:
        print("%-3s %-62s %4d/%-4d seeded=%-4d credited=%s"
              % (x["id"], x["what"][:62], x["agreements_clean"], x["total"],
                 x["agreements_with_seed"], x["credited"]))


if __name__ == "__main__":
    main()
