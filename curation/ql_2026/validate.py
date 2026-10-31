#!/usr/bin/env python3
"""Independent check of the QL_2026 curation. Re-reads the OUTPUT files and the
figures' own printed summaries; does not import the builders.
"""
import csv
import json
import math
import os
import random
import statistics as st

from rdkit import Chem, RDLogger
from rdkit.Chem import Crippen, Descriptors as D, rdMolDescriptors as rd

RDLogger.DisableLog("rdApp.*")
HERE = os.path.dirname(os.path.abspath(__file__))

fails, notes = [], []


def check(ok, msg):
    (notes if ok else fails).append(("PASS" if ok else "FAIL") + "  " + msg)


# Derived, not a literal: 67 curation columns + 5 recomputed descriptors + the
# quality columns in provenance/schema_quality.json. A frozen count here is the
# error class F-002 was logged for -- it fired on a correct file the moment the
# schema legitimately grew, which trains whoever hits it to edit the number.
NCOL = 72 + len(json.load(open(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "..",
    "provenance", "schema_quality.json")))["columns"])


def main():
    rows = list(csv.DictReader(open(os.path.join(HERE, "11_QL_2026_rows.csv"))))
    enum = list(csv.DictReader(open(os.path.join(HERE, "09_enumerated_344.csv"))))
    heat = {r["lipid_name"]: r for r in
            csv.DictReader(open(os.path.join(HERE, "01_fig2b_320.csv")))}
    b2c = list(csv.DictReader(open(os.path.join(HERE, "07_fig2c_bars.csv"))))
    b3c = list(csv.DictReader(open(os.path.join(HERE, "05_fig3c_24.csv"))))
    bars = json.load(open(os.path.join(HERE, "04_fig2def_qc.json")))

    check(len(rows) == 344, f"row count {len(rows)} (expect 344 = 320 + 24)")
    check(len(rows[0]) == NCOL, f"column count {len(rows[0])} (expect {NCOL})")
    check(len(heat) == 320, f"Fig 2b cells {len(heat)}")
    check(sum(1 for r in b3c if r["is_benchmark"] == "0") == 24,
          "Fig 3c library-2 bars = 24")

    # z-scores
    z = [float(r["Experiment_value"]) for r in rows]
    check(abs(st.mean(z)) < 1e-9 and abs(st.pstdev(z) - 1) < 1e-9,
          f"pooled z mean {st.mean(z):.1e} sd {st.pstdev(z):.12f}")
    lib = {r["lipid_name"]: r["library"] for r in enum}
    check(abs(st.mean([float(r["Experiment_value"]) for r in rows
                       if lib[r["IL_name"]] == "1"])) > 0.05,
          "libraries are NOT standardised separately (that would force both "
          "block means to 0 and delete the between-library difference)")

    # composition
    s = [sum(float(r[f"{p}_molratio"]) for p in ("IL", "HL", "CHL", "PEG")) for r in rows]
    check(max(abs(x - 100) for x in s) < 1e-9, "molar percents sum to 100")

    # structures
    smis = [r["IL_SMILES"] for r in rows]
    bad = []
    for sm in set(smis):
        m = Chem.MolFromSmiles(sm)
        if m is None:
            bad.append("unparseable")
        elif Chem.MolToSmiles(m) != sm:
            bad.append("non-canonical")
        elif len(Chem.GetMolFrags(m)) != 1:
            bad.append("disconnected")
        elif sum(a.GetNumRadicalElectrons() for a in m.GetAtoms()):
            bad.append("radical")
    check(not bad, f"{len(set(smis))} unique structures parse, canonical, "
                   f"single-fragment, radical-free {set(bad) if bad else ''}")
    check(len(set(smis)) == 344, f"all 344 structures distinct ({len(set(smis))})")

    # one arm per primary amine, re-derived from the product
    arms = {"70": 1, "75": 1, "80": 1, "400": 1, "113": 2, "114": 2,
            "110": 3, "T3A": 3}
    wrong = 0
    for e in enum:
        n_amide = len(Chem.MolFromSmiles(e["IL_SMILES"]).GetSubstructMatches(
            Chem.MolFromSmarts("[NX3][CX3]=[OX1]")))
        want = arms[e["amine"]]
        # each arm adds one secondary amide plus two imide C=O in the ring
        if n_amide != want * 3 and e["bl"] != "L2":
            wrong += 1
        if e["bl"] == "L2" and n_amide != want:
            wrong += 1
    check(wrong == 0, f"amide count matches one arm per primary amine "
                      f"({wrong} mismatches)")

    # descriptors
    mism = 0
    for r in random.Random(7).sample(rows, 120):
        m = Chem.MolFromSmiles(r["IL_SMILES"])
        for k, v in (("Heavy.Atoms", m.GetNumHeavyAtoms()),
                     ("Nitrogen.Count", sum(1 for a in m.GetAtoms()
                                            if a.GetSymbol() == "N")),
                     ("LogP", round(Crippen.MolLogP(m), 6)),
                     ("ExactMolWt", round(D.ExactMolWt(m), 6))):
            got = float(r[k])
            if abs(got - float(v)) > 1e-5:
                mism += 1
    check(mism == 0, f"descriptors agree with their SMILES on 120 sampled rows ({mism})")

    # figure cross-checks
    plotted = [r["lipid_name"] for r in b2c if r["is_benchmark"] == "0"]
    agree = sum(1 for n in plotted if float(heat[n]["RLU"]) > 2e4)
    check(agree == 27 and len(plotted) == 28,
          f"Fig 2c set vs Fig 2b heatmap: {agree}/{len(plotted)} agree "
          f"(the exception, T3A-L8-C14, is flagged)")

    ratios = []
    for n in plotted:
        if heat[n]["saturated"] == "0" and n != "T3A-L8-C14":
            ratios.append(float([r for r in b2c if r["lipid_name"] == n][0]["RLU"])
                          / float(heat[n]["RLU"]))
    check(0.95 < st.median(ratios) < 1.10,
          f"Fig 2b vs Fig 2c agree on shared lipids, median ratio "
          f"{st.median(ratios):.3f} over n={len(ratios)}")

    check(bars["max_snap_error_pct"] < 0.6,
          f"Fig 2d/2e/2f bars land on exact k/40 and k/64, max error "
          f"{bars['max_snap_error_pct']}%")
    check(bars["implied_total_hits_from_2d"] == bars["implied_total_hits_from_2e"],
          f"printed 2d and 2e agree on the total "
          f"({bars['implied_total_hits_from_2d']})")

    # Fig 3c must reproduce THREE printed statistics. An earlier version
    # asserted only 2C10 and 2C12 and excused 2C11 as "a boundary bar"; 2C11 was
    # exactly the datum that discriminates, and excusing it let a miscalibrated
    # y-axis through. All three are now required.
    lib3 = [r for r in b3c if r["is_benchmark"] == "0"]
    hits3 = [r["lipid_name"] for r in lib3 if float(r["RLU"]) > 7e5]
    by_tail = {t: sum(1 for h in hits3 if h.endswith(t))
               for t in ("2C10", "2C11", "2C12")}
    check(by_tail == {"2C10": 3, "2C11": 1, "2C12": 1},
          f"Fig 3d printed hit counts reproduced exactly: {by_tail}")
    by_bl = {f"L{i}": sum(1 for h in hits3 if f"-L{i}-" in h) for i in range(1, 9)}
    want_bl = {"L1": 1, "L2": 0, "L3": 0, "L4": 0, "L5": 2, "L6": 1, "L7": 1, "L8": 0}
    check(by_bl == want_bl, f"Fig 3e printed per-linker rates reproduced "
                            f"(incl. L8 = 0): {by_bl}")
    n1e5 = sum(1 for r in lib3 if float(r["RLU"]) > 1e5)
    check(n1e5 == 17, f"Results sentence '17 lipids ... greater than 1e5' "
                      f"reproduced: {n1e5}")

    bm2 = {r["lipid_name"]: float(r["RLU"]) for r in b2c if r["is_benchmark"] == "1"}
    bm3 = {r["lipid_name"]: float(r["RLU"]) for r in b3c if r["is_benchmark"] == "1"}
    rr = [bm3[k] / bm2[k] for k in bm2]
    check(st.pstdev(rr) < 0.1 and 1.0 < st.mean(rr) < 1.3,
          f"the 4 shared benchmarks give a consistent scalar offset "
          f"{st.mean(rr):.3f}x (sd {st.pstdev(rr):.3f}), applied to library 2")

    # the biology must survive the scaling: library 2 is ~40x library 1
    import math as _m
    gm = {L: 10 ** st.mean([_m.log10(float(r["RLU"])) for r in
                            ({"1": [x for x in heat.values()]}.get(L) or [])])
          for L in ()}
    z1 = st.mean([float(r["Experiment_value"]) for r in rows
                  if lib[r["IL_name"]] == "1"])
    z2 = st.mean([float(r["Experiment_value"]) for r in rows
                  if lib[r["IL_name"]] == "2"])
    check(z2 - z1 > 1.5,
          f"six-tail library 2 still sits far above library 1 after z-scoring "
          f"(mean z {z2:.2f} vs {z1:.2f}) - the paper's central finding is "
          f"preserved, not standardised away")

    # The conflict row must carry the heatmap value (7,931), not Fig 2c's
    # 58,094. Checked by comparing its z against lipids of known value rather
    # than against the median, which is dragged up by 113 floored cells.
    zof = {r["IL_name"]: float(r["Experiment_value"]) for r in rows}
    zc = zof["T3A-L8-C14"]
    peers_low = [zof[n] for n, h in heat.items()
                 if 7000 < float(h["RLU"]) < 9000 and n in zof]
    peers_high = [zof[n] for n, h in heat.items()
                  if h["saturated"] == "1" and n not in
                  {r["lipid_name"] for r in b2c}]
    check(abs(zc - st.mean(peers_low)) < 0.05 and
          (not peers_high or zc < min(peers_high) - 0.5),
          f"T3A-L8-C14 carries the HEATMAP value as documented: z={zc:.3f} "
          f"matches lipids at ~8e3 RLU (mean z {st.mean(peers_low):.3f}), not "
          f"the saturated band (z>={min(peers_high):.3f})")

    out = {"failures": fails, "passes": notes}
    json.dump(out, open(os.path.join(HERE, "13_validation.json"), "w"), indent=1)
    for line in notes + fails:
        print(" ", line)
    print(f"\n{len(fails)} failures, {len(notes)} passes")
    return 1 if fails else 0


if __name__ == "__main__":
    raise SystemExit(main())
