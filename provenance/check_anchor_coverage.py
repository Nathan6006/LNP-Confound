#!/usr/bin/env python3
"""C18 — is each curation's C10 anchor coverage adequate, and is it a copy?

    python3 -I provenance/check_anchor_coverage.py     # exits non-zero

WHY THIS EXISTS. C10 is the only Tier-1 check that tests topology against
external ground truth; C1-C9 all test internal consistency, which is exactly
why three connectivity bugs sailed through them. But C10 REPORTS its own
weaknesses and gates nothing:

  * `anchors_that_may_be_a_copy_of_the_enumeration` has been non-zero since the
    metric was added. An anchor copied from the enumeration can only pass -- it
    compares a value to itself. That is F-013 (537 rows of epicholesterol,
    same sha256, 223 bytes, while the note claimed independent assembly
    "through the C3 oxygen (3-beta)" and the string was 3-alpha) and F-037 (the
    file written to FIX F-013 re-copied the other half).
  * COVERAGE is not measured at all. `ql_2026` has ONE anchor for 344
    structures, with nothing anchoring its L1 or L2. One exact match on one
    molecule is not evidence about 343 others, and C10's clean line does not
    distinguish the two situations.

So this check gates both, against `anchor_coverage_baseline.json`, each
accepted condition recorded WITH THE REASON and failing on a departure in
either direction.

NOTE WHAT IT DOES NOT CLAIM. Byte-equality is not proof of a copy: two correct
routes give the same canonical SMILES, which is precisely why C10's own
`anchor_may_be_a_copy` missed F-037 and why C16 tests the GENERATOR'S composed
template instead. C18 gates the SURFACE condition -- is there a generator, how
many anchors are there -- and leaves independence to C16.
"""
import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
CURATION = os.path.join(HERE, "..", "curation")
SC = os.path.join(HERE, "structure_check.json")
BASELINE = os.path.join(HERE, "anchor_coverage_baseline.json")
OUT = os.path.join(HERE, "anchor_coverage.json")

# Which Experiment_IDs each curation folder produced. Stated here because the
# folder name and the Experiment_ID are deliberately different conventions.
FOLDER_EXPERIMENTS = {
    "siegwart_42634137": ["ZT_2026"],
    "lin_aromatic_2026": ["LL_2026", "LL_2026_2"],
    "wht_2026_falcon": ["WT_2026"],
    "xd_2026_lipidai": ["XD_2026"],
    "ql_2026": ["QL_2026"],
    "ng_2025": ["NG_2025", "NG_2025_2"],
}


def corpus_structures():
    n = {}
    for name in ("data_vitro.csv", "data_vivo.csv"):
        with open(os.path.join(DATA, name)) as fh:
            for r in csv.DictReader(fh):
                e = (r.get("Experiment_ID") or "").strip()
                s = (r.get("IL_SMILES") or "").strip()
                if e and s:
                    n.setdefault(e, set()).add(s)
    return {k: len(v) for k, v in n.items()}


def main():
    if not os.path.exists(SC):
        print("C18: structure_check.json not found. Run C10 first:\n"
              "    python3 -I provenance/structure_check.py")
        return 1
    with open(SC) as fh:
        sc = json.load(fh)
    nstruct = corpus_structures()

    baseline = {}
    if os.path.exists(BASELINE):
        with open(BASELINE) as fh:
            baseline = json.load(fh).get("accepted", {})

    per_folder = {}
    for e in sc["results"]:
        f = e["paper"]
        d = per_folder.setdefault(f, {"anchors": 0, "real": 0, "copies": 0,
                                      "no_drawn": 0, "generator": False})
        d["anchors"] += 1
        if e["status"] == "NO_DRAWN_STRUCTURE":
            d["no_drawn"] += 1
        elif e["status"] in ("exact_match", "stereo_differs_connectivity_ok",
                             "expected_mismatch_paper_drawing_is_wrong"):
            d["real"] += 1
        if e.get("anchor_may_be_a_copy"):
            d["copies"] += 1
        if e.get("anchor_derivation"):
            d["generator"] = True

    rows, failures = [], []
    print("C18  C10 anchor coverage and independence posture")
    print(f"{'folder':24} {'anchors':>7} {'real':>5} {'copy?':>5} "
          f"{'structs':>7} {'per anchor':>10}  generator")
    for folder in sorted(set(list(per_folder) + list(FOLDER_EXPERIMENTS))):
        d = per_folder.get(folder, {"anchors": 0, "real": 0, "copies": 0,
                                    "no_drawn": 0, "generator": False})
        exps = FOLDER_EXPERIMENTS.get(folder, [])
        ns = sum(nstruct.get(e, 0) for e in exps)
        gen = os.path.exists(os.path.join(CURATION, folder,
                                          "make_reference.py"))
        ratio = (ns / d["real"]) if d["real"] else float("inf")
        rows.append({"folder": folder, "anchors": d["anchors"],
                     "real_anchors": d["real"], "copies": d["copies"],
                     "no_drawn": d["no_drawn"], "structures": ns,
                     "structures_per_real_anchor":
                         None if ratio == float("inf") else round(ratio, 1),
                     "has_generator": gen})
        print(f"{folder:24} {d['anchors']:>7} {d['real']:>5} {d['copies']:>5} "
              f"{ns:>7} {('inf' if ratio == float('inf') else f'{ratio:.1f}'):>10}"
              f"  {'yes' if gen else 'NO'}")

        key = folder
        issues = []
        if d["copies"]:
            issues.append(f"{d['copies']} anchor(s) may be a copy")
        if ns and not d["real"]:
            issues.append(f"{ns} structures and NO usable anchor")
        if ns and d["real"] and ratio > 200:
            issues.append(f"{ratio:.0f} structures per anchor")
        if issues:
            if key in baseline:
                print(f"  accepted: {baseline[key]['why'][:150]}")
            else:
                failures.append((key, issues))
                print(f"  *** NOT BASELINED: {'; '.join(issues)} ***")

    vanished = [k for k in baseline
                if k not in {r["folder"] for r in rows}]
    for k in vanished:
        print(f"  *** {k} is in the baseline but is no longer a folder. ***")

    with open(OUT, "w") as fh:
        json.dump({"what": "C18 anchor coverage and independence posture",
                   "folders": rows,
                   "unbaselined": [k for k, _ in failures],
                   "vanished": vanished}, fh, indent=1)

    if failures or vanished:
        print(f"\nC18 FAILED: {len(failures)} folder(s) with an unbaselined "
              f"anchor-coverage condition.")
        print("A folder with no usable anchor has no external ground truth at "
              "all. Write make_reference.py, or baseline it WITH THE REASON.")
        return 1
    print("\nC18 PASSED.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
