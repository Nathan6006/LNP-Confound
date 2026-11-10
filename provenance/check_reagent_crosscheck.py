#!/usr/bin/env python3
"""C17 — one reagent name, one structure, across every experiment.

    python3 -I provenance/check_reagent_crosscheck.py     # exits non-zero on a
                                                          # gated contradiction

WHY THIS EXISTS. Error class 46 -- "a reagent whose structure already exists,
independently curated, elsewhere in the corpus, and was never consulted" -- was
logged on 2026-11-09 and never implemented.
The corpus is its own best external anchor: a SECOND INDEPENDENT CURATION OF
THE SAME REAGENT is the cheapest C10 there is, and it is already paid for.

It would have caught F-031 for free. `ZT_2026` read `OF-2` off a figure by eye
with all four tails one carbon short -- C17 where the reagent is C18 -- while
`LL_2025` held the same reagent correctly as `OF-02`. Two curations, one
reagent, two structures, and nothing compared them.

WHAT IT MUST NOT DO is report the 288 generic combinatorial labels as errors.
`A1_B1_C1` is in both `JC_2023` and `YX_2024` and they are unrelated libraries
that share a labelling scheme. `CLAUDE.md` §3 is explicit that this is not an
error and that one must never key on `IL_name` alone across experiments. So:

  GATED     a name in the C11 registry, or carried in a non-IL component column
            (HL / CHL / PEG / fifthcomponent). Excipients are catalogue
            reagents, so a collision there is a real contradiction.
  REPORTED  a collision in `IL_name` only, where the generic labels live.
            Printed with its structures for a human; never fails the build.

The baseline is `reagent_crosscheck_baseline.json`, same posture as
`audit_baseline.json`: each accepted contradiction is recorded WITH THE REASON,
and a departure in EITHER direction fails. A known contradiction vanishing
without anyone editing the baseline is how a check that stopped looking
announces itself.

AND THE RULE THAT GOVERNS A FIX. Decide WHICH HALF is wrong before changing
anything. In `AC_2025` the NAMES were wrong and the structures were right --
the COMET source stores SMILES and labels with no names at all, so the values
belonged to the structures, and repairing it the other way would have moved 419
measured values onto the wrong molecules.
"""
import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
LIB = os.path.join(HERE, "..", "curation", "lib")
sys.path.insert(0, LIB)

import reagents                                            # noqa: E402

DATA = os.path.join(HERE, "..", "data")
BASELINE = os.path.join(HERE, "reagent_crosscheck_baseline.json")
OUT = os.path.join(HERE, "reagent_crosscheck.json")


def load_rows():
    rows = []
    for name in ("data_vitro.csv", "data_vivo.csv"):
        with open(os.path.join(DATA, name)) as fh:
            rows += list(csv.DictReader(fh))
    return rows


def main():
    rows = load_rows()
    # Skeleton InChIKey: connectivity. A stereo comparison is run too, because
    # C11 compared skeleton blocks only and was blind to epicholesterol twice
    # (F-013) -- but it is REPORTED, not gated, until a baseline exists for it.
    skel = reagents.cross_check(rows, stereo=False)
    ster = reagents.cross_check(rows, stereo=True)

    print(f"C17  reagent name -> structure, across experiments")
    print(f"     {len(rows)} rows, {skel['n_names']} (name, column) pairs, "
          f"{skel['n_multi_experiment']} carried by >1 experiment")
    print(f"     GATED contradictions:    {len(skel['gated'])}")
    print(f"     reported (generic IL):   {len(skel['reported'])}")
    stereo_only = len(ster["gated"]) - len(skel["gated"])
    print(f"     stereo-only differences: {stereo_only} "
          "(reported, not gated -- no baseline yet)")

    baseline = {}
    if os.path.exists(BASELINE):
        with open(BASELINE) as fh:
            baseline = json.load(fh).get("accepted", {})

    failures = []
    for e in skel["gated"]:
        key = f"{e['prefix']}_name:{e['name']}"
        print()
        print(reagents.explain(e))
        if key in baseline:
            print(f"  ACCEPTED in the baseline: {baseline[key]['why']}")
        else:
            failures.append(key)
            print("  *** NOT IN THE BASELINE -- this fails the check. ***")

    vanished = [k for k in baseline
                if k not in {f"{e['prefix']}_name:{e['name']}"
                             for e in skel["gated"]}]
    for k in vanished:
        print()
        print(f"  *** {k} is in the baseline but no longer contradicts. ***")
        print("      A known condition vanishing without anyone editing the "
              "baseline is how a check that stopped looking announces itself. "
              "If it was repaired, remove it from the baseline IN THE SAME "
              "COMMIT as the repair and say so in CHANGELOG.md.")

    with open(OUT, "w") as fh:
        json.dump({"what": "C17 cross-experiment reagent name -> structure",
                   "gated": skel["gated"], "reported": skel["reported"],
                   "stereo_gated": ster["gated"],
                   "n_names": skel["n_names"],
                   "n_multi_experiment": skel["n_multi_experiment"],
                   "unbaselined": failures, "vanished": vanished}, fh,
                  indent=1)
    print(f"\nwrote {os.path.relpath(OUT, os.path.join(HERE, '..'))}")

    if failures or vanished:
        print(f"\nC17 FAILED: {len(failures)} unbaselined contradiction(s), "
              f"{len(vanished)} vanished.")
        return 1
    print("\nC17 PASSED.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
