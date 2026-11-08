#!/usr/bin/env python3
"""C15 — every number the docs ASSERT, measured against the live corpus.

WHY THIS EXISTS. On 2026-11-07 a sweep of the numbers in `CLAUDE.md` found
three that were wrong and nothing that had noticed:

  * the F-001 flag count was stated as 1,152, the docs told the reader to
    re-measure it with a metric that returns 1,652, and the correct figure for
    what the prose describes is 1,157 (F-024);
  * the chance baseline was `1/31 vitro, 1/17 vivo`, which matched neither the
    live corpus (S = 42 / 19) nor the corpus the frozen pipeline ran on
    (S = 34 / 17) -- it was the S from before a July re-trim (F-026);
  * a trailing-space condition was documented as affecting "several" rows when
    it was 915 cells (F-022).

C1-C14 all check the data against itself or against a paper. None of them check
the data against our own prose, and the prose is what the next session reads
first and what the manuscript's Methods are written from. A stale number in
`CLAUDE.md` is not cosmetic: §6 exists because this project's bugs survive by
being described correctly somewhere nobody re-derives.

The claims live in `docs_claims.json`, each with the file and section that makes
it, so a failure names the prose to fix rather than just a number. Add an entry
whenever a doc states a number about the corpus. A claim with no entry here is a
claim nothing is checking.

    python3 -I provenance/check_docs.py          # exits non-zero on any mismatch

Note what this does NOT do: it does not parse the Markdown. Parsing prose for
numbers is brittle and fails open, which is the worst direction. The expected
value is transcribed here by hand, deliberately -- so updating a doc number
means updating this file too, and the pair either agrees or the build fails.
"""
import collections
import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))


def rows(name):
    with open(os.path.join(ROOT, "data", name)) as f:
        return list(csv.DictReader(f))


def n_names_multi_structure_cross_experiment(rs):
    """Names carrying DIFFERENT structures in DIFFERENT experiments.

    Deliberately not "names with more than one structure" -- that is 292 and is
    what error class 30 counts. This is 288. The two were conflated during the
    sweep that produced this file, by me, which is exactly why the distinction
    is spelled out in a function instead of inlined in a comprehension.
    """
    by = collections.defaultdict(lambda: collections.defaultdict(set))
    for r in rs:
        by[r["IL_name"]][r["Experiment_ID"]].add(r["IL_SMILES"])
    return sum(1 for per_exp in by.values()
               if len(per_exp) > 1
               and len({s for ss in per_exp.values() for s in ss}) > 1)


def n_registered_reagents():
    """Entries in the C11 registry itself, not in the audit's output file."""
    sys.path.insert(0, HERE)
    try:
        import reagent_registry as rr
    finally:
        sys.path.pop(0)
    for attr in ("REGISTRY", "REAGENTS", "NAMED"):
        if hasattr(rr, attr):
            return len(getattr(rr, attr))
    raise SystemExit("check_docs: cannot find the registry dict in "
                     "reagent_registry.py -- update n_registered_reagents()")


def main():
    spec = json.load(open(os.path.join(HERE, "docs_claims.json")))
    env = {
        "V": rows("data_vitro.csv"),
        "W": rows("data_vivo.csv"),
        "O": rows("Original_LNPDB.csv"),
        "A": json.load(open(os.path.join(HERE, "audit_findings.json"))),
        "S": json.load(open(os.path.join(HERE, "schema_quality.json")))["columns"],
        "n_names_multi_structure_cross_experiment":
            n_names_multi_structure_cross_experiment,
        "n_registered_reagents": n_registered_reagents,
    }

    bad = []
    print("C15 docs-vs-data")
    for c in spec["claims"]:
        try:
            got = eval(c["expr"], {"__builtins__": __builtins__}, env)  # noqa: S307
        except Exception as e:                       # a claim that cannot be
            got = f"ERROR: {type(e).__name__}: {e}"  # measured is a failure too
        ok = got == c["expect"]
        if ok:
            print(f"  OK    {c['id']:42} = {got}")
        else:
            bad.append((c, got))
            print(f"  FAIL  {c['id']:42} docs say {c['expect']}, "
                  f"corpus says {got}")

    print()
    if bad:
        print(f"{len(bad)} of {len(spec['claims'])} documented numbers DISAGREE "
              f"with the corpus:\n")
        for c, got in bad:
            print(f"  * {c['id']}: docs say {c['expect']}, corpus says {got}")
            print(f"    claimed in: {c['where']}")
            print(f"    measured by: {c['expr']}\n")
        print("Fix the prose AND this file together, or establish which is "
              "right first -- a doc and a checker that disagree are two claims, "
              "not one.")
        return 1

    print(f"DOCS MATCH DATA -- {len(spec['claims'])} claims checked")
    return 0


if __name__ == "__main__":
    sys.exit(main())
