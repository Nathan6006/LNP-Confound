#!/usr/bin/env python3
"""Derive `Experiment_ID -> lab` and write provenance/lab_map.json.

WHY THIS EXISTS (policy P4, 2026-11-10). `Experiment_ID` is a SCREEN key, not a
laboratory key, and the project had been reading it as a laboratory key -- the
manuscript title says "Cross-laboratory prediction". Those are different
groupings and the difference is large: 43 in-vitro + 20 in-vivo experiment-rows
come from only ~23 distinct laboratories, and the pair carrying 57% of
Result 3's pooled concordance points (`AA_2008` <-> `KW_2014`) is ONE lab.

So the laboratory has to be its own column, derived mechanically from a stored
record rather than remembered.

METHOD, and its limits, stated because this feeds a manuscript claim:

  lab := the LAST AUTHOR of the paper, from the NCBI esummary record cached in
  provenance/esummary_cache.json (refresh with fetch_esummary.py).

Last author is the field's convention for the principal investigator, and it is
mechanical, offline and reproducible, which an editorial judgement would not be.
It is NOT the same as the corresponding author, and it is wrong in at least
these ways, all of which are left visible rather than patched:

  * A paper with co-senior authors has one last author. `LNPDB`'s own paper has
    four starred corresponding authors.
  * Scientific lineage is not captured. `NC_2024` is Whitehead KA, who trained
    in Anderson DG's laboratory; `KW_2014` is Whitehead as FIRST author in
    Anderson's lab. Those are two laboratories by this rule and one tradition by
    any other reading. The cross-lab/same-lab split should be read with that in
    mind, and the manuscript should say so.
  * An author's initials vary between records. NCBI returned both
    "Anderson D" (PMID 36997680) and "Anderson DG" for the same person, which
    would have silently split one laboratory in two -- exactly the P5
    mis-grouping hazard. ALIASES below is the only hand-maintained part of this
    file, and every entry needs a reason.

Run: python3 -I provenance/build_lab_map.py
"""
import csv, json, os, sys, collections

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

# One person, two strings in NCBI's own records. Each entry must name the PMIDs
# that disagree, so a future reader can re-check rather than trust this list.
ALIASES = {
    # PMID 36997680 returns "Anderson D"; PMIDs 18438401, 20080679, 24969323,
    # 31570898, 37127709, 38082180 and others return "Anderson DG". Same person
    # (Daniel G. Anderson, MIT). Without this, BL_2023 would be its own lab and
    # the Anderson group would look like 8 labs instead of 9 experiments.
    "Anderson D": "Anderson DG",
}


def main():
    cache = json.load(open(os.path.join(HERE, "esummary_cache.json")))["result"]
    man = list(csv.DictReader(open(os.path.join(ROOT, "data",
                                                "experiment_manifest.csv"))))
    out, missing = {}, []
    for x in man:
        eid = x["Experiment_ID"]
        pmid = str(x["Publication_PMID"]).split(".")[0]
        rec = cache.get(pmid)
        if not rec:
            missing.append((eid, pmid))
            continue
        raw = (rec.get("lastauthor") or "").strip()
        if not raw:
            missing.append((eid, pmid))
            continue
        lab = ALIASES.get(raw, raw)
        prev = out.get(eid)
        if prev and prev["lab"] != lab:
            # One Experiment_ID appears once per pool; both rows must agree.
            print(f"CONFLICT: {eid} resolves to {prev['lab']!r} and {lab!r}",
                  file=sys.stderr)
            return 1
        out[eid] = {"lab": lab, "lab_raw": raw, "pmid": pmid,
                    "first_author": (rec.get("sortfirstauthor") or "").strip(),
                    "pubdate": (rec.get("pubdate") or "").strip()}

    if missing:
        # A blank lab would silently become its own laboratory of one, which is
        # the same failure as the alias problem in the other direction.
        print("REFUSING to write a partial map. No esummary record for:",
              file=sys.stderr)
        for eid, pmid in missing:
            print(f"  {eid} (PMID {pmid})", file=sys.stderr)
        print("Run: python3 -I provenance/fetch_esummary.py", file=sys.stderr)
        return 1

    labs = collections.Counter(v["lab"] for v in out.values())
    doc = {
        "_what": "Experiment_ID -> laboratory, derived from the LAST AUTHOR of "
                 "each paper's NCBI esummary record. Written by "
                 "provenance/build_lab_map.py; do not hand-edit.",
        "_method": "lab = esummary 'lastauthor', with provenance/build_lab_map.py"
                   "::ALIASES merging initials-variants of one person.",
        "_limits": "Last author is not always the corresponding author, and "
                   "scientific lineage is not captured (Whitehead KA trained in "
                   "Anderson DG's lab). See the module docstring.",
        "_counts": {"experiments": len(out), "distinct_labs": len(labs)},
        "_labs_with_more_than_one_experiment": {
            k: sorted(e for e, v in out.items() if v["lab"] == k)
            for k, n in labs.most_common() if n > 1},
        "map": dict(sorted(out.items())),
    }
    p = os.path.join(HERE, "lab_map.json")
    json.dump(doc, open(p, "w"), indent=1)
    print(f"wrote {p}: {len(out)} experiments, {len(labs)} distinct labs")
    for k, n in labs.most_common(8):
        print(f"  {k:24s} {n:2d}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
