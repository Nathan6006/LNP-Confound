#!/usr/bin/env python3
"""Fetch the NCBI esummary cache that build_manifest.py cross-checks against.

WHY THIS EXISTS. `build_manifest.py` read its cache from a path inside a Claude
session scratchpad:

    ESUM = "/private/tmp/claude-501/-Users-nathanliu-Documents-LNP-confound/
            <session-uuid>/scratchpad/esum.json"

Scratchpads are session-scoped and are cleared. By 2026-11-07 that file was
gone, so `build_manifest.py` raised FileNotFoundError and
`data/experiment_manifest.csv` -- which is tracked in git -- could no longer be
regenerated from anything in the repository. The manifest had become an
artifact with no reproducible lineage, which is the one thing this directory
exists to prevent (review finding F-023).

This is the second instance of the class. `lit_review/build_data.py` was
unrunnable for the same reason on 2026-10-30, and its inputs were recovered
into `lit_review/provenance/`. The rule that follows: **a script in this
repository may not read an input from a session scratchpad.** Scratch is for
intermediates nobody needs twice; anything a tracked output depends on gets
written into the repo.

The cache is a *cross-check*, not a source: the corpus's own PMID and DOI are
authoritative, and this is what they are verified against. It is cached rather
than fetched per run so the manifest build stays offline and deterministic, and
so a run years from now reproduces today's answer instead of silently picking up
a publisher's later metadata correction.

    python3 -I provenance/fetch_esummary.py          # refresh the cache

Writes provenance/esummary_cache.json. Re-run it when a new experiment with a
new PMID is merged; `build_manifest.py` notes "PMID <n> not checked against
NCBI" for anything missing rather than failing, so a stale cache degrades
loudly but safely.
"""
import csv
import json
import os
import sys
import time
import urllib.parse
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))
OUT = os.path.join(HERE, "esummary_cache.json")

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
# NCBI asks automated clients to identify themselves with a tool name. It also
# invites an email address; we deliberately send a generic client name and no
# address -- the user's email is not ours to hand to a third-party service.
TOOL = "lnp-confound-manifest"
BATCH = 100


def corpus_pmids():
    """Every distinct PMID in the emitted corpus, as a sorted list of strings."""
    pmids = set()
    for name in ("data_vitro.csv", "data_vivo.csv"):
        with open(os.path.join(ROOT, "data", name)) as f:
            for r in csv.DictReader(f):
                p = str(r.get("Publication_PMID", "")).strip()
                if p.endswith(".0"):
                    p = p[:-2]
                if p and p.lower() not in ("nan", "na", "none"):
                    pmids.add(p)
    return sorted(pmids, key=int)


def fetch(ids):
    q = urllib.parse.urlencode(
        {"db": "pubmed", "id": ",".join(ids), "retmode": "json", "tool": TOOL})
    with urllib.request.urlopen(f"{EUTILS}?{q}", timeout=60) as r:
        return json.load(r)


def main():
    ids = corpus_pmids()
    print(f"{len(ids)} distinct PMIDs in the corpus")

    merged = {"uids": []}
    for i in range(0, len(ids), BATCH):
        chunk = ids[i:i + BATCH]
        res = fetch(chunk).get("result", {})
        got = res.get("uids", [])
        for u in got:
            merged[u] = res[u]
        merged["uids"].extend(got)
        print(f"  fetched {len(got)}/{len(chunk)}")
        if i + BATCH < len(ids):
            time.sleep(0.4)          # NCBI: <=3 requests/sec without an API key

    missing = [p for p in ids if p not in merged]
    json.dump({"result": merged}, open(OUT, "w"), indent=1)
    print(f"wrote {OUT}  ({len(merged['uids'])} records, "
          f"{os.path.getsize(OUT) / 1024:.1f} kB)")
    if missing:
        print(f"NOT returned by NCBI ({len(missing)}): {missing}", file=sys.stderr)
        # Not fatal: build_manifest.py records an unchecked PMID as a note. A
        # PMID the corpus holds that PubMed does not know is itself a finding
        # worth seeing rather than an error worth crashing on.


if __name__ == "__main__":
    main()
