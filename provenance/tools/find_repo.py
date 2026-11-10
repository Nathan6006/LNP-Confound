#!/usr/bin/env python3
"""Tier-3 of the intake gate: find a paper's data on GitHub.

    python3 -I find_repo.py --doi 10.1002/adhm.202400225 --pmid 38888972 \
                            --title "Branched-tail lipid nanoparticles" \
                            --author Whitehead

Why this belongs in stage 00: a lab that deposits its screen on GitHub has
given us the values as MACHINE-READABLE DATA, which removes the whole
digitisation error class. `WT_2026` came from the authors' own repository and
has zero digitisation defects; the four curations that digitised produced
eight value-changing findings between them.

Searches, in decreasing precision:
  1. code search for the DOI string (a repo that cites its own paper)
  2. code search for the PMID
  3. repo search on title keywords, optionally narrowed by author/lab name
Then lists any data-shaped files (.csv/.xlsx/.json/.smi/.sdf/.txt) in each hit.

Uses the `gh` CLI, which carries the user's existing auth. Read-only.
"""
import argparse, json, subprocess, sys

DATA_EXT = (".csv", ".xlsx", ".xls", ".json", ".smi", ".sdf", ".tsv", ".txt", ".parquet")


def gh(path, **params):
    args = ["gh", "api", "-X", "GET", path]
    for k, v in params.items():
        args += ["-f", f"{k}={v}"]
    r = subprocess.run(args, capture_output=True, text=True, timeout=90)
    if r.returncode != 0:
        return {"_err": r.stderr.strip()[:160]}
    try:
        return json.loads(r.stdout)
    except Exception as e:
        return {"_err": str(e)[:120]}


def code_search(q):
    d = gh("search/code", q=q, per_page="20")
    if "_err" in d:
        return [], d["_err"]
    return [{"repo": i["repository"]["full_name"], "path": i["path"]}
            for i in d.get("items", [])], None


def repo_search(q):
    d = gh("search/repositories", q=q, per_page="12", sort="stars")
    if "_err" in d:
        return [], d["_err"]
    return [{"repo": i["full_name"], "desc": (i.get("description") or "")[:90],
             "stars": i.get("stargazers_count", 0),
             "pushed": (i.get("pushed_at") or "")[:10]}
            for i in d.get("items", [])], None


def tree(repo):
    d = gh(f"repos/{repo}/git/trees/HEAD", recursive="1")
    if "_err" in d:
        return []
    return [{"path": t["path"], "size": t.get("size", 0)}
            for t in d.get("tree", [])
            if t["type"] == "blob" and t["path"].lower().endswith(DATA_EXT)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--doi"); ap.add_argument("--pmid")
    ap.add_argument("--title"); ap.add_argument("--author")
    ap.add_argument("--extra")
    a = ap.parse_args()

    out = {"queries": [], "candidates": {}}
    hits = []
    for label, q in (("doi", f'"{a.doi}"' if a.doi else None),
                     ("pmid", f'"{a.pmid}"' if a.pmid else None),
                     ("extra", a.extra)):
        if not q:
            continue
        items, err = code_search(q)
        out["queries"].append({"kind": "code", "q": q,
                               "n": len(items), "err": err})
        hits += [i["repo"] for i in items]

    # Repo search: a long title finds nothing (GitHub ANDs the terms and repo
    # names are short). Validated against a repo we already knew existed --
    # "FALCON lipid nanoparticle optimization Mao" returned 0 while a 2-3 word
    # query finds it. So try several short queries, most specific first.
    if a.title:
        words = [w for w in a.title.replace("-", " ").split() if len(w) > 3]
        tries = []
        if a.author:
            tries.append(f"{a.author} lipid nanoparticle")
            tries.append(f"{a.author} LNP")
        tries.append(" ".join(words[:3]))
        tries.append(" ".join(words[:2]) + " lipid")
        seen_q = set()
        out["repo_hits"] = []
        for q in tries:
            q = q.strip()
            if not q or q in seen_q:
                continue
            seen_q.add(q)
            items, err = repo_search(q)
            out["queries"].append({"kind": "repo", "q": q,
                                   "n": len(items), "err": err})
            out["repo_hits"] += items
            hits += [i["repo"] for i in items]

    # Relevance filter. A bare DOI/PMID code-search hit is often an unrelated
    # repo that happens to cite a paper, so require the repo name or its data
    # file names to look like this field before reporting it as a candidate.
    import re
    REL = re.compile(r"lipid|lnp|mrna|nanopart|lipidoid|ionizable|transfect|"
                     r"delivery|sirna|formulation", re.I)
    for r in dict.fromkeys(hits):
        files = tree(r)
        if not files:
            continue
        relevant = REL.search(r) or any(REL.search(f["path"]) for f in files)
        if not relevant:
            out.setdefault("filtered_out", []).append(r)
            continue
        out["candidates"][r] = sorted(files, key=lambda f: -f["size"])[:15]
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    sys.exit(main())
