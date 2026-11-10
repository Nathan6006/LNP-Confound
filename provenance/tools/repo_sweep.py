#!/usr/bin/env python3
"""Run find_repo.py over every queued paper. Throttled for GitHub code search
(10 req/min authenticated), so this is slow by design -- run it in background."""
import json, subprocess, sys, time, os
SP = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TOOL = os.path.join(SP, "tools", "find_repo.py")
q = json.load(open("provenance/run_queue.json"))["queue"]
lit = {str(r["pmid"]): r for r in json.load(open("lit_review/data.json"))["rows"]}
out = {}
for x in q:
    if x.get("status") == "duplicate_skip":
        continue
    r = lit.get(x["pmid"], {})
    args = [sys.executable, "-I", TOOL]
    if x.get("doi"): args += ["--doi", x["doi"]]
    if x.get("pmid"): args += ["--pmid", x["pmid"]]
    if x.get("title"): args += ["--title", x["title"]]
    if r.get("last"): args += ["--author", str(r["last"])]
    try:
        p = subprocess.run(args, capture_output=True, text=True, timeout=400)
        d = json.loads(p.stdout) if p.stdout.strip() else {"_err": p.stderr[:120]}
    except Exception as e:
        d = {"_err": str(e)[:120]}
    cands = d.get("candidates", {})
    out[x["folder"]] = d
    if cands:
        print(f"HIT  {x['folder']:<24} {list(cands)}", flush=True)
    else:
        print(f"  -  {x['folder']:<24}", flush=True)
    json.dump(out, open(os.path.join(SP, "repo_sweep.json"), "w"), indent=1)
    time.sleep(8)
print("done")
