#!/usr/bin/env python3
"""Prove curation/lib/emit.py's two refusals fire on THIS curation's data.

A guarantee that is never exercised is a comment. Each seed runs in its own
process (seed_emit_one.py) and the returncode is read, not inferred.
"""
import json, os, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__))
ONE = os.path.join(HERE, "seed_emit_one.py")
SEEDS = [
    ("control", "the real per-item callables", 0),
    ("literal_confidence", "confidence='high' as a LITERAL (error class 29 / "
     "F-011: a literal relabelled 24 'ASSUMED' lipids as the project's "
     "highest confidence)", 1),
    ("caveat_without_note", "confidence='medium' on an 18PG row with an empty "
     "curation_note (F-024: found 6 such rows on its first run)", 1),
    ("bad_vocabulary", "confidence='very high', outside "
     "provenance/schema_quality.json's allowed set", 1),
    ("missing_record", "an item with no enumeration record behind its flag", 1),
]
out = []
for name, what, expect in SEEDS:
    r = subprocess.run([sys.executable, "-I", ONE, name],
                       capture_output=True, text=True)
    out.append({"seed": name, "what": what, "expected_returncode": expect,
                "returncode": r.returncode,
                "raised": r.returncode != 0,
                "CREDITED": r.returncode == expect,
                "message_head": (r.stdout + r.stderr).strip().splitlines()[-1][:160]
                if (r.stdout + r.stderr).strip() else ""})
res = {"seeds": out, "n": len(out),
       "n_credited": sum(1 for s in out if s["CREDITED"]),
       "all_credited": all(s["CREDITED"] for s in out)}
json.dump(res, open(os.path.join(HERE, "..", "04_emit_seed_results.json"), "w"),
          indent=1)
print(json.dumps(res, indent=1))
sys.exit(0 if res["all_credited"] else 1)
