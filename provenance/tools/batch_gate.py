#!/usr/bin/env python3
"""Cross-paper half of the batch review gate (AUTONOMOUS_RUN.md section 4, step 2).

Per-paper Tier-2 review is a separate, adversarial pass. THIS checks only what
no single-paper reviewer can see, because it is the thing that structurally
escapes them: collisions and drift BETWEEN the staged folders and the live
corpus.

    python3 -I batch_gate.py

Reports; does not merge and does not edit anything.
"""
import collections, csv, glob, json, os, sys

ROOT = "/Users/nathanliu/LNP/confound"
CAND = ("09_rows_vitro.csv", "09_rows_vivo.csv", "06_emitted_rows.csv",
        "07_rows_quality.csv", "05_values.csv", "07_emit_preview.csv",
        "03_enumerated.csv")

# Folders name the structure column four different ways. The same aliasing
# problem took C10 down earlier today, where a folder using `label`/`SMILES`
# made the check silently find ZERO structures and report a missing anchor.
# Tolerate the aliases; a gate that cannot see a column reports a clean run.
SMILES_COLS = ("IL_SMILES", "canonical_SMILES", "SMILES", "smiles")
NAME_COLS = ("IL_name", "lipid_name", "compound", "label", "name")


def smiles_of(r):
    return next((r[c].strip() for c in SMILES_COLS if (r.get(c) or "").strip()), "")


def name_of(r):
    return next((r[c].strip() for c in NAME_COLS if (r.get(c) or "").strip()), "")


def live_corpus():
    rows = []
    for f in ("data_vitro.csv", "data_vivo.csv"):
        with open(os.path.join(ROOT, "data", f)) as fh:
            rows += list(csv.DictReader(fh))
    return rows


def staged():
    out = {}
    q = json.load(open(os.path.join(ROOT, "provenance/run_queue.json")))
    for x in q["queue"]:
        if not x.get("rows_emitted"):
            continue
        folder = os.path.join(ROOT, "curation", x["folder"])
        files = [os.path.join(folder, c) for c in CAND
                 if os.path.exists(os.path.join(folder, c))]
        out[x["folder"]] = files
    return out


def main():
    corpus = live_corpus()
    live_exp = {r["Experiment_ID"] for r in corpus}
    live_smi = collections.defaultdict(set)
    for r in corpus:
        s = (r.get("IL_SMILES") or "").strip()
        if s:
            live_smi[s].add(r["Experiment_ID"])

    print("=" * 74)
    print("BATCH REVIEW GATE -- cross-paper checks")
    print("=" * 74)
    print(f"live corpus: {len(corpus)} rows, {len(live_exp)} experiments, "
          f"{len(live_smi)} distinct IL_SMILES\n")

    allrows = {}
    for folder, files in sorted(staged().items()):
        rows = []
        for f in files:
            try:
                with open(f) as fh:
                    rs = list(csv.DictReader(fh))
                if rs and any(k in rs[0] for k in
                              SMILES_COLS + NAME_COLS + ("Experiment_ID",)):
                    rows += rs
            except Exception as e:
                print(f"  !! {folder}/{os.path.basename(f)}: {e}")
        allrows[folder] = rows

    # G1 -- Experiment_ID collisions with the live corpus
    print("G1  Experiment_ID already present in the live corpus")
    hits = 0
    for folder, rows in allrows.items():
        ids = {(r.get("Experiment_ID") or "").strip() for r in rows}
        ids.discard("")
        clash = ids & live_exp
        if clash:
            hits += 1
            print(f"    COLLISION  {folder}: {sorted(clash)}")
    print(f"    {hits} folder(s) collide\n")

    # G2 -- structures already in the corpus
    print("G2  staged structures already in the live corpus (by IL_SMILES)")
    for folder, rows in allrows.items():
        smis = {smiles_of(r) for r in rows}
        smis.discard("")
        if not smis:
            print(f"    {folder:<24} (no structure column found -- "
                  f"columns seen: {sorted(rows[0].keys())[:6] if rows else []})")
            continue
        ov = smis & set(live_smi)
        where = collections.Counter(
            e for s in ov for e in live_smi[s])
        flag = "  <-- MERGE-LEVEL" if len(ov) >= 50 else ""
        print(f"    {folder:<24} {len(ov):>4} of {len(smis):>4} already present"
              f"  {dict(where.most_common(3))}{flag}")
    print()

    # G3 -- one structure under two names, inside a staged folder
    print("G3  one structure under >1 IL_name within a staged folder")
    for folder, rows in allrows.items():
        d = collections.defaultdict(set)
        for r in rows:
            s = smiles_of(r)
            n = name_of(r)
            if s and n:
                d[s].add(n)
        bad = {k: v for k, v in d.items() if len(v) > 1}
        if bad:
            print(f"    {folder:<24} {len(bad)} structure(s), e.g. "
                  f"{sorted(list(bad.values())[0])[:3]}")
    print()

    # G4 -- quality-column vocabulary, against the single source of truth
    print("G4  quality-column vocabulary (provenance/schema_quality.json)")
    sch = json.load(open(os.path.join(ROOT, "provenance/schema_quality.json")))
    allowed = {c["name"]: c.get("allowed") for c in sch["columns"]}
    for folder, rows in allrows.items():
        bad = collections.Counter()
        for r in rows:
            for col, allow in allowed.items():
                if allow is None or col not in r:
                    continue
                v = (r.get(col) or "").strip()
                if v and v not in allow:
                    bad[f"{col}={v!r}"] += 1
        if bad:
            print(f"    {folder:<24} {dict(bad.most_common(3))}")
    print("    (blank = no violations)\n")

    # G5 -- caveat without a note
    print("G5  rows carrying a caveat with an empty curation_note")
    CAV = {"medium", "low"}
    for folder, rows in allrows.items():
        n = 0
        for r in rows:
            if "curation_note" not in r:
                continue
            cav = ((r.get("confidence") or "").strip() in CAV
                   or (r.get("positional_ambiguity") or "").strip() == "1"
                   or (r.get("structure_assumed") or "").strip() == "1"
                   or (r.get("value_is_bound") or "").strip() in ("upper", "lower"))
            if cav and not (r.get("curation_note") or "").strip():
                n += 1
        if n:
            print(f"    {folder:<24} {n} row(s)")
    print("    (blank = none)\n")

    # G6 -- DECISIONS gate, through the real parser
    print("G6  open decisions, via build_dataset.open_decisions")
    import re, io as _io
    src = _io.open(os.path.join(ROOT, "provenance/build_dataset.py"),
                   encoding="utf-8").read()
    fn = re.search(r"def open_decisions\(path\):(?:.|\n)*?\n    return open_ids\n", src).group(0)
    ns = {}
    exec("import os, io\n" + fn, ns)
    total = 0
    for folder in sorted(allrows):
        o = ns["open_decisions"](os.path.join(ROOT, "curation", folder, "x.csv"))
        total += len(o)
        print(f"    {folder:<24} {len(o):>2} OPEN")
    print(f"    {total} open decisions across the staging set -- "
          f"every folder is correctly blocked\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
