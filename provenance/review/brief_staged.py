#!/usr/bin/env python3
"""Generate a review packet for a STAGED curation folder — live facts, not prose.

`brief.py` reads `data/`, which only works for merged papers. The 2026-11-10 run
staged ~7,200 rows across 27 folders that are deliberately NOT merged (every one
carries OPEN decisions, so `build_dataset.py` refuses them). This reads each
folder's own emitted-rows file instead, and recomputes everything the reviewer
would otherwise have to take on trust.

Why generated and not written: a prose brief goes stale at the next edit, and a
stale brief is worse than none — it makes a reviewer confidently verify a number
that is no longer in the data. If a packet disagrees with the repository, the
packet is wrong, and saying so is a valid finding.

    python3 -I provenance/review/brief_staged.py <folder> [<folder>...]
"""
import collections
import csv
import glob
import json
import os
import sys

from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors as D

RDLogger.DisableLog("rdApp.*")
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))


def key(smi):
    m = Chem.MolFromSmiles(smi) if smi else None
    return Chem.MolToInchiKey(m) if m else None


def corpus_keys():
    ks = {}
    for fn in ("data_vitro.csv", "data_vivo.csv"):
        p = os.path.join(ROOT, "data", fn)
        for r in csv.DictReader(open(p)):
            k = key(r.get("IL_SMILES"))
            if k:
                ks.setdefault(k, []).append((r["Experiment_ID"],
                                             (r.get("IL_name") or "").strip()))
    return ks


def emitted_file(folder, want):
    """The file whose row count equals rows_emitted — by SCHEMA, not filename.

    Folders named their output `06_rows.csv`, `06_emitted_rows.csv`,
    `07_<EID>_rows.csv`, `09_rows_vitro.csv` or just `05_values.csv`. A scan
    keyed on one filename found 14 of 27 and reported a lipid count wrong by a
    factor of fourteen.
    """
    cands = []
    for p in set(glob.glob(f"{ROOT}/curation/{folder}/*.csv")) | \
             set(glob.glob(f"{ROOT}/curation/{folder}/**/*.csv", recursive=True)):
        try:
            with open(p) as fh:
                cols = next(csv.reader(fh), None)
            if not cols:
                continue
            cs = set(cols)
            if "IL_SMILES" not in cs or not any(
                    c.startswith("Experiment_value") for c in cs):
                continue
            n = sum(1 for _ in csv.DictReader(open(p)))
            if n:
                cands.append((p, n))
        except Exception:
            continue
    if not cands:
        return None, 0, cands
    exact = [c for c in cands if c[1] == want]
    pick = exact[0] if exact else max(cands, key=lambda t: t[1])
    return pick[0], pick[1], cands


def packet(folder, qrow, ck):
    want = qrow.get("rows_emitted")
    path, n, cands = emitted_file(folder, want)
    L = []
    A = L.append
    A(f"# Review packet — `{folder}`")
    A("")
    A(f"Generated from live files. Queue says **{want}** emitted rows; "
      f"status `{qrow.get('status')}`; proposed `Experiment_ID` "
      f"**{qrow.get('proposed_experiment_id') or 'none'}**; PMID "
      f"{qrow.get('pmid')}.")
    A("")
    if not path:
        A("**NO SCHEMA-VALID ROW FILE FOUND** (needs `IL_SMILES` and a column "
          "starting `Experiment_value`). That is itself a finding: this folder "
          "cannot be merged as it stands.")
        A("")
        A("Files that look row-like but lack one of those columns:")
        for p in sorted(glob.glob(f"{ROOT}/curation/{folder}/*.csv")):
            with open(p) as fh:
                cols = next(csv.reader(fh), [])
            A(f"  - `{os.path.basename(p)}` — {len(cols)} cols, "
              f"IL_SMILES={'IL_SMILES' in cols}, "
              f"value={[c for c in cols if c.startswith('Experiment_value')]}")
        return "\n".join(L)

    rel = os.path.relpath(path, ROOT)
    A(f"## Rows read from `{rel}` — {n} rows"
      + ("" if n == want else f"  **(queue says {want} — DISCREPANCY)**"))
    if len(cands) > 1:
        A("")
        A("Other schema-valid files in this folder (a reviewer should know "
          "which one the merge would take):")
        for p, c in sorted(cands, key=lambda t: -t[1]):
            A(f"  - `{os.path.relpath(p, ROOT)}` — {c} rows")
    rows = list(csv.DictReader(open(path)))
    cols = list(rows[0]) if rows else []
    A("")
    A(f"**{len(cols)} columns.** Corpus schema is 79; this file has "
      f"{len(set(cols) & set(next(csv.reader(open(os.path.join(ROOT,'data','data_vitro.csv'))))))} "
      f"of them.")

    # structures
    ks = collections.Counter()
    bad = 0
    for r in rows:
        k = key(r.get("IL_SMILES"))
        if k:
            ks[k] += 1
        elif (r.get("IL_SMILES") or "").strip():
            bad += 1
    A("")
    A("## Structures")
    A(f"- distinct full InChIKeys: **{len(ks)}**")
    A(f"- unparseable `IL_SMILES`: **{bad}**")
    already = {k: v for k, v in ks.items() if k in ck}
    A(f"- already in the corpus: **{len(already)}** of {len(ks)}")
    if already:
        A("  | InChIKey | rows here | in corpus as |")
        A("  |---|---|---|")
        for k in list(already)[:8]:
            where = sorted({f"{e}:{nm}" for e, nm in ck[k]})[:3]
            A(f"  | `{k}` | {ks[k]} | {', '.join(where)} |")
    A(f"- **NEW to the corpus: {len(ks) - len(already)}**")

    # values
    A("")
    A("## Values")
    vcol = next((c for c in cols if c == "Experiment_value"), None) or \
        next((c for c in cols if c.startswith("Experiment_value")), None)
    vals = []
    for r in rows:
        try:
            vals.append(float(r[vcol]))
        except Exception:
            pass
    if vals:
        mean = sum(vals) / len(vals)
        sd = (sum((v - mean) ** 2 for v in vals) / len(vals)) ** 0.5
        A(f"- `{vcol}`: n={len(vals)}, mean **{mean:.6g}**, sd(pop) "
          f"**{sd:.6g}**, min {min(vals):.4g}, max {max(vals):.4g}")
        A("  - *mean 0 / sd 1 is a TAUTOLOGY after z-scoring — it is not "
          "evidence of anything.* What matters is whether the z was taken over "
          "the right block (error class 22).")
    for c in ("Experiment_ID", "Model_type", "Model_target", "Cargo_type",
              "Experiment_method", "Experiment_batching", "value_normalization_unit"):
        if c in cols:
            cc = collections.Counter((r.get(c) or "").strip() for r in rows)
            A(f"- `{c}`: {dict(cc.most_common(6))}")

    # quality flags
    A("")
    A("## Quality columns (schema_quality.json)")
    for c in ("confidence", "positional_ambiguity", "structure_assumed",
              "value_is_bound", "curation_note"):
        if c in cols:
            cc = collections.Counter((r.get(c) or "").strip() for r in rows)
            nonblank = sum(v for k, v in cc.items() if k)
            A(f"- `{c}`: {nonblank}/{len(rows)} populated — "
              f"{dict(cc.most_common(5))}")
        else:
            A(f"- `{c}`: **COLUMN ABSENT**")

    # folder inventory + decisions
    A("")
    A("## Folder inventory")
    files = sorted(glob.glob(f"{ROOT}/curation/{folder}/**/*", recursive=True))
    for p in files:
        if os.path.isfile(p) and not any(
                s in p for s in ("/source", "__pycache__", "/figs/", "/renders/")):
            A(f"  - `{os.path.relpath(p, ROOT)}` ({os.path.getsize(p)} B)")
    dp = f"{ROOT}/curation/{folder}/DECISIONS.md"
    if os.path.exists(dp):
        txt = open(dp).read()
        opens = txt.count("STATUS: OPEN")
        A("")
        A(f"## DECISIONS.md — **{opens} OPEN**")
        for line in txt.splitlines():
            if line.startswith("### D-"):
                A(f"  - {line[4:].strip()}")
    cp = f"{ROOT}/curation/{folder}/08_contradictions.csv"
    if os.path.exists(cp):
        cr = list(csv.DictReader(open(cp)))
        A("")
        A(f"## 08_contradictions.csv — {len(cr)} rows")
    return "\n".join(L)


def main():
    q = {x["folder"]: x for x in json.load(
        open(os.path.join(ROOT, "provenance", "run_queue.json")))["queue"]}
    ck = corpus_keys()
    outdir = os.path.join(HERE, "round_2026-11-10")
    os.makedirs(outdir, exist_ok=True)
    for folder in sys.argv[1:]:
        if folder not in q:
            print(f"  {folder}: NOT IN QUEUE — refusing", file=sys.stderr)
            continue
        txt = packet(folder, q[folder], ck)
        p = os.path.join(outdir, f"packet_{folder}.md")
        open(p, "w").write(txt + "\n")
        print(f"  wrote {os.path.relpath(p, ROOT)} ({len(txt)} chars)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
