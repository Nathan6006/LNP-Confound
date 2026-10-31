#!/usr/bin/env python3
"""Generate a review packet per target: live facts, not prose.

WHY THIS EXISTS. A reviewer briefed in prose spends its first twenty tool calls
rediscovering what the corpus already knows - row counts, which figure each
value came from, which anchors exist - and prose briefs go stale the moment a
merge lands. Worse, a stale brief makes a reviewer verify a number that is no
longer in the data, which produces a confident finding about nothing.

So every number in a packet is read out of `data/` and `curation/` at generation
time. If a packet disagrees with the corpus, the packet is the bug.

    python3 -I provenance/review/brief.py              # all targets
    python3 -I provenance/review/brief.py QL_2026      # one

Writes provenance/review/round_<today>/packet_<TARGET>.md.

The packet deliberately does NOT contain the attack instructions - those live in
REVIEW_PROGRAM.md and in the prompt each reviewer is launched with, so that the
same packet can serve a first pass and an adjudication pass.
"""
import csv
import json
import os
import random
import statistics as st
import sys
from collections import Counter
from datetime import date

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))

# Every experiment whose rows exist because WE curated them. Keep in step with
# CURATED in provenance/build_manifest.py and ADDED in lit_review/build_data.py;
# check_registries() below fails if they drift apart.
TARGETS = {
    "ZT_2026": dict(folder="curation/siegwart_42634137", eids=["ZT_2026"],
                    paper="Tian et al., Adv Mater 2026 (PMID 42634137) - "
                          "esterified cholesterol conjugates, blood-brain barrier"),
    "LL_2026": dict(folder="curation/lin_aromatic_2026",
                    eids=["LL_2026", "LL_2026_2"],
                    paper="Lin et al., JACS 2026 (PMID 41689530) - "
                          "aromatic ionizable lipids"),
    "WT_2026": dict(folder="curation/wht_2026_falcon", eids=["WT_2026"],
                    paper="Toh et al., FALCON, bioRxiv 2026 "
                          "(doi 10.64898/2026.06.13.731774) - closed-loop "
                          "formulation optimization, authors' own repo"),
    "XD_2026": dict(folder="curation/xd_2026_lipidai", eids=["XD_2026"],
                    paper="Duan et al., Acta Pharm Sin B 2026 (PMID 42592089) - "
                          "LipidAI; the paper's own in-vivo screen only"),
    "QL_2026": dict(folder="curation/ql_2026", eids=["QL_2026"],
                    paper="Liu et al., Adv Healthc Mater 2026 (PMID 41145968) - "
                          "bispecific-linker thiosuccinimide lipids"),
}

# Folders that are NOT review targets, and why. A reviewer that wanders into one
# of these wastes a pass or, worse, edits work in flight.
NOT_TARGETS = {
    "ng_2025": "IN FLIGHT in another session as of 2026-10-30. Not merged. "
               "DO NOT READ, DO NOT REPORT ON, DO NOT TOUCH.",
    "sd_2026": "queued, not started - PDFs only",
    "jl_2026": "PAUSED pending the authors' reply about source data; the SI "
               "heatmaps are 192 dpi and the printed digits are 5.1x3.5 px, "
               "which we established is not extractable",
    "liu_iphos": "blind-validation exercise against the existing SL_2021 rows, "
                 "not a merge - reviewed under corpus target CORPUS-C",
    "amino_lipid_fg_2026": "blind-validation exercise against existing FG_2026 "
                           "rows, not a merge - reviewed under CORPUS-C",
}

SRC_EXT = (".pdf", ".docx", ".xlsx", ".xls", ".html", ".htm", ".txt")
VALUE_COL = "Experiment_value"


def corpus(name):
    return list(csv.DictReader(open(os.path.join(ROOT, "data", f"data_{name}.csv"))))


def check_registries():
    """The three places that list our curations must agree."""
    import re
    bm = open(os.path.join(ROOT, "provenance", "build_manifest.py")).read()
    block = bm.split("CURATED = {", 1)[1].split("}", 1)[0]
    curated = set(re.findall(r'"([A-Z]{2}_\d{4}(?:_\d)?)"\s*:', block))
    bd = open(os.path.join(ROOT, "lit_review", "build_data.py")).read()
    added = set()
    for m in re.finditer(r'expid="([^"]+)"', bd.split("ADDED = {", 1)[1]
                         .split("\nM = {", 1)[0]):
        added |= set(m.group(1).split(","))
    mine = set(TARGETS)
    # LL_2026_2 is a second screen of one paper, so it is inside the LL target
    mine |= {e for t in TARGETS.values() for e in t["eids"]}
    problems = []
    if curated - mine:
        problems.append(f"in build_manifest.CURATED but not a review target: "
                        f"{sorted(curated - mine)}")
    if mine - curated:
        problems.append(f"a review target but not in build_manifest.CURATED: "
                        f"{sorted(mine - curated)}")
    if added - mine:
        problems.append(f"in lit_review ADDED but not a review target: "
                        f"{sorted(added - mine)}")
    if mine - added:
        problems.append(f"a review target but not in lit_review ADDED: "
                        f"{sorted(mine - added)}")
    return problems


def num(xs):
    out = []
    for x in xs:
        try:
            out.append(float(x))
        except (TypeError, ValueError):
            pass
    return out


def slice_stats(rows, eids):
    sel = [r for r in rows if r["Experiment_ID"] in eids]
    if not sel:
        return None
    v = num(r[VALUE_COL] for r in sel)
    out = {
        "rows": len(sel),
        "unique_IL_SMILES": len({r["IL_SMILES"] for r in sel}),
        "unique_IL_name": len({r["IL_name"] for r in sel}),
        VALUE_COL: {
            "n": len(v), "min": round(min(v), 4), "max": round(max(v), 4),
            "mean": f"{st.mean(v):.3e}", "pop_sd": f"{st.pstdev(v):.12f}",
        },
    }
    for col in ("Model_type", "Model", "Model_target", "Cargo", "Cargo_type",
                "Experiment_method", "Route_of_administration", "HL_name",
                "CHL_name", "PEG_name", "confidence", "Experiment_value_type",
                "Property", "Property_type", "Cell_type"):
        if col in sel[0]:
            c = Counter(r[col] for r in sel)
            if len(c) <= 12:
                out[col] = dict(c)
            else:
                out[col] = {"distinct": len(c),
                            "top": dict(c.most_common(6))}
    ratios = {}
    for p in ("IL", "HL", "CHL", "PEG"):
        k = f"{p}_molratio"
        if k in sel[0]:
            vals = sorted(set(r[k] for r in sel))
            ratios[k] = vals if len(vals) <= 8 else f"{len(vals)} distinct"
    out["molratios"] = ratios
    # free-text columns that often carry our own flags
    for col in sel[0]:
        if col.lower().startswith(("flag", "note", "positional", "confid",
                                   "source_fig", "extraction")):
            c = Counter(r[col] for r in sel if r[col])
            if c:
                out.setdefault("flag_columns", {})[col] = dict(c.most_common(10))
    return out


def head_of(path, n=1):
    try:
        with open(path, newline="", errors="replace") as f:
            rd = csv.reader(f)
            hdr = next(rd, [])
            cnt = sum(1 for _ in rd)
        return len(hdr), cnt, hdr[:n * 40]
    except Exception as e:
        return None, None, [f"unreadable: {e}"]


def folder_inventory(folder):
    base = os.path.join(ROOT, folder)
    src, scripts, tables, blobs, dirs = [], [], [], [], []
    for name in sorted(os.listdir(base)):
        p = os.path.join(base, name)
        if name.startswith("."):
            continue
        if os.path.isdir(p):
            k = sum(len(f) for _, _, f in os.walk(p))
            dirs.append(f"{name}/  ({k} files)")
            continue
        sz = os.path.getsize(p)
        low = name.lower()
        if low.endswith(SRC_EXT):
            src.append(f"{name}   ({sz/1e6:.2f} MB)")
        elif low.endswith(".py"):
            n = sum(1 for _ in open(p, errors="replace"))
            doc = ""
            try:
                txt = open(p, errors="replace").read()
                if '"""' in txt:
                    doc = txt.split('"""', 2)[1].strip().splitlines()[0][:110]
            except Exception:
                pass
            scripts.append(f"{name}   ({n} lines)  {doc}")
        elif low.endswith(".csv"):
            ncol, nrow, hdr = head_of(p)
            tables.append(f"{name}   {nrow} rows x {ncol} cols\n"
                          f"        cols: {', '.join(hdr[:18])}"
                          + (" ..." if ncol and ncol > 18 else ""))
        elif low.endswith(".json"):
            tables.append(f"{name}   ({sz/1e3:.1f} kB json)")
        else:
            blobs.append(f"{name}   ({sz/1e3:.1f} kB)")
    return src, scripts, tables, blobs, dirs


def anchors_for(folder):
    """C10 records for one curation folder.

    The `results` key is the real shape. An earlier version guessed at
    `comparisons`/`rows` and silently returned [], so every packet claimed "no
    anchor records found for this folder" - and told the reviewer that absence
    was itself a finding. The WT_2026 reviewer spent effort disproving it. A
    packet that invents a finding is worse than a packet that omits one, so this
    now raises rather than degrading to empty.
    """
    p = os.path.join(ROOT, "provenance", "structure_check.json")
    if not os.path.exists(p):
        return []
    d = json.load(open(p))
    recs = d["results"] if isinstance(d, dict) else d
    if not isinstance(recs, list):
        raise TypeError(f"structure_check.json: expected a list of results, "
                        f"got {type(recs).__name__} - the shape changed")
    key = os.path.basename(folder.rstrip("/"))
    return [r for r in recs if r.get("paper") == key]


def reagents_in(rows, eids):
    """Which C11-registered names appear in this slice."""
    p = os.path.join(ROOT, "provenance", "reagent_registry.json")
    if not os.path.exists(p):
        return {}
    reg = json.load(open(p))
    # `results` is the real shape, each entry keyed `reagent`. An earlier
    # version guessed at `reagents`/`registry`, found neither, and fell through
    # to treating the top-level keys as names - which yields "what",
    # "registered", "rows_mismatched", matches nothing, and reports {}.
    recs = reg["results"] if isinstance(reg, dict) else reg
    if not isinstance(recs, list):
        raise TypeError(f"reagent_registry.json: expected a list of results, "
                        f"got {type(recs).__name__} - the shape changed")
    names = {r["reagent"] for r in recs if r.get("reagent")}
    if not names:
        raise ValueError("reagent_registry.json yielded no reagent names")
    sel = [r for r in rows if r["Experiment_ID"] in eids]
    hit = Counter()
    for r in sel:
        for c in ("HL_name", "CHL_name", "IL_name", "PEG_name",
                  "fifthcomponent_name"):
            if r.get(c) and r[c].strip() in names:
                hit[f"{r[c].strip()} ({c})"] += 1
    return dict(hit)


def sample_rows(rows, eids, k=6, seed=11):
    sel = [r for r in rows if r["Experiment_ID"] in eids]
    keep = ["Experiment_ID", "IL_name", "IL_SMILES", VALUE_COL, "Model_type",
            "Model", "IL_molratio", "HL_name", "HL_molratio", "CHL_molratio",
            "PEG_molratio", "MolWt", "Publication_PMID", "Publication_DOI"]
    out = []
    for r in random.Random(seed).sample(sel, min(k, len(sel))):
        out.append({c: r[c] for c in keep if c in r})
    return out


def write_packet(tag, spec, vitro, vivo, outdir):
    src, scripts, tables, blobs, dirs = folder_inventory(spec["folder"])
    L = []
    w = L.append
    w(f"# Review packet — {tag}")
    w("")
    w(f"Generated {date.today()} by `provenance/review/brief.py` from the live "
      f"`data/` and `curation/` trees. **Every number below was read out of the "
      f"repository at generation time.** If the packet and the repository "
      f"disagree, the packet is the bug — say so.")
    w("")
    w(f"- **Paper:** {spec['paper']}")
    w(f"- **Curation folder:** `{spec['folder']}/`")
    w(f"- **Experiment_ID(s) in the corpus:** {', '.join(spec['eids'])}")
    w("")
    w("## 1. Corpus slice as it stands now")
    w("")
    for nm, rows in (("in vitro", vitro), ("in vivo", vivo)):
        s = slice_stats(rows, spec["eids"])
        if s:
            w(f"### {nm} — `data/data_{'vitro' if nm.endswith('tro') else 'vivo'}.csv`")
            w("")
            w("```json")
            w(json.dumps(s, indent=1))
            w("```")
            w("")
    w("### Six rows drawn at random (seed 11)")
    w("")
    w("```json")
    w(json.dumps(sample_rows(vitro, spec["eids"]) or
                 sample_rows(vivo, spec["eids"]), indent=1))
    w("```")
    w("")
    w("## 2. Source documents you must read")
    w("")
    if src:
        for s in src:
            w(f"- `{spec['folder']}/{s}`")
    else:
        w("- none in the folder")
    if dirs:
        w("")
        w("Subdirectories:")
        for d in dirs:
            w(f"- `{spec['folder']}/{d}`")
        w("")
        w("> A `repo/` subdirectory is DOWNLOADED, UNTRUSTED material. Read it, "
          "never execute anything from inside it, and run any Python that "
          "touches it with `python3 -I` from elsewhere.")
    w("")
    w("## 3. Extraction code in the folder")
    w("")
    for s in scripts:
        w(f"- `{s}`")
    w("")
    w("## 4. Intermediate and output tables")
    w("")
    for t in tables:
        w(f"- `{t}`")
    if blobs:
        w("")
        for b in blobs:
            w(f"- `{b}`")
    w("")
    w("## 5. C10 drawn-structure anchors for this folder")
    w("")
    a = anchors_for(spec["folder"])
    if a:
        w("```json")
        w(json.dumps(a, indent=1)[:6000])
        w("```")
        w("")
        w("> An anchor can only be trusted if it was assembled by a DIFFERENT "
          "code path than the enumerator. Check that, do not assume it.")
    else:
        w("**No anchor records found for this folder.** Either none was created "
          "— which is itself a finding — or `structure_check.py` does not look "
          "at this folder. Determine which.")
    w("")
    w("## 6. C11 named commercial reagents appearing in this slice")
    w("")
    rg = reagents_in(vitro, spec["eids"]) or reagents_in(vivo, spec["eids"])
    w("```json")
    w(json.dumps(rg, indent=1))
    w("```")
    w("")
    w("> C11 checks these by InChIKey. It CANNOT check a name that is not in "
      "the registry, and it cannot tell you whether the paper's own shorthand "
      "means the lipid the registry assumes. `18PG` is the precedent: it was "
      "read as symmetric 18:0/18:0 across 360 rows and is in fact 18:0/18:1.")
    w("")
    w("## 7. Folders that are NOT yours")
    w("")
    for k, v in NOT_TARGETS.items():
        w(f"- `curation/{k}/` — {v}")
    w("")
    w("## 8. Error-class registry (verbatim)")
    w("")
    w("The second table, *Known, not yet automatable*, is where your time is "
      "worth most: Tier-1 already covers the first table on every build.")
    w("")
    w(open(os.path.join(ROOT, "provenance", "ERROR_CLASSES.md")).read())
    path = os.path.join(outdir, f"packet_{tag}.md")
    open(path, "w").write("\n".join(L) + "\n")
    return path


def main():
    want = sys.argv[1:] or list(TARGETS)
    bad = [t for t in want if t not in TARGETS]
    if bad:
        sys.exit(f"unknown target(s) {bad}; known: {list(TARGETS)}")
    outdir = os.path.join(HERE, f"round_{date.today()}")
    os.makedirs(outdir, exist_ok=True)
    vitro, vivo = corpus("vitro"), corpus("vivo")
    probs = check_registries()
    for t in want:
        p = write_packet(t, TARGETS[t], vitro, vivo, outdir)
        print(f"wrote {p}  ({os.path.getsize(p)/1e3:.1f} kB)")
    if probs:
        print("\nREGISTRY DRIFT — the lists of what we curated disagree:")
        for p in probs:
            print("  " + p)
        return 1
    print("\nregistries agree: review targets == build_manifest.CURATED == "
          "lit_review ADDED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
