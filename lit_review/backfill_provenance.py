#!/usr/bin/env python3
"""Backfill missing Publication_PMID and add a canonical Publication_DOI column.

Two gaps this closes:
  1. Four Bowen-Li-lab sources (FG_2026, MO_2026, YX_2026, YX_2026_SDL) carry no
     Publication_PMID - they are identifiable only from free-text Publication_link,
     which is why a PMID-keyed duplicate check silently missed them.
  2. No file carries a parsed DOI, so cross-referencing against a literature search
     requires regex-ing Publication_link every time.

Usage:
    python3 -I backfill_provenance.py --dry-run     # report only (default)
    python3 -I backfill_provenance.py --apply       # write, after .bak backup

Writes <file>.bak_<timestamp> before touching anything. Adds Publication_DOI as the
column immediately after Publication_link. Never alters Experiment_value or any
other analysis column.
"""
import argparse, json, os, re, shutil, sys, time
import pandas as pd

ROOT = "/Users/nathanliu/LNP"
TARGETS = [
    f"{ROOT}/confound/data/LNPDB_vitro.csv",
    f"{ROOT}/confound/data/LNPDB_vivo.csv",
    f"{ROOT}/LNP_ML/new_data/LNPDB_vitro_del_processed.csv",
]

# Experiment_ID -> (PMID, DOI). Identified in the 2026-10-25 literature search:
# the corpus stores only a Publication_link for these, with Publication_PMID blank.
MISSING_PMID = {
    "FG_2026":     ("41922838", "10.1038/s41563-026-02555-0"),   # Gong, Nat Mater 2026
    "MO_2026":     ("42050330", "10.1038/s41587-026-03109-0"),   # Zhou, MOLEA, Nat Biotechnol 2026
    "YX_2026":     ("41742414", "10.1016/j.cell.2026.01.012"),   # Xu, LUMI-lab, Cell 2026
    "YX_2026_SDL": ("41742414", "10.1016/j.cell.2026.01.012"),   # same paper, HBE/self-driving-lab arm
}

DOI_RE = re.compile(r"10\.\d{4,9}/[^\s\"'<>)\],;]+", re.I)


def load_pmid_doi_map():
    """PMID -> DOI, from the PubMed metadata fetched during the literature search."""
    p = os.path.join(os.path.dirname(os.path.abspath(__file__)), "corpus_titles.json")
    if not os.path.exists(p):
        return {}
    return {k: v["doi"] for k, v in json.load(open(p)).items() if v.get("doi")}


def clean_doi(s):
    if not isinstance(s, str):
        return ""
    m = DOI_RE.search(s)
    return m.group(0).rstrip(".").lower() if m else ""


def process(path, pmid2doi, apply_changes):
    name = os.path.basename(path)
    if not os.path.exists(path):
        print(f"  !! missing, skipped: {path}")
        return
    df = pd.read_csv(path, low_memory=False)
    if "Experiment_ID" not in df.columns:
        print(f"  !! no Experiment_ID, skipped: {name}")
        return
    print(f"\n{name}  ({len(df):,} rows, {len(df.columns)} cols)")

    report = {"pmid_filled": {}, "doi_from_link": 0, "doi_from_pmid": 0, "doi_blank": 0}

    # ---- 1. backfill Publication_PMID ---------------------------------------
    if "Publication_PMID" not in df.columns:
        print("   - no Publication_PMID column; not creating one")
    else:
        # cast first: the raw column is float64 wherever PMIDs are missing, so a
        # direct string assignment trips pandas' incompatible-dtype deprecation.
        df["Publication_PMID"] = df["Publication_PMID"].astype(object)
        for eid, (pmid, _doi) in MISSING_PMID.items():
            mask = (df["Experiment_ID"] == eid) & df["Publication_PMID"].isna()
            n = int(mask.sum())
            if n:
                df.loc[mask, "Publication_PMID"] = pmid
                report["pmid_filled"][eid] = (pmid, n)
        # keep the column integral: nullable Int64 avoids the float64/int64 drift
        # already documented in CLAUDE.md's 2026-07-25 changelog entry.
        df["Publication_PMID"] = pd.to_numeric(df["Publication_PMID"], errors="coerce").astype("Int64")

    # ---- 2. derive Publication_DOI -----------------------------------------
    link = df["Publication_link"] if "Publication_link" in df.columns else pd.Series([""] * len(df))
    doi = link.map(clean_doi)
    report["doi_from_link"] = int((doi != "").sum())

    # fall back to the PMID->DOI lookup where the link held no parsable DOI
    known = dict(pmid2doi)
    for _eid, (pmid, d) in MISSING_PMID.items():
        known[pmid] = d
    if "Publication_PMID" in df.columns:
        pm = df["Publication_PMID"].astype("string").str.replace(r"\.0$", "", regex=True)
        fill = pm.map(lambda x: known.get(x, "").lower() if isinstance(x, str) else "")
        need = doi == ""
        doi = doi.where(~need, fill)
        report["doi_from_pmid"] = int((need & (doi != "")).sum())
    report["doi_blank"] = int((doi == "").sum())

    if "Publication_DOI" in df.columns:
        df["Publication_DOI"] = doi
        print("   - Publication_DOI already existed; refreshed in place")
    else:
        cols = list(df.columns)
        at = cols.index("Publication_link") + 1 if "Publication_link" in cols else len(cols)
        df.insert(at, "Publication_DOI", doi)

    for eid, (pmid, n) in report["pmid_filled"].items():
        print(f"   + Publication_PMID {pmid} -> {eid} ({n:,} rows)")
    if not report["pmid_filled"]:
        print("   = no blank Publication_PMID rows matched the backfill table")
    print(f"   + Publication_DOI: {report['doi_from_link']:,} parsed from link, "
          f"{report['doi_from_pmid']:,} via PMID lookup, {report['doi_blank']:,} left blank")

    uniq = sorted(set(d for d in doi if d))
    print(f"   = {len(uniq)} distinct DOIs across "
          f"{df['Experiment_ID'].nunique()} Experiment_IDs")

    if apply_changes:
        bak = f"{path}.bak_{time.strftime('%Y%m%d-%H%M%S')}"
        shutil.copy2(path, bak)
        df.to_csv(path, index=False)
        print(f"   * written; backup at {os.path.basename(bak)}")
    else:
        print("   * dry run, nothing written")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--apply", action="store_true", help="write changes (default: dry run)")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    if a.apply and a.dry_run:
        sys.exit("pick one of --apply / --dry-run")
    pmid2doi = load_pmid_doi_map()
    print(f"PMID->DOI lookup entries: {len(pmid2doi)}")
    for t in TARGETS:
        process(t, pmid2doi, a.apply)
    print("\nRe-run the frozen pipeline (run_all.py) if any downstream stage reads these columns.")


if __name__ == "__main__":
    main()
