#!/usr/bin/env python3
"""Build a per-experiment manifest of the LNP corpus (read-only on the corpus)."""
import json
import sys
import pandas as pd

DATA = "/Users/nathanliu/LNP/confound/data"
OUT = f"{DATA}/experiment_manifest.csv"
ESUM = "/private/tmp/claude-501/-Users-nathanliu-Documents-LNP-confound/ebd1fdd9-2e2e-4284-8eb1-27c940f1a93a/scratchpad/esum.json"

# ---- authoritative "curated by us" map: ADDED dict in lit_review/build_data.py
# cross-checked against the folders present in confound/curation/
CURATED = {
    "ZT_2026":   "curation/siegwart_42634137",
    "LL_2026":   "curation/lin_aromatic_2026",
    "LL_2026_2": "curation/lin_aromatic_2026",
    "WT_2026":   "curation/wht_2026_falcon",
    "XD_2026":   "curation/xd_2026_lipidai",
    "QL_2026":   "curation/ql_2026",
}

# PMIDs externally resolved / corrected against NCBI E-utilities this session.
# key: Experiment_ID -> (field, value, reason)
EXTERNAL = {
    # corpus PMID 29786478 is a 2018 Int J Stroke consensus statement
    # (doi 10.1177/1747493018778713) - an upstream mis-key. esearch on the
    # corpus DOI 10.1016/j.bioactmat.2024.05.012[doi] -> 38779592
    # (Bae SH, Bioact Mater 2024), which matches the SB_2024 initials/year.
    "SB_2024": ("PMID", "38779592",
                "PMID corrected externally: data PMID 29786478 is an unrelated "
                "2018 Int J Stroke paper; esearch on the data DOI gave 38779592 "
                "(Bae SH, Bioact Mater 2024)"),
    # corpus DOI carries a bioRxiv version suffix 'v1'; canonical DOI per
    # esummary on PMID 40060499 is without it.
    "XH_2025": ("DOI", "10.1101/2025.02.25.640222",
                "DOI normalised externally: data value '10.1101/2025.02.25.640222v1' "
                "carries a bioRxiv version suffix; esummary on PMID 40060499 "
                "gives the canonical DOI"),
}

# Folders present in confound/curation/ that are NOT in the ADDED dict: these
# rows came from LNPDB upstream, but local re-extraction / validation material
# exists. curated_by_us stays 'no' (nothing of ours was merged for them).
CURATION_ONLY = {
    "SL_2021": "these rows are upstream LNPDB, but curation/liu_iphos/ holds an "
               "independent re-extraction of the same paper (572 rows, "
               "Experiment_ID SL_2021_reextract) that was NOT merged into the corpus",
    "FG_2026": "these rows are upstream LNPDB, but curation/amino_lipid_fg_2026/ "
               "holds the paper PDF + SI and a QC file for it (no rows merged)",
}

PIPE_COLS = ["Model_type", "Model_target", "Route_of_administration",
             "Cargo", "Cargo_type", "Experiment_method", "HL_name"]


def pipe(series):
    vals = sorted({str(v).strip() for v in series
                   if str(v).strip() not in ("", "nan", "NA", "None")})
    return "|".join(vals)


def one(series, label, eid):
    """Single distinct non-blank value, or '' - warn if several."""
    vals = sorted({str(v).strip() for v in series
                   if str(v).strip() not in ("", "nan", "NA", "None")})
    if len(vals) > 1:
        print(f"  WARN {eid}: {len(vals)} distinct {label}: {vals}", file=sys.stderr)
    return vals[0] if vals else ""


# ---- NCBI esummary cache (PMID -> doi, journal, year, first author)
ncbi = {}
res = json.load(open(ESUM))["result"]
for u in res["uids"]:
    r = res[u]
    ncbi[u] = dict(
        doi="".join(a["value"] for a in r.get("articleids", []) if a["idtype"] == "doi"),
        source=r.get("source", ""), pubdate=r.get("pubdate", ""),
        author=r.get("sortfirstauthor", ""), title=r.get("title", ""))
# the externally resolved SB_2024 record, verified in-session
ncbi["38779592"] = dict(doi="10.1016/j.bioactmat.2024.05.012", source="Bioact Mater",
                        pubdate="2024 Aug", author="Bae SH",
                        title="A lipid nanoparticle platform incorporating trehalose glycolipid")

frames = {c: pd.read_csv(f"{DATA}/data_{c}.csv", low_memory=False)
          for c in ("vitro", "vivo")}

# experiments present in both corpora
both = set(frames["vitro"].Experiment_ID.unique()) & set(frames["vivo"].Experiment_ID.unique())

rows = []
for corpus in ("vitro", "vivo"):
    df = frames[corpus]
    for eid, g in df.groupby("Experiment_ID", sort=True):
        notes = []

        # ---------- identifiers ----------
        pmid_data = one(g.Publication_PMID.astype(str).str.replace(r"\.0$", "", regex=True),
                        "PMID", eid)
        doi_data = one(g.Publication_DOI, "DOI", eid)
        link = one(g.Publication_link, "link", eid)
        pmid, doi = pmid_data, doi_data
        pmid_src = "data" if pmid_data else "missing"
        doi_src = "data" if doi_data else "missing"

        ext = EXTERNAL.get(eid)
        if ext:
            field, value, reason = ext
            if field == "PMID":
                pmid, pmid_src = value, "resolved externally"
            else:
                doi, doi_src = value, "resolved externally"
            notes.append(reason)

        # cross-check the surviving pair against NCBI; lead the note with a
        # human-searchable citation (first author / journal / year)
        rec = ncbi.get(pmid)
        if rec is None:
            notes.append(f"PMID {pmid} not checked against NCBI")
        else:
            cite = f"{rec['author']}, {rec['source']} {rec['pubdate'].split()[0]}"
            notes.insert(0, cite)
            if rec["doi"].lower() != doi.lower():
                notes.append(f"PMID/DOI cross-check FAILED: NCBI reports doi "
                             f"{rec['doi']} for PMID {pmid} (data DOI {doi_data})")
            elif ext:
                notes.append("PMID/DOI pair now consistent with NCBI esummary")
            else:
                notes.append("PMID and DOI both came from the data and were verified "
                             "against NCBI esummary (PMID->DOI matches)")

        # link sanity: does the DOI suffix appear in the link? (publisher URLs
        # routinely drop the 10.xxxx/ registrant prefix, so compare the suffix)
        suffix = doi.split("/", 1)[1].lower() if "/" in doi else doi.lower()
        if doi and link.startswith("http") and suffix not in link.lower():
            if "/pii/" in link.lower():
                pass  # ScienceDirect PII links never carry the DOI; not a mismatch
            else:
                notes.append("WARNING: Publication_link points at a different article "
                             "than the DOI - trust the PMID/DOI, not the link")
        elif doi and not link.startswith("http"):
            notes.append(f"Publication_link is not a URL ({link!r})")

        # curation/ cross-check: material exists for these but was NOT merged,
        # so curated_by_us stays 'no' per the ADDED dict in build_data.py
        extra = CURATION_ONLY.get(eid)
        if extra:
            notes.append(extra)

        # ---------- composition / stats ----------
        v = pd.to_numeric(g.Experiment_value, errors="coerce")
        mw = pd.to_numeric(g["Molecular.Weight"], errors="coerce")

        if eid in both:
            notes.append(f"{eid} appears in BOTH corpora (vitro and vivo); this row is "
                         f"the {corpus} arm only - there is a second row for the other arm")

        rows.append({
            "Experiment_ID": eid,
            "corpus": corpus,
            "n_rows": len(g),
            "n_unique_lipids": g.IL_SMILES.dropna().nunique(),
            "Publication_PMID": pmid,
            "Publication_DOI": doi,
            "Publication_link": link,
            **{c: pipe(g[c]) for c in PIPE_COLS},
            "n_distinct_formulations": g.Formulation_ID.dropna().nunique(),
            "value_min": round(v.min(), 4), "value_max": round(v.max(), 4),
            "value_mean": round(v.mean(), 4) + 0.0,
            "value_sd": round(v.std(ddof=1), 4) if v.notna().sum() > 1 else "",
            "MolWt_min": round(mw.min(), 2), "MolWt_max": round(mw.max(), 2),
            "curated_by_us": "yes" if eid in CURATED else "no",
            "curation_folder": CURATED.get(eid, ""),
            "notes": "; ".join(notes),
        })

COLS = ["Experiment_ID", "corpus", "n_rows", "n_unique_lipids",
        "Publication_PMID", "Publication_DOI", "Publication_link",
        *PIPE_COLS[:-1], "HL_name", "n_distinct_formulations",
        "value_min", "value_max", "value_mean", "value_sd",
        "MolWt_min", "MolWt_max", "curated_by_us", "curation_folder", "notes"]
# HL_name already in PIPE_COLS[:-1]? no - PIPE_COLS ends with HL_name, so rebuild cleanly
COLS = ["Experiment_ID", "corpus", "n_rows", "n_unique_lipids",
        "Publication_PMID", "Publication_DOI", "Publication_link",
        "Model_type", "Model_target", "Route_of_administration",
        "Cargo", "Cargo_type", "Experiment_method", "HL_name",
        "n_distinct_formulations", "value_min", "value_max", "value_mean",
        "value_sd", "MolWt_min", "MolWt_max", "curated_by_us",
        "curation_folder", "notes"]

m = pd.DataFrame(rows)[COLS].sort_values(["corpus", "Experiment_ID"],
                                         kind="mergesort").reset_index(drop=True)
m.to_csv(OUT, index=False)

# ---------------------------------------------------------------- verification
n_vitro_exp = frames["vitro"].Experiment_ID.nunique()
n_vivo_exp = frames["vivo"].Experiment_ID.nunique()
print(f"distinct Experiment_ID: vitro={n_vitro_exp}  vivo={n_vivo_exp}  "
      f"expected rows={n_vitro_exp + n_vivo_exp}  manifest rows={len(m)}")
assert len(m) == n_vitro_exp + n_vivo_exp, "ROW COUNT MISMATCH"

# The assertion that MEANS something is that the manifest's per-experiment
# counts add up to the corpus files themselves, so the reference is the live
# file length, not a literal. An earlier version hardcoded 25,997 / 1,768, which
# is a number that goes stale on every merge: it fired as a DISCREPANCY the
# moment QL_2026 landed, reporting a correct manifest as broken.
#
# BASELINE is kept only as a tripwire for a corpus change nobody intended. It is
# "the size at the last review", so a mismatch prints and does NOT fail - update
# it deliberately when a merge is signed off.
BASELINE = {"vitro": 28885, "vivo": 1780}   # reviewed 2026-10-30
for corpus in ("vitro", "vivo"):
    got = int(m.loc[m.corpus == corpus, "n_rows"].sum())
    live = len(frames[corpus])
    print(f"sum n_rows {corpus}: {got} (corpus file has {live}) "
          f"{'OK' if got == live else '*** DISCREPANCY ***'}")
    assert got == live, (f"manifest n_rows sum {got} != {corpus} file rows "
                         f"{live} - the manifest does not describe the corpus")
    if live != BASELINE[corpus]:
        print(f"  NOTE  {corpus} is {live} rows, baseline {BASELINE[corpus]} "
              f"({live - BASELINE[corpus]:+d}). Intended? Update BASELINE.")

print(f"with PMID: {(m.Publication_PMID != '').sum()}/{len(m)}   "
      f"with DOI: {(m.Publication_DOI != '').sum()}/{len(m)}")
miss = m[(m.Publication_PMID == "") | (m.Publication_DOI == "")]
print("missing:", "none" if miss.empty else
      list(zip(miss.Experiment_ID, miss.corpus)))
print(f"curated_by_us=yes: {(m.curated_by_us=='yes').sum()} rows -> "
      f"{sorted(set(m.loc[m.curated_by_us=='yes','Experiment_ID']))}")
print(f"both-corpora experiments ({len(both)}): {sorted(both)}")
print(f"wrote {OUT}")
