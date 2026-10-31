# Literature sweep provenance

Everything a reviewer needs to reproduce the candidate-source search behind
`../build_data.py` and `../data.json`. These files were rescued from the session
scratchpad on 2026-10-27; they are the only record of the exact queries used.

## Search parameters

| | |
|---|---|
Database | PubMed, via NCBI E-utilities (`esearch` / `efetch` / `esummary`) |
Sweep date | 2026-10-25 |
Queries | **20** (14 in `queries_batch1.txt`, 6 in `queries_batch2.txt`) |
Field tags | `[tiab]` (title/abstract) throughout; one query is date-restricted with `[PDAT]` |
Hits returned | 1,997 PMIDs with duplicates across queries (`pmids_all.txt`) |
Unique PMIDs | 1,341 (`pmids_uniq.txt`) |
Records parsed | 1,996 (`recs3.json`, kept in scratch — regenerate with `parse2.py`) |
Corpus at sweep time | `corpus_pmids_at_sweep.json`, the Experiment_ID → PMID map used for duplicate checking |

`meta.swept` in `data.json` reports 1,996, which is the count of successfully parsed
records, not the 1,341 unique PMIDs. The gap is because `pmids_all.txt` retains
cross-query duplicates; the two numbers measure different things and both are kept here
so the discrepancy is visible rather than hidden.

## Files

- `queries_batch1.txt` — the 14 queries of the first sweep, one per line, verbatim
- `queries_batch2.txt` — the 6 follow-up queries, verbatim
- `pmids_all.txt` / `pmids_uniq.txt` — raw and de-duplicated hit lists
- `corpus_pmids_at_sweep.json` — what was already curated when the sweep ran
- `parse.py`, `parse2.py` — PubMed XML parsers. `parse2.py` is the correct one: it reads
  DOIs only from `./PubmedData/ArticleIdList`. `parse.py` walked every `ArticleId` in the
  record, which reached into `ReferenceList` and returned the last *reference's* DOI. It
  is kept only so the bug is documented.

## Known limitations of the sweep, stated plainly

1. **PubMed only.** No Web of Science, Scopus, Google Scholar, or chemistry-specific
   indexes. ChemRxiv in particular is not in PubMed or Europe PMC, so chemistry preprints
   are under-represented.
2. **`[tiab]` only.** A paper whose library is described solely in the full text or SI,
   and never in the title or abstract, is invisible to these queries. This is a real
   miss mode: Qiaobing Xu's libraries routinely omit the size from the abstract, which is
   why that lab surfaced only through citation mining of an external corpus (LipidDB),
   not through these queries.
3. **Library sizes were auto-extracted from abstracts first**, and that regex produced
   false positives from reagent names (`ALC-0315` → 315, `SM-102` → 102, `PEG2000` →
   2000, `C12-200` → 200, `C14-306` → 306, `YK-407` → 407). Every size reported in
   `data.json` was subsequently re-derived by hand or from full text; see the `SIZES`
   table in `build_data.py` for the per-row provenance (`abstract` / `fulltext` /
   `download` / unknown).
4. **English-language bias** inherent to the index.

## Reproducing

```bash
# per query, collect PMIDs
esearch -db pubmed -query "<line from queries_batch*.txt>" | efetch -format uid
```
or through the REST API as the sweep did:
`https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&retmax=500&term=<urlencoded query>`

Counts will drift upward from 2026-10-25 as new records are indexed.
