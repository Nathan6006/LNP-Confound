# provenance/tools — acquisition and triage helpers

**These live in the repo on purpose.** They were written in a session
scratchpad on 2026-11-10 and moved here the same day, because `CLAUDE.md` §7
already states the rule twice over, from two incidents:

> **No script in this repo may READ an input from a session scratchpad.**
> Scratchpads are session-scoped and get cleared, so a tracked output that
> depends on one has no recoverable lineage.

These became load-bearing — every extraction brief calls them — while still
sitting in a directory that disappears. A curation agent reported exactly the
predicted symptom: *"no `find_repo.py` exists anywhere in the repo (`find`
returns nothing)"*, so it did its GitHub tier by hand and said so. Another
brief omitted the absolute path and the agent could not locate the tool at all.
Third instance of the class (after `lit_review/build_data.py` and
`build_manifest.py`, F-023).

| tool | what it does |
|---|---|
| `find_repo.py` | tier 3 of the intake gate: find a paper's deposited data on GitHub. Code search on DOI/PMID, then SHORT repo-name queries — a full paper title returns 0 because GitHub ANDs the terms. Relevance-filters the hits. |
| `pmc_fetch.py` | PMC full text (JATS) + supplementary files. **The working asset route is `articles/instance/<NUMERIC id>/bin/<name>`**; `articles/PMC<id>/bin/` 404s. Figures come back at 600 dpi. |
| `find_screens.py` | flags spreadsheet sheets that look like a LIBRARY SCREEN — many distinct compound-like labels, no time axis. Reports EVERY qualifying column, because one sheet routinely holds several panels side by side. |
| `triage.py` | stage-00 triage of a PDF/DOCX: pages with a text layer, distinct compound labels, numeric tokens, table-like pages. |
| `fetch_paper.py` | enumerate and download a publisher's assets. Shells out to `curl` — `urllib` is fingerprinted as a bot where curl is not. |
| `batch_gate.py` | the CROSS-PAPER half of the batch review gate: experiment-ID collisions, structures already in the corpus, one structure under two names, quality vocabulary, caveat-without-note, open decisions. |
| `check_state.py` | end-of-turn guard for the autonomous run. Fails if the run is IDLE while papers are ready. |
| `receiver.py` | localhost receiver so the browser can hand large files to disk without the bytes or the session cookies passing through a model's context. |
| `peek_xlsx.py`, `dump_sheet.py`, `repo_sweep.py` | small inspection helpers. |

## Two rules these tools exist to enforce

1. **A downloaded file is not a file until its magic bytes say so.** A size
   floor is the wrong test: a 16 KB Source Data `.xlsx` is legitimate, and
   HTML error pages have been written as `.pdf` and `.jpg` three separate
   times today (ACS 404s, PMC error pages, a 404 saved as a 60 KB "PDF").
2. **A scanner that cannot see a column reports a clean run.** `find_screens`
   hid a 12-lipid panel in column W of an 89-column sheet; `batch_gate` first
   reported "no structure column" for 4 of 10 folders and checked zero
   structures; C10 found zero structures in a folder using `label`/`SMILES`.
   Every one printed a reassuring line.
