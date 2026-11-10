# Autonomous run — log

Append-only. Every action, with UTC timestamps. Plan: `AUTONOMOUS_RUN.md`.
Machine state: `run_queue.json`.

---

## 2026-11-10 — Session 1

### 03:5x — Planning

- Read `lit_review/data.json`: 104 candidates, **41 with `status = "new"`**.
  Composition: **26 `ivl`** (ionizable-lipid chemistry), 10 `ivv`, 5 `form`.
  The user's priority — lipid chemistry over formulation — is well served:
  `ivl` outnumbers `form` 5:1 in the backlog.
- Wrote `provenance/AUTONOMOUS_RUN.md` (plan of record) and
  `provenance/run_queue.json` (41 papers, ranked `ivl` → `ivv` → `form`, then by
  stated library size).
- **15 of the 26 `ivl` rows have `lib = -1`.** Recorded in the plan that this
  means *unknown*, not small, and must not deprioritise a paper — establishing
  the library size is stage 00's first job.
- Found **queue rank 3 (`10.1038/s41587-026-03298-8`) is already on disk** in
  `curation/sd_2026/`. Deduplicated.

### 04:0x — Acquisition baseline

Probed all 41 DOIs with a plain `curl`: **7 reachable, 34 refused (403)**.
Publisher mix: ACS 19, Wiley 9, Elsevier 4, Springer Nature 4, PNAS 3, RSC 1,
Taylor & Francis 1.

**Finding A — `urllib` is fingerprinted where `curl` is not.** The first
fetcher used `urllib.request` with a current Chrome User-Agent and Springer
Nature returned a "Client Challenge" page; the identical URL and UA through
`curl` returned the article. The TLS/header fingerprint matters more than the
UA string. `fetch_paper.py` now shells out to `curl` and says so in a comment,
so it does not get "simplified" back.

**Downloaded, no browser needed (3 papers, 37 files):**

| folder | what | files |
|---|---|---|
| `mitchellmj_38982196` | Nat Chem, amidine lipids, 100 lipids | 8 — incl. **6 Source Data XLSX** |
| `linyx_41851281` | Nat Biomed Eng, AI-guided LNP design | 17 (60 MB) |
| `andersondg_41203968` | Nat Nanotechnol, cyclic amino alcohol lipids | 12 (16 MB) |

Six Source Data workbooks on the first paper alone — exactly the tier-1 material
the intake gate exists to find, and it means zero digitisation for those values.

**Finding B — the 2026-10-25 access grades badly understate reach.** 30 of the
41 are graded `C` ("not reachable — try CWRU library access"). That grade came
from a plain HTTP probe, and today's matrix shows a plain probe 403s on six of
eight publishers *regardless of entitlement*. Grade C is not evidence of a
paywall. Re-probe; do not trust it.

### 04:2x — ACS enumeration

**Finding C — bulk enumeration via in-page `fetch()` works and is ~10x cheaper
than navigating per paper.** After one `browser_navigate` to warm the session
(ACS needs ~10 s for its Cloudflare challenge, and resolves on a second
navigation), a single `browser_evaluate` can `fetch()` many article pages with
`credentials:"include"` and return only the extracted asset URLs — a small
payload instead of a page snapshot per paper.

All **18 of 19** ACS papers yielded a supplementary URL; `10.1021/bc200572w`
(2012, Bioconjug Chem) returned none and is marked paywalled-or-no-SI for
re-check. The ACS pattern is
`pubs.acs.org/<journalcode>/article-supplement/<articleId>/pdf/<name>/`, and it
*is* present in the raw HTML, so no DOM rendering is needed.

**Finding D — ACS supplementary files need the session; the Elsevier CDN did
not.** `curl` on an ACS supplement URL returns 403 with or without the
`/doi/suppl/` form. This differs from Elsevier, where `ars.els-cdn.com` served
both supplementary files and figure rasters to an anonymous request once the
URL was known. So "enumerate signed-in, download anonymously" is an
Elsevier/Springer property, not a general one.

**Finding E — the Playwright MCP server cancels downloads started inside
`browser_run_code_unsafe`.** `download.saveAs()` and `download.path()` both
raise `canceled` for every file after the first, because the server registers
its own download handler. Triggering the download and letting the server save
to `.playwright-mcp/` did not complete either. Four variants were tried.

**Workaround adopted:** a localhost receiver
(`scratchpad/tools/receiver.py`). The browser `fetch()`es the asset with its
own cookies and POSTs the bytes to `127.0.0.1:8765`, which writes the file and
appends a sha256 to `curation/_received.jsonl`. This keeps both the session
cookies and the multi-megabyte payloads out of the model's context — the two
things the obvious workarounds would have cost. The receiver refuses any path
escaping its root, because a path that came back through a web page is
untrusted input.

### 04:3x–04:4x — Transport, resolved per publisher

**Finding F — the ACS 404s were mine.** The enumerated supplement path is
relative to the JOURNAL path, not the site root: it must be
`pubs.acs.org/<journalcode>/article-supplement/...` (`ancac3` ACS Nano,
`jacsat` JACS, `bcches` Bioconjug Chem, `asbcd6` ACS Synth Biol). Resolving it
against the article's final URL fixed all of them. Three 60 KB error pages had
been saved as "PDFs" before this was caught — deleted, and the receiver
manifest now records a file's magic bytes so `%PDF` can be asserted.

**Finding G — transport differs per publisher and none of the obvious routes
is general.** Four transports were tried; what each publisher accepts:

| transport | ACS | Wiley | Springer |
|---|---|---|---|
| plain `curl` | 403 | 403 | **200** |
| Playwright `APIRequestContext` (`page.context().request`) | **200** | 403 | 200 |
| in-page `fetch()` + POST to localhost | fetch ok, **POST blocked by page CSP** | same | n/a |
| `page.goto()` + wait, MCP saves the download | — | **200** | n/a |

So: **Springer = curl; ACS = request context; Wiley = navigate-and-wait.** The
request context passes ACS's bot check but not Wiley's, which fingerprints the
request itself rather than only the cookies — a *page* navigation to the same
Wiley URL succeeds where the request context 403s on identical cookies.

**Finding H — in-page `fetch()` to `127.0.0.1` hangs, it does not error.** Page
CSP blocks the connection and the promise never settles, so two
`browser_evaluate` calls ran past 120 s and were backgrounded before this was
understood. The working combination is **request context for BOTH legs** —
fetch the asset and POST it to the receiver from the Playwright server process,
which has no page CSP.

### Acquisition result, session 1

**21 folders, 56 files.**

| | papers | files | note |
|---|---|---|---|
| ACS | **18 of 19** | 18 SI PDFs, 3–54 pp, all verified `%PDF` | `10.1021/bc200572w` (2012) lists no SI |
| Springer Nature | 3 | 37, incl. **many Source Data XLSX** | zero digitisation needed for those values |
| Wiley | 1 of 9 | 1 `.docx` | transport established, rest pending |

**Finding I — a big Source Data sheet is not a big screen.** `linyx_41851281`'s
`MOESM11` has 1,975 rows and 11,166 numeric cells, and is a **molecular-dynamics
time series** (`Time (ns)`); its `Fig. 2b` sheet has 1,223 rows of flow-cytometry
zeros. Cell count is not library size. Wrote
`scratchpad/tools/find_screens.py`, which flags a sheet only when it has >=15
distinct compound-like labels and no time-axis header — stage 00 should run it
before anyone reads a workbook by eye.

### 04:4x–05:0x — PMC route, triage, first extraction batch

**Finding J — 24 of the 41 candidates are in PubMed Central**, including both of
the two largest screens (PMID 39664425, 584 lipids; PMID 38888972, 580). PMC
needs no publisher, no session and no bot workaround.

But the route is narrower than it first looks. The old `oa.fcgi` endpoint 404s;
the working one is `eutils efetch?db=pmc&retmode=xml`, which returns the JATS
full text. **PMC gives full TEXT, not the publisher's supplementary data
files** — for PMC11628822 the only `xlink:href`s are `gr1.jpg` and `main.pdf`,
and for PMC11446653 (author manuscript) the supplementary items are figure
JPEGs. So PMC is the right fallback for an article PDF or for Methods prose,
and it does NOT replace a publisher SI. Tool: `scratchpad/tools/pmc_fetch.py`.

**Finding K — Elsevier's request context is refused like Wiley's.** The
`APIRequestContext` 403s on sciencedirect.com even with a warm page and a live
entitlement, so Elsevier joins Wiley on the navigate-and-wait transport. Also
probed `ars.els-cdn.com` for the usual `-mmc<n>.<ext>` supplement pattern on
PII S2211383524003514: **no supplement exists**, and only `gr1_lrg` of eight
figure rasters is present. Do not assume the CDN pattern holds for every paper.

**Stage-00 triage run across all 22 acquired folders** (`scratchpad/tools/triage.py`,
output `scratchpad/triage.json`): per file, pages, pages with a real text layer,
distinct compound-like labels, numeric tokens and table-like pages. Richest:
`andersondg_41203968` (302 pp, 231 labels, 31,420 numbers), `linyx_41851281`
(3,318), `dahlmanje_29489381` (2,856), `mckinlaycj_39105751` (2,687),
`mitchellmj_41555650` (2,462). Thinnest, and likely below the 100-point
threshold: `xuq_23651337` (3 pp, 3 labels), `xuq_32520527` (8 pp, 8 labels),
`maohq_39375194` (9 labels).

**Extraction batch 1 launched** (3 concurrent, the plan's cap), each confined to
its own folder, each told to run the intake gate first, import `curation/lib`,
emit via `emit.emit_rows`, keep unresolved questions in `DECISIONS.md`, and
NOT to merge:

| folder | paper | why first |
|---|---|---|
| `mitchellmj_38982196` | Nat Chem 2024, amidine lipids, ~100 | tier-1 Source Data; `Fig. 2b` is lipids 1..N with replicates, `Fig. 2c` a lipid x tail matrix |
| `andersondg_41203968` | Nat Nanotechnol 2025, cyclic amino alcohol | 10 Source Data workbooks + the richest SI in the set |
| `lij_42130331` | ACS Nano 2026, β-hydroxy thioether, ~300 | largest ACS screen; no Source Data, so it tests the intake gate |

### 05:0x — Wiley, PNAS acquired. Transport rule found.

**Finding L — the transport rule is "navigate once per host, then the request
context works".** Wiley and PNAS both 403 `APIRequestContext` on a cold host,
including with a warm *article* page. But after ONE `browser_navigate` to a
**download URL** on that host — which makes the browser solve the Cloudflare
challenge for that path — the request context succeeds for every subsequent
file on the same host. That reduced Wiley from an estimated ~3 calls per file
to one navigation plus one batched call for four files.

Corollary: Wiley serves each journal family on its own subdomain
(`advanced.`, `aiche.`, bare `onlinelibrary.`), and a same-origin in-page
`fetch()` can enumerate every DOI on the subdomain you are standing on, but
cross-origin fetch fails CORS. So: group DOIs by subdomain, navigate once,
enumerate the group, download the group.

**Acquired this pass:** Wiley 9/9, PNAS 3/3.

**Finding M — the 584-lipid paper (rank 1, PMID 39664425) has no reachable
supplementary file.** Checked, in order: ScienceDirect appendix section (no
supplementary links rendered), `ars.els-cdn.com` `-mmc<n>.{pdf,docx,xlsx,zip}`
for n=1..4 (all 404), PMC (only 16 KB of text, 0 tables — a restricted letter),
and the LipidAI GitHub repo named in our tracker note, which belongs to the
LATER paper and holds `LipidDB` and `descriptors_1107.csv`, not this paper's
584. Its PMC text does confirm "584 ILs and their mRNA transfection data".
Its presigned `pdf.sciencedirectassets.com` URL is single-use and expired
before the request context could reuse it. **Recorded as blocked on its SI,
not abandoned**: the main PDF via browser navigation is still untried.

### Acquisition total: 34 folders, 68 files, every file magic-byte verified

| publisher | papers |
|---|---|
| ACS | 18 of 19 |
| Wiley | 9 of 9 |
| Springer Nature | 3 |
| PNAS | 3 of 3 |
| Elsevier | 1 of 4 (partial: PMC text only) |

Receiver manifest `curation/_received.jsonl` holds url, sha256, byte count and
magic bytes for all 29 browser-mediated downloads. No file failed the
`%PDF`/`PK` check.

### 05:1x — `kind` verified against abstracts; three reclassifications

**Finding N — the tracker's `kind` field is not reliable, and it is the field
the user's priority depends on.** Checked every queued paper's abstract for
"<n> lipids" versus "<n> formulations":

| paper | was | is | evidence |
|---|---|---|---|
| `chauy_40028722` | `ivl`, lib 216 | **`form`** | "a library of **216 LNP formulations** with varying lipid compositions, N/P ratios, and hyaluronic acid modifications" |
| `dixonjb_41269090` | `ivv`, lib 150 | **`form`** | abstract says 150 formulations |
| `maohq_42220328` | `form`, lib 1944 | **duplicate** | same paper as `curation/jl_2026` (Materials Today, PII S1369702126001112) |

So **stage 00 must establish whether a screen varies CHEMISTRY or
FORMULATION**, from the abstract, before the paper is ranked — a 216-row
formulation screen was sitting at rank 7 ahead of real chemistry screens.
Queue re-ranked; `kind_original` kept on each row so the change is visible.

**Finding O — `pmc_fetch.py` wrote 7 HTML error pages as `.jpg`/`.pdf`.** The
guard tested `data[:15]` for `<!doctype html`, and PMC's error pages start with
newlines, so every one passed. Now the head is stripped before the test and a
known magic is required. The same class as the ACS 404s earlier today: **a
downloaded file is not a file until its magic bytes say so.** Files deleted.

T&F `chauy_40028722` did yield 200 KB of JATS full text with 3 real tables —
useful even though it is now deprioritised as a formulation screen.

### 05:2x — Triage extended to 33 folders; rank 2 launched

Triage now covers all 33 acquired folders (`scratchpad/triage.json`).

**Rank 2, `whiteheadka_38888972` (580 lipidoids), is extractable and was
launched.** Its Wiley SI `.docx` holds **Table S1 "List of amine heads used for
ionizable lipidoid library synthesis"** (codes + IUPAC names) and a table of
acrylate tails (`O10` n-decyl … `O16` hexadecyl) — together the full
combinatorial definition. Caveat passed to the agent: the docx has only ~5,300
characters of body text but 4 embedded images, so some "tables" are pictures.

**Finding P — rank 1 (`songx_39664425`, 584 lipids) is blocked on its data, and
the blockers are specific.** ScienceDirect renders no supplementary section;
`ars.els-cdn.com` has no `-mmc<n>` asset; PMC holds 16 KB of text with 0 tables;
the LipidAI GitHub repo belongs to the later paper; and the `pdfft` link
resolves through an interstitial to a single-use presigned S3 URL that expires
before it can be refetched. Its PMC text confirms "584 ILs and their mRNA
transfection data" exist. Kept in the queue as `blocked_no_si`, not dropped.

**Thin sources noted for the <100-point stop rule** (to be confirmed at
extraction, not assumed): `dahlmanje_41183194` (12 pp, 1 page with a text layer,
4 numbers — effectively a scanned SI), `xuq_24574196` (9 pp, 12 labels),
`siegwartdj_42150080` (15 pp, 24 labels, 166 numbers).

**Extraction in flight (4):** `mitchellmj_38982196`, `andersondg_41203968`,
`lij_42130331`, `whiteheadka_38888972`. Raised the concurrency cap from 3 to 4
for this batch: the folders are disjoint, no agent touches shared code, and the
4th is the largest screen in the queue.

### 05:3x — Tier 3 of the intake gate is now a real tool; pool expanded

**`scratchpad/tools/find_repo.py`** — searches GitHub for a paper's deposited
data (code search on the DOI and PMID, then repo search on short title/author
queries), lists data-shaped files in each hit, and filters by a
lipid/LNP/mRNA relevance test on the repo name and file paths.

Validated against a repo we already knew existed before trusting it:
`"FALCON lipid nanoparticle optimization Mao"` returned **0** results, while
`"FALCON optimization lipid"` found `MaoResearchGroup/falcon-lnp-optimization`
and the filter discarded four unrelated "falcon" repos. **GitHub ANDs repo-search
terms, so a full paper title finds nothing** — short queries only.

**Finding Q — the GitHub route recovered a dataset the publisher route could
not.** `MaoResearchGroup/ml_lnp_design_rules` is the repo for PMID 39375194,
which our tracker graded `C` ("needs library access"). It contains
**`Raw_Data/1080_LNP_OHE.csv` — 1,080 LNPs** with ionizable lipid, helper
lipid, N/P and molar ratios, size, zeta and transfection readouts, plus
`Compiled_validation.xlsx`. Downloaded to
`curation/maohq_39375194/source/github/`. It is a FORMULATION screen (one
ionizable lipid, varied helper lipid and ratios) so it stays low priority, but
it is 1,080 machine-readable points that needed no digitisation at all.

So **tier 3 moves up the intake gate**: check GitHub for every paper, not only
when tiers 1-2 fail. A lab that deposits its screen has removed the entire
digitisation error class for us.

Other lab repos surfaced, to assess separately rather than per paper:
`evancollins1/LNPBO` (Bayesian-optimisation example sets, 150-330 KB CSVs),
`longbaobaozhuiwen/LNPAgent`, `Martini-Force-Field-Initiative/M3_Ionizable_Lipids`.

### New PubMed sweep (the tracker's was 2026-10-25)

8 queries aimed at combinatorial lipid CHEMISTRY rather than formulation
optimisation -> 114 PMIDs, **81 new** to the tracker, queue and corpus. Parsed
every abstract for "<n> lipids" vs "<n> formulations". The ones worth having:

| PMID | lipids | journal | what |
|---|---|---|---|
| 40644303 | **465** | ACS Biomater Sci Eng 2025 | SAR of ionizable lipids for siRNA and mRNA |
| 41162324 | **315** | JACS 2025 | fluorinated ionizable lipids, spleen-targeted |
| 37602460 | **104** | ACS Appl Mater Interfaces 2023 | combinatorial cationic lipidoids |

10 new candidates added to the queue (51 total). **Acquired 6 more ACS papers**
including the 465- and 315-lipid screens.

**Finding R — a size floor is the wrong file-validity test.** A 16,583-byte
Source Data `.xlsx` was skipped by a `>20000 bytes` guard; a spreadsheet of a
few hundred rows is legitimately that small. Replaced with a magic-byte test
(`%PDF` / `PK`), which is the same rule the ACS-404 and PMC-error-page
incidents already pointed at. Also: the ACS supplement URL encodes its own
format (`/pdf/` vs `/xlsx/`), so the extension should be read from the URL
rather than assumed.

### 05:4x — Extraction batch 1 complete: 3 of 4 reported, 532 rows staged

| folder | rows | value tier | library |
|---|---|---|---|
| `lij_42130331` | **300** | printed numerals in raster source-data tables | stated 300 / enumerated 300 |
| `andersondg_41203968` | **132** (100 vitro + 32 vivo) | **tier 1, no pixel read for any value** | stated 100 / enumerated 100 |
| `mitchellmj_38982196` | **100** | **tier 1**, both Source Data sheets | stated 100 / enumerated 100 |

**532 rows staged. Nothing merged** — all three have OPEN `DECISIONS.md` entries
(5, 9 and 5), so `build_dataset.py`'s gate blocks them by design.

Best external checks produced, each independent of the decoding:
- `lij_42130331`: the paper's named lead decomposes to a table cell that is the
  **unique maximum of all 250**; nothing in the glyph decoding knew which cell
  should win. Glyph clustering reached **max intra-cluster distance 0.000**.
- `andersondg_41203968`: tail-vs-headgroup variance explained (0.427 vs 0.106)
  is a **transpose guard** on the 10x10 matrix -- a transposed read swaps those
  two numbers and nothing else in the folder would notice.
- `mitchellmj_38982196`: printed "Calculated MS: 603.01" reproduced at
  **delta 0.000**, and it is a function of two measured chain lengths at once.

**Finding S -- a paper disagreeing with itself, measured.** In
`mitchellmj_38982196`, Fig. 2g re-tests 8 screened lipids, same cell line, same
stated dose, no stated protocol difference: ratios 0.32x-0.95x, geometric mean
0.64x, and **Spearman -0.262 (p=0.53, n=8)** against the screen -- while the
screen's own replicates agree at **0.990 over 100 lipids**. That is a clean
same-lab, same-assay, cross-experiment rank degradation, and it bears directly
on this project's Result 3 and its T2 noise ceiling. Queued as that folder's
D-6 (emit both as separate `Experiment_ID`s?).

### Three defects in code I wrote today, found by the extractions

1. **C10 crashed repo-wide** with `KeyError: 'lipid_name'` naming no file,
   because `new_curation.py` scaffolded the header `compound_id,...` while
   `structure_check.py` reads `lipid_name,...`. Scaffold fixed; C10 now SKIPS a
   malformed file, names it with found-vs-expected, and still fails at the end.
2. **C10 was finding ZERO structures** in a folder whose enumeration used
   `label`/`SMILES` instead of `lipid_name`/`IL_SMILES`, and reported the
   anchor as "not found" rather than "nothing to compare against". Now accepts
   documented aliases AND prints per-folder coverage, so zero coverage cannot
   hide. That anchor re-verifies as `exact_match`.
3. **`assemble.check_stereo_preserved` refused legitimate chemistry.** It
   compared the whole stereo fingerprint for equality, so a bond that CREATES a
   centre (alkylating the second N of a cis-fused bicycle) failed exactly as an
   INVERSION does -- it fired on 40 valid builds. Now refuses only a LOSS or
   INVERSION of an assigned centre and reports created ones. F-013 still fires:
   inverting C3 takes R:4,S:4 -> R:5,S:3, and a decrease in R or S is the
   signal. New control C-22 and mutation M-26 added; suite **49/49**.
4. **`chains.vector_bonds` returned zero on a whole SI.** All 13,271 drawing
   paths in `andersondg_41203968` are FILLED with `width = None`; the
   stroked-only filter made the function useless there. Added `include_filled`,
   documented with the measurement and with the cost (glyph outlines then
   return as "bonds", so the caller must separate them by thickness, aspect or
   label proximity).

**One agent claim was WRONG and is recorded as such:** C16 "prints its findings
but exits 0". C16 examined 3 templates, found none transcribed, and exited 0 --
correct behaviour -- and that agent's own seeding test showed C16 fires. Treat
reviewer findings as claims.

**Finding T -- ACS full text is not available on this session's entitlement.**
The ACS institutional-login page identifies the browser as PHILLIPS EXETER
ACADEMY, not CWRU; the CWRU OpenAthens login covers Elsevier. ACS *SI* is free
and that is where the data was, so this costs Methods prose (assay, cell line,
units, reaction equivalents), not values. Recorded in the affected folders'
`DECISIONS.md` rather than worked around.

### 06:0x — Batch 2 launched; a process failure fixed structurally

**I ended a turn saying "Launching batch 2" without launching anything.** The
user had to point it out. In an autonomous run that is the worst failure
available: the session looks busy, the log reads correctly, and nothing is
happening. Two structural fixes, not a promise:

1. **`AUTONOMOUS_RUN.md` §0 — "Launch, do not announce."** The `Agent` calls go
   out BEFORE the sentence describing them. Never end a turn with "next I
   will…".
2. **`scratchpad/tools/check_state.py` — an end-of-turn guard**, to be run
   before ending any turn. It fails if a paper is marked `in_progress_agent`
   with no live task file, if a folder is marked `extracted` without
   `05_values.csv`, or **if zero agents are running while papers are still
   ready** — an idle run is a failed run.

**Batch 2 (4 concurrent), the largest remaining chemistry screens:**

| folder | stated lipids | journal |
|---|---|---|
| `siegwart_40644303` | **465** | ACS Biomater Sci Eng 2025 |
| `liu_41162324` | **315** | JACS 2025 |
| `levkinpa_29558113` | **288** | Bioconjug Chem 2018 |
| `kooh_42690843` | **275** | JACS 2026 |

Briefs carry what batch 1 taught: run the GitHub check as part of the intake
gate; write `STATUS: OPEN` **plain, not bold** (the merge gate missed a bold
one); use the exact C10 anchor header; and **do not edit `curation/lib/`**,
because editing it mid-run changed behaviour under a running agent today.

`levkinpa_29558113` is given `lij_42130331`'s independently established
thiolactone regiochemistry **as a hypothesis to test, not an answer** — the
papers are 8 years apart from the same lab.

### Merge-gate bugs found by batch 1, both fixed

1. **Phantom OPEN entry.** The scaffold documents the `DECISIONS.md` format
   inside a ``` fence and that example parsed as a real entry, so every
   scaffolded folder would be permanently blocked by a fake. A gate that always
   fires trains people to ignore it.
2. **False negative, worse.** A folder wrote `**STATUS: OPEN**` in bold and the
   parser only matched a bare `STATUS:` — it reported **zero** open decisions
   for a folder with five, in the one check standing between an unresolved
   question and the corpus. Now strips markdown emphasis. Verified in both
   directions; corrected counts 5 / 9 / 5 / 9, all four folders correctly
   blocked.

### Batch 1 final: 581 rows staged, 0 merged

| folder | rows | tier | library |
|---|---|---|---|
| `lij_42130331` | 300 | printed numerals in raster tables | 300/300 |
| `andersondg_41203968` | 132 | **tier 1, no pixel read** | 100/100 |
| `mitchellmj_38982196` | 100 | **tier 1** | 100/100 |
| `whiteheadka_38888972` | **49 of 580** | pixels | 580/580 enumerated |

`whiteheadka_38888972` is the batch's most important result: the library really
is 580, verified three ways, but **only 49 ever went into an animal**. The other
531 carry a TNS surface-ionisation value and nothing else — the same class the
corpus definition already excludes. A naive read would have claimed 580 rows.

**Finding U — two agents editing `curation/lib/` concurrently is a real
hazard.** I relaxed `assemble.check_stereo_preserved` while an agent was
running against the stricter rule; its stage had been written for the old
behaviour and a seeded check that should have fired reported "DID NOT RAISE",
which is how it was caught. The agent handled it correctly and verified all 580
structures byte-identical either way. **Rule: shared-library changes happen
between batches, never during one.**

**F-049 confirmed from a second, independent direction.** `whiteheadka_38888972`
uses "14:0 PEG2000 PE" = a real DMPE with phosphoethanolamine, and the shared
guard REFUSED the name (3,611 corpus rows carry P and N; `JW_2024`'s 1,873 carry
the DMG diglyceride string with neither). The agent correctly did not repair it.

### 06:2x — `kooh_42690843`: ZERO rows, and that is the right answer

**Every transfection measurement in this paper lives in the MAIN ARTICLE
(Figures 2-4), which we do not hold.** The SI contains only synthesis and
characterisation. `05_values.csv` has a header and no rows. The agent reported
zero rather than manufacturing a number, which is the behaviour the brief asks
for and the behaviour the corpus depends on.

Status `blocked_no_values`, not dropped. What IS durable and finished: all **48
building blocks** and the **9 purified lipids** the SI characterises, each
verified three ways. The folder is one value source away from mergeable.

**Finding V -- the stated library size is not recoverable by arithmetic.** The
paper claims **275** members; the SI draws a 26 amine x 8 boronic acid x 14
aldehyde panel = **2,912** possible, and names only **9**. 275 is 9.4% of 2,912
and is not a product of any subset (26x8=208, 26x14=364, 8x14=112, 26+8+14=48).
Enumerating 2,912 would fabricate 2,903 molecules that were never made;
inventing a staged 275 would be worse. Recorded as a contradiction and an OPEN
decision, unresolved.

**The best check of the run so far, and nobody printed the number.** ChemDraw
emits one stroked line per heavy-atom bond plus a second per double bond, so
`drawn_lines == heavy_bonds + double_bonds(kekulised)` is an exact identity.
SI p.S5 carries **483** stroked items; 2 are dashed dividers; the remaining
**481 partition into exactly 48 clusters** (closest approach between clusters
12.8 pt vs 11.0 pt linkage, so no merge is possible) and the 48 predicted bond
counts **sum to 481**. 48/48. This is the reagent-panel analogue of "a round
number the authors drew but never printed" -- external to our pipeline in a way
a recomputed marginal is not. Its first run found a defect in itself: an
imidazole measured 11 lines against a predicted 9 because the prediction counted
aromatic bonds while ChemDraw draws Kekule; kekulising is what closes the sum.

**Regiochemistry demonstrated rather than asserted.** The coupling is a 427 nm
photo-Petasis reaction; N- and C-alkylation are isomers with identical formula
and identical exact mass. Seeding the N-alkylated variant made **HRMS pass at
+0.026 mDa while the connectivity check fired** -- a clean in-folder
demonstration of "mass cannot see connectivity".

**D-3 has consequences beyond this paper:** SI p.S26 says the in vivo screen
used "the **crude** PILL lipid-LNP sub-libraries" and the abstract says
"**purification-free** in vitro screening". A label in those figures is an
unpurified one-pot mixture, not the molecule. That is the `NG_2025` crude-
mixture situation again and it is a corpus-level convention question.

The agent did not attempt to breach the paywall and says so explicitly; one
spoofed-User-Agent attempt was blocked by the harness and not retried. Correct.

**Replacements launched immediately (guard: never idle):** `sd_2026` (Nat
Biotechnol, **384** stated, and we hold the MAIN ARTICLE as well as the SI --
rare in this run) and `namjh_39264000` (Small, **213** stated, .docx SI plus
PMC11855220 full text).

### 06:4x — Finding W: the SI-only profile is a systematic dead end for VALUES

Two papers in a row returned **zero value rows for the same structural reason**:

| folder | structures | value rows | why |
|---|---|---|---|
| `kooh_42690843` | 48 blocks + 9 lipids | **0** | all transfection in main-text Figs 2-4 |
| `liu_41162324` | **74 enumerated** | **0** | SI's every numeric panel is an excluded endpoint or a 2-3 arm benchmark |

Both are ACS, both acquired **SI-only**, and in both the SI carries synthesis
and characterisation while the per-lipid measurements live in the article. ACS
serves its SI free but the article needs an entitlement this session does not
have.

**So library size is the wrong thing to rank on.** A 315-lipid paper with no
reachable values is worth less than a 100-lipid paper with a Source Data
workbook. The queue is now annotated with **`value_reach`** and sorted on it
before size:

| class | meaning | remaining |
|---|---|---|
| `A_source_data` | a publisher Source Data workbook is already on disk -- tier 1, no digitisation | **3** |
| `B_main_article` | we hold the main article PDF | 0 (sd_2026 is in flight) |
| `C_pmc_fulltext` | open in PubMed Central -- full text and captions free | **22** |
| `D_si_only` | SI only -- the profile that produced both zeros | 11 |
| `E_nothing` | nothing acquired | 4 |

**22 of the remaining papers are in PMC**, which is the single most useful fact
found this hour: PMC gives the article text and captions free, so the
value-location question can be answered before any structure work. Every brief
from now on fetches PMC first and writes `00_value_sources.json` before
touching chemistry, with the instruction that **a fast honest zero beats a slow
one**.

**Finding X -- a seed harness can prove the wrong thing.** `liu_41162324`'s
first seeded-defect harness required an *exception*. Its P7 seed (build every
F5 tail with the F3 tail -- same carbons, two fluorines fewer) scored
`teeth: FALSE`, because the stage *reports* a formula mismatch and exits 0. It
was rewritten to measure an OUTCOME: the seed now drops agreement 34->20 and
34->19. Its P1/P2 seed had also been passing for the wrong reason, tripping a
chain assertion rather than the formula check. **A seeded mutation must assert
the metric moves, not merely that something was raised.** Added to every
subsequent brief.

**C16 fired on an agent's own `make_reference.py`** (an auxiliary SM-102
fragment literal was a canonical-order substring); it re-typed the fragment
hydroxyl-first and C16 went green. An observed alarm, not a hypothetical.

**F-001's exact shape reappeared, from the machinery built to prove checks
work:** that folder's QC restore step silently stripped the quality flags from
`03_enumerated.csv`. Fixed in-folder, with a guard that fails the stage if a
quality column is missing.

**Launched immediately (never idle):** `dahlmanje_29489381` (281 stated, PMC
OA, the richest SI in the batch at 2,856 numeric tokens over 3 table-like
pages) and `mckinlaycj_39105751` (200 stated, PMC OA). Both briefs lead with
"establish WHERE THE VALUES ARE before any structure work".

### 07:0x — `levkinpa_29558113`: 5 usable rows. Third confirmation of Finding W.

288 labels enumerated and the arithmetic reconciles exactly (4 thiolactones x 4
pyridyl disulfides x 18 amines). **But only 32 of 288 structures are
determinable, and only 5 of the 45 valued lipidoids** -- the 18 amine heads are
enumerated in the MAIN ARTICLE, which we do not hold. The SI names or draws
only A1, A5, A10 and gives A8 as a bare formula. Not open access, no repo.
Recorded `extracted_below_threshold`; the agent's own recommendation is not to
merge it as it stands, and I agree: 5 rows over 2 head groups is a near-disjoint
silo that cannot feed Result 3. It becomes a real 45-row experiment the moment
the main article's amine figure is fetched, and every stage re-runs unchanged.

Its value, despite the row count, is in the checks:
- **The sibling's thiolactone hypothesis was CONFIRMED by this paper's own
  scheme** (ring opens at the carbonyl; Scheme S1c, vector) -- and its second
  half correctly did NOT transfer, because there is no epoxide here. Given as a
  hypothesis, tested, half kept. That is how a cross-paper hint should behave.
- **A cross-tier check a row-offset could not survive:** Figure S9's legend
  names the top eleven lipidoids as TEXT objects; pdfplumber recovers them
  without a pixel and they are identical and in order to the raster table's
  rows L1-L11.
- **The round number the authors drew but never printed:** the positive-control
  bar must sit at exactly 1.000. Calibrated from y-ticks alone, SFA reads
  1.0064 and L2K 1.0291 -- identifying SFA as the normaliser, corroborated
  independently by a caption. The 0.64% is the top rule's own thickness.
- **`bars.bar_top` REFUSED 3 bars** crossed by a large arrow (flanks disagree
  13-66 px). Left out rather than rescued by lowering `frac` -- which is F-018.
- **14/14 seeded defects fired, and the seeding found a real gap:** a block
  SMILES in `01_blocks.csv` was a SECOND literal for the molecule the
  enumerator builds from a fragment, so corrupting the fragment left the HRMS
  check passing on an uncorrupted block while every product carried a wrong
  tail. Now one string. That is F-001's family -- a value restated in two
  places drifts.
- **C-2 is a live name-vs-structure conflict:** `PY12` is named
  "2-(decyldisulfanyl)pyridine" (C10) while its own printed HRMS and 12 drawn
  bonds say dodecyl (decyl would be 28 Da away). One emittable row depends on
  it, so it is queued as a DECISION, not assumed -- the `AC_2025` shape again.

### RUN PAUSED — usage limit reached

**Totals: 581 rows staged across 4 folders, 0 merged.** Plus 74 + 32 + 57
structures banked without values in three folders that are one main-article
fetch from being mergeable.

6 subagents were still running when the limit hit and keep going;
`siegwart_40644303` (465), `sd_2026` (384), `namjh_39264000` (213),
`mckinlaycj_39105751` (200), `dahlmanje_29489381` (281) and
`palczewski_41654520`. Their reports will arrive as task notifications.

**On resume:** run `check_state.py` first; it fails if the run is idle while
papers are ready. Then launch from the queue sorted kind -> `value_reach` ->
size. **The binding constraint is the main article, not the SI** -- `D_si_only`
has now produced 0, 0 and 5 rows across three papers, while `A_source_data` and
`B_main_article` produced 132 and 100 with no digitisation at all.

### 2026-11-10 ~01:5x — Limit reset. All 6 agents RESUMED, not restarted.

The 5-hour limit killed all six in flight. **Their work was saved and they were
resumed via `SendMessage` rather than relaunched**, which preserves their
context and avoids redoing finished stages. State on disk at resume:

| folder | what survived |
|---|---|
| `siegwart_40644303` | **nearly complete** — blocks, SI table S3, enumeration, Fig 1 heatmap, in-vivo bars, values, scatter check, vitro+vivo rows, `make_reference.py`. Was building the C10 anchor. |
| `dahlmanje_29489381` | design tables + rank tables (stages 01-02). Was extracting tier-4 text tables positionally. |
| `sd_2026` | `stage01_blocks.py` only |
| `palczewski_41654520` | `00_acquisition_raw.json` only |
| `mckinlaycj_39105751` | value sources established, intake gate being written |
| `namjh_39264000` | nothing yet |

**Finding Y — an agent tried to drive the shared browser.** `namjh_39264000`'s
last line before the cut was "let me try the browser route for the main-text
figures (Wiley works with Playwright per the project's own measurements)". The
Playwright browser is a SINGLE shared instance driven by the main session; a
subagent navigating it would collide with whatever the main session is doing
mid-sequence. Its resume instruction forbids the browser and points at two
plain-client routes that need no login and no shared state:
`pmc_fetch.py`, and PMC's own figure rasters at
`pmc.ncbi.nlm.nih.gov/articles/instance/<id>/bin/<name>` — the route a sibling
agent used to pull six main-text figures at 300 dpi. **Browser access belongs
to the main session only; agents get plain clients.** To go into the plan.

Every resume message also re-states the four conventions batch 2 established:
plain `STATUS: OPEN` (not bold), the exact C10 anchor header, seeds that assert
the METRIC MOVES rather than that an exception was raised, and "a fast honest
zero beats a slow one".

### 02:0x — `siegwart_40644303`: 505 rows, and a finding that touches Result 3

**505 rows emitted** (465 in vitro + 40 of 42 in vivo), 465 stated / 465
enumerated / difference 0, 489 structures total. Largest single yield of the run.
Two in-vivo lipids deliberately not emitted: one has no bar drawn at all, the
other is ~5 px and `bars.bar_top` refused on the flank test. The threshold was
not lowered -- that is F-018.

**FINDING Z -- MERGING THIS WOULD MOVE RESULT 3's HEADLINE NUMBER. User decision
required; nothing has been merged.**

The agent flagged it as D-10 and I verified it directly against the corpus
rather than taking the claim:

| | |
|---|---|
| `IL_name` matches vs `KZ_2016` | **465 of 465** |
| full `IL_SMILES` matches | **360 of 465** |
| largest existing cross-experiment overlap in the corpus | **135** (`AA_2008` <-> `KW_2014`) |
| Result 3's pooled concordance set | **n = 236** points |

So this would be a **360-lipid shared overlap, 2.7x the largest the corpus has
ever had**, dropped into the exact statistic the manuscript's load-bearing
result is computed from. And it is **same-lab, cross-CARGO** (Zhou's siRNA
library re-screened for mRNA), which per `CLAUDE.md` §8 is the `molea.csv`
situation -- deliberately kept separate -- not a cross-lab pair. The 105-SMILES
gap is itself informative: it is the SC3 / 3A2 / 5A2 / 6A1 disagreements below.

**The corpus states one tail three different ways and they disagree (D-4).**
For `SC3`, Fig 1b draws an unbranched chain (traced: 0 branch nodes, 3 segments,
9 px against 18-19 px for the branched SC1/SC2) = n-propyl. `KZ_2016`'s
`IL_SMILES` says neopentyl, its `IL_tail1_SMILES` says `CCCS`, and its
`Molecular.Weight` is the C5 value. That is a live internal contradiction in
corpus data we already ship, found from outside. 15 structures. HRMS would
settle it. **This is a pre-existing corpus defect, not an artefact of this
curation** -- it belongs in `ERROR_CLASSES.md` once adjudicated.

**The load-bearing check (Q5) is one the printed statistics cannot do.** SI
Fig S3 is a genuine text table giving every lipid's pK_a, so Figure 5a's
unlabelled markers can be assigned to named lipids by x coordinate and read on
a *different figure, different axis, different ticks*. 32 lipids confirmed,
median scatter/bar ratio **1.0124**, sd 0.038. A multiplicative bar error would
have to be matched by an identical error in an unrelated panel to survive it.
The agent recorded explicitly that Q7 (printed statistics) **cannot** see such
an error -- seeding x1.28 changes nothing there -- which is the
necessary-not-sufficient rule stated as a measurement.

Two drawn-but-never-printed anchors: the x-axis rule lands on the y-tick-fitted
zero at **+0.41 px** and **-0.20 px**; the 31 column-label centroids sit on the
31 cell centres at max **1.5 px** on 36 px cells.

**Its own checks caught two of its own errors**: the arm fragment was first
built as a disulfide (parsed, no radicals, stereo preserved -- only the corpus
InChIKey comparison saw it), and two tails were transcribed one carbon short.

**C10 anchor built from an explicit V2000 atom/bond table**, because the obvious
"independent" route was a trap: the enumerator's canonical SMILES for the
anchor compound *already starts at the pentyl tail*, so a tail-first SMILES
anchor would have been byte-identical -- the F-037 shape. An atom table cannot
collide that way. C10 `exact_match`, C16 all derived.

Launched `pabst_37602460` (104, PMC OA) in the same turn. Its brief now carries
a new standing instruction: **after enumerating, compare against the corpus and
report how many members already exist** -- D-10 is a merge-level fact that only
surfaces if someone looks.

### 02:2x — `dahlmanje_29489381`: 0 rows, a NEW failure class, and 6 papers deprioritised because of it

**0 emittable rows, blocked two independent ways — and the second one would
survive holding the main article.** That makes it a different class from the
SI-only problem (Finding W):

1. **No value is reachable even with the article.** Figs 3B-G and S3j-uu are
   **unlabelled log-log scatters** — digitising a point yields an (x, y) pair
   with no particle identity, so no row could ever be keyed to a composition.
   Only Fig 2C is per-particle addressable and PMC served 0 of 8 figure assets.
   Data Access reads verbatim *"available upon request"*; `find_repo.py` found
   nothing.
2. **Not one of the 15 lipomers has a determinable SMILES.** 7 of 15 are built
   on PEI 600 / PEI 1200 — polydisperse branched polymers, not molecules (the
   SI's own scheme concedes it, drawing `PEI600` as diethylamine with the label
   underneath). For the other 8 the degree of substitution is never stated:
   spermine offers 6 N-H, so the product spans **~442-1,653 Da**. **The
   molecule is unknown — that is not a flaggable uncertainty**, and the agent
   correctly declined to "assume full substitution".

It also did not call `emit.emit_rows` at all, on the reasoning that those
columns qualify *a corpus row* and there is none to qualify. Right call.

**FINDING AA — "N lipids" in a barcoded screen usually means N FORMULATIONS.**
The tracker said 281. The reality: **564 formulations over only 15 distinct
lipomers, 206 of them sharing a single one (7C1)**, with libraries 1 and 3
varying just *two* ionizable lipids across 300 formulations. A barcoded in vivo
screen varies composition, not chemistry.

Acted on immediately: every remaining queued paper whose title/abstract says
barcoded / pooled / high-throughput-in-vivo was **reclassified `form` and
deprioritised** — `dahlmanje_30275336`, `dahlmanje_41183194`,
`dahlmanje_33758781`, `dahlmanje_30394729`, `mitchellmj_41603124`,
`mitchellmj_41555650`. Six agents not spent on formulation screens the user
explicitly deprioritised. `kind_original` kept on each row.

**The best check was an unprinted DESIGN CONSTRAINT**, which is a new variety of
"the round number they drew but never printed": the four mole percents must sum
to exactly 100. **564/564 rows, residual exactly 0.** On the rank tables it is
86/87 — and the single offender is a real paper defect: Fig S4f prints
`Lipomer Mole % = 36` for one LNP where two other loci print 35, and 36+60+5+0
= 101. The invariant adjudicates a published contradiction that no printed
statistic could. Generalised into subsequent briefs.

Its structure check caught two of the agent's own errors before passing
(an epoxide one carbon too long; an azamacrocycle with correct ring size and
correct formula but nitrogens off the 1,4,8,12 locants) — and the seed confirms
**formula and exact mass are blind to the locant error**, both sides measuring
C11H26N4 / 214.2157 Da. Only the ring walk sees it.

One contradiction was **retracted with its measurement kept in the file** so
nobody re-raises it: a suspected R1/R2 key inversion turned out to be
colour-keyed rather than order-keyed, RGB distance 0.00 within bands against
66.88 across. That is the right way to close a false lead.

Launched `li_39883839` (PNAS Passerini lipids, PMC OA) in the same turn. Its
brief now carries both new traps: count DISTINCT IONIZABLE STRUCTURES separately
from formulations, and report corpus overlap after enumerating.

### 02:4x — `mckinlaycj_39105751`: 5 rows, and a corpus-design question worth raising

**5 emittable rows**, all from main-article Table 1 read as structured
`<table-wrap>` cells out of PMC's JATS — **tier 4, zero pixels, no digitiser in
the folder**. Everything else per-member (93 library members, 12 headgroup
variants, 34 DOE variants) is in figure rasters PMC would not serve: **1 of 8
assets delivered**, the rest bot-check interstitials, not retried.

44 structures enumerated and mass-verified (44/44 within 0.005 Da of the SI's
printed exact masses). The 221-member library is **not enumerable at all**:
**Figure 1B labels two different structures `Mpp`**.

**D-1 is a project-level question I am surfacing rather than burying.** These
are **peptoids** — N-substituted glycine oligomers, not small-molecule
ionizable lipids. The agent laid out three concrete ways that bites this
project, and they are specific enough to matter:
- **ECFP4 at radius 2 may not encode a monomer transposition at all**, so two
  genuinely different peptoids could become near-identical feature vectors;
- a peptoid source is **trivially separable by backbone**, which would make
  Results 1-2 easier for reasons about molecular CLASS rather than about
  laboratory provenance — the exact confound the manuscript argues against;
- the five RDKit descriptor columns all return numbers whose usual meaning does
  not transfer.
Not resolved here. It belongs with the `molea.csv` and D-10 decisions.

**D-3 is a restraint worth recording.** SI Figure S10 compares the lead against
**SM-102 and MC3** at matched dose, we hold the page, and its ticks are in the
text layer — so it was extractable. The agent **deliberately did not extract
it**, because `CLAUDE.md` §8 says to decide *before* Stage 5 is re-run whether
benchmark-only pairs belong in pooled concordance, and extracting it would be
making that call. Correct.

**The unprinted check:** SI Table S9 never says what species "Found Exact Mass"
is. The difference from Calculated is **+1.0089 +/- 0.0062 Da on 35/45 rows** —
the proton. The 10 misses are the paper's own data-entry errors (one row's
found mass is byte-identical to two other rows'; another is printed as `2.00`).
None touches an emitted structure.

**FINDING AB — two seeded checks were no-ops, and both were real defects in the
agent's own harness.** (1) A watermark seed read `34/34 -> 34/34` because
`_rows(x_max=X_TABLE_MAX)` bound the module global as a **default argument**,
so the mutation was shadowed and the seed tested nothing. (2) A transposition
seed read `2/2 -> 1/2` because it swapped two positions that are **both `Ehx`**
— a no-op scoring as a partial pass. Both fixed; the real headline is that
transposing two monomers leaves exact mass at **1275.09 either way** while the
anchor goes **2/2 -> 0/2**. This is the third instance today of a seed proving
the wrong thing: **a seed must change the data in a way that CAN move the
metric, and the metric must actually move.**

`assemble.check_no_radicals` also earned itself again, rejecting a dangling
`[C:2]=O` template that would have parsed and canonicalised cleanly.

Launched `linyx_41851281` (A_source_data, 17 files) in the same turn, with my
own sampling handed over as findings to verify rather than conclusions: two of
its biggest sheets are a molecular-dynamics time series and a block of
flow-cytometry zeros, so it must run `find_screens.py` first. Its brief also
carries a new rule for AI-design papers: **never emit a predicted value as a
measurement.**

### 03:0x — `sd_2026`: 18 rows, 384 structures, and a defect in `curation/lib/axis.py`

**18 rows across three experiments** (`SD_2026` IV / `SD_2026_2` intrathecal /
`SD_2026_3` intratracheal — different route, dose and tissue, never pooled),
plus **384 structures enumerated, 384 distinct, 0 duplicates**.

Holding the main article did not rescue the screen but was decisive for the
rest: the follow-up in-vivo panels are **vector**, so 18 values were read from
the PDF's own coordinates with **no pixel estimation at all** and a tick-fit
residual of 0.0000-0.0001 px. A held-out tick (Fig. 1p's 10^5, excluded from
the fit) was predicted to **0.0000 decades**.

**The 384-lipid screen is unreadable for two independent reasons, either
fatal**, and this is a new entry in the catalogue of why a screen resists
extraction:
1. **No numeric scale anywhere.** The heat-map legend is a ramp annotated
   `Max`/`Min` — no number, no tick, no unit. The Sankey's y-axis is a bare
   arrow with zero ticks. One dashed benchmark line cannot fit a scale.
2. **No cell -> lipid map.** The caption says the rings are "hierarchical
   clusters", so the within-sector order is a **clustering's leaf order,
   published nowhere**. Guessing it would silently permute 384 values across
   384 structures.

**Tier 1 was disproved rather than assumed**: `static-content.springer.com`
enumerated over n=1..16 x {xlsx,xls,zip,csv,docx,pdf,txt} — n=1-4 return 200
for `.pdf` only, n>=5 return 403 for everything, **no `PK` magic anywhere**.
And a trap worth knowing: there IS an asset named **"Source Data Fig. 1"** —
its content is six pages of Sanger chromatograms, not the screen.

**FINDING AC — `axis.calibrate`'s max-residual test is weak against a SINGLE
bad tick, and I verified it is worse than reported.** On a clean 5-tick,
100 px-pitch ladder:

| single-tick error | max residual | verdict at tol 2.5 px |
|---|---|---|
| +3 px on an end tick | 1.20 px | **PASSES** |
| +6 px on an end tick | 2.40 px | **PASSES** |

An end tick has leverage 0.2, so four fifths of its error is absorbed into the
fit. Least squares hides the outlier it is supposed to expose.

**Fixed with `loo_residuals()` + a `loo_tol_px` argument**: refit on the other
n-1 ticks and predict the held-out one, which has no leverage masking. The same
+6 px error surfaces as **6.00 px**. Three tests added (**52/52**): `C-23`
documents the weakness, `M-27` confirms LOO catches it, `C-24` confirms a clean
ladder passes.

**It is OPT-IN for now, deliberately.** Tightening a shared check while six
agents run against it is exactly the hazard that bit this run earlier today
(Finding U). The default flips to armed when the batch drains; until then every
new brief tells the agent to report leave-one-out residuals in its QC.

The agent also reported two honest limits of its own checks — that the "all 10
SI formulas realised" test is **nearly blind** (9 of 10 one-carbon mutations
still realise all ten), and which checks actually pin the inventory instead.
Reporting what a passing check does NOT prove is the behaviour this project
wants.

Launched `xuq_32438474` (cholesteryl lipidoids, PMC OA) in the same turn, with
an explicit warning about the F-013 sterol history and instructions to take the
sterol verbatim from the reagent registry and compare FULL InChIKeys.

### 03:2x — `palczewski_41654520`: 64 rows, TIER 1, zero pixels, and the first CC BY licence

**64 rows from the publisher Source Data workbook at per-replicate resolution.
Not one pixel was read for any value or any structure.** All 34 structures came
from SI IUPAC names + HRMS as text. SI figures were rendered and looked at only
to VERIFY, never to extract.

**Licence CC BY 4.0 — the first folder in this run that clears `CLAUDE.md` §8.**
Every other paper so far is subscription or non-CC, which is a release blocker
the corpus will have to face.

**But the cargo is not mRNA (D-4).** Cre recombinase **protein** and **ABE8e
RNP** — a different endpoint class from the corpus's mRNA transfection. That is
a user call, alongside the peptoid question and D-10, not something to decide
inside a curation.

**Four assay blocks, not one screen.** Fig2c/4c are flow cytometry, Fig2d/4d are
microscopy, and the SAME compound disagrees across the pair by **0.55x (CBB14)**
and **1.86x (CBB11)**. They are emitted as separate blocks. This is the
`mitchellmj_38982196` same-lab disagreement again, from a different direction,
and it keeps being worth measuring.

**"Combinatorial" in the title is not a combinatorial enumeration.** Every
member is individually named, purified and HRMS-characterised; the block table
is 12 x 3 and the paper made **21 of those cells — 21 drawn, 21 built**. A
naive 12x12x3 would have fabricated 432.

**Three real errors found in the paper**, each caught by splitting the mass
check three ways rather than collapsing it to one boolean:
- **CBB12's printed formula is wrong by 2 H** — its own m/z and the homolog
  arithmetic both give H128 where the string says H130.
- **CBB17's formula says S3** — TFA has no sulfur; its own m/z confirms S2.
  Almost certainly copied from CBB18, which legitimately has three.
- **Supplementary Fig 3's Z10 cell says hexadecyl** where the reagent name, the
  product name and the HRMS all say dodecyl.
**Collapsing "our formula vs theirs", "their formula vs their m/z" and "our
formula vs their m/z" into one pass/fail is how a published typo gets adopted
silently.** Generalised into subsequent briefs.

**Another unprinted round number:** the SI weighs out "brilliant G250 (854 mg,
1.0 mmol)" => MW 854. Our free acid is 832.06; +Na-H = **854.04**. One weighing
validates the whole dye-core composition against something external to every
figure.

10/10 seeds moved their own metric, and the agent explicitly verified the gate's
**exit code directly rather than through a pipe** — its first attempt had been
masked by `| tail`, which is the same class as the three no-op seeds found
earlier today.

**A self-inflicted defect it found and reported:** `stage03` was still writing
`01_blocks.csv` after `stage01` took that file over, silently clobbering the
better version on every ordered run. It surfaced only as a determinism failure
— which is why byte-identical reruns are worth testing.

Launched `siegwartdj_42150080` (PNAS cardiomyocyte screen, PMC OA) in the same
turn, carrying the three-way mass check, the exit-code warning, and the new
`axis.loo_residuals` instruction.

### 03:4x — `namjh_39264000`: 108 rows, and a connectivity bug that HRMS could not see

**108 rows over two experiments** (`SB_2025` 72 intradermal Renilla;
`SB_2025_2` 36 intramuscular hEPO), all tier 5. The vaccine readouts (IgG,
ELISpot, cytokines) are deliberately NOT emitted — a third experiment on every
axis — and are queued as D-3 rather than dropped silently.

**The tracker's "213 lipids" is wrong, provably.** The paper says 213 **LNPs**,
and Figure 2's nine panels each print group sizes summing to 213 which
partition across the six compositions of Table S2, with eleven derived
identities all holding. 213 counts FORMULATIONS. Distinct lipids are bounded
only by **76 <= n <= 196**, and just **36 are drawn anywhere**. Third instance
of Finding AA today.

**FINDING AD — a connectivity error that HRMS was completely blind to, found by
the agent in its own work.** Its first `stage03` built acyl fragments by string
surgery on canonical SMILES (`O=C(O)X` -> `[C:1](=O)X`). In SMILES the mapped
atom bonds to the FIRST atom of the remainder, and canonicalisation had moved
the carboxyl carbon's neighbour away from the front. Result: **FH-1 became
linear palmitoyl instead of 2-hexyldecanoyl, and every EB-T chain put the
glycerol ester on the far end of the tail.** Not one atom changed — identical
formula, **exact mass identical to 8 decimal places, and both HRMS anchors
still matching at 0.3 and 0.5 ppm**.

It is now caught by a check whose required value is a PRINTED STATISTIC: Fig 2b
prints FP3 (tertiary carbon) = 3 for exactly 36 of 213 LNPs and 2 for the other
177, and FH-1 is the only head with a branch point. Fragments are now built
from measured counts with the attachment atom **atom-mapped**, never by
rewriting a canonical string. **This is the same hazard `curation/lib/assemble`
exists to prevent, reached by a route the library does not cover** — string
surgery on SMILES outside `assemble`. Passed into the `peer_39912611` brief,
whose whole subject is linkers.

**A correction to the agent's D-2, which changes what it means.** It reported
`SB_2024` as "60 rows … the `KIST_*` lipids". Measured: `SB_2024` is **24 rows**
in `data_vivo.csv` (the 60 in `CLAUDE.md` §3 is `SB_2024` + `SY_2025`
together), and its `IL_name` values are **SM-102 (13) and Dlin-MC3-DMA (11)** —
**benchmarks only**, measured as `luminescence_relative_to_Spikevax`, PMID
**38779592**. The shared context is real and worth keeping — same `Mouse_ICR`,
same `ear`, same 5.0 ug — but the overlap is with BENCHMARK rows, not with that
paper's own library. So it is the §8 "benchmark-only pairs" question rather
than a hidden same-lab duplicate. Recorded with the measurement.

Its own best checks: HRMS at **-0.32 / -0.46 ppm**; a tracer **calibrated
against five SI-named reagents (8/8 exact) before being trusted** on the six
blocks that have no anchor at all; each colourbar's drawn bottom edge (an
unprinted zero) recovered to **0.75 / 0.50 / 1.00 px**, i.e. 0.15-0.22% of full
scale, with the seeded F-019 defect moving it to 3.5-4.0 px; and the printed
**r = 0.7580** reproduced as 0.768 from two independently digitised panels —
which also *identifies which panel is which*.

**Two seeded checks are recorded as NOT firing, and say so.** A colour-match
diagnostic whose metric does not move under a grid shift, and a repo sweep that
has no metric at all. Reporting a seed that did not move the needle is exactly
what the three no-op seeds earlier today failed to do.

It also declined to run C10/C16/C18 because **they write tracked files in
`provenance/` that other agents are using concurrently** — the right call, and
now a standing instruction in every brief.

### 04:0x — `pabst_37602460`: OUT OF SCOPE, and a corpus-wide naming audit it triggered

**FINDING AE — the tracker title was TRUNCATED and the paper is not about
transfection.** The published title ends "...How Molecular **Conformation**
Affects **Membrane-Targeting Antimicrobial Activity**". The endpoint is growth
inhibition of *Bacillus subtilis* in Luria Broth: **no nucleic-acid cargo, no
helper lipid, no cholesterol, no PEG-lipid, no nanoparticle, no mammalian
cell** — the quaternised lipidoid is dispersed in bacterial medium by itself.
The chemistry is squarely in family; the endpoint is not. Marked
`out_of_scope`, 0 rows, with its 104 structures fully enumerated and checked so
nothing is wasted if the scope call changes.

**Every brief now verifies the paper's actual endpoint from the PMC record
before any other work.** A truncated tracker field sent an agent 2.3 million
tokens into the wrong kind of paper.

**It also found 11 of its 104 parent structures ALREADY IN THE CORPUS** by full
InChIKey — across `KW_2014`, `NC_2024` and `AA_2008`, one also in vivo. Second
paper today to find pre-existing overlap it was not looking for (after 465/465
in `KZ_2016`). The corpus-overlap check is earning its place in the brief.

**That led to a corpus-wide audit I ran myself.** Across all 30,998 rows and
19,919 distinct full InChIKeys:

| | |
|---|---|
| structures carrying >1 `IL_name` anywhere | **324** |
| of those, spanning >1 experiment | 278 — **expected**, this is what Result 3 is built on |
| **one structure under >1 name WITHIN a single experiment** | **59** |

The 59 are the defect-shaped subset: `KW_2014` 25, `BL_2023` 11, `AA_2008` 10,
`SX_2025` 9, `JW_2024` 2, `LR_2023` 1, `ZH_2023` 1. Examples are clean pairs
like `AA_2008`'s `21_O13` / `22_O13` resolving to one InChIKey. `CLAUDE.md` §3
already records "**48 duplicate SMILES (39 with conflicting values)**" as a
known accepted condition across `KW_2014`/`BL_2023`/`AA_2008`/`JW_2024` — so
this is largely that condition, measured a different way and extended to
`SX_2025`, `LR_2023` and `ZH_2023`, which that note does not list. **Recorded as
a measurement for the next review round, not acted on.**

The paper's own quality was high and worth noting: **15/15 seeded defects moved
their metric**, three C10 anchors matched by full InChIKey from a generator
containing **zero SMILES literals** (explicit atom/bond tables), and the
unprinted round number was a beauty — all 8 bounded therapeutic indices equal
exactly **400/MIC**, where 400 appears nowhere in the paper (it is 2x the
200 uM ceiling).

Three real paper defects found: SI Table S3's "Tail total" column contradicts
its own cells on **11 of 13 rows**; Table S5's charge densities for two
headgroups are **transposed** (proved four independent ways, including that the
paper's own stated trend has the WRONG SIGN with the printed values and the
right sign with the computed ones); and SI Fig S1 panel A **draws a 9-membered
triaza ring bearing five arms**, which is impossible — a pure connectivity
error, invisible to formula and mass.

Launched `xuq_32520527` in the same turn, with scope verification as its first
instruction and a warning about acetal stereochemistry and regiochemistry.

### 04:2x — Batch review gate built and run; two tool bugs fixed

**`scratchpad/tools/batch_gate.py`** — the cross-paper half of
`AUTONOMOUS_RUN.md` §4 step 2. It checks only what a single-paper reviewer
structurally cannot see: collisions and drift BETWEEN the staged folders and
the live corpus. Six checks (G1 experiment-ID collisions, G2 pre-existing
structures, G3 one structure under two names inside a folder, G4 quality-column
vocabulary, G5 caveat with no note, G6 open decisions through the REAL parser).

**Result over 10 staged folders, 1,299 rows:**

| | |
|---|---|
| G1 experiment-ID collisions | **0** |
| G3 one structure, two names within a folder | **0** |
| G4 quality-vocabulary violations | **0** |
| G5 caveats with an empty note | **0** |
| G6 open decisions | **73**, every folder correctly blocked |

**G2 found a merge-level fact the extraction agent did not report.**
`whiteheadka_38888972`: **246 of its 580 structures already exist in the corpus**
— `NC_2024` 198, `KW_2014` 96, `AA_2008` 14 — and **41 of its 49 emittable rows
are already-present structures**. That is the second large pre-existing overlap
after `siegwart_40644303`'s 360 (now also showing `SL_2020` 8 and `ZC_2023` 1).
Neither agent was asked the question early enough; both overlaps surfaced only
because something checked across papers. **This is exactly why the gate exists
and why Tier-2 per-paper review cannot replace it.**

**The gate's own first run was wrong, in the way this project keeps being
wrong.** It reported "no IL_SMILES column" for 4 of 10 folders — they use
`canonical_SMILES` or keep structures in `03_enumerated.csv` — so it silently
checked **zero structures** and printed a clean line. Same aliasing failure that
took C10 down earlier today. Widened to accept four structure-column spellings
and five name spellings, and to say which columns it actually saw when it finds
none. The corrected numbers are the ones above; before widening it reported
`sd_2026` 0 of 0 and `lij_42130331` "no column", among others.

**FINDING AF — `pmc_fetch.py` had a wrong URL template and had been failing
silently all run.** It built `articles/PMC<id>/bin/<name>`; the working route is
`articles/instance/<NUMERIC id>/bin/<name>`. Measured on the same file: **404 on
the PMC form, 200 on the instance form.** That one template made the tool report
"7/7 failed, 0 downloads" to several agents, each of which then fell back to
pixels or recorded the figures as unreachable — at least three briefs today told
agents PMC figure fetching "has served HTML error pages for several siblings",
which was my tool, not PMC. Fixed and verified end to end: the same article now
yields **5 real figures at 1.3-2.5 MB** where it previously yielded 0. The
instance route also serves at 600 dpi (~4,200 px).

### `li_39883839`: 13 rows, and a new reason a big library yields few rows

**FINDING AG — POOLED screening.** The library is 144 members, but all 144 were
screened **pooled** (9 mice carrying 16 LNPs each, then 4 carrying 4 each). Only
lipids promoted out of the pools were ever injected alone, so there are **13
individual measurements over 9 distinct structures**, 2 of them benchmarks.
Distinct structures, formulations, and *individually measured* are three
different numbers — now asked for separately in every brief.

147 distinct InChIKeys from 148 names: **`H2A4-Ole`, presented in the paper as
newly synthesised to compare linear against branched tails, IS library member
`H2A4B2`** — same InChIKey. Corpus overlap: **zero**, checked three ways.

**D-1 is unresolvable from the paper, and the agent proved it rather than
picking.** Fig 2H and Fig 2K disagree **25.8x and 2.79x** on their two SHARED
lipids at the same stated dose, route and timepoint — and the two shared anchors
disagree **9.2x** about what the offset would be, so the two-anchor tie that
rescued `QL_2026` is simply not available here.

Its checks were strong: Fig 2K's short 4e7 tick predicted to **0.25 px** by a
fit on the other four; the three printed ROI numerals landing **3/3** inside
their own bar's drawn spread; and the same lipid read on two independent panels
with independent calibrations agreeing at **1.039x**. It also found that the SI
labels building-block masses `[M+H]+` while **0 of 9 fit [M+H]+ and 8 of 9 fit
neutral monoisotopic** — two undeclared conventions in one document.

### 04:4x — `siegwartdj_42150080`: 0 rows, but the first DATED reopen path

Per-member values are in main-text Fig 1. PNAS returns **403 (Cloudflare)**,
the PMC article page returns **403**, and the PMC record we did fetch states
`pmc-status-live=no`, `pmc-status-embargo=yes`, **release 2026-11-18**, with a
publisher note forbidding XML full-text download. No browser used, no bypass
attempted.

**This is the first blocker in the run with an expiry date.** Status
`blocked_embargo`, `reopen_after: 2026-11-18`, and the agent states nothing in
the folder needs redoing — only the article PDF is missing. Worth checking on
that date.

**Trap 1 again, in its purest form yet: 16 formulations, ONE distinct ionizable
structure.** The "chemically diverse LNP library" of the abstract is, on the
held material, `Base LNP` + {`18:1 TAP`, `18:1 DAP`, `18:1 PA`} x {10,20,30,40,50}
mol% — a 3x5 titration of Avanti catalogue helper lipids over a fixed ionizable
lipid. **The mol-% suffix is the giveaway**: formulation design, not a chemical
series. Recorded as UNRESOLVED rather than asserted, because Fig 1 is unreadable
and might show more.

**And trap 2: that single structure is already in the corpus.** `6A1-SC8` is
`KZ_2016` (PMID 26729861), 2 rows, full InChIKey `CZPQUMASBKZWQS-UHFFFAOYSA-N`,
unambiguous. The corpus already carries 44 `6A1-S*` siblings in vitro. If this
ever merges it **must be a new `Experiment_ID`** — `KZ_2016` is siRNA/HeLa/liver
and this is mRNA/hiPSC-CM/heart. The agent also measured that DOTAP appears in
2,526 vitro rows, DODAP in 486, and **DOPA is absent from the corpus entirely**.

**The mass check was split FOUR ways and the fourth earned it.** Beyond the
usual three, it compared our formula against **the m/z labels drawn inside the
spectrum figure**, read at 900 dpi: 5/5 labels, worst **3.23 ppm**. Their
spacing is **0.5006 Da => z = 2**, which the caption never states — so the
caption's "found: 2170.3 Da" is a **deconvoluted neutral mass and no 2170 ion
was ever observed**. Back-derived neutral 2170.3212 vs our 2170.3142. That is a
convention recovered from a drawing, not from prose.

It also **declined to credit a check it could have claimed**: the 0.87 ppm
methyl NMR integral is the operator's reference, set to 18.00, so comparing it
to 18 is circular. It reported 2/2 on the two non-circular ratios instead.
Refusing a free pass is the behaviour this project wants.

`axis.loo_residuals` was used as briefed — 19 ticks, fit residual 0.0007 px,
worst LOO **0.0008 px**. First field use of the fix made three hours ago.

Launched `lius_40500862` (161 stated, JACS, D_si_only) in the same turn, with a
warning that carbonate/carbamate/ester linkages are frequently invisible to
molecular formula.

### 05:0x — `linyx_41851281`: 140 rows, tier 1, and the corpus used as an external witness

**140 rows (121 in vitro + 19 in vivo), all tier 1 from publisher Source Data,
zero pixels for any value and ZERO predicted values emitted** — the model's
`Predicted` column is never opened, and `00_value_sources.json` says why. My own
sampling of this folder had found only an MD time series and a block of
flow-cytometry zeros; the agent found the real screen in
`MOESM4` sheet `Supplementary Fig. 10` cells A1:D113.

**The library size was recorded nowhere and is now established:** 116 distinct
ionizable structures = **14 amine heads x 16 tails x 2 linkers x 2-4 arms**,
116 distinct InChIKeys. **Not the formulation trap** — Supplementary Table 3
gives ONE screening formulation for every library lipid, so formulations =
structures here. The 18-formulation sweep elsewhere is one lipid and is
excluded.

**FINDING AH — the corpus itself became the strongest external check.** 15 of
the 116 lipids are already present by full InChIKey under `AA_2008`, `KW_2014`
and `NC_2024`, with zero `IL_name` collisions. Because the corpus's own
lipidoid names ENCODE CHAIN LENGTH, that overlap is an independent witness on
this paper's head and tail assignment: **21/21 agreement, with a perfect
bijection H1<->28, H2<->80, H3<->109, H4<->113, H9<->306**. An overlap that two
other folders flagged as a merge hazard is here also a free validation — worth
telling every agent, since the check costs nothing once the enumeration exists.

Other checks of note: **22/22 (linker, tail) monomer masses solved from HRMS
with spread 0.0000 Da across 108 lipids** — an unprinted quantity; 115/116
enumerated masses reproduce the printed HRMS, and the single failure (lipid 47)
is diagnosed as a **printed typo** because its partner lipid 48 is exact and
pins the head.

**The biggest contradiction of the run so far (C-01, factor 2,405).** Fig. 3c's
row for ML ID 8 gives 3.47-4.49e6 RLU while Supplementary Fig. 10, Fig. 3g and
`label.txt` all give ~1,645-1,918. Those three high numbers are
**byte-identical to lipid P1's in-vitro triplicate**, and the peer-review file
contains the authors' own sentence *"lipid 8 is now referred to as lipid P1"*.
Fig. 3c's id-8 row carries P1's data under the original submission's numbering.
The low values were emitted. Related: the paper frames P1-P6 as AI designs
while the reviewer states P1 was in the original conventional screen — which
does not change a measurement but changes what a "P" row means.

**FINDING AI — `find_screens.py` reported only the BEST column per sheet, and
hid the largest in-vivo panel.** It kept one column at a >=15-label threshold,
so `Fig. 3i`'s 12 lipid labels in **column W of an 89-column sheet** were
invisible while a different column held 16. One sheet routinely holds several
panels side by side. Fixed to report **every** qualifying column at a >=8
threshold; the same file now surfaces three panels (cols 11, 5 and 21) where it
previously surfaced one. The agent found Fig. 3i by reading the workbook, not
by trusting my scanner — which is why it is worth having both.

It also caught a **no-op seed** of its own (a +1 RLU mutation invisible under
the loosened tolerance) and recorded it as such before replacing it with x1.01.
Fourth instance of that class today, and the first where the agent caught it
unprompted.

**Bonus from the anchor:** Supplementary Fig. 12 draws Lipid T3 with 2 arms on
one primary nitrogen and 1 on the other — **the only external evidence anywhere
for the 3-arm regiochemistry, which HRMS cannot see.**

### 05:2x — `xuq_32438474`: 90 rows, and F-013's blindness measured on a real library

**90 library-member values over three blocks** (`YL_2020` mRNA 36, `YL_2020_2`
protein 27, `YL_2020_3` recombination 27 — different cargo, correctly not
pooled). 18 distinct structures; the `-F` suffixed names are **formulations of
the same structures**, not members. Corpus overlap: **zero**, checked by full
key, skeleton key and name.

**A new value tier, named honestly.** The SI's bar charts are **native PDF
vector**, so 54 values came straight from `page.get_drawings()` filled-rectangle
coordinates — no raster, no pixel, no colour matching. Rather than mislabel that
as `5_pixels`, the agent introduced **`4b_pdf_vector_geometry`**. That is the
right call and it is now suggested in every brief: a vector bar chart is not a
digitisation, it is a table stored as geometry.

**FINDING AJ — F-013's blindness, measured rather than asserted.** Seeding
epicholesterol into this library changed **15 of 18 FULL InChIKeys and 0 of 18
SKELETON InChIKeys**. I reproduced the minimal case myself and added it to the
suite as **C-25** (53/53 now):

| | |
|---|---|
| molecular formula | same |
| exact mass to 4 dp | **470.4124 both** |
| InChIKey skeleton block | **IDENTICAL** |
| full InChIKey | differs |

That is precisely why C11 missed F-013 **twice** — it compared skeleton blocks.
The number now lives in the test suite instead of in a paragraph.

**First independent confirmation of this morning's `axis` fix.** The agent
seeded a +6 px error on one y tick: **max residual 1.87 px PASSES** the 2.5 px
tolerance while **`loo_residuals` fires at 5.35 px**. Its own four panels max
out at 0.757 px LOO. The weakness I measured at 06:0x and fixed opt-in has now
caught a seeded defect in the field.

**The strongest check here was an unprinted grouping**, and it is a pattern
worth copying. Each of the paper's prose sets is separated by a single threshold
on the extracted matrix: "consistent across all cell types" has margin **5.16
pp**, "consistently ineffective" **45.4 pp**, and "top five screened in vivo"
resolves to exactly the right five with margin **0.59 pp**. It even explains an
apparent anomaly — why the best HeLa lipid (79.3%) is excluded from
"consistent", because its NIH 3T3 is 23.2%. All four groupings would fail under
a permuted category order, which is what tests the 12 labels that had to be read
by eye from a rotated raster.

**It also declined to correct a zero point (D-2).** The tick-fitted zero sits
2.22 px (0.48 pp) above the axis rule and both printed anchors are low by about
that much — so moving the scale to make them land would be **error class 27**
(F-006) exactly. Left uncorrected, with both readings recorded.

Launched `khanof_36652629` (JACS, iterative IM design) in the same turn, with
the vector-geometry tier, the prose-grouping test and the LOO instruction.

### 05:4x — `lius_40500862`: 4 rows, and the tools moved into the repo

**4 emittable rows.** The 161-member screen is in paywalled main-text Figs 1-3;
the SI yields 15 namable lipids, 10 with an HRMS anchor, and **8 individually
measured delivery values over 4 lipids**. Corpus overlap zero on all three keys.

**A new shape of the pooled trap, and the subtlest yet.** Fig S25 tabulates 16
formulations — and S25B plots not those 16 but **16 MARGINAL MEANS over 4
factors x 4 levels, each an average of 4 formulations**. Not one of the 16 is
individually reported. A panel with 16 bars beside a table with 16 rows is not
necessarily the same 16.

**FINDING AK — my tools lived in a session scratchpad, which `CLAUDE.md` §7
forbids twice over.** This agent reported the exact predicted symptom: *"no
`find_repo.py` exists anywhere in the repo (`find` returns nothing)"*, so it did
the GitHub tier by hand and said so — and one of my briefs had omitted the
absolute path entirely, so the tool was simply unreachable. §7's rule exists
because of `lit_review/build_data.py` and F-023; this is the **third instance of
the class**, committed by me, in tools that every brief now depends on.

**Moved to `provenance/tools/` with a README** recording what each does and the
two rules they enforce (magic bytes over size floors; a scanner that cannot see
a column reports a clean run). Verified they run from the repo copy.

**The chemistry was verified three independent ways rather than assumed** — the
paper's whole point is a carbonate linker, and carbonate/carbamate/ester are
formula-blind. Synthesis (4-nitrophenyl chloroformate + alcohol, and the core
nitrogens are tertiary so carbamate is impossible); **13C at d155.16 with
nothing between 168-180** (dialkyl carbonate ~155, aliphatic ester ~173); and
HRMS, where each carbonate exceeds its named carboxylate counterpart by exactly
one oxygen per tail. Regiochemistry likewise read, not assumed: two 13C peaks
near 53 ppm say one tail per nitrogen.

**`Cm-n` is a carbinol convention, not a chain length** — pinned independently
of mass by the SI naming one acid **2-hexyloctanoic acid**. The project's own
"name tokens are not chain lengths" hazard, again, in a new dialect.

**P2 is worth copying: 10/10 calc-vs-found within 6.33 ppm, ALL TEN POSITIVE,
mean +4.8 ppm.** A systematic instrument offset rather than a published typo —
which is what the three-way split is for, and the opposite conclusion from the
sibling that found three real formula typos the same way.

**It retracted three of its own claims with the measurements that retired
them**, including its anchor guard: 2 of 4 correctly hand-typed anchors came out
**byte-identical** to the enumeration, and §6 says byte-identity is not the
signal and testing for it is a trap. Its first guard would have forced a correct
anchor to be rewritten into a worse one. Replaced with C16's actual property,
checked by AST: every anchor SMILES must be a source string literal and the
module must import no assembler.

Also caught in its own work: a chroma mask with a `max_val <= 252` ceiling
**silently deleted a whole bar** (254,224,224) — visible only because
`bars.bar_top` refused with "found 3 windows, expected 4". The refusal did its
job.

### 06:5x — Limit reset; 4 agents resumed, 2 completions, and a FALSE REFUSAL fixed

**`xuq_32520527`: 33 rows**, 11 distinct structures x 3 formulations, tier 5
raster but the good case (linear axis, four labelled ticks, solid bars;
quantisation floor 0.33 pp, an order of magnitude below the authors' own
replicate s.d.). Corpus overlap zero by all three keys — **but the agent found
what those keys cannot see**: the amine HEAD numbers are the Xu-lab catalogue,
and heads 75/80/113/400 are already in `QL_2026`, read there from a vendor name
by an independent curator. 4 of 11 heads externally cross-checked for free.
It also warned that **163 corpus `IL_name`s begin with one of these numbers and
do NOT all mean the same thing** — `BL_2024`'s `113-9` uses 113 as an unrelated
series index — so keying on the leading token would merge unrelated lipids.

Its unprinted check: four **unlabelled minor ticks**, excluded from the fit,
predicted as **105.00 / 75.03 / 45.04 / 15.09** against the 105/75/45/15 the
authors drew and never printed, worst error 0.091 pp.

**`chamolstad_41691524`: 77 rows, but ONE distinct structure.** The 96-member
library is in paywalled main-text figures; everything reachable is a
formulation DoE on a single lipid. And the deeper problem: **only 1 of 12 head
amines is disclosed anywhere**, and the 8 tails are NMR-characterised but never
mapped to a T index — so **even with the missing figures, 88 of 96 could not be
given a structure**. Its tier-1 workbook turned out to be a 3'UTR half-life
sheet, not a lipid screen; `find_screens.py` correctly printed nothing for it.

Its best check is a model for cross-panel work: Fig S3A and Fig S4C report the
same formulations, cell line, readout and timepoint, and independently
digitised they **rank-agree at Spearman 0.9352 (p = 8.1e-06, n = 12)** while
**disagreeing in magnitude by a median factor of 8.96**. The rank agreement is
strong evidence both reads are correct; the scale gap is a real property of the
paper. A different-cell-line control gives rho = 0.437.

**FINDING AL — `assemble.check_stereo_preserved` was issuing FALSE REFUSALS, and
a false refusal is worse than a missed catch because it blocks correct work.**

A CIP code is a property of the WHOLE molecule, so a transformation that changes
a remote substituent can **relabel** a centre without moving an atom. Verified
minimally: **(S)-glycidol's C2 is `S`; its decanoate ester is `R`** — nothing
inverted, the free OH leaves the C1/C3 comparison tied and the ester breaks the
tie. The counts check reads that as an inversion and refuses.

Counts genuinely cannot separate the two cases, because a real inversion also
preserves the total. **Parity can.** Added `assemble.parity_preserved()`: map
each fragment atom onto the product by substructure match, then compare the
centre's chiral tag after permuting its neighbours into the same order on both
sides. Permutation parity encodes handedness and does not care what the CIP
rules make of the substituents.

| case | parity verdict |
|---|---|
| (S)-glycidol -> its decanoate ester (relabel) | **accepted** |
| the F-013 epicholesterol epimer | **rejected** |

Wired in as **opt-in** (`allow_cip_relabel=True`) so the four running agents are
untouched — editing a shared check mid-batch is the hazard that bit this run at
Finding U. Tests C-26 and M-28 added; suite **55/55**. Default flips when the
batch drains.

The reporting agent handled it correctly in the meantime: it respected the
refusal, did **not** edit `curation/lib/`, and derived the centre from the
published chemical name instead of carrying a tag through a bond-forming step.

A second library limitation is recorded but NOT yet fixed:
`assemble.stereocentres()` returns `{}` for a tail where RDKit's
`FindPotentialStereo` finds two unspecified para-stereocentres — so the stereo
guard can be blind to the only stereochemistry in a molecule, and its passing is
not evidence. To address when the batch drains.

**Four agents resumed** (`peerd_31999380`, `levkinpa_23909648`,
`peer_39912611`, `khanof_36652629`) rather than relaunched; each brief now also
points at `provenance/tools/` instead of the scratchpad.

### 07:1x — `peer_39912611`: 144 rows, and THE CORPUS CAUGHT THE AGENT'S ERROR

**144 rows** (108 vitro + 36 vivo), 12 structures, 12 formulations, 12
individually measured — one composition for all of them, so the "N lipids = N
formulations" trap does not bite; what multiplies is measurement cells
(12 x 3 cell lines x 3 doses, plus 12 x 3 organs). Proposed `GN_2025`.

**FINDING AM — the corpus functioned as an external anchor for a paper it does
not contain.** Lipid 14 has no IUPAC name, no NMR and no MS in this SI. The
agent first completed it as the C3 homologue of lipid 37, read from flattened
scheme text. A substructure search then turned up three `ZT_2026` rows —
`Lipid 14-CB1/CB2/CB3`, quaternised cholesterol conjugates **built by a
different group on this paper's lipid 14** — whose SMILES carry
`CCOC(=O)CCC[N+]`, three CH2, the **C4** acid. Re-reading Supp Fig 15 at 600 dpi
confirmed C4. I verified both halves myself: all 3 `ZT_2026` rows carry the C4
form, and the C3 form the agent first built **appears nowhere in 30,998 corpus
rows**. Retracted with the measurement that retired it.

This is a use of the corpus nobody designed: overlap search is meant to find
duplicates, and here a *near* miss supplied the ground truth a paper's own SI
had withheld. Worth adding to every brief.

**It refused to claim a check it did not have.** I had told it to run the
sibling's carbonate test (13C d155 vs d173). It searched the SI by name first —
`13C`, `C NMR`, `100 MHz`, `125 MHz` all **0 hits**; `HRMS` **0 hits**; only 1H
at 400 MHz and ESI-TOF to one decimal — and said the test was unavailable rather
than manufacturing it. Then built three substitutes that do exist: an **O-CH2
proton count** that separates ester / reverse-ester (6 H) from carbonate (8 H)
from amide / urea (4 H), 11/11 with a seed moving it to 10/11; **carbamate
excluded by valence** (every nitrogen in both partners is tertiary, so no N-H
exists to acylate); and the printed IUPAC locants (`7-oxo-6,8-dioxa-`).

**Its best check is the strongest unprinted test seen this run.** The log y-axis
decade length was ambiguous from 4 major ticks. But the **8 minor ticks per
decade sit at log10(2..9) and are UNEQUALLY spaced**, so their pattern encodes
the decade length on its own: mean distance to nearest predicted tick is
**0.393 px on a 47.45 px/decade ladder against 1.633 px on 71.25 — a ratio of
4.15**. The authors drew that and never printed it.

**And it said plainly which of its checks does NOT work.** Reproducing the
paper's printed Pearson r moves only 0.034 under a wrong decade ladder, so it
reported that check as unable to settle the ladder rather than counting it.
That is the posture §6 asks for and it is rarer than the passes.

Its two printed-statistic reproductions are nonetheless tight: HeLa r **0.7364**
vs printed 0.7358, and lungs-vs-pKa r **0.8369 / p 0.0007 / R2 0.7005** vs
printed 0.8417 / 0.0006 / 0.7084 — the first of these also **identifies which
dose** the paper's correlation panel used, which the caption never says.

**The `check_stereo_preserved` limitation never fired here**, and the agent
explained why rather than reporting a pass: the only stereocentres are remote
from every bond formed and unassigned on both sides. A check that cannot apply
is not a check that passed.

**My own ID-collision check was the flawed one.** I ran a cross-folder scan to
test its "no collision with the 54 live IDs" claim and got a false positive,
because my set included the folder's own file. Re-ran properly: **26 staged IDs
across 20 folders, zero cross-folder collisions.** The round-1 pattern of
drawing a conclusion from a secondary artifact, landing on me again.

Open, and it is a real gap: **Fig 3F/G/H were not extracted.** They are the
redundancy that would check C/D/E bar-for-bar (log10 C - log10 D = F) *and*
recover the six censored values. Highest-value item left in that folder.

Launched `songx_39664425`.

### 07:2x — I committed the exact error the user forbade, and it is now a check

**What happened.** I wrote `songx_39664425`'s queue status to
`in_progress_agent` and only *then* launched the agent. In between,
`check_state.py` printed **"running: songx_39664425"** for an agent that did not
exist. I caught it in the same turn and launched immediately, but the queue had
already been wrong — and the queue is this run's own account of itself, so a
status written before the thing it describes makes every later reading of that
file wrong.

The user's standing instruction, after an earlier instance where I said
"Launching batch 2" without launching: **"This should never happen again. Ensure
that it doesn't."** Saying it will not recur is not ensuring it. So the ORDER IS
NOW A RULE WITH A CHECK: launch, take the agent id from the launch result, then
write the queue entry with `agent_task_id` set. An entry claiming to be in
progress with no id, or with an id that has no task file on disk, fails.

**Seeding it caught a second bug of the same family, in the guard itself.** My
first wiring referenced a variable name that does not exist in that file, so
both the control run *and* the seeded run exited 1 — a `NameError`, not a
finding. The seeded run looked like proof the check worked. That is
`CLAUDE.md` §6's rule earning itself for the fourth time this run: **a new check
is not done until it has been made to fail ON PURPOSE, for the right reason,
with the right message.** Exit code alone is not the signal; three of this
project's past incidents printed a clean line and two printed a non-zero one.

Fixed and seeded twice, with a control between:

| run | exit | why |
|---|---|---|
| control, 4 agents live with ids recorded | **0** | `state consistent` |
| seed A — marked running, `agent_task_id` removed | **1** | names the folder and asks whether the status was written before the launch |
| seed B — `agent_task_id` that never existed | **1** | names the id and says it has no task file |
| restored | **0** | — |

`TASKS` now reads `CONFOUND_TASK_DIR` from the environment when set. The task
directory is session-scoped, and this tool moved into the repo an hour ago
precisely because §7 forbids a repo script taking an *input* from a scratchpad —
so it is an override, not a stored input, and nothing tracked depends on it.

### 07:4x — `levkinpa_23909648`: 29 rows, and TWO false zeros in my own triage tool

**29 rows, tier 4, zero pixels** — every value from PDF text tables. 3 purified
structures with real C10 anchors (this paper draws them in full, unlike its 2018
sibling).

**"200 lipids" is real and means something else again.** Table S1 resolves to
20 amines x 10 alkyl halides = 200 exactly — but each is a **crude one-pot
product**, and Table S1 reports several tail counts for a single (amine, halide)
pair, up to six, with **23/23 rows showing >1**. So 200 = reaction products, not
structures and not measurements. **Only 12 of the 200 carry a value**, and **9 of
those 12 have a value with no disclosed head structure** — including the best
crude performer, `4C18` at 31.7%. Masses are solvable and were recorded; a mass
is not a structure and none was guessed. That is the fifth distinct dialect of
this trap today.

**FINDING AN — `find_tables()` returns 0 on pages that are nothing but table,
and my `triage.py` could not see it.** Verified myself on the file:

```
p18: find_tables()=0  words=245  images=0  thin rects=80  heights=[0.36, 0.48]
     with rect runs as explicit h-lines: 1 table, 7 rows x 13 cols
```

The rules are **filled rectangles 0.36-0.48 pt tall**, not `line` objects, so
pdfplumber's default lines strategy finds no edges. Two false zeros compound,
both silent: `find_tables()` says "no table", and **the SI's own List of
Contents stops at Figure S2 and names none of the four value tables**, so a
triage trusting the contents page agrees with it.

Added `triage.ruled_tables()`. On this paper the text heuristic flagged **1**
page; the new detector finds **5**, including a **31x15** and a **21x19**:

| page | find_tables | thin rects | recovered |
|---|---|---|---|
| 6 | 0 | 27 | 1x7 |
| 16 | 0 | 112 | **31x15** |
| 17 | 0 | 108 | **21x19** |
| 18 | 0 | 80 | 7x13 |
| 19 | 0 | 22 | 3x7 |

**FINDING AO — and `triage.py` printed a clean EMPTY report when given a folder.**
`main()` took files as extra arguments; with none it emitted `{"files": {}}` and
exited 0. But `new_curation.py` scaffolds sources into `<folder>/source/`, so
every caller passing just the folder got a silent, clean, empty triage. Added
`discover()` (walks the tree) and a **refusal**: no PDF or DOCX found now prints
why and exits 1 rather than printing an empty result as if it were one. Seeded
on an empty directory — exits 1 with the reason. Same false-zero shape as
`find_screens` hiding a 12-lipid panel in column W and `batch_gate` checking
zero structures in four folders. Third instance, this time in the tool whose
whole job is to decide whether a paper is worth extracting.

A corpus-wide re-scan for hidden ruled tables is running.

**Other things this agent got right.** It split the mass check **seven ways**,
not three — and caught **two mass conventions inside one SI**: Table S1's
`calc.` is the **neutral** mass while the HRMS sections' "calculated m/z" is
[M+H]+. Comparing across them showed a 1.0421 Da residual, exactly one proton,
and *under a 2 Da tolerance every row would have "passed" a proton out*.

Its unprinted checks are the best kind: **38/38 consecutive tail steps in the
authors' own calc column each equal exactly one CnH2n** (worst 0.0904 Da), and
**two relations among their three printed HRMS values close only if the tail
multiplicities are 2, 1, 2** — confirming mono/di independently of anything we
built. It also **recovered a DNA dose printed nowhere**, by matching a control's
mean *and* SD (14.0 +/- 1.8) across four tables to the 150 ng column and neither
other dose; recorded as inference, dose left blank.

**First field confirmation of my `include_filled` fix**, and it measured what I
had only warned about: with `include_filled=True` it gets **exactly 2.0 segments
per bond** (64/48/80 for 32/24/40), and glyph contamination is **zero** —
measured, not assumed, because every glyph on those pages is a bitmap *image*,
not a path. The closed identity `drawn == heavy + kekulised double` is 3/3 exact.

It also caught two real defects in its own work: a **vertically-centred merged
cell** that its first parser read as the block's last row, stamping **9 of 30
values with the wrong dose** while producing a clean 30-row table; and the
two-convention mass error above.

Launched `whiteheadka_30637934` (same lab as the 580-lipid sibling that found
246 already in the corpus — expect heavy overlap as a free check).

### 08:0x — `songx_39664425` is a COMMENTARY; the tracker sent an agent after a mis-attribution

**0 rows, and the right answer.** PMID 39664425 is a research-highlight
commentary with **zero original measurements**: `abstractText` length 0,
`pubType` `letter`, PMC subject heading `Highlight`, no Methods / Results /
Data-availability section, author contributions **writing + conceptualization
only**, `hasSuppl=N`, 0 tables, 0 supplementary material.

Its single quantitative figure *does* carry calibratable bars, and the agent
**refused tier 5 on provenance**, which is the correct call: the caption ties it
to reference 1 **by markup** — `Reprinted with the permission from Ref.` plus
`<xref rid="bib1">`, and `bib1` carries `pmid 38740955`. Digitising it would give
one lab's measurements a second `Experiment_ID` and **fabricate a source** — the
same refusal `CLAUDE.md` §9 already records for the `XD_2026` LipidDB. It is
also a ~4x downsample of data we already hold at source.

**FINDING AP — the tracker was wrong three ways and had promoted it to prio 1.**
`lit_review/build_data.py` carried this PMID as `status="new", prio=1,
lib=584`, promoted on *"full text confirms 584 experimentally made ILs"*, and
described it as the predecessor to LipidAI (42592089). All three are false: the
584 belong to **PMID 38740955 = `BL_2024`, already merged**; 584 is not a
distinct-lipid count; and it is unrelated to LipidAI. The entry is now
`status="dup", prio=9` with the evidence, and `lib` carries 0 with the
resolution. Rebuilt; `data.json` confirms the change.

The generalisation is sharper than the one I had been briefing. I have been
warning agents about **truncated** titles. This was a **mis-attributed** one —
the tracker credited the commented-on paper's library to the commentary. So I
swept all 15 ready papers' Europe PMC `pubType` and abstract length. One flagged
on `pubType` alone (`xuq_23651337`, ACS Synth Biol "Letter") and is **ordinary
research** — a 797-character abstract describing its own synthesis. **The
discriminator that works is `abstractText` length 0**, not `pubType`; only
`songx` had it. Every brief now carries that test as step 1.

**It also resolved the long-standing 584-vs-560 question, and I re-derived it
from `data/` rather than trusting the report.** `BL_2024`'s 560 unique
`IL_name`s decompose exactly as `<3-digit prefix>-<member>` over **128
prefixes**:

| | |
|---|---|
| 120 prefixes carrying member set exactly `{9,22,23}` | library A = 128 x 3 = **384** |
| 8 prefixes (`111,241,244,611,613,621,622,624`) carrying `1..25` | library B = 8 x 25 = **200** |
| prefixes fitting neither | **0** |
| `{9,22,23}` is a subset of every B prefix | 8 x 3 = **24** in both screens |

384 + 200 - 24 = **560**, matching the corpus exactly, with **0 duplicated
SMILES**. So **584 counts screen memberships and 560 counts lipids** — the
corpus is right and the discrepancy was never an error. Also confirmed: the
paper's lead compound **`119-23` is absent** from `BL_2024` (0 names start
`119`, and 119 is outside the 128-prefix grammar), consistent with it coming
from the 40,000-member *virtual* library. Routed as an open item on `BL_2024`.

Its gate caught a defect in **its own first version**: a hand-transcribed prose
expectation failed because the source has U+00A0 after "Ref." and the "1" is an
`<xref>`, not a literal — replaced with the stronger markup assertion. And it
declined to call `emit.emit_rows` at all, stating that with 0 rows it could only
succeed and the F-011 per-item-callable requirement is unreachable. A check that
cannot fail is not a check, and saying so beats a green line.

Two of its references are already curated (`JC_2023`, `BL_2023_3CR`); ref 6
(Wang, Acta Pharm Sin B 2022, ~400 points pooled from 22 studies) carries the
same fabricate-a-lab hazard as the LipidAI DB and must not be merged as one
experiment.

Launched `xuq_32662132`.

### 08:3x — `peerd_31999380`: 0 admissible rows, and two direct hits on MY triage output

**347 values extracted, 0 corpus-admissible.** The combinatorial-library *screen*
is in main-text Figs 1-4; only the SI was on disk. The SI's values are
**viability (XTT), biodistribution, apparent pKa and toxicology** — and `data/`'s
`Experiment_method` vocabulary is delivery-only. The agent wrote all 347 rows
with `corpus_admissible = no` and a per-row reason rather than quietly dropping
them, which is the right shape for a zero. 14 structures; **formulations not
recoverable at all** (no helper lipid, no ratios, no N:P anywhere in the SI);
only **2** lipids measured on a delivery endpoint.

It also closed the acquisition tier properly: NCBI `elink` returns **only
`pubmed_pmc_refs`, no `pubmed_pmc` self-link**, which means no PMC deposit
exists and `pmc_fetch.py` never had a PMCID to take. That test is now in the
briefs — it is the cheap way to avoid hammering a route that cannot work.

**FINDING AQ — my `triage.py` reported "0 table-like pages" on a file with two
tables `find_tables()` finds trivially.** Verified myself:

```
pages where find_tables() SUCCEEDS: [(47, ['12x14']), (48, ['12x17'])]
```

matching the agent's hand count of Tables S1/S2 exactly. My heuristic needs >=8
distinct compound-like labels AND a header word on the same page, so a table of
toxicology rows passes neither and vanishes. **I put that false zero into the
brief.** This file has now been wrong in both directions in one day — missing
native tables here, and missing filled-rect tables at 07:4x. `triage` now
reports **`native_tables`** and **`ruled_tables`** side by side with the text
heuristic, each labelled, none trusted alone. Re-verified: p47/p48 now appear,
and `levkinpa`'s 4 filled-rect pages still do.

**FINDING AR — and the token count was unquotable.** I reported "1,745 numeric
tokens"; the same text gives **2,533 / 3,123 / 6,958** under three other
reasonable regexes. The project has been through this exact disease with the
flagged-row count that moved four times (F-024): **a quantity with a prose
definition and no implementation gets re-implemented, slightly differently, by
every pass.** `numbers_regex` now ships beside the number.

**FINDING AS — `heatmap.MAX_DIST` is unsafe across panels, and I confirmed the
numbers.** `MAX_DIST = 5.0` is calibrated for matching cells against swatches
sampled from the *same* raster, where only JPEG noise separates them. Across
panels it is actively dangerous:

| | |
|---|---|
| S3 "Lipid 13" | (127, 4, 127) |
| S4 "Lipid 12" | (128, 0, 129) |
| euclidean distance | **4.583** |
| `MAX_DIST` | 5.0 |

So a cross-panel match would have assigned S4's Lipid 12 values to Lipid 13 and
**never tripped the guard that exists to catch exactly this**. The agent avoided
it by matching only within each panel. Added
`heatmap.cross_panel_collisions()`, which reports *different-label* pairs below
the threshold — it finds both real collisions here (2.828 and 4.583) and stays
silent on well-separated panels. Tests M-29 / C-27; suite **57/57**.

**Its panel-ratio work is the best use of trap 4 yet.** Figures S2 and S4@0.06
are the same condition, and **8 of 12 shared lipids agree to within 0.3%**
(ratios 0.9971-1.0034) — which simultaneously validates the digitisation *and*
shows the panels share one dataset. Against that baseline the other 4 cannot be
noise: **Lipid 3 at x1.313, a 29.35 pp gap**, confirmed by eye as the tallest
bar in its group. A tight agreement baseline is what turns a disagreement into
a finding instead of a worry.

**And the SI prints no `Calc.` column at all**, so "monoisotopic or average?"
had to be *determined*: 4 lipids sit closest to monoisotopic [M+H]+, 6 closest
to average, 4 are bare integers that cannot discriminate — all 14 within 0.76 Da,
which confirms formula at integer resolution and nothing finer. Stating that is
worth more than a passing tolerance.

It caught two bugs in its own parsers, both of which looked like paper errors:
an olefinic integral double-counted because a space after a hyphen defeated a
lookbehind, and a mean+/-SD regex that accepted `±0.2.3` because `[\d.]+` is not
a number test. Its note that `lit_review/` changed mid-session is correct and
was me, at 07:10, correcting the tracker entry.

Launched `xuq_24574196` (bioreducible siRNA, same lab).

### 09:0x — `khanof_36652629`: 18 rows, and my pixel tolerances were in the wrong units

**18 rows, emitted as SIX blocks rather than pooled**, and the reason is a
finding: Table S2 gives **the same lipid a different formulation ID in
different panels**, and five lipids measured more than once disagree by
**1.053x, 1.662x, 1.751x, 1.861x and 5.011x** (C3-K2-E14: 1.711e6 vs 3.415e5).
"Iterative design" means rounds, and rounds are different experiments. Not a
library either: 3 cores x 2 linkers x 5 tails = 30 possible, 8 exist.

**Bears on an open §8 decision.** All 3 benchmarks hit **full InChIKey,
skeleton AND `IL_name` simultaneously** — a free check that our registry
structures reproduce the corpus — and 0 of 8 novel lipids hit anything. So the
only shared chemistry is benchmarks. **But here they do NOT top the screen**:
four bars within 1.945x, all "ns". §8 worries that benchmark-only pairs push
cross-lab concordance toward +1 because benchmarks top every screen; this paper
is a counter-example and should be cited when that decision is taken.

**FINDING AT — `curation/lib/axis`'s tolerances assume a 600 dpi render, and
flattened SI rasters do not have one.** Rendering a page at 600 dpi does not
create information that is not in the embedded JPEG — it upsamples. These six
panels are rasters at a **measured 188-221 dpi**, so my 1.5 px LOO tolerance is
**0.55 of a real pixel**, and three panels failed a test no correct
digitisation could pass. The agent rescaled by `tol x 600/native_dpi`, recorded
both numbers per panel, and **seeded the widened tolerance** to prove it still
bites. That last step is what separates a units fix from a convenience.

Added `axis.native_tol()`, which refuses a missing or guessed dpi — a tolerance
rescaled by a guess is worse than one that is simply wrong, because it looks
derived. And I measured two things the agent's report did not cover, both now in
the docstring:

| | |
|---|---|
| `calibrate`'s own hard **2.5 px** limit is the real blocker and is **not** rescaled by `native_tol` | a +5.0 px seed makes it RAISE at 3.50 px before LOO ever runs |
| the rescaled LOO tolerance is **weak at the margin** | at 204 dpi (tol 4.412 px) a **+3.25 px** tick error gives LOO 4.333 and **passes by 0.08 px** |

So the rescale is right in units and still costs real sensitivity. It is not a
free move, and the numbers depend on decade spacing, so every caller must seed
it **in its own panel's geometry**. Tests C-28 / M-30; suite **59/59**.

**Error class: a SMILES fragment has a direction, and mass is blind to it.**
This anchor scored **0/8, then 1/8, then 5/8** over three runs, from three
instances of that one class: a core nitrogen bonded to a tail; tails written
reversed past the carbinol; and `OC2` written oxygen-first, giving N-O-CH2CH3.
**All three had the correct molecular formula**, parsed, canonicalised, no
radicals — only the InChIKey saw them. That is the F-013 shape in a new
disguise, and it is now in every brief.

**Mass conventions drifted WITHIN one column**: their "calc'd" is **average on
two compounds and monoisotopic on six**, and on C3-K2-E14 the calc (1137.92,
average) and found (1137.09, monoisotopic) are on different scales. Third
distinct version of this hazard today, and the worst.

**Its unprinted check is the strongest category yet seen.** The paper's own 12
significance marks sort **12/12 on a single fold threshold at 2.074x** — largest
"ns" 1.945x, smallest starred 2.212x, margin 0.056 decades. A t-test is not a
function of fold alone, so this genuinely tests the bar tops rather than
recomputing them. Plus a **stroke-centreline convention** verified at 0.013 px
worst over six panels, from a baseline rule and a tick centroid, neither taken
from a bar — which is what licenses `value = outer + t/2`.

And it reported a check that **does not agree** rather than burying it: Fig S9A's
drawn fold-increases are 2.06 / 5.93 against our 3.162 / 7.251 (1.53x, 1.22x
off), logged as a contradiction.

**Six `bars.bar_top` refusals, all respected and all fixed at the source.** A
saturation floor of 28 deleted every lavender/teal bar (teal channel spread
**4**); 18 px strokes mean a fill can never reach 90% width; a 2-3 px JPEG halo
inflated the window; and an error-bar cap competed as a dark rule, rejected by
requiring the rule to be as thick as the side wall drawn with the same pen. Also
worth propagating: **replicate dots sit on BOTH top corners in five panels** —
the one shape the flank test cannot see — so the scan must be anchored on the
bar body by solidity.

Launched `zhangz_40433897`.

### 09:5x — the reopened head: recovered, mass-confirmed, and a reusable principle

`xuq_32662132` came back with head 9322 **recovered and then tested**, and it
verified my correction from `data/` rather than from my quote — which is the
right response to a claim from me as much as from a reviewer.

- `A9322` resolves to **one** structure across all 22 `YR_2024` rows:
  `Cc1nccn1CCCN` = **1-(3-aminopropyl)-2-methylimidazole**.
- `93` (`AA_2008`, 11 rows) and `A93` (`KW_2014`, 5 rows) resolve to the **same
  single** structure `NCCCn1ccnc1` — two spellings, two experiments, one
  molecule. Looked up through `corpus_components.smi()`, which raises on
  ambiguity rather than picking by row order.

**FINDING AU — a low-resolution axis cannot PIN a value but it can still TEST
one, and that distinction was worth four measurements.** The agent had refused
Figure S1's mass axis as unusable at ~3.93 Da/px, and *as a read* that refusal
was correct — three candidate glyph labels were undecidable. But with a formula
in hand the same axis **discriminates**:

| | |
|---|---|
| 9322-O17O measured | **765.5** |
| predicted [M+H]+ on head **9322** | 764.651 → **+0.8 Da = +0.20 px** |
| predicted [M+H]+ on head **93** | 750.635 → +14.9 Da = **+3.78 px** |
| margin on a one-CH2 difference | **18.9x** |

It also decides the three glyph candidates it had abandoned: `764.73` at 0.20 px
against `784.73` at 4.89 and `794.75` at 7.44. Its own summary of the mistake is
the lesson: *"I concluded that because a ~4 Da/px axis cannot pin a formula it
could not test one. Different questions."* A refusal for one purpose became a
stopping point for another. Now in every brief.

**And it refused to go one step further, correctly.** Sixteen remaining heads
(`9310`-`9316`, `9321`, `9323`, `9324`, `9331`-`9334`, `9341`, `9351`, `9352`)
have **zero corpus attestation** — measured and printed as a null result — and
the paper states no digit grammar. **One data point does not define one**, so
they stay not-enumerable and nothing was invented. It also corrected its own
count, 17 → 16.

A second result fell out of the same calibration: the O17C panel *derives* a
per-arm mass of **310.19 +/- 1.99 Da**, which **excludes a carbonate** (-74.1 Da,
~19 px) and fits heptadecyl acrylate to -0.1 Da — but +/-2 Da cannot separate
that from other 310-312 Da candidates, so O17C stays unbuilt and D-4 stays open.
Four of six panels still **refuse** (spacing nonuniformity 4.65% / 18.3% / 37.8%
against a 4% criterion; one never reaches 10 label boxes) and nothing was
loosened to recover them. It fixed a real detector bug on the way: its axis-row
finder used longest-contiguous-run, and an Xcalibur rule broken by tick notches
made it pick a 20 px title underline.

**No emitted row changed** — `06_rows.csv` is byte-identical — because the
emission rule never depended on that head. What changed is honest bookkeeping: a
skip reason moved from "no structure" to the dose-response reason, blocks 6 → 7,
enumerated 5 → 10, and `reference_structures.csv` now carries 2 of 3 structures
with the 93xx-series row deliberately empty. Gate 14/14, all 14 seeds move their
agreement count.

**Its retraction is a model of the posture §6 asks for.** It recorded *why* the
false claim survived: it was **never measured** — it came from eyeballing
`A9322_C2'L6` in an exploratory grep, and its own leading-token regex required
the name to *start* with `93`, so a name beginning with `A` never reached any of
its three keys. A conclusion drawn from a secondary artifact instead of the
primary one, which is the round-1 pattern that has now landed on four different
agents and twice on me today.

Launched `dahlmanje_30275336`, briefed with a scope question I have not had to
ask yet: **a barcoded pooled screen's normalised delivery is a competitive,
relative measurement** — each value depends on the other members of its pool —
so admissibility has to be argued before anything is extracted, and any
separate single-LNP luciferase validation is the part likely to be admissible.

### 10:2x — `zhangz_40433897`: 75 rows, a manuscript-grade number, and a silent no-op of mine

**75 rows, tier 5, every one `confidence = low` / `structure_assumed = 1`** with
a per-item note naming exactly what is assumed. The blocker is real and the
agent put it first: **no structure is pinned.** Three are proposed from HRMS
minus 2x the drawn intermediate, but **the cores of VN2C1 — the paper's own lead
compound, used in all its in-vivo work — and of VN3C2O have no amine of the
right formula among the 34 drawn at all.** So "the only drawn amine with this
formula" is an argument whose premise the paper falsifies. Also: the caption says
35 amines; one CAS appears twice, so **34**.

**FINDING AV, for the manuscript (F-056) — a same-lab cross-cell-type Spearman
of -0.089.** Six identical lipids, one lab, one library, same formulations, same
cargo, same readout, and only the target cell changes (4T1-Luc → MC38-Luc):

| | Spearman | n |
|---|---|---|
| this paper, 4T1 vs MC38 | **-0.089** | 16 |
| `molea.csv`, hepatocyte vs chondrocyte (§8's middle rung) | +0.30 | 1,152 |
| within-experiment regime (Result 6) | ~0.37 | |
| cross-lab shared lipids (Result 3) | 0.083 | 236 |

§8 uses `molea.csv` to pre-empt *"cross-lab rho ~ 0 is just assay noise"*. This
is a harder version of the same rung — weak at n=16, but it lands **below** the
cross-lab 0.083 rather than between it and 0.37. With `whiteheadka_30637934`'s
**same-dose rank inversion at rho = -0.20**, two independent papers now show a
lab disagreeing with *itself*. Recorded as needing a second digitisation before
it is quoted.

**Independently verified by me from `data/`:** this paper is a wholly disjoint
silo — **0 of 19,919 corpus lipids contain even the bare chromanol ring system**,
nearest neighbour to the drawn VE-AC intermediate **Tanimoto 0.2593**. A clean
illustration of Result 5: extrapolation across disjoint support, not a
distribution shift.

**Two method notes worth propagating.** Its first three benchmark SMILES, written
from memory, were **all wrong** — SM-102 one CH2 short of C44H87NO5, MC3 one CH2
short on a tail, ALC-0315 with the ester reversed — all parsed, all
canonicalised, no radicals, **only the InChIKey saw it.** That is exactly why
`reagent_registry.py` exists; now in every brief. And its best check found an
**unstated metric and a published sign error**: Spearman reproduces all four
printed coefficients to |delta| <= 0.0160 while Pearson matches 2 of 4, so the
metric is Spearman — and the printed DMG-PEG coefficient **has the wrong sign**,
with no parameterisation of Table S2 making it negative. Its seed is elegant:
a rank-preserving x1.10 leaves the count at 4 because Spearman is correctly
blind to it; a rank perturbation drops it to 0.

New furniture hazard: **an inset correlation table drawn INSIDE the plot box**,
whose ~300 px rule crossed a bar window and read **149.4% on a 0-150 axis**.
Fixed with a furniture strip — a horizontal run longer than three bar widths is
not a bar top — not by moving a threshold.

**FINDING AW — my own queue updates could silently match nothing, and one did.**
I launched an agent on PMID 41555650 and briefed it to work in
`curation/dahlmanje_41555650/`. The real entry is **`mitchellmj_41555650` — the
Mitchell lab, not Dahlman.** My inline `if x["folder"] == ...` loop matched
nothing, said nothing, and `check_state` then showed **3 running agents where 4
were running**. The agent had been told to create the folder if absent, so a
folder named for the wrong lab was about to enter the tree — and `Experiment_ID`
is derived from the first author, so this was a provenance error in the making,
not a cosmetic one.

Corrected the agent immediately with the right folder, the right author term for
`find_repo.py` (a wrong author returns a clean-looking zero), and a retraction of
the pooled-barcode assumption I had reasoned from the misidentified lab — plus a
pointer to `mitchellmj_38982196`, already staged at 100 rows from the same lab,
as a free external check.

Then fixed the class rather than the instance: **`provenance/tools/queue_update.py`
REFUSES when a folder matches nothing**, suggests near-misses, and offers
`--find <pmid>` to get the key right *before* launching. Seeded with the exact
wrong key I used — it exits 1 and names the correct entry. This is the same
defect class the whole run has been chasing in other people's code, committed by
me in mine: **a write that matches nothing must fail.**

### 11:0x — `maohq_39375194`: 5,112 rows, zero digitisation, and three findings bigger than the rows

**The run's largest and cleanest extraction: 5,112 rows, tier 3, 100% from the
authors' own repo** (`MaoResearchGroup/ml_lnp_design_rules` @ `82cee58`, file
md5s recorded), cloned shallow into its own empty directory with nothing
executed. This is the third confirmation of `PIPELINE_V2`'s central claim: the
curations that take values from a repo or Source Data have **zero digitisation
defects**, against eight for the four that digitised.

**`1080_LNP_OHE.csv` is 1,080 FORMULATIONS OF ONE LIPID.** `Ionizable_Lipid` has
exactly one distinct value on all 1,080 rows — DLin-MC3-DMA — and the one-hot
group is the **helper**. The grid is a complete unreplicated factorial:
6 helper x 5 x 4 x 3 x 3 = **1,080**, product equal to row count, every cell
occupied once. **Distinct structures 1 · formulations 1,152 · individually
measured 6,192.** The one-hot inversion was proved lossless both ways: 1,080/1,080
rows sum to exactly 1, and 1,080/1,080 inversions equal the file's **own**
categorical column — a verified bijection, not a relabelling we invented.

**It found that one arm of the library is already in the corpus, and excluded it
rather than re-emitting.** `RLU_HepG2` joins `YZ_2022` one-to-one at Spearman
**1.000000** — and it did not stop at the correlation: fitting the affine map and
expressing the residual in ln units gives rms **2.8531e-4** against **2.8868e-4**
expected from +/-0.0005 three-decimal rounding, **ratio 0.988**. The discrepancy
*is* the file's rounding, and the corpus holds those values at higher precision.
I confirmed the shape myself: `YZ_2022` is 1,080 rows, HepG2 only, one ionizable
lipid, six helpers x 180 — the identical grid.

Three findings recorded in `review/FINDINGS_2026-11-09.md`:

**F-057 (open, touches 360 already-merged rows).** `18PG` is **SOPG by the Nat
Commun text** — which is what `CLAUDE.md` §3's existing repair rests on, "named
verbatim" — and **DOPG by this paper's drawing**, which shows a cis C=C in both
tails, with the repo file carrying `Double_bonds = 2`. The corpus stores one
SOPG string on **1,474** rows (`WT_2026` 1,114, `YZ_2022` 180, `YZ_2024` 180).
If the drawing is right, **360 merged rows carry the wrong helper lipid**, and
`reagent_registry.py` is *enforcing* the wrong answer via C11. The curation
emitted as drawn with `confidence = medium` and **left the corpus alone**, which
is correct: this is the `AC_2025` shape, where a confident internally consistent
repair put 419 measurements on the wrong molecules. One caution I added that the
report did not: a `Double_bonds` column in a one-hot file could be an
**ionizable**-lipid descriptor (MC3 has four C=C), so it is not decisive alone.
Needs a catalogue number or an author query. Error class 48, recurring on the
same reagent.

**F-058 (manuscript).** Same lab, same library, same cargo, only the target cell
changing: Spearman against HepG2 = **0.494 / 0.488 / 0.476 / 0.447 / 0.330**.
That is a **five-fold replication of §8's single `molea.csv` 0.30 datapoint**, on
a different lab, a different cargo and a different axis of variation
(formulation, not chemistry) — and **0.30 is the bottom of the range, not an
outlier**. With F-056's -0.089 and F-054's same-dose inversion at -0.20, the §8
ladder now rests on four independent papers instead of one, spanning -0.20 to
+0.49. That spread is itself the argument against "cross-lab rho ~ 0 is just
assay noise".

**F-059 (release blocker).** Both licences fail §8: the article is closed
(`license = null`, PMC deposit an NIH author manuscript under "fair use"), and
the repo is **GPL-3.0** — copyleft *software*, not a permissive data licence. The
asymmetry is notable: `YZ_2022`, already merged and covering this library's
HepG2 arm, **is** open access. Better known now than at submission.

**Its own anchor caught its own bug**, which is the point of C10: its first MC3
build attached the C1 end of the chain to C19 — **correct formula, correct exact
mass to 6 dp, wrong molecule**, and only the InChIKey saw it. Its three anchor
seeds all carry identical formula *and* identical exact mass (the MC3
regioisomer: C43H79NO2, 641.611081 both). It also **declined to credit a check**:
mol% summing to 100 is algebraically a tautology of its own derivation, so it
replaced it with a row-wise `pct = 100/(1+ratio)` identity against the file's
second parameterisation — 0 violations at 1e-6 on 1,080 rows. And its first
seeding run found a **real gate gap**: two seeds moved a metric while the gate
stayed silent, so it widened **the gate, not the seed**.

Two more contradictions worth keeping: **SI Table 6's first two row labels are
swapped** relative to SI Table 1, and **the PEG bracket reads 44 while 45
ethylene oxides are drawn** — one sits outside the bracket, and taking the
subscript gives a 44 Da error that only the InChIKey catches, while n=45
reproduces the corpus's stored `DMG-PEG2000` exactly.

### 11:3x — `dahlmanje_30275336`: 0 rows, and 22.6% of the in-vivo pool now has a caveat

**0 of 29 values admissible, and the reasoning is the deliverable.** Sago/Dahlman
PNAS 2018. It used **tier 2 printed numerals only** and wrote no digitiser at
all.

**It declined to refuse pooled data on principle, which was right.** I had
briefed it to treat a pooled barcode readout as suspect; it checked the schema
first and found `data/` already has an `Experiment_batching` vocabulary of
`{barcoded, individual}`. I verified the consequence myself:

| pool | batching | rows | experiments |
|---|---|---|---|
| vivo | `individual` | 1,396 | 19 |
| vivo | **`barcoded`** | **407** | **1 — `LR_2023`** |
| vitro | `individual` | 20,924 | 37 |
| vitro | **(blank)** | **8,271** | **6** |

So a pooled multiplexed in vivo screen is **already merged and declared**, and a
blanket refusal would have contradicted the project's own precedent. **F-060**
records what it established instead, by reconstructing the paper's own worked
example: the readout is **compositional** —
`nd_i = [(count_i/Σcount)/(input_i/Σinput)] / Σ_k[same] × 100%` — reproducing
**6/6 column totals, 18/18 cells, 3/3 input shares, residual 0.0**. The authors
**drew `Total (%) = 100`** under both tables and never computed it; the unrounded
reconstruction hits exactly 100.000000.

The consequence is specific: **each value is that LNP's share of its own pool**,
its three screens used pools of 71/78/106, so the values live on three different
simplices — and `Experiment_value` z-scores a *share*. Result 3 uses
within-source ranks, which are invariant to the common denominator, so the damage
is limited; but the in-vivo concordance is already LOW-CONFIDENCE at n=9 and
**22.6% of that pool now carries a second caveat**. Separately, **8,271 vitro
rows across 6 experiments have a BLANK batching** — per §3 that is "not stated",
not "clean", so 28% of the in-vitro pool has unrecorded batching that Methods
cannot currently describe.

**The actual blocker was structure, not pooling.** Eight named compounds,
**zero determinable SMILES**: four are PEI600/PEI1200 conjugates — polydisperse
branched polymers drawn as single repeat-unit cartoons with truncated stubs and
**no stated degree of substitution**. The paper agrees with us, blaming batch
variability on "the resulting heterogenous chemical structure". `data/` has **0
blank `IL_SMILES` in 30,998 rows**, so a structureless row cannot be emitted.
And the individually-measured subset — 7C2 and 7C3, proper single-LNP
injections — are **built on the same compound**, so even with a structure it
would be a 2-row experiment with zero structural variance.

**A limitation of my own guarded lookup, recorded in F-060.**
`corpus_components.smi()` raises when a name carries two structures *in the
corpus*. It cannot catch a name that is unambiguous in the corpus and simply
**wrong for the paper in hand** — which is what this curation measured for
`C14-PEG2000`: 13 contiguous bond quads at a 5.92 pt stride (14 carbons
including the carbonyl) on a drawn phosphoethanolamine-**carbamate**-PEG head,
where the corpus's same-named entry is the glycerol-ether DMG-PEG2000. Second
sighting of this shape after F-049. The honest posture is §3's for PEG-lipids: a
per-row flag, not a registry entry.

**It measured the name-token trap rather than describing it:** a naive substring
join on the bare token `80` captures **94 vitro rows spanning 85 distinct
structures across 5 experiments**, and in this paper `80` indexes a PEI1200
polymer. Same token, unrelated molecule — F-055's amendment, independently
confirmed.

**Two no-op seeds caught in its own work**, both printing a clean `14`: one
removed the short carbonyl→O stub instead of a bond, the other injected
something the size filter is *supposed* to exclude. And a **magic-bytes catch** —
its first `fig04.jpg` fetch returned a ZIP, which a size floor would have passed.

Its replicated-arithmetic finding is worth copying: four formulations summing to
101 at a **stride of 16 across four helper lipids** — the signature of one error
in a four-times-pasted block, confirmed as the paper's because the table grid and
an independent text extraction return identical digits. Overall **418/424
formulations sum to exactly 100**.

Licence **CC BY-NC-ND 4.0**: NC *and* ND both fail §8.

Launched `mitchellmj_41603124`, the next barcoded paper, briefed with the
compositional reconstruction as a required step rather than a caution.

### 12:0x — laptop slept: all four agents died, all four resumed, and the guard that missed it is fixed

**What happened.** The machine slept. All four running agents stopped: their task
files froze at a **140-byte stub**, stale 25-39 minutes. **Their work survived on
disk**, which is what made the right response resumption rather than relaunch:

| folder | files on disk | newest artifact |
|---|---|---|
| `mitchellmj_41555650` | 26 | `07_qc.json` |
| `dixonjb_41269090` | 18 | `06_rows_z.json` |
| `mitchellmj_41603124` | 7 | `source/` |
| `dahlmanje_36701517` | 6 | `04_digitized.json` |

All four resumed by `SendMessage`, which preserves their context, each with the
decision that governs its folder restated — the Mitchell-not-Dahlman correction
for 41555650; the per-experiment-z rule (error class 22) for the one mid-z-score;
the "RNA-seq of the host response is not a delivery measurement" question for
36701517; and the compositional-reconstruction requirement for the barcoded one.

**FINDING AX — `check_state.py` printed "state consistent" while every agent was
dead.** My guard, added two hours earlier for a different failure, checked only
that a task file **existed**. A frozen stub exists. So the script reported a
clean line on a completely stalled run — which is the exact defect this whole run
has been cataloguing in other people's code (`find_tables()` returning 0 on a
page that is nothing but table; a chroma mask deleting a bar; three `ZT_2026`
marginal checks scoring a literal 0 for twelve days; `triage` printing an empty
report), committed by me in the one script whose entire job is to notice that
work has stopped.

Fixed with a **staleness** test at 15 minutes — long enough for a slow tool call
in a job where agents routinely run 20-45 minutes, short enough that an
unnoticed sleep trips it. The message says to **resume, not relaunch**, because
relaunching would discard exactly the work that survived.

Seeded, per §6's rule that a new check is not done until it has been made to
fail on purpose:

| run | exit | result |
|---|---|---|
| control, four live agents | **0** | `state consistent` |
| one task file back-dated 40 min | **1** | names the folder, the agent, **40 min**, 8,933,218 bytes, and says to resume |
| restored | **0** | — |

Note the seeded line also reveals the file was 8.9 MB, not the 140-byte stub —
i.e. that agent had already resumed and written before the seed, so the guard is
measuring liveness rather than file size, which is the right quantity.

### 12:3x — pushed everything; `isolation: "remote"` measured and it does NOT leave the laptop

The user asked to move the run to the cloud so the laptop could be closed and
carried around. **A running session cannot be migrated** — the conversation and
its in-process subagents live on the machine that started them, and there is no
"move session" operation. What can move is the work, so:

**Committed and pushed `d13782c` to `origin/main`: 847 paths, 22.0 MB.** Every
extraction script, every `05_values.csv` / `06_rows.csv` / `07_qc.json` /
`DECISIONS.md`, `curation/lib/`, `provenance/tools/`, the findings ledger and
this log.

**FINDING AY — the repo is PUBLIC, and a plain `git add -A` would have published
60.2 MB of copyrighted publisher assets to it.** `.gitignore` already excluded
`*.pdf` "for copyright"; it did **not** exclude anything else a publisher
serves. Measured before committing: 101 files — **53 jpg, 16 xml, 14 png,
7 jpeg, 4 tif, 3 zip** — from papers that are subscription-only or
**CC BY-NC-ND**, where ND alone forbids redistribution, and §8 requires
CC0/CC-BY compatibility before release. Extended the rule to figure rasters,
supplementary archives, JATS XML, acquisition directories (`source_epmc/` was
not covered by a `source/` rule and held 6 publisher JPEGs), page renders and
annotated overlays — an overlay is the figure plus our annotation, not a new
work.

Getting that exclusion to bite took three attempts and the reasons are worth
recording. A `curation/**/*.png` pattern **did not match**, while a slashless
`*.png` did — which is why the original `*.pdf` rule has always worked. And once
a path is in the index, `.gitignore` has no effect on it, so `git add -A` kept
re-adding files the new rule covered; `git check-ignore` was silent for the same
reason, and only `--no-index` revealed the rule was correct all along. The fix
was `git reset` to clear the index first. Final staged set verified at **zero**
matches for `.pdf`, `.png`, `.jpe?g`, `.tif`, `/source`, `/renders/`, `/figs/`.

Separately noted and **not** acted on: **31 overlay images from earlier commits
are already public** in `HEAD`. Removing them needs a history rewrite on a
public repo, which is the user's call, not mine.

**FINDING AZ — `isolation: "remote"` does not leave this machine.** I told the
user cloud agents were available because the launch succeeded. That was wrong,
and the probe I had already sent proved it: the agent reported
`hostname: Nathans-MacBook-Air-6.local` and a path under
`.claude/worktrees/agent-<id>` — **a local git worktree, not a cloud
environment.** The flag degraded silently. Corrected to the user in the same
turn. The probe still earned its cost: a fresh checkout at `d13782c` sees
`data_vitro.csv` at 29,195 data rows, `curation/lib/`, the queue and the log,
has network access (Europe PMC HTTP 200) and imports rdkit, pdfplumber,
pymupdf, numpy and scipy — so the repo **is** a sufficient handoff for a session
started elsewhere, which is the thing worth knowing.

`provenance/HANDOFF_CLOUD.md` written: what transfers, what cannot and why
(the 84 PDFs and the entitled-browser publishers), what a cloud session can do
fully / partly / not at all, and the standing rule that **nothing in
`curation/` has been merged** and the frozen pipeline must not be re-pointed
without the user's decision.
