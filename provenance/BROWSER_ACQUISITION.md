# Browser acquisition — measured publisher behaviour, 2026-11-09

Tested with the Playwright MCP browser against the publishers the corpus
actually draws on. **The headline is not what item 9 of the planning document
(superseded by `PIPELINE_V2.md`)
assumed.**

## The assumption that was wrong

The plan said the browser's job was to get past **institutional SSO** to reach
paywalled supplementary material. Measured, that is not the binding constraint.
The binding constraint is **bot-blocking**, and it bites on the *article page*
long before any paywall does.

Two separate things were conflated and they behave differently:

| | |
|---|---|
| **Bot-blocking** | Cloudflare / Akamai refusing a non-browser client. Hits `curl` on 6 of 8 endpoints tested. Nothing to do with subscriptions. |
| **Paywalling** | the subscription question. Only matters once the page loads. |

## The matrix

`curl` used a current desktop Chrome user-agent. Playwright is the MCP browser,
not logged in to anything.

| publisher | corpus papers | `curl` | Playwright | note |
|---|---|---|---|---|
| Springer Nature | **71** | **200** | 200 | no block, no login |
| `static-content.springer.com` (Source Data) | — | **200** | 200 | **open to a plain request** |
| Elsevier / ScienceDirect | 19 | 403 | **200 after a homepage visit** | login needed to ENUMERATE; see below |
| `ars.els-cdn.com` (Elsevier assets) | — | **200** | 200 | **open — no session needed** |
| Wiley | 12 | 403 | **200** | browser wins |
| PNAS | 9 | 403 | **200** | browser wins |
| ACS | 6 | 403 | **200** | resolved on the 2nd navigation |
| RSC | 4 | 403 | **200** | Cloudflare challenge self-resolved in ~8 s |
| bioRxiv | 2 | 429 | not tested | rate-limited |
| PubMed (`pubmed.ncbi.nlm.nih.gov`) | — | — | **reCAPTCHA** | use the E-utilities API instead |

## What this changes

1. **The browser earns its place, but for a different reason.** It unblocks
   Wiley, PNAS, ACS and RSC — 31 papers' worth of publisher — where a plain
   request gets 403. That is the actual value, and no login was involved in any
   of it.

2. **For Springer Nature, which is 71 of the corpus, the browser is only a
   URL enumerator.** The article page and every `static-content.springer.com`
   asset return 200 to a plain request. The browser is needed to *read the
   anchor hrefs off the page*; the download itself needs nothing. So the
   autonomous stage-00 path for the largest publisher is: resolve the DOI in
   the browser, enumerate `MOESM` links, then fetch with an ordinary client.

3. **Elsevier is not a hard block — I was wrong about this, and the correction
   matters.** My first conclusion was that ScienceDirect refuses the automated
   browser outright, because the article URL returned `403 Just a moment...`
   and did not yield to a 20-second wait. That was an artifact of **how I hit
   it**: a cold request straight to an article route, with no session.

   Visiting `https://www.sciencedirect.com/` first returns 200 and sets a
   session; the *same* article URL then loads immediately. So the rule for a
   fetcher is **warm up on the publisher's homepage, then request the article**
   — and the general lesson is that a 403 on a deep link is not evidence about
   the site, only about the request. Both of my Elsevier attempts had been deep
   links.

   What Elsevier *does* require the institutional login for is the content
   itself: the article resolves to `/science/article/**abs**/pii/...`, offers
   "Purchase PDF", and **exposes no supplementary links at all** to an
   unentitled session. On Nature and Wiley the supplementary assets are listed
   whether or not you are entitled; on Elsevier they are not, so without a
   login you cannot even enumerate what exists. 19 papers.

4. **Never scrape PubMed.** It serves a reCAPTCHA to the automated browser. The
   repo already does the right thing with `provenance/fetch_esummary.py` over
   the E-utilities API, and that must stay the route for PMIDs and DOIs.

5. **Cloudflare needs a wait, and the wait is the whole trick.** RSC went from
   `403 Just a moment...` to the real article in ~8 s with no interaction, and
   ACS resolved on a second navigation. A stage-00 fetcher that treats the
   first 403 as final will report "paywalled" for sites that would have let it
   in. Navigate, wait, re-read, and only then conclude.

## The login route, mapped as far as it goes without credentials

For `curation/jl_2026` (PII `S1369702126001112`, *Materials Today*):

1. `https://www.sciencedirect.com/` — warm up the session.
2. The article URL — loads, `/abs/`, "Purchase PDF", no supplementary links.
3. "Access through your organization" ->
   `sciencedirect.com/user/institution/login?targetUrl=...`
4. Elsevier's "Find your organization" -> type the **organization name**
   (never the user's email, which is private) -> "Case Western Reserve
   University Libraries" is the first match.
5. -> `login.openathens.net/auth/case.edu/...` -> "Sign in at Case Western
   Reserve University"
6. -> `https://login.case.edu/cas/login?...` — the CWRU CAS form, CWRU ID +
   passphrase. **This is where automation stops.**

Two cookie-consent overlays intercept clicks on the way and must be dismissed
first: Elsevier's OneTrust banner (no "reject all" — decline Performance
Cookies in the preference centre, then "Confirm my choices") and OpenAthens'
own (has "Accept necessary cookies"). Both were set to the minimal option.

Note the CAS page's own banner: *"This URL must BEGIN WITH EXACTLY:
https://login.case.edu/ If it does NOT, do not enter anything."* Check that
before typing, every time.

## What the institutional login actually bought — tested 2026-11-09

Signed in through CWRU, on `curation/jl_2026` (PII `S1369702126001112`):

| | unentitled | entitled |
|---|---|---|
| article URL | `/science/article/**abs**/pii/...` | `/science/article/pii/...` |
| header | "Purchase PDF" | "View PDF", *"Brought to you by: Case Western Reserve University Libraries"* |
| supplementary links | **none listed at all** | Appendix A, "Supplementary Data 1" |

Then the part that matters operationally. Having obtained the asset URLs,
**both the supplementary PDF and the figure images download with NO session**:

```
ars.els-cdn.com/.../S1369702126001112-mmc1.pdf       200   3,265,424 bytes
ars.els-cdn.com/.../S1369702126001112-gr2_lrg.jpg    200     628,906 bytes
sciencedirect.com/.../pdfft?md5=...-main.pdf         403   (needs the session)
```

So on Elsevier the login is required to **enumerate** the supplementary assets
and to fetch the **article PDF**, and not to fetch the supplementary files
themselves. That is the opposite shape from what "paywall" suggests, and it
means a signed-in enumeration pass can hand its URL list to an ordinary client.

**For this particular paper the login found nothing new**: the sole
supplementary item is `mmc1.pdf`, and `curation/jl_2026/` already holds it at
exactly 3,265,424 bytes. Worth stating plainly — the honest result of the first
real Elsevier login is that there was nothing more to get.

### The genuinely useful find: publisher figure images, open on the CDN

Elsevier publishes each figure as a standalone high-resolution JPEG at a
predictable URL:

```
https://ars.els-cdn.com/content/image/1-s2.0-<PII>-gr<N>_lrg.jpg
```

Measured on this paper, against the alternative of rendering the PDF page:

| | pixels |
|---|---|
| `gr1_lrg.jpg`, `gr3_lrg.jpg` (the figure alone) | **3398 x 4336**, **3227 x 4336** |
| `gr2_lrg.jpg`, `gr5_lrg.jpg` | **3549 x 3530**, 3549 x 3534 |
| a 600 dpi render of the whole A4 page | 4961 x 6615, of which a half-page figure is ~2500 x 3300 |

So the CDN figure is roughly 1.4x the linear resolution our 600 dpi render
gives for the same figure, and it is the publisher's own raster rather than a
re-rasterisation of an embedded image. `CLAUDE.md` section 5 says "resolution
matters more than cleverness: glyph clustering that was hopeless at 300 dpi
became exact at 600." This is the next step up that same ladder, it is free,
and it needs no login.

**Check these URLs before rendering any Elsevier figure.** The same should be
checked for other publishers; Nature's equivalents are under
`media.springernature.com`.

## Operational rules, unchanged from the plan

These still hold and are now load-bearing rather than precautionary:

- **Credentials are never handled here.** When Elsevier is attempted, the user
  signs in by hand once into a persistent profile directory; scripts drive that
  profile and read cookies they never wrote. No password in a file, an env var,
  a prompt, a script or a log. On expiry, fail with "log in again".
- **The profile directory lives outside the repo** and is gitignored by path.
- **One DOI at a time, on request.** No crawling, no citation-graph following,
  no speculative prefetch, no retry storms. Automated bulk retrieval is what
  gets an institution's whole IP range blocked, and the subscription is the
  university's.
- **Downloads are untrusted data**: own new empty directory, scripts elsewhere,
  paths passed as arguments, readers run under `python3 -I`.
- **Record the acquisition**, with URL, sha256, bytes, content type and UTC
  timestamp — see `curation/amino_lipid_fg_2026/00_acquisition.json`.
- **Record the licence at acquisition time.** `CLAUDE.md` §8 needs every
  constituent source CC0/CC-BY-compatible before release. Not yet done for the
  FG_2026 files, and noted as outstanding in that file.

## Housekeeping

The Playwright MCP server writes page snapshots and console logs into
`.playwright-mcp/` **in whatever directory it was launched from**, which is the
repo root. Added to `.gitignore` 2026-11-09. They are session artifacts and
nothing reads them.

## First real use: the FG_2026 blocker

`curation/amino_lipid_fg_2026/07_qc.json` has recorded a BLOCKER since the
folder was created: the 960 per-lipid values, the 16 aldehyde tails and the 2
bis-isonitrile cores are in neither uploaded PDF, and the paper says "Source
data are provided with this paper".

All five Source Data workbooks were retrieved (10.1038/s41563-026-02555-0,
PMID 41922838). **The blocker is only partially relieved and must not be
recorded as resolved.**

- **Supplied:** Fig. 1d gives all **30** amino-acid head-group blocks with
  normalised values (Cha = 1.000, Asn = 0.404, Lys = 0.294, …). Fig. 1e gives
  **20** named members with values and confirms the naming scheme is
  `<core><AminoAcid>-<tail index>` over cores `CH` and `LH`.
- **Not supplied:** the 960. Only **20** distinct member names appear anywhere
  across all five workbooks — 2 cores, 15 of 30 amino acids, 9 of 16 tail
  indices. 1,615 numeric cells in total; the largest single sheet is Fig. 5h at
  608. The 16 aldehydes and the 2 cores are absent.

So this folder can now do value-level validation on 20 of 960 members plus the
full block ranking, which checks the naming scheme and the block ordering. It
still cannot do the 960-row recall comparison it exists for. Detail and hashes
in `curation/amino_lipid_fg_2026/00_acquisition.json`.
