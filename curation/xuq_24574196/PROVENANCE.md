# PROVENANCE — curation/xuq_24574196

Wang M, Alberti K, Varone A, Pouli D, Georgakoudi I, Xu Q.
**Enhanced Intracellular siRNA Delivery Using Bioreducible Lipid-Like Nanoparticles.**
*Adv. Healthcare Mater.* 2014. PMID **24574196**, DOI 10.1002/adhm.201400039.

**ENDPOINT: siRNA-mediated gene knockdown** — GFP knockdown read by flow
cytometry, and Plk-1 knockdown read by Western blot. That is a nucleic-acid
delivery endpoint and is admissible in kind. The paper is original research,
not a commentary: `abstractText` is 493 characters and the SI carries a full
Materials and Methods with cell densities, doses, N/P ratios and timepoints.

**ROWS EMITTED: 0.** The Supporting Information is all this folder holds, and
it contains no delivery measurement. Every knockdown panel is in the main
article, which is closed access and was not acquired. This is an honest zero,
not an unfinished extraction.

---

## 1. What the SI actually contains

| SI locus | content | admissible |
|---|---|---|
| pp. 2–6 | Materials and Methods | metadata only |
| p. 3 | ¹H NMR + ESI-MS of **1-O16B** and **1-O16** | structural anchors |
| p. 7, Fig S1 | the two NMR spectra as rasters | no |
| p. 8, Fig S2 | TEM, before/after GSH, no numeric scale | no |
| p. 9, Fig S3 | AlamarBlue **cell viability**, 6 bars, 1–32 nM | **no — viability is not a delivery endpoint** |

Fig S3 was deliberately **not digitised**. `data/`'s `Experiment_method`
vocabulary is delivery-only; digitising a viability panel would manufacture
rows the corpus cannot accept.

## 2. Acquisition, and why the values are missing

Measured, not assumed:

- NCBI `elink` pubmed→pmc returns **`pubmed_pmc_refs` only** (31 citing
  papers), **no `pubmed_pmc` self-link** → no PMC deposit, no PMCID.
- Europe PMC: `inEPMC=N`, `isOpenAccess=N`, `hasSuppl=N`; `fullTextXML`,
  `supplementaryFiles` and `fullTextRepo` all **404**.
- OpenAlex: `is_oa=false`, `oa_status="closed"`,
  `any_repository_has_fulltext=false`, `oa_url=null`.
- Wiley landing page via plain `curl`: HTTP 200 serving a Cloudflare
  interstitial (`<title>Just a moment...</title>`), **0 assets enumerable**.
- `find_repo.py` on DOI + PMID + author: **0 relevant repositories**.

Browser acquisition was out of scope for this run. **No paywall or bot-check
bypass was attempted.**

## 3. Structures — what is pinned and what is not

Four lipidoids enumerated in `03_enumerated.csv`, built through
`curation/lib/assemble` (atom-mapped `attach` + `esterify`, no string surgery,
no `RemoveAtom` near a stereocentre — there are no stereocentres in this
family at all).

The obvious build order was **refused by the library and the refusal was
respected**: esterifying the mapped `[CH2:2]` fragment first leaves it
valence-incomplete, and `check_no_radicals` fired with
`1 radical atom(s): [(0, 'C', 1)]`. The route was reordered
(attach → esterify → attach → esterify) so the mapped carbon is always an
`attach` attachment atom. Nothing was loosened.

| lipidoid | status | confidence |
|---|---|---|
| **1-O16B** | pinned by `C39H78N2O4S4` **and** the 40 H bulk integral | high |
| **1-O16** | pinned by `C43H86N2O4` **and** the 56 H bulk integral | high |
| 1-O14B | name only; S-alkyl extrapolated | medium, `structure_assumed=1` |
| 1-O18B | name only; S-alkyl extrapolated | medium, `structure_assumed=1` |

**`O14` and `O18` are not enumerated** — the nonbioreducible C14/C18
counterparts are never named anywhere in this SI, and writing them would be
inventing library members. **Amine 2 is not enumerated** — SI p.2 says "amine 1
or 2" and never says what amine 2 is.

### The two traps this paper sets, both measured

**A name token is not a chain length.** `O16B` is **not** a C16 alkyl. Its
S-alkyl is **dodecyl**, and two independent printed facts say so: the formula
`C39H78N2O4S4` (a C16 alkyl gives `C47H94N2O4S4`, eight carbons too many) and
the bulk-methylene integral of **40 H** (a C16 alkyl predicts 56 H). The token
counts **chain heavy atoms**: 2 C + 2 S + 12 C = 16 — which is also the only
reading under which the paper's own claim that O16B and O16 have "similar
chemical structures" is true.

**A leading number is not an identity.** The `1` in `1-O16` is a local label.
It is **not** the corpus's amine 1: `AA_2008`'s `1_O15` is
`C40H79NO5` — one nitrogen, five oxygens — against the SI's printed
`C43H86N2O4`. Keying on the number token would have put every value in this
paper on a molecule the paper's own mass spectrum excludes.

## 4. Corpus overlap

Four keys over 30,998 corpus rows / 19,919 unique `IL_SMILES`.

- **Exact full-InChIKey: 0 of 4.** Skeleton InChIKey: **0 of 4.**
- **Near miss, and it is the valuable result.** `AA_2008` (also in `KW_2014`)
  carries an `80_O<n>` series — `80_O10/12/13/15/18` — whose head is
  **InChIKey-identical to ours** after both are truncated to the O-methyl
  ester (3 of 3 probes). So **Xu's "amine 1" is the Anderson/Akinc catalogue's
  amine 80**, N,N-dimethyl-1,3-propanediamine. That head was derived here from
  NMR and a printed formula and is corroborated by a structure an independent
  curator drew years earlier from a different paper. `80_O16` is **absent** from
  the corpus, so `1-O16` really is new; its nearest neighbours are `80_O15`
  (−1 CH₂ per tail) and `80_O18` (+2).
- Same amine again in **`QL_2026`** (`80-L2-C*`, head `C(=O)NCCCN(C)C`), from
  this same lab — a third, independent sighting of catalogue number 80.
- The disulfide **arm** is new: 410 corpus structures contain a disulfide and
  98 contain a `C(=O)OCCSS` motif, but those are `ZC_2023`'s symmetric
  bis-ester linkers (`…OCCSSCCOC(=O)…`), not an alkyl-terminated
  `C(=O)OCCSS-C₁₂H₂₅`.

Zero exact overlap **costs** something: no corpus value cross-checks ours, so
the mass, NMR and anchor evidence carries the whole external burden. The head
identity above is the one place that burden is shared.

## 5. Checks — `07_qc.json`

Every check reports its agreement count **and** the count it scores on a
seeded defect. A check is credited only where the two differ.

| check | real | seeded | teeth |
|---|---|---|---|
| P1 built formula vs printed | **2/2** | 1 | yes |
| P2 **unprinted**: Σ printed integrals = formula H (78, 86) | **2/2** | 1 | yes |
| P3 bulk-CH₂ integral predicted from structure (40, 56) | **2/2** | 1 | yes |
| P4 ESI-MS mass convention | **2/2** | 0 | yes |
| P5 **unprinted**: name token = chain heavy atoms | **4/4** | 3 | yes |
| P6 **external**: head ≡ corpus amine 80 | **4/4** | 2 | yes |
| P7 name-token trap: corpus amine 1 excluded by formula | **1/1** | 0 | yes |
| P8 **unprinted**: `80_O*` series steps are exact CₙH₂ₙ | **4/4** | 1 | yes |
| P10 anchor vs enumeration, full InChIKey | **2/2** | 0 | yes |
| P9 parse / canonical round-trip / radicals / stereo | **4/4** | n/a | see below |

P9's seeded failure is historical rather than synthetic: `curation/lib`
refused a radical at build time on this paper's first run.

**P4 — which mass convention.** The SI prints no `Calc.` column header, so the
convention was determined from the data. `Calcd` is the **neutral average
mass**: residuals **0.031** and **0.029** Da, against 0.19/0.47 for
[M+H]⁺ monoisotopic, 0.54/0.82 for neutral monoisotopic and 1.04/0.98 for
[M+H]⁺ average. The margin over the next-best candidate is 0.16 and 0.44 Da.
But an ESI-MS peak is necessarily an **ion**, so `Calcd` and `Found` are not
the same quantity (`Found` − neutral average = +0.369 and +0.529 Da). **This
instrument confirms the molecular formula at roughly unit resolution and
nothing finer.** No sub-Dalton claim may be made from it.

**My own parser bug, caught and fixed.** P10's independence test first scanned
`make_reference.py` for the substrings `import assemble` / `esterify` and
reported the module unclean — because those words appear in its *docstring*,
where they are prose about what it deliberately does not use. Rewritten to
parse the import statements with `ast`. A substring scan cannot tell a module
name from a sentence about one.

## 6. Anchor posture (recorded by hand — C10/C16/C18 were not run)

`provenance/structure_check.py`, `check_anchor_independence.py` and
`check_anchor_coverage.py` write tracked files other agents read and were
**not** run, per this run's brief.

- `reference_structures.csv` holds **2 anchors** (`1-O16B`, `1-O16`) and **2
  deliberately empty rows** (`1-O14B`, `1-O18B`), so the absence is visible.
- **Independence is claimed by provenance, not by bytes.** `make_reference.py`
  imports `csv`, `os`, `sys` and RDKit and **nothing else** — verified by an
  `ast` walk of its import statements, not by grep. It imports no
  `curation/lib` module and not the enumerator. Every fragment is a source
  string literal carrying the NMR shift and integral that fixes it, and the
  traversal is **head-first**, the reverse of the enumerator's.
- Both anchors are **byte-identical** to the enumeration after
  canonicalisation. **That is expected and is not evidence of copying** — a
  correct re-derivation canonicalises to the same string as a correct
  enumeration, which is exactly why C10's own `anchor_may_be_a_copy` missed
  F-037. InChIKeys: `KEEDANUIQHUIOP-UHFFFAOYSA-N` (1-O16B),
  `ULZCWFRFMLJWOX-UHFFFAOYSA-N` (1-O16).
- **The gate was proved to fire.** `python3 -I make_reference.py` exits **0**;
  `python3 -I make_reference.py --seed`, which reads the token `16` as a C16
  S-alkyl, exits **1** with
  `ANCHOR GATE FAILED: 1-O16B: anchor formula C47H94N2O4S4 != printed C39H78N2O4S4`.
  Exit codes read from `returncode`, not from a pipe.
- **Weaker ground truth than usual, stated plainly:** this SI draws **no
  chemical structure anywhere**. The anchor rests on the printed ¹H NMR
  assignment list, which is prose the authors wrote, not a rendering they had
  to make. See `DECISIONS.md` D-5.

## 7. Licence

The SI carries *"Copyright WILEY-VCH Verlag GmbH & Co. KGaA, 69469 Weinheim,
Germany, 2014."* and OpenAlex reports the article `closed`. **This source is
not CC0/CC-BY-compatible**, which `CLAUDE.md` §8 requires before corpus
release. Moot while this folder emits zero rows; blocking the moment it does not.

## 8. Not merged

`data/` was not touched. `DECISIONS.md` carries five `OPEN` entries, so the
merge gate will refuse this folder — correctly.

## 9. What could NOT be verified

- Amine 2, and therefore every lipidoid built on it.
- `O14B`/`O18B` S-alkyl lengths — a one-point extrapolation. The primary
  reference (Ryu *et al.*, *Langmuir* 2010, **26**, 7086) would settle it and
  was not acquired.
- The N/P ratio of the Fig S3 viability experiment: 2:1 (Methods) vs 5:1
  (caption), with no third statement.
- The library's true size. The SI never says how many lipidoids were made.
  **Distinct structures: 4 enumerable** (of an unstated total).
  **Formulations: not recoverable per lipid** — the screen used bare crude
  lipidoid, the imaging work a 16:4:1 w/w lipidoid/cholesterol/DOPE film, and
  no per-lipid formulation table exists anywhere in the SI.
  **Individually measured: unknown, and 0 obtainable.**
