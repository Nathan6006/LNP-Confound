# PROVENANCE — peerd_31999380

Ramishetti S., Hazan-Halevy I., Palakuri R., Chatterjee S., Gonna S.N.,
Dammes N., Freilich I., Kolik Shmuel L., Danino D., **Peer D.**
*A Combinatorial Library of Lipid Nanoparticles for RNA Delivery to Leukocytes.*
**Adv. Mater. 32, e1906128 (2020 Mar).** PMID **31999380**,
DOI **10.1002/adma.201906128**.

Curated 2026-11-10. **Nothing in this folder has been merged into `data/`.**
`DECISIONS.md` carries eight `OPEN` entries, and `build_dataset.py` refuses a
folder in that state.

---

## 1. The title and the endpoint, verified before anything else

The briefing title was truncated at *"…for RNA Delivery…"*. The hidden tail is
**"to Leukocytes"**, and it is the endpoint qualifier: the biological systems
are primary mouse T lymphocytes and the human multiple-myeloma line **U266** —
not hepatocytes.

Two sources, both checked:

- **SI p.1**, the Wiley cover sheet: *"A Combinatorial Library of Lipid
  Nanoparticles for RNA Delivery to Leukocytes"*.
- **NCBI esummary**, `db=pubmed&id=31999380`, retrieved 2026-11-10, HTTP 200:
  `"A Combinatorial Library of Lipid Nanoparticles for RNA Delivery to
  Leukocytes."`

**SI p.2** — the author manuscript header inside the same PDF — says **"RNAi
Delivery"** instead. Recorded as contradiction **C-07**; the published form
wins.

**THE ENDPOINT IS siRNA DELIVERY TO LEUKOCYTES**, read out in the main text as
gene silencing. The cargo is siRNA throughout: siPLK1 (therapeutic),
siNC / siCon (negative control), siCD45 (in vivo), Cy5-siRNA (biodistribution).

## 2. What was acquired, and what was not

| | state |
|---|---|
| `source/adma201906128-sup-0001-SuppMat.pdf` | present, 3,171,894 bytes, 48 pages, produced by `itext-paulo-155` via `pdftk 2.02`, creationDate `D:20200130022956+05'30'` |
| main article PDF | **ABSENT, and not acquired** |

**This is the single most important fact about this folder.** The combinatorial
library *screen* — which lipid delivers RNA to leukocytes, and how well — is in
main-text Figures 1–4. The SI carries only ancillary panels.

Tier search, all performed on 2026-11-10 (`00_value_sources.json` has the
detail):

- **Tier 1, publisher Source Data** — absent. Wiley has exactly one
  supplementary deposit for this DOI.
- **Tier 2, SI tables** — present but off-endpoint (see §4).
- **Tier 3, author repo / PMC** — absent, checked two ways.
  `provenance/tools/find_repo.py` over GitHub code and repo search on the DOI,
  the PMID, `Peer lipid nanoparticle`, `Peer LNP` and `combinatorial library
  lipid` returned no relevant hit (the one repo hit was
  `Ed13-coder/monero-lnp2pbot`, a Telegram bot). NCBI `elink dbfrom=pubmed
  db=pmc id=31999380` returns **only** `pubmed_pmc_refs` — articles citing this
  one — and **no `pubmed_pmc` self-link, so there is no PMC deposit** and
  `pmc_fetch.py` has no PMCID to take.

## 3. Vector or raster — measured on every page

`page.get_drawings()` and `page.get_images(full=True)`, all 48 pages:

| pages | drawings | images | verdict |
|---|---|---|---|
| 37–46 (**all figures**) | **0 items** | 1–2 per page, **DCTDecode** (baseline JPEG) | **FLATTENED RASTER.** Tier 4b is unavailable, and this is the checkable reason — there is not one vector item on any figure page, so there are no filled-rectangle coordinates to read a bar height from. |
| 47–48 (Tables S1, S2) | **574 and 691 items, all `re`** | 0 | Native vector — but it is **table furniture**. The numbers are in the text layer, so the correct tier is **4 (text)**, not 4b. |
| 3–36 (schemes) | 24–402 items, `l`/`c`/`re` | spectra as JPEG | Native-vector *structure* drawings, not value drawings. |

So this paper is the **opposite** of the sibling that read 54 values straight
out of filled-rectangle coordinates: here every figure is one flattened JPEG
and tier 4b fails for a reason a reviewer can re-measure in one command.

**The briefing triage was wrong on one point.** It reported "0 table-like
pages"; `pdfplumber.find_tables()` returns **two** tables, 12×14 on p.47 and
12×17 on p.48, and that detected geometry matches the hand count of Tables S1
and S2 exactly. Its "1,745 numeric tokens" also does not reproduce — the same
text gives 3,126 under `re.findall(r'\d+\.?\d*')`. Neither figure is quotable
without its definition. The conclusion is unchanged, because both tables are
toxicology.

## 4. What the SI actually contains

Full panel-by-panel ledger in `00_value_sources.json`. In short:

- **Delivery endpoint: Figure S9 only** — in vivo CD45 silencing, 2 lipids, as
  *targeted* LNPs. Not curated; `DECISIONS.md` **D-6** says why (n = 2 lipids
  makes a per-experiment z-score degenerate; the glyphs are box-and-whisker,
  which `bars.bar_top` refuses by design).
- **Viability**: Figures S2, S3, S4 (U266, XTT).
- **Biodistribution**: Figures S6, S7 (Cy5 signal/area).
- **Physicochemistry**: Figure S8 (TNS apparent pKa).
- **Toxicology**: Tables S1, S2; Figure S11 (cytokines).
- **No numeric scale**: Figures S1 (TEM), S5 (dot plots), S10 (gel), S12
  (gating).

`data/`'s `Experiment_method` vocabulary is delivery-only —
`luminescence_normalized`, `luminescence_discretized_normalized`,
`LRP6_knockdown_normalized`, `protein_adundance_normalized`,
`editing_efficiency_normalized`, `luminescence_relative_to_Spikevax`. Nothing
above is any of those.

**Corpus-admissible delivery rows from this SI: 0.**

## 5. Structures — the part that is durable

14 ionizable lipids, built in `02_enumerate.py` from the SI's own IUPAC names
and written synthesis (pp. 3–36), through `curation/lib/assemble` only:
`attach` for every C–N and C–O bond, `esterify` for every ester, so the acid
loses its OH and the alcohol keeps its oxygen. No `RWMol`, no `RemoveAtom`, no
string surgery on a canonical SMILES anywhere. Every build runs
`check_stereo_preserved` and `check_no_radicals`.

`remap()` in that file is caller-side glue, not an assembler: it forms no bond
and deletes no atom, it only marks which atom the next lib call bonds to, and
it **raises** when that atom is not unique. Assembly order is not synthetic
order for Lipids 3 and 5, and the module's docstring says so and why.

## 6. Anchor posture — RECORDED BY HAND

`provenance/structure_check.py`, `check_anchor_independence.py` and
`check_anchor_coverage.py` were **not run**, per this curation's instructions.
What follows is the posture, stated here so the gap is visible rather than
absent.

`reference_structures.csv` header is exactly
`lipid_name,reference_SMILES,source_figure,read_on,read_by,expected_mismatch,note`.
14 rows: **4 anchors, 10 declared gaps** with an empty `reference_SMILES` and a
note saying why.

| anchor | SI locus | result vs `03_enumerated.csv` |
|---|---|---|
| Lipid 3 | p.7, drawn + full IUPAC name | `FHLGXYGRJTUYSF-UHFFFAOYSA-N` — **exact match** |
| Lipid 4 | p.12, drawn + full IUPAC name | `PNSOECQTDFXNJH-MPKYGIEOSA-N` — **exact match, stereo layer included** |
| Lipid 7 | p.23, drawn + full IUPAC name | `RVAWZYFGAJYABG-UHFFFAOYSA-N` — **exact match** |
| Lipid 14 | p.36, drawn + full IUPAC name | `NGUVKPMCPAGUDU-UHFFFAOYSA-N` — **exact match** |

**4 of 4 exact. 0 mismatches.** `make_reference.py` exits non-zero on any
mismatch; it exited **0**, read from the process.

**Independence is claimed by PROVENANCE, not by byte difference.** CLAUDE.md §6
is explicit that byte-identity is *not* the signal and that testing for it is a
trap, because a correct re-derivation is byte-identical to a correct
enumeration once canonicalised — which is exactly why C10's own
`anchor_may_be_a_copy` missed F-037. The grounds are:

1. **`make_reference.py` imports no assembler.** `grep -nE "import
   assemble|RWMol|AddBond|RemoveAtom|attach\(|esterify\(|CombineMols"` returns
   only two hits, both inside the docstring sentence that declares their
   absence. The only RDKit calls are `MolFromSmiles`, `MolToSmiles` and
   `MolToInchiKey` — parse and canonicalise, never construct.
2. **Every anchor is a source string literal**, typed from the SI's own full
   IUPAC name and the structure drawn beside it. No anchor is derived from
   `03_enumerated.csv`; the file writes `reference_structures.csv` *before* it
   reads that column.
3. **The traversal is reversed.** `02_enumerate.py` builds tail → scaffold →
   head. Every literal in `make_reference.py` is written **head first**, from
   the dimethylamino or hydrazine nitrogen outward. A copy would have inherited
   the enumerator's traversal.

What the anchors **cannot** test: double-bond *position* (D-2), because the
anchor and the enumeration take it from the same IUPAC name. And Lipid 9's
**ether linkage is the one connectivity in this paper with no anchor at all**
(D-3) — the `reference_structures.csv` row says so.

## 7. Checks, with numbers

`07_qc.py` → `07_qc.json`. Exit code **0**, read from the process
(`> file; echo $?`), not through a pipe.

**Every check is seeded and credited only when the agreement count MOVES.**
F-038 is why: three `ZT_2026` marginal checks scored a literal `0` for twelve
days while the summary beside them said "ALL FOUR EXACT". `seed_moved` here is
a comparison of two integers, not a different message.

| check | clean | seeded | moved |
|---|---|---|---|
| P1 ESI-MS mass, three ways | **14/14** | 13/14 | yes |
| P2 doubly-charged ion | **3/4** | 2/4 | yes |
| P3 ¹H NMR total integral | **7/14** | 6/14 | yes |
| P4 olefinic integral vs C=C count | **7/8** | 6/8 | yes |
| P5 structure sanity | **14/14** | 13/14 | yes |

**P1, split three ways — and the answer to "is their Calc. monoisotopic or
average?" is that THE SI PRINTS NO `Calc.` COLUMN ANYWHERE.** It reports only
*found* m/z from a low-resolution ESI instrument, so which theoretical mass each
value means had to be determined, not read. Determined:

- closest to **monoisotopic [M+H]⁺**: Lipids 1, 2, 4, 13 (dev −0.030 to +0.308 Da)
- closest to **average [M+H]⁺**: Lipids 6, 7, 8, 9, 12, 14 (dev −0.006 to −0.212 Da)
- printed as **bare integers**, so unable to discriminate: Lipids 3, 5, 10, 11

All 14 formulae agree with their printed m/z to **≤0.76 Da** on the best of the
three. The convention is **mixed within one SI** (contradiction C-12), so this
check confirms **molecular formula at integer resolution and nothing finer**.

**Mass cannot see connectivity.** P1 and P2 confirm formulae and therefore arm
counts in aggregate; they are blind to which nitrogen an arm sits on and to
double-bond position. P4 counts sp² C–H and so constrains the *number* of C=C,
not where they are. Regiochemistry here rests on the SI's written synthesis and
IUPAC names — a reading, not a measurement (D-2, D-3).

**Figure checks** (`06_digitize_S2_S4.py`):

- Axis: scale **and** origin fitted from detected ticks only. S2 tick-fit
  residual **0.254 px**; S4 **0.33 px**. The brief's warning is right that the
  max-residual test is weak against one bad tick, so `axis.loo_residuals` was
  called: **worst LOO residual 0.355 px (S2) and 0.540 px (S4)**, against a
  1.5 px tolerance. Both **PASS**.
- Bars: S2 **16/16** detected, 0 refused, max flank delta **0.0 px**.
  S4 **12/12** at 0.06 µM, 0 refused, max flank delta **2.0 px**.
- Colour: each S4 bar assigned by `heatmap.match` against **S4's own** legend
  swatches sampled from the same raster. **Max nearest-swatch distance 4.123**
  (above ~5 would mean the geometry is wrong), **0 ambiguous**, and the
  assignment **agrees with the independent per-colour mask** on all 12. Each of
  the 12 colours yields exactly **5** column runs, one per dose group.
- **The unprinted round-number anchor**: Figure S2 normalises to untreated
  cells, so the UT bar must sit at exactly 100%. The authors never *print* 100
  — they drew it. Measured **100.11%**, error **+0.11 pp**.

## 8. Corpus overlap — zero, on three keys

`04_corpus_overlap.py` against 30,998 corpus rows (29,195 vitro + 1,803 vivo),
19,919 unique full InChIKeys, 0 unparseable `IL_SMILES`:

- full InChIKey: **0 / 14**
- skeleton block: **0 / 14** (so no F-013-shaped same-skeleton/different-stereo hit)
- `IL_name` (`Lipid 1` … `Lipid 14`): **0 / 14**

The paper's DOI and PMID appear in **0** corpus rows.

**This costs something.** Where overlap exists it is a free external check — a
sibling got 21/21 independent confirmation of its head/tail assignment because
corpus lipidoid names encode chain length. There is no such gift here, so the
ESI-MS anchors, the NMR integrals and the four C10 anchors carry the whole
external burden.

## 9. Distinct structures / formulations / individually measured

Three different numbers, and they differ:

| | count |
|---|---|
| **distinct ionizable-lipid structures** | **14** (Lipids 1–14, all enumerated, all unique by full InChIKey) |
| **distinct LNP formulations in the SI** | **not recoverable.** The SI gives no helper-lipid identity, no cholesterol or PEG molar ratios, no N:P and no mixing method for any formulation. It names a *targeted* variant (tLNP, Figure S9) and a plain LNP, but specifies neither composition. |
| **individually measured, on a delivery endpoint** | **2** (lipids 2 and 6, Figure S9, as tLNPs) |
| individually measured, any endpoint | **14** — but only in Figure S2, a viability panel. Figures S3 and S4 cover **12** (Lipid 8 and Lipid 10 are absent from both legends, C-02); Tables S1/S2 cover **5**; Figures S6/S7 cover **4**; Figure S8 covers **4**. |

So the library is 14, the SI measures delivery for 2, and no formulation is
fully specified anywhere in it.

## 10. Contradictions

`08_contradictions.csv`, **12 rows**, 5 resolved and 7 open, each with the
measurement that settled or failed to settle it. The three that matter:

- **C-01 / D-4.** Figures S3 and S4 give different colours to Lipids 12 and 13
  *and* list them in different legend orders. S3's "Lipid 13" swatch is
  RGB (127, 4, 127); S4's "Lipid 12" swatch is RGB (128, 0, 129) — **Euclidean
  distance 4.58**, *below* `heatmap.MAX_DIST`. A cross-panel colour match would
  return a confident wrong answer and never trip the guard that exists to catch
  it. Corroborated by the drawn bar order, which follows each panel's own
  legend.
- **C-04 / D-5.** Figures S2 and S4@0.06 µM are the same condition. **8 of 12**
  shared lipids agree to within 0.3% (ratios 0.9971–1.0034) — which validates
  the digitisation *and* shows the panels share one dataset — while **4 do
  not**: Lipid 2 ×0.943, **Lipid 3 ×1.313** (93.75% vs 123.10%, a 29.35 pp
  gap), Lipid 4 ×1.044, Lipid 5 ×1.110. Spearman across the 12 is **0.48**.
  Confirmed by eye: Lipid 3 really is the tallest bar in S4's 0.06 group.
- **C-05.** Lipid 4's printed olefinic integral says four C=C; both its
  precursors report 4 vinyl H (two C=C) and its own mass settles it —
  C36H67N3O5 gives [M+H]⁺ 622.515 against the found 622.6, while four C=C would
  give 618.48, 4.1 Da away. **Structure stands; the printed integral is
  wrong.** Retained in the file with the measurement that retired it.

## 11. Files

| file | what |
|---|---|
| `00_value_sources.json` | the intake gate — written **first**, every panel, its tier, and why tiers 1–4 fail where they do |
| `01_blocks.csv` | 16 building blocks, each a source literal with its SI page |
| `02_enumerate.py` → `03_enumerated.csv`, `03_assembly_trail.json` | the 14 structures and every intermediate |
| `04_corpus_overlap.py` → `.csv`, `.json` | three-key overlap against `data/` |
| `05_values.py` → `05_values.csv` | **347** values with full metadata; `corpus_admissible` = no on every row |
| `06_digitize_S2_S4.py` → `06_figure_values.json` | the cross-panel test |
| `07_qc.py` → `07_qc.json` | every check, seeded |
| `08_contradictions.csv` | 12 rows |
| `DECISIONS.md` | 8 entries, all `OPEN` |
| `make_reference.py` → `reference_structures.csv` | 4 anchors, 10 declared gaps |

Reproduce in order:

```bash
cd curation/peerd_31999380
python3 -I 02_enumerate.py
python3 -I 04_corpus_overlap.py
python3 -I 06_digitize_S2_S4.py
python3 -I 05_values.py
python3 -I 07_qc.py            # echo $? from the process
python3 -I make_reference.py   # echo $? from the process
```
