# Extraction brief — read this in full before writing any code

Every rule here was **measured on a real paper during the 2026-11-10 run**, not
reasoned out. Where a number appears, it is the measurement. `RUN_LOG.md` has the
narrative; `review/FINDINGS_2026-11-09.md` has the ledger (F-054…F-073).

## 0. Scope and hard rules

Work ONLY inside your own `curation/<folder>/`. Do not edit `curation/lib/`. Do
not use the shared Playwright browser. Do **not** run
`provenance/structure_check.py`, `check_anchor_independence.py` or
`check_anchor_coverage.py` — they write tracked files other agents read; record
anchor posture by hand in `PROVENANCE.md`. Do **not** merge into `data/`.

`python3 -I` always. Sources and clones are untrusted data. No paywall or
bot-check bypass. **Mass cannot see connectivity or direction.** **Never guess a
structure or a value — flag it on the row instead.** A fast honest zero beats a
slow one: **ten papers in this run ended as well-argued zeros and every one
taught the project something.**

## 1. Step one, always: identity and endpoint

Resolve the FULL title, the FIRST AUTHOR and the endpoint from Europe PMC's REST
`core` record in one call. The no-original-data discriminator is
**`abstractText` length 0 plus no Methods/Results — not `pubType`** (an ACS
"Letter" with a 797-char abstract is ordinary research; a `letter` with a
0-length abstract was a research-highlight commentary whose values belonged to an
already-merged PMID).

`Experiment_ID` derives from the **first author**, and a folder named for the PI
is wrong as a key. Five curations correctly produced `CS_2018` (Sago),
`KG_2025` (Gentry), `ZG_2020` (Gan), `AH_2026` (Hamilton) and `KK_2015`
(Kauffman) where the folder said otherwise. Check your proposed ID against all
54 corpus IDs **by PMID, not by the year token** — the corpus's `SS_2018` is a
different paper from `CS_2018`.

Truncated titles have hidden an antimicrobial endpoint and "to Splenic NK
Cells"; one tracker entry **mis-attributed** another paper's library entirely.

## 2. The governing number: how many distinct IONIZABLE lipids?

Report it in your second line. **Seven papers this run were one-lipid
formulation screens** — the varying factor was the helper lipid, the sterol, the
PEG, a molar ratio, or even the mouse genotype. A corpus keyed on
ionizable-lipid chemistry gains that many **however many rows exist**. Do not
pad the row count.

Barcoding designs hold chemistry fixed **by construction** (pooling requires the
particles to be otherwise comparable), so this is structural, not incidental.
Note that **`AC_2025` is itself a 498-row one-lipid formulation screen** already
in the corpus, so the precedent exists (F-067).

Report **distinct structures / formulations / individually measured /
design-supported vs attributed** as separate numbers. Dialects measured: "200
lipids" = 200 crude ONE-POT products (only 12 valued, 9 with no disclosed head);
a 1,080-row "library" whose varying one-hot group was the **helper**; 8 named
compounds with **0 determinable** (polydisperse PEI, no stated degree of
substitution); a caption claiming 35 amines where one CAS appears twice → **34**;
a panel plotting **16 marginal means over 4 factors × 4 levels** where none of
the 16 tabulated formulations was individually reported; designs supporting 146,
309, 1,744 and ~2,700 per-LNP values that attributed **5, 0, 0 and 0**.
"Not determinable" and "not recoverable" are correct answers. **If a factorial is
claimed, prove it** (4×2×8=64; 10×3×4=120; 18×2×4=144; 5×6=30=row count;
4×2×4×3×3×3=864).

## 3. Endpoint admissibility

`data/`'s `Experiment_method` is **delivery-only**. Admissible: luminescence,
protein abundance (EPO/hFIX ELISA, %GFP+, Thy1.1+), editing efficiency, TIDE
indel %, protein-silencing MFI, RT-qPCR knockdown, and **Cre→tdTomato
recombination** (downstream of translation).

**NOT delivery:** viability, particle biodistribution, pKa, size, PDI,
encapsulation efficiency, **cellular uptake**, cytokine panels, body weight,
tumour volume, scRNA-seq of the host response, and **cargo abundance measured
before translation** (F-071 — one paper's authors said verbatim they sampled
early because "peak mRNA translation occurs relatively late").

**Pooled/barcoded is NOT itself disqualifying** — `LR_2023` contributes **407 of
1,803 in-vivo rows** as `barcoded`. But the precedent is narrower: it carries a
**translated** readout. The line for Methods: *barcode recovered from cells
sorted on a translated reporter is delivery; barcode recovered from tissue is
biodistribution.*

If pooled, test whether the readout is **compositional**: reconstruct the
normalisation from the paper's own worked example and report the residual.
Achieved: `nd_i = [(c_i/Σc)/(in_i/Σin)]/Σ_k[same] × 100%` at **residual 0.0,
21/21 and 18/18 cells exact** — and that curation **said the 0.0 was against
integers the authors rounded**, giving pre-rounding discrepancies of 0.299/0.291
pp against 0.25 expected, and a **±0.5 pp sensitivity limit**. Two others
reported the residual **not computable** because the paper publishes two
different formulas and no worked example — say that rather than invent one.

Put `corpus_admissible` yes/no **with a reason on every row** of
`05_values.csv`, counts summing. `06_rows.csv` header-only if zero, so the zero
is visible.

## 4. Acquisition

Write `00_value_sources.json` FIRST, each tier a checkable sentence. Tier 1
publisher Source Data / SI workbook · 2 SI numeric table · 3 author repo · 4 PDF
text/vector · **4b `4b_pdf_vector_geometry`** · 5 raster pixels.

**Tier 1 and 3 are worth real effort:** the curations that took values from a
repo or Source Data file have **zero digitisation defects between them, against
eight for the four that digitised.**

- `find_repo.py --pmid … --title "<short words>" --author <surname>` — run for
  **both** first author and PI; a wrong author term returns a clean-looking zero.
- **Grep the full JATS AND the SI** for `GSE\d+`/SRP/PRJNA/SRX/SRR/SRA/
  accession/deposit/github/zenodo/dryad/figshare/osf.io and **report a NULL
  RESULT EXPLICITLY with the character counts searched** (one reported 0 hits
  across 1,895,389 chars). A GEO/SRA accession is **not** a delivery value.
- `pmc_fetch.py <PMCID> <outdir>` — **first check `elink` for a real
  `pubmed_pmc` SELF-link**; `pubmed_pmc_refs` alone means no deposit exists.
  PMC's `/bin/` route serves **figures** but returns ~21.5 kB **HTML error pages
  for supplements**.
- Europe PMC `/<PMCID>/supplementaryFiles` is **per-PMCID**: of three probes, one
  gave a valid 5.4 MB ZIP, one HTTP 500, one a 297-byte XML stub. Judge by
  **magic bytes AND size** — `504b0304` for a real archive.
- **Elsevier CDN, anonymous and reliable**: DOI → `api.crossref.org/works/<doi>`
  → `alternative-id`'s `S`-prefixed entry is the PII →
  `https://ars.els-cdn.com/content/image/1-s2.0-<PII>-mmc<N>.<ext>`. This
  recovered three SIs including a tier-1 Source Data workbook.
- **Elsevier high-res figures**: `…-gr<N>_lrg.jpg` gave one curation 2371×2857
  against 790×952 — **9× the pixels**.
- **Validate every download by MAGIC BYTES.** A `.jpg` fetch returned a ZIP; two
  HTML error pages were correctly rejected. A size floor passes all of these.

**F-065 — the session scratchpad is SHARED, and a directory NAME is not a safe
container.** A glob copy put **17 figure files from other papers** into one
folder's `source/figs/`; another agent's clone landed **inside a pre-existing
`repo/` directory holding a different paper's project, with its own `.git`**. A
file-count or size check passes both. A digitiser pointed at a foreign figure
yields values that are internally consistent, geometrically valid and
**attributed to the wrong paper** — the `AC_2025` class. So: never
`cp <scratchpad>/*`; use a **uniquely named, verified-empty** directory; and
attest the result against something only the right source can produce (one
sibling matched its SI's **sha256 to the file inside Europe PMC's own ZIP**;
another matched an ACS watermark path to the manuscript id inside its own DOI).

## 5. Corpus overlap — join on InChIKey, and every control must assert its premise

**F-070 and its amendment.** One sweep returned a clean 0 because its SMARTS
demanded a ring CH₂ and could not match a molecule known to be present. Another
then caught **its own control**: a hand-written sterol SMARTS returned **6**
structures and reported `passes: True` **while not matching the registry
cholesterol at all**.

> Every corpus query carries a positive control, **and the control must assert
> that it matches a specific named molecule you know is present.** A non-zero
> count is not evidence the control works.

Working controls, with numbers: `DSPC` name scan **7,038** rows; InChIKey
`DLin-MC3-DMA` **3,098** across **two spellings**; cholesterol **30,397**;
SMARTS tertiary amine **18,472**; cKK-E12 name **216** → InChIKey **241**.

**A structure can sit under several names.** cKK-E12 is `CKK-E12` (`AC_2025`,
216), **`other`** (`JW_2024`, 19) and **`MD-1`** (`BL_2023`, 6 vivo) — one
InChIKey. **324 of 19,919 structures carry more than one `IL_name`; 292 names
carry more than one structure.** C11 pins name→structure; **nothing checks the
reverse** (F-063).

**The in-vivo pool has exactly TWO structures shared across experiments** —
`DLin-MC3-DMA` (6 sources) and `SM-102` (4). Any in-vivo lipid already present
elsewhere is **rare and load-bearing**; report it prominently. Organ caveat:
SM-102's four in-vivo sources are *ear* and *whole_body* only, so a new liver
source raises the count while adding **zero organ-matched pairs**.

Run a **substructure sweep labelled a SUPERSET** (measured: 19 loose → 0 strict;
23 → 2 → 0 exact; 1,535 → 1; 1,247/18,720 and 1,584 both correctly called
chemotypes).

**Token-pattern test, mandatory.** The amine catalogue IS shared across labs
(`93` = 1-(3-aminopropyl)imidazole in `AA_2008`/`KW_2014`; `9322` = its 2-methyl
analogue in `YR_2024`; `80` in three experiments) and has pinned structures
papers never drew — one got **43/43 by retro-Michaelisation**. **But it is
sparse and non-contiguous (255 values; only 7 of 1–16 present), so contiguous
1..N numbering means LOCAL INDICES.** Measured damage from naive substring
joins: `80` → 94 rows/85 structures; `7C1` → 36/35 with **none** the intended
lipid; `200` → 526 rows/24 structures/≥3 chemotypes; bare `9` → **5,590
rows/4,817 structures**; bare `1` → 20,837. And the **wrong-namespace** case:
`C14`/`C18` naming a *PEG tail* in one paper and *lipidoid tails* in the corpus;
`A-9…A-17` being *helper* lipids. Also a false **zero**: **`MC-3`, a paper's own
spelling, matches 0 corpus rows for a lipid present 3,098 times.**

Prefer a **structural attestation**. Use `curation/corpus_components.py`'s
guarded `load`/`smi`, which **RAISES on an ambiguous name** — but it cannot catch
a name unambiguous in the corpus and simply **wrong for your paper**
(`C14-PEG2000`, drawn as a phosphoethanolamine-carbamate-PEG, stored as the
glycerol-ether). **Never conclude from an exploratory grep**: one agent called 22
rows "an unrelated index with no imidazole" when all 22 carried the head and
belonged to a lab it had misnamed, because its regex required the name to *start*
with `93`.

**Report whether benchmarks top the screen.** Measured: three papers no (four
bars within 1.945× all "ns"; a novel lipid at **1.411× SM-102**; LNP7 beating MC3
**2.39×**/**4.06×**), one yes outright (SM-102 first in all six organs, the lead
winning only a **13.6×** selectivity ratio), one found its own marked "lead"
ranks **2 of 22**, two use no benchmark at all (moot — say so), and one had **4
of 4 shared lipids be universal controls**, which is §8's hazard in its purest
form. A "winner" selected into a second screen tops it **tautologically**.

**Check whether an arm of your screen is already merged.** One found its column
joined an existing experiment at Spearman **1.000000**, then showed the ln-unit
residual rms was **2.8531e-4 against 2.8868e-4 expected from three-decimal
rounding** — the discrepancy *was* the rounding — and excluded that arm.
Duplicating a merged experiment under a new `Experiment_ID` would manufacture
exactly the source-separability signal this project measures.

## 6. Reading values

**A column header is a claim, not a measurement.** One workbook column labelled
"Fluorescence intensity" **held relative expression** (values 0.06–0.83 against
2.5e4–2.3e7 in the real fluorescence columns, and one of them was the abstract's
own printed figure). The corpus's own `van.der.Waals.Molecular.Volume` never held
a van der Waals volume. **Check a column's values against its name.**

**Read the primary object, not a proxy.** Counting `<row>` tags in an xlsx gave
51 and 49 where reading the cells gave **25 data rows**. Three tag-counting
methods, three wrong answers.

**`find_screens.py` silence is not evidence** — it printed nothing for a 72-row
validation workbook and for a 12 kB tier-1 Source Data file, both found by
dumping the sheet XML directly. Its sibling `triage.py`'s `table_pages` is a
**TEXT heuristic wrong in BOTH directions**: read `native_tables` and
`ruled_tables`, verify by hand. `find_tables()` returns 0 on pages that are
nothing but table when the rules are **filled rects 0.36–0.48 pt tall** (re-read
with `horizontal_strategy="explicit"` and the rect runs as
`explicit_horizontal_lines`, requiring real filled **and numeric** cells — a
recovered grid of empty cells is not a table). An SI value table can be a pure
**BITMAP** with ~300 chars of page text, or **three full data tables embedded as
rasters that look like plots**. One SI's contents page named none of its four
value tables; another's TOC was **cyclically permuted by one** against its own
captions. Cite by caption, match by content, say which you used.

**Mass conventions drift, sometimes within one column** — `[M+Na]+` ×9 and
`[M+H]+` ×1 in one column; "calc'd" average on two compounds and monoisotopic on
six; a neutral "calc." against [M+H]⁺ giving a **1.0421 Da residual, exactly one
proton** (under a 2 Da tolerance every row "passes" a proton out); and **no
`Calc.` column at all**, where the convention was determined from the data with a
stated margin over next-best. Same hazard for **pKa**: one SI's two tables were
LNP *surface* pKa by TNS and *lipidoid* pKa by titration, differing by up to
**1.14 pH** and disagreeing on direction.

**A name token is not a chain length** — `C24` = two C4 tails; `Cm-n` = a
carbinol position; `O16B`'s "16" counted chain **heavy atoms** (2 C + 2 S + 12
C); `17D9.2` summed to **18.993** carbons; `A-N` carried **N+1** acyl carbons
(10/10); one key read `C12 = C12H25` while the arm is CH₂–CH(OH)–R so "C12"
carries **14**; `7C1` packs three numbers, none a chain length. **Derivable
method:** within a series holding the amine fixed and walking the tail,
**dC/d(tail) equals the arm count** — measured 5.0,5.0,5.0,5.0 for a 5-arm series
and 1.0,1.0,1.0 for a 1-arm one.

**If two sources share a member, measure the ratio.** Measured: 8 of 12 agreeing
within **0.3%** (which validated the digitisation *and* proved shared data,
making the other 4 non-noise at up to ×1.313); a within-paper same-dose **rank
inversion at ρ = −0.20**; a same-lab cross-cell-type **−0.089**; a five-fold
same-lab replication at **0.33–0.49**; round repeats at 1.053×–5.011×;
**14/14 inside the PBS band the figure itself draws**, joined **on label not
index**; and two panels of one condition disagreeing in **opposite directions**
(×1.129, ×0.703), moving the paper's own effect size from 1.1818 to 1.8984.
**Never average or rescale two disagreeing sources** — emit separately or not at
all. If a shared value is blocked by a refusal, report the ratio **not
computable** rather than using the better-looking flank.

**Do not colour-match ACROSS panels** — two panels gave different colours *and*
different legend orders to the same two lipids at distance **4.583 < `MAX_DIST`
5.0**. Call `heatmap.cross_panel_collisions(a, b)`; match within a panel only. A
heat map with no numeric scale is not extractable; one **with** a scale is still
unusable if its columns are unlabelled.

## 7. Digitisation, only if you must

**Check `page.get_drawings()` first.** One curation read a whole paper at tier
**4b with no page rendered**: 239 path groups separating bar body / outline /
error bar / axis by **exact stroke width** (0.4596 / 0.2757 / 1.1949 pt), so
error class 32 **cannot arise** — nothing is thresholded, so no threshold can be
wrong — and it set its own tolerance at **0.01 pt** because a vector read has no
detector noise. Another found exactly **2 vector paths per page, both the
publisher watermark**, which is the checkable reason tier 4b fails. Native vector
may be **table furniture** or **scheme drawings** rather than data: if the
numbers are in the text layer, that is tier **4**, not 4b.

**`chains.vector_bonds` returns zero where bonds are FILLED POLYGONS with
`width=None`** (measured three times; 603 of 915 paths on one page) — pass
`include_filled=True`, separate glyph outlines from bonds (by stroke thickness,
aspect ratio or label proximity), and prefer the closed identity
`drawn == heavy_bonds + kekulised_double_bonds` (closed 3/3 and 48/48).

**If raster, MEASURE NATIVE DPI** (`px / rect_pt × 72`; measured **137–355.7**).
Use `axis.native_tol(axis.LOO_TOL_PX, native_dpi=<measured>)`, which **refuses a
guessed dpi**. A standalone JPEG with no page rect has no derivable dpi — say so.
Two limits measured by the orchestrator: **`calibrate`'s hard 2.5 px threshold is
NOT rescaled by that helper** and raises first; and the rescaled tolerance is
**weak at the margin** (at 204 dpi a +3.25 px tick error passes by 0.08 px).
**Seed the widening in your own geometry**, report both figures and the dpi,
never widen silently — **and if the strict 1.5 px gate already passes, say the
widening was not relied on** (six curations did exactly that). State precision
honestly: 18.42 px/decade means 1 px = factor **1.133**, so a 2.5 px tolerance
would permit a factor-**1.37** error and pass.

**Respect every `bars.bar_top` refusal and fix it at the source**, never by
loosening a threshold. Diagnosed: replicate rings straddling a top corner; a
2-row right-flank rule; a saturation floor of 28 deleting every teal bar
(channel spread **4**); 18 px strokes meaning a fill can never reach 90% width; a
2–3 px JPEG halo causing 9 spurious refusals; **legend swatches inside the plot
box**; a grey mask swallowing the axis rule; a chroma `max_val <= 252` ceiling
**deleting a whole bar** (254,224,224); an **inset table drawn INSIDE the plot
box** whose ~300 px rule read **149.4% on a 0–150 axis** (fixed with a furniture
strip — *a horizontal run longer than three bar widths is not a bar top*); an
**error-bar stem down the bar centre** biasing **36 of 43 "clean" reads LOW by up
to 12.5×**, every one moving upward when fixed, the signature of occlusion not
nudging; **replicate dots spanning a full bar width** making one bar read
**97.73 instead of 77.84** until a solidity walk up from the baseline was added;
and **parallelogram bars** (~2 px shear) where a de-skew was **tried and
rejected** because neither the spine slope nor the bar-edge slope validated
itself. **Replicate dots sit on BOTH top corners** in some panels — the one shape
the flank test cannot see. Prefer reading bars by **colour identity** where the
figure allows it.

## 8. Structures

Use `curation/lib/` — `axis`, `bars`, `heatmap`, `chains`, `assemble`,
`reagents`, `emit`, `errors`. Import it; never write your own geometry,
assembly or emit code.

**Build fragments with atom-mapped attachment points, through `assemble`.** **A
SMILES fragment has a DIRECTION that mass cannot see** — this bit five agents,
every version carrying the correct formula and exact mass, parsed and
canonicalised, **caught only by InChIKey**: an anchor scoring 0/8 then 1/8 then
5/8 from three instances of one class; a chain's C1 end on the wrong carbon; DSPC
built as **a ketone bearing an ether** at the right formula `C44H88NO8P` and mass
`789.624755`; a Markovnikov regiochemistry at **Δ 0.000000 Da**; and a tail built
by string concatenation that became **a ketone bearing a free alcohol instead of
an acrylate on all 11 rows — two lines below its own comment warning against
string surgery**, caught only by printing the column.

**Do not hand-write a sterol.** One agent's cholesterol came out a
**diastereomer with 4 of 8 centres inverted** — identical formula, identical
exact mass to six decimals, **identical skeleton InChIKey block**, i.e. F-013
reproduced from scratch — and it correctly **discarded the anchor rather than
flipping tokens until it matched the registry**. Another **declared its registry
cholesterol pin NOT independent** rather than counting it as an anchor. Do both.

**Never write a benchmark SMILES from memory** — one agent's SM-102, MC3 and
ALC-0315 were **all three wrong** (one CH₂ short, one CH₂ short on a tail, one
ester reversed), all parsed and canonicalised; only the InChIKey saw it. Use
`provenance/reagent_registry.py`, and note **not every excipient is in it** (no
18:1 PA entry; the only PA is `14PA` = DMPA). A drawn **PEG bracket can disagree
with its own subscript** (44 printed, 45 drawn, one outside the bracket; the
subscript gives a 44 Da error only the InChIKey catches).

Respect `assemble`'s refusals and reorder the route (one fixed a radical refusal
by reordering to attach→esterify→attach→esterify). **If achiral, say the stereo
assertions are VACUOUS, not passing**; assert double-bond geometry explicitly,
since CIP counting is blind to it. `assemble.parity_preserved(product, fragment)`
compares **spatial parity** rather than CIP counts, for a transformation that
merely relabels a centre.

## 9. The C10 anchor

Byte-identity is **NOT** the signal and testing for it is a trap (F-037) — a
correct re-derivation is byte-identical to a correct enumeration once
canonicalised.

Verify independence by **walking your module's own AST at run time, inspecting
`Import`/`ImportFrom`/`Call`/`Attribute` nodes only** so the ban list cannot
match itself (one agent's guard failed on its own module for exactly that
reason; another's substring scan flagged its own **docstring**). Add a
**template-provenance guard**: no literal of ≥25 characters may be a substring of
any enumerated SMILES. **Make both guards fail on purpose** in a child process,
with a clean control that must come back clean. Prefer a traversal **opposite**
to your enumerator's, and pick anchors of **different chain lengths** — one
length cannot catch an off-by-one that scales with n.

**Prove the gate tests connectivity** by seeding a defect with **identical
formula and identical exact mass**: achieved with a reversed ester (Δ 0.000000),
an MC3 regioisomer (641.611081 both), a Markovnikov opening (299.318815 both), a
1,3-diacyl-2-phospho regioisomer, and a cholesterol epimer where formula, exact
mass **and the skeleton InChIKey block** all stayed identical.

**If the paper draws nothing, store empty `reference_SMILES` rows so the absence
is visible** — and look for an intermediate it DOES draw in full, especially one
appearing more than once per member (one anchored a tocopheryl-acrylate that
appears twice in every library member). **If the ionizable lipid is undrawn,
anchor the excipients instead** — three curations got real C10 results that way.
**Against an empty anchor the guards are VACUOUS, not passing**; saying "C10
clean" there is F-013/F-037 exactly.

## 10. Checks and seeds

**A CHECK THAT CANNOT FAIL IS NOT A CHECK, AND A SEED THAT ONLY RAISES IS NOT
CREDITED.** Credit a seed only when the **agreement count MOVES**. Shapes caught:

- a check scoring **11/11 clean and 11/11 seeded** because it measured spread
  about the mean of its own set (shift-invariant);
- a seed exiting non-zero while **moving nothing** because `calibrate` refused
  first — correctly split into a subtle seed that passes calibrate and fails LOO,
  **plus a gross one credited explicitly by the refusal**, recorded as a distinct
  route;
- two seeds moving a metric while **the gate stayed silent** — that agent widened
  **the gate, not the seed**;
- no-op seeds that removed the wrong object, injected something a size filter is
  **meant** to exclude, or rolled the whole image (moving ticks and bars
  together);
- a seed exiting from a `TypeError` **without reaching its metric** — that
  harness now **refuses to credit a seed whose metric line is absent**. Adopt
  that. One agent then hit it twice more and **retargeted the seeds rather than
  weakening the guard**, crediting the guard with a seed of its own;
- a `spine_erased` seed that was a **no-op until its bounds came from the image**
  instead of a constant;
- a **1e-9 round-trip bound that could only ever fail**, replaced by observed
  **2.0164e-5** against **2.0233e-5** predicted from rounding alone, ratio
  0.9966;
- a **tautological check**, refused on those grounds and replaced by a row-wise
  identity against the file's *second* parameterisation (0 violations at 1e-6 on
  1,080 rows). **`mol% sums to 100` is a tautology where cholesterol is the
  stated balance** — say so; where the four percents are printed independently it
  is a real transcription check (measured 248/248, 158/160 with two at **stride
  4** in a twice-pasted block, 418/424 at stride 16, and **46/47 catching one
  real error**).

**Report your measured sensitivity limit** honestly — examples: ±0.5 pp; a +1 px
tick shift moving no value by >0.25% recorded as *harmless, not credited*; a
full-width marker **leaking** at 18 px where the widest real cap was 0.647 of the
bar against a 0.90 gate (margin 0.253, a stated limitation not a live defect);
and a builder that **cannot see double-bond position** (Δ9→Δ6 passes at Δmass
0.000000, caught only by InChIKey).

**`python3 … | tail` returns tail's exit code** — F-012's shape; read verdicts
from `returncode`, never stdout or a pipeline, **one child process per seed**.

**A SEED HARNESS MUST WRITE TO A TEMP DIRECTORY, NEVER THE CURATION FOLDER
(F-069).** One harness pinned each mutated copy's `HERE` to the real folder so it
could read inputs, and the copies wrote their **outputs** there — leaving
`reference_structures.csv` holding a **seeded regioisomer anchor** (Δmass
**0.000000 Da**) written by a run that exited non-zero, **so the non-zero exit
concealed the write**. A file-count or exit-code check passes; only a cross-run
hash comparison sees it. Use a separate `OUTDIR` and assert every tracked
artifact is byte-identical across two clean runs **after** the seed suite runs.

**Look for a round number the authors DREW or IMPLIED but never printed** — the
best checks in this run were all of this kind: `Total (%) = 100` drawn under a
table and never computed (residual 0.0); an **average column equal to the mean of
its own three organ columns to 1.67e-9 on 22/22**, validating all 88 transcribed
numerals because one misread digit breaks it; **pie labels pinning a unique
composition row and recovering an unlabelled slice**; a formulation matching the
published **Onpattro ratio 50 : 38.5 : 10 : 1.5 to |Δ| 0.0 mol%**; **"2×; 160
LNPs tested in vivo" drawn and never computed**, where 120+40 = 160; a printed
yield range **"77–92%" that reproduces only with the disputed compound
included**, settling the paper's own off-by-one; a **plot box exactly 16 bins
wide** (199.6795 vs 199.6799 pt, **0.00048 pt**); a **frame top on an unlabelled
110% level**; **minor ticks predicted from a fit that EXCLUDED them** (passed as
`forbid_target`, making class 27 a type error); **8 minor ticks per decade at
log10(2..9), unequally spaced, pinning the decade length at a 4.15× margin**; a
dashed reference rule landing on **1.000**; an untreated bar that must be exactly
100% (read **100.11%**); the **zero tick on the drawn baseline to 0.000 px**;
**11 x-label right edges at tick + 2.994 pt with 0.002 pt spread**; **Σ printed
NMR integrals = the formula's H count exactly**; consecutive analogues **exactly
one CH₂ apart to ≤0.05 mDa** where the paper never says so; **38/38 tail steps
each exactly one CₙH₂ₙ**; **two relations among three printed HRMS values that
close only if the tail multiplicities are 2,1,2**; the paper's own **12
significance marks sorting 12/12 on a single fold threshold** (largest "ns"
1.945×, smallest starred 2.212× — a t-test is not a function of fold alone); a
panel's **red-vs-black colouring separated by one threshold** revealing an
unstated 50% cut-off; **four printed correlations reproduced to ≤0.016 by
Spearman but not Pearson**, identifying the unstated metric *and* exposing a
printed **sign error**; and an **unprinted dose recovered by matching a
control's mean AND SD across four tables**.

**Prose groupings are a strong unprinted test**: show each qualitative set the
paper names is separated by a single threshold on your matrix, with the margin
(one got 5.16 pp, 45.4 pp and 0.59 pp, and explained an apparent anomaly).

**A LOW-RESOLUTION AXIS CANNOT PIN A VALUE BUT CAN STILL TEST ONE.** One agent
refused a ~4 Da/px mass axis as unusable — correct as a *read* — then used it to
**discriminate** candidates once it had a formula: one CH₂ apart predicted +0.8
vs +14.9 Da, measured **+0.20 px against +3.78 px, an 18.9× margin**. Pinning and
testing are different questions; do not let a refusal for one become a stopping
point for the other.

## 11. Deliverables, in order

1. `00_value_sources.json` (**first**). 2. `01_blocks.csv` / `03_enumerated.csv`
+ `02_corpus_overlap.json`, columns `IL_SMILES` / `IL_name`. 3. `05_values.csv` —
every value with source location, readout/units **as printed**,
cargo/cell/organ/route/dose/timepoint/assay, pool identity and size if pooled,
and `corpus_admissible` + reason. 4. `06_rows.csv` for the admissible ones
(header-only if zero; read the 79-column header **live** from
`data_vitro.csv`/`data_vivo.csv`), flags **on the row** via `emit.emit_rows` with
**per-item callables** — a hard-coded literal `confidence` must raise (F-011), a
caveat with no note must raise (F-024). **One z per experiment**, pooling cell
types and organs (`value_normalization_unit = experiment`); a per-cell,
per-organ or per-cohort z is **error class 22** and once deleted a paper's own
published selectivity, and once would have forced a paper's entire finding to
zero. Carry design differences in `Model_target`, not in separate
standardisation. 5. `07_qc.json` — every check with **numbers**.
6. `08_contradictions.csv`, **keeping retracted suspicions with the measurement
that retired them**. 7. `DECISIONS.md` — `### D-n` then a plain `STATUS: OPEN`
line (**not bold**; bold made the merge gate report zero open decisions for a
folder with five). 8. `make_reference.py` + `reference_structures.csv`, header
exactly `lipid_name,reference_SMILES,source_figure,read_on,read_by,expected_mismatch,note`.

**Record the LICENCE** for the article and any repo, **off the Europe PMC record,
not the publisher landing page**. §8 needs CC0/CC-BY. Measured this run: three
cleanly **CC BY**; two **CC BY-NC-ND** (NC *and* ND each fail); one failing twice
(closed article + **GPL-3.0** repo — copyleft *software*, not a data licence);
three NIH author manuscripts under "text mining … fair use", which is **not** CC;
and one All Rights Reserved.

## 12. Report back

THE ENDPOINT (first line) · **how many distinct IONIZABLE lipids this paper
contributes** (second) · where the values are and the tier · whether the readout
requires expression, and if pooled whether it is compositional with the residual
or an explicit "not computable" and why · distinct / formulations /
individually-measured / design-supported vs attributed · corpus overlap
including your positive controls, any already-merged arm, any in-vivo shared
lipid, and the token-pattern test · rows emittable · every check with numbers ·
licences · DECISIONS · what you could NOT verify.
