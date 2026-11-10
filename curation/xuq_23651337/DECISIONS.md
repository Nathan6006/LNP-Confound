# DECISIONS — `xuq_23651337`

Wang M, Sun S, Alberti KA, Xu Q. *A combinatorial library of unsaturated
lipidoids for efficient intracellular gene delivery.* ACS Synth Biol 2012;1(9).
PMID 23651337, DOI 10.1021/sb300023h.

The pipeline detects that a source is silent; it does not get to decide what
the source meant. Everything below is queued, not resolved.

**This folder emits zero rows**, so none of these decisions currently blocks a
number. They block any *future* attempt to curate this paper once the main
article is in hand — which is the only circumstance in which it becomes
curatable at all.

---

### D-1 — This paper cannot be curated from the material in hand. Acquire the main article, or drop it.

Three independent reasons, each on its own sufficient:

1. **No admissible value exists in the file we have.** The 3-page SI holds
   exactly three panels: TEM micrographs (no number), an ethidium-bromide
   exclusion curve (nucleic-acid binding, not delivery) and a 16-bar cell
   viability chart (cytotoxicity, not delivery). `05_values.csv` records 18
   rows, `corpus_admissible = no` on all 18.
2. **The delivery values are unobtainable.** They are in the main article.
   ACS subscription; `elink` returns only `pubmed_pmc_refs` (17 citing papers)
   with **no `pubmed_pmc` self-link**; `esearch db=pmc` on the DOI returns
   **count=0**; `find_repo.py` finds no author repository; the repo-wide search
   finds no copy. Details and the exact queries in `00_acquisition.json`.
3. **No structure exists either.** Zero of the 16 lipidoids are enumerable:
   the SI names, draws and formulates no amine. A row would carry a value and
   no molecule.

STATUS: OPEN

The decision needed is whether to obtain the main article through the
institutional route that already supplied the SI (its watermark shows it came
from a subscription session), or to close this paper as not-curatable and mark
it `cite` rather than `added` in `lit_review`. Nothing else in this folder can
move until that is settled.

---

### D-2 — Which of `16I` and `16II` is the unsaturated lipidoid

`Figure S2`'s caption lists `16I/DNA and 16II/DNA`; the same figure's legend
lists `Unsaturated Lipidoid` then `Saturated Lipidoid`. Nothing states that
the two orders correspond, and the SI defines neither numeral.

What is certain: lipidoid 16 exists in two forms and the numeral distinguishes
the oleyl series from the n-octadecyl series.

Two inferences are available and they point **opposite ways**:

* `Figure S1` shows **16II** with both DNA and mRNA, which looks like the lead
  compound, and the abstract says the unsaturated lipidoids are superior →
  16II = unsaturated.
* The same abstract ascribes the difference to *"the large, tightly packed
  lipoplexes of **saturated** lipidoids"*, which is a reason to show a TEM of
  the saturated one → 16II = saturated.

Not assigned. Guessing would attach a label to the wrong series, which is the
`AC_2025` failure mode in miniature (there, the first confident repair would
have put 419 measured values on the wrong molecules).

STATUS: OPEN

Resolvable from the main article's nomenclature paragraph, which the SI does
not reproduce.

---

### D-3 — Is the saturated comparator an acryl*amide* or an acryl*ate*?

The abstract names the unsaturated tail explicitly (`oleyl acrylamide`) but
names the saturated comparator only as *"saturated lipidoids with n-octadecyl
tails"*. `01_blocks.csv` therefore carries `octadecyl_acrylamide` at
`confidence = medium` with `structure_assumed` in substance: the **linkage** is
parallel construction, not a printed fact.

The corpus supplies real but *indirect* support: `AA_2008`'s `*_N18` family is
44 lipids (43 distinct structures) that are bis-Michael adducts of
N-octadecylacrylamide, and **43/43 of them liberate exactly the structure this
folder built from the name** (`02_corpus_overlap.json` → `d1`). That attests
*the structure for the name*. It does **not** attest that this paper used the
amide rather than the ester — `AA_2008` runs both `N` (amide) and `O` (ester)
families over the same amines.

STATUS: OPEN

Resolvable from the main article's synthesis description. Until then the
`medium` confidence stands and no row depends on it.

---

### D-4 — The SI's licence is unrecorded, and `isOpenAccess = N`

`00_acquisition.json` records the gap: no landing page was visited, so no
licence string was captured, and Europe PMC reports the article as not open
access. `CLAUDE.md` §8 requires every constituent source to be
CC0/CC-BY-compatible before corpus release.

STATUS: OPEN

This is a release-blocking question for the corpus, not a curation question for
this folder, and it would apply to any future merge of this paper.

---

### D-5 — Figure S3's 16 viability values: keep the measurement, or delete it?

The panel is digitised exactly (eight checks, all passing — three-way bar-top
agreement to 0.00269 pt, 6/6 held-out minor ticks, the plot box exactly 16 bins
wide to 0.00048 pt) and **none of it is admissible**. The values sit in
`05_values.csv` with `corpus_admissible = no` and a reason.

Arguments for keeping them: they are what pins the library size at 16; they are
the evidence that "we looked" is a measurement rather than an assertion; and
they would be the natural cross-check if the main article is ever obtained
(a viability panel and a transfection panel over the same 16 lipidoids at the
same dose is a ratio test).

Argument for deleting them: a non-admissible value in a curation folder has
been mistaken for an admissible one before.

STATUS: OPEN

Recommendation: keep, because `corpus_admissible` is an explicit per-row field
and `06_rows.csv` is independently empty. But it is the user's call whether
non-delivery panels are digitised at all.
