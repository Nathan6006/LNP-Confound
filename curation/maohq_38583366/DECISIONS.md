# DECISIONS — maohq_38583366 (Zhu et al., *Biomaterials* 2024, PMID 38583366)

The pipeline detects that a source is silent; it does not get to decide what
the source meant. Everything below is queued, not resolved.

**Nothing has been merged into `data/`. `06_rows.csv` is header-only.**

### D-1 — The in-vitro arm is ALREADY IN THE CORPUS. Do not emit it.

Supplementary Tables S1–S6 print 192 `RLU in HepG2` values (6 helpers × top 32).
Every one of them is already in `data/data_vitro.csv` as `YZ_2022`
(PMID 35879315), which this paper cites as ref [16] for that screen:

* corpus `YZ_2022` joins the Mao lab's deposited 1,080-row table on
  **composition only** (helper + four component percents to 1 dp + N:P), **1,080
  of 1,080 rows**, Spearman(`Experiment_value`, `RLU_HepG2`) = **1.000000**;
* the ln-unit residual rms is **2.8531e-4** against **2.8868e-4** predicted from
  three-decimal rounding alone (ratio 0.9884) — the discrepancy *is* the
  rounding in the deposit, which is why an exact-equality join fails and is not
  the right test;
* the paper's 192 values are present in that column as a multiset **192/192**,
  and are exactly each helper's top 32.

Emitting them under a new `Experiment_ID` would duplicate a merged experiment
and manufacture exactly the source-separability signal this project measures.

STATUS: RESOLVED — not emitted.

### D-2 — Should the 20 individually measured in-vivo values be emitted at all?

Fig 4B gives 20 **individually measured** liver fluxes (intraduodenal, 100 µg
pDNA), digitised here and reported in `05_values.csv`. They are a translated
luciferase readout, so they are admissible in kind. Two reasons not to emit,
and one reason someone might want to:

* **Against.** The paper contributes **one** ionizable lipid, DLin-MC3-DMA,
  already in the corpus on 3,098 rows. Within this source the 20 values differ
  only by *formulation*, so the within-source normalised rank that Result 3 and
  Stage 5 consume would be a formulation rank of a single benchmark molecule,
  not a chemistry rank. `CLAUDE.md` §8 already flags benchmark-only pairs as a
  hazard for cross-lab concordance; this is that hazard in its purest form,
  because the benchmark is the *only* lipid.
* **Against.** Every value is tier 5 (pixels). No Source Data file exists
  (mmc2–mmc4 are 404 on the Elsevier CDN in four extensions each), there is no
  repository, and no article PDF could be obtained.
* **For.** The in-vivo pool has only **two** structures shared across
  experiments. DLin-MC3-DMA is one of them, present in six in-vivo sources, and
  **three of those already carry liver** (`NG_2025_2`, `ZC_2023`, `ZH_2023`). A
  seventh liver source would add organ-matched pairs, which is rare and
  load-bearing — unlike SM-102, whose four sources are ear and whole-body only.
  The precedent also exists: `AC_2025` is itself a 498-row one-lipid
  formulation screen already in the corpus.

STATUS: OPEN

### D-3 — `Experiment_ID` collision. `YZ_2024` is taken.

The first author is **Zhu Yining**, so the key derives to `YZ_2024` — and
`YZ_2024` already exists in `data_vitro.csv` as a *different* paper,
**PMID 38082180** (mRNA, 1,080 rows). `YZ_2022` is PMID 35879315, also Zhu Y.
The documented `_2` suffix is for a *second screen from one paper*, not for a
*different paper* whose initials and year collide, so there is no convention to
apply. This must be settled before anything from this paper is ever emitted.

STATUS: OPEN

### D-4 — `18PG` (F-057). This paper cannot settle it.

The complete and verbatim statement is *"a mixture of helper lipids (DOTAP,
DDAB, DOPE, DSPC, 14PA, 18PG) from Avanti"*. No chemical name, no catalogue
number, no acyl specification; no helper lipid is drawn anywhere in the article
or the SI. A search for `stearoyl`, `oleoyl`, `glycero`, `phospho`, `18:0`,
`18:1`, `DOPG`, `SOPG`, `catalog` and `Cat` over **60,875** characters of JATS
full text and **17,600** characters of SI text returns **zero hits for every
one**. The paper supports neither SOPG nor DOPG.

`data/` is unchanged. Recorded so the next reader does not re-open this paper
hoping for an answer.

STATUS: RESOLVED (as a null result) — no change to `data/`.

### D-5 — Fig 4B's sixth named formulation does not clear the stated threshold.

`08_contradictions.csv` C-4. The paper's six named formulations are exactly our
six largest bars (6/6), but `DI-α3` reads 4.47 × 10⁴ against the stated
"exceeding 5 × 10⁴", a 10.6% shortfall — about 13 px on a panel whose worst
held-out tick residual is 0.48 px. We do not adjust a read to meet a sentence.
Nothing depends on it here because nothing is emitted.

STATUS: OPEN
