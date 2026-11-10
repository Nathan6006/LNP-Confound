# Triage of the four hand-downloaded papers — one extracted, three below the bar

The user set the threshold: **extract only if the paper yields more than ~75
datapoints.** All four files arrived intact; naming differed only in case
(`suppmat` vs `SuppMat`), which is cosmetic — the real names are recorded in
`curation/SOURCES.md`. Every number below was **measured from the files**, not
inferred from the abstract.

| folder | library | measured yield | verdict |
|---|---|---|---|
| **`liu_42585110`** | **40-member RAIL** | **~41 in vitro + ~82 in vivo ≈ 123** | **EXTRACTING** |
| `leejy_42610355` | 19 ionizable lipids | ~19 + 3 benchmarks ≈ **22** | below bar |
| `cheng_41656866` | **10** PILs | ≈ 10 × a few organs | below bar |
| `xuq_22148515` | combinatorial, crude → refined | ~1 attributable lead | below bar |

## 1. `liu_42585110` — EXTRACTED

Fig 1 caption, verbatim: *"(a) Luciferase expression in DC2.4 cells … with
FLuc mRNA-loaded Dlin-MC3 LNPs and **40 RAIL-LNP formulations**. (c)
Quantification of luciferase expression in the **liver and spleen**."* The text
describes *"a modular **40-member** library"* from *"systematic variation of the
ionizable scaffold, linker architecture, and retinoid-derived hydrophobic
domain"*.

So **40 distinct structures**, not a formulation screen — the first paper in
some time that is not. Fig 1d is spleen-to-liver **ratios**, which §8 already
rules non-poolable (`mTP_Difference`), so only the absolute 1c values are
candidates. Retinoic-acid-derived chemistry should have **zero** corpus
analogue.

## 2. `leejy_42610355` — below the bar, but the closest call

*"A library of **19** ionizable lipids was synthesized in an efficient
three-step protocol"*, screened in Jurkat T cells for eGFP. **One condition**:
`40 : 48 : 10 : 2` molar ratio, **N/P 8**, **100 ng mRNA per formulation** —
stated as held constant across the library and the benchmarks. So ~19 members
plus SM-102 / ALC-0315 / MC3 ≈ **22 datapoints**.

**This one is worth a second look despite the count**, and the reason should be
on the record:

- it is **CC BY** — the **only one of the four whose licence passes §8**, where
  the other three are Wiley/ACS with no CC terms;
- **19 novel structures** is more chemical support than most of this run's
  large folders supplied: seven papers were one-lipid formulation screens
  contributing **one** structure each, and the single largest extraction
  (5,112 rows) also contributed **one**;
- **Jurkat T cells** would be a new environment — worth checking against the
  corpus's 12 cell types.

By rows it is small; by the quantity this project actually measures — distinct
chemistry on which two labs can be compared — it is worth more than its row
count. **Flagged for the user, not overridden.**

## 3. `cheng_41656866` — below the bar

*"we synthesized a library of **10** structurally diverse PILs"*, and again
*"a systematic evaluation … using a library of 10 structurally diverse LNPs"*.
Figs 4-6 are all single-lead work on `B12-a13Dab4` (comparison with SM-102,
liver gene editing, an HT-1 rescue model). Only **3** distinct compound labels
appear in the main text. Ten structures across a few organs does not reach 75.

## 4. `xuq_22148515` — below the bar, and the 2012 paper has a further problem

*"A combinatorial library of lipidoids was constructed"*, screened **crude** in
HeLa and then refined — but only **one** lipidoid label (`14O-87`) appears in
the whole PDF, and the figures are the only value source (PMC declares
`hasSuppl=N`; its JATS has 0 `<table-wrap>`). This is the "crude one-pot
product" shape measured on `levkinpa_23909648`, where 200 reaction products
yielded 12 valued members and **9 of those 12 had no disclosed head structure**.
Low yield and high structural uncertainty together.

## Note on the files themselves

All four validated by magic bytes. Two points worth keeping:

- `liu_42585110`'s SI has **`<w:tbl>` = 0 but 83 images** in `word/media/` —
  the EMF-table trap `goulart_41138890` hit, where a tag count reported one
  table for a document holding four, recovered from 694/534/1,571 glyph
  records. Flagged to the agent.
- `cheng_41656866`'s 44 MB SI has **1,196 images across 47 pages** and only
  34,509 text characters — i.e. it is almost entirely raster, which is the
  expensive tier even if the count had cleared.
