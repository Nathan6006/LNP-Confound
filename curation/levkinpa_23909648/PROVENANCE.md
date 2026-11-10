# PROVENANCE — LL_2013

Li, Wang, Wu, Davidson, Levkin. *Combinatorial synthesis and high-throughput
screening of alkyl amines for nonviral gene delivery.* Bioconjugate Chemistry
**2013**, 24(9), 1543–1551. PMID 23909648, DOI 10.1021/bc400158w.

Curated 2026-11-10. **NOT MERGED into `data/`** — see `DECISIONS.md`, which
carries eight `OPEN` entries; `build_dataset.py` refuses such a folder, and
D-1 and D-2 block a merge outright.

---

## 1. The endpoint, verified before any extraction work

**Average transfection efficiency, printed as "the percentage of transfected
cells", with an SD. Cargo: plasmid DNA. Lipoplex with DOPE. In vitro.**

The tracker's title was truncated ("…Alkylated…"), and a sibling in this run
burned its budget on a paper whose truncation hid an antimicrobial endpoint.
Checked first, via Europe PMC `EXT_ID:23909648`: the full title ends "…alkyl
amines for nonviral gene delivery". Then confirmed against the SI's own text
layer — `transfect` 4 hits (the four table captions), `plasmid` 2, and
`luciferase` / `luminesc` / `RLU` / `GFP` / `siRNA` / `knockdown` **0 hits
each**. This is a flow-style percentage, not luminescence and not knockdown;
it is stored exactly as printed and is **not** log-transformed, converted or
z-scored by this curation.

## 2. Where the values are

**Tables S3, S4, S5 and S6, SI pp. 18–19 — all four in the PDF's text layer.**

**Tier 4 on every value; zero pixels, zero digitisation.** No axis, no colour
ramp, no bar, so none of the geometric error classes (20, 27, 32, 33, 35) can
apply here.

Two things worth the next reader's attention:

- **The SI's own List of Contents (p.2) omits all four tables.** It stops at
  Figure S2. Every one of the 87 delivery values lives in a table the contents
  page does not mention (`08_contradictions.csv` C-4). A triage that trusted
  the contents page would have reported this paper as having no endpoint.
- **`pdfplumber.find_tables()` returns 0 on pages 4, 5, 18 and 19** — *not*
  because the tables are images. Every cell is a text object; the rules are
  drawn as filled rects 0.36–0.48 pt tall rather than `line` objects, so the
  default `lines` strategy finds no edges. The fix is to hand it the rect runs
  as explicit column edges, not to reach for pixels.

## 3. What the paper actually contains — the three numbers

| | count | what it is |
|---|---|---|
| **Distinct structures** | **3** determined | `1C12-D`, `1C16-S`, `1C16-D` — column-purified, full IUPAC names, 1H/13C NMR, HRMS, and drawn in full. |
| **Formulations** | **9** fully specified | 3 lipidoid:DOPE ratios × 3 plasmid-DNA doses (Table S5). Plus 3 ratios at an unprinted dose (S4), 6 blend compositions (S6), and 1 unstated condition (S3). |
| **Individually measured** | **87** printed values | over **19** distinct test articles (12 crude lipidoids, 3 purified compounds, 4 binary blends) plus **6** L2K controls. |

**The library is 200 members and that number is confirmed — but it does not
mean 200 structures.** Table S1 resolves 20 amines × 10 alkyl halides = 200,
matching the tracker exactly (gate `Q-2`). Each of the 200 is a **crude
one-pot reaction product**: Table S1 reports *several different tail counts
for a single (amine, alkyl halide) pair*, up to six, and **all 23 of its rows
show more than one** (gate `Q-3`). The paper says so itself — "The 1C12-D
product was isolated from the crude lipidoid mixture by column
chromatography". So the screened material is a mono/di/tri/… mixture in an
unreported ratio.

**Only 12 of the 200 members carry a value anywhere in this SI.** The other
188 were screened (that is the paper's Figure 1) and are not reported here.

## 4. The structure blocker

**Nine of the twelve valued lipidoids have no amine structure in this SI.**
The 20 heads are bare integers in Table S1's `Amine` column; the SI names
exactly one (amine 1 = 1-(2-aminoethyl)pyrrolidine, p.3) and draws no amine
panel (`Scheme` 0 hits, `structure` 0 hits over the full text layer). This is
the same shape of blocker the sibling `levkinpa_29558113` hit, and it is why
that folder yielded 5 rows.

Each unknown amine's neutral mass *is* solvable from Table S1 and those masses
are recorded in `01_blocks.json` — but **a mass is not a structure** and no
structure was guessed from one. `DECISIONS.md` D-1.

## 5. Rows emitted

**29**, all from Tables S5 and S6, all on the three purified compounds
(27 = 3 compounds × 3 DOPE ratios × 3 doses, plus 2 unblended Table S6 rows).
58 values were skipped, each with its own reason on the row in
`06_skipped.csv`: 48 crude/unknown-head, 6 L2K (a proprietary mixture with no
defined structure), 4 binary blends (no single `IL_SMILES`). 29 + 58 = 87,
gated at `Q-10`.

Quality columns written by `curation/lib/emit.emit_rows`, the single writer.
`confidence` is **derived**, not literal — `_confidence()` scores each
structure against three named witnesses (printed IUPAC name, HRMS within 5
ppm, 1H integral sum within 3 H) and returns the note with it; a literal would
be refused, which is F-011's lesson. All 29 come out `high`.
**`rows_with_a_caveat` = 1**: the single `value_is_bound = upper` row, Table
S5's `1C16-S:DOPE=1:2` at 150 ng, printed `0 ± 0`, i.e. a non-detect whose
true value is at or below the detection limit — and which cannot be
log-transformed at all (D-5).

`value_normalization_unit` and `Experiment_value_within_model` were left
**unassessed** on purpose: no standardisation was applied, and writing a unit
would assert one that was never performed.

## 6. Corpus overlap — reported, not acted on

Against the live corpus (29,195 vitro + 1,803 vivo, 19,919 unique
`IL_SMILES`): **0 of 3 match by full InChIKey, 0 by skeleton block, 0
`IL_name` collisions.** So, unlike four siblings this run, nothing here is
already present.

But the **amine-1 head group appears on 638 corpus rows across 5 experiments**
(`KZ_2016`, `YX_2024`, `YR_2024`, `XH_2025`, `YX_2026_SDL`), always carrying
far more elaborate ester/amide tails. These three compounds are the
unelaborated version of a head the corpus already spans. Bearing on D-8: this
experiment would contribute **nothing** to Result 3.

## 7. The anchor — posture recorded BY HAND

`provenance/structure_check.py`, `check_anchor_independence.py` and
`check_anchor_coverage.py` were **not run** (out of scope for this task), so
C10/C16/C18 have not passed on this folder. What follows is the hand record.

**The paper draws all three purified compounds in full**, SI pp. 7, 10 and 13,
each above a printed `Chemical Formula` / `Exact Mass` / `Molecular Weight`
block. So an anchor is possible for every structure this folder emits, and
`reference_structures.csv` carries **three anchors and no empty rows**.

**Independence.** There are three structure-building routes in this folder and
no two share a template:

| route | method | template for 1C12-D |
|---|---|---|
| `stage03_enumerate.py` | RDKit fragment assembly, `curation/lib/assemble.attach`, atom-mapped | head `C1CCN(C1)CC[NH2:1]` + `[CH2:2]`+`C`×11, one attach per tail |
| `stage07_qc.py` | string composition, **tail-first** | `f"{alk}N({alk})CCN1CCCC1"` |
| `make_reference.py` | string composition, **head-first**, written from the drawing | `C1CCCN1CCN(CCCCCCCCCCCC)CCCCCCCCCCCC` |

The anchor generator's composed template is **byte-different** from the
enumeration's canonical SMILES
(`CCCCCCCCCCCCN(CCCCCCCCCCCC)CCN1CCCC1`) for all three compounds, and all
three canonicalise to the same molecule. Gate `Q-13` checks **the template**,
not the emitted anchor, because byte-identity of the *emitted* anchor is not
the signal — a correct re-derivation canonicalises to the same string as a
correct enumeration, which is exactly why C10's own `anchor_may_be_a_copy`
missed F-037. Seed `S26` replays F-037 by making the template a byte copy; the
emitted anchor is unchanged and still correct, and `Q-13` still fires (6 → 4).

**No stereocentres anywhere** in this chemistry (asserted per structure), so
the F-013 chiral-tag class cannot arise here. The connectivity class still
can, which is why §8's geometric check exists.

## 8. Checks — 22 gates, 26 seeds

`python3 -I stage07_qc.py` → **22/22 pass, exit 0**.
`python3 -I 10_seed_checks.py` → **26/26 seeds credited, 22/22 gates
exercised, exit 0**, with an unmutated control copy that still passes.

A seed is credited only if the gate command **exits non-zero** — read from
`subprocess.returncode`, never through a pipe — **and** the named gate's
agreement count **moves**. Every numeric seed is multiplicative (×1.01 or
more), because a sibling's +1-unit seed was absorbed by a loose tolerance and
proved nothing.

The mass evidence is split **seven** ways, and each leg can fail alone:

| leg | what it compares | result |
|---|---|---|
| `Q-5a` | our [M+H]+ vs **their calculated** HRMS m/z | −0.053 / −0.072 / −0.006 ppm |
| `Q-5b` | **their calculated vs their own observed** | −0.100 / −0.100 / +0.340 mDa |
| `Q-5c` | our [M+H]+ vs **their observed** | −0.275 / −0.367 / +0.597 ppm |
| `Q-5d` | **unprinted relations among their three HRMS values** | −0.021 and −0.021 mDa |
| `Q-5e` | our neutral masses vs **Table S1's calc column** | ≤ 0.046 Da |
| `Q-5f` | our formula vs the **drawn** formula | 3/3 identical |
| `Q-5g` | our neutral exact mass vs the **drawn** Exact Mass | ≤ 0.0005 mDa |

**Quantities the authors never printed** (CLAUDE.md §5 prefers these to any
recomputed marginal):

- **`Q-4`** — the per-tail increment inside their own Table S1 calc column
  must equal exactly one CₙH₂ₙ. **38/38** consecutive steps hold, worst
  residual 0.0904 Da against a 1-decimal column. Pure chemistry; nothing of
  ours enters.
- **`Q-5d`** — two differences among their three printed HRMS values:
  1C16-D − 1C16-S = 1 × C16H32 (224.25038 vs 224.25040) and
  1C16-D − 1C12-D = 2 × C4H8 (112.12518 vs 112.12520). **These only close if
  the tail multiplicities are 2, 1, 2**, so they confirm the mono/di
  assignment independently of any structure we built.
- **`Q-7`** — the DNA dose of Figures 2A, 3 and 8, which is printed nowhere.
  Table S5's L2K control is 6.9 ± 1.7 / 14.0 ± 1.8 / 22.4 ± 4.0 at 75/150/200
  ng (three distinct values); Tables S3, S4 and S6 each print **14.0 ± 1.8**,
  matching the 150 ng block in mean *and* SD and matching neither other dose.
  Recorded as an **inference** and queued as D-4; the dose is left blank on
  the rows.

**`Q-12` is the only check here that is not arithmetic over a printed
number.** Everything else would still agree if the picture showed something
else. The drawings are vector, with bonds exported as **filled outlines**, so
`chains.vector_bonds` in its default stroked mode correctly returns **0** —
the library refusing rather than inventing. In `include_filled=True` mode it
returns **exactly 2.0 segments per bond** on all three pages (64/48/80), being
the two long edges of each outline. The library's docstring warns that filled
mode also returns glyph outlines as bonds; **no separation was needed here,
and that is measured rather than assumed** — every glyph on these pages
(the `N`, the `H`, the compound label, the whole formula block) is a bitmap
*image*, not a path, which is also why none of that text is in the text layer.
The exact 2:1 ratio is itself the evidence that nothing leaked in.

The closed identity, as a sibling used it today: ChemDraw emits one path per
heavy-atom bond plus one per double bond (kekulised). This chemistry is fully
saturated, so n_double = 0 and the identity reduces to drawn = heavy bonds —
but the term is computed, not dropped. **32 / 24 / 40 drawn bonds against 32 /
24 / 40 expected, 3/3 exact.** A C12 tail mistaken for C14 moves this by 2; a
mono- mistaken for a di-alkylamine moves it by 16.

## 9. What could NOT be verified

- **The cell line, the reporter gene and the assay timepoint.** Absent from
  the SI; `Model_target` and `Model_type` are written `NOT STATED IN THE SI`
  and were not inferred from the lab's other papers (D-2).
- **The amine heads 2–20** (D-1) — only their masses are recoverable.
- **The composition of every crude library member** — the mono:di ratio is
  never reported, so no crude entry can be given a structure.
- **Which of Tables S5 and S6 is right** where they disagree 5.0× on `1C16-S`
  (D-6).
- **Table S6's lipidoid:DOPE ratio** — not printed.
- **The licence** — ACS, not CC-BY (D-7).
- **C10 / C16 / C18 have not been run** on this folder; §7 is a hand record.

## 10. Reproduce

```bash
cd curation/levkinpa_23909648
python3 -I stage01_blocks.py        # Table S1: library shape, masses
python3 -I stage02_drawn.py         # drawn-bond geometry
python3 -I stage03_enumerate.py     # the 3 structures
python3 -I stage04_corpus_overlap.py
python3 -I stage05_values.py        # Tables S3-S6 -> 05_values.csv
python3 -I stage06_emit.py          # join + quality columns
python3 -I make_reference.py        # C10 anchors
python3 -I stage07_qc.py            # THE GATE: 22/22, exit 0
python3 -I 10_seed_checks.py        # 26/26 seeds credited, exit 0
```
