# DECISIONS — peerd_31999380

Ramishetti, Hazan-Halevy, Palakuri, … Peer.
*A Combinatorial Library of Lipid Nanoparticles for RNA Delivery to Leukocytes.*
Adv. Mater. **32**, e1906128 (2020). PMID 31999380, DOI 10.1002/adma.201906128.

`build_dataset.py` refuses to merge a folder with an `OPEN` entry. Nothing here
is resolved alone: PIPELINE_V2 §3 lists these as the classes that get queued,
and `AC_2025` is why — its first repair was confident, internally consistent,
and would have put 419 measured values on the wrong molecules.

---

### D-1. Does this paper enter the corpus at all, given that it yields zero admissible delivery values?

The 14 ionizable-lipid structures are solid: every one is built from the SI's
own IUPAC name and written synthesis, all 14 agree with their ESI-MS to ≤0.76
Da, 4 of 4 C10 anchors match exactly, and the corpus has none of them on any of
three keys. The **values** are the problem. This SI's only delivery-endpoint
panel is Figure S9, and it covers 2 lipids as antibody-targeted LNPs (D-6). The
library screen is in main-text Figures 1–4, which were not acquired.

So the choice is between (a) leaving the folder unmerged until the main article
is acquired, (b) merging the structures with no `Experiment_value`, which the
79-column schema has no row shape for, and (c) acquiring the main text and
curating the screen properly.

Recommendation: **(a), then (c)**. A structure with no measurement adds a row
the corpus cannot use and a lipid the chance baseline has to count.

STATUS: OPEN

---

### D-2. Double-bond POSITION in Lipids 4 and 12 rests on the SI's prose and is not independently testable here.

Lipid 4 uses *cis*-2-nonen-1-ol; Lipid 12 uses *cis*-3-nonen-1-ol. The two
arms are **C15H27O2 in both cases** — identical formula, identical exact mass,
identical degree of unsaturation. Mass cannot see it. The olefinic integral
(P4) counts C=C bonds and is blind to where they sit. The C10 anchor cannot
separate them either, because the anchor and the enumeration both take the
position from the *same* IUPAC name, so agreement between them is not
independent evidence.

What does separate them is the ester OCH₂ **chemical shift**, which is where the
allylic/homoallylic difference actually shows up. Measured off the SI text, all
five spectra in the two series:

| compound | ester OCH₂ | multiplicity |
|---|---|---|
| Lipid 4 Boc-precursor (ESI 609) | δ **4.61**, 4H | d |
| Lipid 4 free hydrazine (ESI 509.5) | δ **4.61**, 4H | d |
| Lipid 4 (ESI 622.6) | δ **4.61**, 4H | t *(see below)* |
| Lipid 12 precursor (ESI 538.7) | δ **4.07**, 4H | t |
| Lipid 12 (ESI 651.8) | δ **4.01**, 4H | t |

A 0.60 ppm gap, consistent across three independent spectra on one side and two
on the other, is the allylic ester (`-O-CH₂-CH=CH-`) sitting downfield of the
homoallylic one (`-O-CH₂-CH₂-CH=CH-`). Both of Lipid 4's precursors also print
the **doublet** an allylic OCH₂ must give, against the triplet seen throughout
the Lipid 12 series.

Lipid 4's own line prints "t" where its two precursors print "d". That is the
*third* defect in that one NMR line — C-05 is the second (its olefinic integral
says four C=C where mass and both precursors say two) — so the line is
unreliable, not the structure.

This is still a reading of five numbers, by me, and it has not had a second
pass.

Recommendation: accept the assignment, record it as `confidence = high` with a
`curation_note` naming the δ 4.61 / δ 4.01 evidence, and have a reviewer
confirm the five shifts before merge.

STATUS: OPEN

---

### D-3. Lipid 9's printed ¹H NMR is not obviously the spectrum of the ether product.

Lipid 9 is the one ether in the library:
(linoleyl)₂N-CH₂CH₂-O-(CH₂)₄-NMe₂, made by displacing EA-1's mesylate with
4-(dimethylamino)butan-1-olate (SI p.27). Its mass fits — found 658.0 against
an average [M+H]⁺ of 658.18 — and its total integral is exactly right at 84/84.

But the dimethylamino singlet is in the wrong place. Every other
dimethylamino-bearing lipid in this SI reports it near **δ 2.2 (6H, s)**.
Lipid 9 reports **δ 3.36 (6H, s)** and has no singlet at 2.2 at all. A 6H
singlet that far downfield is what a *quaternised* or strongly deshielded
N-methyl pair looks like — and this is the one synthesis in the paper that
exposes a tertiary amine to an alkyl mesylate, where over-alkylation to a
quaternary ammonium is the obvious side reaction.

The ether linkage is also the **only connectivity in this paper with no C10
anchor** (`reference_structures.csv` says so on the Lipid 9 row).

Recommendation: emit the neutral ether as drawn and described, with
`confidence = medium`, `structure_assumed = 1` and a `curation_note` that names
the 3.36 ppm singlet. Do **not** emit a quaternary ammonium — the mass argues
against it and the paper never claims one. Flag, do not guess.

STATUS: OPEN

---

### D-4. Figures S3 and S4 disagree about which lipid a colour means, and the gap is too small for any colour check to catch.

Contradiction C-01, measured: S3's "Lipid 13" swatch is RGB (127, 4, 127) and
S4's "Lipid 12" swatch is RGB (128, 0, 129). **Euclidean distance 4.58**, below
`heatmap.MAX_DIST`. The two legends also list the right column in different
orders (…11, **12, 13**, 14 vs …11, **13, 12**, 14), and S4's drawn bar order
follows S4's legend, so the plotted data order differs too.

This curation never matched across panels, so nothing we produced is affected.
It matters for anyone who later reads S3 — including whoever curates the main
text, where the same palette may recur.

Recommendation: record it in `ERROR_CLASSES.md`. The generalisable lesson is
that a legend is **per panel**, and that a colour LUT carried between panels can
be wrong by less than the distance threshold that exists to detect exactly that.

STATUS: OPEN

---

### D-5. Figures S2 and S4 report the same condition and disagree about four of twelve lipids.

Contradiction C-04. Eight of the twelve shared lipids agree to within 0.3%
(ratios 0.9971–1.0034), which simultaneously validates the digitisation and
shows the two panels largely share one dataset. Against that baseline the other
four cannot be noise: Lipid 2 ×0.943, Lipid 3 ×**1.313** (93.75% vs 123.10%, a
29.35 pp gap, confirmed by eye as the tallest bar in S4's 0.06 µM group),
Lipid 4 ×1.044, Lipid 5 ×1.110. Spearman across the twelve is only 0.48.

Neither panel is corpus-admissible, so nothing hangs on resolving it. It is
recorded because it is evidence about **this lab's figure preparation**, and
that bears on how much weight the main-text screen should carry when it is
curated.

Recommendation: no action on the data. Note it in the folder and raise it if
the main-text panels show the same pattern.

STATUS: OPEN

---

### D-6. Is Figure S9 — the only delivery readout in the SI — worth curating?

Figure S9 is in vivo CD45 silencing in CD4⁺ and CD8⁺ T lymphocytes, spleen and
lymph node, for lipid 2 and lipid 6 **targeted** LNPs vs Mock. It is in kind the
same quantity as the corpus's `LRP6_knockdown_normalized` and
`protein_adundance_normalized`. Three things stop it:

1. **n = 2 ionizable lipids.** A per-experiment z-score over two points is ±1 by
   construction and carries no rank information — it would add an
   `Experiment_ID` to the chance baseline while contributing nothing to any
   shared-lipid comparison.
2. **tLNPs.** Both are antibody-decorated. The 79-column schema represents a
   targeting moiety only through `fifthcomponent_*`, and the paper gives no
   molar ratio for it.
3. **Box-and-whisker, not bars.** `bars.bar_top` requires a full-width top rule
   with both flanks agreeing, down to a baseline; a box's lower edge floats, so
   it refuses. Measuring a box top would need bespoke geometry, which this
   curation was instructed not to write — and which is error class 49, the
   reason `curation/lib` exists.

Recommendation: **do not digitise S9.** If the main article is acquired, the
untargeted library screen in Figures 1–4 is the thing to curate, and S9 becomes
a supporting detail.

STATUS: OPEN

---

### D-7. Two numbers printed on one line give Lipid 13 two different masses.

Contradiction C-06: `880.12 ([M+1]+)` implies M = 879.11 while `441.75
([M/2+1]+)` on the same line implies M = 881.49 — the paper disagrees with
itself by **2.37 Da**, with no structure of ours involved. The other three
doubly-charged ions in the SI (Lipids 2, 9, 12) agree with their own structures
to 0.03, 0.16 and 0.33 Da, so the method is normally sound.

Our C54H106N2O6 matches the singly-charged ion (monoisotopic [M+H]⁺ 879.812 vs
880.12 printed) and is independently supported by the C10 anchors for Lipid 7
(tridecan-7-yl decanoate arms) and Lipid 14 (the N-CH₂CH₂-O head ester), which
between them cover every motif in Lipid 13.

Recommendation: keep the structure, anchored on the singly-charged ion; carry
C-06 forward so the discrepancy is visible to anyone who re-derives the mass.

STATUS: OPEN

---

### D-8. Two cells in Table S1 are not usable numbers, and no reading of them is defensible.

- Lipid 10 / 2 h haematocrit: **`51.300 ±0.2.3`** — "0.2.3" is not a number.
  "0.23", "0.2" and "2.3" are all plausible and the SI gives no way to choose.
- Lipid 8 / 2 h platelets: **`1020.2 ±0.18`** — the other ten platelet SDs in
  that column run 72.42 to 393.2. The haematocrit cell in the *same row* also
  ends in 0.18, which looks like a copy, but that is a hypothesis about a
  spreadsheet, not a reading of the paper.

Both are emitted as `value_raw_text` with an `extraction_note`, never as parsed
numbers. Neither is corpus-admissible in any case (toxicology).

Recommendation: leave both unrepaired. An implausible number is not a licence
to invent a plausible one.

STATUS: OPEN
