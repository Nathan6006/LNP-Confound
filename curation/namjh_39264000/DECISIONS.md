# DECISIONS — namjh_39264000

Questions this curation must NOT answer on its own. `build_dataset.py` refuses
to merge this folder while any entry below is `OPEN`.

WHY THIS FILE GATES THE MERGE. The pipeline is good at detecting that a source
is silent. It is not entitled to resolve it. In `AC_2025` the first repair was
confident, internally consistent, and would have assigned **419 measured values
to the wrong molecules** — the names were wrong and the structures were right,
and nothing in the data said which. So the rule is: detect, queue, and stop.

## Format

```
### D-n  <one line>
STATUS:    OPEN | RESOLVED <date> by <who>
AMBIGUITY: what the source does not say.
CANDIDATES: each reading, and what it implies for the emitted rows.
EVIDENCE NEEDED: what would settle it.
DECISION:  (when resolved) what was chosen and why. Keep the entry.
```

## Entries

### D-1  The paper's headline is 213 LNPs; only 36 ionizable lipids are drawn anywhere.
STATUS: OPEN
AMBIGUITY: The tracker records this paper as **213 lipids**. The paper says
"213 LNPs", built from "14 types of FHs" and "14 types of EB–T chains". Those
14×14 blocks are drawn in no figure and named in no procedure: Figure S1 is a
generic scaffold with R1 and R2, and Figure 3a draws 6 heads, 3 bridges and 2
tails — the 36-member subset and nothing more. Figure 2's nine printed group
sizes partition 213 across the six compositions of Table S2 (76 / 25 / 26 / 76
/ 3 / 7, all nine partitions reconciling), which establishes that **213 counts
formulations, not molecules**. The number of distinct ionizable lipids in the
study is stated nowhere and is bounded only by 76 ≤ n ≤ 196.
CANDIDATES:
 (a) Emit the 36 structures with the 108 measurements that exist for them, and
     record that 177 of the 213 LNPs are unrecoverable. **This is what the
     folder currently emits.**
 (b) Emit nothing, on the grounds that a 36-of-213 slice of a screen is a
     biased subset — these six heads were chosen by the authors *because* the
     ML ranked their substructures highly, so the slice is enriched for
     actives and is not a random sample of the screen.
 (c) Request the full 213-row table from the corresponding author (the Data
     Availability Statement invites exactly this) and curate it properly.
EVIDENCE NEEDED: either the 213-row dataset, or a decision on whether a
selection-biased 36-member slice belongs in a corpus whose central claim is
about cross-laboratory rank agreement. Candidate (b) is not obviously wrong:
Result 3 compares *shared* lipids across sources, and a slice chosen for being
good performers is precisely the kind of subset that can inflate a
cross-source correlation for reasons unrelated to transferable SAR — the same
hazard already flagged in CLAUDE.md §8 for benchmark-only pairs.

### D-2  The `Experiment_ID`, and its collision risk with `SB_2024`.
STATUS: OPEN
AMBIGUITY: The corpus key is first-author initials plus year. The first author
is **Seo-Hyeon Bae**, giving `SB_2025`, which is what is emitted. But
`SB_2024` is already in the corpus — and it is not merely a name clash. Its 60
rows are `Mouse_ICR`, `Model_target = ear`, `Route = intradermal`,
`Dose = 5.0 µg`, `Mixing_method = microfluidics`, `Aqueous_buffer = citrate`:
**the same assay, the same strain, the same tissue and the same dose as this
paper's I.D. screen.** `SB_2024` is PMID 38779592 (Bioactive Materials 2024),
whose lipids are the `KIST_*` series.
CANDIDATES:
 (a) `SB_2025` / `SB_2025_2`, as emitted.
 (b) `JN_2025` after the corresponding author (Jae-Hwan Nam), to avoid the
     visual collision.
 (c) Something that marks the relationship explicitly, if the two papers share
     a laboratory.
EVIDENCE NEEDED: whether PMID 38779592 and this paper share an author or an
animal facility. **This matters well beyond naming.** The project's central
claim is that cross-laboratory prediction is confound-limited, and Result 3
rests on how identical lipids rank across *independent* sources. Two corpus
"sources" that are in fact one laboratory running one protocol are not two
environments, and neither the source-classification result nor the
shared-lipid concordance should treat them as such. This curation is not
entitled to decide that.

### D-3  The vaccine readouts are not emitted. Should they be?
STATUS: OPEN
AMBIGUITY: This is a vaccine paper and its primary biological claims rest on
total IgG, IgG1/IgG2a, IFN-γ ELISpot, TNF-α and IL-2 (Fig 4f, 4g, Fig S4).
None of it is emitted. The two emitted experiments are reporter-protein
expression: Renilla luciferase in the ear after an I.D. injection, and human
EPO in serum after an I.M. injection.
CANDIDATES:
 (a) Keep them out. They are a third experiment on every axis that matters:
     different cargo (HPV E6/E7 mRNA, not a reporter), different schedule
     (prime/boost a week apart, read a week after the boost, versus a single
     injection read at 6 h), a different and internally contradicted mouse
     strain (08_contradictions.csv row 6), and units that are antibody titre
     and spot-forming cells rather than a normalised signal. The corpus's
     `Experiment_method` vocabulary has nothing for either.
 (b) Emit them as a third `Experiment_ID` with a new `Experiment_method`.
     An antibody titre is a delivery-adjacent but not a delivery measurement,
     and pooling it with luminescence would put two different physical
     quantities on one z-scale.
EVIDENCE NEEDED: a project-level decision on whether immunogenicity endpoints
belong in a delivery corpus at all. Note that the panels also cover only a
subset of the 36 and are grouped rather than per-lipid in places.

### D-4  Table S2's "Mass ratio" cannot be mapped onto `IL_to_nucleicacid_massratio`.
STATUS: OPEN
AMBIGUITY: Table S2 gives a "Mass ratio" of 1 or 0.5 per LNP type, and Figure
2l confirms those are the only two values across all 213 LNPs. The corpus
column `IL_to_nucleicacid_massratio` holds an actual lipid-to-nucleic-acid mass
ratio — `SB_2024`, the most comparable source, carries 12.63 and 6.31 there.
A value of 1 or 0.5 is not that quantity in those units, and the paper never
defines what it is a ratio of.
CANDIDATES:
 (a) Leave it blank. **This is what is emitted**, on the F-004 precedent: three
     mutually inconsistent candidates for `N:P` were blanked rather than
     reconstructed.
 (b) Store 1.0 and 0.5 literally, accepting that the column then carries two
     different quantities across sources.
 (c) Derive it from N/P and the molar ratios, which requires an assumed mRNA
     nucleotide mass and an assumed ionizable-lipid charge fraction at pH 4 —
     i.e. two assumptions to recover one undefined number.
EVIDENCE NEEDED: a definition of "mass ratio" from the authors, or the
project's rule for a column whose units differ between sources. `N/P` is
unaffected and IS emitted (12 for type 1, 6 for type 4) because the corpus's
`IL_to_nucleicacid_chargeratio` is the same quantity.

### D-5  The source licence is CC BY-NC-**ND**, which is not CC-BY-compatible.
STATUS: OPEN
AMBIGUITY: The PMC record's `<permissions>` block states
`https://creativecommons.org/licenses/by-nc-nd/4.0/` — "no modifications or
adaptations are made". CLAUDE.md §8 requires every constituent source to be
CC0/CC-BY-compatible before corpus release. Digitising a figure into a value
table and z-scoring it is, on the strictest reading, an adaptation; the NC
clause also bites on any commercial redistribution of the corpus.
CANDIDATES:
 (a) Include the rows and carry the licence restriction on the source list.
 (b) Include only facts (numbers), on the view that data are not copyrightable
     and the ND clause attaches to the expression, not to the measurements.
 (c) Exclude the source from any released corpus while keeping it for
     internal analysis.
EVIDENCE NEEDED: the project's position on NC/ND sources generally. This is not
specific to this paper — it needs deciding once, for all Wiley and Elsevier
sources in the corpus.

### D-6  `Model_target` for the I.M. hEPO experiment: `serum` or `muscle`?
STATUS: OPEN
AMBIGUITY: hEPO is expressed at the injection site (muscle) and measured by
ELISA in **serum** (Methods, Serum Analysis). `serum` is emitted, because it is
where the measurement was made. The corpus has no `serum` target; it has 497
in-vivo rows at `muscle`, and it is not clear from here whether those are
expression-site or sampling-site labels.
CANDIDATES: (a) `serum`, as emitted; (b) `muscle`, matching the existing 497
rows; (c) a split convention, recorded once for the whole corpus.
EVIDENCE NEEDED: what the existing 497 `muscle` rows actually measured.

### D-7  Are Fig S2c and Fig S2e two formulations within one experiment, or two experiments?
STATUS: OPEN
AMBIGUITY: They are pooled into one `Experiment_ID` (`SB_2025`, 72 rows) and
z-scored together, so the roughly two-fold difference in their maxima survives
— deliberately, because F-016 is 5,040 rows whose between-assay offset was
deleted by standardising sub-blocks separately. But the Figure S2 caption says
each result was normalised to "the values of the Nil group included in **the
experiment in question**", which is the paper's own word for treating them as
separate experiments with separate controls.
CANDIDATES:
 (a) One experiment, 72 rows, pooled z. **Emitted.** `value_normalization_unit
     = experiment`. Preserves the type-1 vs type-4 offset.
 (b) Two experiments, 36 rows each. Deletes the offset, and makes the LNP-type
     comparison — which is one of the paper's own findings — unrecoverable
     from the corpus.
 (c) One experiment with `value_normalization_unit = experiment_x_model_type`,
     which the schema reserves for exactly this situation but describes as an
     error class (22) that is only legitimate where the offset is
     unrecoverable. Here it is recoverable, because both panels are already
     fold change over their own Nil.
EVIDENCE NEEDED: whether the type-1 and type-4 panels share a Nil control. If
they do, (a) is right without qualification. If they do not, the offset is a
comparison between two normalisations and (c) becomes arguable.

### D-8  36 of 36 structures are emitted; only 2 of 36 have an external anchor.
STATUS: OPEN
AMBIGUITY: 24 of the 36 members contain at least one block (FH-1, FH-2, FH-4,
FH-6, EB-6, EB-7) that the paper names in no procedure, no reagent list, no NMR
assignment and no mass spectrum. Their only evidence is the Figure 3a drawing.
They are emitted at `confidence = medium` with the reason on the row. The
tracer that measured them reproduces all five SI-named blocks exactly
(2-hexyl-1-decanol → 16 C, 2-octyl-1-dodecanol → 20 C, glutaric anhydride → 5
backbone C, 3,5-dihydroxybenzoic acid → locants [3,5], benzoyl →
unsubstituted), which is the strongest statement available and is still not an
anchor for the blocks themselves.
CANDIDATES:
 (a) Emit all 36 at `medium`, as now.
 (b) Emit only the 12 members whose every block is SI-named (FH-3 and FH-5 ×
     EB-5 × T-1/T-2), i.e. 4 lipids and 12 rows.
EVIDENCE NEEDED: the project's threshold for a structure whose only evidence is
a drawing that a calibrated tracer read. CLAUDE.md §6 lists "whether a
`confidence = low` library belongs in the corpus" as a queued decision; this is
the `medium` case of the same question.
