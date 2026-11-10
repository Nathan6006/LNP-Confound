# DECISIONS -- pabst_37602460

Questions this curation must NOT answer on its own. `build_dataset.py` refuses
to merge this folder while any entry below is `OPEN`.

WHY THIS FILE GATES THE MERGE. The pipeline is good at detecting that a source
is silent. It is not entitled to resolve it. In `AC_2025` the first repair was
confident, internally consistent, and would have assigned **419 measured values
to the wrong molecules** -- the names were wrong and the structures were right,
and nothing in the data said which. So the rule is: detect, queue, and stop.

## Format

```
### D-1  <one line>
STATUS:    OPEN | RESOLVED <date> by <who>
AMBIGUITY: what the source does not say.
CANDIDATES: each reading, and what it implies for the emitted rows.
EVIDENCE NEEDED: what would settle it.
DECISION:  (when resolved) what was chosen and why. Keep the entry.
```

## Entries

### D-1  The endpoint is antibacterial MIC, not delivery. Does anything here belong in the corpus?
STATUS: OPEN
AMBIGUITY: This paper is not a transfection screen. Its measured quantity is
growth inhibition of *Bacillus subtilis* in Luria Broth, read as OD420-580 over
24 h, with "active" meaning >99.9% growth inhibition. There is no nucleic-acid
cargo, no helper lipid, no cholesterol, no PEG-lipid, no nanoparticle and no
mammalian cell line anywhere in the assay: the quaternised lipidoid is dispersed
in bacterial growth medium by itself. The corpus's in-vitro pool is a
delivery-only pool by construction -- CLAUDE.md's 2026-07-23 entry records 90
`hemolysis_percent` rows and 77 out-of-domain cell-line rows being trimmed out
of it precisely to keep it so. An antibacterial MIC has no `Model_type`,
`Model_target` or `Cargo` value in that schema that is not an invention.
CANDIDATES:
 (a) Do not merge. The chemistry is in-family (polyamine + acrylate lipidoids;
     eleven members are already in the corpus, see D-2) but the readout is a
     different biological question, and pooling it would put a bacterial
     growth-inhibition ordinal on the same `Experiment_value` scale as
     luciferase transfection. Costs nothing that the paper's own data can give,
     since none of it is a delivery measurement.
 (b) Merge as a separate, clearly-labelled non-delivery experiment with its own
     `Model_type`. This creates the corpus's first non-delivery `Experiment_ID`
     and every downstream analysis (Stages 3-8) would then be pooling two
     incommensurable endpoints unless each stage is taught to exclude it.
 (c) Keep it standalone, as `molea.csv` is kept (CLAUDE.md section 8): extracted
     and documented in this folder, never merged, available as a
     *cross-endpoint* comparison point -- eleven of these lipidoids have a
     transfection value in `KW_2014`/`AA_2008`/`NC_2024` AND an antibacterial
     score here, which is a same-chemistry, different-endpoint concordance test
     of exactly the kind Result 3 is built from.
EVIDENCE NEEDED: the user's call on whether the in-vitro pool is
delivery-only by definition. Nothing in the paper can settle it; this is a
corpus-scope decision, not a reading of the source. Note that (c) is the only
option that uses the eleven-member overlap, and that the overlap is a fact
about the corpus, not about this paper.

### D-2  Which structure is "the lipid" -- the Michael adduct, or the quaternised cation that was tested?
STATUS: OPEN
AMBIGUITY: Two different molecules have a claim to the row. The Michael adduct
(polyamine + 4/5/6 acrylates) is the classic lipidoid and is what the corpus
already stores for eleven of these members. The species actually dosed is that
adduct after exhaustive methylation with CH3I: every nitrogen bears a methyl and
a permanent positive charge, with an iodide counterion each. Figure 1A draws the
quaternised product in full and the text says the methylation was done to make
the compounds dispersible; the unmethylated adducts "were poorly soluble in a
variety of aqueous solutions, including bacterial growth media" and were never
assayed.
CANDIDATES:
 (a) Emit the parent adduct as `IL_SMILES`. Then eleven rows carry a FULL
     InChIKey already present under `KW_2014` (8), `AA_2008` (2) and `NC_2024`
     (1) -- `2N14`, `2N16`, `3N13B`, `3N14`, `3N16`, `3NB14`, `3NB16`, `4N14`,
     `4N16`, `4NC14`, `4NC16`. They would become "shared lipids" for Result 3
     while carrying a value measured on a *chemically different* compound (the
     quaternised salt). That is the `AC_2025` failure mode in a new costume: a
     measured value attached to the wrong molecule.
 (b) Emit the quaternised cation. Chemically honest -- it is what was dosed --
     and the corpus overlap drops to zero (checked: 0 of 104 methylated cations
     match any corpus SMILES by full or skeleton InChIKey). Costs: the corpus
     has no convention for a tetracationic permanent quaternary ammonium with
     counterions, and `num_permanent_cationic_N` would be 2-4 where every
     existing lipidoid row is 0.
 (c) Emit both, in separate columns, and let the merge decide later. This is
     what `03_enumerated.csv` does now -- `IL_SMILES` is the cation and
     `parent_SMILES` is the adduct -- but a corpus row has one structure column.
EVIDENCE NEEDED: the user's call. Also worth noting: the paper reports NO
analysis quantifying methylation (no NMR integral for the N-CH3, no elemental
analysis, no titration); SI Table S1's `IN/IA` ratios quantify tail
substitution only. So "every nitrogen is methylated" rests on Figure 1A's
drawing alone, which is recorded in 08_contradictions.csv C-6.

### D-3  If merged, which quantity is `Experiment_value` -- the 0-6 ordinal or the MIC?
STATUS: OPEN
AMBIGUITY: The paper has two value layers with different coverage and different
scales. SI Table S3 gives a 0-6 point score for all 104 members; the score is an
ordinal encoding of a 3-point concentration ladder (25/50/100 uM) crossed with
"delayed" vs "no growth", defined in SI Table S2. Main-text Table 1 gives a true
MIC in uM, but only for 27 members, and the two are consistent (score 6 <-> MIC
<= 25 uM, 27/27, 07_qc.json P7).
CANDIDATES:
 (a) The 104 scores. Full coverage; but it is a discretised ordinal, so the
     corpus convention applies -- CLAUDE.md section 3: "for a discretized
     source the quantity z-scored is the BIN ORDINAL, not a midpoint", as for
     `SL_2021`. Scores 0 and 6 are both censored (0 = no inhibition at any
     tested concentration; 6 = active at the lowest tested concentration), so
     47 rows would need `value_is_bound`.
 (b) The 27 MICs. A real concentration, log-scalable, but it discards 77 of the
     104 members and is itself censored at both ends (4NC8's 100 uM is the top
     of its ladder; several 2 uM values sit at the bottom of the 1-12.5 uM
     range).
 (c) Both, as two `Experiment_ID`s, which is how a second screen in one paper is
     handled elsewhere (`LL_2026` / `LL_2026_2`). But they are not independent
     screens -- the 27 are a subset of the 104, re-measured.
EVIDENCE NEEDED: the user's call, and it only arises if D-1 resolves to merge.
The censoring counts are computable from 05_values.csv either way.

### D-4  SI Table S5's charge densities for 3NB and 4NC are transposed. Repair which half?
STATUS: OPEN
AMBIGUITY: The printed charge density for 3NB is 0.0200 and for 4NC is 0.0274;
the values computed from the structures Figure 1A draws are 0.027353 (3NB, TREN,
4 N / 146.23 Da) and 0.019967 (4NC, cyclam, 4 N / 200.32 Da). Each printed value
is the other's computed value to four decimals. Either the two table entries
were transposed, or one of the two structures is not what Figure 1A appears to
draw.
CANDIDATES:
 (a) The TABLE is wrong (two entries transposed). Supported by: 6 of 8 rows are
     exact under the same definition; the paper's own stated trend ("lower
     charge densities ... achieving the highest scores") has the wrong sign with
     the printed values (Spearman +0.259) and the right sign with the computed
     ones (-0.247); and our 4NC parent structures match `KW_2014`'s `A138_C10`
     and `A138_C12` by full InChIKey, independently confirming cyclam.
 (b) A STRUCTURE is wrong -- e.g. 4NC is cyclen (12-ane-N4, 172.27 Da) rather
     than cyclam. Rejected on the drawing (the ring alternates 3- and 2-carbon
     spacers, which cyclen does not) and on the arithmetic (4/172.27 = 0.0232,
     matching neither printed value), but recorded because this is exactly the
     `AC_2025` question and that one was got wrong on the first attempt.
EVIDENCE NEEDED: nothing from us changes either way -- no descriptor from
Table S4 or S5 is emitted. The decision matters only if the paper's Figure 1E,
or its charge-density trend, is ever cited. Flagged rather than fixed.

### D-5  Three tail groups are commercial isomer mixtures. Emit one representative, or flag and stop?
STATUS: OPEN
AMBIGUITY: `13B` is isodecyl acrylate, sold as a mixture of branched C10
isomers; Figure 1A draws one representative (8-methylnonyl) and the paper gives
no NMR, HRMS or supplier isomer statement. `9C` is the tricyclodecanyl
(hydrogenated dicyclopentadienyl) ester, whose commercial form is a mixture of
8-yl/9-yl and endo/exo isomers; the drawing fixes the cage and puts the
attachment on the two-carbon bridge next to a bridgehead (the 8-yl) but cannot
distinguish the rest. `8CB` is isobornyl acrylate, constitutionally fixed by the
drawn cage and by its 208.3 Da entry in Table S4, but with no stated
stereochemistry and commercial bornyl isomer present. Together these are 24 of
the 104 members.
CANDIDATES:
 (a) Emit the drawn representative with `structure_assumed = 1` and a note
     naming the feature. This is what `03_enumerated.csv` records now, and it
     matches how `ZT_2026` treats the same iso-decyl acrylate ("branched isomer
     mixture", `confidence = medium`) -- CLAUDE.md section 6 already lists "the
     exact branch of iso-decyl acrylate" as unverifiable from the paper alone.
 (b) Drop the 24 affected members.
EVIDENCE NEEDED: the user's call, consistent with how `ZT_2026`'s isodecyl rows
were settled. Only arises if D-1 resolves to merge.
