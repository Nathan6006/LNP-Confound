# DECISIONS — XZ_2020 (Zhao et al. 2020, PMID 32662132)

This folder is **not merge-ready, deliberately.** `build_dataset.py` refuses a
folder with an `OPEN` entry, and every entry below is open. Each is a place
where the source is silent or self-contradictory and the pipeline does not get
to decide what it meant.

---

### D-1 Does a single-lipid formulation DoE belong in the corpus at all?

`06_rows.csv` holds 13 rows and **one** unique ionizable structure, 93-O17S, at
twelve excipient ratios plus the lipidoid-only formulation. There is precedent
(`ZT_2026` is a 486-formulation DoE on one lipid), but the consequence is sharp:
the per-experiment z-score of a one-lipid block encodes formulation, not
chemistry, so the block can contribute nothing to shared-lipid concordance
(Result 3) and adds one more near-disjoint source to a corpus whose problem is
near-disjoint sources. Merging it would add an `Experiment_ID` whose entire
chemical content is one molecule.

STATUS: OPEN

### D-2 Should the four Figure S3 %GFP values for 93-O17S be emitted too?

They are the **only printed, exactly transcribed** delivery values in the whole
document (2.68 / 4.16 / 3.19 / 7.05 % GFP+ at 5 / 15 / 30 / 60 ug/mL) and they
cost no geometry at all. They are also a dose response of one lipid, so four
rows would share one `(Experiment_ID, IL_name)` key. Options: emit all four as a
second experiment `XZ_2020_2`; emit only the 30 ug/mL point as the
representative; or leave them in `05_values.csv`, which is what this run did.

STATUS: OPEN

### D-3 Is `O16B` the same tail as `O16S-S`, and `N16B` the same as `N16S-S`?

`O16B` and `N16B` appear only as Figure S1 panel titles; `O16S-S` and `N16S-S`
appear only in the text and in Figures S4/S5. The chain-atom counts match and
"B" plausibly abbreviates bisulfide.

**Re-attempted 2026-11-10 with head 9322 in hand and still not settled.** The
predictions now exist — 9322-O16S-S `[M+H]+` **804.487**, 9322-N16S-S
**802.519** — but both panels **refused** the tick-label calibration (O16B
spacing nonuniformity 18.3 %, N16B never reaches 10 label boxes, against a 4 %
criterion that was **not** loosened). And even a clean calibration could not
decide it: the two differ by **1.97 Da = 0.47 px** on a 4.2 Da/px axis. Neither
is built.

STATUS: OPEN

### D-4 What does the "C" in `O17C` designate?

The name-token rule (5/5 on the tails the SI synthesises) fixes the chain at 17
heavy atoms, and nothing fixes the linker. A plain C17 alkyl chain, a carbonate
and a carbamate all satisfy the name. There is no 13C spectrum anywhere in the
SI — searched for by name, not assumed absent.

**New evidence 2026-11-10, enough to exclude one candidate and not enough to
choose.** With head 9322 known, the 9322-O17C panel calibrates and gives a
**derived per-arm mass of 310.19 ± 1.99 Da**. That **excludes a carbonate**
(the C22H40O5 variant is 74.1 Da away, ~19 px) and is consistent with
heptadecyl acrylate (C20H38O2, 310.287 — **0.1 Da** away). But ±2 Da cannot
separate that from other 310–312 Da candidates, so the linker is still not
fixed. `O17C` is not built.

STATUS: OPEN

### D-5 Should this paper be curated at all without the main text?

Everything the paper is *about* — the ex-vivo screen of 18 amine heads against
up to 8 tails, and the in-vivo lead — is plotted only in Figures 1-4, which are
behind the Wiley subscription wall with no deposited author manuscript and no
open copy (`00_acquisition.json`). What is reachable is 2 lipidoids, 1 of them
structureless. The alternative to merging 13 degenerate rows is to park the
folder until the article can be obtained through a licensed route and then
curate the screen properly.

STATUS: OPEN

### D-6 Is a catalogue number from *other* papers an acceptable structure source?

Head 93's identity comes entirely from outside this paper: `AA_2008` and
`KW_2014` corpus rows naming the exact 3-aminopropylimidazole structure
`93_O12 ... 93_O18`, `93_N12 ... 93_N18`, `A93_C12`, and `curation/xuq_32520527`,
which read head 93 off a drawn figure in a third Xu-lab paper. Three independent
records agree. But the inference "this paper's 93 is that 93" is a convention of
the Xu/Anderson amine catalogue, not a statement this paper makes — and the run
brief's own warning is that a leading number token is not an identity (22 corpus
rows starting `9322` are an unrelated index). Every emitted row carries
`confidence = medium` and a note saying exactly this.

STATUS: OPEN

### D-7 The four 9322-O17S measurements — now with a molecule.

Figure S3 gives 9322-O17S at 8.13 / 6.83 / 11.1 / 0.74 % GFP+ — printed and
exact. **This entry previously said the head was unobtainable; that was my
error** (`08_retracted_suspicions.csv` R-06). Head 9322 is
1-(3-aminopropyl)-2-methylimidazole, attested by 22 `YR_2024` corpus rows under
the guarded lookup and confirmed against this paper's own Figure S1 at 0.20 px.

So the question is no longer "structureless row or discard" but the same one as
D-2: these four values are a **dose response of one lipid**, and emitting them
would put four rows on one `(Experiment_ID, IL_name)` key. If D-2 is resolved in
favour of emitting, this resolves with it and adds **8 rows** on **2** distinct
structures — which would be the only thing in this folder that puts more than
one chemistry in one block.

STATUS: OPEN

### D-8 Weight ratios into a mole-percent schema.

Figure S2a's categories are **w/w** ratios (lipidoid : cholesterol : DOPE :
DSPE-PEG). `HL_molratio`, `CHL_molratio` and `PEG_molratio` are mole percents.
Cholesterol and DOPE could be converted; DSPE-PEG cannot, because its PEG block
length is a distribution and no single molar mass is right — which is why
`reagent_registry.py` deliberately holds no PEG-lipid. This run left all three
blank and put the w/w ratio in `curation_note`. A partial conversion would make
two of four components comparable to the rest of the corpus and the third
silently absent, which is worse than blank.

STATUS: OPEN

### D-9 Licence not captured.

`CLAUDE.md` section 8 requires every constituent source to be CC0/CC-BY-
compatible before corpus release. The Wiley landing page could not be read
without the browser this run is forbidden to use, so the SI's licence string is
unrecorded.

STATUS: OPEN
