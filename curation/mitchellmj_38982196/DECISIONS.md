# DECISIONS — `mitchellmj_38982196` (Wang et al., Nat Chem 2024, PMID 38982196)

Things this paper does not settle. Per `PIPELINE_V2.md` §3 these are queued, not
resolved here: the pipeline detects that a source is silent, it does not get to
decide what the source meant. `build_dataset.py` refuses to merge a folder with
an `OPEN` entry, which is the intended state of this folder.

The reason this list is short is that one SI panel did a lot of work.
**Supplementary Fig. 1a draws the model T-MCR product `19T-O12` in full, with an
arm on each of amine 19's two primary nitrogens and none on its tertiary
N-methyl, and prints the mass spectrum that confirms it (401.129, `[M+2H]2+`;
our two-arm structure gives 401.17, the one-arm alternative 230.4 or 459.8).**
Before that panel was found, "how many arms does a polyamine get" was the whole
curation's blocker. What remains open is only what that one drawing cannot
reach.

---

### D-1
**STATUS: RESOLVED FROM THE SOURCE — recorded here because it was the central question and a reader should see how it was settled, not just the answer.**

**AMBIGUITY:** How many amidine-thioether-ester arms does an amine with more
than one primary nitrogen receive? Fig. 1a's scheme shows `R1–NH2` producing one
product and the naming convention `<amine>T-O<tail>` carries no stoichiometric
suffix, so the figure that defines the library is silent.

**RESOLVED BY:** Supplementary Fig. 1a (SI p.3). It draws `19T-O12` as the
bis-adduct and prints `19 (1 equiv.)  Traut's reagent (3 equiv.)  O12 (2.5
equiv.)` beneath the reagents, plus a mass spectrum peak labelled
`401.129  100.00%  [M+2H]2+`. Reproduced as check **P1b**: our two-arm structure
is 800.32 Da average, `[M+2H]2+` = 401.17 (0.04 from the printed value, against
an instrument that reports 603.73 for a 603.01 compound elsewhere in the same
SI); the one-arm alternative is 458.8 Da, i.e. 459.8 at 1+ or 230.4 at 2+, more
than 50 Da away in either charge state.

**EMITTED:** one arm per primary nitrogen, for every amine.

---

### D-2
**STATUS: OPEN**

**AMBIGUITY:** Amines **18** and **24** each carry a SECONDARY N–H in addition to
their two primary nitrogens. 2-Iminothiolane is used as a primary-amine reagent,
but nothing in this paper states that its secondary amines were left alone, and
the only product the paper draws (`19T-O12`, from amine 19) has no secondary
N–H, so it cannot speak to this.

**CANDIDATES:**
1. Two arms, on the two primary nitrogens only; the secondary N–H untouched.
   This is what 03_enumerate.py emits, with `structure_assumed = 1` and
   `positional_ambiguity = 1` on all 8 rows, and the reason written into
   `curation_note`.
2. Three arms, the secondary N–H also acylated. With 3 equiv Traut's reagent
   this is not excluded by stoichiometry.
3. A mixture of (1) and (2), which is what an unpurified one-pot product would
   actually be.

**EVIDENCE NEEDED:** the main article's Methods (equivalents and time for the
library reactions — not obtained, see `00_value_sources.json`), or an MS of any
`18T-*` or `24T-*` member. Supplementary Fig. 5 purified five lipids and says
"target products were confirmed by mass spectrometry" but prints no m/z, and
none of the five is `18T-*` or `24T-*`.

**AFFECTS:** 8 of 100 rows (amines 18 and 24 × 4 tails).

---

### D-3
**STATUS: OPEN**

**AMBIGUITY:** Amines **20** and **22** have three primary nitrogens and amine
**25** has four. The only stoichiometry the paper prints anywhere is Supplementary
Fig. 1a's model reaction: 1 equiv amine, **3 equiv** Traut's reagent, **2.5
equiv** acrylate. If the library was run at those equivalents for every amine,
then for amine 25 that is 0.75 equiv of Traut's and 0.63 equiv of acrylate *per
primary nitrogen* — below what exhaustive substitution needs — and for amines 20
and 22 it is 1.0 and 0.83.

**CANDIDATES:**
1. Exhaustive substitution of every primary nitrogen, because that is the
   declared design (`25 x 4 = 100`, one compound per amine × tail pair, one name
   per cell). This is what 03_enumerate.py emits, at `confidence = low` with
   `structure_assumed = 1` and the reason in `curation_note`.
2. A distribution of partially substituted species whose mean substitution is
   set by the equivalents, with the stored name denoting the intended product.
3. The equivalents were scaled per primary nitrogen for the library and the
   model reaction's numbers apply only to amine 19.

**EVIDENCE NEEDED:** the main article's Methods, or an MS of any `20T-*`,
`22T-*` or `25T-*` member. Note the direction of the risk: under candidate 2 the
STORED MOLECULAR WEIGHT of a `25T-*` lipid would be up to ~1.0 kDa too high,
which no check in this folder can see, because every check here is either a
count, an index or a mass of a *different* compound.

**AFFECTS:** 12 of 100 rows (amines 20, 22, 25 × 4 tails).

**DO NOT resolve by picking the chemically tidier option.** `AC_2025`'s first
repair was confident, internally consistent, and put 419 measured values on the
wrong molecules.

---

### D-4
**STATUS: OPEN**

**AMBIGUITY:** Every value in `05_values.csv` was measured on an **unpurified**
one-pot reaction mixture. Supplementary Fig. 5 says so in as many words
("Selected AID-lipids were purified… and target products were confirmed by mass
spectrometry"), and the screen it compares against is the unpurified material.
The corpus has no column that says "this value belongs to a crude mixture whose
major component is the stored SMILES".

**CANDIDATES:**
1. Emit as-is. The paper's own purified-vs-unpurified comparison (5 lipids, i.m.
   in vivo, n=2) spans 0.81×–1.23× and does not separate them, so the crude
   material is a reasonable proxy for the pure compound *in this assay*.
2. Set `structure_assumed = 1` on all 100 rows with a note saying the measured
   species is a crude mixture. This is literally true but would make the flag
   uninformative — it would no longer distinguish the 20 rows where something
   specific about the connectivity is assumed.
3. Add a corpus-level property rather than a per-row flag, since it is a property
   of the whole experiment.

**EVIDENCE NEEDED:** a project decision on where "the screened material was
crude" belongs. This is not unique to this paper — `NG_2025`'s Mannich products
are the same situation and were handled with `positional_ambiguity`, which does
not fit here because nothing is positionally ambiguous; the question is purity,
not regiochemistry.

**EMITTED MEANWHILE:** candidate 1, with the fact recorded in
`08_contradictions.csv` rather than on the row.

---

### D-5
**STATUS: OPEN**

**AMBIGUITY:** The licence. `CLAUDE.md` §8 requires every constituent source to
be CC0/CC-BY-compatible before the corpus is released. Europe PMC reports
`isOpenAccess = "N"` and `license = null` for PMID 38982196; the PMC copy
(PMC11446653) is an NIH author manuscript. The licence line on the Nature
landing page was not read because `www.nature.com` redirects to `idp.nature.com`
without a session.

**CANDIDATES:** (a) the Source Data files are separately licensed, as some
Springer Nature supplementary material is; (b) the article is subscription-only
and the data cannot be redistributed; (c) an NIH public-access licence applies to
the author manuscript but not to the publisher's Source Data.

**EVIDENCE NEEDED:** the licence string from the article landing page, or
Springer Nature's terms for `static-content.springer.com` ESM. Every value in
this curation came from those ESM files.

---

### D-6
**STATUS: OPEN**

**AMBIGUITY:** Which experiment is the corpus row? The screen
(Fig. 2b/2c, 100 lipids, n=2) and the follow-up (Fig. 2g, 8 lipids, n=2, same
cell line, same stated 15 ng dose) disagree: geometric mean ratio 0.64×, and the
eight shared lipids rank at **Spearman −0.262** between them while the screen's
own replicates agree at **0.990**. See `08_contradictions.csv`, row 1.

**CANDIDATES:**
1. Emit the screen only, as one experiment of 100 rows. This is what
   `05_values.csv` does, because it is the only panel that covers the library.
2. Emit both as two separate `Experiment_ID`s, which would give the corpus a
   rare same-lab repeat pair — directly relevant to Result 3 and to the
   assay-noise ceiling in T2.
3. Emit the screen and record the follow-up as a reproducibility datum outside
   the corpus.

**EVIDENCE NEEDED:** a project decision. Note that candidate 2 is the only one
that puts the −0.262 into the analysis rather than into a note, and that an
8-lipid pair is small — the same objection T8 raises about the in-vivo pool.
