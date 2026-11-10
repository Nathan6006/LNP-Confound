# DECISIONS — `mitchellmj_41555650` (Han et al., ACS Nano 2026, PMID 41555650)

Queued, never resolved alone. `build_dataset.py` refuses a folder with an `OPEN`
entry, which is the point: none of these should be settled by the pipeline.

---

### D-1. The headline screen is a POOLED mRNA-BARCODE experiment. Is any of it corpus-admissible?

**What the readout is, in the paper's own words.** The abstract: *"We utilize a
high-throughput in vivo screening approach based on mRNA barcoding to study a
large library of LNPs."* Table S1's caption: a greyed LNP *"was not included in
the **pooled administration** for high throughput in vivo screening due to not
meeting the selection criteria (< 200 nm in size, < 0.3 PDI, > 50% encapsulation
efficiency, and > 15 ng/uL mRNA)."* So library members are formulated with
distinct barcoded mRNAs, mixed, injected together, and delivery is read as
barcode abundance by deep sequencing.

**Pool size, measured (check P14).** 98 b-LNP numbers were assigned; 8 of them
(19, 41, 55, 56, 57, 67, 93, 95) never appear in Tables S1/S2 at all; 93 rows are
tabulated (90 library + 3 benchmarks); **12 are greyed** (1, 2, 7, 20, 23, 33,
36, 51, 60, 63, 73, 74). **78 b-LNPs went into the pooled administration.**
Greying is a colour property invisible to text extraction, so this looked "not
determinable" until the character colours were read; `non_stroking_color` has
exactly two values on those pages, black and 191/255 grey. It is independently
confirmed: applying the caption's four numeric thresholds to Table S2's printed
numbers reproduces the grey/black split on **90 of 90 rows, 0 disagreements**.

**MY FIRST VERDICT HERE WAS WRONG AND IS CORRECTED.** I wrote that a pooled,
competitive barcode value is not commensurate with `data/`. Measured against the
corpus, that is false: `data/data_vivo.csv` carries an **`Experiment_batching`**
column whose vocabulary is `{barcoded, individual}`, and **`LR_2023` already
contributes 407 barcoded in-vivo rows** (383 unique ionizable lipids, liver,
`protein_adundance_normalized`). So **pooling alone is not disqualifying** — the
corpus's convention is to accept barcoded data and *declare it on the row*. The
competitive-measurement concern is real and is what the column exists to record,
not a reason to exclude.

**So the binding constraint is ACQUISITION, not design.** The screen's values are
in main-article Figure 2, which was not acquired (D-3). Nothing about the pooled
design rules it out, and this folder should be re-opened when the article is in
hand.

**What remains genuinely unsettled** is the normalisation. A sibling curation
reconstructed a barcoded screen's normalisation from the paper's own worked
example as `nd_i = [(count_i/Σcount)/(input_i/Σinput)] / Σ_k[same] × 100%`,
reproducing 6/6 column totals and 18/18 cells at residual 0.0, which established
that each value is a **share of its own pool**. **That reconstruction could not
be attempted here**: the worked example, the normalisation statement and the
values are all in the unacquired main article. The residual is therefore
**not computed**, and is recorded as not done rather than quietly skipped (D-2).

**Emitted today:** 0 pooled values (none acquired), and 48 values from the
*separate single-LNP validation*, which carry `Experiment_batching = individual`.

STATUS: OPEN

---

### D-2. The pooled screen's normalisation residual was NOT computed.

For a pooled barcode screen the cheap analogue of "a round number the authors
drew but never printed" is that normalised delivery across one pool must sum to
the stated normalisation (1, or 100 %). **That test was not run.** The pool's
values, the normalisation formula and any worked example are all in the
unacquired main article (D-3). Recorded as not done, with the reason, rather
than omitted.

What *was* available and was tested instead: the **composition** of the pool
(78 members, P14, 90/90 agreement with the stated criteria) and the **per-animal
percentage identity** of the single-LNP validation, where SI Fig S1A's organ
percentages must be `100 × S2A_organ / Σ S2A_organs` — checked at median 6.2 %
over 46 values (P11).

STATUS: OPEN

---

### D-3. The main article was never acquired, and that is the binding constraint.

ACS (`pubs.acs.org`) returns a Cloudflare interstitial to a plain client;
`provenance/tools/fetch_paper.py` read the title `"Just a moment..."` and
enumerated 0 assets. `BROWSER_ACQUISITION.md` records that ACS needs a browser
and no login, but this run was instructed not to use the shared Playwright
browser, and a further curl probe of the ACS supplementary URLs was refused by
the sandbox as bot-check circumvention. **No bypass was attempted.**

Consequently: main-article Figures 1–5 are unseen. That removes the barcode
screen, Figure 3 (the single-LNP validation as the authors present it), Figure 4
(flow cytometry) and any structure the paper draws. The 48 emitted values come
from the SI's own Figure S2A, which carries the same validation experiment.

**What a browser run would add:** Figure 3's panels (which Table S3's CVs would
then validate dimensionlessly), the pooled screen, and a drawn structure for the
C10 anchor. This folder is worth re-opening once the article is in hand.

STATUS: OPEN

---

### D-4. `Experiment_ID` is proposed, not reserved.

`EH_2026`, from first author **Emily L. Han**. Not checked for collision against
a reservation list, and not merged. Note the lab is **Mitchell**, not Dahlman —
the brief that commissioned this extraction misattributed it, and the folder was
renamed from `dahlmanje_41555650` to `mitchellmj_41555650`.

STATUS: OPEN

---

### D-5. Table S1's four molar figures do not sum to 100 on 75 of 93 rows.

Measured (check P7): the sums take six distinct values — 100, 102.5, 110, 112.5,
120, 122.5. Only **18 of 93** sum to 100: the fifteen b-LNPs at the base
(35 % IL, 2.5 % PEG) setting, plus all three benchmarks (MC3 50/10/38.5/1.5,
ALC-0315 46.3/9.4/42.7/1.6, SM-102 50/10/38.5/1.5 — each exactly 100).

The pattern is unambiguous: the authors varied IL (35/45/55) and PEG (2.5/5)
while leaving DOPE at 16 and cholesterol at 46.5 **unchanged**, so the printed
quartet is a set of *parts* rather than mole percentages for every row except
the base one. Two readings follow and they are not reconcilable from what is
printed:

* **(a)** the true mol % are the printed figures renormalised by their sum
  (b-LNP 14 → 40.9 / 14.5 / 42.3 / 2.3), or
* **(b)** the printed figures are the intended mol % and DOPE/cholesterol were
  reduced to accommodate them, with the table not updated — in which case the
  true helper and sterol percentages are unrecoverable.

**Nothing was guessed.** The emitted rows store the printed figures verbatim in
`IL_molratio` / `HL_molratio` / `CHL_molratio` / `PEG_molratio`, carry
`molratio_sum`, and every row whose sum is not 100 says so in its
`curation_note`. 30 of the 48 emitted rows are affected (b-LNPs 14, 70, 96, 97,
98). This is the same shape as F-004's `N:P`: several mutually inconsistent
candidates, so the honest act is to record the printed value and the
inconsistency, not to pick.

STATUS: OPEN

---

### D-6. Licence is not release-compatible.

ACS Nano, `isOpenAccess = N`, Europe PMC `license = null`, the only full-text URL
marked *"Subscription required"*. There is no author repository (0 hits). So this
source is **not CC0/CC-BY**, and `CLAUDE.md` §8 requires constituent-source
licence compatibility before corpus release. Structures derived from a
systematic name and values read off a figure are facts rather than expression,
but that is a judgement for the release decision, not for a curation agent.

STATUS: OPEN

---

### D-7. `data/` carries a conflicting structure for `C12-200`, and we do not edit `data/`.

Measured in `02_corpus_overlap.json` and `08_contradictions.csv`: `KL_2010`'s
`200C8/200C10/200C12/200C14` are **C20–C26, three nitrogens, ONE
2-hydroxyalkyl arm and two N-propyl groups**, whereas this paper's `C12-200`
(and `AC_2025`'s, on 498 rows, and `JW_2024`'s) is **C70H145N5O5, five
nitrogens, five arms**. `KL_2010` is Love et al. 2010 — the original epoxide
lipidoid library this amine numbering comes from — and its `114C12`/`114C14`
match our `C12-114`/`C14-114` on the **full InChIKey**, so the experiment is not
wholesale wrong; its 200-series is.

This is a finding ABOUT `data/`, raised by this curation and **not acted on**.
Deciding which half to repair is exactly the `AC_2025` situation where the first
confident repair would have put 419 measurements on the wrong molecules.

STATUS: OPEN

---

### D-8. Three `Model_target` values would be new to the corpus.

`data/data_vivo.csv`'s `Model_target` vocabulary is `bone_marrow, brain, ear,
liver, lung_epithelium, muscle, spleen, whole_body`. The emitted rows use
`brain`, `liver` and `spleen`, which exist, and **`heart`, `kidney` and `lung`,
which do not**. `lung` is deliberately not mapped onto the existing
`lung_epithelium`: that value names an epithelial compartment, and SI Fig S2A
measures whole-organ flux from an excised lung. Every affected row carries
`Model_target_is_new_in_corpus = yes` (24 of 48). Adding a `Model_target` value
is a corpus-vocabulary change, not a curation call.

STATUS: OPEN
