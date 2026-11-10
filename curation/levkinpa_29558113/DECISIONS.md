# DECISIONS — PL_2018 (Molla & Levkin, Bioconjug Chem 2018, PMID 29558113)

Queued, not resolved. `build_dataset.py` refuses a folder with an `OPEN` entry,
which is the point: every one of these would change what enters `data/`, and
PIPELINE_V2.md §3 is explicit that the pipeline detects that a source is silent
and does not get to decide what the source meant.

---

### D-1
STATUS: OPEN

AMBIGUITY: The library's 18 amine head groups A1–A18 are enumerated in the MAIN
ARTICLE, which this curation does not hold. The SI names or draws only A1, A5
and A10 and gives a molecular formula for A8. Table S2's 45 valued lipidoids use
13 distinct amines (A1, A3, A4, A5, A6, A7, A8, A9, A11, A13, A14, A15, A17), so
**40 of the 45 rows have a measured value and no structure**. Europe PMC,
2026-11-10: `isOpenAccess=N`, `inEPMC=N`, `pmcid=null`, "Subscription required";
no author repository (`find_repo.py`, 0 hits, though two GitHub *code* queries
were rate-limited and are recorded as incomplete rather than negative).

CANDIDATES:
 (a) Acquire the main article (ACS, browser, no login per
     `provenance/BROWSER_ACQUISITION.md`), read the amine panel, re-run stages
     02–06, and merge 45 rows.
 (b) Merge only the 5 structurally determined rows (L1, L11, L31, L35, L44) as
     `PL_2018`.
 (c) Do not merge. A 5-lipid experiment spanning 2 head groups adds one more
     near-disjoint silo and cannot contribute to Result 3 (shared-lipid
     concordance) at all.

EVIDENCE NEEDED: The main article's amine figure. Nothing else resolves this —
no amount of further work on the SI will, because the information is not in it.

---

### D-2
STATUS: OPEN

AMBIGUITY: A8 is determined only as the molecular formula **C7H17N3**, solved
from two MALDI formulae (T18U-PY12-A8 = C41H81N4O2S2; T14-PY12-A8 =
C37H75N4O2S2; their difference C4H6 is exactly T18U − T14, which is the
consistency cross-check). A8 carries the paper's **second and third best**
lipidoids (L2 at 89 %TE, L3 at 87 %TE), so this is not a marginal gap.

CANDIDATES: 1-(2-aminoethyl)-4-methylpiperazine; N-(3-aminopropyl)-N,N'-
dimethylethylenediamine; 3-(4-methylpiperazin-1-yl)propylamine minus a CH2;
and other C7H17N3 triamines. A formula is not a structure and nothing in the SI
chooses between them.

EVIDENCE NEEDED: the main article's amine panel (same as D-1), or the 1H NMR of
T14-PY12-A8 in Figure S12 read quantitatively — its assignment letters are
printed on the spectrum and might fix the head group's topology, but that is a
separate piece of work and was not attempted here.

---

### D-3
STATUS: OPEN

AMBIGUITY: A10 is drawn (Figure S1 top scheme and the panel d inset) but never
named. The only geometric measurement the raster supports — the drawn polygon's
width-to-vertical-side ratio — is **inconclusive** once its tolerance is swept:
it returns hexagon at some tolerances and heptagon at others (1.618 / 1.732 /
1.802 are the 5-, 6- and 7-gon predictions; `01_amines.json`). 1-ethylpiperazine
is the obvious chemical candidate and is recorded as a *candidate*, not as a
structure.

CANDIDATES: (a) accept 1-ethylpiperazine with `confidence = low`; (b) leave A10
unenumerated, as now.

EVIDENCE NEEDED: the main article's amine panel, or a higher-resolution copy of
Figure S1 than the 834×1189 JPEG the SI embeds. **Nothing currently depends on
it: A10 appears in no row of Table S2.**

---

### D-4
STATUS: OPEN

AMBIGUITY: How to put these values on the corpus's scale. `Experiment_value` is
a per-experiment z-score of log10 (CLAUDE.md §3), but Table S2 reports only the
**45 of 288** lipidoids at or above 40 %TE. Every moment computed from what we
hold is a moment of a **truncated top tail**, not of the experiment. The mean
and sd of the 243 unreported lipidoids are unknown and unknowable from this SI.
`06_rows.csv` therefore leaves `Experiment_value` blank and
`value_normalization_unit` empty ("not stated"), and stores the raw integer
percent in `value_raw_percent_TE`.

CANDIDATES:
 (a) z-score over the 45 reported rows and flag every row with a
     `curation_note` saying the block is a truncated top tail;
 (b) store the raw percentage and let the corpus carry a non-z-scored block, as
     `SB_2024`/`SY_2025` already do in vivo (CLAUDE.md §8, itself still open);
 (c) treat the 243 unreported lipidoids as left-censored at 40 %TE and z-score
     the full 288 under an explicit censoring model.

EVIDENCE NEEDED: a project-level call, not more data. Note that (a) and (b) give
different answers for the same rows and the difference is not small: the
reported range is 40–95 on a scale whose true floor is 0.

---

### D-5
STATUS: OPEN

AMBIGUITY: Three further transfection panels exist in this SI and none is on the
same footing as Table S2 — Figure S6 (HEK293T but a **10% eGFP / 90% LacZ**
payload, relative to SFA), Figure S8 (**HepG2**, % transfection efficiency,
top eleven lipidoids), and Figure S3 (the same experiment as Table S2 but
expressed **relative to SFA**). None was digitised (`00_value_sources.json`
VS-3, VS-4), because with 40 of 45 head groups unknown they would add unusable
rows.

CANDIDATES: (a) digitise Figure S8 to give the paper a second `Model_type`;
(b) digitise Figure S6 as a separate experiment with its own payload;
(c) leave both, as now.

EVIDENCE NEEDED: whether D-1 resolves in favour of acquiring the main article.
If it does not, this is moot.

---

### D-6
STATUS: OPEN

AMBIGUITY: Licence. CLAUDE.md §8 requires every constituent source to be
CC0/CC-BY-compatible before corpus release. This is an ACS subscription article
(Europe PMC `isOpenAccess=N`) and its SI carries ACS copyright; the PDF held
here was downloaded through an institutional session (the per-page watermark
names the institution). No CC licence string was captured from the landing page
because this curation did not visit it.

CANDIDATES: (a) treat extracted numeric values and structures as facts, not a
derivative work, and merge; (b) exclude PL_2018 from any released corpus while
keeping it for internal analysis; (c) seek permission.

EVIDENCE NEEDED: the landing page's rights statement, and a project-level
position on (a) that covers the other ACS-sourced experiments already in the
corpus.

---

### D-7
STATUS: OPEN

AMBIGUITY: `T16-PY18-A13` appears **three times** in Table S2 — L17 at 57, L18
at 57, L26 at 50. 45 rows carry 43 distinct names while the caption says 45
lipidoids. Either two rows duplicate one measurement, or one or two of the three
names are mistyped.

CANDIDATES: (a) keep all three rows as printed (current state); (b) drop L18 as
a duplicated row; (c) treat one or more names as typos and reassign.

EVIDENCE NEEDED: the main article's Figure of the full screen, or the authors.
Not repaired here: reassigning a name on a guess is exactly how `AC_2025`'s
first repair would have moved 419 measured values onto the wrong molecules.
None of the three is emittable today (A13 is undetermined), so nothing in
`data/` turns on it yet.

---

### D-8
STATUS: OPEN

AMBIGUITY: `PY12`'s name and its own data disagree. SI p.5 calls it
"2-(**decyl**disulfanyl) pyridine" (C10); its printed HRMS (C17H29NS2, m/z
312.1821) and the twelve stroked bonds Scheme S1b draws both say **dodecyl**
(C12). We emitted dodecyl (08_contradictions.csv C-2).

CANDIDATES: (a) dodecyl — the mass, the drawing and the label `PY12` all agree,
and a decyl compound's [M+H]+ would be 284.151, 28 Da away; (b) decyl, treating
the name as authoritative and the mass and drawing as both wrong.

EVIDENCE NEEDED: the main article's experimental section, or the reagent's
catalogue number. **One emittable row depends on this** — L11, `T16-PY12-A1`,
68 %TE — so it is not cosmetic. The AC_2025 lesson is that the confident repair
can be the wrong half; here three independent lines point one way and one name
points the other, which is why (a) was taken, but the call is recorded rather
than assumed.

---

### D-9
STATUS: OPEN

AMBIGUITY: Whether `A5`'s two N-ethyl branches are good enough to emit. A5's
three-carbon linker IS measured (7 of 8 lag settings agree; a two-carbon linker
would make A5 identical to A1), but its two N-substituents could not be measured
— the Figure S1 raster's red ink drops out entirely at the far end of one bond,
and Douglas-Peucker returns 4 vertices for a one-corner path at every eps from
1.2 to 5.0. They are READ as ethyls. The two affected rows (L31, L44) are
emitted with `confidence = medium`, `structure_assumed = 1` and a
`curation_note` naming the feature.

CANDIDATES: (a) emit as now; (b) withhold L31 and L44 until the main article
names A5, leaving 3 emittable rows instead of 5.

EVIDENCE NEEDED: the main article's amine panel (D-1), or a higher-resolution
Figure S1.
