# DECISIONS — mckinlaycj_39105751 (Webster, Peck, … Zuckermann, Deutsch, McKinlay, *ACS Nano* 2024, PMID 39105751)

Questions this curation must NOT answer on its own. `build_dataset.py` refuses
to merge this folder while any entry below is `OPEN`.

The pipeline is good at detecting that a source is silent. It is not entitled
to resolve it. In `AC_2025` the first repair was confident, internally
consistent, and would have assigned **419 measured values to the wrong
molecules**. So: detect, queue, stop.

---

### D-1. Does a peptoid belong in a corpus of small-molecule ionizable lipids?

STATUS: OPEN

AMBIGUITY: Every other ionizable lipid in `data_vitro.csv` / `data_vivo.csv` is
a small molecule — a head group, a linker and two to four tails, typically
300–1,400 Da, and fully described by a single SMILES whose *connectivity* is the
chemistry. These are **peptoids**: N-substituted glycine oligomers of 4–9
residues, 824–1,670 Da, built on a solid support, and defined by their monomer
**order**. The paper itself makes the point that order is a separate axis —
it synthesised peptoid **266** specifically as a composition-matched sequence
isomer of **236** to test it. The two molecules share a molecular formula, an
exact mass of 1275.09 and (with 237, 259 and 264) a calculated-mass cell in SI
Table S9, and they are five different molecules.

This matters to the project's own machinery, not just to taxonomy:

  (a) **Featurisation.** §7 pins source classification and Tanimoto NN on
      Morgan/ECFP radius 2, 2048 bits. A 7-residue peptoid is a long chain of
      near-identical `-N(R)-CH2-C(=O)-` units; ECFP4 on these will be dominated
      by backbone amide environments shared by every member, and the bits that
      distinguish 236 from 266 — a monomer **transposition** — may not exist at
      radius 2 at all. If peptoids are entered and ECFP cannot separate 236 from
      266, the corpus silently contains near-duplicate features on genuinely
      different molecules.
  (b) **Result 1 and Result 2.** "Ionizable-lipid chemistry predicts its source
      of origin at near-ceiling accuracy" and "ionizable lipids are
      ~source-exclusive". A peptoid source is trivially separable from every
      other source by backbone alone. Adding one makes Result 1 *easier* in a
      way that is about molecular class, not about laboratory provenance — which
      is the distinction the whole paper rests on.
  (c) **The five descriptor columns.** `Num_carbon_in_tail`,
      `num_protonatable_nitrogens`, `num_permanent_cationic_N`,
      `formal_net_charge`, `num_unsaturated_cc_bonds` are RDKit-computed and
      will all *return* a number here. `num_protonatable_nitrogens` will count
      every backbone amide nitrogen's neighbours in a molecule whose measured
      particle pKa is 5–7 and whose predicted molecular pKa the paper puts at
      5.5–5.7 — i.e. the descriptor's usual meaning does not transfer.

CANDIDATES:
  (a) admit them as ordinary rows, one SMILES each, and accept (a)–(c);
  (b) admit them with a new `molecular_class` column (small_molecule / oligomer)
      so analyses can include or exclude them deliberately — note this is a
      **schema change**, which is a bigger decision than this paper;
  (c) keep the folder curated and unmerged, as a standalone set, the way
      `molea.csv` and `XD_2026`'s LipidDB are kept (§8);
  (d) drop the paper.

EVIDENCE NEEDED: none from the paper — the paper is unambiguous about what
these molecules are. This is a project-level scope call about what the corpus
is a corpus *of*, and it should be made once, for peptoids and any future
polymer/oligomer screen, not per paper. Worth deciding before Result 1 and
Result 2 are rerun, because the answer changes their denominator.

---

### D-2. Is a 5-row experiment admissible, and against what is it z-scored?

STATUS: OPEN

AMBIGUITY: Exactly five per-member delivery values exist in everything we hold
(main-article Table 1; see `00_value_sources.json` for the exhaustive
inventory). The corpus convention is that `Experiment_value` is a
**per-experiment z-score of log10**, population sd. Over n = 5 that is a
well-defined but very thin standardisation: mean(log10) = 8.8391, population
sd = 0.7426, and the five z-values are determined by four degrees of freedom.
Worse, the five are **not a random sample of the screen** — four of them are
the campaign's top performers and the fifth (226) is the deliberately chosen
pre-optimisation control, so the spread is a selected contrast, not an
experiment's distribution. The 34-member experiment these five were lifted out
of is in Figure 5D, which we do not hold.

A z-score over a hand-picked best-and-worst pair of groups is not comparable to
a z-score over a screen, and Result 3 and Result 4 both consume
`Experiment_value` as if it were.

CANDIDATES:
  (a) merge the 5 rows with `value_normalization_unit = experiment` as now, and
      caveat;
  (b) merge them with the raw p/s retained and `Experiment_value` left blank,
      the way `SB_2024`/`SY_2025` are flagged-not-fixed (§3);
  (c) set a minimum experiment size for admission to the corpus (the in-vivo
      pool's current smallest is `XD_2026` at 12 rows) and exclude this one;
  (d) do not merge until the Figure 5D raster is obtainable, at which point the
      experiment is 34 rows and the standardisation is honest.

EVIDENCE NEEDED: a project rule on minimum experiment size, and a ruling on
whether a *selected* subset may be z-scored as though it were a screen. Both
are general; this paper is only the first place they bite.

---

### D-3. Should SI Figure S10 (236 vs SM-102 vs MC3) be digitised?

STATUS: OPEN

AMBIGUITY: SI Figure S10 compares Nutshell 236 against **SM-102** and
**DLin-MC3-DMA** at a matched 0.75 mg/kg IV dose, 6 h, in Balb/c — a linear
axis with ticks at 1–5 × 10^11 that are present in the page's text layer, so
`curation/lib` `axis` + `bars` could calibrate it from its own furniture with no
estimation step. We hold this page.

That is exactly the shape §8 warns about. CLAUDE.md §8: the in-vivo pool's only
cross-source chemistry is MC3 and SM-102, and *"a source whose only shared
lipids are benchmarks that top every screen can push cross-lab concordance
toward +1 for reasons that have nothing to do with transferable SAR"*, with an
explicit instruction to **decide before Stage 5 is re-run whether benchmark-only
pairs belong in the pooled concordance at all**. Digitising this panel would
hand Result 3 a new benchmark-only pair. Extracting it *is* making that call, so
this curation did not extract it.

It would also be a **third** experiment: 0.75 mg/kg here, 0.125 mg/kg for the
Table 1 five, 2.5 µg for SI Figure S7B. It cannot join the emitted rows.

CANDIDATES:
  (a) leave it unextracted (what this folder does);
  (b) extract it as its own 3-row in-vivo experiment and let the §8 decision
      handle it downstream;
  (c) extract it but mark the two benchmark rows so Stage 5 can exclude
      benchmark-only pairs — depends on the §8 ruling existing first.

EVIDENCE NEEDED: the §8 ruling on benchmark-only pairs. Nothing from the paper.

---

### D-4. Source licences were not captured.

STATUS: OPEN

AMBIGUITY: §8 requires every constituent source to be CC0/CC-BY-compatible
before corpus release, and PIPELINE_V2 §3 says to capture the licence string off
the landing page at acquisition because it is far cheaper then. Neither was
captured here. The PMC record for PMC11342374 is open access and the JATS XML we
parse Table 1 out of came from PMC; the SI PDF came from `pubs.acs.org` and its
terms were not read. ACS open-access articles carry several different licences
(CC-BY, CC-BY-NC-ND, and ACS AuthorChoice), and CC-BY-NC-ND would be a problem
for redistribution of derived rows.

CANDIDATES:
  (a) read the licence line off the PMC landing page and record it, then
      proceed;
  (b) treat PMC open-access as sufficient for the five values (which come from
      the PMC XML) and the SI as reference-only;
  (c) hold the folder until the licence is confirmed.

EVIDENCE NEEDED: one line from the PMC or ACS landing page. Not fetched here —
the only automated requests this folder made to those hosts returned a bot-check
interstitial, and per CLAUDE.md's non-negotiables no attempt was made to work
around it.

---

### D-5. `Model_target` for the Table 1 values: liver, or whole animal?

STATUS: OPEN

AMBIGUITY: Table 1's column header is `luciferase expression (p/s)` and names no
organ. The sentence citing it pairs it with Figure 5D, captioned *"Total Fluc
bioluminescence in balb/c liver"*, which argues liver. But the paper's own
`Total Flux (p/s)` axis elsewhere (SI Figures S7B, S10) is whole-animal, and the
Experimental Methods describe both a whole-animal image and a separate ex-vivo
organ image without saying which Table 1 reports. `05_values.csv` carries both
readings in `organ_as_printed` and `organ_inferred` and does not choose.

CANDIDATES:
  (a) `Model_target = liver`, on the Figure 5D pairing;
  (b) `Model_target = whole_body`, on the column's own wording;
  (c) leave `Model_target` blank with the note, as F-004 blanked `N:P` on 129
      `WT_2026` rows when three candidate values were all unsupported.

EVIDENCE NEEDED: Figure 5D's y-axis label, which would say which quantity is
plotted. It is in the main-article raster we could not obtain.
