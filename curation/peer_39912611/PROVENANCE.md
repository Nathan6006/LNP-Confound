# PROVENANCE -- peer_39912611

Naidu GS, Rampado R, Sharma P, Ezra A, Kundoor GR, Breier D, Peer D.
*Ionizable Lipids with Optimized Linkers Enable Lung-Specific, Lipid
Nanoparticle-Mediated mRNA Delivery for Treatment of Metastatic Lung Cancer.*
ACS Nano 2025. PMID **39912611**, PMCID **PMC11841047**,
DOI **10.1021/acsnano.4c18636**. Licence **CC BY 4.0**, taken from the
article's own `<permissions>` block, not assumed.

**Nothing in this folder has been merged into `data/`.** No file outside
`curation/peer_39912611/` was created or modified.

---

## 1. Where the files came from, and which routes failed

`00_acquisition.json` carries the sha256, byte count, URL and licence of all
nine inputs, plus the six acquisition routes that did NOT work and why. The two
that matter:

* **PMC's `/articles/<id>/bin/<name>` route served an HTML error page for all
  ten declared assets**, reproducing what sibling curations in this run
  reported.
* **PMC's article-PDF route is behind a proof-of-work bot check**
  (`POW_CHALLENGE`, `POW_DIFFICULTY=4`, cookie `cloudpmc-viewer-pow`). It was
  not bypassed, so the fallback the task suggested -- the article PDF's own
  embedded image streams -- is unavailable.

The figures were obtained instead from **Europe PMC's
`europepmc.org/api/fulltextRepo`**, which serves them to a plain client with no
session. The parameter `mimeType` is mandatory: without it the same URL returns
HTTP 500 with a 55-byte JSON body.

`pubs.acs.org`'s larger rasters were **not attempted**:
`provenance/BROWSER_ACQUISITION.md` measures ACS at 403 to a plain client and
200 only from a browser on a second navigation, and this run is restricted to
plain clients. **700 px is therefore the resolution limit for Figures 2 and 3,
and it is a constraint, not an oversight.** Every value in `05_values.csv`
carries it: one pixel is 0.020 decades in Figure 2A, 0.040 in 2B/2C, 0.021 in
3C, 0.019 in 3D and 0.014 in 3E.

## 2. What the paper is, in corpus terms

**12 distinct ionizable lipid structures** -- lipid 14 and lipids 32-42 -- in
**one** LNP composition each (40 : 48.5 : 10 : 1.5 ionizable : cholesterol :
DSPC : PEG-DMG2000, N/P 9). So 12 lipids = 12 formulations = 12 LNPs; the
"N lipids means N formulations" trap does not bite here. What multiplies is the
*values*: 3 cell lines x 3 doses in vitro (108) and 3 organs in vivo (36).

**Corpus overlap: none** (machine-readable in `02_corpus_overlap.json`). All 30,998 rows of `data/data_vitro.csv` and
`data/data_vivo.csv` were parsed and compared on FULL InChIKey, on the InChIKey
skeleton block, and on `IL_name`. Zero matches on each. A substructure search
for this paper's bis(8-((2-hexyldecanoyl)oxy)octyl)amine core found **3 rows in
`ZT_2026`** -- `Lipid 14-CB1`, `-CB2`, `-CB3` in `data/data_vivo.csv`, which
are quaternised cholesterol conjugates built on this paper's lipid 14 by a
different group. They are different molecules, so they are not an overlap; they
are an independent corroboration of lipid 14's scaffold, and they are what
caught our own first misreading of it (see `08_contradictions.csv` row 1).

## 3. Anchor posture -- verified BY HAND

`provenance/structure_check.py`, `check_anchor_independence.py` and
`check_anchor_coverage.py` were **not run**: they write tracked files in
`provenance/` that other agents are using concurrently. The posture is
therefore argued here and the evidence is in the folder.

* `reference_structures.csv` holds **12 anchors plus one explicitly empty
  row**, with the exact header the task specifies.
* `make_reference.py` builds them. It imports **no** `curation/lib` module,
  constructs **no** `RWMol`, calls **no** `esterify` or `attach`, and never
  rewrites a canonical SMILES. It emits characters, walking each lipid's
  PRINTED IUPAC NAME from the core nitrogen outward, one atom token per
  locant, and it builds the tail arm by counting out "octane-8,1-diyl" and
  "2-hexyldecanoate" rather than reading `stage01`'s `ARM` constant. The two
  routes share no function and no constant.
* **Byte identity is not the signal, and testing for it is the F-037 trap.**
  All 12 anchors are byte-identical to the enumeration once canonicalised,
  which is what two correct routes SHOULD produce. What was verified instead is
  that the route has teeth: the token list for lipid 34 was changed from
  `C,C,N,C(=O),...` (amide) to `C,C,C(=O),N,...` (reverse amide) in a copy
  under the scratchpad, and **the agreement count moved 12 -> 11 and the exit
  code went 0 -> 1**. The seeded anchor's InChIKey came out equal to lipid
  36's, which is precisely the isomer confusion a mass check cannot see.
* The one real gap is recorded as an empty row: **lipid 14 has no synthesis
  write-up, no IUPAC name, no NMR and no ESI-MS in this SI** -- the authors
  cite their own earlier work -- so its anchor comes from the drawn scheme
  (Supplementary Figure 15) rather than from a name.

## 4. Why the linker directions are trustworthy

The point of this paper is the linker, and its linker pairs are **isomeric**:
lipid 14/32 are both C56H110N2O6, lipid 34/36 both C56H111N3O5 (the SI prints
[M+1]+ 907.2 and 907.3), lipid 37/38 both C58H113N3O6 (949.2 and 949.2), lipid
40/42 both C58H114N4O5 (948.2 and 948.3). **No mass, formula or HRMS check can
tell any of these pairs apart.** The direction is carried by three independent
symbolic statements, all read at 600 dpi:

1. the full IUPAC name printed as each synthesis heading (SI pp. 21-41);
2. the reagent pair drawn in Supplementary Figure 15 (SI p. 17);
3. **which intermediate the write-up starts from** -- alcohol **2**, acid **3**
   or amine **4** -- which fixes the direction on its own, because the acid can
   only contribute the carbonyl and the alcohol can only contribute the oxygen.

Check `S2` in `07_qc.json` asserts the direction per lipid with a SMARTS and
was seeded with the lipid-34-as-lipid-36 swap; its agreement count moved
12 -> 11 on a defect that leaves the exact mass identical to three decimals.

## 4b. What this SI does NOT contain, searched rather than assumed

The linker is this paper's subject, so the usual linker tests were looked for
by name across all 42 SI pages before anything was claimed:

* **No ¹³C NMR.** `13C`, `C NMR`, `100 MHz`, `125 MHz`: **zero hits**. The
  "dialkyl carbonate at δ155 with nothing at 168-180" test is **unavailable
  here**, and asserting it would be inventing evidence.
* **No HRMS.** `HRMS`: **zero hits**. The General Information says only
  "Mass spectra were recorded on an ESI-TOF mass spectrometer", and every
  printed m/z is given to **one decimal** -- so there is no ppm-level anchor
  either, and the best mass agreement available is the ±0.04 Da systematic
  offset in check `S3`.

The ¹H equivalent that *is* available is check **`S5`**. An O-CH₂ of an ester
or carbonate sits at δ 4.0-4.25. Every member carries two tail esters, a fixed
4 H at δ 4.05. An ester or reverse-ester adds **one** more O-CH₂ (2 H); a
**carbonate adds two** (4 H), because both of its oxygens are alkylated; an
amide, urea or reverse-amide adds **none** and instead shows an N-CH₂ at
δ 3.3-3.5 and an N-H. Observed vs predicted, **11/11**: lipids 33 and 39
(carbonate) 8 H, lipids 32/37/38 (ester and reverse-ester) 6 H, lipids
34/35/36/40/41/42 4 H. Seeding lipid 33 with lipid 32's reverse-ester
structure moves the agreement count **11 -> 10**.

Carbamate is excluded by **valence, not preference**: the CDI acylation that
makes lipids 33 and 39 is between alcohol **2** and a head ALCOHOL, and every
nitrogen in both partners is tertiary -- the core N carries two C8 arms plus
the hydroxyethyl, the head N carries two methyls or is a piperazine N. There is
no N-H anywhere to acylate.

## 4c. A library limitation this curation did not hit

`curation/lib/assemble.check_stereo_preserved` compares CIP-code *counts*, so a
transformation that merely RELABELS a centre reads as an inversion and is
refused. **That refusal never fired here**, and it could not have: the only
stereocentres in these molecules are the two 2-hexyldecanoate C2 atoms, they
are remote from every bond this build forms, and they are carried through
unassigned on both sides (DECISIONS D-5). Nothing was worked around and nothing
in `curation/lib` was edited.

## 5. Checks that are external to our pipeline

The strongest are the ones computed through a quantity we did not digitise.
Supplementary Figure 2's legend prints a pKa for all 12 lipids **as vector
text** (tier 4, exact), and the authors then printed Pearson coefficients
against it:

| what | ours | printed | delta |
|---|---|---|---|
| HeLa luminescence vs pKa, middle dose | r = 0.7364, p = 0.0063 | 0.7358, 0.0064 | **0.0006** |
| HepG2 luminescence vs pKa, middle dose | r = 0.7191, p = 0.0084 | 0.7226, 0.0079 | 0.0035 |
| lung radiance vs pKa, raw scale | r = 0.8369, p = 0.0007, R² = 0.7005 | 0.8417, 0.0006, 0.7084 | **0.0048** |

Each was seeded, and each is credited only because the seed MOVED ITS
AGREEMENT COUNT. Rotating the lipid labels by one bar group takes the HeLa
coefficient from 0.7364 to **0.5106**, a move of 0.225, and the check's hit
count from 2/2 to 0/2.

The decade ladder of panel 3C needed its own evidence, and one of these checks
turned out NOT to supply it. Panel 3C's bottom axis label is eight pixels of
JPEG and reads equally as 10^5 or 10^6. Re-reading it as 10^6 -- three decades
over the same four ticks -- moves the lung Pearson coefficient by only
**0.034**, from 0.8369 to 0.8708, because a monotone rescaling barely perturbs
a correlation. **So `F3_2` is consistent with the ladder but does not settle
it, and saying otherwise would have been exactly the kind of clean line this
project's §6 warns about.** Two geometric checks settle it instead:

* `F3_3`: the four major ticks are equally spaced at 47.45 px (gaps 47.50,
  47.00, 48.00). A 10^8..10^6 reading would need one of them unlabelled.
* `F3_4`, the decisive one: a GraphPad log axis draws **eight minor ticks per
  decade at log10(2..9)**, which are UNEQUALLY spaced and therefore encode the
  decade length rather than merely being consistent with it. Seven are
  resolvable between the top two majors (two of the eight merge at 2.4 px
  apart). Mean distance to the nearest predicted position is **0.39 px** on
  the 47.45 px ladder and **1.63 px** on the 71.25 px alternative, a factor of
  4.2. The authors drew this and never printed it.

## 6. What could NOT be verified

* **The in-vitro dose series.** Methods and Figure 2 disagree irreconcilably
  (DECISIONS D-1).
* **Six values are upper bounds**, not measurements (D-3).
* **The Untreated control bars** in Figure 3C/D/E are hollow and were refused
  (D-4).
* **Figure 3F/G/H were not extracted.** They are the redundancy that would both
  check C/D/E bar-for-bar (log10 C - log10 D must equal F) and recover the six
  censored values. Writing that second digitiser was out of scope here; it is
  the highest-value thing left to do in this folder.
* **Lipid 14's structure rests on a drawing plus an outside corroboration**, not
  on this paper's own characterisation.
* **Stereochemistry is absent by decision** (D-5), because the source states
  none.
