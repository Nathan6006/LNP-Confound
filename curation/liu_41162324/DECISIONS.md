# DECISIONS — curation/liu_41162324

Liu et al., *Fluorinated Ionizable Lipids for Efficient Spleen-Targeted mRNA
Delivery in Cancer Immunotherapy*, J Am Chem Soc 147(45), 2025.
PMID 41162324 · DOI 10.1021/jacs.5c13982.

The pipeline detects that a source is silent or self-contradictory. It does not
get to decide what the source meant. Every entry below is OPEN and the merge
gate refuses this folder until a human closes it. Nothing here was resolved by
the curation.

Supporting numbers: `07_qc.json`. Both loci of every disagreement:
`08_contradictions.csv`. Where the values are and are not: `00_value_sources.json`.

---

### D-1. The abstract says 315 lipids; the SI's own roster is 74.

STATUS: OPEN

AMBIGUITY: The brief states the abstract claims **315** lipids. The SI, which is
the only material held, is self-consistent at **74** across four independent
loci: 42 numbered synthesis headings in section 3.2 plus 32 SP-series members
named in the NMR/MS figure captions; the NMR figure captions S12–S29; the
ESI-MS figure captions S30–S47; and the x-axis of Figures S48–S51, whose tick
marks were counted mechanically at **75** (SM102 + 74 FILs, pitch CV 0.028,
seeded and confirmed to drop to 74 when one tick is erased). The strings
`315`, `library`, `libraries`, `combinatorial`, `total of` and `screen` occur
**0 times** in the 100-page SI. 315 is not 74, 75, 2×74, or any arithmetic this
curation can construct from the blocks (4 cores × 10 linkers × 4 fluoroalkyls =
160 is the full cross-product, of which 74 were made).

CANDIDATES:
 (a) 315 counts something other than distinct ionizable lipids — e.g. LNP
     formulations across cell types or doses, or a virtual/screened set;
 (b) 315 counts a larger synthesised library of which the SI characterises 74;
 (c) the abstract's number is an error.

EVIDENCE NEEDED: the article text. Specifically the sentence in which 315
appears and the main-text figure that shows the screen — if that figure draws
315 columns, (b); if it draws 74 and 315 is a formulation count, (a). Until
then `03_enumerated.csv` holds 74 and makes no claim about 315.

---

### D-2. The SP series' lipid tail: 2-heptylnonan-1-ol or 2-heptyldecan-1-ol?

STATUS: OPEN

AMBIGUITY: This decides the structure of **32 of the 74 lipids** — the entire SP
series. SI 3.2.2 is three sentences for the whole sub-library and says
"2-Heptyldecan-1-ol was used as the starting material" (C17). Every other
locus says C16: the Materials list (SI p.3) stocks `2-heptylnonan-1-ol (97%)`
and does not stock 2-heptyldecan-1-ol; SI 3.1.47's heading names the product
`2-heptylnonyl 8-((2-hydroxyethyl)amino)octanoate`; its body charges
`2-heptylnonan-1-ol (4.8 g, 20.0 mmol)`; and its stated product is named the
same way. **No SP member has a molecular formula or an m/z anywhere in the SI**,
so nothing the authors measured discriminates the two.

This is not a cosmetic difference. 2-heptyldecan-1-ol is C17H36O, isomeric with
the SS series' heptadecan-9-ol, so candidate (b) would make **every SP lipid an
exact isomer of its SS counterpart** — identical formula, identical exact mass,
and a corpus in which `IL_name` and structure diverge in the way CLAUDE.md §3
warns about. Candidate (a) makes them one CH2 lighter and plainly distinct.

CANDIDATES:
 (a) 2-heptylnonyl (C16) — built here; 3 loci, including the purchase record;
 (b) 2-heptyldecyl (C17) — 1 locus, the SP section's own sentence.

EVIDENCE NEEDED: any SP-series mass. The article may print one for the SP lead
compound; failing that, the ESI-MS rasters in Figures S38–S44 carry peak
labels that were not read in this curation (they are spectra, and reading a
peak label off a 1152-px raster is tier 5 for a number that tier 4 should
supply). Reading one SPC3F1 peak would settle it, and the two candidates
are well separated: **[M+H]+ 702.568 for (a)** (C40H76FNO7) against
**716.584 for (b)** (C41H78FNO7). Note that (b) is *exactly* the mass the
SI prints for SSC3F1 (716.58) -- so under (b) an SP spectrum and its SS
counterpart would be indistinguishable, which is a further reason the
answer must come from a locus that is not a mass.
Until then all 32 rows carry `confidence = low`, `structure_assumed = 1` and a
note naming both loci.

---

### D-3. Eight lipids where the SI's formula AND its own m/z agree with each other and disagree with the structure the SI describes.

STATUS: OPEN

AMBIGUITY: This is the `AC_2025` shape — a name and a number disagree and
nothing in the data says which half to repair. In `AC_2025` the confident,
internally consistent first repair would have put 419 measured values on the
wrong molecules.

Two groups, both at exactly one CH2 (14.016 Da):

 * **The P series (SSCPF1/F3/F5/F7, and by the series rule SPCPF1/F3/F5/F7).**
   The SI prints C46 formulas and, for PF3/PF5/PF7, [M+H]+ values that match
   those formulas to ≤0.02 Da. Against that: the reagent is named
   "…-3-oxopropoxy)**benzyl** 6-bromohexanoate"; the intermediate at SI p.11 is
   "3-(4-(((6-bromohexanoyl)oxy)**methyl**)phenoxy)propanoic acid" with
   ¹H NMR 5.04 (s, 2H); SSCPF1 and SSCPF5 each show their own 2H singlet near
   5.0 ppm; and the Materials list stocks **4-hydroxybenzyl alcohol**. A
   structure matching the printed C46 would be a *phenyl* ester with no
   benzylic CH2 and therefore no such singlet.
 * **SSCB4F7.** Printed C45H80F7NO7 and [M+H]+ 880.59 agree to 0.001 Da.
   Against that: its own reagent sentence names a **nonanoate** (C9 → C46), and
   SSCB4F1/F3/F5 are all C46 and all match our build exactly.

CANDIDATES:
 (a) the printed formula and mass are wrong (both derived from one wrong
     formula) and the described structure is right — what is built here;
 (b) the described structure is wrong and the mass is right — which for the P
     series requires a phenyl ester that contradicts four other loci;
 (c) hold the eight rows out of any merge.

EVIDENCE NEEDED: the ESI-MS spectra in Figures S37 (SSCPF3, SSCPF5) and S35
(SSCB4F7), whose peak labels are drawn and were not read here; or ¹⁹F/¹³C data
the SI does not provide; or the article's own structural figure. All eight rows
carry `confidence = low`, `structure_assumed = 1` and a note naming both loci.

---

### D-4. This paper supplies ZERO corpus rows from the material held. Does it enter the corpus at all?

STATUS: OPEN

AMBIGUITY: `00_value_sources.json` enumerates all sixteen numeric SI panels.
Every one is either an endpoint the corpus excludes (size, zeta, PDI, viability,
serum chemistry, histology, antibody titre) or a 2–3 arm in-vivo comparison in
which one or two arms are the SM102 benchmark. **The 74-member screen's
transfection values exist only in main-text figures, which were not supplied**,
and the paper is not open access (Europe PMC `isOpenAccess: N`, `inEPMC: N`,
`hasSuppl: N`; a plain request to pubs.acs.org returns HTTP 403 behind a
Cloudflare interstitial). So this folder currently holds 74 structures with
their flags and no `Experiment_value`.

CANDIDATES:
 (a) acquire the article PDF and extract the screen, then merge;
 (b) merge the structures with blank values, as a chemistry-only contribution;
 (c) leave the folder unmerged until (a) is possible.

EVIDENCE NEEDED: a decision on (a) — which requires a browser-based fetch of
pubs.acs.org (`BROWSER_ACQUISITION.md` records ACS as browser-but-no-login) or
an institutional copy. (b) is noted as a candidate only; a row with a blank
`Experiment_value` was the reason 328 in-vitro rows were excluded at Stage 1 of
the frozen pipeline, so it is not obviously admissible.

---

### D-5. If the in-vivo SI panels are ever digitised, they are benchmark-only comparisons.

STATUS: OPEN

AMBIGUITY: Figures S52, S53 and S58 quantify luciferase for SSC6F5 against
SM102 (and against SM102/18PA, a SORT formulation of the same ionizable lipid).
That is one novel lipid versus one or two benchmark arms. CLAUDE.md §8 already
has this question open for the in-vivo pool: "a source whose only shared lipids
are benchmarks that top every screen can push cross-lab concordance toward +1
for reasons that have nothing to do with transferable SAR", and asks that it be
decided before Stage 5 is re-run.

CANDIDATES:
 (a) exclude benchmark-only panels from any future merge of this paper;
 (b) include them and flag the rows;
 (c) defer to whatever is decided for the corpus-wide question in CLAUDE.md §8.

EVIDENCE NEEDED: the corpus-level decision in CLAUDE.md §8. Nothing specific to
this paper is missing; this entry exists so that a later merge of this folder
cannot slip the rows in without the question being asked.
