# DECISIONS — goulart_41138890 (PMID 41138890, proposed `MR_2025`)

The pipeline detects that a source is silent; it does not get to decide what the
source meant. Everything below that is `OPEN` blocks the merge gate, which is
the intended behaviour, not a failure of this curation.

Format note: the `STATUS:` line is deliberately plain, not bold — bold made the
merge gate report zero open decisions for a folder with five.

---

### D-1 — Should `MR_2025` be merged at all? All three of its ionizable lipids are universal controls already in the corpus.

The paper contributes **three** distinct ionizable lipids — C12-200, SM-102,
ALC-0315 — and **zero** of them are new. Measured by InChIKey join:

| lipid | vitro experiments today | vivo experiments today |
|---|---|---|
| C12-200 | `AC_2025` (498, name `C12-200`), `JW_2024` (20, name `other`) | none |
| SM-102 | `AC_2025`, `WT_2026` (2,295 rows) | `SB_2024`, `SS_2018`, `SY_2025`, `XD_2026` |
| ALC-0315 | `AC_2025`, `LZ_2025` (312 rows) | `XD_2026` |

Merging would take each from 2 vitro sources to 3. That is exactly the hazard
`CLAUDE.md` §8 already has open for the in-vivo pool: *"a source whose only
shared lipids are benchmarks that top every screen can push cross-lab
concordance toward +1 for reasons that have nothing to do with transferable
SAR."* Here it is **3 of 3**, the purest form of it.

Two things cut the other way, and they are why this is a question and not a
refusal:

1. **The benchmarks do NOT top this screen.** C12-200 beats SM-102 by
   **2.2375×** and ALC-0315 by **3.2630×** at 72 h. A benchmark-only source
   that *disagrees* with the "controls win" prior is informative about
   cross-lab rank disagreement rather than degenerate.
2. **The cell type is new.** `primary_rat_cortical_neuron` appears nowhere in
   the corpus's 17 in-vitro `Model_type` values, and no corpus row of any kind
   is neuronal (`Model_target = brain` exists only on 51 `ZT_2026` **vivo**
   rows). This is the first neuronal in-vitro environment.

But see D-6: the ranking is confounded.

STATUS: OPEN

---

### D-2 — `Experiment_method = GFP_fluorescence_normalized` is a new token in a three-value vocabulary.

The in-vitro corpus uses exactly three method strings:
`luminescence_normalized` (24,072), `luminescence_discretized_normalized`
(5,011), `LRP6_knockdown_normalized` (112). This readout is GFP fluorescence
intensity per area from pDNA, which is none of them.

There is precedent for a bespoke token — `LL_2012`'s 112 rows are pDNA/GFP in
HEK293T under `LRP6_knockdown_normalized` — so the vocabulary is already
per-paper rather than controlled. Adding a fourth token is cheap; silently
reusing `luminescence_normalized` for a fluorescence readout would not be.

STATUS: OPEN

---

### D-3 — `Model_type = primary_rat_cortical_neuron`: a new cell type, and the second rat cell in the corpus.

The only other rat-derived line is `INS-1E` (123 rows). Everything else is
human or mouse. A primary rat culture is also not a line: it is a mixed
preparation, which the paper itself quantifies in Fig S3 (NeuN+, GFAP+, DAPI)
and in Fig 4C — where the GFAP+/DAPI+ (astrocyte) bar is the **tallest of the
three**, above both NeuN+ and total GFP+. The percentages were not digitised
(that panel supplies no corpus value), but the ordering is unambiguous on the
rendered figure, and the paper marks both comparisons significant. So the
`Model_type` string would say "neuron" for a culture whose most-transfected
population is measurably not neuronal.

STATUS: OPEN

---

### D-4 — `PEG_SMILES` is blank because the corpus name resolves to the wrong molecule.

Figure 1A draws `C14PEG2k` as a 1,2-diacyl-*sn*-glycero-3-phospho-ethanolamine
*N*-carbamoyl-mPEG45 (ammonium salt): a phosphate, an `-NH-C(=O)O-` carbamate,
and the subscript **45** printed inside the bracket. The corpus name
`C14-PEG2000` (`BL_2023_3CR`, 480 rows) resolves **unambiguously** to InChIKey
`NLDCRRFPGOSAHV`, formula `C122H242O50` — a glycerol ether with no N and no P.

So `corpus_components.smi()` would **not** raise: the name is unambiguous in
the corpus and simply wrong for this paper. `provenance/reagent_registry.py`
has no PEG-lipid by design (the block length is a distribution). The acyl
chains are drawn as zig-zags and would have to be counted off a raster, so no
structure is written.

This is the eighth sighting of the shape behind the open **F-049**, and this
curation measured it directly: `DMPE-PEG2000` itself carries **two** structures
in the corpus — the real phosphoethanolamine `OJTFVTWLLKMOJL` on 3,611 rows
across 8 experiments, and the glycerol ether `NLDCRRFPGOSAHV` on 1,873 rows in
`JW_2024` alone.

STATUS: OPEN

---

### D-5 — `Dose_ug_nucleicacid` left blank for Figure 4G.

Neighbouring panels use three different doses — Fig S5 says 0.025 µg/well,
Fig S6 says 0.01 µg/well, Fig 3B is a dose-response over 0.005–0.8 µg/well —
and Fig 4G's own dose appears in no source we hold (the article is paywalled;
the SI caption for Fig 4 does not exist because Fig 4 is a main-text figure).
Reconstructing it from Fig 3B's optimum would be a guess dressed as a
derivation.

Blanked, following the precedent the user set on **F-004** (`WT_2026`'s N:P:
three mutually inconsistent candidates, blanked rather than chosen).

STATUS: SETTLED — blank, by F-004's precedent.

---

### D-6 — The three-lipid ranking is confounded by formulation.

The formulation used in Fig 4G (IL 35 / DOPE 30.5 / cholesterol 32 /
C14PEG2k 2.5) is `LNP4`, which was selected as the winner of an 18-member
screen **run entirely on C12-200**. SM-102 and ALC-0315 are then evaluated in a
formulation optimised for a different ionizable lipid.

So "C12-200 beats both FDA benchmarks 2.24× and 3.26×" is not a clean chemistry
comparison, and this is a milder form of the brief's *"a winner selected into a
second screen tops it tautologically"*: the winner here is the formulation, not
the lipid, but the lipid rides along with it.

The confound is written onto all three emitted rows in `curation_note`. Whether
a corpus row should carry a value produced under it is the question.

STATUS: OPEN

---

### D-7 — Only the 72 h timepoint is emitted, of three measured.

Fig 4G reports 24 h, 48 h and 72 h for each lipid. All nine values are recorded
in `05_values.csv`; only the three at 72 h are emitted.

Emitting all nine would put three conflicting `Experiment_value`s on one SMILES
inside one `Experiment_ID` — the shape of the corpus's 48 duplicate-SMILES
groups with conflicting values, which `CLAUDE.md` §3 already lists as an
unresolved problem. 72 h is the paper's terminal point and where its own
statistics and its stated conclusion sit.

STATUS: SETTLED — 72 h only; the other six values are recorded, not emitted.

---

### D-8 — The in-vivo and in-vitro arms are NOT pooled, and the in-vivo arm is not emitted.

The title promises two pools and they are genuinely different experiments:
in vivo is a pooled b-DNA barcode screen read by NGS from five organs, in vitro
is per-well GFP imaging in primary neurons. They are not one `Experiment_ID`
and no z is taken across them.

The in-vivo arm is not emitted at all, for a reason that is not about pooling:
**barcode recovered from tissue is biodistribution, not delivery** (brief §3).
The in-vivo *luciferase* arm (Fig 6C, brain Total Flux) IS a translated readout
and would be admissible — but LNP2/LNP3/LNP4 are all C12-200 at three DOPE
contents, so it carries one ionizable lipid and no chemistry varies.

STATUS: SETTLED — separate experiments, in-vivo arm not emitted, reasons above.
