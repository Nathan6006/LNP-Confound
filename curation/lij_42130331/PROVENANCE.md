# lij_42130331 -- provenance

**Status: EXTRACTED 2026-11-10. NOT merged; `DECISIONS.md` has five OPEN
entries and `build_dataset.py` refuses a folder with one.**

| | |
|---|---|
| Paper | Guo, Zhang, Yang, Liu, Sun, Luo, Wang, Meng, Pan, Wang, Yang, Sun, Huang, Gong, Li. *beta-Hydroxy Thioether-Derived Ionizable Lipids for Spleen-Tropic mRNA Delivery and In Vivo Chimeric Antigen Receptor T Cell Engineering*, ACS Nano 2026, 20(20):14602-14614 |
| DOI | 10.1021/acsnano.6c00706 |
| PMID | 42130331 (supplied by the caller; confirm with `provenance/fetch_esummary.py` before merge) |
| Experiment_ID | `QG_2026` -- first author Qiang Guo, matching the convention of `ZT_2026` (Tian), `LL_2026` (Lin), `XD_2026` (Duan), `QL_2026` (Liu). No collision with any existing id. |
| Rows | 300 emitted (`07_QG_2026_rows.csv`), 50 from Library 1 + 250 from Library 2 |
| Value source | SI Tables S1 and S2 -- the authors' own source data. See `00_value_sources.json`. |
| Licence | (c) 2026 American Chemical Society, subscription-only, NOT CC0/CC-BY. See `DECISIONS.md` D-5. |

## Where the values come from

The paper's own source-data tables, captioned verbatim "Source data for the
heatmap in Fig. 1C" (Table S1, 50 values) and "...Fig. 1E" (Table S2, 250
values). They are published as **flattened 1-bit raster images inside the SI
PDF**, so the tier is `5_pixels` -- but it is the most favourable kind: the
pixels encode printed numerals, not a colour ramp or a bar height. There is no
axis, no origin, no scale factor and no colormap anywhere in the read, so error
classes 16, 20, 23, 27, 32, 33 and 43 are **structurally inapplicable** rather
than merely checked.

Tiers 1-4 were each ruled out by an observation a reviewer can repeat:

| tier | why it failed |
|---|---|
| 1 publisher Source Data | The article page (loaded 2026-11-10, HTTP 200 on the second navigation) exposes exactly ONE supplementary asset, `/ancac3/article-supplement/5105918/pdf/nn6c00706_si_001/`, a PDF. Enumerated by querying every `<a href>` for `suppl_file\|supporting\|si_<n>\|.xlsx\|.zip\|.csv\|.xls`; the two hits are that same PDF under two labels. |
| 2 SI tables | The SI *contains* the tables and they *are* the source data -- but as images. `pdfplumber.extract_text()` on pages 38 and 39 returns 194 characters each, which is the caption plus the download watermark; `find_tables()` returns 0. |
| 3 author repository | Europe PMC core record for PMID 42130331 (2026-11-10): `isOpenAccess=N`, `inEPMC=N`, `inPMC=N`, `hasSuppl=N`, one full-text URL marked "Subscription required". No repository, accession or deposit. |
| 4 PDF text table | `page.get_drawings()` returns **0 paths** on every page of the file; each table page holds exactly one image object (1422x193 and 1421x1126, 1-bit PNG). |

The decode does not use OCR. CLAUDE.md section 5: "OCR is not good enough.
Tesseract misread 1.2x10^5 as 109." Instead every connected component in every
cell is normalised and clustered with complete linkage: **2,243 glyphs fall
into 45 clusters whose maximum intra-cluster distance is 0.000** -- byte-identical
bitmaps -- with the nearest two clusters 2.236 apart. The 45 cluster
representatives were labelled once by eye off a montage that shows each at its
native pixel size magnified 10x, and the labels are keyed by the sha256 of the
representative's bitmap, so a label cannot drift onto a different shape.

28 of the 45 clusters are two or three characters that touch at this
resolution. Labelling those by eye would put one person on both sides of the
comparison, so `stage04_tables.py --verify` re-derives each from the
single-character templates alone, by a dynamic programme that never sees the
eye label: **25 of 28 reconstruct with a residual of exactly 0 pixels.** The
three that do not -- `AA`, `OG`, `3A` -- are the only clusters containing A, O
or G, which occur nowhere in isolation in either table, so the template
alphabet cannot spell them; their residuals are 90-159 px against 0 for every
confirmed one. All three appear only in NAME cells and each is printed
independently in Fig. S1B, Fig. S2 and Fig. S3.

## The paper's own cross-checks

**Passing all of them is necessary and not sufficient** -- `QL_2026` reproduced
every printed hit rate on bars up to 28% too high. Here the values are not
measured from geometry at all, so the risk is not a multiplicative bias but a
mis-mapping of cells to compounds, and the checks are chosen for that.

| check | result |
|---|---|
| Stated library size (abstract: "300 distinct ionizable lipids") | 10 amines x 5 epoxides = 50, x 5 acyl chlorides = 250, total **300**, all with distinct canonical SMILES |
| Values present for every enumerated structure | 300 of 300, with none left over in either direction |
| **The lead compound is the maximum.** The paper names `113-AA-C8C14` throughout. That name decomposes to Table S2's column `113-AA-C8` and row `Cl14`, and nothing in the decoding knows which cell should win | It is the **unique maximum of all 250** (9.8E+04; the runner-up, 9.3E+04, is the same column's `Cl12`) |
| Seeded failure of that check: `ROW_TO_ACYL` rotated by one row | REFUSED -- "the maximum is 113-AA-C8C6". The check has teeth. |
| Seeded failure of the glyph-label gate: one cluster's label deleted | REFUSED, exit 1, nothing written |
| Seeded failure of C16: `TEMPLATE_113` replaced by the enumeration's canonical SMILES | FLAGGED -- "1 TEMPLATE(S) APPEAR TO BE TRANSCRIBED FROM THE ENUMERATION" |

There is no printed statistic in the SI to reproduce -- no hit rate, no stated
fold-change over the library, no marginal. The lead-compound-is-the-maximum
test is the one external anchor this paper offers, and it is a good one
precisely because it is a claim made in words and tested against geometry.

## Chain lengths

The vector route does not exist here (`get_drawings()` = 0 paths), so two
independent RASTER routes are used and must agree: Douglas-Peucker vertices via
`chains.raster_chain_length`, and arc length divided by the bond length
measured from each region's own longest chain with the label clearance added
back. A third route (per-column ink-centroid turns) is recorded but not used
for agreement, because several arms contain near-vertical and near-horizontal
bonds for which `y(x)` is not single-valued.

**Six chain lengths rest on one automated route plus a human reading a 2000 dpi
render**, and that is error class 48 exposure recorded rather than hidden:
amine 75's two N-ethyls, amine 96's propylene, one of amine 400's two
propylenes, one of T3A's three arms, and the branched 2-ethylhexyl of OG8
(where both automated routes follow the longest path and miss the ethyl branch
entirely). The renders are committed in `00_pages/zoom_*.png`. For T3A the
other two arms measure 3 and 3 with both routes agreeing.

**Name tokens are not chain lengths, and here the mismatch is the chemistry.**
`C8` is 1,2-epoxyoctane and leaves a SIX-carbon R on the carbinol carbon --
measured as 6, both routes agreeing. `OG8` is 2-ethylhexyl glycidyl ether, 11
carbons in the reagent and 8 in the group the token names. The acyl tokens do
match and are not a measurement at all: Fig. S2 *prints* n = 3, 5, 7, 9, 11
against Cy = C6, C8, C10, C12, C14.

## Contradictions

See `08_contradictions.csv`: six. The one that would have cost values is the
amine column order -- Table S1's header puts **400 before 96** while Table S2's
block headers and Fig. S3's axis groups both put **96 before 400**. Every
column is labelled in full and the decoder keys on the label, so nothing moved;
a positional read would have swapped ten Library-1 values between two amines
with no downstream symptom.

## C10 anchors

`make_reference.py` BUILDS two anchors -- `113-AA-C8C14` (the lead, two
AA arms) and `96-AA-C8C6` (the two-arm case, where a secondary N-H opens an
epoxide directly and the arm has neither sulfur nor amide). Both match their
enumerated structure on the **full** InChIKey, stereo block included.

The generator shares as little as possible with the enumerator: it contains no
`RWMol`, no `AddBond`, no `RemoveAtom` and no atom map, where the enumerator is
built entirely on those; it composes one SMILES string centre-out from the
tertiary nitrogen, where the enumerator hangs arms off a mapped amine and
canonicalises starting at an acyl tail; the acetamide is spelled `NC(=O)C` here
and `NC(C)=O` there. C16 reports both templates as derived, and was made to
fail on purpose first.

The emitted anchor strings ARE byte-identical to the enumeration. That is
expected and is not the signal: a correct re-derivation canonicalises to the
same string, which is exactly why C10's own `anchor_may_be_a_copy` missed F-037
and why C16 tests the generator's composed template instead.

**What the anchors do not prove.** They re-derive the drawing. They cannot test
the drawing, because this paper publishes nothing against which the drawn
connectivity could be checked -- no HRMS, no NMR table, no SMILES/InChI, no
structure file. Mass cannot see connectivity, and here there is not even a
mass. That is `DECISIONS.md` D-4.
