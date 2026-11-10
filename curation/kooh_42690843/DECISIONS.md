# DECISIONS — kooh_42690843 (Chun, Nam, … Koo, *JACS* 2026, PMID 42690843)

Questions this curation must NOT answer on its own. `build_dataset.py` refuses
to merge this folder while any entry below is `OPEN`.

The pipeline is good at detecting that a source is silent. It is not entitled
to resolve it. In `AC_2025` the first repair was confident, internally
consistent, and would have assigned **419 measured values to the wrong
molecules**. So: detect, queue, stop.

---

### D-1. Which 275 of the 2,912 possible combinations are the library?

STATUS: OPEN

AMBIGUITY: The abstract states that **275** α-branched amine ionizable lipids
were generated. SI p.S5 draws **26 amines × 8 boronic acids × 14 aldehydes =
2,912** possible products. The SI names **9** members (verified by script:
`stage00_scan_for_values.py` finds 9 distinct `AxByCz` labels across all 49
pages) and never states the library's composition. 275 is 9.4 % of 2,912 and is
not a product of any subset of the panel sizes — 26×8 = 208, 26×14 = 364,
8×14 = 112, 26+8+14 = 48 — so no arithmetic recovers it. The design is almost
certainly staged (one block varied at a time), which is consistent with the SI's
own nine being exactly `A2/A18 × B1/B3/B8 × C1` plus `A18B8 × C8/C9/C14`, but
"almost certainly" is a guess and a guess here invents up to 266 molecules.

CANDIDATES:
  (a) take the member list from main-article Figure 2, which plots them;
  (b) treat only the 9 SI-characterised lipids as the library and say so;
  (c) reconstruct a staged design from the SI's own nine and the dashed
      grouping on p.S5 (A1–A15 monofunctional, A16–A24 bifunctional,
      A25–A26 trifunctional) — **rejected by this curation as fabrication**;
  (d) drop the paper.

EVIDENCE NEEDED: the main-article PDF, specifically Figure 2's axis labels.
Nothing else in the public record carries the member list — the paper has no
data-availability repository, is not open access, and is not in PMC.

---

### D-2. Does a folder with zero values merge at all?

STATUS: OPEN

AMBIGUITY: `05_values.csv` is empty by necessity, not by failure (see
`00_value_sources.json`). The nine structures are fully verified and their
quality columns are already written (`07_rows_quality.csv`), but a corpus row
needs an `Experiment_value` and there is none. The corpus has no precedent for
a structure-only contribution.

CANDIDATES:
  (a) hold the folder unmerged until the main article is acquired — the
      structures and checks are durable and nothing is lost;
  (b) merge structures with blank values, which would put 9 rows into
      `data_vitro.csv` that no analysis can use and that every value-based
      check would have to special-case;
  (c) keep the folder as a standalone reference set, the way `XD_2026`'s
      LipidDB was kept out of the corpus.

EVIDENCE NEEDED: a decision from the project owner; and, for (a), the main
article PDF.

---

### D-3. The screened material is CRUDE, not purified. Is a crude one-pot mixture a corpus "lipid"?

STATUS: OPEN

AMBIGUITY: SI p.S26 (Figure S4) says the in vivo screen evaluated "the **crude**
PILL lipid-LNP sub-libraries", and p.S27 (Figure S5) distinguishes those from
"**purified** PILL ionizable lipids". The abstract likewise says the platform
"was directly integrated with **purification-free** in vitro screening". So the
value attached to a label like `A18B8C14` in Figures 2–3 is the activity of an
unpurified reaction mixture whose nominal major component is that structure;
only Figure 4 and the nine SI compounds are purified material. The corpus
convention treats `IL_SMILES` as *the molecule that was tested*. For the
screening figures that is an approximation, and the SI quantifies the problem:
p.S4 Figure S2 exists solely to show that a DCM precipitation step removes
"bulk insoluble impurities", i.e. the crude contains a lot of something else.
Only `A18B8C14` has a stated purity anywhere (HPLC 96.00 %, p.S19).

CANDIDATES:
  (a) admit crude-screen rows with `structure_assumed = 1` and a note saying
      the measurement is of an unpurified mixture;
  (b) admit only the purified compounds (Figure 4 / the nine);
  (c) admit crude rows with no flag — **rejected**: it would be the first time
      the corpus silently stored a mixture as a molecule.

EVIDENCE NEEDED: the project owner's call on whether purification-free
screening data enters the corpus at all; this question is not specific to this
paper and will recur (it is the same class as the `NG_2025` crude Mannich
mixtures already carrying `positional_ambiguity = 1`).

---

### D-4. Double-bond geometry on C4, C5 and C7.

STATUS: OPEN

AMBIGUITY: The bond-line check in `stage02_blocks.py` verifies every block's
topology against the drawing bond-for-bond, but two parallel lines are two
parallel lines — it is blind to E/Z. `C4` ((2E,6Z)-nona-2,6-dienal), `C5`
((2E,4E)-deca-2,4-dienal) and `C7` (farnesal) carry their geometry **only**
from a human read of the 600 dpi render in `00_pages/`, and none of the three is
pinned by a reagent amount (they are not used in any of the nine SI compounds,
so no mass evidence exists for them). The corpus already holds 212 rows with an
unassigned stereogenic C=C, so storing them unassigned is a precedent, not a
novelty.

CANDIDATES:
  (a) store as read, flagged `structure_assumed = 1` with a note naming the
      double bond;
  (b) store with the alkene geometry stripped (unassigned), which is honest and
      lossy;
  (c) resolve from the commercial catalogue entry for each aldehyde, which is
      external evidence this curation did not consult.

EVIDENCE NEEDED: only matters if D-1 resolves such that a C4/C5/C7 lipid is
emitted. Defer until then.

---

### D-5. Licence. This source is not known to be CC0/CC-BY compatible.

STATUS: OPEN

AMBIGUITY: `CLAUDE.md` §8 requires every constituent source to be
CC0/CC-BY-compatible before corpus release. Europe PMC reports
`isOpenAccess = N`, `license = null`, `availability = "Subscription required"`,
and the landing page could not be loaded to capture the licence string (HTTP
403, Cloudflare). The SI PDF we hold is watermarked as served through an
institutional subscription session.

CANDIDATES: (a) capture the licence string from the landing page and decide;
(b) exclude the paper from any released corpus while keeping it for internal
analysis; (c) contact the authors.

EVIDENCE NEEDED: the article landing page, or the ACS AuthorChoice status.

---

### D-6. `Experiment_ID`, and how many experiments this paper is.

STATUS: OPEN

AMBIGUITY: `HK_2026` is proposed (Heebeom Koo, corresponding author), but
`HK` is not currently taken and the first author is Simin Chun with Hyorim Nam
as co-equal first — the corpus convention (`INITIALS_YEAR`) does not say whose
initials. Separately, this paper contains at least three distinguishable
screens that the corpus would have to split or pool: an **in vitro** HeLa
luciferase screen of the crude library (Figure 2), an **in vivo** ex vivo organ
screen of crude sub-libraries (Figure 3), and an **in vivo** validation of
purified lipids (Figure 4, n = 3 and n = 4 in different panels per SI p.S27).
Figures 3 and 4 are different material (crude vs purified) at the same dose and
endpoint, which is the configuration that produced `LL_2026` / `LL_2026_2`.

CANDIDATES: (a) `HK_2026` vitro + `HK_2026` vivo, Figures 3 and 4 pooled;
(b) `HK_2026` (Fig 2, vitro) / `HK_2026` (Fig 3, vivo, crude) /
`HK_2026_2` (Fig 4, vivo, purified); (c) initials from the first author.

EVIDENCE NEEDED: the main-article figures, and the project owner's convention
call on the initials.

---

### D-7. `A2B8C1`'s reported HRMS is −123 ppm from its own calcd.

STATUS: OPEN

AMBIGUITY: SI p.S14 prints `[M+H]+ calcd for C22H49N2, 341.3890; found,
341.3471` — a −41.9 mDa, −122.7 ppm error. The other eight lipids in the same
SI span −7.1 to +6.1 ppm. Our structure reproduces the **calcd** to +0.03 mDa
and is independently confirmed by the compound name
("N1-(heptadecan-9-yl)-N3,N3-dimethylpropane-1,3-diamine") and by the drawn
bond-line count (23 lines, measured). So the structure is not in doubt; what is
in doubt is whether the instrument reading was transcribed correctly, and
therefore whether this compound's identity was actually confirmed by mass.

CANDIDATES: (a) emit with `confidence = high` and a note recording the
discrepancy; (b) emit with `confidence = medium` because the only external
mass confirmation for this one compound is self-inconsistent; (c) treat the
printed `found` as a typo for 341.3871 and say so.

EVIDENCE NEEDED: the raw HRMS trace on SI p.S42, which is a raster image; a
reviewer could read the labelled peak off it at 600 dpi. Not done here because
it is an isolated reading and reading it would make this curation both the
source and the judge of its own correction.
