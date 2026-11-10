# DECISIONS — curation/songx_39664425

PMID 39664425 · PMC11628822 · doi 10.1016/j.apsb.2024.08.032
*Pioneering integration of combinatorial chemistry and machine learning to accelerate
the development of tailored LNPs for mRNA delivery* · Acta Pharm Sin B 2024
He X, Wang P, Qing L, Song X (Sichuan University)

Outcome: **0 rows emittable.** The article is a Highlight commentary.

---

### D-1
STATUS: SETTLED

**Do not digitise Figure 1 panel 4, although its axes are calibratable.**

Panel ④ carries four real bar charts: total flux in muscle (p s⁻¹, linear axis to
5×10¹⁰, 119-23 vs MC3 vs SM102, P = 0.0024 and P = 0.0367 printed) and %tdTomato⁺ in
liver, spleen and lung across 8–9 immune/stromal subsets with per-comparison P values
printed. Tier 5 is technically reachable.

It is refused on provenance. The caption attributes the panel to reference 1 *by
markup* — `Reprinted with the permission from Ref.` followed by
`<xref rid="bib1">`, and `bib1` carries `pub-id pmid 38740955`. That paper is already
in the corpus as **`BL_2024`** (560 vitro + 51 vivo). Emitting these values under a new
`Experiment_ID` would give one laboratory's measurements two source keys, which is
precisely what `CLAUDE.md` §9 refused for the `XD_2026` LipidDB ("one `Experiment_ID`
would fabricate a lab"). Because `Experiment_ID` is the environment key for every
result in this project, a fabricated second source would directly corrupt Results 1–4.

It is also strictly worse data: the served raster is 740×917 px with EXIF recording a
2960×3666 original, i.e. a ~4× downsample of a figure whose primary version is already
curated.

### D-2
STATUS: SETTLED

**No `make_reference.py`; `reference_structures.csv` carries one empty-anchor row.**

A C10 anchor tests an enumerated structure against the structure the paper *draws*.
This paper enumerates nothing and draws nothing of its own — every structure in
Figure 1 is inside the licensed reprint. An anchor built from it would be an anchor for
`BL_2024`, written in the wrong folder.

A generator script here could only succeed, and a code path that can only pass is the
exact anti-pattern `CLAUDE.md` §6 names ("my first two versions of C16 itself reported
'all derived' while examining nothing"). So none was written, and the *absence* of an
anchor is recorded explicitly in `reference_structures.csv` with
`expected_mismatch = no_anchor_possible` so the gap is visible rather than silent.

Anchor posture is recorded by hand in `PROVENANCE.md`, as instructed.
`provenance/structure_check.py`, `check_anchor_independence.py` and
`check_anchor_coverage.py` were **not** run — they write tracked files other agents read.

### D-3
STATUS: OPEN

**Is `BL_2024` missing the lead compound 119-23?**

Both the commentary text and Figure 1 identify `119-23` as the ML-designed lead of
PMID 38740955, the compound that beat MC3 and SM102. No `IL_name` beginning `119`
exists anywhere in `data_vitro.csv` or `data_vivo.csv`, and prefix `119` is outside
`BL_2024`'s 128-prefix grammar — consistent with 119-23 coming from the 40,000-member
*virtual* library and being synthesised separately for validation, rather than
belonging to either synthesised screen.

This belongs to the `BL_2024` curation, not to this folder, and nothing outside this
folder was touched. Handed to the caller to route. Note that this curation cannot
supply the value: the only number for 119-23 reachable here is inside the reprint
refused under D-1.

### D-4
STATUS: OPEN

**The `lit_review` tracker entry for this PMID is wrong and should be corrected.**

`lit_review/build_data.py` records this PMID as `status="new"`, `prio=1`, promoted
2026-10-26 with the note *"full text confirms 584 experimentally made ILs"*, a
`lib` size of 584, and a claim that it is *"the predecessor to PMID 42592089"*
whose data may sit in the LipidAI repo.

All of that is wrong, three ways: the 584 ILs belong to PMID 38740955 (already curated
as `BL_2024`), 584 is a count of screen memberships rather than distinct lipids
(see C-4 — the distinct count is 560), and the paper is unrelated to LipidAI.
Recommended status: **`cite`**, or `dup` against `BL_2024`.

Not applied — `lit_review/` is outside this folder and the brief is explicit that only
this directory may be written.

### D-5
STATUS: SETTLED

**`emit.emit_rows` was not called.**

With zero rows the call could only succeed, so it would test nothing, and the
per-item-callable requirement (a hard-coded literal `confidence` must raise, F-011) is
unreachable here. Stated in `07_qc.json` rather than papered over with a vacuous call
that would have produced a clean line while examining nothing.
