# DECISIONS — siegwartdj_42150080

Sun, Pien, Xiao, … Olson, Siegwart. *A human cardiomyocyte screen identifies
optimized lipid nanoparticles for in vivo cardiac gene editing.* PNAS 2026,
123(21) e2528477123. PMID 42150080, PMC13214034.

`build_dataset.py` refuses a folder with an `OPEN` entry. Every entry below is
open, and the folder emits nothing, so there is nothing to block — but the
entries are what a later session needs in order to decide whether to reopen it.

---

### D-1. The per-member screen values are not in any held file. 0 rows emitted.

The abstract says the lead was found "by screening a chemically diverse LNP
library". Those per-member values are in the main article's Figure 1. The main
article is not held and was not obtainable:

- `www.pnas.org/doi/10.1073/pnas.2528477123` → HTTP 403 (Cloudflare "Just a
  moment…") to curl with a browser UA and to WebFetch.
- `pmc.ncbi.nlm.nih.gov/articles/PMC13214034/` → HTTP 403; the PMC record we
  did fetch says `pmc-status-live=no`, `pmc-status-embargo=yes`, release
  **2026-11-18**, and carries the comment *"The publisher of this article does
  not allow downloading of the full text in XML form."*
- The shared Playwright browser belongs to the main session; the task brief
  forbids its use here. No paywall or bot-check bypass was attempted.
- Tier 1 (publisher source data): PMC declares exactly one supplementary item,
  the `.sapp.pdf` we hold; the article's further `sm*` supplements are MOVIES
  (the JATS badge subject is `video`). Crossref registers no components.
- Tier 3 (author repo): `find_repo.py` → 0 repo hits, 0 surviving code hits.

Everything inside the held SI is a non-corpus endpoint (diameter, PDI,
viability), a single-formulation experiment (D-3), or a qualitative image.

**Recommendation: leave this paper uncurated until 2026-11-18**, then retry PMC
for the article PDF and re-run the intake gate. Expect Figure 1 at tier
`5_pixels`. Nothing in this folder needs redoing.

**If it is ever merged, it must be a NEW `Experiment_ID`, never pooled with
`KZ_2016`.** `KZ_2016` is siRNA in HeLa and in mouse liver; this paper is mRNA
in hiPSC-CMs and mouse heart. Different cargo and different assay (brief's trap
4).

STATUS: OPEN

---

### D-2. `18:1 TAP` / `18:1 DAP` / `18:1 PA` are mapped to DOTAP / DODAP / DOPA by catalogue convention, not by anything the SI says.

The Fig S2/S3 axes print the Avanti shorthand only. The expansion to
1,2-dioleoyl-3-trimethylammonium-propane, -3-dimethylammonium-propane and
sn-glycero-3-phosphate is the catalogue meaning of that shorthand, and it is
what makes the charge-class column in `01_blocks.csv` mean anything — but it is
an inference. It is recorded in `01_blocks.py` as `HELPER`, flagged there as
"NOT read from the paper", so a reviewer can object to the mapping itself
rather than to a number derived from it.

The main text's Methods would settle it. Not held (D-1).

Note for context, not as evidence: the corpus already carries DOTAP in
`HL_name` (2,526 vitro / 38 vivo rows) and DODAP in `fifthcomponent_name` (486
vitro / 51 vivo rows — the `ZT_2026` Siegwart DoE curation). DOPA is absent
from the corpus entirely.

STATUS: OPEN

---

### D-3. Fig S6B is a real delivery endpoint with one lipid. Emit it, or not?

SI p.7 Fig S6B: IVIS average radiance (p/s/cm²/sr) of thymus, lung, heart,
diaphragm and liver, PBS vs 18:1 TAP10 LNP-Luc, intrathoracic route, 6 h,
1 mg/kg total Luc mRNA, n = 2 (treatment group; the PBS group's n is not
stated — see `08_contradictions.csv`). The y-axis carries printed linear ticks
(0, 1×10⁷ … 4×10⁷), so it IS calibratable at tier `5_pixels`.

Argument against emitting, which is why nothing was emitted: one ionizable
lipid across five organs is a biodistribution block, not a screen. The
resulting `Experiment_ID` would hold a single unique `IL_SMILES`. It can
contribute nothing to shared-lipid concordance (Result 3) and nothing to LOSO
(Result 4), and the per-experiment z-score convention would be computing a
z-score over five organs of the same molecule.

The same argument covers Fig S7B (%tdTomato⁺ by cardiac cell population, one
formulation, one route).

Not my call to make alone — queued, per `PIPELINE_V2.md` §3.

STATUS: OPEN

---

### D-4. This paper draws no chemical structure anywhere, so there is no C10 anchor. Is the mass/NMR anchor an acceptable substitute?

`reference_structures.csv` carries a row for `6A1-SC8` with an **empty**
`reference_SMILES`, so the absence is visible rather than silently missing, as
`PIPELINE_V2.md` requires.

What the SI gives instead is external and strong, but it is a *composition*
anchor, not a *topology* anchor. From `07_qc.json`:

| check | what it compares | result |
|---|---|---|
| q1 | our formula ← corpus 6A1-SC8, vs the caption's printed string | 5/5 elements, `C110H203N5O24S6` |
| q2 | the caption's string vs the caption's own printed m/z | Δ 0.0042 Da (rounding) |
| q3 | our formula vs the caption's printed "found" | Δ 0.0142 Da |
| q4 | our formula vs the ladder DRAWN inside Fig S1B | 5/5 labels, worst 3.23 ppm, z=2 derived from 0.5006 Da spacing |
| q6 | integral ratios vs the structure's own H counts | 2/2 within 5% (0.79%, 3.5%) |

**Mass cannot see connectivity.** q1–q4 fix the molecular formula; q6 fixes one
ethylene-glycol diester and one octyl thioether per arm, and the S count fixes
6 arms. None of them can see which nitrogen of the pentamine core carries which
arm. A C10 anchor for the 6A1 series would have to come from the *synthesis
route* in the original paper (Zhou et al. 2016, PMID 26729861 = `KZ_2016`),
which is where the corpus structure came from in the first place — so it would
not be independent of the thing it checks.

STATUS: OPEN

---

### D-5. Licence is CC BY-NC-ND 4.0, which is not CC0/CC-BY-compatible.

Read two ways on 2026-11-10: the PMC JATS `<ali:license_ref
content-type="ccbyncndlicense">` and the Crossref record's single `license`
entry, both `https://creativecommons.org/licenses/by-nc-nd/4.0/`,
`content-version` "vor", `delay-in-days` 0.

`CLAUDE.md` §8 requires every constituent source's licence to be
CC0/CC-BY-compatible before corpus release. **NC** and **ND** are both
incompatible with that requirement as written. This does not block curation —
extracted factual values are not the licensed expression — but it does bear on
the Data & code availability section, and it should be decided for the whole
corpus rather than per paper.

STATUS: OPEN

---

### D-6. "6A1-SC8 is the ionizable lipid of all 16 formulations" is an inference.

What the held SI actually says: Fig S1 characterises 6A1-SC8 by ¹H-NMR and
LC/MS, and it is the only lipid synthesised or characterised anywhere in the
SI. Figs S2/S3 name a "Base LNP" and fifteen variants of it. The SI never
states the Base LNP's composition — that is in the main text's Methods, not
held (D-1).

`01_blocks.csv` therefore writes `6A1-SC8 (inferred, see note)` in the
`ionizable_lipid` column rather than a bare name. If the inference is wrong,
the only thing that changes in this folder is that column; no value depends on
it, because no value was extracted.

STATUS: OPEN
