# Validation review — F-013 (ZT_2026 steroid C3 epimer)

Reviewer: validation session, 2026-10-30. Read-only except this file.
Scope: `curation/siegwart_42634137/` (Tian et al., *Adv Mater* 2026, PMID 42634137),
`curation/build_new_curation.py`, `data/data_vitro.csv`, `data/data_vivo.csv`,
`provenance/structure_check.py`, `provenance/reagent_registry.{py,json}`.

---

## 1. VERDICT

**Correct but incomplete.** It is the right half: the paper uses ordinary
3beta-cholesterol, never epicholesterol, and nothing in the main text, SI,
figures or NMR compares the two epimers — so rewriting all 51 structures to
3beta destroys no real distinction. The proposed `_cholesteryl_ester` rewrite is
chemically and computationally correct: I rebuilt it independently and its output
is identical, stereocentre for stereocentre, to a plain-SMILES ground truth that
never touches `RWMol`, for all three linkers, and end-to-end through
`quaternize()` all 36 enumerated lipids then liberate `HVYWMOMLDIMFJA-DPAQBDIFSA-N`.

It is incomplete for one reason: **the repair as proposed touches only
`curation/siegwart_42634137/enumerate.py`, which accounts for 36 of the 51 unique
structures (522 of the 537 rows). The other 15 unique structures (15 rows) — every
`MC3-CB*`, `SM-102-CB*`, `OF-2-CB*`, `TT3-CB*`, `Lipid 14-CB*` — are built by a
verbatim second copy of the same bug at `curation/build_new_curation.py:179-190`,
which the proposal does not mention.** Fixing only `enumerate.py` leaves 15
epicholesterol structures in the corpus and, worse, leaves them there *after* a
fix has been announced. Both copies must change in the same commit.

Two smaller caveats, neither fatal: the proposed snippet drops the existing
`if oh is None: return None, None` guard (bare `next(...)` raises
`StopIteration`), and the repair does not address the fact that `quaternize()`
bonds the amine to *atom index 0 of the canonical SMILES* of the returned
fragment — a silent assumption that happens to hold (verified below) but is the
same species of implicit-index assumption that produced this defect.

---

## 2. EVIDENCE FROM THE PAPER

### 2.1 Is there any epi / 3alpha / allo sterol anywhere? No.

Exhaustive search of both text layers (main PDF 15 pp., SI 23 pp.) for
`epi`, `allo`, `3-alpha`, `3α`, `3β`, `3-beta`, `cholestan`, `CAS`:

- SI: **zero hits** for any of these terms.
- Main text: the only hits are `allowing` (p.5) and `allow` (p.12) — ordinary English.

Every one of the ~45 occurrences of the word in the main text is plain
"cholesterol" / "cholesteryl" / "cholesterol-based" / "cholesterol tails".
There is no epimer arm, no epicholesterol control, and no comparison of C3
configurations anywhere in this paper. **Rewriting 3alpha -> 3beta destroys no
distinction the authors drew.**

### 2.2 Reagent provenance

Main text, §4.1 Materials (p.11, right column):

> "All chemicals and reagents were of analytical grade and used without further
> purification. Amines, **cholesterol**, Evans blue, Dulbecco's modified
> phosphate buffered saline (PBS), RPMI-1640 medium, fetal bovine serum (FBS),
> and trypsin-EDTA (0.25%) were obtained from **Sigma-Aldrich**."

No CAS number, no catalogue number, no stereochemical qualifier. "Cholesterol"
from Sigma-Aldrich is CAS 57-88-5, cholest-5-en-3beta-ol; epicholesterol is a
separate, far more expensive catalogue item (CAS 474-77-1) that is always named
as such. This is supporting, not decisive, evidence — the decisive evidence is
below. Note also that §4.2/§4.3 describe only the Michael addition and the
quaternization; **the esterification that makes B1/B2/B3 is not described
anywhere**, so the structures must be read from the drawings and the NMR.

### 2.3 The structure drawings — decisive

**Figure 1b, page 3, bottom band ("Cholesterol tails"), B1 / B2 / B3.**
Rendered from the PDF at 900–1800 dpi (the figure is a flattened raster, image
xref 82, 1848x2452 px; no vector text layer).

All three are drawn as the standard cholesteryl-ester depiction:
Br-(CH2)n-C(=O)-O-C3, Delta-5,6 alkene, C18 and C19 angular methyls as **solid
bold wedges**, H-9/H-14/H-17 and the C20 methyl as **hashes**, H-8 as a bold
wedge. In all three, **the C3-O bond is a solid bold wedge, narrow at C3 and
widening toward O** — i.e. cofacial with the two bold angular methyls.

The wedge-sense relation is the robust reading, because it survives any flip or
rotation of the whole drawing: *in a flat steroid depiction the 3beta-oxygen is
drawn with the same wedge sense as C18/C19; the 3alpha-oxygen with the opposite
sense.* I confirmed this against RDKit renderings of both epimers of cholesteryl
3-bromopropanoate (`B1_3beta.png`, `B1_3alpha.png` in scratch): RDKit laid both
out in a vertically mirrored orientation, where the 3beta isomer draws C3-O
hashed *together with* hashed C18/C19, and the 3alpha isomer draws C3-O **bold
against hashed** C18/C19. The paper's relation (all three bold) is the 3beta one.

**Figure 2a, page 5, the full C6O2B2 structure.** The sterol is drawn in a
different orientation (ring A upper-left). Same relation: C3-O solid wedge,
C19 methyl solid wedge, H-9 hashed. 3beta again.

While in Fig 2a I also re-read the rest of the lead compound, since the C10
anchor claims to come from it (§6): it draws **five** dodecyl-acrylate arms
(2 on the left tertiary N, 1 on the middle N, 2 on the far-right N), the N+ on a
**piperazine ring nitrogen**, and the cholesteryl **pentanoate** (5-carbon)
linker on that N+. The quaternized ring N is specifically the one whose
exocyclic branch leads to 2 further nitrogens bearing 3 ester arms — and
`_quaternization_site` picks exactly that nitrogen (idx 57; the other ring N,
idx 60, leads to 1 nitrogen / 2 arms). So the **connectivity** of the stored
C6O2B2 matches Fig 2a; only the stereocentre is wrong.

### 2.4 1H NMR — C3-H shift and multiplicity

From the SI figure captions (text layer, exact quotes):

| cmpd | SI fig | olefinic H-6 | **C3-H** | C18-Me |
|---|---|---|---|---|
| B1 | S4, p.5 | 5.38 (dd, J = 5.3, 2.0 Hz, 1H) | **4.71 – 4.63 (m, 1H)** | 0.68 (s, 3H) |
| B2 | S5, p.6 | 5.47 – 5.25 (m, 1H) | **4.66 - 4.58 (m, 1H)** | 0.68 (s, 3H) |
| B3 | S6, p.7 | 5.48 – 5.24 (m, 1H) | **4.61 (m, 1H)** | 0.67 (s, 3H) |
| C6O2B2 | S7, p.8 | 5.30 (d, J = 5.0 Hz, 1H) | 4.54 – 4.51 (m, 1H) | 0.61 (s, 3H) |

For B1 and B2 the reported range is 0.08 ppm wide = **~32 Hz at 400 MHz**. That
is the signature of the *axial* H-3alpha of a **3beta**-acyloxy sterol, which is a
tt with J ~ 11.5, 11.5, 4.5, 4.5 Hz (total width ~32 Hz). The 3alpha epimer puts
H-3 equatorial: a narrow multiplet, W ~ 8-12 Hz (~0.03 ppm), and ~0.4 ppm further
downfield (epicholesteryl acetate H-3 ~ 5.0 ppm vs cholesteryl acetate 4.60).
Both the shift and the width say 3beta.

(The C6O2B2 entry, 4.54–4.51, is narrower; that peak sits on the shoulder of the
huge aliphatic envelope in a 160H-per-molecule spectrum and I would not use it
either way. The three clean B-intermediate spectra are the ones that carry the
argument.)

### 2.5 13C NMR — independently measured from the spectrum, decisive

The SI's typed 13C list for C6O2B2 jumps straight from **122.61 to 64.94** — the
C3-O carbon is simply missing from the transcription. So I measured it off the
plotted spectrum (Figure S7, SI p.8) directly.

Method: rendered SI p.8 at 600 dpi, extracted the column-wise peak-height profile
of the 13C trace, located peak centroids, and calibrated the ppm axis on two
widely separated *assigned* peaks from the authors' own list (173.09 and 122.61),
giving 14.045 px/ppm. Calibration then reproduces the rest of the authors' list
to within ~0.3 ppm: observed 77.5 (CDCl3), **64.94**, 56.75/56.25, 50.34/49.84,
42.45, 39.5/38.3 — against reported 64.94/64.78, 56.69/56.13, 50.03/49.56,
42.30, 39.73/39.51/38.16.

**There is a clean, isolated, intense peak at 73.9 +/- 0.3 ppm** (pixel 2163;
height comparable to the 173 ppm ester carbonyls), sitting just to the right of
the CDCl3 triplet and unlabelled by the authors.

- 3beta-cholesteryl ester C-3: **73.7 – 74.3 ppm**  <- observed
- 3alpha (epicholesteryl) ester C-3: **70.5 – 71.0 ppm** — would fall at pixel
  ~2205. **Nothing is there.**

This is a direct, quantitative measurement from the paper's own raw data, on the
lead compound itself, and it independently confirms the wedges.

### 2.6 Supporting structural reads (confirm the rest of the tail assignment)

- B1: 3.58 (t, J = 6.9, 2H) = CH2Br; 2.89 (t, J = 6.9, 2H) = CH2C(=O) — consistent
  with cholesteryl **3-bromopropanoate**, n_link = 3.
- B2 / B3: CH2Br multiplets at 3.68–3.26 / 3.46 and CH2CO at 2.32 / 2.42–2.18
  (4H, overlapping the sterol C4-H2) — consistent with the 5- and 7-carbon
  bromoalkanoates. Linker lengths as coded.

**Conclusion of §2: certain. The paper's sterol is 3beta-cholesterol. The stored
3alpha structures are wrong. This is the right half to repair.**

---

## 3. CODE REVIEW

All results below are from my own rdkit runs (scratch dir
`/private/tmp/claude-501/-Users-nathanliu-Documents-LNP-confound/ebd1fdd9-2e2e-4284-8eb1-27c940f1a93a/scratchpad/v013/`).

### 3.1 The mechanism is exactly as diagnosed

Minimal reproduction, nothing to do with steroids:

```
start                     : C[C@H](O)CC   -> CC[C@H](C)O     CIP S
RWMol: AddAtom(F), AddBond(C2,F), RemoveAtom(O)
result                    : CC[C@@H](C)F  CIP R      <-- inverted
plain-SMILES truth        : CC[C@H](C)F              CIP S
```

`RWMol.RemoveAtom` deletes the neighbour from the middle of the stereocentre's
bond list and does not re-parity the chiral tag. Confirmed.

### 3.2 `CHOL` itself is fine

`CC(C)CCC[C@@H](C)[C@H]1CC[C@H]2[C@@H]3CC=C4C[C@@H](O)CC[C@]4(C)[C@H]3CC[C@]12C`
-> **HVYWMOMLDIMFJA-DPAQBDIFSA-N**, genuine cholesterol. The input was never the
problem; the surgery is.

### 3.3 Old vs proposed vs an independent ground truth

Ground truth built by *plain string substitution*, never touching `RWMol`:
`CHOL.replace("(O)", "(OC(=O)" + "C"*(n-1) + ")", 1)`.

| n_link | OLD == TRUTH | **NEW == TRUTH** | OLD C3 CIP | **NEW C3 CIP** | OLD cleaved sterol | **NEW cleaved sterol** |
|---|---|---|---|---|---|---|
| 3 | False | **True** | R | **S** | …-VEIPTCAHSA-N (epi) | **…-DPAQBDIFSA-N** |
| 5 | False | **True** | R | **S** | …-VEIPTCAHSA-N (epi) | **…-DPAQBDIFSA-N** |
| 7 | False | **True** | R | **S** | …-VEIPTCAHSA-N (epi) | **…-DPAQBDIFSA-N** |

Not just the C3 label: the **full 8-element CIP label vector** of NEW equals
TRUTH's for all three linkers, and differs from OLD's only at C3. So the
Delta-5 alkene and the other seven stereocentres (C8, C9, C10, C13, C14, C17,
C20) are untouched. Stereo-stripped canonical SMILES of OLD and NEW are
identical — the change really is stereochemistry only.

### 3.4 Index safety — `acyl_c` does not shift

| n | acid SMILES | OH idx | OH's heavy neighbours | acyl_c | after RemoveAtom: acyl_c still a carbonyl C? | indices below OH unchanged? |
|---|---|---|---|---|---|---|
| 3 | `CCC(=O)O` | 4 | `[2]` | 2 | yes | yes |
| 5 | `CCCCC(=O)O` | 6 | `[4]` | 4 | yes | yes |
| 7 | `CCCCCCC(=O)O` | 8 | `[6]` | 6 | yes | yes |

The acid's -OH is always the **last** atom, so removing it cannot renumber
anything below it, and it has exactly one heavy neighbour, so `GetNeighbors()[0]`
is unambiguous. `off = frag.GetNumAtoms()` is correctly read *after* the removal.
Safe.

Also safe but worth stating: the acyl carbon is momentarily under-valent between
`RemoveAtom` and `AddBond`. Because the snippet does not sanitize in between and
sanitizes only at the end, the H count comes out right — formulas are
C30H50O2 / C32H54O2 / C34H58O2, identical to TRUTH. **Do not insert a
`SanitizeMol` between those two lines.**

### 3.5 Regiochemistry — acyl is on the C3 oxygen, not a side-chain methyl

The earlier bug in this function. Verified it stays fixed, with a full-skeleton
SMARTS on the assembled C6O2B2 (not just a substructure of the ester):

```
[N+;R]CCCCC(=O)O[CH;R]1CC[C]2(C)C(=CC[CH]3[CH]4CC[CH]([CH](C)CCCC(C)C)
[C]4(C)CC[CH]32)C1          -> MATCH (True), C3 CIP = S
```

i.e. piperazine N+ -> four CH2 -> carbonyl -> **ring** C3 oxygen -> intact
cholest-5-ene with its isooctyl chain terminating in `CCCC(C)C`. Correct.

### 3.6 `quaternize()` introduces no further stereo damage — but relies on an implicit index

`quaternize` does `InsertMol` + `AddBond` + `SetFormalCharge` + `SetNoImplicit`
and **no** `RemoveAtom`/`RemoveBond`, so it cannot invert anything. Confirmed
end-to-end over the whole 6x2x3 grid:

```
OLD (current repo): 36 lipids, liberated sterol keys = {HVYWMOMLDIMFJA-VEIPTCAHSA-N}
NEW (proposed)    : 36 lipids, liberated sterol keys = {HVYWMOMLDIMFJA-DPAQBDIFSA-N}
```

The radical guard still passes on every one. Also checked the one place stereo
could be lost rather than inverted: `OF-2-CB1/2/3` carry eight STEREOZ alkenes
in the corpus today, so the insertion path preserves bond stereo.

**Latent hazard, flagged not blocking.** `quaternize` does
`off = rw.GetNumAtoms(); rw.InsertMol(frag); rw.AddBond(tidx, off, ...)` — it
bonds the nitrogen to **atom index 0 of the fragment**, which is whatever atom
RDKit's canonicalisation happened to put first in the string
`_cholesteryl_ester` returns. I verified it is currently the chain terminus for
all three linkers:

| n | atom 0 | degree | nH | bonds to carbonyl C | expected |
|---|---|---|---|---|---|
| 3 | C | 1 | 3 | 2 | 2 |
| 5 | C | 1 | 3 | 4 | 4 |
| 7 | C | 1 | 3 | 6 | 6 |

So the proposed repair is safe today. But this is an unstated invariant of
exactly the kind that already shipped a bug in this function. Recommend the
repair also tag the attachment atom explicitly — return the fragment with an
atom map (`[CH3:1]...`) or a dummy, and have `quaternize` look it up — rather
than trusting canonical ordering.

### 3.7 Minor

- The proposed snippet replaces `next((...), None)` + `if oh is None: return None, None`
  with a bare `next(...)`, which raises `StopIteration` instead of returning the
  sentinel. `CHOL` is a module constant so it cannot fire, but keep the guard.
- Keep the trailing `Chem.SanitizeMol` try/except and the `(smiles, 0)` return
  shape; `quaternize` unpacks two values.

**Confidence: certain** (everything in §3 is a reproduced measurement, not a reading).

---

## 4. SECOND COPIES

Repo-wide grep for `RemoveAtom`, `RWMol`, `InsertMol`, `RemoveBond` and for any
hard-coded steroid SMILES (`CC=C4C`).

### 4.1 The live second copy — must be fixed in the same commit

**`/Users/nathanliu/LNP/confound/curation/build_new_curation.py:164-218`**, function
`quaternize(core_smiles, n_link)`:

- **line 177-178** — its own hard-coded `chol = "CC(C)CCC[C@@H](C)..."` (the
  string is correct 3beta cholesterol).
- **line 179** — `rwf = Chem.RWMol(Chem.MolFromSmiles("C" * (n_link - 1) + "C(=O)O"))`
- **line 188** — `rwf.InsertMol(cm)`
- **line 189** — `rwf.AddBond(ester_o, off0 + c3, Chem.BondType.SINGLE)`
- **line 190** — `rwf.RemoveAtom(off0 + oh.GetIdx())`   <- **the identical defect**

Its comment block (lines 174-176) even recites the *previous* bug in this
function as if it were the current state of the art. The inversion is not
mentioned.

**Blast radius.** This function produces, via the `COMMERCIAL_CORES` /
`PARENT` branch at lines 303-336, the 15 commercial-core conjugates:

```
MC3-CB1/2/3  SM-102-CB1/2/3  OF-2-CB1/2/3  TT3-CB1/2/3  Lipid 14-CB1/2/3
```

Corpus split, measured:

| builder | unique structures | corpus rows |
|---|---|---|
| `siegwart_42634137/enumerate.py` | 36 | 522 (486 vitro + 36 vivo) |
| `build_new_curation.py` (this copy) | 15 | 15 (vivo) |
| **total** | **51** | **537** |

Independent confirmation over the live corpus: I cleaved the ring-carbon ester of
every steroid-bearing SMILES in `data/data_vitro.csv` + `data/data_vivo.csv` and
assigned CIP at C3.

```
('CHL_SMILES', PMID 18438401, C3 = S)  30064 rows   <- free cholesterol, correct
('CHL_SMILES', PMID 32080183, C3 = S)     13 rows   <- free cholesterol, correct
('IL_SMILES',  PMID 42634137, C3 = R)    537 rows   <- ZT_2026, ALL WRONG
```

All 51/51 ZT_2026 unique structures liberate `HVYWMOMLDIMFJA-VEIPTCAHSA-N`.
**Nowhere else in the corpus is any C3 = R.** The two builders are the only
sources and both must change.

### 4.2 Hard-coded sterol SMILES — all correct, none 3alpha

| file:line | context | verdict |
|---|---|---|
| `curation/siegwart_42634137/enumerate.py:57` | `CHOL` constant | 3beta, correct |
| `curation/build_new_curation.py:177` | local `chol` in `quaternize` | 3beta, correct |
| `curation/ng_2025/stage07_emit.py:29` | `CHOL` constant | 3beta string, identical to the above (grep only — that folder is owned by another session and I did not read or run it) |
| `provenance/reagent_registry.json:62,65` | `Cholesterol` registry entry + variant | 3beta, correct; InChIKey pin `HVYWMOMLDIMFJA-DPAQBDIFSA-N` |
| `provenance/reagent_registry.py:64` | registry source | 3beta, correct |

No 3alpha string is hard-coded anywhere. Every stored epimer in the repo is a
*product* of the two `RemoveAtom` sites, not a typo.

### 4.3 Other `RemoveAtom` / `RemoveBond` sites — latent, not currently firing

| file:line | call | assessment |
|---|---|---|
| `curation/ql_2026/enumerate_library.py:112` | `m.RemoveAtom(dummy)` after `AddBond(n_idx, anchor)` | **Same pattern.** If `anchor` were ever a stereocentre the tag would invert. `09_enumerated_344.csv`: 344 rows, **0** `@`, 0 `/C=C` — no stereo to damage today. Latent. |
| `curation/lin_aromatic_2026/enumerate_library.py:98` | `frag.RemoveBond(o, ch2)` (epoxide opening) | Same class; `03_enumerated.csv` (840) and `13_enumerated_147.csv` (147): 0 `@`, 0 `/C=C`. Latent. Its `add_acrylate` also demotes a C=C to single, which would strand a stale bond-stereo directive — the three `/C=C\` strings in `blocks.csv` are commercial helper lipids (DOPE/DOPS/DOPC), not IL blocks, so nothing is affected. |
| `curation/liu_iphos/enumerate_library.py:133-138` | `RWMol` + `InsertMol` only, no removal | safe; `03_enumerated.csv` 572 rows, 0 `@`. |
| `curation/ng_2025/stage08_reference.py:72` | `m.RemoveAtom(off + star[0])` | **Flagged for that session, not reviewed here** — out of scope by instruction. The grep line shows a dummy-atom removal immediately after fragment assembly, i.e. the same shape as `ql_2026:112`. Worth their checking, especially as `ng_2025/stage07_emit.py:29` carries a cholesterol constant. |

**Confidence: certain** for §4.1-4.2 (measured over the live corpus), **probable**
for the ng_2025 line (grep only, file not opened).

---

## 5. THE PERMANENT CHECK

Proposed: after regenerating, cleave the steroid ester of every `ZT_2026`
`IL_SMILES` and assert the liberated alcohol's InChIKey ==
`HVYWMOMLDIMFJA-DPAQBDIFSA-N`, taken from `provenance/reagent_registry.json`
under `Cholesterol`.

### Would it have caught this on day one? Yes.

It tests precisely the invariant that was violated, at the layer where it was
violated, against a value that is pinned independently of the generator. I ran
exactly this check as written and it separates OLD from NEW cleanly for all 51
structures and both builders. It is the right check.

### Why nothing existing caught it (this is the useful part)

- **`reagent_registry.py` cannot, by construction.** Two reasons, both structural:
  (a) it compares **InChIKey *skeleton* only** — line 101, `…split("-")[0]`, with
  the docstring rationale at line 20, *"so stereo descriptors do not cause
  noise."* An epimer shares its skeleton block exactly, so the registry is
  blind to this entire error class by design. (b) It only inspects the
  `IL/HL/CHL/PEG/fifthcomponent_name` + `_SMILES` pairs; **all 537 ZT_2026 rows
  have an empty `CHL` column** (verified — the paper's whole claim is
  cholesterol-free LNPs), so the sterol here lives *inside* `IL_SMILES` where the
  registry never looks. The proposed check fixes both gaps at once. It must not
  inherit the skeleton-only habit: **compare the full InChIKey.**
- **`structure_check.py` cannot**, because its only anchor for this paper is a
  byte-copy of the thing it checks (§6) — *and* because even a correct anchor
  would only have raised `stereo_differs_connectivity_ok`, a non-failing flag
  (lines 128-132, 157). The stereo class is currently un-gated repo-wide.
- **`audit.py` C9** checks only *unassigned C=C* bond stereo
  (`FindPotentialStereoBonds`), never tetrahedral centres.
- Mass / formula / descriptors are blind, as the proposal already says.

### How to break it

1. **Vacuous pass.** If the row selector finds nothing — paper key renamed,
   `Publication_PMID` reformatted, rows temporarily dropped — the `all(...)`
   is trivially true and the check reports green. *Fix: assert the counts.
   `537` rows, `51` unique `IL_SMILES`, `>= 1` cleavage site per structure.
   A structure with zero ester cleavage sites must FAIL, not be skipped.*
2. **Empty-key pass.** `Chem.MolToInchiKey` returns `''` on InChI failure, which
   does happen on large permanently-charged species. An implementation shaped
   `if key and key != EXPECTED: fail` passes silently. *Fix: `''`/`None` is a
   failure.*
3. **Data-only "repair".** Someone string-edits the CSV instead of the
   generators; the check goes green while both `RemoveAtom` sites stay broken,
   and the next regeneration re-introduces 537 wrong rows. *Fix: pair it with a
   regenerate-and-diff determinism test, so corpus and generator cannot drift.*
4. **Scope.** As written it is ZT_2026-specific. The same bug in the next paper
   with a sterol is uncovered. *Fix: make it corpus-wide and column-agnostic —
   the query I ran in §4.1 takes seconds over 30,614 rows: any molecule
   containing the cholest-5-ene skeleton, in any `*_SMILES` column, must present
   C3 as `S` (free sterol or ester alike).*

### What it misses even when hardened

It validates **one fragment**. It says nothing about the amine head, the number
of acrylate arms, the quaternization site, the linker length, or the alkyl
identity — the other four things this curation had to infer. It is a narrow,
cheap, permanent guard on one reagent, and should be described as that rather
than as structural validation.

One genuine strength worth keeping: because it runs on `data/*.csv` rather than
on `03_enumerated.csv`, it covers **both** builders, including the
`build_new_curation.py` copy that the repair proposal overlooks. Run it there.

**Recommended additional check, at the function level.** The corpus check is
downstream; add a three-line unit test that is upstream and is a true oracle:

```python
assert _cholesteryl_ester(n)[0] == Chem.MolToSmiles(Chem.MolFromSmiles(
    CHOL.replace("(O)", "(OC(=O)" + "C"*(n-1) + ")", 1)))   # n in (3, 5, 7)
```

This is the comparison I used in §3.3. It is independent of RDKit's molecule
surgery entirely, it pins regiochemistry and all eight stereocentres at once,
and it fails at the function rather than 537 rows later.

**Confidence: certain** on the diagnosis of why the existing checks missed it
(read from the code); **probable** that the four break modes are the main ones —
I constructed them by inspection, not by fuzzing an implementation that does not
exist yet.

---

## 6. THE C10 ANCHOR

### The finding is confirmed, and it is worse than "compared to itself"

```
reference_structures.csv  reference_SMILES   sha256 a08277aa15758b07…  (223 bytes)
03_enumerated.csv         C6O2B2 IL_SMILES   sha256 a08277aa15758b07…  (223 bytes)
byte-identical: True
```

And the note attached to it reads:

> "Connectivity read off the drawn structure and assembled from fragments,
> **independently of enumerate.py**: … cholesteryl 5-carbon ester **through the
> C3 oxygen (3-beta)** … Cholesterol stereo descriptors are the standard ones,
> not re-read from the figure."

The string it annotates is 3**alpha**. So the note does not merely fail to be
checkable — it asserts, in prose, the opposite of what the stored string says,
and it was the single anchor relied on for this paper. `structure_check.py`'s
docstring (lines 44-51) anticipates the copy failure and concludes it "is NOT
checked automatically and cannot be." That is right in general and wrong in this
instance; see the guard below.

Separately: even a correct 3beta anchor would not have **failed** this. The
epimer changes only the second InChIKey block, which `structure_check.py` routes
to `stereo_differs_connectivity_ok` — reported, never failing (lines 128-132).
Fixing the anchor without fixing that policy leaves the defect visible but
un-gated.

### What an honestly independent anchor should be, and where to read it

**Read the steroid half from Figure 1b, page 3 — not from Figure 2a.** Fig 1b's
"Cholesterol tails" band draws B1/B2/B3 standalone at roughly four times the
scale of the sterol in Fig 2a, with every wedge and hash resolvable at 900 dpi
(I read all three; §2.3). Fig 2a (page 5) remains the right source for the
**amine half** — the five arms, the piperazine N+, and which ring nitrogen bears
the tail — because it is the only place the whole lead compound is drawn.

So: `source_figure` should read **"Fig 1b p.3 (steroid half) + Fig 2a p.5 (amine
half)"**, not "Fig 2a" alone. Two figures, two halves, both re-readable.

**Derivation.** Assemble the reference from two sources, neither of which is
`03_enumerated.csv`:

1. **sterol half** — take `registry_smiles` verbatim from
   `provenance/reagent_registry.json` -> `Cholesterol`
   (`…[C@@H](O)…`, pinned `HVYWMOMLDIMFJA-DPAQBDIFSA-N`). This makes the note's
   claim — "stereo descriptors are the standard ones, not re-read from the
   figure" — *true and mechanically enforced* instead of an unverified
   assertion beside a 3alpha string.
2. **amine half** — type the quaternary-ammonium polyamine with its five
   dodecyl-acrylate arms by hand from Fig 2a.
3. **join** by plain SMILES substitution, no `RWMol` anywhere:
   `CHOL.replace("(O)", "(OC(=O)CCCC[N+]1(…)…)", 1)`.

That is a genuinely different code path from the enumerator (which does fragment
surgery), it is the same oracle I used in §3.3, and it would have failed OLD
loudly.

**Record the observables, not just the conclusion**, in the `note` column, so the
next reviewer can re-check the reading without re-deriving it:
C3-O bold wedge cofacial with the bold C18/C19 methyls in Fig 1b and Fig 2a;
1H 4.66–4.58 (m, ~32 Hz, axial H-3) for B2; 13C 73.9 ppm measured off Fig S7.

**Add one automatic guard** to `structure_check.py`. The docstring is right that
an honestly re-derived reference canonicalises to the same string as the
enumerated one, so a post-canonicalisation comparison proves nothing. But the
**raw stored strings** are comparable, and that is cheap and informative: a
human typing a 223-character lipid from a figure essentially never reproduces
RDKit's canonical output character-for-character, including ring-closure
numbering and atom order. So:

> if `reference_SMILES` as stored is byte-identical to `IL_SMILES` as stored
> **and** is already in RDKit canonical form, emit `ANCHOR_MAY_BE_A_COPY`.

It is a heuristic and will occasionally be wrong (someone may legitimately paste
a canonicalised re-derivation), which is why it should be a loud status rather
than a hard failure — but it fires on today's file and would have fired on day
one.

Finally, **make reagent-derived stereo a hard gate**: keep the soft
`stereo_differs_connectivity_ok` policy for hand-typed reference stereo, and put
the sterol's stereochemistry under the §5 registry check instead, where the
expected value comes from a pin rather than from someone's typing. That is the
right division: hand-typed stereo stays advisory; registry-pinned stereo fails
the build.

**Confidence: certain** on the byte-identity and on the stereo-soft-flag policy
(both read directly from the files); **probable** that the two-figure split is
the best anchor design — it is my judgement after reading both figures, not a
measurement.

---

## 7. CONFIDENCE SUMMARY

| § | Claim | Confidence |
|---|---|---|
| 1 | Right half; repair correct; **incomplete — second builder untouched** | **certain** |
| 2.1 | No epi / 3alpha / allo sterol anywhere in paper or SI | **certain** (exhaustive text search of both PDFs) |
| 2.2 | Sigma-Aldrich "cholesterol", no CAS given; B1-B3 synthesis not described | **certain** (quoted); the inference Sigma "cholesterol" = 57-88-5 is **probable** |
| 2.3 | Fig 1b (p.3) and Fig 2a (p.5) both draw C3-O cofacial with C18/C19 = 3beta | **certain** (read at 900-2400 dpi, cross-checked against rendered reference epimers) |
| 2.3 | Fig 2a draws 5 arms, N+ on the piperazine N matching `_quaternization_site` | **probable** (read from a raster figure; the arm count is unambiguous, the ring-N assignment I traced by branch topology) |
| 2.4 | 1H C3-H at 4.58-4.71, ~32 Hz wide = axial H-3 of a 3beta ester | **certain** on the quoted values; **probable** on the ~0.4 ppm epimer separation (literature recall, not re-measured) |
| 2.5 | 13C peak at 73.9 +/- 0.3 ppm = 3beta ester C-3; nothing at ~70.9 | **probable→certain**: the peak and the calibration are my own measurement and reproduce 8 of the authors' assigned shifts to ~0.3 ppm; the 73.7/70.6 epimer reference values are literature recall |
| 3 | Proposed code == independent plain-SMILES truth, all 8 stereocentres, n=3/5/7; regiochemistry on C3-O; `quaternize` adds no damage | **certain** (reproduced) |
| 3.6 | `quaternize`'s atom-index-0 assumption holds today but is implicit | **certain** that it holds; **probable** that it is worth changing |
| 4.1 | `build_new_curation.py:179-190` is a live second copy; 15 structures / 15 rows; 537/537 corpus rows carry C3 = R; no other C3 = R in the corpus | **certain** (measured over `data/*.csv`) |
| 4.3 | ql_2026 / lin_aromatic sites are the same pattern but have no stereo today | **certain**; ng_2025 line **probable** (grep only, by instruction) |
| 5 | Check would have caught it; registry is blind by design (skeleton-only + empty CHL column); four break modes | **certain** on why existing checks missed it; **probable** that the break list is complete |
| 6 | Anchor is a byte-copy with a note that contradicts its own string; stereo would have been soft-flagged anyway | **certain** (sha256 + code read) |
| 6 | Recommended anchor design | **probable** (judgement) |
