# PROVENANCE — `dahlmanje_33758781`

Gan Z, Lokugamage MP, Hatit MZC, Loughrey D, Paunovska K, Sato M, Cristian A,
**Dahlman JE**. *Nanoparticles containing constrained phospholipids deliver mRNA
to liver immune cells in vivo without targeting ligands.*
**Bioeng Transl Med 2020;5(3):e10161. PMID 33758781 · PMCID PMC7974401 ·
DOI 10.1002/btm2.10161.**

**First author is Gan Z**, not Dahlman — the folder is named after the senior
author. Prospective `Experiment_ID` would be **`ZG_2020`**; no such id exists in
either corpus file (43 vitro + 20 vivo ids enumerated 2026-11-10).

**Outcome: ZERO rows emitted.** `06_rows.csv` is header-only (79 columns, read
live from `data/data_vivo.csv`) so the zero is visible. This folder must **not**
be added to `build_dataset.py`'s `MERGE`, and `DECISIONS.md` carries five `OPEN`
entries that would block it anyway.

---

## Licence

| | |
|---|---|
| Article | **CC BY 4.0** — Europe PMC `core` record `license` field reads `cc by`; the JATS `<permissions>` block cites `creativecommons.org/licenses/by/4.0/`. §8-compatible. |
| Repository / data | **None exists.** No deposit, no code release. |

Read off the Europe PMC record, not the publisher landing page.

---

## Acquisition, and the F-065 check

Nothing was copied out of the session scratchpad. Every source file was fetched
directly into `source/`, or extracted in-process from a ZIP fetched there, with
a `^[A-Za-z0-9._-]+$` name guard confining writes to one directory.

| file | origin | bytes | sha256 (head) |
|---|---|---|---|
| `source/btm210161-sup-0001-AppendixS1.pdf` | **already in the folder at task start** | 3,324,917 | `24d89f54d064e46b…` |
| `source/pmc/PMC7974401.jats.xml` | `provenance/tools/pmc_fetch.py PMC7974401` (efetch) | 89,915 | — |
| `source/epmc_suppl/suppl.zip` | `ebi.ac.uk/europepmc/webservices/rest/PMC7974401/supplementaryFiles` | 3,257,992 | `31afaacbc877af4b…` |
| `source/epmc_suppl/extracted/*` | unzipped from the above, 19 entries | — | each hashed in `07_qc.json` |
| `source/article_epmc.pdf` | `europepmc.org/articles/PMC7974401?pdf=render` | 5,626 | **not a PDF** — Cloudflare interstitial, kept unused, **not bypassed** |

**The pre-existing Appendix S1 was verified, not trusted.** Its sha256 is
byte-identical to `BTM2-5-e10161-s001.pdf` inside the Europe PMC
`supplementaryFiles` ZIP for PMC7974401, which we downloaded independently. That
proves the file belongs to *this* PMCID — the check a file-count or size
comparison would have passed. `pmc_fetch.py`'s own binary downloads all returned
~21.5 kB HTML error pages and were correctly refused by its magic-byte
validation; only the XML landed.

---

## Anchor posture — recorded BY HAND

`provenance/structure_check.py`, `provenance/check_anchor_independence.py` and
`provenance/check_anchor_coverage.py` were **not run**: they write tracked files
other agents read. What follows is the posture they would record.

`reference_structures.csv` has the exact header
`lipid_name,reference_SMILES,source_figure,read_on,read_by,expected_mismatch,note`
and **7 rows**:

| lipid | `reference_SMILES` | posture |
|---|---|---|
| `A-9` | present | **DERIVED.** Independent generator, InChIKey-identical to the enumeration. |
| `A-17` | present | **DERIVED.** Same, deliberately a different chain length (C18 acyl vs A-9's C10). |
| `Cholesterol` | present | **REGISTRY PIN, NOT INDEPENDENT — declared as such on the row.** Taken verbatim from `provenance/reagent_registry.py` so the provenance claim is mechanically true. Fig 2a's drawn skeleton agrees (four fused rings, 3-OH, 5-ene, iso-octyl tail); its **eight stereocentres are not resolvable at 709 px** and are inherited from the registry, so they are *not* evidence from this paper. |
| `cKK-E15` | **empty** | Drawn in Fig 2a but the tail label's subscripts are 2–3 px tall; see `DECISIONS.md` D-2. |
| `beta-Sitosterol` | **empty** | The C24 ethyl creates a stereocentre a 709 px drawing cannot resolve. Not in the registry. **Discarded rather than hand-written** — hand-writing a sterol is how a sibling produced a 4-of-8-inverted diastereomer today with identical formula, identical exact mass and an identical skeleton InChIKey block. |
| `20alpha-OH Cholesterol` | **empty** | The 20-hydroxyl's configuration is the name's whole content and is not resolvable. Not in the registry. |
| `C14PEG2000` | **empty** | **Not anchorable in principle** — the PEG block length is a distribution, which is why the registry deliberately holds no PEG-lipid. |

### How independence is actually established

Byte-identity of the emitted anchor is **not** the signal, and testing for it is
a trap (F-037): a correct re-derivation canonicalises to the same string as a
correct enumeration. Both our anchors *are* byte-identical to the enumeration,
and that fact carries no information either way. Independence is tested on the
**generator**, by two guards that run every time `make_reference.py` runs:

- **Guard 1 (AST, at run time).** Walks `make_reference.py`'s own syntax tree
  and refuses any RDKit editing primitive (`RWMol`, `AddBond`, `RemoveAtom`,
  `CombineMols`, `RunReactants`, …) or an `assemble` import. It inspects
  `Call`/`Import` nodes **only, never string constants**, so the ban list
  cannot match itself — the failure mode a sibling's substring guard hit today.
- **Guard 2 (template provenance).** No string literal of ≥25 characters in the
  module may be a substring of any enumerated SMILES, or vice versa. This is
  F-037's defect stated directly, and it is the part that byte-comparison of
  the *output* structurally cannot see.

The generator differs from the enumerator in route as well as in text: the
enumerator starts from intermediate 4 and calls `assemble.esterify`
(combine → `AddBond` → `RemoveAtom`), while the anchor is written as plain text
in the **opposite traversal** — adamantyl cage first, then acetyl ester,
glycerol, fatty acyl, phosphocholine last — off Figure 1a/1b.

The stereocentre is **selected, not guessed**: `_pick_R` builds both tags,
labels each with RDKit's modern CIP labeller, and returns the one labelled `R`,
`R` being the descriptor the paper prints for intermediate 4 and all ten
products. The ranking `OH > CH2-O-P > CH2-O-acyl > H` is derived in the
enumerator's own comments and is identical at all three stages of the route, so
no CIP relabel can occur along it. The adamantyl cage's printed `(3R,5R,7R)`
descriptors are **vacuous** — bridgeheads are not stereogenic, RDKit finds one
centre per molecule, and the paper itself writes them lower-case in some
entries.

### The guards were made to fail on purpose

`python3 -I make_reference.py --selftest` — three mutated copies, each run in
its **own child process**, verdict from `returncode`, and **credit refused to
any seed whose metric line is absent**:

| seed | what it injects | result |
|---|---|---|
| `G1_banned_call_injected` | a `Chem.RWMol(...)` call inside `build()` | **FIRED** — `guard_ast_violations=1 [('call','RWMol',187)]` |
| `G2_template_copied_from_enumeration` | replaces the template with the enumeration's own canonical SMILES for A-9 — **F-037's exact defect** | **FIRED** — `guard_template_provenance_violations=9` |
| `G3_regioisomer_identical_mass` | swaps the two acyl groups on the glycerol: **identical formula, exact mass identical to 0.000000 Da** | **FIRED** — `inchikey_agreements=0/2` |

`G3` is the one that matters: it proves the anchor gate tests **connectivity**
and not composition. On the first run all three exited non-zero *before reaching
a metric* (the copies could not find `03_enumerated.json`), and the harness
correctly refused to credit them — **the harness was fixed, not the seeds**.

---

## What was checked, with numbers

Full detail in `07_qc.json`; `stage07_qc.py` reads every figure out of the file
that computed it rather than transcribing it.

| check | result |
|---|---|
| HRMS exact mass, our build vs printed `calcd` | **10/10**, max abs delta **0.04 mDa** |
| HRMS adduct formula | **9/10**; A-16's printed formula omits its own `[M+Na]+` sodium (`X-3`) |
| Adduct conventions in one column | `[M+Na]+` ×9, `[M+H]+` ×1 — both reproduce (trap 2) |
| Fig 1c neutral formulas (a second, independently typeset source) | **10/10 exact** |
| Fig 1c molecular weight | median printed offset −0.435 Da; max residual after offset **0.0** Da once A-13 is excused; A-13's own column steps are 13.0 and 15.0 against CH₂ = 14.027, summing to exactly 28.0 — a one-digit typo (`X-4`) |
| Fig 1c LogP monotone in chain length | **8/8** increments positive; A-17-2Z < A-17 as two *cis* bonds require |
| Name token vs chain length | **10/10** — `A-N` is the acyl chain minus one carbon (A-9 = decanoyl C10, A-17 = stearoyl C18) |
| Fig 2b molar ratios sum to 100 | **4/4**, 0 offenders |
| Factorial 10 × 3 × 4 | **= 120**, the number the paper states twice; 120 − 109 = 11 dropped at QC, never identified |
| Fig 3h enriched-bar count by colour identity | **5 teal**, matching the caption's printed "five"; 5 + 4 negative + 1 ≈ zero = **10** |
| Fig 3h rank | tallest teal bar is **A-11**, matching the text's "most highly enriched" |
| **Printed yield range 77–92%** | reproduces **only with A-9 included** (without it: 77–89); 92% is A-9's yield and no other compound's — the authors' own number settles their "A10-17-2Z" off-by-one (`X-2`) |
| Seeded mutations | **14 seeds + 3 controls, 0 uncredited**; 13 fired, 1 blind as predicted |
| Anchor guards | **3 seeds, 0 uncredited**, all fired |

**Measured blind spot.** Seed `S5` swaps the adamantylacetyl and fatty acyl on
all ten compounds and the HRMS gate stays **green**: formula identical, exact
mass identical to **0.000000 Da**. That is why `X-1`'s regiochemistry is settled
from the synthesis route and the uniform ¹H NMR pattern, and never from mass.

**Measured sensitivity limit.** Seed `S10` shifts every Fig 3h bar slot by half
the 22.8 px pitch and breaks both the sign count and the rank — so the slot
positions are load-bearing at that scale and are *detected*, not assumed.

---

## What was deliberately NOT done

- **No digitiser was run, anywhere.** The only obtainable figures are Europe PMC
  web JPEGs at 709 px on the long edge with **no page rect**, so native dpi is
  not derivable and `axis.native_tol()` would refuse them. Every pixel-derived
  quantity in this folder is a **count, a sum, a sign or a rank** — never a
  value. No `axis.calibrate`, no `bars.bar_top`, no `heatmap` read.
- **No browser.** The publisher's full-resolution figures (Wiley) need one, and
  Europe PMC's `?pdf=render` returns a Cloudflare interstitial. Not bypassed.
- **No tracked file outside this folder was written or read-modified.**
- **`provenance/structure_check.py`, `check_anchor_independence.py` and
  `check_anchor_coverage.py` were not run** — posture recorded above by hand.

## What could not be verified

1. **`cKK-E15`'s tail length** (`D-2`). The paper draws it, cites nothing for
   it, and abbreviates the tails to an illegible `CnH2n+1`. Left as an
   inference, with an **empty** SMILES everywhere it would otherwise be stored.
2. **The `Normalized Delivery` formula** (`X-5`). The Results and the Methods
   give two *different* normalisations and the paper publishes no worked
   example and no counts, so **no reconstruction residual can be reported**.
   The family is decidable from the axes: both centre on 1, where a plain
   composition over 109 LNPs would centre on 1/109 = 0.00917, ~94–109× lower.
3. **LNP1's sterol and molar ratio** (`X-6`). Never stated anywhere.
4. **Which 11 of 120 formulations failed QC.** Shown only as a change in a pie
   chart, so the 109 members cannot be enumerated.
5. **Supplementary Figures S1–S3.** Cited in the main text, absent from
   Appendix S1 and from every copy obtainable.
6. **Cholesterol's eight stereocentres, and both other sterols** — see the
   anchor table.

---

## One defect found in this folder's own tooling

A determinism check — run every stage from clean twice and compare hashes —
showed three outputs differing between runs. Cause: **the seed harnesses pinned
the mutated copies' `HERE` to the real folder** (so a copy in a temp directory
could still find its inputs), and the copies therefore wrote their *outputs*
there too. `reference_structures.csv` was left holding seed **G3's regioisomer
anchors** — structurally wrong molecules, in the file whose entire job is to be
the external ground truth, written by a run that had exited non-zero.

Fixed by giving every seeded stage an `OUTDIR` distinct from `HERE`; the
harnesses now redirect `OUTDIR` to their temp directory while leaving `HERE`
pinned for reads. Re-verified: the 14 + 3 seed verdicts are unchanged, and ten
outputs are byte-identical across three independent runs **including runs
performed after both harnesses execute**.

Worth recording as a class: *a mutation harness that shares an output path with
the thing it mutates will leave the seeded defect behind in a tracked file, and
the seed's own non-zero exit is what hides it* — the run looks like a correct
refusal. Only a cross-run hash comparison sees it. A file-count or
"did the script exit 0" check would have passed.
