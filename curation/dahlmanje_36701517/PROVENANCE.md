# PROVENANCE — dahlmanje_36701517

Radmand et al. 2023, PMID 36701517. **Licence: CC BY** (Europe PMC core record,
`license: "cc by"`). Nothing outside this folder was written.

## Anchor posture — recorded by hand

`provenance/structure_check.py` (C10), `check_anchor_independence.py` (C16) and
`check_anchor_coverage.py` (C18) were **not run**: they write tracked files that
other agents read. What they would report is stated here instead, with the
numbers that back it.

### What the paper draws, and what it does not

| lipid | drawn? | where | `reference_SMILES` |
|---|---|---|---|
| **cKK-E12** (the only ionizable lipid, all 144 formulations) | **NO — nowhere** | named once in the SI Materials with a vendor and catalogue number (Organix Inc., O-8744) | **EMPTY, on purpose** |
| 18:0 DDAB (Cat-LNP helper) | yes, **twice** | Fig S2a row C3; redrawn in Fig 3a | built |
| 18:0 PC / DSPC (Neu-LNP helper) | yes, **twice** | Fig S2b row N7; redrawn in Fig S15a | built |
| 18:1 PA (An-LNP helper) | yes, **twice** | Fig S2c row A1; redrawn in Fig S15c | built |

C18 would record this folder as having **one anchor gap**, and the gap names
work that cannot be done: there is no drawn structure for the ionizable lipid to
anchor against. The empty row is in `reference_structures.csv` so the absence is
visible rather than silently filled.

### Independence — by provenance, executed, not asserted

`make_reference.py` walks **its own AST at run time** before building anything
(`_ast_guard`). It refuses if the module imports `assemble`, `reagents`,
`corpus_components`, `reagent_registry` or `structure_check`; if it calls
`RWMol`, `AddBond`, `RemoveAtom`, `CombineMols`, `EditableMol`,
`ReplaceSubstructs`, `esterify` or `attach`; or if any **executable** string
literal names a registry or corpus file. Docstrings and the guard's own
blocklist are exempt, which is why the guard can describe what it forbids.

Clean run: `forbidden_calls_found 0`, `forbidden_imports_found 0`,
**1776 AST nodes walked**, exit 0.

**The guard was made to fail on purpose.** `--seed-guard` re-parses the module
with `import assemble` prepended; the guard fires, naming that import, and exits
**1**. Byte identity of the emitted anchor is deliberately **not** tested —
F-037 established that a correct re-derivation is byte-identical to a correct
enumeration once canonicalised, which is exactly why C10's own
`anchor_may_be_a_copy` missed it.

### The gate tests connectivity — seeded with an isomeric defect

`--seed-ester` writes the DSPC sn-1 acyl fragment in its other direction
(`C(=O)CCC…` where `CCC…C(=O)` belongs). The result:

| | clean | seeded |
|---|---|---|
| formula | C44H88NO8P | **C44H88NO8P** (same) |
| exact mass | 789.624755 | **789.624755** (Δ 4.06e-7 Da, floating point only) |
| parses, canonicalises, 0 radicals | yes | yes |
| skeleton InChIKey | `NRJAVPSFFCBXDT` | `YXPWLPQHOHNCMY` (**different**) |
| connectivity (esters, ketones, dialkyl ethers, quaternary N) | (2, 0, 0, 1) | **(1, 0, 1, 1)** |
| verdict | pass | **REFUSED, exit 1** |

This is not a hypothetical seed. **The first version of this file made exactly
that mistake by accident**, reported `C44H88NO8P` at `789.624755` — the right
formula and the right mass for DSPC — and was a ketone bearing an ether. It was
caught only by comparing InChIKeys. `_direction_checks` is the permanent guard:
*a SMILES fragment has a direction that mass cannot see*.

Verdicts above are read from `returncode`, never from stdout and never through
a pipe (`python3 … | tail` returns tail's exit code — F-012's shape).

### Why string composition is safe here specifically

`CLAUDE.md` mandates fragment assembly because string concatenation once fused
piperazine to cholesterol through a **ring-closure digit collision**. All three
molecules built here are **acyclic** — `CalcNumRings == 0` is asserted on each —
so no ring-closure digit exists in any fragment and that collision cannot arise.
The opposite hazard (`RWMol.RemoveAtom` dropping an atom out of a stereocentre's
neighbour list, F-013) also cannot arise, because no atom is ever removed.

Stereo assertions: 18:0 DDAB has **0 stereocentres**, so its stereo assertion is
**VACUOUS, not passing**. DSPC and 18:1 PA each have 2 assigned centres and 0
unassigned. 18:1 PA's two C=C bonds are asserted explicitly as `STEREOCIS`
(cis-9, as drawn) — CIP centre counting is blind to double-bond geometry, so it
is checked separately.

### C10-equivalent comparison, made locally

| anchor | registry entry | skeleton | full InChIKey | note |
|---|---|---|---|---|
| 18:0 DDAB | `DDAB` | match | **match** | — |
| 18:0 PC (DSPC) | `DSPC` | match | **match** | — |
| 18:1 PA | *(none)* | — | — | the registry's only PA is `14PA` = DMPA, 14:0/14:0 — a different molecule. **An-LNP's helper lipid is the one excipient in this paper with no second attestation.** |
| cKK-E12 | *(none)* | — | — | empty anchor; the paper draws nothing |

## Tier posture (full detail in `00_value_sources.json`)

Everything numeric in this paper is **tier 5, pixels**, because tiers 1–4 do not
exist: the publisher deposited exactly one file (a PDF, confirmed against the
Europe PMC `supplementaryFiles` zip, which is the publisher's own deposit), that
PDF contains no table object and no text-layer numerals inside any figure, no
repository or accession is named anywhere, and `find_repo.py` found none.

Native resolution is **measured**, not assumed: every figure page is a flattened
raster at **219.7–220.2 dpi** (`px / rect_pt × 72`) with **0 vector paths**, the
one exception being the three red lead-LNP boxes of Fig S6, which are genuine
vector strokes of colour exactly (1.0, 0.0, 0.0) at 1.25 pt.

`axis.native_tol(1.5, native_dpi=220)` = **4.09 px** was available and **was not
used**. Every panel is gated on the **unrescaled 1.5 px**, and the worst
leave-one-out tick residual across all 21 sub-panels is **0.282 px**. Reporting
a pass against the tighter gate when the looser one was available is the
stronger statement; nothing was widened.

Main-text figure rasters (PMC's 650 × 935 px JPEGs) were used for **labels and
printed numerals only, never for a digitised value** — they are standalone JPEGs
with no page rect, so their native dpi cannot be measured by the px/pt rule and
`axis.native_tol` refuses a guessed dpi. Nothing is lost: every value in
main-text Fig 3b is also in Fig S13a, raw and with its PBS control, at a
measured 220 dpi.
