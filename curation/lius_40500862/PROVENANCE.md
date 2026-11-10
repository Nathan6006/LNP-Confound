# PROVENANCE — lius_40500862

Wu S.‡, Shi L.‡, Su K., Lin L., Yan X., Zhang Y., Gu T., Wang Z., Xu T.,
**Liu S.*** — *Carbonate-Bearing Ionizable Lipids for mRNA Delivery to Splenic
NK Cells*, J. Am. Chem. Soc. **2025**, 147, 28665–28673.
DOI `10.1021/jacs.5c01389` · PMID `40500862` · Zhejiang University.

**The tracker title is truncated.** It reads "…for mRNA Delivery to Spleen…";
the published title ends "**to Splenic NK Cells**". The endpoint was verified
against SI p.1, the Crossref record and the Europe PMC record before any
chemistry was attempted.

---

## 1. What we hold

| | |
|---|---|
| file | `source/ja5c01389_si_001.pdf` |
| bytes | 4,493,154 |
| sha256 | `668f5fc7ce7bf4dc1b1dbae27c8d6189e0ff70b0b56bc0b170a219bd2d31c161` |
| pages | 54 — Experimental Procedures 2–12, Figures S1–S37 13–53, References 54 |
| article | **NOT HELD.** Europe PMC `isOpenAccess=N`, `inPMC=N`, one full-text URL marked "Subscription required". No paywall or bot-check bypass was attempted. |
| licence | ACS standard (`10.15223/policy-029`, `-037`, `-045`, all `stm-asf`). **Not CC0/CC-BY** — see `DECISIONS.md` D-6. |

## 2. Anchor posture — recorded by hand

`provenance/structure_check.py` (C10), `check_anchor_independence.py` (C16) and
`check_anchor_coverage.py` (C18) were **not run**: they write tracked files that
other agents are using. This section is the hand record they would otherwise
produce. Nothing in `data/`, `provenance/` or `curation/lib/` was modified by
this curation.

**Anchor source.** Figure S18 (SI p.30) draws four lipids in full:
`A5-CE-C10`, `A5-C10`, `A5-CE-C7-6`, `A5-C7-6` — which are exactly the four
lipids this folder proposes to emit. Read at 200 dpi on 2026-11-10.

**Coverage.** 4 anchored / 15 enumerated. The other 11 split into two groups,
both declared rather than left silent:

- 2 extra-tail controls and 4 further `A5-CE-Cm` members that carry HRMS and
  ¹H NMR but no full drawing — structure is externally anchored by mass, not by
  a drawing, so no `reference_SMILES` row is claimed for them;
- 5 lipids named only on Figure S5's x axis, which carry **an explicit empty
  `reference_SMILES`** in `reference_structures.csv` so the gap is visible
  (`DECISIONS.md` D-5).

**Independence, and what is actually asserted.** The enumerator
(`stage02_structures.py`) builds **core-first**: it parses an atom-mapped
tetraol core and calls `assemble.attach` once per tail, so core atoms are
numbered before any tail exists. The anchors (`make_reference.py`) are
**hand-transcribed string literals**, typed tail-first in one left-to-right
pass while reading Figure S18, with no fragments, no atom maps, no RDKit
editing, and no import of `assemble` or of the enumerator.

Result: **4/4 anchors match the enumeration by full InChIKey and by molecular
formula.**

**2 of the 4 hand-typed templates are byte-identical to the enumerator's
canonical output** (`A5-CE-C7-6`, `A5-C7-6`), and that is *not* treated as a
failure. `CLAUDE.md` §6 states plainly that byte-identity of the emitted anchor
is not the signal and that testing for it is a trap, because a correct
re-derivation canonicalises to the same string as a correct enumeration —
which is exactly why C10's own `anchor_may_be_a_copy` missed F-037. An earlier
version of `make_reference.py` raised on byte-identity; it would have forced a
correct anchor to be rewritten into a worse one until the bytes differed. That
guard was removed and replaced with C16's actual property, tested on this
file's own AST: **every anchor SMILES must be a source string literal, and the
module must import no assembler.** Seeded by turning one literal into a
`.replace()` call; the script exits 1. Both states were verified from the
process exit code directly.

## 3. Lineage

```
ja5c01389_si_001.pdf
  ├─ pages 13-54 : one flattened JPEG each (0 vector drawing items)
  │    └─ figs/*.jpeg          native-resolution extraction, no resampling
  ├─ stage01_values.py         -> 01_digitised.json   (23 points, 5 panels)
  ├─ stage02_structures.py     -> 01_blocks.csv, 03_enumerated.csv,
  │                               02_structure_qc.json   (15 structures)
  ├─ (corpus scan)             -> 04_corpus_overlap.json
  ├─ stage03_emit.py           -> 05_values.csv        (quality cols via emit)
  ├─ stage07_checks.py         -> 07_qc.json           THE GATE, exits non-zero
  ├─ stage07_seed_test.py      -> 10/10 checks credited
  └─ make_reference.py         -> reference_structures.csv
```

`figs/` holds images extracted from the SI at their embedded resolution rather
than page renders: PIPELINE_V2 §3 says to prefer a publisher raster to a
re-rasterised page, and the same logic applies within the PDF. Rendering page
36 at 600 dpi would upsample a 627 × 463 JPEG and add no information.

## 4. Nothing was merged

`data/` is untouched. `05_values.csv` marks 4 rows
`emitted_as_corpus_candidate=True` and 19 held. Eight `DECISIONS.md` entries
are `OPEN`, so `build_dataset.py` would refuse this folder as it stands, which
is the intended state.
