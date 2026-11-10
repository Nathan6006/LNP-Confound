# PROVENANCE — chamolstad_41691524

Kim JG, Xu J, Lee D, … Mitchell MJ, Cha-Molstad H.
**Comprehensive Engineering of Ionizable Lipid Nanoparticles and mRNA Elements
for Next-Generation Vaccines.** *ACS Nano* 2026.
PMID **41691524** · DOI **10.1021/acsnano.5c19036**.

Title verified against the Europe PMC REST `core` record for `EXT_ID:41691524`,
retrieved 2026-11-10 UTC. That record also gives `isOpenAccess: N` and
`license: None`.

## 1. The endpoint

**mRNA delivery.** Firefly-luciferase mRNA transfection *in vitro* (DC2.4,
C2C12), whole-animal bioluminescence *in vivo*, and SARS-CoV-2 RBD vaccine
immunogenicity. Abstract: *"From a 96-member biodegradable ionizable lipid
library, stereodefined H9T6 emerged as the lead ionizable lipid… achieved more
than 3.5-fold higher dendritic-cell transfection than SM-102."*

The paper contains a **second, separate screen** of **mRNA elements** —
>630,000 3′ UTR variants by half-life and >68,000 5′ UTR sequences by ribosome
loading. That varies the cargo, not the lipid, and is a different experiment.
Nothing from it is emitted.

## 2. Sources on disk

| file | sha256 | bytes |
|---|---|---|
| `source/nn5c19036_si_001.xlsx` | `3cd5f22e5d1ae6c7a39d026d05c5ab58591d50ce481522901baa9f2bce5fe5a1` | 16,583 |
| `source/nn5c19036_si_002.pdf` | `baa616a0faa6c2d0762dcc226b49e1893c3f2fb55e1bf0a4e833771424ab5c13` | 4,092,216 |

Both were already present when this curation started; neither was fetched here,
so no acquisition URL or landing-page licence string is recorded. **The licence
is unrecorded and the Europe PMC record says the article is not open access** —
that must be settled before any corpus release (`CLAUDE.md` §8).

The main-text PDF is **not** on disk and was not obtained: ACS is paywalled and
no paywall or bot-check bypass is permitted.

## 3. Where the values came from — the one-line version

The tier-1 publisher Source Data workbook is **not a lipid screen**. It is
11 rows × 21 columns of 3′ UTR decay data (9 selected UTRs, CPM at 0/1/2/4/8 h,
four rate constants, eleven half-life replicates). `find_screens.py` printed
nothing for it, which is the right answer.

Full tier-by-tier accounting, including every panel deliberately NOT used and
why, is in **`00_value_sources.json`**, written before any chemistry.

Emitted values are **tier 5, pixels**, from Figure S4A and S4C on SI p. 11.

## 4. Anchor posture (recorded by hand)

`provenance/structure_check.py`, `check_anchor_independence.py` and
`check_anchor_coverage.py` were **not run** — they write tracked files that
other agents are using. This section stands in for their output.

**`reference_structures.csv` holds one row, H9T6, with a non-empty
`reference_SMILES` and an empty `expected_mismatch`.**

**The paper draws no structure of H9T6 anywhere in the SI.** The anchor is
therefore derived from the SI's **written preparation** (p. 5: tris(2-amino-
ethyl)amine + 5.7 equiv of the tail; p. 4: 0.95 × N equiv with N = 2 ×
(primary amines) = 6), not from a drawing. Stated plainly: **this anchor cannot
catch an error that the written procedure also contains.** It tests the
enumerator against the paper's words.

**Independence of the generator, which is what C16 tests** — not byte-identity
of the emitted anchor, because a correct re-derivation is byte-identical to a
correct enumeration once canonicalised (F-037):

| | `02_enumerate.py` (enumerator) | `make_reference.py` (anchor) |
|---|---|---|
| bond formation | `RWMol.RemoveBond` on the oxirane, then `assemble.attach` with atom maps, six times | **none** — plain `str.replace` on a template |
| RDKit calls | `RWMol`, `AddBond`, `RemoveBond`, `SanitizeMol`, `assemble.esterify`, `assemble.attach` | `MolFromSmiles`, `MolToSmiles`, `AssignStereochemistry`, `MolToInchiKey` only |
| stereocentre | taken from a CIP selection over both chiral tags on a **flat ester constitution** | written once, in the arm template, **from the nitrogen outward** (opposite traversal) |
| site selection | SMARTS `[NX3;H1,H2;!$(N-C=O)]`, "first remaining N-H", six passes | six literal `%` marks in the core template |

Mechanically checked with `grep -cF`: the generator's composed template
`N(CCN(%)%)(CCN(%)%)CCN(%)%` and its arm `%C[C@H](O)COC(=O)CCCCCCCCC` occur
**0 times** in `02_enumerate.py`. `RWMol` and `attach` appear in
`make_reference.py` only inside its docstring, describing the other path.

**Result: FULL InChIKey match, stereo layer included.**
`UMBWDSGQAXLQTI-DVIOPHQQSA-N` from both. All six CIP codes are `S`.
The canonical SMILES are byte-identical (169 chars) — expected, and not the
signal.

No other structure in this folder has an anchor. The 4- and 5-tailed LC-MS
congeners in `03_enumerated.csv` are positional-isomer mixtures with no
determined constitution and are **not emitted**; H9T2 is not emitted (D-6).

## 5. Corpus overlap — reported, not acted on

Checked against the live corpus (`data/data_vitro.csv`, 29,195 rows, 18,720
unique parseable `IL_SMILES`; `data/data_vivo.csv`, 1,803 rows, 1,506 unique):

* **full InChIKey**: 0 hits for H9T6, for its 4- and 5-tailed congeners, or for
  the (S)-glycidyl decanoate tail.
* **InChIKey skeleton block**: 0 hits, so this is not a stereochemistry-only
  miss.
* **`IL_name`**: 0 rows named `H9T6` or `H9T2`; 0 rows anywhere in either pool
  whose `IL_name` matches `H<digits>T<digits>`.

Overlap is **zero**, in both directions, on both keys. No `Experiment_ID` to
report.

## 6. What was not done

* `data/` is untouched. Nothing merged.
* `curation/lib/` is untouched.
* The shared Playwright browser was not used; the only network call was a plain
  `curl` to the public Europe PMC REST endpoint.
* `provenance/` is untouched.

## 7. One finding about our own tooling

`assemble.esterify` **refused a correct build**. Esterifying (S)-glycidol with
decanoic acid relabels C2 from `S` to `R` with **no spatial change**: glycidol's
free OH leaves the C1/C3 CIP comparison tied at the first sphere, and the ester
breaks that tie. `check_stereo_preserved` compares CIP **counts**, so it reads
the relabel as an inversion and raises.

The refusal was respected, not worked around. The tail's stereocentre is instead
derived from the compound's **published name** — SI p. 4 names it
"(S)-(–)-Glycidyl decanoate" — by building the flat constitution with
`assemble.esterify` on achiral glycidol and then selecting, of the two chiral
tags, the one RDKit's CIP implementation labels `S`. No tag is carried through a
bond-forming step at all.

`curation/lib/assemble.py` was **not edited**. The behaviour is reported here so
someone who owns that file can decide whether `check_stereo_preserved` should
compare a numbering-independent parity instead of CIP codes; a CIP code is a
property of the whole molecule, so it is not conserved by a transformation that
conserves geometry.
