# PROVENANCE — `KG_2025`, and the C10 anchor posture

Written by hand, because this curation was run under instructions **not** to
execute `provenance/structure_check.py`, `check_anchor_independence.py` or
`check_anchor_coverage.py` — those write tracked files other agents read. So
C10, C16 and C18 have **not** been run against this folder. What follows is the
posture a reviewer should check when they are.

---

## 1. What is anchored, and what is not

`reference_structures.csv` has ten rows. **Two carry a `reference_SMILES`; eight
are deliberately empty.**

| lipid | anchor | why |
|---|---|---|
| `C16 Galactosyl(a) Ceramide` | **YES**, from Fig. 1A | drawn in full, with stereochemistry |
| `C8 Galactosyl(a) Ceramide` | **YES**, from Fig. 1A | drawn in full, with stereochemistry |
| `SM-102` | empty | named in a box in Fig. 1A; **no structure is drawn anywhere in the article or SI** |
| `cKK-E12` | empty | same |
| `Cholesterol`, `DSPC`, `ALC-0159`, `DMG-PEG`, `C14-PEG2000`, `C18-PEG2000` | empty | named only |

**The ionizable lipids cannot be anchored by this paper.** C10 compares an
enumerated structure against the structure the paper DRAWS, and this paper draws
no ionizable lipid at all. Storing an empty `reference_SMILES` with a note is
the prescribed response, so that the *absence* of an anchor is visible rather
than implied (`PIPELINE_V2` §"The anchor").

So the **excipients are anchored instead** — which is also the right choice on
the merits: the two glycolipids are the paper's novel components and are
**absent from the corpus under every name tried** (`C8GaCer`, `C8GalCer`,
`GalCer`, `alpha-GalCer`, `C16GaCer`, `C16GalCer`), while SM-102 and cKK-E12 are
already in `data/` on 2,312 and 241 rows respectively.

When C18 (`check_anchor_coverage.py`) is run, this folder will report **8
components with no anchor**. Every one of them is a component the paper never
draws, and each says so in its `note` column. That is a true gap in coverage and
should be baselined as one, not explained away.

---

## 2. The independence claim, stated precisely

F-013 was an anchor that was a byte copy of the enumeration. F-037 was the file
written to *fix* F-013, whose amine half was still the enumeration's canonical
SMILES with the sterol swapped out. The lesson the project drew: **byte-identity
of the emitted anchor is not the signal, and testing for it is a trap**, because
a correct re-derivation is byte-identical to a correct enumeration once
canonicalised.

So here is exactly what is and is not independent in this folder. **I am not
claiming two independent derivations of the glycolipid SMILES string.**

### What IS independent: the inputs

| route | input | what it establishes |
|---|---|---|
| **pixels** — `stage04_chains.py` | Fig. 1A, the drawn structure | the N-acyl chain measures **16 carbons**: three methods (centreline turning points, mid-line crossings, apex vertices) all give 15 bonds, +1 for the carbonyl. `chains.check_name_token` confirms the token "C16" means the acyl chain — the `KS_2024` `C24` trap, tested and passed. |
| **the name** — `make_reference.py` | the string "C16 Galactosyl(α) Ceramide" | N-hexadecanoyl-1-O-(α-D-galactopyranosyl)-D-erythro-sphingosine, built from fragments whose stereodescriptors are checked against their IUPAC names |

Those two share no input: one reads a raster, the other reads a reagent name.
A drawing showing 14 carbons would have produced a mismatch.

### What is NOT independent

`stage05_emit.py` **reads** the SMILES out of `reference_structures.csv` rather
than constructing a second string. That is deliberate and is stated in the
emitter's own docstring. Writing a second derivation here would not have been
independent — it would have been me transcribing one decision twice, which is
precisely the F-037 failure mode dressed up as a check. **A reviewer should
treat the glycolipid's SMILES as having ONE source**, with the acyl length and
the α/4E features independently corroborated from pixels.

### The AST guard, and its own seed

`make_reference.assert_independent()` walks **this module's own AST at run time**
and refuses if it imports, reads or names the emitter. It is not a substring
scan — a scan would flag the module's own docstring, which names
`stage05_emit` several times.

**The guard is itself seeded.** My first two versions of C16 reported "all
derived" while examining nothing, so:

```
python3 -I curation/dahlmanje_41183194/make_reference.py . --prove-guard
  -> GUARD SEED: raised as required -- make_reference: NOT INDEPENDENT --
     this module reaches the emitter: ['import stage05_emit'].
  -> METRIC guard_fires=True
```

The seed injects a real `import` into the tree the guard examines. A guard
looking at nothing would not raise and the seed run would exit 0, which the
caller treats as a failure.

---

## 3. How the anchor was built, and the refusal that shaped it

**`assemble.attach` REFUSED the obvious route and the route was changed, not the
check.** Two refusals, in order:

1. Passing the intact carboxylic acid to `attach` for the amide bond raised
   `sanitize failed: Explicit valence for atom # 15 C, 5`. Fixed at the source
   by writing the acyl fragments already dehydrated — the acid loses its OH,
   which is `esterify`'s rule applied by hand because the partner is an amine.
2. The glycosidic bond then raised
   `check_stereo_preserved: the product's stereocentres are {R:6, S:1, ?:1}, the
   fragments' sum is {R:4, S:3, ?:0}`. This is a **CIP relabel, not an
   inversion**: forming a glycoside turns the anomeric OH into an OR and changes
   its substituent priorities, so two centres are relabelled S→R with no atom
   moving. `attach` adds a bond and removes nothing, so the F-013 `RemoveAtom`
   mechanism cannot arise in it — but it does not expose `allow_cip_relabel`,
   and editing a shared library mid-batch is itself a documented hazard.

So the molecule is written head-first from the two verified fragment strings,
each kept in **its own verified atom order**, and spatial parity is asserted
instead:

```
parity_preserved vs alpha-D-galactopyranose: True   centres_compared 5, flipped []
parity_preserved vs D-erythro-sphingosine:   True   centres_compared 2, flipped []
```

`parity_preserved` is numbering- and substituent-independent, so it is *stronger*
here than the fingerprint count: it compares centre by centre rather than
counting labels.

### The stereo assertions are NOT vacuous

Both fragments are chiral, and each is checked against a **discriminating
negative control** that differs at exactly one centre and shares its formula and
exact mass — the test the project demands after the epicholesterol episode,
where formula, exact mass **and the skeleton InChIKey block** all stayed
identical:

| fragment | IUPAC | control | shares formula | shares exact mass |
|---|---|---|---|---|
| α-D-galactopyranose | (2S,3R,4S,5R,6R)-6-(hydroxymethyl)oxane-2,3,4,5-tetrol | the **β anomer** | yes | yes |
| D-erythro-sphingosine | (2S,3R,4E)-2-aminooctadec-4-ene-1,3-diol | the **C2 epimer** | yes | yes |

Each control is asserted to produce a *different* CIP set; if it did not, the
check could not fail and `_verify_fragment_stereo()` raises saying so.

**A direction bug was caught this way.** The first version of the sphingosine
fragment copied `[C@H]` tokens from the verified string into a differently
ordered one (`[CH2:2](O)[C@H]...` rather than `OC[C@H]...`). A SMILES fragment
has a direction mass cannot see: the reordering silently produced the C2 epimer,
and it showed up as the stereocentre counts disagreeing, not as any formula or
mass difference.

### The string-composition hazard is assertion-gated

Fragment assembly was mandated after string concatenation fused piperazine to
cholesterol through a ring-closure digit collision. Because this build composes
strings, `build()` **asserts** that the sphingosine and acyl fragments contain
no ring-closure digit at all and that the sugar carries exactly one, so only the
sugar's `1` exists in the product. That assertion fired on the first run — the
acyl fragments still carried an atom map `[C:1]` — and the fragments were
rewritten without maps.

### Results

| | C16 Galactosyl(α) Ceramide | C8 Galactosyl(α) Ceramide |
|---|---|---|
| InChIKey | `VJLLLMIZEJJZTE-KAQXTUTFSA-N` | `LCEXEEHGNKGJES-MKWYEPMISA-N` |
| formula | C40H77NO8 | C32H61NO8 |
| exact mass | 699.564918 | 587.439718 |
| stereocentres, all assigned | 7 | 7 |

Mass difference **112.125200 Da = C8H16 exactly** — the octanoyl → hexadecanoyl
difference, a closed identity across the two anchors that no single count gives.

---

## 4. What a reviewer should do when the gates are run

1. **C10** will compare `reference_structures.csv` against the emitted
   `HL_SMILES`. For the glycolipid these come from the same string, so C10 will
   report `exact_match` while comparing a value to itself. **That is the F-013
   shape and it should be recorded as such**, not read as corroboration. §2
   above states what the real corroboration is.
2. **C10 should also flag `anchor_may_be_a_copy`** unless it credits the
   presence of `make_reference.py`. The generator exists and does build the
   anchor; what it does not do is build it twice.
3. **C16** works on the generator's composed template. The template here is
   `STEREO_SPEC[...]["smiles"]` plus `HEXADECANOYL`/`OCTANOYL`, none of which
   appears anywhere in `stage05_emit.py` — the emitter has no SMILES literal for
   the glycolipid at all. C16 should pass, and its own guard has been seeded
   separately (§2).
4. **C18** will record 8 components with no anchor. All 8 are undrawn. Baseline
   them with that reason.
5. **C11 / C17** (one reagent name, one structure): this curation introduces two
   reagent names the registry does not have. If either glycolipid is ever added
   to `provenance/reagent_registry.py`, note that the sphingoid base length is
   from the chemical name and not measured (DECISIONS D-6), so it is a weaker
   authority than the registry's other entries.
