# mckinlaycj_39105751 — Webster … McKinlay, *ACS Nano* 2024, PMID 39105751

**Not merged.** `DECISIONS.md` has five `OPEN` entries and `build_dataset.py`
refuses a folder with any; confirmed by running its `open_decisions()` against
this folder, which returns all five.

| | |
|---|---|
| DOI | 10.1021/acsnano.4c05513 |
| PMCID | PMC11342374 (open access — the source of every emitted value) |
| Value source | see `00_value_sources.json` |
| Library | 221 clustered screening peptoids (SI Table S1) + 45 optimisation peptoids 222–266 (SI Table S9) |
| Enumerated | **44** structures (222–266 less peptoid 231) |
| **Rows with a value** | **5** |

## Where the values come from

**Main-article Table 1, and nowhere else.** It is carried in the PMC JATS record
as structured `<table-wrap>` markup, so the five values are read from declared
table cells — **tier 4, zero pixels, no digitiser in this folder**.

Everything else is an acquisition gap, not a digitisation failure. The
per-member delivery values for the 93 screened library members (Figure 2C), the
12 headgroup variants (Figure 3B), the 34 DOE lipid-block variants (Figure 5D)
and the 4 aRSV formulations (Figure 6A) all live in main-article figure rasters.
PMC served **one** of eight (Figure 1, 600 dpi, the monomer pool); the other
seven returned a bot-check interstitial on a single plain request each and were
**not retried** — CLAUDE.md's non-negotiables forbid working around a bot check.

The SI's nine tables were inventoried exhaustively (`00_value_sources.json`,
`si_tables_inventory`) and carry **particle properties only** — size, PDI, %
encapsulation, zeta, pKa, hemolysis — none of which is a delivery endpoint under
the corpus definition. `stage00_scan_for_values.py` re-derives that negative
claim in one command: of 46 SI pages, the only two carrying both expression
vocabulary and a flux-magnitude number are pages 20 and 24, and on both the
number is an **axis tick label**, not a table cell keyed to a peptoid.

## These are peptoids, and order is the whole point

N-substituted glycine oligomers, 4–9 residues, built C→N on Rink amide resin,
so the C-terminus is a primary carboxamide and the N-terminus keeps its
hydrogen. The one trap: **Dea and Eae are secondary amines**, so their backbone
nitrogen carries two arms and no hydrogen. Getting that wrong is exactly one NH
(15.011 Da) and it is what the mass check caught on the first run — peptoid 227
came out at 1077.95 against the printed 1062.94.

**A formula cannot see monomer order.** Peptoids 236, 237, 259, 264 and 266 all
have calculated exact mass 1275.09 and are five different molecules; 266 was
synthesised *deliberately* as 236's composition-matched sequence isomer. So the
exact-mass check, strong as it is, confirms composition and nothing else. Order
is confirmed twice independently: SI Table S4 reprints all 34 sequences
(34/34 agree with Table S9), and SI Figure S7A **draws 236 and 266 in full**.

## The paper's own cross-checks, and one it never printed

| check | result |
|---|---|
| enumerated mass vs SI Table S9's printed **Calculated Exact Mass** | **44/44** within 0.005 Da |
| SI Table S4's three printed columns (lipid monomers, lipid carbons, branched/unsaturated) recomputed | **34/34** on all three |
| SI Table S4 vs SI Table S9 sequences | **34/34** identical |
| Table 1 values parsed with an explicit `<sup>` exponent | **5/5** |
| Dod and Ehx chain lengths from the PDF's **own stroked bond coordinates** | 13 and 9 merged vertices → **12 C** and **8 C** |

**Better than a printed statistic — a regularity the authors drew and never
printed.** SI Table S9 has a `Found Exact Mass` column and never says what
species was found. Across the table the difference from `Calculated` is
**+1.0089 ± 0.0062 Da on 35 of 45 rows** — the proton. That detector is external
to our pipeline, and the ten rows that miss it turn out to be **the paper's own
data-entry errors**: peptoid 253's found mass is byte-identical to 245's and
255's, and 262's is printed as `2.00`. None touches a structure we emit, because
every structure uses the `Calculated` column.

Passing all of these is **necessary and not sufficient**. This folder avoids the
digitisation error classes by not using the panels, not by being careful with
them.

## Contradictions

`08_contradictions.csv`, eleven of them. The ones that cost something: Figure 1B
labels **two different structures `Mpp`**, which is the single largest reason the
221-member screening library is not enumerated here; Figure 1B prints `C12H23`
for dodecyl (C12H25); and the paper's stated optimum is **55** total lipid
carbons in the text but **52.5 ± 8.8** in Table 1 — a difference that flips the
same paragraph's claim that peptoid 245's 44 carbons are "just outside of the
predicted range" (43.7 is the lower bound, so 44 is inside).

## C10 anchor

`make_reference.py` rebuilds **236 and 266** and both reproduce the enumeration
exactly. Independence is on four axes and all four are executed, not asserted:

- **different input** — the monomer order comes from the fill colour and x-centre
  of SI page 20's own rounded side-chain boxes, sorted; the enumerator reads SI
  Table S9's text layer;
- **different traversal** — C→N here, N→C there;
- **different construction** — plain SMILES string composition, no `assemble`,
  no `RWMol`, no atom maps;
- **`blocks.py` is not imported** by `make_reference.py`.

`reference_structures.csv` also carries **empty** anchor rows for 226, 245 and
263, because the paper draws them nowhere and that absence should be visible.

## Seeded defects

`07_qc.json`, seven of them, and each is credited only if an **agreement count
moves** — not merely if an exception is raised. Two seeds failed on their first
run and both were real defects in this folder:

1. the watermark seed reported `34/34 → 34/34` because `_rows(x_max=X_TABLE_MAX)`
   bound the module global as a **default argument**, so the mutation was
   shadowed and the seed tested nothing;
2. the order-transposition seed reported `2/2 → 1/2` because it swapped
   positions 1 and 2 of peptoid 266, which are **both Ehx** — half the seed was
   a no-op that read as a partial pass.

The headline seed: transposing two monomers leaves the exact mass at **1275.09
either way** and the C10 anchor still goes `2/2 → 0/2`. That is the one check in
the folder that can see a sequence error.
