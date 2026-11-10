# PROVENANCE — `curation/mitchellmj_41555650`

Han, E. L. *et al.* **High-Throughput In Vivo Screening Identifies Structural
Factors Driving mRNA Lipid Nanoparticle Delivery to the Brain.** *ACS Nano*
(2026). PMID **41555650**, DOI **10.1021/acsnano.5c19138**.
Lab: **Michael J. Mitchell**, University of Pennsylvania.

> The extraction brief named this a **Dahlman**-lab paper and asked for the
> folder `curation/dahlmanje_41555650`. It is a **Mitchell**-lab paper. The
> empty wrong-name folder was removed; all work is here. `find_repo.py` was run
> with `--author Mitchell`, not the brief's `--author Dahlman`, because a wrong
> author term returns a zero that looks like a clean negative.

## Sources

| | |
|---|---|
| Main article PDF | **NOT ACQUIRED.** ACS Cloudflare; shared browser forbidden for this run; a curl probe of the ACS supplementary routes was refused by the sandbox as bot-check circumvention and **no bypass was attempted**. |
| SI PDF | `source/nn5c19138_si_001.pdf`, 9,108,406 bytes, 28 pages, `%PDF-1.3` verified by magic bytes. sha256 `c9206958b9debb76fdbb63d81b10e1d6777fa5aa7c2f3f58c616a9509f8e1ba5`, md5 `c4cbe865e34aa6d9e3d5a5d3c64802da`. Pre-staged before this run; page footers show it was fetched through an entitled institutional session that is not ours. |
| PMC deposit | **None.** NCBI `elink` returns only `pubmed_pmc_refs` (4 citing ids), no `pubmed_pmc` SELF-link. No PMCID exists to take, so `pmc_fetch.py` was not run. |
| Author repository | **None.** `find_repo.py`: DOI and PMID code search 0 hits; `Mitchell LNP`, `mRNA barcoding brain`, `mRNA barcoding lipid` 0 hits each; 12 irrelevant candidates filtered. |
| Licence | ACS Nano, `isOpenAccess = N`, Europe PMC `license = null`, full-text availability *"Subscription required"*. **Not CC0/CC-BY.** See `DECISIONS.md` D-6. |

## C10 anchor posture — recorded BY HAND

`provenance/structure_check.py`, `check_anchor_independence.py` and
`check_anchor_coverage.py` were **not run**: they write tracked files other
agents read, and this run was instructed not to touch them. The posture below is
therefore stated here and is reproducible from `make_reference.py`, which exits
non-zero if the anchors disagree.

**Coverage: 2 of 15 lipids carry an anchor.** `reference_structures.csv` holds a
row for all 15; the other 13 carry an **empty `reference_SMILES`** with a note,
so the absence is visible rather than silent.

**The paper draws no full structure anywhere we hold.** The SI's structural
content is systematic IUPAC names (pp. 2–4) and LC-MS traces (pp. 4–12); any
drawn structure is in the unacquired main article. So the anchors are built from
the **systematic names**, which pin every bond, and not from a drawing.

**Independence is a property of the route, and it is executed, not asserted.**
`make_reference.py::audit_independence` walks **its own AST at run time** and
confirms the module imports no assembler and calls none of `RWMol`, `AddBond`,
`RemoveAtom`, `CombineMols`, `esterify`, `attach`, `AddAtom`,
`ReplaceSubstructs`, `ReactionFromSmarts`, `RunReactants`. The AST is walked
rather than the file substring-scanned, because a substring scan flags the words
in the module's own docstring.

| | enumerator (`03_enumerate.py`) | anchor (`make_reference.py`) |
|---|---|---|
| start point | the amine core | the far methyl of a lipid tail |
| direction | head → tail | **tail → head** |
| bond formation | `assemble.attach`, one arm per N-H | none — one hand-written SMILES literal per anchor |
| shares with the other | — | no fragment string, no core string, no bond-forming code |

Anchors are **different chain lengths and different cores** so one systematic
slip cannot pass both: `C10-306` (C10 tails, 4 arms, acyclic triamine) and
`C14-200` (C14 tails, 5 arms, piperazine pentamine).

**Byte-identity is not the signal (F-037).** The canonical forms are identical to
the enumeration, which is what two correct derivations of the same molecule must
produce. What was authored is the written string, recorded verbatim in
`reference_report.json` alongside the canonical form.

**The gate was made to fail on purpose.** Seed: one arm's hydroxyl moved from C2
to C3 on `C10-306` — **identical molecular formula `C47H99N3O4` and exact-mass
delta 0.000000**. The InChIKey comparison flips `True → False` and the seeded
process exits 1 (`07_seed_results.json`, check 1). Mass and formula are blind to
it; connectivity is not.

## Pooled screen — admissibility, corrected

My first verdict was that a pooled barcode value is not commensurate with
`data/`. **That was wrong**, and the correction is recorded in `DECISIONS.md`
D-1 and `08_contradictions.csv`. `data/data_vivo.csv` has an
**`Experiment_batching`** column, vocabulary `{barcoded, individual}`, and
**`LR_2023` contributes 407 barcoded in-vivo rows** over 383 unique ionizable
lipids. Pooling alone is not disqualifying; the corpus declares it on the row.
The screen is excluded here purely because **it was not acquired**. The 48
emitted rows carry `Experiment_batching = individual`.

The pool itself is pinned: **78 of 98 b-LNP numbers** went into the pooled
administration (12 greyed, 8 never tabulated), recovered from character colour
and confirmed against the caption's four stated thresholds on **90/90 rows, 0
disagreements** (P14). The screen's **normalisation residual was not computed** —
the formula and values are in the unacquired article (D-2).

## Geometry, if you are checking the digitisation

Both digitised panels are **embedded rasters**, confirmed by
`page.get_drawings() == 0` and exactly one image on each of SI pp. 13–21.
Native resolution measured as `px / rect_pt × 72`:

| panel | native dpi | `LOO_TOL_PX` @600 dpi | rescaled to native | measured LOO worst |
|---|---|---|---|---|
| Fig S2A (p. 14) | **355.7** | 1.5 | **2.530** | **0.370 px** |
| Fig S1A (p. 13) | **220.2** | 1.5 | **4.088** | **0.500 px** |

Both pass the **strict** 600 dpi tolerance, so the widened tolerance is reported
but **not relied on**. It was seeded anyway, in this panel's own geometry
(`07_seed_checks.py` check 8): a 5.0 px tick displacement is caught at the
native tolerance; 3.25 px is the known weak margin and is reported, not hidden.

Scale and origin come only from the seven detected decade tick marks. A round
number the authors drew and never printed: on **both** panels the lowest decade
tick sits on the drawn baseline to **0.000 px**.
