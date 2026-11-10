# _scaffold_tmp_kj -- provenance

**Status: SCAFFOLDED 2026-11-10. Nothing extracted yet.**

| | |
|---|---|
| Paper | TODO: authors, journal, year |
| DOI | TODO |
| PMID | TODO -- and a paper may be recorded as unindexed ONLY with a dated query result in `provenance/esummary_cache.json`. "It is a preprint" is not a reason: PubMed has indexed bioRxiv deposits since 2020, and that false claim has already been made twice (F-027, and `QL_2026` before it). |
| Experiment_ID | TODO |
| Rows | TODO |
| Value source | see `00_value_sources.json` |

## Where the values come from

Filled from `00_value_sources.json`. Any panel recorded as `5_pixels` must say
here why tiers 1-4 failed, in a sentence a reviewer can check.

## The paper's own cross-checks

List every statistic the paper derived from a figure, and whether we reproduce
it. **Passing all of them is necessary and not sufficient** -- `QL_2026`
reproduced every printed hit rate on bars up to 28% too high, and `XD_2026`
landed 12/12 inside its stated band on tops 0-21 px wrong. Geometry has to be
checked geometrically.

Better than a printed statistic: a round number the authors DREW and never
printed -- a dashed threshold, an axis limit at a round value, a colorbar's own
tick notch. Those are external to our pipeline in a way a recomputed marginal
is not.

## Contradictions

See `08_contradictions.csv`. Every paper curated so far has contradicted
itself; `NG_2025` alone yielded seven. An empty file means nobody looked.

## C10 anchor

`make_reference.py` BUILDS the anchor; `reference_structures.csv` stores what it
built. The generator must share as little as possible with the enumerator --
C16 tests the generator's composed template, because byte-identity of the
emitted anchor is not the signal (a correct re-derivation is byte-identical to a
correct enumeration once canonicalised, which is exactly why C10's own
`anchor_may_be_a_copy` missed F-037).
