#!/usr/bin/env python3
"""ONE guarded way to look a component's structure up in the live corpus.

WHY THIS FILE EXISTS. Several emitters needed the SMILES for `DOPE`,
`Cholesterol`, `DMG-PEG2000` and similar, and each grew its own version of:

    comp.setdefault(name, smiles)      # first row wins

over `data_vitro.csv`, which is written in alphabetical experiment order. So the
structure a name resolved to was decided by **which experiment sorts first**.
That is error class 30, and it has already cost this project once: `MC3-CB1/2/3`
were built on the KC2 dioxolane because `DLin-MC3-DMA` was ambiguous in the
corpus -- 2,190 rows of the correct ester plus 203 dioxolane rows from
`AC_2025`, which sorts before everything -- and the lookup silently returned the
wrong one (F-014). Same formula, same exact mass to four decimals, invisible to
every composition check.

Round 1 hardened exactly one caller, `build_new_curation.py`, and left the
others. Round 2 found the same construct still live in three more emitters, with
`DLin-MC3-DMA` among the names one of them resolves -- unambiguous *today* only
because an unrelated name repair landed, which is the "fixed by luck of row
order, not by the guard" case that file's own docstring warns about (F-044).

So the fix is a shared, refusing loader rather than a fourth copy of the guard.
**`smi()` raises on an ambiguous name.** It does not pick, warn, or prefer: a
name that means two molecules is a question the caller has to answer, and
answering it by row order is how a measurement ends up on the wrong structure.

    from corpus_components import load, smi
    comp, ambiguous = load(corpus_rows)
    dope = smi(comp, ambiguous, "DOPE")

Prefer `provenance/reagent_registry.py` over this for a catalogue reagent whose
structure is fixed and public: that file states the structure from the chemical
name, so it does not inherit whatever the corpus happens to contain. Use this
only for components that must come from the corpus.
"""
PREFIXES = ("IL", "HL", "CHL", "PEG", "fifthcomponent")


def load(corpus, prefixes=PREFIXES):
    """name -> set of distinct SMILES, over every component column.

    Returns (resolved, ambiguous): `resolved` maps each unambiguous name to its
    single SMILES, `ambiguous` maps each multi-structure name to the full set.
    Nothing is dropped and nothing is chosen -- the caller is told.
    """
    seen = {}
    for r in corpus:
        for pre in prefixes:
            n = (r.get(f"{pre}_name") or "").strip()
            s = (r.get(f"{pre}_SMILES") or "").strip()
            if n and s:
                seen.setdefault(n, set()).add(s)
    resolved = {n: next(iter(v)) for n, v in seen.items() if len(v) == 1}
    ambiguous = {n: sorted(v) for n, v in seen.items() if len(v) > 1}
    return resolved, ambiguous


def smi(resolved, ambiguous, name):
    """The one SMILES for `name`, or raise. Never guesses."""
    n = (name or "").strip()
    if n in ambiguous:
        raise SystemExit(
            f"corpus_components.smi({name!r}): AMBIGUOUS -- this name carries "
            f"{len(ambiguous[n])} different structures in the corpus:\n  "
            + "\n  ".join(ambiguous[n])
            + "\n\nRefusing to pick one by row order. That is error class 30 and "
              "it put 419 measurements on the wrong molecules once already "
              "(F-014). Resolve it: add the reagent to "
              "provenance/reagent_registry.py and take the structure from there, "
              "or patch the experiment whose rows are wrong.")
    if n not in resolved:
        raise SystemExit(
            f"corpus_components.smi({name!r}): not found in any component "
            f"column of the corpus. Check the spelling against the "
            f"{len(resolved)} names that are present.")
    return resolved[n]
