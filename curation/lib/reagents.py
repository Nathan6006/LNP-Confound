#!/usr/bin/env python3
"""Reagent name -> structure, refusing on ambiguity, plus the cross-experiment
check that error class 46 asks for and nothing has implemented.

PRECEDENCE, which matters more than the code here:

  1. provenance/reagent_registry.py  for a catalogue reagent whose structure is
     fixed and public. That file writes the structure FROM THE CHEMICAL NAME,
     so it does not inherit whatever the corpus happens to contain. 15 reagents
     are pinned. Prefer it always.
  2. curation/corpus_components.py   for a component that must come from the
     corpus. Raises on an ambiguous name; never resolves by row order.
  3. nothing else. A name typed into a dict in a curation script is how `18PG`
     became symmetric 18:0/18:0 when it is SOPG 18:0/18:1.

THE NEW PART: cross_check(). Error class 46 -- "a reagent whose structure
already exists, independently curated, elsewhere in the corpus, and was never
consulted" -- was logged on 2026-11-09 and never implemented. For every reagent name carried by MORE THAN ONE experiment,
compare the skeleton InChIKeys. The corpus is its own best external anchor: a
second independent curation of the same reagent is the cheapest C10 there is,
and it is already paid for.

It would have caught F-031 for free. `ZT_2026` read `OF-2` off a figure by eye
with ALL FOUR TAILS one carbon short (C17 for C18), while `LL_2025` held the
same reagent correctly as `OF-02`. Two curations, one reagent, two structures,
and nothing compared them.

WHAT IT MUST NOT DO is report the 288 generic combinatorial labels as errors.
`A1_B1_C1` appears in both `JC_2023` and `YX_2024` and they are unrelated
libraries that happen to share a labelling scheme -- CLAUDE.md section 3 is
explicit that this is not an error and that one must never key on `IL_name`
alone across experiments. So `cross_check` SEPARATES the two populations:

  GATED    a name in the registry, or carried in a non-IL component column
           (HL/CHL/PEG/fifthcomponent). Excipients are catalogue reagents, so a
           collision there is a real contradiction.
  REPORTED a collision in IL_name only, which is where the generic labels live.
           Reported with its structures for a human to read, never failed on.

`alias` normalisation reuses the registry's own ALIASES, so `Dlin-MC3-DMA` and
`DLin-MC3-DMA` are one name here as they are there.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
_PROV = os.path.join(HERE, "..", "..", "provenance")
_CUR = os.path.join(HERE, "..")
for _p in (_PROV, _CUR):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from rdkit import Chem                                    # noqa: E402

from errors import AmbiguityError                         # noqa: E402

try:
    from reagent_registry import REGISTRY, ALIASES, PREFIXES   # noqa: E402
except Exception:                                          # pragma: no cover
    REGISTRY, ALIASES = {}, {}
    PREFIXES = ("IL", "HL", "CHL", "PEG", "fifthcomponent")

try:
    import corpus_components                               # noqa: E402
except Exception:                                          # pragma: no cover
    corpus_components = None

GATED_PREFIXES = tuple(p for p in PREFIXES if p != "IL")


def normalise(name):
    n = (name or "").strip()
    return ALIASES.get(n, n)


def skeleton(smiles):
    """The InChIKey skeleton block, or None. Connectivity without stereo."""
    m = Chem.MolFromSmiles(smiles or "")
    if m is None:
        return None
    return Chem.MolToInchiKey(m).split("-")[0]


def full_key(smiles):
    """The FULL InChIKey. Stereo included -- the block a skeleton comparison
    cannot see, which is why C11 was blind to epicholesterol twice (F-013)."""
    m = Chem.MolFromSmiles(smiles or "")
    if m is None:
        return None
    return Chem.MolToInchiKey(m)


def registered(name):
    """The registry's SMILES for a catalogue reagent, or None if not pinned."""
    n = normalise(name)
    hit = REGISTRY.get(n)
    return hit[0] if hit else None


def smi(name, corpus_rows=None):
    """The one structure for `name`, by the precedence above, or raise.

    Tries the registry first. Falls back to the corpus only if `corpus_rows` is
    given, and then through corpus_components.smi, which refuses on ambiguity.
    """
    n = normalise(name)
    pinned = registered(n)
    if pinned:
        return pinned
    if corpus_rows is None:
        raise AmbiguityError(
            f"reagents.smi({name!r}): not in provenance/reagent_registry.py "
            f"({len(REGISTRY)} reagents pinned) and no corpus rows were "
            "supplied to fall back to.\n"
            "If this is a catalogue reagent, ADD IT TO THE REGISTRY, writing "
            "the structure from its chemical name -- that is what makes C11 "
            "cover it. Do not type a SMILES into a curation script: `18PG` "
            "became symmetric 18:0/18:0 that way when it is SOPG 18:0/18:1.")
    if corpus_components is None:                          # pragma: no cover
        raise AmbiguityError("reagents.smi: corpus_components unavailable.")
    resolved, ambiguous = corpus_components.load(corpus_rows)
    return corpus_components.smi(resolved, ambiguous, n)


def collect(corpus_rows, prefixes=PREFIXES):
    """(name, prefix) -> {experiment_id: set(smiles)} over the corpus."""
    out = {}
    for r in corpus_rows:
        exp = (r.get("Experiment_ID") or "").strip()
        for pre in prefixes:
            n = normalise(r.get(f"{pre}_name"))
            s = (r.get(f"{pre}_SMILES") or "").strip()
            if n and s and exp:
                out.setdefault((n, pre), {}).setdefault(exp, set()).add(s)
    return out


def cross_check(corpus_rows, gated_prefixes=GATED_PREFIXES, stereo=False):
    """Compare each multi-experiment reagent name's structures across experiments.

    stereo=False compares skeleton InChIKeys (connectivity). stereo=True
    compares full InChIKeys, which is what catches an inverted centre -- C11
    compared skeleton blocks only and was blind to epicholesterol twice.

    Returns dict(gated=[...], reported=[...], n_names=int, n_multi=int).
    Each entry: name, prefix, experiments, keys {key: [experiments]}, smiles.
    The CALLER decides the exit code; this function never exits.
    """
    keyfn = full_key if stereo else skeleton
    seen = collect(corpus_rows)
    gated, reported = [], []
    n_multi = 0
    for (name, pre), byexp in sorted(seen.items()):
        if len(byexp) < 2:
            continue
        n_multi += 1
        keys = {}
        for exp, smis in byexp.items():
            for s in smis:
                k = keyfn(s)
                keys.setdefault(k or "UNPARSEABLE", set()).add(exp)
        if len(keys) < 2:
            continue
        entry = {
            "name": name, "prefix": pre,
            "n_experiments": len(byexp),
            "n_distinct_structures": len(keys),
            "keys": {k: sorted(v) for k, v in sorted(keys.items())},
            "registered": registered(name) is not None,
            "comparison": "full_inchikey" if stereo else "skeleton_inchikey",
        }
        if pre in gated_prefixes or entry["registered"]:
            gated.append(entry)
        else:
            reported.append(entry)
    return {"gated": gated, "reported": reported,
            "n_names": len(seen), "n_multi_experiment": n_multi}


def explain(entry):
    """One paragraph a reviewer can act on, for a cross_check entry."""
    lines = [
        f"{entry['name']} ({entry['prefix']}_name) carries "
        f"{entry['n_distinct_structures']} different structures across "
        f"{entry['n_experiments']} experiments, compared by "
        f"{entry['comparison']}:"]
    for k, exps in entry["keys"].items():
        lines.append(f"    {k}  <- {', '.join(exps)}")
    if entry["prefix"] != "IL" or entry["registered"]:
        lines.append(
            "  This is GATED. An excipient or a registered catalogue reagent "
            "has one structure, so two curations disagree and one of them is "
            "wrong. Decide WHICH before changing anything: in AC_2025 the "
            "NAMES were wrong and the structures were right, and repairing it "
            "the other way would have moved 419 measured values onto the wrong "
            "molecules.")
    else:
        lines.append(
            "  This is REPORTED, not gated. Generic combinatorial labels are "
            "reused by unrelated libraries -- 288 IL names legitimately carry "
            "different structures in different experiments (CLAUDE.md section "
            "3). Read it, do not key on IL_name across experiments, and only "
            "escalate if these two papers really do mean the same compound.")
    return "\n".join(lines)
