# DECISIONS -- _scaffold_tmp_kj

Questions this curation must NOT answer on its own. `build_dataset.py` refuses
to merge this folder while any entry below is `OPEN`.

WHY THIS FILE GATES THE MERGE. The pipeline is good at detecting that a source
is silent. It is not entitled to resolve it. In `AC_2025` the first repair was
confident, internally consistent, and would have assigned **419 measured values
to the wrong molecules** -- the names were wrong and the structures were right,
and nothing in the data said which. So the rule is: detect, queue, and stop.

Entries that belong here:

- which half to repair when a NAME and a STRUCTURE disagree;
- whether a `confidence = low` library belongs in the corpus at all;
- a value with several mutually inconsistent candidates (F-004's `N:P` had
  three -- 7.199, 16.766, 6.0 -- and was blanked, on the user's decision);
- anything that would change a published number.

## Format

```
### D-1  <one line>
STATUS:    OPEN | RESOLVED <date> by <who>
AMBIGUITY: what the source does not say.
CANDIDATES: each reading, and what it implies for the emitted rows.
EVIDENCE NEEDED: what would settle it.
DECISION:  (when resolved) what was chosen and why. Keep the entry.
```

## Entries

### D-1  placeholder -- delete when the first real decision is recorded
STATUS:    RESOLVED 2026-11-10 by scaffold
AMBIGUITY: none; this entry exists so the file is never empty by accident.
DECISION:  n/a.
