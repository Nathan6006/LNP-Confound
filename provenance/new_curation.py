#!/usr/bin/env python3
"""Scaffold a new curation folder: stage stubs, intake gate, anchor, decisions.

    python3 -I provenance/new_curation.py <folder-name> <paper.pdf> [si.pdf ...]

WHY THIS EXISTS. curation/PROTOCOL.md already specifies a file-based, resumable
stage sequence, and it works. What was missing is that nothing GENERATED it, so
every paper's extraction code was written from scratch -- 50 scripts, 9,606
lines across 7 folders, of which 20 files and 3,853 lines read pixels or build
structures. That is how the same geometric defect (error class 20) came to have
four independent instances in four separately written files, the fourth of them
in a SIBLING FILE of a folder where the class had already been found and fixed
(F-033).

So this script writes stubs that import from curation/lib/ and refuse to be
skipped. It creates:

  00_pages/                 renders and pdftotext output land here
  00_value_sources.json     THE INTAKE GATE. Every panel or table that supplies
                            a value names the tier it came from, and anything
                            answering "pixels" must say why tiers 1-4 failed.
  00_acquisition.json       where each input file came from, with sha256
  DECISIONS.md              the questions this curation must NOT answer alone
  stage01..stage07 stubs    importing curation/lib, each raising until written
  make_reference.py         the C10 anchor GENERATOR -- not a stored string
  reference_structures.csv  with an empty anchor row, so a missing anchor is
                            visible rather than absent
  PROVENANCE.md             a skeleton whose claims C15 will check

It does not write extraction logic. A paper's own geometry constants, panel
bounds and library rules stay in the folder, next to the figure, where a
reviewer reads them.
"""
import csv
import hashlib
import json
import os
import shutil
import sys
from datetime import datetime, timezone

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
CURATION = os.path.join(ROOT, "curation")

VALUE_TIERS = [
    ["1_publisher_source_data", "a Source Data file published with the paper"],
    ["2_si_tables", "an xlsx/csv/docx table in the supplementary information"],
    ["3_author_repository", "the authors' own repo, Zenodo, Dryad, figshare"],
    ["4_pdf_text_table", "a real text table inside the PDF (pdfplumber)"],
    ["5_pixels", "digitised from a figure. THE LAST RESORT."],
]

STAGES = [
    ("stage01_structures.py", "read the paper's scheme; one structure per "
                              "library member", ["chains", "assemble"]),
    ("stage02_blocks.py", "the building blocks, from the figure that draws "
                          "them", ["chains", "reagents"]),
    ("stage03_enumerate.py", "combine blocks into the full library",
     ["assemble", "reagents"]),
    ("stage04_values.py", "the measured values, from the tier "
                          "00_value_sources.json names", ["axis", "bars",
                                                          "heatmap"]),
    ("stage05_formulation.py", "formulation, dose, route, cell line, metadata",
     []),
    ("stage06_verify.py", "the paper's OWN cross-checks, measured by script",
     ["axis"]),
    ("stage07_emit.py", "emit rows in the 79-column schema", ["emit"]),
]

HEADER = '''#!/usr/bin/env python3
"""{title}

GENERATED STUB -- replace this docstring with what this stage actually does,
and say where each number came from.

Run from this folder:   python3 -I {fname}

Imports from curation/lib, which holds ONE audited implementation per hazard
class. Do not copy a geometry helper in here: {classes}
"""
import json
import os
import sys

LIB = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib")
sys.path.insert(0, LIB)
{imports}

# ---------------------------------------------------------------------------
# HAZARD: THE SESSION SCRATCHPAD IS SHARED BY EVERY AGENT IN THE RUN.
#
# Measured 2026-11-10 (`mitchellmj_41603124`): a glob copy out of the scratchpad
# put **17 figure files belonging to OTHER PAPERS** into that curation's
# `source/figs/`. It was caught only because the agent listed the directory
# immediately after the copy and re-hashed the survivors against its own
# download record (6/6 identical). **A file-count or a size check would have
# passed**, and a digitiser pointed at a foreign figure produces values that are
# internally consistent, geometrically valid and attributed to the wrong paper.
#
# That is the worst failure shape this project has: not a wrong number, but a
# right number on the wrong molecule, which is the `AC_2025` class (419
# measurements, F-014).
#
# RULES:
#   * never `cp <scratchpad>/*` into a curation folder -- copy named files only;
#   * after any copy into `source/`, list the directory and assert every file is
#     one you fetched, by name AND by hash against your own download record;
#   * prefer fetching straight into `<folder>/source/` so the scratchpad is
#     never in the path at all.
# ---------------------------------------------------------------------------

def main():
    raise SystemExit(
        "{fname} is still the generated stub. Write it, or delete it and "
        "record in DECISIONS.md why this stage does not apply to this paper. "
        "An empty stage that silently returns is how three ZT_2026 marginal "
        "checks scored a literal 0 for twelve days while the summary beside "
        "them said ALL FOUR EXACT (F-038).")


if __name__ == "__main__":
    main()
'''

DECISIONS = '''# DECISIONS -- {folder}

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
STATUS:    RESOLVED {date} by scaffold
AMBIGUITY: none; this entry exists so the file is never empty by accident.
DECISION:  n/a.
'''

PROVENANCE = '''# {folder} -- provenance

**Status: SCAFFOLDED {date}. Nothing extracted yet.**

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
'''

MAKE_REF = '''#!/usr/bin/env python3
"""BUILD the C10 anchor for @@FOLDER@@. Do not store a string -- derive it.

    python3 -I make_reference.py

C10 compares each enumerated structure against the structure the paper DRAWS.
It is the only Tier-1 check that tests topology against external ground truth.
It only tests anything if the anchor was built by a DIFFERENT CODE PATH than the
enumerator.

That independence has been asserted twice and false twice:

  F-013  the stored anchor was a BYTE COPY of the enumeration -- same sha256,
         223 bytes -- while its note claimed independent fragment assembly
         "through the C3 oxygen (3-beta)" and the string was 3-alpha. C10
         reported `exact_match` while comparing a value to itself, on the
         flagship curation, for the one molecule it exists to protect. 537 rows
         of epicholesterol.
  F-037  the file written in round 1 TO FIX THAT had an AMINE_HALF template
         that WAS the enumeration's canonical SMILES with the sterol swapped
         for `%` and one ring digit renumbered -- byte-identical over 141
         characters. So the repair cured the half the bug was in and re-copied
         the other half.

PREFER A ROUTE THAT SHARES AS LITTLE AS POSSIBLE with the enumerator. The
working replacement uses plain string substitution with no RWMol / AddBond /
RemoveAtom anywhere, and takes the sterol VERBATIM from reagent_registry.json
so the provenance claim is mechanically true.

If the paper draws no full structure, write a row with an empty
reference_SMILES and a note saying so, so the ABSENCE of an anchor is visible.
"""
import csv
import os
import sys

LIB = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib")
sys.path.insert(0, LIB)

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   "reference_structures.csv")


def build():
    """Return [(compound_id, reference_SMILES, source, note)].

    `source` names the FIGURE OR SCHEME the structure was read from, e.g.
    "Fig 2a, lead compound, read head-first". `note` says what route built it.
    """
    raise SystemExit(
        "make_reference.py is still the generated stub. Read the structure the "
        "paper DRAWS -- usually the lead compound or the final product of a "
        "'precise synthesis' scheme -- and build it here by a route that "
        "shares as little as possible with the enumerator.")


if __name__ == "__main__":
    rows = build()
    with open(OUT, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["lipid_name", "reference_SMILES", "source_figure",
                    "read_on", "read_by", "expected_mismatch", "note"])
        w.writerows(rows)
    print(f"wrote {len(rows)} anchor(s) to {OUT}")
'''


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main(argv):
    if len(argv) < 3:
        print(__doc__)
        return 2
    folder = argv[1]
    inputs = argv[2:]
    dest = os.path.join(CURATION, folder)
    if os.path.exists(dest):
        print(f"REFUSING: {dest} already exists. Scaffolding over a live "
              "curation would overwrite its stage files. Pick another name, "
              "or work in the folder that is there.")
        return 1
    for p in inputs:
        if not os.path.exists(p):
            print(f"REFUSING: input not found: {p}")
            return 1

    now = datetime.now(timezone.utc).strftime("%Y-%m-%d")
    os.makedirs(os.path.join(dest, "00_pages"))
    os.makedirs(os.path.join(dest, "source"))

    # --- 00_acquisition.json -----------------------------------------
    acq = {"_what": "Where each input file came from, with its hash. A "
                    "provenance claim that is not executed is F-045; this file "
                    "is the executed half.",
           "scaffolded": now, "files": []}
    for p in inputs:
        tgt = os.path.join(dest, "source", os.path.basename(p))
        shutil.copy2(p, tgt)
        acq["files"].append({
            "path": os.path.relpath(tgt, ROOT),
            "sha256": sha256(tgt),
            "bytes": os.path.getsize(tgt),
            "retrieved_from": "TODO: the URL or DOI this came from",
            "retrieved_utc": "TODO",
            "licence": "TODO: the licence string on the landing page. "
                       "CLAUDE.md s8 requires every constituent source to be "
                       "CC0/CC-BY-compatible before corpus release, and that "
                       "is far cheaper to record now than to reconstruct "
                       "across 60-odd papers later.",
        })
    with open(os.path.join(dest, "00_acquisition.json"), "w") as fh:
        json.dump(acq, fh, indent=1)

    # --- 00_value_sources.json: THE INTAKE GATE ----------------------
    vs = {
        "_what": "THE INTAKE GATE. One entry per panel or table that supplies "
                 "a VALUE. Fill it BEFORE writing any digitiser.",
        "_why": "The two curations that took their values from a publisher "
                "Source Data file (NG_2025) and from the authors' own "
                "repository (WT_2026) have ZERO digitisation defects between "
                "them. The four that digitised have eight value-changing "
                "findings: F-006, F-008, F-017, F-018, F-019, F-020, F-033, "
                "F-040.",
        "_tiers": {k: v for k, v in VALUE_TIERS},
        "_rule": "A panel recorded as 5_pixels MUST fill "
                 "why_tiers_1_to_4_failed with a sentence a reviewer can "
                 "check. 'Not available' is not such a sentence.",
        "_also": "Once the numbers are in hand from tiers 1-3, use the figure "
                 "to VERIFY rather than to extract. NG_2025's stage06_heatmap "
                 "mapped all 310 drawn cells onto workbook values at rho >= "
                 "0.998 per panel, which independently confirmed block order, "
                 "row order and panel assignment -- strictly stronger than "
                 "anything available when the figure IS the source.",
        "panels": [{
            "panel": "TODO e.g. 'Fig 2b heatmap'",
            "supplies": "TODO what values, how many",
            "tier": "TODO one of " + ", ".join(k for k, _ in VALUE_TIERS),
            "locator": "TODO file name, sheet, URL, or page",
            "why_tiers_1_to_4_failed": "",
            "verified_against": "",
        }],
    }
    with open(os.path.join(dest, "00_value_sources.json"), "w") as fh:
        json.dump(vs, fh, indent=1)

    # --- 08_contradictions.csv ---------------------------------------
    with open(os.path.join(dest, "08_contradictions.csv"), "w",
              newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["claim", "locus_a", "locus_b", "what_we_emitted", "why"])

    # --- stage stubs -------------------------------------------------
    for fname, title, mods in STAGES:
        imports = ("import " + ", ".join(mods) + "  # noqa: F401\n"
                   if mods else "")
        classes = (", ".join(mods) + " already hold the audited versions."
                   if mods else "nothing geometric belongs in this stage.")
        with open(os.path.join(dest, fname), "w") as fh:
            fh.write(HEADER.format(title=title, fname=fname,
                                   imports=imports, classes=classes))

    with open(os.path.join(dest, "make_reference.py"), "w") as fh:
        fh.write(MAKE_REF.replace("@@FOLDER@@", folder))
    with open(os.path.join(dest, "reference_structures.csv"), "w",
              newline="") as fh:
        w = csv.writer(fh)
        # THIS HEADER IS C10's CONTRACT -- provenance/structure_check.py
        # reads these exact field names. The scaffold shipped
        # `compound_id,reference_SMILES,source,note` on 2026-11-10 and one
        # scaffolded folder then crashed C10 for the WHOLE repo with a
        # KeyError that named no file. Do not rename these.
        w.writerow(["lipid_name", "reference_SMILES", "source_figure",
                    "read_on", "read_by", "expected_mismatch", "note"])
        w.writerow(["", "", "", "", "", "", "SCAFFOLD: no anchor yet. C10 cannot test "
                                "topology against external ground truth until "
                                "make_reference.py builds one."])
    with open(os.path.join(dest, "DECISIONS.md"), "w") as fh:
        fh.write(DECISIONS.format(folder=folder, date=now))
    with open(os.path.join(dest, "PROVENANCE.md"), "w") as fh:
        fh.write(PROVENANCE.format(folder=folder, date=now))

    print(f"scaffolded curation/{folder}/")
    for f in sorted(os.listdir(dest)):
        print(f"  {f}")
    print("\nNEXT: fill 00_value_sources.json BEFORE writing any digitiser. "
          "Exhaust tiers 1-4 first; the two papers that never digitised a "
          "value have zero digitisation defects between them.")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
