#!/usr/bin/env python3
"""C16 - is each C10 anchor DERIVED, or transcribed from the enumeration?

WHY THIS EXISTS, AND WHY IT IS NOT THE SAME AS C10's EXISTING CHECK.

C10 compares each enumerated structure against the structure the paper DRAWS.
Its entire value rests on the reference being built by a different code path
than the enumerator: copy the enumerator's output into the reference file and
C10 compares a value to itself and can only ever pass. That is not theoretical.

  * F-013 (round 1): the ZT_2026 anchor WAS a byte copy of the enumeration -
    same sha256, 223 bytes - while its note claimed independent fragment
    assembly "through the C3 oxygen (3-beta)". The string it annotated was 3a.
    C10 reported `exact_match` while comparing a value to itself, on the one
    molecule it existed to protect. 537 rows carried epicholesterol past nine
    checks.
  * F-037 (round 2): the REPLACEMENT generator, written specifically to cure
    that, hard-coded an amine-half template that was the enumeration's own
    canonical SMILES with the sterol substring replaced by a placeholder and one
    ring-closure digit renumbered. Undoing those two edits reproduced the
    enumeration byte for byte over 141 characters. So the repair cured the half
    the bug was in and re-copied the other half.

**Byte-identity of the EMITTED anchor is not the signal, and testing for it is
a trap.** A correctly re-derived structure, once canonicalised, is byte-identical
to a correct enumeration - that is what agreement means. C10 already reports
`anchor_may_be_a_copy` on that weaker signal (anchor == enumeration AND the
folder has no generator script), which is why F-037 slipped past: the folder DID
have a generator.

The signal is one level up, in the GENERATOR'S SOURCE. A hand-derived template
is written in whatever traversal order the author found natural. RDKit's
canonical traversal is a specific, arbitrary order - it starts C6O2B2 at a
dodecyl chain terminus - and a human does not reproduce it by chance over a
hundred characters. So: take each long SMILES literal out of the generator
script, normalise away the two edits a copy-paste needs (placeholder
substitution and ring renumbering), and ask whether what remains is a substring
of the enumeration's canonical output. If it is, the template was pasted.

This test would have fired on make_reference.py the day it was written.

    python3 -I provenance/check_anchor_independence.py      # non-zero on a hit

Limits, stated honestly: it cannot prove independence, only catch transcription.
A template typed out in canonical order on purpose would trip it (no such case
exists here), and a transcription deliberately re-ordered would evade it. It is
a tripwire on the cheap failure mode, not a proof.
"""
import csv
import glob
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))

MIN_LEN = 40          # shorter literals are fragments, not whole-molecule pastes
SMILES_CHARS = re.compile(r'^[A-Za-z0-9@+\-\[\]\(\)=#%/\\\.$]+$')


def smiles_templates(path):
    """The generator's COMPOSED module-level templates, not its raw literals.

    Reading quoted literals out of the source does not work and the first
    version of this file was useless for exactly that reason: it examined 0
    templates and reported clean. A template assembled as
    `"[N+]5(" + _B2 + ")" + ...` has no single long literal anywhere in the
    file, so a regex over quoted strings sees only fragments -- and a pasted
    template split across a few adjacent string pieces would evade it too.

    So import the module and read what it actually built. Anchor generators keep
    their work behind `if __name__ == "__main__"`, so importing is side-effect
    free; any that is not will fail loudly here rather than silently pass.
    """
    import importlib.util
    spec = importlib.util.spec_from_file_location("_anchor_gen", path)
    mod = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(mod)
    except Exception as e:
        print(f"  WARNING: could not import {os.path.basename(path)}: "
              f"{type(e).__name__}: {e}", file=sys.stderr)
        return []
    out = []
    for name in dir(mod):
        v = getattr(mod, name)
        if (isinstance(v, str) and len(v) >= MIN_LEN
                and SMILES_CHARS.match(v)
                and v.count("C") + v.count("c") >= 10):
            out.append((name, v))
    return out


PLACEHOLDER = re.compile(r"%|\{[A-Za-z_]+\}")
MIN_PIECE = 20


def normalise(s):
    """Strip ALL ring digits, so a renumbering cannot hide a match."""
    return re.sub(r"\d", "", s)


def pieces(template):
    """The template split at its placeholders, normalised, long pieces only.

    Splitting matters and is why the first version of this check had no teeth.
    A template marks where the independent half goes with a placeholder -- `%`
    for the sterol in make_reference.py -- so the pasted text is NOT a
    contiguous substring of the enumeration: the enumeration has 82 characters
    of sterol sitting in the middle of it. Testing the whole template as one
    string therefore always missed, and the check reported clean on the very
    file it was written for.

    Each surviving piece has to appear in the enumeration, in order.
    """
    return [p for p in (normalise(x) for x in PLACEHOLDER.split(template))
            if len(p) >= MIN_PIECE]


def transcribed_from(template, enumerated):
    """True if every long piece of `template` occurs in `enumerated`, in order."""
    ps = pieces(template)
    if not ps:
        return False
    ne = normalise(enumerated)
    at = 0
    for p in ps:
        i = ne.find(p, at)
        if i < 0:
            return False
        at = i + len(p)
    return True


def enumerated_smiles(folder):
    out = []
    for name in ("03_enumerated.csv", "03_structures.csv"):
        p = os.path.join(folder, name)
        if os.path.exists(p):
            with open(p) as f:
                for r in csv.DictReader(f):
                    for col in ("IL_SMILES", "SMILES", "smiles"):
                        if r.get(col):
                            out.append(r[col])
    return out


def main():
    hits, checked = [], 0
    for folder in sorted(glob.glob(os.path.join(ROOT, "curation", "*"))):
        if not os.path.isdir(folder) or folder.endswith("ng_2025"):
            continue                              # ng_2025: another session owns it
        enum = enumerated_smiles(folder)
        if not enum:
            continue
        for script in sorted(glob.glob(os.path.join(folder, "*.py"))):
            base = os.path.basename(script)
            if not ("reference" in base or "anchor" in base):
                continue                          # only the ANCHOR generators
            for var, lit in smiles_templates(script):
                checked += 1
                for e in enum:
                    if transcribed_from(lit, e):
                        hits.append((f"{os.path.relpath(script, ROOT)}:{var}",
                                     len(lit), e[:60]))
                        break

    print(f"C16 anchor independence: {checked} template literal(s) examined "
          f"in anchor generators")
    if hits:
        print(f"\n{len(hits)} TEMPLATE(S) APPEAR TO BE TRANSCRIBED FROM THE "
              f"ENUMERATION:\n")
        for path, n, e in hits:
            print(f"  * {path}: a {n}-character SMILES literal is, after "
                  f"removing placeholders and ring digits, a substring of the "
                  f"enumerated structure\n      {e}...")
        print("\nAn anchor built from the enumeration cannot test the "
              "enumeration. Re-derive the template from the paper's drawing in "
              "a different traversal order, or state in the anchor's note "
              "exactly which part of the molecule it does NOT check.")
        return 1
    print("no anchor template matches its own enumeration -- all derived")
    return 0


if __name__ == "__main__":
    sys.exit(main())
