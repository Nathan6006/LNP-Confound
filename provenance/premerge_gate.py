#!/usr/bin/env python3
"""THE PRE-MERGE GATE. Every staged curation row is checked HERE before it can
reach data/.

WHY IT EXISTS. `merge_curations()` routed each row with one expression:

    (vitro if r["Model"] == "in_vitro" else vivo).append(row)

so a row whose `Model` is blank, misspelled, or absent does not fail -- it goes
to **data_vivo.csv**. Measured 2026-11-10 on the staged run: `lij_42130331` has
**300 rows with a blank `Model`**, which would have entered the in-vivo pool and
become **16.6% of it** (1,803 rows today) without one line of output saying so.
The same sweep found **14 of 26** staged files off the 79-column schema.

That is this project's recurring shape -- a defect that produces a clean report
-- and the merge was the one place with no check at all. `open_decisions()`
gates on a FOLDER's unresolved questions; nothing gated the ROWS.

The design rules, each paid for by an earlier finding:

  * **Report every violation, not the first.** A gate that stops at violation
    one turns a single fix-and-rerun cycle into fourteen.
  * **Never canonicalise before testing.** C7's blank check stripped before
    testing and so could not see 915 padded cells (F-022, error class 37). The
    checks here see the bytes as emitted.
  * **An empty input is a failure, not a pass.** `find_screens`, `batch_gate`
    and the first two versions of C16 all reported clean while examining
    nothing. `check()` refuses a file with zero rows.
  * **The vocabulary comes from schema_quality.json**, which is the single
    source of truth; it is not restated here (F-001).

Run standalone over the staged folders:

    python3 -I provenance/premerge_gate.py                 # every curation/*/
    python3 -I provenance/premerge_gate.py <file.csv> ...  # named files
"""
import csv, glob, json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

QUALITY_SCHEMA = json.load(open(os.path.join(HERE, "schema_quality.json")))
# {column: set(allowed)} for the enum-valued quality columns only. `allowed:
# null` means free text or float and is deliberately not constrained.
VOCAB = {c["name"]: set(c["allowed"])
         for c in QUALITY_SCHEMA["columns"] if c.get("allowed")}

# Routing is the whole point of this file: these two strings, exactly, and
# nothing else. Not a regex, not a startswith, not a case fold -- a row that
# does not say which pool it belongs in has not been curated, it has been
# guessed at.
MODELS = {"in_vitro", "in_vivo"}

# Columns without which a row cannot be interpreted at all.
REQUIRED = ["Experiment_ID", "Model", "IL_SMILES", "Experiment_value"]


def _num(x):
    try:
        float(x)
        return True
    except (TypeError, ValueError):
        return False


def check(path, known_experiments=None, parse_smiles=True, in_merge=False):
    """Return a list of violation strings for one emitted-rows CSV.

    `in_merge` says whether this file is actually listed in
    build_dataset.py's MERGE. A file that is not yet in MERGE is STAGED: its
    violations are things to clear before it may be added, not things blocking
    today's build. Saying otherwise would make the gate cry wolf on 30 folders
    that nobody has proposed merging, and a gate that always fires trains
    people to ignore it -- the same reasoning that fixed the phantom OPEN entry
    in `open_decisions()`.
    """
    v = []
    name = os.path.relpath(path, ROOT)
    try:
        with open(path, encoding="utf-8", errors="replace", newline="") as f:
            rows = list(csv.DictReader(f))
    except Exception as e:
        return [f"{name}: UNREADABLE ({type(e).__name__}: {e})"]

    if not rows:
        # A zero-row file passing silently is how a scanner that cannot see its
        # subject announces success -- but a DELIBERATE zero is also a real
        # result this run produced 10 of, and flagging those as violations is
        # noise. So it blocks only if someone has actually put it in MERGE.
        if in_merge:
            return [f"{name}: ZERO ROWS but it is listed in MERGE. An empty "
                    f"file must not be a merge input: either it yields rows or "
                    f"it comes out of MERGE."]
        return []

    cols = set(rows[0].keys())
    for c in REQUIRED:
        if c not in cols:
            v.append(f"{name}: missing required column {c!r}")
    if any(c for c in cols if c is None):
        v.append(f"{name}: ragged header -- a row has more fields than the header")

    bad_model, blank_model = {}, 0
    blank_val, bad_val, bad_smi, blank_smi = 0, 0, [], 0
    bad_vocab = {}
    eids = set()

    for i, r in enumerate(rows, start=2):   # line 1 is the header
        m = r.get("Model")
        if m is None or str(m).strip() == "":
            blank_model += 1
        elif m not in MODELS:
            bad_model.setdefault(repr(m), []).append(i)

        ev = r.get("Experiment_value")
        if ev is None or str(ev).strip() == "":
            blank_val += 1
        elif not _num(ev):
            bad_val += 1

        s = r.get("IL_SMILES")
        if s is None or str(s).strip() == "":
            blank_smi += 1
        elif parse_smiles:
            from rdkit import Chem, RDLogger
            RDLogger.DisableLog("rdApp.*")
            if Chem.MolFromSmiles(str(s)) is None:
                bad_smi.append(i)

        if r.get("Experiment_ID"):
            eids.add(r["Experiment_ID"])

        for col, allowed in VOCAB.items():
            if col not in r:
                continue
            val = r[col]
            if val is None:
                continue
            if str(val) not in allowed:
                bad_vocab.setdefault((col, repr(val)), 0)
                bad_vocab[(col, repr(val))] += 1

    if blank_model:
        v.append(f"{name}: {blank_model} row(s) have a BLANK `Model`. These "
                 f"would be routed to data_vivo.csv by the else-branch and "
                 f"enter the in-vivo pool silently.")
    for m, lines in bad_model.items():
        v.append(f"{name}: {len(lines)} row(s) have Model={m}, which is neither "
                 f"'in_vitro' nor 'in_vivo' (first at line {lines[0]}). Routed "
                 f"to vivo by the else-branch.")
    if blank_val:
        v.append(f"{name}: {blank_val} row(s) have a BLANK Experiment_value. "
                 f"The corpus column is empty, so this curation has not "
                 f"standardised its values yet -- check whether they are "
                 f"sitting in a raw column such as "
                 f"`Experiment_value_raw_as_printed` instead.")
    if bad_val:
        v.append(f"{name}: {bad_val} row(s) have a NON-NUMERIC "
                 f"Experiment_value (present but unparseable as a float)")
    if blank_smi:
        v.append(f"{name}: {blank_smi} row(s) have a blank IL_SMILES")
    if bad_smi:
        v.append(f"{name}: {len(bad_smi)} row(s) have an unparseable IL_SMILES "
                 f"(first at line {bad_smi[0]})")
    for (col, val), n in sorted(bad_vocab.items()):
        v.append(f"{name}: {col}={val} on {n} row(s) is not in "
                 f"schema_quality.json's vocabulary {sorted(VOCAB[col])}. "
                 f"Policy P5: a genuinely new value may be added BY HAND after "
                 f"review; a script may not add one, and two different things "
                 f"must never be merged under one token.")

    if known_experiments:
        mine = os.path.basename(os.path.dirname(os.path.abspath(path)))
        for e in sorted(eids):
            owner = known_experiments.get(e)
            # Two files in ONE folder sharing an Experiment_ID is normal and
            # correct: `siegwart_40644303` emits 09_rows_vitro.csv and
            # 09_rows_vivo.csv for one screen, exactly as BL_2023 and KZ_2016
            # already span both pools in data/. Only a tag reused by a
            # DIFFERENT PAPER is a collision. The first version of this check
            # reported three such pairs as collisions -- a false positive in
            # the gate, found by reading its own output.
            if owner and owner[0] != mine:
                v.append(f"{name}: Experiment_ID {e!r} is already used by "
                         f"{owner[1]}. Policy P4(b): resolve by publication date, "
                         f"earlier takes 'a' and later takes 'b' "
                         f"(e.g. {e}a / {e}b). A bare colliding tag is never "
                         f"emitted.")
    return v


def discover():
    """Emitted-rows files across the staged curations.

    A file counts if it has BOTH an `Experiment_value` and an `IL_SMILES`
    column -- by SCHEMA, not by filename. Keying on the name `06_rows.csv` is
    exactly the mistake F-076 made: 14 of 29 folders use another name, and the
    scan declared ~200 verified rows unmergeable because of it.
    """
    out = []
    for p in sorted(glob.glob(os.path.join(ROOT, "curation", "*", "*.csv"))):
        try:
            with open(p, encoding="utf-8", errors="replace", newline="") as f:
                hdr = next(csv.reader(f), [])
        except Exception:
            continue
        if "Experiment_value" in hdr and any("SMILES" in (c or "") for c in hdr):
            out.append(p)
    return out


def merge_paths():
    """The files build_dataset.py would actually read, resolved to real paths."""
    try:
        import importlib.util
        spec = importlib.util.spec_from_file_location(
            "_bd", os.path.join(HERE, "build_dataset.py"))
        # Importing build_dataset executes module-level work; read MERGE out of
        # the source instead so this gate stays side-effect free.
    except Exception:
        pass
    out = set()
    src = open(os.path.join(HERE, "build_dataset.py"), encoding="utf-8").read()
    i = src.find("MERGE = [")
    if i < 0:
        return out
    block = src[i:src.find("]", i) + 1]
    for frag in block.split('"'):
        if frag.endswith(".csv"):
            cand = os.path.normpath(os.path.join(HERE, "..", frag)) \
                if not os.path.isabs(frag) else frag
            out.add(os.path.normpath(cand))
            out.add(os.path.normpath(os.path.join(ROOT, frag)))
    return out


def main():
    paths = sys.argv[1:] or discover()
    if not paths:
        print("NO EMITTED-ROWS FILE FOUND. Refusing to print an empty report "
              "as if it were a clean one.")
        return 1
    known = {}
    inm = merge_paths()
    problems, staged = [], []
    for p in paths:
        is_m = os.path.normpath(os.path.abspath(p)) in inm
        v = check(p, known_experiments=known, in_merge=is_m)
        (problems if is_m else staged).extend(v)
        try:
            with open(p, encoding="utf-8", errors="replace", newline="") as f:
                for r in csv.DictReader(f):
                    if r.get("Experiment_ID"):
                        known.setdefault(r["Experiment_ID"],
                                         (os.path.basename(os.path.dirname(
                                             os.path.abspath(p))),
                                          os.path.relpath(p, ROOT)))
        except Exception:
            pass
    nm = sum(1 for p in paths
             if os.path.normpath(os.path.abspath(p)) in inm)
    print(f"pre-merge gate: {len(paths)} emitted-rows file(s) examined "
          f"({nm} in MERGE, {len(paths) - nm} staged)")
    if staged:
        print(f"\n{len(staged)} ISSUE(S) in STAGED files -- not blocking "
              f"today's build, but each must be cleared before that file may "
              f"be added to MERGE:\n")
        for x in staged:
            print(f"  - {x}")
    if problems:
        print(f"\n{len(problems)} VIOLATION(S) in files that ARE in MERGE -- "
              f"the build must refuse these:\n")
        for x in problems:
            print(f"  * {x}")
        return 1
    print("\nno violations in any file listed in MERGE")
    return 0


if __name__ == "__main__":
    sys.exit(main())
