#!/usr/bin/env python3
"""Update run_queue.json entries BY FOLDER, and REFUSE when a folder is unknown.

WHY THIS EXISTS. Every queue update in this run was an inline loop of the shape

    for x in state["queue"]:
        if x["folder"] == "some_folder":
            x["status"] = "..."

which does nothing at all when the folder name is wrong, and says nothing about
it. On 2026-11-10 I launched an agent against PMID 41555650, briefed it to work
in `curation/dahlmanje_41555650/`, and updated that key. The real entry is
`mitchellmj_41555650` -- **the Mitchell lab, not Dahlman** -- so the update was a
silent no-op and `check_state.py` showed 3 agents where 4 were running. Worse,
the agent had been told to create the folder if absent, so a folder named for
the wrong lab would have appeared in the tree, and `Experiment_ID` is derived
from the first author.

A no-op that reports success is the same defect class as the whole run has been
chasing in other people's code: `find_tables()` returning 0 on a page that is
nothing but table, `find_screens` hiding a column, `batch_gate` checking zero
structures, `triage` printing an empty report, three ZT_2026 marginal checks
scoring a literal 0 for twelve days. **A write that matches nothing must fail.**

    python3 -I provenance/tools/queue_update.py <folder> key=value [key=value...]
    python3 -I provenance/tools/queue_update.py --find 41555650

Values are parsed as JSON when possible, else kept as strings. `--find <text>`
searches folder and pmid and prints matches, which is the cheap way to get the
key right BEFORE launching, not after.
"""
import json
import os
import sys

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                     "..", ".."))
QUEUE = os.path.join(ROOT, "provenance", "run_queue.json")


def load():
    with open(QUEUE) as f:
        return json.load(f)


def find(state, text):
    t = str(text).lower()
    return [x for x in state["queue"]
            if t in x["folder"].lower() or t in str(x.get("pmid", "")).lower()]


def main(argv):
    state = load()
    if not argv:
        print(__doc__)
        return 2
    if argv[0] == "--find":
        if len(argv) < 2:
            print("--find needs a search term")
            return 2
        hits = find(state, argv[1])
        if not hits:
            print(f"no queue entry matches {argv[1]!r}")
            return 1
        for x in hits:
            print(f"  folder={x['folder']!r}  pmid={x.get('pmid')}  "
                  f"status={x.get('status')}  rows={x.get('rows_emitted')}  "
                  f"task={x.get('agent_task_id')}")
        return 0

    folder, assignments = argv[0], argv[1:]
    hits = [x for x in state["queue"] if x["folder"] == folder]
    if not hits:
        near = find(state, folder.split("_")[-1])
        print(f"REFUSING: no queue entry has folder == {folder!r}.\n"
              f"An update that matches nothing is a silent no-op, and this "
              f"exact mistake made check_state report 3 running agents where 4 "
              f"were running, with a folder named for the wrong lab.", file=sys.stderr)
        if near:
            print("\nDid you mean one of these?", file=sys.stderr)
            for x in near:
                print(f"  {x['folder']!r}  (pmid {x.get('pmid')})", file=sys.stderr)
        else:
            print("\nNothing similar either. Run with --find <pmid> first.",
                  file=sys.stderr)
        return 1
    if len(hits) > 1:
        print(f"REFUSING: {len(hits)} entries share folder {folder!r}. "
              f"The queue should have one entry per folder.", file=sys.stderr)
        return 1

    entry, changed = hits[0], []
    for a in assignments:
        if "=" not in a:
            print(f"REFUSING: {a!r} is not key=value", file=sys.stderr)
            return 1
        k, v = a.split("=", 1)
        try:
            v = json.loads(v)
        except Exception:
            pass
        changed.append(f"{k}: {entry.get(k)!r} -> {v!r}")
        entry[k] = v
    import datetime
    state["updated_utc"] = datetime.datetime.now(
        datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    with open(QUEUE, "w") as f:
        json.dump(state, f, indent=1)
    print(f"{folder}:")
    for c in changed:
        print(f"  {c}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
