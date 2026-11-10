#!/usr/bin/env python3
"""END-OF-TURN GUARD. Run this before ending ANY turn in the autonomous run.

WHY IT EXISTS. On 2026-11-10 I ended a turn with the sentence "Launching batch
2 on the three largest remaining chemistry screens" and did not launch
anything. The user had to notice and say so. In an autonomous run that is the
worst possible failure: the session looks busy, the log reads correctly, and no
work is happening.

A promise in prose is not a launched agent. This script compares what the queue
CLAIMS against what is actually on disk, and prints a single verdict line.

    python3 -I check_state.py

Rules enforced:
  1. Every queue entry marked `in_progress_agent` must have a live task file.
  2. If nothing is in progress and the queue still has `ready_to_extract`
     papers, that is an IDLE RUN -- the run is supposed to be continuous.
  3. Every folder claimed `extracted` must actually have 05_values.csv.
"""
import glob, json, os, sys, time

ROOT = "/Users/nathanliu/LNP/confound"
# An agent writing nothing for this long is treated as dead.
# Agents here routinely run 20-45 min, so this must exceed a
# slow tool call but stay under a sleep that goes unnoticed.
STALE_SECONDS = 15 * 60
# Session-scoped, so it is an OVERRIDE, not a stored input: this tool now lives
# in the repo (CLAUDE.md §7 forbids a repo script READING an input from a
# scratchpad). Nothing tracked depends on it -- it only makes the agent-liveness
# cross-check possible in the session that is running the batch.
TASKS = os.environ.get("CONFOUND_TASK_DIR") or (
    "/private/tmp/claude-501/-Users-nathanliu-LNP-confound/"
    "c9d997e8-37d9-4383-83e3-9eec1d2b80dd/tasks")


def claimed_running_without_an_agent(state, task_dir):
    """A queue entry marked in_progress MUST name the agent that is doing it.

    WHY. On 2026-11-10 I set a folder's status to `in_progress_agent` and only
    then launched the agent -- and in between, this very script printed
    "running: songx_39664425" for an agent that did not exist. The queue is
    the run's own account of itself, so a status written before the thing it
    describes makes every later reading of that file wrong, and the failure
    reads as a clean line. The user's standing instruction after an earlier
    instance of exactly this was: "This should never happen again. Ensure that
    it doesn't."

    So the ORDER IS NOW A RULE, enforced here rather than remembered:
    launch the agent, take its id from the launch result, THEN write the queue
    entry with `agent_task_id` set. An entry claiming to be in progress with no
    id, or with an id that has no task file on disk, is a lie about the state of
    the run and fails.
    """
    bad = []
    for x in state.get("queue", []):
        if x.get("status") != "in_progress_agent":
            continue
        tid = x.get("agent_task_id")
        if not tid:
            bad.append((x["folder"], "marked in_progress with no agent_task_id "
                                     "-- was it written BEFORE the launch?"))
        else:
            f = os.path.join(task_dir, f"{tid}.output")
            if not os.path.exists(f):
                bad.append((x["folder"],
                            f"agent_task_id {tid} has no task file"))
                continue
            # A STALE TASK FILE IS NOT A LIVE AGENT. On 2026-11-10 the laptop
            # slept; all four agents died, their task files froze at a 140-byte
            # stub, and this script printed "state consistent" because it only
            # checked that the file EXISTED. That is the same defect the whole
            # run has been chasing -- a check that cannot see the condition it
            # is for reports a clean line. Their WORK survived on disk, so the
            # fix is to notice and resume, not to relaunch.
            try:
                age = time.time() - os.path.getmtime(f)
                size = os.path.getsize(f)
            except OSError:
                continue
            if age > STALE_SECONDS:
                bad.append((x["folder"],
                            f"agent {tid} has not written for {int(age/60)} min "
                            f"({size} bytes) -- probably dead. RESUME it with "
                            f"SendMessage (its work is on disk); do not "
                            f"relaunch, which would discard that work."))
    return bad

def main():
    q = json.load(open(os.path.join(ROOT, "provenance/run_queue.json")))
    rows = q["queue"]
    inprog = [x for x in rows if x.get("status") == "in_progress_agent"]
    ready = [x for x in rows if x.get("status") == "ready_to_extract"]
    done = [x for x in rows if x.get("status") == "extracted"]

    # agent task files touched in the last 6 hours = plausibly this run's
    live = []
    for p in glob.glob(os.path.join(TASKS, "*.output")):
        try:
            if time.time() - os.path.getmtime(os.path.realpath(p)) < 6 * 3600:
                live.append(os.path.basename(p).split(".")[0])
        except OSError:
            pass

    problems = []
    for x in done:
        f = os.path.join(ROOT, "curation", x["folder"], "05_values.csv")
        if not os.path.exists(f):
            problems.append(f"{x['folder']} is marked extracted but has no "
                            f"05_values.csv")
    if not inprog and ready:
        problems.append(f"IDLE: 0 agents in progress but {len(ready)} papers "
                        f"are ready to extract. The run is meant to be "
                        f"continuous -- launch, do not announce.")

    print(f"queue: {len(rows)} papers | in_progress {len(inprog)} | "
          f"ready {len(ready)} | extracted {len(done)} | "
          f"recent agent task files {len(live)}")
    for x in inprog:
        print(f"  running: {x['folder']}")
    if problems:
        print("\nSTATE PROBLEMS:")
        for p in problems:
            print(f"  * {p}")
        return 1
    orphans = claimed_running_without_an_agent(q, TASKS)
    if orphans:
        print("\nQUEUE CLAIMS WORK THAT HAS NO AGENT:")
        for folder, why in orphans:
            print(f"  * {folder}: {why}")
        print("\nLaunch first, read the agent id from the launch result, THEN "
              "write the queue entry. A status written before the launch makes "
              "every later reading of run_queue.json wrong.")
        return 1
    print("state consistent")
    return 0


if __name__ == "__main__":
    sys.exit(main())
