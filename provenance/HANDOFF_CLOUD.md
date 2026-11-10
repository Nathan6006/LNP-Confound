# Moving this run to the cloud — what transfers, what cannot, and how

Written 2026-11-10 when the laptop needed to be free. Read with
`AUTONOMOUS_RUN.md` (the plan), `RUN_LOG.md` (what happened, in order) and
`run_queue.json` (machine-readable state).

## The honest answer first

**A running Claude Code session cannot be migrated.** There is no "move this
session" operation: the conversation, its context and its in-process subagents
live on the machine that started them. What can move is **the work**, and after
commit `d13782c` all of it is on `origin/main`. A cloud session is therefore a
*fresh* session that picks the run up from the repo — not this one relocated.

## What is on the remote, and what deliberately is not

Pushed: **847 paths, 22 MB** — every extraction script, every `05_values.csv` /
`06_rows.csv`, every `07_qc.json`, every `DECISIONS.md`, `curation/lib/`,
`provenance/tools/`, the findings ledger and this run's log.

**Not pushed, on purpose:** the 84 publisher PDFs (368 MB) and 101 publisher
figure rasters / supplementary archives / JATS XML (60.2 MB), plus page renders
and annotated overlays. `.gitignore` already excluded `*.pdf` "for copyright";
this run extended it to every other asset a publisher serves, after measuring
that a plain `git add -A` would have published **60.2 MB in 101 files to a
PUBLIC repo** from papers that are subscription-only or **CC BY-NC-ND** — where
ND alone forbids redistribution. `CLAUDE.md` §8 requires CC0/CC-BY
compatibility before release, and several curations this run recorded their
licence as failing it.

**Consequence, stated plainly: a cloud session cannot re-read those PDFs.** It
has the values we extracted from them and the scripts that did it, but not the
pages themselves.

## So what can a cloud session actually do

**Fully, with no local dependency:**
- adjudicate the open findings (F-054 … F-061) — they are corpus questions,
  answered against `data/`, which is tracked
- run Tier-1 (`audit.py`, `check_docs.py`, C16/C17/C18) and `test_lib.py`
- decide the merge for the 23 folders that carry rows, and run
  `build_dataset.py` — every folder still has OPEN decisions, so the builder
  refuses them until a human resolves each one, by design
- the §8 release questions: licences, the excipient mapping, the batching blanks
- re-run the frozen pipeline (the open decision in §8), and the manuscript work

**Partly:**
- **new extractions of OPEN-ACCESS papers.** `provenance/tools/pmc_fetch.py`,
  `find_repo.py` and the Europe PMC `api/fulltextRepo` route all work from any
  machine with network access. The run's largest and cleanest result (5,112
  rows, zero digitisation) came from a GitHub deposit and would reproduce in the
  cloud exactly.

**Not at all:**
- any paper whose values need one of the 84 local PDFs
- any paper behind a publisher that requires the entitled browser session —
  `BROWSER_ACQUISITION.md` records which ones, measured. ACS is **not** covered
  by the institution in use, and Wiley/Elsevier/PNAS needed a browser warm-up.

## Two ways to run it, and when each is right

**(a) A fresh cloud session on the repo** — best when the laptop must be free.
Open `claude.ai/code` on `Nathan6006/LNP-Confound` and tell it to read
`provenance/AUTONOMOUS_RUN.md`, then `RUN_LOG.md` from the end backwards, then
`run_queue.json`. It inherits the rules, the tooling and the findings. It will
not inherit this session's conversation, so anything decided here but *not*
written to a file is lost — which is the reason `RUN_LOG.md` exists and is
written after every single paper.

**(b) Remote agents from a local session — MEASURED, AND IT DOES NOT WORK HERE.**
`isolation: "remote"` is documented as putting the *agent* in a remote cloud
environment. Probed on 2026-11-10 with a throwaway agent told to report its own
environment. It came back from:

```
hostname: Nathans-MacBook-Air-6.local
pwd:      /Users/nathanliu/LNP/confound/.claude/worktrees/agent-<id>
```

— i.e. **the flag silently degraded to a LOCAL git worktree on the same
laptop.** Nothing went to any cloud. The probe is still useful (it confirms a
fresh checkout at `d13782c` sees `data/data_vitro.csv` at 29,195 data rows,
`curation/lib/`, the queue, the log, network access to Europe PMC at HTTP 200,
and rdkit/pdfplumber/pymupdf/numpy/scipy all importing) — but as a route off
this machine it is a dead end, and an agent launched that way dies with the lid
exactly like an ordinary one.

**So there is no cloud path from inside this session.** Option (a) — a fresh
session started from the browser — is the only one.

## Keeping a LOCAL run alive across a closed lid

The failure already happened once this run, at 12:0x: the machine slept, all
four agents died, their task files froze at a 140-byte stub, and `check_state.py`
reported "state consistent" because it only checked that a task file *existed*.
The work survived on disk and all four were resumed. Two things came out of it:

- `check_state.py` now has a **staleness** test (15 min) that says *resume, not
  relaunch* — relaunching would discard the surviving work. Seeded: control 0,
  a 40-minute-stale agent exits 1 naming the folder and the age.
- To prevent it rather than detect it:
  `caffeinate -dimsu -w $$` for the session's lifetime, or run the session under
  `caffeinate -s`. On a Mac with the lid closed this still requires AC power;
  without it, macOS sleeps regardless and agents will die again. **The
  resumption path is the reliable half; caffeinate is the optimistic half.**

## The one rule a cloud session must not lose

Nothing in `curation/` has been merged. `data/` is unchanged at 29,195 x 79 /
1,803 x 79, and the frozen `results/` still belong to the 19,028-row corpus they
were computed on. **Do not merge, and do not re-point the frozen pipeline,
without the user's decision** — `CLAUDE.md` §2 and §7 both say so, and the row
counts move every reported number.
