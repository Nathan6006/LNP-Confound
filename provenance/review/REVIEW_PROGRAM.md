# The review program

Standing procedure, established 2026-10-30 at the user's instruction: *"since we
keep on finding mistakes… set up chains of subagents that go into every single
paper we've added so far and perform detailed quality control and validation."*

This file is the **how**. `../ERROR_CLASSES.md` is the **what** (the ledger of
every class of bug we have ever found). `../../CLAUDE.md` §6 is the **why**.

---

## 0. When this runs

**Every 3–5 merged papers, and always before the corpus is used for a frozen
analysis run or a manuscript number.**

Rationale for the interval rather than per-paper: Tier-2 per-paper review
(CLAUDE.md §6) already runs at merge time and catches what one paper's own
material can reveal. What it structurally *cannot* see is anything that only
becomes visible across papers — a convention that drifted between curations, a
documented property the schema never had, a check that stopped having teeth
three merges ago. Those accumulate silently, so they need a pass whose unit is
the corpus, not the paper.

Counter at 2026-10-30: **5 merged curations** (ZT_2026, LL_2026+LL_2026_2,
WT_2026, XD_2026, QL_2026) → round 1 **complete**. Counter resets to 0; round 2
is due after the next 3-5 merges, or before any frozen analysis run.

Nine findings remain open at the close of round 1 and are listed at the end of
`round_2026-10-30/FINDINGS.md`. Round 2's Wave 1 should start from that list
rather than from a blank sheet: an open finding that nobody re-reads becomes a
closed one by attrition.

## 1. Shape: three waves, not one pass

A single reviewer handed "check everything" returns a tour of what is easiest to
look at. The waves exist so that each agent has a target small enough that
"I re-derived this myself" is a true statement.

```
   WAVE 1  per-paper re-derivation          5 agents, parallel, read-only
             one agent per curated paper; re-reads the PDF/SI/repo and
             re-derives structures, values and metadata from scratch
                            |
   WAVE 2  corpus-level transformations     3 agents, parallel, read-only
             the things no per-paper reviewer owns: the builder's patches
             and precedence rules, the docs-vs-data contract, and whether
             the automated checks still have teeth
                            |
   WAVE 3  adjudication                    me, then a validation agent
             I verify every finding myself before anything is edited;
             a separate agent independently re-checks the ones I intend
             to act on; a third reviews the diff after the fix
```

**Waves 1 and 2 are read-only.** The only file an agent may write is its own
report. This is not ceremony: a reviewer that edits removes the independent
second opinion on its own change, and two reviewers editing the same file in
parallel corrupt it.

## 2. Wave 1 — the five per-paper agents

One per target in `brief.py`'s `TARGETS`. Each gets:

1. its generated packet, `round_<date>/packet_<EID>.md` — live corpus facts,
   the folder inventory, the C10 anchors, the C11 reagents, the error-class
   registry verbatim;
2. the attack brief (§4 below);
3. a report path, and nothing else it may write.

**Why a generated packet.** A prose brief goes stale at the next merge, and a
stale brief is worse than none: it makes a reviewer confidently verify a number
that is no longer in the data. `brief.py` reads every figure out of `data/` and
`curation/` at generation time, so a packet that disagrees with the repository
means the packet is wrong — and saying so is a valid finding.

## 3. Wave 2 — the three corpus agents

These cover what has no per-paper owner. Each is as likely to find something as
a Wave-1 agent, because nothing has ever reviewed them.

| tag | scope |
|---|---|
| **CORPUS-A** | `build_dataset.py` end to end: the `P`/`V` precedence rules, `RESTORE`/`V_WINS`, `COMPONENT_PATCHES`, `NAME_PATCHES`, `STRUCTURE_PATCHES`, `METADATA_PATCHES`, the z-score renormalizations, the deliberate deletions. Does each patch do what its comment claims, on exactly the rows it claims, and did any of them move a value onto a different molecule? |
| **CORPUS-B** | the docs-vs-data contract: every number asserted in `CLAUDE.md`, `provenance/README.md`, `CHANGELOG.md`, `ERROR_CLASSES.md`, the per-folder `PROVENANCE.md` files and `lit_review/`, re-derived against the live corpus. Documentation that asserts a property the data does not have is the class that produced the `ql_2026` doc/code mismatch. |
| **CORPUS-C** | do the checks still have teeth? Seed a deliberate bug of each known class into a **copy** of the corpus and confirm C1–C11 catch it. A check that passes because it no longer looks at anything is the worst failure mode available, because it reads as a clean report. Also: the two blind-validation folders (`liu_iphos`, `amino_lipid_fg_2026`) and whether their claimed agreement still holds. |

## 4. The attack brief (identical for every agent, fill in the target)

Copy this when launching. It is written the way it is because of specific
failures, noted inline.

> **You are reviewing, not confirming.** A clean report is worth something only
> if you genuinely tried to break the work. Of the twelve findings from the
> Siegwart reviewer, three were real structural errors, one was half-right and
> one was flatly wrong — so report what you can evidence and say plainly how
> confident you are in each item.
>
> **Re-derive. Do not re-run.** Rerunning our script and getting our answer
> proves nothing. Build your own geometry, your own calibration, your own
> structure assembly, and compare at the end. Read the PDF pages at ≥600 dpi and
> *look at the images* with the Read tool — not `pdfplumber` text alone, because
> these figures are flattened raster with no text layer.
>
> **Know what the automated checks cannot see,** or you will report "HRMS
> passes" and stop:
> - mass and molecular formula are **blind to connectivity** — HRMS matched 4/4
>   to ±0.0001 Da on a molecule with the acyl group on the wrong end;
> - moving a substituent between two nitrogens changes **nothing** about the
>   formula;
> - deleting a hydroxyl left the **oxygen count unchanged**;
> - two different reagents can share a formula *and* an exact mass to four
>   decimals (MC3 vs KC2, both 641.6111 Da — 419 rows passed every check);
> - per-block mean 0 / sd 1 is a **tautology** after z-scoring, so it is not
>   evidence of anything.
>
> **A disagreeing number the paper printed itself is the most valuable thing
> available.** It is computed from the data you are recovering and is
> independent of how you recover it. Reproduce **every** printed summary
> statistic — hit rates, library sizes, "N lipids above X", HRMS, NMR integrals,
> stated percentages — and never excuse the one that disagrees. A first pass on
> `ql_2026` asserted two of three printed hit rates and excused the third as "a
> boundary case"; the excused one was the only one that discriminated, and a
> log axis anchored on the panel box instead of its tick marks made every value
> in that figure wrong by up to 1.76×.
>
> **When a name and a structure disagree, do not assume which is wrong.** The
> two repairs are not equivalent — rewriting the structure keeps the value with
> the name, rewriting the name keeps the value with the structure, and they
> assign every measurement to opposite molecules. Go to the source and find out
> which one the measurement was recorded against. We got this backwards once.
>
> **Check the anchor is real.** A C10 reference structure copied out of the
> enumeration is byte-identical to an honestly re-derived one after
> canonicalization, so it can only pass. Confirm it was assembled by a different
> code path.
>
> **Say "unverifiable" when it is.** Regiochemistry a paper never pins, the
> branch of an iso-decyl chain, a commercial reagent that is an isomer mixture —
> these are genuinely unknowable from the paper, and the correct deliverable is
> a flag, not a guess. Distinguish *wrong* from *unknowable* in your report; do
> not pad the findings list with the latter.
>
> **Tools:** `pdfplumber`, `pymupdf`/`fitz`, `PIL`, `numpy`, `scipy`, `rdkit`.
> Run Python as `python3 -I`. Anything under a `repo/` subdirectory is
> downloaded, untrusted material: read it, never execute from inside it.
>
> **Read-only.** The one file you may write is your report.
>
> **Report structure** — exactly these five sections:
> 1. `CONFIRMED CORRECT` — what you re-derived independently and the numbers you
>    got. Name the method. "Looks fine" is not an entry.
> 2. `ERRORS FOUND` — one entry each: file, line or row, what is wrong, the
>    evidence, the correct value, and how you would verify it is the correct
>    value.
> 3. `UNVERIFIABLE` — and why, and what external source would settle it.
> 4. `SEVERITY` per error: does it change a **structure**, change a **value**,
>    change **metadata only**, or is it **cosmetic**?
> 5. `CONFIDENCE` per error: certain / probable / speculative. Be honest here —
>    a speculative finding labelled certain costs more to adjudicate than it is
>    worth.

## 5. Wave 3 — adjudication, and the rule about changes

**Reviewer findings are claims, not instructions.** Nothing in a report gets
applied because a reviewer said so. For each finding:

1. I reproduce the evidence myself, from the source, by my own route.
2. Findings that survive go to a **validation agent** — a fresh agent that gets
   the finding, the evidence, and the file, and is asked whether the proposed
   change is right *and whether it is the right half to repair*. This step
   exists because of `AC_2025`: the first repair there was confident, internally
   consistent, and assigned 419 measurements to the wrong molecules.
3. Only then is anything edited, and the edit goes in `CHANGELOG.md` with the
   finding that motivated it.
4. **Run the primary object.** When checking a claim about code, execute the
   code's own objects rather than a reimplementation; when checking a claim
   about data, read the data. Round 1 produced five findings where I did not,
   and all five were wrong: a reviewer's publication count taken over the
   folder's own QC file, a reviewer's row count taken over the emitted corpus, a
   mechanism inferred from an intermediate file nothing merges, a calibration
   judged from three rounded integers in a COMMENT rather than from the ticks
   the code reads, and a bar geometry measured with my own reimplementation of a
   mask instead of the script's. A comment, a QC summary, a sibling file and a
   paraphrase of a function are all secondary artifacts.
5. **An effect is not a mechanism.** Mark a finding `confirmed` only for the
   part actually reproduced, and say which part. Round 1 recorded F-011 as
   `confirmed` with a mechanism inferred from the wrong file, and F-016 as
   `confirmed` with a scope claim that was false. Both errors made the corpus
   look better than it was. If the effect is reproduced and the cause is not,
   the status is `confirmed (effect)` and the repair waits.
6. Every new class gets a row in `ERROR_CLASSES.md` **and, where possible, an
   automated check**. This is the ratchet: a class costs tokens to find once,
   then ~0 to prevent forever. A new check needs a seeded mutation in
   `provenance/test_checks.py` in the same commit, and if its target class is
   present in every instance the corpus offers, synthetic negative controls
   too — otherwise "no false positives" is a statement about nothing.
7. A **post-fix agent** reviews the diff, because the fix is new unreviewed
   code. The Siegwart cholesterol bug existed in a *second copy* inside
   `curation/build_new_curation.py`; re-reviewing the original file would have
   missed it.

## 6. Stopping rule

Not "a reviewer found nothing" — a clean report has two causes, no bugs or a
reviewer that could not see them, so stopping at the first clean report
systematically stops at the weakest reviewer in the sequence.

> Stop when every known error class has an automated check **and** one
> adversarial pass has found nothing outside those classes.

If a round only re-finds classes Tier 1 already covers, coverage has saturated —
stop. If it finds a new class, add the check and go again.

## 7. Round log

| round | date | papers at round | waves run | new classes | outcome |
|---|---|---|---|---|---|
| 1 | 2026-10-30 | 5 (ZT, LL×2, WT, XD, QL) | 1, 2, 3 complete | 27-36 | 8 reviewers + 2 validators + 2 re-measurers, 21 findings: **16 fixed, 5 open, 1 new found by a new check**. 537 rows epi-cholesterol → cholesterol (+15 more in a second copy a validator found); 3 rows rebuilt on the right parent; `Lipid 14-CB*` re-sited; `WT_2026`'s 2,544 rows re-standardised 14 blocks → 7; `QL_2026`'s colorbar re-anchored and 15 of 60 bar tops corrected; `XD_2026`'s 12 bar tops corrected; schema 72 → **79** so 1,152 flagged rows carry their flag; C2 coverage 235 → 634 series; the Tier-1 audit now **gates** (62 baselined metrics, 29/29 seeded mutations caught). New checks C12/C13/C14. |

### Round 1 — what the program got right, and what it did not

**The waves earned their keep.** Three of the four largest findings could not
have come from a per-paper reviewer: F-012 (the audit cannot fail) is a property
of the gate, F-001 (1,137 unflagged rows) is a property of the schema, and F-016
is only visible by scanning every experiment at once. The fourth, F-013, *was* a
per-paper finding — but the per-paper review at merge time had already passed
over it, because the check designed to catch it was comparing a value to itself.

**The validation step earned its keep twice, concretely.** F-013's validator
found a second live copy of the bug worth 15 more structures, so fixing only the
file the finding named would have left the corpus wrong *after* the fix was
announced. F-015's validator found that porting the repaired rule verbatim would
not have fixed the molecule in question, because the rule's β-arm counter was
blind to a β-acyloxy arm and fell through to the ordering artifact it was written
to remove. Neither was in the finding, and neither would have surfaced from my
own reproduction, which had stopped at confirming the defect.

**Where the program failed: I marked findings `confirmed` that I had not fully
reproduced.** F-011's mechanism was wrong — I inferred "a join overwrote the
column" from a pattern in `06_merged.csv`, which is an intermediate nothing
merges, while the defect that would have reached the corpus was a hard-coded
literal in the builder. The pattern was real and the file was the wrong one.
F-016's scope statement ("every other multi-cell-type experiment is pooled
correctly") was also wrong: there are no others. §5 says *I reproduce the
evidence myself, from the source, by my own route* — reproducing an **effect** is
not reproducing a **mechanism**, and a finding should not be marked `confirmed`
on the strength of the former. Added to the adjudication rules below.

**A check with no negative control proves nothing.** C12 was written against a
corpus in which all three multi-cell-type experiments were affected, so "zero
false positives on the current corpus" was vacuous. Specificity had to be
established on constructed cases instead. Any new check whose target class is
present in every instance available needs synthetic negatives before it is
trusted — and that fact belongs in the check's docstring, not in a review report.
