# Session prompt — regenerate both campaigns from the released code

Paste this to start a **fresh session whose only job is to run the campaigns**. It is a
companion to `PROMPT.md`, not a replacement: read `code/README.md` first, then
`PROMPT.md` for the artefact's context, then this file for the run.

Delete this file before the release tag — it describes one night's work, not the package.

> `PROMPT.md` is stale in one section. "What the package is" still claims
> `placeholders.jl` is kept deliberately and the package is ~4 400 lines. Both were
> true on 23 Sep and are false now: `placeholders.jl`, `cutloop.jl` and
> `solverutils.jl` are deleted, `src/` is 2 818 lines, and no reference to
> `Equivalent`, `cc`, `bus_split` or `subbus` survives. `README.md` is current;
> trust it over `PROMPT.md` on anything about the code surface.

---

## Why this run exists

Every result on disk before it carries `git_sha 853ae50` — code that has since lost
35% of `src/`. The mathematics did not change and the acceptance test still passes
(IEEE-57 → 7.382, IEEE-118 → 5.29, both secure), but the **counters** did: removing the
`subbus` index changed a JuMP container from dense to sparse, which changed variable
ordering, which changed branching. IEEE-118 REF now takes **101 LP solves and 108
Benders iterations** where the recorded campaign says 79 and 87. IEEE-57 is unchanged
at 25 and 41.

So a reader running the released code would not reproduce the released per-cell
counters. This run regenerates every cell from the shipped working tree, and in doing
so establishes the new counter baseline.

## Hard constraints

1. **`docs/data` is read-only. Never write to it.** It holds the reference campaign the
   paper reports. New results go to `results/` only. If the new numbers deviate beyond
   noise, *report it and stop* — changing the reference data is a decision for Benoit,
   taken in discussion, not a step in this run.
2. **Nothing else runs on this machine while a campaign is in flight.** Wall time is a
   reported quantity: IEEE-118 REF measures 17 s alone and 27 s against a second job.
   That includes your own diagnostics — no test suites, no `julia -e` probes, no builds.
   The live watcher is the one permitted exception (it is a 20 s `sleep` loop).
3. **Exactly one live watcher.** Check for a stale one before starting, and kill it:
   `pgrep -fl live_report.jl`. Two watchers racing on `results/LIVE.md` has happened.
4. **No source edits during the run.** A campaign that spans a code change produces a
   mixed-`git_sha` result set, which is the exact failure this run exists to repair.
   Queue anything you notice; do not fix it now.
5. **Per-cell cap ≤ 1 hour.** Already true: `cap = 900` on Gurobi, `CAP_SCIP = 3000` on
   SCIP. Any individual diagnostic run you launch afterwards obeys the same bound.
6. **Do not commit, and do not `git add`.** `logs/*.log` are tracked and would bury the
   real diff in churn; a `logs/` gitignore is pending and is Benoit's call.
7. **Report what happened, including what failed.** A stage that dies must never be
   reported as finished — that has happened once here, and `set -euo pipefail` plus
   `failed_count` exist because of it.

## Complementary runs — these go FIRST

Two short runs produce claims the paper is waiting on, and one of them needs a source
edit. Both reasons put them **before** the campaigns: constraint 4 forbids editing source
mid-campaign, and constraint 2 forbids running anything alongside one. Do them in order,
one at a time, nothing else on the machine, then start the campaigns and leave them.

### Before either: a way to force a branch position

The released code has no way to pin a branch open or closed — the mechanism was stripped
in the extraction. Add the smallest thing that works: an optional dict of branch label to
position, applied in `src/master.jl` next to `warmstart_openings!`, defaulting to empty so
every existing call is untouched. **Commit nothing.** Make this edit before the campaigns
start, so that both stages run against one unchanging tree.

### Run 1 — 40-41 and 40-42 must be open

**Claim under test:** on the IEEE-118 reference instance every feasible plan opens both
branches, so de-energization is *necessary* here rather than merely permitted. That is a
sharper motivation than anything now in the paper.

- Reference instance and operating point of §5.1; keep the heuristic warm start and
  contingency embedding.
- Force 40-41 **closed** ($u_e = 1$). Set `H = 0` — which in `master.jl:51` means the
  Hamming constraint is never added, so the search is unrestricted, not frozen — and
  **relax the switchable set to all branches**. Cap 2 h.
- **Expected: infeasible.** Then repeat forcing 40-42 closed instead.
- **If run 1 exceeds 15 min, skip run 2.** One branch proves the point, and the paper
  sentence must then name only the branch that has a run behind it.
- If either run returns *feasible*, stop and report. The paper sentence cannot be written
  and the claim in the plan is wrong — do not soften it to fit.
- Fallback if unrestricted proves too slow: `H = 6`, switchable set still relaxed. Record
  which configuration produced the result; the paper must quote that one.

### Run 2 — IEEE-118 has no feasible base case at TLF = 100%

**Claim under test:** a sentence already standing in `results.tex` §5.4, currently
unverified. It explains why the experiments use TLF = 150% — the systems are stressed to
be a demanding but solvable challenge, not tuned until they solve.

- IEEE-118, `contingencies = []`, TLF = 100%. **Expected: infeasible.**
- If it is *feasible*, say so and stop: the manuscript sentence is then false and must be
  cut before the paper goes anywhere near a reviewer.

### For both

Store the results under `results/` like any other run, reference each script from
`README.md`, and say in the README what outcome each one demonstrates — a reviewer should
be able to read the conclusion without running anything. Gurobi only.

---

## The run

From `/Users/benoitjeanson/vsCode/TUD/IJEPES2026-OTSD/code`:

```bash
# 0. Nothing already running, and a clean tree apart from the known edits.
pgrep -fl 'live_report.jl|run_ablation.jl|julia' ; git status --short | head -40

# 1. The live view, first, one instance, in the background.
nohup julia --project=. experiments/live_report.jl --watch \
  > logs/$(date +%Y%m%d_%H%M%S)_live_report.log 2>&1 &

# 2. The campaigns. Archives the superseded cells to /tmp, then Gurobi, then SCIP,
#    stopping on the first stage that fails.
nohup ./run_night.sh > logs/$(date +%Y%m%d_%H%M%S)_night.log 2>&1 &
```

`run_night.sh` moves `results/*/` to `/tmp/otsd_results_superseded_<stamp>` before
anything else. That is deliberate and must not be skipped: a campaign skips cells
already on disk, so leaving them would silently preserve the superseded numbers.

Follow progress in `results/LIVE.md` (the full matrix, the in-flight cell, a
config-aware ETA, and a `vs reference` column read from `docs/data/index.json`) and in
`logs/stages.log` (one line per stage boundary).

## What to expect

| stage | cells | cap | expected |
|---|---|---|---|
| gurobi | 54 = 42 ablation (7 configs × 2 systems × 3 seeds) + 12 p_c | 900 s | ≈ 2.1 h |
| scip | 48 = 36 ablation (6 configs × 2 systems × 3 seeds) + 12 p_c | 3000 s | ≈ 3.5 h |

≈ 5.6 h in total. SCIP omits NO-LOCALSEARCH on purpose — it never converges, so its
cells would burn the cap and show nothing the Gurobi campaign does not.

Anchors that must hold, or something is wrong:

- IEEE-57 REF → 7.382, secure; IEEE-118 REF → 5.29, secure, 8 branches open.
- IEEE-57 REF counters 25 LP solves / 41 Benders iterations (unchanged by the cut).
- p_c: the PC-UNIF arm must reproduce REF exactly — it is the cheapest check that the
  weighting plumbing left the baseline alone. Expect IEEE-118 5.29 unweighted /
  6.649 weighted against PC-LEN's 5.43 / 6.611, and IEEE-57 7.382 / 3.011 in both arms.
- The three `118 NO-EMBED` SCIP cells, which died last night on a `SingularException`
  in `secured_dcpf`, must now finish. That is the one behavioural fix in this run:
  `_solve_angles` raises `DisconnectedTopology`, and `src/scip.jl` catches it inside
  the handler instead of letting it unwind through SCIP's C frame.
- Both stages exit 0, `failed_count == 0`, and `logs/stages.log` ends in `all done`.

A cell that reaches its cap is telling you something is wrong, not that it needed the
time: the slowest cell ever measured here is IEEE-118 NO-SCREEN on SCIP at 980 s.

## When it finishes

1. **Record the new counter baseline** — LP solves and Benders iterations per cell for
   seed 0, both systems, all configs, both backends. This replaces the void one.
2. **Compare against `docs/data/index.json`** (read-only): objective, secure, openings,
   wall time, counters. Say plainly which cells agree, which differ, and by how much.
   Objectives are the claim; counters are allowed to move, wall times are noisy.
3. **Write the outcome into the paper's revision plan**, §9 of
   `~/vsCode/TUD/papers/IJEPES2026OTSD/review/revision_plan.md`. A ⚠️ block there
   currently records that the IEEE-118 baseline is void and that this run supplies the
   replacement — resolve it with the measured numbers.
4. `results/PC.md` is written by the p_c report; re-render from disk with
   `julia --project=. experiments/run_pc.jl --report` if needed (it runs nothing).
5. Leave `results/` as it is and report. Do not commit.

## If a stage fails

The script stops at the first failure, which is intended. Then:

- read the tail of the stage log named in `logs/stages.log` — never `cat` a campaign log;
- the campaign resumes from what is on disk, so re-running the stage retries only the
  cells that are missing;
- known traps, all fixed but worth recognising: a Julia exception unwinding through
  SCIP's C frame kills the process a cell or two later (`_separate_count` returns `-1`
  instead); a `result.json` read mid-write raises `UnexpectedEOF` (results are written
  to `.tmp` and `mv`'d, and the watcher skips unparseable files); a licence check-out
  per model is slow and can fail under load on a token licence (Gurobi shares one env).

## Complementary runs

> **Benoit fills this in.** Anything added here runs **strictly after** both stages
> finish, one at a time, never overlapping a campaign or each other — the wall-time
> constraint above applies to these too. Each entry should say what to run, what result
> would confirm it, and what to do if it doesn't.

- [ ] …
