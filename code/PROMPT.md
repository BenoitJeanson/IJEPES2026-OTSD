# Session prompt — the OTSD reproducibility artefact

Paste this to start a session on the code, the solvers, or the experiments.

---

You are working on the **reproducibility artefact** for the journal paper
IJEPES-D-26-03028, *"Improved Transmission Switching with De-energization via a Benders
Decomposition with Closed-Form Loss-of-Load Cuts"* (Jeanson, Tanneau, Tindemans),
currently under major revision.

Working directory: `/Users/benoitjeanson/vsCode/TUD/IJEPES2026-OTSD` — a **public**
GitHub repo (`BenoitJeanson/IJEPES2026-OTSD`, Apache 2.0) that already served the
interactive campaign viewer at `docs/`. The package lives in `code/`.

Read `code/README.md` first, then this file.

## Why it lives here and not with the paper

The manuscript repo (`papers/IJEPES2026OTSD`) carries `review/` — reviewer reports and
a revision plan containing candid assessments and reviewers' names. It can never be
published, and a Zenodo DOI would archive it permanently. This repo was already public
and already held the campaign data. One repo, one DOI, code beside the data it
produced. (Decision D22 in the paper's `review/revision_plan.md`; §9 there covers the
artefact from the paper's side.)

## What the package is

A **frozen one-way extraction** from the research monorepo `~/vsCode/TUD/tnr`:
~4 400 lines against 13 000, 16 dependencies against 28. Removed: bus splitting
(`MagnitudeSearch`), Dantzig–Wolfe, `ZoneScreen`, network equivalents,
`ReduceViolations`, the Makie stack, `PowerSystems`, the expanding-H search, the light
master, and ~250 lines of unreachable callback helpers.

**Kept deliberately:** empty `SubstationConfs` and `EquivalentSet` in
`placeholders.jl`. Every model builder takes them and every use is guarded by an
emptiness test; removing the arguments means editing model-building code, which would
risk the reproduction for no gain. The feeder/bus indirection *was* stripped, because
its empty path is a plain equality.

tnr is untouched and will not be refactored to depend on this. Divergence is accepted.

## The acceptance test — the thing that matters

`Pkg.test()` runs it, and nothing ships unless it passes:

| system | objective | operating point |
|---|---|---|
| IEEE-57 | **7.382** | H=3, d_viol=2, d_sol=0, k=3, seed 0 |
| IEEE-118 | **5.29** | same |

Exact reproduction, both secure, openings identical to the recorded campaign. The warm
start is the expert heuristic's output from the companion paper, recorded as data in
`experiments/instances.jl` — the heuristic itself is not in this package, and starting
anywhere else changes the numbers.

When refactoring, also check the **internal counters** (LP solves, Benders iterations),
not just the objective: equality there is what proves a refactor changed nothing.

## Solver backends — this was settled the hard way

The algorithm needs cuts injected at integer-feasible nodes inside one branch-and-cut
tree. Three backends exist in `backend.jl`; `supports_lazy` distinguishes them.

- **`GurobiBackend`** — lazy-constraint callback. The published configuration.
- **`SCIPBackend`** — **the same algorithm, no licence.** SCIP has no lazy callback but
  has the mechanism that generalises one: a *constraint handler* (`scip.jl`). `check`
  runs separation dry and answers on the count; `enforce_lp_sol` runs it for real and
  reports `SCIP_CONSADDED`. Verified: 7.382 (65 s) and 5.29 (297 s), openings identical
  to Gurobi.
- **`HiGHSBackend`** — a *different* algorithm: classical solve-separate-resolve
  (`cutloop.jl`), because HiGHS genuinely cannot take a lazy constraint from Julia.

**The user's instruction: the artefact must run the same algorithm as Gurobi, because
anything else would have to be explained in the paper and that is unwanted.** SCIP
satisfies this. **HiGHS and `cutloop.jl` are therefore slated for removal** — that work
is not yet done. Check with the user before deleting, then remove the backend, the
loop, the violation-filter branch that exists only for it, and the README section.

### Why HiGHS cannot, definitively

HiGHS 1.15.1 has nine callback types (upstream `docs/src/callbacks.md`); lazy
constraints are not among them. `HighsCallbackDataOut` gives the incumbent and a
read-only cut pool; `HighsCallbackDataIn` accepts only an interrupt flag and a user
solution. No C function adds a row. `kHighsCallbackMipDefineLazyConstraints` exists in
the header as an undocumented placeholder with no data channel. Also checked: **GLPK**
supports lazy constraints but its MIP solver would collapse on this master (and it is
GPL); **SCIP's `MOI.LazyConstraintCallback` is unsupported** — the constraint handler
is the right layer, which is one level below where a first look stops.

### Four things in the SCIP path are load-bearing and non-obvious

1. **Variable locks in `lock`.** Without them SCIP's dual presolve fixes the branch
   variables and returns a *wrong answer silently* — the toy model returned 0.0 instead
   of 12.
2. **Both `misc/allowstrongdualreds` and `misc/allowweakdualreds` off.** SCIP split the
   old single flag the docs still name.
3. **The handler cannot be called `benders`** — SCIP's own native Benders framework
   already registers that name. Ours is `otsdbenders`.
4. **Subproblem LPs run on HiGHS** (`lp_backend` in `backend.jl`). The feasibility cut
   is read off the Farkas dual of an infeasible LP, and SCIP does not expose one
   through MathOptInterface. HiGHS does — but only with **presolve off**, since a
   presolve-proved infeasibility leaves no basis to read a certificate from.

## Experiments

```bash
julia --project=. experiments/run_ablation.jl gurobi   # 42 cells: 7 configs x 2 systems x 3 seeds
julia --project=. experiments/live_report.jl --watch    # regenerates results/LIVE.md every 20 s
```

Results are one `result.json` per run under `results/`; a cell already on disk is
skipped, so campaigns resume after interruption.

**Run experiments strictly serially.** Wall time is a reported quantity and contention
distorts it — 118 REF measured 17 s alone and 27 s contended. A whole night was lost
once to a chained script that dead-locked: it polled `pgrep -f "run_ablation.jl gurobi"`,
which matched the shell that had *written* the script. Sequence stages inside one
script; never poll for a process name.

**Open:** the campaign definition still names HiGHS cells and should move to SCIP once
that backend is confirmed as the licence-free one. A clean serial Gurobi campaign was
running when this prompt was written.

## What the ablation shows (contended run, indicative)

Every configuration reaches the reference objective except NO-LOCALSEARCH. On IEEE-118:
screening is worth ~180x in LP solves, contingency embedding ~10x in wall time, the
closed-form cuts ~15x in LP solves, cut inheritance is real but modest. NO-LOCALSEARCH
is the only one that costs *quality*, and its three seeds straddle the reference on
IEEE-57 (7.331 / 8.532 / 7.361 against 7.382), none converged — which is the concrete
evidence that the outer loop is heuristic.

## Working style

Discuss before implementing anything substantial. Verify claims against the code or a
run rather than asserting them — several conclusions here reversed on inspection
(HiGHS "might" support callbacks; SCIP "doesn't" support lazy constraints; a
one-seed reading of NO-LOCALSEARCH). Test the mechanism on a toy model before wiring it
into the master; that is how each of the four SCIP gotchas was found cheaply.
