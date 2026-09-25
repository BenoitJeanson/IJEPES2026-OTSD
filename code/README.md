# OTSD — the code behind the paper

Optimal Transmission Switching with De-energization, solved by a Benders decomposition
with closed-form loss-of-load cuts. This is the implementation that produced every
number in the paper, reduced to what the experiments need.

```julia
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=.

julia> using OTSD
julia> r = solve_otsd("case118"; tlf = 1.5, H = 3, d_viol = 2, d_sol = 0, k = 3)
julia> r.objective, r.secure
```

`solve_otsd` starts from the fully closed network unless given a `warm_start`. The
published runs start from the expert heuristic of the companion paper; its output is
recorded in [`experiments/instances.jl`](experiments/instances.jl), because the
heuristic itself is not part of this package.

## Reproducing the paper

```julia
include("experiments/instances.jl")
inst, p = IEEE118, PAPER_OPERATING_POINT
r = solve_otsd(inst.case; tlf = inst.tlf, H = p.H, d_viol = p.d_viol, d_sol = p.d_sol,
               k = p.k, warm_start = inst.warm_start, embed = inst.embedded, seed = p.seed)
r.objective   # 5.29, per unit on a 100 MVA base
```

IEEE-57 at its own operating point returns 7.382. Both are checked by `Pkg.test()`.

The component ablation:

```bash
julia --project=. experiments/run_ablation.jl gurobi   # 42 cells: 7 configs × 2 systems × 3 seeds
julia --project=. experiments/run_ablation.jl scip     # the same, without a licence
```

The `p_c` sensitivity study, which asks whether a non-uniform contingency weighting
changes the outcome:

```bash
julia --project=. experiments/run_pc.jl          # 12 runs, then results/PC.md
julia --project=. experiments/run_pc.jl --report # re-render from disk
```

Results land one JSON per run under `results/`. A run already on disk is skipped, so
an interrupted campaign resumes where it stopped. Run campaigns **serially**: wall
time is one of the reported quantities, and contention distorts it — IEEE-118 at the
reference measured 17 s alone and 27 s against a second job.

## Two solvers, one algorithm

| | `GurobiBackend` | `SCIPBackend` |
|---|---|---|
| licence | commercial | none |
| cuts enter | from a lazy-constraint callback | from a constraint handler |
| | at integer-feasible nodes, inside one branch-and-cut tree | the same |
| reproduces the published numbers | yes, exactly | yes, exactly |

The decomposition wants one thing its solver may not offer: the ability to add a cut
the moment an integer-feasible topology appears, so that the whole search fits in a
single branch-and-cut tree. Gurobi provides it directly. SCIP has no lazy-constraint
callback — and `MOI.LazyConstraintCallback` is unsupported there, which is where a
first look stops — but it has the mechanism a lazy callback is a special case of: a
*constraint handler*, one level below MathOptInterface.

A handler is asked two questions during the search: `check`, is this candidate
solution acceptable, and `enforce_lp_sol`, it is not, so do something about it. That
is the Benders contract. `check` runs the separation without adding anything and
answers on the count; `enforce_lp_sol` runs it for real and reports `SCIP_CONSADDED`.
See [`src/scip.jl`](src/scip.jl).

Both backends drive the *same* separation routine (`_separate!` in
[`src/master.jl`](src/master.jl)), build cuts from the same `cut_constraint`
definitions, and receive them at the same points in the same kind of search. This is
the same algorithm on either solver, not a portable approximation of it. At the
published operating point it reproduces exactly on both systems — 7.382 and 5.29,
both secure, openings identical to the Gurobi run.

SCIP's MIP search is sequential, so expect it to be slower: measured serially here,
42 s against 12 s on IEEE-57 and 276 s against 19 s on IEEE-118.

What the two do *not* share is how a tie is broken. The objective does not price an
opening that carries no flow, so a topology can be optimal in several ways, and the
two solvers enumerate a degenerate restricted problem in different orders. In the
ablation configurations this shows: every one reaches 7.382, secure, but four of the
six reach it with one or two extra branches open. The reference configuration —
the one the paper reports — lands on the same openings on both solvers.

One detail that is not bookkeeping: the handler declares **variable locks** for every
branch variable. Without them SCIP's dual presolve concludes that nothing constrains
those variables, fixes them at a bound, and returns a wrong answer *silently*. The two
`misc/allow{strong,weak}dualreds` flags go off for the same reason.

### HiGHS is here, but not as a master

The feasibility cut is read off the Farkas dual of an infeasible subproblem LP, and
SCIP does not expose one through MathOptInterface. HiGHS does, so a SCIP master pairs
with HiGHS subproblems (`lp_backend` in [`src/backend.jl`](src/backend.jl)) — with
**presolve off**, since a presolve-proved infeasibility leaves no basis to read a
certificate from. Both are open source, and the pairing is invisible to the
algorithm: the cuts, and the order they are generated in, are unchanged.

HiGHS cannot host the master. `kHighsCallbackMipDefineLazyConstraints` exists in the
enum as an undocumented placeholder, but `HighsCallbackDataIn` carries no channel for
a new row and no C function adds one. `HiGHSLPBackend` is therefore not a master
backend, and the master rejects it rather than silently solving without cuts.

## What is here, and what is not

About 3 700 lines. The research monorepo this was extracted from carries roughly
13 000, and the difference is all work that belongs to other papers: substation
reconfiguration, a Dantzig–Wolfe decomposition, a magnitude search, network
equivalents, and the plotting stack.

| file | |
|---|---|
| `grid.jl`, `case.jl` | the graph model and PGLib loading |
| `dcpf.jl`, `pocket.jl` | DC power flow, the security analysis, pockets and bridges |
| `blocks.jl` | the model: variables, KCL, Ohm, limits, energization |
| `master.jl` | the master problem and the separation routine |
| `subproblem.jl` | the contingency LP and its duals |
| `cutpool.jl`, `cutsink.jl` | the two cut families, and how a cut is delivered |
| `scip.jl` | the constraint handler that takes the place of the lazy callback |
| `sbs.jl`, `phase.jl`, `session.jl` | the switchable branch set and the local search |
| `backend.jl` | what differs between solvers |

Nothing here is carried for another line of work. Earlier versions of this package
kept substation configurations and network equivalents as empty types that every model
builder accepted and every use guarded with an emptiness test; they are gone, along
with the `subbus` index on the energization variables, which existed only so that a
bus could be split. What remains is the model the paper describes, with no argument
that the published runs never set.

### The network has no parallel circuits

A branch is keyed on its pair of buses, so where PGLib lists parallel circuits on a
corridor the corridor becomes a single branch: IEEE-57 goes from 80 branch records to
78 branches, IEEE-118 from 186 to 179. This is deliberate, and it is discussed in the
paper. Collapsing a double circuit can only increase the number of branches whose
outage severs a pocket, never reduce it.

## Licence

Apache 2.0. See [`../LICENSE`](../LICENSE).
