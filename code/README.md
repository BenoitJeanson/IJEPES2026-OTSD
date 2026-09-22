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
julia --project=. experiments/run_ablation.jl gurobi
```

Results land one JSON per run under `results/`. A run already on disk is skipped, so
an interrupted campaign resumes where it stopped.

## Two solvers

| | `GurobiBackend` | `HiGHSBackend` |
|---|---|---|
| licence | commercial | none |
| cuts enter | from a lazy-constraint callback, inside one branch-and-cut tree | between master solves |
| reproduces the published numbers | yes, exactly | no — see below |

The decomposition wants one thing its solver may not offer: the ability to add a cut
the moment an integer-feasible topology appears. Gurobi provides it, and the whole
search then fits in a single tree. HiGHS does not — `kHighsCallbackMipDefineLazyConstraints`
exists in the enum, but `HighsCallbackDataIn` carries no channel for a new row and the
Julia wrapper exposes none. The HiGHS backend therefore solves the master to
optimality, separates against its solution, adds whatever that solution violates, and
solves again, until a round finds nothing to add.

Both backends drive the *same* separation routine (`_separate!` in
[`src/master.jl`](src/master.jl)) and build cuts from the same `cut_constraint`
definitions, so the cuts cannot drift apart. What differs is when they arrive, and
that changes the search: expect different topologies where the restricted problem has
ties, and expect the licence-free path to be substantially slower, because every
round discards the tree and re-solves from scratch.

One consequence worth stating plainly: a cut that is valid but not *violated* at the
current point makes no progress. Inside a tree the solver absorbs it harmlessly; in a
loop it must be filtered out, or every round re-adds the same inert cuts and the loop
cannot tell convergence from deadlock. `cut_violation` in
[`src/cutpool.jl`](src/cutpool.jl) is that filter.

## What is here, and what is not

About 4 400 lines. The research monorepo this was extracted from carries roughly
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
| `cutloop.jl` | Benders without callbacks |
| `sbs.jl`, `phase.jl`, `session.jl` | the switchable branch set and the local search |
| `backend.jl` | what differs between solvers |

Two type families in `placeholders.jl` are defined but never populated:
substation configurations and network equivalents. The model builders accept them
because they are shared with those other lines of work, and every use is guarded by
an emptiness test. They are kept rather than removed because taking the arguments out
means editing model-building code, and the value of this package is that its model is
the published one.

### The network has no parallel circuits

A branch is keyed on its pair of buses, so where PGLib lists parallel circuits on a
corridor the corridor becomes a single branch: IEEE-57 goes from 80 branch records to
78 branches, IEEE-118 from 186 to 179. This is deliberate, and it is discussed in the
paper. Collapsing a double circuit can only increase the number of branches whose
outage severs a pocket, never reduce it.

## Licence

Apache 2.0. See [`../LICENSE`](../LICENSE).
