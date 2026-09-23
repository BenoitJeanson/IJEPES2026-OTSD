"""
    OTSD

Optimal Transmission Switching with De-energization: a Benders decomposition with
closed-form loss-of-load cuts.

This package is the reproducibility artefact for the paper. It is a frozen extraction
of the research code that produced the published numbers, reduced to what the
experiments need — the graph model, the DC security analysis, the master problem, the
two cut families, and the local search around them.

Two solver backends are available, described in `src/backend.jl`:

  * `GurobiBackend` — the algorithm as published: one branch-and-cut tree, cuts
    injected from a lazy-constraint callback on every incumbent.
  * `HiGHSBackend`  — a licence-free path. HiGHS cannot accept lazy constraints from
    Julia, so the master is re-solved between rounds of cut generation.

Both drive the same separation routine, so the cuts are identical by construction;
only the moment at which they enter the master differs.

Entry point: [`solve_otsd`](@ref).
"""
module OTSD

using Dates
using Logging
using LinearAlgebra
using Random
using Serialization
using SparseArrays

using Graphs
using MetaGraphsNext
using JuMP
using PGLib
using Printf

using Gurobi
using HiGHS
using SCIP


# ── Graph model ───────────────────────────────────────────────────────────────
export ELabel, VLabel, PGLibtograph, scale_branch_limits!, balance!
export getbridges, edge_distance_map, hamming, connectedcomponent, neighbor_labels
export incident, incident_signed, opposite, openbranchesset, sub_graph
export labels, edge_labels

# ── Cases ─────────────────────────────────────────────────────────────────────
export GridCase, ElementaryCase, RichCase, load_case

# ── Security analysis ─────────────────────────────────────────────────────────
export reduced_sa, secured_dcpf, create_bridge_to_pocket, violating_contingencies

# ── Benders ───────────────────────────────────────────────────────────────────
export BendersCut, CutRecord, mastercutpool, solve_benders_phase
export CutSink, LazySink, DirectSink, benders_cut_loop!

# ── Local search ──────────────────────────────────────────────────────────────
export SessionState, run_benders_iterations!, sa_induced_followed, extend_sbs_by_hops

# ── Backends ──────────────────────────────────────────────────────────────────
export Backend, GurobiBackend, SCIPBackend, HiGHSBackend, solve_otsd
export backend_name, solver_version, supports_lazy

# Gurobi is held in a single environment for the life of the session: a licence
# check-out per model is slow and, on a token licence, can fail under load.
const GRB_ENV_REF = Ref{Gurobi.Env}()

function __init__()
    try
        GRB_ENV_REF[] = Gurobi.Env(Dict{String,Any}())
    catch e
        @warn "No Gurobi environment; the HiGHS backend remains available." exception = e
    end
end

include("grid.jl")
include("gridcase.jl")
include("elementarycase.jl")
include("placeholders.jl")
include("pocket.jl")
include("case.jl")

include("dcpf.jl")

include("backend.jl")
include("cutsink.jl")
include("blocks.jl")
include("benders_commons.jl")
include("cutpool.jl")
include("cutloop.jl")
include("scip.jl")
include("subproblem.jl")
include("master.jl")

include("sbs.jl")
include("phase.jl")
include("session.jl")
include("solve.jl")

end # module
