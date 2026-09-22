# ── Solver backends ───────────────────────────────────────────────────────────
#
# The decomposition needs one thing from its solver that is not portable: a way to
# add a cut once an integer-feasible topology appears. Gurobi offers it as a
# lazy-constraint callback, so the whole search fits in one branch-and-cut tree.
# HiGHS does not — `kHighsCallbackMipDefineLazyConstraints` exists in the enum, but
# `HighsCallbackDataIn` carries no channel for a new row and the Julia wrapper
# exposes none. On HiGHS the master is therefore re-solved between rounds of cut
# generation.
#
# The separation itself is shared: see `separate_cuts` in `master.jl`. Only the
# moment at which a cut reaches the master differs between the two.

"""
    Backend

How the master problem is solved. Concrete backends are [`GurobiBackend`](@ref) and
[`HiGHSBackend`](@ref); [`supports_lazy`](@ref) distinguishes them.
"""
abstract type Backend end

"""
    GurobiBackend(; threads = 4)

The algorithm as published: a single branch-and-cut tree with cuts injected from a
lazy-constraint callback. Requires a Gurobi licence.
"""
Base.@kwdef struct GurobiBackend <: Backend
    threads::Int = 4
end

"""
    HiGHSBackend(; threads = 4, max_rounds = 200)

A licence-free path. Cuts are added between master solves rather than inside the
tree, so the master is re-solved once per round until a round produces no cut. The
cuts are the same; the search that finds them is not, and neither are the runtimes.
`max_rounds` bounds the loop.
"""
Base.@kwdef struct HiGHSBackend <: Backend
    threads::Int = 4
    max_rounds::Int = 200
end

"Can this backend accept a cut from inside the branch-and-cut tree?"
supports_lazy(::GurobiBackend) = true
supports_lazy(::HiGHSBackend) = false

"Short identifier used in run tags and result records."
backend_name(::GurobiBackend) = "gurobi"
backend_name(::HiGHSBackend) = "highs"

# ── Model construction ────────────────────────────────────────────────────────

"""
    new_model(backend; log_path = "") -> JuMP.Model

An empty model on `backend`, silent on the console and logging to `log_path` when one
is given. Gurobi uses a direct model over a shared environment: a licence check-out
per model is slow, and on a token licence it can fail under load.
"""
function new_model(b::GurobiBackend; log_path::String = "")
    isassigned(GRB_ENV_REF) || (GRB_ENV_REF[] = Gurobi.Env(Dict{String,Any}()))
    m = direct_model(Gurobi.Optimizer(GRB_ENV_REF[]))
    set_optimizer_attribute(m, "LogToConsole", 0)
    set_optimizer_attribute(m, "Threads", b.threads)
    isempty(log_path) || set_optimizer_attribute(m, "LogFile", log_path)
    m
end

function new_model(b::HiGHSBackend; log_path::String = "")
    m = direct_model(HiGHS.Optimizer())
    set_optimizer_attribute(m, "output_flag", false)
    set_optimizer_attribute(m, "threads", b.threads)
    isempty(log_path) || set_optimizer_attribute(m, "log_file", log_path)
    m
end

# ── Attributes that both solvers have under different names ───────────────────

set_seed!(::GurobiBackend, m, seed::Int) = set_attribute(m, "Seed", seed)
set_seed!(::HiGHSBackend, m, seed::Int) = set_attribute(m, "random_seed", seed)

set_timeout!(::GurobiBackend, m, seconds::Real) = set_attribute(m, "TimeLimit", Float64(seconds))
set_timeout!(::HiGHSBackend, m, seconds::Real) = set_attribute(m, "time_limit", Float64(seconds))

"Emphasis on finding good incumbents early. Gurobi only; HiGHS has no equivalent."
set_mip_focus!(::GurobiBackend, m, focus::Int) = set_attribute(m, "MIPFocus", focus)
set_mip_focus!(::HiGHSBackend, _, _) = nothing

"""
    request_infeasibility_certificate!(backend, m)

Ask for a Farkas certificate when the LP turns out infeasible — the feasibility cut
is read off it. Gurobi must be told in advance; HiGHS returns a dual ray whenever it
has one, so there is nothing to set.
"""
request_infeasibility_certificate!(::GurobiBackend, m) = set_optimizer_attribute(m, "InfUnbdInfo", 1)
# HiGHS returns a dual ray only from the simplex solve itself: if presolve proves
# infeasibility first there is no basis to read a certificate from, and the dual
# objective comes back empty. Presolve is therefore off on the subproblem LPs.
request_infeasibility_certificate!(::HiGHSBackend, m) = set_optimizer_attribute(m, "presolve", "off")

"Announce that the model will receive lazy constraints. Gurobi requires this up front."
enable_lazy!(m) = MOI.set(m, MOI.RawOptimizerAttribute("LazyConstraints"), 1)

"""
    solver_version(backend) -> String

For the run manifest, so a result can be tied to the solver that produced it.
"""
function solver_version(::GurobiBackend)
    major, minor, technical = Ref{Cint}(), Ref{Cint}(), Ref{Cint}()
    Gurobi.GRBversion(major, minor, technical)
    "$(major[]).$(minor[]).$(technical[])"
end

solver_version(::HiGHSBackend) = unsafe_string(HiGHS.Highs_version())

"""
    default_backend() -> Backend

`GurobiBackend` when a licence is present, `HiGHSBackend` otherwise, so the examples
run for a reader without a commercial solver.
"""
default_backend() = isassigned(GRB_ENV_REF) ? GurobiBackend() : HiGHSBackend()
