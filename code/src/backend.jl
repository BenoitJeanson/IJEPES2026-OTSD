# ── Solver backends ───────────────────────────────────────────────────────────
#
# The decomposition needs one thing from its solver that is not portable: a way to
# add a cut once an integer-feasible topology appears, so that the whole search fits
# in a single branch-and-cut tree. Gurobi offers it directly, as a lazy-constraint
# callback. SCIP has no such callback but has the mechanism a callback is a special
# case of — a constraint handler (`scip.jl`) — and reaches the same place.
#
# Both master backends therefore run the *same* algorithm. The separation itself is
# shared code either way (`_separate!` in `master.jl`); only the interface through
# which a cut reaches the solver differs.
#
# HiGHS appears here in one role only: it solves the contingency LPs under a SCIP
# master, because SCIP does not expose a Farkas certificate. It cannot host the
# master — `kHighsCallbackMipDefineLazyConstraints` exists in the enum, but
# `HighsCallbackDataIn` carries no channel for a new row and the Julia wrapper
# exposes none.

"""
    Backend

Which solver runs a problem. The master backends are [`GurobiBackend`](@ref) and
[`SCIPBackend`](@ref), and both take cuts inside one branch-and-cut tree
([`supports_lazy`](@ref)). [`HiGHSLPBackend`](@ref) is not a master backend: it
solves the contingency LPs for a SCIP master.
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
    HiGHSLPBackend(; threads = 4)

Not a master backend. HiGHS solves the contingency LPs under a [`SCIPBackend`](@ref)
master, because the feasibility cut is read off the Farkas dual of an infeasible LP
and SCIP does not expose one through MathOptInterface. Passing this to the master
is rejected — HiGHS cannot take a cut from inside a branch-and-cut tree.
"""
Base.@kwdef struct HiGHSLPBackend <: Backend
    threads::Int = 4
end

"""
    SCIPBackend(; threads = 1)

The published algorithm without a licence. SCIP has no lazy-constraint callback, but
it has the mechanism that generalises one — a constraint handler — so cuts still
enter at integer-feasible nodes inside a single branch-and-cut tree. See `scip.jl`.
SCIP's MIP search is sequential; `threads` reaches its LP solves only.
"""
Base.@kwdef struct SCIPBackend <: Backend
    threads::Int = 1
end

"""
    supports_lazy(backend) -> Bool

Can this backend take a cut at an integer-feasible node, inside one tree? True of
every backend that may host the master; the master rejects one for which it is false.
"""
supports_lazy(::GurobiBackend) = true
supports_lazy(::SCIPBackend) = true
supports_lazy(::HiGHSLPBackend) = false

"Short identifier used in run tags and result records."
backend_name(::GurobiBackend) = "gurobi"
backend_name(::SCIPBackend) = "scip"
backend_name(::HiGHSLPBackend) = "highs-lp"

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

function new_model(b::SCIPBackend; log_path::String = "")
    m = direct_model(SCIP.Optimizer())
    set_optimizer_attribute(m, "display/verblevel", 0)
    set_optimizer_attribute(m, "lp/threads", b.threads)
    m
end

function new_model(b::HiGHSLPBackend; log_path::String = "")
    m = direct_model(HiGHS.Optimizer())
    set_optimizer_attribute(m, "output_flag", false)
    set_optimizer_attribute(m, "threads", b.threads)
    isempty(log_path) || set_optimizer_attribute(m, "log_file", log_path)
    m
end

# ── Attributes the two masters have under different names ─────────────────────
#
# Set on the master only, so there is nothing to define for the LP backend.

set_seed!(::GurobiBackend, m, seed::Int) = set_attribute(m, "Seed", seed)
set_seed!(::SCIPBackend, m, seed::Int) = set_attribute(m, "randomization/randomseedshift", seed)

set_timeout!(::GurobiBackend, m, seconds::Real) = set_attribute(m, "TimeLimit", Float64(seconds))
set_timeout!(::SCIPBackend, m, seconds::Real) = set_attribute(m, "limits/time", Float64(seconds))

"""
    lp_backend(backend) -> Backend

Which solver solves the contingency subproblems. They are pure LPs, and the
feasibility cut is read off the dual ray of an infeasible one, so the only
requirement is that the solver hands back a Farkas certificate.

Gurobi does. SCIP does not expose one through MathOptInterface, so a SCIP master
pairs with HiGHS subproblems — both open source, and invisible to the algorithm: the
cuts and the order they are generated in are unchanged.
"""
lp_backend(b::GurobiBackend) = b
lp_backend(b::HiGHSLPBackend) = b
lp_backend(b::SCIPBackend) = HiGHSLPBackend(threads = b.threads)

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
# SCIP's subproblem LPs are solved through the same path; it returns a dual ray
# without being asked, as HiGHS does.
request_infeasibility_certificate!(::SCIPBackend, m) = nothing
request_infeasibility_certificate!(::HiGHSLPBackend, m) = set_optimizer_attribute(m, "presolve", "off")

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

solver_version(::SCIPBackend) = string(SCIP.SCIPmajorVersion(), ".", SCIP.SCIPminorVersion(),
                                       ".", SCIP.SCIPtechVersion())
solver_version(::HiGHSLPBackend) = unsafe_string(HiGHS.Highs_version())

"""
    default_backend() -> Backend

`GurobiBackend` when a licence is present, `SCIPBackend` otherwise, so the examples
run for a reader without a commercial solver — and run the same algorithm.
"""
default_backend() = isassigned(GRB_ENV_REF) ? GurobiBackend() : SCIPBackend()
