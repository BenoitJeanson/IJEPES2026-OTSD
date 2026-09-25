# ── SCIP: the same algorithm, without a licence ───────────────────────────────
#
# Gurobi takes a cut from a lazy-constraint callback. SCIP has no such callback, but
# it has the more general mechanism the callback is a special case of: a *constraint
# handler*. A handler is asked two questions during the search —
#
#   check(sol)          is this candidate solution acceptable?
#   enforce_lp_sol()    it is not; do something about it
#
# — which is exactly the Benders contract. `check` runs the separation without adding
# anything and answers on the count; `enforce_lp_sol` runs it for real and reports
# `SCIP_CONSADDED`. Cuts therefore enter at integer-feasible nodes inside one
# branch-and-cut tree, as they do on Gurobi, and the algorithm is the same one.
#
# `lock` is not bookkeeping. A handler that declares no variable locks tells SCIP's
# dual presolve that nothing constrains those variables, and it will fix them at a
# bound — silently returning a wrong answer rather than an error. Every variable the
# handler can constrain is locked in both directions.

"""
    ConshdlrSink(m, o, sol, added, dry)

A [`CutSink`](@ref) reading the solution SCIP is currently asking about. During
enforcement `sol` is `C_NULL`, which SCIP reads as "the solution under enforcement";
during `check` it is the candidate being judged.
"""
struct ConshdlrSink <: AccumulatingSink
    m::Model
    o::SCIP.Optimizer
    sol::Ptr{SCIP.SCIP_SOL}
    added::Base.RefValue{Int}
    dry::Bool
end

ConshdlrSink(m::Model, o::SCIP.Optimizer, sol = C_NULL; dry::Bool = false) =
    ConshdlrSink(m, o, sol, Ref(0), dry)

solution_value(s::ConshdlrSink, var) =
    SCIP.SCIPgetSolVal(s.o, s.sol, SCIP.var(s.o, JuMP.index(var)))

"""
    BendersConshdlr

The handler SCIP consults. `separate!` is the solver-independent separation routine
of `mastercutpool`; `vars` are the branch variables, locked so that presolve leaves
them alone.
"""
mutable struct BendersConshdlr <: SCIP.AbstractConstraintHandler
    o::SCIP.Optimizer
    m::Model
    separate!::Function
    vars::Vector{MOI.VariableIndex}
end

"""
    _separate_count(ch, sol; dry) -> Int

Run separation against `sol`, adding cuts unless `dry`; return how many were violated.

Nothing may throw out of here. SCIP calls the handler from its own C frame, and a
Julia exception unwinding through it is undefined behaviour — in practice it corrupts
SCIP's state and the process dies a cell or two later, which is how a campaign came to
report "all done" over a stage that had crashed. A candidate we cannot evaluate is
reported as unacceptable instead, which is both safe and true: `-1` is returned so the
caller reads it as violated, and SCIP rejects the solution or branches elsewhere.

`DisconnectedTopology` is the expected case — SCIP asks about candidates that have not
yet cleared `base_connectivity!`. Anything else is a real bug, so it is logged once
rather than swallowed silently.
"""
function _separate_count(ch::BendersConshdlr, sol; dry::Bool)
    sink = ConshdlrSink(ch.m, ch.o, sol; dry)
    try
        ch.separate!(sink)
    catch e
        e isa DisconnectedTopology || @error "separation failed on a SCIP candidate" exception = (e, catch_backtrace())
        return -1
    end
    cuts_added(sink)
end

function SCIP.check(ch::BendersConshdlr, constraints, sol, checkintegrality,
                    checklprows, printreason, completely)
    _separate_count(ch, sol; dry = true) == 0 ? SCIP.SCIP_FEASIBLE : SCIP.SCIP_INFEASIBLE
end

function _enforce(ch::BendersConshdlr)
    n = _separate_count(ch, C_NULL; dry = false)
    # `CONSADDED` promises SCIP a constraint it can make progress on. A failed
    # separation added none, so it must report plain infeasibility instead and let
    # SCIP branch — claiming otherwise sends it looking for progress that never comes.
    n < 0 && return SCIP.SCIP_INFEASIBLE
    n == 0 ? SCIP.SCIP_FEASIBLE : SCIP.SCIP_CONSADDED
end

SCIP.enforce_lp_sol(ch::BendersConshdlr, constraints, nusefulconss, solinfeasible) =
    _enforce(ch)

SCIP.enforce_pseudo_sol(ch::BendersConshdlr, constraints, nusefulconss,
                        solinfeasible, objinfeasible) = _enforce(ch)

function SCIP.lock(ch::BendersConshdlr, constraint, locktype, nlockspos, nlocksneg)
    for vi in ch.vars
        var_ = SCIP.var(ch.o, vi)
        var_ == C_NULL && continue
        SCIP.@SCIP_CALL SCIP.SCIPaddVarLocksType(
            ch.o, var_, locktype, nlockspos + nlocksneg, nlockspos + nlocksneg)
    end
end

"""
    register_benders_handler!(m, g, separate!)

Attach the separation routine to `m` as a SCIP constraint handler. The branch
variables `v` are the ones a cut can constrain, so they are the ones to lock.
"""
function register_benders_handler!(m::Model, g::MetaGraph, separate!::Function)
    o = JuMP.unsafe_backend(m)::SCIP.Optimizer
    vars = [JuMP.index(m[:v][br...]) for br in edge_labels(g)]
    # SCIP ships its own Benders framework, whose handler is already named "benders";
    # ours needs a different name. And with `needs_constraints = false` the handler is
    # consulted on every solution, which SCIP's dual reductions do not account for --
    # the documented companion to the variable locks in `lock` below. (SCIP split the
    # old `misc/allowdualreds` into a strong and a weak flag; both must go.)
    set_optimizer_attribute(m, "misc/allowstrongdualreds", false)
    set_optimizer_attribute(m, "misc/allowweakdualreds", false)
    ch = BendersConshdlr(o, m, separate!, vars)
    SCIP.include_conshdlr(o, ch;
        name = "otsdbenders",
        description = "OTSD Benders cuts, separated at integer-feasible nodes",
        needs_constraints = false,
        enforce_priority = -15,
        check_priority = -7_000_000)
    ch
end
