struct OBendersCut #Optimality cuts
    edges::Vector{ELabel}
    contingency::ELabel
    v0::Vector{Bool}
    rc::Vector{Float64}
    var::Symbol
    objval::Float64
end

struct FBendersCut # Feasibility cuts
    edges::Vector{ELabel}
    contingency::ELabel
    v0::Vector{Bool}
    rc::Vector{Float64}
    FBendersCut(edges, contingency, v0, rc) = new(edges, contingency, v0, rc) # in case, objval is known = 1
    FBendersCut(edges, contingency, v0, rc, objval) = new(edges, contingency, v0, rc ./ objval) # normalize Farkas ray
end

struct GenericCut
    edges::Vector{ELabel}
    v0::Vector{Bool}
    rc::Vector{Float64}
    val::Float64
end

const BendersCut = Union{OBendersCut,FBendersCut,GenericCut}

"""Cut record for persistence and later reuse across phases.

`required_openings` lists branches that must be switchable (not frozen) for the cut
to be valid. Empty for regular Benders cuts; non-empty for bypass (GenericCut) cuts,
where it holds the inner-pocket branches the cut prescribes to open.
"""
struct CutRecord
    cut::BendersCut
    iteration::Int
    hamming::Int
    openbranches::Set{ELabel}
    required_openings::Set{ELabel}
end
CutRecord(cut, iter, ham, open) = CutRecord(cut, iter, ham, open, Set{ELabel}())

# ── Delivering a cut ──────────────────────────────────────────────────────────
#
# `cut_constraint` is the single definition of what each cut family asserts; the
# sink (see `cutsink.jl`) decides whether it is submitted to a callback or added to
# the model. The two delivery paths therefore cannot drift apart.

# ── What each cut family asserts ──────────────────────────────────────────────

cut_constraint(m::Model, bc::OBendersCut) = @build_constraint(
    bc.objval + sum(bc.rc[i] * (m[:v][br...] - bc.v0[i]) for (i, br) in enumerate(bc.edges)) ≤ m[bc.var][bc.contingency...])

cut_constraint(m::Model, bc::FBendersCut) = @build_constraint(
    1 + sum(bc.rc[i] * (m[:v][br...] - bc.v0[i]) for (i, br) in enumerate(bc.edges)) ≤ 0)

cut_constraint(m::Model, bc::GenericCut) = @build_constraint(
    sum(bc.rc[i] * (m[:v][br...]) for (i, br) in enumerate(bc.edges)) ≤ bc.val)

# ── Is the cut violated here? ─────────────────────────────────────────────────
#
# A Benders cut is valid everywhere, but only a *violated* one makes progress. Inside
# a branch-and-cut tree that distinction is the solver's problem: a lazy constraint
# that already holds is simply absorbed. A cut loop has to make it explicitly —
# otherwise every round adds the same inert cuts, the master never moves, and the
# loop cannot tell convergence from deadlock. The margin by which each family is
# violated is written once, here, and used by both the loop and its stopping test.

"How much `bc` is violated at the solution `s` reads. Non-positive means it holds."
function cut_violation(s::CutSink, bc::OBendersCut)
    lhs = bc.objval + sum(bc.rc[i] * (solution_value(s, s.m[:v][br...]) - bc.v0[i])
                          for (i, br) in enumerate(bc.edges))
    lhs - solution_value(s, s.m[bc.var][bc.contingency...])
end

function cut_violation(s::CutSink, bc::FBendersCut)
    1 + sum(bc.rc[i] * (solution_value(s, s.m[:v][br...]) - bc.v0[i])
            for (i, br) in enumerate(bc.edges))
end

cut_violation(s::CutSink, bc::GenericCut) =
    sum(bc.rc[i] * solution_value(s, s.m[:v][br...]) for (i, br) in enumerate(bc.edges)) - bc.val

"Below this, a cut is treated as already satisfied."
const CUT_VIOLATION_TOL = 1e-6

# ── Delivery ──────────────────────────────────────────────────────────────────

# In the tree, hand every cut to the solver: one that already holds costs nothing.
apply_cut!(s::LazySink, bc) = MOI.submit(s.m, MOI.LazyConstraint(s.cb_data), cut_constraint(s.m, bc))

# Everywhere else, take only what the current point violates. `dry` sinks count
# without adding, which is how a solution is judged acceptable without changing the
# model underneath the solver.
function apply_cut!(s::AccumulatingSink, bc)
    cut_violation(s, bc) > CUT_VIOLATION_TOL || return nothing
    is_dry(s) || add_cut_to_model!(s, cut_constraint(s.m, bc))
    s.added[] += 1
    nothing
end

add_cut_to_model!(s::AccumulatingSink, con) = add_constraint(s.m, con)


function create_pklostload_optimality_cut(edg::Vector{ELabel}, br::ELabel, pk::Pocket)
    pk.d ≤ 0 && return nothing
    OBendersCut(edg, br,
        [br ∉ pk.branches for br in edg],
        [br in pk.branches ? -pk.d : 0 for br in edg],
        :lostload, pk.d)
end

"""
    extract_cuts_for_phase(records) -> Vector{BendersCut}

Select cuts to inherit into the next Benders phase.
Regular cuts (OBendersCut / FBendersCut) are always included.
Bypass cuts (GenericCut) are excluded: they are conditional on the pocket
forming and would cause infeasibility as hard constraints. They are stored in
records solely for SBS extension via `required_branches_from_cutpool`.
"""
function extract_cuts_for_phase(records::Vector{CutRecord})::Vector{BendersCut}
    BendersCut[rec.cut for rec in records if !isa(rec.cut, GenericCut)]
end

"""Add an inherited cut as a hard constraint (not lazy)."""
function add_inherited_cut!(m::Model, bc::OBendersCut)
    # Guard: skip :s_flows cuts if variable doesn't exist (thermal mode without slack)
    if !haskey(object_dictionary(m), bc.var)
        return
    end
    con = @build_constraint(
        bc.objval + sum(bc.rc[i] * (m[:v][br...] - bc.v0[i]) for (i, br) in enumerate(bc.edges)) ≤ m[bc.var][bc.contingency...])
    add_constraint(m, con)
end

function add_inherited_cut!(m::Model, bc::FBendersCut)
    con = @build_constraint(
        1 + sum(bc.rc[i] * (m[:v][br...] - bc.v0[i]) for (i, br) in enumerate(bc.edges)) ≤ 0)
    add_constraint(m, con)
end

function add_inherited_cut!(m::Model, bc::GenericCut)
    con = @build_constraint(
        sum(bc.rc[i] * (m[:v][br...]) for (i, br) in enumerate(bc.edges)) ≤ bc.val)
    add_constraint(m, con)
end

function add_inherited_cuts!(m::Model, cuts::Vector{BendersCut})
    for cut in cuts
        add_inherited_cut!(m, cut)
    end
end