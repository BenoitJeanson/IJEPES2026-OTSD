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

struct SubpbBCPool
    contingency::ELabel
    OBC::Vector{OBendersCut}
    FBC::Vector{FBendersCut}
end
SubpbBCPool(contingency) = SubpbBCPool(contingency, OBendersCut[], FBendersCut[])

struct BCPool
    edges::Vector{ELabel}
    subpbbcp::Dict{ELabel,SubpbBCPool}
    appliedFcuts::Vector{Int}
end

BCPool(g::MetaGraph) = BCPool(collect(edge_labels(g)), Dict{ELabel,SubpbBCPool}(), [0])

add_cut!(bcpool::BCPool, bc::OBendersCut) =
    push!(get!(bcpool.subpbbcp, bc.contingency, SubpbBCPool(bc.contingency)).OBC, bc)

add_cut!(bcpool::BCPool, bc::FBendersCut) =
    push!(get!(bcpool.subpbbcp, bc.contingency, SubpbBCPool(bc.contingency)).FBC, bc)

add_cuts!(bcpool::BCPool, bcs::Vector) =
    foreach(bc -> add_cut!(bcpool, bc), bcs)

count_fcuts(bcpool::BCPool) =
    sum(length(subpb.FBC) for subpb in values(bcpool.subpbbcp))

function evaluate_benders_cut(bc::OBendersCut, v, subpbvalue::Float64)
    return bc.objval + sum(bc.rc .* (v .- bc.v0)) - subpbvalue
end
function evaluate_benders_cut(bc::FBendersCut, v)
    return 1 + sum(bc.rc .* (v .- bc.v0))
end

function evaluate_benders_cut(bc::OBendersCut, openbranches::Set{ELabel}, subpbvalue::Float64)
    v = [br ∉ openbranches for br in bc.edges]
    return evaluate_benders_cut(bc, v, subpbvalue)
end

function evaluate_benders_cut(bc::FBendersCut, openbranches::Set{ELabel})
    v = [br ∉ openbranches for br in bc.edges]
    return evaluate_benders_cut(bc, v)
end

function get_activated_ocuts(bcpool::BCPool, contingency::ELabel, v, subpbvalue::Float64)
    ocuts = get(bcpool.subpbbcp, contingency, SubpbBCPool(contingency, OBendersCut[], FBendersCut[])).OBC
    return Dict(ocut => val for ocut in ocuts
                for val = evaluate_benders_cut(ocut, v, subpbvalue)
                if val > 0)
end

function get_activated_fcuts(bcpool::BCPool, contingency::ELabel, v)
    fcuts = get(bcpool.subpbbcp, contingency, SubpbBCPool(contingency, OBendersCut[], FBendersCut[])).FBC
    return Dict(ocut => val for ocut in fcuts
                for val = evaluate_benders_cut(ocut, v)
                if val > 0)
end

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

# In the loop, add only what the current solution violates, and count it: that count
# is what tells the loop whether the round achieved anything.
function apply_cut!(s::DirectSink, bc)
    cut_violation(s, bc) > CUT_VIOLATION_TOL || return nothing
    add_constraint(s.m, cut_constraint(s.m, bc))
    s.added[] += 1
    nothing
end

apply_cuts!(s::CutSink, bcs) = foreach(bc -> apply_cut!(s, bc), bcs)

function apply_best_feasibility_cut!(m, bcpool::BCPool, sink::CutSink, v_0, iteration=0)
    feascutfound = false
    t1 = now()
    for contingency in keys(bcpool.subpbbcp)
        fcuts = get_activated_fcuts(bcpool, contingency, v_0)
        isempty(fcuts) && continue
        feascutfound = true
        _, cut = findmax(fcuts)
        bcpool.appliedFcuts[1] += 1
        apply_cut!(sink, cut)
    end
    t2 = now()
    if feascutfound
        return true
    end
    return false
end


function create_pklostload_optimality_cut(edg::Vector{ELabel}, br::ELabel, pk::Pocket)
    pk.d ≤ 0 && return nothing
    OBendersCut(edg, br,
        [br ∉ pk.branches for br in edg],
        [br in pk.branches ? -pk.d : 0 for br in edg],
        :lostload, pk.d)
end

function create_pklostload_optimality_cuts(bridge_to_pocket::Dict{ELabel,Pocket})
    cuts = OBendersCut[]
    for (br, pk) in bridge_to_pocket
        pk.d ≤ 0 && continue
        branches = union(pk.branches, [br])
        for contingency in branches
            pkbranches = setdiff(branches, [contingency])
            push!(cuts, create_pklostload_optimality_cut(g, contingency, Pocket(pk.buses, pkbranches, pk.innerbranches, pk.d)))
        end
    end
    cuts
end

function create_extra_pkloadloss_optimality_cuts(g::MetaGraph, br::ELabel, pk::Pocket)::Union{Vector{OBendersCut},Nothing}
    pk.d ≤ 0 && return OBendersCut[]
    branches = union(pk.branches, [br])
    cuts = OBendersCut[]
    for contingency in branches
        pkbranches = setdiff(branches, [contingency])
        push!(cuts, create_pklostload_optimality_cut(g, contingency, Pocket(g, pk.buses, pkbranches, pk.d)))
    end
    cuts
end

function create_extra_pkloadloss_optimality_cuts!(bcpool::BCPool, g::MetaGraph, br::ELabel, pk::Pocket)
    add_cuts!(bcpool, create_extra_pkloadloss_optimality_cuts(g, br, pk))
end

function create_extra_feasibility_cuts_from_zeros(edg, openbranches, contingency, reduced_cost, dual_val; atol=1e-6, nb_zeros_min=0)
    cuts = FBendersCut[]
    _null_coefs = Set(br for br in edg
                      if (isapprox(reduced_cost[br...] / dual_val, 0; atol=atol)) && (br ≠ contingency))
    length(_null_coefs) ≤ nb_zeros_min && return cuts
    _openbranches = Set([openbranches; contingency])

    open_wo_null = Set(br for br in _openbranches if !(br in _null_coefs))
    open_or_null = Set(br for br in edg if (br in _null_coefs || br in _openbranches))

    nb = 0
    for nonnull in open_wo_null
        nb += 1
        cut_openbranches = setdiff(open_wo_null, [nonnull])
        push!(cuts,
            FBendersCut(edg, contingency,
                [br ∉ cut_openbranches for br in edg],
                [br in cut_openbranches ? -1 :
                 br ∉ open_or_null ? 1 : 0
                 for br in edg]
            ))
    end
    cuts
end

function create_extra_feasibility_cuts_from_zeros!(bcpool::BCPool, openbranches, contingency, reduced_cost, dual_val; atol=1e-6, nb_zeros_min=0)
    f_cuts = create_extra_feasibility_cuts_from_zeros(bcpool.edges, openbranches, contingency, reduced_cost, dual_val; atol=atol, nb_zeros_min=nb_zeros_min)
    add_cuts!(bcpool, f_cuts)
end

function create_sflow_optimality_cuts(ec::ElementaryCase, openbranches, contingency::ELabel, vb::ELabel, overflow::Float64, bridge_to_pocket::Dict{ELabel,Pocket}, feasibility::Bool=false)
    g = ec.g
    overflow ≤ 0 && return nothing
    cuts = feasibility ? FBendersCut[] : OBendersCut[]

    if vb in keys(bridge_to_pocket)
        # if vb is a bridge, this means that contingency is also a bridge, and both pockets are not compatible together
        perimeter = collect(union(bridge_to_pocket[vb].branches, bridge_to_pocket[contingency].branches))
        nb_branches = length(perimeter)
        if feasibility
            push!(cuts, FBendersCut(collect(perimeter), vb, zeros(Int, nb_branches), -overflow * ones(Int, nb_branches), overflow))
        else
            push!(cuts, OBendersCut(collect(perimeter), vb, zeros(Int, nb_branches), -overflow * ones(Int, nb_branches), :s_flows, overflow))
        end
        return cuts
    end

    br_to_pk_with_ctg = create_bridge_to_pocket(ec, Set([openbranches; contingency]))
    if vb in keys(br_to_pk_with_ctg)
        # if vb becomes a bridge because of contingency, contingency is in the perimeter of the pocket, and any topology that result in that perimeter after contingency results in that same overflow, if the inside remains connected (or at least if the topology inside is the same).
        # in that case contingency is surely not a bridge
        pk = br_to_pk_with_ctg[vb]
        sb = collect(setdiff(union(pk.innerbranches, pk.branches, [vb]), [contingency]))
        @info "$openbranches\n$contingency -> $vb \n$sb"
        @info pk
        # sb = collect(edge_labels(g))
        v0 = [br ∉ openbranches for br in sb]
        sol = [("41", "56"), ("42", "56"), ("24", "25"), ("7", "29"), ("21", "22"), ("56", "57"), ("3", "4")]
        fcuts = FBendersCut(collect(sb), vb, v0, overflow * (2 * v0 .- 1), overflow)
        ebc = evaluate_benders_cut(fcuts, Set(sol))
        if ebc > 0
            @warn "cut: $fcuts\n$(ebc)"
        else
            @info "cut: $fcuts\n$(ebc)"
        end

        if feasibility
            push!(cuts, FBendersCut(collect(sb), vb, v0, overflow * (2 * v0 .- 1), overflow))
        else
            push!(cuts, OBendersCut(collect(sb), vb, v0, overflow * (2 * v0 .- 1), :s_flows, overflow))
        end
        return cuts
    end

    # (i) if vb ∉ pockets everything in pockets is unsensitive for vb, if vb ∈ pockets, everything outside the smallest of them is unsensitive => put no coef to branches insensitive
    # (ii) the sflow is ≥ overflow with fixing the non unsensitive branches to their v0.
    vb_pk = [pk for pk in values(bridge_to_pocket) if vb in pk.innerbranches]
    sensitive_branches = Set{ELabel}()
    if isempty(vb_pk)
        innerbranches = union(pk.innerbranches for pk in vb_pk)
        sensitive_branches = setdiff(edge_labels(g), innerbranches) # all branches not in any pocket are sensitive, perimeters remain sensitive
    else
        pk = argmin(pk -> length(pk.buses), vb_pk) # take the smallest pocket
        sensitive_branches = union(pk.branches, pk.innerbranches) # are sensitive, the branches that makes the pocket, and the inner branches
    end
    setdiff!(sensitive_branches, [contingency]) # per definition the contingency is not sensitive
    contingency in keys(bridge_to_pocket) && union!(sensitive_branches, bridge_to_pocket[contingency].branches)

    sb = collect(sensitive_branches)
    v0 = [br ∉ openbranches for br in sb]
    if feasibility
        push!(cuts, FBendersCut(collect(sb), vb, v0, overflow * (2 * v0 .- 1), overflow))
    else
        push!(cuts, OBendersCut(collect(sb), vb, v0, overflow * (2 * v0 .- 1), :s_flows, overflow))
    end
    return cuts
end

"""Save cut records to a .jls file."""
function save_cutpool(filename::String, records::Vector{CutRecord})
    open(filename, "w") do io
        Serialization.serialize(io, records)
    end
end

"""Load cut records from a .jls file."""
function load_cutpool(filename::String)::Vector{CutRecord}
    isfile(filename) || return CutRecord[]
    open(filename, "r") do io
        Serialization.deserialize(io)
    end
end

"""Save initial cuts (BendersCut vector) to a .jls file for phase reproducibility."""
function save_initial_cuts(filename::String, cuts::AbstractVector{BendersCut})
    open(filename, "w") do io
        Serialization.serialize(io, cuts)
    end
end

"""Load initial cuts from a .jls file."""
function load_initial_cuts(filename::String)::Vector{BendersCut}
    isfile(filename) || return BendersCut[]
    open(filename, "r") do io
        Serialization.deserialize(io)
    end
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

"""
    required_branches_from_cutpool(records) -> Set{ELabel}

Union of all `required_openings` across cut records. Used to extend the SBS so
that deferred bypass cuts become applicable in the next phase.
"""
function required_branches_from_cutpool(records::Vector{CutRecord})::Set{ELabel}
    result = Set{ELabel}()
    for rec in records
        union!(result, rec.required_openings)
    end
    result
end

"""Filter cuts to keep only those binding at the optimal solution v*."""
function filter_active_cuts(cuts::Vector{BendersCut}, openbranches::Set{ELabel}; atol::Float64=1e-4)::Vector{BendersCut}
    active = BendersCut[]
    for cut in cuts
        if isa(cut, FBendersCut)
            val = evaluate_benders_cut(cut, openbranches)
            # Keep if nearly binding or violated
            (val ≥ -atol) && push!(active, cut)
        elseif isa(cut, OBendersCut)
            # For OBendersCut without subproblem value, use heuristic: keep if objval term is reasonable
            # (conservative: keep most cuts, rely on LP solver to ignore non-binding ones)
            push!(active, cut)
        end
    end
    active
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