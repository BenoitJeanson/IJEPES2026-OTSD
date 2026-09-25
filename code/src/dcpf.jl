struct SA_result
    gc::GridCase
    bus_orig::VLabel
    contingencies::Dict{ELabel,Int}
    branches::Dict{ELabel,Int}
    flows::Matrix{Float64}
    black_buses::Dict{ELabel,Set{VLabel}}
end

Base.show(io::IO, sr::SA_result) = print(io, "SA_result($(sr.gc), $(sr.bus_orig), $(length(sr.contingencies)) contingencies, $(length(sr.branches)) branches)")

function setflows!(g::MetaGraph, flows::Vector{Float64})
    for (i, e) in enumerate(edge_labels(g))
        g[e...].p = flows[i]
    end
end


function setflows!(ec::ElementaryCase, flows::Vector{Float64})
    setflows!(ec.g, flows)
end



"""
    DisconnectedTopology(outages)

The reduced susceptance matrix is singular: with `outages` open, some buses are no
longer reachable from the reference bus, so their angles are not determined.

A master that enforces `base_connectivity!` never produces such a topology, and the
Gurobi path never raises this. SCIP's constraint handler is consulted on candidate
solutions before they have cleared every other handler, so it can be asked about one —
and an exception crossing back into SCIP's C frame is undefined behaviour, which is
why `scip.jl` catches this rather than letting it escape.
"""
struct DisconnectedTopology <: Exception
    outages::Set{ELabel}
end

Base.showerror(io::IO, e::DisconnectedTopology) =
    print(io, "DisconnectedTopology: the network is not connected to the reference bus ",
              "with ", length(e.outages), " branches open")

"Solve for the bus angles, reporting a disconnected topology as such."
function _solve_angles(A, b, outages)
    try
        A \ b
    catch e
        e isa LinearAlgebra.SingularException || rethrow()
        throw(DisconnectedTopology(outages))
    end
end

function secured_dcpf(
    ec::ElementaryCase,
    outages::Set{ELabel} = Set{ELabel}(),
    contingencies::Union{Nothing,Vector{ELabel}} = nothing;
    bridge_to_pocket::Union{Nothing,Dict{ELabel,Pocket}} = nothing,
    include_base_case::Bool = true,
)
    _ec = ec
    contingencies0 = contingencies
    _outages = outages

    g = _ec.g
    bus_orig = _ec.bus_orig
    _bridge_to_pocket =
        isnothing(bridge_to_pocket) ? create_bridge_to_pocket(_ec, _outages) :
        bridge_to_pocket
    black_buses = Dict{ELabel,Set{VLabel}}(
        br=>_bridge_to_pocket[br].buses for br in keys(_bridge_to_pocket)
    )

    BASECONTINGENCY = ("", "")
    _contingencies = include_base_case ? ELabel[BASECONTINGENCY] : ELabel[]
    append!(
        _contingencies,
        isnothing(contingencies0) ? collect(edge_labels(g)) : contingencies0,
    )

    orig_id = code_for(g, bus_orig)
    edge_ids = Dict(e => i for (i, e) in enumerate(edge_labels(g)))

    A = incidence_matrix(g; oriented = true)
    for br in _outages
        A[:, edge_ids[br]] .= 0
    end

    D = spdiagm(map(e -> g[e...].b, edge_labels(g)))
    B = A * D * A'

    B_wo_orig = B[1:end .≠ orig_id, 1:end .≠ orig_id]
    p_wo_orig = [g[bus] for bus in labels(g) if bus ≠ bus_orig]

    ϕ_wo_orig_to_ϕ = spzeros(nv(g), nv(g) - 1)
    foreach(i -> ϕ_wo_orig_to_ϕ[i+(i≥orig_id), i] = 1, 1:(nv(g)-1))
    ϕ_to_flows = D * A' * ϕ_wo_orig_to_ϕ

    res_flows = zeros(length(_contingencies), ne(g))

    function _change_b_wo_orig(i, j, δb)
        (i == orig_id || j == orig_id) && return
        B_wo_orig[i-(i>orig_id), j-(j>orig_id)] += δb
    end

    for (i, br) in enumerate(_contingencies)
        br ≠ BASECONTINGENCY && (e_id = edge_ids[br])
        from_id, to_id = 0, 0

        if br ≠ BASECONTINGENCY
            from_id, to_id = code_for(g, from(br)), code_for(g, to(br))

            # remove the influence of the impedance of the open branch (don't do it if already open)
            b = g[br...].b
            if br ∉ _outages
                _change_b_wo_orig(from_id, from_id, -b)
                _change_b_wo_orig(to_id, to_id, -b)
                _change_b_wo_orig(from_id, to_id, +b)
                _change_b_wo_orig(to_id, from_id, +b)
            end
        end

        if br in keys(_bridge_to_pocket)  # A pocket will be deenergized
            pk = _bridge_to_pocket[br]
            pkbuses = Set(pk.buses)
            inbus_ids, pkbus_ids = Int[], Int[]
            for (j, bus) in enumerate(labels(g))
                j == orig_id && continue
                k = j - (j ≥ orig_id)
                if (bus ∉ pkbuses)
                    push!(inbus_ids, k)
                else
                    push!(pkbus_ids, k)
                end
            end

            δp = zeros(nv(g)) # to be added to p_wo_orig

            #balancing inside the area
            imbalance = -sum(p_wo_orig[pkbus_ids]) # the global system is balanced, so the imbalance is minus what is in the pocket
            ingen_ids = [i for i in inbus_ids if p_wo_orig[i] ≤ 0]
            genglob = sum(p_wo_orig[ingen_ids]) + (g[bus_orig] ≤ 0 ? g[bus_orig] : 0)  # if the bus orig is a generator, it was part of the initial balance
            foreach(
                bus_id -> δp[bus_id] = -p_wo_orig[bus_id] * imbalance / genglob,
                ingen_ids,
            )

            ϕ_wo_orig = zeros(nv(g) - 1)
            ϕ_wo_orig[inbus_ids] = _solve_angles(
                B_wo_orig[inbus_ids, inbus_ids], p_wo_orig[inbus_ids] + δp[inbus_ids], _outages)
            flows = ϕ_to_flows * ϕ_wo_orig
            br ≠ BASECONTINGENCY && (flows[e_id] = 0)
            res_flows[i, :] = flows

        else            # no imbalance to handle
            ϕ_wo_orig = _solve_angles(B_wo_orig, p_wo_orig, _outages)
            flows = ϕ_to_flows * ϕ_wo_orig
            br ≠ BASECONTINGENCY && (flows[e_id] = 0)
            res_flows[i, :] = flows
        end

        # restore the influence of the impedance of the open branch for next iterations
        if br ≠ BASECONTINGENCY && br ∉ _outages
            _change_b_wo_orig(from_id, from_id, +b)
            _change_b_wo_orig(to_id, to_id, +b)
            _change_b_wo_orig(from_id, to_id, -b)
            _change_b_wo_orig(to_id, from_id, -b)
        end
    end

    SA_result(
        _ec,
        bus_orig,
        Dict(e => i for (i, e) in enumerate(_contingencies)),
        edge_ids,
        res_flows,
        black_buses,
    )
end

function flow(sr::SA_result, ctg::ELabel, br::ELabel)
    sr.flows[sr.contingencies[ctg], sr.branches[br]]
end 

function violated_branches(sr::SA_result, contingency::ELabel)::Set{ELabel}
    branches = Set{ELabel}()
    g = sr.gc.g
    for br in edge_labels(g)
        if abs(flow(sr, contingency, br)) > g[br...].p_max
            push!(branches, br)
        end
    end
    branches
end

function violated_branches(sr::SA_result)::Dict{ELabel,Set{ELabel}}
    ctg_to_vbrs = Dict{ELabel,Set{ELabel}}()
    for ctg in keys(sr.contingencies)
        vbrs = violated_branches(sr, ctg)
        !isempty(vbrs) && (ctg_to_vbrs[ctg] = vbrs)
    end
    ctg_to_vbrs
end

function max_overload(sr::SA_result, contingency::ELabel, branches::Vector{ELabel})
    maxval, index =
        findmax(br -> abs(flow(sr, contingency, br)) / sr.gc.g[br...].p_max, branches)
    maxval, branches[index]
end

function max_overload(sr::SA_result)
    max_ol = 0.0
    all_branches = collect(edge_labels(sr.gc.g))
    for br in all_branches
        max_ol = maximum(max_ol, max_overload(sr, br, all_branches)[1])
    end
    max_ol
end

function violating_contingencies(sr::SA_result)
    Set(ctg for (ctg, i) in sr.contingencies if !isempty(violated_branches(sr, ctg)))
end

function identify_most_constraining_contingency(sr::SA_result)
    ctg_to_overloads = Dict{ELabel,Set{ELabel}}()
    for ctg in keys(sr.contingencies)
        vb = violated_branches(sr, ctg)
        !isempty(vb) && (ctg_to_overloads[ctg] = vb)
    end
    base_violation = violated_branches(sr, NULLEDGE)
    isempty(ctg_to_overloads) && isempty(base_violation) && return nothing
    sorted_contingencies =
        sort(collect(keys(ctg_to_overloads)), by = ctg -> length(ctg_to_overloads[ctg]))
    return (
        ctg_to_overloads = ctg_to_overloads,
        sorted_contingencies = sorted_contingencies,
        base_violation = base_violation,
    )
end

function identify_most_constraining_contingency(
    rc::RichCase,
    openbranches = Set{ELabel}(),
    contingencies::Union{Nothing,Vector{ELabel}} = nothing,
)
    g = rc.gc.g
    SA_res = secured_dcpf(rc.gc, openbranches, contingencies)
    return identify_most_constraining_contingency(SA_res)
end

function reduced_sa(rc::RichCase, outages::Set{ELabel} = Set{ELabel}())
    sa_res = secured_dcpf(rc.gc, outages)
    vc = violating_contingencies(sa_res)
    Dict(c => violated_branches(sa_res, c) for c in vc)
end

str(rsa::Dict{ELabel,Set{ELabel}}) =
    join(["$(str(branch)) → $(str(rsa[branch]))" for branch in sort(collect(keys(rsa)))], ",   ")