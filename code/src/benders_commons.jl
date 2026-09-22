
function init_cb(cb_data, cb_where::Cint)
    resultP = Ref{Cint}()
    GRBcbget(cb_data, cb_where, GRB_CB_MIPNODE_STATUS, resultP)
    Gurobi.load_callback_variable_primal(cb_data, cb_where)
end

function get_v_0_openbranches(g, cb_data, m)
    v_0 = zeros(Bool, ne(g))
    openbranches = Set{ELabel}()
    for (i, br) in enumerate(edge_labels(g))
        v = round(callback_value(cb_data, m[:v][br...]))
        v_0[i] = isapprox(v, 1.0)
        v == 0.0 && push!(openbranches, br)
    end
    v_0, openbranches
end

function get_cb_value(cb_data, m::Model, fieldname::Symbol, indices)
    return Dict(id => callback_value(cb_data, m[fieldname][id...]) for id in indices)
end

function benders_subpb_res(m::Model)
    if is_solved_and_feasible(m; dual = true)
        return (
            is_feasible = true,
            obj = objective_value(m),
            reduced_cost = reduced_cost.(m[:v]),
            model = m,
        )
    else
        return (
            is_feasible = false,
            dual_obj = dual_objective_value(m),
            reduced_cost = reduced_cost.(m[:v]),
            model = m,
        )
    end
end

function identify_most_often_overloaded(
    rc,
    openbranches,
    contingencies;
    bridge_to_pocket::Union{Nothing,Dict{ELabel,Pocket}} = nothing,
)
    SA_res = secured_dcpf(
        rc.gc,
        Set(openbranches),
        contingencies;
        bridge_to_pocket,
    )
    br_to_overloads = Dict{ELabel,Vector{NamedTuple}}()
    for cbr in contingencies, br in edge_labels(g)
        flow = abs(flow(SA_res, cbr, br))
        limit = rc.gc.g[br...].p_max
        flow ≤ limit && continue
        br_to_overloads[br] = push!(
            get!(br_to_overloads, (contingency = br, ol = flow / limit), NamedTuple[]),
            cbr,
        )
    end
    isempty(br_to_overloads) && return nothing
    most_constrained = findmax(length, br_to_overloads)[2]
    branch_overloads = br_to_overloads[most_constrained]
    ol, big_contingency = findmax(bo -> bo.ol, branch_overloads)
    @info "br_to_overloads: $(br_to_overloads)\nmost_constrained: $most_constrained,\tbig_contingency: $big_contingency,\tol: $ol"
    return (
        branch = most_constrained,
        contingencies = br_to_overloads[most_constrained],
        big_contingency = big_contingency,
    )
end

