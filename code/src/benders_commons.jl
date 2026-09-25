
function init_cb(cb_data, cb_where::Cint)
    resultP = Ref{Cint}()
    GRBcbget(cb_data, cb_where, GRB_CB_MIPNODE_STATUS, resultP)
    Gurobi.load_callback_variable_primal(cb_data, cb_where)
end

function get_v_0_openbranches(g, sink::CutSink, m)
    v_0 = zeros(Bool, ne(g))
    openbranches = Set{ELabel}()
    for (i, br) in enumerate(edge_labels(g))
        v = round(solution_value(sink, m[:v][br...]))
        v_0[i] = isapprox(v, 1.0)
        v == 0.0 && push!(openbranches, br)
    end
    v_0, openbranches
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

