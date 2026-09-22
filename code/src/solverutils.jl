"""
    print_conflict!(model)

Compute and print a conflict for an infeasible `model`.
"""
function print_conflict!(model)
    JuMP.compute_conflict!(model)
    ctypes = list_of_constraint_types(model)
    for (F, S) in ctypes
        cons = all_constraints(model, F, S)
        for i in eachindex(cons)
            isassigned(cons, i) || continue
            con = cons[i]
            cst = MOI.get(model, MOI.ConstraintConflictStatus(), con)
            cst == MOI.IN_CONFLICT && @info con
        end
    end
    return nothing
end
    