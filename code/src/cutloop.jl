# ── Benders without callbacks ─────────────────────────────────────────────────
#
# Gurobi lets a cut be added the moment an integer-feasible topology appears, so the
# whole decomposition fits in one branch-and-cut tree. A solver without that hook
# needs the outer form instead: solve the master to optimality, separate against its
# solution, add whatever cuts came back, and solve again. The loop ends when a round
# separates nothing, at which point the incumbent satisfies every cut the oracle can
# produce — the same termination condition the callback version reaches inside the
# tree.
#
# What differs is cost, not correctness. Each round throws away the tree, so the
# master is re-solved from scratch as many times as there are rounds, and the tree
# never gets to prune with a cut it has not yet been told about. Expect it to be
# slower, and expect a different topology when the restricted problem has ties:
# the two searches visit different incumbents on the way.

"""
    benders_cut_loop!(m, separate!, backend, logfilename, t0, timeout) -> Int

Solve `m` by alternating master solves with separation, returning the number of
rounds used. `separate!` takes a [`CutSink`](@ref) and adds cuts for one topology;
it is the same routine the lazy callback drives.

Stops when a round produces no cut, when the master stops being solvable, when
`backend.max_rounds` is reached, or when `timeout` seconds have elapsed since `t0`.
"""
function benders_cut_loop!(m::Model, separate!::Function, backend::Backend,
                           logfilename::String, t0::DateTime, timeout::Float64)
    for r in 1:backend.max_rounds
        if timeout > 0
            spent = (now() - t0).value / 1000
            remaining = timeout - spent
            if remaining ≤ 1.0
                write_in_logfile(logfilename, "⏱ cut loop out of time after $(r - 1) rounds")
                return r - 1
            end
            set_timeout!(backend, m, remaining)
        end

        optimize!(m)

        # No incumbent to separate against: infeasible, unbounded, or out of time.
        if primal_status(m) != MOI.FEASIBLE_POINT
            write_in_logfile(logfilename,
                "cut loop stopped in round $r: $(termination_status(m))")
            return r - 1
        end

        sink = DirectSink(m)
        separate!(sink)
        added = cuts_added(sink)
        write_in_logfile(logfilename,
            "round $r: obj=$(round_obj(m))  cuts=$added")

        added == 0 && return r
    end

    write_in_logfile(logfilename, "cut loop hit max_rounds = $(backend.max_rounds)")
    backend.max_rounds
end

round_obj(m::Model) = try
    round(objective_value(m); digits = 3)
catch
    NaN
end
