"""
In the main problem:
- only the base case
- the loss function is reducing the lostload
In the subproblems:
- in the first iteration, call to all the contingencies, then
- add one contingency at a time (the most constraining) to the contingency list when the actual contingency list no longer raise any feasibility cut.
"""
function mastercutpool(ec::ElementaryCase, logfilename::String;
    backend::Backend=default_backend(),
    contingencies_in_master::Union{Nothing,Vector{ELabel}}=nothing,
    ov_lostload_coef::Union{Float64,Dict{ELabel,Float64}}=1.,
    contingencies=nothing,
    warmstart_openings::Set{ELabel}=Set{ELabel}(),
    max_hamming::Int=0,
    sbs::Set{ELabel}=Set{ELabel}(),
    bigM_bound_multiplier::Real=2.0,
    initial_cuts::AbstractVector{BendersCut}=BendersCut[],
    seed::Int=0,
    ablate_screening::Bool=false,
    cf_cuts::Symbol=:closedform,
    timeout::Float64=0.0,
)
    g = ec.g; bus_orig = ec.bus_orig
    cf_cuts in (:closedform, :lp_pocket, :lp_free) ||
        error("cf_cuts must be :closedform, :lp_pocket or :lp_free — got :$cf_cuts")

    function _build_model(g::MetaGraph, basecase, master_totalloss_coef, contingencies_in_master, logfilename, bigM_bound_multiplier)

        n_1cases = isnothing(contingencies_in_master) ? [] : branch_to_case.(contingencies_in_master)
        cases = isnothing(contingencies_in_master) ? [BASECASEID] : [BASECASEID; n_1cases]
        m = init_model(backend, logfilename)

        @variable(m, v[busfrom in labels(g), outneighbor_labels(g, busfrom)], Bin)
        @expression(m, w[c in cases, busfrom in labels(g), busto in outneighbor_labels(g, busfrom)], case_to_branch(c) == (busfrom, busto) ? 0 : v[busfrom, busto])

        @constraint(m, [br in edge_labels(g); br ∉ sbs], m[:v][br...] == 1)
        @variable(m, -nv(g) ≤ c_flows[busfrom in labels(g), outneighbor_labels(g, busfrom)] ≤ nv(g))

        @variable(m, flows[cases, busfrom in labels(g), outneighbor_labels(g, busfrom)])
        @variable(m, ϕ[cases, bus in labels(g), 1:1])
        @variable(m, θf[cases, busfrom in labels(g), outneighbor_labels(g, busfrom)])
        @variable(m, θt[cases, busfrom in labels(g), outneighbor_labels(g, busfrom)])
        align_feeder_to_bus_angles!(m, g, cases)

        @variable(m, load[cases, labels(g)] ≥ 0)
        @variable(m, gen[cases, labels(g)] ≥ 0)

        if !isempty(warmstart_openings)
            warmstart_openings!(m, g, warmstart_openings)
            max_hamming ≠ 0 &&
                @constraint(m, allowed_modif, sum(m[:v][bus] for bus in warmstart_openings) +
                                              sum(1 - m[:v][bus] for bus in edge_labels(g) if bus ∉ warmstart_openings) ≤
                                              max_hamming)
        end

        phase_reference!(m, bus_orig, cases)
        bus_KCL!(m, g, bus_orig, cases)
        balance_basecase!(m, g, basecase)
        bigM = compute_tight_bigM(g; bigM_bound_multiplier)
        ohm!(m, g, cases, bigM)
        base_connectivity!(m, g, [bus_orig])

        if !isnothing(contingencies_in_master)
            flowslimits!(m, g, cases)
        else
            flowslimits!(m, g, [BASECASEID])
        end

        if !isnothing(contingencies_in_master)
            create_energization_state_variables!(m, g, n_1cases, false)
            align_feeder_to_bus_energization_state!(m, g, n_1cases)

            @variable(m, σ[n_1cases] ≥ 0)
            balance_N_1cases!(m, g, n_1cases)
            align_N_1_energization!(m, g, [bus_orig], n_1cases)
            # if false
            #     @variable(m, -nv(g) ≤ n_1c_flows[n_1cases, busfrom in labels(g), outneighbor_labels(g, busfrom)] ≤ nv(g))
            #     _fictitious_flows_N_1(m, g, _fixed_buses, n_1cases)
            # end
        end

        @variable(m, lostload[busfrom in labels(g), outneighbor_labels(g, busfrom)] ≥ 0)
        @expression(m, totallostload, sum(m[:lostload]))
        if !isnothing(contingencies_in_master)
            @constraint(m, [c in n_1cases], lostload[case_to_branch(c)...] == sum(max(g[bus], 0) - load[c, bus] for bus in labels(g)))
        end

        # A Dict form weights each contingency's shed load by its own p_c. Only the
        # objective changes:
        # `lostload[br]` is per contingency and the optimality cuts bound it individually
        # (`create_pklostload_optimality_cut`), so every cut stays valid under any weights.
        if isa(master_totalloss_coef, Dict)
            loss_obj = sum(master_totalloss_coef[br...] * m[:lostload][br...] for br in edge_labels(g))
        else
            loss_obj = master_totalloss_coef * m[:totallostload]
        end

        @objective(m, Min, loss_obj)
        m
    end

    # Separation. Given an integer-feasible topology, decide which contingencies need
    # a cut and hand each one to `sink`. Nothing here knows how the topology was
    # reached or how the cut is delivered, so both backends share it verbatim.
    function _separate!(sink::CutSink)
        iteration += 1
        v_0, openbranches = get_v_0_openbranches(g, sink, m)
        t1 = now()
        write_in_logfile(lfn, "Iteration $iteration \t$(canonicalize(t1-t0))\topen branches: $openbranches")


        bridge_to_pocket = create_bridge_to_pocket(ec, openbranches)

        for (br, pk) in bridge_to_pocket
            max_ll[br...] = max(max_ll[br...], pk.d)
        end


        sa_res = secured_dcpf(ec, Set(openbranches), contingencies, bridge_to_pocket=bridge_to_pocket, include_base_case=false)

        v_ctg = violating_contingencies(sa_res)

        vbrs = reduce(union, (vbrs for (ctg, vbrs) in violated_branches(sa_res)), init=Set())
        risk = sum(pk.d for pk in values(bridge_to_pocket))
        ov = risk + 0
            #  sum(ov_n1violation_coef[br] *
            #      (maximum(abs(flow(sa_res, ctg, br)) for ctg in v_ctg) -
            #       g[br...].p_max)
            #      for br in vbrs; init=0.)
        write_in_logfile(lfn, "v_ctg $v_ctg\n")
        if isempty(v_ctg)
            write_in_logfile(logfilename, "\tR=$ov $openbranches,")
            solutions[openbranches] = ov
            risk < best_risk && (best_risk = risk; best_sol = openbranches)
        end




        # Security screening filters which contingencies get a subproblem. Ablating it
        # probes every contingency; `v_ctg` stays the security verdict and still gates
        # the loss-of-load cuts below, so only the filter is removed.
        probed_ctg = ablate_screening ? contingencies : v_ctg
        for contingency in probed_ctg
            contingency in openbranches && continue
            res_subpb = contingency_subproblem(ec, openbranches, contingency, bridge_to_pocket, 0e-5; bigM_bound_multiplier=bigM_bound_multiplier, backend=lp_backend(backend))
            lp_feas_solves += 1
            if res_subpb.is_feasible
                # Expected for every non-violating contingency once screening is ablated.
                if !ablate_screening
                    @warn "Feasible should not happen"
                    @info "openbranches: $openbranches\ncontingency: $contingency"
                end
            else
                fcut = FBendersCut(edg, contingency, v_0, [res_subpb.reduced_cost[br...] for br in edg], res_subpb.dual_obj)
                apply_cut!(sink, fcut)
                push!(collected_cuts, CutRecord(fcut, iteration, hamming(openbranches, warmstart_openings), openbranches))

                appliedcuts += 1
            end
        end

        write_in_logfile(lfn, "applied cuts: $appliedcuts\n")

        # !isempty(v_ctg) && !is_reduceviolations && return

        ll_cuts = if cf_cuts === :closedform
            lostload_optimality_cuts(m, sink, edg, v_ctg, bridge_to_pocket)
        else
            cuts, n_lp = lp_lostload_optimality_cuts(m, sink, ec, edg, v_0, openbranches,
                v_ctg, bridge_to_pocket;
                free_π=(cf_cuts === :lp_free), bigM_bound_multiplier)
            lp_opt_solves += n_lp
            cuts
        end
        for ocut in ll_cuts
            push!(collected_cuts, CutRecord(ocut, iteration, hamming(openbranches, warmstart_openings), openbranches))
        end

    end

    # Gurobi calls this on every integer-feasible incumbent.
    function _callback(cb_data, cb_where::Cint)
        cb_where != GRB_CB_MIPSOL && return
        init_cb(cb_data, cb_where)
        _separate!(LazySink(m, cb_data))
    end

    write_in_logfile(logfilename, "contingencies: $contingencies")
    write_in_logfile(logfilename, "warmstart_openings: $warmstart_openings")
    write_in_logfile(logfilename, "max_hamming: $max_hamming")
    write_in_logfile(logfilename, "sbs: $(length(sbs)) branches")
    write_in_logfile(logfilename, "bigM_bound_multiplier: $bigM_bound_multiplier")
    write_in_logfile(logfilename, "seed: $seed")
    write_in_logfile(logfilename, "ablate_screening: $ablate_screening")
    write_in_logfile(logfilename, "cf_cuts: $cf_cuts")
    write_in_logfile(logfilename, "timeout: $timeout")

    lfn = "$(logfilename)_cb"

    edge_to_id = Dict(br => i for (i, br) in enumerate(edge_labels(g)))

    iteration = 0
    appliedcuts = 0
    lp_feas_solves = 0
    lp_opt_solves = 0
    bypasspockets = Pocket[]
    solutions = Dict{Set{ELabel}, Float64}()
    best_risk = Inf
    best_sol = nothing
    isnothing(contingencies) && (contingencies = collect(edge_labels(g)))
    m = _build_model(g, BASECASEID, ov_lostload_coef, contingencies_in_master, logfilename, bigM_bound_multiplier)

    add_inherited_cuts!(m, initial_cuts)
    if !isempty(initial_cuts)
        write_in_logfile(logfilename, "injected $(length(initial_cuts)) cuts from previous phase")
    end

    # Collect cuts during this solve
    collected_cuts = CutRecord[]

    ρ = Dict(br => 1. for br in edge_labels(g))
    max_ll = Dict(br => 0. for br in edge_labels(g))

    edg = collect(edge_labels(g))
    πstat = Dict(bus => 0 for bus in labels(g))

    # Both master backends inject cuts at integer-feasible nodes inside one tree, so
    # the search is one branch-and-cut tree either way. They differ only in the
    # interface the solver offers for it: a lazy-constraint callback on Gurobi, a
    # constraint handler on SCIP.
    supports_lazy(backend) || throw(ArgumentError(
        "$(backend_name(backend)) cannot take a cut inside the branch-and-cut " *
        "tree; the master needs GurobiBackend or SCIPBackend."))
    if backend isa SCIPBackend
        register_benders_handler!(m, g, _separate!)
    else
        enable_lazy!(m)
        MOI.set(m, Gurobi.CallbackFunction(), _callback)
    end

    t0 = now()

    # set_optimizer_attribute(m, "Cutoff", 0)
    set_seed!(backend, m, seed)
    timeout > 0.0 && set_timeout!(backend, m, timeout)

    optimize!(m)

    has_sol = is_solved_and_feasible(m) || primal_status(m) == MOI.FEASIBLE_POINT
    openings = has_sol ? getopenings(m) : warmstart_openings
    write_in_logfile(logfilename, "#iterations: $iteration\n Hamming distance: $(hamming(openings, warmstart_openings))")
    has_sol && write_in_logfile(logfilename, "#openings: $openings")

    write_in_logfile(logfilename, "collected $(length(collected_cuts)) cuts")

    write_in_logfile(logfilename, "LP subproblem solves: $(lp_feas_solves + lp_opt_solves) (feasibility: $lp_feas_solves, optimality: $lp_opt_solves)")

    (model=m, iterations=iteration, cuts=collected_cuts, best_sol=best_sol, best_risk=best_risk, solutions=solutions,
     lp_feas_solves=lp_feas_solves, lp_opt_solves=lp_opt_solves, lp_solves=lp_feas_solves + lp_opt_solves)
end

function lostload_optimality_cuts(m, sink::CutSink, edg, violating_ctg, bridge_to_pocket, πstat=nothing)
    applied = OBendersCut[]
    for (br, pk) in bridge_to_pocket
        br in violating_ctg && continue
        !isnothing(πstat) && foreach(bus -> πstat[bus] += 1, pk.buses)
        ocut = create_pklostload_optimality_cut(edg, br, pk)
        isnothing(ocut) && continue
        apply_cut!(sink, ocut)
        push!(applied, ocut)
    end
    applied
end

"""
    lp_lostload_optimality_cuts(m, sink, ec, edg, v_0, openbranches, violating_ctg,
                                bridge_to_pocket; free_π=false, ...)

LP counterpart of [`lostload_optimality_cuts`](@ref). Same candidate set — the bridges
of `bridge_to_pocket` that are not already handled by a feasibility cut — but the cut is
read off the duals of `contingency_subproblem` instead of the closed-form cut-set
argument. With `free_π=true` the subproblem drops `pocket_π_to_0!` and constrains `π`
through the generic connectivity block, so no property of the pocket enters the cut.

Returns `(cuts, n_lp)`: the cuts applied, and how many subproblem LPs were solved.
"""
function lp_lostload_optimality_cuts(m, sink::CutSink, ec, edg, v_0, openbranches, violating_ctg,
                                     bridge_to_pocket;
                                     free_π::Bool=false,
                                     bigM_bound_multiplier::Real=2.0,
                                     atol::Float64=1e-6)
    applied = OBendersCut[]
    n_lp = 0
    for (br, pk) in bridge_to_pocket
        br in violating_ctg && continue
        # Same guard as `create_pklostload_optimality_cut`: an empty pocket sheds nothing,
        # so the three cut modes are compared over an identical candidate set.
        pk.d ≤ 0 && continue
        res = contingency_subproblem(ec, openbranches, br, bridge_to_pocket;
            free_π, bigM_bound_multiplier)
        n_lp += 1
        res.is_feasible || continue
        res.obj ≤ atol && continue
        ocut = OBendersCut(edg, br, v_0, [res.reduced_cost[b...] for b in edg], :lostload, res.obj)
        apply_cut!(sink, ocut)
        push!(applied, ocut)
    end
    applied, n_lp
end
