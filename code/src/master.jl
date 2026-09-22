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
    forced_branch_openings::Union{Nothing,Set{ELabel}}=nothing,
    ov_lostload_coef::Union{Float64,Dict{ELabel,Float64}}=1.,
    ov_openingpenalty_coef::Float64=0.,
    ov_n1violation_coef::Union{Float64,Dict{ELabel,Float64},Nothing}=nothing,
    n1violation_multiplier::Float64=0.,
    master_enforce_n1violation::Bool=false,
    analytic_viol::Bool=false,
    dual_viol::Bool=true,
    withcallbacks=true,
    contingencies=nothing,
    warmstart_openings::Set{ELabel}=Set{ELabel}(),
    include_modularity_cut::Bool=false,
    max_hamming::Int=0,
    nb_common_openings::Int=0,
    sbs::Set{ELabel}=Set{ELabel}(),
    forced_open_branches::Vector{ELabel}=ELabel[],
    forced_energized_buses::Vector{String}=String[],
    with_cutpool=true,
    cp_nb_zeros_min=0,
    is_master_π_binary::Bool=false,
    draw_callback::Bool=false,
    detailed_cb_logs::Bool=false,
    bigM_π::Float64=1e2,
    bigM_flows::Float64=1e3,
    tight_bigM::Bool=true,
    θ_max_bigM::Real=π,
    bigM_bound_multiplier::Real=2.0,
    nb_max_openings::Int=0,
    initial_cuts::AbstractVector{BendersCut}=BendersCut[],
    rc=nothing,
    best_known_obj::Float64=Inf,
    Δ_stop::Float64=Inf,
    kernel_branches::Vector{ELabel}=ELabel[],
    min_inside_kernel::Int=0,
    max_outside_kernel::Int=typemax(Int),
    hotkernel_branches::Vector{ELabel}=ELabel[],
    min_inside_hotkernel::Int=0,
    use_bypass_cuts::Bool=false,
    seed::Int=0,
    mip_focus::Int=0,
    save_cuts_to_disk::Bool=false,
    extra_warmstarts::Vector{Vector{ELabel}}=Vector{ELabel}[],
    ablate_screening::Bool=false,
    cf_cuts::Symbol=:closedform,
    timeout::Float64=0.0,
)
    g = ec.g; bus_orig = ec.bus_orig
    cf_cuts in (:closedform, :lp_pocket, :lp_free) ||
        error("cf_cuts must be :closedform, :lp_pocket or :lp_free — got :$cf_cuts")

    function _build_model(g::MetaGraph, basecase, bigM_π, bigM_flows, master_totalloss_coef, ov_openingpenalty_coef, ov_n1violation_coef, contingencies_in_master, nb_max_openings, nb_common_openings, logfilename, tight_bigM, θ_max_bigM, bigM_bound_multiplier)

        n_1cases = isnothing(contingencies_in_master) ? [] : branch_to_case.(contingencies_in_master)
        cases = isnothing(contingencies_in_master) ? [BASECASEID] : [BASECASEID; n_1cases]
        m = init_model(backend, logfilename)

        @variable(m, v[busfrom in labels(g), outneighbor_labels(g, busfrom)], Bin)
        @expression(m, w[c in cases, busfrom in labels(g), busto in outneighbor_labels(g, busfrom)], case_to_branch(c) == (busfrom, busto) ? 0 : v[busfrom, busto])

        @constraint(m, [br in edge_labels(g); br ∉ sbs], m[:v][br...] == 1)
        for br in forced_open_branches
            @constraint(m, m[:v][br...] == 0)
        end
        nb_max_openings > 0 && @constraint(m, sum(1 - m[:v][br...] for br in edge_labels(g)) ≤ nb_max_openings)

        @variable(m, -nv(g) ≤ c_flows[busfrom in labels(g), outneighbor_labels(g, busfrom)] ≤ nv(g))

        @variable(m, flows[cases, busfrom in labels(g), outneighbor_labels(g, busfrom)])
        @variable(m, ϕ[cases, bus in labels(g), 1:1])
        @variable(m, θf[cases, busfrom in labels(g), outneighbor_labels(g, busfrom)])
        @variable(m, θt[cases, busfrom in labels(g), outneighbor_labels(g, busfrom)])
        align_feeder_to_bus_angles!(m, g, cases)

        @variable(m, load[cases, labels(g)] ≥ 0)
        @variable(m, gen[cases, labels(g)] ≥ 0)

        if !isempty(warmstart_openings)
            if !isempty(extra_warmstarts)
                Base.invokelatest(multi_warmstart_openings!, m, g, warmstart_openings, extra_warmstarts)
            else
                warmstart_openings!(m, g, warmstart_openings)
            end
            include_modularity_cut && modularity_cut!(m, g, warmstart_openings)
            max_hamming ≠ 0 &&
                @constraint(m, allowed_modif, sum(m[:v][bus] for bus in warmstart_openings) +
                                              sum(1 - m[:v][bus] for bus in edge_labels(g) if bus ∉ warmstart_openings) ≤
                                              max_hamming)
            nb_common_openings ≠ 0 &&
                @constraint(m, common_openings, sum(1 - m[:v][bus] for bus in edge_labels(g) if bus ∈ warmstart_openings) ≥ nb_common_openings)
        end

        forced_branch_status!(m, g, forced_branch_openings)

        if !isempty(kernel_branches) && (min_inside_kernel > 0 || max_outside_kernel < typemax(Int))
            kernel_set = Set(kernel_branches)
            min_inside_kernel > 0 && @constraint(m, kernel_constraint,
                sum(1 - m[:v][br...] for br in kernel_branches if br in sbs) >= min_inside_kernel)
            max_outside_kernel < typemax(Int) && @constraint(m, outside_kernel_constraint,
                sum(1 - m[:v][br...] for br in edge_labels(g) if br ∉ kernel_set && br in sbs) <= max_outside_kernel)
        end
        if !isempty(hotkernel_branches) && min_inside_hotkernel > 0
            @constraint(m, hotkernel_constraint,
                sum(1 - m[:v][br...] for br in hotkernel_branches if br in sbs) >= min_inside_hotkernel)
        end

        phase_reference!(m, bus_orig, cases)
        bus_KCL!(m, g, bus_orig, cases)
        balance_basecase!(m, g, basecase)
        if tight_bigM
            bigM = compute_tight_bigM(g; θ_max=θ_max_bigM, bigM_bound_multiplier=bigM_bound_multiplier)
            ohm!(m, g, cases, bigM)
            base_connectivity!(m, g, [bus_orig])
        else
            ohm!(m, g, cases, bigM_flows)
            base_connectivity!(m, g, [bus_orig])
        end

        if !isnothing(contingencies_in_master) && !is_reduceviolations
            flowslimits!(m, g, cases)
        else
            flowslimits!(m, g, [BASECASEID])
        end

        if !isnothing(contingencies_in_master)
            create_energization_state_variables!(m, g, n_1cases, SubstationConfs(), is_master_π_binary)
            align_feeder_to_bus_energization_state!(m, g, n_1cases)

            @variable(m, σ[n_1cases] ≥ 0)
            balance_N_1cases!(m, g, n_1cases)
            align_N_1_energization!(m, g, [bus_orig], n_1cases)
            for bus in forced_energized_buses
                @constraint(m, [c in n_1cases], m[:π][c, bus, 1] == 1)
            end

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

        @expression(m, totalopening, sum(1 .- m[:v]))
        if is_reduceviolations
            @variable(m, s_flows[busfrom in labels(g), outneighbor_labels(g, busfrom)] ≥ 0)
            master_enforce_n1violation && @constraint(m, [busfrom in labels(g), busto in outneighbor_labels(g, busfrom)], m[:s_flows][busfrom, busto] == 0)
            @expression(m, totaln1violation, sum(m[:s_flows][:, :]))
            if !isnothing(contingencies_in_master)
                flowslimit_slack!(m, g)
            end
        else
            @expression(m, totaln1violation, 0)
        end
        if isnothing(ov_n1violation_coef)
            viol_n1 = 0
        elseif isa(ov_n1violation_coef, Dict)
            viol_n1 = sum(ov_n1violation_coef[br...] * m[:s_flows][br...] for br in edge_labels(g))
        else
            viol_n1 = ov_n1violation_coef * sum(m[:s_flows][:, :])
        end

        # A Dict form weights each contingency's shed load by its own p_c, the same
        # pattern `ov_n1violation_coef` already uses above. Only the objective changes:
        # `lostload[br]` is per contingency and the optimality cuts bound it individually
        # (`create_pklostload_optimality_cut`), so every cut stays valid under any weights.
        if isa(master_totalloss_coef, Dict)
            loss_obj = sum(master_totalloss_coef[br...] * m[:lostload][br...] for br in edge_labels(g))
        else
            loss_obj = master_totalloss_coef * m[:totallostload]
        end

        @objective(m, Min, loss_obj +
                           ov_openingpenalty_coef * m[:totalopening] +
                           viol_n1)
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

        with_cutpool && apply_best_feasibility_cut!(m, bcpool, sink, v_0, iteration) && return

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

        if isempty(v_ctg) && ov < best_known_obj - Δ_stop
        # if ov < best_known_obj - Δ_stop
            write_in_logfile(logfilename, "⚡ Early stop: ov=$(round(ov;digits=3)) < best_known=$(best_known_obj) - Δ=$(Δ_stop)")
            early_stopped[] = true
            GRBterminate(JuMP.backend(m))
            return
        end

        if draw_callback
            vb_dict = violated_branches(sa_res)
            vb = isempty(vb_dict) ? [] : union(values(vb_dict)...)
            display(draw(rc;
                edge_color=[br in openbranches ? :black : br in vb ? :red : :green for br in edge_labels(rc.gc.g)],
                edge_widths=[br in openbranches ? 6 : br in vb ? 6 : 1 for br in edge_labels(rc.gc.g)]))
        end


        # Security screening filters which contingencies get a subproblem. Ablating it
        # probes every contingency; `v_ctg` stays the security verdict and still gates
        # the loss-of-load cuts below, so only the filter is removed.
        probed_ctg = ablate_screening ? contingencies : v_ctg
        for contingency in probed_ctg
            contingency in openbranches && continue
            if is_reduceviolations
                for vb in violated_branches(sa_res, contingency)
                    sf = abs(flow(sa_res, contingency, vb)) - g[vb...].p_max
                    if dual_viol
                        lp_feas_solves += 1
                        res_subpb = contingency_subproblem(ec, openbranches, contingency, bridge_to_pocket, bigM_π, bigM_flows; reduce_violations=true, monitored_branch=vb, tight_bigM=tight_bigM, θ_max_bigM=θ_max_bigM, bigM_bound_multiplier=bigM_bound_multiplier, backend=backend)
                        if res_subpb.is_feasible
                            ocut = OBendersCut(edg, vb, v_0, [res_subpb.reduced_cost[br...] for br in edg], :s_flows, res_subpb.obj)
                            apply_cut!(sink, ocut)
                            push!(collected_cuts, CutRecord(ocut, iteration, hamming(openbranches, warmstart_openings), openbranches))
                        else
                            # rare: f > α·p_max — disconnected island; fall back to feasibility cut
                            lp_feas_solves += 1
                            res_feas = contingency_subproblem(ec, openbranches, contingency, bridge_to_pocket, bigM_π, bigM_flows; tight_bigM=tight_bigM, θ_max_bigM=θ_max_bigM, bigM_bound_multiplier=bigM_bound_multiplier, backend=backend)
                            if !res_feas.is_feasible
                                fcut = FBendersCut(edg, contingency, v_0, [res_feas.reduced_cost[br...] for br in edg], res_feas.dual_obj)
                                apply_cut!(sink, fcut)
                                push!(collected_cuts, CutRecord(fcut, iteration, hamming(openbranches, warmstart_openings), openbranches))
                            end
                            break
                        end
                    end
                    if analytic_viol
                        ocuts = create_sflow_optimality_cuts(ec, openbranches, contingency, vb, sf, bridge_to_pocket)
                        apply_cuts!(sink, ocuts)
                        write_in_logfile(lfn, "analytic cuts: $(length(ocuts)) for $vb")
                    end
                end
            else
                res_subpb = contingency_subproblem(ec, openbranches, contingency, bridge_to_pocket, bigM_π, bigM_flows, 0e-5; tight_bigM=tight_bigM, θ_max_bigM=θ_max_bigM, bigM_bound_multiplier=bigM_bound_multiplier, backend=backend)
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

                    with_cutpool && create_extra_feasibility_cuts_from_zeros!(bcpool, openbranches, contingency, res_subpb.reduced_cost, res_subpb.dual_obj;
                        nb_zeros_min=cp_nb_zeros_min)
                    appliedcuts += 1
                end
            end
        end

        write_in_logfile(lfn, "applied cuts: $appliedcuts\n")

        # !isempty(v_ctg) && !is_reduceviolations && return

        ll_cuts = if cf_cuts === :closedform
            lostload_optimality_cuts(m, sink, edg, v_ctg, bridge_to_pocket)
        else
            cuts, n_lp = lp_lostload_optimality_cuts(m, sink, ec, edg, v_0, openbranches,
                v_ctg, bridge_to_pocket, bigM_π, bigM_flows;
                free_π=(cf_cuts === :lp_free), tight_bigM, θ_max_bigM, bigM_bound_multiplier)
            lp_opt_solves += n_lp
            cuts
        end
        for ocut in ll_cuts
            push!(collected_cuts, CutRecord(ocut, iteration, hamming(openbranches, warmstart_openings), openbranches))
        end

        if use_bypass_cuts
            bp_cuts = bypass_cut_from_pockets(m, sink, g, openbranches, sbs, bridge_to_pocket, bypasspockets)
            for (cut, req) in bp_cuts
                push!(collected_cuts, CutRecord(cut, iteration, hamming(openbranches, warmstart_openings), openbranches, req))
            end
        end
    end

    # Gurobi calls this on every integer-feasible incumbent.
    function _callback(cb_data, cb_where::Cint)
        cb_where != GRB_CB_MIPSOL && return
        init_cb(cb_data, cb_where)
        _separate!(LazySink(m, cb_data))
    end

    write_in_logfile(logfilename, "withcallbacks: $withcallbacks")
    write_in_logfile(logfilename, "contingencies: $contingencies")
    write_in_logfile(logfilename, "warmstart_openings: $warmstart_openings")
    write_in_logfile(logfilename, "extra_warmstarts: $(length(extra_warmstarts))")
    write_in_logfile(logfilename, "include_modularity_cut: $include_modularity_cut")
    write_in_logfile(logfilename, "max_hamming: $max_hamming")
    write_in_logfile(logfilename, "sbs: $(length(sbs)) branches")
    write_in_logfile(logfilename, "forced_open_branches: $forced_open_branches")
    write_in_logfile(logfilename, "kernel_branches: $kernel_branches")
    write_in_logfile(logfilename, "min_inside_kernel: $min_inside_kernel")
    write_in_logfile(logfilename, "max_outside_kernel: $max_outside_kernel")
    write_in_logfile(logfilename, "hotkernel_branches: $hotkernel_branches")
    write_in_logfile(logfilename, "min_inside_hotkernel: $min_inside_hotkernel")
    write_in_logfile(logfilename, "with_cutpool: $with_cutpool")
    write_in_logfile(logfilename, "cp_nb_zeros_min: $cp_nb_zeros_min")
    write_in_logfile(logfilename, "is_master_π_binary: $is_master_π_binary")
    write_in_logfile(logfilename, "draw_callback: $draw_callback")
    write_in_logfile(logfilename, "bigM_π: $bigM_π")
    write_in_logfile(logfilename, "bigM_flows: $bigM_flows")
    write_in_logfile(logfilename, "tight_bigM: $tight_bigM")
    write_in_logfile(logfilename, "θ_max_bigM: $θ_max_bigM")
    write_in_logfile(logfilename, "bigM_bound_multiplier: $bigM_bound_multiplier")
    write_in_logfile(logfilename, "nb_max_openings: $nb_max_openings")
    write_in_logfile(logfilename, "nb_common_openings: $nb_common_openings")
    write_in_logfile(logfilename, "seed: $seed")
    write_in_logfile(logfilename, "mip_focus: $mip_focus")
    write_in_logfile(logfilename, "ablate_screening: $ablate_screening")
    write_in_logfile(logfilename, "cf_cuts: $cf_cuts")
    write_in_logfile(logfilename, "timeout: $timeout")

    lfn = detailed_cb_logs ? "$(logfilename)_cb" : ""

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
    is_reduceviolations = !isnothing(ov_n1violation_coef) || master_enforce_n1violation
    ablate_screening && is_reduceviolations &&
        error("ablate_screening is only defined for the feasibility-cut path (is_reduceviolations=false)")
    m = _build_model(g, BASECASEID, bigM_π, bigM_flows, ov_lostload_coef, ov_openingpenalty_coef, ov_n1violation_coef, contingencies_in_master, nb_max_openings, nb_common_openings, logfilename, tight_bigM, θ_max_bigM, bigM_bound_multiplier)

    cutdir = joinpath(dirname(Base.active_project()), "tmp")
    if save_cuts_to_disk && !isempty(initial_cuts)
        initfile = joinpath(cutdir, "$(logfilename)_initial_cuts.jls")
        save_initial_cuts(initfile, initial_cuts)
        write_in_logfile(logfilename, "saved $(length(initial_cuts)) initial cuts to $(basename(initfile))")
    end
    add_inherited_cuts!(m, initial_cuts)
    if !isempty(initial_cuts)
        write_in_logfile(logfilename, "injected $(length(initial_cuts)) cuts from previous phase")
    end

    # Collect cuts during this solve
    collected_cuts = CutRecord[]
    early_stopped = Ref(false)

    ρ = Dict(br => 1. for br in edge_labels(g))
    max_ll = Dict(br => 0. for br in edge_labels(g))

    bcpool = BCPool(g)
    edg = collect(edge_labels(g))
    πstat = Dict(bus => 0 for bus in labels(g))

    # A backend without lazy constraints separates between master solves instead;
    # `benders_cut_loop!` below drives the same `_separate!`.
    if withcallbacks && supports_lazy(backend)
        enable_lazy!(m)
        MOI.set(m, Gurobi.CallbackFunction(), _callback)
    end

    t0 = now()

    # set_optimizer_attribute(m, "Cutoff", 0)
    set_seed!(backend, m, seed)
    mip_focus > 0 && set_mip_focus!(backend, m, mip_focus)
    timeout > 0.0 && set_timeout!(backend, m, timeout)

    if withcallbacks && !supports_lazy(backend)
        benders_cut_loop!(m, _separate!, backend, logfilename, t0, timeout)
    else
        optimize!(m)
    end

    has_sol = is_solved_and_feasible(m) || primal_status(m) == MOI.FEASIBLE_POINT
    openings = has_sol ? getopenings(m) : warmstart_openings
    write_in_logfile(logfilename, "#iterations: $iteration\n Hamming distance: $(hamming(openings, warmstart_openings))")
    has_sol && write_in_logfile(logfilename, "#openings: $openings")

    if save_cuts_to_disk && !isempty(collected_cuts)
        cutfile = joinpath(cutdir, "$(logfilename)_cuts.jls")
        save_cutpool(cutfile, collected_cuts)
        write_in_logfile(logfilename, "saved $(length(collected_cuts)) cuts to $(basename(cutfile))")
    end
    write_in_logfile(logfilename, "collected $(length(collected_cuts)) cuts")

    write_in_logfile(logfilename, "LP subproblem solves: $(lp_feas_solves + lp_opt_solves) (feasibility: $lp_feas_solves, optimality: $lp_opt_solves)")

    (model=m, iterations=iteration, cuts=collected_cuts, early_stopped=early_stopped[], best_sol=best_sol, best_risk=best_risk, solutions=solutions,
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
                                bridge_to_pocket, bigM_π, bigM_flows; free_π=false, ...)

LP counterpart of [`lostload_optimality_cuts`](@ref). Same candidate set — the bridges
of `bridge_to_pocket` that are not already handled by a feasibility cut — but the cut is
read off the duals of `contingency_subproblem` instead of the closed-form cut-set
argument. With `free_π=true` the subproblem drops `pocket_π_to_0!` and constrains `π`
through the generic connectivity block, so no property of the pocket enters the cut.

Returns `(cuts, n_lp)`: the cuts applied, and how many subproblem LPs were solved.
"""
function lp_lostload_optimality_cuts(m, sink::CutSink, ec, edg, v_0, openbranches, violating_ctg,
                                     bridge_to_pocket, bigM_π, bigM_flows;
                                     free_π::Bool=false,
                                     tight_bigM::Bool=true,
                                     θ_max_bigM::Real=π,
                                     bigM_bound_multiplier::Real=2.0,
                                     atol::Float64=1e-6)
    applied = OBendersCut[]
    n_lp = 0
    for (br, pk) in bridge_to_pocket
        br in violating_ctg && continue
        # Same guard as `create_pklostload_optimality_cut`: an empty pocket sheds nothing,
        # so the three cut modes are compared over an identical candidate set.
        pk.d ≤ 0 && continue
        res = contingency_subproblem(ec, openbranches, br, bridge_to_pocket, bigM_π, bigM_flows;
            free_π, tight_bigM, θ_max_bigM, bigM_bound_multiplier)
        n_lp += 1
        res.is_feasible || continue
        res.obj ≤ atol && continue
        ocut = OBendersCut(edg, br, v_0, [res.reduced_cost[b...] for b in edg], :lostload, res.obj)
        apply_cut!(sink, ocut)
        push!(applied, ocut)
    end
    applied, n_lp
end

function modularity_cut!(m, g, openbranches)
    @constraint(m, sum(1 - m[:v][bus] for bus in edge_labels(g) if bus ∉ openbranches) ≤
                   (ne(g) - length(openbranches)) * sum(m[:v][bus] for bus in openbranches))
end

"""
    bypass_cut_from_pockets(...) -> Vector{Tuple{GenericCut, Set{ELabel}}}

Generate pocket-topology bypass cuts. Returns ALL cuts (applied and deferred) as
`(cut, required_openings)` pairs so the caller can store them in the cut pool.

- Infeasibility cuts (pocket OTS infeasible): applied immediately; `required_openings=Set()`.
- Normal bypass cuts: applied immediately when `pk_openings ⊆ sbs`;
  deferred (stored but not applied) otherwise. Both cases carry
  `required_openings = pk_openings` so future phases can test SBS compatibility.
"""
function bypass_cut_from_pockets(m, sink::CutSink, g, openbranches, sbs,
                                  bridge_to_pocket, known_pockets::Vector{Pocket})
    all_cuts = Vector{Tuple{GenericCut, Set{ELabel}}}()
    for (bridge, pk) in bridge_to_pocket
        isempty(pk.branches) && continue
        pk in known_pockets && continue
        subbusorig = from(bridge) in pk.buses ? from(bridge) : to(bridge)
        subg = create_pocket_subgraph(ElementaryCase(g, subbusorig), openbranches, pk, bridge)
        m_ots = ots(ElementaryCase(subg, subbusorig))
        perimeter     = collect(pk.branches)
        innerbranches = collect(pk.innerbranches)

        if !is_solved_and_feasible(m_ots)
            # Infeasible pocket: at least one perimeter branch must stay closed.
            # -sum_P v_i ≤ -1  ↔  sum_P v_i ≥ 1
            isempty(perimeter) && continue
            infeas_cut = GenericCut(perimeter, falses(length(perimeter)), fill(-1.0, length(perimeter)), -1.0)
            apply_cut!(sink, infeas_cut)
            push!(all_cuts, (infeas_cut, Set{ELabel}()))
            push!(known_pockets, pk)
            continue
        end

        push!(known_pockets, pk)
        pk_openings = getopenings(m_ots)
        req = pk_openings

        edg = vcat(perimeter, innerbranches)
        v_0 = vcat(zeros(Bool, length(perimeter)), [br in pk_openings ? false : true for br in innerbranches])
        # Card(Innerbranches) * sum_branches v_i - sum_innerbranches |v_i - v_i^0| ≥ 0
        # -Card(Innerbranches) * sum_branches v_i + sum_{innerbranches^0} v_i - sum_{innerbranches^1} v_i ≤ -card(Innerbranches^1)
        rc = vcat(
            -length(innerbranches) .* ones(length(perimeter)),
            [br in pk_openings ? 1 : -1 for br in innerbranches],
        )
        cut = GenericCut(edg, v_0, rc, -length([br for br in innerbranches if br ∉ pk_openings]))

        if !issubset(pk_openings, sbs)
            # Inner topology incompatible with current SBS — defer for future phases.
            push!(all_cuts, (cut, req))
            continue
        end

        apply_cut!(sink, cut)
        push!(all_cuts, (cut, req))
    end
    all_cuts
end