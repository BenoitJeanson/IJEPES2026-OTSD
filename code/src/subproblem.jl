
const ZEROSTAB = 0

function contingency_subproblem(ec::ElementaryCase, outages::Set{ELabel}, contingency::ELabel, bridge_to_pocket::Dict{ELabel,Pocket}, bigM_π::Float64, bigM_flows::Float64, margin=0e-2;
    include_base_connectivity::Bool=false,
    reduce_violations::Bool=false,
    monitored_branch::ELabel=("", ""),
    tight_bigM::Bool=false,
    θ_max_bigM::Real=π,
    bigM_bound_multiplier::Real=2.0,
    free_π::Bool=false,
    backend::Backend=default_backend())
    g = ec.g; bus_orig = ec.bus_orig

    function _build_model(g::MetaGraph, outages::Set{ELabel}, contingency::ELabel, bridge_to_pocket, bigM_π, bigM_flows, tight_bigM, θ_max_bigM, bigM_bound_multiplier)

        case = branch_to_case(contingency)

        m = init_model(backend, "", "In subproblem contingency: $contingency")
        openbranches = outages ∪ Set([contingency])
        @variable(m, v[busfrom in labels(g), busto in outneighbor_labels(g, busfrom)])
        foreach(br -> fix(v[br...], !(br in openbranches)), collect(edge_labels(g)))

        if include_base_connectivity
            @variable(m, -nv(g) ≤ c_flows[busfrom in labels(g), outneighbor_labels(g, busfrom)] ≤ nv(g))
            base_connectivity!(m, g, [bus_orig])
        end

        @variable(m, w[[case], busfrom in labels(g), busto in outneighbor_labels(g, busfrom)])
        @constraint(m, w[case, contingency...] == ZEROSTAB)
        @constraint(m, [br in edge_labels(g); br ≠ contingency], w[case, br...] == v[br...])

        @variable(m, load[[case], labels(g)] ≥ 0)
        @variable(m, gen[[case], labels(g)] ≥ 0)
        @expression(m, lostload, sum(max(g[bus], 0) - load[case, bus] for bus in labels(g)))

        @variable(m, flows[[case], busfrom in labels(g), outneighbor_labels(g, busfrom)])
        @variable(m, ϕ[[case], bus in labels(g), 1:1])
        @variable(m, θf[[case], busfrom in labels(g), outneighbor_labels(g, busfrom)])
        @variable(m, θt[[case], busfrom in labels(g), outneighbor_labels(g, busfrom)])
        align_feeder_to_bus_angles!(m, g, [case])
        create_energization_state_variables!(m, g, [case], SubstationConfs(), false)
        align_feeder_to_bus_energization_state!(m, g, [case])
        @variable(m, σ[[case]])

        phase_reference!(m, bus_orig, [case])
        bus_KCL!(m, g, bus_orig, [case])
        balance_N_1cases!(m, g, [case])
        if tight_bigM
            bigM = compute_tight_bigM(g; θ_max=θ_max_bigM, bigM_bound_multiplier=bigM_bound_multiplier)
            ohm!(m, g, [case], bigM)
        else
            ohm!(m, g, [case], bigM_flows)
        end
        # `free_π` drops the graph-derived de-energization and replaces it with the
        # generic connectivity block, so nothing about the pocket enters the LP.
        if free_π
            n_1_connectivity!(m, g, [bus_orig], [case])
        else
            pocket_π_to_0!(m, [case], bridge_to_pocket)
        end
        align_N_1_energization!(m, g, bus_orig, [case])
        if reduce_violations
            @variable(m, s_flows[[from(monitored_branch)], [to(monitored_branch)]] ≥ 0)
            flowslimit_slack!(m, g, monitored_branch)
            @expression(m, n1violation, sum(m[:s_flows]))
            @objective(m, Min, n1violation)
        else
            flowslimits!(m, g, [case], margin)
            @objective(m, Min, m[:lostload])
        end

        m
    end

    m = _build_model(g, outages, contingency, bridge_to_pocket, bigM_π, bigM_flows, tight_bigM, θ_max_bigM, bigM_bound_multiplier)
    set_optimizer_attribute(m, "InfUnbdInfo", 1)
    optimize!(m)

    return benders_subpb_res(m)
end