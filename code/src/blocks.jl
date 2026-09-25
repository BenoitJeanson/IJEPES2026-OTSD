
"""
    branch_to_case(case::ELabel)

Convert a branch tuple `(from, to)` into the contingency case identifier string.

Returns:
- `String`: identifier of the form `"from-to"`.
"""
branch_to_case(case::ELabel) = case[1] * BASECASEID * case[2]

"""
    case_to_branch(st::String)

Convert a contingency identifier string back into a branch tuple.

Returns:
- `ELabel`: branch tuple `(from, to)`.
"""
function case_to_branch(st::String)
    parts = split(st, BASECASEID)
    length(parts) == 2 || throw(ArgumentError("invalid case identifier: $st"))
    (parts[1], parts[2])
end

"""
    initlogpath(logfilename::String)

Build the absolute path of the log file in the project `tmp/` directory.

Arguments:
- `logfilename`: base log filename (without extension).

Returns:
- `String`: path `tmp/<logfilename>.log`.
"""
function initlogpath(logfilename::String)
    filename = isempty(logfilename) ? "log" : logfilename
    joinpath(dirname(Base.active_project()), "tmp", "$filename.log")
end

"""
    write_in_logfile(logfilename::String, message::String)

Append a message to the configured log file when logging is enabled.

Arguments:
- `logfilename`: logical log file name (without `.log` suffix).
- `message`: line to append.

Returns:
- `Nothing`.

Side effects:
- Appends text to `tmp/<logfilename>.log` when `logfilename != ""`.
"""
function write_in_logfile(logfilename::String, message::String)
    logfilename == "" && return
    open(initlogpath(logfilename), "a") do f
        println(f, message)
    end
end

"""
    init_model(backend::Backend, logfilename::String="", message::String="")

An empty model on `backend`, silent on the console, logging to `logfilename` when one
is given. `message` is written as a header before the solver attaches to the file.
"""
function init_model(backend::Backend, logfilename::String = "", message::String = "")
    log_path = ""
    if logfilename != ""
        log_path = initlogpath(logfilename)
        open(log_path, "a") do f
            println(f, "________________________________________________________________")
            message ≠ "" && print(f, message)
        end
    end
    new_model(backend; log_path)
end

"""
    base_connectivity!(m::Model, g::MetaGraph, fixed_buses::Vector{VLabel}, big_M; ...)

Add base-case connectivity constraints, optionally with bus-split activation logic.

Arguments:
- `m`: JuMP model.
- `g`: network graph.
- `fixed_buses`: buses excluded from unit-flow connectivity equations.

Returns:
- `Nothing`.

Side effects:
- Adds connectivity and branch-flow linking constraints to `m`.
"""
function base_connectivity!(
    m::Model,
    g::MetaGraph,
    fixed_buses::Vector{VLabel},
    allow_branch_openings::Bool = true,
)

    big_M = ne(g) + 1

    for bus in labels(g)
        bus in fixed_buses && continue
        @constraint(
            m,
            sum(m[:c_flows][busfrom, bus] for busfrom in inneighbor_labels(g, bus)) -
            sum(m[:c_flows][bus, busto] for busto in outneighbor_labels(g, bus)) == 1
        )
    end

    if allow_branch_openings
        @constraint(m, m[:c_flows][:, :] .≤ m[:v][:, :] .* big_M)
        @constraint(m, m[:c_flows][:, :] .≥ -m[:v][:, :] .* big_M)
    end
end

"""
    compute_tight_bigM(g; θ_max=π, bigM_bound_multiplier=2.0)

Return per-branch tight big-M values for the `flows!` constraints.
- `bigM_ohm[br]`   = 2·θ_max·b_e             — Ohm's law decoupling when branch is open.
- `bigM_bound[br]` = bigM_bound_multiplier·p_max_e — flow-zeroing when branch is open.

`bigM_bound_multiplier > 1` (default 2) is deliberate: setting it to 1 would make
`|f| ≤ p_max·w` a hard thermal limit for closed branches (w=1), collapsing the slack
variable in `flowslimit_slack!` and causing infeasibility in ρ>0 (reduce-violations) mode.
With multiplier=2 the big-M gives headroom for the LP while the slack constraint
`|f| ≤ p_max + s_flows` remains the actual binding thermal constraint.

Note: bigM_ohm is practically valid without explicit angle bounds — the LP has no incentive
to drive de-energized bus angles to extreme values (no objective term on ϕ; energized buses
anchored by power balance). See EXP-07.
"""
function compute_tight_bigM(
    g::MetaGraph;
    θ_max::Real = π,
    bigM_bound_multiplier::Real = 2.0,
)
    bigM_ohm   = Dict(br => 2.0 * θ_max * g[br...].b for br in edge_labels(g))
    bigM_bound = Dict(br => bigM_bound_multiplier * Float64(g[br...].p_max) for br in edge_labels(g))
    (bigM_ohm = bigM_ohm, bigM_bound = bigM_bound)
end

"""Convenience overload of `base_connectivity!` for one fixed reference bus."""
function base_connectivity!(
    m::Model,
    g::MetaGraph,
    bus_orig::String,
    allow_branch_openings::Bool = true,
)
    base_connectivity!(m, g, [bus_orig], allow_branch_openings)
end

"""
    phase_reference!(m::Model, fixed_buses::Vector{VLabel}, fixed_ϕ::Vector{Float64}, cases::Vector{String})

Fix phase-angle reference values for selected buses over all cases.

Returns:
- `Nothing`.

Side effects:
- Adds equality constraints on `m[:ϕ]`.
"""
function phase_reference!(
    m::Model,
    fixed_buses::Vector{VLabel},
    fixed_ϕ::Vector{Float64},
    cases::Vector{String},
)
    @constraint(
        m,
        [c in cases, fixed_id in 1:length(fixed_buses)],
        m[:ϕ][c, fixed_buses[fixed_id], 1] == fixed_ϕ[fixed_id]
    )
end

"""Convenience overload of `phase_reference!` with zero phase at one bus."""
function phase_reference!(m::Model, bus_orig::VLabel, cases::Vector{String})
    phase_reference!(m, [bus_orig], [0.0], cases)
end

# Tight big-M variant: per-branch M dicts produced by compute_tight_bigM().
# bigM_ohm[br]   bounds the Ohm's law decoupling constraint (active when branch open).
# bigM_bound[br] bounds the flow-zeroing constraint (active when branch open).
function flows!(
    m::Model,
    g::MetaGraph,
    fixed_buses::Vector{VLabel},
    fixed_ϕ::Vector{Float64},
    cases::Vector{String},
    bigM_ohm::Dict,
    bigM_bound::Dict,
)
    @constraint(
        m,
        [c in cases, fixed_id in 1:length(fixed_buses)],
        m[:ϕ][c, fixed_buses[fixed_id]] == fixed_ϕ[fixed_id]
    )

    @constraint(
        m,
        [c in cases, busfrom in labels(g), busto in outneighbor_labels(g, busfrom)],
        m[:flows][c, busfrom, busto] ≤
        m[:w][c, busfrom, busto] * bigM_bound[(busfrom, busto)]
    )
    @constraint(
        m,
        [c in cases, busfrom in labels(g), busto in outneighbor_labels(g, busfrom)],
        m[:flows][c, busfrom, busto] ≥
        -m[:w][c, busfrom, busto] * bigM_bound[(busfrom, busto)]
    )

    @constraint(
        m,
        [c in cases, busfrom in labels(g), busto in outneighbor_labels(g, busfrom)],
        m[:flows][c, busfrom, busto] -
        g[busfrom, busto].b * (m[:ϕ][c, busto] - m[:ϕ][c, busfrom]) ≤
        bigM_ohm[(busfrom, busto)] * (1 - m[:w][c, busfrom, busto])
    )
    @constraint(
        m,
        [c in cases, busfrom in labels(g), busto in outneighbor_labels(g, busfrom)],
        m[:flows][c, busfrom, busto] -
        g[busfrom, busto].b * (m[:ϕ][c, busto] - m[:ϕ][c, busfrom]) ≥
        -bigM_ohm[(busfrom, busto)] * (1 - m[:w][c, busfrom, busto])
    )

    @constraint(
        m,
        [c in cases, bus in labels(g); !(bus in fixed_buses)],
        m[:load][c, bus] - m[:gen][c, bus] ==
        sum(m[:flows][c, busfrom, bus] for busfrom in inneighbor_labels(g, bus)) -
        sum(m[:flows][c, bus, busto] for busto in outneighbor_labels(g, bus))
    )
end

function flows!(
    m::Model,
    g::MetaGraph,
    bus_orig::String,
    cases::Vector{String},
    bigM_ohm::Dict,
    bigM_bound::Dict,
)
    flows!(m, g, [bus_orig], [0.0], cases, bigM_ohm, bigM_bound)
end

# NOTE: body of flows!(…, w::Function) was truncated in the merge conflict.
# Verify this matches the IOS-revisited version.
function flows!(
    m::Model,
    g::MetaGraph,
    fixed_buses::Vector{VLabel},
    fixed_ϕ::Vector{Float64},
    cases::Vector{String},
    w::Function,
)
    @constraint(
        m,
        [c in cases, fixed_id in 1:length(fixed_buses)],
        m[:ϕ][c, fixed_buses[fixed_id]] == fixed_ϕ[fixed_id]
    )
    @constraint(
        m,
        [c in cases, busfrom in labels(g), busto in outneighbor_labels(g, busfrom)],
        m[:flows][c, busfrom, busto] ==
        w(c, busfrom, busto) * g[busfrom, busto].b * (m[:ϕ][c, busto] - m[:ϕ][c, busfrom])
    )
    @constraint(
        m,
        [c in cases, bus in labels(g); !(bus in fixed_buses)],
        m[:load][c, bus] - m[:gen][c, bus] ==
        sum(m[:flows][c, busfrom, bus] for busfrom in inneighbor_labels(g, bus)) -
        sum(m[:flows][c, bus, busto] for busto in outneighbor_labels(g, bus))
    )
end

"""
    bus_KCL!(m::Model, g::MetaGraph, fixed_buses::Vector{VLabel}, cases::Vector{String})

Add standard KCL constraints for all non-fixed buses.

Returns:
- `Nothing`.

Side effects:
- Adds KCL constraints to `m`.
"""
function bus_KCL!(
    m::Model,
    g::MetaGraph,
    fixed_buses::Vector{VLabel},
    cases::Vector{String},
)
    @constraint(
        m,
        [c in cases, bus in labels(g); !(bus in fixed_buses)],
        m[:load][c, bus] - m[:gen][c, bus] ==
        sum(m[:flows][c, busfrom, bus] for busfrom in inneighbor_labels(g, bus)) -
        sum(m[:flows][c, bus, busto] for busto in outneighbor_labels(g, bus))
    )
end

"""Convenience overload of `bus_KCL!` for one fixed bus."""
function bus_KCL!(m::Model, g::MetaGraph, bus_orig::VLabel, cases::Vector{String})
    bus_KCL!(m, g, [bus_orig], cases)
end

"""
    ohm!(m::Model, g::MetaGraph, cases::Vector{String}, big_M)

Add big-M DC Ohm-law constraints and branch-status flow gating.

Returns:
- `Nothing`.

Side effects:
- Adds constraints to `m`.
"""
function ohm!(m::Model, g::MetaGraph, cases::Vector{String}, big_M::NamedTuple)
    for c in cases, br in edge_labels(g)
        bigM_ohm = big_M.bigM_ohm[br]
        bigM_bound = big_M.bigM_bound[br]
        @constraint(m, m[:flows][c, br...] .≤ m[:w][c, br...] .* bigM_bound)
        @constraint(m, m[:flows][c, br...] .≥ -m[:w][c, br...] .* bigM_bound)

        @constraint(
            m,
            m[:flows][c, br...] - sum(g[br...].b * (m[:θt][c, br...] - m[:θf][c, br...])) ≤
            bigM_ohm * (1 - m[:w][c, br...])
        )
        @constraint(
            m,
            m[:flows][c, br...] - sum(g[br...].b * (m[:θt][c, br...] - m[:θf][c, br...])) ≥
            -bigM_ohm * (1 - m[:w][c, br...])
        )
    end
end

"""
    ohm!(m::Model, g::MetaGraph, cases::Vector{String}, w::Function)

Add Ohm-law constraints when branch status is provided by callback/function `w`.

Returns:
- `Nothing`.

Side effects:
- Adds constraints to `m`.
"""
function ohm!(m::Model, g::MetaGraph, cases::Vector{String}, w::Function)
    @constraint(
        m,
        [c in cases, busfrom in labels(g), busto in outneighbor_labels(g, busfrom)],
        m[:flows][c, busfrom, busto] ==
        w(c, busfrom, busto) *
        (g[busfrom, busto].b * (m[:θt][c, busfrom, busto] - m[:θf][c, busfrom, busto]))
    )
end

"""
    flowslimits!(m::Model, g::MetaGraph, cases::Vector{String}, margin=0.; branches=nothing)

Add symmetric thermal flow limits, optionally on a subset of branches.

Returns:
- `Nothing`.

Side effects:
- Adds limit constraints to `m`.
"""
function flowslimits!(
    m::Model,
    g::MetaGraph,
    cases::Vector{String},
    margin = 0.0;
    branches::Union{Nothing,Vector{ELabel}} = nothing,
)
    branchiterator = isnothing(branches) ? edge_labels(g) : branches
    @constraint(
        m,
        flowLimitsP[c in cases, (busfrom, busto) in branchiterator],
        m[:flows][c, busfrom, busto] ≤ g[busfrom, busto].p_max - margin
    )
    @constraint(
        m,
        flowLimitsM[c in cases, (busfrom, busto) in branchiterator],
        -m[:flows][c, busfrom, busto] ≤ g[busfrom, busto].p_max - margin
    )
end

"""
    flowslimit_slack!(m::Model, g::MetaGraph, branches=nothing)

Relax flow limits with nonnegative branch slack variables.

Returns:
- `Nothing`.

Side effects:
- Adds slack-augmented limit constraints to `m`.
"""
function flowslimit_slack!(
    m::Model,
    g::MetaGraph,
    branches::Union{Nothing,Vector{ELabel}} = nothing,
)
    _branches = isnothing(branches) ? edge_labels(g) : branches
    @constraint(
        m,
        [br in _branches],
        m[:flows][:, br...] .≤ g[br...].p_max .+ m[:s_flows][br...]
    )
    @constraint(
        m,
        [br in _branches],
        -m[:flows][:, br...] .≤ g[br...].p_max .+ m[:s_flows][br...]
    )
end

"""Convenience overload of `flowslimit_slack!` for a single branch."""
function flowslimit_slack!(m::Model, g::MetaGraph, branch::ELabel)
    flowslimit_slack!(m, g, [branch])
    # @constraint(m, m[:flows][:, branch...] .≤ g[branch...].p_max + m[:s_flows][branch...])
    # @constraint(m, -m[:flows][:, branch...] .≤ g[branch...].p_max + m[:s_flows][branch...])
    # @constraint(m, m[:s_flows] ≥ 0)
end

"""
    create_energization_state_variables!(m::Model, g::MetaGraph, n_1cases, is_π_binary::Bool; no_deenergization::Bool=false)

Create bus/feeder energization variables (or fixed expressions when de-energization is disabled).

Returns:
- `Nothing`.

Side effects:
- Adds `:π`, `:ψf`, and `:ψt` containers to `m`.
"""
function create_energization_state_variables!(
    m::Model,
    g::MetaGraph,
    n_1cases,
    is_π_binary::Bool;
    no_deenergization::Bool = false,
)
    if no_deenergization
        @expression(
            m,
            π[n_1cases, bus in labels(g)],
            1
        )
        @expression(
            m,
            ψf[n_1cases, busfrom in labels(g), outneighbor_labels(g, busfrom)],
            1
        )
        @expression(
            m,
            ψt[n_1cases, busfrom in labels(g), outneighbor_labels(g, busfrom)],
            1
        )
        return
    end
    if is_π_binary
        @variable(
            m,
            π[n_1cases, bus in labels(g)],
            Bin
        )
        @variable(
            m,
            ψf[n_1cases, busfrom in labels(g), outneighbor_labels(g, busfrom)],
            Bin
        )
        @variable(
            m,
            ψt[n_1cases, busfrom in labels(g), outneighbor_labels(g, busfrom)],
            Bin
        )
    else
        @variable(
            m,
            0 ≤ π[n_1cases, bus in labels(g)] ≤ 1
        )#, Bin)
        @variable(
            m,
            0 ≤ ψf[n_1cases, busfrom in labels(g), outneighbor_labels(g, busfrom)] ≤ 1
        )#, Bin)
        @variable(
            m,
            0 ≤ ψt[n_1cases, busfrom in labels(g), outneighbor_labels(g, busfrom)] ≤ 1
        )#, Bin)
    end
end

# ── Feeder-to-bus alignment ───────────────────────────────────────────────────
#
# Each branch endpoint carries its own angle (`θf`, `θt`) and energization state
# (`ψf`, `ψt`) variable, tied here to the variable of the bus it lands on. The
# indirection exists so that a substation-reconfiguration extension can split a bus
# and send its feeders to different sub-buses; nothing in this package splits a bus,
# so the tie is an equality and the model is the one the paper reports.

"Tie a feeder-level variable to the bus-level variable it belongs to."
_tie_feeder_to_bus!(m::Model, cases, branch, bus, feeder_var::Symbol, bus_var::Symbol) =
    @constraint(m, [c in cases], m[feeder_var][c, branch...] .== m[bus_var][c, bus, 1])

"Tie both endpoint angles of every branch to their bus angle."
function align_feeder_to_bus_angles!(m::Model, g::MetaGraph, cases)
    for busfrom in labels(g), busto in outneighbor_labels(g, busfrom)
        br = (busfrom, busto)
        _tie_feeder_to_bus!(m, cases, br, busfrom, :θf, :ϕ)
        _tie_feeder_to_bus!(m, cases, br, busto, :θt, :ϕ)
    end
end

"Tie both endpoint energization states of every branch to their bus state."
function align_feeder_to_bus_energization_state!(m::Model, g::MetaGraph, n_1cases)
    for busfrom in labels(g), busto in outneighbor_labels(g, busfrom)
        br = (busfrom, busto)
        _tie_feeder_to_bus!(m, n_1cases, br, busfrom, :ψf, :π)
        _tie_feeder_to_bus!(m, n_1cases, br, busto, :ψt, :π)
    end
end

"""
    warmstart_openings!(m::Model, g::MetaGraph, warmstart_openings::Vector{ELabel})

Set start values for branch-opening binary variables.

Returns:
- `Nothing`.

Side effects:
- Mutates optimizer start values in `m`.
"""
function warmstart_openings!(m::Model, g::MetaGraph, warmstart_openings::Set{ELabel})
    isempty(warmstart_openings) && return
    for br in edge_labels(g)
        set_start_value(m[:v][br...], !(br in warmstart_openings))
    end
end

"""
    balance_basecase!(m::Model, g::MetaGraph, basecase::String)

Fix base-case load and generation from net bus injections.

Returns:
- `Nothing`.

Side effects:
- Adds equalities on `m[:load]` and `m[:gen]`.
"""
function balance_basecase!(m::Model, g::MetaGraph, basecase::String)
    for bus in labels(g)
        p = g[bus]
        @constraint(m, m[:load][basecase, bus] == max(0, p))
        @constraint(m, m[:gen][basecase, bus] == -min(0, p))
    end
end

"""
    balance_N_1cases!(m::Model, g::MetaGraph, n_1cases::Vector{String}, σ_max=2.0)

Add N-1 load/generation balance constraints using energization variables and scaling
variable `σ`.

Returns:
- `Nothing`.

Side effects:
- Adds balance and big-M coupling constraints to `m`.
"""
function balance_N_1cases!(
    m::Model,
    g::MetaGraph,
    n_1cases::Vector{String},
    σ_max = 2.0,
)
    @constraint(m, m[:σ][:] .≤ σ_max)
    @constraint(m, [c in n_1cases], sum(m[:load][c, :]) == sum(m[:gen][c, :]))
    for bus in labels(g)
        p = g[bus]
        big_M = σ_max * abs(p)
        for c in n_1cases
            if p < 0
                @constraint(m, m[:load][c, bus] == 0)
                @constraint(
                    m,
                    m[:gen][c, bus] + m[:σ][c] * p ≤ big_M * (1 - m[:π][c, bus])
                )
                @constraint(
                    m,
                    -(m[:gen][c, bus] + m[:σ][c] * p) ≤ big_M * (1 - m[:π][c, bus])
                )
                @constraint(m, m[:gen][c, bus] ≤ big_M * m[:π][c, bus])
                @constraint(m, -m[:gen][c, bus] ≤ big_M * m[:π][c, bus])
            else
                @constraint(m, m[:gen][c, bus] == 0)
                @constraint(m, m[:load][c, bus] == m[:π][c, bus] * p)
            end
        end
    end
end

"""
    align_N_1_energization!(m::Model, g::MetaGraph, fixed_buses::Vector{VLabel}, n_1cases::Vector{String})

Couple endpoint energization variables and force fixed buses energized in N-1 cases.

Returns:
- `Nothing`.

Side effects:
- Adds energization alignment constraints to `m`.
"""
function align_N_1_energization!(
    m::Model,
    g::MetaGraph,
    fixed_buses::Vector{VLabel},
    n_1cases::Vector{String},
)
    @constraint(m, [bus in fixed_buses], m[:π][:, bus] .== 1)
    for busfrom in labels(g), busto in outneighbor_labels(g, busfrom)
        @constraint(
            m,
            [c in n_1cases],
            m[:ψf][c, busfrom, busto] ≥
            m[:ψt][c, busfrom, busto] + m[:w][c, busfrom, busto] - 1
        )
        @constraint(
            m,
            [c in n_1cases],
            m[:ψt][c, busfrom, busto] ≥
            m[:ψf][c, busfrom, busto] + m[:w][c, busfrom, busto] - 1
        )
    end
end

"""Convenience overload of `align_N_1_energization!` for one fixed bus."""
function align_N_1_energization!(
    m::Model,
    g::MetaGraph,
    bus_orig::String,
    n_1cases::Vector{String},
)
    align_N_1_energization!(m, g, [bus_orig], n_1cases)
end

"""
    pocket_π_to_0!(m::Model, n_1cases::Vector{String}, bridge_to_pocket::Dict{ELabel,Pocket})

Add constraints forcing pocket buses to de-energize when linking branches are open.

Returns:
- `Nothing`.

Side effects:
- Adds constraints on `m[:π]` and `m[:w]`.
"""
function pocket_π_to_0!(
    m::Model,
    n_1cases::Vector{String},
    bridge_to_pocket::Dict{ELabel,Pocket},
)
    for case in n_1cases
        cbr = case_to_branch(case)
        !haskey(bridge_to_pocket, cbr) && continue
        pocket = bridge_to_pocket[cbr]
        @constraint(
            m,
            sum(m[:π][case, bus] for bus in pocket.buses) ≤
            length(pocket.buses) * sum(m[:w][case, br...] for br in pocket.branches)
        )
    end
end

"""
    n_1_connectivity!(m::Model, g::MetaGraph, fixed_buses::Vector{VLabel}, n_1cases::Vector{String})

Tie `m[:π]` to topological connectivity in the N-1 cases with a fictitious
single-commodity flow gated by `m[:w]`: every energized bus injects one unit, the
`fixed_buses` absorb, and nothing crosses an open branch. A bus in an island detached
from a fixed bus therefore cannot be energized.

This is the topology-agnostic counterpart of [`pocket_π_to_0!`](@ref), which obtains the
same de-energization from a precomputed bridge/pocket decomposition.

Returns:
- `Nothing`.

Side effects:
- Declares `m[:n_1c_flows]` and adds constraints on `m[:π]`, `m[:w]`.
"""
function n_1_connectivity!(
    m::Model,
    g::MetaGraph,
    fixed_buses::Vector{VLabel},
    n_1cases::Vector{String},
)
    @variable(
        m,
        -nv(g) ≤
        n_1c_flows[n_1cases, busfrom in labels(g), outneighbor_labels(g, busfrom)] ≤
        nv(g)
    )
    big_M = ne(g)
    @constraint(
        m,
        [c in n_1cases, bus in labels(g); !(bus in fixed_buses)],
        m[:π][c, bus] ==
        sum(sign * m[:n_1c_flows][c, br...] for (br, sign) in incident_signed(g, bus))
    )
    @constraint(
        m,
        [c in n_1cases, br in edge_labels(g)],
        m[:n_1c_flows][c, br...] ≤ big_M * m[:w][c, br...]
    )
    @constraint(
        m,
        [c in n_1cases, br in edge_labels(g)],
        m[:n_1c_flows][c, br...] ≥ -big_M * m[:w][c, br...]
    )
end

"""
    violated_branches(model::Model, g::MetaGraph, case::String)

Return branches whose absolute flow exceeds thermal rating in `case`.

Returns:
- `Vector{ELabel}`.
"""
function violated_branches(model::Model, g::MetaGraph, case::String)
    branches = Vector{ELabel}()
    for br in edge_labels(g)
        if abs(value(model[:flows][case, br...])) > g[br...].p_max
            push!(branches, br)
        end
    end
    branches
end

"""
    max_overload(model::Model, g::MetaGraph, case::String, branches::Vector{ELabel})

Return the maximum normalized overload and corresponding branch in a case.

Returns:
- Tuple `(Float64, ELabel)`.
"""
function max_overload(model::Model, g::MetaGraph, case::String, branches::Vector{ELabel})
    maxval, index =
        findmax(br -> abs(value(model[:flows][case, br...])) / g[br...].p_max, branches)
    maxval, branches[index]
end

"""
    getopenings(m::Model)

Return the set of open branches from binary decision values.

Returns:
- `Set{ELabel}`.
"""
getopenings(m::Model) =
    haskey(m, :v) ? Set(br for br in eachindex(m[:v]) if round(value(m[:v][br])) == 0) :
    Set{ELabel}()

