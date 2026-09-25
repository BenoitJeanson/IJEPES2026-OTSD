"""
    solve_benders_phase(rc, label, hamming, warm_openings, initial_cuts, sbs; kwargs...)

Wraps a single `mastercutpool` call. Returns a NamedTuple with solve results.
All keyword arguments (except `logdir` and `timestamp`) are forwarded to `mastercutpool`,
allowing per-experiment overrides of any solver option.

The result always includes:
  label, hamming, obj, iter, time, secure, openings, cuts, sbs_size,

Default solver options match the EXP-12 baseline configuration.
"""
function solve_benders_phase(rc::RichCase, label::String, hamming::Int,
                              warm_openings::Set{ELabel},
                              initial_cuts::Vector{BendersCut},
                              sbs::Set{ELabel};
                              logdir::String,
                              timestamp::String,
                              ov_lostload_coef::Union{Float64,Dict{ELabel,Float64}} = 1.0,
                              bigM_bound_multiplier::Real       = 2.0,
                              timeout::Float64                  = 0.0,
                              kwargs...)
    logdir_abs = joinpath(dirname(Base.active_project()), "tmp", logdir)
    isdir(logdir_abs) || mkpath(logdir_abs)
    logfile  = "$(logdir)/$(timestamp)_$(label)"
    sbs_size = length(sbs)
    @info "  ▷ $label  H=$hamming  warm=$(length(warm_openings))  cuts=$(length(initial_cuts))  |SBS|=$sbs_size"

    t = @elapsed res = mastercutpool(rc.gc, logfile;
        warmstart_openings         = warm_openings,
        max_hamming                = hamming,
        sbs                        = sbs,
        initial_cuts               = initial_cuts,
        ov_lostload_coef           = ov_lostload_coef,
        bigM_bound_multiplier      = bigM_bound_multiplier,
        timeout                    = timeout,
        kwargs...)

    t_r           = round(t; digits=1)
    solved        = is_solved_and_feasible(res.model)
    has_incumbt   = primal_status(res.model) == MOI.FEASIBLE_POINT
    lp_solves     = hasproperty(res, :lp_solves) ? res.lp_solves : 0
    # A run stopped by its own time cap still carries a usable incumbent; keep it rather
    # than reporting Inf. Only reachable when a cap was asked for (`timeout > 0`).
    capped        = timeout > 0.0 && !solved && has_incumbt
    stop_tag      = capped ? " ⏱" : ""

    if solved || capped
        ov     = round(objective_value(res.model); digits=3)
        ops    = getopenings(res.model)
        secure = isempty(reduced_sa(rc, Set(ops)))
        @info "    ✓ obj=$ov  iter=$(res.iterations)  t=$(t_r)s  secure=$secure$stop_tag"
        return (label=label, hamming=hamming, obj=ov, iter=res.iterations,
                time=t_r, secure=secure, openings=ops, cuts=res.cuts,
                sbs_size=sbs_size, solutions=res.solutions,
                lp_solves=lp_solves, capped=capped)
    else
        @error "    ✗ not solved"
        return (label=label, hamming=hamming, obj=Inf, iter=res.iterations,
                time=t_r, secure=false, openings=warm_openings, cuts=res.cuts,
                sbs_size=sbs_size, solutions=Dict{Set{ELabel}, Float64}(),
                lp_solves=lp_solves, capped=false)
    end
end
