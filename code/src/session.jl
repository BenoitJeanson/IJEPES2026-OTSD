# ── TeeLogger: routes every log message to both a stream and a file ──────────
struct TeeLogger <: AbstractLogger
    console::AbstractLogger
    file::AbstractLogger
end
Logging.shouldlog(l::TeeLogger, level, _module, group, id) =
    Logging.shouldlog(l.console, level, _module, group, id) ||
    Logging.shouldlog(l.file, level, _module, group, id)
Logging.min_enabled_level(l::TeeLogger) =
    min(Logging.min_enabled_level(l.console), Logging.min_enabled_level(l.file))
Logging.catch_exceptions(::TeeLogger) = false
function Logging.handle_message(l::TeeLogger, level, message, _module, group, id, file, line; kwargs...)
    Logging.handle_message(l.console, level, message, _module, group, id, file, line; kwargs...)
    Logging.handle_message(l.file, level, message, _module, group, id, file, line; kwargs...)
end

# ── SessionState ──────────────────────────────────────────────────────────────
mutable struct SessionState
    all_records::Vector{CutRecord}
    best_openings::Set{ELabel}
    best_obj::Float64
    accumulated_sbs::Set{ELabel}
    iteration::Int
end

# ── Core Benders iteration loop ───────────────────────────────────────────────
function run_benders_iterations!(state::SessionState, rc::RichCase,
                                  all_branches::Vector{ELabel},
                                  label_prefix::String;
                                  H::Int,
                                  k::Int,
                                  SBS_SA_depth::Int,
                                  seed::Int,
                                  contingencies_in_master::Union{Nothing,Vector{ELabel}},
                                  logdir::String,
                                  timestamp::String,
                                  logdir_abs::String,
                                  warmstart_hop_sbs::Int = 0,
                                  use_cut_inheritance::Bool = true,
                                  timeout::Float64 = 0.0,
                                  solver_kwargs::NamedTuple = NamedTuple(),
)
    edm = warmstart_hop_sbs > 0 ? edge_distance_map(rc.gc.g) : nothing
    results_log = String[]
    iters       = NamedTuple[]
    cumul_time  = 0.0

    for _ in 1:k
        state.iteration += 1
        iter  = state.iteration
        label = "$(label_prefix)_i$(lpad(iter, 3, '0'))"

        initial_cuts  = use_cut_inheritance ? extract_cuts_for_phase(state.all_records) : BendersCut[]

        @info "── $label_prefix  iter $iter  |SBS|=$(length(state.accumulated_sbs))  |cuts|=$(length(initial_cuts))"
        r = solve_benders_phase(rc, label, H, state.best_openings, initial_cuts, state.accumulated_sbs;
            logdir, timestamp,
            contingencies_in_master,
            bigM_bound_multiplier = 2.0,
            seed,
            timeout,
            solver_kwargs...)

        append!(state.all_records, r.cuts)
        cumul_time += r.time

        improved = r.obj < state.best_obj - 1e-4
        if improved
            @info "  ↑ improvement: $(state.best_obj) → $(r.obj)"
            state.best_obj      = r.obj
            state.best_openings = r.openings
        end

        prev_sbs = length(state.accumulated_sbs)
        union!(state.accumulated_sbs, sa_induced_followed(rc, state.best_openings, SBS_SA_depth))
        if warmstart_hop_sbs > 0
            union!(state.accumulated_sbs,
                extend_sbs_by_hops(Set(state.best_openings), all_branches, edm; d=warmstart_hop_sbs))
        end
        @info "  |SBS| $prev_sbs → $(length(state.accumulated_sbs))"

        open(joinpath(logdir_abs, "$(timestamp)_i$(lpad(iter, 3, '0'))_state.jls"), "w") do io
            Serialization.serialize(io, state)
        end

        push!(iters, (iteration=iter, H=H, obj=r.obj, n_openings=length(r.openings),
                      inner_iter=r.iter, time=r.time, cumul_time=round(cumul_time; digits=1),
                      secure=r.secure, sbs_size=length(state.accumulated_sbs), improved=improved,
                      lp_solves=r.lp_solves, capped=r.capped, openings=r.openings))

        push!(results_log,
            "i=$(lpad(iter,3))  H=$H  obj=$(rpad(r.obj,8))  openings=$(lpad(length(r.openings),3))  " *
            "inner_iter=$(lpad(r.iter,5))  t=$(lpad(r.time,7))s  cumul=$(lpad(round(cumul_time; digits=1),8))s  " *
            "secure=$(r.secure)  |SBS|=$(length(state.accumulated_sbs))  impr=$(improved)  " *
            "lp=$(r.lp_solves)$(r.capped ? "  ⏱capped" : "")")
    end

    (results_log=results_log, cumul_time=cumul_time, iters=iters)
end
