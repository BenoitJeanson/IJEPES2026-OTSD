# ── Entry point ───────────────────────────────────────────────────────────────

"""
    solve_otsd(case; kwargs...) -> NamedTuple

Solve OTSD on a PGLib case and return the result of the local search.

    solve_otsd("case118"; tlf = 1.5, H = 3, d_viol = 2, d_sol = 0, k = 3,
               warm_start = …, backend = GurobiBackend())

Keyword arguments:

  * `tlf`        — thermal limit factor (1.5 for IEEE-118, 1.0 for IEEE-57 in the paper).
  * `H`          — Hamming bound on the master. `0` removes the bound.
  * `d_viol`     — hops traced along the flow path from each violation into the SBS.
  * `d_sol`      — topological hops around the incumbent added to the SBS.
  * `k`          — number of local-search phases.
  * `warm_start` — the openings the search starts from. Defaults to the closed network,
                   which is valid but slower than the heuristic seed the paper uses.
  * `embed`      — contingencies embedded in the master from the outset.
  * `inherit`    — carry cuts across phases.
  * `backend`    — [`GurobiBackend`](@ref) or [`SCIPBackend`](@ref); both run the
                   same algorithm, the second without a licence.
  * `seed`       — solver seed.
  * `timeout`    — wall-clock cap in seconds for the whole search; `0` for none.

Any further keyword arguments reach the master problem unchanged; the ablation study
uses `ablate_screening` and `cf_cuts` this way.

Returns `(; objective, openings, secure, phases, wall, lp_solves, benders_iterations)`.
"""
function solve_otsd(case::String;
                    tlf::Union{Nothing,Real} = nothing,
                    H::Int = 3,
                    d_viol::Int = 2,
                    d_sol::Int = 0,
                    k::Int = 3,
                    warm_start::Set{ELabel} = Set{ELabel}(),
                    embed::Union{Nothing,Vector{ELabel}} = nothing,
                    inherit::Bool = true,
                    backend::Backend = default_backend(),
                    seed::Int = 0,
                    timeout::Real = 0.0,
                    logdir::String = "otsd",
                    label::String = "OTSD",
                    solver_kwargs...)
    rc = load_case(case; tlf)
    all_branches = collect(ELabel, edge_labels(rc.gc.g))
    sbs = H > 0 ? sa_induced_followed(rc, warm_start, d_viol) : Set(all_branches)

    state = SessionState(CutRecord[], copy(warm_start), Inf, sbs, 0)
    timestamp = Dates.format(now(), "yyyy-mm-dd_HHMMSS")
    logdir_abs = joinpath(pwd(), "tmp", logdir)
    mkpath(logdir_abs)

    wall = @elapsed itr = run_benders_iterations!(
        state, rc, all_branches, label;
        H, k, SBS_SA_depth = d_viol, seed,
        contingencies_in_master = embed,
        logdir, timestamp, logdir_abs,
        warmstart_hop_sbs = d_sol,
        use_cut_inheritance = inherit,
        timeout = Float64(timeout),
        solver_kwargs = (; backend, solver_kwargs...))

    (; objective = state.best_obj,
       openings = state.best_openings,
       secure = isempty(reduced_sa(rc, state.best_openings)),
       phases = itr.iters,
       wall = round(wall; digits = 1),
       lp_solves = sum(it.lp_solves for it in itr.iters; init = 0),
       benders_iterations = sum(it.inner_iter for it in itr.iters; init = 0))
end
