# ── Component ablation ────────────────────────────────────────────────────────
#
# Five acceleration strategies carry the reported speedup. Each configuration below
# disables exactly one of them, so the cost of removing it can be read off directly.
# Runs are strictly serial on a fixed thread count: wall time is the quantity being
# compared, and thread contention would distort it.
#
# Results are files. Every finished run writes its own `result.json`, so the campaign
# can be interrupted and resumed, and the summary is a pure function of what is on disk.

using OTSD
using Dates, JSON3, Logging, Printf

include(joinpath(@__DIR__, "instances.jl"))

"""
Configurations. `REF` disables nothing; each `NO-*` disables one component.
`embed`, `inherit` and `local_search` act on the search, the rest reach the master.
"""
const CONFIGS = Dict(
    "REF"            => (; note = "nothing disabled — reference"),
    "NO-SCREEN"      => (; ablate_screening = true, note = "1: security screening"),
    "NO-CF-POCKET"   => (; cf_cuts = :lp_pocket,    note = "2: closed-form cuts → LP, pocket π"),
    "NO-CF-FREE"     => (; cf_cuts = :lp_free,      note = "2: closed-form cuts → LP, free π"),
    "NO-EMBED"       => (; embed = false,           note = "3: contingency embedding"),
    "NO-INHERIT"     => (; inherit = false,         note = "4: cut inheritance"),
    "NO-LOCALSEARCH" => (; local_search = false,    note = "5: local search"),
)

const CONFIG_ORDER = ["REF", "NO-SCREEN", "NO-CF-POCKET", "NO-CF-FREE",
                      "NO-EMBED", "NO-INHERIT", "NO-LOCALSEARCH"]

"""
    run_cell(instance, config; backend, seed, cap, outdir) -> Dict

One row of the study. Returns the record it also writes to disk.

Ablating the local search removes both of its knobs at once — the Hamming ball and
the switchable branch set — leaving a single global Benders solve. That is a
different operating point, not `H = 3` with a switch flipped, so the record says so.
"""
function run_cell(inst::Instance, cfgid::String;
                  backend::Backend, seed::Int = 0, cap::Real = 900,
                  outdir::String = "results")
    opts = CONFIGS[cfgid]
    p = PAPER_OPERATING_POINT

    local_search = get(opts, :local_search, true)
    embed = get(opts, :embed, true)
    inherit = get(opts, :inherit, true)
    solver_kwargs = Base.structdiff(opts, NamedTuple{(:note, :embed, :inherit, :local_search)})

    H, k = local_search ? (p.H, p.k) : (0, 1)
    d_viol, d_sol = local_search ? (p.d_viol, p.d_sol) : (0, 0)

    tag = "$(inst.name)_$(cfgid)_$(backend_name(backend))_s$(seed)"
    dir = joinpath(outdir, tag)
    mkpath(dir)

    rec = Dict{String,Any}(
        "instance" => inst.name, "config" => cfgid, "disabled" => opts.note,
        "backend" => backend_name(backend), "solver_version" => solver_version(backend),
        "seed" => seed, "H" => H, "k" => k, "d_viol" => d_viol, "d_sol" => d_sol,
        "cap_seconds" => cap, "julia_version" => string(VERSION),
        "started_at" => string(now()), "status" => "running",
    )

    try
        r = solve_otsd(inst.case;
            tlf = inst.tlf, H, d_viol, d_sol, k,
            warm_start = inst.warm_start,
            embed = embed ? inst.embedded : nothing,
            inherit, backend, seed, timeout = Float64(cap),
            logdir = joinpath("ablation", tag), label = tag,
            solver_kwargs...)

        rec["status"] = isfinite(r.objective) ? "done" : "no solution"
        rec["objective"] = isfinite(r.objective) ? r.objective : nothing
        rec["secure"] = r.secure
        rec["wall_seconds"] = r.wall
        rec["lp_solves"] = r.lp_solves
        rec["benders_iterations"] = r.benders_iterations
        rec["openings"] = ["$(a)-$(b)" for (a, b) in sort(collect(r.openings))]
        rec["phases"] = [Dict("iteration" => it.iteration, "objective" => it.obj,
                              "secure" => it.secure, "sbs_size" => it.sbs_size,
                              "seconds" => it.time, "capped" => it.capped)
                         for it in r.phases]
    catch e
        rec["status"] = "failed"
        rec["error"] = sprint(showerror, e)
        @error "✗ $tag" exception = (e, catch_backtrace())
    end

    rec["finished_at"] = string(now())
    open(io -> JSON3.pretty(io, rec), joinpath(dir, "result.json"), "w")
    rec
end

"""
    run_campaign(cells; outdir) -> Vector{Dict}

Run `cells`, each a `(instance, config, backend, seed, cap)` tuple, one at a time.
A cell whose `result.json` already exists is skipped, so an interrupted campaign
resumes where it stopped.
"""
function run_campaign(cells; outdir::String = "results")
    mkpath(outdir)
    records = Dict{String,Any}[]
    t0 = now()
    for (i, c) in enumerate(cells)
        tag = "$(c.instance.name)_$(c.config)_$(backend_name(c.backend))_s$(c.seed)"
        path = joinpath(outdir, tag, "result.json")
        if isfile(path)
            @info "[$i/$(length(cells))] ↻ $tag (on disk)"
            push!(records, copy(JSON3.read(read(path, String), Dict{String,Any})))
            continue
        end
        @info "[$i/$(length(cells))] ▶ $tag  cap=$(c.cap)s  elapsed=$(canonicalize(now() - t0))"
        push!(records, run_cell(c.instance, c.config;
                                backend = c.backend, seed = c.seed, cap = c.cap, outdir))
        summarise(records)
    end
    records
end

"Print what is on disk so far. Safe to call at any point."
function summarise(records)
    println("\n", "="^94)
    @printf("%-5s %-15s %-8s %-5s %-9s %-8s %-10s %-7s %s\n",
            "sys", "config", "backend", "seed", "objective", "secure", "wall (s)", "LPs", "status")
    println("-"^94)
    for r in records
        @printf("%-5s %-15s %-8s %-5s %-9s %-8s %-10s %-7s %s\n",
                r["instance"], r["config"], r["backend"], r["seed"],
                something(get(r, "objective", nothing), "—"),
                something(get(r, "secure", nothing), "—"),
                something(get(r, "wall_seconds", nothing), "—"),
                something(get(r, "lp_solves", nothing), "—"),
                r["status"])
    end
    println("="^94, "\n")
end

"A cell specification."
cell(inst, config; backend, seed = 0, cap = 900) =
    (; instance = inst, config, backend, seed, cap)
