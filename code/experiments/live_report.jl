# Render campaign progress to a markdown file, from whatever is on disk.
#
#   julia --project=. experiments/live_report.jl           # once
#   julia --project=. experiments/live_report.jl --watch    # every 20 s
#
# Nothing here talks to the running campaign: a finished run is a `result.json`, and
# the run in flight is whichever directory has none yet. So this is safe to run at any
# time, and it cannot disturb the timings it reports. It deliberately does not load
# OTSD — it needs no solver, and a second Gurobi environment would be one more thing
# competing with the campaign it is watching.
#
# The reference campaign in `docs/data` is read, never written. It is the data the
# paper reports; a new run is compared against it, and any decision to change it is
# the user's.

using JSON3, Dates, Printf

const RESULTS = joinpath(@__DIR__, "..", "results")
const LOGS = joinpath(@__DIR__, "..", "logs")
const OUT = joinpath(RESULTS, "LIVE.md")
const REF_INDEX = joinpath(@__DIR__, "..", "..", "docs", "data", "index.json")

const ORDER = ["REF", "NO-SCREEN", "NO-CF-POCKET", "NO-CF-FREE",
               "NO-EMBED", "NO-INHERIT", "NO-LOCALSEARCH", "PC-UNIF", "PC-LEN"]
# The p_c arms run inside the campaign, so they belong in the matrix; they have no
# counterpart in the reference campaign, which only ever ran p_c ≡ 1.
const PC_ARMS = ["PC-UNIF", "PC-LEN"]
const SYSTEMS = ("57", "118")
const SEEDS = 0:2

# What each campaign is supposed to contain, so the matrix can show what has not run
# yet. The SCIP block omits NO-LOCALSEARCH: it never converges, so its cells burn the
# cap outright, and the claim is about the search rather than the solver.
const EXPECTED = Dict("gurobi" => ORDER,
                      "scip" => filter(!=("NO-LOCALSEARCH"), ORDER))

# ── What is on disk ───────────────────────────────────────────────────────────

celldir(tag) = joinpath(RESULTS, tag)
tag_of(sys, cfg, backend, seed) = "$(sys)_$(cfg)_$(backend)_s$(seed)"

"""
Every finished cell on disk. A `result.json` caught mid-write parses as truncated
JSON, so an unreadable file is skipped rather than fatal: the campaign is writing it
right now and the next tick will read it whole.
"""
function records()
    out = []
    for d in sort(readdir(RESULTS))
        f = joinpath(RESULTS, d, "result.json")
        (isdir(joinpath(RESULTS, d)) && isfile(f)) || continue
        try
            push!(out, JSON3.read(read(f, String)))
        catch
            continue
        end
    end
    out
end

"Directories with no result.json yet — at most one is genuinely in flight."
function in_flight()
    dirs = [d for d in readdir(RESULTS)
            if isdir(joinpath(RESULTS, d)) && !isfile(joinpath(RESULTS, d, "result.json"))]
    isempty(dirs) ? nothing : last(sort(dirs, by = d -> mtime(celldir(d))))
end

"""
Newest campaign log for `backend`, so the reader knows where the detail is. Sorted by
mtime rather than by name: logs are named `<timestamp>_<label>.log`, and an older log
without the prefix would otherwise sort last and win.
"""
function newest_log(backend)
    isdir(LOGS) || return nothing
    ls = filter(f -> occursin(backend, f) && endswith(f, ".log"), readdir(LOGS))
    isempty(ls) ? nothing : last(sort(ls, by = f -> mtime(joinpath(LOGS, f))))
end

capped(r) = get(r, :wall_seconds, 0.0) ≥ get(r, :cap_seconds, Inf) - 1

# ── The reference campaign ────────────────────────────────────────────────────

"""
Per `(system, config)`, what the recorded campaign in `docs/data` found at the paper
operating point across its seeds. Read-only: this is the data the paper reports.
"""
function reference()
    isfile(REF_INDEX) || return Dict()
    out = Dict{Tuple{String,String},Vector{Any}}()
    for r in JSON3.read(read(REF_INDEX, String)).runs
        (r.H == 3 && r.d_viol == 2 && r.d_sol == 0) || continue
        sys = r.system == "ieee118" ? "118" : "57"
        push!(get!(out, (sys, r.config), []), r)
    end
    out
end

"""
How a finished cell compares with the reference: the objective must match, and the
LP count is expected to land inside the spread the reference campaign shows across
its seeds. Outside that spread is not a failure — it is the signal worth looking at.
"""
function vs_reference(r, ref)
    rs = get(ref, (string(r.instance), string(r.config)), nothing)
    rs === nothing && return "—"
    obj = get(r, :objective, nothing)
    obj === nothing && return "no objective"
    objok = any(x -> abs(x.obj - obj) < 1e-3, rs)
    lps = [x.lp for x in rs]
    lp = get(r, :lp_solves, nothing)
    inrange = lp !== nothing && minimum(lps) ≤ lp ≤ maximum(lps)
    string(objok ? "obj ✓" : "**obj ✗ $(rs[1].obj)**", " · LP ",
           inrange ? "in" : "**out**", " $(minimum(lps))–$(maximum(lps))")
end

# ── Rendering ─────────────────────────────────────────────────────────────────

fmt(x) = x === nothing ? "—" : x isa AbstractFloat ? string(round(x; digits = 3)) : string(x)

"""
The run matrix: every cell the campaign is meant to produce, and where each one is.
Pending cells are the point — a table of finished runs cannot show what is left.
"""
function matrix(io, backend, done, running)
    cfgs = get(EXPECTED, backend, ORDER)
    println(io, "\n### Matrix — ", uppercase(backend), "\n")
    print(io, "| config |")
    for sys in SYSTEMS, s in SEEDS
        print(io, " $sys·s$s |")
    end
    println(io, "\n|---|", repeat(":-:|", length(SYSTEMS) * length(SEEDS)))
    for cfg in cfgs
        print(io, "| ", cfg, " |")
        for sys in SYSTEMS, s in SEEDS
            t = tag_of(sys, cfg, backend, s)
            mark = if haskey(done, t)
                r = done[t]
                r.status != "done" ? "✗" : capped(r) ? "⏱" : "✓"
            elseif t == running
                "▶"
            else
                "·"
            end
            print(io, " ", mark, " |")
        end
        println(io)
    end
    n = length(cfgs) * length(SYSTEMS) * length(SEEDS)
    pending = [(sys, cfg) for cfg in cfgs, sys in SYSTEMS, s in SEEDS
               if !haskey(done, tag_of(sys, cfg, backend, s))]
    println(io, "\n`✓` done `⏱` hit its cap `▶` running `·` not started `✗` failed — ",
            "**$(n - length(pending))/$n**")
    pending
end

function table(io, rs, backend, ref)
    rows = [r for r in rs if r.backend == backend]
    isempty(rows) && return
    println(io, "\n### Results — ", uppercase(backend), "\n")
    println(io, "| system | config | seed | objective | secure | wall (s) | LP solves | Benders | vs reference |")
    println(io, "|---|---|---:|---:|:---:|---:|---:|---:|---|")
    key(r) = (r.instance == "118" ? 0 : 1,
              something(findfirst(==(r.config), ORDER), 99), r.seed)
    for r in sort(rows, by = key)
        println(io, "| ", r.instance, " | ", r.config, " | ", r.seed,
                " | **", fmt(get(r, :objective, nothing)), "**",
                " | ", get(r, :secure, nothing) === true ? "✓" : "✗",
                " | ", fmt(get(r, :wall_seconds, nothing)), capped(r) ? " ⏱" : "",
                " | ", fmt(get(r, :lp_solves, nothing)),
                " | ", fmt(get(r, :benders_iterations, nothing)),
                " | ", vs_reference(r, ref), " |")
    end
end

"""
    eta(rows, pending) -> String

Estimated time left. Each pending cell is costed from finished cells of the same
`(system, config)` when there are any — seed 0 of every configuration runs first, so
usually there are — and from the campaign mean otherwise. Costing everything at the
campaign mean would badly under-estimate: NO-LOCALSEARCH never converges and always
spends its whole cap, while a REF cell is seconds.
"""
function eta(rows, pending)
    isempty(pending) && return "done"
    walls = Dict{Tuple{String,String},Vector{Float64}}()
    for r in rows
        w = get(r, :wall_seconds, nothing)
        w === nothing && continue
        push!(get!(walls, (string(r.instance), string(r.config)), Float64[]), Float64(w))
    end
    isempty(walls) && return "unknown"
    allw = reduce(vcat, values(walls))
    mean(v) = sum(v) / length(v)
    secs = sum(haskey(walls, k) ? mean(walls[k]) : mean(allw) for k in pending)
    string("~", Int(secs ÷ 3600), "h", lpad(Int(secs % 3600 ÷ 60), 2, '0'))
end

function render()
    rs = records()
    ref = reference()
    running = in_flight()
    done = Dict(string(r.instance, "_", r.config, "_", r.backend, "_s", r.seed) => r for r in rs)

    open(OUT, "w") do io
        println(io, "# Campaign progress\n")
        println(io, "_", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"), "_\n")
        if running !== nothing
            age = round(Int, (time() - mtime(celldir(running))) / 60)
            println(io, "**In flight:** `", running, "` — started ", age, " min ago\n")
        else
            println(io, "**In flight:** nothing\n")
        end

        for b in ("gurobi", "scip")
            any(r -> r.backend == b, rs) || haskey(EXPECTED, b) || continue
            pending = matrix(io, b, done, running)
            rows = [r for r in rs if r.backend == b]
            lg = newest_log(b)
            println(io, "\nremaining ", eta(rows, pending),
                    lg === nothing ? "" : " · log `logs/$lg`")
            table(io, rows, b, ref)
        end

        println(io, "\n---\n")
        println(io, "`vs reference` compares against `docs/data`, the recorded campaign the paper ",
                    "reports: the objective must match, and the LP count is shown against the ",
                    "spread that campaign found across its seeds. Landing outside it is a signal ",
                    "to look, not a failure. **`docs/data` is never written by a campaign.**")
        println(io, "\n⏱ marks a run stopped by its time cap: its objective is the best found, ",
                    "not a converged optimum.")
    end
    OUT
end

if "--watch" in ARGS
    # A watcher that dies takes the only view of a twelve-hour campaign with it, so no
    # single failed render is allowed to end the loop.
    while true
        try
            render()
        catch e
            @warn "render failed; retrying next tick" exception = e
        end
        sleep(20)
    end
else
    println(render())
end
