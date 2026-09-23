# Render campaign progress to a markdown file, from whatever is on disk.
#
#   julia --project=. experiments/live_report.jl           # once
#   julia --project=. experiments/live_report.jl --watch    # every 20 s
#
# Nothing here talks to the running campaign: a finished run is a `result.json`, and
# the run in flight is whichever directory has none yet. So this is safe to run at any
# time, and it cannot disturb the timings it reports.

using JSON3, Dates, Printf

const RESULTS = joinpath(@__DIR__, "..", "results")
const OUT = joinpath(RESULTS, "LIVE.md")
const ORDER = ["REF", "NO-SCREEN", "NO-CF-POCKET", "NO-CF-FREE",
               "NO-EMBED", "NO-INHERIT", "NO-LOCALSEARCH"]

records() = [JSON3.read(read(f, String)) for f in
             sort(filter(isfile, [joinpath(RESULTS, d, "result.json")
                                  for d in readdir(RESULTS) if isdir(joinpath(RESULTS, d))]))]

"Directories with no result.json yet — at most one is genuinely in flight."
function in_flight()
    dirs = [d for d in readdir(RESULTS)
            if isdir(joinpath(RESULTS, d)) && !isfile(joinpath(RESULTS, d, "result.json"))]
    isempty(dirs) ? nothing : last(sort(dirs, by = d -> mtime(joinpath(RESULTS, d))))
end

fmt(x) = x === nothing ? "—" : x isa AbstractFloat ? string(round(x; digits = 3)) : string(x)

function table(io, rs, backend)
    rows = [r for r in rs if r.backend == backend]
    isempty(rows) && return
    println(io, "\n## $(uppercase(backend)) — $(length(rows)) runs\n")
    println(io, "| system | config | seed | objective | secure | wall (s) | LP solves | Benders | status |")
    println(io, "|---|---|---:|---:|:---:|---:|---:|---:|---|")
    key(r) = (r.instance == "118" ? 0 : 1,
              something(findfirst(==(r.config), ORDER), 99), r.seed)
    for r in sort(rows, by = key)
        capped = get(r, :wall_seconds, 0.0) ≥ get(r, :cap_seconds, Inf) - 1
        println(io, "| ", r.instance, " | ", r.config, " | ", r.seed,
                " | **", fmt(get(r, :objective, nothing)), "**",
                " | ", get(r, :secure, nothing) === true ? "✓" : "✗",
                " | ", fmt(get(r, :wall_seconds, nothing)), capped ? " ⏱" : "",
                " | ", fmt(get(r, :lp_solves, nothing)),
                " | ", fmt(get(r, :benders_iterations, nothing)),
                " | ", r.status, " |")
    end
end

function render()
    rs = records()
    open(OUT, "w") do io
        println(io, "# Campaign progress\n")
        println(io, "_", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"), " — ",
                length(rs), " runs complete._\n")
        running = in_flight()
        if running !== nothing
            println(io, "**In flight:** `", running, "`\n")
        end
        for b in ("gurobi", "highs")
            table(io, rs, b)
        end
        println(io, "\n⏱ marks a run stopped by its time cap: its objective is the best found, not a converged optimum.")
    end
    OUT
end

if "--watch" in ARGS
    while true
        render()
        sleep(20)
    end
else
    println(render())
end
