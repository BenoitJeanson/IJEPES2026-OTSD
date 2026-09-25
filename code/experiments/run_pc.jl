# ── p_c sensitivity ───────────────────────────────────────────────────────────
#
# The paper calls `p_c` the probability of occurrence of contingency `c`, but every
# run behind the reported numbers used the uniform value 1 — so the objective is a
# *total* over the single-branch contingencies, not an expectation. This study asks
# the one question the paper can honestly answer: does a non-uniform `p_c` change the
# *outcome*?
#
# It is deliberately not a statement about techno-economic value. Pricing shed load
# against redispatch would need cost data this work does not have; the paper's claim
# is that an algorithm exists for the problem, not that its solutions are economically
# optimal.
#
# The weighting. In practice `p_c` comes from outage statistics and asset models. A
# reproducible stand-in is failure rate proportional to line length, and the DC model's
# only length-like quantity is the series reactance `x = 1/b`:
#
#     w_c = x_c                       for a line
#     w_c = mean(x over lines) / 5    for a transformer
#
# A transformer's reactance is winding impedance, not a length, so feeding it through
# the same rule is meaningless — on IEEE-57 the largest transformer reactance (1.355)
# exceeds every line, which would make it the *most* likely failure. Transformers
# therefore take a flat value, five times less likely than a typical line.
#
# Weights are normalised to mean 1 over the contingency set, so the uniform arm is
# exactly `w ≡ 1` and both arms' objectives stay on the scale the paper reports.
# Normalising to sum 1 would also be defensible but rescales by 1/|C| and makes
# nothing comparable to the published numbers.
#
# What is measured is the topology, not the clock. For each arm's final topology both
# objectives are evaluated — the weighted one the solver minimised, and the unweighted
# total lost load — with `create_bridge_to_pocket`, which is deterministic and needs no
# solver. Each arm should win its own metric; the question is whether the topologies
# differ at all, and what each costs under the other's measure.
#
#   julia --project=. experiments/run_pc.jl            # run, then report
#   julia --project=. experiments/run_pc.jl --report   # re-render from disk only

include(joinpath(@__DIR__, "ablation.jl"))

using PGLib

const PC_SEEDS = 0:2
const PC_CAP = 900
const PC_OUT = joinpath(@__DIR__, "..", "results")

# ── The weighting ─────────────────────────────────────────────────────────────

"""
    is_transformer(case) -> Dict{ELabel,Bool}

Which branches of the graph are transformers. The graph keeps only `b = 1/x` and
`p_max`, so the tap ratio has to be read back from the PGLib case. Parallel circuits
are collapsed into one branch, and a collapsed pair counts as a transformer if either
circuit is one.
"""
function is_transformer(case::String)
    c = pglib(case)
    xf = Dict{ELabel,Bool}()
    for br in values(c["branch"])
        f, t = string(br["f_bus"]), string(br["t_bus"])
        k = parse(Int, f) < parse(Int, t) ? (f, t) : (t, f)
        istr = get(br, "transformer", false) || abs(get(br, "tap", 1.0) - 1.0) > 1e-9
        xf[k] = get(xf, k, false) || istr
    end
    xf
end

"""
    pc_weights(rc, case) -> Dict{ELabel,Float64}

`p_c` proportional to branch length, normalised to mean 1 over the contingency set.
Length is proxied by the series reactance `x = 1/b`; transformers take a flat
`mean(x over lines) / 5`. See the header for why.
"""
function pc_weights(rc::RichCase, case::String)
    g = rc.gc.g
    xf = is_transformer(case)
    els = collect(ELabel, edge_labels(g))
    x = Dict(e => 1 / g[e...].b for e in els)
    lines = [e for e in els if !get(xf, e, false)]
    isempty(lines) && error("no lines in $case — the transformer rule has no anchor")
    x_tr = sum(x[e] for e in lines) / length(lines) / 5
    raw = Dict(e => (get(xf, e, false) ? x_tr : x[e]) for e in els)
    scale = length(els) / sum(values(raw))       # mean 1, so the uniform arm is w ≡ 1
    Dict(e => v * scale for (e, v) in raw)
end

"Total lost load of a topology, summed over every single-branch contingency."
risk_unweighted(rc::RichCase, openings) =
    sum(pk.d for pk in values(create_bridge_to_pocket(rc.gc, Set{ELabel}(openings))); init = 0.0)

"The same, with each contingency weighted by its `p_c`."
risk_weighted(rc::RichCase, openings, w::Dict{ELabel,Float64}) =
    sum(get(w, br, 1.0) * pk.d
        for (br, pk) in create_bridge_to_pocket(rc.gc, Set{ELabel}(openings)); init = 0.0)

# ── The campaign ──────────────────────────────────────────────────────────────
#
# Two arms at the paper's operating point, differing only in `ov_lostload_coef`: the
# uniform arm leaves it at its default of 1.0, the weighted arm passes the Dict. The
# weights are per system, so they travel on the cell rather than in `CONFIGS`.

"The two arms for one instance, as `(config id, options)` pairs."
function pc_arms(inst::Instance)
    rc = load_case(inst.case; tlf = inst.tlf)
    w = pc_weights(rc, inst.case)
    ws = sort(collect(values(w)))
    @info "  $(inst.name): |C| = $(length(w))  w ∈ [$(round(ws[1]; digits=3)), " *
          "$(round(ws[end]; digits=3))]  median $(round(ws[end÷2]; digits=3))"
    [("PC-UNIF", (; note = "p_c uniform (= 1), as the paper reports")),
     ("PC-LEN", (; ov_lostload_coef = w,
                   note = "p_c ∝ branch length (x = 1/b); transformers flat at mean(x_line)/5"))]
end

"""
    pc_cells(; backend, cap, seeds, systems)

The p_c arms as campaign cells, so the study runs inside `run_ablation.jl` and shows
up in the live view alongside the ablation rather than as a separate campaign.
"""
function pc_cells(; backend = GurobiBackend(threads = 4), cap = PC_CAP,
                    seeds = PC_SEEDS, systems = (IEEE57, IEEE118))
    cells = []
    for inst in systems, (cfgid, opts) in pc_arms(inst), seed in seeds
        push!(cells, cell(inst, cfgid; backend, seed, cap, opts))
    end
    cells
end

# ── The comparison ────────────────────────────────────────────────────────────

parse_opening(s::AbstractString) = (a = split(s, "-"); (String(a[1]), String(a[2])))

"One row per (system, arm, seed): the topology reached, and what it costs either way."
function pc_rows(; outdir = PC_OUT, backend = "gurobi")
    rows = []
    for inst in (IEEE57, IEEE118)
        rc = load_case(inst.case; tlf = inst.tlf)
        w = pc_weights(rc, inst.case)
        for arm in ("PC-UNIF", "PC-LEN"), seed in PC_SEEDS
            path = joinpath(outdir, "$(inst.name)_$(arm)_$(backend)_s$(seed)", "result.json")
            isfile(path) || continue
            r = JSON3.read(read(path, String))
            op = Set(parse_opening(o) for o in r.openings)
            push!(rows, (; system = inst.name, arm, seed, status = r.status,
                           openings = op, n_open = length(op),
                           unweighted = risk_unweighted(rc, op),
                           weighted = risk_weighted(rc, op, w),
                           solver_obj = r.objective, wall = r.wall_seconds))
        end
    end
    rows
end

"Branches open in one topology and not the other."
hamming_sets(a, b) = length(symdiff(a, b))

function write_pc_report(rows; path = joinpath(PC_OUT, "PC.md"))
    io = IOBuffer()
    println(io, "# p_c sensitivity — outcome, not performance\n")
    println(io, "Lost load in per unit on a 100 MVA base. `weighted` applies p_c ∝ branch")
    println(io, "length (x = 1/b), transformers flat at mean(x_line)/5, normalised to mean 1 —")
    println(io, "so the uniform arm is exactly p_c ≡ 1 and both columns are on the scale the")
    println(io, "paper reports. Each arm minimises its own column; the question is whether the")
    println(io, "topologies differ at all, and what each costs under the other's measure.\n")

    for sys in unique(r.system for r in rows)
        srows = [r for r in rows if r.system == sys]
        isempty(srows) && continue
        println(io, "## IEEE-$sys\n")
        println(io, "| arm | seed | openings | unweighted | weighted | solver obj | wall (s) |")
        println(io, "|---|---:|---:|---:|---:|---:|---:|")
        for r in sort(srows, by = r -> (r.arm, r.seed))
            println(io, "| ", r.arm, " | ", r.seed, " | ", r.n_open,
                    " | ", round(r.unweighted; digits = 3),
                    " | ", round(r.weighted; digits = 3),
                    " | ", r.solver_obj, " | ", r.wall, " |")
        end

        unif = [r for r in srows if r.arm == "PC-UNIF"]
        len = [r for r in srows if r.arm == "PC-LEN"]
        if !isempty(unif) && !isempty(len)
            d = [hamming_sets(u.openings, l.openings)
                 for u in unif for l in len if u.seed == l.seed]
            println(io, "\nPer-seed topology difference (branches open in one arm only): ",
                    isempty(d) ? "—" : join(d, ", "))
            bu = minimum(r.unweighted for r in unif); bl = minimum(r.unweighted for r in len)
            wu = minimum(r.weighted for r in unif); wl = minimum(r.weighted for r in len)
            println(io, "\nBest of each arm, under both measures:\n")
            println(io, "| | unweighted | weighted |")
            println(io, "|---|---:|---:|")
            println(io, "| PC-UNIF | ", round(bu; digits = 3), " | ", round(wu; digits = 3), " |")
            println(io, "| PC-LEN  | ", round(bl; digits = 3), " | ", round(wl; digits = 3), " |")
            println(io, "\nEach arm should win its own column. ",
                    bu ≤ bl && wl ≤ wu ? "It does here." :
                    "It does not here, which means the outer loop's heuristic search, " *
                    "not the weighting, decided the outcome — see NO-LOCALSEARCH in the ablation.")
        end
        println(io)
    end
    open(f -> write(f, String(take!(io))), path, "w")
    path
end

# ── Entry point ───────────────────────────────────────────────────────────────

if abspath(PROGRAM_FILE) == @__FILE__
    if "--report" ∉ ARGS
        cells = pc_cells()
        @info "p_c sensitivity: $(length(cells)) runs (2 systems × 2 arms × $(length(PC_SEEDS)) seeds)"
        run_campaign(cells; outdir = PC_OUT)
    end
    println(write_pc_report(pc_rows()))
end
