# Overnight campaign: the component ablation on both systems and both backends.
#
#   julia --project=. experiments/run_ablation.jl            # everything below
#   julia --project=. experiments/run_ablation.jl gurobi     # one backend only
#
# Cells already on disk are skipped, so this is safe to re-run after an interruption.

include(joinpath(@__DIR__, "ablation.jl"))

const WHICH = isempty(ARGS) ? "all" : ARGS[1]
const OUT = joinpath(@__DIR__, "..", "results")

# Gurobi: every component, both systems, three seeds. This is the claim the paper
# makes, regenerated from the released code.
gurobi_cells = [
    cell(INSTANCES[sys], cfg; backend = GurobiBackend(threads = 4), seed, cap = 900)
    for sys in ("57", "118"), cfg in CONFIG_ORDER, seed in 0:2
]

# HiGHS: the same on IEEE-57, and on IEEE-118 only the configurations that are cheap
# under Gurobi. Re-solving the master per round costs more than the callback, and a
# licence-free path is a portability claim, not a performance one.
const HIGHS_118 = ["REF", "NO-CF-POCKET", "NO-CF-FREE", "NO-INHERIT"]
highs_cells = vcat(
    [cell(IEEE57, cfg; backend = HiGHSBackend(threads = 4), seed = 0, cap = 900)
     for cfg in CONFIG_ORDER],
    [cell(IEEE118, cfg; backend = HiGHSBackend(threads = 4), seed = 0, cap = 900)
     for cfg in HIGHS_118],
)

cells = WHICH == "gurobi" ? vec(gurobi_cells) :
        WHICH == "highs"  ? highs_cells :
                            vcat(vec(gurobi_cells), highs_cells)

@info "Campaign: $(length(cells)) cells, results in $OUT"
records = run_campaign(cells; outdir = OUT)
summarise(records)
