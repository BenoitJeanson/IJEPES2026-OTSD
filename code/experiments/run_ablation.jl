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

# HiGHS: a portability claim, not a performance one. Re-solving the master once per
# cut round costs far more than the callback, so this is deliberately a handful of
# cells with a generous cap — enough to show the package runs without a commercial
# licence and returns a secure topology, not enough to compare runtimes.
highs_cells = vcat(
    [cell(IEEE57, cfg; backend = HiGHSBackend(threads = 4), seed = 0, cap = 2400)
     for cfg in ("REF", "NO-CF-POCKET", "NO-INHERIT")],
    [cell(IEEE118, "REF"; backend = HiGHSBackend(threads = 4), seed = 0, cap = 2400)],
)

cells = WHICH == "gurobi" ? vec(gurobi_cells) :
        WHICH == "highs"  ? highs_cells :
                            vcat(vec(gurobi_cells), highs_cells)

@info "Campaign: $(length(cells)) cells, results in $OUT"
records = run_campaign(cells; outdir = OUT)
summarise(records)
