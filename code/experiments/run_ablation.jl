# Overnight campaign: the component ablation on both systems and both backends.
#
#   julia --project=. experiments/run_ablation.jl            # everything below
#   julia --project=. experiments/run_ablation.jl gurobi     # one backend only
#   julia --project=. experiments/run_ablation.jl scip       # the licence-free block
#
# Cells already on disk are skipped, so this is safe to re-run after an interruption.

# `run_pc.jl` pulls in `ablation.jl`; its own driver is guarded, so including it
# here brings the p_c arms into the campaign without running them twice.
include(joinpath(@__DIR__, "run_pc.jl"))

const WHICH = isempty(ARGS) ? "all" : ARGS[1]
const OUT = joinpath(@__DIR__, "..", "results")

# No cell may run longer than an hour. Nothing observed comes close — the slowest SCIP
# cell measured is IEEE-118 NO-SCREEN at 980 s — so a cell that reaches this cap is
# telling us something is wrong, not that it needed the time.
const CAP_SCIP = 3000

# Gurobi: every component, both systems, three seeds. This is the claim the paper
# makes, regenerated from the released code.
const GUROBI = GurobiBackend(threads = 4)

gurobi_cells = vcat(
    vec([cell(INSTANCES[sys], cfg; backend = GUROBI, seed, cap = 900)
         for sys in ("57", "118"), cfg in CONFIG_ORDER, seed in 0:2]),
    # The p_c sensitivity arms: the same operating point, weighted and unweighted.
    # The uniform arm reproduces REF, which is the cheapest check that the weighting
    # plumbing left the baseline alone.
    pc_cells(backend = GUROBI, cap = 900),
)

# SCIP: the same algorithm, run without a licence. Because the cuts enter at the same
# places, this block is a reproduction claim and not merely a portability one — the
# objective must match the Gurobi cell above it, phase by phase. The final openings
# need not: the local search warm-starts each phase from the previous one, so a tie
# broken differently in phase 1 legitimately lands phase 2 elsewhere.
#
# NO-LOCALSEARCH is absent on purpose. It never converges — it is the evidence that
# the outer loop is heuristic — so its cells burn the cap outright, and the claim is
# about the search, not the solver. Six capped cells would cost six hours here and
# show nothing the Gurobi campaign does not already show.
#
# Cells are ordered by expected cost, cheapest first. SCIP is 4–15x slower than
# Gurobi and the two IEEE-118 cells at the end may take an hour each, so a campaign
# stopped at a wall-clock budget loses only its most expensive cells — and every cell
# that finished is on disk, because a campaign resumes from what it finds there.
const SCIP_COST_ORDER = [        # measured on IEEE-57, extrapolated on IEEE-118
    ("57", "NO-EMBED"), ("57", "NO-INHERIT"), ("57", "REF"),
    ("57", "NO-CF-POCKET"), ("57", "NO-CF-FREE"), ("57", "NO-SCREEN"),
    ("118", "REF"), ("118", "NO-INHERIT"), ("118", "NO-CF-FREE"),
    ("118", "NO-CF-POCKET"), ("118", "NO-EMBED"), ("118", "NO-SCREEN"),
]

scip_cells = vcat(
    [cell(INSTANCES[sys], cfg; backend = SCIPBackend(), seed, cap = CAP_SCIP)
     for (sys, cfg) in SCIP_COST_ORDER[1:6] for seed in 0:2],      # IEEE-57 ablation
    pc_cells(backend = SCIPBackend(), cap = CAP_SCIP, systems = (IEEE57,)),
    [cell(INSTANCES[sys], cfg; backend = SCIPBackend(), seed, cap = CAP_SCIP)
     for (sys, cfg) in SCIP_COST_ORDER[7:10] for seed in 0:2],     # IEEE-118, the cheap four
    pc_cells(backend = SCIPBackend(), cap = CAP_SCIP, systems = (IEEE118,)),
    [cell(INSTANCES[sys], cfg; backend = SCIPBackend(), seed, cap = CAP_SCIP)
     for (sys, cfg) in SCIP_COST_ORDER[11:12] for seed in 0:2],    # the two expensive ones, last
)

cells = WHICH == "gurobi" ? vec(gurobi_cells) :
        WHICH == "scip"   ? scip_cells :
                            vcat(vec(gurobi_cells), scip_cells)

@info "Campaign: $(length(cells)) cells, results in $OUT"
records = run_campaign(cells; outdir = OUT)
summarise(records)

# Exit nonzero when any cell failed, so a stage that dies cannot be mistaken for one
# that finished. A campaign script that ignores this is how six missing cells came to
# be reported as "all done".
nfailed = failed_count(records)
if nfailed > 0
    @error "$nfailed of $(length(records)) cells did not finish"
    exit(1)
end
