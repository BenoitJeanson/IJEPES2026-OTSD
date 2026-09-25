using OTSD
using Test
using MetaGraphsNext: labels, edge_labels

include(joinpath(@__DIR__, "..", "experiments", "instances.jl"))

@testset "OTSD" begin

    @testset "case loading" begin
        # Parallel circuits collapse: a branch is keyed on its pair of buses.
        # These counts are what every reported figure was computed on.
        rc118 = load_case("case118"; tlf = 1.5)
        rc57 = load_case("case57"; tlf = 1.0)
        @test length(collect(edge_labels(rc118.gc.g))) == 179
        @test length(collect(edge_labels(rc57.gc.g))) == 78

        # The reference bus is the largest net injection, which is what the
        # published runs used.
        @test rc118.gc.bus_orig == "69"
        @test rc57.gc.bus_orig == "8"

        # Balancing leaves no net imbalance.
        @test sum(rc118.gc.g[b] for b in labels(rc118.gc.g)) ≈ 0 atol = 1e-6
    end

    @testset "thermal limit factor scales limits, not injections" begin
        plain = load_case("case57")
        scaled = load_case("case57"; tlf = 2.0)
        br = first(edge_labels(plain.gc.g))
        @test scaled.gc.g[br...].p_max ≈ 2 * plain.gc.g[br...].p_max
        @test scaled.gc.g[br[1]] ≈ plain.gc.g[br[1]]
    end

    @testset "the warm start is not already optimal" begin
        # Guards against a warm start that silently makes the search a no-op.
        for inst in (IEEE57, IEEE118)
            rc = load_case(inst.case; tlf = inst.tlf)
            @test !isempty(inst.warm_start)
            @test all(br -> br in Set(edge_labels(rc.gc.g)), inst.warm_start)
        end
    end

    # The acceptance test. These two numbers are what the paper reports at
    # PAPER_OPERATING_POINT; if either moves, this package no longer reproduces it.
    @testset "reproduces the published objective" begin
        if !isassigned(OTSD.GRB_ENV_REF)
            @info "No Gurobi licence — skipping the reproduction test."
        else
            p = PAPER_OPERATING_POINT
            for inst in (IEEE57, IEEE118)
                r = solve_otsd(inst.case;
                    tlf = inst.tlf, H = p.H, d_viol = p.d_viol, d_sol = p.d_sol, k = p.k,
                    warm_start = inst.warm_start, embed = inst.embedded, seed = p.seed,
                    backend = GurobiBackend(threads = 4),
                    logdir = "test", label = "TEST_$(inst.name)")
                @test r.objective ≈ reference_objective(inst) atol = 1e-3
                @test r.secure
            end
        end
    end

    @testset "backends" begin
        # Both master backends take a cut inside the tree; that is what makes them
        # the same algorithm. HiGHS is here as the subproblem LP solver only, and
        # the master must refuse it rather than silently solving without cuts.
        @test supports_lazy(GurobiBackend())
        @test supports_lazy(SCIPBackend())
        @test !supports_lazy(HiGHSLPBackend())
        @test backend_name(SCIPBackend()) == "scip"
        @test !isempty(solver_version(SCIPBackend()))

        # A SCIP master pairs with HiGHS subproblems: SCIP exposes no Farkas dual.
        @test OTSD.lp_backend(SCIPBackend()) isa HiGHSLPBackend
        @test OTSD.lp_backend(GurobiBackend()) isa GurobiBackend
    end
end
