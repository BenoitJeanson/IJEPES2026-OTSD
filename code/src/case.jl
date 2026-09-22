# ── Loading a benchmark case ──────────────────────────────────────────────────
#
# Both published experiments run on PGLib's IEEE-57 and IEEE-118. The graph keys a
# branch on its pair of buses, so where PGLib lists parallel circuits on a corridor
# the corridor becomes one branch: 80 → 78 on IEEE-57, 186 → 179 on IEEE-118. This is
# deliberate and is discussed in the paper; it can only increase the number of
# branches whose outage severs a pocket, never reduce it.

"""
    load_case(case; tlf = nothing, reference = nothing) -> RichCase

Build the case the experiments run on.

  * `case`      — a PGLib name, e.g. `"case118"`.
  * `tlf`       — thermal limit factor: scales every branch limit uniformly. The paper
                  uses 1.5 for IEEE-118 and 1.0 for IEEE-57, both stressed operating
                  points. `nothing` leaves the limits as PGLib sets them.
  * `reference` — the reference bus. Buses it cannot reach over closed branches are
                  de-energized. Defaults to the largest net injection, which is bus
                  `"69"` on IEEE-118 — the value the published runs used.

Injections are balanced before the limits are scaled, so the two are independent.
"""
function load_case(case::String; tlf::Union{Nothing,Real} = nothing,
                   reference::Union{Nothing,VLabel} = nothing)
    g = PGLibtograph(case)
    balance!(g)
    isnothing(tlf) || scale_branch_limits!(g, tlf)
    ref = something(reference, largest_injection_bus(g))
    RichCase(ElementaryCase(g, ref))
end

"The bus with the most negative net injection — the natural slack of a DC model."
largest_injection_bus(g::MetaGraph) =
    argmin(bus -> g[bus], collect(labels(g)))
