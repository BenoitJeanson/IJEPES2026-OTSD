"""
A `ElementaryCase` is the simplest GridCase: a graph with a reference bus.
`bus_orig` fixes the reference bus — buses not reachable from it through
non-open branches are treated as de-energized.
"""
struct ElementaryCase <: GridCase
    g        ::MetaGraph
    bus_orig ::VLabel
end

Base.copy(ec::ElementaryCase) = ElementaryCase(copy(ec.g), ec.bus_orig)

connectedcomponent(ec::ElementaryCase, outages::Set{ELabel}) =
    connectedcomponent(ec.g, ec.bus_orig, outages)