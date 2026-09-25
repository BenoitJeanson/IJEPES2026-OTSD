struct Pocket
    buses::Set{VLabel}
    branches::Set{ELabel}
    innerbranches::Set{ELabel}
    d::Float64
end

Base.show(io::IO, pk::Pocket) = print(io, "Pocket($(pk.buses), $(pk.branches), $(pk.innerbranches), $(pk.d))")

function Base.show(io::IO, ::MIME"text/plain", pk::Pocket)
    println(io, "Pocket:")
    println(io, "  buses:         ", pk.buses)
    println(io, "  branches:      ", pk.branches)
    println(io, "  innerbranches: ", pk.innerbranches)
    print(io,   "  d:             ", pk.d)
end

function pk_perimeter(buses::Set{VLabel}, outages::Set{ELabel})
    Set{ELabel}(
            outage for outage in outages
            if (from(outage) ∈ buses) ⊻ (to(outage) ∈ buses)
        )
end

function pk_innerbranches(buses::Set{VLabel}, g::MetaGraph)
    ib = Set{ELabel}()
    for bus in buses
        for br in incident(g, bus)
            if opposite(br, bus) in buses && !(br in ib)
                push!(ib, br)
            end
        end
    end
    ib
end

Pocket(g::MetaGraph, buses::Set{VLabel}, outages::Set{ELabel}) = 
    Pocket(buses, pk_perimeter(buses, outages), pk_innerbranches(buses, g), sum(max(g[bus], 0) for bus in buses))

function Pocket(g::MetaGraph, buses::Set{VLabel}, branches::Set{ELabel}, d::Float64)
    pk = Pocket(buses, branches, pk_innerbranches(buses, g), d)
end


"""
    create_bridge_to_pocket(ec::ElementaryCase, outages::Set{ELabel}=Set{ELabel}()) -> Dict{ELabel,Pocket}

Get the buses and branches that are shed in the given outages and bridges.
Returns a dictionary with the bridges as keys and the downstream Pocket as value.
"""
function create_bridge_to_pocket(
    ec::ElementaryCase,
    outages::Set{ELabel} = Set{ELabel}(),
)::Dict{ELabel,Pocket}
    g = ec.g;
    bus_orig = ec.bus_orig

    function r_visit_bridges!(
        bridge_to_buses,
        bus::VLabel,
        outages::Set{ELabel},
        bridges,
        visited_buses = VLabel[],
        crossed_bridges = ELabel[],
    )
        push!(visited_buses, bus)
        for br in incident(g, bus)
            if br in outages
                continue
            end
            other_bus = opposite(br, bus)
            other_bus in visited_buses && continue
            for bridge in crossed_bridges
                push!(bridge_to_buses[bridge], other_bus)
            end
            if br in bridges
                bridge_to_buses[br] = Set{VLabel}([other_bus])
                crossed_bridges2 = copy(crossed_bridges)
                push!(crossed_bridges2, br)
                r_visit_bridges!(
                    bridge_to_buses,
                    other_bus,
                    outages,
                    bridges,
                    visited_buses,
                    crossed_bridges2,
                )
            else
                r_visit_bridges!(
                    bridge_to_buses,
                    other_bus,
                    outages,
                    bridges,
                    visited_buses,
                    crossed_bridges,
                )
            end
        end
    end

    bridge_to_buses = Dict{ELabel,Set{VLabel}}()

    bridges = getbridges(g, outages)
    r_visit_bridges!(bridge_to_buses, bus_orig, outages, bridges)

    bridge_to_pocket = Dict{ELabel,Pocket}()
    for bridge in bridges
        !(bridge in keys(bridge_to_buses)) && continue
        buses = bridge_to_buses[bridge]
        d = sum(max(g[bus], 0) for bus in buses)
        bridge_to_pocket[bridge] = Pocket(g, buses, outages)
    end
    bridge_to_pocket
end
