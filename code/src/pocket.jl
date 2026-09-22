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

function create_systematic_pocket(g::MetaGraph)::Vector{Pocket}
    sp = Pocket[]
    for bus in labels(g)
        max(g[bus], 0) ≠ 0 &&
            push!(sp, Pocket(g, Set([bus]), Set(incident(g, bus)), g[bus]))
    end

    for br in edge_labels(g)
        f, t = from(br), to(br)
        branches = union(Set(incident(g, f)), Set(incident(g, t)))
        filter!(b -> b ≠ br, branches)
        d = max(g[f], 0) + max(g[t], 0)
        d ≠ 0 && push!(sp, Pocket(g, Set{VLabel}([f, t]), branches, d))
    end
    sp
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

"""
Create a pocket that contains the bus_orig and for which the branches are the bridges and borders of the pockets that are not embeded to a bigger one.
"""
function create_bus_orig_pocket(
    ec::ElementaryCase,
    bridge_to_pocket::Dict{ELabel,Pocket},
)::Pocket
    g = ec.g;
    bus_orig = ec.bus_orig
    pocket_borders = Set([br for pk in values(bridge_to_pocket) for br in pk.branches])
    pocket_bridges = Set(keys(bridge_to_pocket))
    union!(pocket_borders, pocket_bridges)

    function _r_visit_graph(bus, buses, branches, innerbranches)
        bus in buses && return 0
        push!(buses, bus)
        d = max(g[bus], 0)
        for br in incident(g, bus)
            if br in pocket_borders
                br ∉ branches && push!(branches, br)
                continue
            end
            br ∉ innerbranches && push!(innerbranches, br)
            d += _r_visit_graph(opposite(br, bus), buses, branches, innerbranches)
        end
        return d
    end
    buses, branches, innerbranches = Set{VLabel}(), Set{ELabel}(), Set{ELabel}()
    d = _r_visit_graph(bus_orig, buses, branches, innerbranches)
    Pocket(buses, branches, innerbranches, d)
end

"""
Input: bus_orig: if pocket is the main pocket, then the bus_orig is the one of the main one else it is the bus of the bridge that is in the pocket.
"""
function create_pocket_subgraph(
    ec::ElementaryCase,
    openbranches::Set{ELabel},
    pk::Pocket,
    bridge::ELabel,
)
    g = ec.g;
    bus_orig = ec.bus_orig
    dcpf_res = dcpf(ec; outages = collect(openbranches))
    h = _initgraph()
    for bus in pk.buses
        h[bus] = g[bus]
    end
    if bridge == case_to_branch(BASECASEID)
        for br in pk.branches
            busin = from(br) in pk.buses ? from(br) : to(br)
            h[busin] += (busin == from(br) ? 1 : -1) * dcpf_res.flows[br...]
        end
    else
        h[bus_orig] += (bus_orig == from(bridge) ? 1 : -1) * dcpf_res.flows[bridge...]
    end

    for br in pk.innerbranches
        h[br...] = g[br...]
    end
    return h
end
