
const ELabel = Tuple{String,String}
const VLabel = String

mutable struct Branch
    b::Float64
    p_max::Float64
    v_nom1::Float64
    v_nom2::Float64
    p::Float64
end

Base.show(io::IO, b::Branch) = print(
    io,
    "Branch(b=$(b.b), p_max=$(b.p_max), v_nom1=$(b.v_nom1), v_nom2=$(b.v_nom2), p=$(b.p))",
)

Branch(b, p_max) = Branch(b, p_max, 0, 0, 0)
Branch(b, p_max, v_nom1, v_nom2) = Branch(b, p_max, v_nom1, v_nom2, 0)
Branch(b::Branch) = Branch(b.b, b.p_max, b.v_nom1, b.v_nom2, b.p)
Branch(b::Branch, p_max) = Branch(b.b, p_max, b.v_nom1, b.v_nom2, b.p)
Branch() = Branch(0, 0, 0, 0, 0)

function e_label_for(g::MetaGraph, e::Edge)
    l1 = label_for(g, src(e))
    l2 = label_for(g, dst(e))
    haskey(g, l1, l2) ? (l1, l2) : (l2, l1)
end

from(e::ELabel) = e[1]
to(e::ELabel) = e[2]

str(e::ELabel) = e[1] * "-" * e[2]
str(branches::Union{AbstractSet{ELabel},AbstractVector{ELabel}}) = join(str.(branches), ", ")

lsrc(g::MetaGraph, e::Graphs.SimpleEdge) = label_for(g, src(e))
ldst(g::MetaGraph, e::Graphs.SimpleEdge) = label_for(g, dst(e))

incident_signed(g::MetaGraph, bus::VLabel) = Iterators.flatten((
    (((busin, bus), 1) for busin in inneighbor_labels(g, bus)),
    (((bus, busout), -1) for busout in outneighbor_labels(g, bus)),
))

incident(g::MetaGraph, bus::VLabel) = Iterators.flatten((
    ((busin, bus) for busin in inneighbor_labels(g, bus)),
    ((bus, busout) for busout in outneighbor_labels(g, bus)),
))

neighbor_labels(g::MetaGraph, bus::VLabel) =
    Iterators.flatten((inneighbor_labels(g, bus), outneighbor_labels(g, bus)))

opposite(e::ELabel, v::VLabel) = v == e[1] ? e[2] : e[1]

openbranchesset(outages::Set{ELabel}, tripping::Union{Nothing,ELabel} = nothing) =
    isnothing(tripping) ? outages : union(outages, Set([tripping]))

function _initgraph(;
    edge_label_builder = (bus1, bus2) -> bus1 < bus2 ? (bus1, bus2) : (bus2, bus1),
)
    MetaGraph(
        DiGraph();
        label_type = String,
        vertex_data_type = Float64,
        edge_data_type = Branch,
        graph_data = (edge_label_builder = edge_label_builder,),
    )
end

function build_simple_grid(; micro = true)
    g = _initgraph()
    g["1"] = micro ? -2 : -3
    for i = 2:(micro ? 3 : 4)
        g["$i"] = 1
    end
    g["1", "2"] = Branch(1, 1)
    g["2", "3"] = Branch(1, 1)
    if micro
        g["1", "3"] = Branch(1, 1)
    else
        g["3", "4"] = Branch(1, 1)
        g["2", "4"] = Branch(1, 1)
        g["1", "4"] = Branch(1, 1)
    end
    g
end

function PGLibtograph(
    case::String;
    buslabels = num -> "$num",
    edge_order = (bus1, bus2) -> parse(Int, bus1) < parse(Int, bus2),
)
    c = pglib(case)

    edge_label_builder(bus1, bus2) = edge_order(bus1, bus2) ? (bus1, bus2) : (bus2, bus1)

    g = _initgraph(; edge_label_builder = edge_label_builder)

    for (label, bus) in c["bus"]
        g[buslabels(bus["bus_i"])] = 0
    end

    for load in values(c["load"])
        g[buslabels(load["load_bus"])] += load["pd"]
    end
    for gen in values(c["gen"])
        g[buslabels(gen["gen_bus"])] -= gen["pg"]
    end

    for br in values(c["branch"])
        br_label = edge_label_builder(buslabels(br["f_bus"]), buslabels(br["t_bus"]))
            # edge_order(buslabels(br["f_bus"]), buslabels(br["t_bus"])) ?
            # (buslabels(br["f_bus"]), buslabels(br["t_bus"])) :
            # (buslabels(br["t_bus"]), buslabels(br["f_bus"]))
        g[br_label...] = Branch(1 / br["br_x"], br["rate_a"])
    end

    to_remove = [bus for bus in labels(g) if isempty(neighbor_labels(g, bus))]
    foreach(bus -> rem_vertex!(g, code_for(g, bus)), to_remove)

    g
end

@enum BalanceType all_non_zero_uniform gen_proportional

function balance!(g::MetaGraph, btype::BalanceType = gen_proportional)
    non_zeros = [label for label in labels(g) if g[label] ≠ 0]
    imbalance = sum(g[label] for label in non_zeros)

    if btype == all_non_zero_uniform
        uniform = imbalance / length(non_zeros)
        for label in non_zeros
            g[label] -= uniform
        end

    elseif btype == gen_proportional
        generators = [label for label in labels(g) if g[label] ≤ 0]
        if isempty(generators)
            println("no generator to balance")
        else
            total_gen = sum(g[label] for label in generators)
            for label in generators
                g[label] -= imbalance * g[label] / total_gen
            end
        end

    else
        println("BALANCE TYPE TO BE IMPLEMENTED")
    end #TODO oter types
end

function check_flow_consistency(g; v::Bool = false)
    flows = Dict(l => g[l...].p for l in edge_labels(g))
    injections = Dict(l => g[l] for l in labels(g))

    # Initialize net flows for each node
    net_flows = Dict{String,Float64}()

    # Update net flows based on branch flows
    for ((from, to), flow) in flows
        net_flows[from] = get(net_flows, from, 0.0) - flow
        net_flows[to] = get(net_flows, to, 0.0) + flow
    end

    # Compare with injections
    consistent = true

    v && println("\nComparison with injections:")
    for (node, injection) in injections
        net_flow = get(net_flows, node, 0.0)

        v && println(
            "Node $node: Injection = $injection, Net Flow = $net_flow, Difference = $(injection - net_flow)",
        )
        if abs(injection - net_flow) > 1e-9  # Tolerance for floating-point comparison
            consistent = false
        end
    end

    if v
        if consistent
            println("\nThe flows and injections are consistent.")
        else
            println("\nThe flows and injections are not consistent.")
        end
    end
    consistent
end

function scale_branch_limits!(g::MetaGraph, ratio)
    for br in edge_labels(g)
        g[br...].p_max *= ratio
    end
end

function create_graph_wo_outages(g::MetaGraph, outages::Set{ELabel})
    h = copy(g)
    for br in outages
        delete!(h, br...)
    end
    h
end

function getbridges(g::MetaGraph, outages::Set{ELabel} = Set{ELabel}())::Vector{ELabel}
    h = create_graph_wo_outages(g, outages)
    _bridges = Graphs.bridges(Graph(h.graph))
    [e_label_for(h, br) for br in _bridges]
end

function connectedcomponent(g::MetaGraph, bus::VLabel, outages::Set{ELabel})

    function _expand!(cc_buses, cc_edges, bus)
        push!(cc_buses, bus)
        if haskey(g, bus)
            for edg in incident(g, bus)
                (edg ∈ outages || edg ∈ cc_edges) && continue
                push!(cc_edges, edg)
                nb = opposite(edg, bus)
                nb ∉ cc_buses && _expand!(cc_buses, cc_edges, nb)
            end
        end
    end

    buses, edges = Set{VLabel}(), Set{ELabel}()
    _expand!(buses, edges, bus)
    return (buses = buses, edges = edges)
end

function connectedbusessets(g::MetaGraph, outages = Set{ELabel}())::Vector{Set{VLabel}}
    remaining = Set(labels(g))
    ccs = Vector{Set{VLabel}}()
    while !isempty(remaining)
        component = Set{VLabel}()
        stack = [first(remaining)]
        while !isempty(stack)
            bus = pop!(stack)
            bus ∈ remaining || continue
            push!(component, bus)
            delete!(remaining, bus)
            for br in incident(g, bus)
                br in outages && continue
                nb = opposite(br, bus)
                nb ∈ remaining && push!(stack, nb)
            end
        end
        push!(ccs, component)
    end
    return ccs
end

function edge_distance_map(g::MetaGraph)
    lg = SimpleGraph(ne(g))

    for (i, e1) in enumerate(edges(g))
        for (j, e2) in enumerate(edges(g))
            # connect if they share a vertex
            if length(intersect([src(e1), dst(e1)], [src(e2), dst(e2)])) > 0
                add_edge!(lg, i, j)
            end
        end
    end

    sp = floyd_warshall_shortest_paths(lg)

    edge_distance_map = Dict{Tuple{ELabel,ELabel},Int64}(
        (ei, ej) => sp.dists[i, j] for (i, ei) in enumerate(edge_labels(g)) for
        (j, ej) in enumerate(edge_labels(g))
    )

    return edge_distance_map
end

"""
For a given edge, returns the neighboring edges reached by graph paths of lengths lower than depth.
"""
function edge_neighborhood(branch::ELabel, g::MetaGraph, depth::Int)
    function _r_neighbors!(visited::Dict{ELabel,Int}, bus::VLabel, actual_depth::Int)
        for br in incident(g, bus)
            br in keys(visited) && visited[br] ≥ actual_depth && continue
            visited[br] = actual_depth
            actual_depth ≥ 1 && _r_neighbors!(visited, opposite(br, bus), actual_depth - 1)
        end
    end

    result = Dict{ELabel,Int}(branch => depth)
    _r_neighbors!(result, from(branch), depth - 1)
    _r_neighbors!(result, to(branch), depth - 1)
    keys(result)
end


"""
In a Dictionary, returns for each edge the neighboring edges reached by graph paths of lengths lower than depth.
"""
function edge_neighborhood(g::MetaGraph, depth::Int)
    function r_visit!(visited, visited_nodes, remaining_depth, g, br)
        push!(visited, br)
        remaining_depth == 0 && return
        for bus in br
            bus in visited_nodes && continue
            push!(visited_nodes, bus)
            for b_r in incident(g, bus)
                b_r in visited && continue
                r_visit!(visited, visited_nodes, remaining_depth - 1, g, b_r)
            end
        end
    end

    neighbors = Dict{ELabel,Vector{ELabel}}()
    for br in edge_labels(g)
        visited = ELabel[]
        visited_nodes = VLabel[]
        r_visit!(visited, visited_nodes, depth, g, br)
        neighbors[br] = visited
    end
    neighbors
end

function hamming(v1, v2)
    s1 = Set(v1)
    s2 = Set(v2)
    length(union(setdiff(s1, s2), setdiff(s2, s1)))
end


function sub_graph(g::MetaGraph, buses_branches::Dict{VLabel,Vector{ELabel}}) #TODO: consider not building remaining graph.

    function _r_add_bus!(sub, buses_to_rm, edges_to_rm, boundariestovisit, g, bus_label)
        if !haskey(sub, bus_label)
            sub[bus_label] = g[bus_label]
            if haskey(buses_branches, bus_label)
                delete!(boundariestovisit, bus_label)
                _edges = buses_branches[bus_label]
            else
                push!(buses_to_rm, bus_label)
                _edges = incident(g, bus_label)
            end
            for e in _edges
                bus2 = opposite(e, bus_label)
                _r_add_bus!(sub, buses_to_rm, edges_to_rm, boundariestovisit, g, bus2)
                sub[e...] = g[e...]
                push!(edges_to_rm, e)
            end
        end
        !isempty(boundariestovisit) && _r_add_bus!(
            sub,
            buses_to_rm,
            edges_to_rm,
            boundariestovisit,
            g,
            first(boundariestovisit),
        )
    end

    sub = _initgraph()
    buses_to_rm = VLabel[]
    edges_to_rm = ELabel[]
    boundariestovisit = Set{VLabel}(keys(buses_branches))
    _r_add_bus!(
        sub,
        buses_to_rm,
        edges_to_rm,
        boundariestovisit,
        g,
        first(boundariestovisit),
    )
    remaining = copy(g)
    foreach(e -> delete!(remaining, e[1], e[2]), edges_to_rm)
    foreach(b -> delete!(remaining, b), buses_to_rm)
    sub, remaining
end