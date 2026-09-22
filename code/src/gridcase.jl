"""
A `GridCase` is an object which a DCPF can be computed on.
"""
abstract type GridCase end

"""
A `RichCase` is a `GridCase` with additional information for drawing.
`trippings`: is a list of trippings for which the order is important.
`bus_orig` is stored in the inner `GridCase` (e.g. `ElementaryCase` or `EquivalentCase`).
"""
struct RichCase{T<:GridCase}
    gc        ::T
    drawparams::Dict{Symbol,Any}
    trippings ::Union{Nothing,Vector{ELabel}}
end

RichCase(gc::T, drawparams::Dict{Symbol,Any}=Dict{Symbol,Any}()) where T<:GridCase =
    RichCase(gc, drawparams, nothing)
