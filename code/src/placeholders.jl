# ── Inert extension points ────────────────────────────────────────────────────
#
# The model builders in `blocks.jl` are shared with two lines of work that are not
# part of this paper: substation reconfiguration (splitting a bus and sending its
# feeders to different sub-buses) and network equivalents (replacing an external
# area by a reduced model). Both enter as arguments that default to empty, and every
# use of them is guarded by an emptiness test, so with the definitions below the
# generated model is exactly the one the published runs solved.
#
# They are defined here, empty and unpopulated, rather than deleted, because
# removing the arguments touches model-building code across `blocks.jl` and the
# value of this package is that its model is the published one. Nothing in the
# package constructs a non-empty value of either type.

# ── Substation configurations ─────────────────────────────────────────────────

"One way of splitting a bus: which branches land on which sub-bus."
struct Conf
    subbus2br::Dict{Int,Vector{ELabel}}
    br2subbus::Dict{ELabel,Int}
end

const Confs = Vector{Conf}

"Candidate splits per bus. Always empty here: no bus is split."
const SubstationConfs = Dict{VLabel,Confs}

"A chosen split per bus. Always empty here."
const SubstationConf = Dict{VLabel,Conf}

subbuses(confs::Confs) = unique(sb for conf in confs for sb in keys(conf.subbus2br))

# ── Network equivalents ───────────────────────────────────────────────────────

"An external area replaced by a reduced model. Never constructed here."
struct Equivalent
    buses::Vector{VLabel}
end

"""
    EquivalentSet()

The equivalents applied to a case. Always empty in this package: both benchmark
systems are modelled in full.
"""
struct EquivalentSet
    eqs::Vector{Equivalent}
    connectors::Dict{VLabel,Float64}
    apply_branch_limits::Vector{Bool}
end

EquivalentSet() = EquivalentSet(Equivalent[], Dict{VLabel,Float64}(), Bool[])
