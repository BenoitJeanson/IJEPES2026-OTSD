"""
    CutSink

Where separated cuts go, and where the incumbent is read from.

Two kinds. [`LazySink`](@ref) hands a cut straight to a solver that is waiting for
one, inside the branch-and-cut tree; a cut that happens to hold already is absorbed
harmlessly, so nothing is filtered. The others accumulate cuts into the model and
must filter: a cut that is valid but not *violated* at the current point changes
nothing, and adding it anyway makes progress indistinguishable from deadlock.
"""
abstract type CutSink end

"""
    AccumulatingSink

A sink that adds cuts to the model rather than answering a solver's request for one.
Cuts are filtered by violation and counted, and that count is what tells the caller
whether anything was achieved.
"""
abstract type AccumulatingSink <: CutSink end

"Cuts submitted as lazy constraints from inside the tree (Gurobi)."
struct LazySink{D} <: CutSink
    m::Model
    cb_data::D
end

"Cuts added to the master between solves (a solver with no callback at all)."
struct DirectSink <: AccumulatingSink
    m::Model
    added::Base.RefValue{Int}
    dry::Bool
end

DirectSink(m::Model; dry::Bool = false) = DirectSink(m, Ref(0), dry)

model(s::CutSink) = s.m

"Value of `var` in the incumbent this sink is separating against."
solution_value(s::LazySink, var) = callback_value(s.cb_data, var)
solution_value(s::DirectSink, var) = value(var)

"How many cuts this sink has taken. Zero means the point was already cut-feasible."
cuts_added(s::AccumulatingSink) = s.added[]

"When true, count violated cuts but do not add them — used to answer \"is this solution acceptable?\"."
is_dry(s::AccumulatingSink) = s.dry
