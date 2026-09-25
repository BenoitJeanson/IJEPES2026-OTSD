"""
    CutSink

Where separated cuts go, and where the incumbent is read from.

Two kinds. [`LazySink`](@ref) hands a cut straight to a solver that is waiting for
one, inside the branch-and-cut tree; a cut that happens to hold already is absorbed
harmlessly, so nothing is filtered. An [`AccumulatingSink`](@ref) adds cuts to the
model itself and must filter: a cut that is valid but not *violated* at the current
point changes nothing, and counting it would make an acceptable solution look
rejected. `ConshdlrSink` in `scip.jl` is the one in use.
"""
abstract type CutSink end

"""
    AccumulatingSink

A sink that adds cuts to the model rather than answering a solver's request for one.
Cuts are filtered by violation and counted, and that count is what tells the caller
whether anything was achieved — for SCIP's constraint handler, whether the solution
it is asking about is acceptable.
"""
abstract type AccumulatingSink <: CutSink end

"Cuts submitted as lazy constraints from inside the tree (Gurobi)."
struct LazySink{D} <: CutSink
    m::Model
    cb_data::D
end

model(s::CutSink) = s.m

"Value of `var` in the incumbent this sink is separating against."
solution_value(s::LazySink, var) = callback_value(s.cb_data, var)

"How many cuts this sink has taken. Zero means the point was already cut-feasible."
cuts_added(s::AccumulatingSink) = s.added[]

"When true, count violated cuts but do not add them — used to answer \"is this solution acceptable?\"."
is_dry(s::AccumulatingSink) = s.dry
