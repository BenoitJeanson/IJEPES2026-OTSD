# ── Where separated cuts go ───────────────────────────────────────────────────
#
# A cut is built once and delivered in one of two ways. Inside a branch-and-cut tree
# it is submitted as a lazy constraint against the callback handle; between master
# solves it is simply added to the model. `cut_constraint` is the single definition
# of what each cut family says, so the two paths cannot drift apart.

"""
    CutSink

Where separated cuts go, and where the incumbent is read from. [`LazySink`](@ref)
talks to a solver callback, [`DirectSink`](@ref) to a solved model.
"""
abstract type CutSink end

"Cuts submitted as lazy constraints from inside the tree (Gurobi)."
struct LazySink{D} <: CutSink
    m::Model
    cb_data::D
end

"Cuts added to the master between solves (HiGHS, and any solver without callbacks)."
struct DirectSink <: CutSink
    m::Model
    added::Base.RefValue{Int}
end

DirectSink(m::Model) = DirectSink(m, Ref(0))

model(s::CutSink) = s.m

"Value of `var` in the incumbent this sink is separating against."
solution_value(s::LazySink, var) = callback_value(s.cb_data, var)
solution_value(s::DirectSink, var) = value(var)

"How many cuts this sink has received. Meaningful for `DirectSink`; drives the loop."
cuts_added(s::DirectSink) = s.added[]

