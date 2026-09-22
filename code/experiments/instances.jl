# ── The two published instances ───────────────────────────────────────────────
#
# Each instance fixes what the paper fixes: the PGLib case, the thermal limit factor,
# the contingencies embedded in the master from the outset, and the warm start.
#
# The warm start is the output of the expert heuristic of the companion paper
# (Jeanson, Tanneau & Tindemans, "Scalable Iterative Algorithm for Solving Optimal
# Transmission Switching with De-energization"). That heuristic is not part of this
# package, so its result is recorded here as data. Every reported figure starts from
# these openings; starting anywhere else changes the search and therefore the numbers.

"The operating point the paper reports: H = 3, d_viol = 2, d_sol = 0, k = 3."
const PAPER_OPERATING_POINT = (; H = 3, d_viol = 2, d_sol = 0, k = 3, seed = 0)

struct Instance
    name::String
    case::String
    tlf::Float64
    embedded::Vector{ELabel}
    warm_start::Set{ELabel}
end

const IEEE118 = Instance(
    "118", "case118", 1.5,
    ELabel[("9", "10"), ("37", "38"), ("65", "66"), ("8", "9"),
           ("38", "65"), ("64", "65"), ("5", "8")],
    Set{ELabel}([("54", "55"), ("54", "56"), ("51", "58"), ("49", "54"), ("49", "51"),
                 ("13", "15"), ("40", "41"), ("40", "42"), ("48", "49"), ("30", "38")]),
)

const IEEE57 = Instance(
    "57", "case57", 1.0,
    ELabel[("35", "36"), ("7", "8"), ("8", "9"), ("6", "8"), ("36", "37"), ("37", "38")],
    Set{ELabel}([("24", "25"), ("23", "24"), ("39", "57"), ("9", "55"),
                 ("4", "5"), ("4", "6"), ("40", "56")]),
)

const INSTANCES = Dict("118" => IEEE118, "57" => IEEE57)

"""
    reference_objective(instance) -> Float64

The objective the paper reports at [`PAPER_OPERATING_POINT`], in per unit on a
100 MVA base. Used as the acceptance test for this package: the Gurobi backend must
reproduce it exactly.
"""
reference_objective(i::Instance) = i.name == "118" ? 5.29 : 7.382
