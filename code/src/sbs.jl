"""Minimum hop distance from branch `br` to any branch in `targets`."""
min_hop(br::ELabel, targets, edm) =
    isempty(targets) ? typemax(Int) : minimum(edm[(br, t)] for t in targets)

extend_sbs_by_hops(sbs, all_branches, edm; d::Int=1) =
    Set(br for br in all_branches if min_hop(br, sbs, edm) ≤ d)

"""
Fixed-core SBS: branches within `d_fixed` hops of violated branches and
violating N-1 contingencies. Computed once at startup.
"""
function compute_fixed_sbs(rc::RichCase, all_branches::Vector{ELabel}, edm; d_fixed::Int=2)
    BASECONTINGENCY = ELabel(("", ""))
    sa0 = secured_dcpf(rc.gc)
    viol_dict = violated_branches(sa0)
    all_viol = union((v for v in values(viol_dict))...)
    viol_ctgs = Set(k for (k, v) in viol_dict if k ≠ BASECONTINGENCY && !isempty(v))
    # targets   = union(all_viol, viol_ctgs)
    targets = all_viol
    @info "Fixed SBS targets: $(length(all_viol)) violated branches + $(length(viol_ctgs)) violating contingencies"
    extend_sbs_by_hops(targets, all_branches, edm; d=d_fixed)
end

"""
Canonical dynamic SBS builder.
Returns fixed_sbs ∪ d_dynamic-hop neighborhood of current_openings.
"""
function build_sbs(all_branches::Vector{ELabel}, fixed_sbs::Set{ELabel},
    current_openings, edm; d_dynamic::Int=3)
    ops_set = Set{ELabel}(current_openings)
    dynamic_sbs = extend_sbs_by_hops(ops_set, all_branches, edm; d=d_dynamic)
    union(fixed_sbs, dynamic_sbs)
end

"""
Solution-induced SBS: for each open branch, close it and run SA; collect induced violations.
SBS = current_openings ∪ d_induced-hop neighborhood of all induced violations.
Returns NamedTuple (sbs, sbs_size, n_viol, per_branch).
"""
function solution_induced_sbs(rc::RichCase, all_branches::Vector{ELabel}, edm,
    openings; d_induced::Int=1)
    all_induced = Set{ELabel}()
    per_branch = Pair{ELabel,Set{ELabel}}[]
    for b in openings
        partial = setdiff(Set{ELabel}(openings), (b,))
        sa_b = secured_dcpf(rc.gc, partial)
        vd = violated_branches(sa_b)
        viols_b = isempty(vd) ? Set{ELabel}() : union(values(vd)...)
        push!(per_branch, b => viols_b)
        union!(all_induced, viols_b)
    end
    hop_nbhd = isempty(all_induced) ? Set{ELabel}() :
               extend_sbs_by_hops(all_induced, all_branches, edm; d=d_induced)
    sbs = union(Set{ELabel}(openings), hop_nbhd)
    (sbs=sbs, sbs_size=length(sbs), n_viol=length(all_induced), per_branch=per_branch)
end


function follow_the_flows(sa_res::SA_result, contingency::ELabel, branch::ELabel, max_hops::Int=typemax(Int))::Set{ELabel}
    function _r_follows!(branches, br, is_forward, remaining_hops)
        remaining_hops == 0 && return
        fl = flow(sa_res, contingency, br)
        b, f = fl > 0 ? (br[1], br[2]) : (br[2], br[1])

        starting_bus = is_forward ? f : b
        for next in incident(g, starting_bus)
            next == br && continue
            flow_f = (next[1] == starting_bus ? 1 : -1) * flow(sa_res, contingency, next)
            ((is_forward && flow_f ≤ 0) || (!is_forward && flow_f ≥ 0)) && continue
            push!(branches, next)
            _r_follows!(branches, next, is_forward, remaining_hops - 1)
        end
    end

    g = sa_res.gc.g
    branches = Set{ELabel}()
    push!(branches, branch)
    _r_follows!(branches, branch, true, max_hops)
    _r_follows!(branches, branch, false, max_hops)
    branches
end

function follow_the_flows(sa_res::SA_result, max_hops::Int=typemax(Int))::Set{ELabel}
    ctg2vbs = violated_branches(sa_res)
    reduce(union, (follow_the_flows(sa_res, ctg, br, max_hops) for (ctg, vbs) in ctg2vbs for br in vbs), init=Set{ELabel}())
end

function sa_induced_followed(rc::RichCase, solution_openings::Set{ELabel}, max_hops::Int=typemax(Int))
    sbs = Set{ELabel}()
    union!(sbs, solution_openings)
    for ctg in solution_openings
        _openings = setdiff(solution_openings, [ctg])
        sa_res2 = secured_dcpf(rc.gc, Set(_openings))
        union!(sbs, follow_the_flows(sa_res2, max_hops))
    end
    sbs
end