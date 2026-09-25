"""Minimum hop distance from branch `br` to any branch in `targets`."""
min_hop(br::ELabel, targets, edm) =
    isempty(targets) ? typemax(Int) : minimum(edm[(br, t)] for t in targets)

extend_sbs_by_hops(sbs, all_branches, edm; d::Int=1) =
    Set(br for br in all_branches if min_hop(br, sbs, edm) ≤ d)


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