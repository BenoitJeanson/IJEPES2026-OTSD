# Extract the authoritative SBS for every phase of every run.
#
# The text logs record |SBS| but never its membership, and the union of branches
# actually opened across a phase's candidates covers only ~39% of it -- so the SBS
# has to come from the serialized solver state, not be inferred.
#
# `i{N}_state.jls` holds `state.accumulated_sbs` AFTER phase N grew it, i.e. the set
# phase N+1 searches in. Phase 1's set is never serialized, so it is recomputed with
# TNROpt's own `sa_induced_followed` -- the same call ablation.jl makes at line 256.

using TNROpt, Serialization, JSON3, Dates
using TNROpt: ELabel, create_case, sa_induced_followed

const ABL = "/Users/benoitjeanson/vsCode/TUD/tnr/tmp/IJEPES/ablation"

case_data(s) = s == "ieee118" ?
    (case="case118", ratio=1.5, heur=Set{ELabel}([("54","55"),("54","56"),("51","58"),
        ("49","54"),("49","51"),("13","15"),("40","41"),("40","42"),("48","49"),("30","38")])) :
    (case="case57", ratio=1.0, heur=Set{ELabel}([("24","25"),("23","24"),("39","57"),
        ("9","55"),("4","5"),("4","6"),("40","56")]))

lab(e) = "$(e[1])-$(e[2])"
const RC = Dict{String,Any}()
rich(sys) = get!(RC, sys) do
    cd = case_data(sys); create_case(cd.case, cd.ratio)
end

out = Dict{String,Any}()
mism = String[]
dirs = sort(filter(d -> !startswith(d, "_") && isdir(joinpath(ABL, d)), readdir(ABL)))

for d in dirs
    rp = joinpath(ABL, d, "result.json")
    isfile(rp) || continue
    res = JSON3.read(read(rp, String))
    man = JSON3.read(read(joinpath(ABL, d, "manifest.json"), String))
    sys = String(res["system"])
    # Files are stamped yyyy-mm-dd_HHMMSS from a clock read just before the manifest
    # is written, so the two can differ by a second; and a directory may hold two
    # campaigns. Take the stamp nearest the manifest's start time.
    t0 = DateTime(first(split(String(man["started_at"]), ".")))
    stamps = unique([m.match for m in eachmatch(r"\d{4}-\d{2}-\d{2}_\d{6}",
                                                join(readdir(joinpath(ABL, d)), " "))])
    stamp = argmin(s -> abs(DateTime(s, "yyyy-mm-dd_HHMMSS") - t0), stamps)

    cd = case_data(sys)
    # phase 1 searches the heuristic seed's SA-induced set
    init = res["config"] == "NO-LOCALSEARCH" ?
        Set(TNROpt.edge_labels(rich(sys).gc.g)) :
        sa_induced_followed(rich(sys), cd.heur, res["d_viol"])
    length(init) == res["init_sbs"] ||
        push!(mism, "$d init $(length(init)) vs $(res["init_sbs"])")

    used = Dict{String,Any}("1" => sort(lab.(collect(init))))
    after = Dict{String,Any}()
    for it in res["iterations"]
        n = it["iteration"]
        cands = filter(f -> startswith(basename(f), stamp),
                       filter(f -> occursin(Regex("_i0*$(n)_state\\.jls\$"), f),
                              readdir(joinpath(ABL, d), join=true)))
        if isempty(cands)
            push!(mism, "$d i$n: no state.jls for stamp $stamp")
            continue
        end
        st = open(deserialize, cands[end])
        s = sort(lab.(collect(st.accumulated_sbs)))
        after[string(n)] = s
        used[string(n + 1)] = s
        length(s) == it["sbs_size"] ||
            push!(mism, "$d i$n after $(length(s)) vs $(it["sbs_size"])")
    end
    out[d] = Dict("used" => used, "after" => after)
    print(".")
end

println("\nruns: ", length(out), "   mismatches vs logged |SBS|: ", length(mism))
for m in first(mism, 10); println("  ", m); end
open("/private/tmp/claude-501/-Users-benoitjeanson-vsCode-TUD-tnr/3d72863e-a524-46cb-a45c-fd004df95463/scratchpad/sbs.json", "w") do io
    JSON3.write(io, out)
end
