using CSV
using DataFrames

mgx_genus = CSV.read("paper-taxonomic-levels/resonance_mgx_genus.csv", DataFrame)
dada2_genus = CSV.read("paper-taxonomic-levels/resonance_dada2_genera.csv", DataFrame)

dif = setdiff(replace.(dada2_genus.genus, Ref("g__"=>"")), mgx_genus.taxname)
uni = union(replace.(dada2_genus.genus, Ref("g__"=>"")), mgx_genus.taxname)


mpa_genera = String[]
for line in eachline("paper-taxonomic-levels/mpa_v30_CHOCOPhlAn_201901_marker_info.txt")
    for m in eachmatch(r"g__(\w+)", line)
        push!(mpa_genera, m.captures[1])
    end
end
mpa_genera

println.(setdiff(dif, mpa_genera))

open("paper-taxonomic-levels/genus_union.txt", "w") do io
    for g in uni
        endswith(g, "aceae") && continue
        occursin(r"\d", g) && continue
        println(io, g)
    end
end