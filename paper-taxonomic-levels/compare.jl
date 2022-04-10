using CSV
using DataFrames
using Statistics
using PrettyTables

mgx_family = CSV.read("paper-taxonomic-levels/resonance_mgx_family.csv", DataFrame)
mgx_genus = CSV.read("paper-taxonomic-levels/resonance_mgx_genus.csv", DataFrame)
mgx_species = CSV.read("paper-taxonomic-levels/resonance_mgx_species.csv", DataFrame)

dada2_family = CSV.read("paper-taxonomic-levels/resonance_dada2_families.csv", DataFrame)
dada2_genus = CSV.read("paper-taxonomic-levels/resonance_dada2_genera.csv", DataFrame)
dada2_species = CSV.read("paper-taxonomic-levels/resonance_dada2_species.csv", DataFrame)

summary = DataFrame(method =repeat(["mgx", "amplicon"], outer=3), 
                    level  = repeat(["family", "genus", "species"], inner=2),
                    median = zeros(6),
                    mean   = zeros(6),
                    stdev  = zeros(6))

function unclassified_summary(df, col)
    uncl = [x for x in filter(col => ==("UNCLASSIFIED"), df)[1, 2:end]]
    return (median(1 .- uncl), mean(1 .- uncl), std(1 .- uncl))
end


summary[1, 3:end] .= unclassified_summary(mgx_family, :taxname)
summary[2, 3:end] .= unclassified_summary(dada2_family, :family)
summary[3, 3:end] .= unclassified_summary(mgx_genus, :taxname)
summary[4, 3:end] .= unclassified_summary(dada2_genus, :genus)
summary[5, 3:end] .= unclassified_summary(mgx_species, :taxname)
summary[6, 3:end] .= unclassified_summary(dada2_species, :species)
