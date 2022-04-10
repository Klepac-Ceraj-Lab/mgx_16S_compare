# ---
# author:
# - Kevin Bonham
# - PhD
# title: RESONANCE amplicon updated
# ---
#
# # RESONANCE amplicon updated
#
# ## Running the pipeline
#
# Using QIIME.jl
#
# ```sh
# qiimejl -n test -i v4v5_seqs -o ./ --threads 16 -v
# ```
#
# ### Massage dataframe

using CSV
using DataFrames
using Chain
using MGXAmplicon: norm_taxon


features = CSV.read("/augusta/staff/kevin/echo/echo_amplicon/pipeline/dada2/feature-table.tsv", DataFrame, header=2, delim='\t')
taxa = CSV.read("/augusta/staff/kevin/echo/echo_amplicon/pipeline/taxa/taxonomy.tsv", DataFrame, delim='\t')

rename!(features, "#OTU ID"=>"feature")
rename!(taxa, "Feature ID"=>"feature")

labeled = leftjoin(select(taxa, ["feature", "Taxon"]), features, on="feature")
# some columns don't have any taxonomic assigment, remove them
labeled = labeled[!, [!(eltype(labeled[!, col]) <: Union{Missing, <:Number}) || 
                        sum(labeled[!, col]) > 0 
                        for col in names(labeled)]]

# convert to relative abundance
select!(labeled, [:Taxon => (x -> identity(x)), 
                 [C => c -> c ./ sum(c) for C in names(labeled, r"^[CM]")]...],
                 renamecols=false)
any(isnan, [sum(labeled[:, i]) for i in 2:ncol(labeled)-2])

# Now, we need to massage the labels a bit,
# since we have things like `s__Prevotella_sp.` listed as a species,
# and things like ` g__[Ruminococcus]_gnavus_group` not resolved to species.

transform!(labeled, :Taxon => ByRow(x-> norm_taxon(x; splitter=';', kind=:amplicon)) => [:taxlevel, :taxname])

labeled.family = [row.taxlevel in (:genus, :species, :family) ? String(split(row.taxname, ';')[5]) : "UNCLASSIFIED" for row in eachrow(labeled)]
labeled.genus = [row.taxlevel in (:genus, :species) ? String(split(row.taxname, ';')[6]) : "UNCLASSIFIED" for row in eachrow(labeled)]
labeled.species = [row.taxlevel != :species ? "UNCLASSIFIED" : String(split(row.taxname, ';')[7]) for row in eachrow(labeled)]

fa = @chain labeled begin
    groupby(:family)
    combine([n => sum for n in names(labeled, r"^[CM]")], renamecols=false)
end

ge = @chain labeled begin
    groupby(:genus)
    combine([n => sum for n in names(labeled, r"^[CM]")], renamecols=false)
end


sp = @chain labeled begin
    groupby(:species)
    combine([n => sum for n in names(labeled, r"^[CM]")], renamecols=false)
end

CSV.write("/home/kevin/repos/danielle-thesis/paper-taxonomic-levels/resonance_dada2_families.csv", fa)
CSV.write("/home/kevin/repos/danielle-thesis/paper-taxonomic-levels/resonance_dada2_genera.csv", ge)
CSV.write("/home/kevin/repos/danielle-thesis/paper-taxonomic-levels/resonance_dada2_species.csv", sp)

# ### Get some statistics

using Statistics

sp_totals = sum.(eachcol(filter(:species=> !=("UNCLASSIFIED"), sp)[!, r"^[CM]"]));
mean(sp_totals)
median(sp_totals)

ge_totals = sum.(eachcol(filter(:genus => !=("UNCLASSIFIED"), ge)[!, r"^[CM]"]));
mean(ge_totals)
median(ge_totals)
