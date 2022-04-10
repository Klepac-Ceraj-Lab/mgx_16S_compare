# ---
# author = "Kevin Bonham, PhD"
# title = "MGX taxonomic profiles"
# ---
#
# # MGX taxonomic profiles
#
# Begining from MetaPhlAn (v3) tables for each sample,
# we want to generate tables that include only species or genera,
# which sum to 1 (including unidentified). 
# This is easier than with QIIME, since in MetaPhlAn tables,
# each taxonomic level already sums to 100,
# so it's just a matter of filtering.
## load packages
using DataFrames
using CSV
using Microbiome
using Chain
using MGXAmplicon: norm_taxon

inpath = "/augusta/echo/analysis/diabimmune/"
profiles = String[]

## Recurssively walk directory, find files in `metaphlan/` output directories,
## and add them to the list of profiles
for (root, dirs, files) in walkdir(inpath)
    occursin("metaphlan", root) || continue
    filter!(f-> occursin("profile", f), files)
    append!(profiles, joinpath.(Ref(root), files))
end

@info first(profiles, 5)

# Now that we have all of the file paths, we will load each file
# and append it to a DataFrame in tidy format, with columns for
# `sample`, `taxname`, `taxlevel`, and `abundance`:

profilesdf = DataFrame()

for p in profiles
    profile = CSV.read(p, DataFrame, delim='\t', header=["taxon", "taxid", "abundance", "adtl"], datarow=5)
    sample = replace(splitext(basename(p))[1], r"_profile"=>"")
    transform!(profile, :taxon => ByRow(norm_taxon) => [:taxlevel, :taxname])
    profile[!, :sample] .= sample
    append!(profilesdf, select(profile, ["sample", "taxname", "taxlevel", "abundance"]))
end

profilesdf

# Now, we need to fileter on the desired taxonomic level

species = @chain profilesdf begin
    filter(:taxlevel => ==(:species), _)
    select(Not(:taxlevel))
    unique([:taxname, :sample])
    unstack(:sample, :abundance)
end

genus = @chain profilesdf begin
    filter(:taxlevel => ==(:genus), _)
    select(Not(:taxlevel))
    unique([:taxname, :sample])
    unstack(:sample, :abundance)
end

family = @chain profilesdf begin
    filter(:taxlevel => ==(:family), _)
    select(Not(:taxlevel))
    unique([:taxname, :sample])
    unstack(:sample, :abundance)
end

species[!, r"^G"] .= coalesce.(species[!, r"^G"] ./ 100, 0.)
genus[!, r"^G"] .= coalesce.(genus[!, r"^G"] ./ 100, 0.)
family[!, r"^G"] .= coalesce.(family[!, r"^G"] ./ 100, 0.)

println.(species.taxname);

#- 

println.(genus.taxname);

#- 

println.(family.taxname);

# Many of these taxon names are not _actually_ classified to the level suggested.
# eg `Desulfovibrionaceae_unclassified` is listed in genera,
# but only has a family classifier.
# Here, we'll sum up a bunch of these
# and put them in a single "UNCLASSIFIED" row.

unclassified_patterns = [
    r"^Candidatus",
    r"_unclassified",
    r"[Ii]ncertae_[Ss]edis",
    r"_sp_",
    r"_CAG_",
    r"\d+"
]

@chain species begin
    filter(row-> any(p-> occursin(p, row.taxname), unclassified_patterns), _)
    DataFrame("taxname"=>["UNCLASSIFIED"], (n => sum(_[!,n]) for n in names(_, r"^G"))...)
    append!(species, _)
    filter!(row-> !any(p-> occursin(p, row.taxname), unclassified_patterns), _)
end

@assert all(col-> isapprox(sum(col), 1, atol=1e-5), eachcol(species[!, r"^G"]))

@chain genus begin
    filter(row-> any(p-> occursin(p, row.taxname), unclassified_patterns), _)
    DataFrame("taxname"=>["UNCLASSIFIED"], (n => sum(_[!,n]) for n in names(_, r"^G"))...)
    append!(genus, _)
    filter!(row-> !any(p-> occursin(p, row.taxname), unclassified_patterns), _)
end

@assert all(col-> isapprox(sum(col), 1, atol=1e-5), eachcol(genus[!, r"^G"]))


@chain family begin
    filter(row-> any(p-> occursin(p, row.taxname), unclassified_patterns), _)
    DataFrame("taxname"=>["UNCLASSIFIED"], (n => sum(_[!,n]) for n in names(_, r"^G"))...)
    append!(family, _)
    filter!(row-> !any(p-> occursin(p, row.taxname), unclassified_patterns), _)
end

@assert all(col-> isapprox(sum(col), 1, atol=1e-5), eachcol(family[!, r"^G"]))
#-


CSV.write("paper-diabimmune-method-comparison/diabimmune_mgx_species.csv", species)
CSV.write("paper-diabimmune-method-comparison/diabimmune_mgx_genus.csv", genus)
CSV.write("paper-diabimmune-method-comparison/diabimmune_mgx_family.csv", family)
