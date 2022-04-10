# ---
# author:
# - Kevin Bonham
# - PhD
# title: DIABIMMUNE re-run
# ---

# # DIABIMMUNE re-run

# ## 16S files

# ### Downloading files

# [DIABIMMUNE](https://diabimmune.broadinstitute.org/diabimmune/three-country-cohort/resources/16s-sequence-data)

# 1.  Command from website doesn\'t work
#     `wget -r -np -nd https://pubs.broadinstitute.org/diabimmune/data/10`
#     1.  This downloads a file called `10` that is actually html
#     2.  The html file has urls inside it
#     3.  I got those out with a series of `sed` commands
        
#         ``` {.bash org-language="sh"}
#         sed -E "s/<a href='([^']+)'/\1\n/g" 10 > urls.txt
#         sed -E 's/^.+(https.+)/\1/' urls.txt| grep https > urls1.txt
#         mv urls1.txt urls.txt
#         ```

#         Trying to download from these urls, it says they\'ve moved, so
#         need to change `pubs.broadinstitute.org/diabimmune` to
#         `diabimmune.broadinsitute.org/diabimmune`

#         ``` {.bash org-language="sh"}
#         sed -i 's/pubs/diabimmune/' urls.txt
#         ```

#     4.  Then downloaded in parallel with `xargs`

#         ``` {.bash org-language="sh"}
#         cat urls.txt | xargs -n1 -P8 curl -LO
#         ```

# ### Running the pipeline

# Using QIIME.jl

# ```sh
#  qiimejl --threads 16 -v basic -n pipeline -i rawfastq -o ./ --input-format CasavaOneEightLanelessPerSampleDirFmt --fwd-primer GTGCCAGCMGCCGCGGTAA --rev-primer GGACTACHVGGGTWTCTAAT --fwd-trunc 150 --rev-trunc 150
# ```

# ### Massage dataframe

using CSV
using DataFrames
using Chain
using MGXAmplicon: norm_taxon


features = CSV.read("/augusta/staff/kevin/echo/diabimmune_amplicon/pipeline/dada2/feature-table.tsv", DataFrame, header=2, delim='\t')
taxa = CSV.read("/augusta/staff/kevin/echo/diabimmune_amplicon/pipeline/taxa/taxonomy.tsv", DataFrame, delim='\t')

rename!(features, "#OTU ID"=>"feature")
rename!(taxa, "Feature ID"=>"feature")

labeled = leftjoin(select(taxa, ["feature", "Taxon"]), features, on="feature")
# some columns don't have any taxonomic assigment, remove them
labeled = labeled[!, [!(eltype(labeled[!, col]) <: Union{Missing, <:Number}) || 
                        sum(labeled[!, col]) > 0 
                        for col in names(labeled)]]

# convert to relative abundance
select!(labeled, [:Taxon => (x -> identity(x)), 
                 [C => c -> c ./ sum(c) for C in names(labeled, r"^G")]...],
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
    combine([n => sum for n in names(labeled, r"^G")], renamecols=false)
end

ge = @chain labeled begin
    groupby(:genus)
    combine([n => sum for n in names(labeled, r"^G")], renamecols=false)
end


sp = @chain labeled begin
    groupby(:species)
    combine([n => sum for n in names(labeled, r"^G")], renamecols=false)
end


CSV.write("/home/kevin/repos/danielle-thesis/diabimmune/karalia_dada2_family.csv", fa)
CSV.write("/home/kevin/repos/danielle-thesis/diabimmune/karalia_dada2_genera.csv", ge)
CSV.write("/home/kevin/repos/danielle-thesis/diabimmune/karalia_dada2_species.csv", sp)

# ### Get some statistics

using Statistics

sp_totals = sum.(eachcol(filter(:species=> !=("UNCLASSIFIED"), sp)[!, r"^G"]));
mean(sp_totals)
median(sp_totals)

ge_totals = sum.(eachcol(filter(:genus => !=("UNCLASSIFIED"), ge)[!, r"^G"]));
mean(ge_totals)
median(ge_totals)

# ## Metagenomes


# ### Similar to above, need to parse download

# ``` {.shell}
# wget -r -np -nd https://pubs.broadinstitute.org/diabimmune/data/16
# sed -E "s/<a href='([^']+)'/\1\n/g" 16 > urls.txt
# sed -E 's/^.+(https.+)/\1/' urls.txt| grep https > urls1.txt
# mv urls1.txt urls.txt
# cat urls.txt | xargs -n1 -P8 curl -LO
# ```

# ### Environment setup
#
# #### Load singularity module
#
# ``` shell
# module add singularity/3.7.0
# ```

# #### Install biobakery workflows container

# ``` bash
# singularity build biobakery.simg docker://biobakery/workflows
# ```

# #### Bind paths for access

# - By default, singularity only has access to `$HOME`, `/tmp`, and `$PWD`
# - Can add paths via `--bind /path/on/syestem:/path/in/container`
# - Can also add to `$SINGULARITY_BINDPATH` env variable

#   eg.
  
#   ``` bash
#   export SINGULARITY_BINDPATH=/pool001/vklepacc:/pool,/nobackup1/vklepacc:/nobackup
#   export BIOBAKERY_WORKFLOWS_DATABASES=/pool/databases
#   ```

# #### Install databases

# ``` shell
# singularity exec ~/software/biobakery/biobakery.simg kneaddata_database --download human_genome bowtie2 /pool/databases/kneaddata
# singularity exec ~/software/biobakery/biobakery.simg kneaddata_database --download ribosomal_RNA bowtie2 /pool/databases/kneaddata
# singularity exec ~/software/biobakery/biobakery.simg humann_databases --download chocophlan full /pool/databases/humann
# singularity exec ~/software/biobakery/biobakery.simg humann_databases --download uniref uniref90_diamond /pool/databases/humann
# singularity exec ~/software/biobakery/biobakery.simg humann_databases --download utility_mapping full /pool/databases/humann
# ```

# NOTE: metaphlan database is already installed. Others are too big

# Singularity is unable to update the `humann_config`, so these options need to be manually included in the workflow run

# #### Run the workflows

# ``` shell

# ```


