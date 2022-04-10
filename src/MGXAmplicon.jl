module MGXAmplicon

export
    norm_taxon


const taxlevels = (
    d = :domain,
    k = :kingodom,
    p = :phylum,
    c = :class,
    o = :order,
    f = :family,
    g = :genus,
    s = :species,
    t = :subspecies
)

# Many taxon names are not _actually_ classified to the level suggested.
# eg `Desulfovibrionaceae_unclassified` is listed in genera,
# but only has a family classifier.

unclassified_patterns = (
    r"^Candidatus",
    r"unclassified",
    r"uncultured",
    r"[Ii]ncertae_[Ss]edis",
    r"Clade_I+",
    r"_sp_",
    r"_CAG_",
    r"\d+",
    r"_bacterium$",
    r"[Uu]nknown",
    r"\-",
    r"Lineage",
    r"g__\w+aceae"
)

special_cases = Dict(
    "Eubacterium_coprostanoligenes_group"=> 
        "d__Bacteria;p__Firmucutes;c__Clostridia;o__Eubacteriales;f__Eubacteriaceae;g__Eubacterium;s__Eubacterium_coprostanoligenes",
    "f__Chloroplast"=>
        "d__Bacteria;p__Cyanobacteria;c__Cyanobacteriia;o__Chloroplast",
    "Mitochondria"=>
        "d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rickettsiales",
    "f__Peptostreptococcales-Tissierellales"=>
        "d__Bacteria;p__Firmicutes;c__Clostridia;o__Peptostreptococcales-Tissierellales",
    "Clostridium_methylpentosum_group"=>
        "d__Bacteria;p__Firmucutes;c__Clostridia;o__Eubacteriales;f__Oscillospiraceae;g__Clostridium;s__Clostridium_methylpentosum",
    "d__Eukyota"=>
        "d__Eukyota",
    "Alcanivoracaceae1"=>
        "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Oceanospirillales;f__Alcanivoracaceae;g__Alcanivorax",
    "Faecalitalea_sp."=>
        "d__Bacteria; p__Firmicutes; c__Bacilli; o__Erysipelotrichales; f__Erysipelotrichaceae;g__Faecalitalea",
)

function norm_taxon(tax::AbstractString; splitter='|', kind=:mgx)
    taxsplit = split(replace(tax, r"[\[\]]"=> ""), splitter)
    if length(taxsplit) == 1 && any(==(lowercase(tax)), ("unassigned", "unclassified", "unknown")) 
        return (taxlevel=:UNKNOWN, taxname=tax)
    end
    kind == :amplicon ? norm_amp_taxon(taxsplit) : norm_taxon(taxsplit)
end

function norm_taxon(taxvec::AbstractVector)
    isempty(taxvec) && throw(ArgumentError("Something went wrong parsing taxon: Vector is empty!"))
    tax = strip(pop!(taxvec))
    
    taxletter = Symbol(first(tax))
    !in(taxletter, keys(taxlevels)) && return norm_taxon(taxvec)
    
    taxlevel = taxlevels[taxletter]
    taxname = replace(tax, r"[kpcofgs]__"=> "")
    return any(p-> occursin(p, taxname), unclassified_patterns) ?
                norm_taxon(taxvec) :
                (; taxlevel, taxname)
end

function norm_amp_taxon(taxvec::AbstractVector)
    taxvec = strip.(taxvec)
    taxstring = join(taxvec, ";")
    tax = last(taxvec)
    any(p-> occursin(p, replace(tax, r"^[dkpcofgs]__"=>"")), unclassified_patterns) && return norm_amp_taxon(taxvec[1:end-1])
    
    for (key, value) in special_cases
        (key == taxstring || taxstring == value) && continue
        occursin(key, taxstring) && return norm_amp_taxon(split(value, ';'))
    end


    if length(taxvec) >= 7
        sp = taxvec[7]
        # species that end with `sp` or don't start with a capital letter
        # are not acutally classified to the species level
        if any(pat-> occursin(pat, sp), [r"_sp.$", r"^s__[^A-Z]"]) || endswith(sp, "_x")
            # recurse up to genus
            return norm_amp_taxon(taxvec[1:6])
        end
    end
    
    if length(taxvec) >= 6
        ge = taxvec[6]
        ge = replace(ge, r"\[|\]"=>"")
        !occursin(r"^g__[A-Z]", ge) && return norm_amp_taxon(taxvec[1:5])
        if occursin(r"_group$", ge) || occursin(r"_sensu_stricto", ge)
            mg = match(r"g__([A-Z][a-z]+_[a-z]+)_group$", ge)
            mss = match(r"g__([A-Z][a-z]+)_sensu_stricto", ge)
            if !isnothing(mg)
                taxvec[6] = replace(ge, r"_[a-z]+_group$"=>"")
                push!(taxvec, "s__$(mg.captures[1])")
                return norm_amp_taxon(taxvec)
            elseif !isnothing(mss)
                taxvec[6] = "g__$(mss.captures[1])"
                return norm_amp_taxon(taxvec)
            else
                ge = replace(ge, r"_[a-z]+_group$"=>"")
            end
        end
        taxvec[6] = ge
    end
    
    if length(taxvec) >= 5
        fa = taxvec[5]
        if fa == "f__Coriobacteriales_Incertae_Sedis"
            taxvec[5] = "Eggerthellaceae"
            taxvec[4] = "Eggerthellales"
            return norm_amp_taxon(taxvec)
        end

        fa = first(split(replace(fa, r"^f__"=>""), "_"))
        if !endswith(fa, "aceae") || fa == replace(taxvec[4], "^o__"=> "")
            return norm_amp_taxon(taxvec[1:4])
        else
            taxvec[5] = "f__" * fa
        end
    end

    taxlevel = taxlevels[Symbol(first(last(taxvec)))]
    taxname = join(taxvec, ";")
    
    return (; taxlevel, taxname)
end


end # module MGXAmplicon