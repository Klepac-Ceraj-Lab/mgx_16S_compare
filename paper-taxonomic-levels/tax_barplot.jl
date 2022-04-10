using CairoMakie
using Clustering
using Distances
using CSV
using DataFrames


u1dm = braycurtis(u1_spec)
middm = braycurtis(mid_spec)
o2dm = braycurtis(o2_spec)

u1clust = hclust(u1dm, linkage=:complete, branchorder=:optimal)
midclust = hclust(middm, linkage=:complete, branchorder=:optimal)
o2clust = hclust(o2dm, linkage=:complete, branchorder=:optimal)


function topx(cp, n=10)
    totals = vec(featuretotals(cp))
    rows = partialsortperm(totals, 1:n, rev=true)
    top = cp[rows, :]
    other = sum(abundances(cp[collect(1:nfeatures(cp))[Not(rows)], :]), dims=1)

    return (names = vcat(featurenames(top), ["other"]), abundances = vcat(abundances(top), other))
end

allnames = setdiff(union([topx(cp).names for cp in (u1_spec, mid_spec, o2_spec)]...), Set(["other"]))
name_dict = Dict(n => colormap[i] for (i, n) in enumerate(allnames))
name_dict["other"] = ColorSchemes.Greys_3.colors[1]

function plottopn!(fig, layout, axrow, cp, clust, n, title, markerspecies = String[]; kwargs...)
    topnames, topabund = topx(cp, n)
    topabund = topabund[:, clust.order]
    ax = layout[axrow,1] = Axis(fig; title)
    append!(markerspecies, [topnames[i] for i in 1:n+1])
    
    for i in (n+1):-1:1
        v = vec(sum(topabund[1:i, :], dims=1))
        barplot!(ax, 1:nsamples(cp), v, color=name_dict[topnames[i]])
    end

    tightlimits!(ax)
    hidexdecorations!(ax)
    fig, markerspecies
end

sublayout = GridLayout()
figure1[1:2, 4:6] = sublayout
(_, markerspecies) = plottopn!(figure1, sublayout, 1, u1_spec, u1clust, 10, "Top 10 species, kids under 1 yo")
plottopn!(figure1, sublayout, 2, mid_spec, midclust, 10, "Top 10 species, kids 1-2 yo", markerspecies)
plottopn!(figure1, sublayout, 3, o2_spec, o2clust, 10, "Top 10 species, kids over 2 yo", markerspecies)
unique!(markerspecies)
sort!(markerspecies)
fig1c_leg = sublayout[4, 1] = Legend(figure1, [MarkerElement(color = name_dict[m], marker = :rect, strokecolor = :black) for m in markerspecies], [m for m in markerspecies], tellheight=true, nbanks = 2)
