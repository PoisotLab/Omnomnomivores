using CairoMakie
using GLMakie
using Makie.Colors
using Plots

# get global extrema
extremas = map(extrema, meta_comm)
global_min = minimum(t->first(t), extremas)
global_max = maximum(t->last(t), extremas)
# these limits have to be shared by the maps and the colorbar
clims = (global_min, global_max)

species_col = fill("", species_richness,)

for i in axes(species_col,1)
   if trophic_level[i] == 1
      species_col[i] = "green"
   elseif trophic_level[i] == 2
      species_col[i] = "blue"
   else 
      species_col[i] = "red"
   end
end

# Some diagnostic plots

#burnin community
fig = Figure()
axs = [
    Axis(fig[1, 1],
    xlabel = "Environment value",
    ylabel = "Abundance"),
    Axis(fig[1, 2],
    xlabel = "Environment value",
    ylabel = "Abundance"),
    Axis(fig[1, 3],
    xlabel = "Environment value",
    ylabel = "Abundance"),
    Axis(fig[2, 1:3],
    xlabel = "Generation",
    ylabel = "Abundance"),
    Axis(fig[3, 1:3],
    xlabel = "Generation",
    ylabel = "Species Richness"),
]
for species in axes(metacommunity_burnin, 3)
    tl = trophic_level[species]
    scatter!(axs[tl], vec(environment_burnin), vec(metacommunity_burnin[:, :, species, end]))
end

abund = dropdims(mapslices(sum, metacommunity_burnin; dims = (1, 2)); dims = (1, 2))
for species in axes(abund, 1)
    lines!(axs[4], abund[species, 1:end], color = species_col[species])
end

abund[findall(abund .> 0.0), 1] .= 1.0
lines!(axs[5], vec(sum(abund[1:40,:], dims = 1)), color = "green", label = "plant")
lines!(axs[5], vec(sum(abund[41:63,:], dims = 1)), color = "blue", label = "herbivore")
lines!(axs[5], vec(sum(abund[64:end,:], dims = 1)), color = "red", label = "carnivore")
axislegend()

current_figure()

save("figures/diagnostics_burnin.png", fig)

# 'heated' community
fig = Figure()
axs = [
    Axis(fig[1, 1],
    xlabel = "Environment value",
    ylabel = "Abundance"),
    Axis(fig[1, 2],
    xlabel = "Environment value",
    ylabel = "Abundance"),
    Axis(fig[1, 3],
    xlabel = "Environment value",
    ylabel = "Abundance"),
    Axis(fig[2, 1:3],
    xlabel = "Generation",
    ylabel = "Abundance"),
    Axis(fig[3, 1:3],
    xlabel = "Generation",
    ylabel = "Species Richness"),
]
for species in axes(metacommunity, 3)
    tl = trophic_level[species]
    scatter!(axs[tl], vec(environment), vec(metacommunity[:, :, species, end]))
end

abund = dropdims(mapslices(sum, metacommunity; dims = (1, 2)); dims = (1, 2))
for species in axes(abund, 1)
    lines!(axs[4], abund[species, 1:end], color = species_col[species])
end

abund[findall(abund .> 0.0), 1] .= 1.0
lines!(axs[5], vec(sum(abund[1:40,:], dims = 1)), color = "green", label = "plant")
lines!(axs[5], vec(sum(abund[41:63,:], dims = 1)), color = "blue", label = "herbivore")
lines!(axs[5], vec(sum(abund[64:end,:], dims = 1)), color = "red", label = "carnivore")
axislegend()

current_figure()

save("figures/diagnostics.png", fig)