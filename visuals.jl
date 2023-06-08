using GLMakie
using ColorSchemes

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
    ylims!(axs[tl], -0.001, 0.001)
end

abund = dropdims(mapslices(sum, metacommunity; dims = (1, 2)); dims = (1, 2))
for species in axes(abund, 1)
    lines!(axs[4], abund[species, 1:end], color = species_col[species])
    ylims!(axs[4], -0.001, 0.01)
end

abund[findall(abund .> 0.0), 1] .= 1.0
lines!(axs[5], vec(sum(abund[1:40,:], dims = 1)), color = "green", label = "plant")
lines!(axs[5], vec(sum(abund[41:63,:], dims = 1)), color = "blue", label = "herbivore")
lines!(axs[5], vec(sum(abund[64:end,:], dims = 1)), color = "red", label = "carnivore")
axislegend()

current_figure()

save("figures/diagnostics.png", fig)