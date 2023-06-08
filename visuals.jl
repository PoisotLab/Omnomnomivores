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
    Axis(fig[2, 1],
    xlabel = "Environment value",
    ylabel = "Abundance"),
    Axis(fig[2, 2],
    xlabel = "Generation",
    ylabel = "Abundance"),
]
for species in axes(metacommunity, 3)
    tl = trophic_level[species]
    scatter!(axs[tl], vec(environment), vec(metacommunity[:, :, species, end]))
    ylims!(axs[tl], -0.001, 0.001)
end

abund = dropdims(mapslices(sum, metacommunity; dims = (1, 2)); dims = (1, 2))
for species in axes(abund, 1)
    lines!(axs[end], abund[species, 1:end], color = species_col[species])
    ylims!(axs[end], -0.001, 0.01)
end

current_figure()

