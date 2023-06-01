using CairoMakie
using GLMakie
using Makie.Colors
using Plots

# get global extrema
extremas = map(extrema, metacommunity)
global_min = minimum(t->first(t), extremas)
global_max = maximum(t->last(t), extremas)
# these limits have to be shared by the maps and the colorbar
clims = (global_min, global_max)

anim = @animate for i=1:generations
   Plots.plot(Plots.heatmap(environment, title="Environment"),
   Plots.heatmap(metacommunity[:, :, 1, i], title="Plant", clims = clims),
   Plots.heatmap(metacommunity[:, :, 55, i], title="Herbivore", clims = clims),
   Plots.heatmap(metacommunity[:, :, 78, i], title="Carnivore", clims = clims),
   layout=(2,2))
end
gif(anim, "figures/all_test.gif", fps = 1)