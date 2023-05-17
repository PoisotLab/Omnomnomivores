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

anim = @animate for i=1:100
   plot(heatmap(environment_value, title="Environment"),
   heatmap(meta_comm[:, :, 1, i], title="Plant (Eopt 0.25)", clims = (global_min, global_max)),
   heatmap(meta_comm[:, :, 6, i], title="Herbivore (Eopt 1.0)", clims = (global_min, global_max)),
   heatmap(meta_comm[:, :, 7, i], title="Carnivore (Eopt 0.5)", clims = (global_min, global_max)),
   layout=(2,2))
end
gif(anim, "figures/all_test.gif", fps = 1)