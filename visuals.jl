using Plots

heatmap(environment_value)

heatmap(meta_comm[:, :, 1, 100])
heatmap(meta_comm[:, :, 6, 1])
heatmap(meta_comm[:, :, 7, 1])

anim = @animate for i=1:100
   plot(heatmap(meta_comm[:, :, 2, i], title="Plant"),
   heatmap(meta_comm[:, :, 6, i], title="Herbivore"),
   heatmap(meta_comm[:, :, 7, i], title="Carnivore"),
   layout=(3,1))
end
gif(anim, "figures/all_test.gif", fps = 1)