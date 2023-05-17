using Plots

heatmap(environment_value)

heatmap(meta_comm[:, :, 1, 100])
heatmap(meta_comm[:, :, 5, 1])
heatmap(meta_comm[:, :, 7, 1])

anim = @animate for i=1:100
    heatmap(meta_comm[:, :, 1, i])
end
gif(anim, "figures/herbivire_test.gif", fps = 15)