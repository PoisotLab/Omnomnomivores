# # Step 2 - boundary analysis
#
# Now we can start playing around with boundaries - so lets grab 
# `Spatial Boundaries`

using EcologicalNetworks
using SpatialBoundaries
using SpeciesDistributionToolkit

# > metacommuinty[M1, M2, S, generation, connectivity]

# ## Step 2 - Create some Networks

# have matrix that is the netwrok matrix for each cell, so 80X80 for the 20x20

N = zeros(Bool, (species_richness, species_richness))
interaction_networks = [copy(N) for x in 1:first(landscape_size), y in 1:last(landscape_size)]

# lets make this one landscape network matrix into an array for each level of
# landscpae connectivity

landscape_networks = [deepcopy(interaction_networks) for x in eachindex(c)]

# make interaction netwrok based on interaction strength as well as abundance
# (this is binary)

for k in eachindex(c)
    for i in axes(interaction_strength, 1)
        for j in axes(interaction_strength, 2)
            if !iszero(interaction_strength[i,j])
                for x in axes(interaction_networks, 1)
                    for y in axes(interaction_networks, 2)
                        if (metacommunity[x, y, i, end, k] > 0) && (metacommunity[x, y, j, end, k] > 0)
                            landscape_networks[k][x, y][i, j] = 1
                        end
                    end
                end
            end
        end
    end
end

# lets calculate the connectance for each landscape patch as a 'network measure'

# new matrix to store netwrok measure
network_measure = fill(zeros((landscape_size)), length(c))

for k in eachindex(c)
    for x in axes(interaction_networks, 1)
        for y in axes(interaction_networks, 2)
            N = UnipartiteNetwork(landscape_networks[k][x,y])
            network_measure[k][x,y] = connectance(simplify(N))
        end
    end
end


# ## Step 1 - Boundary time

# make some 'container' matrices

rates = fill(zeros((2,2)), length(c))
directions = fill(zeros((2,2)), length(c))
candidate_boundaries = fill(zeros((2,2)), length(c))

for i in eachindex(c)
    wombled_layers = wombling(network_measure[i])
    rates[i] = wombled_layers.m
    directions[i] = wombled_layers.θ
end

# visuals

fig = Figure()
axs = [
    Axis(fig[1, 1];
        title = ""),
    Axis(fig[1, 2];
        title = "Environment"),
    Axis(fig[1, 3];
        title = ""),
    Axis(fig[2, 1];
        title = ""),
    Axis(fig[2, 2];
        title = "Rate of Change"),
    Axis(fig[2, 3];
        title = ""),
    Axis(fig[3, 1];
        title = ""),
    Axis(fig[3, 2];
        title = "Direction of Change"),
    Axis(fig[3, 3];
        title = ""),
]
#colsize!(fig.layout, 1, Aspect(1, 1))

heatmap!(axs[1], landscape_connectivity[:, :, 1])
heatmap!(axs[2], landscape_connectivity[:, :, 2])
heatmap!(axs[3], landscape_connectivity[:, :, 3])
heatmap!(axs[4], rates[1])
heatmap!(axs[5], rates[2])
heatmap!(axs[6], rates[3])
heatmap!(axs[7], directions[1], colormap=:romaO, colorrange=(0., 360.))
heatmap!(axs[8], directions[2], colormap=:romaO, colorrange=(0., 360.))
heatmap!(axs[9], directions[3], colormap=:romaO, colorrange=(0., 360.))

current_figure()
save("figures/heatmaps.png", fig)