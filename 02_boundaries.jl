# # Step 2 - boundary analysis
#
# Now we can start playing around with boundaries - so lets grab 
# `Spatial Boundaries`

using EcologicalNetworks
using SpatialBoundaries
using SpeciesDistributionToolkit

# ## Step 2.1 - Some wrangling first

# Lets start with making our networks though...

# > metacommuinty[M1, M2, S, generation]

# have matrix that is the netwrok matrix for each cell, so 80X80 for the 20x20

N = zeros(Bool, (species_richness, species_richness))
interaction_networks = [copy(N) for x in 1:first(landscape_size), y in 1:last(landscape_size)]

# make interaction netwrok based on interaction strength as well as abundance
# (this is binary)

for i in axes(interaction_strength, 1)
    for j in axes(interaction_strength, 2)
        if !iszero(interaction_strength[i,j])
            for x in axes(interaction_networks, 1)
                for y in axes(interaction_networks, 2)
                    if (metacommunity[x, y, i, end] > 0) && (metacommunity[x, y, j, end] > 0)
                        interaction_networks[x, y][i, j] = 1
                    end
                end
            end
        end
    end
end

# lets calculate the connectance for each landscape patch

# new matrix to store netwrok measure

network_measure = fill(0.0, (landscape_size))

for x in axes(interaction_networks, 1)
    for y in axes(interaction_networks, 2)
        N = UnipartiteNetwork(interaction_networks[x,y])
        network_measure[x,y] = connectance(simplify(N))
    end
end

# species richness

metacommunity2 = metacommunity[:,:,:,end]
metacommunity2[findall(metacommunity[:,:,:,end] .> 0.0), 1] .= 1

species_richness = fill(0, (landscape_size))

for x in axes(interaction_networks, 1)
    for y in axes(interaction_networks, 2)
        species_richness[x,y] = sum(vec(metacommunity2[x,y,:]))
    end
end

# ## Step 2.2 - Boundary time

# lets concatinate the three variable to make looping through them easier

L = [environment_heating[:, :, end], species_richness, network_measure]
rates = fill(zeros(Float16, (2,2)), length(L))
directions = fill(zeros(Float16, (2,2)), length(L))
candidate_boundaries = fill(zeros(Float16, (2,2)), length(L))

for i in 1:3
    wombled_layers = wombling(L[i])
    rates[i] = wombled_layers.m
    directions[i] = wombled_layers.θ

    # calculate boundaries
    rate, direction = SimpleSDMPredictor(wombled_layers)
    b = similar(rate)
    b.grid[boundaries(wombled_layers, 0.1; ignorezero = true)] .= 1.0
    candidate_boundaries[i] = b.grid

    # only keep upper (0.95) and lower (0.05) quantiles
    quants = quantile(rates[i], [0.05, 0.90])
    # we will replace non quantile values with zero
    #rates[i][findall(t ->  quants[1] < t < quants[2], rates[i])] .= 0
    rates[i][findall(t ->  t < quants[2], rates[i])] .= 0
   
end

# visuals

fig = Figure()
axs = [
    Axis(fig[1, 1];
        title = "Environment"),
    Axis(fig[1, 2];
        title = "Richness"),
    Axis(fig[1, 3];
        title = "Connectance"),
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

heatmap!(axs[1], environment_heating[:, :, end])
heatmap!(axs[2], species_richness)
heatmap!(axs[3], network_measure)
heatmap!(axs[4], rates[1], colormap=[:transparent, :green], colorrange=(minimum(rates[1]), maximum(rates[1])))
heatmap!(axs[5], rates[2], colormap=[:transparent, :green], colorrange=(minimum(rates[2]), maximum(rates[2])))
heatmap!(axs[6], rates[3], colormap=[:transparent, :green], colorrange=(minimum(rates[3]), maximum(rates[3])))
heatmap!(axs[7], directions[1], colormap=:romaO, colorrange=(0., 360.))
heatmap!(axs[8], directions[2], colormap=:romaO, colorrange=(0., 360.))
heatmap!(axs[9], directions[3], colormap=:romaO, colorrange=(0., 360.))

current_figure()
save("figures/heatmaps.png", fig)