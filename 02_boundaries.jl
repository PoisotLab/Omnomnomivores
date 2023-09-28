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

interaction_networks =
    fill(zeros(Bool, (species_richness, species_richness)), (landscape_size))

# make interaction netwrok based on interaction strength as well as abundance
# (this is binary)

for x in axes(interaction_networks, 1)
    for y in axes(interaction_networks, 2)
        for i in axes(interaction_strength, 1)
            for j in axes(interaction_strength, 2)
                if interaction_strength[i, j] !=0 && metacommunity[x, y, i, end] > 0 &&
                   metacommunity[x, y, j, end] > 0
                    interaction_networks[x, y][i, j] = 1
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
candidate_boundaries = fill(zeros(Float16, (2,2)), length(L))

for i in 1:3
    wombled_layers = wombling(L[i])
    rate, direction = SimpleSDMPredictor(wombled_layers)
    b = similar(rate)
    b.grid[boundaries(wombled_layers, 0.1; ignorezero = true)] .= 1.0
    candidate_boundaries[i] = b.grid
    
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
        title = "Environental Boundaries"),
    Axis(fig[2, 2];
        title = "Richness Boundaries"),
    Axis(fig[2, 3];
        title = "Connectance Boundaries"),
]
#colsize!(fig.layout, 1, Aspect(1, 1))

heatmap!(axs[1], environment_heating[:, :, end])
heatmap!(axs[2], species_richness)
heatmap!(axs[3], network_measure)
heatmap!(axs[4], candidate_boundaries[1], colormap=[:transparent, :green])
heatmap!(axs[5], candidate_boundaries[2], colormap=[:transparent, :green])
heatmap!(axs[6], candidate_boundaries[3], colormap=[:transparent, :green])

current_figure()
save("figures/heatmaps.png", fig)