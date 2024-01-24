# # Step 2 - boundary analysis
#
# Now we can start playing around with boundaries - so lets grab 
# `Spatial Boundaries`

using EcologicalNetworks
using SpatialBoundaries
using SpeciesDistributionToolkit

# > metacommuinty[M1, M2, S, generation, connectivity]

# ## Step 2 - Create some Networks

# have matrix that is the network matrix for each cell, so 80X80 for the 20x20

N = zeros(Bool, (species_richness, species_richness))
interaction_networks = [deepcopy(N) for x in 1:first(landscape_size), y in 1:last(landscape_size)]

# lets make this one landscape network matrix into an array for each level of
# landscpae connectivity

landscape_networks = [deepcopy(interaction_networks) for x in eachindex(c)]

# make interaction network based on interaction strength as well as abundance
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

# species richness

metacommunity2 = metacommunity[:,:,:,end,:]
metacommunity2[findall(metacommunity[:,:,:,end,:] .> 0.0), 1] .= 1

richness_landscape = fill(0.0, (landscape_size..., length(c)))

for k in eachindex(c) 
    for x in axes(interaction_networks, 1)
        for y in axes(interaction_networks, 2)
            richness_landscape[x,y,k] = sum(vec(metacommunity2[x,y,:,k]))
        end
    end
end

# lets calculate the connectance for each landscape patch as a 'network measure'

# new matrix to store network measure
network_measure = fill(0.0, (landscape_size..., length(c)))

for k in eachindex(c)
    for x in axes(interaction_networks, 1)
        for y in axes(interaction_networks, 2)
            N = UnipartiteNetwork(landscape_networks[k][x,y])
            network_measure[x,y,k] = connectance(simplify(N))
        end
    end
end

# ## Step 1 - Boundary time

# make some 'container' matrices

L = [environment_heating[:, :, end, :], richness_landscape, network_measure]
rates = fill(0.0, (landscape_size .- 1)..., length(L), length(c))
directions = fill(0.0, (landscape_size .- 1)..., length(L), length(c))
candidate_boundaries = fill(0.0, (landscape_size .- 1)..., length(L), length(c))

for l in eachindex(L)
    for k in eachindex(c)
        wombled_layers = wombling(L[l][:,:,k])
        rates[:,:,l,k] = wombled_layers.m
        directions[:,:,l,k] = wombled_layers.θ
    end
end

# visuals

fig = Figure()
axs = [
    Axis(fig[1, 1]),
    Axis(fig[1, 2];
        title = "Rate of Change (Environment"),
    Axis(fig[1, 3]),
    Axis(fig[2, 1]),
    Axis(fig[2, 2];
        title = "Rate of Change (Richness)"),
    Axis(fig[2, 3]),
    Axis(fig[3, 1];
        xlabel = "low connectivity"),
    Axis(fig[3, 2];
        title = "Rate of Change (Connectence)",
        xlabel = "mid connectivity"),
    Axis(fig[3, 3];
        xlabel = "high connectivity"),
]
#colsize!(fig.layout, 1, Aspect(1, 1))

heatmap!(axs[1], rates[:,:,1,1])
heatmap!(axs[2], rates[:,:,1,2])
heatmap!(axs[3], rates[:,:,1,3])
heatmap!(axs[4], rates[:,:,2,1])
heatmap!(axs[5], rates[:,:,2,2])
heatmap!(axs[6], rates[:,:,2,3])
heatmap!(axs[7], rates[:,:,3,1])
heatmap!(axs[8], rates[:,:,3,2])
heatmap!(axs[9], rates[:,:,3,3])

current_figure()
save("figures/heatmaps.png", fig)