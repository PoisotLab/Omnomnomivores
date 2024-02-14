# # Step 1 - creating communities in the environment
#
# This step is actually two steps since we will first be creating a 'burn-in'
# phase where the landscape is uniform for a set number of generations and then
# we will incrementally begin to heat/cool the different environmental patches
# until we reach the 'final state'. 

# ## Dependencies

# Nothing fancy

using Distributions
using CairoMakie
using NeutralLandscapes
using Random
using ProgressMeter

# ## Functionality
#
# Load the functions we need from the lib folder

include("lib/01_species_creation.jl")
include("lib/02_model_internals.jl")
include("lib/03_model_specifications.jl")

# ## Initiation

Random.seed!(66)

# Let's specify parameters for the model

landscape_size = (20, 20)
species_richness = 80
generations = 10
landscape_connectivity = range(start = 0.1; stop = 0.99, length = 5)
landscape = zeros(Float64, (landscape_size..., generations, length(landscape_connectivity)))

# Now we can begin to create the species metadata which we will populate in a
# bit(you can see a more detailed breakdown of these in
# `lib/01_species_creation.jl`)

interaction_strength = zeros(Float64, (species_richness, species_richness))
trophic_level = zeros(Int8, species_richness)
environmental_optimum = zeros(Float64, species_richness, generations, length(landscape_connectivity))
dispersal_decay = zeros(Float64, species_richness)
dispersal_rate = zeros(Float64, species_richness)


# We can populate the landscape

set_landscape!(generations, landscape_connectivity, landscape, landscape_size, species_richness)

# Now we can populate the species metadata

set_trophic_levels!(trophic_level)
set_interaction_strength!(interaction_strength; trophic_level)
set_dispersal_rate!(dispersal_rate)
set_dispersal_decay!(dispersal_decay; trophic_level)

# environmental optimum is a bit more 'tricksy' since we have different landscapes
# so we can just loop through the respective landscapes

# 🐛 I hate this but it will have to do for now... I think we can fix this by
# creating an alternive function state for environmental_optimum!

for i in 1:generations
    for k in eachindex(landscape_connectivity)
        environmental_optimum[:, i, k] = 
        set_environmental_optimum!(environmental_optimum[:, i, k],
                                landscape[:, :, i, k],
                                trophic_level,
    )
    end
end


# Now we can create the matrix that will store each timedtamp of community.
# We will also create an initial abundance for all species for the first
# timestamp.

metacommunity = fill(0.0, (landscape_size..., species_richness, generations, length(landscape_connectivity)));
metacommunity[:, :, :, 1, :] .= 0.1;
# Now we can simply run the model

for k in eachindex(c)

    #run the simulation
    for g in 1:generations
        _meta_comm = metacommunity[:, :, :, :, k]
        _eopt = environmental_optimum[:, g, k]
        simulate!(
            _meta_comm,
            dispersal_rate,
            dispersal_decay,
            landscape[:, :, g, k],
            _eopt,
            interaction_strength,
        )
        metacommunity[:, :, :, :, k] = _meta_comm
    end   
end

# ### Diagnostics

using CairoMakie
using Makie.Colors

fig = Figure()
axs = [
    Axis(fig[1, 1];
        xlabel = "Environment value",
        ylabel = "Abundance"),
    Axis(fig[1, 2];
        xlabel = "Environment value",
        ylabel = "Abundance"),
    Axis(fig[1, 3];
        xlabel = "Environment value",
        ylabel = "Abundance"),
    Axis(fig[2, 1:3];
        xlabel = "Generation",
        ylabel = "Abundance"),
    Axis(fig[3, 1:3];
        xlabel = "Generation",
        ylabel = "Species Richness"),
]
for species in axes(metacommunity, 3)
    tl = trophic_level[species]
    scatter!(
        axs[tl],
        vec(landscape[:, :, end, 2]),
        vec(metacommunity[:, :, species, end, 1]),
    )
end

abund = dropdims(mapslices(sum, metacommunity[:,:,:,:,2]; dims = (1, 2)); dims = (1, 2))
for species in axes(abund, 1)
    lines!(axs[4], abund[species, 1:end]; color = species_col[species])
end

abund[findall(abund .> 0.0), 1] .= 1.0
lines!(axs[5], vec(sum(abund[1:40, :, 1]; dims = 1)); color = "green", label = "plant")
lines!(axs[5], vec(sum(abund[41:63, :, 1]; dims = 1)); color = "blue", label = "herbivore")
lines!(axs[5], vec(sum(abund[64:end, :, 1]; dims = 1)); color = "red", label = "carnivore")
axislegend()

current_figure()

save("figures/diagnostics.png", fig)

# This is good - we can see different species peaking at different environmental
# optima as well as a few species going extinct. Oh and somewhat stable community
# dynamics I guess...

