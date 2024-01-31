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
using Makie.Colors
using NeutralLandscapes
using Random
using ProgressMeter

# ## Functionality
#
# Load the functions we need from the lib folder

include("lib/01_species_creation.jl")
include("lib/02_model_internals.jl")

# ## Initiation

Random.seed!(66) # the time has come

# First we will specify the size of our landscape as well as the community 

landscape_size = (20, 20)
species_richness = 80

# Now we can begin to create the species metadata which we will populate in a
# bit(you can see a more detailed breakdown of these in
# `lib/01_species_creation.jl`)

interaction_strength = zeros(Float64, (species_richness, species_richness))
trophic_level = zeros(Int8, species_richness)
environmental_optimum = zeros(Float64, species_richness)
dispersal_decay = zeros(Float64, species_richness)
dispersal_rate = zeros(Float64, species_richness)

# ## Burn-in
#
# Let's specify the number of burn-in generations

generations = 200

# For the burn-in we want to keep the landscape uniform so we will populate the
# landscape with an environmental value of 10.0

environment_burnin = zeros(Float64, landscape_size)
environment_burnin .= 10.0

# Now we can populate the species metadata

set_trophic_levels!(trophic_level)
set_interaction_strength!(interaction_strength; trophic_level)
set_environmental_optimum!(environmental_optimum, environment_burnin, trophic_level)
set_dispersal_rate!(dispersal_rate)
set_dispersal_decay!(dispersal_decay; trophic_level)

# Now we can create the matrix that will store each timedtamp of our burn-in
# community. We will also create an initial abundance for all species for the
# first timestamp.

metacommunity_burnin = fill(0.0, (landscape_size..., species_richness, generations))
metacommunity_burnin[:, :, :, 1] .= 0.1

# Now we can simply run the model

simulate!(
    metacommunity_burnin,
    dispersal_rate,
    dispersal_decay,
    environment_burnin,
    environmental_optimum,
    interaction_strength,
)

# ### Diagnostics

# We can compile some visuals to see how the community is changing over time.

## this is for some colour allocation
species_col = fill("", species_richness)
for i in axes(species_col, 1)
    if trophic_level[i] == 1
        species_col[i] = "green"
    elseif trophic_level[i] == 2
        species_col[i] = "blue"
    else
        species_col[i] = "red"
    end
end

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
for species in axes(metacommunity_burnin, 3)
    tl = trophic_level[species]
    scatter!(
        axs[tl],
        vec(environment_burnin),
        vec(metacommunity_burnin[:, :, species, end]),
    )
end

abund = dropdims(mapslices(sum, metacommunity_burnin; dims = (1, 2)); dims = (1, 2))
for species in axes(abund, 1)
    lines!(axs[4], abund[species, 1:end]; color = species_col[species])
end

abund[findall(abund .> 0.0), 1] .= 1.0
lines!(axs[5], vec(sum(abund[1:40, :]; dims = 1)); color = "green", label = "plant")
lines!(axs[5], vec(sum(abund[41:63, :]; dims = 1)); color = "blue", label = "herbivore")
lines!(axs[5], vec(sum(abund[64:end, :]; dims = 1)); color = "red", label = "carnivore")
axislegend()

current_figure()

save("figures/diagnostics_burnin.png", fig)

# The good news here is that we can see some sort of population dynamics
# emerging when we look at the abundance curves.

# ## Community generation
#
# Now that we have created a somewhat 'stocahstic' community it is time to
# introduce some environmental variability - we will do this using
# `NeutralLandscapes.jl`. But before we do that lets specify our 'heating'
# generations time. 

# > Note because of how the simulation function is currently set up we will be
# > setting two generation time variables. `generations_heating` is the actual
# > number of generations and 'generations' is simply to make the simulate!
# > function go brrr.

generations_heating = 500
generations = 2

# because we are 'ramping' the environment (i.e. simulating a gradual change) we
# are creating a matrix that will keep each environmental layer so we can sample
# it for each simulation run. 

# but we also want to create landscapes with varying spatial-auto correlation so
# lets also create those

# specify connectivity values
c = [0, 0.5, 0.99]
# create empty landscape matrix
landscape_connectivity = zeros(Float64, (landscape_size..., length(c)))

# populate with environmental values
for i in eachindex(c)
    landscape_connectivity[:, :, i] = rand(DiamondSquare(c[i]), landscape_size)
end

environment_heating = fill(0.0, (landscape_size..., generations_heating))

# TODO #13 The plan here is to have a series of generations with set values
# (proofing); then we do the increase in temperature (baking); and then some
# number of generations with the final values (resting). This should be handled
# by a single function, and we can simply use the baking update schedule to
# manifest the environmental values at each timestep. The switch between
# logisting and linear can be handled by passing a function as an argument
# somewhere.

environment_heating[:, :, 1] = environment_burnin
environment_heating[:, :, end] = rand(DiamondSquare(), landscape_size) .* species_richness

# now we can calculate the magnitude of change for each timestep to have the
# landscape reach its 'final state'

heating_step = (environment_heating[:, :, end] - environment_burnin) ./ generations_heating

# now we can sequentially add this value to each intermediate landscape state.
# To make the environmental change a bit more gradual we can add a logisitc
# 'tweak' to the environmental change.

for e in 2:(generations_heating - 1)
    environment_heating[:, :, e] = (environment_heating[:, :, e - 1] + heating_step)
end

for i in 1:generations_heating
    environment_heating[:, :, i] =
        environment_heating[:, :, i] * (1 / (1 + exp(-(i / generations_heating))))
end

# we also need to reassign the envirnomental optima of all species (recall these
# were optimised to the uniform landscape)

set_environmental_optimum!(
    environmental_optimum,
    environment_heating[:, :, end],
    trophic_level,
)

# create the new metacommuity matrix and assign the first timestep as the
# abundance values of the final timestep of the burnin community.

metacommunity = fill(0.0, (landscape_size..., species_richness, generations_heating))
metacommunity[:, :, :, 1] .= metacommunity_burnin[:, :, :, end]

# run the simulations

for g in 1:(generations_heating - 1)
    _meta_comm = metacommunity[:, :, :, g:(g + 1)]
    _eopt = set_environmental_optimum!(
        environmental_optimum,
        environment_heating[:, :, g + 1],
        trophic_level,
    )
    simulate!(
        _meta_comm,
        dispersal_rate,
        dispersal_decay,
        environment_heating[:, :, g + 1],
        _eopt,
        interaction_strength,
    )
    metacommunity[:, :, :, g:(g + 1)] = _meta_comm
end

# ### Diagnostics
#
# We can repeat the same series of plots we used to look at the burn-in
# community.

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
        vec(environment_heating[:, :, end]),
        vec(metacommunity[:, :, species, end]),
    )
end

abund = dropdims(mapslices(sum, metacommunity; dims = (1, 2)); dims = (1, 2))
for species in axes(abund, 1)
    lines!(axs[4], abund[species, 1:end]; color = species_col[species])
end

abund[findall(abund .> 0.0), 1] .= 1.0
lines!(axs[5], vec(sum(abund[1:40, :]; dims = 1)); color = "green", label = "plant")
lines!(axs[5], vec(sum(abund[41:63, :]; dims = 1)); color = "blue", label = "herbivore")
lines!(axs[5], vec(sum(abund[64:end, :]; dims = 1)); color = "red", label = "carnivore")
axislegend()

current_figure()

save("figures/diagnostics.png", fig)

# This is good - we can see different species peaking at different environmental
# optima as well as a few species going extinct. Oh and somewhat stable community
# dynamics I guess...

