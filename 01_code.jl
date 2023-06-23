using Distributions
using NeutralLandscapes
using SpatialBoundaries
using Random
using ProgressMeter

landscape_size = (20, 20)
species_richness = 80

# Load the functions we need from the lib folder
include("lib/01_species_creation.jl")
include("lib/02_model_internals.jl")

# Prepare burn in run
generations = 200

interaction_strength = zeros(Float64, (species_richness, species_richness))
trophic_level = zeros(Int8, species_richness)
environmental_optimum = zeros(Float64, species_richness)
dispersal_decay = zeros(Float64, species_richness)
dispersal_rate = zeros(Float64, species_richness)

environment_burnin = zeros(Float64, landscape_size)
environment_burnin .= 10.0

# Initial values for the species
set_trophic_levels!(trophic_level)
set_interaction_strength!(interaction_strength; trophic_level)
set_environmental_optimum!(environmental_optimum, environment_burnin, trophic_level)
set_dispersal_rate!(dispersal_rate)
set_dispersal_decay!(dispersal_decay; trophic_level)

# Set an initial metacommunity object, only the first timestep is set to 10.0
metacommunity_burnin = fill(0.0, (landscape_size..., species_richness, generations))
metacommunity_burnin[:, :, :, 1] .= 0.1

# Run the model
simulate!(
    metacommunity_burnin,
    dispersal_rate,
    dispersal_decay,
    environment_burnin,
    environmental_optimum,
    interaction_strength,
)

# set up for 'heating' run

# set heating generation timestep
generations_heating = 200
generations = 2

# create 'ramped' environmental change landscape
environment_heating = fill(0.0, (landscape_size..., generations_heating))
environment_heating[:,:,1] = environment_burnin
environment_heating[:,:,end] = rand(DiamondSquare(), landscape_size) .* species_richness
# final - initial environment
heating_step = (environment_heating[:,:,end] - environment_burnin)./generations_heating
# add said quotient to each patch sequentially
for e in 2:(generations_heating - 1)
    environment_heating[:,:,e] = environment_heating[:,:,e-1] + heating_step
end

# re-assign new environmental optimum
set_environmental_optimum!(environmental_optimum, environment_heating[:,:,end], trophic_level)

# create new metacommunity matrix
metacommunity = fill(0.0, (landscape_size..., species_richness, generations_heating))
# assign T1 as final community from burnin
metacommunity[:, :, :, 1] .= metacommunity_burnin[:, :, :, end]

# run simulations

for g in 1:(generations_heating - 1)
_meta_comm = metacommunity[:, :, :, g:(g + 1)]
_eopt = set_environmental_optimum!(environmental_optimum, environment_heating[:,:,g+1], trophic_level)
    simulate!(
    _meta_comm,
    dispersal_rate,
    dispersal_decay,
    environment_heating[:,:,g+1],
    _eopt,
    interaction_strength,
    )
metacommunity[:, :, :, g:(g + 1)] = _meta_comm
end