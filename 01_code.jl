using Distributions
using NeutralLandscapes
using SpatialBoundaries
using Random
using ProgressMeter

landscape_size = (20, 20)
species_richness = 80
generations = 200

# Load the functions we need from the lib folder
include("lib/01_species_creation.jl")
include("lib/02_model_internals.jl")

# Prepare the run
interaction_strength = zeros(Float64, (species_richness, species_richness))
trophic_level = zeros(Int8, species_richness)
environmental_optimum = zeros(Float64, species_richness)
dispersal_decay = zeros(Float64, species_richness)
dispersal_rate = zeros(Float64, species_richness)

environment = rand(DiamondSquare(), landscape_size) .* species_richness

# Initial values for the species
set_trophic_levels!(trophic_level)
set_interaction_strength!(interaction_strength; trophic_level)
set_environmental_optimum!(environmental_optimum, environment, trophic_level)
set_dispersal_rate!(dispersal_rate)
set_dispersal_decay!(dispersal_decay; trophic_level)

# Set an initial metacommunity object, only the first timestep for plants is set to 10.0
metacommunity = fill(0.0, (landscape_size..., species_richness, generations))
metacommunity[:, :, findall(trophic_level .< 2), 1] .= 15.0
metacommunity[:, :, findall(trophic_level .== 2), 50] .= 7.0
metacommunity[:, :, findall(trophic_level .== 3), 100] .= 1.0

# Run the model
simulate!(
    metacommunity,
    dispersal_rate,
    dispersal_decay,
    environment,
    environmental_optimum,
    interaction_strength,
)