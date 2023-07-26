# # Step 1 - creating communities in the environment
#
# This step is actually two steps since we will first be creating a 'burn-in'
# phase where the landscape is uniform for a set number of generations and then
# we will incrementally begin to heat/cool the different environmental patches
# until we reach the 'final state'. 

# ## Dependencies

# Nothing fancy

using Distributions
using NeutralLandscapes
using Random
using ProgressMeter

# ## Functionality
#
# Load the functions we need from the lib folder

include("lib/01_species_creation.jl")
include("lib/02_model_internals.jl")

# ## Initiation
#
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

# For the burn-in we want to keep the landscape uniform so we will use an empty
# landscape

environment_burnin = zeros(Float64, landscape_size)

# Now we can populate the species metadata

set_trophic_levels!(trophic_level)
set_interaction_strength!(interaction_strength; trophic_level)
set_environmental_optimum!(environmental_optimum, environment_burnin, trophic_level)
set_dispersal_rate!(dispersal_rate)
set_dispersal_decay!(dispersal_decay; trophic_level)

# Now we can create the matrix that will store each timedtamp of our burnin
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

# ## Community generation
#
# Now that we have created a somewhat 'stocahstic' community it is time to
# introduce some environmental variability - we will do this using
# `NeutralLandscapes.jl`. But before we do that lets specify our 'heating'
# generations time. 

# > Note because of how the simulation function is currently set up we will be
# > setting two generation time variables. `generations_heating` is the actual
# > number of generations and 'generations' is simly to make the simulate!
# > function go brrr.

generations_heating = 200
generations = 2

# becuase we are 'ramping' the environment (i.e. simulating a gradual change) we
# are creating a matrix that will keep each environmental layer so we can sample
# it for each simulation run. 

environment_heating = fill(0.0, (landscape_size..., generations_heating))

# then we can set the inital and final environemntal state

environment_heating[:,:,1] = environment_burnin
environment_heating[:,:,end] = rand(DiamondSquare(), landscape_size) .* species_richness

# now we can calculate the magnitude of change for each timestep to have the
# landscape reach its 'final state'

heating_step = (environment_heating[:,:,end] - environment_burnin)./generations_heating

# now we can sequentially add this value to each intermediate landscape state

for e in 2:(generations_heating - 1)
    environment_heating[:,:,e] = environment_heating[:,:,e-1] + heating_step
end

# **TODO** to make the environmental change a bit more gradual we can log
# transform the values

# we also need to reassign the envirnomental optima of all species (recall these
# were optimised to the uniform landscape)

set_environmental_optimum!(environmental_optimum, environment_heating[:,:,end], trophic_level)

# create the new metacommuity matrix and assign the first timestep as the
# ubandance values of the final timestep of the burnin community.

metacommunity = fill(0.0, (landscape_size..., species_richness, generations_heating))
metacommunity[:, :, :, 1] .= metacommunity_burnin[:, :, :, end]

# run the simulations

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