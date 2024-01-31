using Distributions
using CairoMakie
using Makie.Colors
using NeutralLandscapes
using Random
using ProgressMeter

Random.seed!(66) # the time has come

include("lib/01_species_creation.jl")
include("lib/02_model_internals.jl")

landscape_size = (20, 20)
species_richness = 80

mutable struct OmnomnomSimulation
    landscape::Matrix{Float64}
    proofing_value::Float64
    proofing::Int64
    baking::Int64
    cooling::Int64
end

mutable struct OmnomnomCommunity
    interaction_strength::Matrix{Float64}
    trophic_level::Vector{UInt8}
    environmental_optimum::Vector{Float64}
    dispersal_decay::Vector{Float64}
    dispersal_rate::Vector{Float64}
end

function OmnomnomCommunity(S::Int64)
    return OmnomnomCommunity(
        zeros(Float64, (S, S)),
        zeros(UInt8, S),
        zeros(Float64, S),
        zeros(Float64, S),
        zeros(Float64, S)
    )
end

comm = OmnomnomCommunity(20)
sim = OmnomnomSimulation(rand(10, 10), 10.0, 10, 20, 10)

function setup!(comm::OmnomnomCommunity, sim::OmnomnomSimulation; plants=0.5, herbivores=0.3, carnivores=0.2, mean_dispersal_rate=0.1,)
    set_trophic_levels!(
        comm.trophic_level;
        plants = plants,
        herbivores = herbivores,
        carnivores = carnivores,
    )
    set_interaction_strength!(comm.interaction_strength, comm.trophic_level)
    set_environmental_optimum!(comm.environmental_optimum, sim.landscape, comm.trophic_level)
    set_dispersal_rate!(comm.dispersal_rate, comm.trophic_level; mean_dispersal_rate=mean_dispersal_rate)
    set_dispersal_decay!(comm.dispersal_decay, comm.trophic_level)
    sim.landscape .*= length(comm.trophic_level)
    return comm, sim
end

setup!(comm, sim)

function simulate(comm::OmnomnomCommunity, sim::OmnomnomSimulation)
    S = length(comm.trophic_level)
    L = size(sim.landscape)
    runtime = sim.proofing + sim.baking + sim.cooling
    tracker = zeros(Float64, (L..., S, runtime))
    # Main simulation loop
    for generation in 1:runtime
        if generation <= sim.proofing
            # Burn-in phase
            this_landscape = fill(sim.proofing_value, L)
        elseif sim.proofing < generation <= (sim.proofing + sim.baking)
            # Progressive warmup phase
        else
            # Cooling phase
            this_landscape = sim.landscape
        end
    end
    # End
    return tracker
end

simulate(comm, sim)