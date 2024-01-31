using Distributions
using CairoMakie
using Makie.Colors
using NeutralLandscapes
using Random
using ProgressMeter

#Random.seed!(66) # the time has come

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
        zeros(Float64, S),
    )
end

comm = OmnomnomCommunity(80)
sim = OmnomnomSimulation(rand(15, 15), 10.0, 200, 200, 200)

function setup!(
    comm::OmnomnomCommunity,
    sim::OmnomnomSimulation;
    plants = 0.5,
    herbivores = 0.3,
    carnivores = 0.2,
    mean_dispersal_rate = 0.1,
)
    sim.landscape .*= length(comm.trophic_level)
    set_trophic_levels!(
        comm.trophic_level;
        plants = plants,
        herbivores = herbivores,
        carnivores = carnivores,
    )
    set_interaction_strength!(comm.interaction_strength, comm.trophic_level)
    set_environmental_optimum!(
        comm.environmental_optimum,
        sim.landscape,
        comm.trophic_level,
    )
    set_dispersal_rate!(
        comm.dispersal_rate,
        comm.trophic_level;
        mean_dispersal_rate = mean_dispersal_rate,
    )
    set_dispersal_decay!(comm.dispersal_decay, comm.trophic_level)
    return comm, sim
end

setup!(comm, sim)

function RichardsCurve(f; M = 1.0, α = 1.08, c = 205.0)
    return M / (1 + α^(c * (0.5 - f)))
end

function simulate(
    comm::OmnomnomCommunity,
    sim::OmnomnomSimulation;
    schedule = (x) -> RichardsCurve(x),
)
    S = length(comm.trophic_level)
    L = size(sim.landscape)
    runtime = sim.proofing + sim.baking + sim.cooling
    tracker = zeros(Float64, (L..., S, runtime+1))
    tracker[:, :, :, 1] .= 0.1
    L0 = fill(sim.proofing_value, L)
    O0 = fill(sim.proofing_value, S)
    # Main simulation loop
    @showprogress for t in 1:runtime
        if t <= sim.proofing
            # Burn-in phase
            this_landscape = L0
            this_optimum = O0
        elseif sim.proofing < t <= (sim.proofing + sim.baking)
            # Progressive warmup phase
            f = (t - sim.proofing) / sim.baking # proportion of simulation done for logistic warmup
            this_landscape = L0 .+ (sim.landscape .- L0) .* schedule(f)
            this_optimum = O0 .+ (comm.environmental_optimum .- O0) .* schedule(f)
        else
            # Cooling phase
            this_landscape = sim.landscape
            this_optimum = comm.environmental_optimum
        end
        Threads.@threads for s in 1:S
            rate_of_increase = comm.trophic_level[s] == 0x01 ? 1e-1 : -1e-2
            for x in axes(this_landscape, 1)
                for y in axes(this_landscape, 2)
                    if tracker[x, y, s, t] > 0.0
                        patch_location = (x, y)
                        environment = _environmental_effect(
                            tracker,
                            s,
                            patch_location,
                            t,
                            this_landscape,
                            this_optimum,
                        )
                        interaction = _interaction_effect(
                            tracker,
                            s,
                            patch_location,
                            t,
                            comm.interaction_strength,
                        )
                        growth_rate = exp(rate_of_increase + environment + interaction)
                        tracker[x, y, s, t + 1] +=
                            tracker[x, y, s, t] * growth_rate

                        immigration!(
                            tracker,
                            s,
                            patch_location,
                            t,
                            comm.dispersal_rate[s],
                            comm.dispersal_decay[s],
                        )
                        tracker[x, y, s, t + 1] -=
                            tracker[x, y, s, t] *
                            comm.dispersal_rate[s]
                    end
                end
            end
        end
    end
    # End
    return tracker
end

metacommunity = simulate(comm, sim)

## this is for some colour allocation
species_col = fill("", length(comm.trophic_level))
for i in axes(species_col, 1)
    if comm.trophic_level[i] == 1
        species_col[i] = "green"
    elseif comm.trophic_level[i] == 2
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
for species in axes(metacommunity, 3)
    tl = comm.trophic_level[species]
    scatter!(
        axs[tl],
       # vec(comm.environmental_optimum),
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