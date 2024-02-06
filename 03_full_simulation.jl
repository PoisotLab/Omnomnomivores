using Distributions
using CairoMakie
using Makie.Colors
using NeutralLandscapes
using Random
using ProgressMeter
using Loess

#Random.seed!(66) # the time has come

include("lib/01_species_creation.jl")
include("lib/02_model_internals.jl")

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

function setup!(
    comm::OmnomnomCommunity,
    sim::OmnomnomSimulation;
    plants = 0.5,
    herbivores = 0.3,
    carnivores = 0.2,
    mean_dispersal_rate = 0.1,
)
    plants, herbivores, carnivores = (plants, herbivores, carnivores) ./ (plants + herbivores + carnivores)
    sim.landscape .*= 0.5length(comm.trophic_level) # 'environmental amplitude'
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

EnvironmentalChange(x; b=1.0) = x^b / (x^b + (1-x)^b)

function simulate(
    comm::OmnomnomCommunity,
    sim::OmnomnomSimulation;
    schedule = (x) -> EnvironmentalChange(x; b=1.0), h=300.0, σ=50.0
)
    S = length(comm.trophic_level)
    L = size(sim.landscape)
    runtime = sim.proofing + sim.baking + sim.cooling
    tracker = zeros(Float64, (L..., S, runtime+1))
    tracker[:, :, :, 1] .= 1e-2
    # Pre-allocate the landscape and optimum
    ℒ = zeros(Float64, (L..., runtime+1))
    𝒪 = zeros(Float64, (S, runtime+1))
    L0 = fill(sim.proofing_value, L)
    O0 = fill(sim.proofing_value, S)
    ΔL = sim.landscape .- L0
    ΔO = comm.environmental_optimum .- O0
    ℒ[:,:,1] = L0
    𝒪[:,1] = O0
    for t in 1:runtime
        if t <= sim.proofing
            # Burn-in phase
            ℒ[:,:,t+1] = L0
            𝒪[:,t+1] = O0
        elseif sim.proofing < t <= (sim.proofing + sim.baking)
            # Progressive warmup phase
            f = (t - sim.proofing) / sim.baking # proportion of simulation done for logistic warmup
            ℒ[:,:,t+1] = L0 .+ ΔL .* schedule(f)
            𝒪[:,t+1] = O0 .+ ΔO .* schedule(f)
        else
            # Cooling phase
            ℒ[:,:,t+1] = sim.landscape
            𝒪[:,t+1] = comm.environmental_optimum
        end
    end
    # Main simulation loop
    @showprogress for t in 1:runtime
        Threads.@threads for s in 1:S
            rate_of_increase = comm.trophic_level[s] == 0x01 ? 1e-1 : -1e-3
            # Environmental effect - the maximum value is the ACTUAL max, not the OBSERVED max
            Δ = ((ℒ[:,:,t] .- 𝒪[s,t]).^2.0) ./ (2*σ*σ)
            A = (1/(σ*sqrt(2π)) * h) .* (exp.(-Δ) .- 1.0)
            
            for x in axes(ℒ[:,:,t], 1)
                for y in axes(ℒ[:,:,t], 2)
                    if tracker[x, y, s, t] > 0.0
                        patch_location = (x, y)
                        #=environment = _environmental_effect(
                            tracker,
                            s,
                            patch_location,
                            t,
                            ℒ[:,:,t],
                            𝒪[:,t],
                        )=#
                        interaction = _interaction_effect(
                            tracker,
                            s,
                            patch_location,
                            t,
                            comm.interaction_strength,
                        )
                        growth_rate = exp(rate_of_increase + A[x,y] + interaction)
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
        oops = findall(x -> x <= 1e-3, tracker[:,:,:,t])
        if !isempty(oops)
            for oop in oops
                tracker[oop.I..., t+1] = 0.0
            end
        end
    end
    # End
    return tracker, ℒ, 𝒪
end

# Simulation starts here
landscape_size = (20, 10)
species_richness = 80
landscape = rand(DiamondSquare(0.99), landscape_size)
comm = OmnomnomCommunity(species_richness)
sim = OmnomnomSimulation(landscape, 0.5species_richness, 500, 2000, 1000)
setup!(comm, sim; plants=5, herbivores=3, carnivores=2)

schedule = (x) -> EnvironmentalChange(x; b=1)
metacommunity, environment, optima = simulate(comm, sim; schedule=schedule, h=300.0, σ=50.0)

## this is for some colour allocation

palette = (plant=colorant"#00798c", herbivore=colorant"#edae49", carnivore=colorant"#d1495b")

species_col = fill(colorant"#ffffff", length(comm.trophic_level))
for i in axes(species_col, 1)
    if comm.trophic_level[i] == 1
        species_col[i] = palette.plant
    elseif comm.trophic_level[i] == 2
        species_col[i] = palette.herbivore
    else
        species_col[i] = palette.carnivore
    end
end

fig = Figure(; resolution=(1000, 900))
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
    if sum(metacommunity[:, :, species, end]) > 0.0
        X = vec(sim.landscape)
        Xr = range(extrema(X)...; step=1.0)
        Y = vec(metacommunity[:, :, species, end])
        Yr = predict(loess(X, Y, span=0.4), Xr)
        lines!(
            axs[tl], 
            Xr,
            Yr;
            color = palette[Int64(tl)]
        )
        scatter!(
            axs[tl],
            X,
            Y,
            color = palette[Int64(tl)],
            alpha=0.6,
            markersize=2
        )
    end
end

abund = dropdims(mapslices(sum, metacommunity; dims = (1, 2)); dims = (1, 2))./prod(landscape_size)
vlines!(axs[4], sim.proofing, color=:black)
vlines!(axs[4], sim.proofing+sim.baking, color=:black)
for species in axes(abund, 1)
    lines!(axs[4], abund[species, 1:end]; color = species_col[species], alpha=0.5)
end

ylims!(axs[1], low=0.0)
ylims!(axs[2], low=0.0)
ylims!(axs[3], low=0.0)

xlims!(axs[1], extrema(sim.landscape))
xlims!(axs[2], extrema(sim.landscape))
xlims!(axs[3], extrema(sim.landscape))

ylims!(axs[4], (0, maximum(abund)*1.1))

abund[findall(abund .> 0.0), 1] .= 1.0

plt_idx = findall(comm.trophic_level .== 0x01)
hrb_idx = findall(comm.trophic_level .== 0x02)
crn_idx = findall(comm.trophic_level .== 0x03)

n_plt = vec(sum(abund[plt_idx, :]; dims = 1))
n_hrb = n_plt .+ vec(sum(abund[hrb_idx, :]; dims = 1))
n_crn = n_hrb .+ vec(sum(abund[crn_idx, :]; dims = 1))

band!(axs[5], 1:length(n_plt), 0.0, n_crn; color = palette.carnivore, label = "carnivore")
band!(axs[5], 1:length(n_plt), 0.0, n_hrb; color = palette.herbivore, label = "herbivore")
band!(axs[5], 1:length(n_plt), 0.0, n_plt; color=palette.plant, label = "plant")
axislegend(; position=:lb, nbanks=3)

vlines!(axs[5], sim.proofing, color=:black)
vlines!(axs[5], sim.proofing+sim.baking, color=:black)

xlims!(axs[4], (0, size(metacommunity, 4)))
xlims!(axs[5], (0, size(metacommunity, 4)))
ylims!(axs[5], (0, species_richness))

current_figure()

#=
vec(mapslices(maximum, environment, dims=(1,2))) |> scatter
vec(mapslices(minimum, environment, dims=(1,2))) |> scatter!
vec(mapslices(median, environment, dims=(1,2))) |> scatter!
current_figure()
=#