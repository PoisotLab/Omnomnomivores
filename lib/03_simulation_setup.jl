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
    env_range = 0.5
)
    plants, herbivores, carnivores = (plants, herbivores, carnivores) ./ (plants + herbivores + carnivores)
    sim.landscape .*= 0.5*env_range # 'environmental amplitude'
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
