"""
    _immigration

Calculates the number of individuals of the current species that will immigrate
(𝐼) to the current landscape patch from all other landscape patches.
"""
function _immigration(
    abundances,
    location,
    rate,
    decay,
)
    d = zeros(Float64, size(abundances))
    for i in axes(abundances, 1)
        for j in axes(abundances, 2)
            d[i, j] = sqrt(sum(((i, j) .- location) .^ 2.0))
        end
    end
    d /= sum(d)
    k = exp.(-decay .* d .* 250.0)
    k[location...] = 0.0
    return rate * sum(k .* abundances)
end

## Environmental effect term

"""
    _environmental_effect

Calculates the impact of environment (𝐴) on the abundance of the current
species in the current landscape patch. The scaling parameter and σ can be
specified but default to h = 300, σ = 50.
"""
function _environmental_effect(
    metacommunity,
    species,
    patch,
    generation,
    landscape::Matrix{Float64},
    environmental_optimum::Vector{Float64};
    h = 300.0,
    σ = 50.0,
)
    Δ = landscape[patch...] - environmental_optimum[species]
    ξ = 2σ^2.0
    modifier = exp(-(Δ^2.0) / (ξ))
    return h * modifier - h
end

"""
    _interaction_effect

Calcualtes the per capita effect (𝐵) of all species on the abundance of the
current species based on interaction strength and current abundance.
"""
function _interaction_effect(
    metacommunity,
    species,
    patch,
    generation,
    interaction_strength,
)
    others = metacommunity[patch..., :, generation]
    return sum(interaction_strength[:, species] .* others)
end

"""
    metacommunity_model

TODO
"""
function simulate!(
    metacommunity::Array{Float64, 4},
    dispersal_rate::Vector{Float64},
    dispersal_decay::Vector{Float64},
    landscape::Matrix{Float64},
    environmental_optimum::Vector{Float64},
    interaction_strength::Matrix{Float64},
)
    @showprogress for generation in 1:(generations - 1)
        for species in axes(metacommunity, 3)
            rate_of_increase = trophic_level[species] == 1 ? 1e-1 : -1e-3
            abundances = view(metacommunity, :, :, species, generation)
            for x in axes(metacommunity, 1)
                for y in axes(metacommunity, 2)
                    patch_location = (x, y)
                    current_abundance = abundances[x, y]
                    environment = _environmental_effect(
                        metacommunity,
                        species,
                        patch_location,
                        generation,
                        landscape,
                        environmental_optimum,
                    )
                    immigration = _immigration(
                        abundances,
                        patch_location,
                        dispersal_rate[species],
                        dispersal_decay[species],
                    )
                    interaction = _interaction_effect(
                        metacommunity,
                        species,
                        patch_location,
                        generation,
                        interaction_strength,
                    )
                    emmigration = current_abundance * dispersal_rate[species]
                    growth_rate = exp(rate_of_increase + environment + interaction)
                    metacommunity[x, y, species, generation + 1] =
                        metacommunity[x, y, species, generation] * growth_rate
                    metacommunity[x, y, species, generation + 1] += immigration
                    metacommunity[x, y, species, generation + 1] -= emmigration
                    if metacommunity[x, y, species, generation + 1] <= 1e-3
                        metacommunity[x, y, species, generation + 1] = 0.0
                    end
                    if isinf(metacommunity[x, y, species, generation + 1])
                        metacommunity[x, y, species, generation + 1] = 0.0
                    end
                end
            end
        end
    end
    return metacommunity
end
