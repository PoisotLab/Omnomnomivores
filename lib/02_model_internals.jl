"""
    immigration!

Calculates the number of individuals of the current species that will immigrate
(𝐼) to the current landscape patch from all other landscape patches.
"""
function immigration!(
    metacommunity,
    species,
    patch,
    generation,
    rate,
    decay,
)
    abundances = view(metacommunity, :, :, species, generation)
    d = zeros(Float64, size(abundances))
    for i in axes(abundances, 1)
        for j in axes(abundances, 2)
            d[i, j] = sqrt(sum(((i, j) .- patch) .^ 2.0))
        end
    end
    k = exp.(-decay .* d)
    k[patch...] = 0.0
    k /= sum(k)
    metacommunity[:, :, species, generation + 1] .+= (rate * abundances[patch...]) .* k
    return nothing
end

## Environmental effect term

"""
    _environmental_effect

Calculates the impact of environment (𝐴) on the abundance of the current
species in the current landscape patch. The scaling parameter and σ can be
specified but default to h = 4 000, σ = 5.
"""
function _environmental_effect(
    metacommunity,
    species,
    patch,
    generation,
    landscape::Matrix{Float64},
    environmental_optimum::Vector{Float64};
    h = 4000.0,
    σ = 5.0,
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
            rate_of_increase = trophic_level[species] == 1 ? 1e-1 : -1e-2
            for x in axes(metacommunity, 1)
                for y in axes(metacommunity, 2)
                    if metacommunity[x, y, species, generation] > 0.0
                        patch_location = (x, y)
                        environment = _environmental_effect(
                            metacommunity,
                            species,
                            patch_location,
                            generation,
                            landscape,
                            environmental_optimum,
                        )
                        interaction = _interaction_effect(
                            metacommunity,
                            species,
                            patch_location,
                            generation,
                            interaction_strength,
                        )

                        growth_rate = exp(rate_of_increase + environment + interaction)
                        metacommunity[x, y, species, generation + 1] +=
                            metacommunity[x, y, species, generation] * growth_rate

                        immigration!(
                            metacommunity,
                            species,
                            patch_location,
                            generation,
                            dispersal_rate[species],
                            dispersal_decay[species],
                        )
                        metacommunity[x, y, species, generation + 1] -=
                            metacommunity[x, y, species, generation] *
                            dispersal_rate[species]
                    end
                end
            end
        end
        bzzt_wrongo =
            findall(v -> (v <= 1e-3) | isinf(v), metacommunity[:, :, :, generation + 1])
        for bztwrg in bzzt_wrongo
            metacommunity[bztwrg.I..., generation + 1] = 0.0
        end
    end
    return metacommunity
end
