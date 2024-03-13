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
