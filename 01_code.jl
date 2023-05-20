using Distributions
using NeutralLandscapes
using SpatialBoundaries
using Random

"""
    set_trophic_levels!

Assigns the trophic levels for all species in the community. The number for each
trophic level (1 = plant, 2 = herbivore, 3 = carnivore) can be specified but
default to 5:3:2 (plant:herbivore:carnivore)
"""
function set_trophic_levels!(
    trophic_level::Vector{Int8};
    plants::Float64 = 0.5,
    herbivores::Float64 = 0.3,
    carnivores::Float64 = 0.2,
)
    _p, _h, _ = [plants, herbivores, carnivores] ./ (plants + herbivores + carnivores)
    n_plants = ceil(Int, length(trophic_level) * _p)
    n_herbivores = ceil(Int, length(trophic_level) * _h)
    n_carnivores = length(trophic_level) - (n_plants + n_herbivores)
    trophic_level[1:n_plants] .= Int8(1)
    trophic_level[(n_plants + 1):(n_plants + n_herbivores)] .= Int8(2)
    trophic_level[(end - n_carnivores):end] .= Int8(3)
    return trophic_level
end

"""
    set_dispersal_decay!

Determines the strength of exponetial decrease in dispersal distance (𝐿) for
all species. The value is drawn from a normal distribution based on the trophic
level of the species. Normal distributionas are charactersied as follows: plants
μ = 0.3, herbivores μ = 0.2, and carnivores μ = 0.1. In all instances σ = μ/4.
"""
function set_dispersal_decay!(dispersal_decay::Vector{Float64}; trophic_level)
    p = [(L, 0.25L) for L in [0.3, 0.2, 0.1]]
    for s in axes(trophic_level, 1)
        dispersal_decay[s] = rand(Normal(p[trophic_level[s]]...))
    end
    return dispersal_decay
end

"""
    set_interaction_strength!

Determines the interaction strength between species by drawing from a uniform
distribution, which is determined by the trophic level of the two interacting
species. The ranges of the distributions are as follows: plant-plant (-0.1,
0.0), herbivore-plant (-0.3, 0.0), plant-herbivore (0.0, 0.1),
carnivore-herbivore (-0.1, 0.0), herbivore-carnivore (0.0, 0.08), all other
combinations are set to zero.

Interaction strength go FROM ROW, TO COLUMN
"""
function set_interaction_strength!(interaction_strength::Matrix{Float64}; trophic_level)
    plant_plant = Uniform(-0.1, 0.0)
    herb_plant = Uniform(-0.3, 0.0)
    plant_herb = Uniform(0.0, 0.1)
    pred_herb = Uniform(-0.1, 0.0)
    herb_pred = Uniform(0.0, 0.08)
    M = [
        plant_plant plant_herb nothing
        herb_plant nothing herb_pred
        nothing pred_herb nothing
    ]
    S = length(trophic_level)
    for (i, tl_i) in enumerate(trophic_level)
        for (j, tl_j) in enumerate(trophic_level)
            interaction_strength[i, j] =
                isnothing(M[tl_i, tl_j]) ? 0.0 : rand(M[tl_i, tl_j])
        end
    end
    return interaction_strength ./ (0.33S)
end

"""
    set_dispersal_rate!

Determines the proportion of individuals that will migrate for each species by
drawing from a normal distribution. The μ can be specified but defaults to 0.25,
σ = μ/4.
"""
function set_dispersal_rate!(
    dispersal_rate::Vector{Float64};
    mean_dispersal_rate::Float64 = 0.1,
)
    for s in axes(trophic_level, 1)
        dispersal_rate[s] = rand(Normal((mean_dispersal_rate, mean_dispersal_rate / 4)...))
    end
    return dispersal_rate
end

"""
    set_environmental_optimum!

Determines the environmental optimum for all species. Optima are equally
distributed across the entire environmental range. This is done seperately for
each trophic level.
"""
function set_environmental_optimum!(
    environmental_optimum::Vector{Float64},
    environment_value::Matrix{Float64},
    trophic_level::Vector{Int8},
)
    env_range = extrema(environment_value)
    for tl in sort(unique(trophic_level))
        sp_at_level = findall(==(tl), trophic_level)
        if ~isempty(sp_at_level)
            environmental_optimum[sp_at_level] .=
                LinRange(env_range..., length(sp_at_level))
        end
    end
    return environmental_optimum
end

## Immigration term
# 🐛 I cant actually conceptualise distance and there must be a 🐛 here...

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
    immigration_counter = zero(eltype(abundances))
    for i in axes(abundances, 1)
        for j in axes(abundances, 2)
            d = sqrt(sum(((i, j) .- location) .^ 2.0))
            immigration_counter += exp(-decay * d) * (rate * abundances[i, j])
        end
    end
    return immigration_counter
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
    h = 1.0,
    σ = 2.0,
)
    Δ = landscape[patch...] - environmental_optimum[species]
    ξ = 2σ^2.0
    modifier = exp(-(Δ^2.0) / (ξ))
    return h * (modifier - 1)
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
    for generation in 1:(generations - 1)
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
                    immigration = emmigration = 0.0
                    growth_rate = exp(rate_of_increase + environment + interaction)
                    migration_change = immigration - emmigration
                    metacommunity[x, y, species, generation + 1] =
                        metacommunity[x, y, species, generation] * growth_rate +
                        migration_change
                end
            end
        end
    end
    return metacommunity
end

## 'Workflow'

landscape_size = (5, 7)
species_richness = 12
generations = 5

interaction_strength = zeros(Float64, (species_richness, species_richness))
trophic_level = zeros(Int8, species_richness)
environmental_optimum = zeros(Float64, species_richness)
dispersal_decay = zeros(Float64, species_richness)
dispersal_rate = zeros(Float64, species_richness)

environment = rand(DiamondSquare(), landscape_size) .* 2.0 .- 1.0

set_trophic_levels!(trophic_level)
set_interaction_strength!(interaction_strength; trophic_level)
set_environmental_optimum!(environmental_optimum, environment, trophic_level)
set_dispersal_rate!(dispersal_rate)
set_dispersal_decay!(dispersal_decay; trophic_level)

# Set an initial metaco object
metacommunity = fill(10.0, (landscape_size..., species_richness, generations))

simulate!(
    metacommunity,
    dispersal_rate,
    dispersal_decay,
    environment,
    environmental_optimum,
    interaction_strength,
)
