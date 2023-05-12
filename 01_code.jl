## packages here
using Distributions
using NeutralLandscapes
using OffsetArrays #(at least at some point we'll be brave and try 0 indexing)
using Plots
using SpatialBoundaries
using Random

## Parameters and data framework

t = 5 # number of timestamps
_landscape_size = (20, 20)
_species_richness = 8 # number of species

current_community = fill(
    10.0, (
        _landscape_size...,
        _species_richness,
    ))

interaction_strength = zeros(Float64, (_species_richness, _species_richness))
trophic_level = zeros(Int8, _species_richness)
environmental_optimum = zeros(Float64, _species_richness)
dispersal_decay = zeros(Float64, _species_richness)
dispersal_rate = zeros(Float64, _species_richness)

environment_value = rand(EdgeGradient(), _landscape_size) #❗
patch_position = CartesianIndices((1:_landscape_size[1], 1:_landscape_size[2]))

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
    p = [(0.3, 0.3 / 4), (0.2, 0.2 / 4), (0.1, 0.1 / 4)]
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
"""
function set_interaction_strength!(interaction_strength::Matrix{Float64}; trophic_level)
    plant_plant = Uniform(-0.1, 0.0)
    herb_plant = Uniform(-0.3, 0.0)
    plant_herb = Uniform(0.0, 0.1)
    pred_herb = Uniform(-0.1, 0.0)
    herb_pred = Uniform(0.0, 0.08)
    for i in axes(interaction_strength, 2)
        for j in axes(interaction_strength, 2)
            if (trophic_level[i] == 1 && trophic_level[j] == 1)
                interaction_strength[i, j] =
                    rand(plant_plant, 1)[1] / 0.33 * length(trophic_level)
            elseif (trophic_level[i] == 2 && trophic_level[j] == 1)
                interaction_strength[i, j] =
                    rand(herb_plant, 1)[1] / 0.33 * length(trophic_level)
            elseif (trophic_level[i] == 1 && trophic_level[j] == 2)
                interaction_strength[i, j] =
                    rand(plant_herb, 1)[1] / 0.33 * length(trophic_level)
            elseif (trophic_level[i] == 3 && trophic_level[j] == 2)
                interaction_strength[i, j] =
                    rand(pred_herb, 1)[1] / 0.33 * length(trophic_level)
            elseif (trophic_level[i] == 2 && trophic_level[j] == 3)
                interaction_strength[i, j] =
                    rand(herb_pred, 1)[1] / 0.33 * length(trophic_level)
            else
                interaction_strength[i, j] = 0.0
            end
        end
    end
end

"""
    set_dispersal_rate!

Determines the proportion of individuals that will migrate for each species by
drawing from a normal distribution. The μ can be specified but defaults to 0.25,
σ = μ/4.
"""
function set_dispersal_rate!(
    dispersal_rate::Vector{Float64};
    mean_dispersal_rate::Float64 = 0.25,
)
    for s in axes(trophic_level, 1)
        dispersal_rate[s] = rand(Normal((mean_dispersal_rate, mean_dispersal_rate / 4)...))
    end
    return dispersal_rate
end

## Environmental optimum
#❗ how optima are distribted is sus...
# also could maybe be optimised??

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
    env_range = maximum(environment_value) - minimum(environment_value)
    n_plants = count(==(1), trophic_level)
    n_herbivores = count(==(2), trophic_level)
    n_carnivores = count(==(3), trophic_level)
    for i in 1:n_plants
        environmental_optimum[i] = env_range / n_plants * i
    end
    for j in 1:n_herbivores
        environmental_optimum[(n_plants + j)] =
            env_range / n_herbivores * j
    end
    for k in 1:n_carnivores
        environmental_optimum[n_plants + n_herbivores + k] = env_range / n_carnivores * k
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
    community_abundance,
    species_id,
    dispersal_rate::Vector{Float64},
    dispersal_decay::Vector{Float64},
    patch_location,
    patch_position::CartesianIndices,
)
    return sum(
        dispersal_rate[species_id] * community_abundance[i, j] *
        exp(
            -dispersal_decay[species_id] * sqrt(
                sum(
                    (
                        (patch_position[patch_location[1], patch_location[2]].I) .-
                        (patch_position[i, j].I)
                    ) .^ 2.0,
                ),
            ),
        ) for
        i in axes(patch_position, 1), j in axes(patch_position, 1)
    )
end

## Environmental effect term

"""
    _environmental_effect

Calculates the impact of environment (𝐴) on the abundance of the current
species in the current landscape patch. The scaling parameter and σ can be
specified but default to h = 300, σ = 50.
"""
function _environmental_effect(
    patch_location,
    species_id,
    environment_value::Matrix{Float64},
    environmental_optimum::Vector{Float64};
    h = 300,
    σ = 50,
)
    return h - (
        h * exp(
            -(
                environment_value[patch_location[1], patch_location[2]] -
                environmental_optimum[species_id]
            )^2 /
            (2σ^2),
        )
    )
end

## Interaction effect term

"""
    _interaction_effect

Calcualtes the per capita effect (𝐵) of all species on the abundance of the
current species based on interaction strength and current abundance.
"""
function _interaction_effect(
    patch_location,
    species_id,
    current_community,
    interaction_strength,
)
    return sum(
        interaction_strength[species_id, n] *
        current_community[patch_location[1], patch_location[2], n] for
        n in axes(current_community, 3)
    )
end

## Full metacommunity model

"""
    metacommunity_model

TODO
"""
function metacommunity_model(
    current_community,
    dispersal_rate::Vector{Float64},
    dispersal_decay::Vector{Float64},
    patch_position,
    environment_value::Matrix{Float64},
    environmental_optimum::Vector{Float64},
    interaction_strength::Matrix{Float64};
    rate_of_increase::Float64 = 0.05,
)
    _next_community = similar(current_community)
    for i in axes(current_community, 3)
        community_abundance = current_community[:, :, i]
        for j in axes(current_community, 2), k in axes(current_community, 2)
            patch_location = [j, k]
            species_id = i
            current_abundance = current_community[j, k, i]
            environment = _environmental_effect(
                patch_location,
                species_id,
                environment_value,
                environmental_optimum,
            )
            immigration = _immigration(
                community_abundance,
                species_id,
                dispersal_rate,
                dispersal_decay,
                patch_location,
                patch_position,
            )
            interaction = _interaction_effect(
                patch_location,
                species_id,
                current_community,
                interaction_strength,
            )
            emmigration = current_abundance * dispersal_rate[i]
            new_abundance =
                current_abundance * exp(rate_of_increase + interaction + environment) +
                immigration - emmigration
            _next_community[j, k, i] = new_abundance # this is extra but makes for clearer reading?
        end
    end
    return _next_community
end

## 'Workflow'

set_trophic_levels!(trophic_level)
set_interaction_strength!(interaction_strength; trophic_level)
set_environmental_optimum!(environmental_optimum, environment_value, trophic_level)
set_dispersal_rate!(dispersal_rate)
set_dispersal_decay!(dispersal_decay; trophic_level)

metacommunity_model(
    current_community,
    dispersal_rate,
    dispersal_decay,
    patch_position,
    environment_value,
    environmental_optimum,
    interaction_strength,
)
