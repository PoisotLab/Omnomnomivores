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
        dispersal_rate[s] =
            rand(Truncated(Normal(mean_dispersal_rate, 0.25mean_dispersal_rate), 0.0, 1.0))
    end
    return dispersal_rate
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
function set_interaction_strength!(interaction_strength::Matrix{Float64}, trophic_level::Vector{Int8})
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
            if i == j
                interaction_strength[i, j] = tl_i == 1 ? -0.2 : -0.15
            end
        end
    end
    interaction_strength ./= (0.33 * S)
    return interaction_strength
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
