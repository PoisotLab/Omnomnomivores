## Method with matrix

patch_position = CartesianIndices((1:_landscape_size[1], 1:_landscape_size[2]))
patch_distance = zeros(Float64, (prod(_landscape_size), prod(_landscape_size)))
_patch_distance = zeros(Float64, _landscape_size)
for x in axes(patch_distance, 1)
    for i in axes(patch_position, 1), j in axes(patch_position, 2)
        _patch_distance[i, j] =
            sqrt(sum(((patch_position[x].I) .- (patch_position[i, j].I)) .^ 2.0))
    end
    patch_distance[x, :] = vec(_patch_distance)
end
patch_distance

function _immigration(
    community_abundance,
    species_id,
    dispersal_rate::Vector{Float64},
    dispersal_decay::Vector{Float64},
    patch_distance::Matrix{Float64},
)
    for i in axes(patch_distance, 1) # 🐛 it is this indexing...
        _comm_vector = vec(community_abundance)
        return sum(
            dispersal_rate[species_id] * _comm_vector[l] *
            exp(-dispersal_decay[species_id] * patch_distance[i, l]) for
            l in axes(patch_distance, 1)
        )
    end
end

## Method with in situ calculations

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

@time _immigration(
    community_abundance,
    species_id,
    dispersal_rate,
    dispersal_decay,
    patch_location,
    patch_position,
)

@time _immigration(
    community_abundance,
    species_id,
    dispersal_rate,
    dispersal_decay,
    patch_distance,
)

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


@time metacommunity_model(
    current_community,
    dispersal_rate,
    dispersal_decay,
    patch_distance,
    environment_value,
    environmental_optimum,
    interaction_strength,
)

@time metacommunity_model(
    current_community,
    dispersal_rate,
    dispersal_decay,
    patch_position,
    environment_value,
    environmental_optimum,
    interaction_strength,
)
