"""
create_landscape

This function creates the various environental layers for each generation
"""
function set_landscape!(
    generations,
    landscape_connectivity,
    landscape,
    landscape_size,
    species_richness
)
    # set generations for different phases
    proofing = 1:round(Int32, generations*0.25)
    cooling = round(Int32, generations*0.75):generations
    baking = proofing[end]+1:cooling[1]

    # set the proofing landscapes to fixed value
    for i in proofing
        landscape[:, :, i, :] .= 10.0
    end

    # set the cooling landscape
    for i in cooling
        for k in eachindex(landscape_connectivity)
            landscape[:, :, i, k] .= rand(DiamondSquare(landscape_connectivity[k]), landscape_size) .* species_richness
        end
    end

    # get the interval of change for the landscape during the baking period
    heating_step = zeros(Float64, (landscape_size..., length(landscape_connectivity)))
    for k in eachindex(landscape_connectivity)
        heating_step[:, :, k] = (landscape[:, :, cooling[1], k] - landscape[:, :, proofing[1], k])/length(baking)
    end

    # set the baking landscape values
    for j in proofing[end]:cooling[1]-1
        for k in eachindex(landscape_connectivity)
            landscape[:, :, j+1, k] = (landscape[:, :, j, k] + heating_step[:, :, k]) * (1 / (1 + exp(-(j / length(baking)))))
        end
    end

    return landscape
end