## packages here
using Distributions
using NeutralLandscapes
using OffsetArrays #(at least at some point we'll be brave and try 0 indexing)
using Plots
using SpatialBoundaries
using Random


## Parameters

t = 5 # number of timestamps
M = 10.0 # number of habitat patches (uniform for now)
_landscape_size = (20, 20)
_species_richness = 8 # number of species
Ci = 0.05 # rate of increase
h = 300 # scaling param
σ = 50 # std dev relating to env effect
μa = 0.25 # mean dispersal rate
σa = μa*0.62 # std dev of dispersal rate

# _THE_ ArrayTM


current_community = zeros(
    Float64,
    (
        _landscape_size...,
        _species_richness
    )
)

interaction_strength = zeros(Float64, (_species_richness, _species_richness))
trophic_level = zeros(Int8, _species_richness)
environmental_optimum = zeros(Float64, _species_richness)
dispersal_decay = zeros(Float64, _species_richness)
dispersal_rate = zeros(Float64, _species_richness)

function set_trophic_levels!(trophic_level::Vector{Int8}; plants::Float64=0.5, herbivores::Float64=0.3, carnivores::Float64=0.2)
    _p, _h, _ = [plants, herbivores, carnivores] ./ (plants + herbivores + carnivores)
    n_plants = ceil(Int, length(trophic_level)*_p)
    n_herbivores = ceil(Int, length(trophic_level)*_h)
    n_carnivores = length(trophic_level) - (n_plants + n_herbivores)
    trophic_level[1:n_plants] .= Int8(1)
    trophic_level[(n_plants+1):(n_plants+n_herbivores)] .= Int8(2)
    trophic_level[(end-n_carnivores):end] .= Int8(3)
    return trophic_level
end

function set_dispersal_decay!(dispersal_decay::Vector{Float64}; trophic_level)
    p = [(0.3, 0.3/4), (0.2, 0.2/4), (0.1, 0.1/4)]
    for s in axes(trophic_level, 1)
        dispersal_decay[s] = rand(Normal(p[trophic_level[s]]...))
    end
    return dispersal_decay
end

patch_position =  CartesianIndices((1:_landscape_size[1], 1:_landscape_size[2]))
D = zeros(Float64, (prod(_landscape_size), prod(_landscape_size)))
for i in axes(patch_position, 1)
    for j in axes(patch_position, 2)
        D[i,j] = sqrt(sum(((patch_position[i].I) .- (patch_position[j].I)).^2.0))
    end
end

function set_interaction_strength!(interaction_strength::Matrix{Float64}; trophic_level)
    plant_plant = Uniform(-0.1, 0.0)
    herb_plant = Uniform(-0.3, 0.0)
    plant_herb = Uniform(0.0, 0.1)
    pred_herb = Uniform(-0.1, 0.0)
    herb_pred = Uniform(0.0, 0.08)
    for i in axes(interaction_strength, 1)
        for j in axes(interaction_strength, 2)
            if (trophic_level[i] == 1 && trophic_level[j] == 1)
                interaction_strength[i,j] = rand(plant_plant,1)/0.33*length(trophic_level)
            elseif(trophic_level[i] == 2 && trophic_level[j] == 1)
                interaction_strength[i,j] = rand(herb_plant,1)/0.33*length(trophic_level)
            elseif(trophic_level[i] == 1 && trophic_level[j] == 2)
                interaction_strength[i,j] = rand(plant_herb,1)/0.33*length(trophic_level)
            elseif(trophic_level[i] == 3 && trophic_level[j] == 2)
                interaction_strength[i,j] = rand(pred_herb,1)/0.33*length(trophic_level)
            elseif(trophic_level[i] == 2 && trophic_level[j] == 3)
                interaction_strength[i,j] = rand(herb_pred,1)/0.33*length(trophic_level)
            else
                interaction_strength[i,j] = 0.0
            end 
        end
    return interaction_strength
end

function set_dispersal_rate!(dispersal_rate::Vector{Float64}; mean_dispersal_rate::Float64=0.25)
    for s in axes(trophic_level, 1)
        dispersal_rate[s] = rand(Normal((mean_dispersal_rate, mean_dispersal_rate/4)...))
    end
    return dispersal_rate
end

_next_community = similar(current_community)

Random.seed!(66) # Execute order 66


## Environmental optima

# normal distribution equally distributed across trophic levels
# setting to 0 for now in Sm

# use the max and min env variables from M to get range
# Neutral Landscapes means  range is 0 - 1.0
# for each trophic level (1-3) divide range by number species per level (n) and assign to spp
# collect(range(0, 1, length = n))
# ❗ might want to create a env _range_ centered around the optima - need to think about the σ of these though


## Waves hand and things happen (but only for one timestamp)

function whaveter!(next, current; )
end

# function _Immigration
# function _Environmental_effect

A1 = copy(A)

for i in 1:S
    a = rand(Normal(μa, σa))
    for j in 1:M
        # growth = ...
        # trophic = ...
        # dispersal = ...
        A1[j, i+3] = A[j, i+3]exp(Ci+sum(B[i,n]*A[1, n+3] for n in 1:S)+(h-h*exp(-((A[j,3]-Sm[i,3])^2/2σ^2))))+sum(a*A[l, i+3]^(-Sm[i,4]*sqrt((A[j,1]-A[l,1])^2+(A[j,2]-A[l,2])^2)) for l in 1:M)-A[j, i+3]a
    end
end