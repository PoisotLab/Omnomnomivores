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
S = 8 # number of species
Ci = 0.05 # rate of increase
h = 300 # scaling param
σ = 50 # std dev relating to env effect
μa = 2.5
σa = μa*0.25

# _THE_ ArrayTM

# we'll get back to this
# the idea though? - a M x S x t miltidim array
#A = [x + y + z  for x in 1:M, y ∈ 1:S+3, z = 1:t+1]

# timestamp 1
A = [repeat([1.0], M) 1.0:M repeat([0.0], M) fill(10.0, (M, S))]
# ❗ re-index to -2? (guessing there is a zero...)

# Species metadata
# for now only an id, trophic level, env optima (set to zero for now)

# ❗ TODO trophic level can be more representative of real world ratios
Sm = [1:S rand(1:3,S) repeat([0], S)]

## Per capita effect (rules)

# empty interaction matrix
B = zeros(Float64, S, S)

Random.seed!(66) # Execute order 66

plant_plant = Uniform(-0.1, 0.0)
herb_plant = Uniform(-0.3, 0.0)
plant_herb = Uniform(0.0, 0.1)
pred_herb = Uniform(-0.1, 0.0)
herb_pred = Uniform(0.0, 0.08)

# determine 'interaction strength'
# ❗ I think this could be a 'wide table' that appends to species metadata but also maybe not...
for i in 1:S
    for j in 1:S
        if (Sm[i,2] == 1 && Sm[j,2] == 1)
            B[i,j] = rand(plant_plant,1)[1]/0.33*S
        elseif(Sm[i,2] == 2 && Sm[j,2] == 1)
            B[i,j] = rand(herb_plant,1)[1]/0.33*S
        elseif(Sm[i,2] == 1 && Sm[j,2] == 2)
            B[i,j] = rand(plant_herb,1)[1]/0.33*S
        elseif(Sm[i,2] == 3 && Sm[j,2] == 2)
            B[i,j] = rand(pred_herb,1)[1]/0.33*S
        elseif(Sm[i,2] == 2 && Sm[j,2] == 3)
            B[i,j] = rand(herb_pred,1)[1]/0.33*S
        else
            B[i,j] = 0.0
        end 
    end
end

## Environmental optima

# normal distribution equally distributed across trophic levels
# setting to 0 for now in Sm

# use the max and min env variables from M to get range
# Neutral Landscapes means  range is 0 - 1.0
# for each trophic level (1-3) divide range by number species per level (n) and assign to spp
# collect(range(0, 1, length = n))
# ❗ might want to create a env _range_ centered around the optima - need to think about the σ of these though

## Waves hand and things happen (but only for one timestamp)

A1 = copy(A)

for i in 1:S
    a = rand(Normal(μa, σa))[1]
    for j in 1:M
        A1[j, i+3] = A[j, i+3]exp(Ci+sum(B[i,n]*A[1, n+3] for n in 1:S)+(h-h*exp(-((A[j,3]-Sm[i,3])^2/2σ^2))))-A[j, i+3]a
    end
end

## e.g. 'sketch'

ΔN(Nt,Pt) = λ*Nt*exp((1-Nt/K)-a*Pt)
ΔP(Nt,Pt) = n*Nt*(1-exp(-a*Pt))

ΔN(Nt, Pt)
ΔP(Nt, Pt)

t = 9e2

P = zeros(t)
N = zeros(t)
N[1] = 25.0
P[1] = 10.0

for i in 2:length(P)
    j = i-1
    N[i] = ΔN(N[j], P[j])
    P[i] = ΔP(N[j], P[j])
end

p1 = plot(eachindex(N), N, lab="Host", c=:black)
plot!(p1, eachindex(P), P, lab="Parasitoid", c=:red)
xaxis!(p1, "Time")
yaxis!(p1, "Abundances")

p2 = plot(N, P, leg=false, c=:grey)
scatter!(p2, [N[1]], [P[1]], c=:black)
xaxis!(p2, (0, maximum(N)), "Host")
yaxis!(p2, (0, maximum(P)), "Parasitoid")

plot(p1, p2)
