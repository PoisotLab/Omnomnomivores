## packages here
using Distributions
using NeutralLandscapes
using Plots
using Random

## Assets to store things to

# number of timestamps
t = 5

# create the matrices for abundance
X = Array{Any}(undef, t)

## Parameters

# ❗ TODO this is a uniform, linear landscape with the same 'E' value will vary later
M = fill(0, (5, t))
S = 8 # number of species
Ci = 0.05 # rate of increase
h = 300 # scaling param
σ = 50 # std dev relating to env effect
X[1] = fill(10, (S, S)) # abundance of spp in patch at t = 1 (there is value in zero indexing and I hate that Phil is right...)

# Species metadata
# for now only an id and trophic level

# ❗ TODO trophic level can be more representative of real world ratios
Sm = [1:S rand(1:3,S)]

## Per capita effect (rules)

# empty interaction matrix
B = zeros(Float64, S, S)

Random.seed!(66) #Setting the seed

plant_plant = Uniform(-0.1, 0.0)
herb_plant = Uniform(-0.3, 0.0)
plant_herb = Uniform(0.0, 0.1)
pred_herb = Uniform(-0.1, 0.0)
herb_pred = Uniform(0.0, 0.08)

# determine 'interaction strength'
# ❗ TODO actually populate the entire (1/2) matrix (bonus make matrix sparse)
for i in 1:S-1
    j = i+1
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




# e.g. 'sketch'

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
