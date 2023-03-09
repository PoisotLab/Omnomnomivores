# packages here
using Distributions
using NeutralLandscapes
using Plots
using Random


# parameters (using Thompson for now)

M = 200 # number of patches
S = 80 # number of species
Ci = 0.05 # rate of increase
h = 300 # scaling param
σ = 50 # std dev relating to env effect
X[1,1] = 10 # abundance of spp in patch

# Per capita effect (rules)

## drawn from nor


Random.seed!(123) # Setting the seed
d = Normal(μ=0.16, σ=0.05)
n=rand(d,1000)



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
