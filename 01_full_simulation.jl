using Distributions
using CairoMakie
using Makie.Colors
using NeutralLandscapes
using Random
using ProgressMeter
using Loess

#Random.seed!(66) # the time has come

include("lib/01_species_creation.jl")
include("lib/02_model_internals.jl")
include("lib/03_simulation_setup.jl")


landscape_size = (26, 26)
species_richness = 100
#connectivity = collect(0.1:0.2:0.99)
connectivity = [0.1, 0.5, 0.99]
metacommunity = []
environment = []

for i in eachindex(connectivity)

    landscape = rand(DiamondSquare(connectivity[i]), landscape_size)
    comm = OmnomnomCommunity(species_richness)
    sim = OmnomnomSimulation(landscape, 0.5species_richness, 500, 1000, 500)
    
    setup!(comm, sim; plants=5, herbivores=3, carnivores=2)
    
    schedule = (x) -> EnvironmentalChange(x; b=1)
    metacomm, env, optima = simulate(comm, sim; schedule=schedule, h=30.0, σ=50.0)
    
    push!(metacommunity, metacomm)
    push!(environment, env)

end



## this is for some colour allocation

palette = (plant=colorant"#00798c", herbivore=colorant"#edae49", carnivore=colorant"#d1495b")

species_col = fill(colorant"#ffffff", length(comm.trophic_level))
for i in axes(species_col, 1)
    if comm.trophic_level[i] == 1
        species_col[i] = palette.plant
    elseif comm.trophic_level[i] == 2
        species_col[i] = palette.herbivore
    else
        species_col[i] = palette.carnivore
    end
end

fig = Figure(; resolution=(1000, 900))
axs = [
    Axis(fig[1, 1];
        xlabel = "Environment value",
        ylabel = "Abundance"),
    Axis(fig[1, 2];
        xlabel = "Environment value",
        ylabel = "Abundance"),
    Axis(fig[1, 3];
        xlabel = "Environment value",
        ylabel = "Abundance"),
    Axis(fig[2, 1:3];
        xlabel = "Generation",
        ylabel = "Abundance"),
    Axis(fig[3, 1:3];
        xlabel = "Generation",
        ylabel = "Species Richness"),
]
for species in axes(metacommunity, 3)
    tl = comm.trophic_level[species]
    if sum(metacommunity[:, :, species, end]) > 0.0
        X = vec(sim.landscape)
        Xr = range(extrema(X)...; step=1.0)
        Y = vec(metacommunity[:, :, species, end])
        Yr = predict(loess(X, Y, span=0.4), Xr)
        lines!(
            axs[tl], 
            Xr,
            Yr;
            color = palette[Int64(tl)]
        )
        scatter!(
            axs[tl],
            X,
            Y,
            color = palette[Int64(tl)],
            alpha=0.6,
            markersize=2
        )
    end
end

abund = dropdims(mapslices(sum, metacommunity; dims = (1, 2)); dims = (1, 2))./prod(landscape_size)
vlines!(axs[4], sim.proofing, color=:black)
vlines!(axs[4], sim.proofing+sim.baking, color=:black)
for species in axes(abund, 1)
    lines!(axs[4], abund[species, 1:end]; color = species_col[species], alpha=0.5)
end

ylims!(axs[1], low=0.0)
ylims!(axs[2], low=0.0)
ylims!(axs[3], low=0.0)

xlims!(axs[1], extrema(sim.landscape))
xlims!(axs[2], extrema(sim.landscape))
xlims!(axs[3], extrema(sim.landscape))

ylims!(axs[4], (0, maximum(abund)*1.1))

abund[findall(abund .> 0.0), 1] .= 1.0

plt_idx = findall(comm.trophic_level .== 0x01)
hrb_idx = findall(comm.trophic_level .== 0x02)
crn_idx = findall(comm.trophic_level .== 0x03)

n_plt = vec(sum(abund[plt_idx, :]; dims = 1))
n_hrb = n_plt .+ vec(sum(abund[hrb_idx, :]; dims = 1))
n_crn = n_hrb .+ vec(sum(abund[crn_idx, :]; dims = 1))

band!(axs[5], 1:length(n_plt), 0.0, n_crn; color = palette.carnivore, label = "carnivore")
band!(axs[5], 1:length(n_plt), 0.0, n_hrb; color = palette.herbivore, label = "herbivore")
band!(axs[5], 1:length(n_plt), 0.0, n_plt; color=palette.plant, label = "plant")
axislegend(; position=:lb, nbanks=3)

vlines!(axs[5], sim.proofing, color=:black)
vlines!(axs[5], sim.proofing+sim.baking, color=:black)

xlims!(axs[4], (0, size(metacommunity, 4)))
xlims!(axs[5], (0, size(metacommunity, 4)))
ylims!(axs[5], (0, species_richness))

current_figure()
save("figures/diagnostics.png", fig)

#=
vec(mapslices(maximum, environment, dims=(1,2))) |> scatter
vec(mapslices(minimum, environment, dims=(1,2))) |> scatter!
vec(mapslices(median, environment, dims=(1,2))) |> scatter!
current_figure()
=#