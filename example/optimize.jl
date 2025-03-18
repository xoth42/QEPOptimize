using Revise
includet("../src/QEPO.jl") # TODO once we have an actual library this will be just a `using` statement

# Revise.track("src/evolve.jl")
using .QEPO

using .QEPO: initialize_pop!,step!,Performance # TODO export these

ENV["JULIA_DEBUG"] = "debug"
# ENV["JULIA_REVISE_POLL"] = 1
pop = Population()

initialize_pop!(pop) # TODO do a 3-register run

# TODO evolve for 100 generations

# ppf is 'purified_pairs_fidelity'
bestPerfs_ppf::Vector{Float64} =Float64[]
worstPerfs_ppf::Vector{Float64} =Float64[]
 
for generation in 1:100
    step!(pop)
    push!(bestPerfs_ppf,pop.individuals[1].performance.purified_pairs_fidelity)
    push!(worstPerfs_ppf,pop.individuals[100].performance.purified_pairs_fidelity)
end


# TODO plot history of the evolution (best, worst, etc fidelities over generations)

using CairoMakie
fig = Figure()
ax = Axis(fig[1,1])
ax.title = "Purified pairs probability per 100 generations" 
gens = 1:100
points = Point2f.(gens,bestPerfs_ppf)
lines!(ax,gens,bestPerfs_ppf)
display(fig)

# TODO plot F_in vs F_out and F_in vs P for the best circuits
