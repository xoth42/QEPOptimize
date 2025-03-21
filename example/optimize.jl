using Revise
includet("../src/QEPO.jl") # TODO once we have an actual library this will be just a `using` statement

# Revise.track("src/evolve.jl")
using .QEPO

using .QEPO: initialize_pop!,step!,Performance, NetworkFidelity # TODO export these

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
# Best 
points = Point2f.(gens,bestPerfs_ppf)
lines!(ax,gens,bestPerfs_ppf)
# worst
points = Point2f.(gens,worstPerfs_ppf)
lines!(ax,gens,worstPerfs_ppf)
display(fig)

# TODO plot F_in vs F_out and F_in vs P for the best circuits



# Function to optimize given the F_in
function optimize_f_in(f_in;steps=100)
    this_pop = Population()
    noises = [NetworkFidelity(f_in)]
    initialize_pop!(this_pop;noises=noises)
    for _ in 1:steps
        step!(this_pop;num_simulations=500,noises=noises)
    end
    return this_pop.individuals[1] # returns the best individual
end




f_ins = LinRange(0,1,10)
# individuals = [optimize_f_in(f_in;steps=500) for f_in in f_ins]

using OhMyThreads: tmap

individuals = tmap(optimize_f_in, f_ins)

f_outs = [indiv.performance.logical_qubit_fidelity for indiv in individuals]
p_outs = [indiv.performance.purified_pairs_fidelity for indiv in individuals]

# the graphs side by side
fig = Figure()
ax1 = Axis(fig[1,1])
ax1.title = "F_in vs F_out"
points = Point2f.(f_ins,f_outs)
lines!(ax1,f_ins,f_outs)

# f_in by P 
ax2 = Axis(fig[1,2])
ax2.title = "F_in vs P"
points = Point2f.(f_ins,p_outs)
lines!(ax2,f_ins,p_outs)

# add a simple y = x line 
lines!(ax1, 0:0.1:1, 0:0.1:1, color="black")
lines!(ax2, 0:0.1:1, 0:0.1:1, color="black")

display(fig)