using QEPOptimize

using QEPOptimize: initialize_pop! # TODO export these

##

pop = Population()

initialize_pop!(pop) # TODO do a 3-register run

# TODO evolve for 100 generations

# TODO plot history of the evolution (best, worst, etc fidelities over generations)

# TODO plot F_in vs F_out and F_in vs P for the best circuits
