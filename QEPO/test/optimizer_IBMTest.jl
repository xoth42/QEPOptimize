
using Revise

using QEPO.Configurable

config = Configuration()

using QEPO.Optimizer
pop = Population()
config.num_simulations = 10000

initialize_pop_with_constraints!(pop,config)
run_with_constraints_history!(pop,config)

using QEPO.Visualizer: display_top_circuits, plot_performance_metrics

display_top_circuits(pop.individuals,10)

plot_performance_metrics(pop)