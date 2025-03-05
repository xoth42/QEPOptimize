using Revise

using QEPO.Visualizer: display_top_circuits, plot_performance_metrics

using QEPO.Configurable
dataPath = "QEPO/data/ibm_sherbrooke_calibrations_2024-10-09.csv"
valid_qubits::Array{Int} = ([ 43, 44, 45, 46, 47,  48, 49, 50])
hw_cfg = HardwareConfiguration(dataPath,valid_qubits)
adv_cfg = AdvancedConfiguration()
adv_cfg.population_size = 100
adv_cfg.children_per_pair = 2
adv_cfg.starting_pop_multiplier = 400
adv_cfg.starting_ops = 5

config = Configuration(hw_cfg,adv_cfg)

config.raw_bell_pairs = 3
config.num_registers = 3
config.max_ops = 10
# config
using QEPO.Optimizer
pop = Population()
config.num_simulations = 1000
config.max_gen = 10
config.advanced_config.communication_fidelity_in = .5
# create a test thread data object
# td = ThreadData()
# config.optimize_for = average_marginal_fidelity
config.optimize_for = purified_pairs_fidelity

# config.advanced_config.p_add_operation = 0.5
# initialize_pop_with_constraints!(pop,config) 
# debug
# @run run_with_constraints_history!(pop,config)

run_with_constraints_history!(pop,config)

# sort!(pop.individuals)

# config.optimize_for
display_top_circuits(pop.individuals,3)
# pop.individuals[1]
# config.advanced_config.communication_fidelity_in
# # plot_performance_metrics(pop)
# config.advanced_config.population_size

# print(pop.individuals[1].ops)
# using QEPO.Optimizer: calculate_performance!
# calculate_performance!(pop.individuals[1],20,1,4,config.optimize_for,config.advanced_config)
# for i in 1:20
#     calculate_performance!(pop.individuals[i],20,1,4,config.optimize_for,config.advanced_config)

# end