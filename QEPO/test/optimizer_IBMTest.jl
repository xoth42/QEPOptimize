using Revise

using QEPO.Visualizer: display_top_circuits, plot_performance_metrics

using QEPO.Configurable
dataPath = "QEPO/data/ibm_sherbrooke_calibrations_2024-10-09.csv"
valid_qubits::Array{Int} = ([ 43, 44, 45, 46, 47,  48, 49, 50])
hw_cfg = HardwareConfiguration(dataPath,valid_qubits)
adv_cfg = AdvancedConfiguration()
config = Configuration(hw_cfg,adv_cfg)
config
using QEPO.Optimizer
pop = Population()
config.num_simulations = 100
config.max_gen = 20
# initialize_pop_with_constraints!(pop,config) 
run_with_constraints_history!(pop,config)
config.optimize_for
display_top_circuits(pop.individuals,5)
pop.individuals[1]
config.advanced_config.communication_fidelity_in
# plot_performance_metrics(pop)
config.advanced_config.population_size

print(pop.individuals[1].ops)
using QEPO.Optimizer: calculate_performance!
calculate_performance!(pop.individuals[1],20,1,4,config.optimize_for,config.advanced_config)
for i in 1:20
    calculate_performance!(pop.individuals[i],20,1,4,config.optimize_for,config.advanced_config)

end