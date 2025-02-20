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

# create a test thread data object
# td = ThreadData()

# initialize_pop_with_constraints!(pop,config) 
run_with_constraints_history!(pop,config)

display_top_circuits(pop.individuals,5)
pop.individuals[1]
# plot_performance_metrics(pop)