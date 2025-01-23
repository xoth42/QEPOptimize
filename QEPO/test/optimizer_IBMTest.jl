# include("../src/Optimizer.jl")
# include("../src/Configurable.jl")

using Revise
using QEPO.Configurable


# Let's break down the key components of the simulation process, which can be summarized as 5 steps 

################################################
##     choose SC quantum computer and path    ##
################################################

# The first step is to choose the quantum computer we want to use for our simulations. 
# In this example, we’re working with IBM’s superconducting quantum computers, here i showcased two examples 
# Monte Carlo Simulation for calculating the performance of circuits 
num_simulations = 10000

# read the calibration data and Define the valid qubits

# Ex1: 27-qubit ibmq_kolkata
# ibmqc = "ibmq_kolkata"
# calibration_data = read_calibration_data("ibmq_kolkata_calibrations_2024-02-02.csv")
# measurement_duration = 1216e-9
# gate_times = 533e-9
# valid_qubits = ([ 4, 7, 10, 12, 15, 18])                  # select path (6 qubits) with high fidelity 
# valid_qubits = ([ 1, 4, 7, 10, 12, 15, 18, 21])           # select path (8 qubits) with high fidelity 

# Ex2: 127-qubit ibmq_sherbrooke
ibmqc = "ibmq_sherbrooke"
calib_path = "QEPO/test/ibm_sherbrooke_calibrations_2024-10-09.csv"

calibration_data = read_calibration_data(calib_path)

t1 = 286e-6
t2 = 251e-6
gate_times = 533e-9
measurement_duration = 640e-9
valid_qubits = ([ 43, 44, 45, 46, 47,  48, 49, 50])

# For this simulation, we select 6-10 qubits from ibmq_sherbrooke.
# Calibration data for ibmq_sherbrooke is read from a file, and we focus on valid qubits that exhibit the lowest error rates.

# noise_model = create_T1_T2_noise_model(calibration_data, valid_qubits, measurement_duration, gate_times)

##################################
## initialization of population ##
##################################

# Next, we initialize the population of quantum circuits for optimization. Here’s what we define:
n = 6                                # number of raw bell pairs
k = 1                                # number of purified pairs
r = 4                                # number of registers (length of the valid_qubits)
Configurable.CostFunction(0)
optimize_for = CostFunction(0)       # 0: logical_qubit_fidelity  1: purified_pairs_fidelity 2 : average_marginal_fidelity
code_distance = 3                    # for logical_qubit_fidelity
f_in = 0.90                          # fidelity that being generated at communication qubits 
                                     # f_in = 0.9 would be appropriate

########################################
##  parameters for genetic algorithm  ##
########################################

starting_ops = 17                    # should be adapted to r  
max_ops = 17                         # should be appropriate!
population_size = 20                 # target population after distillation 
starting_pop_multiplier = 200        # multiplier 20 or 200
max_gen = 20

pairs = 20
children_per_pair = 3
mutants_per_individual_per_type = 5

p_lose_operation = 0.9
p_add_operation = 0.7
p_swap_operations = 0.8
p_mutate_operations = 0.8
individuals = []
selection_history = Dict()


hardware_config = HardwareConfiguration(calib_path,calibration_data,valid_qubits)
advanced_config = AdvancedConfiguration(code_distance,f_in,population_size,starting_pop_multiplier,starting_ops,pairs,children_per_pair,mutants_per_individual_per_type,p_lose_operation,p_add_operation,p_swap_operations,p_mutate_operations)
config = Configuration(num_simulations,n,k,r,optimize_for,max_gen,max_ops,hardware_config,advanced_config)
############################
##  optimization process  ##
############################

# Once the population is initialized, the optimization process begins. 
using QEPO.Optimizer
using Revise

population::Population = Population(individuals,selection_history)

run_with_constraints_history!(population,config)
# initialize_pop_with_constraints!(population,config)
# population = Population_hardware(n, k, r, code_distance, optimize_for, f_in, population_size, starting_pop_multiplier, max_gen, max_ops, starting_ops, pairs, children_per_pair, mutants_per_individual_per_type, p_lose_operation, p_add_operation, p_swap_operations, p_mutate_operations, individuals, selection_history, num_simulations)

# We run the genetic algorithm with constraints derived from the calibration data and the valid qubit paths we defined earlier.

# run_with_constraints!(population, calibration_data, valid_qubits)
#run_with_constraints_history!(population, calibration_data, valid_qubits)

using QEPO.Configurable
# user
config = DEFAULT_CONFIGURATION()

using QEPO.Optimizer
pop = Population()
initialize_pop_with_constraints!(pop,config)
run_with_constraints_history!(pop,config)
