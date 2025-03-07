module Optimizer

# TODO:  add communication between qubits and hardware-specific errors from csv data to calculate_performance! Not the 'example error rates' 

# TODO update exports, clean up what actually needs to be exported and what can remain internal to this file
# TODO: improve clarity of use and document how to use everything all together
# TODO: provide and document full usage of optimizer.jl within QEPO.jl
# TODO document Individual,population, and docs for all functions
# TODO: update abstractions and data types
# TODO: rewrite gain_op_with_constraints

"""
    Organize the noise in the ciruit. Some noise is in mutate, calculate performance, etc. 

    Abstract away noise in initial state and in circuit

    Thermal relaxation notes for later:
    # add thermal relaxation noise T1 T2 into circuits
    # TODO: incorporate given hardware noise? 
    # t1_avg, t2_avg, gate_times = 286e-6, 251e-6, 533e-9  # example for testing: average T1, T2 and t on 'ibmq_sherbrooke'
    # λ₁, λ₂ = thermal_relaxation_error_rate(t1_avg, t2_avg, gate_times)
    # noisy_ops = add_thermal_relaxation_noise(indiv.ops, λ₁, λ₂)

"""


"""
    Optimizer.jl - Genetic optimizer from YipiaoWu/QuantumHardware

    ***SETTING UP PARAMETERS
    In order to run the optimizer, these parameters are needed:

    1.  General optimizer parameters 
        1a. Most important bits 
        Stored in the 'Configuration' struct

        num_simulations::Int        
        raw_bell_pairs::Int         # (n) incoming bell pairs
        purified_pairs::Int         # (k) outgoing bell pairs
        num_registers::Int          # amount of registers
        optimize_for::CostFunction  # optimization goal (see @enum CostFunction)
        max_gen::Int                # TODO more info
        max_ops::Int                # Limits the number of operations in each individual circuit


        1b.  Generic optimizer internals and quantum ciruit specifications
        Stored in the 'AdvancedConfiguration' struct

        code_distance::Int                              # for logical_qubit_fidelity
        communication_fidelity_in::Float64              # fidelity that being generated at communication qubits. 0.9 would be appropriate (f_in)
        population_size::Int                            #  target number of individuals in the population for each generation after initialization
        starting_pop_multiplier::Int                    #  A multiplier used to determine the size of the initial population
        starting_ops::Int                               # initial number of operations in each individual circuit of the population
        pairs::Int                                      # Number of parent pairs selected for breeding new individuals
        children_per_pair::Int                          # Number of offspring produced by each pair of parents
        mutants_per_individual_per_type::Int            # Number of mutations applied to each individual for each mutation type
        # Probabilities for various operations to be done during optimization:
        p_lose_operation::Float64
        p_add_operation::Float64
        p_swap_operation::Float64
        p_mutate_operation::Float64
        
        
        1c. Hardware specifications -> for IBM calibration data only, at the moment
        Stored in the 'HardwareConfiguration' struct

        calibration_data_path::String # Path to the csv file containing IBM calibration data
        calibration_data::Dict{Any,Any}
        valid_qubits::Array{Int}

    2. Formatting
        All of these parameters are stored like so, with the HardwareConfiguration and AdvancedConfiguration going inside of the Configuration struct, under the names 
        hardware_config and advancd_config, respectively. 

    ***RUNNING THE OPTIMIZER 
    Once the config is set up and stored in a variable, create your population that will be used for the optimizer. For example this makes an empty population:
            population = Population()
    Now, you can call these functions:
            1. run_with_constraints!(population::Population, config::Configuration)
                - Runs the basic simulation and updates your population with the results
            2. run_with_constraints_history!(population::Population, config::Configuration)
                - Runs the simulation and updates your population with added data tracking, and plots the fidelities to see the optimizer's progress once finished running.

    The functions with ! are meant to change every part of the population, not the config.
"""


# TODO: change these to an enum
hist_list = ["manual", "survivor", "random", "child", "drop_m", "gain_m", "swap_m", "ops_m"]

using Plots
using Base.Threads #  for multithreading
using Random
using Statistics
using LinearAlgebra
using BPGates # efficient representation of purification circuits
using QEPO.Configurable # import Configuration type from Configurable.jl
using QuantumClifford # general-purpose tools for Clifford circuits
using QuantumClifford.Experimental.NoisyCircuits
using BPGates: T1NoiseOp, T2NoiseOp
using OhMyThreads: tmap

# import Base: sort! # imports necessary for defining new methods of functions defined in Base

export generate_noisy_BellSwap_ops_for_individual, long_range_entanglement_generation!, Population, Performance, Individual, generate_valid_pairs,NoisyBellSwap, initialize_pop_with_constraints!,run_with_constraints_history!, ThreadData

### Genetic optimizer setup

"""A convenient structure to store various purification performance metrics."""
struct Performance
    """a list of probabilities as such [probability for no errors, probability for one Bell pair to be eroneous, probability for two Bell pairs to be eroneous, ..., probability for all Bell pairs to be eroneous]"""
    error_probabilities::Vector{Float64}
    """probability for no errors"""
    purified_pairs_fidelity::Float64
    """(TODO this is complicated and we are not using it right now) the fidelity of the logical qubit, as measured after teleportation with the purified bell pairs followed by error correction"""
    logical_qubit_fidelity::Float64
    """the average marginal fidelity of all purified pairs (note: correlation of errors is ignored here):
    Fidelity of i-th Bell pair is `Fᵢ = ⟨A|Trᵢ(ρ)|A⟩` where `Trᵢ` is "partial trace for all subspaces except the i-th one".
    The average marginal fidelity is `mean([F₁, F₂, ... Fₖ])`."""
    average_marginal_fidelity::Float64
    """the proportion of runs of a given protocol that do not have detected errors (i.e. Alice and Bob do not measure any errors)"""
    success_probability::Float64
end

# Define PauliNoiseBellGate and others to be droppable/not
is_droppable(::Any) = false
is_droppable(::PauliNoiseBellGate) = true
is_droppable(::NoisyBellMeasureNoisyReset) = true



mutable struct Individual
    history::String
    ops::Vector{Any}      # A vector containing a sequence of quantum operations that make up the individual's circuit
    performance::Performance
    fitness::Float64
    Individual() = new("", [], Performance(Float64[], 0.0, 0.0, 0.0, 0.0), 0.0)
    Individual(history::String) = new(history, [], Performance(Float64[], 0.0, 0.0, 0.0, 0.0), 0.0)
    Individual(history::String, ops::Vector{Any}, performance::Performance, fitness::Float64) = new(history, ops, performance, fitness)
end

mutable struct Population
    individuals::Vector{Individual}
    selection_history::Dict{String,Vector{Int64}} # Keeps track of the selection history for different types of individuals (e.g., survivors, mutants)
    # help in understanding which types of individuals contribute most to the improvement of the population
    Population() = new([], Dict{String, Vector{Int64}}())
    Population(individuals, selection_hist) = new(individuals,selection_hist)

end



# Create as struct for each thread's data used in calculate_performance
mutable struct ThreadData
    count_success::Int64
    counts_marginals::Vector{Int64}
    counts_nb_errors::Vector{Int64}
    err_count::Int64
end

### Quantum Circuit and qubit setup

"""convert two gate error rate to pauli noise representation"""
function p2_to_pauli(p2::Float64)
    px = py = pz = (1 - p2) / 4
    return px, py, pz
end

"""
    f_in_to_pauli(f_in::Float64)

Converts the f_in parameter to Pauli X, Y, and Z noise,
used in `calculate_performance!` to set up the initial noise.
"""
function f_in_to_pauli(f_in::Float64)
    px = py = pz = (1 - f_in) / 3
    return px, py, pz
end


""" Convert p2 of two Bell pairs to Pauli representation """
function map_p2s_to_pauli(p2_A::Float64, p2_B::Float64)

    px_A = py_A = pz_A = p2_A/ 4        # pauli noise on bell pair A (two qubit)
    px_B = py_B = pz_B = p2_B/ 4        # pauli noise on bell pair B (two qubit)

    px = px_A + px_B - px_A * px_B
    py = py_A + py_B - py_A * py_B
    pz = pz_A + pz_B - pz_A * pz_B

    return px, py, pz
end


""" Convert pswap of two Bell pairs to Pauli representation """
function map_pswap_to_pauli(p2_A::Float64, p2_B::Float64)

    # calculate effective error probabilities of each swap gate
    # px_swap = 1 - (1-px)^3 -（1-px)*px^2 = 3*px - 4*px^2 + 2*px^3 ≈ 3*px
    # py_swap = 1 - (1-py)^3 -（1-py)*px^2 = 3*py - 4*py^2 + 2*py^3 ≈ 3*py
    # pz_swap = 1 - (1-pz)^3 -（1-pz)*px^2 = 3*pz - 4*pz^2 + 2*pz^3 ≈ 3*pz

    # infidities of swap gate is as 3 times as fidilities of cnot gate
    # For simplicity, higher-order terms can be ignored because p2 are small
    px_A = py_A = pz_A = p2_A/ 4 *3
    px_B = py_B = pz_B = p2_B/ 4 *3

    px = px_A + px_B - px_A * px_B
    py = py_A + py_B - py_A * py_B
    pz = pz_A + pz_B - pz_A * pz_B

    return px, py, pz
end

""" readout error probabilities of a Bell pair """
function pair_readout_error(eta_A::Float64, eta_B::Float64)

    eta = eta_A + eta_B - eta_A * eta_B    # eta = 1- (1-eta_A)(1-eta_B) = 1 - (1- eta_A - eta_B - eta_A * eta_B)

    return eta
end

""" Function to add thermal relaxation T1 and T2 noise to each gate in the circuit, including idle times """
# retrieve average t1 and t2 time
function avg_t1_t2(calibration_data, valid_qubits)
    t1_times = Float64[]
    t2_times = Float64[]
    for q in valid_qubits
        t1 = calibration_data[q][:T1]
        t2 = calibration_data[q][:T2]
        push!(t1_times, t1)
        push!(t2_times, t2)
    end
    t1_avg = mean(t1_times)
    t2_avg = mean(t2_times)
    return t1_avg, t2_avg
end

function thermal_relaxation_error_rate(t1, t2, gate_time) # experimental private function for internal use
    λ₁ = 1 - exp(-gate_time/t1)
    t_ϕ = t1*t2 / (2*t1 - t2)
    λ₂ = 1 - exp(-gate_time/t_ϕ)
    return λ₁, λ₂
end

# TODO this function does not add thermal noise to the "wait" times in between gates on the qubits that are not acted upon on the current timestep
function add_thermal_relaxation_noise(circuit, λ₁, λ₂) # experimental private function for intenral use
    max_steps = 100
    max_pairs = 4
    time_table = falses(max_steps, max_pairs)   # tracks the activity of each qubit over time
    thermal_noisy_circuit = []
    step = 1
    for gate in circuit
        # move to next time step if the current gate cannot be added at this step
        if !can_fill_table(time_table, gate, step)         # check if a gate can be added to the current time step
            add_idle_noise!(thermal_noisy_circuit, time_table, step, λ₁, λ₂)  # Add noise to idle qubits at each time step.
            step = step +1  # Move to the next step
        end
        push!(thermal_noisy_circuit, gate)                # Apply the original gate
        push_noise!(thermal_noisy_circuit, gate, λ₁, λ₂)  # Apply noise for active qubits in the gate
        fill_table!(time_table, gate, step)  # Mark qubits as active for this gate
    end

    # a workaround for slow performance -- this is a bad piece of code -- the proper way to fix this is to hook into the `compactify_circuit` capabilities of QuantumClifford or to just not use Julia ;)
    # UType = Union{PauliNoiseBellGate{CNOTPerm}, T1NoiseOp, PauliNoiseBellGate{BellSwap}, T2NoiseOp, NoisyBellMeasureNoisyReset}
    # UType = QuantumOperation
    # thermal_noisy_circuit = convert(Vector{UType}, thermal_noisy_circuit);
    return thermal_noisy_circuit
end

### helper functions
# Checks if a gate can be scheduled at the given step without conflicts
function can_fill_table(time_table, gate, step)
    qubits = get_qubits_involved(gate)
    for q in qubits
        if time_table[step, q]    # if is true, it means that qubit q is already engaged in another operation at this time step, so the gate cannot be scheduled
            return false  # Conflict: qubit is already active at this time step
        end
    end
    return true  # No conflicts; gate can be scheduled
end

# mark the qubits involved in the current gate as "active"
function fill_table!(time_table, gate, step)
    qubits = get_qubits_involved(gate)
    for q in qubits
        time_table[step, q] = true  # Mark qubit as active
    end
end

# retrieve qubits involved in a gate
get_qubits_involved(gate::PauliNoiseBellGate) = [gate.g.idx1, gate.g.idx2]
get_qubits_involved(gate::NoisyBellMeasureNoisyReset) = [gate.m.sidx]
get_qubits_involved(gate) = []

# Multimethod to apply noise for different gate types (in order to be extendable)
function push_noise!(circuit, gate::PauliNoiseBellGate{T}, λ₁, λ₂) where T
    push!(circuit, T1NoiseOp(gate.g.idx1, λ₁))
    push!(circuit, T2NoiseOp(gate.g.idx1, λ₂))
    push!(circuit, T1NoiseOp(gate.g.idx2, λ₁))
    push!(circuit, T2NoiseOp(gate.g.idx2, λ₂))
end

push_noise!(circuit, gate::NoisyBellMeasureNoisyReset, λ₁, λ₂) = []  # No thermal relaxation added to measurement
push_noise!(circuit, any::Any, λ₁, λ₂) = [] 
function add_idle_noise!(circuit, time_table, step, λ₁, λ₂)
    for q in 1:size(time_table, 2)
        if !time_table[step, q]  # If the qubit is idle at this time step
            push!(circuit, T1NoiseOp(q, λ₁))
            push!(circuit, T2NoiseOp(q, λ₂))
        end
    end
end


"""
    reset_population!(population,population_size::Int,starting_pop_multiplier::Int)

    Initialize individuals with empty operations and default performance 
"""
function reset_population!(population,population_size::Int,starting_pop_multiplier::Int)
    reset_selection_history!(population)
    population.individuals=[
        Individual("random")
        for _ in 1:population_size * starting_pop_multiplier
    ]
end

### Setup algorithms


"""
    calculate_performance!(indiv::Individual, num_simulations::Int, purified_pairs::Int,number_registers::Int, optimize_for::CostFunction, code_distance::Int,communication_fidelity_in::Float64)


    updates the individual's performance and fidelity, and returns it
        
    simulates quantum operations using Monte Carlo trajectories to evaluate the performance of a given quantum purification circuit

    Uses a the supplied initial state
"""
function calculate_performance!(indiv::Individual, num_simulations::Int, purified_pairs::Int,number_registers::Int, optimize_for::CostFunction, code_distance::Int,communication_fidelity_in::Float64)
    count_success = 0
    counts_marginals = zeros(Int,purified_pairs) # an array to find F₁, F₂, …, Fₖ (tracks how often each purified bell pair is in the desired state)
    # probability of having 0, 1, 2, ... k 'erroneous' BP's
    counts_nb_errors = zeros(Int,purified_pairs+1) # an array to find P₀, P₁, …, Pₖ -- Careful with indexing it!
    # Threads.@threads for _ in 1:num_simulations # TODO from Stefan: this is a good place for threads
    initial_noise_circuit = [PauliNoiseOp(i, f_in_to_pauli(communication_fidelity_in)...) for i in 1:number_registers]
    # start with modeling 'clean' operations on noisy circuits

    for _ in 1:num_simulations
        
        # 1. randomize initial state
        initial_noisy_state, res = mctrajectory!(BellState(number_registers), initial_noise_circuit) # network noise model circuit
        purified_state, res = mctrajectory!(initial_noisy_state, indiv.ops) 
        # If the circuit execution was 'successful'
        if res == continue_stat
            count_success += 1
            err_count = 0
            for i in 1:purified_pairs # for each purified pair
                if purified_state.phases[2i-1] || purified_state.phases[2i] # checks whether an error has occurred based on binary representation in BPGates
                    err_count += 1
                else
                    counts_marginals[i] += 1 # tracks the i'th purified pair is in the desired state
                end
            end
            # For Julia indexing; index 1 corresponds to 0 errors, index 2 to 1 error, ...
            counts_nb_errors[err_count+1] += 1
        end
    end

        # count_success = 0 
        # counts_marginals = zeros(Int,purified_pairs) # an array to find F₁, F₂, …, Fₖ (tracks how often each purified bell pair is in the desired state)
        # # probability of having 0, 1, 2, ... k 'erroneous' BP's
        # counts_nb_errors = zeros(Int,purified_pairs+1) # an array to find P₀, P₁, …, Pₖ -- Careful with indexing it!

        # Threads.@threads for _ in 1:num_simulations # TODO from Stefan: this is a good place for threads
        # for _ in 1:num_simulations

        # Multithreading 
        # TODO: redo this implementation, so that it uses less space. right now it creates 2 n-vectors for each thread, which may be a lot
        # Get amount of threads
        # threads = Threads.nthreads()
    
        # # Fill a vector with empty structs for each thread
        # threads_data = fill(ThreadData(
        #     0,                              # count successes
        #     zeros(Int, purified_pairs),     # count marginals
        #     zeros(Int, purified_pairs+1),   # counts_nb_errors
        #     0),                             # err couunt
        # threads)
        
        # # Figure out how many simulations each thread will run
        # simulations_per_thread = div(num_simulations, threads)
        # # Loop through all of the data containers (eg, threads)
        # Threads.@threads for data in threads_data
        #     # Loop through all of the designated simulations for this thread
        #     for _ in 1:simulations_per_thread
        #         res_state, res = mctrajectory!(copy(initial_state), indiv.ops) # Simulates the operations on the initial bell state

        #         # If the circuit execution was 'successful'
        #         if res == continue_stat
        #             data.count_success += 1
        #             data.err_count = 0
        #             for i in 1:purified_pairs # for each purified pair
        #                 if res_state.phases[2i-1] || res_state.phases[2i] # checks whether an error has occurred based on binary representation in BPGates
        #                     data.err_count += 1
        #                 else
        #                     data.counts_marginals[i] += 1 # tracks the i'th purified pair is in the desired state
        #                 end
        #             end
        #             # For Julia indexing; index 1 corresponds to 0 errors, index 2 to 1 error, ...
        #             data.counts_nb_errors[data.err_count+1] += 1
        #         end
        #     end
        # end

        # # Combine the results from all threads
        # count_success =    sum([thread.count_success    for thread in threads_data])
        # counts_marginals = sum([thread.counts_marginals for thread in threads_data])
        # counts_nb_errors = sum([thread.counts_nb_errors for thread in threads_data])

        

    if count_success == 0
        @warn "No successful simulations; marginals and error probabilities will be undefined."
    
    end
    
    p_success = count_success    / num_simulations # proportion of successful simulations
    marginals = counts_marginals / count_success # marginal fidelities of individual purified pairs
    err_probs = counts_nb_errors / count_success # Distribution of errors across simulations : an array containing in each index i, how many errors occurred in (i-1)-bell-pair

    correctable_errors = div(code_distance  - 1, 2) # Maximum number of correctable errors based on code distance after teleportation
    indiv_logical_qubit_fidelity = sum(err_probs[1:min(end, correctable_errors+1)]) # Calculates the logical qubit fidelity by summing the probabilities of correctable errors

    # TODO: find out why logical_qubit_fidelity is being set to 1
    
    # Apply performance data to individual
    indiv.performance =  Performance(err_probs, err_probs[1],indiv_logical_qubit_fidelity, mean(marginals), p_success)

    
    # # debug 
    # if (indiv.performance.logical_qubit_fidelity == 1)
    #     @warn  "logical qubit fidelity is 1. p_success: $count_success, marginals: $marginals, err_probs: $err_probs, correctable_errors: $correctable_errors"
    # else 
    #     @warn "valid lqf. p_success: $count_success, marginals: $marginals, err_probs: $err_probs, correctable_errors: $correctable_errors"
    # end
    # Sets the fitness value based on the optimization goal
    if optimize_for == logical_qubit_fidelity
        indiv.fitness =  indiv.performance.logical_qubit_fidelity
    elseif optimize_for == purified_pairs_fidelity
        indiv.fitness =  indiv.performance.purified_pairs_fidelity
    elseif optimize_for == average_marginal_fidelity
        indiv.fitness =  indiv.performance.average_marginal_fidelity
    elseif optimize_for == success_probability
        indiv.fitness =  indiv.performance.success_probability
    else 
        indiv.fitness = 0.0
    end

    if count_success <= 0
        indiv.fitness = 0.0
        indiv.performance =  Performance(err_probs, 0,0, 0, 0)
        # should the rest of the individual data be set to zero? Otherwise, the data is skewed. Ideally this individual should be deleted
    end


    return indiv.performance

end


##################################################
###          initilize polulation              ###
##################################################

""" helper function to to map valid qubits to pairs """
function map_valid_qubits_to_pairs(valid_pair, valid_pairs)
    # create the dictionary where we are looking for the key
    valid_pair_to_num = Dict{Tuple{Int, Int}, Int}()
    for (i, pair) in enumerate(valid_pairs)
        valid_pair_to_num[pair] = i
        valid_pair_to_num[reverse(pair)] = i          # Ensures bidirectional access to pairs, treating (a, b) and (b, a) equivalently
    end
    return get(valid_pair_to_num, valid_pair, 0)       # retrieves the index of valid_pair from valid_pair_to_num, or 0 if valid_pair is not found
end

function map_num_to_valid_qubits(num::Int, valid_pairs)
    if num > 0 && num <= length(valid_pairs)
        return valid_pairs[num]
    end
    return (0, 0)
end


function generate_noisy_BellSwap_ops_for_individual(num_registers,valid_pairs,calibration_data)::Vector{Any}
    ##### Should num_gates:
    # num_gates = rand(1:get_starting_ops(advanced_config) - 1)  # Randomly determine the number of gates to include in each individual's circuit)
    ##### Be used here ??
    
    """ Create a sequence of BellSwap gate that will effectively move the qubits from the lowest index to the highest index in a structured manner"""
    swap_gates = []
    # for i in population.r:-1:2              # The outer loop ensures that for each qubit register from the second one to the topmost one

    for i in 2:1:num_registers
        for j in num_registers:-1:i          # The inner loop performs the swaps for the current qubit register i with all registers down to i-1
            push!(swap_gates, BellSwap(j, j-1))
        end
    end

    """ Wrap the noise of each BellSwap gate in a NoiseBellSwap based on the hardware's error rates """
    noisy_BellSwap = [
        begin
            pair_A, pair_B = map_num_to_valid_qubits(swap.idx1, valid_pairs), map_num_to_valid_qubits(swap.idx2, valid_pairs)
            # Are you indexing calibration_data with heterogenous keys? -- this is not a good idea in general
            # add two-qubit gate error rate from calibration data
            p2_A = calibration_data[pair_A].two_qubit_error
            p2_B = calibration_data[pair_B].two_qubit_error
            NoisyBellSwap(swap.idx1, swap.idx2, map_pswap_to_pauli(p2_A::Float64, p2_B::Float64)... )
        end
        for swap in swap_gates
    ]

    return noisy_BellSwap
end

""" Generate pairs in the order specified """
function generate_valid_pairs(valid_qubits)
    # Memoize this function so if called again for the same input, it will not need to recalculate
    memo = Dict{Array{Int, 1}, Array{Tuple{Int, Int}, 1}}()
    if haskey(memo, valid_qubits)
        return memo[valid_qubits]
    
    else
        valid_pairs = [(valid_qubits[i], valid_qubits[i + 1]) for i in 1:2:length(valid_qubits)-1]
        memo[valid_qubits] = valid_pairs
        return valid_pairs
    end
end

"""
    long_range_entanglement_generation!(population::Population,config::Configuration)

    For each individual in the population, set their operations to noisy BellSwaps. 

    The sequence of BellSwap gates that will effectively move the qubits from the lowest index to the highest index in a structured manner, and wraps the gates in noisy bellswaps.
"""
function long_range_entanglement_generation!(population::Population,config::Configuration)
    reset_population!(population,config.advanced_config.population_size,config.advanced_config.starting_pop_multiplier,)

    """ Generate pairs in the order specified """
    valid_pairs = generate_valid_pairs(config.hardware_config.valid_qubits)

    """ Each individual in the population is processed in parallel to speed up the initialization process """
    Threads.@threads for indiv in population.individuals


        # num_gates = rand(1:get_starting_ops(advanced_config) - 1)  # Randomly determine the number of gates to include in each individual's circuit
        
        noisy_BellSwap = generate_noisy_BellSwap_ops_for_individual(config.num_registers, valid_pairs,config.hardware_config.calibration_data)

        indiv.ops = noisy_BellSwap

    end

end

"""
    initialize_pop_with_constraints!(population::Population, config::Configuration)

 initilize a polulation of quantum circuits with constraints about gate connectivity  and ibm noise (cnot, measurement) 

It is always called with sort and cull after, so I am adding the sort and cull call to the end of this method.
"""
function initialize_pop_with_constraints!(population::Population, config::Configuration)
    # TODO: instead of passing entire config, only pass what is needed.
    # Attributes needed:
    #  population, calibration_data, valid qubits, num registers, purified pairs, starting ops, anything needed for reset_population!, communication_fidelity_in
    calibration_data = config.hardware_config.calibration_data
    valid_qubits = config.hardware_config.valid_qubits
    # For calibration_data, there at the *very least* MUST be some data validation so that we know type information
    # on calibration_data. Otherwise, we cannot take advantage of any of the type-based performance advantages that Julia offers.
    # Yipiao: Ensure calibration data is validated and properly typed?
    num_registers = config.num_registers
    purified_pairs = config.purified_pairs
    starting_ops = config.advanced_config.starting_ops

    reset_population!(population,config.advanced_config.population_size,config.advanced_config.starting_pop_multiplier)
   

    """ Generate pairs defining which qubits can interact (i.e., nearest neighbors) """
    valid_pairs = generate_valid_pairs(valid_qubits)
    # Creates pairs of consecutive qubits (e.g., (q0, q1), (q2, q3), (q4, q5))

    """ Each individual in the population is processed in parallel to speed up the initialization process """
    Threads.@threads for indiv in population.individuals
        ###### num_gates-> this is not used, but defined here and in similar locations. Is this intended to be used ?? TODO: find out from Yipiao/Stefan
        # num_gates = rand(1:get_starting_ops(advanced_config) - 1)  # Randomly determine the number of gates to include in each individual's circuit
       
        """ Create a sequence of BellSwap gate that will effectively move the qubits from the lowest index to the highest index in a structured manner, and wraps the gates in noisy bellswaps."""
        noisy_BellSwap = generate_noisy_BellSwap_ops_for_individual(num_registers,valid_pairs,calibration_data)

        """ Create Random CNOT Gates with Nearest Neighbor Constraints """
        random_gates::Vector{Any} = []
        # num_gates = rand(1:starting_ops-1)             # strategy 1: randomly determines the number of gates to be included in each individual's circuit.
        num_gates = rand(starting_ops*0.5 : starting_ops*0.7)  # strategy 2: control the ratio between the number of cnotperm gate and measurements

        for _ in 1:num_gates

            pair_num_A = rand(1:num_registers)

            neighbors = []
            if pair_num_A > 1
                push!(neighbors, pair_num_A - 1)
            end
            if pair_num_A < num_registers
                push!(neighbors, pair_num_A + 1)
            end

            pair_num_B = neighbors[rand(1:length(neighbors))]

            """ create bi-CNOT gate between two Bell pairs (source pair and scarifice pair)"""
            push!(random_gates, rand(CNOTPerm, (pair_num_A, pair_num_B)...))
        end

        """" Wrap the noise of each CNOTPerm gate in a PauliNoiseBellGate based on the hardware's error rates """
        noisy_random_gates = [
            begin
                pair_A, pair_B = map_num_to_valid_qubits(g.idx1, valid_pairs), map_num_to_valid_qubits(g.idx2, valid_pairs)
                # add two-qubit gate error rate from calibration data
                p2_A = calibration_data[pair_A][:two_qubit_error]
                p2_B = calibration_data[pair_B][:two_qubit_error]
                PauliNoiseBellGate(g, map_p2s_to_pauli(p2_A::Float64, p2_B::Float64)...)
            end
            for g in random_gates
        ]

        """" Generate Random Noisy Measurement Operations """
        random_measurements = [
            begin

                q = rand(purified_pairs+1 : num_registers)  # Randomly select a source register

                q_pair = map_num_to_valid_qubits(q, valid_pairs)

                # Check if the register is the bottom-most pair. If so, create a raw Bell pair with f_in
                if q == num_registers
                    # So this is effectively considering the phase damping that occurs only during measurement.
                    # Definitely need to consider two qubit gate durations for high depth circuits.
                    NoisyBellMeasureNoisyReset(
                        rand(BellMeasure, q),
                        pair_readout_error(calibration_data[q_pair[1]][:readout_error], calibration_data[q_pair[2]][:readout_error]),   #❓change and see difference
                        f_in_to_pauli(config.advanced_config.communication_fidelity_in)...                              # only consider F_in, not consider T1, T2
                        )
                else
                # Other registers are NoisyBellMeasureNoReset
                    NoisyBellMeasureNoisyReset(
                        rand(BellMeasure, q),
                        pair_readout_error(calibration_data[q_pair[1]][:readout_error], calibration_data[q_pair[2]][:readout_error]),
                        0.25, 0.25, 0.25                     # Trick of NoRest: create a maximally mixed state with px = py = pz = pI = 1/4
                    )
                end

            end
            for _ in 1:(starting_ops - num_gates)
        ]


        # TODO: delete these comments below

        # """ Get T1, T2 noise parameters from calibration data for the specific qubits """
        # noise_params = [
        #     begin
        #         a, b = valid_pair
        #         # t, T1, and T2 aren't defined here?-- they're symbols
        #         map_phase_damping_to_pauli( population.f_in,
        #             # What's t here?-- probably the time it takes to execute a single measurement
        #             calibration_data[a][:t], calibration_data[a][:T1], calibration_data[a][:T2],
        #             calibration_data[b][:t], calibration_data[b][:T1], calibration_data[b][:T2]
        #         )
        #         # +
        #         # map_amplitude_damping_to_pauli(
        #         #     calibration_data[a][:t], calibration_data[a][:T1],
        #         #     calibration_data[b][:t], calibration_data[b][:T1]
        #         # )
        #     end
        #     for valid_pair in valid_pairs
        # ]

        """ Keep BellSwap gates ordered while Randomizing the other operations to create a diverse set of quantum circuits """
        all_ops = vcat(noisy_random_gates, random_measurements)
        # noisy_random_gates_measurements = all_ops[length(noisy_BellSwap)+1:end]
        # shuffled_ops = vcat(noisy_BellSwap, noisy_random_gates_measurements[randperm(length(noisy_random_gates_measurements))])   
        shuffled_ops = vcat(all_ops[randperm(length(all_ops))])

        # indiv.ops =  convert(Vector{Union{PauliNoiseBellGate{CNOTPerm}, NoisyBellMeasureNoisyReset, PauliNoiseBellGate{BellSwap}}}, shuffled_ops)  # Converts the operations into a vector of gate types

        # indiv.ops =  convert(QuantumOperation, shuffled_ops)  # Converts the 
        indiv.ops = shuffled_ops
        # operations into a vector of gate types
        # TODO from Stefan: the convert above should not be necessary -- probably there is something else in the code that makes things messy if this is needed
    end

    simulate_and_sort!(population,config.num_simulations,config.purified_pairs,config.num_registers,config.optimize_for,config.advanced_config)
    cull!(population,config.advanced_config.population_size)

end


"""
    reset_selection_history!(population::Population)

    Reset the selection history for the population
"""
function reset_selection_history!(population::Population)
    # TODO from Stefan: this list gets repeated frequently, probably it makes sense to put it in a "convenience" global variable (and as mentioned elsewhere, probably a list of symbols, not a list of strings)
    for hist in hist_list 
        population.selection_history[hist] = Vector{Int64}()
    end

end


"""
    update_selection_history!(population::Population)

    Update the selection history for the individuals in the population
"""
function update_selection_history!(population::Population)
    for hist in hist_list
        push!(population.selection_history[hist], reduce(+, [1 for indiv in population.individuals if indiv.history  == hist], init=0)) # TODO from Stefan: `sum` is the same as `reduce(+)`
    end
end


##################################################
###          optimized process                 ###
##################################################

"""
    run_with_constraints!(population::Population, config::Configuration)

    Execution of a Genetic Algorithm Designed to Evolve a Population of Quantum Circuits 
"""
function run_with_constraints!(population::Population, config::Configuration) 
    # TODO from Stefan: the fact that there is `run` but also `run_with_constraints` kinda sounds like this can be written more neatly and simply if we use "multiple dispatch"  

    initialize_pop_with_constraints!(population, config)

    # Sort and cull are now called automatically with the initialize_pop_with_constraints call.
    # simulate_and_sort!(population,config)  # Evaluate the performance and Sorts the individuals by fitness, placing the best-performing individuals at the top
    # cull!(population,config.advanced_config.population_size)  # Removes excess individuals to maintain the target population size

    for _ = 1:config.max_gen

        # Produce the next generation of individuals by performing selection, crossover, mutation, and other genetic operations
        step_with_constraints!(population,config.max_ops,config.hardware_config.valid_qubits,config.purified_pairs,config.num_registers,config.hardware_config.calibration_data,config.num_simulations,config.optimize_for,config.advanced_config)
        update_selection_history!(population)
       
        # Calculate performance for each individual in parallel using OhMyThreads.jl
        tmap(indiv -> calculate_performance!(indiv,
            config.num_simulations,
            config.purified_pairs,
            config.num_registers,
            config.optimize_for,
            config.advanced_config.code_distance,
            config.advanced_config.communication_fidelity_in,
            ), population.individuals)
    end

end

# alternative one when printing out the history
"""
    run_with_constraints_history!(population::Population, config::Configuration)

    Includes a plot for the history of the optimization for single F_in
TBW
"""
function run_with_constraints_history!(population::Population, config::Configuration)
    max_gen = config.max_gen
    
    max_purified_fidelities_gen_dic = zeros(max_gen,1)
    min_purified_fidelities_gen_dic = zeros(max_gen,1)
    ave_purified_fidelities_gen_dic = zeros(max_gen,1)

    initialize_pop_with_constraints!(population, config)
    # simulate_and_sort!(population,config)
    # cull!(population,population_size)

    for i = 1:config.max_gen
        # Mutate
        step_with_constraints!(population,config.max_ops,config.hardware_config.valid_qubits,config.purified_pairs,config.num_registers,config.hardware_config.calibration_data,config.num_simulations,config.optimize_for,config.advanced_config)
        update_selection_history!(population)

        # Calculate performance for each individual in parallel using OhMyThreads.jl
        performances = tmap(indiv -> calculate_performance!(indiv,
            config.num_simulations,
            config.purified_pairs,
            config.num_registers,
            config.optimize_for,
            config.advanced_config.code_distance,
            config.advanced_config.communication_fidelity_in
            ), population.individuals)
        
        purified_fidelities = [perf.purified_pairs_fidelity for perf in performances]
        # Filter out NaN values
        valid_fidelities = filter(!isnan, purified_fidelities)
        
        max_purified_fidelities_gen_dic[i]= maximum(valid_fidelities)
        min_purified_fidelities_gen_dic[i]= minimum(valid_fidelities)
        ave_purified_fidelities_gen_dic[i]= mean(valid_fidelities)
        
        # check running progress 
        println("Running process: generation ",i)

    end

    # plot the history of the optimization for single F_in
    default(fontfamily="Times")
    generations = 1:max_gen
    plot(generations, max_purified_fidelities_gen_dic, label="Best purified pairs fidelity",lw=2,ylim=(0.7,1.01))
    plot!(generations, min_purified_fidelities_gen_dic, label="Worst purified pairs fidelity",lw=2,ylim=(0.7,1.01))
    plot!(generations, ave_purified_fidelities_gen_dic, label="Average purified pairs fidelity",lw=2,ylim=(0.7,1.01))
    xlabel!("Generation")
    ylabel!("Fitness")
    title!("Optimization History")

end


"""
    add_mutations!(individuals::Vector{Individual}, valid_qubits::Array{Int}, purified_pairs::Int, num_registers::Int,max_ops::Int,calibration_data::Dict{Any,Any}, adv_config::AdvancedConfiguration)

    Adds mutated individuals to the given individual vector, the amount mutated is decided by the probabilities of operations in the AdvancedConfiguration part of the Configuration.

    TODO: rewrite for efficiency
"""
function add_mutations!(individuals::Vector{Individual}, valid_qubits::Array{Int}, purified_pairs::Int, num_registers::Int,max_ops::Int,calibration_data::Dict{Any,Any}, adv_config::AdvancedConfiguration)
    # TODO: add 'traits' for gates that the mutations use

    # introduce diversity, enabling the population to explore new solutions
    # Apply Mutation Operations with Constraints (drop, gain, swap, mutate)

    # Store the children so this can use multithreading
    indivContainers = [[] for _ in 1:adv_config.population_size]

    # For every individual, up to the pop size limit
    Threads.@threads for i in 1:adv_config.population_size
        indiv = individuals[i]
        # For every mutation per individual, up to the limit
        for _ in 1:adv_config.mutants_per_individual_per_type
            drop_op::Bool = rand() < adv_config.p_lose_operation  && length(indiv.ops) > 0
            gain_op::Bool = rand() < adv_config.p_add_operation  && length(indiv.ops) < max_ops 
            swap_op::Bool = rand() < adv_config.p_swap_operation  && length(indiv.ops) > 0
            mutate_op::Bool = rand() < adv_config.p_mutate_operation && length(indiv.ops) > 0 

            new_individuals::Int = drop_op + gain_op + swap_op + mutate_op
            if new_individuals > 0
                j::UInt = 1
                new_individuals_vec::Array{Union{Nothing,Individual}} = [nothing for _ in 1:new_individuals]

                if drop_op
                    new_individuals_vec[j] = drop_op_with_constraints(indiv)
                    j+=1
                end
                if gain_op
                    new_individuals_vec[j] = gain_op_with_constraints(indiv, calibration_data, valid_qubits, purified_pairs,num_registers,adv_config.communication_fidelity_in) 
                    j+=1
                end
                if swap_op
                    new_individuals_vec[j] = swap_op_with_constraints(indiv)
                    j+=1
                end
                if mutate_op
                    new_individuals_vec[j] = mutate_with_constraints(indiv)
                end
 
                append!(indivContainers[i], new_individuals_vec)
            end
        end
    end

    ## add all children back to the individuals vector
    append!(individuals, vcat(indivContainers...))
end

"""
    step_with_constraints!(population::Population, max_ops::Int,  valid_qubits::Array{Int},purified_pairs::Int,num_registers::Int,calibration_data::Dict{Any,Any},advanced_config::AdvancedConfiguration)

    Important: This calls sort and cull.

    Execute one generation step of the genetic algorithm

    Each generation, the population is optimized by applying mutation, selection, and crossover operations
    while respecting constraints (such as gate connectivity and noise calibration)
"""
function step_with_constraints!(population::Population, max_ops::Int,  valid_qubits::Array{Int},purified_pairs::Int,num_registers::Int,calibration_data::Dict{Any,Any}, num_simulations::Int, optimize_for::CostFunction, advanced_config::AdvancedConfiguration)
    # Mark existing individuals as survivors
    # Survivors ensure that some individuals are carried over unchanged, maintaining good solutions
    for indiv in population.individuals
        indiv.history =  "survivor"
    end

    # Select pairs of parents randomly from the population
    parents = [(rand(population.individuals), rand(population.individuals)) for _ = 1:advanced_config.pairs]
    # Generate children from each pair of individuals
    for (p1, p2) in parents # every step, adds vertically to individual vector. May be a more efficient way TODO 
        population.individuals = vcat(population.individuals, [new_child(p1, p2,max_ops,num_registers) for j = 1:advanced_config.children_per_pair])
    end

    add_mutations!(population.individuals,valid_qubits,purified_pairs,num_registers,max_ops,calibration_data,advanced_config)
    
    # Sort the population by fitness and cull the excess individuals to maintain the population size
    simulate_and_sort!(population,num_simulations,purified_pairs,num_registers,optimize_for,advanced_config)
    cull!(population,advanced_config.population_size)

end

function sort_pop!(population::Population)
    # update the population with the sorted vector of individuals by fitness
    population.individuals = sort(population.individuals, by=indiv -> indiv.fitness, rev=true)
end

"""
    simulate_and_sort!(population::Population,num_simulations::Int,purified_pairs::Int,num_registers::Int,optimize_for::CostFunction,advanced_config::AdvancedConfiguration)

"" Evaluate and Sort the individuals in descending order of fitness 
"""
function simulate_and_sort!(population::Population,num_simulations::Int,purified_pairs::Int,num_registers::Int,optimize_for::CostFunction,advanced_config::AdvancedConfiguration)
    # calculate and update each individual's performance
    Threads.@threads for indiv in population.individuals
        calculate_performance!(indiv,
            num_simulations,
            purified_pairs,
            num_registers, 
            optimize_for, 
            advanced_config.code_distance,
            advanced_config.communication_fidelity_in) 
    end
    sort_pop!(population)
end

"""
    cull!(population::Population,population_size::Int)

    Reduce the population size to the target population_size 
"""
function cull!(population::Population,population_size::Int)
    population.individuals = population.individuals[1:population_size]
end


##########################################
##  (1) Parent Selection and Crossover   ##
##########################################

"""creates a new child individual from two parent individuals by combining their operations"""
function new_child(indiv::Individual, indiv2::Individual, max_ops::Int,num_registers::Int)::Individual
    if length(indiv2.ops) == 0
        return deepcopy(indiv) # No crossover if one of the parents has no operations
    elseif length(indiv.ops) == 0
        return deepcopy(indiv2) 
    end
    new_indiv = deepcopy(indiv)

    ## filiter out the BellSwap gate
    # method 1
    ops1 = filter(op -> op isa PauliNoiseBellGate{CNOTPerm} || op isa NoisyBellMeasureNoisyReset, indiv.ops)
    ops2 = filter(op -> op isa PauliNoiseBellGate{CNOTPerm} || op isa NoisyBellMeasureNoisyReset, indiv2.ops)

    # method 2
    # i = indiv.r*(indiv.r-1)/2         # number of NoisyBellSwap gate is r*(r-1)/2
    # ops1 = indiv.ops[i+1:end]
    # ops1 = indiv.ops[i+1:end]

    # explore new circuit configurations
    # """With a 50% probability, the operations list of the first parent is reversed"""
    if rand() < 0.5
        ops1 = ops1[end:-1:1]
    end

    # """With a 50% probability, the operations list of the second parent is reversed"""
    if rand() < 0.5
        ops2 = ops2[end:-1:1]
    end

    # Potential problem: This selection may not be optimal because some valuable operations could be at the end of the parents' lists,
    # but the child might miss them due to the biased selection from the beginning only.

    # Randomly selects how many operations to take from each parent
    sample1 = 1:min(length(ops1), max_ops)
    if length(sample1) == 0
        return deepcopy(indiv2)
    end   
    num_ops1 = rand(sample1)
    
    sample2 = 1:min(length(ops2), max_ops - num_ops1)
    if length(sample2) == 0
        return deepcopy(indiv)
    end
    num_ops2 = rand(sample2)

    # Combining the selected operations from both parents
    # TODO: fix types here?
    new_indiv.ops = vcat(
        # indiv.ops[1:num_registers*(num_registers-1)÷2],  # Include BellSwap gates
        ops1[1:num_ops1],
        ops2[1:num_ops2]
    )[1:min(end, max_ops)]  # Ensure the total number of operations does not exceed max_ops
    new_indiv.history =  "child"  # indicating that it was created through crossover

    return new_indiv

end


##########################################
##      (2) methods of mutations        ##
##########################################


"""
    swap_op_with_constraints(indiv::Individual)

    Swap two randomly selected operations, that can be cnot or measurement. Create new individual with these new ops and return it.
"""
function swap_op_with_constraints(indiv::Individual)::Individual
    if length(indiv.ops) < 2
        return indiv
    end
    new_indiv = deepcopy(indiv)
    ops = indiv.ops 

    """ Randomly select two positions (cnot or measurement) """
    sample = [i for i in 1:length(ops) if isa(ops[i],PauliNoiseBellGate{CNOTPerm} ) || isa(ops[i], NoisyBellMeasureNoisyReset)]

    if length(sample) < 2
        return indiv
    end
    
    ind1 = rand(sample)

    # do not swap with the same operation
    sample2 = [i for i in 1:length(ops) if i != ind1 && (isa(ops[i],PauliNoiseBellGate{CNOTPerm} ) || isa(ops[i], NoisyBellMeasureNoisyReset))]

    if length(sample2) < 1
        return indiv
    end
    ind2 = rand(sample2)

    """ Get the operations at the selected positions """
    op1, op2 = ops[ind1], ops[ind2]

    """ Swap the operations """
    
    new_indiv.ops[ind1] = op2
    new_indiv.ops[ind2] = op1
    new_indiv.history = "swap_m"

    return new_indiv
end

"""
    drop_op_with_constraints(indiv::Individual)::Individual

    mutating methods with gate connectivity constraints

    TODO
    make trait checking function
        "is_droppable"
        default is true, subtypes can make it false

    Define is_droppable for any object
    So you can use gates from other libraries etc

"""
function drop_op_with_constraints(indiv::Individual)::Individual
    new_indiv = deepcopy(indiv)
    # Filter the indices of operations that can be dropped
    drop_indices = [i for i in 1:length(new_indiv.ops) if is_droppable(new_indiv.ops[i])]

    # If there are no droppable operations, return the individual as is
    if  isempty(drop_indices)
        return new_indiv
    else
        # Randomly select and delete one of the droppable operations
        deleteat!(new_indiv.ops, rand(drop_indices))
        new_indiv.history = "drop_m"
        return new_indiv
    end
end

"""
    gain_op_with_constraints(indiv::Individual,  calibration_data::Dict, valid_qubits::Array{Int},purified_pairs,num_registers,f_in)::Individual

    Add a new operation to an individual copy, and return it
    
"""
function gain_op_with_constraints(indiv::Individual,  calibration_data::Dict, valid_qubits::Array{Int},purified_pairs,num_registers,f_in)::Individual
    new_indiv = deepcopy(indiv)

    """ Generate pairs in the order specified """
    valid_pairs = generate_valid_pairs(valid_qubits)
    valid_pair_A = valid_pairs[rand(1:length(valid_pairs))]

    """ Ensure the selected neighbor B is different from A """
    neighbors = []
    index_A = map_valid_qubits_to_pairs(valid_pair_A, valid_pairs )
    if index_A > 1
        push!(neighbors, valid_pairs[index_A - 1])
    end
    if index_A < length(valid_pairs)
        push!(neighbors, valid_pairs[index_A + 1])
    end
    valid_pair_B = neighbors[rand(1:length(neighbors))]

    # TODO: check if these numbers should be magic numbers or not?
    """ Create a random operation, with a 70% chance of being a PauliNoiseBellGate
        and a 30% chance of being a NoisyBellMeasureNoisyReset """
    rand_op = if rand() < 0.7

        pair_num_A = map_valid_qubits_to_pairs(valid_pair_A, valid_pairs)
        pair_num_B = map_valid_qubits_to_pairs(valid_pair_B, valid_pairs)
        random_gate = rand(CNOTPerm, (pair_num_A, pair_num_B)...)

        qubit_pair = map_num_to_valid_qubits(random_gate.idx1, valid_pairs)
        p2= calibration_data[qubit_pair][:two_qubit_error]
        PauliNoiseBellGate(random_gate, p2/3, p2/3, p2/3)

    else
        # q = rand(1:indiv.r)  # Randomly select a register
        q = rand(purified_pairs+1 : num_registers)          # debug
        q_pair = map_num_to_valid_qubits(q, valid_pairs)

        # Check if the register is the bottom-most pair. If so, create a raw Bell pair
        if q == num_registers
            # Definitely need to consider two qubit gate durations for high depth circuits.
            NoisyBellMeasureNoisyReset(
                rand(BellMeasure, q),
                pair_readout_error(calibration_data[q_pair[1]][:readout_error], calibration_data[q_pair[2]][:readout_error]),   #❓change and see difference
                f_in_to_pauli(f_in)...                              # only consider F_in, not consider T1, T2
                )
        else
        # Other registers are NoisyBellMeasureNoReset
            NoisyBellMeasureNoisyReset(
                rand(BellMeasure, q),
                pair_readout_error(calibration_data[q_pair[1]][:readout_error], calibration_data[q_pair[2]][:readout_error]),
                0.25, 0.25, 0.25                     # Trick of NoReset: create a maximally mixed state
            )
        end


        # TODO: can delete this commented out code?
        # """ Get T1 T2 noise parameters from calibration data for the specific qubits """
        # a, b = q_pair
        # noise_params = map_phase_damping_to_pauli(indiv.f_in,
        #     calibration_data[a][:t], calibration_data[a][:T1], calibration_data[a][:T2],
        #     calibration_data[b][:t], calibration_data[b][:T1], calibration_data[b][:T2]
        # )

        # if q == population.r
        #     NoisyBellMeasureNoisyReset(
        #         rand(BellMeasure, q),
        #         1 - pair_readout_error(calibration_data[q_pair[1]][:readout_error], calibration_data[q_pair[2]][:readout_error]),
        #         noise_params[1], noise_params[2], noise_params[3]

        #     )
        # else
        #     NoisyBellMeasureNoisyReset(
        #         rand(BellMeasure, q),
        #         1 - pair_readout_error(calibration_data[q_pair[1]][:readout_error], calibration_data[q_pair[2]][:readout_error]),
        #         # 0.25, 0.25, 0.25
        #         noise_params[1], noise_params[2], noise_params[3]
        #     )
        # end
    end


    if isempty(new_indiv.ops)
        push!(new_indiv.ops, rand_op)
    else
        position = 1
        if length(new_indiv.ops) == 2
            # edge case, ignore bellswaps and just add the gate somewhere
            position = rand((1,2))
        else
            """ Insert the random operation at a random position after BellSwap gates """
            # TODO: remove consideration for BellSwaps (to be removed)

            # Search through the ops to find the last consecutive noisybellswaps
            i = 1
            while i <= length(new_indiv.ops) && isa(new_indiv.ops[i], BellSwap)
                i += 1
            end
            # i = num_registers*(num_registers-1)/2         # number of NoisyBellSwap gate is r*(r-1)/2 // not necessarlly true, especially after running long simulations where noisybellswaps are dropped
            # position = rand(i + 1 : length(new_indiv.ops)) |> Int   # Ensures integer index
            position = rand(i: length(new_indiv.ops)) |> Int   # Ensures integer index

        end
        insert!(new_indiv.ops, position, rand_op)
    end

    new_indiv.history = "gain_m"
    return new_indiv
end

"""
    mutate_with_constraints(indiv::Individual)::Individual

    applying the appropriate mutation function to each of its gates, returns new individual
"""
function mutate_with_constraints(indiv::Individual)::Individual
    if length(indiv.ops) == 0
        return indiv
    end
    new_indiv = deepcopy(indiv)
    # Apply mutation only to PauliNoiseBellGate{CNOTPerm} and NoisyBellMeasureNoisyReset operations
    new_indiv.ops = [isa(gate, PauliNoiseBellGate{CNOTPerm}) || isa(gate, NoisyBellMeasureNoisyReset) ? mutate(gate) : gate for gate in new_indiv.ops]
    new_indiv.history =  "ops_m"
    return new_indiv
end

"""
    mutate(gate::NoisyBellMeasureNoisyReset)

 The measurement component (X,Y,Z) of the gate is randomized while keeping the other parameters (p, px, py, pz) the same 
"""
function mutate(gate::NoisyBellMeasureNoisyReset)
    return NoisyBellMeasureNoisyReset(rand(BellMeasure, gate.m.sidx), gate.p, gate.px, gate.py, gate.pz)
end

"""
    mutate(gate::PauliNoiseBellGate)

 The permutation component of the gate is randomized while keeping the noise parameters (px, py, pz) the same 
"""
function mutate(gate::PauliNoiseBellGate)
    return PauliNoiseBellGate(rand(CNOTPerm, gate.g.idx1, gate.g.idx2), gate.px, gate.py, gate.pz)
end


#######################################################################################
##     renoise bell pair for different input fidilities when generating dataframe    ##
#######################################################################################

"""
    refresh_noise(n::PauliNoiseBellGate, f_in::Float64)

    reapply the noise parameters without considering p2
    No Change to Noise Parameters
TBW

"""
function refresh_noise(n::PauliNoiseBellGate, f_in::Float64)
    return PauliNoiseBellGate(n.g, n.px, n.py, n.pz)
end

"""
    refresh_noise(n::NoisyBellMeasureNoisyReset, f_in::Float64)

    The noise parameters are updated based on f_in
"""
function refresh_noise(n::NoisyBellMeasureNoisyReset, f_in::Float64)
    px, py, pz = f_in_to_pauli(f_in)
    return NoisyBellMeasureNoisyReset(n.m, n.p, px, py, pz)
end

"""
    refresh_noise(indiv::Individual, f_in::Float64)

    Reset and return an individual's performance and fitness, and refresh the noise of their operations. 
"""
function refresh_noise(indiv::Individual, f_in::Float64)
    new_indiv = deepcopy(indiv)
    new_indiv.ops = [refresh_noise(op, f_in) for op in indiv.ops]
    new_indiv.performance = Performance([], 0, 0, 0, 0)
    new_indiv.fitness = 0.0
    return new_indiv
end

# TODO: get docs
"""
    total_raw_pairs(ops,num_registers,purified_pairs)

TBW
"""
function total_raw_pairs(ops,num_registers,purified_pairs)
    total = num_registers
    last_ops_reg = Set(1:num_registers)
    # TODO: fix the algorithm to not need to use deleition to count
    for op in reverse(ops)
        if isa(op, NoisyBellMeasureNoisyReset)
            t = op.m.sidx
            if t in last_ops_reg
                delete!(last_ops_reg, t)
                if t < purified_pairs
                    total += 1
                end
            else
                total += 1
            end
        elseif isa(op, PauliNoiseBellGate{CNOTPerm})
            for t in [op.g.idx1, op.g.idx2]
                delete!(last_ops_reg, t)
            end
        elseif isa(op, BellSwap)
            for t in [op.sidx1, op.sidx2]
                delete!(last_ops_reg, t)
            end
        end
    end

    return total
end

# TODO: check changed hash functionality
function Base.hash(indiv::Individual)
    return hash((Individual, [hash(op) for op in indiv.ops]))
end

function Base.hash(g::CNOTPerm)
    return hash((CNOTPerm, g.single1, g.single2, g.idx1, g.idx2))
end

function Base.hash(n::PauliNoiseBellGate{CNOTPerm})
    return hash((PauliNoiseBellGate{CNOTPerm}, n.px, n.py, n.pz))
end

function Base.hash(m::BellMeasure)
    return hash((BellMeasure, m.midx, m.sidx))
end

function Base.hash(n::NoisyBellMeasureNoisyReset)
    return hash((NoisyBellMeasureNoisyReset, hash(n.m), n.p, n.px, n.py, n.pz))
end


end