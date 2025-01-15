# TODO: default configurations
# TODO: methods for creating configurations
# TODO: delete NOTES
# TODO: document all functions, create docs
# TODO: store default configs somewhere else?
# TODO: Make sure these work well as the defaults


"""NOTES: 
# function abc(a::AdvancedConfiguration)

# abstract type AbstractConfiguration

# function optimizer_config(conf::AbstractConfiguration)
#     conf.prop # this is bad
#     getprop(conf) # this will compose well
#     return ...
# end

# function mutate(individual)
#     Individual([mutate(g) for g in individual.gates])
# end

# function mutate(m::AbstractMeasurement)
#     typeof(m)(affectedqubit(m)+1)
# end
"""

module Configurable

    using CSV 
    using DataFrames

    """What the optimizer uses to define success """
    @enum  CostFunction logical_qubit_fidelity purified_pairs_fidelity average_marginal_fidelity
    
    export CostFunction
    export AbstractConfiguration, AbstractHardwareConfiguration, AbstractAdvancedConfiguration
    export HardwareConfiguration, AdvancedConfiguration, Configuration
    export get_num_simulations, get_raw_bell_pairs, get_purified_pairs, get_num_registers, get_optimize_for, get_max_gen, get_max_ops, get_hardware_config, get_advanced_config
    export set_num_simulations!, set_raw_bell_pairs!, set_purified_pairs!, set_num_registers!, set_optimize_for!, set_max_gen!, set_max_ops!, set_hardware_config!, set_advanced_config!
    export get_calibration_data_path, get_calibration_data, get_valid_qubits
    export set_calibration_data_path!, set_calibration_data!, set_valid_qubits!
    export get_population_size, get_starting_pop_multiplier, get_starting_ops, get_pairs, get_children_per_pair, get_mutants_per_individual_per_type, get_p_lose_operation, get_p_add_operation, get_p_swap_operations, get_p_mutate_operations
    export set_population_size!, set_starting_pop_multiplier!, set_starting_ops!, set_pairs!, set_children_per_pair!, set_mutants_per_individual_per_type!, set_p_lose_operation!, set_p_add_operation!, set_p_swap_operations!, set_p_mutate_operations!
    export get_code_distance, get_communication_fidelity_in, set_code_distance!, set_communication_fidelity_in!
    export DEFAULT_HARDWARE_CONFIG, DEFAULT_ADVANCED_CONFIG, DEFAULT_CONFIGURATION, read_calibration_data, update_hardware_data!


    """Type for the main configuration"""
    abstract type AbstractConfiguration end
    abstract type AbstractHardwareConfiguration end
    abstract type AbstractAdvancedConfiguration end

    """
        HardwareConfiguration

    Configuration specific to given hardware, used in the circuit simulation.

    # Fields
    - `calibration_data_path::String`: Path to the CSV file containing IBM calibration data.
    - `calibration_data::Dict{Any,Any}`: Dictionary containing the calibration data.
    - `valid_qubits::Array{Int}`: Array of integers representing the valid qubits."""
    mutable struct HardwareConfiguration <: AbstractHardwareConfiguration
        calibration_data_path::String # Path to the csv file containing IBM calibration data
        calibration_data::Dict{Any,Any}
        valid_qubits::Array{Int}
        # additional overrides?
    end

    """
        AdvancedConfiguration

    Additional optimizer configurations. Normally can be left alone.

    # Fields
    - `code_distance::Int`: For the Genetic optimizer, used for logical_qubit_fidelity.
    - `communication_fidelity_in::Float64`: Fidelity that is being generated at communication qubits. 0.9 would be appropriate (f_in).
    - `population_size::Int`: Target number of individuals in the population for each generation after initialization.
    - `starting_pop_multiplier::Int`: A multiplier used to determine the size of the initial population.
    - `starting_ops::Int`: Initial number of operations in each individual circuit of the population.
    - `pairs::Int`: Number of parent pairs selected for breeding new individuals.
    - `children_per_pair::Int`: Number of offspring produced by each pair of parents.
    - `mutants_per_individual_per_type::Int`: Number of mutations applied to each individual for each mutation type.
    - `p_lose_operation::Float64`: Probability of losing an operation during optimization.
    - `p_add_operation::Float64`: Probability of adding an operation during optimization.
    - `p_swap_operations::Float64`: Probability of swapping operations during optimization.
    - `p_mutate_operations::Float64`: Probability of mutating operations during optimization."""
    mutable struct AdvancedConfiguration <: AbstractAdvancedConfiguration
        # For the Genetic optimizer:
        # TODO: pair code_distance with Configuration.CostFunction?
        code_distance::Int                  # for logical_qubit_fidelity
        communication_fidelity_in::Float64                         # fidelity that being generated at communication qubits. 0.9 would be appropriate (f_in)
        population_size::Int                            #  target number of individuals in the population for each generation after initialization
        starting_pop_multiplier::Int                    #  A multiplier used to determine the size of the initial population
        starting_ops::Int                               # initial number of operations in each individual circuit of the population
        pairs::Int                                      # Number of parent pairs selected for breeding new individuals
        children_per_pair::Int                          # Number of offspring produced by each pair of parents
        mutants_per_individual_per_type::Int            # Number of mutations applied to each individual for each mutation type
        # Probabilities for various operations to be done during optimization:
        p_lose_operation::Float64
        p_add_operation::Float64
        p_swap_operations::Float64
        p_mutate_operations::Float64
    end

    """
        Configuration

    All configuration needed for running optimizer. Hardware configuration and advanced configuration need not be touched by the user in most scenarios.

    # Fields
    - `num_simulations::Int`: Number of simulations that the optimizer will run.
    - `raw_bell_pairs::Int`: Number of incoming bell pairs (n).
    - `purified_pairs::Int`: Number of outgoing bell pairs (k).
    - `num_registers::Int`: Amount of registers.
    - `optimize_for::CostFunction`: Optimization goal (see `@enum CostFunction`).
    - `max_gen::Int`: Maximum number of generations (TODO: more info needed).
    - `max_ops::Int`: Limits the number of operations in each individual circuit.
    - `hardware_config::AbstractHardwareConfiguration`: Hardware specifications, error rates.
    - `advanced_config::AbstractAdvancedConfiguration`: Additional optimizer parameters that more or less can be left alone."""
    mutable struct Configuration <: AbstractConfiguration
        # Common user defined/edited parameters
        num_simulations::Int        # number of simulations that the optimizer will run
        raw_bell_pairs::Int         # (n) incoming bell pairs
        purified_pairs::Int         # (k) outgoing bell pairs
        num_registers::Int          # amount of registers
        optimize_for::CostFunction  # optimization goal (see @enum CostFunction)
        max_gen::Int                # TODO more info
        max_ops::Int                # Limits the number of operations in each individual circuit
        hardware_config::AbstractHardwareConfiguration # for hardware specifications, error rates
        advanced_config::AbstractAdvancedConfiguration # for additional optimizer parameters that more or less can be left alone
    end



    """
    #######
    DEFAULT CONFIGS 
    #######
    """

    """
        DEFAULT_HARDWARE_CONFIG()::AbstractHardwareConfiguration

    Returns the default hardware configuration for the system. This includes the path to the calibration data and the calibration data itself.

    # Returns
    - `AbstractHardwareConfiguration`: The default hardware configuration.
    """
    function DEFAULT_HARDWARE_CONFIG()::AbstractHardwareConfiguration
        #=
        calibration_data_path = "ibm_sherbrooke_calibrations_2024-10-09.csv",
        calibration_data = read_calibration_data("ibm_sherbrooke_calibrations_2024-10-09.csv")
        valid_qubits = ([ 43, 44, 45, 46, 47,  48, 49, 50])

        =#
        return HardwareConfiguration("ibm_sherbrooke_calibrations_2024-10-09.csv", read_calibration_data("ibm_sherbrooke_calibrations_2024-10-09.csv"),([ 43, 44, 45, 46, 47,  48, 49, 50])
        )
    end

    """
        DEFAULT_ADVANCED_CONFIG()::AbstractAdvancedConfiguration

    Returns the default advanced configuration for the system. This includes various parameters for operations, pairs, children per pair, mutants per individual per type, and probabilities for different operations.

    # Returns
    - `AbstractAdvancedConfiguration`: The default advanced configuration.
    """
    function DEFAULT_ADVANCED_CONFIG()::AbstractAdvancedConfiguration
        #=
        code_distance = 3                    # for logical_qubit_fidelity
        communication_fidelity_in = 0.90                          # fidelity that being generated at communication qubits 
        population_size = 20                 # target population after distillation 
    starting_pop_multiplier = 200        # multiplier 20 or 200
        starting_ops = 17,
        
        pairs = 20,
        children_per_pair = 3,
        mutants_per_individual_per_type = 5,
        p_lose_operation = 0.9,
        p_add_operation = 0.7,
        p_swap_operations = 0.8,
        p_mutate_operations = 0.8
        =#
        return AdvancedConfiguration(3,.9,20,200,17,20,3,5,0.9,0.7,0.8,0.8)
    end



    """
        DEFAULT_CONFIGURATION(hardware_conf::AbstractHardwareConfiguration, advanced_conf::AbstractAdvancedConfiguration)::AbstractConfiguration

    Returns the default configuration for the system using the provided hardware and advanced configurations. This includes the number of simulations, raw and purified bell pairs, number of registers, optimization target, maximum generations, and maximum operations.

    # Arguments
    - `hardware_conf::AbstractHardwareConfiguration`: The hardware configuration to use.
    - `advanced_conf::AbstractAdvancedConfiguration`: The advanced configuration to use.

    # Returns
    - `AbstractConfiguration`: The default configuration.
    """
    function DEFAULT_CONFIGURATION(hardware_conf::AbstractHardwareConfiguration,advanced_conf::AbstractAdvancedConfiguration)::AbstractConfiguration
        #=
        num_simulations = 1000,
        raw_bell_pairs = 6,
        purified_pairs = 1,
        num_registers = 4,
        optimize_for = CostFunction.purified_pairs_fidelity,
        max_gen = 20,
        max_ops = 17,
        hardware_config = hardware_conf,
        advanced_config = advanced_conf
        =#

        return Configuration(1000,6, 1,4,CostFunction(0),20,17,hardware_conf,advanced_conf)
    end


    """
        DEFAULT_CONFIGURATION()::AbstractConfiguration

    Returns the default configuration for the system using the default hardware and advanced configurations.

    # Returns
    - `AbstractConfiguration`: The default configuration.
    """
    function DEFAULT_CONFIGURATION()::AbstractConfiguration
        return DEFAULT_CONFIGURATION(DEFAULT_HARDWARE_CONFIG(),DEFAULT_ADVANCED_CONFIG())
    end




    """
    read the calibration data of IBM quantum computer:
    Reads the calibration data of an IBM quantum computer from a CSV file.
    It parses the T1, T2, readout error, and two-qubit gate error and time data.
    Returns a dictionary containing the calibration data.
    """
    function read_calibration_data(filename::String)::Dict{Any,Any}

        # Read the CSV file into a DataFrame
        df = CSV.read(filename, DataFrame)

        # Initialize an empty dictionary to store the calibration data
        calibration_data::Dict{Any,Any} = Dict()
        # t = 400e-9                                        # TODO Two-qubit gate time (simplified)

        # Process the single qubit data
        for row in eachrow(df)
            qubit = row["Qubit"]
            T1 = row["T1 (us)"] * 1e-6
            T2 = row["T2 (us)"] * 1e-6
            readout_error = row["Readout assignment error "]
            calibration_data[qubit] = (T1=T1, T2=T2, readout_error=readout_error)
        end

        # Process the two-qubit data
        for row in eachrow(df)
            cnot_errors_str = row["CNOT error "]
            cnot_errors_list = split(cnot_errors_str, ';')
            cnot_time_str = row["Gate time (ns)"]
            cnot_time_list = split(cnot_errors_str, ';')

            for cnot_error in cnot_errors_list
                qubit_pair, cx_error_rates = split(cnot_error, ':')
                qubits = Tuple(map(x -> parse(Int, x), split(qubit_pair, '_')))
                cnot_error_rate = parse(Float64, cx_error_rates)
                calibration_data[qubits] = (two_qubit_error=cnot_error_rate,)
            end

            # Add when considering T1 T2 noise model
            for cnot_time in cnot_time_list
                qubit_pair, cx_gate_time = split(cnot_time, ':')
                qubits = Tuple(map(x -> parse(Int, x), split(qubit_pair, '_')))
                cnot_gate_time = parse(Float64, cx_gate_time)* 1e-9 # convert to seconds
                # Combining Two-Qubit Error and Gate Time
                existing_data = calibration_data[qubits]
                calibration_data[qubits] = merge(existing_data, (gate_time=cnot_gate_time,))
            end

        end

        return calibration_data
    end

    """
        update_hardware_data!(config::Configuration)

    Updates the hardware configuration data within the given `config` object. 
    This function retrieves the current hardware configuration, reads the 
    calibration data from the appropriate path, and sets the updated hardware 
    configuration back into the `config` object.

    # Arguments
    - `config::Configuration`: The configuration object to be updated.

    """
    function update_hardware_data!(config::Configuration)
        hardware_config::HardwareConfiguration = get_hardware_config(config)
        set_calibration_data!(hardware_config, read_calibration_data(get_calibration_data_path(hardware_config)))
        set_hardware_config!(config,hardware_config)
    end



    ##### Configuration struct

    # Getters

    """
        get_num_simulations(config::Configuration)::Int

    Returns the number of simulations from the given `Configuration` object.
    """
    function get_num_simulations(config::Configuration)::Int
        return config.num_simulations
    end

    """
        get_raw_bell_pairs(config::Configuration)::Int

    Returns the number of raw Bell pairs from the given `Configuration` object.
    """
    function get_raw_bell_pairs(config::Configuration)::Int
        return config.raw_bell_pairs
    end

    """
        get_purified_pairs(config::Configuration)::Int

    Returns the number of purified pairs from the given `Configuration` object.
    """
    function get_purified_pairs(config::Configuration)::Int
        return config.purified_pairs
    end

    """
        get_num_registers(config::Configuration)::Int

    Returns the number of registers from the given `Configuration` object.
    """
    function get_num_registers(config::Configuration)::Int
        return config.num_registers
    end

    """
        get_optimize_for(config::Configuration)::CostFunction

    Returns the cost function to optimize for from the given `Configuration` object.
    """
    function get_optimize_for(config::Configuration)::CostFunction
        return config.optimize_for
    end

    """
        get_max_gen(config::Configuration)::Int

    Returns the maximum number of generations from the given `Configuration` object.
    """
    function get_max_gen(config::Configuration)::Int
        return config.max_gen
    end

    """
        get_max_ops(config::Configuration)::Int

    Returns the maximum number of operations from the given `Configuration` object.
    """
    function get_max_ops(config::Configuration)::Int
        return config.max_ops
    end

    """
        get_hardware_config(config::Configuration)::HardwareConfiguration

    Returns the hardware configuration from the given `Configuration` object.
    """
    function get_hardware_config(config::Configuration)::HardwareConfiguration
        return config.hardware_config
    end

    """
        get_advanced_config(config::Configuration)::AdvancedConfiguration

    Returns the advanced configuration from the given `Configuration` object.
    """
    function get_advanced_config(config::Configuration)::AdvancedConfiguration
        return config.advanced_config
    end

    # Setters

    """
        set_num_simulations!(config::Configuration, value::Int)

    Sets the number of simulations in the given `Configuration` object to `value`.
    """
    function set_num_simulations!(config::Configuration, value::Int)
        config.num_simulations = value
    end

    """
        set_raw_bell_pairs!(config::Configuration, value::Int)

    Sets the number of raw Bell pairs in the given `Configuration` object to `value`.
    """
    function set_raw_bell_pairs!(config::Configuration, value::Int)
        config.raw_bell_pairs = value
    end

    """
        set_purified_pairs!(config::Configuration, value::Int)

    Sets the number of purified pairs in the given `Configuration` object to `value`.
    """
    function set_purified_pairs!(config::Configuration, value::Int)
        config.purified_pairs = value
    end

    """
        set_num_registers!(config::Configuration, value::Int)

    Sets the number of registers in the given `Configuration` object to `value`.
    """
    function set_num_registers!(config::Configuration, value::Int)
        config.num_registers = value
    end

    """
        set_optimize_for!(config::Configuration, value::CostFunction)

    Sets the cost function to optimize for in the given `Configuration` object to `value`.
    """
    function set_optimize_for!(config::Configuration, value::CostFunction)
        config.optimize_for = value
    end

    """
        set_max_gen!(config::Configuration, value::Int)

    Sets the maximum number of generations in the given `Configuration` object to `value`.
    """
    function set_max_gen!(config::Configuration, value::Int)
        config.max_gen = value
    end

    """
        set_max_ops!(config::Configuration, value::Int)

    Sets the maximum number of operations in the given `Configuration` object to `value`.
    """
    function set_max_ops!(config::Configuration, value::Int)
        config.max_ops = value
    end

    """
        set_hardware_config!(config::Configuration, value::HardwareConfiguration)

    Sets the hardware configuration in the given `Configuration` object to `value`.
    """
    function set_hardware_config!(config::Configuration, value::HardwareConfiguration)
        config.hardware_config = value
    end

    """
        set_advanced_config!(config::Configuration, value::AdvancedConfiguration)

    Sets the advanced configuration in the given `Configuration` object to `value`.
    """
    function set_advanced_config!(config::Configuration, value::AdvancedConfiguration)
        config.advanced_config = value
    end

    #### HardwareConfiguration struct

    # Getters

    """
        get_calibration_data_path(config::HardwareConfiguration)::String

    Returns the calibration data path from the given `HardwareConfiguration` object.
    """
    function get_calibration_data_path(config::HardwareConfiguration)::String
        return config.calibration_data_path
    end

    """
        get_calibration_data(config::HardwareConfiguration)::Dict{Any,Any}

    Returns the calibration data from the given `HardwareConfiguration` object.
    """
    function get_calibration_data(config::HardwareConfiguration)::Dict{Any,Any}
        return config.calibration_data
    end

    """
        get_valid_qubits(config::HardwareConfiguration)::Array{Int}

    Returns the array of valid qubits from the given HardwareConfiguration object.
    """
    function get_valid_qubits(config::HardwareConfiguration)::Array{Int}
        return config.valid_qubits
    end

    # Setters

    """
        set_calibration_data_path!(config::HardwareConfiguration, value::String)

    Sets the calibration data path in the given `HardwareConfiguration` object to `value`.
    """
    function set_calibration_data_path!(config::HardwareConfiguration, value::String)
        config.calibration_data_path = value
    end

    """
        set_calibration_data!(config::HardwareConfiguration, value::Dict{Any,Any})

    Sets the calibration data in the given `HardwareConfiguration` object to `value`.
    """
    function set_calibration_data!(config::HardwareConfiguration, value::Dict{Any,Any})
        config.calibration_data = value
    end

    """
        set_valid_qubits!(config::HardwareConfiguration, value::Array{Int})

    Sets the array of valid qubits in the given `HardwareConfiguration` object to `value`.
    """
    function set_valid_qubits!(config::HardwareConfiguration, value::Array{Int})
        config.valid_qubits = value
    end

    ### AdvancedConfiguration struct

    # Getters

    """
        get_population_size(config::AdvancedConfiguration)::Int

    Returns the population size from the given `AdvancedConfiguration` object.
    """
    function get_population_size(config::AdvancedConfiguration)::Int
        return config.population_size
    end

    """
        get_starting_pop_multiplier(config::AdvancedConfiguration)::Int

    Returns the starting population multiplier from the given `AdvancedConfiguration` object.
    """
    function get_starting_pop_multiplier(config::AdvancedConfiguration)::Int
        return config.starting_pop_multiplier
    end

    """
        get_starting_ops(config::AdvancedConfiguration)::Int

    Returns the starting operations from the given `AdvancedConfiguration` object.
    """
    function get_starting_ops(config::AdvancedConfiguration)::Int
        return config.starting_ops
    end

    """
        get_pairs(config::AdvancedConfiguration)::Int

    Returns the number of pairs from the given `AdvancedConfiguration` object.
    """
    function get_pairs(config::AdvancedConfiguration)::Int
        return config.pairs
    end

    """
        get_children_per_pair(config::AdvancedConfiguration)::Int

    Returns the number of children per pair from the given `AdvancedConfiguration` object.
    """
    function get_children_per_pair(config::AdvancedConfiguration)::Int
        return config.children_per_pair
    end

    """
        get_mutants_per_individual_per_type(config::AdvancedConfiguration)::Int

    Returns the number of mutants per individual per type from the given `AdvancedConfiguration` object.
    """
    function get_mutants_per_individual_per_type(config::AdvancedConfiguration)::Int
        return config.mutants_per_individual_per_type
    end

    """
        get_p_lose_operation(config::AdvancedConfiguration)::Float64

    Returns the probability of losing an operation from the given `AdvancedConfiguration` object.
    """
    function get_p_lose_operation(config::AdvancedConfiguration)::Float64
        return config.p_lose_operation
    end

    """
        get_p_add_operation(config::AdvancedConfiguration)::Float64

    Returns the probability of adding an operation from the given `AdvancedConfiguration` object.
    """
    function get_p_add_operation(config::AdvancedConfiguration)::Float64
        return config.p_add_operation
    end

    """
        get_p_swap_operations(config::AdvancedConfiguration)::Float64

    Returns the probability of swapping operations from the given `AdvancedConfiguration` object.
    """
    function get_p_swap_operations(config::AdvancedConfiguration)::Float64
        return config.p_swap_operations
    end

    """
        get_p_mutate_operations(config::AdvancedConfiguration)::Float64

    Returns the probability of mutating operations from the given `AdvancedConfiguration` object.
    """
    function get_p_mutate_operations(config::AdvancedConfiguration)::Float64
        return config.p_mutate_operations
    end

    # Setters

    """
        set_population_size!(config::AdvancedConfiguration, value::Int)

    Sets the population size in the given `AdvancedConfiguration` object to `value`.
    """
    function set_population_size!(config::AdvancedConfiguration, value::Int)
        config.population_size = value
    end

    """
        set_starting_pop_multiplier!(config::AdvancedConfiguration, value::Int)

    Sets the starting population multiplier in the given `AdvancedConfiguration` object to `value`.
    """
    function set_starting_pop_multiplier!(config::AdvancedConfiguration, value::Int)
        config.starting_pop_multiplier = value
    end

    """
        set_starting_ops!(config::AdvancedConfiguration, value::Int)

    Sets the starting operations in the given `AdvancedConfiguration` object to `value`.
    """
    function set_starting_ops!(config::AdvancedConfiguration, value::Int)
        config.starting_ops = value
    end

    """
        set_pairs!(config::AdvancedConfiguration, value::Int)

    Sets the number of pairs in the given `AdvancedConfiguration` object to `value`.
    """
    function set_pairs!(config::AdvancedConfiguration, value::Int)
        config.pairs = value
    end

    """
        set_children_per_pair!(config::AdvancedConfiguration, value::Int)

    Sets the number of children per pair in the given `AdvancedConfiguration` object to `value`.
    """
    function set_children_per_pair!(config::AdvancedConfiguration, value::Int)
        config.children_per_pair = value
    end

    """
        set_mutants_per_individual_per_type!(config::AdvancedConfiguration, value::Int)

    Sets the number of mutants per individual per type in the given `AdvancedConfiguration` object to `value`.
    """
    function set_mutants_per_individual_per_type!(config::AdvancedConfiguration, value::Int)
        config.mutants_per_individual_per_type = value
    end

    """
        set_p_lose_operation!(config::AdvancedConfiguration, value::Float64)

    Sets the probability of losing an operation in the given `AdvancedConfiguration` object to `value`.
    """
    function set_p_lose_operation!(config::AdvancedConfiguration, value::Float64)
        config.p_lose_operation = value
    end

    """
        set_p_add_operation!(config::AdvancedConfiguration, value::Float64)

    Sets the probability of adding an operation in the given `AdvancedConfiguration` object to `value`.
    """
    function set_p_add_operation!(config::AdvancedConfiguration, value::Float64)
        config.p_add_operation = value
    end

    """
        set_p_swap_operations!(config::AdvancedConfiguration, value::Float64)

    Sets the probability of swapping operations in the given `AdvancedConfiguration` object to `value`.
    """
    function set_p_swap_operations!(config::AdvancedConfiguration, value::Float64)
        config.p_swap_operations = value
    end

    """
        set_p_mutate_operations!(config::AdvancedConfiguration, value::Float64)

    Sets the probability of mutating operations in the given `AdvancedConfiguration` object to `value`.
    """
    function set_p_mutate_operations!(config::AdvancedConfiguration, value::Float64)
        config.p_mutate_operations = value
    end
    """
        get_code_distance(config::AdvancedConfiguration)::Int

    Returns the code distance from the given `AdvancedConfiguration` object.
    """
    function get_code_distance(config::AdvancedConfiguration)::Int
        return config.code_distance
    end

    """
        get_communication_fidelity_in(config::AdvancedConfiguration)::Float64

    Returns the communication fidelity in from the given `AdvancedConfiguration` object.
    """
    function get_communication_fidelity_in(config::AdvancedConfiguration)::Float64
        return config.communication_fidelity_in
    end

    """
        set_code_distance!(config::AdvancedConfiguration, value::Int)

    Sets the code distance in the given `AdvancedConfiguration` object to `value`.
    """
    function set_code_distance!(config::AdvancedConfiguration, value::Int)
        config.code_distance = value
    end

    """
        set_communication_fidelity_in!(config::AdvancedConfiguration, value::Float64)

    Sets the communication fidelity in the given `AdvancedConfiguration` object to `value`.
    """
    function set_communication_fidelity_in!(config::AdvancedConfiguration, value::Float64)
        config.communication_fidelity_in = value
    end
end