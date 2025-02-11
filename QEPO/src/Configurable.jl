
# TODO: methods for creating configurations
# TODO: delete NOTES?
# TODO: document all functions, create docs
# TODO: Make sure these work well as the defaults

"""NOTES: 
DON'T HESITATE TO EMAIL!
silly question, annoyances
https://modernjuliaworkflows.org/
https://www.youtube.com/watch?v=kc9HwsxE1OY


Priority:

    1.
    Done - Minimize the amount of state that you are moving around
    Make struct internals more standard
    Done - Restructure getters/setters -> 
    List what API is expected by optimizer
    1.1 Mutate APi

        Use any type for operations, get the set of functions needed for the optimizer to run, use multiple dispatch


    What does each gate need etc (mutate, droppable, etc) -> any new gate will work for optimizer

    Don't worry about error messages for methods etc (for now)


    1.2 Preliminary testing
        Have list of maual gates (~10)=> have code to evaulate quality of circuit


    2.
    Constructors & defaults
    @kwdef <--

    3.
    Documentation macro
    https://github.com/JuliaDocs/DocStringExtensions.jl





Later:
    3.1 Testing 




    Test macros 
    @test ?
    @testItems -> https://github.com/julia-vscode/TestItems.jl
    https://github.com/QuantumSavory/QuantumSymbolics.jl/pull/75

    4. Set up docs 
    4.1 Doc tests?


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
    
    export CostFunction, logical_qubit_fidelity, purified_pairs_fidelity, average_marginal_fidelity
    export AbstractConfiguration, AbstractHardwareConfiguration, AbstractAdvancedConfiguration
    export HardwareConfiguration, AdvancedConfiguration, Configuration
   
    export read_calibration_data, update_hardware_data!


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
        - `valid_qubits::Array{Int}`: Array of integers representing the valid qubits.
    """
    mutable struct HardwareConfiguration <: AbstractHardwareConfiguration
        calibration_data_path::String # Path to the csv file containing IBM calibration data
        calibration_data::Dict{Any,Any}
        valid_qubits::Array{Int}
        """
        ##Default constructor
        calibration_data_path = "ibm_sherbrooke_calibrations_2024-10-09.csv",
        calibration_data = read_calibration_data("ibm_sherbrooke_calibrations_2024-10-09.csv")
        valid_qubits = ([ 43, 44, 45, 46, 47,  48, 49, 50])"""
        HardwareConfiguration() = new("QEPO/data/ibm_sherbrooke_calibrations_2024-10-09.csv", read_calibration_data("QEPO/data/ibm_sherbrooke_calibrations_2024-10-09.csv"),([ 43, 44, 45, 46, 47,  48, 49, 50]))
        HardwareConfiguration(dataPath::String,valid_qubits::Array{Int}) = new(dataPath,read_calibration_data(dataPath),valid_qubits)
        
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
        - `p_swap_operation::Float64`: Probability of swapping operations during optimization.
        - `p_mutate_operation::Float64`: Probability of mutating operations during optimization.
    """
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
        p_swap_operation::Float64
        p_mutate_operation::Float64
        """
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
        p_swap_operation = 0.8,
        p_mutate_operation = 0.8"""
        AdvancedConfiguration() = new(3,.9,20,200,17,20,3,5,0.9,0.7,0.8,0.8)
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
        - `advanced_config::AbstractAdvancedConfiguration`: Additional optimizer parameters that more or less can be left alone.
    """
    mutable struct Configuration <: AbstractConfiguration
        """ Common user defined/edited parameters"""
        """number of simulations that the optimizer will run"""
        num_simulations::Int        

        raw_bell_pairs::Int         # (n) incoming bell pairs
        purified_pairs::Int         # (k) outgoing bell pairs
        num_registers::Int          # amount of registers
        optimize_for::CostFunction  # optimization goal (see @enum CostFunction)
        max_gen::Int                # TODO more info
        max_ops::Int                # Limits the number of operations in each individual circuit
        hardware_config::AbstractHardwareConfiguration # for hardware specifications, error rates
        advanced_config::AbstractAdvancedConfiguration # for additional optimizer parameters that more or less can be left alone
       """ 
        num_simulations = 1000,
        raw_bell_pairs = 6,
        purified_pairs = 1,
        num_registers = 4,
        optimize_for = CostFunction.purified_pairs_fidelity,
        max_gen = 20,
        max_ops = 17,
        hardware_config = hardware_conf,
        advanced_config = advanced_conf"""
        Configuration() =  new(1000,6, 1,4,CostFunction(1),20,17,HardwareConfiguration(),AdvancedConfiguration())
        Configuration(hardware_config::HardwareConfiguration,adv_config::AdvancedConfiguration) = new(1000,6, 1,4,CostFunction(1),20,17,hardware_config,adv_config)
    end

    """
    read the calibration data of IBM quantum computer:
    Reads the calibration data of an IBM quantum computer from a CSV file.
    It parses the T1, T2, readout error, and two-qubit gate error and time data.
    Returns a dictionary containing the calibration data."""
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
        config.hardware_config.calibration_data = read_calibration_data(config.hardware_config.calibration_data_path)
    end


end