module Configurable
# TODO: default configurations
# TODO: methods for creating configurations

using CSV 
using DataFrames


export AbstractAdvancedConfiguration, AbstractConfiguration, AbstractHardwareConfiguration, HardwareConfiguration, AdvancedConfiguration, HardwareConfiguration, DEFAULT_ADVANCED_CONFIG, DEFAULT_CONFIGURATION, DEFAULT_HARDWARE_CONFIG

"""What the optimizer uses to define success """
@enum  CostFunction logical_qubit_fidelity purified_pairs_fidelity average_marginal_fidelity

"""Type for the main configuration"""
abstract type AbstractConfiguration end
abstract type AbstractHardwareConfiguration end
abstract type AbstractAdvancedConfiguration end

"""Config specific to given hardware, used in the circuit simulation."""
struct HardwareConfiguration <: AbstractHardwareConfiguration
    calibration_data_path::String # Path to the csv file containing IBM calibration data
    calibration_data::Dict{Any,Any}
    # additional overrides?
end

"""Additional optimizer configurations. Normally can be left alone."""
struct AdvancedConfiguration <: AbstractAdvancedConfiguration
    # For the Genetic optimizer:
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


"""All configuration needed for running optimizer. Hardware configuration and advanced configuration need not be touched by the user in most scenarios. """
struct Configuration <: AbstractConfiguration
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



"""DEFAULT CONFIGS 
TODO: store somewhere else
Make sure these work well as the defaults
"""
function DEFAULT_HARDWARE_CONFIG()::AbstractHardwareConfiguration
    return HardwareConfiguration(
    calibration_data_path = "ibm_sherbrooke_calibrations_2024-10-09.csv",
    calibration_data = read_calibration_data("ibm_sherbrooke_calibrations_2024-10-09.csv")
)
end

function DEFAULT_ADVANCED_CONFIG()::AbstractAdvancedConfiguration
    return AdvancedConfiguration(
    starting_ops = 17,
    pairs = 20,
    children_per_pair = 3,
    mutants_per_individual_per_type = 5,
    p_lose_operation = 0.9,
    p_add_operation = 0.7,
    p_swap_operations = 0.8,
    p_mutate_operations = 0.8
)
end



function DEFAULT_CONFIGURATION(hardware_conf::AbstractHardwareConfiguration,advanced_conf::AbstractAdvancedConfiguration)::AbstractConfiguration
    return Configuration(
    num_simulations = 1000,
    raw_bell_pairs = 6,
    purified_pairs = 1,
    num_registers = 4,
    optimize_for = CostFunction.purified_pairs_fidelity,
    max_gen = 20,
    max_ops = 17,
    hardware_config = hardware_conf,
    advanced_config = advanced_conf
)
end


function DEFAULT_CONFIGURATION()::AbstractConfiguration
    return DEFAULT_CONFIGURATION(DEFAULT_HARDWARE_CONFIG(),DEFAULT_ADVANCED_CONFIG())
end


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

function update_hardware_data!(config::Configuration)
    hardware_config::HardwareConfiguration = get_hardware_config(config)
    set_calibration_data!(hardware_config, read_calibration_data(get_calibration_data_path(hardware_config)))
    set_hardware_config!(config,hardware_config)
end



##### Configuration struct
# Getters
function get_num_simulations(config::Configuration)::Int
    return config.num_simulations
end

function get_raw_bell_pairs(config::Configuration)::Int
    return config.raw_bell_pairs
end

function get_purified_pairs(config::Configuration)::Int
    return config.purified_pairs
end

function get_num_registers(config::Configuration)::Int
    return config.num_registers
end

function get_optimize_for(config::Configuration)::CostFunction
    return config.optimize_for
end

function get_max_gen(config::Configuration)::Int
    return config.max_gen
end

function get_max_ops(config::Configuration)::Int
    return config.max_ops
end

function get_hardware_config(config::Configuration)::HardwareConfiguration
    return config.hardware_config
end

function get_advanced_config(config::Configuration)::AdvancedConfiguration
    return config.advanced_config
end

# Setters
function set_num_simulations!(config::Configuration, value::Int)::Void
    config.num_simulations = value
end

function set_raw_bell_pairs!(config::Configuration, value::Int)::Void
    config.raw_bell_pairs = value
end

function set_purified_pairs!(config::Configuration, value::Int)::Void
    config.purified_pairs = value
end

function set_num_registers!(config::Configuration, value::Int)::Void
    config.num_registers = value
end

function set_optimize_for!(config::Configuration, value::CostFunction)::Void
    config.optimize_for = value
end

function set_max_gen!(config::Configuration, value::Int)::Void
    config.max_gen = value
end

function set_max_ops!(config::Configuration, value::Int)::Void
    config.max_ops = value
end

function set_hardware_config!(config::Configuration, value::HardwareConfiguration)::Void
    config.hardware_config = value
end

function set_advanced_config!(config::Configuration, value::AdvancedConfiguration)::Void
    config.advanced_config = value
end

#### HardwareConfiguration struct

# Getters
function get_calibration_data_path(config::HardwareConfiguration)::String
    return config.calibration_data_path
end

function get_calibration_data(config::HardwareConfiguration)::Dict{Any,Any}
    return config.calibration_data
end

# Setters
function set_calibration_data_path!(config::HardwareConfiguration, value::String)::Void
    config.calibration_data_path = value
end

function set_calibration_data!(config::HardwareConfiguration, value::Dict{Any,Any})::Void
    config.calibration_data = value
end

### AdvancedConfiguration struct

# Getters
function get_starting_ops(config::AdvancedConfiguration)::Int
    return config.starting_ops
end

function get_pairs(config::AdvancedConfiguration)::Int
    return config.pairs
end

function get_children_per_pair(config::AdvancedConfiguration)::Int
    return config.children_per_pair
end

function get_mutants_per_individual_per_type(config::AdvancedConfiguration)::Int
    return config.mutants_per_individual_per_type
end

function get_p_lose_operation(config::AdvancedConfiguration)::Float64
    return config.p_lose_operation
end

function get_p_add_operation(config::AdvancedConfiguration)::Float64
    return config.p_add_operation
end

function get_p_swap_operations(config::AdvancedConfiguration)::Float64
    return config.p_swap_operations
end

function get_p_mutate_operations(config::AdvancedConfiguration)::Float64
    return config.p_mutate_operations
end

# Setters
function set_starting_ops!(config::AdvancedConfiguration, value::Int)::Void
    config.starting_ops = value
end

function set_pairs!(config::AdvancedConfiguration, value::Int)::Void
    config.pairs = value
end

function set_children_per_pair!(config::AdvancedConfiguration, value::Int)::Void
    config.children_per_pair = value
end

function set_mutants_per_individual_per_type!(config::AdvancedConfiguration, value::Int)::Void
    config.mutants_per_individual_per_type = value
end

function set_p_lose_operation!(config::AdvancedConfiguration, value::Float64)::Void
    config.p_lose_operation = value
end

function set_p_add_operation!(config::AdvancedConfiguration, value::Float64)::Void
    config.p_add_operation = value
end

function set_p_swap_operations!(config::AdvancedConfiguration, value::Float64)::Void
    config.p_swap_operations = value
end

function set_p_mutate_operations!(config::AdvancedConfiguration, value::Float64)::Void
    config.p_mutate_operations = value
end

end