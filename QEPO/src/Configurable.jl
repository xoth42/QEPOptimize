module Configurable
# TODO: default configurations
# TODO: methods for creating configurations

using CSV 
using DataFrames

"""What the optimizer uses to define success """
@enum  CostFunction logical_qubit_fidelity purified_pairs_fidelity average_marginal_fidelity

"""All configuration needed for running optimizer. Hardware configuration and advanced configuration need not be touched by the user in most scenarios. """
struct Configuration
    # Common user defined/edited parameters
    num_simulations::Int        # number of simulations that the optimizer will run
    raw_bell_pairs::Int         # (n) incoming bell pairs
    purified_pairs::Int         # (k) outgoing bell pairs
    num_registers::Int          # amount of registers
    optimize_for::CostFunction  # optimization goal (see @enum CostFunction)
    max_gen::Int                # TODO more info
    max_ops::Int                # Limits the number of operations in each individual circuit
    hardware_config::HardwareConfiguration # for hardware specifications, error rates
    advanced_config::AdvancedConfiguration # for additional optimizer parameters that more or less can be left alone
end
"""Config specific to given hardware, used in the circuit simulation."""
struct HardwareConfiguration
    calibration_data_path::String # Path to the csv file containing IBM calibration data
    calibration_data::Dict{Any,Any}
    # additional overrides?
end

"""Additional optimizer configurations. Normally can be left alone."""
struct AdvancedConfiguration
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


"""
read the calibration data of IBM quantum computer:
Reads the calibration data of an IBM quantum computer from a CSV file.
It parses the T1, T2, readout error, and two-qubit gate error and time data.
Returns a dictionary containing the calibration data.
"""
function read_calibration_data(filename::String)

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



end
