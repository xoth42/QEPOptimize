using Test

# include("test_Configurable.jl")
# include("test_Optimizer.jl")

# using Test
using CSV
using DataFrames
using QEPO.Configurable

# Mock data for testing
mock_calibration_data = """
Qubit,T1 (us),T2 (us),Readout assignment error ,CNOT error ,Gate time (ns)
0,50,60,0.01,"0_1:0.02","0_1:200"
1,55,65,0.02,"1_0:0.03","1_0:210"
"""

# Write mock data to a temporary file
mock_calibration_file = tempname() * ".csv"
open(mock_calibration_file, "w") do f
    write(f, mock_calibration_data)
end

# Create a mock HardwareConfiguration
mock_hardware_config = HardwareConfiguration(mock_calibration_file,Dict(),([1,2]))
set_calibration_data_path!(mock_hardware_config,mock_calibration_file)
mock_config = DEFAULT_CONFIGURATION(mock_hardware_config,DEFAULT_ADVANCED_CONFIG())

# Define a small tolerance for floating-point comparisons
tolerance = 1e-5

# Test the update_hardware_data! function
@testset "update_hardware_data!" begin
    Configurable.update_hardware_data!(mock_config)
    hardware_config = Configurable.get_hardware_config(mock_config)
    calibration_data = Configurable.get_calibration_data(hardware_config)
    
    @test isapprox(calibration_data[0].T1, 50e-6, atol=tolerance)
    @test isapprox(calibration_data[0].T2, 60e-6, atol=tolerance)
    @test isapprox(calibration_data[0].readout_error, 0.01, atol=tolerance)
    
    @test isapprox(calibration_data[1].T1, 55e-6, atol=tolerance)
    @test isapprox(calibration_data[1].T2, 65e-6, atol=tolerance)
    @test isapprox(calibration_data[1].readout_error, 0.02, atol=tolerance)
    
    @test isapprox(calibration_data[(0, 1)].two_qubit_error, 0.02, atol=tolerance)
    @test isapprox(calibration_data[(0, 1)].gate_time, 200e-9, atol=tolerance)
    
    @test isapprox(calibration_data[(1, 0)].two_qubit_error, 0.03, atol=tolerance)
    @test isapprox(calibration_data[(1, 0)].gate_time, 210e-9, atol=tolerance)
end

@testset "QEPO Library Tests" begin
    @testset "Test Case 1" begin
        @test 1 + 1 == 2
    end

    @testset "Test Case 2" begin
        @test 2 * 2 == 4
    end

    @testset "Test Case 3" begin
        @test 3 - 1 == 2
    end
end


using Test
# include("../src/Optimizer.jl")
using QEPO.Optimizer
using QEPO.Configurable
# using QEPO.Optimizer: generate_noisy_BellSwap_ops_for_individual, long_range_entanglement_generation!, NoisyBellSwap, Population, get_individuals, get_ops

# Mock data for testing
mutable struct MockCalibrationData
    two_qubit_error::Float64
end

mutable struct MockConfig
    num_registers::Int
    starting_ops::Int
    valid_qubits::Vector{Int}
    calibration_data::Dict{Tuple{Int, Int}, MockCalibrationData}
end

function get_num_registers(config::MockConfig)
    return config.num_registers
end

function get_starting_ops(config::MockConfig)
    return config.starting_ops
end

function get_valid_qubits(config::MockConfig)
    return config.valid_qubits
end

function get_calibration_data(config::MockConfig)
    return config.calibration_data
end



# Test generate_noisy_BellSwap_ops_for_individual
# @testset "generate_noisy_BellSwap_ops_for_individual" begin
    valid_qubits = ([1,2,3,4])
    num_registers = 4
    valid_pairs = generate_valid_pairs(valid_qubits)
    starting_ops = 5

    calibration_data = Dict(
        (0,0) => MockCalibrationData(0.01),
        (1, 2) => MockCalibrationData(0.01),
        (2, 1) => MockCalibrationData(0.02),
        (1, 0) => MockCalibrationData(0.03),
        (0, 1) => MockCalibrationData(0.03),
        (3,4)=>MockCalibrationData(0.01),
        (2,3)=>MockCalibrationData(0.01),
    )

    noisy_ops = generate_noisy_BellSwap_ops_for_individual(num_registers, valid_pairs, starting_ops,calibration_data)

    @test length(noisy_ops) > 0
    noisy_ops
# end
# Test long_range_entanglement_generation!
# @testset "long_range_entanglement_generation!" begin
    population = Population([], Dict())
    # config = MockConfig(
    #     4,
    #     5,
       
        
    # )
    valid_qubits =  [1, 2, 3, 4]
    calibration_data = Dict(
        (0,0) => MockCalibrationData(0.01),
        (1, 2) => MockCalibrationData(0.01),
        (2, 1) => MockCalibrationData(0.02),
        (1, 0) => MockCalibrationData(0.03),
        (0, 1) => MockCalibrationData(0.03),
        (3,4)=>MockCalibrationData(0.01),
        (2,3)=>MockCalibrationData(0.01),
    )
    hardware_conf = HardwareConfiguration("",calibration_data,valid_qubits)
    adv_conf = DEFAULT_ADVANCED_CONFIG()
    set_starting_ops!(adv_conf,5)

    config = DEFAULT_CONFIGURATION(hardware_conf,adv_conf)
    set_num_registers!(config,4)
    long_range_entanglement_generation!(population, config)
    indiv = get_individuals(population)
    get_ops(indiv[1])
    @test length(get_individuals(population)) > 0
    @test all(indiv -> length(get_ops(indiv)) > 0, get_individuals(population))
# end