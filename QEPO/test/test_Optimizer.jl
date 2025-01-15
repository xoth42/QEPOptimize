using Test
# include("../src/Optimizer.jl")
using QEPO.Optimizer
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
@testset "generate_noisy_BellSwap_ops_for_individual" begin
    num_registers = 4
    valid_pairs = [(1, 2), (2, 3), (3, 4)]
    starting_ops = 5
    calibration_data = Dict(
        (1, 2) => MockCalibrationData(0.01),
        (2, 3) => MockCalibrationData(0.02),
        (3, 4) => MockCalibrationData(0.03)
    )

    noisy_ops = generate_noisy_BellSwap_ops_for_individual(num_registers, valid_pairs, starting_ops)

    @test length(noisy_ops) > 0
    @test all(op -> isa(op, NoisyBellSwap), noisy_ops)
end

# Test long_range_entanglement_generation!
@testset "long_range_entanglement_generation!" begin
    population = Population([], Dict())
    config = MockConfig(
        4,
        5,
        [1, 2, 3, 4],
        Dict(
            (1, 2) => MockCalibrationData(0.01),
            (2, 3) => MockCalibrationData(0.02),
            (3, 4) => MockCalibrationData(0.03)
        )
    )

    long_range_entanglement_generation!(population, config)

    @test length(get_individuals(population)) > 0
    @test all(indiv -> length(get_ops(indiv)) > 0, get_individuals(population))
end