using Test
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
mock_hardware_config = HardwareConfiguration(mock_calibration_file, Dict())
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