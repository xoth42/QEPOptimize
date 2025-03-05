using Test
# include("../src/Optimizer.jl")
using Revise
using QEPO.Optimizer:generate_noisy_BellSwap_ops_for_individual,generate_valid_pairs, calculate_performance!, push_noise!, Individual
using QEPO.Configurable
using BPGates
using Quantikz
# Basic config
config = Configuration()
# Valid qubits([ 43, 44, 45, 46, 47, 48, 49, 50])
# raw_bell_pairs = 6
# purified_pairs = 1
# num_registers = 4

# Simple version of gates (noise = 0)
# start with entanglement purification bellswaps
simpleGates = [BPGates.PauliNoiseBellGate{BPGates.BellSwap}(BPGates.BellSwap(4, 3), 0.0, 0.0, 0.0),BPGates.PauliNoiseBellGate{BPGates.BellSwap}(BPGates.BellSwap(3, 2), 0.0, 0.0, 0.0),BPGates.PauliNoiseBellGate{BPGates.BellSwap}(BPGates.BellSwap(2, 1), 0.0, 0.0, 0.0), BPGates.PauliNoiseBellGate{BPGates.BellSwap}(BPGates.BellSwap(4, 3), 0.0, 0.0, 0.0), BPGates.PauliNoiseBellGate{BPGates.BellSwap}(BPGates.BellSwap(3, 2), 0.0, 0.0, 0.0), BPGates.PauliNoiseBellGate{BPGates.BellSwap}(BPGates.BellSwap(4, 3), 0.0, 0.0, 0.0),
 
# now 2 cnot perms 
BPGates.CNOTPerm(6, 2, 3, 4), BPGates.CNOTPerm(6, 2, 1, 2),

# now 2 bellmeasures
BPGates.BellMeasure(1, 4),BPGates.BellMeasure(1, 2)]

# displaycircuit(simpleGates)


valid_pairs = generate_valid_pairs(config.hardware_config.valid_qubits)
# Normal list of gates
# start with entanglement purification bellswaps, using the generator here
manGates = generate_noisy_BellSwap_ops_for_individual(config.num_registers,valid_pairs,config.hardware_config.calibration_data)

# add in 2 CNOTPerms 1->2, 3->4, with the same permutation index (for now)
# Ignore noisy gates on these (for now)
manGates = vcat(manGates,[BPGates.CNOTPerm(6, 2, 3, 4),BPGates.CNOTPerm(6, 2, 1, 2)])

# add 2 Bellmeasures in x on register 2 and 4, ignoring noise for now.
manGates = vcat(manGates,[BPGates.BellMeasure(1, 4),BPGates.BellMeasure(1, 2)])

# displaycircuit(manGates)

# Calculate performance
function return_performance(ops)
    indiv = Individual()
    indiv.ops = ops
    calculate_performance!(indiv,
        config.num_simulations,
        config.purified_pairs,
        config.num_registers,
        config.optimize_for,
        config.advanced_config)
    return indiv.performance
end

function log_performance(perf)
    println("Error probabilities: ", perf.error_probabilities)
    println("Purified pairs fidelity: ", perf.purified_pairs_fidelity)
    println("Logical qubit fidelity: ", perf.logical_qubit_fidelity)
    println("Average marginal fidelity: ", perf.average_marginal_fidelity)
end
# Results
manPerf = return_performance(manGates)
println("Manual gates (bellswaps with noise)")
log_performance(manPerf)
# Basic generated gates from short simulation

genGates = [

    BPGates.PauliNoiseBellGate{BPGates.BellSwap}(BPGates.BellSwap(4, 3), 0.008382719809533715, 0.008382719809533715, 0.008382719809533715),BPGates.PauliNoiseBellGate{BPGates.BellSwap}(BPGates.BellSwap(3, 2), 0.008854846264460238, 0.008854846264460238, 0.008854846264460238),BPGates.PauliNoiseBellGate{BPGates.BellSwap}(BPGates.BellSwap(2, 1), 0.009904959580707464, 0.009904959580707464, 0.009904959580707464), BPGates.PauliNoiseBellGate{BPGates.BellSwap}(BPGates.BellSwap(4, 3), 0.008382719809533715, 0.008382719809533715, 0.008382719809533715), BPGates.PauliNoiseBellGate{BPGates.BellSwap}(BPGates.BellSwap(3, 2), 0.008854846264460238, 0.008854846264460238, 0.008854846264460238), BPGates.PauliNoiseBellGate{BPGates.BellSwap}(BPGates.BellSwap(4, 3), 0.008382719809533715, 0.008382719809533715, 0.008382719809533715), 
    
    BPGates.PauliNoiseBellGate{BPGates.CNOTPerm}(BPGates.CNOTPerm(6, 2, 3, 4), 0.0027975865743341466, 0.0027975865743341466, 0.0027975865743341466), BPGates.NoisyBellMeasureNoisyReset(BPGates.BellMeasure(2, 4), 0.02531267, 0.033333333333333326, 0.033333333333333326, 0.033333333333333326), BPGates.NoisyBellMeasureNoisyReset(BPGates.BellMeasure(2, 3), 0.019603000000000002, 0.25, 0.25, 0.25), BPGates.NoisyBellMeasureNoisyReset(BPGates.BellMeasure(2, 3), 0.019603000000000002, 0.25, 0.25, 0.25), BPGates.PauliNoiseBellGate{BPGates.CNOTPerm}(BPGates.CNOTPerm(6, 3, 1, 2), 0.0033067465712277007, 0.0033067465712277007, 0.0033067465712277007), BPGates.NoisyBellMeasureNoisyReset(BPGates.BellMeasure(2, 4), 0.02531267, 0.033333333333333326, 0.033333333333333326, 0.033333333333333326), BPGates.NoisyBellMeasureNoisyReset(BPGates.BellMeasure(1, 3), 0.019603000000000002, 0.25, 0.25, 0.25), BPGates.PauliNoiseBellGate{BPGates.CNOTPerm}(BPGates.CNOTPerm(4, 5, 1, 2), 0.0033067465712277007, 0.0033067465712277007, 0.0033067465712277007), BPGates.NoisyBellMeasureNoisyReset(BPGates.BellMeasure(2, 2), 0.0146464, 0.25, 0.25, 0.25), BPGates.NoisyBellMeasureNoisyReset(BPGates.BellMeasure(2, 2), 0.0146464, 0.25, 0.25, 0.25), BPGates.PauliNoiseBellGate{BPGates.CNOTPerm}(BPGates.CNOTPerm(6, 2, 3, 4), 0.0027975865743341466, 0.0027975865743341466, 0.0027975865743341466)]


# Now trying the genGates

genPerf = return_performance(genGates)
println("Generated gates")
log_performance(genPerf)

# simple gates
simplePerf = return_performance(simpleGates)
println("Simple gates (no noise)")
log_performance(simplePerf)


