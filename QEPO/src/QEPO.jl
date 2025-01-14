module QEPO

"""
Main Tasks:
Configurable.jl 
QEPO.jl


TODO:Reformat code to follow Julia guidelines"""



using Base.Threads #  for multithreading
using Random
using Statistics
using LinearAlgebra
using BPGates # efficient representation of purification circuits
using QuantumClifford # general-purpose tools for Clifford circuits
using QuantumClifford.Experimental.NoisyCircuits
using BPGates: T1NoiseOp, T2NoiseOp
using OhMyThreads: tmap

import Base: sort! # imports necessary for defining new methods of functions defined in Base


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

# abstract type QuantumOperations 
# TODO document this
mutable struct Individual
    history::String
    k::Int                                                                           # number of purified pairs
    r::Int                                                                           # number of registers
    f_in::Float64
    code_distance::Int
    ops::Vector{Union{PauliNoiseBellGate{CNOTPerm},NoisyBellMeasureNoisyReset, PauliNoiseBellGate{BellSwap}}}      # A vector containing a sequence of quantum operations that make up the individual's circuit
    performance::Performance
    fitness::Float64
    optimize_for::CostFunction
end

end # module QEPO
