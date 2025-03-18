module QEPO

using BPGates
using BPGates: mctrajectory!, continue_stat # TODO these should be exported by default

using Statistics: mean

using Random: randperm

export Individual, calculate_performance!, f_in_to_pauli, NetworkFidelity, NetworkPauliNoise

include("noises.jl") # TODO (low priority) this should be upstreamed to BPGates
include("evolution_optimizer.jl")

end
