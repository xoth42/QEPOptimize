"perform a mutation"
function mutate end

function mutate(op) # a default no-op "mutation" in case one is not implemented
    return op
end

function mutate(indiv::Individual)
    if length(indiv.ops) == 0
        return indiv
    end
    mutated_ops = mutate.(indiv.ops)
    return Individual(:opsmutate, mutated_ops)
end

function mutate(gate::BellMeasure)
    return rand(BellMeasure, gate.m.sidx) # TODO (low priority) this `rand` is a "pun"; should be changed to use a keyword argument to specify affected qubit, but that is a breaking change in BPGates.jl
end

function mutate(gate::CNOTPerm)
    return rand(CNOTPerm, gate.g.idx1, gate.g.idx2) # TODO (low priority) this `rand` is a "pun"; should be changed to use a keyword argument to specify affected qubit, but that is a breaking change in BPGates.jl
end
