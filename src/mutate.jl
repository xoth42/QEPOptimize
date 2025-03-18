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


is_droppable(::Any) = false
is_droppable(::CNOTPerm) = true
is_droppable(::BellMeasure) = true

"make a new individual with randomly deleted operations"
function drop_op(indiv::Individual)
    # Filter the indices of operations that can be dropped
    drop_indices = [i for (i,op) in pairs(indiv.ops) if is_droppable(op)]

    if  isempty(drop_indices)
        # If there are no droppable operations, return the individual as is
        return copy(indiv)
    else
        # Randomly select and delete one of the droppable operations
        new_ops = copy(indiv.ops)
        deleteat!(new_ops, rand(drop_indices))
        return Individual(:drop, new_ops)
    end
end
