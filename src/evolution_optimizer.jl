# TODO (low priority) this would be a great place to use an Enum or, even better, an algebraic data type (ADT)
const HISTORIES = [:manual, :survivor, :random, :child, :drop_m, :gain_m, :swap_m, :ops_m]

"An individual (circuit) in the population we are evolving"
mutable struct Individual
    "How did this individual come to be (conventionally a symbol from the `HISTORIES` list)"
    history::Symbol
    "A vector containing a sequence of quantum operations that make up the individual's circuit"
    ops::Vector{Any}
    "A variety of performance estimates"
    performance::Performance
    "Overall fitness, derived from `performance` and some configuration options given to the performance evaluator/optimizer"
    fitness::Float64
    # TODO (low priority) inner constructor that checks that `history` is a member of HISTORIES (will be unnecessary if we used an enum of an ADT)
end

Individual() = Individual(:manual)
Individual(history::Symbol) = Individual(history, [])
Individual(history::Symbol, ops) = Individual(history, ops, Performance(Float64[], 0.0, 0.0, 0.0, 0.0), 0.0)
