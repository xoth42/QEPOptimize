
"""
Main Tasks:
Configurable.jl 
QEPO.jl
Visualizer.jl


TODO:Reformat code to follow Julia guidelines
See TODOs in Configurable.jl, Optimizer.jl
"""
module QEPO

include("Configurable.jl")
include("Optimizer.jl")
include("Visualizer.jl")
# using QEPO.Configurable
# using QEPO.Optimizer

# function run_optimizer(config::AbstractConfiguration)::Performance
#     ### Run default optimization
#     # Initialize genetic optimizer 
#     advanced_config = get_advanced_config(config)
#     hardware_config = get_hardware_config(config)
#     individuals = []
#     selection_history = Dict()
#     population = Population_hardware(
#         get_raw_bell_pairs(config),
#         get_purified_pairs(config),
#         get_num_registers(config),
#         get_code_distance(advanced_config),
#         get_optimize_for(config),
#         get_communication_fidelity_in(advanced_config),

#     )

# end


end # module QEPO
